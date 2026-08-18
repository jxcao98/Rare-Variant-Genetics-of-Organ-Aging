#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(GSVA)
  library(limma)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(paste(
    "Usage: run_gsva_limma.R TISSUE_MANIFEST.tsv",
    "CANDIDATE_SETS.tsv OUTPUT_DIR"
  ))
}

manifest_path <- normalizePath(args[[1]], mustWork = TRUE)
manifest_root <- dirname(manifest_path)
candidate_path <- args[[2]]
output_dir <- args[[3]]
scores_dir <- file.path(output_dir, "scores")
dir.create(scores_dir, recursive = TRUE, showWarnings = FALSE)

resolve_manifest_path <- function(path) {
  if (grepl("^([A-Za-z]:[\\\\/]|/)", path)) path else file.path(manifest_root, path)
}

require_columns <- function(table, columns, label) {
  missing <- setdiff(columns, names(table))
  if (length(missing) > 0) {
    stop(label, " is missing columns: ", paste(missing, collapse = ", "))
  }
}

manifest <- fread(manifest_path)
candidates <- fread(candidate_path)
require_columns(
  manifest,
  c("tissue", "tissue_key", "expression_file", "metadata_file"),
  "Tissue manifest"
)
require_columns(candidates, c("gene", "set"), "Candidate-set table")

candidates[, gene := sub("\\..*$", "", as.character(gene))]
candidates <- unique(candidates[!is.na(gene) & gene != "" & !is.na(set) & set != ""])
gene_sets <- split(candidates$gene, candidates$set)

analyze_tissue <- function(row) {
  tissue <- row$tissue[[1]]
  tissue_key <- row$tissue_key[[1]]
  expression_path <- resolve_manifest_path(row$expression_file[[1]])
  metadata_path <- resolve_manifest_path(row$metadata_file[[1]])

  expression_table <- fread(expression_path)
  metadata <- fread(metadata_path)
  require_columns(metadata, c("sample_id", "age", "sex", "rin", "ischemic_time"), tissue)
  if (ncol(expression_table) < 3) {
    message("Skipping ", tissue, ": expression matrix has fewer than two samples")
    return(NULL)
  }

  genes <- sub("\\..*$", "", as.character(expression_table[[1]]))
  expression <- as.matrix(expression_table[, -1, with = FALSE])
  storage.mode(expression) <- "double"
  rownames(expression) <- genes

  common_samples <- intersect(colnames(expression), metadata$sample_id)
  expression <- expression[, common_samples, drop = FALSE]
  metadata <- metadata[match(common_samples, sample_id)]
  metadata[, age_group := fifelse(age <= 39, "Young", fifelse(age >= 60, "Old", NA_character_))]
  metadata <- metadata[
    !is.na(age_group) & complete.cases(sex, rin, ischemic_time)
  ]

  group_sizes <- table(metadata$age_group)
  if (!all(c("Young", "Old") %in% names(group_sizes)) || min(group_sizes) < 5) {
    message("Skipping ", tissue, ": fewer than five complete samples in an age group")
    return(NULL)
  }
  if (uniqueN(metadata$sex) < 2 || stats::sd(metadata$rin) == 0 ||
      stats::sd(metadata$ischemic_time) == 0) {
    message("Skipping ", tissue, ": a required covariate has no variation")
    return(NULL)
  }

  present_sets <- lapply(gene_sets, intersect, y = rownames(expression))
  present_sets <- present_sets[lengths(present_sets) >= 2]
  if (length(present_sets) == 0) {
    message("Skipping ", tissue, ": no candidate gene set has at least two expressed genes")
    return(NULL)
  }

  gsva_parameters <- gsvaParam(
    exprData = expression,
    geneSets = present_sets,
    kcdf = "Gaussian",
    maxDiff = TRUE
  )
  scores <- gsva(gsva_parameters, verbose = FALSE)
  scores_table <- as.data.table(scores, keep.rownames = "gene_set")
  fwrite(
    scores_table,
    file.path(scores_dir, paste0(tissue_key, ".gsva.tsv.gz")),
    sep = "\t"
  )

  analysis_samples <- metadata$sample_id
  analysis_scores <- scores[, analysis_samples, drop = FALSE]
  model_data <- data.frame(
    age_group = factor(metadata$age_group, levels = c("Young", "Old")),
    sex = factor(metadata$sex),
    rin = as.numeric(metadata$rin),
    ischemic_time = as.numeric(metadata$ischemic_time)
  )
  design <- model.matrix(~ 0 + age_group + sex + rin + ischemic_time, data = model_data)
  if (!is.fullrank(design)) {
    stop("Design matrix is not full rank for tissue: ", tissue)
  }

  contrast <- makeContrasts(
    old_vs_young = age_groupOld - age_groupYoung,
    levels = design
  )
  fit <- lmFit(analysis_scores, design)
  fit <- contrasts.fit(fit, contrast)
  fit <- eBayes(fit)
  results <- as.data.table(topTable(fit, coef = "old_vs_young", number = Inf),
                           keep.rownames = "gene_set")
  results[, `:=`(
    tissue = tissue,
    tissue_key = tissue_key,
    comparison = "young_le39_vs_old_ge60",
    n_young = unname(group_sizes[["Young"]]),
    n_old = unname(group_sizes[["Old"]])
  )]
  setcolorder(
    results,
    c("tissue", "tissue_key", "gene_set", "comparison", "n_young", "n_old")
  )
  results
}

results <- rbindlist(
  lapply(seq_len(nrow(manifest)), function(index) analyze_tissue(manifest[index])),
  fill = TRUE
)
if (nrow(results) == 0) {
  stop("No tissue produced a valid GSVA/limma comparison")
}

results[, fdr := p.adjust(P.Value, method = "BH")]
setorder(results, fdr, tissue, gene_set)
fwrite(results, file.path(output_dir, "gsva_limma_results.tsv"), sep = "\t")
