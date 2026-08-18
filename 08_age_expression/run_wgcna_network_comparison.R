#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(WGCNA)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4 || length(args) > 7) {
  stop(paste(
    "Usage: run_wgcna_network_comparison.R WHOLE_BLOOD_EXPRESSION.tsv[.gz]",
    paste(
      "METADATA.tsv CANDIDATE_SETS.tsv OUTPUT_DIR",
      "[MODULE_ANNOTATIONS.tsv|-] [PERMUTATIONS] [THREADS]"
    )
  ))
}

expression_path <- args[[1]]
metadata_path <- args[[2]]
candidate_path <- args[[3]]
output_dir <- args[[4]]
module_annotation_path <- if (length(args) >= 5 && args[[5]] != "-") args[[5]] else NA_character_
permutations <- if (length(args) >= 6) as.integer(args[[6]]) else 10000L
threads <- if (length(args) == 7) as.integer(args[[7]]) else 1L
if (is.na(permutations) || permutations < 1 || is.na(threads) || threads < 1) {
  stop("PERMUTATIONS and THREADS must be positive integers")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(66)
allowWGCNAThreads(nThreads = threads)

require_columns <- function(table, columns, label) {
  missing <- setdiff(columns, names(table))
  if (length(missing) > 0) {
    stop(label, " is missing columns: ", paste(missing, collapse = ", "))
  }
}

expression_table <- fread(expression_path)
metadata <- fread(metadata_path)
candidates <- fread(candidate_path)
require_columns(metadata, c("sample_id", "age", "sex", "rin", "ischemic_time"), "Metadata")
require_columns(candidates, c("gene", "set"), "Candidate-set table")

genes <- sub("\\..*$", "", as.character(expression_table[[1]]))
expression <- as.matrix(expression_table[, -1, with = FALSE])
storage.mode(expression) <- "double"
rownames(expression) <- genes
candidates[, gene := sub("\\..*$", "", as.character(gene))]
candidates <- unique(candidates[!is.na(gene) & gene != "" & !is.na(set) & set != ""])

common_samples <- intersect(colnames(expression), metadata$sample_id)
metadata <- metadata[match(common_samples, sample_id)]
metadata[, age_group := fifelse(age <= 39, "Young", fifelse(age >= 60, "Old", NA_character_))]
metadata <- metadata[
  !is.na(age_group) & complete.cases(age, sex, rin, ischemic_time)
]
expression <- expression[, metadata$sample_id, drop = FALSE]

group_sizes <- table(metadata$age_group)
if (!all(c("Young", "Old") %in% names(group_sizes)) || min(group_sizes) < 20) {
  stop("At least 20 complete whole-blood samples are required in each age group")
}

gene_mad <- apply(expression, 1, stats::mad, na.rm = TRUE)
mad_cutoff <- unname(stats::quantile(gene_mad, probs = 0.25, na.rm = TRUE))
keep_genes <- is.finite(gene_mad) & gene_mad > 0 & gene_mad >= mad_cutoff
expression <- expression[keep_genes, , drop = FALSE]
if (nrow(expression) < 30) {
  stop("Fewer than 30 genes remain after retaining the top 75% by MAD")
}

select_power <- function(dat_expr, label) {
  powers <- c(1:10, seq(12, 30, by = 2))
  fit <- pickSoftThreshold(
    dat_expr,
    RsquaredCut = 0.85,
    powerVector = powers,
    corFnc = "bicor",
    corOptions = list(use = "p", maxPOutliers = 0.1),
    networkType = "signed",
    verbose = 0
  )
  fit_indices <- as.data.table(fit$fitIndices)
  eligible <- fit_indices[is.finite(SFT.R.sq) & SFT.R.sq >= 0.85]
  if (nrow(eligible) == 0) {
    stop("No soft-thresholding power reached scale-free fit >= 0.85 for ", label)
  }
  list(power = min(eligible$Power), fit_indices = fit_indices)
}

build_network <- function(group_expression, group_metadata, label) {
  if (uniqueN(group_metadata$sex) < 2 || stats::sd(group_metadata$rin) == 0 ||
      stats::sd(group_metadata$ischemic_time) == 0) {
    stop("A required covariate has no variation in the ", label, " group")
  }

  dat_expr_raw <- as.data.frame(t(group_expression[, group_metadata$sample_id, drop = FALSE]))
  rownames(dat_expr_raw) <- group_metadata$sample_id
  covariates <- data.frame(
    sex = factor(group_metadata$sex),
    rin = as.numeric(group_metadata$rin),
    ischemic_time = as.numeric(group_metadata$ischemic_time)
  )
  removed_covariates <- model.matrix(~ sex + rin + ischemic_time, data = covariates)
  removed_covariates <- removed_covariates[, -1, drop = FALSE]
  retained_covariates <- data.frame(age = as.numeric(group_metadata$age))

  adjusted <- empiricalBayesLM(
    data = dat_expr_raw,
    removedCovariates = removed_covariates,
    retainedCovariates = retained_covariates,
    automaticWeights = "bicov",
    aw.maxPOutliers = 0.1,
    verbose = 0
  )$adjustedData

  quality <- goodSamplesGenes(adjusted, verbose = 0)
  adjusted <- adjusted[quality$goodSamples, quality$goodGenes, drop = FALSE]
  if (nrow(adjusted) < 20 || ncol(adjusted) < 30) {
    stop("Insufficient samples or genes after WGCNA quality checks for ", label)
  }

  power_selection <- select_power(adjusted, label)
  network <- blockwiseModules(
    adjusted,
    power = power_selection$power,
    networkType = "signed",
    TOMType = "signed",
    corType = "bicor",
    maxPOutliers = 0.1,
    minModuleSize = 30,
    mergeCutHeight = 0.25,
    numericLabels = FALSE,
    nThreads = threads,
    randomSeed = 66,
    verbose = 0
  )

  module_colors <- as.character(network$colors)
  names(module_colors) <- colnames(adjusted)
  list(
    label = label,
    dat_expr = adjusted,
    module_colors = module_colors,
    module_eigengenes = orderMEs(network$MEs),
    power = power_selection$power,
    soft_threshold_fit = power_selection$fit_indices
  )
}

young_metadata <- metadata[age_group == "Young"]
old_metadata <- metadata[age_group == "Old"]
young_network_path <- file.path(output_dir, "young_wgcna_network.rds")
old_network_path <- file.path(output_dir, "old_wgcna_network.rds")

if (!is.na(module_annotation_path) && file.exists(young_network_path) &&
    file.exists(old_network_path)) {
  message("Loading cached Young and Old networks for the annotated comparison")
  young_network <- readRDS(young_network_path)
  old_network <- readRDS(old_network_path)
} else {
  young_network <- build_network(expression, young_metadata, "Young")
  old_network <- build_network(expression, old_metadata, "Old")
  saveRDS(young_network, young_network_path)
  saveRDS(old_network, old_network_path)
}

network_metrics <- function(network) {
  kme <- signedKME(
    network$dat_expr,
    network$module_eigengenes,
    corFnc = "bicor",
    corOptions = "use = 'p', maxPOutliers = 0.1"
  )
  assigned_kme <- vapply(seq_len(nrow(kme)), function(index) {
    column <- paste0("kME", network$module_colors[[index]])
    if (column %in% colnames(kme)) kme[index, column] else NA_real_
  }, numeric(1))
  data.table(
    gene = colnames(network$dat_expr),
    module = unname(network$module_colors),
    kME = assigned_kme
  )
}

young_metrics <- network_metrics(young_network)
old_metrics <- network_metrics(old_network)

# Independent WGCNA color labels have no cross-network identity. For example,
# Young "blue" and Old "turquoise" may represent the same biological program,
# while two modules both named "blue" may represent different programs. The
# workflow therefore pauses after exporting memberships. Annotate each module
# from its enrichment results, fill biological_module in the generated template,
# and rerun this script with the completed annotation table as argument five.
fwrite(
  young_metrics,
  file.path(output_dir, "young_module_membership.tsv"),
  sep = "\t"
)
fwrite(
  old_metrics,
  file.path(output_dir, "old_module_membership.tsv"),
  sep = "\t"
)

module_annotation_template <- rbindlist(list(
  young_metrics[, .(module_size = .N), by = .(network = "Young", module_color = module)],
  old_metrics[, .(module_size = .N), by = .(network = "Old", module_color = module)]
))
module_annotation_template[, biological_module := ""]
setorder(module_annotation_template, network, -module_size, module_color)
fwrite(
  module_annotation_template,
  file.path(output_dir, "module_annotation_template.tsv"),
  sep = "\t"
)
fwrite(
  rbindlist(list(
    cbind(group = "Young", young_network$soft_threshold_fit),
    cbind(group = "Old", old_network$soft_threshold_fit)
  ), fill = TRUE),
  file.path(output_dir, "soft_threshold_selection.tsv"),
  sep = "\t"
)

if (is.na(module_annotation_path)) {
  message(paste(
    "PAUSED: module memberships and module_annotation_template.tsv were written to",
    output_dir
  ))
  message(paste(
    "Determine Young and Old module functions from enrichment results, fill",
    "biological_module, and rerun with that table as argument five."
  ))
  quit(save = "no", status = 0)
}

module_annotations <- fread(module_annotation_path)
require_columns(
  module_annotations,
  c("network", "module_color", "biological_module"),
  "Module-annotation table"
)
module_annotations[, `:=`(
  network = as.character(network),
  module_color = as.character(module_color),
  biological_module = trimws(as.character(biological_module))
)]
module_annotations <- unique(
  module_annotations[, .(network, module_color, biological_module)]
)
if (module_annotations[, anyDuplicated(paste(network, module_color, sep = "\r"))] > 0) {
  stop("Each network/module_color pair must have exactly one biological_module label")
}

required_annotations <- module_annotation_template[module_color != "grey"]
annotation_check <- merge(
  required_annotations[, .(network, module_color)],
  module_annotations,
  by = c("network", "module_color"),
  all.x = TRUE
)
missing_annotations <- annotation_check[
  is.na(biological_module) | biological_module == ""
]
if (nrow(missing_annotations) > 0) {
  stop(
    "Missing biological module annotations for: ",
    paste(
      paste(missing_annotations$network, missing_annotations$module_color, sep = "/"),
      collapse = ", "
    )
  )
}

annotate_metrics <- function(metrics, network_name) {
  mapping <- module_annotations[
    network == network_name,
    .(module = module_color, biological_module)
  ]
  annotated <- merge(metrics, mapping, by = "module", all.x = TRUE, sort = FALSE)
  annotated[module == "grey", biological_module := NA_character_]
  annotated
}

young_metrics <- annotate_metrics(young_metrics, "Young")
old_metrics <- annotate_metrics(old_metrics, "Old")
gene_metrics <- merge(
  young_metrics,
  old_metrics,
  by = "gene",
  suffixes = c("_young", "_old")
)
gene_metrics[, `:=`(
  delta_kME = kME_old - kME_young,
  absolute_delta_kME = abs(kME_old - kME_young),
  module_changed = fifelse(
    is.na(biological_module_young) | is.na(biological_module_old),
    NA,
    biological_module_young != biological_module_old
  )
)]
setorder(gene_metrics, -absolute_delta_kME)

candidate_metrics <- merge(candidates, gene_metrics, by = "gene")
setorder(candidate_metrics, set, -absolute_delta_kME)

kme_set_tests <- rbindlist(lapply(sort(unique(candidates$set)), function(set_name) {
  finite_metrics <- gene_metrics[is.finite(absolute_delta_kME)]
  set_genes <- intersect(candidates[set == set_name, gene], finite_metrics$gene)
  if (length(set_genes) == 0) return(NULL)
  observed <- mean(finite_metrics[gene %in% set_genes, absolute_delta_kME])
  background <- finite_metrics$absolute_delta_kME
  null <- replicate(
    permutations,
    mean(sample(background, length(set_genes), replace = FALSE))
  )
  data.table(
    set = set_name,
    genes = length(set_genes),
    mean_absolute_delta_kME = observed,
    p_value = (sum(null >= observed) + 1) / (permutations + 1)
  )
}))
kme_set_tests[, fdr := p.adjust(p_value, method = "BH")]
setorder(kme_set_tests, fdr, set)

candidate_pairs <- rbindlist(lapply(sort(unique(candidates$set)), function(set_name) {
  set_genes <- sort(intersect(candidates[set == set_name, gene], gene_metrics$gene))
  if (length(set_genes) < 2) return(NULL)
  pairs <- as.data.table(t(combn(set_genes, 2)))
  setnames(pairs, c("gene1", "gene2"))
  pairs[, set := set_name]
  setcolorder(pairs, c("set", "gene1", "gene2"))
  pairs
}))
if (nrow(candidate_pairs) == 0) {
  stop("No candidate gene set contains at least two genes in both networks")
}

sample_random_pairs <- function(available_genes, number) {
  first <- sample.int(length(available_genes), number, replace = TRUE)
  second <- sample.int(length(available_genes), number, replace = TRUE)
  identical_index <- first == second
  while (any(identical_index)) {
    second[identical_index] <- sample.int(
      length(available_genes), sum(identical_index), replace = TRUE
    )
    identical_index <- first == second
  }
  low <- pmin(first, second)
  high <- pmax(first, second)
  data.table(gene1 = available_genes[low], gene2 = available_genes[high])
}

common_network_genes <- intersect(
  colnames(young_network$dat_expr),
  colnames(old_network$dat_expr)
)
null_pairs <- sample_random_pairs(common_network_genes, permutations)
unique_candidate_pairs <- unique(candidate_pairs[, .(gene1, gene2)])
needed_pairs <- unique(rbind(unique_candidate_pairs, null_pairs))
needed_pairs[, pair_id := paste(gene1, gene2, sep = "\r")]

adjacency_summary <- function(network, common_genes, pairs) {
  adjacency_matrix <- adjacency(
    network$dat_expr[, common_genes, drop = FALSE],
    type = "signed",
    power = network$power,
    corFnc = "bicor",
    corOptions = list(use = "p", maxPOutliers = 0.1)
  )

  total <- 0
  total_squares <- 0
  count <- 0
  if (nrow(adjacency_matrix) >= 2) {
    for (column in 2:ncol(adjacency_matrix)) {
      values <- adjacency_matrix[seq_len(column - 1), column]
      finite <- values[is.finite(values)]
      total <- total + sum(finite)
      total_squares <- total_squares + sum(finite * finite)
      count <- count + length(finite)
    }
  }
  edge_mean <- total / count
  edge_sd <- sqrt((total_squares - (total * total / count)) / (count - 1))
  if (!is.finite(edge_sd) || edge_sd == 0) stop("Adjacency variance is zero")

  row_index <- match(pairs$gene1, common_genes)
  column_index <- match(pairs$gene2, common_genes)
  selected <- adjacency_matrix[cbind(row_index, column_index)]
  names(selected) <- pairs$pair_id
  rm(adjacency_matrix)
  invisible(gc())
  list(values = selected, mean = edge_mean, sd = edge_sd)
}

young_adjacency <- adjacency_summary(young_network, common_network_genes, needed_pairs)
old_adjacency <- adjacency_summary(old_network, common_network_genes, needed_pairs)
needed_pairs[, `:=`(
  adjacency_young = unname(young_adjacency$values[pair_id]),
  adjacency_old = unname(old_adjacency$values[pair_id])
)]
needed_pairs[, `:=`(
  z_young = (adjacency_young - young_adjacency$mean) / young_adjacency$sd,
  z_old = (adjacency_old - old_adjacency$mean) / old_adjacency$sd
)]
needed_pairs[, z_difference := z_old - z_young]

null_ids <- paste(null_pairs$gene1, null_pairs$gene2, sep = "\r")
null_distribution <- needed_pairs$z_difference[match(null_ids, needed_pairs$pair_id)]
pair_statistics <- merge(unique_candidate_pairs, needed_pairs, by = c("gene1", "gene2"))
pair_statistics[, p_value := vapply(z_difference, function(observed) {
  (sum(abs(null_distribution) >= abs(observed), na.rm = TRUE) + 1) /
    (sum(is.finite(null_distribution)) + 1)
}, numeric(1))]
pair_statistics[, fdr := p.adjust(p_value, method = "BH")]

edge_results <- merge(candidate_pairs, pair_statistics, by = c("gene1", "gene2"))
edge_results[, absolute_z_difference := abs(z_difference)]
setorder(edge_results, fdr, set, -absolute_z_difference)

fwrite(gene_metrics, file.path(output_dir, "all_gene_network_metrics.tsv"), sep = "\t")
fwrite(candidate_metrics, file.path(output_dir, "candidate_gene_network_metrics.tsv"), sep = "\t")
fwrite(kme_set_tests, file.path(output_dir, "candidate_set_kme_permutation.tsv"), sep = "\t")
fwrite(edge_results, file.path(output_dir, "candidate_edge_rewiring.tsv"), sep = "\t")
top_candidate_genes <- gene_metrics[gene %in% unique(candidates$gene)]
setorder(top_candidate_genes, -absolute_delta_kME)
fwrite(
  head(top_candidate_genes, 20),
  file.path(output_dir, "top20_candidate_kme_changes.tsv"),
  sep = "\t"
)
