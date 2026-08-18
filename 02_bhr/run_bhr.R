#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(bhr)
  library(data.table)
  library(dplyr)
})

parse_options <- function(args) {
  result <- list()
  index <- 1
  while (index <= length(args)) {
    key <- sub("^--", "", args[[index]])
    if (index == length(args) || startsWith(args[[index + 1]], "--")) {
      result[[key]] <- TRUE
      index <- index + 1
    } else {
      result[[key]] <- args[[index + 1]]
      index <- index + 2
    }
  }
  result
}

options <- parse_options(commandArgs(trailingOnly = TRUE))
required_options <- c("manifest", "baseline", "output-dir")
missing_options <- required_options[!required_options %in% names(options)]
if (length(missing_options) > 0) {
  stop("Missing options: ", paste(missing_options, collapse = ", "))
}

output_dir <- options[["output-dir"]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

manifest <- fread(options$manifest)
required_manifest <- c("trait", "variant_class", "path")
if (!all(required_manifest %in% names(manifest))) {
  stop("The manifest must contain: ", paste(required_manifest, collapse = ", "))
}
baseline <- fread(options$baseline)

maf_bins <- list(
  ultra_rare = c(0, 1e-5),
  rare = c(1e-5, 1e-4)
)

read_sumstats <- function(path) {
  frame <- fread(path)
  required <- c("gene", "chromosome", "gene_position", "N", "beta", "se", "AF", "phenotype_key")
  if (!all(required %in% names(frame))) {
    stop(path, " is missing BHR columns: ", paste(setdiff(required, names(frame)), collapse = ", "))
  }
  frame
}

flatten_component <- function(component) {
  values <- unlist(component, recursive = TRUE, use.names = TRUE)
  data.frame(metric = names(values), value = as.numeric(values), stringsAsFactors = FALSE)
}

safe_bhr <- function(...) {
  tryCatch(BHR(...), error = function(error) {
    warning(conditionMessage(error))
    NULL
  })
}

heritability_results <- list()
aggregate_results <- list()
result_index <- 0

for (row_index in seq_len(nrow(manifest))) {
  trait <- manifest$trait[[row_index]]
  variant_class <- manifest$variant_class[[row_index]]
  sumstats <- read_sumstats(manifest$path[[row_index]])
  nonempty_bins <- list()

  for (bin_name in names(maf_bins)) {
    bounds <- maf_bins[[bin_name]]
    selected <- sumstats %>% filter(AF >= bounds[[1]], AF < bounds[[2]])
    if (nrow(selected) == 0) next
    nonempty_bins[[bin_name]] <- selected

    fit <- safe_bhr(
      trait1_sumstats = selected,
      mode = "univariate",
      annotations = list(baseline)
    )
    if (is.null(fit)) next

    result_index <- result_index + 1
    heritability_results[[result_index]] <- flatten_component(fit$mixed_model$heritabilities) %>%
      mutate(trait = trait, variant_class = variant_class, maf_bin = bin_name, .before = 1)
  }

  if (length(nonempty_bins) >= 2) {
    fit <- safe_bhr(
      ss_list = nonempty_bins,
      trait_list = list(trait),
      mode = "aggregate",
      annotations = list(baseline)
    )
    if (!is.null(fit)) {
      aggregate_results[[length(aggregate_results) + 1]] <- flatten_component(fit) %>%
        mutate(trait = trait, variant_class = variant_class, .before = 1)
    }
  }
}

if (length(heritability_results) > 0) {
  fwrite(bind_rows(heritability_results), file.path(output_dir, "heritability_by_maf.tsv"), sep = "\t")
}
if (length(aggregate_results) > 0) {
  fwrite(bind_rows(aggregate_results), file.path(output_dir, "heritability_aggregate.tsv"), sep = "\t")
}

if (!is.null(options[["tissue-annotations"]])) {
  tissues <- fread(options[["tissue-annotations"]])
  if (!"gene" %in% names(tissues)) stop("Tissue annotations require a gene column")
  tissue_names <- setdiff(names(tissues), "gene")
  enrichment_results <- list()

  lof_manifest <- manifest %>% filter(variant_class == "LoF")
  for (row_index in seq_len(nrow(lof_manifest))) {
    trait <- lof_manifest$trait[[row_index]]
    sumstats <- read_sumstats(lof_manifest$path[[row_index]]) %>% filter(AF < 1e-5)
    for (tissue in tissue_names) {
      annotation <- tissues %>% select(gene, all_of(tissue))
      fit <- safe_bhr(
        trait1_sumstats = sumstats,
        mode = "univariate",
        annotations = list(baseline, annotation)
      )
      if (is.null(fit)) next
      enrichment_results[[length(enrichment_results) + 1]] <-
        flatten_component(fit$mixed_model$enrichments) %>%
        mutate(trait = trait, tissue = tissue, .before = 1)
    }
  }
  if (length(enrichment_results) > 0) {
    fwrite(bind_rows(enrichment_results), file.path(output_dir, "tissue_enrichment.tsv"), sep = "\t")
  }
}

lof_manifest <- manifest %>% filter(variant_class == "LoF")
trait_names <- unique(lof_manifest$trait)
if (length(trait_names) >= 2) {
  pairs <- combn(trait_names, 2, simplify = FALSE)
  correlation_results <- lapply(pairs, function(pair) {
    path1 <- lof_manifest %>% filter(trait == pair[[1]]) %>% slice(1) %>% pull(path)
    path2 <- lof_manifest %>% filter(trait == pair[[2]]) %>% slice(1) %>% pull(path)
    ss1 <- read_sumstats(path1) %>% filter(AF < 1e-5)
    ss2 <- read_sumstats(path2) %>% filter(AF < 1e-5)
    fit <- safe_bhr(
      trait1_sumstats = ss1,
      trait2_sumstats = ss2,
      mode = "bivariate",
      annotations = list(baseline)
    )
    if (is.null(fit) || is.null(fit$rg)) return(NULL)
    data.frame(
      trait1 = pair[[1]],
      trait2 = pair[[2]],
      rg = fit$rg$rg_mixed,
      rg_se = fit$rg$rg_mixed_se
    )
  }) %>% bind_rows()
  if (nrow(correlation_results) > 0) {
    correlation_results <- correlation_results %>%
      mutate(p_value = 2 * pnorm(-abs(rg / rg_se)), fdr = p.adjust(p_value, method = "BH"))
    fwrite(correlation_results, file.path(output_dir, "genetic_correlations.tsv"), sep = "\t")
  }
}

if (!is.null(options[["bag-matrix"]])) {
  bag <- fread(options[["bag-matrix"]])
  traits <- intersect(trait_names, names(bag))
  if (length(traits) >= 2) {
    pairs <- combn(traits, 2, simplify = FALSE)
    phenotypic <- lapply(pairs, function(pair) {
      complete <- complete.cases(bag[[pair[[1]]]], bag[[pair[[2]]]])
      test <- cor.test(bag[[pair[[1]]]][complete], bag[[pair[[2]]]][complete], method = "pearson")
      data.frame(
        trait1 = pair[[1]], trait2 = pair[[2]], n = sum(complete),
        correlation = unname(test$estimate), p_value = test$p.value
      )
    }) %>% bind_rows() %>%
      mutate(fdr = p.adjust(p_value, method = "BH"))
    fwrite(phenotypic, file.path(output_dir, "phenotypic_correlations.tsv"), sep = "\t")
  }
}
