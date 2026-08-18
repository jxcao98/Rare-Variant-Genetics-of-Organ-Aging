#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: summarize_saige_results.R MANIFEST.tsv OUTPUT_DIR")
}
manifest <- fread(args[[1]])
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!all(c("trait", "path") %in% names(manifest))) {
  stop("Manifest requires trait and path columns")
}

read_result <- function(trait, path) {
  frame <- fread(path)
  required <- c("Region", "Group", "max_MAF", "MAC", "Pvalue")
  if (!all(required %in% names(frame))) {
    stop(path, " is missing columns: ", paste(setdiff(required, names(frame)), collapse = ", "))
  }
  frame %>% mutate(trait = trait, source_file = basename(path), .before = 1)
}

all_results <- lapply(seq_len(nrow(manifest)), function(index) {
  read_result(manifest$trait[[index]], manifest$path[[index]])
}) %>% bind_rows()

allowed_groups <- c("LoF", "LoF;DMis_Revel75;DMis_Revel50", "LoF:DMis_Revel75:DMis_Revel50")
tested <- all_results %>%
  filter(MAC >= 10, Group %in% allowed_groups) %>%
  mutate(
    mask = if_else(Group == "LoF", "LoF", "LoF+DMis"),
    Pvalue = pmin(pmax(Pvalue, .Machine$double.xmin), 1),
    fdr_all_tests = p.adjust(Pvalue, method = "BH")
  )

beta_candidates <- intersect(c("BETA_Burden", "BETA"), names(tested))
se_candidates <- intersect(c("SE_Burden", "SE"), names(tested))
if (length(beta_candidates) == 0 || length(se_candidates) == 0) {
  tested$burden_beta <- NA_real_
  tested$burden_se <- NA_real_
} else {
  beta_column <- beta_candidates[[1]]
  se_column <- se_candidates[[1]]
  tested <- tested %>%
    mutate(burden_beta = .data[[beta_column]], burden_se = .data[[se_column]])
}

best <- tested %>%
  group_by(trait, Region) %>%
  slice_min(Pvalue, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    exome_wide_significant = Pvalue < 2.5e-6,
    candidate = exome_wide_significant & fdr_all_tests < 0.05
  )

candidate <- best %>% filter(candidate)
pleiotropy <- candidate %>%
  group_by(Region) %>%
  summarise(
    n_traits = n_distinct(trait),
    traits = paste(sort(unique(trait)), collapse = ";"),
    effect_pattern = case_when(
      all(burden_beta > 0, na.rm = TRUE) ~ "positive",
      all(burden_beta < 0, na.rm = TRUE) ~ "negative",
      any(burden_beta > 0, na.rm = TRUE) & any(burden_beta < 0, na.rm = TRUE) ~ "mixed",
      TRUE ~ "unavailable"
    ),
    .groups = "drop"
  ) %>%
  filter(n_traits >= 2) %>%
  arrange(desc(n_traits), Region)

inflation <- tested %>%
  group_by(trait, mask, max_MAF) %>%
  summarise(
    tests = n(),
    lambda_gc = median(qchisq(1 - Pvalue, df = 1), na.rm = TRUE) / qchisq(0.5, df = 1),
    .groups = "drop"
  )

fwrite(tested, file.path(output_dir, "all_gene_tests.tsv"), sep = "\t")
fwrite(best, file.path(output_dir, "best_test_per_gene_trait.tsv"), sep = "\t")
fwrite(candidate, file.path(output_dir, "candidate_gene_trait_associations.tsv"), sep = "\t")
fwrite(pleiotropy, file.path(output_dir, "pleiotropic_genes.tsv"), sep = "\t")
fwrite(inflation, file.path(output_dir, "genomic_inflation.tsv"), sep = "\t")
