#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(TwoSampleMR)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: run_two_sample_mr.R MANIFEST.tsv OUTPUT_DIR")
}
manifest <- fread(args[[1]])
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required_manifest <- c(
  "pair_id", "protein", "trait", "exposure_gwas", "outcome_gwas", "instruments",
  "exposure_sample_set", "outcome_sample_set"
)
if (!all(required_manifest %in% names(manifest))) {
  stop("MR manifest is missing: ", paste(setdiff(required_manifest, names(manifest)), collapse = ", "))
}

format_summary <- function(frame, type) {
  required <- c(
    "variant_id", "beta", "se", "effect_allele", "other_allele",
    "eaf", "p_value", "chromosome", "position"
  )
  if (!all(required %in% names(frame))) {
    stop(type, " GWAS is missing: ", paste(setdiff(required, names(frame)), collapse = ", "))
  }
  TwoSampleMR::format_data(
    frame,
    type = type,
    snp_col = "variant_id",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    eaf_col = "eaf",
    pval_col = "p_value",
    chr_col = "chromosome",
    pos_col = "position"
  )
}

for (index in seq_len(nrow(manifest))) {
  record <- manifest[index, ]
  if (record$exposure_sample_set[[1]] == record$outcome_sample_set[[1]]) {
    stop("Exposure and outcome sample-set labels are identical for ", record$pair_id[[1]])
  }

  instruments <- fread(record$instruments[[1]], header = FALSE)$V1
  exposure_raw <- fread(record$exposure_gwas[[1]]) %>%
    filter(variant_id %in% instruments) %>%
    mutate(f_statistic = (beta / se)^2) %>%
    filter(is.finite(f_statistic), f_statistic >= 10)
  if (nrow(exposure_raw) < 3) {
    warning("Fewer than three strong instruments for ", record$pair_id[[1]])
    next
  }
  outcome_raw <- fread(record$outcome_gwas[[1]]) %>%
    filter(variant_id %in% exposure_raw$variant_id)

  exposure <- format_summary(exposure_raw, "exposure")
  exposure$exposure <- record$protein[[1]]
  outcome <- format_summary(outcome_raw, "outcome")
  outcome$outcome <- record$trait[[1]]
  harmonised <- harmonise_data(exposure, outcome, action = 2) %>% filter(mr_keep)
  if (nrow(harmonised) < 3) {
    warning("Fewer than three harmonised instruments for ", record$pair_id[[1]])
    next
  }

  run_dir <- file.path(output_dir, record$pair_id[[1]])
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  methods <- c(
    "mr_ivw",
    "mr_egger_regression",
    "mr_weighted_median",
    "mr_weighted_mode"
  )
  mr_results <- mr(harmonised, method_list = methods)
  heterogeneity <- mr_heterogeneity(harmonised)
  pleiotropy <- mr_pleiotropy_test(harmonised)
  leave_one_out <- mr_leaveoneout(harmonised)

  fwrite(mr_results, file.path(run_dir, "mr_results.tsv"), sep = "\t")
  fwrite(heterogeneity, file.path(run_dir, "heterogeneity.tsv"), sep = "\t")
  fwrite(pleiotropy, file.path(run_dir, "horizontal_pleiotropy.tsv"), sep = "\t")
  fwrite(leave_one_out, file.path(run_dir, "leave_one_out.tsv"), sep = "\t")
  fwrite(harmonised, file.path(run_dir, "harmonised_instruments.tsv"), sep = "\t")
}

