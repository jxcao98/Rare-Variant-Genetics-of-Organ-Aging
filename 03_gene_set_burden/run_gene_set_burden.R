#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(broom)
  library(data.table)
  library(dplyr)
})

parse_options <- function(args) {
  result <- list()
  index <- 1
  while (index <= length(args)) {
    key <- sub("^--", "", args[[index]])
    result[[key]] <- args[[index + 1]]
    index <- index + 2
  }
  result
}

options <- parse_options(commandArgs(trailingOnly = TRUE))
required <- c("phenotypes", "burdens", "trait-manifest", "output")
missing <- required[!required %in% names(options)]
if (length(missing) > 0) stop("Missing options: ", paste(missing, collapse = ", "))

phenotypes <- fread(options$phenotypes)
burdens <- fread(options$burdens)
traits <- fread(options[["trait-manifest"]])

required_burden <- c("sample_id", "gene_set", "variant_class", "burden")
if (!all(required_burden %in% names(burdens))) {
  stop("Burden input requires: ", paste(required_burden, collapse = ", "))
}
if (!all(c("trait", "raw_bag", "int_bag") %in% names(traits))) {
  stop("Trait manifest requires trait, raw_bag, and int_bag columns")
}
if (!"sample_id" %in% names(phenotypes)) stop("Phenotypes require a sample_id column")

allowed_gene_sets <- c("Exome", "pLI", "AgingAtlas")
allowed_classes <- c("LoF", "DMis", "Synonymous")
burdens <- burdens %>%
  filter(gene_set %in% allowed_gene_sets, variant_class %in% allowed_classes)

covariates <- c("Sex", "Age", "Age2", "AssessmentCenter", "TDI", paste0("PC", 1:20))
missing_covariates <- setdiff(covariates, names(phenotypes))
if (length(missing_covariates) > 0) {
  stop("Phenotype input is missing covariates: ", paste(missing_covariates, collapse = ", "))
}
phenotypes <- phenotypes %>%
  mutate(Sex = factor(Sex), AssessmentCenter = factor(AssessmentCenter))

fit_burden <- function(data, outcome) {
  model_terms <- c("burden", covariates)
  formula <- reformulate(model_terms, response = outcome)
  fit <- lm(formula, data = data)
  broom::tidy(fit) %>% filter(term == "burden") %>% mutate(n = nobs(fit))
}

results <- list()
result_index <- 0
for (trait_index in seq_len(nrow(traits))) {
  trait <- traits$trait[[trait_index]]
  outcomes <- c(raw = traits$raw_bag[[trait_index]], inverse_normal = traits$int_bag[[trait_index]])
  missing_outcomes <- setdiff(outcomes, names(phenotypes))
  if (length(missing_outcomes) > 0) {
    stop("Missing BAG columns for ", trait, ": ", paste(missing_outcomes, collapse = ", "))
  }

  for (gene_set in allowed_gene_sets) {
    for (variant_class in allowed_classes) {
      burden_subset <- burdens %>%
        filter(.data$gene_set == gene_set, .data$variant_class == variant_class) %>%
        select(sample_id, burden)
      if (nrow(burden_subset) == 0) next
      analysis <- inner_join(phenotypes, burden_subset, by = "sample_id")

      for (scale_name in names(outcomes)) {
        result_index <- result_index + 1
        results[[result_index]] <- fit_burden(analysis, outcomes[[scale_name]]) %>%
          mutate(
            trait = trait,
            bag_scale = scale_name,
            gene_set = gene_set,
            variant_class = variant_class,
            .before = 1
          )
      }
    }
  }
}

final <- bind_rows(results) %>%
  group_by(bag_scale) %>%
  mutate(fdr = p.adjust(p.value, method = "BH")) %>%
  ungroup() %>%
  rename(beta = estimate, standard_error = std.error, p_value = p.value)

fwrite(final, options$output, sep = "\t")
