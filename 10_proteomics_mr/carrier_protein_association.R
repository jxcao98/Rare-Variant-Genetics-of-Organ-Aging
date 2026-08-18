#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(broom)
  library(data.table)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: carrier_protein_association.R MAPPING.tsv CARRIERS.tsv PROTEOMICS.tsv OUTPUT.tsv")
}

mapping <- fread(args[[1]])
carriers <- fread(args[[2]])
proteomics <- fread(args[[3]])

required_mapping <- c("gene", "protein", "variant_class", "max_maf")
required_carriers <- c("sample_id", "gene", "variant_class", "max_maf", "burden")
required_covariates <- c(
  "sample_id", "Age", "Sex", "FastingTime", "CollectionToAssayInterval",
  "AssessmentCenter", "Platform", paste0("PC", 1:20)
)
if (!all(required_mapping %in% names(mapping))) stop("Mapping is missing required columns")
if (!all(required_carriers %in% names(carriers))) stop("Carrier table is missing required columns")
if (!all(required_covariates %in% names(proteomics))) {
  stop("Proteomics table is missing: ", paste(setdiff(required_covariates, names(proteomics)), collapse = ", "))
}

proteomics <- proteomics %>%
  mutate(
    Age2 = Age^2,
    Sex = factor(Sex),
    AssessmentCenter = factor(AssessmentCenter),
    Platform = factor(Platform)
  )
mapping <- mapping %>% distinct(gene, protein, variant_class, max_maf)
carriers <- carriers %>% mutate(max_maf = as.numeric(max_maf))
mapping <- mapping %>% mutate(max_maf = as.numeric(max_maf))
duplicate_masks <- mapping %>% count(gene, protein) %>% filter(n > 1)
if (nrow(duplicate_masks) > 0) {
  stop("Mapping must select exactly one rare-variant mask for each gene-protein pair")
}

covariates <- c(
  "Age", "Age2", "Sex", "FastingTime", "CollectionToAssayInterval",
  "AssessmentCenter", "Platform", paste0("PC", 1:20)
)

results <- lapply(seq_len(nrow(mapping)), function(index) {
  record <- mapping[index, ]
  protein <- record$protein[[1]]
  if (!protein %in% names(proteomics)) {
    warning("Protein column not found: ", protein)
    return(NULL)
  }
  carrier_subset <- carriers %>%
    filter(
      gene == record$gene[[1]],
      variant_class == record$variant_class[[1]],
      max_maf == record$max_maf[[1]]
    ) %>%
    transmute(sample_id, carrier = as.integer(burden > 0))
  if (nrow(carrier_subset) == 0) return(NULL)

  analysis <- proteomics %>%
    select(all_of(c(required_covariates, "Age2", protein))) %>%
    inner_join(carrier_subset, by = "sample_id") %>%
    mutate(protein_value = .data[[protein]])
  formula <- reformulate(c("carrier", covariates), response = "protein_value")
  fit <- lm(formula, data = analysis)
  coefficient <- broom::tidy(fit) %>% filter(term == "carrier")
  if (nrow(coefficient) == 0) return(NULL)

  coefficient %>%
    transmute(
      gene = record$gene[[1]],
      protein = protein,
      variant_class = record$variant_class[[1]],
      max_maf = record$max_maf[[1]],
      n = nobs(fit),
      carriers = sum(model.frame(fit)$carrier == 1),
      beta = estimate,
      standard_error = std.error,
      p_value = p.value
    )
}) %>% bind_rows() %>%
  mutate(fdr = p.adjust(p_value, method = "BH")) %>%
  arrange(fdr, p_value)

fwrite(results, args[[4]], sep = "\t")
