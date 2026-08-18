#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: summarize_phewas.R PHEWAS.tsv CANDIDATE_GENE_TRAIT.tsv OUTPUT_DIR")
}

phewas <- fread(args[[1]])
candidates <- fread(args[[2]])
output_dir <- args[[3]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required_phewas <- c(
  "GeneName", "Ancestry", "Collapsing model", "Phenotype",
  "Phenotypic category", "P value"
)
if (!all(required_phewas %in% names(phewas))) {
  stop("PheWAS input is missing: ", paste(setdiff(required_phewas, names(phewas)), collapse = ", "))
}
if (!all(c("trait", "gene") %in% names(candidates))) {
  stop("Candidate input requires trait and gene columns")
}

filtered <- phewas %>%
  transmute(
    gene = GeneName,
    ancestry = Ancestry,
    model = `Collapsing model`,
    phenotype = Phenotype,
    disease_category = `Phenotypic category`,
    p_value = `P value`
  ) %>%
  filter(
    ancestry == "European",
    model %in% c("ptv", "ptvraredmg"),
    str_starts(phenotype, "Union#"),
    !str_detect(phenotype, "Block"),
    p_value < 2.5e-6
  ) %>%
  group_by(gene, phenotype) %>%
  slice_min(p_value, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    disease_category = str_squish(
      str_remove(disease_category, regex("^Chapter\\s+[IVXLCDM]+\\s*", ignore_case = TRUE))
    )
  )

candidate_diseases <- candidates %>%
  distinct(trait, gene) %>%
  inner_join(filtered, by = "gene")

summary <- candidate_diseases %>%
  group_by(trait, disease_category) %>%
  summarise(
    associations = n_distinct(paste(gene, phenotype, sep = "::")),
    genes = n_distinct(gene),
    .groups = "drop"
  ) %>%
  arrange(trait, desc(associations), disease_category)

fwrite(candidate_diseases, file.path(output_dir, "candidate_gene_disease_associations.tsv"), sep = "\t")
fwrite(summary, file.path(output_dir, "organ_disease_category_counts.tsv"), sep = "\t")

