#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(GOSemSim)
  library(org.Hs.eg.db)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3 || length(args) > 4) {
  stop("Usage: go_semantic_similarity.R CANDIDATE_SETS.tsv AGING_ATLAS.tsv OUTPUT.tsv [KEYTYPE]")
}
keytype <- ifelse(length(args) == 4, args[[4]], "ENSEMBL")

candidates <- fread(args[[1]])
aging <- fread(args[[2]])
if (!all(c("gene", "set") %in% names(candidates))) {
  stop("Candidate input requires gene and set columns")
}
if (!all(c("gene", "hallmark") %in% names(aging))) {
  stop("Aging Atlas input requires gene and hallmark columns")
}

candidates <- candidates %>% filter(!is.na(gene), gene != "") %>% distinct(gene, set)
aging <- aging %>% filter(!is.na(gene), gene != "") %>% distinct(gene, hallmark)

sasp <- "senescenceassociated_secretory_phenotype"
nfkb <- "NF_kB_related_gene"
chronic <- aging %>%
  filter(hallmark %in% c(sasp, nfkb)) %>%
  transmute(gene, hallmark = "Chronic_Inflammation")
aging <- aging %>%
  filter(!hallmark %in% c(sasp, nfkb)) %>%
  bind_rows(chronic) %>%
  bind_rows(aging %>% distinct(gene) %>% mutate(hallmark = "Total_Aging")) %>%
  distinct(gene, hallmark)

candidate_sets <- split(candidates$gene, candidates$set) %>% lapply(unique)
aging_sets <- split(aging$gene, aging$hallmark) %>% lapply(unique)
names(candidate_sets) <- paste0("candidate::", names(candidate_sets))
names(aging_sets) <- paste0("aging::", names(aging_sets))

semantic_data <- GOSemSim::godata(
  OrgDb = org.Hs.eg.db,
  keytype = keytype,
  ont = "BP",
  computeIC = FALSE
)
similarity <- GOSemSim::mclusterSim(
  c(candidate_sets, aging_sets),
  semData = semantic_data,
  measure = "Wang",
  combine = "BMA"
)

candidate_names <- names(candidate_sets)
aging_names <- names(aging_sets)
result <- rbindlist(lapply(candidate_names, function(candidate_name) {
  rbindlist(lapply(aging_names, function(aging_name) {
    candidate_genes <- candidate_sets[[candidate_name]]
    aging_genes <- aging_sets[[aging_name]]
    data.frame(
      candidate_set = sub("^candidate::", "", candidate_name),
      aging_set = sub("^aging::", "", aging_name),
      similarity = similarity[candidate_name, aging_name],
      overlapping_genes = length(intersect(candidate_genes, aging_genes)),
      candidate_genes = length(candidate_genes),
      aging_genes = length(aging_genes)
    )
  }))
}))

fwrite(result, args[[3]], sep = "\t")

