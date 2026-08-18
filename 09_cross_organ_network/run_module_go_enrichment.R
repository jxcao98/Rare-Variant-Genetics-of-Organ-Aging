#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(gprofiler2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: run_module_go_enrichment.R CYTOSCAPE_MCL_MODULES.tsv OUTPUT_DIR")
}

modules <- fread(args[[1]])
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

rename_alias <- function(table, target, aliases, label) {
  if (target %in% names(table)) return(invisible(NULL))
  normalized_names <- tolower(gsub("[^a-z0-9]", "", names(table)))
  normalized_aliases <- tolower(gsub("[^a-z0-9]", "", aliases))
  matched <- which(normalized_names %in% normalized_aliases)
  if (length(matched) != 1) {
    stop(label, " requires a column for ", target)
  }
  setnames(table, names(table)[matched], target)
}

rename_alias(modules, "gene", c("gene", "protein name"), "Cytoscape module export")
rename_alias(modules, "module", c("module", "cluster number"), "Cytoscape module export")
modules <- unique(modules[!is.na(gene) & gene != "" & !is.na(module), .(gene, module)])

module_gene_sets <- split(unique(modules$gene), modules$module)
enrichment <- gost(
  query = module_gene_sets,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = TRUE,
  user_threshold = 0.05,
  correction_method = "g_SCS",
  sources = c("GO:BP", "GO:MF")
)

if (is.null(enrichment) || is.null(enrichment$result) || nrow(enrichment$result) == 0) {
  empty <- data.table(
    module = character(),
    source = character(),
    term_id = character(),
    term_name = character(),
    adjusted_p_value = numeric(),
    module_size = integer(),
    overlap_size = integer(),
    overlapping_genes = character()
  )
  fwrite(empty, file.path(output_dir, "module_go_enrichment.tsv"), sep = "\t")
  fwrite(empty, file.path(output_dir, "module_top_go_terms.tsv"), sep = "\t")
  quit(save = "no", status = 0)
}

results <- enrichment$result %>%
  transmute(
    module = as.character(query),
    source,
    term_id,
    term_name,
    adjusted_p_value = p_value,
    module_size = query_size,
    overlap_size = intersection_size,
    overlapping_genes = vapply(
      intersection,
      function(genes) paste(genes, collapse = ","),
      character(1)
    )
  ) %>%
  arrange(as.numeric(module), source, adjusted_p_value)

top_results <- results %>%
  filter(adjusted_p_value < 0.05) %>%
  group_by(module, source) %>%
  slice_min(adjusted_p_value, n = 5, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(as.numeric(module), source, adjusted_p_value)

fwrite(results, file.path(output_dir, "module_go_enrichment.tsv"), sep = "\t")
fwrite(top_results, file.path(output_dir, "module_top_go_terms.tsv"), sep = "\t")
