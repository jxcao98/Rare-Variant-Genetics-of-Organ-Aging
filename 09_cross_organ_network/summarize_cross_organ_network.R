#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop(paste(
    "Usage: summarize_cross_organ_network.R CYTOSCAPE_MCL_MODULES.tsv STRING_EXPORT.tsv",
    "CANDIDATE_SETS.tsv MODULE_GO.tsv CHEA3.tsv OUTPUT_DIR"
  ))
}

modules <- fread(args[[1]])
string_edges <- fread(args[[2]])
candidates <- fread(args[[3]])
go_results <- fread(args[[4]])
chea3 <- fread(args[[5]])
output_dir <- args[[6]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

require_columns <- function(table, columns, label) {
  missing <- setdiff(columns, names(table))
  if (length(missing) > 0) {
    stop(label, " is missing columns: ", paste(missing, collapse = ", "))
  }
}

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
rename_alias(string_edges, "source", c("source", "#node1", "node1"), "STRING export")
rename_alias(string_edges, "target", c("target", "node2"), "STRING export")
rename_alias(string_edges, "score", c("score", "combined_score"), "STRING export")
rename_alias(go_results, "module", c("module", "module_name", "query"), "GO results")
rename_alias(
  go_results,
  "adjusted_p_value",
  c("adjusted_p_value", "p_value"),
  "GO results"
)
rename_alias(chea3, "module", c("module", "cluster_id"), "ChEA3 results")
rename_alias(chea3, "MeanRank", c("MeanRank", "mean rank"), "ChEA3 results")
rename_alias(
  chea3,
  "overlapping_genes",
  c("overlapping_genes", "Overlapping_Genes"),
  "ChEA3 results"
)

modules <- unique(modules[!is.na(gene) & gene != "" & !is.na(module), .(gene, module)])
string_edges <- unique(
  string_edges[
    !is.na(source) & source != "" & !is.na(target) & target != "" & !is.na(score),
    .(source, target, score)
  ]
)
require_columns(modules, c("gene", "module"), "Cytoscape module export")
require_columns(string_edges, c("source", "target", "score"), "STRING export")
require_columns(candidates, c("gene", "set"), "Candidate-set table")
require_columns(
  go_results,
  c("module", "source", "term_id", "term_name", "adjusted_p_value"),
  "GO results"
)
require_columns(
  chea3,
  c("module", "TF", "MeanRank", "overlapping_genes"),
  "ChEA3 results"
)

modules[, module := as.character(module)]
go_results[, module := as.character(module)]
chea3[, `:=`(module = as.character(module), MeanRank = as.numeric(MeanRank))]
candidates <- unique(candidates[!is.na(gene) & gene != "" & !is.na(set) & set != ""])

pleiotropic <- candidates %>%
  distinct(gene, set) %>%
  count(gene, name = "trait_count") %>%
  mutate(is_pleiotropic = trait_count >= 2)
candidate_genes <- unique(candidates$gene)
pleiotropic_genes <- pleiotropic$gene[pleiotropic$is_pleiotropic]

module_membership <- modules %>%
  inner_join(candidates, by = "gene") %>%
  distinct(module, gene, set)
module_summary <- module_membership %>%
  group_by(module) %>%
  summarise(
    genes = n_distinct(gene),
    represented_traits = n_distinct(set),
    traits = paste(sort(unique(set)), collapse = ";"),
    .groups = "drop"
  ) %>%
  arrange(as.numeric(module))

pleiotropy_tests <- bind_rows(lapply(sort(unique(modules$module)), function(module_id) {
  in_module <- intersect(modules$gene[modules$module == module_id], candidate_genes)
  out_module <- setdiff(candidate_genes, in_module)
  contingency <- matrix(
    c(
      length(intersect(in_module, pleiotropic_genes)),
      length(setdiff(in_module, pleiotropic_genes)),
      length(intersect(out_module, pleiotropic_genes)),
      length(setdiff(out_module, pleiotropic_genes))
    ),
    nrow = 2,
    byrow = TRUE
  )
  test <- fisher.test(contingency, alternative = "greater")
  data.frame(
    module = module_id,
    module_genes = length(in_module),
    module_pleiotropic_genes = length(intersect(in_module, pleiotropic_genes)),
    odds_ratio = if (length(test$estimate) == 0) NA_real_ else unname(test$estimate),
    p_value = test$p.value
  )
})) %>%
  mutate(fdr = p.adjust(p_value, method = "BH")) %>%
  arrange(fdr, as.numeric(module))

top_go <- go_results %>%
  filter(adjusted_p_value < 0.05, source %in% c("GO:BP", "GO:MF")) %>%
  group_by(module, source) %>%
  slice_min(adjusted_p_value, n = 5, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(as.numeric(module), source, adjusted_p_value)

top_tf <- chea3 %>%
  filter(is.finite(MeanRank)) %>%
  group_by(module) %>%
  slice_min(MeanRank, n = 5, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(as.numeric(module), MeanRank)

gene_sets <- candidates %>%
  group_by(gene) %>%
  summarise(sets = paste(sort(unique(set)), collapse = ";"), .groups = "drop")
gene_nodes <- modules %>%
  left_join(gene_sets, by = "gene") %>%
  transmute(
    node_id = gene,
    label = gene,
    node_type = "Gene",
    module,
    sets = ifelse(is.na(sets), "", sets)
  )
tf_nodes <- top_tf %>%
  distinct(module, TF) %>%
  transmute(
    node_id = paste("TF", module, TF, sep = "::"),
    label = TF,
    node_type = "TF",
    module,
    sets = ""
  )
cytoscape_nodes <- bind_rows(gene_nodes, tf_nodes)

module_genes <- unique(modules$gene)
ppi_edges <- string_edges %>%
  filter(source %in% module_genes, target %in% module_genes) %>%
  transmute(
    source,
    target,
    edge_type = "PPI",
    score = as.numeric(score),
    mean_rank = NA_real_
  )

regulatory_edges <- bind_rows(lapply(seq_len(nrow(top_tf)), function(index) {
  row <- top_tf[index, ]
  targets <- trimws(unlist(strsplit(as.character(row$overlapping_genes), "[,;]")))
  targets <- intersect(targets[targets != ""], modules$gene[modules$module == row$module])
  if (length(targets) == 0) return(NULL)
  data.frame(
    source = paste("TF", row$module, row$TF, sep = "::"),
    target = targets,
    edge_type = "Regulation",
    score = NA_real_,
    mean_rank = row$MeanRank
  )
}))
cytoscape_edges <- bind_rows(ppi_edges, regulatory_edges)

fwrite(module_summary, file.path(output_dir, "module_summary.tsv"), sep = "\t")
fwrite(pleiotropy_tests, file.path(output_dir, "module_pleiotropy_enrichment.tsv"), sep = "\t")
fwrite(top_go, file.path(output_dir, "module_top_go_terms.tsv"), sep = "\t")
fwrite(top_tf, file.path(output_dir, "module_top_transcription_factors.tsv"), sep = "\t")
fwrite(cytoscape_nodes, file.path(output_dir, "cytoscape_nodes.tsv"), sep = "\t")
fwrite(cytoscape_edges, file.path(output_dir, "cytoscape_edges.tsv"), sep = "\t")
