#!/usr/bin/env python3
"""Test candidate-gene enrichment across Tabula Sapiens broad cell classes."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", required=True, type=Path)
    parser.add_argument("--candidate-sets", required=True, type=Path)
    parser.add_argument("--protein-coding-genes", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--group-column", default="broad_cell_class")
    parser.add_argument("--gene-id-column")
    parser.add_argument("--protein-coding-column", default="gene")
    parser.add_argument("--min-cells-per-class", type=int, default=2)
    parser.add_argument("--expected-classes", type=int, default=40)
    parser.add_argument("--log2fc-threshold", type=float, default=0.58)
    parser.add_argument("--fdr-threshold", type=float, default=0.05)
    parser.add_argument("--strip-gene-version", action="store_true")
    return parser.parse_args()


def require_columns(frame: pd.DataFrame, columns: set[str], label: str) -> None:
    missing = columns.difference(frame.columns)
    if missing:
        raise ValueError(f"{label} is missing columns: {', '.join(sorted(missing))}")


def normalize_gene_ids(values: pd.Series, strip_version: bool) -> pd.Series:
    normalized = values.astype("string").str.strip()
    if strip_version:
        normalized = normalized.str.replace(r"\.[0-9]+$", "", regex=True)
    return normalized


def load_candidate_sets(path: Path, strip_version: bool) -> dict[str, set[str]]:
    frame = pd.read_csv(path, sep="\t", dtype="string")
    require_columns(frame, {"gene", "set"}, "Candidate-set input")
    frame = frame[["gene", "set"]].dropna()
    frame["gene"] = normalize_gene_ids(frame["gene"], strip_version)
    frame["set"] = frame["set"].str.strip()
    frame = frame.loc[frame["gene"].ne("") & frame["set"].ne("")].drop_duplicates()
    if frame.empty:
        raise ValueError("Candidate-set input contains no valid gene-set assignments")
    return {
        name: set(group["gene"])
        for name, group in frame.groupby("set", sort=True)
    }


def load_protein_coding_genes(
    path: Path, column: str, strip_version: bool
) -> set[str]:
    frame = pd.read_csv(path, sep="\t", dtype="string")
    require_columns(frame, {column}, "Protein-coding gene input")
    genes = normalize_gene_ids(frame[column].dropna(), strip_version)
    result = set(genes.loc[genes.ne("")])
    if not result:
        raise ValueError("Protein-coding gene input contains no valid genes")
    return result


def benjamini_hochberg(p_values: pd.Series) -> np.ndarray:
    values = p_values.to_numpy(dtype=float)
    if values.size == 0:
        return np.array([], dtype=float)
    if not np.isfinite(values).all() or ((values < 0) | (values > 1)).any():
        raise ValueError("P values must be finite and between zero and one")
    order = np.argsort(values)
    ranked = values[order] * len(values) / np.arange(1, len(values) + 1)
    adjusted_sorted = np.clip(np.minimum.accumulate(ranked[::-1])[::-1], None, 1.0)
    adjusted = np.empty_like(adjusted_sorted)
    adjusted[order] = adjusted_sorted
    return adjusted


def run_differential_expression(
    args: argparse.Namespace, protein_coding_genes: set[str]
) -> tuple[pd.DataFrame, pd.DataFrame, set[str]]:
    try:
        import scanpy as sc
    except ImportError as error:
        raise RuntimeError("scanpy is required for Tabula Sapiens analysis") from error

    backed = sc.read_h5ad(args.h5ad, backed="r")
    try:
        if args.group_column not in backed.obs.columns:
            raise ValueError(f"AnnData obs is missing {args.group_column}")
        if args.gene_id_column is None:
            gene_ids = pd.Series(backed.var_names, index=backed.var_names, dtype="string")
        else:
            if args.gene_id_column not in backed.var.columns:
                raise ValueError(f"AnnData var is missing {args.gene_id_column}")
            gene_ids = pd.Series(
                backed.var[args.gene_id_column].to_numpy(),
                index=backed.var_names,
                dtype="string",
            )
        gene_ids = normalize_gene_ids(gene_ids, args.strip_gene_version)
        selected = gene_ids.isin(protein_coding_genes)
        selected_ids = gene_ids.loc[selected]
        if selected_ids.empty:
            raise ValueError("No protein-coding genes matched the AnnData variables")
        duplicates = selected_ids.loc[selected_ids.duplicated(keep=False)].dropna().unique()
        if len(duplicates):
            preview = ", ".join(map(str, duplicates[:5]))
            raise ValueError(f"AnnData gene identifiers are not unique: {preview}")

        observed_groups = backed.obs[args.group_column]
        observed = observed_groups.notna().to_numpy()
        if not observed.any():
            raise ValueError(f"AnnData obs column {args.group_column} contains no labels")
        adata = backed[observed, selected.to_numpy()].to_memory()
        adata.var_names = pd.Index(selected_ids.astype(str), name="gene")
    finally:
        backed.file.close()

    adata.obs[args.group_column] = (
        adata.obs[args.group_column]
        .astype("string")
        .str.strip()
        .replace("", pd.NA)
        .astype("category")
    )
    adata = adata[adata.obs[args.group_column].notna()].copy()
    adata.obs[args.group_column] = adata.obs[args.group_column].cat.remove_unused_categories()
    class_counts = (
        adata.obs[args.group_column]
        .value_counts(sort=False)
        .rename_axis("cell_class")
        .reset_index(name="cell_count")
        .sort_values("cell_class")
        .reset_index(drop=True)
    )
    if args.expected_classes and len(class_counts) != args.expected_classes:
        raise ValueError(
            f"Expected {args.expected_classes} broad cell classes, found {len(class_counts)}"
        )
    small_classes = class_counts.loc[
        class_counts["cell_count"].lt(args.min_cells_per_class), "cell_class"
    ]
    if not small_classes.empty:
        raise ValueError(
            "Cell classes below --min-cells-per-class: " + ", ".join(small_classes.astype(str))
        )

    sc.pp.filter_genes(adata, min_cells=1)
    if adata.n_vars == 0:
        raise ValueError("No protein-coding genes were detected in the selected cells")
    universe = set(adata.var_names.astype(str))
    sc.tl.rank_genes_groups(
        adata,
        groupby=args.group_column,
        reference="rest",
        method="t-test",
        corr_method="benjamini-hochberg",
        use_raw=False,
        n_genes=adata.n_vars,
    )
    results = sc.get.rank_genes_groups_df(adata, group=None).rename(
        columns={
            "group": "cell_class",
            "names": "gene",
            "scores": "t_statistic",
            "logfoldchanges": "log2_fold_change",
            "pvals": "p_value",
            "pvals_adj": "fdr",
        }
    )
    result_columns = [
        "cell_class",
        "gene",
        "t_statistic",
        "log2_fold_change",
        "p_value",
        "fdr",
    ]
    require_columns(results, set(result_columns), "Differential-expression result")
    results = results[result_columns].copy()
    results["cell_class"] = results["cell_class"].astype(str)
    results["gene"] = results["gene"].astype(str)
    for column in ["t_statistic", "log2_fold_change", "p_value", "fdr"]:
        results[column] = pd.to_numeric(results[column], errors="raise")
    results = results.sort_values(
        ["cell_class", "fdr", "p_value", "gene"], kind="stable"
    ).reset_index(drop=True)
    return results, class_counts, universe


def calculate_enrichment(
    markers: pd.DataFrame,
    cell_classes: list[str],
    candidate_sets: dict[str, set[str]],
    universe: set[str],
    fdr_threshold: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    marker_sets = {
        name: set(group["gene"]).intersection(universe)
        for name, group in markers.groupby("cell_class", sort=False)
    }
    records = []
    coverage = []
    for set_name, genes in sorted(candidate_sets.items()):
        tested_candidates = genes.intersection(universe)
        coverage.append(
            {
                "candidate_set": set_name,
                "input_gene_count": len(genes),
                "tested_gene_count": len(tested_candidates),
                "unmapped_gene_count": len(genes.difference(universe)),
            }
        )
        if not tested_candidates:
            raise ValueError(f"Candidate set {set_name} has no genes in the analysis universe")
        for cell_class in cell_classes:
            marker_set = marker_sets.get(cell_class, set())
            overlap = tested_candidates.intersection(marker_set)
            a = len(overlap)
            b = len(marker_set) - a
            c = len(tested_candidates) - a
            d = len(universe) - a - b - c
            if min(a, b, c, d) < 0:
                raise ValueError("Invalid Fisher exact-test contingency table")
            odds_ratio, p_value = fisher_exact([[a, b], [c, d]], alternative="greater")
            records.append(
                {
                    "candidate_set": set_name,
                    "cell_class": cell_class,
                    "universe_gene_count": len(universe),
                    "candidate_gene_count": len(tested_candidates),
                    "marker_gene_count": len(marker_set),
                    "overlap_gene_count": a,
                    "overlap_genes": ",".join(sorted(overlap)),
                    "odds_ratio": float(odds_ratio),
                    "p_value": float(p_value),
                }
            )

    results = pd.DataFrame(records)
    results["fdr"] = results.groupby("candidate_set", sort=False)["p_value"].transform(
        benjamini_hochberg
    )
    results["significant"] = results["fdr"].lt(fdr_threshold)
    results = results.sort_values(
        ["candidate_set", "fdr", "p_value", "cell_class"], kind="stable"
    ).reset_index(drop=True)
    return results, pd.DataFrame(coverage).sort_values("candidate_set").reset_index(drop=True)


def validate_args(args: argparse.Namespace) -> None:
    if args.min_cells_per_class < 2:
        raise ValueError("--min-cells-per-class must be at least 2")
    if args.expected_classes < 0:
        raise ValueError("--expected-classes cannot be negative")
    if args.log2fc_threshold < 0:
        raise ValueError("--log2fc-threshold cannot be negative")
    if not 0 < args.fdr_threshold < 1:
        raise ValueError("--fdr-threshold must be between zero and one")


def main() -> None:
    args = parse_args()
    validate_args(args)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    candidate_sets = load_candidate_sets(args.candidate_sets, args.strip_gene_version)
    protein_coding_genes = load_protein_coding_genes(
        args.protein_coding_genes,
        args.protein_coding_column,
        args.strip_gene_version,
    )
    differential_expression, class_counts, universe = run_differential_expression(
        args, protein_coding_genes
    )
    markers = differential_expression.loc[
        differential_expression["fdr"].lt(args.fdr_threshold)
        & differential_expression["log2_fold_change"].gt(args.log2fc_threshold)
    ].copy()
    enrichment, coverage = calculate_enrichment(
        markers,
        class_counts["cell_class"].astype(str).tolist(),
        candidate_sets,
        universe,
        args.fdr_threshold,
    )

    differential_expression.to_csv(
        args.output_dir / "broad_cell_class_differential_expression.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    markers.to_csv(
        args.output_dir / "broad_cell_class_marker_genes.tsv", sep="\t", index=False
    )
    enrichment.to_csv(
        args.output_dir / "broad_cell_class_enrichment.tsv", sep="\t", index=False
    )
    class_counts.to_csv(
        args.output_dir / "broad_cell_class_counts.tsv", sep="\t", index=False
    )
    coverage.to_csv(args.output_dir / "candidate_set_coverage.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
