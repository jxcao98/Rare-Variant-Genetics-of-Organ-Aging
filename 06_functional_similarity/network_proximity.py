#!/usr/bin/env python3
"""Test Human Gene Connectome proximity between candidate and Aging Atlas gene sets."""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


SASP = "senescenceassociated_secretory_phenotype"
NFKB = "NF_kB_related_gene"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-sets", required=True, type=Path)
    parser.add_argument("--aging-atlas", required=True, type=Path)
    parser.add_argument("--distances", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--permutations", type=int, default=100_000)
    parser.add_argument("--seed", type=int, default=66)
    parser.add_argument("--chunk-size", type=int, default=5_000)
    return parser.parse_args()


def load_sets(candidate_path: Path, aging_path: Path) -> tuple[dict[str, set[str]], dict[str, set[str]]]:
    candidates = pd.read_csv(candidate_path, sep="\t")
    aging = pd.read_csv(aging_path, sep="\t")
    if not {"gene", "set"}.issubset(candidates.columns):
        raise ValueError("Candidate input requires gene and set columns")
    if not {"gene", "hallmark"}.issubset(aging.columns):
        raise ValueError("Aging Atlas input requires gene and hallmark columns")

    candidate_sets = {
        name: set(group["gene"].dropna().astype(str))
        for name, group in candidates.groupby("set", sort=True)
    }
    aging_sets = {
        name: set(group["gene"].dropna().astype(str))
        for name, group in aging.groupby("hallmark", sort=True)
        if name not in {SASP, NFKB}
    }
    chronic = set(aging.loc[aging["hallmark"].isin([SASP, NFKB]), "gene"].dropna().astype(str))
    if chronic:
        aging_sets["Chronic_Inflammation"] = chronic
    aging_sets["Total_Aging"] = set(aging["gene"].dropna().astype(str))
    return candidate_sets, aging_sets


def stream_distance_means(
    path: Path, aging_sets: dict[str, set[str]]
) -> tuple[list[str], np.ndarray, list[str]]:
    hallmark_names = sorted(aging_sets)
    hallmark_index = {name: index for index, name in enumerate(hallmark_names)}
    gene_memberships: dict[str, list[int]] = defaultdict(list)
    for hallmark, genes in aging_sets.items():
        index = hallmark_index[hallmark]
        for gene in genes:
            gene_memberships[gene].append(index)

    sums: dict[str, np.ndarray] = {}
    counts: dict[str, np.ndarray] = {}
    background: set[str] = set()

    def update(target: str, memberships: list[int], distance: float) -> None:
        if target not in sums:
            sums[target] = np.zeros(len(hallmark_names), dtype=float)
            counts[target] = np.zeros(len(hallmark_names), dtype=np.int32)
        sums[target][memberships] += distance
        counts[target][memberships] += 1

    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.split()
            if len(fields) < 3:
                continue
            gene1, gene2 = fields[0], fields[1]
            try:
                distance = float(fields[2])
            except ValueError:
                if line_number == 1:
                    continue
                raise
            background.update((gene1, gene2))
            if gene1 in gene_memberships:
                update(gene2, gene_memberships[gene1], distance)
            if gene2 in gene_memberships:
                update(gene1, gene_memberships[gene2], distance)

    for gene, memberships in gene_memberships.items():
        background.add(gene)
        update(gene, memberships, 0.0)

    genes = sorted(background)
    means = np.full((len(genes), len(hallmark_names)), np.nan, dtype=float)
    for row, gene in enumerate(genes):
        if gene not in sums:
            continue
        valid = counts[gene] > 0
        means[row, valid] = sums[gene][valid] / counts[gene][valid]
    return genes, means, hallmark_names


def permutation_p_value(
    values: np.ndarray,
    observed: float,
    set_size: int,
    permutations: int,
    rng: np.random.Generator,
    chunk_size: int,
) -> float:
    extreme = 0
    completed = 0
    while completed < permutations:
        batch = min(chunk_size, permutations - completed)
        draws = np.vstack([rng.choice(len(values), size=set_size, replace=False) for _ in range(batch)])
        null_means = values[draws].mean(axis=1)
        extreme += int(np.count_nonzero(null_means <= observed))
        completed += batch
    return (extreme + 1) / (permutations + 1)


def benjamini_hochberg(p_values: pd.Series) -> np.ndarray:
    values = p_values.to_numpy(dtype=float)
    order = np.argsort(values)
    ranked = values[order] * len(values) / np.arange(1, len(values) + 1)
    adjusted_sorted = np.clip(np.minimum.accumulate(ranked[::-1])[::-1], None, 1.0)
    adjusted = np.empty_like(adjusted_sorted)
    adjusted[order] = adjusted_sorted
    return adjusted


def main() -> None:
    args = parse_args()
    candidate_sets, aging_sets = load_sets(args.candidate_sets, args.aging_atlas)
    genes, distance_means, hallmark_names = stream_distance_means(args.distances, aging_sets)
    gene_index = {gene: index for index, gene in enumerate(genes)}
    rng = np.random.default_rng(args.seed)

    records = []
    for candidate_name, candidate_genes in candidate_sets.items():
        candidate_rows = [gene_index[gene] for gene in candidate_genes if gene in gene_index]
        for hallmark_column, hallmark_name in enumerate(hallmark_names):
            background_values = distance_means[:, hallmark_column]
            background_values = background_values[np.isfinite(background_values)]
            observed_values = distance_means[candidate_rows, hallmark_column]
            observed_values = observed_values[np.isfinite(observed_values)]
            if len(observed_values) == 0 or len(background_values) < len(observed_values):
                continue
            observed = float(observed_values.mean())
            p_value = permutation_p_value(
                background_values,
                observed,
                len(observed_values),
                args.permutations,
                rng,
                args.chunk_size,
            )
            records.append(
                {
                    "candidate_set": candidate_name,
                    "aging_set": hallmark_name,
                    "mean_network_distance": observed,
                    "p_value": p_value,
                    "permutations": args.permutations,
                    "candidate_genes_tested": len(observed_values),
                    "candidate_genes_total": len(candidate_genes),
                    "aging_genes": len(aging_sets[hallmark_name]),
                    "overlapping_genes": len(candidate_genes & aging_sets[hallmark_name]),
                }
            )

    result = pd.DataFrame(records)
    if not result.empty:
        result["fdr"] = benjamini_hochberg(result["p_value"])
        result = result.sort_values(["p_value", "candidate_set", "aging_set"])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
