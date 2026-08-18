#!/usr/bin/env python3
"""Combine MR results with rare-variant and carrier-protein evidence."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mr-dir", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--rare-evidence", required=True, type=Path)
    parser.add_argument("--protein-evidence", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def benjamini_hochberg(values: pd.Series) -> np.ndarray:
    p_values = values.to_numpy(dtype=float)
    order = np.argsort(p_values)
    ranked = p_values[order] * len(p_values) / np.arange(1, len(p_values) + 1)
    adjusted_sorted = np.clip(np.minimum.accumulate(ranked[::-1])[::-1], None, 1.0)
    adjusted = np.empty_like(adjusted_sorted)
    adjusted[order] = adjusted_sorted
    return adjusted


def read_optional(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def main() -> None:
    args = parse_args()
    manifest = pd.read_csv(args.manifest, sep="\t")
    required_manifest = {"pair_id", "protein", "trait"}
    if not required_manifest.issubset(manifest.columns):
        raise ValueError("Manifest requires pair_id, protein, and trait columns")

    records = []
    method_records = []
    for row in manifest.itertuples(index=False):
        run_dir = args.mr_dir / row.pair_id
        mr = read_optional(run_dir / "mr_results.tsv")
        if mr.empty:
            continue
        mr["pair_id"] = row.pair_id
        mr["protein"] = row.protein
        mr["trait"] = row.trait
        method_records.append(mr)

        ivw = mr.loc[mr["method"].str.fullmatch("Inverse variance weighted", na=False)]
        if ivw.empty:
            ivw = mr.loc[mr["method"].str.contains("Inverse variance weighted", na=False)]
        if ivw.empty:
            continue
        ivw = ivw.iloc[0]
        sensitivity_methods = {"MR Egger", "Weighted median", "Weighted mode"}
        sensitivity = mr.loc[mr["method"].isin(sensitivity_methods) & mr["b"].notna()]
        concordant = bool(
            sensitivity_methods.issubset(set(sensitivity["method"]))
            and np.all(np.sign(sensitivity["b"].to_numpy(dtype=float)) == np.sign(float(ivw["b"])))
        )

        pleiotropy = read_optional(run_dir / "horizontal_pleiotropy.tsv")
        pleiotropy_p = float(pleiotropy["pval"].iloc[0]) if not pleiotropy.empty else np.nan
        heterogeneity = read_optional(run_dir / "heterogeneity.tsv")
        ivw_heterogeneity = heterogeneity.loc[
            heterogeneity["method"].str.contains("Inverse variance weighted", na=False)
        ] if not heterogeneity.empty else pd.DataFrame()
        heterogeneity_p = (
            float(ivw_heterogeneity["Q_pval"].iloc[0]) if not ivw_heterogeneity.empty else np.nan
        )
        records.append(
            {
                "pair_id": row.pair_id,
                "protein": row.protein,
                "trait": row.trait,
                "instruments": int(ivw["nsnp"]),
                "ivw_beta": float(ivw["b"]),
                "ivw_standard_error": float(ivw["se"]),
                "ivw_p_value": float(ivw["pval"]),
                "sensitivity_directions_concordant": concordant,
                "egger_intercept_p_value": pleiotropy_p,
                "cochran_q_p_value": heterogeneity_p,
            }
        )

    summary = pd.DataFrame(records)
    if summary.empty:
        raise ValueError("No complete IVW results were found")
    summary["ivw_fdr"] = benjamini_hochberg(summary["ivw_p_value"])
    summary["credible_causal_association"] = (
        summary["ivw_fdr"].lt(0.05)
        & summary["sensitivity_directions_concordant"]
        & summary["egger_intercept_p_value"].ge(0.05)
    )

    rare = pd.read_csv(args.rare_evidence, sep="\t")
    protein = pd.read_csv(args.protein_evidence, sep="\t")
    if not {"protein", "trait", "rare_beta", "rare_p_value"}.issubset(rare.columns):
        raise ValueError("Rare evidence requires protein, trait, rare_beta, and rare_p_value")
    if not {"protein", "beta", "p_value", "fdr"}.issubset(protein.columns):
        raise ValueError("Protein evidence requires protein, beta, p_value, and fdr")
    protein = protein.rename(
        columns={"beta": "carrier_protein_beta", "p_value": "carrier_protein_p_value", "fdr": "carrier_protein_fdr"}
    )
    summary = summary.merge(rare, on=["protein", "trait"], how="left")
    summary = summary.merge(
        protein[["protein", "carrier_protein_beta", "carrier_protein_p_value", "carrier_protein_fdr"]],
        on="protein",
        how="left",
    )
    summary["evidence_direction"] = (
        np.sign(summary["rare_beta"]).astype("Int64").astype(str)
        + "/"
        + np.sign(summary["carrier_protein_beta"]).astype("Int64").astype(str)
        + "/"
        + np.sign(summary["ivw_beta"]).astype("Int64").astype(str)
    )
    summary = summary.sort_values(["ivw_fdr", "ivw_p_value", "protein", "trait"])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.output, sep="\t", index=False)

    if method_records:
        pd.concat(method_records, ignore_index=True).to_csv(
            args.output.with_name(args.output.stem + "_all_methods.tsv"), sep="\t", index=False
        )


if __name__ == "__main__":
    main()
