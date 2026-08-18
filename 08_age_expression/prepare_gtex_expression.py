#!/usr/bin/env python3
"""Prepare tissue-level GTEx v8 expression and metadata for aging analyses."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd


SAMPLE_COLUMNS = {
    "SAMPID": "sample_id",
    "SMTSD": "tissue",
    "SMGEBTCH": "batch",
    "SMCENTER": "center",
    "SMRIN": "rin",
    "SMTSISCH": "ischemic_time",
}
SUBJECT_COLUMNS = {
    "SUBJID": "subject_id",
    "AGE": "age",
    "SEX": "sex",
    "DTHHRDY": "death_cause",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Apply the GTEx v8 expression filters described in the manuscript and "
            "write one expression matrix and metadata table per tissue."
        )
    )
    parser.add_argument("sample_attributes", type=Path)
    parser.add_argument("subject_phenotypes", type=Path)
    parser.add_argument("gene_tpm_gct", type=Path)
    parser.add_argument(
        "protein_coding_genes",
        type=Path,
        help="Table containing an Ensembl gene identifier column named gene or gene_id",
    )
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--rin-min", type=float, default=6.0)
    parser.add_argument("--tpm-min", type=float, default=0.1)
    parser.add_argument("--minimum-expressed-fraction", type=float, default=0.5)
    return parser.parse_args()


def read_table(path: Path) -> pd.DataFrame:
    suffixes = "".join(path.suffixes).lower()
    separator = "," if suffixes.endswith(".csv") else "\t"
    return pd.read_csv(path, sep=separator, low_memory=False)


def normalize_gene_id(value: object) -> str:
    gene = str(value).strip()
    if "|" in gene:
        ensembl_tokens = [token for token in gene.split("|") if token.startswith("ENSG")]
        gene = ensembl_tokens[0] if ensembl_tokens else gene.split("|")[0]
    return gene.split(".")[0]


def age_midpoint(value: object) -> float:
    if pd.isna(value):
        return np.nan
    if isinstance(value, (int, float, np.integer, np.floating)):
        return float(value)
    match = re.fullmatch(r"\s*(\d+)\s*-\s*(\d+)\s*", str(value))
    if match:
        return (float(match.group(1)) + float(match.group(2))) / 2.0
    try:
        return float(value)
    except ValueError:
        return np.nan


def tissue_key(tissue: str) -> str:
    key = re.sub(r"[^A-Za-z0-9]+", "_", tissue).strip("_")
    if not key:
        raise ValueError(f"Cannot create a file-safe name for tissue: {tissue!r}")
    return key


def load_protein_coding_ids(path: Path) -> set[str]:
    table = read_table(path)
    normalized = {column.lower(): column for column in table.columns}
    column = normalized.get("gene_id") or normalized.get("gene")
    if column is None:
        raise ValueError("Protein-coding gene table requires a gene or gene_id column")
    return {
        normalize_gene_id(value)
        for value in table[column].dropna()
        if normalize_gene_id(value)
    }


def load_metadata(sample_path: Path, subject_path: Path, rin_min: float) -> pd.DataFrame:
    samples = pd.read_csv(sample_path, sep="\t", low_memory=False)
    subjects = pd.read_csv(subject_path, sep="\t", low_memory=False)

    missing_sample = set(SAMPLE_COLUMNS).difference(samples.columns)
    missing_subject = set(SUBJECT_COLUMNS).difference(subjects.columns)
    if missing_sample or missing_subject:
        raise ValueError(
            "Missing GTEx metadata columns: "
            + ", ".join(sorted(missing_sample | missing_subject))
        )

    samples = samples.loc[:, list(SAMPLE_COLUMNS)].rename(columns=SAMPLE_COLUMNS)
    subjects = subjects.loc[:, list(SUBJECT_COLUMNS)].rename(columns=SUBJECT_COLUMNS)
    samples["subject_id"] = samples["sample_id"].str.split("-").str[:2].str.join("-")
    metadata = samples.merge(subjects, on="subject_id", how="left", validate="many_to_one")

    metadata["age"] = metadata["age"].map(age_midpoint)
    metadata["sex"] = metadata["sex"].map({1: "Male", 2: "Female"}).fillna(
        metadata["sex"].astype("string")
    )
    for column in ("rin", "ischemic_time"):
        metadata[column] = pd.to_numeric(metadata[column], errors="coerce")

    metadata = metadata.dropna(subset=["sample_id", "tissue", "age", "rin"])
    metadata = metadata.loc[metadata["rin"] >= rin_min].copy()
    return metadata.drop_duplicates("sample_id")


def load_expression(path: Path, protein_coding_ids: set[str]) -> pd.DataFrame:
    expression = pd.read_csv(path, sep="\t", skiprows=2, low_memory=False)
    required = {"Name", "Description"}
    if not required.issubset(expression.columns):
        raise ValueError("GTEx GCT file requires Name and Description columns")

    expression["gene"] = expression["Name"].map(normalize_gene_id)
    expression = expression.loc[expression["gene"].isin(protein_coding_ids)]
    expression = expression.drop(columns=["Name", "Description"])
    expression = expression.set_index("gene")
    expression = expression.apply(pd.to_numeric, errors="coerce")
    if expression.index.has_duplicates:
        expression = expression.groupby(level=0, sort=False).mean()
    return expression


def main() -> None:
    args = parse_args()
    if not 0 < args.minimum_expressed_fraction <= 1:
        raise ValueError("--minimum-expressed-fraction must be in (0, 1]")

    expression_dir = args.output_dir / "expression"
    metadata_dir = args.output_dir / "metadata"
    expression_dir.mkdir(parents=True, exist_ok=True)
    metadata_dir.mkdir(parents=True, exist_ok=True)

    protein_coding_ids = load_protein_coding_ids(args.protein_coding_genes)
    metadata = load_metadata(args.sample_attributes, args.subject_phenotypes, args.rin_min)
    expression = load_expression(args.gene_tpm_gct, protein_coding_ids)

    available_samples = [sample for sample in metadata["sample_id"] if sample in expression.columns]
    metadata = metadata.set_index("sample_id").loc[available_samples].reset_index()
    expression = expression.loc[:, available_samples]

    manifest_rows: list[dict[str, object]] = []
    for tissue in sorted(metadata["tissue"].unique()):
        tissue_metadata = metadata.loc[metadata["tissue"] == tissue].copy()
        samples = tissue_metadata["sample_id"].tolist()
        tissue_tpm = expression.loc[:, samples]
        expressed_fraction = tissue_tpm.ge(args.tpm_min).mean(axis=1)
        tissue_tpm = tissue_tpm.loc[expressed_fraction >= args.minimum_expressed_fraction]
        tissue_expression = np.log2(tissue_tpm + 1.0)

        key = tissue_key(tissue)
        expression_relative = Path("expression") / f"{key}.expression.tsv.gz"
        metadata_relative = Path("metadata") / f"{key}.metadata.tsv"
        tissue_expression.rename_axis("gene").to_csv(
            args.output_dir / expression_relative, sep="\t", compression="gzip"
        )
        tissue_metadata.to_csv(args.output_dir / metadata_relative, sep="\t", index=False)

        manifest_rows.append(
            {
                "tissue": tissue,
                "tissue_key": key,
                "samples": len(samples),
                "genes": tissue_expression.shape[0],
                "expression_file": expression_relative.as_posix(),
                "metadata_file": metadata_relative.as_posix(),
            }
        )

    if not manifest_rows:
        raise RuntimeError("No GTEx samples remained after metadata and RIN filtering")
    pd.DataFrame(manifest_rows).to_csv(
        args.output_dir / "tissue_manifest.tsv", sep="\t", index=False
    )


if __name__ == "__main__":
    main()
