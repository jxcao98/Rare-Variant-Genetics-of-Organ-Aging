#!/usr/bin/env python3
"""Merge PLINK 2 association results with rare-variant annotations for BHR."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


VALID_CLASSES = ("LoF", "DMis", "Synonymous")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trait", required=True)
    parser.add_argument("--association-dir", required=True, type=Path)
    parser.add_argument("--association-pattern", default="*.glm.linear")
    parser.add_argument("--variant-annotation", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def load_associations(directory: Path, pattern: str) -> pd.DataFrame:
    paths = sorted(directory.glob(pattern))
    if not paths:
        raise FileNotFoundError(f"No association files matched {directory / pattern}")

    frames = []
    required = {"ID", "OBS_CT", "A1_FREQ", "BETA", "SE", "P"}
    for path in paths:
        frame = pd.read_csv(path, sep="\t", low_memory=False)
        missing = required.difference(frame.columns)
        if missing:
            raise ValueError(f"{path} is missing columns: {', '.join(sorted(missing))}")
        frames.append(frame.loc[frame["P"].notna(), sorted(required)])
    return pd.concat(frames, ignore_index=True)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    associations = load_associations(args.association_dir, args.association_pattern)
    associations = associations.rename(
        columns={"ID": "variant_id", "OBS_CT": "N", "BETA": "beta", "SE": "se"}
    )
    associations["A1_FREQ"] = pd.to_numeric(associations["A1_FREQ"], errors="coerce")
    associations["N"] = pd.to_numeric(associations["N"], errors="coerce")
    associations["sample_maf"] = np.minimum(
        associations["A1_FREQ"], 1.0 - associations["A1_FREQ"]
    )
    associations["MAC"] = 2.0 * associations["N"] * associations["sample_maf"]

    annotation = pd.read_csv(args.variant_annotation, sep="\t", low_memory=False)
    required_annotation = {
        "variant_id",
        "gene",
        "chromosome",
        "gene_position",
        "variant_class",
        "AF",
    }
    missing = required_annotation.difference(annotation.columns)
    if missing:
        raise ValueError(
            f"Variant annotation is missing columns: {', '.join(sorted(missing))}"
        )
    annotation = annotation[
        ["variant_id", "gene", "chromosome", "gene_position", "variant_class", "AF"]
    ].drop_duplicates()
    annotation["AF"] = pd.to_numeric(annotation["AF"], errors="coerce")
    annotation["AF"] = np.minimum(annotation["AF"], 1.0 - annotation["AF"])

    merged = associations.merge(annotation, on="variant_id", how="inner", validate="many_to_many")
    merged = merged.dropna(subset=["gene", "AF", "beta", "se", "N"])
    merged["phenotype_key"] = args.trait

    output_columns = [
        "variant_id",
        "gene",
        "chromosome",
        "gene_position",
        "N",
        "beta",
        "se",
        "AF",
        "A1_FREQ",
        "MAC",
        "phenotype_key",
    ]
    manifest_rows = []
    for variant_class in VALID_CLASSES:
        subset = merged.loc[merged["variant_class"].eq(variant_class), output_columns].copy()
        output_path = args.output_dir / f"{args.trait}.{variant_class}.tsv"
        subset.to_csv(output_path, sep="\t", index=False)
        manifest_rows.append(
            {
                "trait": args.trait,
                "variant_class": variant_class,
                "path": str(output_path.resolve()),
                "variant_count": len(subset),
            }
        )

    pd.DataFrame(manifest_rows).to_csv(
        args.output_dir / f"{args.trait}.manifest.tsv", sep="\t", index=False
    )


if __name__ == "__main__":
    main()
