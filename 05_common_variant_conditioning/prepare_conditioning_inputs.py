#!/usr/bin/env python3
"""Map candidate genes to published GWAS loci and prepare lead-SNP covariates."""

from __future__ import annotations

import argparse
from pathlib import Path
from urllib.parse import quote

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", required=True, type=Path)
    parser.add_argument("--gwas-loci", required=True, type=Path)
    parser.add_argument("--lead-snp-dosages", required=True, type=Path)
    parser.add_argument("--phenotype-dir", required=True, type=Path)
    parser.add_argument("--group-file-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--genome-build", required=True, choices=("GRCh37", "GRCh38"))
    return parser.parse_args()


def require_columns(frame: pd.DataFrame, columns: set[str], label: str) -> None:
    missing = columns.difference(frame.columns)
    if missing:
        raise ValueError(f"{label} is missing columns: {', '.join(sorted(missing))}")


def map_overlaps(candidates: pd.DataFrame, loci: pd.DataFrame) -> pd.DataFrame:
    records = []
    grouped_loci = {
        key: group for key, group in loci.groupby(["trait", "chromosome"], sort=False)
    }
    for candidate in candidates.itertuples(index=False):
        subset = grouped_loci.get((candidate.trait, candidate.chromosome))
        if subset is None:
            continue
        overlap = subset.loc[
            subset["start"].le(candidate.end) & subset["end"].ge(candidate.start)
        ]
        for locus in overlap.itertuples(index=False):
            records.append(
                {
                    "trait": candidate.trait,
                    "gene": candidate.gene,
                    "chromosome": candidate.chromosome,
                    "gene_start": candidate.start,
                    "gene_end": candidate.end,
                    "locus_id": locus.locus_id,
                    "lead_snp": locus.lead_snp,
                }
            )
    return pd.DataFrame(records)


def extract_group_lines(source: Path, gene: str) -> list[str]:
    matches = []
    with source.open("r", encoding="utf-8") as handle:
        for line in handle:
            first_field = line.split("\t", maxsplit=1)[0]
            if first_field.split("|", maxsplit=1)[0] == gene:
                matches.append(line)
    return matches


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    candidates = pd.read_csv(args.candidates, sep="\t")
    loci = pd.read_csv(args.gwas_loci, sep="\t")
    candidate_columns = {"trait", "gene", "chromosome", "start", "end", "genome_build"}
    locus_columns = {
        "trait",
        "locus_id",
        "chromosome",
        "start",
        "end",
        "lead_snp",
        "genome_build",
    }
    require_columns(candidates, candidate_columns, "Candidate table")
    require_columns(loci, locus_columns, "GWAS-locus table")

    if not candidates["genome_build"].eq(args.genome_build).all():
        raise ValueError("Candidate coordinates do not all match --genome-build")
    if not loci["genome_build"].eq(args.genome_build).all():
        raise ValueError("GWAS-locus coordinates do not all match --genome-build")

    mapping = map_overlaps(candidates, loci)
    if mapping.empty:
        raise ValueError("No candidate gene overlapped a trait-matched GWAS locus")
    mapping = mapping.drop_duplicates(["trait", "gene", "lead_snp"])

    dosages = pd.read_csv(args.lead_snp_dosages, sep="\t", dtype={"sample_id": str})
    require_columns(dosages, {"sample_id"}, "Dosage table")
    job_rows = []
    for row in mapping.itertuples(index=False):
        if row.lead_snp not in dosages.columns:
            raise ValueError(f"Lead-SNP dosage is missing for {row.lead_snp}")

        phenotype_path = args.phenotype_dir / f"{row.trait}.tsv"
        phenotype = pd.read_csv(phenotype_path, sep="\t", dtype={"sample_id": str})
        require_columns(phenotype, {"sample_id"}, str(phenotype_path))
        merged = phenotype.merge(
            dosages[["sample_id", row.lead_snp]], on="sample_id", how="left", validate="one_to_one"
        )
        mean_dosage = merged[row.lead_snp].mean(skipna=True)
        if pd.isna(mean_dosage):
            raise ValueError(f"All dosages are missing for {row.lead_snp}")
        merged["LeadSNP_Dosage"] = merged[row.lead_snp].fillna(mean_dosage)
        merged = merged.drop(columns=row.lead_snp)

        lead_snp_dir = quote(str(row.lead_snp), safe="")
        run_dir = args.output_dir / row.trait / row.gene / lead_snp_dir
        run_dir.mkdir(parents=True, exist_ok=True)
        conditioned_phenotype = run_dir / "phenotype_with_lead_snp.tsv"
        merged.to_csv(conditioned_phenotype, sep="\t", index=False)

        chromosome = str(row.chromosome).removeprefix("chr")
        source_group = args.group_file_dir / f"chr{chromosome}.group.txt"
        lines = extract_group_lines(source_group, row.gene)
        if not lines:
            raise ValueError(f"No SAIGE group records found for {row.gene} in {source_group}")
        conditioned_group = run_dir / "gene.group.txt"
        conditioned_group.write_text("".join(lines), encoding="utf-8")

        job_rows.append(
            {
                **row._asdict(),
                "phenotype_file": str(conditioned_phenotype.resolve()),
                "group_file": str(conditioned_group.resolve()),
            }
        )

    mapping.to_csv(args.output_dir / "candidate_gwas_mapping.tsv", sep="\t", index=False)
    pd.DataFrame(job_rows).to_csv(
        args.output_dir / "conditioning_jobs_local.tsv", sep="\t", index=False
    )


if __name__ == "__main__":
    main()
