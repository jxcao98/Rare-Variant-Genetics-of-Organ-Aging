#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: $0 MANIFEST.tsv SAIGE_OR_PLINK_IMAGE_FILE DX_DESTINATION" >&2
  echo "Manifest columns: trait, chrom, phenotype_column, phenotype_file, pgen, pvar, psam" >&2
}

if [[ $# -ne 3 ]]; then
  usage
  exit 2
fi

manifest=$1
image_file=$2
destination=${3%/}

command -v dx >/dev/null 2>&1 || { echo "The DNAnexus CLI (dx) is required." >&2; exit 1; }
[[ -f "$manifest" ]] || { echo "Manifest not found: $manifest" >&2; exit 1; }

covariates="Sex,Age,Age2,AssessmentCenter,TDI"
for pc in $(seq 1 20); do
  covariates+=",PC${pc}"
done

object_name() {
  dx describe "$1" --name
}

tail -n +2 "$manifest" | while IFS=$'\t' read -r trait chrom phenotype_column phenotype_file pgen pvar psam; do
  [[ -n "$trait" ]] || continue

  phenotype_name=$(object_name "$phenotype_file")
  pgen_name=$(object_name "$pgen")
  pfile_prefix=${pgen_name%.pgen}
  output_prefix="${trait}.chr${chrom}"

  analysis_command="plink2 --pfile '${pfile_prefix}' \
    --pheno '${phenotype_name}' --pheno-name '${phenotype_column}' \
    --covar '${phenotype_name}' --covar-name '${covariates}' \
    --covar-variance-standardize --glm hide-covar no-x-sex \
    --no-input-missing-phenotype --out '${output_prefix}'"

  dx run swiss-army-knife \
    -iin="$phenotype_file" \
    -iin="$pgen" \
    -iin="$pvar" \
    -iin="$psam" \
    -iimage_file="$image_file" \
    -icmd="$analysis_command" \
    --name="BHR summary statistics: ${trait} chr${chrom}" \
    --destination="${destination}/${trait}/" \
    --brief -y
done

