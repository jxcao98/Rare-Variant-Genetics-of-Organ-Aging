#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat >&2 <<'EOF'
Usage:
  run_saige_gene_dnanexus.sh step1 STEP1_MANIFEST.tsv SAIGE_IMAGE_FILE DX_DESTINATION
  run_saige_gene_dnanexus.sh step2 STEP2_MANIFEST.tsv SAIGE_IMAGE_FILE DX_DESTINATION

Step 1 columns:
  trait phenotype_column phenotype_file bed bim fam sparse_grm sparse_grm_ids

Step 2 columns:
  trait chrom model_file variance_ratio_file bgen bgi sample group_file sparse_grm sparse_grm_ids
EOF
}

if [[ $# -ne 4 ]]; then
  usage
  exit 2
fi

stage=$1
manifest=$2
image_file=$3
destination=${4%/}

command -v dx >/dev/null 2>&1 || { echo "The DNAnexus CLI (dx) is required." >&2; exit 1; }
[[ -f "$manifest" ]] || { echo "Manifest not found: $manifest" >&2; exit 1; }

object_name() {
  dx describe "$1" --name
}

covariates="Sex,Age,Age2,AssessmentCenter,TDI"
for pc in $(seq 1 20); do
  covariates+=",PC${pc}"
done

if [[ "$stage" == "step1" ]]; then
  tail -n +2 "$manifest" | while IFS=$'\t' read -r trait phenotype_column phenotype_file bed bim fam sparse_grm sparse_grm_ids; do
    [[ -n "$trait" ]] || continue
    phenotype_name=$(object_name "$phenotype_file")
    bed_name=$(object_name "$bed")
    plink_prefix=${bed_name%.bed}
    grm_name=$(object_name "$sparse_grm")
    grm_ids_name=$(object_name "$sparse_grm_ids")
    output_prefix="${trait}.saige_step1"

    analysis_command="step1_fitNULLGLMM.R \
      --sparseGRMFile='${grm_name}' \
      --sparseGRMSampleIDFile='${grm_ids_name}' \
      --useSparseGRMtoFitNULL=TRUE \
      --plinkFile='${plink_prefix}' \
      --phenoFile='${phenotype_name}' \
      --traitType=quantitative --invNormalize=TRUE \
      --phenoCol='${phenotype_column}' \
      --covarColList='${covariates}' \
      --qCovarColList='Sex,AssessmentCenter' \
      --sampleIDColinphenoFile=sample_id \
      --isCateVarianceRatio=TRUE \
      --IsOverwriteVarianceRatioFile=TRUE \
      --outputPrefix='${output_prefix}'"

    dx run swiss-army-knife \
      -iin="$phenotype_file" -iin="$bed" -iin="$bim" -iin="$fam" \
      -iin="$sparse_grm" -iin="$sparse_grm_ids" \
      -iimage_file="$image_file" -icmd="$analysis_command" \
      --name="SAIGE-GENE+ Step 1: ${trait}" \
      --destination="${destination}/step1/${trait}/" --brief -y
  done
elif [[ "$stage" == "step2" ]]; then
  annotations="LoF,LoF:DMis_Revel75:DMis_Revel50"
  maf_thresholds="0.00001,0.0001,0.001,0.01"

  tail -n +2 "$manifest" | while IFS=$'\t' read -r trait chrom model_file variance_ratio_file bgen bgi sample group_file sparse_grm sparse_grm_ids; do
    [[ -n "$trait" ]] || continue
    model_name=$(object_name "$model_file")
    variance_name=$(object_name "$variance_ratio_file")
    bgen_name=$(object_name "$bgen")
    bgi_name=$(object_name "$bgi")
    sample_name=$(object_name "$sample")
    group_name=$(object_name "$group_file")
    grm_name=$(object_name "$sparse_grm")
    grm_ids_name=$(object_name "$sparse_grm_ids")
    output_file="${trait}.chr${chrom}.saige_gene.tsv"

    analysis_command="step2_SPAtests.R \
      --bgenFile='${bgen_name}' --bgenFileIndex='${bgi_name}' --sampleFile='${sample_name}' \
      --AlleleOrder=ref-first --minMAF=0 --minMAC=1 --is_imputed_data=FALSE --maxMissing=0.05 \
      --LOCO=FALSE --sparseGRMFile='${grm_name}' --sparseGRMSampleIDFile='${grm_ids_name}' \
      --GMMATmodelFile='${model_name}' --varianceRatioFile='${variance_name}' \
      --groupFile='${group_name}' --annotation_in_groupTest='${annotations}' \
      --maxMAF_in_groupTest='${maf_thresholds}' \
      --is_output_markerList_in_groupTest=TRUE --is_output_moreDetails=TRUE \
      --IsOutputBETASEinBurdenTest=TRUE --is_overwrite_output=TRUE \
      --SAIGEOutputFile='${output_file}'"

    dx run swiss-army-knife \
      -iin="$model_file" -iin="$variance_ratio_file" \
      -iin="$bgen" -iin="$bgi" -iin="$sample" -iin="$group_file" \
      -iin="$sparse_grm" -iin="$sparse_grm_ids" \
      -iimage_file="$image_file" -icmd="$analysis_command" \
      --name="SAIGE-GENE+ Step 2: ${trait} chr${chrom}" \
      --destination="${destination}/step2/${trait}/" --brief -y
  done
else
  usage
  exit 2
fi

