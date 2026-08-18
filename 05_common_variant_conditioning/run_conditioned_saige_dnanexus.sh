#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "Usage: $0 MANIFEST.tsv SAIGE_IMAGE_FILE DX_DESTINATION" >&2
  echo "Manifest columns: trait gene lead_snp chrom phenotype_file group_file bed bim fam bgen bgi sample sparse_grm sparse_grm_ids" >&2
  exit 2
fi

manifest=$1
image_file=$2
destination=${3%/}
command -v dx >/dev/null 2>&1 || { echo "The DNAnexus CLI (dx) is required." >&2; exit 1; }

object_name() {
  dx describe "$1" --name
}

covariates="Sex,Age,Age2,AssessmentCenter,TDI,LeadSNP_Dosage"
for pc in $(seq 1 20); do
  covariates+=",PC${pc}"
done

tail -n +2 "$manifest" | while IFS=$'\t' read -r trait gene lead_snp chrom phenotype_file group_file bed bim fam bgen bgi sample sparse_grm sparse_grm_ids; do
  [[ -n "$trait" ]] || continue
  safe_lead_snp=${lead_snp//[^[:alnum:]_.-]/_}
  phenotype_name=$(object_name "$phenotype_file")
  group_name=$(object_name "$group_file")
  bed_name=$(object_name "$bed")
  plink_prefix=${bed_name%.bed}
  bgen_name=$(object_name "$bgen")
  bgi_name=$(object_name "$bgi")
  sample_name=$(object_name "$sample")
  grm_name=$(object_name "$sparse_grm")
  grm_ids_name=$(object_name "$sparse_grm_ids")
  prefix="${trait}.${gene}.${safe_lead_snp}.conditioned"

  analysis_command="step1_fitNULLGLMM.R \
    --sparseGRMFile='${grm_name}' --sparseGRMSampleIDFile='${grm_ids_name}' \
    --useSparseGRMtoFitNULL=TRUE --plinkFile='${plink_prefix}' \
    --phenoFile='${phenotype_name}' --traitType=quantitative --invNormalize=TRUE \
    --phenoCol='${trait}' --covarColList='${covariates}' \
    --qCovarColList='Sex,AssessmentCenter' --sampleIDColinphenoFile=sample_id \
    --isCateVarianceRatio=TRUE --IsOverwriteVarianceRatioFile=TRUE --outputPrefix='${prefix}' && \
    step2_SPAtests.R --bgenFile='${bgen_name}' --bgenFileIndex='${bgi_name}' \
    --sampleFile='${sample_name}' --AlleleOrder=ref-first --minMAF=0 --minMAC=1 \
    --is_imputed_data=FALSE --maxMissing=0.05 --LOCO=FALSE \
    --sparseGRMFile='${grm_name}' --sparseGRMSampleIDFile='${grm_ids_name}' \
    --GMMATmodelFile='${prefix}.rda' --varianceRatioFile='${prefix}.varianceRatio.txt' \
    --groupFile='${group_name}' --annotation_in_groupTest='LoF,LoF:DMis_Revel75:DMis_Revel50' \
    --maxMAF_in_groupTest='0.00001,0.0001,0.001,0.01' \
    --IsOutputBETASEinBurdenTest=TRUE --is_output_moreDetails=TRUE \
    --SAIGEOutputFile='${prefix}.saige_gene.tsv'"

  dx run swiss-army-knife \
    -iin="$phenotype_file" -iin="$group_file" \
    -iin="$bed" -iin="$bim" -iin="$fam" \
    -iin="$bgen" -iin="$bgi" -iin="$sample" \
    -iin="$sparse_grm" -iin="$sparse_grm_ids" \
    -iimage_file="$image_file" -icmd="$analysis_command" \
    --name="Conditioned SAIGE-GENE+: ${trait} ${gene} ${lead_snp}" \
    --destination="${destination}/${trait}/${gene}/${safe_lead_snp}/" --brief -y
done
