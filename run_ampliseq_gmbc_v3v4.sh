#!/usr/bin/env bash
#SBATCH --job-name=gmbc_v3v4_ampliseq
#SBATCH --output=gmbc_v3v4_ampliseq_%A.out
#SBATCH --error=gmbc_v3v4_ampliseq_%A.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=128GB
#SBATCH --time=48:00:00

set -euo pipefail

module load gcc12-env
module load nextflow
module load singularity

# V3-V4 primers
FW_PRIMER="CCTACGGGAGGCAGCAG"
RV_PRIMER="GGACTACHVGGGTWTCTAAT"

nextflow run nf-core/ampliseq \
  -r 2.16.1 \
  -profile ccga_cau \
  -c ampliseq_custom.config \
  -resume \
  --input samplesheet.tsv \
  --outdir results_filtered \
  --FW_primer "${FW_PRIMER}" \
  --RV_primer "${RV_PRIMER}" \
  --dada_assign_taxlevels Kingdom,Phylum,Class,Order,Family,Genus,Species \
  --ignore_failed_filtering \
  --trunclenf 260 \
  --trunclenr 200 \
  --max_ee 2 \
  --ignore_empty_input_files \
  --dada_ref_taxonomy gtdb \
  --min_samples 2 \
  --min_frequency 10 \
  --skip_diversity_indices \
  --metadata metadata.tsv