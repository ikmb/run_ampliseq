#!/usr/bin/env bash
#SBATCH --job-name=gmbc_18s_ampliseq
#SBATCH --output=gmbc_18s_ampliseq_%A.out
#SBATCH --error=gmbc_18s_ampliseq_%A.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=128GB
#SBATCH --time=48:00:00

set -euo pipefail

module load gcc12-env
module load nextflow
module load singularity

# 18S primers
FW_PRIMER="TTAAARVGYTCGTAGTYG"
RV_PRIMER="CCGTCAATTHCTTYAART"

nextflow run nf-core/ampliseq \
  -r 2.16.1 \
  -profile ccga_cau \
  -c ampliseq_custom_18s.config \
  -resume \
  --input samplesheet.tsv \
  --outdir results_filtered_18s \
  --FW_primer "${FW_PRIMER}" \
  --RV_primer "${RV_PRIMER}" \
  --skip_cutadapt \
  --dada_assign_taxlevels Domain,Supergroup,Division,Subdivision,Class,Order,Family,Genus,Species \
  --ignore_failed_filtering \
  --trunclenf 250 \
  --trunclenr 210 \
  --max_ee 2 \
  --skip_dada_addspecies \
  --ignore_empty_input_files \
  --db pr2 \
  --min_samples 2 \
  --min_frequency 10 \
  --skip_diversity_indices \
  --metadata metadata.tsv