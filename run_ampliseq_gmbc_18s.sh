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

nextflow run nf-core/ampliseq \
  -r 2.16.1 \
  -profile ccga_cau \
  -c ampliseq_custom_18s.config \
  -params-file params_18s.yml \
  -resume