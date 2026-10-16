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

# nf-schema requires Java 17+; gcc12-env/nextflow modules pull in OpenJDK 16
export JAVA_HOME=/zfshome/sw/spack/spack0.18.1/usr/opt/spack/linux-rocky8-x86_64/gcc-12.1.0/openjdk-21.0.3_9-yx33rojs7zmwjdm4guhkcihsqjfglque
export PATH="${JAVA_HOME}/bin:${PATH}"

nextflow run nf-core/ampliseq \
  -r 2.16.1 \
  -profile ccga_cau \
  -c ampliseq_custom_18s.config \
  -params-file params_18s_gmbc.yml \
  -resume