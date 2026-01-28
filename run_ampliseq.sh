#!/bin/bash
#SBATCH --output=ampliseq_pipeline_%A.out  # Standard output file
#SBATCH --error=ampliseq_pipeline_%A.err   # Standard error file
#SBATCH --cpus-per-task=16                 # Number of CPUs
#SBATCH --mem=128GB                        # Memory allocation
#SBATCH --time=48:00:00                    # Maximum runtime

# Load required modules
module load gcc12-env
module load nextflow
module load singularity

show_help() {
  cat << EOF
Usage: ${0##*/} --input <samplesheet.tsv> --outdir <results_dir> [options]

Required arguments:
  --input FILE                    Path to the nf-core/ampliseq samplesheet
  --outdir DIR                    Output directory for results

Optional arguments:
  --primers v3v4|v1v2|its2|archaea|18s                           Primer profile (default: v3v4)
  --db silva|gtdb|pr2|greengenes2|rdp|unite-fungi=10.0           Taxonomy database to use (default: rdp)
  --single_end                                                   Process only forward reads (18S only)
  --multiple_sequencing_runs                                     Enables handling of multiple sequencing runs
  -resume                                                        Resume a previously failed/interrupted run
  -h, --help                                                     Show this help message and exit

Example:
  ./run_ampliseq.sh --input samplesheet.tsv --outdir results --primers its2 --db unite --multiple_sequencing_runs -resume
  ./run_ampliseq.sh --input samplesheet.tsv --outdir results --primers 18s --db pr2 --single_end
EOF
}

# ---------------------
# Default values
# ---------------------
DB="rdp"
PRIMER_PROFILE="v3v4"
SINGLE_END=false
NEXTFLOW_ARGS=()

# ---------------------
# Parse named arguments
# ---------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --input) INPUT_FILE="$2"; shift 2 ;;
    --outdir) OUTPUT_DIR="$2"; shift 2 ;;
    --db) DB="$2"; shift 2 ;;
    --primers) PRIMER_PROFILE="$2"; shift 2 ;;
    --single_end) SINGLE_END=true; shift ;;
    -resume|--multiple_sequencing_runs) NEXTFLOW_ARGS+=("$1"); shift ;;
    -h|--help) show_help; exit 0 ;;
    *) echo "Unknown argument: $1"; show_help; exit 1 ;;
  esac
done

# Check for required named arguments
if [[ -z "$INPUT_FILE" || -z "$OUTPUT_DIR" ]]; then
  echo "ERROR: --input and --outdir are required."
  echo
  show_help
  exit 1
fi

# ---------------------
# Set primer sequences
# ---------------------
if [[ "$PRIMER_PROFILE" != "its2" ]]; then
  case "$PRIMER_PROFILE" in
    v3v4)
      FW_PRIMER="CCTACGGGAGGCAGCAG"
      RV_PRIMER="GGACTACHVGGGTWTCTAAT"
      ;;
    archaea)
      FW_PRIMER="CAGCMGCCGCGGTAA"
      RV_PRIMER="GGACTACVSGGGTATCTAAT"
      ;;
    18s)
      FW_PRIMER="TTAAARVGYTCGTAGTYG"
      RV_PRIMER="CCGTCAATTHCTTYAART"
      ;;
    *)
      echo "ERROR: Unknown --primers option: $PRIMER_PROFILE"
      echo "       Valid options: v3v4, v1v2, its2, archaea, 18s"
      exit 1
      ;;
  esac
fi

# ---------------------
# Run nf-core/ampliseq
# ---------------------
CMD=(
  nextflow run nf-core/ampliseq
  -r 2.14.0
  -profile ccga_cau
  --input "$INPUT_FILE"
  --ignore_failed_filtering
  --dada_ref_taxonomy "$DB"
  --dada_assign_taxlevels Kingdom,Phylum,Class,Order,Family,Genus,Species
  --outdir "$OUTPUT_DIR"
  --skip_qiime
  --skip_phyloseq
  --skip_tse
  --skip_qiime_downstream
  --skip_fastqc
  --skip_multiqc
  --skip_report
  --ignore_empty_input_files
  --ignore_failed_filtering
)

# ITS2-specific adjustments
if [[ "$PRIMER_PROFILE" == "its2" ]]; then
  CMD+=(
    --illumina_pe_its
    --skip_cutadapt
    --trunclenf 230
    --trunclenr 150
    --max_ee 4
    --skip_dada_addspecies
    -c /work_ikmb/ikmb_repository/shared/microbiome/RUN_AMPLISEQ/ampliseq_custom.config
  )
elif [[ "$PRIMER_PROFILE" == "v1v2" ]]; then
  CMD+=(
    --skip_cutadapt
  )
elif [[ "$PRIMER_PROFILE" == "v3v4" ]]; then
  CMD+=(
    --trunclenf 260
    --trunclenr 200
    --max_ee 2
    --FW_primer "$FW_PRIMER"
    --RV_primer "$RV_PRIMER"
  )
elif [[ "$PRIMER_PROFILE" == "18s" ]]; then
  CMD+=(
    --skip_cutadapt
    --dada_assign_taxlevels Domain,Supergroup,Division,Subdivision,Class,Order,Family,Genus,Species
    --max_ee 2
    --skip_dada_addspecies
  )
  
  # Adjust truncation and config based on single-end vs paired-end
  if [[ "$SINGLE_END" == true ]]; then
    CMD+=(
      --single_end 
      --trunclenf 250
      -c /work_ikmb/sukmb662/test_16s/ampliseq_custom_18s.config
    )
  else
    CMD+=(
      --trunclenf 250 
      --trunclenr 210
      -c /work_ikmb/sukmb662/test_16s/ampliseq_custom_18s_no_tryRC.config
    )
  fi
else
  CMD+=(
    --FW_primer "$FW_PRIMER"
    --RV_primer "$RV_PRIMER"
  )
fi

# Append any user-defined flags like -resume
CMD+=("${NEXTFLOW_ARGS[@]}")

# Execute command
"${CMD[@]}"