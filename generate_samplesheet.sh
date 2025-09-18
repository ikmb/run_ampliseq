#!/bin/bash

show_help() {
  cat << EOF
Usage: ${0##*/} --input_dir DIR [--output FILE]

Generates a samplesheet.tsv file for use with nf-core/ampliseq.

Options:
  --input_dir DIR     Path to the folder containing the raw FASTQ files (required)
  --output FILE       Name of the output samplesheet file (default: samplesheet.tsv)
  -h, --help          Show this help message and exit

Notes:
  - sampleID will be prefixed with a 'b' and '-' replaced with '_'
  - run ID is inferred from the second underscore-separated part before a dash

Example:
  /work_ikmb/ikmb_repository/shared/microbiome/RUN_AMPLISEQ/generate_samplesheet.sh --input_dir /path/to/fastq --output samplesheet.tsv
EOF
}

# Defaults
INPUT_DIR=""
OUTPUT="samplesheet.tsv"

# Parse arguments
while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --input_dir) INPUT_DIR="$2"; shift 2 ;;
    --output) OUTPUT="$2"; shift 2 ;;
    -h|--help) show_help; exit 0 ;;
    *) echo "Unknown option: $1" >&2; show_help; exit 1 ;;
  esac
done

# Check required input
if [[ -z "$INPUT_DIR" ]]; then
  echo "Error: --input_dir is required" >&2
  show_help
  exit 1
fi

# Check directory
if [ ! -d "$INPUT_DIR" ]; then
  echo "Error: Directory '$INPUT_DIR' not found." >&2
  exit 1
fi

# Create samplesheet
echo -e "sampleID\tforwardReads\treverseReads\trun" > "$OUTPUT"

for forward_file in "$INPUT_DIR"/*_R1_001.fastq.gz; do
  [ -e "$forward_file" ] || continue
  base_name=$(basename "$forward_file" "_R1_001.fastq.gz")
  reverse_file="${INPUT_DIR}/${base_name}_R2_001.fastq.gz"
  run=$(echo "$base_name" | sed -E 's/.*_([^-]+)-.*/\1/')
  sample_id="X$(echo "$base_name" | grep -oE '[0-9]{2}[A-Za-z]{3}[0-9]{3}-DL[0-9]{3}' | sed 's/-/_/')"
  echo -e "${sample_id}\t${forward_file}\t${reverse_file}\t${run}" >> "$OUTPUT"
done

echo "Samplesheet generated: $OUTPUT"