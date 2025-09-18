#!/bin/bash

set -euo pipefail

# Function to show help message
show_help() {
  cat << EOF
Usage: ./prepare_seqtab.sh

This script should be executed inside the base directory of a DADA2 run,
where the 'work/' directory is located.

It performs the following:
  - Finds and copies the first *.seqtab.rds file to ./seqtab.Rds
  - Cleans the rownames of seqtab.Rds to standardized format (e.g. 25Jun798_DL123)
  - Saves the result as seqtab.cleaned.Rds
  - Removes the work/ directory to save space

Requirements:
  - R/4.3.3 module with stringr package installed
  - You must run this script from the directory that contains 'work/'

EOF
}

# Check for --help flag
if [[ "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then
  show_help
  exit 0
fi

echo "Searching and copying .seqtab.rds file..."
SEQTAB_FILE=$(find work/ -name "*.seqtab.rds" | head -n 1)

if [[ -z "$SEQTAB_FILE" ]]; then
  echo ".seqtab.rds file found in work/ directory. Aborting."
  exit 1
fi

cp "$SEQTAB_FILE" ./seqtab.Rds

# Load modules
module load gcc12-env
module load R/4.3.3

echo "Cleaning sample names in seqtab.Rds..."

module load gcc12-env
module load R/4.3.3

Rscript --vanilla -e '
seqtable <- readRDS("seqtab.Rds")
rownames(seqtable) <- sub(".*?(\\d{2}[A-Za-z]{3}\\d{2,3}-DL\\d+).*", "\\1", rownames(seqtable))
saveRDS(seqtable, "seqtab.Rds")
cat("Done. Example rownames:\n")
print(head(rownames(seqtable)))
'
# Clean up work directory
echo "Removing work/ directory..."
rm -rf work/
echo "Cleanup complete."