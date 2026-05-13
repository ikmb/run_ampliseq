**BIP Platform wrapper scripts for standardized execution of [nf-core/ampliseq](https://nf-co.re/ampliseq)**  
This repository provides helper scripts and usage instructions to simplify processing of raw amplicon data into cleaned, standardized ASV tables.  
It ensures consistent sample naming, reproducible parameterization, and automated cleanup of intermediate files.

---

## 📂 Scripts included

| Script                      | Purpose                                                                          |
| --------------------------- | -------------------------------------------------------------------------------- |
| `generate_samplesheet.sh`   | Create a standardized `samplesheet.tsv` from raw FASTQ filenames                 |
| `run_ampliseq.sh`           | Launch the `nf-core/ampliseq` pipeline with chosen primers and database          |
| `run_ampliseq_gmbc_v3v4.sh` | Run GMBC V3-V4 settings with QIIME2 branch enabled for phylogenetic tree outputs |
| `run_ampliseq_gmbc_18s.sh`  | Run GMBC 18S settings with QIIME2 branch enabled for phylogenetic tree outputs   |
| `prepare_seqtab.sh`         | Clean and rename the final DADA2 output, remove intermediate files               |

These scripts are also available at:  
`/work_ikmb/ikmb_repository/shared/microbiome/RUN_AMPLISEQ/`

---

## How to Run

### 1. Prepare Working Directory

```bash
cd /work_ikmb/sukmb276/Microbiome/clean_data_from_dada2/Runs_v.1.10_Fungi
mkdir XXMonthXXX
cd XXMonthXXX
```

### 2. Generate Samplesheet

Use the helper script to create a standardized `samplesheet.tsv` based on raw FASTQ filenames:

```bash
/work_ikmb/ikmb_repository/shared/microbiome/RUN_AMPLISEQ/generate_samplesheet.sh \
  --input_dir /path/to/raw_data
```

### 3. Run the Pipeline

We recommend running the pipeline inside a tmux session to avoid job interruptions:

```bash
tmux new-session

# Example 1: ITS2 fungi (paired-end)
/work_ikmb/ikmb_repository/shared/microbiome/RUN_AMPLISEQ/run_ampliseq.sh \
  --input samplesheet.tsv \
  --outdir results \
  --primers its2 \
  --db unite

# Example 2: 18S eukaryotes (paired-end)
/work_ikmb/ikmb_repository/shared/microbiome/RUN_AMPLISEQ/run_ampliseq.sh \
  --input samplesheet.tsv \
  --outdir results \
  --primers 18s \
  --db pr2

# Example 3: 18S eukaryotes (single-end, forward only)
/work_ikmb/ikmb_repository/shared/microbiome/RUN_AMPLISEQ/run_ampliseq.sh \
  --input samplesheet.tsv \
  --outdir results \
  --primers 18s \
  --db pr2 \
  --single_end
```

### 3b. GMBC V3-V4 run (with QIIME2 branch and phylogenetic tree)

For GMBC V3-V4 processing, use the dedicated script in this repository:

```bash
bash run_ampliseq_gmbc_v3v4.sh
```

This script uses:

- `nf-core/ampliseq` `2.16.1`
- V3-V4 primers (`CCTACGGGAGGCAGCAG` / `GGACTACHVGGGTWTCTAAT`)
- `--db gtdb`
- `--outdir results_filtered`
- Filtering settings (`--min_samples 2`, `--min_frequency 10`)

Why this is needed for GMBC:

- The generic wrapper may skip downstream QIIME2 steps in some modes.
- The GMBC script keeps the QIIME2 branch active (no `--skip_qiime*` flags), which is required to produce phylogeny outputs.
- This is the expected path when you need phylogenetic tree generation for downstream analyses.

Input files expected in the working directory:

- `samplesheet.tsv`
- `metadata.tsv`
- `custom.config` (already provided in this repository and points to `ampliseq_custom.config`)

### 3c. GMBC 18S run (with QIIME2 branch and phylogenetic tree)

For GMBC 18S processing, use the dedicated script in this repository:

```bash
bash run_ampliseq_gmbc_18s.sh
```

This script uses:

- `nf-core/ampliseq` `2.16.1`
- 18S primers (`TTAAARVGYTCGTAGTYG` / `CCGTCAATTHCTTYAART`)
- `--db pr2`
- `--outdir results_filtered_18s`
- Filtering settings (`--min_samples 2`, `--min_frequency 10`)

Why this is needed for GMBC:

- It mirrors the GMBC V3-V4 dedicated execution model for consistency.
- It keeps the QIIME2 branch active (no `--skip_qiime*` flags), which is required for phylogeny outputs.
- It applies 18S-specific taxonomy levels and uses the 18S custom config.

Input files expected in the working directory:

- `samplesheet.tsv`
- `metadata.tsv`

### 4. Post-process Output (optional)

After the pipeline finishes, clean the output and extract seqtab with:

```bash
cd /work_ikmb/sukmb276/Microbiome/clean_data_from_dada2/Runs_v.1.10_Fungi/XXMonthXXX

/work_ikmb/ikmb_repository/shared/microbiome/RUN_AMPLISEQ/prepare_seqtab.sh
```

## Other Use Cases

### Single-end mode (forward reads only)

For 18S rRNA analysis with single-end reads:

```bash
/work_ikmb/ikmb_repository/shared/microbiome/RUN_AMPLISEQ/run_ampliseq.sh \
  --input samplesheet.tsv \
  --outdir results \
  --primers 18s \
  --db pr2 \
  --single_end
```

**Note:** For 18S analysis:

- **Paired-end mode**: Uses `ampliseq_custom_18s.config` (with `tryRC=TRUE` for proper taxonomic assignment)
- **Single-end mode** (default): Uses `ampliseq_custom_18s_no_tryRC.config` (without `tryRC=TRUE`)

### Resume after failure

```bash
/work_ikmb/ikmb_repository/shared/microbiome/RUN_AMPLISEQ/run_ampliseq.sh \
  --input samplesheet.tsv \
  --outdir results \
  --primers its2 \
  --db unite \
  --resume
```

### Multiple sequencing runs

```bash
/work_ikmb/ikmb_repository/shared/microbiome/RUN_AMPLISEQ/run_ampliseq.sh \
  --input samplesheet.tsv \
  --outdir results \
  --primers its2 \
  --db unite \
  --multiple_sequencing_runs
```

## Help

To explore script options:

```bash
/work_ikmb/ikmb_repository/shared/microbiome/RUN_AMPLISEQ/run_ampliseq.sh --help
/work_ikmb/ikmb_repository/shared/microbiome/RUN_AMPLISEQ/generate_samplesheet.sh --help
/work_ikmb/ikmb_repository/shared/microbiome/RUN_AMPLISEQ/prepare_seqtab.sh --help
```
