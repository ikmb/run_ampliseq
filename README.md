**BIP Platform wrapper scripts for standardized execution of [nf-core/ampliseq](https://nf-co.re/ampliseq)**  
This repository provides helper scripts and usage instructions to simplify processing of raw amplicon data into cleaned, standardized ASV tables.  
It ensures consistent sample naming, reproducible parameterization, and automated cleanup of intermediate files.

---

## 📂 Scripts included

| Script | Purpose |
|--------|---------|
| `generate_samplesheet.sh` | Create a standardized `samplesheet.tsv` from raw FASTQ filenames |
| `run_ampliseq.sh`         | Launch the `nf-core/ampliseq` pipeline with chosen primers and database |
| `prepare_seqtab.sh`       | Clean and rename the final DADA2 output, remove intermediate files |

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

###  2. Generate Samplesheet

Use the helper script to create a standardized `samplesheet.tsv` based on raw FASTQ filenames:

```bash
/work_ikmb/ikmb_repository/shared/microbiome/RUN_AMPLISEQ/generate_samplesheet.sh \
  --input_dir /path/to/raw_data
```

###  3. Run the Pipeline
We recommend running the pipeline inside a tmux session to avoid job interruptions:
```bash
tmux new-session

/work_ikmb/ikmb_repository/shared/microbiome/RUN_AMPLISEQ/run_ampliseq.sh \
  --input samplesheet.tsv \
  --outdir results \
  --primers its2 \
  --db unite
```

###  4. Post-process Output (optional)
After the pipeline finishes, clean the output and extract seqtab with:
```bash
cd /work_ikmb/sukmb276/Microbiome/clean_data_from_dada2/Runs_v.1.10_Fungi/XXMonthXXX

/work_ikmb/ikmb_repository/shared/microbiome/RUN_AMPLISEQ/prepare_seqtab.sh
```

## Other Use Cases

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
