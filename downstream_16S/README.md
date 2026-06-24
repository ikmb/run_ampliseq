# Downstream 16S Analysis Pipeline

A modular collection of R scripts for comprehensive downstream analysis of 16S amplicon sequencing data. Designed to be student-friendly with clean, minimal code.

## Overview

This pipeline provides a standardized workflow for:

- **Alpha diversity** analysis with statistical testing (Wilcoxon/Kruskal-Wallis/ANOVA)
- **Beta diversity** analysis with ordination (PCoA/NMDS), PERMANOVA, and convex hulls
- **Relative abundance** visualization at different taxonomic levels with faceting support
- **Differential abundance analysis** at ASV level using MaAslin2

All analyses are configurable through a single configuration file and can be run together via the master script or individually.

## Key Features

- 🎨 **Custom color palettes** for publication-ready plots
- 📊 **Statistical captions** automatically added to plots
- 🔬 **ASV-level differential abundance** with abundance filtering
- 🎯 **Convex hulls and centroids** in beta diversity plots
- 📑 **Faceting support** for stratified relative abundance plots
- ✨ **Clean, minimal code** designed for students
- 📦 **Example phyloseq preparation** script included

## Directory Structure

```
downstream_16S/
├── 00_master_script.R                      # Main script to run all analyses
├── 01_alpha_diversity.R                    # Alpha diversity analysis
├── 02_beta_diversity.R                     # Beta diversity analysis
├── 03_relative_abundance.R                 # Relative abundance plots
├── 04_differential_abundance_maaslin2.R    # DAA with MaAslin2 at ASV level
├── prepare_phyloseq.R                      # Example: prepare phyloseq from nf-core/ampliseq
├── config_template.R                       # Configuration template
├── config.R                                # Your configuration (create this)
└── README.md                               # This file
```

## Quick Start

### 1. Install Required Packages

```r
# Install from CRAN
install.packages(c("dplyr", "tidyr", "ggplot2", "vegan",
                   "ape", "FSA", "scales"))

# Install from Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("phyloseq", "microbiome", "Maaslin2", "genefilter"))

# Optional: Install picante for Faith's Phylogenetic Diversity
install.packages("picante")
```

### 2. Prepare Configuration File

```bash
# Copy the template
cp config_template.R config.R

# Edit config.R with your data paths and analysis parameters
```

### 3. Configure Your Analysis

Edit `config.R` to specify:

```r
# Input files
INPUT_ASV <- "path/to/asv_table.tsv"
INPUT_TAXONOMY <- "path/to/taxonomy.tsv"
INPUT_METADATA <- "path/to/metadata.csv"

# Or use existing phyloseq object
INPUT_PHYLOSEQ <- "path/to/phyloseq.rds"

# Output directory
OUTPUT_DIR <- "results"

# Analysis parameters
ALPHA_GROUP_VAR <- "treatment"
ALPHA_COVARIATES <- c("age", "sex")  # NULL for no covariates
BETA_COLOR_VAR <- "treatment"
BETA_PERMANOVA_VARS <- c("treatment", "age", "sex")
RELAB_TAX_LEVEL <- "Genus"
RELAB_FACET_VAR <- NULL  # or "timepoint" for faceting
DAA_FIXED_EFFECTS <- c("treatment")

# Custom color palette (Paul Tol's muted palette)
COLOR_PALETTE <- c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933",
                   "#DDCC77", "#CC6677", "#882255", "#AA4499", "#661100")
```

### 4. Run the Analysis

#### Option A: Run all analyses together

```r
source("00_master_script.R")
```

This will:

1. Load your data
2. Perform preprocessing and filtering
3. Run all downstream analyses
4. Generate plots and tables
5. Create a summary report

#### Option B: Run analyses individually

```r
# First, prepare the data
source("00_master_script.R")  # This creates ps_filtered

# Then run individual analyses
source("01_alpha_diversity.R")
alpha_results <- run_alpha_diversity(ps_filtered, config, OUTPUT_DIR)

source("02_beta_diversity.R")
beta_results <- run_beta_diversity(ps_filtered, config, OUTPUT_DIR)

source("03_relative_abundance.R")
relab_results <- run_relative_abundance(ps_filtered, config, OUTPUT_DIR)

source("04_differential_abundance_maaslin2.R")
maaslin2_results <- run_maaslin2_daa(ps_filtered, config, OUTPUT_DIR)
```

## Preparing Phyloseq Object

The `prepare_phyloseq.R` script provides an example of creating a phyloseq object from nf-core/ampliseq output:

```r
# Example structure
source("prepare_phyloseq.R")

# Input files from nf-core/ampliseq
feature_table <- "../ampliseq_output/dada2/feature-table.tsv"
taxonomy_file <- "../ampliseq_output/dada2/ASV_tax.gtdb.tsv"
species_taxonomy <- "../ampliseq_output/dada2/ASV_tax_species.gtdb.tsv"
tree_file <- "../ampliseq_output/dada2/tree.nwk"
metadata_file <- "metadata.csv"

# Creates phyloseq object with:
# - ASV intersection (filtered feature table matched with taxonomy)
# - ASV renaming (ASV_1, ASV_2, ...)
# - Metadata integration with type handling
# - Phylogenetic tree with matched tip labels
```

## Input Data Format

### ASV/Feature Table (TSV)

```
ASV_ID    Sample1    Sample2    Sample3    ...
ASV_001   100        150        80         ...
ASV_002   50         30         120        ...
```

### Taxonomy Table (TSV)

```
ASV_ID      Kingdom    Phylum         Class          Order          Family         Genus          Species_exact    confidence    sequence
ASV_001     Bacteria   Firmicutes     Bacilli        Lactobacillales Lactobacillaceae Lactobacillus  L. acidophilus   0.95         ATCG...
ASV_002     Bacteria   Bacteroidetes  Bacteroidia    Bacteroidales  Bacteroidaceae  Bacteroides    B. fragilis      0.89         GCTA...
```

Required columns: `ASV_ID`, `Kingdom`, `Phylum`, `Class`, `Order`, `Family`, `Genus`, `Species_exact`, `confidence`, `sequence`

### Metadata Table (CSV or TSV)

```
SampleID,treatment,age,sex,timepoint
Sample1,control,25,F,baseline
Sample2,treated,30,M,baseline
Sample3,control,28,F,week4
```

First column should be sample IDs (or include a `SampleID` column)

## Output Structure

Results are organized by analysis type:

```
results/
├── analysis_summary.txt                     # Overall summary report
├── phyloseq_raw.rds                        # Raw phyloseq object
├── phyloseq_filtered.rds                   # Preprocessed phyloseq object
├── alpha_diversity/
│   ├── alpha_diversity_indices.csv         # All indices for all samples
│   ├── alpha_diversity_all_indices.png     # Comparison of all indices
│   ├── alpha_diversity_shannon.png         # Plot with statistical caption
│   └── alpha_diversity_shannon_stats.txt   # Statistical test results
├── beta_diversity/
│   ├── beta_diversity_plot.png             # Ordination with hulls/centroids
│   ├── permanova_results.csv               # PERMANOVA table (marginal tests)
│   ├── distance_matrix_bray.rds            # Distance matrix
│   ├── ordination_PCoA.rds                 # Ordination object
│   └── beta_diversity_summary.txt          # Summary statistics
├── relative_abundance/
│   ├── relative_abundance_Genus.csv        # Long format table
│   ├── relative_abundance_Genus_wide.csv   # Wide format table
│   └── relative_abundance_Genus.png        # Stacked bar plot (with facets if configured)
└── differential_abundance_maaslin2/
    ├── maaslin2_significant_results.csv    # Significant ASVs (q < 0.25)
    ├── maaslin2_all_results.csv            # All tested ASVs
    ├── maaslin2_summary.txt                # Summary text
    └── maaslin2_output/                    # Full MaAslin2 output
```

## Configuration Options

### Preprocessing

- `PERFORM_RAREFACTION`: Should rarefaction be performed?
- `RAREFACTION_DEPTH`: Rarefaction depth (NULL for automatic)
- `MIN_COUNT_PER_SAMPLE`: Minimum count threshold
- `MIN_SAMPLES`: Minimum sample prevalence
- `MIN_TAXONOMY_CONFIDENCE`: Confidence threshold (0-1 scale, per-taxonomic-level)

### Alpha Diversity

- `ALPHA_INDICES`: Which indices to calculate (observed, shannon, gini_simpson, pielou, fisher, coverage, faith_pd)
- `ALPHA_GROUP_VAR`: Grouping variable for comparisons
- `ALPHA_METRIC`: Main metric for detailed analysis
- `ALPHA_COVARIATES`: Additional covariates for ANOVA (NULL for Wilcoxon/Kruskal-Wallis)
- `ALPHA_REFERENCE_GROUP`: Reference group for pairwise tests

### Beta Diversity

- `BETA_ORD_METHOD`: "PCoA" or "NMDS"
- `BETA_DISTANCE`: Distance metric ("bray", "jaccard", "unifrac", "wunifrac")
- `BETA_PERMANOVA_VARS`: Variables for marginal PERMANOVA (by="margin")
- `BETA_COLOR_VAR`: Variable for coloring ordination
- `BETA_SHAPE_VAR`: Variable for point shapes (or NULL)
- `BETA_VIS_STYLE`: Visualization style (hulls and centroids automatically added)

### Relative Abundance

- `RELAB_TAX_LEVEL`: Taxonomic level ("Phylum", "Genus", etc.)
- `RELAB_GROUP_VAR`: Grouping variable (use "SampleID" for individual samples)
- `RELAB_FACET_VAR`: Variable for faceting (NULL for no facets)
- `RELAB_THRESHOLD`: Abundance threshold for display (0.01 = 1%)
- `RELAB_SORT_TAXON`: Optional taxon to sort by
- `MIN_TAXONOMY_CONFIDENCE`: Applied before agglomeration

### Differential Abundance (MaAslin2)

- `DAA_FIXED_EFFECTS`: Main variables of interest
- `DAA_RANDOM_EFFECTS`: Random effects (for repeated measures)
- `DAA_COVARIATES`: Additional covariates to adjust for
- `DAA_REFERENCE`: Reference levels (e.g., `c(treatment = "control")`)
- `DAA_NORMALIZATION`: Normalization method ("TSS", "CLR", "CSS", "NONE", "TMM")
- `DAA_ANALYSIS_METHOD`: Statistical method ("LM", "CPLM", "ZICP", "NEGBIN", "ZINB")
- `DAA_MIN_ABUNDANCE`: Minimum reads per ASV (e.g., 4)
- `DAA_MIN_PREVALENCE`: Minimum fraction of samples (e.g., 0.1 = 10%)
- `DAA_MAX_SIGNIFICANCE`: Q-value threshold (e.g., 0.25)

### Plotting

- `PLOT_WIDTH`, `PLOT_HEIGHT`, `PLOT_DPI`: Plot dimensions
- `COLOR_PALETTE`: Vector of hex colors for custom palette

## Statistical Methods

### Alpha Diversity

- **Two groups (no covariates)**: Wilcoxon rank-sum test
- **Multiple groups (no covariates)**: Kruskal-Wallis test + Dunn's post-hoc (if significant)
- **With covariates**: Type II ANOVA (linear model)
- Statistical results shown in plot captions

### Beta Diversity

- **Ordination**: PCoA or NMDS on distance matrices (Bray-Curtis, UniFrac, etc.)
- **Testing**: PERMANOVA with marginal tests (`adonis2(..., by="margin")`)
  - Marginal tests show unique contribution of each variable after controlling for others
- **Visualization**: Convex hulls for group boundaries, centroids with triangles
- PERMANOVA results shown in plot captions with R² and p-values

### Relative Abundance

- **Abundance filtering**: ≥4 reads in ≥10% of samples (configurable)
- **Confidence filtering**: Applied before taxonomic agglomeration
- **Unknown taxa**: Low-confidence taxa combined and placed at bottom of bars
- **Faceting**: Optional stratification by additional variable

### Differential Abundance (MaAslin2)

- **Analysis level**: ASV-level (not genus-level agglomeration)
- **Abundance filtering**: ≥4 reads in ≥10% of samples (default)
- **Normalization**: Total Sum Scaling (TSS) by default
- **Transformation**: LOG transformation
- **Method**: Linear Models (LM)
- **Multiple testing**: Benjamini-Hochberg FDR correction
- **Significance**: q-value < 0.25 (default)
- **ASV labeling**: "ASV_X_GenusName" format for interpretability

## Tips and Best Practices

### Rarefaction

- Rarefaction is **optional** but recommended for alpha/beta diversity
- For differential abundance, use raw counts (rarefaction not needed)
- Check sample depth distribution before choosing rarefaction depth

### Filtering

- Remove low-abundance/low-prevalence taxa to reduce noise
- Typical thresholds: min 5 counts in at least 2-5 samples
- Adjust based on your dataset size

### Statistical Power

- For differential abundance, ensure adequate sample size (n ≥ 10 per group)
- Consider batch effects and confounders
- Use random effects for repeated measures or paired designs

### Visualization

- Adjust plot dimensions and DPI for publication quality
- Custom color palette using Paul Tol's muted qualitative scheme (color-blind friendly)
- Statistical captions automatically added to plots
- Convex hulls and centroids enhance beta diversity visualization
- Use faceting (`RELAB_FACET_VAR`) to stratify relative abundance by additional variables
- Unknown taxa automatically placed at bottom of stacked bars
- Customize further by modifying the plotting functions

## Troubleshooting

### "Package not found" error

Install missing packages:

```r
BiocManager::install("package_name")
```

### "Phyloseq object not found" error

Run the master script first to create `ps_filtered`:

```r
source("00_master_script.R")
```

### Memory issues with large datasets

- Increase R memory limit with `memory.limit(size = 16000)` (Windows)
- Run analyses on filtered/rarefied data to reduce dataset size
- Process samples in batches if needed

### No significant results in DAA

- Check if your effect size is realistic
- Increase sample size if possible
- Adjust prevalence/abundance thresholds
- Consider using different normalization methods

## Citation

If you use this pipeline, please cite the relevant tools:

- **phyloseq**: McMurdie PJ, Holmes S (2013). phyloseq: An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data. PLoS ONE 8(4): e61217.
- **MaAslin2**: Mallick H, et al. (2021). Multivariable association discovery in population-scale meta-omics studies. PLoS Computational Biology 17(11): e1009442.
- **vegan**: Oksanen J, et al. (2020). vegan: Community Ecology Package. R package version 2.5-7.
- **microbiome**: Lahti L, Shetty S (2017). microbiome R package.

## Contact

For questions, issues, or contributions, please contact the pipeline maintainer.

## License

This pipeline is provided as-is for research purposes.
