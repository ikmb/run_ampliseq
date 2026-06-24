################################################################################
# Configuration File for 16S Downstream Analysis
# 
# This file contains all user-defined parameters for downstream 16S analysis.
# Copy this file to 'config.R' and modify the values according to your data.
################################################################################

# =============================================================================
# INPUT FILES
# =============================================================================

# Option 1: Provide individual files
INPUT_ASV <- "path/to/asv_table.tsv"           # ASV/feature table (TSV format)
INPUT_TAXONOMY <- "path/to/taxonomy.tsv"       # Taxonomy table (TSV format)
INPUT_METADATA <- "path/to/metadata.csv"       # Metadata/sample info (CSV or TSV)

# Option 2: Provide phyloseq object (if you already have one)
INPUT_PHYLOSEQ <- NULL                          # Path to phyloseq .rds file or NULL

# =============================================================================
# OUTPUT DIRECTORY
# =============================================================================

OUTPUT_DIR <- "results"                         # Directory for all output files

# =============================================================================
# PREPROCESSING OPTIONS
# =============================================================================

# Rarefaction settings
PERFORM_RAREFACTION <- TRUE                     # Should rarefaction be performed?
RAREFACTION_DEPTH <- NULL                       # NULL for automatic (minimum depth), or specify a number

# Filtering settings
MIN_COUNT_PER_SAMPLE <- 1                       # Minimum count per sample for a taxon
MIN_SAMPLES <- 5                                # Minimum number of samples a taxon must appear in
MIN_TAXONOMY_CONFIDENCE <- 70                   # Minimum confidence threshold (0-100)

# =============================================================================
# ALPHA DIVERSITY SETTINGS
# =============================================================================

# Indices to calculate (available: "observed", "diversity_shannon", 
# "diversity_gini_simpson", "evenness_pielou", "dominance_core_abundance", 
# "diversity_fisher")
ALPHA_INDICES <- c("observed", "diversity_shannon", "diversity_gini_simpson",
                   "evenness_pielou", "diversity_fisher")

# Grouping variable(s) for statistical testing
ALPHA_GROUP_VAR <- "treatment"                  # Column name in metadata

# Metric for detailed analysis
ALPHA_METRIC <- "diversity_shannon"             # Main metric for plots and tests

# Covariates for multivariable model (leave as NULL for simple comparison)
ALPHA_COVARIATES <- NULL                        # c("age", "sex") or NULL

# Reference group (only used if ALPHA_COVARIATES is NULL)
ALPHA_REFERENCE_GROUP <- NULL                   # Specific group name or NULL for first level

# =============================================================================
# BETA DIVERSITY SETTINGS
# =============================================================================

# Ordination method
BETA_ORD_METHOD <- "PCoA"                       # "PCoA" or "NMDS"

# Distance metric
BETA_DISTANCE <- "bray"                         # "bray", "jaccard", "unifrac", etc.

# PERMANOVA formula variables
BETA_PERMANOVA_VARS <- c("treatment")           # Variables for PERMANOVA testing

# Visualization
BETA_COLOR_VAR <- "treatment"                   # Variable for coloring points
BETA_SHAPE_VAR <- NULL                          # Variable for point shapes (or NULL)
BETA_VIS_STYLE <- "Basic"                       # "Basic", "Centroids", or "Convex Hulls"

# =============================================================================
# RELATIVE ABUNDANCE SETTINGS
# =============================================================================

# Taxonomic level for agglomeration
RELAB_TAX_LEVEL <- "Genus"                      # "Phylum", "Class", "Order", "Family", "Genus", "Species"

# Grouping variable
RELAB_GROUP_VAR <- "treatment"                  # Column name in metadata for grouping

# Abundance threshold for "Other" category
RELAB_THRESHOLD <- 1                            # Taxa below this % grouped as "< 1%"

# Sorting (optional: specific taxon name to sort by, or NULL)
RELAB_SORT_TAXON <- NULL                        # Taxon name or NULL

# =============================================================================
# DIFFERENTIAL ABUNDANCE ANALYSIS (DAA) SETTINGS
# =============================================================================

# MaAslin2/3 settings
DAA_TAX_LEVEL <- "Genus"                        # Taxonomic level for DAA

# Fixed effects (required)
DAA_FIXED_EFFECTS <- c("treatment")             # Main variables of interest

# Reference level(s) for fixed effects (optional, named vector)
DAA_REFERENCE <- c(treatment = "control")       # e.g., c(treatment = "control", time = "baseline")

# Random effects (optional, for repeated measures or batch effects)
DAA_RANDOM_EFFECTS <- NULL                      # c("subject_id") or NULL

# Covariates to adjust for
DAA_COVARIATES <- NULL                          # c("age", "sex") or NULL

# Normalization method
DAA_NORMALIZATION <- "TSS"                      # "TSS", "CLR", "CSS", "NONE", "TMM"

# Analysis method (for MaAslin2)
DAA_ANALYSIS_METHOD <- "LM"                     # "LM", "CPLM", "ZICP", "NEGBIN", "ZINB"

# Statistical significance thresholds
DAA_MIN_ABUNDANCE <- 0.0                        # Minimum abundance threshold
DAA_MIN_PREVALENCE <- 0.1                       # Minimum prevalence (fraction of samples)
DAA_MAX_SIGNIFICANCE <- 0.05                    # Q-value threshold for significance

# =============================================================================
# PLOTTING OPTIONS
# =============================================================================

# Plot dimensions
PLOT_WIDTH <- 10                                # Width in inches
PLOT_HEIGHT <- 6                                # Height in inches
PLOT_DPI <- 300                                 # DPI for PNG exports

# Color palette (color-blind friendly)
COLOR_PALETTE <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                   "#0072B2", "#D55E00", "#CC79A7", "#999999")

# =============================================================================
# COMPUTATIONAL OPTIONS
# =============================================================================

# Random seed for reproducibility
RANDOM_SEED <- 42

# Number of cores for parallel processing (where applicable)
N_CORES <- 4

################################################################################
# END OF CONFIGURATION
################################################################################
