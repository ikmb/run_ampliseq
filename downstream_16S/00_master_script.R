################################################################################
# 00_master_script.R
# 
# Master script to run all downstream 16S analyses
# This script coordinates data loading, preprocessing, and all downstream analyses
################################################################################

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(speedyseq)
})

# =============================================================================
# LOAD CONFIGURATION
# =============================================================================

cat("\n=== 16S Downstream Analysis Pipeline ===\n\n")

if (!file.exists("config.R")) {
  stop("Configuration file 'config.R' not found. Please copy 'config_template.R' to 'config.R' and configure it.")
}

cat("Loading configuration...\n")
source("config.R")

# Set random seed for reproducibility
set.seed(RANDOM_SEED)

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
cat("Output directory:", OUTPUT_DIR, "\n")


# =============================================================================
# LOAD DATA
# =============================================================================

cat("\n=== Loading Data ===\n\n")

if (!is.null(INPUT_PHYLOSEQ) && file.exists(INPUT_PHYLOSEQ)) {
  cat("Loading phyloseq object from:", INPUT_PHYLOSEQ, "\n")
  ps <- readRDS(INPUT_PHYLOSEQ)
  cat("Loaded successfully\n")
  
} else {
  # Load individual files and construct phyloseq object
  message("Loading ASV table, taxonomy, and metadata...")
  
  # Load ASV table
  if (!file.exists(INPUT_ASV)) {
    stop("ASV table file not found: ", INPUT_ASV)
  }
  asv <- read.table(INPUT_ASV, header = TRUE, row.names = 1, 
                   sep = "\t", check.names = FALSE)
  message(paste0("  ASV table: ", nrow(asv), " ASVs, ", ncol(asv), " samples"))
  
  # Load taxonomy
  if (!file.exists(INPUT_TAXONOMY)) {
    stop("Taxonomy file not found: ", INPUT_TAXONOMY)
  }
  tax <- read.table(INPUT_TAXONOMY, header = TRUE, sep = "\t", 
                   stringsAsFactors = FALSE)
  
  # Set ASV_ID as row names
  if ("ASV_ID" %in% colnames(tax)) {
    rownames(tax) <- tax$ASV_ID
    tax$ASV_ID <- NULL
  }
  message(paste0("  Taxonomy table: ", nrow(tax), " entries"))
  
  # Load metadata
  if (!file.exists(INPUT_METADATA)) {
    stop("Metadata file not found: ", INPUT_METADATA)
  }
  
  if (grepl("\\.csv$", INPUT_METADATA)) {
    meta <- read.csv(INPUT_METADATA, header = TRUE, stringsAsFactors = FALSE)
  } else {
    meta <- read.table(INPUT_METADATA, header = TRUE, sep = "\t", 
                      stringsAsFactors = FALSE)
  }
  
  # Set sample IDs as row names (assuming first column or a column named "SampleID")
  if ("SampleID" %in% colnames(meta)) {
    rownames(meta) <- meta$SampleID
  } else {
    rownames(meta) <- meta[, 1]
    meta <- meta[, -1, drop = FALSE]
  }
  message(paste0("  Metadata: ", nrow(meta), " samples"))
  
  # Create phyloseq object
  message("Creating phyloseq object...")
  ps <- phyloseq(
    otu_table(as.matrix(asv), taxa_are_rows = TRUE),
    tax_table(as.matrix(tax)),
    sample_data(meta)
  )
  
  # Save phyloseq object
  ps_file <- file.path(OUTPUT_DIR, "phyloseq_raw.rds")
  saveRDS(ps, ps_file)
  message(paste0("Saved raw phyloseq object: ", ps_file))
}

# Print phyloseq summary
message("\nPhyloseq object summary:")
print(ps)

cat("\n=== Preprocessing ===\n\n")

# Check if phyloseq object has tree
if (!is.null(phy_tree(ps, errorIfNULL = FALSE))) {
  cat("Phylogenetic tree found:", ntaxa(phy_tree(ps)), "tips\n")
} else if (!is.null(TREE_FILE) && file.exists(TREE_FILE)) {
  cat("Loading phylogenetic tree from file...\n")
  library(ape)
  tree <- read.tree(TREE_FILE)
  phy_tree(ps) <- tree
  cat("Added tree to phyloseq object\n")
} else {
  cat("No phylogenetic tree available (Faith's PD and UniFrac will be skipped)\n")
}

# Filter low-abundance/low-prevalence taxa
cat("\nFiltering taxa...")
keep_taxa <- phyloseq::genefilter_sample(
  ps,
  filterfun_sample(function(x) x >= MIN_COUNT_PER_SAMPLE),
  A = MIN_SAMPLES
)
ps_filtered <- phyloseq::prune_taxa(keep_taxa, ps)
cat(" Kept", ntaxa(ps_filtered), "of", ntaxa(ps), "taxa\n")

# Rarefaction (if enabled)
if (PERFORM_RAREFACTION) {
  cat("\nRarefying to depth:", RAREFACTION_DEPTH, "...")
  
  # Check if any samples will be removed
  low_samples <- sample_sums(ps_filtered) < RAREFACTION_DEPTH
  if (any(low_samples)) {
    n_low <- sum(low_samples)
    cat("\n  Warning:", n_low, "sample(s) will be removed due to low depth\n")
  }
  
  # Rarefy
  ps_filtered <- phyloseq::rarefy_even_depth(
    ps_filtered, 
    sample.size = RAREFACTION_DEPTH,
    rngseed = RANDOM_SEED,
    replace = FALSE,
    trimOTUs = TRUE,
    verbose = TRUE
  )
  
  message(paste0("  After rarefaction: ", nsamples(ps_filtered), " samples, ", 
                ntaxa(ps_filtered), " taxa"))
}

# Save preprocessed phyloseq object
ps_filtered_file <- file.path(OUTPUT_DIR, "phyloseq_filtered.rds")
saveRDS(ps_filtered, ps_filtered_file)
message(paste0("\nSaved preprocessed phyloseq object: ", ps_filtered_file))


# =============================================================================
# RUN ANALYSES
# =============================================================================

# Create a config list to pass to analysis functions
config <- list(
  # Alpha diversity
  ALPHA_INDICES = ALPHA_INDICES,
  ALPHA_GROUP_VAR = ALPHA_GROUP_VAR,
  ALPHA_METRIC = ALPHA_METRIC,
  ALPHA_COVARIATES = ALPHA_COVARIATES,
  ALPHA_REFERENCE_GROUP = ALPHA_REFERENCE_GROUP,
  
  # Beta diversity
  BETA_DISTANCE = BETA_DISTANCE,
  BETA_ORD_METHOD = BETA_ORD_METHOD,
  BETA_PERMANOVA_VARS = BETA_PERMANOVA_VARS,
  BETA_COLOR_VAR = BETA_COLOR_VAR,
  BETA_SHAPE_VAR = BETA_SHAPE_VAR,
  BETA_VIS_STYLE = BETA_VIS_STYLE,
  
  # Relative abundance
  RELAB_TAX_LEVEL = RELAB_TAX_LEVEL,
  RELAB_GROUP_VAR = RELAB_GROUP_VAR,
  RELAB_FACET_VAR = RELAB_FACET_VAR,
  RELAB_THRESHOLD = RELAB_THRESHOLD,
  RELAB_SORT_TAXON = RELAB_SORT_TAXON,
  MIN_TAXONOMY_CONFIDENCE = MIN_TAXONOMY_CONFIDENCE,
  
  # Differential abundance
  DAA_TAX_LEVEL = DAA_TAX_LEVEL,
  DAA_FIXED_EFFECTS = DAA_FIXED_EFFECTS,
  DAA_RANDOM_EFFECTS = DAA_RANDOM_EFFECTS,
  DAA_COVARIATES = DAA_COVARIATES,
  DAA_REFERENCE = DAA_REFERENCE,
  DAA_NORMALIZATION = DAA_NORMALIZATION,
  DAA_ANALYSIS_METHOD = DAA_ANALYSIS_METHOD,
  DAA_MIN_ABUNDANCE = DAA_MIN_ABUNDANCE,
  DAA_MIN_PREVALENCE = DAA_MIN_PREVALENCE,
  DAA_MAX_SIGNIFICANCE = DAA_MAX_SIGNIFICANCE,
  
  # Plotting
  PLOT_WIDTH = PLOT_WIDTH,
  PLOT_HEIGHT = PLOT_HEIGHT,
  PLOT_DPI = PLOT_DPI,
  COLOR_PALETTE = COLOR_PALETTE
)

# Run analyses
message("\n========================================")
cat("\nRunning downstream analyses...\n")

# 1. Alpha diversity
cat("\n1. Alpha diversity...")
source("01_alpha_diversity.R")
alpha_results <- run_alpha_diversity(ps_filtered, config, OUTPUT_DIR)

# 2. Beta diversity
cat("\n2. Beta diversity...")
source("02_beta_diversity.R")
beta_results <- run_beta_diversity(ps_filtered, config, OUTPUT_DIR)

# 3. Relative abundance
cat("\n3. Relative abundance...")
source("03_relative_abundance.R")
relab_results <- run_relative_abundance(ps_filtered, config, OUTPUT_DIR)

# 4. Differential abundance - MaAslin2
cat("\n4. Differential abundance (MaAslin2)...\n")
tryCatch({
  source("04_differential_abundance_maaslin2.R")
  maaslin2_results <- run_maaslin2_daa(ps_filtered, config, OUTPUT_DIR)
}, error = function(e) {
  cat("Error in MaAslin2:", e$message, "\n")
  cat("Make sure MaAslin2 is installed: BiocManager::install('Maaslin2')\n")
})


# =============================================================================
# GENERATE SUMMARY REPORT
# =============================================================================

cat("\n=== Generating Summary Report ===\n\n")

summary_file <- file.path(OUTPUT_DIR, "analysis_summary.txt")
sink(summary_file)

cat("========================================\n")
cat("16S Downstream Analysis Summary\n")
cat("========================================\n\n")
cat(paste0("Date: ", Sys.Date(), "\n"))
cat(paste0("Output directory: ", OUTPUT_DIR, "\n\n"))

cat("DATA SUMMARY\n")
cat("------------\n")
cat(paste0("Total samples: ", nsamples(ps_filtered), "\n"))
cat(paste0("Total taxa: ", ntaxa(ps_filtered), "\n"))
cat(paste0("Rarefaction performed: ", PERFORM_RAREFACTION, "\n"))
if (PERFORM_RAREFACTION) {
  cat(paste0("Rarefaction depth: ", RAREFACTION_DEPTH, "\n"))
}
cat("\n")

cat("ANALYSES PERFORMED\n")
cat("------------------\n")
cat("✓ Alpha diversity\n")
cat(paste0("  - Indices: ", paste(ALPHA_INDICES, collapse = ", "), "\n"))
cat(paste0("  - Grouping: ", ALPHA_GROUP_VAR, "\n\n"))

cat("✓ Beta diversity\n")
cat(paste0("  - Method: ", BETA_ORD_METHOD, " (", BETA_DISTANCE, " distance)\n"))
cat(paste0("  - PERMANOVA variables: ", paste(BETA_PERMANOVA_VARS, collapse = ", "), "\n\n"))

cat("✓ Relative abundance\n")
cat(paste0("  - Taxonomic level: ", RELAB_TAX_LEVEL, "\n"))
cat(paste0("  - Grouping: ", RELAB_GROUP_VAR, "\n\n"))

cat("✓ Differential abundance (MaAslin2)\n")
cat("  - Analysis level: ASV (filtered by abundance)\n")
cat(paste0("  - Fixed effects: ", paste(DAA_FIXED_EFFECTS, collapse = ", "), "\n"))
if (!is.null(DAA_RANDOM_EFFECTS)) {
  cat(paste0("  - Random effects: ", paste(DAA_RANDOM_EFFECTS, collapse = ", "), "\n"))
}
cat("\n")

cat("OUTPUT FILES\n")
cat("------------\n")
cat("Results are organized in subdirectories:\n")
cat("  - alpha_diversity/\n")
cat("  - beta_diversity/\n")
cat("  - relative_abundance/\n")
cat("  - differential_abundance_maaslin2/\n\n")

cat("KEY FILES\n")
cat("---------\n")
cat("- phyloseq_filtered.rds : Preprocessed phyloseq object\n")
cat("- alpha_diversity/alpha_diversity_indices.csv\n")
cat("- beta_diversity/permanova_results.csv\n")
cat("- relative_abundance/relative_abundance_*.csv\n")
cat("- differential_abundance_maaslin2/maaslin2_significant_results.csv\n\n")

cat("========================================\n")
cat("Analysis Complete!\n")
cat("========================================\n")

sink()

message(paste0("Saved: ", summary_file))

message("\n========================================")
message("ALL ANALYSES COMPLETE!")
message("========================================\n")
message(paste0("Results saved to: ", OUTPUT_DIR))
message("\nTo view the summary report:")
message(paste0("  cat ", summary_file))
message("\nTo explore results:")
message(paste0("  cd ", OUTPUT_DIR))
message("  ls -R\n")
