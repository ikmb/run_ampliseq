################################################################################
# Differential Abundance Analysis with MaAslin2
################################################################################

library(phyloseq)
library(Maaslin2)
library(dplyr)
library(ggplot2)
library(genefilter)

# Prepare data for MaAslin2 (ASV-level analysis)
prepare_maaslin2_data <- function(ps, tax_level = "Genus", min_reads = 4, min_prevalence = 0.1) {
  
  # Get taxonomy table
  tax_df <- as.data.frame(tax_table(ps)) %>%
    mutate(across(everything(), as.character))
  
  # Apply abundance filtering at ASV level
  min_samples <- ceiling(nsamples(ps) * min_prevalence)
  abundance_filter <- genefilter_sample(
    ps,
    filterfun_sample(function(x) x >= min_reads),
    A = min_samples
  )
  
  # Filter taxa with valid genus names
  if (tax_level == "Genus") {
    genus <- tax_df[, "Genus"]
    genus_filter <- !is.na(genus) & trimws(genus) != "" & genus != "Unknown"
  } else {
    genus_filter <- rep(TRUE, nrow(tax_df))
  }
  
  # Combine filters
  combined_filter <- abundance_filter & genus_filter
  ps_filt <- prune_taxa(combined_filter, ps)
  
  cat("Kept", ntaxa(ps_filt), "of", ntaxa(ps), "ASVs after filtering\n")
  
  # Get filtered taxonomy and strip GTDB prefixes
  tax_df_filt <- as.data.frame(tax_table(ps_filt)) %>%
    mutate(across(everything(), as.character))
  
  for (col in colnames(tax_df_filt)) {
    tax_df_filt[[col]] <- gsub("^[a-z]__", "", tax_df_filt[[col]])
  }
  
  # Label ASVs as "ASV_X_Genus"
  asv_labels <- paste(taxa_names(ps_filt), tax_df_filt[[tax_level]], sep = "_")
  asv_labels <- make.unique(asv_labels, sep = ".")
  taxa_names(ps_filt) <- asv_labels
  
  # Extract abundance table (samples as rows)
  otu_df <- as.data.frame(otu_table(ps_filt))
  if (taxa_are_rows(ps_filt)) {
    otu_df <- t(otu_df)
  }
  otu_df <- round(otu_df)  # Round to integers
  
  # Extract metadata
  meta_df <- data.frame(sample_data(ps_filt), stringsAsFactors = FALSE, 
                       row.names = rownames(sample_data(ps_filt)))
  
  return(list(abundance = otu_df, metadata = meta_df))
}

# Create summary plot
create_summary_plot <- function(results_df, output_dir, max_significance = 0.05, 
                               width = 10, height = 8, dpi = 300) {
  
  sig_results <- results_df %>%
    filter(qval < max_significance) %>%
    arrange(qval) %>%
    head(20)
  
  if (nrow(sig_results) == 0) {
    cat("No significant results to plot\n")
    return(NULL)
  }
  
  p <- ggplot(sig_results, aes(x = coef, y = reorder(feature, coef), color = metadata)) +
    geom_point(size = 3) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme_classic() +
    labs(x = "Coefficient", y = "Feature", color = "Variable",
         title = "Top Significant Associations")
  
  plot_file <- file.path(output_dir, "maaslin2_summary_plot.png")
  ggsave(plot_file, plot = p, width = width, height = height, dpi = dpi)
  
  return(p)
}

# Main analysis function
run_maaslin2_daa <- function(ps, config, output_dir) {
  cat("\n=== MaAslin2 Differential Abundance Analysis ===\n\n")
  
  maaslin_dir <- file.path(output_dir, "differential_abundance_maaslin2")
  dir.create(maaslin_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Prepare data (ASV-level)
  data_list <- prepare_maaslin2_data(ps, config$DAA_TAX_LEVEL, min_reads = 4, min_prevalence = 0.1)
  
  # Combine fixed effects and covariates
  all_fixed_effects <- c(config$DAA_FIXED_EFFECTS, config$DAA_COVARIATES)
  
  # Set reference levels
  metadata <- data_list$metadata
  if (!is.null(config$DAA_REFERENCE)) {
    for (var in names(config$DAA_REFERENCE)) {
      if (var %in% colnames(metadata)) {
        metadata[[var]] <- relevel(factor(metadata[[var]]), ref = config$DAA_REFERENCE[var])
      }
    }
  }
  
  # Run MaAslin2
  cat("Running MaAslin2...\n")
  maaslin_results <- Maaslin2(
    input_data = data_list$abundance,
    input_metadata = metadata,
    output = file.path(maaslin_dir, "maaslin2_output"),
    fixed_effects = all_fixed_effects,
    random_effects = config$DAA_RANDOM_EFFECTS,
    normalization = config$DAA_NORMALIZATION,
    transform = "LOG",
    analysis_method = "LM",
    min_abundance = config$DAA_MIN_ABUNDANCE,
    min_prevalence = config$DAA_MIN_PREVALENCE,
    max_significance = config$DAA_MAX_SIGNIFICANCE,
    correction = "BH",
    standardize = TRUE,
    cores = 1
  )
  
  # Read and save results
  results_file <- file.path(maaslin_dir, "maaslin2_output", "all_results.tsv")
  
  if (file.exists(results_file)) {
    all_results <- read.table(results_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    sig_results <- all_results %>% filter(qval < config$DAA_MAX_SIGNIFICANCE)
    
    # Save cleaned results
    write.csv(sig_results, file.path(maaslin_dir, "maaslin2_significant_results.csv"), row.names = FALSE)
    write.csv(all_results, file.path(maaslin_dir, "maaslin2_all_results.csv"), row.names = FALSE)
    
    cat("Found", nrow(sig_results), "significant associations\n")
    
    # Create summary plot
    if (nrow(sig_results) > 0) {
      create_summary_plot(sig_results, maaslin_dir, config$DAA_MAX_SIGNIFICANCE,
                         config$PLOT_WIDTH, config$PLOT_HEIGHT, config$PLOT_DPI)
    }
    
    # Save summary
    summary_file <- file.path(maaslin_dir, "maaslin2_summary.txt")
    sink(summary_file)
    cat("MaAslin2 Differential Abundance Analysis\n\n")
    cat("Analysis level: ASV (filtered by abundance and genus confidence)\n")
    cat("Fixed effects:", paste(all_fixed_effects, collapse = ", "), "\n")
    cat("Normalization:", config$DAA_NORMALIZATION, "\n")
    cat("Total associations tested:", nrow(all_results), "\n")
    cat("Significant associations:", nrow(sig_results), "\n\n")
    
    if (nrow(sig_results) > 0) {
      cat("Top 10 significant associations:\n")
      print(head(sig_results[, c("feature", "metadata", "coef", "pval", "qval")], 10))
    }
    sink()
    cat("Saved:", summary_file, "\n")
  }
  
  cat("\nMaAslin2 analysis complete!\n")
  return(maaslin_results)
}
