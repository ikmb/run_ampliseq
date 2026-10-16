################################################################################
# Relative Abundance Analysis
################################################################################

library(phyloseq)
library(dplyr)
library(tidyr)
library(ggplot2)

# Agglomerate to taxonomic level (keeps unknown taxa)
agglomerate_taxa <- function(ps, tax_level, min_conf = 0.5) {
  
  # Get taxonomy table
  tax_df <- as.data.frame(tax_table(ps))
  
  # Mark low-confidence taxa as "Unknown" before agglomeration
  if (tax_level %in% colnames(tax_df)) {
    conf_col <- paste0("confidence_", tolower(tax_level))
    if (conf_col %in% colnames(tax_df)) {
      conf_values <- as.numeric(tax_df[[conf_col]])
      low_conf_idx <- is.na(conf_values) | conf_values < min_conf
      tax_df[[tax_level]][low_conf_idx] <- "Unknown"
      
      # Update the taxonomy table in phyloseq
      tax_table(ps) <- as.matrix(tax_df)
    }
    
    # Also mark empty or NA taxonomy as "Unknown"
    empty_idx <- is.na(tax_df[[tax_level]]) | tax_df[[tax_level]] == ""
    tax_df[[tax_level]][empty_idx] <- "Unknown"
    tax_table(ps) <- as.matrix(tax_df)
  }
  
  # Agglomerate - NArm = FALSE keeps taxa with "Unknown"
  ps_agg <- tax_glom(ps, taxrank = tax_level, NArm = FALSE)
  return(ps_agg)
}

# Calculate relative abundance
calculate_relative_abundance <- function(ps, tax_level) {
  
  # Transform to relative abundance (as proportion 0-1)
  ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
  
  # Get taxonomy and abundance
  tax_df <- as.data.frame(tax_table(ps_rel))
  otu_df <- as.data.frame(otu_table(ps_rel))
  
  # Strip GTDB prefixes
  for (col in colnames(tax_df)) {
    tax_df[[col]] <- gsub("^[a-z]__", "", tax_df[[col]])
  }
  
  # Get taxon names (confidence filtering was already done in agglomeration)
  if (tax_level %in% colnames(tax_df)) {
    taxon_names <- tax_df[[tax_level]]
  } else {
    taxon_names <- taxa_names(ps_rel)
  }
  
  # Add taxon names to otu_df
  otu_df$Taxon <- taxon_names
  
  # Convert to long format
  otu_long <- otu_df %>%
    pivot_longer(-Taxon, names_to = "SampleID", values_to = "Abundance")
  
  # Merge with metadata
  meta <- data.frame(sample_data(ps_rel), stringsAsFactors = FALSE)
  meta$SampleID <- rownames(meta)
  
  relab_df <- left_join(otu_long, meta, by = "SampleID")
  
  return(relab_df)
}

# Plot stacked bar chart
plot_relative_abundance <- function(relab_df, group_var, output_file, 
                                   width = 12, height = 6, dpi = 300, facet_var = NULL) {
  
  # Define color palette (Paul Tol's muted qualitative palette + extensions)
  sophisticated_colors <- c(
    "#332288", "#88CCEE", "#44AA99", "#117733", "#999933",
    "#DDCC77", "#CC6677", "#882255", "#AA4499", "#661100",
    "#6699CC", "#AA4466", "#4477AA", "#228833", "#CCBB44",
    "#EE6677", "#AA3377", "#BBBBBB", "#999999", "#666666"
  )
  
  # If showing individual samples, use raw data; otherwise aggregate
  if (group_var == "SampleID") {
    plot_data <- relab_df
  } else {
    # Prepare grouping variables for aggregation
    group_vars <- group_var
    if (!is.null(facet_var)) {
      group_vars <- c(group_var, facet_var)
    }
    
    # Summarize by group and taxon
    plot_data <- relab_df %>%
      group_by(across(all_of(c(group_vars, "Taxon")))) %>%
      summarize(Abundance = mean(Abundance), .groups = "drop")
  }
  
  # Reorder taxa: put Unknown at the bottom (first level)
  taxa_levels <- unique(plot_data$Taxon)
  if ("Unknown" %in% taxa_levels) {
    # Unknown first, then all others
    taxa_levels <- c("Unknown", setdiff(taxa_levels, "Unknown"))
  }
  plot_data$Taxon <- factor(plot_data$Taxon, levels = taxa_levels)
  
  # Create plot
  p <- ggplot(plot_data, aes_string(x = group_var, y = "Abundance", fill = "Taxon")) +
    geom_bar(stat = "identity", position = "stack") +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = sophisticated_colors) +
    theme_classic() +
    theme(legend.position = "right") +
    labs(x = NULL, y = "Relative Abundance", fill = "Taxon")
  
  # Add faceting if specified
  if (!is.null(facet_var)) {
    p <- p + facet_grid(as.formula(paste("~", facet_var)), scales = "free_x", space = "free_x") +
      theme(axis.text.x = element_blank(),  # Hide x-axis labels when showing many samples
            axis.ticks.x = element_blank(),
            strip.text = element_text(size = 12, face = "bold"))
  } else {
    p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  
  ggsave(output_file, plot = p, width = width, height = height, dpi = dpi)
  return(p)
}

# Main analysis function
run_relative_abundance <- function(ps, config, output_dir) {
  cat("\n=== Relative Abundance Analysis ===\n\n")
  
  relab_dir <- file.path(output_dir, "relative_abundance")
  dir.create(relab_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Apply abundance filter: taxa must have >= 4 reads in >= 10% of samples
  filter_taxa <- genefilter_sample(ps, filterfun_sample(function(x) x >= 4), 
                                   A = 0.1 * nsamples(ps))
  ps_filtered <- prune_taxa(filter_taxa, ps)
  cat("After abundance filtering:", ntaxa(ps_filtered), "taxa (removed", 
      ntaxa(ps) - ntaxa(ps_filtered), "rare taxa)\n")
  
  # Agglomerate to taxonomic level
  ps_agg <- agglomerate_taxa(ps_filtered, config$RELAB_TAX_LEVEL, config$MIN_TAXONOMY_CONFIDENCE)
  cat("Agglomerated to", config$RELAB_TAX_LEVEL, "level:", ntaxa(ps_agg), "taxa\n")
  
  # Calculate relative abundance
  relab_df <- calculate_relative_abundance(ps_agg, config$RELAB_TAX_LEVEL)
  
  # Save table
  relab_file <- file.path(relab_dir, paste0("relative_abundance_", config$RELAB_TAX_LEVEL, ".csv"))
  write.csv(relab_df, relab_file, row.names = FALSE)
  cat("Saved:", relab_file, "\n")
  
  # Plot
  plot_file <- file.path(relab_dir, paste0("relative_abundance_", config$RELAB_TAX_LEVEL, "_barplot.png"))
  plot_relative_abundance(relab_df, config$RELAB_GROUP_VAR, plot_file,
                         config$PLOT_WIDTH, config$PLOT_HEIGHT, config$PLOT_DPI,
                         config$RELAB_FACET_VAR)
  
  cat("\nRelative abundance analysis complete!\n")
  return(relab_df)
}
