################################################################################
# Beta Diversity Analysis
################################################################################

library(phyloseq)
library(vegan)
library(ape)
library(dplyr)
library(tidyr)
library(ggplot2)

# Helper function to compute convex hull
compute_hulls <- function(df, x_var = "Axis.1", y_var = "Axis.2", group_var) {
  df %>%
    group_by(!!sym(group_var)) %>%
    slice(chull(!!sym(x_var), !!sym(y_var))) %>%
    ungroup()
}

# Calculate distance matrix
calculate_distance <- function(ps, method) {
  
  # Check if method requires tree
  if ((method == "unifrac" || method == "wunifrac") && is.null(phy_tree(ps, errorIfNULL = FALSE))) {
    cat("Warning: UniFrac requires phylogenetic tree - skipping\n")
    return(NULL)
  }
  
  # Calculate distance
  dist_matrix <- phyloseq::distance(ps, method = method)
  return(dist_matrix)
}

# Perform ordination
perform_ordination <- function(ps, dist_matrix, method = "PCoA") {
  
  if (method == "PCoA") {
    ord <- ordinate(ps, method = "PCoA", distance = dist_matrix)
  } else if (method == "NMDS") {
    ord <- ordinate(ps, method = "NMDS", distance = dist_matrix)
  }
  
  return(ord)
}

# Plot ordination with hulls and centroids
plot_ordination_custom <- function(ps, ord, color_var, shape_var = NULL, output_file, 
                                   width = 8, height = 6, dpi = 300, color_palette = NULL,
                                   permanova_results = NULL) {
  
  # Base ordination plot
  p <- plot_ordination(ps, ord, color = color_var, shape = shape_var)
  
  # Extract plot data
  plot_data <- p$data
  
  # Calculate centroids for each group
  if (!is.null(color_var) && color_var %in% colnames(plot_data)) {
    centroids <- plot_data %>%
      group_by(!!sym(color_var)) %>%
      summarize(
        Axis.1 = mean(Axis.1, na.rm = TRUE),
        Axis.2 = mean(Axis.2, na.rm = TRUE),
        .groups = "drop"
      )
    
    # Calculate convex hulls
    hulls <- compute_hulls(plot_data, "Axis.1", "Axis.2", color_var)
  }
  
  # Build plot with hulls and centroids
  p <- p +
    geom_polygon(data = hulls, 
                 aes_string(x = "Axis.1", y = "Axis.2", 
                           fill = color_var, color = color_var),
                 alpha = 0.1, show.legend = FALSE) +
    geom_point(size = 3, alpha = 0.7) +
    geom_point(data = centroids, 
               aes_string(x = "Axis.1", y = "Axis.2", color = color_var),
               size = 6, shape = 17, stroke = 1.5) +
    theme_classic() +
    theme(legend.position = "right")
  
  # Apply custom color palette if provided
  if (!is.null(color_palette)) {
    p <- p +
      scale_color_manual(values = color_palette) +
      scale_fill_manual(values = color_palette)
  }
  
  # Add PERMANOVA caption
  if (!is.null(permanova_results)) {
    caption_texts <- c()
    for (var in names(permanova_results)) {
      perm_res <- permanova_results[[var]]
      if (!is.null(perm_res)) {
        r2 <- perm_res$R2[1]
        p_val <- perm_res$`Pr(>F)`[1]
        
        p_text <- if (p_val < 0.001) {
          "< 0.001"
        } else {
          paste0("= ", round(p_val, 3))
        }
        sig <- if (p_val < 0.05) " *" else ""
        caption_texts <- c(caption_texts, paste0(var, " R²=", round(r2, 3), ", p", p_text, sig))
      }
    }
    if (length(caption_texts) > 0) {
      p <- p + labs(caption = paste0("PERMANOVA: ", paste(caption_texts, collapse = "; "))) +
        theme(plot.caption = element_text(hjust = 0, size = 9, face = "italic"))
    }
  }
  
  ggsave(output_file, plot = p, width = width, height = height, dpi = dpi)
  return(p)
}

# PERMANOVA test
run_permanova <- function(dist_matrix, metadata, variables, nperm = 999) {
  
  # Build formula with all variables for marginal tests (Type II)
  formula_str <- paste("dist_matrix ~", paste(variables, collapse = " + "))
  
  # Run adonis2 with by="margin" to get unique contribution of each variable
  perm_result <- adonis2(as.formula(formula_str), data = metadata, 
                         permutations = nperm, by = "margin")
  
  # Extract results for each variable
  results_list <- list()
  for (var in variables) {
    if (var %in% rownames(perm_result)) {
      # Create a single-row data frame for this variable
      var_result <- data.frame(
        R2 = perm_result$R2[which(rownames(perm_result) == var)],
        `Pr(>F)` = perm_result$`Pr(>F)`[which(rownames(perm_result) == var)],
        check.names = FALSE
      )
      results_list[[var]] <- var_result
    }
  }
  
  # Store full model for reference
  results_list$full_model <- perm_result
  
  return(results_list)
}

# Main analysis function
run_beta_diversity <- function(ps, config, output_dir) {
  cat("\n=== Beta Diversity Analysis ===\n\n")
  
  beta_dir <- file.path(output_dir, "beta_diversity")
  dir.create(beta_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Get metadata
  metadata <- data.frame(sample_data(ps), stringsAsFactors = FALSE, 
                        row.names = rownames(sample_data(ps)))
  
  permanova_results_all <- list()
  
  # Loop through distance methods
  for (distance in config$BETA_DISTANCE) {
    cat("Processing", distance, "distance...\n")
    
    # Calculate distance
    dist_matrix <- calculate_distance(ps, distance)
    if (is.null(dist_matrix)) next
    
    # Save distance matrix
    dist_file <- file.path(beta_dir, paste0("distance_matrix_", distance, ".rds"))
    saveRDS(dist_matrix, dist_file)
    
    # PERMANOVA - run first to include in plots
    permanova_results <- run_permanova(dist_matrix, metadata, config$BETA_PERMANOVA_VARS)
    permanova_results_all[[distance]] <- permanova_results
    
    # Ordination
    for (ord_method in config$BETA_ORD_METHOD) {
      ord <- perform_ordination(ps, dist_matrix, ord_method)
      
      # Save ordination
      ord_file <- file.path(beta_dir, paste0("ordination_", ord_method, "_", distance, ".rds"))
      saveRDS(ord, ord_file)
      
      # Plot with PERMANOVA results
      plot_file <- file.path(beta_dir, paste0("ordination_", ord_method, "_", distance, ".png"))
      plot_ordination_custom(ps, ord, config$BETA_COLOR_VAR, config$BETA_SHAPE_VAR, 
                           plot_file, config$PLOT_WIDTH, config$PLOT_HEIGHT, config$PLOT_DPI,
                           config$COLOR_PALETTE, permanova_results)
    }
  }
  
  # Save PERMANOVA results
  perm_file <- file.path(beta_dir, "permanova_results.txt")
  sink(perm_file)
  cat("PERMANOVA Results\n\n")
  for (distance in names(permanova_results_all)) {
    cat("Distance method:", distance, "\n")
    cat("======================\n\n")
    for (var in names(permanova_results_all[[distance]])) {
      cat("Variable:", var, "\n")
      print(permanova_results_all[[distance]][[var]])
      cat("\n")
    }
  }
  sink()
  cat("Saved:", perm_file, "\n")
  
  cat("\nBeta diversity analysis complete!\n")
  return(permanova_results_all)
}
