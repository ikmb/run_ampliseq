################################################################################
# Alpha Diversity Analysis
################################################################################

library(phyloseq)
library(microbiome)
library(dplyr)
library(tidyr)
library(ggplot2)
library(car)
library(FSA)
library(picante)

# Calculate alpha diversity indices
calculate_alpha_diversity <- function(ps, indices) {
  
  # Separate faith_pd from other indices
  has_faith <- "faith_pd" %in% indices
  other_indices <- setdiff(indices, "faith_pd")
  
  # Calculate standard indices
  if (length(other_indices) > 0) {
    alpha_df <- microbiome::alpha(ps, index = other_indices)
  } else {
    alpha_df <- data.frame(row.names = sample_names(ps))
  }
  
  # Calculate Faith's PD if tree available
  if (has_faith && !is.null(phy_tree(ps, errorIfNULL = FALSE))) {
    otu_tab <- as(otu_table(ps), "matrix")
    if (!taxa_are_rows(ps)) otu_tab <- t(otu_tab)
    
    otu_tab_t <- t(otu_tab > 0)
    tree <- phy_tree(ps)
    faith_values <- picante::pd(otu_tab_t, tree, include.root = FALSE)
    alpha_df$faith_pd <- faith_values$PD[match(rownames(alpha_df), rownames(faith_values))]
  }
  
  # Add sample IDs
  alpha_df$SampleID <- as.character(rownames(alpha_df))
  
  # Merge with metadata
  meta <- data.frame(sample_data(ps), stringsAsFactors = FALSE)
  meta$SampleID <- as.character(rownames(meta))
  alpha_merged <- left_join(alpha_df, meta, by = "SampleID")
  
  return(alpha_merged)
}

# Plot all indices
plot_all_alpha_indices <- function(alpha_df, group_var, output_file, width = 10, height = 6, dpi = 300, color_palette = NULL) {
  
  index_cols <- c("observed", "diversity_shannon", "diversity_gini_simpson",
                  "evenness_pielou", "diversity_fisher", "faith_pd")
  index_cols <- intersect(index_cols, colnames(alpha_df))
  
  alpha_long <- alpha_df %>%
    select(SampleID, all_of(group_var), all_of(index_cols)) %>%
    pivot_longer(cols = all_of(index_cols), names_to = "Index", values_to = "Value")
  
  p <- ggplot(alpha_long, aes_string(x = group_var, y = "Value", fill = group_var)) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
    facet_wrap(~ Index, scales = "free_y", ncol = 3) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text = element_text(size = 10, face = "bold"),
          legend.position = "bottom") +
    labs(x = NULL, y = "Index Value", fill = group_var)
  
  # Apply custom color palette if provided
  if (!is.null(color_palette)) {
    p <- p + scale_fill_manual(values = color_palette)
  }
  
  ggsave(output_file, plot = p, width = width, height = height, dpi = dpi)
  return(p)
}

# Plot single metric
plot_alpha_metric <- function(alpha_df, metric, group_var, output_file, 
                              width = 8, height = 6, dpi = 300, color_palette = NULL, 
                              test_results = NULL) {
  
  p <- ggplot(alpha_df, aes_string(x = group_var, y = metric, fill = group_var)) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
    labs(x = NULL, y = gsub("_", " ", tools::toTitleCase(metric)))
  
  # Apply custom color palette
  if (!is.null(color_palette)) {
    p <- p + scale_fill_manual(values = color_palette)
  }
  
  # Add statistical test caption
  if (!is.null(test_results)) {
    caption_text <- ""
    if (!is.null(test_results$p_value)) {
      # Simple test (Kruskal-Wallis or Wilcoxon)
      p_val_text <- if (test_results$p_value < 0.001) {
        "p < 0.001"
      } else {
        paste0("p = ", round(test_results$p_value, 3))
      }
      caption_text <- paste0(test_results$test, " test: ", p_val_text)
      if (test_results$p_value < 0.05) {
        caption_text <- paste0(caption_text, " *")
      }
    } else if (!is.null(test_results$anova)) {
      # ANOVA with covariates - show all variables
      anova_df <- as.data.frame(test_results$anova)
      p_values <- anova_df[, "Pr(>F)"]
      var_names <- rownames(anova_df)
      
      # Format each variable's p-value (excluding Residuals row)
      var_texts <- c()
      for (i in seq_along(var_names)) {
        if (!is.na(p_values[i]) && var_names[i] != "Residuals") {
          p_text <- if (p_values[i] < 0.001) {
            "< 0.001"
          } else {
            paste0("= ", round(p_values[i], 3))
          }
          sig <- if (p_values[i] < 0.05) " *" else ""
          var_texts <- c(var_texts, paste0(var_names[i], " p", p_text, sig))
        }
      }
      caption_text <- paste0("ANOVA: ", paste(var_texts, collapse = "; "))
    }
    
    if (nchar(caption_text) > 0) {
      p <- p + labs(caption = caption_text) +
        theme(plot.caption = element_text(hjust = 0, size = 9, face = "italic"))
    }
  }
  
  ggsave(output_file, plot = p, width = width, height = height, dpi = dpi)
  return(p)
}

# Statistical testing
test_alpha_diversity <- function(alpha_df, metric, group_var, covariates = NULL) {
  
  # Remove missing values
  test_df <- alpha_df %>%
    filter(!is.na(.data[[metric]]), !is.na(.data[[group_var]]))
  
  if (!is.null(covariates)) {
    test_df <- test_df %>% filter(if_all(all_of(covariates), ~ !is.na(.)))
  }
  
  results <- list()
  
  # Simple comparison without covariates
  if (is.null(covariates)) {
    n_groups <- length(unique(test_df[[group_var]]))
    
    if (n_groups == 2) {
      wtest <- wilcox.test(as.formula(paste(metric, "~", group_var)), data = test_df)
      results$test <- "Wilcoxon"
      results$statistic <- wtest$statistic
      results$p_value <- wtest$p.value
    } else if (n_groups > 2) {
      kwtest <- kruskal.test(as.formula(paste(metric, "~", group_var)), data = test_df)
      results$test <- "Kruskal-Wallis"
      results$statistic <- kwtest$statistic
      results$p_value <- kwtest$p.value
      
      if (kwtest$p.value < 0.05) {
        dunn_test <- FSA::dunnTest(as.formula(paste(metric, "~", group_var)), 
                                   data = test_df, method = "bh")
        results$posthoc <- dunn_test$res
      }
    }
  } else {
    # Multivariable model with covariates
    formula_str <- paste(metric, "~", paste(c(group_var, covariates), collapse = " + "))
    lm_model <- lm(as.formula(formula_str), data = test_df)
    anova_res <- car::Anova(lm_model, type = 2)
    
    results$test <- "ANOVA"
    results$model <- lm_model
    results$anova <- anova_res
    results$summary <- summary(lm_model)
  }
  
  return(results)
}

# Main analysis function
run_alpha_diversity <- function(ps, config, output_dir) {
  cat("\n=== Alpha Diversity Analysis ===\n\n")
  
  alpha_dir <- file.path(output_dir, "alpha_diversity")
  dir.create(alpha_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Calculate indices
  alpha_df <- calculate_alpha_diversity(ps, config$ALPHA_INDICES)
  
  # Save results
  alpha_file <- file.path(alpha_dir, "alpha_diversity_indices.csv")
  write.csv(alpha_df, alpha_file, row.names = FALSE)
  cat("Saved:", alpha_file, "\n")
  
  # Plot all indices
  plot_all_alpha_indices(
    alpha_df, config$ALPHA_GROUP_VAR,
    file.path(alpha_dir, "alpha_diversity_all_indices.png"),
    config$PLOT_WIDTH, config$PLOT_HEIGHT, config$PLOT_DPI,
    config$COLOR_PALETTE
  )
  
  # Analyze specific metric
  if (!is.null(config$ALPHA_METRIC)) {
    # Run statistical test first
    test_results <- test_alpha_diversity(
      alpha_df, config$ALPHA_METRIC, config$ALPHA_GROUP_VAR, config$ALPHA_COVARIATES
    )
    
    # Plot with test results
    plot_alpha_metric(
      alpha_df, config$ALPHA_METRIC, config$ALPHA_GROUP_VAR,
      file.path(alpha_dir, paste0("alpha_diversity_", config$ALPHA_METRIC, ".png")),
      config$PLOT_WIDTH * 0.8, config$PLOT_HEIGHT, config$PLOT_DPI,
      config$COLOR_PALETTE,
      test_results
    )
    
    # Save stats
    results_file <- file.path(alpha_dir, paste0("alpha_diversity_", config$ALPHA_METRIC, "_stats.txt"))
    sink(results_file)
    cat("Alpha Diversity Statistical Analysis:", config$ALPHA_METRIC, "\n\n")
    cat("Test:", test_results$test, "\n")
    
    if (!is.null(test_results$p_value)) {
      cat("Statistic:", round(test_results$statistic, 4), "\n")
      cat("P-value:", signif(test_results$p_value, 3), "\n")
      cat("Significant:", ifelse(test_results$p_value < 0.05, "YES", "NO"), "\n\n")
    }
    
    if (!is.null(test_results$posthoc)) {
      cat("\nPost-hoc comparisons:\n")
      print(test_results$posthoc)
    }
    
    if (!is.null(test_results$anova)) {
      cat("\nANOVA results:\n")
      print(test_results$anova)
      cat("\nModel summary:\n")
      print(test_results$summary)
    }
    sink()
    cat("Saved:", results_file, "\n")
  }
  
  cat("\nAlpha diversity analysis complete!\n")
  return(alpha_df)
}
