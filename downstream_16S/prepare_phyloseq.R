################################################################################
# prepare_phyloseq.R
#
# Prepare phyloseq object from nf-core/ampliseq outputs
# Handles intersection of feature table with taxonomy and ASV renaming
################################################################################

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  library(ape)
})

#' Prepare Phyloseq Object from ampliseq outputs
#'
#' @param feature_table_file Path to feature-table.tsv (from QIIME2)
#' @param tax_genus_file Path to ASV_tax.gtdb.tsv (taxonomy to genus)
#' @param tax_species_file Path to ASV_tax_species.gtdb.tsv (taxonomy to species)
#' @param metadata_file Path to metadata file
#' @param tree_file Path to tree.nwk file
#' @param sample_prefix Prefix to remove from sample names (e.g., "S_")
#' @param output_file Output RDS file path
#' @return Phyloseq object
prepare_phyloseq <- function(feature_table_file,
                            tax_genus_file,
                            tax_species_file,
                            metadata_file,
                            tree_file,
                            sample_prefix = "S_",
                            output_file = "phyloseq_object.rds") {
  
  message("========================================")
  message("Preparing Phyloseq Object from ampliseq")
  message("========================================\n")
  
  # ===========================
  # 1. Load feature table
  # ===========================
  message("Loading feature table...")
  # Read feature table, skipping first comment line
  feature_lines <- readLines(feature_table_file)
  header_line <- grep("^#OTU ID", feature_lines)[1]
  
  feature_df <- read.table(
    feature_table_file,
    sep = "\t",
    header = TRUE,
    skip = header_line - 1,
    comment.char = "",
    check.names = FALSE,
    row.names = 1
  )
  
  # Remove '#' from first column name if present
  colnames(feature_df)[1] <- gsub("^#", "", colnames(feature_df)[1])
  
  # Get ASV IDs that passed filtering
  kept_asvs <- rownames(feature_df)
  message(paste0("  Feature table has ", length(kept_asvs), " ASVs"))
  
  # Remove sample prefix
  if (!is.null(sample_prefix) && sample_prefix != "") {
    colnames(feature_df) <- gsub(paste0("^", sample_prefix), "", colnames(feature_df))
  }
  message(paste0("  Feature table has ", ncol(feature_df), " samples"))
  
  # ===========================
  # 2. Load and filter taxonomy
  # ===========================
  message("\nLoading taxonomy tables...")
  
  # Load genus-level taxonomy
  tax_genus <- read.table(
    tax_genus_file,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  message(paste0("  Genus taxonomy has ", nrow(tax_genus), " ASVs"))
  
  # Load species-level taxonomy
  tax_species <- read.table(
    tax_species_file,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  message(paste0("  Species taxonomy has ", nrow(tax_species), " ASVs"))
  
  # Filter to keep only ASVs in feature table
  tax_genus_filt <- tax_genus[tax_genus$ASV_ID %in% kept_asvs, ]
  tax_species_filt <- tax_species[tax_species$ASV_ID %in% kept_asvs, ]
  
  message(paste0("  After intersection: ", nrow(tax_genus_filt), " ASVs"))
  
  # Check if all feature table ASVs have taxonomy
  missing_tax <- setdiff(kept_asvs, tax_genus_filt$ASV_ID)
  if (length(missing_tax) > 0) {
    warning(paste0("  WARNING: ", length(missing_tax), " ASVs in feature table have no taxonomy"))
  }
  
  # ===========================
  # 3. Create ASV mapping (hash to ASV_X)
  # ===========================
  message("\nCreating ASV renaming scheme...")
  asv_mapping <- data.frame(
    hash_id = kept_asvs,
    new_id = paste0("ASV_", seq_along(kept_asvs)),
    stringsAsFactors = FALSE
  )
  message(paste0("  Renamed ", nrow(asv_mapping), " ASVs: ASV_1 to ASV_", nrow(asv_mapping)))
  
  # Apply renaming to feature table
  rownames(feature_df) <- asv_mapping$new_id[match(rownames(feature_df), asv_mapping$hash_id)]
  
  # Apply renaming to taxonomy tables
  tax_genus_filt$ASV_ID_original <- tax_genus_filt$ASV_ID
  tax_genus_filt$ASV_ID <- asv_mapping$new_id[match(tax_genus_filt$ASV_ID, asv_mapping$hash_id)]
  
  tax_species_filt$ASV_ID_original <- tax_species_filt$ASV_ID
  tax_species_filt$ASV_ID <- asv_mapping$new_id[match(tax_species_filt$ASV_ID, asv_mapping$hash_id)]
  
  # ===========================
  # 4. Create combined taxonomy table
  # ===========================
  message("\nCreating combined taxonomy table...")
  
  # Use genus-level taxonomy as base (confidence for genus)
  tax_combined <- tax_genus_filt %>%
    select(ASV_ID, Kingdom, Phylum, Class, Order, Family, Genus, confidence_genus = confidence)
  
  # Add species information from species-level taxonomy
  tax_species_select <- tax_species_filt %>%
    select(ASV_ID, Species, Species_exact, confidence_species = confidence)
  
  tax_combined <- left_join(tax_combined, tax_species_select, by = "ASV_ID")
  
  # Remove GTDB prefixes
  tax_cols <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  for (col in tax_cols) {
    if (col %in% colnames(tax_combined)) {
      tax_combined[[col]] <- gsub("^[a-z]__", "", tax_combined[[col]])
    }
  }
  
  # Create final taxonomy table for phyloseq
  tax_final <- tax_combined %>%
    select(ASV_ID, Kingdom, Phylum, Class, Order, Family, Genus, Species, 
           confidence_genus, confidence_species) %>%
    as.data.frame()
  
  rownames(tax_final) <- tax_final$ASV_ID
  tax_final$ASV_ID <- NULL
  
  # ===========================
  # 5. Load metadata
  # ===========================
  message("\nLoading metadata...")
  metadata <- read.csv(metadata_file, row.names = 1, check.names = FALSE)
  
  # Rename columns to match expected names
  if ("bristol_stool_scale" %in% colnames(metadata)) {
    colnames(metadata)[colnames(metadata) == "bristol_stool_scale"] <- "BSS"
    message("  Renamed 'bristol_stool_scale' to 'BSS'")
  }
  
  # Remove sample prefix from rownames if present
  if (!is.null(sample_prefix) && sample_prefix != "") {
    rownames(metadata) <- gsub(paste0("^", sample_prefix), "", rownames(metadata))
  }
  
  # Keep only samples in feature table
  common_samples <- intersect(rownames(metadata), colnames(feature_df))
  metadata <- metadata[common_samples, , drop = FALSE]
  feature_df <- feature_df[, common_samples, drop = FALSE]
  
  # Filter out samples with missing key variables (BSS, age, sex, locality)
  key_vars <- c("BSS", "age", "sex", "locality")
  available_vars <- key_vars[key_vars %in% colnames(metadata)]
  
  if (length(available_vars) > 0) {
    complete_samples <- complete.cases(metadata[, available_vars, drop = FALSE])
    n_removed <- sum(!complete_samples)
    
    if (n_removed > 0) {
      message(paste0("  Removing ", n_removed, " samples with missing values in: ", 
                     paste(available_vars, collapse = ", ")))
      metadata <- metadata[complete_samples, , drop = FALSE]
      feature_df <- feature_df[, rownames(metadata), drop = FALSE]
    }
  }
  
  message(paste0("  Metadata has ", nrow(metadata), " samples"))
  message(paste0("  ", ncol(metadata), " metadata variables"))
  
  # ===========================
  # 6. Load tree
  # ===========================
  message("\nLoading phylogenetic tree...")
  tree <- read.tree(tree_file)
  message(paste0("  Tree has ", length(tree$tip.label), " tips"))
  
  # Rename tree tips
  tree$tip.label <- asv_mapping$new_id[match(tree$tip.label, asv_mapping$hash_id)]
  
  # Check tree matches ASVs
  tree_in_data <- sum(tree$tip.label %in% rownames(feature_df))
  message(paste0("  ", tree_in_data, " tree tips match ASVs in data"))
  
  # ===========================
  # 7. Create phyloseq object
  # ===========================
  message("\nCreating phyloseq object...")
  
  # Create phyloseq components
  otu_mat <- as.matrix(feature_df)
  tax_mat <- as.matrix(tax_final)
  
  OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)
  TAX <- tax_table(tax_mat)
  META <- sample_data(metadata)
  TREE <- phy_tree(tree)
  
  # Create phyloseq object
  ps <- phyloseq(OTU, TAX, META, TREE)
  
  message("\nPhyloseq object summary:")
  print(ps)
  
  # ===========================
  # 8. Save
  # ===========================
  message(paste0("\nSaving phyloseq object to: ", output_file))
  saveRDS(ps, output_file)
  
  # Also save ASV mapping
  mapping_file <- gsub("\\.rds$", "_asv_mapping.tsv", output_file)
  write.table(asv_mapping, mapping_file, sep = "\t", row.names = FALSE, quote = FALSE)
  message(paste0("Saved ASV mapping to: ", mapping_file))
  
  message("\n========================================")
  message("Phyloseq preparation complete!")
  message("========================================\n")
  
  return(ps)
}


# ===========================
# Run if called directly
# ===========================
if (sys.nframe() == 0) {
  # Example usage
  ps <- prepare_phyloseq(
    feature_table_file = "../argentina/data/feature-table.tsv",
    tax_genus_file = "../argentina/data/ASV_tax.gtdb.tsv",
    tax_species_file = "../argentina/data/ASV_tax_species.gtdb.tsv",
    metadata_file = "../argentina/data/metadata.csv",
    tree_file = "../argentina/data/tree.nwk",
    sample_prefix = "S_",
    output_file = "../argentina/phyloseq_object_ampliseq.rds"
  )
}
