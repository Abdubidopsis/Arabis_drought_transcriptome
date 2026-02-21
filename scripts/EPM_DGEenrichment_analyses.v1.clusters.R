# Load necessary libraries
library(tidyverse)
# library(tidygraph) # Not strictly needed if only using igraph directly later
library(dplyr)
library(tidyr)
library(igraph)
library(ggraph)
library(ggplot2)
library(ggrepel)
library(topGO)
library(org.At.tair.db) # Make sure this is installed: BiocManager::install("org.At.tair.db")
library(viridisLite)
library(RColorBrewer)
library(pheatmap)
library(stringr) # Added for str_remove
library(tibble) # For column_to_rownames

#########################
# --- Define Dynamic Parameters ---
species_tag <- "Asag_to_Asag"      # e.g., "Asag", "Athaliana"
epm_tag <- "SW"
condition_tag <- "wilting" # e.g., "recovery", "stress_treatment"
output_dir <- "../studies"  # Directory to save results
p_value_threshold <- 0.00001 # Significance threshold for filtering GO terms
z_threshold <- 3          # Significance threshold for filtering Clusters in DGE enrichment comparison vs NS
cluster_file_k <- 25      # Specify k for cluster file name

##################################################################################
# Set the working directory (ensure this path is correct for your system)
# setwd("/home/ibg-4/Desktop/arab_env_2024i/workflows") # Uncomment if needed
# If running interactively, you might not need setwd() if paths are relative to project root

# Define the file paths (ensure these paths are correct)
file_paths_file_EPM <- c(
  "../studies/AnemAsag_EPMclu/AsagSW_E-04_2025d_gene_none-q1q9_SZ.csv")
#  "~/Desktop/Rhome/creome/assays/tools/protocols/moca/results/Anem_SWX0.8dP7K_Asag/annotation/annotated_motifs_q1q9.csv")
file_paths_file_Q0S <- c(
  "../studies/AnemAsag_dge_quadrants_khan/Asag-D0vW-dge1_reg0.csv")
file_paths_file_ortho <- c(
  "../studies/AnemAsag_dge_quadrants_khan/arabis_sag_hlx2024h.pep.tsv")
# *** NEW: Path to the EPM cluster assignment file ***
file_path_clusters <- paste0("~/Desktop/arab_env_2024i/studies/AnemAsag_EPM_set_analyses_K5K9/EPM_cluster_assignments_k", cluster_file_k, ".csv") # Adjust path as needed


# --- File Reading ---
cat("Reading input files...\n")
tryCatch({
  epm_df <- read.csv(file_paths_file_EPM)
  Q0S_df <- read.csv(file_paths_file_Q0S)
  ortho_df <- read.csv(file_paths_file_ortho, sep = "\t")
  # *** NEW: Read cluster assignment file ***
  cluster_df <- read.csv(file_path_clusters)
  cat("Files read successfully.\n")
}, error = function(e) {
  stop("Error reading input files: ", e$message)
})

###################################################################################
# --- Data Preparation: EPM, Clusters, and DGE Merging ---
cat("Preparing EPM, Cluster, and DGE data...\n")

# *** NEW: Clean Cluster Data ***
# Keep only the first 17 characters of epm and ensure cluster is factor/character
cluster_df_cleaned <- cluster_df %>%
  mutate(
    epm_short = substr(epm, 1, 17), # Extract first 17 characters
    cluster = as.factor(cluster)     # Treat cluster ID as a factor
  ) %>%
  dplyr::select(epm_short, cluster) %>%
  distinct() # Ensure unique short EPM to cluster mapping

cat("Cleaned cluster assignments (head):\n")
print(head(cluster_df_cleaned))

# Select EPM and gene_id, shorten EPM name
epm_df0 <- epm_df %>%
  dplyr::select(epm, gene_id) %>%
  mutate(epm_short = substr(epm, 1, 17)) # Extract first 17 characters

# *** NEW: Merge EPM data with Clusters ***
epm_cluster_df <- epm_df0 %>%
  left_join(cluster_df_cleaned, by = "epm_short") %>%
  filter(!is.na(cluster)) # Keep only EPMs that were assigned to a cluster

# Check how many EPMs were lost due to missing cluster assignment
original_epms <- length(unique(epm_df0$epm_short))
epms_with_clusters <- length(unique(epm_cluster_df$epm_short))
cat("Original unique short EPMs:", original_epms, "\n")
cat("Unique short EPMs with cluster assignments:", epms_with_clusters, "\n")
if (original_epms > epms_with_clusters) {
  cat("WARNING:", original_epms - epms_with_clusters, "EPMs did not have a cluster assignment and were removed.\n")
}

head(epm_cluster_df)
head(Q0S_df)
# Merge clustered EPMs with DGE data
merged_df_Q0S_clustered <- epm_cluster_df %>%
  inner_join(Q0S_df, by = c("gene_id" = "gene_id"), relationship = "many-to-many") %>%
  # Ensure regulation column is character and handle potential NAs/empty strings
  mutate(regulation = as.character(regulation)) %>%
  filter(!is.na(regulation), regulation != "", regulation %in% c("UP", "DO", "NS")) %>%
  # Add a direction column consistent with later GO analysis (UP/DO)
  mutate(direction = ifelse(regulation %in% c("UP", "DO"), regulation, NA_character_))

cat("Head of merged Cluster-DGE data:\n")
# Select relevant columns for head display
print(head(dplyr::select(merged_df_Q0S_clustered, cluster, epm, gene_id, regulation, direction)))

################################################################################
# --- Cluster Enrichment Analysis (Based on DGE Ratios vs NS) ---
# Purpose: Identify CLUSTERS showing significant differences in the proportion of
#          UP- or DO-regulated genes compared to NS-regulated genes.
#          The genes from these significant CLUSTERS will be used for GO analysis.
cat("\n--- Starting CLUSTER Enrichment Analysis (based on DGE Ratios vs NS) ---\n")

# Step 1: Calculate counts per CLUSTER for EACH regulation status (Counts of unique genes)
cluster_gene_counts_by_regulation <- merged_df_Q0S_clustered %>%
  group_by(cluster, regulation) %>%
  summarize(
    unique_genes = n_distinct(gene_id), # Unique genes in this Cluster/regulation combo
    .groups = 'drop'
  )

# Step 1b: Calculate total unique genes per cluster
total_unique_genes_per_cluster <- merged_df_Q0S_clustered %>%
  group_by(cluster) %>%
  summarise(total_unique_genes = n_distinct(gene_id), .groups = 'drop')

# Step 2: Pivot unique gene counts to wide format
cluster_counts_wide <- cluster_gene_counts_by_regulation %>%
  pivot_wider(names_from = regulation,
              values_from = unique_genes,
              names_prefix = "unique_genes_", # Add prefix for clarity
              values_fill = 0) # Fill missing regulation statuses with 0 genes

# Step 2b: Join total unique genes per cluster
# Ensure cluster column types match if necessary before joining
cluster_enrich <- cluster_counts_wide %>%
  left_join(total_unique_genes_per_cluster, by = "cluster") %>%
  # Ensure total_unique_genes is not NA (e.g., if a cluster had 0 genes overall)
  # Also handle cases where UP/DO/NS columns might be missing after pivot if no genes existed
  mutate(
    unique_genes_UP = if("unique_genes_UP" %in% names(.)) unique_genes_UP else 0,
    unique_genes_DO = if("unique_genes_DO" %in% names(.)) unique_genes_DO else 0,
    unique_genes_NS = if("unique_genes_NS" %in% names(.)) unique_genes_NS else 0,
    total_unique_genes = replace_na(total_unique_genes, 0)
  )

cat("Counts per Cluster (UP, DO, NS) and Total:\n")
print(head(cluster_enrich))
# print(colnames(cluster_enrich)) # Optional: Check column names for debugging

# Step 3: Calculate Rates, Variances, and Z-Scores
cluster_enrich <- cluster_enrich %>%
  mutate(
    # Calculate rates using the pivoted counts and total unique genes
    rate_UP = ifelse(total_unique_genes > 0, unique_genes_UP / total_unique_genes, 0),
    rate_DO = ifelse(total_unique_genes > 0, unique_genes_DO / total_unique_genes, 0),
    rate_NS = ifelse(total_unique_genes > 0, unique_genes_NS / total_unique_genes, 0),
    
    # Now calculate variances using the just-calculated rates and total_unique_genes
    var_rate_UP = ifelse(total_unique_genes > 0, rate_UP * (1 - rate_UP) / total_unique_genes, 0),
    var_rate_DO = ifelse(total_unique_genes > 0, rate_DO * (1 - rate_DO) / total_unique_genes, 0),
    var_rate_NS = ifelse(total_unique_genes > 0, rate_NS * (1 - rate_NS) / total_unique_genes, 0),
    
    # Calculate Z-scores, ensuring variances are positive for sqrt
    z_UP = ifelse(total_unique_genes > 0 & (var_rate_UP + var_rate_NS) > 1e-9,
                  (rate_UP - rate_NS) / sqrt(var_rate_UP + var_rate_NS),
                  NA_real_),
    z_DO = ifelse(total_unique_genes > 0 & (var_rate_DO + var_rate_NS) > 1e-9,
                  (rate_DO - rate_NS) / sqrt(var_rate_DO + var_rate_NS),
                  NA_real_)
  ) %>%
  # Flag significance using the Z-threshold
  mutate(
    sig_UP_Z = ifelse(!is.na(z_UP) & abs(z_UP) > z_threshold, TRUE, FALSE),
    sig_DO_Z = ifelse(!is.na(z_DO) & abs(z_DO) > z_threshold, TRUE, FALSE)
  )


print("Data with Rates, Z-Scores (DGE Rates vs NS background):")
print(dplyr::select(cluster_enrich, cluster, starts_with("unique_genes_"), total_unique_genes, starts_with("rate_"), z_UP, sig_UP_Z, z_DO, sig_DO_Z) %>% head())


# Step 4: Create Grouping Variable for Plotting (Based on Z-Score Significance vs NS)
# (This step should now work correctly as z_UP, z_DO etc. are calculated)
cluster_enrich <- cluster_enrich %>%
  mutate(
    sig_UP_Z_plot = ifelse(is.na(sig_UP_Z), FALSE, sig_UP_Z),
    sig_DO_Z_plot = ifelse(is.na(sig_DO_Z), FALSE, sig_DO_Z),
    sig_combined_z = case_when(
      sig_UP_Z_plot & sig_DO_Z_plot ~ "Both Significant",
      sig_UP_Z_plot                ~ "UP Significant",
      sig_DO_Z_plot                ~ "DO Significant",
      TRUE                         ~ "Not Significant"
    ),
    cluster_group_plot = ifelse(sig_UP_Z_plot | sig_DO_Z_plot, as.character(cluster), "Not Significant")
  )

# Define colors for plot - adjust as needed
plot_colors <- c("UP Significant" = "red", "DO Significant" = "blue",
                 "Both Significant" = "purple", "Not Significant" = "grey")

# (Inside the DGE Cluster Enrichment section, Step 5)

# Step 5: Create the Z-score Plot (DGE Rates vs NS)
cluster_enrich_plot <- cluster_enrich %>%
  mutate(z_UP_plot = replace_na(z_UP, 0),
         z_DO_plot = replace_na(z_DO, 0))

if(nrow(cluster_enrich_plot) > 0) {
  p_dge_enrich <- ggplot(cluster_enrich_plot, aes(x = z_UP_plot, y = z_DO_plot, color = sig_combined_z)) +
    geom_point(aes(size = total_unique_genes), alpha = 0.7) +
    # *** MODIFIED: Increased text size ***
    scale_size_continuous(name = "Total Unique Genes", range = c(3, 10)) + # Slightly larger points
    geom_text_repel(
      data = filter(cluster_enrich_plot, sig_combined_z != "Not Significant"),
      aes(label = cluster),
      # *** MODIFIED: Increased text size ***
      size = 4, color = "black", # Label significant points in black (Increased size)
      box.padding = 0.6, point.padding = 0.7, segment.color = "grey50",
      max.overlaps = 30 # Increased max.overlaps
    ) +
    geom_hline(yintercept = c(-z_threshold, z_threshold), linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = c(-z_threshold, z_threshold), linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = 0, linetype = "solid", linewidth = 0.4, color = "black") + # Slightly thicker line
    geom_vline(xintercept = 0, linetype = "solid", linewidth = 0.4, color = "black") + # Slightly thicker line
    scale_color_manual(values = plot_colors, name = "Significance vs NS") +
    labs(title = "Cluster Enrichment based on DGE Rates (UP/DO vs NS)",
         subtitle = paste("Z-Scores comparing gene rate in UP/DO vs NS. Threshold = |", z_threshold, "|"),
         x = "Z-Score (UP Rate vs NS Rate)",
         y = "Z-Score (DO Rate vs NS Rate)") +
    # *** MODIFIED: Increased base font size ***
    theme_minimal(base_size = 14) + # Increased base font size for theme elements
    theme(panel.grid.minor = element_blank(),
          legend.position = "bottom",
          plot.title = element_text(size = rel(1.1)), # Relative size increase for title
          plot.subtitle = element_text(size = rel(1.0)),
          axis.title = element_text(size = rel(1.0)),
          legend.title = element_text(size = rel(1.0)),
          legend.text = element_text(size = rel(0.9))
    )
  
  print(p_dge_enrich)
  
  # Save plot
  if (!dir.exists(output_dir)) {dir.create(output_dir, recursive = TRUE)}
  dge_enrich_plot_file <- file.path(output_dir, paste0("CLUSTER_DGErate_enrichment_", species_tag, "_", epm_tag, "_", condition_tag, "_ZscorePlot_vs_NS_k", cluster_file_k, ".png"))
  # *** MODIFIED: Increased DPI and adjusted dimensions ***
  ggsave(filename = dge_enrich_plot_file, plot = p_dge_enrich, width = 8, height = 7.5, units = "in", dpi = 600)
  cat("Saved DGE rate enrichment plot (High Res) to:", dge_enrich_plot_file, "\n")
  
  # *** NEW: Optional PDF output ***
  dge_enrich_plot_file_pdf <- file.path(output_dir, paste0("CLUSTER_DGErate_enrichment_", species_tag, "_", epm_tag, "_", condition_tag, "_ZscorePlot_vs_NS_k", cluster_file_k, ".pdf"))
  ggsave(filename = dge_enrich_plot_file_pdf, plot = p_dge_enrich, width = 8, height = 7.5, units = "in", device = cairo_pdf)
  cat("Saved DGE rate enrichment plot (PDF) to:", dge_enrich_plot_file_pdf, "\n")
  
  
} else {
  print("No data available to generate the DGE rate Z-score plot.")
}
# Step 6: Select and Save Genes from Significantly Enriched CLUSTERS (Based on DGE Ratios vs NS)
clusters_enriched_by_dge <- cluster_enrich %>%
  filter(sig_UP_Z | sig_DO_Z) %>% # Clusters where UP or DO rate is significantly different from NS
  distinct(cluster, .keep_all = TRUE)

# Get the subset of original DGE results (UP/DO genes only) for these significant CLUSTERS
genes_from_dge_enriched_clusters <- merged_df_Q0S_clustered %>%
  # *** MODIFIED: Filter by significant clusters ***
  filter(cluster %in% clusters_enriched_by_dge$cluster,
         !is.na(direction))             # Only UP or DO regulated genes

# *** MODIFIED: Select cluster, gene_id, direction ***
cluster_DGEenriched_subset_for_GO <- genes_from_dge_enriched_clusters %>%
  dplyr::select(cluster, gene_id, direction) # Select relevant columns for GO analysis input

# Print summaries
# *** MODIFIED: Print significant clusters ***
cat("\nSignificant Clusters (based on DGE rate vs NS):\n")
if(nrow(clusters_enriched_by_dge)>0) print(dplyr::select(clusters_enriched_by_dge, cluster, z_UP, sig_UP_Z, z_DO, sig_DO_Z)) else print("None found.")
cat("Total genes from these significant clusters passed to GO analysis:", nrow(cluster_DGEenriched_subset_for_GO), "\n")
cat("Head of genes passed to GO analysis:\n")
if(nrow(cluster_DGEenriched_subset_for_GO)>0) print(head(cluster_DGEenriched_subset_for_GO)) else print("None")

# Save intermediate results
# *** MODIFIED: Filenames ***
dge_enrich_full_stats_file <- file.path(output_dir, paste0("clusters_DGErate_enrichment_", species_tag, "_",epm_tag, "_", condition_tag, "_full_stats_vs_NS_k", cluster_file_k, ".csv"))
dge_enrich_sig_clusters_file <- file.path(output_dir, paste0("clusters_DGErate_enrichment_", species_tag, "_",epm_tag, "_", condition_tag, "_significant_vs_NS_k", cluster_file_k, ".csv"))
dge_enrich_gene_list_file <- file.path(output_dir, paste0("genes_from_DGErate_enriched_clusters_", species_tag,epm_tag, "_", "_", condition_tag, "_for_GO_k", cluster_file_k, ".csv"))

write.csv(cluster_enrich, file = dge_enrich_full_stats_file, row.names = FALSE)
write.csv(clusters_enriched_by_dge, file = dge_enrich_sig_clusters_file, row.names = FALSE)
write.csv(cluster_DGEenriched_subset_for_GO, file = dge_enrich_gene_list_file, row.names = FALSE)
cat("Saved DGE rate enrichment results and gene list for GO.\n")

cat("--- Finished CLUSTER Enrichment Analysis (based on DGE Ratios vs NS) ---\n")

#################################################################################
# --- Ortholog Mapping ---
cat("\n--- Preparing Ortholog Mapping ---\n")

# Rename ortho_df columns for clarity
new_column_names <- c("Orthogroup", "Species", "target_gene_id", "At_orthologs")
if (length(new_column_names) == ncol(ortho_df)) {
  colnames(ortho_df) <- new_column_names
} else {
  warning(paste("Ortholog file column count mismatch. Expected", length(new_column_names), "got", ncol(ortho_df)))
}

# Filter for Arabidopsis and create long format
ortho_long_df <- ortho_df %>%
  filter(Species == "Arabidopsis_thaliana.TAIR10.59.pep.all") %>% # Ensure this species name is exact
  filter(!is.na(target_gene_id) & target_gene_id != "",
         !is.na(At_orthologs) & At_orthologs != "") %>% # Filter out empty strings before separating
  separate_rows(target_gene_id, sep = ",\\s*") %>%
  separate_rows(At_orthologs, sep = ",\\s*") %>%
  dplyr::select(Orthogroup, target_gene_id, At_orthologs) # Keep relevant columns

# Clean ortholog data: remove version suffixes and duplicates
ortho_cleaned_df <- ortho_long_df %>%
  dplyr::select(target_gene_id, At_orthologs) %>%
  mutate(across(c(target_gene_id, At_orthologs), ~ str_remove(.x, "\\.\\d+$"))) %>%
  filter(target_gene_id != "", At_orthologs != "") %>% # Filter again after removing suffix
  distinct()

cat("Cleaned ortholog mapping prepared.\n")
print(head(ortho_cleaned_df))

##################################################################################
# --- Data Preparation for GO Analysis ---
cat("\n--- Preparing Data for GO Analysis ---\n")

# Merge gene list (from DGE-rate-enriched CLUSTERS) with cleaned orthologs
# *** MODIFIED: Use the cluster-based gene list ***
merged_data_for_GO <- cluster_DGEenriched_subset_for_GO %>%
  left_join(ortho_cleaned_df,
            by = c("gene_id" = "target_gene_id"),
            relationship = "many-to-many") %>%
  filter(!is.na(At_orthologs) & At_orthologs != "") # Keep only rows with a matched Arabidopsis ortholog

cat("Head of data prepared for GO analysis (genes from DGE-enriched CLUSTERS + orthologs):\n")
# *** MODIFIED: Show cluster column ***
print(head(dplyr::select(merged_data_for_GO, cluster, gene_id, direction, At_orthologs)))

# Define the universe of genes (all unique Arabidopsis orthologs involved)
# Using all orthologs from the initial merged_df_Q0S_clustered provides a broader background
all_orthologs_in_data <- merged_df_Q0S_clustered %>%
  # Make sure cluster info is present if needed, but for universe it's just the genes
  inner_join(ortho_cleaned_df, by = c("gene_id" = "target_gene_id")) %>%
  pull(At_orthologs) %>%
  unique() %>%
  na.omit() %>%
  .[. != ""]

universe_genes <- all_orthologs_in_data # Using broader background
cat("\nTotal unique Arabidopsis orthologs in GO Universe:", length(universe_genes), "\n")

# Get list of unique CLUSTER conditions to test (from the subsetted data)
# *** MODIFIED: Get unique clusters ***
cluster_list_for_GO <- unique(merged_data_for_GO$cluster)
cat("Cluster conditions to be tested in GO:", paste(sort(cluster_list_for_GO), collapse=", "), "\n\n")

##################################################################################
# --- GO Enrichment Function ---
# *** MODIFIED: Function adapted for clusters ***
run_go_enrichment_cluster <- function(data, cluster_conditions, universe, direction_filter,
                                      ontology = "BP", node_size = 10,
                                      algorithm = "elim", statistic = "fisher", top_nodes = 50) {
  
  go_results_list <- list()
  cat("Running GO Enrichment for Direction:", direction_filter, "\n")
  
  # *** MODIFIED: Loop through clusters ***
  for (cl in cluster_conditions) {
    cat("  Processing Cluster:", as.character(cl), "...")
    # Subset genes of interest for this CLUSTER and direction
    genes_interest_cl <- data %>%
      # *** MODIFIED: Filter by cluster ***
      filter(cluster == cl, direction == direction_filter) %>%
      pull(At_orthologs) %>%
      unique()
    
    # Skip if no genes found for this cluster/direction
    if (length(genes_interest_cl) == 0) {
      cat(" No genes found. Skipping.\n")
      next
    }
    # Ensure genes of interest are within the defined universe
    genes_interest_cl <- intersect(genes_interest_cl, universe)
    if (length(genes_interest_cl) == 0) {
      cat(" No genes found *within universe*. Skipping.\n")
      next
    }
    cat(" Found", length(genes_interest_cl), "genes within universe.")
    
    # Create named factor geneList required by topGO
    geneList <- factor(as.integer(universe %in% genes_interest_cl))
    names(geneList) <- universe
    
    # Create topGO data object
    suppressMessages({ # Suppress numerous messages from topGO
      GOdata <- tryCatch(new("topGOdata",
                             # *** MODIFIED: Description ***
                             description = paste("GO Enrichment for Cluster", cl, direction_filter),
                             ontology = ontology,
                             allGenes = geneList,
                             geneSel = function(x) x == 1,
                             nodeSize = node_size,
                             mapping = "org.At.tair.db",
                             annot = annFUN.org), # Requires org.At.tair.db
                         error = function(err) {
                           cat(" ERROR creating topGOdata for Cluster", cl, ":", err$message, "\n")
                           return(NULL)
                         })
    })
    
    if(is.null(GOdata)) next # Skip if GOdata creation failed
    
    # Run enrichment test
    resultFisher <- tryCatch(runTest(GOdata, algorithm = algorithm, statistic = statistic),
                             error = function(err){
                               cat(" ERROR running runTest for Cluster", cl, ":", err$message, "\n")
                               return(NULL)
                             })
    
    if(is.null(resultFisher)) next # Skip if runTest failed
    
    # Generate results table
    goEnrichment <- tryCatch(GenTable(GOdata, Fisher = resultFisher, orderBy = "Fisher", topNodes = top_nodes),
                             error = function(err){
                               cat(" ERROR running GenTable for Cluster", cl, ":", err$message, "\n")
                               return(NULL)
                             })
    
    if(is.null(goEnrichment) || nrow(goEnrichment) == 0) {
      cat (" No significant terms found or error in GenTable.\n")
      next # Skip if GenTable failed or returned empty
    }
    
    # *** MODIFIED: Add Cluster identifier and store ***
    goEnrichment$Cluster <- cl # Add cluster ID
    go_results_list[[as.character(cl)]] <- goEnrichment # Use cluster ID as list name
    cat(" Found", nrow(goEnrichment), "enriched terms.\n")
    
  } # End loop over cluster_conditions
  
  # Combine results from the list into a single dataframe
  combined_results <- bind_rows(go_results_list)
  # Ensure Fisher is numeric, handle potential "< 1e-30" strings
  if (nrow(combined_results) > 0 && "Fisher" %in% names(combined_results)) {
    combined_results <- combined_results %>%
      mutate(Fisher = suppressWarnings(as.numeric(str_replace(Fisher, "<\\s*", ""))))
  }
  # Ensure Cluster column is of the correct type (factor or character as needed)
  if(nrow(combined_results) > 0 && "Cluster" %in% names(combined_results)) {
    # Convert Cluster back to factor if it was numeric in the list names
    combined_results <- combined_results %>% mutate(Cluster = factor(Cluster, levels = levels(cluster_conditions)))
  }
  
  
  cat("Finished GO Enrichment for Direction:", direction_filter, "\n\n")
  return(combined_results)
}

##################################################################################
# --- Run GO Enrichment for UP and DO Regulated Genes (per Cluster) ---
cat("--- Running GO Enrichment per Cluster ---\n")
# *** MODIFIED: Run on cluster data and use cluster function ***
go_results_cluster_UP <- run_go_enrichment_cluster(merged_data_for_GO, cluster_list_for_GO, universe_genes, "UP")
go_results_cluster_DO <- run_go_enrichment_cluster(merged_data_for_GO, cluster_list_for_GO, universe_genes, "DO")

# --- Save Raw GO Results ---
# *** MODIFIED: Filenames ***
go_up_raw_file <- file.path(output_dir, paste0("GO_cluster_enrichment_raw_", species_tag, "_", condition_tag, "_UP_k", cluster_file_k, ".csv"))
go_do_raw_file <- file.path(output_dir, paste0("GO_cluster_enrichment_raw_", species_tag, "_", condition_tag, "_DO_k", cluster_file_k, ".csv"))
write.csv(go_results_cluster_UP, file = go_up_raw_file, row.names = FALSE)
write.csv(go_results_cluster_DO, file = go_do_raw_file, row.names = FALSE)
cat("Saved raw GO results to:\n", go_up_raw_file, "\n", go_do_raw_file, "\n\n")

##################################################################################
# --- Analysis 1: Identify Significantly Enriched GO Terms (Per CLUSTER) ---
# *** MODIFIED: Analysis per Cluster ***
cat("--- Analysis 1: Identifying Significant GO Terms per Cluster (p <", p_value_threshold, ") ---\n")

filter_significant_go <- function(go_results, p_thresh) {
  if (is.null(go_results) || nrow(go_results) == 0) return(tibble()) # Return empty tibble if no data
  go_results %>%
    filter(!is.na(Fisher) & Fisher < p_thresh) %>%
    mutate(minusLog10Fisher = -log10(Fisher + .Machine$double.eps)) %>%
    # *** MODIFIED: Arrange by Cluster ***
    arrange(Cluster, Fisher)
}

# *** MODIFIED: Use cluster results ***
significant_go_UP <- filter_significant_go(go_results_cluster_UP, p_value_threshold)
significant_go_DO <- filter_significant_go(go_results_cluster_DO, p_value_threshold)

cat("Total significant GO terms found (UP):", nrow(significant_go_UP), "\n")
cat("Total significant GO terms found (DO):", nrow(significant_go_DO), "\n")

# Summarize counts per CLUSTER
# *** MODIFIED: Count by Cluster ***
summary_sig_go_UP <- significant_go_UP %>% count(Cluster, name = "Significant_GO_Count_UP")
summary_sig_go_DO <- significant_go_DO %>% count(Cluster, name = "Significant_GO_Count_DO")

cat("Summary of significant GO term counts per Cluster (UP):\n")
print(summary_sig_go_UP)
cat("Summary of significant GO term counts per Cluster (DO):\n")
print(summary_sig_go_DO)

# Save the lists of significant GO terms
# *** MODIFIED: Filenames ***
sig_go_up_file <- file.path(output_dir, paste0("GO_cluster_enrichment_significant_", species_tag, "_", condition_tag, "_UP_k", cluster_file_k, ".csv"))
sig_go_do_file <- file.path(output_dir, paste0("GO_cluster_enrichment_significant_", species_tag, "_", condition_tag, "_DO_k", cluster_file_k, ".csv"))
write.csv(significant_go_UP, file = sig_go_up_file, row.names = FALSE)
write.csv(significant_go_DO, file = sig_go_do_file, row.names = FALSE)
cat("Saved significant GO terms lists to:\n", sig_go_up_file, "\n", sig_go_do_file, "\n\n")

cat("--- Finished Analysis 1 ---\n")

##################################################################################
# --- Analysis 2: Identify CLUSTERS Enriched for Significant GO Terms ---
# *** MODIFIED: Analysis of Clusters ***
cat("\n--- Analysis 2: Identifying CLUSTERS Enriched for Significant GO Terms ---\n")

# *** MODIFIED: Function calculates enrichment for Clusters ***
calculate_cluster_go_enrichment <- function(significant_go_results, direction_label) {
  if (is.null(significant_go_results) || nrow(significant_go_results) == 0) {
    cat("No significant GO terms found for", direction_label, ", skipping Cluster enrichment calculation.\n")
    return(tibble(Cluster = factor(), Significant_GO_Count = integer(), Mean_NegLog10Pval = double()))
  }
  
  cluster_go_summary <- significant_go_results %>%
    # *** MODIFIED: Group by Cluster ***
    group_by(Cluster) %>%
    summarise(
      Significant_GO_Count = n(),
      Mean_NegLog10Pval = mean(minusLog10Fisher, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(Significant_GO_Count), desc(Mean_NegLog10Pval)) # Rank by count, then p-value score
  
  # Add rank
  cluster_go_summary <- cluster_go_summary %>%
    mutate(Rank = row_number())
  
  cat("Clusters ranked by enrichment of significant GO terms (", direction_label, "):\n")
  print(cluster_go_summary)
  return(cluster_go_summary)
}

# *** MODIFIED: Use significant cluster results ***
cluster_go_enriched_UP <- calculate_cluster_go_enrichment(significant_go_UP, "UP")
cluster_go_enriched_DO <- calculate_cluster_go_enrichment(significant_go_DO, "DO")

# Save Cluster enrichment results based on GO
# *** MODIFIED: Filenames ***
cluster_go_enrich_up_file <- file.path(output_dir, paste0("Cluster_enrichment_by_GO_", species_tag, "_", condition_tag, "_UP_k", cluster_file_k, ".csv"))
cluster_go_enrich_do_file <- file.path(output_dir, paste0("Cluster_enrichment_by_GO_", species_tag, "_", condition_tag, "_DO_k", cluster_file_k, ".csv"))
write.csv(cluster_go_enriched_UP, file = cluster_go_enrich_up_file, row.names = FALSE)
write.csv(cluster_go_enriched_DO, file = cluster_go_enrich_do_file, row.names = FALSE)
cat("Saved Cluster enrichment rankings based on GO results to:\n", cluster_go_enrich_up_file, "\n", cluster_go_enrich_do_file, "\n\n")

cat("--- Finished Analysis 2 ---\n")

##################################################################################
# --- Network Building Function ---
# *** MODIFIED: Function uses Cluster info for edges ***
build_go_network_cluster <- function(go_results, p_val_thresh) {
  cat("Building GO network based on Clusters...\n")
  
  # Ensure input is a tibble for more consistent dplyr behavior
  # Use the *significant* GO results as input for the network
  go_results_tbl <- as_tibble(go_results) %>%
    filter(!is.na(Fisher) & Fisher < p_val_thresh) %>% # Ensure only significant terms are used
    mutate(minusLog10Fisher = -log10(Fisher + .Machine$double.eps))
  
  if(nrow(go_results_tbl) == 0) {
    cat("  No significant GO terms found below threshold", p_val_thresh, ". Cannot build network.\n")
    return(NULL)
  }
  cat("  Found", nrow(go_results_tbl), "significant term entries for network building.\n")
  
  # *** MODIFIED: Check for Cluster column ***
  required_cols <- c("Cluster", "Term", "minusLog10Fisher", "Annotated")
  if (!all(required_cols %in% names(go_results_tbl))) {
    stop("Missing required columns in input data for build_go_network_cluster: ",
         paste(setdiff(required_cols, names(go_results_tbl)), collapse=", "))
  }
  filtered_go <- go_results_tbl %>%
    mutate(Term = as.character(Term),
           # *** MODIFIED: Ensure Cluster is character or factor for grouping ***
           Cluster = as.character(Cluster), # Convert factor to character for grouping/edges
           Annotated = as.integer(Annotated)) %>% # Ensure Annotated is integer
    filter(!is.na(Term) & Term != "", !is.na(Cluster) & Cluster != "", !is.na(Annotated))
  
  
  # Build edge list based on shared CLUSTERS
  edge_list <- filtered_go %>%
    # *** MODIFIED: Group by Cluster ***
    group_by(Cluster) %>%
    filter(n() > 1) %>%
    summarise(
      pairs = list(combn(Term, 2, simplify = FALSE)), # Term is already character
      weight = mean(minusLog10Fisher, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(lengths(pairs) > 0) %>% # Ensure combn produced pairs
    unnest(pairs) %>%
    filter(lengths(pairs) == 2) %>%
    mutate(
      Term1 = sapply(pairs, `[`, 1),
      Term2 = sapply(pairs, `[`, 2)
    ) %>%
    # *** MODIFIED: Select Cluster as edge attribute ***
    dplyr::select(from = Term1, to = Term2, Cluster, weight)
  
  if(nrow(edge_list) == 0) {
    cat("  No edges could be formed (no Cluster groups with >1 significant term).\n")
  } else {
    cat("  Built", nrow(edge_list), "edges.\n")
  }
  
  # Build node table (remains the same - based on unique GO Terms)
  node_df <- filtered_go %>%
    group_by(Term) %>%
    arrange(desc(Annotated)) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    mutate(name = trimws(Term)) %>% # Ensure 'name' column exists
    dplyr::select(name, Annotated)
  
  if(nrow(node_df) == 0 && nrow(edge_list) == 0) {
    cat("  No nodes could be determined. Cannot build network.\n")
    return(NULL)
  } else if (nrow(node_df) == 0 && nrow(edge_list) > 0) {
    all_terms_in_edges <- unique(c(edge_list$from, edge_list$to))
    node_df <- tibble(name = all_terms_in_edges, Annotated = 1)
    cat("  Created node list from edges as initial node list was empty.\n")
  }
  
  # Ensure all nodes from edge list are in node list
  if(nrow(edge_list) > 0) {
    all_terms_in_edges <- unique(c(edge_list$from, edge_list$to))
    missing_nodes <- setdiff(all_terms_in_edges, node_df$name)
    if (length(missing_nodes) > 0) {
      dummy_nodes <- tibble(name = missing_nodes, Annotated = 1)
      node_df <- bind_rows(node_df, dummy_nodes)
      cat("  Added", length(missing_nodes), "dummy nodes present in edges but not initially in node list.\n")
    }
  }
  
  # Ensure node names are unique
  node_df <- node_df %>% distinct(name, .keep_all = TRUE)
  cat("  Built node table with", nrow(node_df), "unique terms.\n")
  
  # Create igraph object
  if (nrow(node_df) > 0) {
    if (nrow(edge_list) == 0) {
      cat("  Creating graph with nodes but no edges.\n")
      graph_obj <- make_empty_graph(n = nrow(node_df), directed = FALSE)
      V(graph_obj)$name <- node_df$name
      V(graph_obj)$Annotated <- node_df$Annotated
    } else {
      # *** MODIFIED: Ensure Cluster is factor for coloring if needed ***
      edge_list <- edge_list %>% mutate(Cluster = factor(Cluster))
      graph_obj <- graph_from_data_frame(d = edge_list, vertices = node_df, directed = FALSE)
    }
    cat("Finished building GO network based on Clusters.\n\n")
    return(graph_obj)
  } else {
    cat("  Cannot create graph: No nodes were finalized.\n\n")
    return(NULL)
  }
}

##################################################################################
# --- Build Networks for UP and DO (Based on Clusters) ---
cat("--- Building Cluster-based GO Networks ---\n")
# *** MODIFIED: Use cluster results and cluster network function ***
cat("Processing UP regulated genes...\n")
# Check structure of input
# str(go_results_cluster_UP)
g_UP <- build_go_network_cluster(go_results_cluster_UP, p_value_threshold)

cat("\nProcessing DOWN regulated genes...\n")
# Check structure of input
# str(go_results_cluster_DO)
g_DO <- build_go_network_cluster(go_results_cluster_DO, p_value_threshold)

# Optional: Inspect graphs
# print(g_UP)
# print(g_DO)

##################################################################################
# --- Network Visualization Prep ---
cat("\n--- Preparing Network Visualization ---\n")

# *** MODIFIED: Find unique CLUSTERS in graphs ***
unique_clusters_in_graphs <- character(0)
if (!is.null(g_UP) && "Cluster" %in% edge_attr_names(g_UP)) {
  unique_clusters_in_graphs <- union(unique_clusters_in_graphs, unique(na.omit(E(g_UP)$Cluster)))
}
if (!is.null(g_DO) && "Cluster" %in% edge_attr_names(g_DO)) {
  unique_clusters_in_graphs <- union(unique_clusters_in_graphs, unique(na.omit(E(g_DO)$Cluster)))
}
# Ensure they are characters for naming colors
unique_clusters_in_graphs <- as.character(unique_clusters_in_graphs)


custom_colors <- NULL
if (length(unique_clusters_in_graphs) > 0) {
  tryCatch({
    # Use viridis or other qualitative palette for potentially many clusters
    if(length(unique_clusters_in_graphs) > 8){ # Use viridis if many clusters
      custom_colors <- viridis(length(unique_clusters_in_graphs))
    } else if (length(unique_clusters_in_graphs) > 1) { # Use Brewer for fewer
      # Use a qualitative palette suitable for factors
      qual_palette <- brewer.pal.info[brewer.pal.info$category == 'qual',]
      # Pick a palette with enough colors, e.g., Set1, Set2, Paired
      palette_name <- rownames(qual_palette)[qual_palette$maxcolors >= length(unique_clusters_in_graphs)][1]
      if(is.na(palette_name)) palette_name <- "Set1" # Fallback
      custom_colors <- brewer.pal(length(unique_clusters_in_graphs), palette_name)
    } else { # Only one cluster
      custom_colors <- "blue"
    }
    names(custom_colors) <- unique_clusters_in_graphs
    cat("Defined colors for Clusters:", paste(names(custom_colors), collapse=", "), "\n")
    # print(custom_colors) # Can be long if many clusters
  }, error = function(e) {
    cat("Error generating colors:", e$message, "\nAttempting fallback.\n")
    # Fallback
    custom_colors <- setNames(rainbow(length(unique_clusters_in_graphs)), unique_clusters_in_graphs)
    cat("Using fallback colors for Clusters:", paste(names(custom_colors), collapse=", "), "\n")
    # print(custom_colors)
  })
} else {
  cat("No Clusters found in significant graph edges to define colors for.\n")
}

##################################################################################
# --- Network Plotting Function (Node size by Degree) ---
# *** MODIFIED: Function plots cluster-based network ***
# (Replace the entire plot_go_network_cluster function)

plot_go_network_cluster <- function(graph_obj, direction_label, color_palette,
                                    sp_tag, cond_tag, clust_k, out_dir, layout_type = "kk") {
  
  if (is.null(graph_obj) || length(V(graph_obj)) == 0) {
    cat("Skipping network plot for", direction_label, "- No graph data.\n\n")
    return(NULL)
  }
  cat("Plotting CLUSTER network (Node Size by Degree) for:", direction_label, "\n")
  
  plot_title <- paste("GO Term Network (Edges by Cluster, Size by Degree) -", direction_label, "Genes")
  # *** MODIFIED: Filename base ***
  filename_base <- file.path(out_dir, paste0("Network_Cluster_", sp_tag, "_", cond_tag, "_", direction_label, "_node_size_degree_k", clust_k))
  filename_png <- paste0(filename_base, ".png")
  filename_pdf <- paste0(filename_base, ".pdf") # For optional PDF output
  
  # Calculate node degree
  if (length(V(graph_obj)) > 0) {
    node_degrees <- igraph::degree(graph_obj, mode = "all")
    graph_obj <- set_vertex_attr(graph_obj, "degree", value = node_degrees)
  } else {
    graph_obj <- set_vertex_attr(graph_obj, "degree", value = 0)
  }
  
  has_edges <- length(E(graph_obj)) > 0
  has_cluster_attr <- if(has_edges) "Cluster" %in% edge_attr_names(graph_obj) else FALSE
  
  # Determine node label size based on number of nodes (optional adjustment)
  num_nodes <- length(V(graph_obj))
  # *** MODIFIED: Larger base node label size, slight decrease for very large networks ***
  node_label_size <- max(2.5, 5 - num_nodes * 0.005) # Adjust base size (5) and scaling factor as needed
  
  p <- ggraph(graph_obj, layout = layout_type)
  
  if(has_edges && has_cluster_attr) {
    # *** MODIFIED: Slightly thicker edges ***
    p <- p + geom_edge_link(aes(color = Cluster), linewidth = 1.2, alpha = 0.6)
  } else if (has_edges) {
    p <- p + geom_edge_link(linewidth = 1.2, alpha = 0.6, color="grey")
  }
  
  p <- p + geom_node_point(aes(size = degree), color = "grey40", alpha=0.8) +
    # *** MODIFIED: Significantly increased node text size ***
    geom_node_text(aes(label = name),
                   size = node_label_size, # Use dynamic or fixed larger size
                   color = "black", repel = TRUE, max.overlaps = 25, # Increased max.overlaps further
                   bg.color = "white", bg.r = 0.1) +
    # *** MODIFIED: Adjusted node size range ***
    scale_size_continuous(name = "Node Degree", range = c(4, 15)) + # Increased range
    # *** MODIFIED: Increased base theme size ***
    theme_graph(base_family = "sans", base_size = 14) + # Increased base size
    labs(title = plot_title, edge_color = "Cluster") +
    # *** MODIFIED: Adjust legend text size if needed ***
    theme(
      plot.title = element_text(size = rel(1.1)),
      legend.title = element_text(size = rel(1.0)),
      legend.text = element_text(size = rel(0.9)),
      legend.position = "right" # Adjust legend position if needed
    )
  
  
  if (has_cluster_attr && !is.null(color_palette)) {
    p <- p + scale_edge_color_manual(values = color_palette, na.value = "grey")
  } else if (has_cluster_attr) {
    p <- p + scale_edge_color_viridis_d(na.value = "grey")
  }
  
  
  # Display and Save
  tryCatch({
    print(p)
    cat("Plot display attempted.\n")
    # *** MODIFIED: Increased DPI and adjusted dimensions ***
    ggsave(filename_png, plot = p, width = 14, height = 12, units = "in", dpi = 600, bg = "white") # Larger dimensions
    cat("Saved cluster network plot (High Res PNG) to:", filename_png, "\n")
    
    # *** NEW: Optional PDF output ***
    ggsave(filename_pdf, plot = p, width = 14, height = 12, units = "in", device = cairo_pdf)
    cat("Saved cluster network plot (PDF) to:", filename_pdf, "\n\n")
    
  }, error = function(e) {
    cat("ERROR displaying/saving plot:", e$message, "\n")
  })
  
  return(p)
}
cat("Cluster Network plot function redefined with larger fonts and higher resolution.\n") # Update message

##################################################################################
# --- Plot Networks ---
# Ensure plotting device is active if running non-interactively or after dev.off()
# if(is.null(dev.list())) png(file.path(output_dir,"dummy_plot_start.png")); dev.off() # Opens a device if none are open

cat("\n--- Plotting Cluster-based GO Networks ---\n")
# *** MODIFIED: Use cluster plot function and pass cluster k ***
p_network_UP <- plot_go_network_cluster(g_UP, "UP", custom_colors, species_tag, condition_tag, cluster_file_k, output_dir)
p_network_DO <- plot_go_network_cluster(g_DO, "DO", custom_colors, species_tag, condition_tag, cluster_file_k, output_dir)

##################################################################################
# --- Heatmap Preparation and Plotting ---
cat("\n--- Preparing and Plotting Cluster Heatmaps ---\n")

# Heatmap Data Preparation Function (uses significant GO terms per cluster)
# *** MODIFIED: Prepares data based on Clusters ***
prepare_heatmap_data_cluster <- function(significant_go_results, p_val_thresh) { # Takes significant results now
  cat("Preparing cluster heatmap data...\n")
  
  if(is.null(significant_go_results) || nrow(significant_go_results) == 0) {
    cat("  No significant GO terms provided. Cannot create heatmap data.\n")
    return(NULL)
  }
  
  # Data is already filtered and has minusLog10Fisher
  heatmap_data_wide <- significant_go_results %>%
    # *** MODIFIED: Select Cluster, Term, value ***
    dplyr::select(Cluster, Term, minusLog10Fisher) %>%
    # Ensure unique Term/Cluster pairs - take max p-value if duplicated
    # *** MODIFIED: Group by Cluster, Term ***
    group_by(Cluster, Term) %>%
    summarise(minusLog10Fisher = max(minusLog10Fisher, na.rm = TRUE), .groups = "drop") %>%
    # Handle potential Inf values if p-value was exactly 0 + epsilon
    mutate(minusLog10Fisher = ifelse(is.infinite(minusLog10Fisher),
                                     # Find max non-infinite value or use 0 if none exist
                                     max(c(0, .$minusLog10Fisher[!is.infinite(.$minusLog10Fisher)]), na.rm=TRUE) + 1,
                                     minusLog10Fisher)) %>%
    # *** MODIFIED: Pivot wider with Cluster as rows ***
    pivot_wider(names_from = Term, values_from = minusLog10Fisher, values_fill = 0)
  
  cat("Finished preparing cluster heatmap data matrix. Dimensions:", dim(heatmap_data_wide), "\n\n")
  return(heatmap_data_wide) # Return wide format for pheatmap
}

# Heatmap Plotting Function (pheatmap)
# *** MODIFIED: Plots cluster-based heatmap ***
# (Replace the entire plot_go_heatmap_pheatmap_cluster function)

plot_go_heatmap_pheatmap_cluster <- function(heatmap_matrix_data, direction_label, color_vec,
                                             sp_tag, cond_tag, clust_k, out_dir,
                                             cluster_rows = TRUE, cluster_cols = TRUE,
                                             dist_method = "euclidean", clust_method = "ward.D2",
                                             # *** MODIFIED: Increased base fontsize ***
                                             fontsize = 10, angle_col = 45) {
  
  if (is.null(heatmap_matrix_data) || nrow(heatmap_matrix_data) < 2 || ncol(heatmap_matrix_data) < 2) {
    cat("Skipping pheatmap plot for", direction_label, "- Insufficient rows/columns (need >=2 of each).\n\n")
    return(NULL)
  }
  cat("Preparing cluster pheatmap with dendrograms for:", direction_label, "\n")
  
  heatmap_matrix <- heatmap_matrix_data %>%
    mutate(Cluster = as.character(Cluster)) %>%
    tibble::column_to_rownames("Cluster") %>%
    as.matrix()
  
  can_cluster_rows <- cluster_rows && nrow(heatmap_matrix) >= 2 && !any(apply(heatmap_matrix, 1, function(r) length(unique(na.omit(r))) <= 1))
  can_cluster_cols <- cluster_cols && ncol(heatmap_matrix) >= 2 && !any(apply(heatmap_matrix, 2, function(c) length(unique(na.omit(c))) <= 1))
  if (!can_cluster_rows && cluster_rows) cat("  Disabling row clustering due to insufficient rows or constant values.\n")
  if (!can_cluster_cols && cluster_cols) cat("  Disabling column clustering due to insufficient columns or constant values.\n")
  
  
  plot_title <- paste("Clustered Heatmap of GO Enrichment (-log10 P) by Cluster -", direction_label, "Genes")
  # *** MODIFIED: Filename base ***
  filename_base <- file.path(out_dir, paste0("pHeatmap_Clustered_Cluster_", sp_tag, "_", cond_tag, "_", direction_label, "_k", clust_k))
  filename_png <- paste0(filename_base, ".png")
  filename_pdf <- paste0(filename_base, ".pdf") # For optional PDF output
  
  
  # *** MODIFIED: Dynamic font sizes with larger minimums and adjusted plot dimensions ***
  # Adjust base font size calculation if needed, ensure minimums are larger
  row_fontsize <- max(6, fontsize - nrow(heatmap_matrix) * 0.03) # Larger minimum, adjust scaling
  col_fontsize <- max(6, fontsize - ncol(heatmap_matrix) * 0.02) # Larger minimum, adjust scaling
  # Increase base width/height further if labels are long or numerous
  plot_width <- max(8, ncol(heatmap_matrix) * 0.2 + 5)  # Increased base and scaling factor
  plot_height <- max(8, nrow(heatmap_matrix) * 0.2 + 5) # Increased base and scaling factor
  
  cat("  Generating pheatmap plot...\n")
  cat("  Dimensions (w x h):", plot_width, "x", plot_height, "inches\n")
  cat("  Row Fontsize:", row_fontsize, " Col Fontsize:", col_fontsize, "\n")
  
  plot_object <- NULL
  tryCatch({
    # Use a file device explicitly for better control, especially for PDF
    # PNG Saving
    png(filename_png, width = plot_width, height = plot_height, units = "in", res = 600) # Increased resolution
    pheatmap(
      heatmap_matrix,
      color = colorRampPalette(color_vec)(50),
      border_color = "grey90",
      cluster_rows = can_cluster_rows,
      cluster_cols = can_cluster_cols,
      clustering_distance_rows = dist_method,
      clustering_distance_cols = dist_method,
      clustering_method = clust_method,
      main = plot_title,
      fontsize = fontsize, # Base fontsize for title etc.
      fontsize_row = row_fontsize,
      fontsize_col = col_fontsize,
      angle_col = angle_col,
      silent = FALSE # Plot to the device
    )
    dev.off() # Close PNG device
    cat("Saved cluster pheatmap plot (High Res PNG) to:", filename_png, "\n")
    
    # *** NEW: Optional PDF output ***
    # Use cairo_pdf for better embedding and handling of fonts/symbols
    pdf(filename_pdf, width = plot_width, height = plot_height) # Use default PDF device or cairo_pdf if installed/preferred
    # cairo_pdf(filename_pdf, width = plot_width, height = plot_height) # Alternative
    pheatmap(
      heatmap_matrix,
      color = colorRampPalette(color_vec)(50),
      border_color = "grey90",
      cluster_rows = can_cluster_rows,
      cluster_cols = can_cluster_cols,
      clustering_distance_rows = dist_method,
      clustering_distance_cols = dist_method,
      clustering_method = clust_method,
      main = plot_title,
      fontsize = fontsize,
      fontsize_row = row_fontsize,
      fontsize_col = col_fontsize,
      angle_col = angle_col,
      silent = FALSE # Plot to the device
    )
    dev.off() # Close PDF device
    cat("Saved cluster pheatmap plot (PDF) to:", filename_pdf, "\n\n")
    
  }, error = function(e) {
    cat("ERROR generating or saving cluster pheatmap:", e$message, "\n")
    # Make sure device is off if error occurred during plotting
    if(dev.cur() != 1) dev.off()
  })
  
  # Return NULL as pheatmap plots directly to device or file here
  return(invisible(NULL))
}
cat("Cluster Heatmap preparation and plotting functions redefined with larger fonts and higher resolution.\n") # Update message

# Prepare heatmap data using the significant GO results per cluster
# *** MODIFIED: Use cluster heatmap prep function ***
heatmap_data_UP <- prepare_heatmap_data_cluster(significant_go_UP, p_value_threshold)
heatmap_data_DO <- prepare_heatmap_data_cluster(significant_go_DO, p_value_threshold)

# Plot Heatmaps
# Define color gradients - use more steps for smoother gradient
up_heatmap_colors <- c("white", "orange", "darkorange2")
do_heatmap_colors <- c("white", "lightblue", "cyan4")

# *** MODIFIED: Use cluster plot function and pass cluster k ***
p_heatmap_UP <- plot_go_heatmap_pheatmap_cluster(heatmap_data_UP, "UP", up_heatmap_colors, species_tag, condition_tag, cluster_file_k, output_dir)
p_heatmap_DO <- plot_go_heatmap_pheatmap_cluster(heatmap_data_DO, "DO", do_heatmap_colors, species_tag, condition_tag, cluster_file_k, output_dir)

##################################################################################
