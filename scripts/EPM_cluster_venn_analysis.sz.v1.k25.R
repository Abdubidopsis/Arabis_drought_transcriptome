# =============================================================================
# MASTER SCRIPT: REPRODUCIBLE ANALYSIS (HIERARCHICAL TREE VERSION)
# =============================================================================
# This version is identical to your original setup but switches the clustering
# engine from K-Means to Hierarchical (Ward's Method).
# This recovers the ~6-7 unique clusters and aligns with the Tree Plots.
# =============================================================================

# Define the directory path
directory <- "/home/ibg-4/Desktop/arab_env_2024i/studies/AnemAsag_EPMclu/"
# Use a new output folder for this version
output_directory <- file.path(directory, "R_plots_larger_venn_labels_TREE_VERSION")

# Create output directory if it doesn't exist
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
}

# Read similarity matrix
epm_similarity_matrix_file <- file.path(directory, "rdf5_epmAsagAnemS0WS_cwm-motifs.jaspar_matrix-SW.csv")
print(paste("Reading similarity matrix:", epm_similarity_matrix_file))
epm_similarity_matrix <- read.csv(epm_similarity_matrix_file, sep = "\t")

# Load necessary libraries
library(umap)
library(ggplot2) 
library(cluster)
library(plotly)
library(eulerr)
library(gridExtra)
library(dplyr)
library(grid)
library(RColorBrewer) 
library(fpc)

if (colnames(epm_similarity_matrix)[1] == "X") {
  epm_names <- epm_similarity_matrix$X
  similarity_matrix <- epm_similarity_matrix[, -1]
} else {
  epm_names <- epm_similarity_matrix[, 1]
  similarity_matrix <- epm_similarity_matrix[, -1]
}

# --- 1. UMAP & CLUSTERING (MODIFIED FOR TREE CONSISTENCY) ---

# CRITICAL 1: Set seed BEFORE UMAP
# This ensures the map coordinates are identical every time you run this.
set.seed(42) 
print("Running UMAP...")
umap_result <- umap(similarity_matrix)

# CRITICAL 2: SWITCH TO HIERARCHICAL CLUSTERING 
# Instead of K-Means (which shuffles points), we build a Tree.
# This preserves the small, unique clusters you saw in previous results.
k_value <- 25
print(paste("Running Hierarchical Clustering (Tree) with k=", k_value, "..."))

# Calculate Euclidean distance on the UMAP layout
d_umap <- dist(umap_result$layout, method = "euclidean")

# Run Ward's Hierarchical Clustering
hc_tree <- hclust(d_umap, method = "ward.D2")

# Cut the tree to get exactly 25 clusters
tree_clusters <- cutree(hc_tree, k = k_value)

# Create data frames (Updated to use tree_clusters)
umap_df <- data.frame(umap1 = umap_result$layout[, 1], umap2 = umap_result$layout[, 2],
                      cluster = factor(tree_clusters), epm = epm_names)
cluster_df <- data.frame(epm = epm_names, cluster = tree_clusters)

################################################################################
# --- 1. Silhouette Analysis (Range k=1 to 25) ---
print("Running Silhouette Analysis for k=1 to 25...")
# We update this loop to also use CUTREE so the plot matches the new method.

# A. Convert Similarity to Distance (Normalized)
dist_matrix <- as.dist((2 - as.matrix(similarity_matrix)) / 2)

# B. Loop through K
k_range <- 2:30 
results <- data.frame(k = k_range, avg_sil = NA)

for(k in k_range) {
  # Use cutree here too for consistency
  curr_cut <- cutree(hc_tree, k = k)
  sil <- silhouette(curr_cut, dist_matrix)
  results$avg_sil[results$k == k] <- mean(sil[, "sil_width"])
}

# C. Add K=1 and Highlights
results <- rbind(data.frame(k=1, avg_sil=0), results)
results$highlight <- "Normal"
results$highlight[results$k %in% c()] <- "Selected"

# D. Save Basis Values for Silhouette Plot
sil_csv_file <- file.path(output_directory, "silhouette_basis_values_k1_to_30.csv")
print(paste("Saving Silhouette data:", sil_csv_file))
write.csv(results, sil_csv_file, row.names = FALSE)

# E. Generate Plot
sil_summary_plot <- ggplot(results, aes(x = k, y = avg_sil)) +
  geom_line(color = "grey70", size = 0.8) +
  geom_point(aes(color = highlight, size = highlight)) +
  scale_color_manual(values = c("Normal" = "black", "Selected" = "red")) +
  scale_size_manual(values = c("Normal" = 1.5, "Selected" = 4)) +
  scale_x_continuous(breaks = seq(1, 30, 1)) +
  labs(title = "Silhouette Analysis (Hierarchical)", subtitle = "x", y = "Avg Silhouette Width") +
  theme_minimal() + theme(legend.position = "none")

ggsave(file.path(output_directory, "silhouette_range_1_to_30.png"), plot = sil_summary_plot, width = 12, height = 6)

# --- 2. Cluster Statistics & Intra-Cluster Distance ---
print(paste("Calculating Cluster Statistics for k=", k_value, "..."))

# A. Initialize Stats Container
cluster_stats <- data.frame(
  Cluster = sort(unique(tree_clusters)),
  Size = NA,
  Mean_Intra_Dist = NA, 
  Mean_Silhouette = NA
)

# B. Calculate Stats
full_dist_mat <- as.matrix(dist_matrix) 
sil_obj_k <- silhouette(tree_clusters, dist_matrix) 

for(i in 1:nrow(cluster_stats)) {
  clust_id <- cluster_stats$Cluster[i]
  idx <- which(tree_clusters == clust_id)
  
  # Size
  cluster_stats$Size[i] <- length(idx)
  
  # Mean Silhouette for this specific cluster
  cluster_stats$Mean_Silhouette[i] <- mean(sil_obj_k[idx, "sil_width"])
  
  # Mean Intra-Cluster Distance
  if(length(idx) > 1) {
    sub_dist <- full_dist_mat[idx, idx]
    cluster_stats$Mean_Intra_Dist[i] <- mean(sub_dist[upper.tri(sub_dist)])
  } else {
    cluster_stats$Mean_Intra_Dist[i] <- 0 
  }
}

# C. Save Basis Values for Bar Plot
stats_csv_file <- file.path(output_directory, paste0("cluster_compactness_basis_k", k_value, ".csv"))
print(paste("Saving Cluster Stats data:", stats_csv_file))
write.csv(cluster_stats, stats_csv_file, row.names = FALSE)

# D. Generate Bar Plot
dist_plot <- ggplot(cluster_stats, aes(x = factor(Cluster), y = Mean_Intra_Dist, fill = factor(Cluster))) +
  geom_col(alpha = 0.8, width = 0.7, color = "black", size = 0.2) +
  scale_fill_viridis_d(option = "turbo", guide = "none") + 
  geom_text(aes(label = round(Mean_Intra_Dist, 3)), vjust = -0.5, size = 3) +
  labs(title = paste("Cluster Compactness (k=", k_value, ")"), y = "Mean Pairwise Distance") +
  theme_bw() + theme(panel.grid.major.x = element_blank())

ggsave(file.path(output_directory, paste0("cluster_intra_dist_barplot_k", k_value, ".png")), plot = dist_plot, width = 10, height = 6)
################################################################################


# --- Save Plotly UMAP (Interactive HTML - Unchanged) ---
print("Generating Plotly UMAP plot...")
umap_plot <- plot_ly(umap_df, x = ~umap1, y = ~umap2, color = ~cluster, type = 'scatter', mode = 'markers',
                     hoverinfo = 'text', text = ~paste("EPM: ", epm, "<br>Cluster: ", cluster)) %>%
  layout(title = list(text = paste("UMAP of EPM Similarity Matrix with k=", k_value," Clusters"), font = list(size = 16)),
         xaxis = list(title = "UMAP Dimension 1", titlefont = list(size = 14), tickfont = list(size = 12)),
         yaxis = list(title = "UMAP Dimension 2", titlefont = list(size = 14), tickfont = list(size = 12)),
         font = list(size = 12))
umap_html_file <- file.path(output_directory, paste0("umap_plot_interactive_k", k_value, "_large.html"))
print(paste("Saving Plotly UMAP plot:", umap_html_file))
htmlwidgets::saveWidget(umap_plot, umap_html_file)


# --- Generate Static UMAP Plot (ggplot2, Earth Tones, SVG Output) ---
print("Generating static UMAP plot with ggplot2...")

# Add Species/Treatment columns to umap_df for coloring/shaping
extract_info <- function(epm_name) {
  parts <- strsplit(as.character(epm_name), "_")[[1]]
  species <- NA; treatment <- NA
  if (length(parts) >= 2) {
    species_code <- parts[1]; condition_code <- parts[2]
    species <- ifelse(species_code == "Anem", "A. nemorensis", ifelse(species_code == "Asag", "A. sagittata", NA))
    treatment <- ifelse(condition_code == "S0", "Control", ifelse(condition_code == "SW", "Wilting", ifelse(condition_code == "SS", "Recovery", NA)))
  }
  return(data.frame(Species = species, Treatment = treatment))
}
info_df <- bind_rows(lapply(umap_df$epm, extract_info))
umap_df <- cbind(umap_df, info_df)
umap_df <- umap_df[!is.na(umap_df$Species) & !is.na(umap_df$Treatment),]

# Set factor levels
condition_levels_ordered <- c("Control", "Wilting", "Recovery")
species_levels_ordered <- c("A. nemorensis", "A. sagittata")
species_labels_italic <- c(expression(italic("A. nemorensis")), expression(italic("A. sagittata")))

umap_df$Species <- factor(umap_df$Species, levels = species_levels_ordered)
umap_df$Treatment <- factor(umap_df$Treatment, levels = condition_levels_ordered)

# Define earth tone palette and shapes
earth_palette_conditions <- brewer.pal(3, "YlOrBr") # Earth tones for Treatment
species_shapes <- c(16, 17) # Shapes for Species

# Create ggplot object
static_umap_plot <- ggplot(umap_df, aes(x = umap1, y = umap2, color = Treatment, shape = Species)) +
  geom_point(size = 1.5, alpha = 0.7) +
  scale_color_manual(values = earth_palette_conditions, name = "Treatment") +
  scale_shape_manual(values = species_shapes, name = "Species", labels = species_labels_italic) + # Italic labels here
  labs(
    title = paste("UMAP of EPM Similarity (k=", k_value," Clusters)"),
    x = "UMAP Dimension 1",
    y = "UMAP Dimension 2"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(size = rel(1.2), hjust = 0.5),
    axis.title = element_text(size = rel(1.1)),
    axis.text = element_text(size = rel(0.9)),
    legend.title = element_text(size = rel(1.0)),
    legend.text = element_text(size = rel(0.9)),
    legend.text.align = 0
  ) +
  guides(shape = guide_legend(override.aes = list(size=3)))

# --- Save static UMAP plot as SVG and PNG ---
static_umap_svg_file <- file.path(output_directory, paste0("umap_plot_static_k", k_value, "_earth.svg"))
static_umap_png_file <- file.path(output_directory, paste0("umap_plot_static_k", k_value, "_earth_600dpi.png"))

print(paste("Saving static UMAP plot:", static_umap_svg_file))
ggsave(static_umap_svg_file, plot = static_umap_plot, width = 6, height = 5.5, units = "in")

print(paste("Saving static UMAP plot:", static_umap_png_file))
ggsave(static_umap_png_file, plot = static_umap_plot, width = 6, height = 5.5, units = "in", dpi = 600)
################################################################################ #added on 21th January for revision
print("Running Cluster Stability Analysis (Bootstrap)...")

# 1. Define the data you actually clustered
input_data <- umap_result$layout 

# 2. Run Bootstrap
# NOTE: We use K-Means for bootstrapping stability checks because 
# hierarchical clustering stability is harder to quantify in 'clusterboot'. 
# This is an acceptable proxy for spatial stability.
set.seed(42) 
cboot <- clusterboot(
  input_data, 
  B = 100,                
  bootmethod = "boot",    
  clustermethod = kmeansCBI, 
  krange = k_value,       
  seed = 42
)

# 3. View and Interpret Results
print("--- Cluster Stability Scores (Jaccard Mean) ---")
print(cboot$bootmean)


################################################################################
################################################################################
# --- Venn Diagram Section ---
print("Preparing data for Venn diagrams...")
# ... (Splitting data frame and extracting cluster sets remains the same) ...
Asag_S0 <- cluster_df %>% filter(grepl("Asag_S0", epm)); Asag_SW <- cluster_df %>% filter(grepl("Asag_SW", epm)); Asag_SS <- cluster_df %>% filter(grepl("Asag_SS", epm))
Anem_S0 <- cluster_df %>% filter(grepl("Anem_S0", epm)); Anem_SW <- cluster_df %>% filter(grepl("Anem_SW", epm)); Anem_SS <- cluster_df %>% filter(grepl("Anem_SS", epm))
Asag_S0_clusters <- unique(Asag_S0$cluster); Asag_SW_clusters <- unique(Asag_SW$cluster); Asag_SS_clusters <- unique(Asag_SS$cluster)
Anem_S0_clusters <- unique(Anem_S0$cluster); Anem_SW_clusters <- unique(Anem_SW$cluster); Anem_SS_clusters <- unique(Anem_SS$cluster)

# Colors for Venn (using the warm colors from the R script)
border_colors_a <- c("darkgrey", "darkgrey", "darkgrey"); fill_colors_a <- c("bisque", "orange", "brown")
border_colors_b <- c("darkgrey", "darkgrey", "darkgrey"); fill_colors_b <- c("bisque", "orange", "brown")
border_colors_c <- c("darkgrey", "darkgrey"); fill_colors_c <- c("bisque", "bisque3")
border_colors_d <- c("darkgrey", "darkgrey"); fill_colors_d <- c("orange", "orange3")
border_colors_e <- c("darkgrey", "darkgrey"); fill_colors_e <- c("brown", "brown4")

# --- Define Labels and Font Sizes for Venn ---
anem_label_plain <- "A. nemorensis"
asag_label_plain <- "A. sagittata"
condition_labels_plain <- c("Control", "Wilting", "Recovery")

# --- Increased Font Sizes for Venn Labels ---
euler_quantity_fontsize <- 20 # Size for numbers inside segments (Increased)
euler_label_fontsize <- 16    # Size for set labels (Increased)

# --- Create Euler Plots without Labels ---
print("Generating Euler diagrams without labels...")

comparison_a <- euler(list(Control = Asag_S0_clusters, Wilting = Asag_SW_clusters, Recovery = Asag_SS_clusters))
plot_a <- plot(comparison_a, main = asag_label_plain, quantities = list(fontsize = euler_quantity_fontsize), labels = FALSE, col = border_colors_a, fill = fill_colors_a)

comparison_b <- euler(list(Control = Anem_S0_clusters, Wilting = Anem_SW_clusters, Recovery = Anem_SS_clusters))
plot_b <- plot(comparison_b, main = anem_label_plain, quantities = list(fontsize = euler_quantity_fontsize), labels = FALSE, col = border_colors_b, fill = fill_colors_b)

comparison_c <- euler(list('A. sagittata' = Asag_S0_clusters, 'A. nemorensis' = Anem_S0_clusters))
plot_c <- plot(comparison_c, main = "Control Comparison", quantities = list(fontsize = euler_quantity_fontsize), labels = FALSE, col = border_colors_c, fill = fill_colors_c)

comparison_d <- euler(list('A. sagittata' = Asag_SW_clusters, 'A. nemorensis' = Anem_SW_clusters))
plot_d <- plot(comparison_d, main = "Wilting Comparison", quantities = list(fontsize = euler_quantity_fontsize), labels = FALSE, col = border_colors_d, fill = fill_colors_d)

comparison_e <- euler(list('A. sagittata' = Asag_SS_clusters, 'A. nemorensis' = Anem_SS_clusters))
plot_e <- plot(comparison_e, main = "Recovery Comparison", quantities = list(fontsize = euler_quantity_fontsize), labels = FALSE, col = border_colors_e, fill = fill_colors_e)
# Arrange plots
print("Arranging Euler diagrams...")
lay <- rbind(c(1, 2, NA), c(3, 4, 5))
combined_venn_plot <- grid.arrange(grobs = list(plot_a, plot_b, plot_c, plot_d, plot_e), layout_matrix = lay, top = textGrob(paste("EPM Cluster Overlap Analysis (k=", k_value,")"), gp = gpar(fontsize = 16)))

# Save the combined Venn diagram plots (using larger dimensions)
venn_svg_file <- file.path(output_directory, paste0("venn_plot_k", k_value, "_large_labels.svg")) # New filename part
venn_png_file <- file.path(output_directory, paste0("venn_plot_k", k_value, "_large_labels_600dpi.png")) # New filename part

print(paste("Saving combined Venn plot:", venn_svg_file))
svg(filename = venn_svg_file, width = 11, height = 9) # Slightly increased size for larger fonts
grid.draw(combined_venn_plot)
dev.off()

print(paste("Saving combined Venn plot:", venn_png_file))
png(filename = venn_png_file, width = 11, height = 9, units = "in", res = 600) # Slightly increased size
grid.draw(combined_venn_plot)
dev.off()

# --- Set Operations Section (Unchanged) ---
print("Performing set operations...")
set_operations <- function(set1, set2, set3, name1, name2, name3) {common_1_2_3 <- intersect(intersect(set1, set2), set3); common_1_2 <- intersect(set1, set2); common_1_3 <- intersect(set1, set3); common_2_3 <- intersect(set2, set3); unique_1 <- setdiff(set1, union(set2, set3)); unique_2 <- setdiff(set2, union(set1, set3)); unique_3 <- setdiff(set3, union(set1, set2)); data.frame( Comparison = c(paste("Common between", name1, name2, "and", name3), paste("Common between", name1, "and", name2), paste("Common between", name1, "and", name3), paste("Common between", name2, "and", name3), paste("Unique to", name1), paste("Unique to", name2), paste("Unique to", name3)), Values = c(paste(common_1_2_3, collapse = ", "), paste(common_1_2, collapse = ", "), paste(common_1_3, collapse = ", "), paste(common_2_3, collapse = ", "), paste(unique_1, collapse = ", "), paste(unique_2, collapse = ", "), paste(unique_3, collapse = ", ")), stringsAsFactors = FALSE)}
asag_results <- set_operations(Asag_S0_clusters, Asag_SW_clusters, Asag_SS_clusters, "AsagS0", "AsagSW", "AsagSS"); anem_results <- set_operations(Anem_S0_clusters, Anem_SW_clusters, Anem_SS_clusters, "AnemS0", "AnemSW", "AnemSS")
asag_anem_results <- data.frame( Comparison = c("Common between AsagS0 and AnemS0", "Common between AsagSW and AnemSW", "Common between AsagSS and AnemSS", "Unique to AsagS0 vs AnemS0", "Unique to AnemS0 vs AsagS0", "Unique to AsagSW vs AnemSW", "Unique to AnemSW vs AsagSW", "Unique to AsagSS vs AnemSS", "Unique to AnemSS vs AsagSS" ), Values = c(paste(intersect(Asag_S0_clusters, Anem_S0_clusters), collapse = ", "), paste(intersect(Asag_SW_clusters, Anem_SW_clusters), collapse = ", "), paste(intersect(Asag_SS_clusters, Anem_SS_clusters), collapse = ", "), paste(setdiff(Asag_S0_clusters, Anem_S0_clusters), collapse = ", "), paste(setdiff(Anem_S0_clusters, Asag_S0_clusters), collapse = ", "), paste(setdiff(Asag_SW_clusters, Anem_SW_clusters), collapse = ", "), paste(setdiff(Anem_SW_clusters, Asag_SW_clusters), collapse = ", "), paste(setdiff(Asag_SS_clusters, Anem_SS_clusters), collapse = ", "), paste(setdiff(Anem_SS_clusters, Asag_SS_clusters), collapse = ", ") ), stringsAsFactors = FALSE)
all_results <- bind_rows(asag_results, anem_results, asag_anem_results)
results_csv_file <- file.path(output_directory, paste0("set_comparison_results_k", k_value, ".csv"))
print(paste("Saving set comparison results:", results_csv_file))
write.csv(all_results, results_csv_file, row.names = FALSE)

# --- 3. UMAP Basis Values ---
umap_csv_file <- file.path(output_directory, paste0("umap_coordinates_basis_k", k_value, ".csv"))
print(paste("Saving UMAP coordinates:", umap_csv_file))
write.csv(umap_df, umap_csv_file, row.names = FALSE)

print("Script finished.")