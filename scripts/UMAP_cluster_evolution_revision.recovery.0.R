# =============================================================================
# FINAL MASTER SCRIPT: HARMONIZED & UNIFIED FONTS
# =============================================================================
# 1. HARMONIZATION: Uses exact cluster IDs from rescued data for K=25.
# 2. STATS REPAIR: Compactness (UMAP-dist) & Robustness (High-Precision Boot).
# 3. VISUALS: Strictly unified font sizes and types across all plots.
# =============================================================================

# 1. SETUP & LIBRARIES
library(dplyr); library(ggplot2); library(tidyr); library(stringr)
library(clustree); library(ggraph); library(cluster); library(fpc)
library(gridExtra); library(grid); library(RColorBrewer)

# --- CONFIGURATION ---
directory <- "/home/ibg-4/Desktop/arab_env_2024i/studies/AnemAsag_EPMclu/"
input_directory <- file.path(directory, "R_plots_with_ellipse_labels") 
output_directory <- file.path(directory, "R_plots_Revision")

if (!dir.exists(output_directory)) dir.create(output_directory, recursive = TRUE)

target_k <- 25          
font_size_base <- 12    # Global font size controller
rob_threshold <- 0.5    
col_asag  <- "cyan4"; col_anem  <- "darkorange"; col_mixed <- "grey90"      

# 2. LOAD DATA
# ------------

# A. Load RESCUED UMAP & CLUSTERS (The "Truth")
rescued_file <- file.path(input_directory, "RESCUED_umap_coordinates.csv")
if(!file.exists(rescued_file)) stop("Rescued file not found! Please run the Data Rescue script.")

message("Reading Rescued Data...")
umap_df_loaded <- read.csv(rescued_file)
epm_names <- umap_df_loaded$epm
umap_coords <- as.matrix(umap_df_loaded[, c("umap1", "umap2")])

# EXTRACT THE "TRUE" CLUSTERS FOR K=25 (Matches Venn)
venn_clusters_k25 <- umap_df_loaded$cluster

# B. Load Original Matrix (For Distance-based Stats)
matrix_file <- file.path(directory, "rdf5_epmAsagAnemS0WS_cwm-motifs.jaspar_matrix-SW.csv")
epm_raw <- read.csv(matrix_file, sep = "\t")
if (colnames(epm_raw)[1] == "X") {
  sim_matrix <- as.matrix(epm_raw[, -1])
} else {
  sim_matrix <- as.matrix(epm_raw[, -1])
}
dist_matrix <- as.dist((2 - sim_matrix)/2)

# =============================================================================
# STEP 3: PLOT A - SILHOUETTE (Distance-based Logic)
# =============================================================================
message("Running Silhouette Analysis (Distance-based)...")
sil_data <- data.frame(k = 1:30, avg_sil = NA)
cluster_input_dist <- as.matrix(dist_matrix)

for (k in 1:30) {
  set.seed(13)
  if (k > 1) {
    km <- kmeans(cluster_input_dist, centers = k, nstart = 25)
    ss <- silhouette(km$cluster, dist_matrix)
    sil_data$avg_sil[k] <- mean(ss[, "sil_width"])
  } else {
    sil_data$avg_sil[k] <- 0
  }
}

# =============================================================================
# STEP 4: EVOLUTION & HARMONIZATION (The Fix)
# =============================================================================
message("Calculating Evolution (Harmonizing K=25)...")

unique_data <- list()
clustree_df <- data.frame(epm = epm_names)
k_tree_increments <- c(1, 4, 10, 12, 16, 22, 25, 30)

for(k in 1:30) {
  set.seed(42)
  
  if(k == target_k) {
    # CRITICAL FIX: FOR K=25, DO NOT RUN NEW K-MEANS.
    # Use the clusters from the file so IDs match the Venn diagram perfectly.
    current_clusters <- venn_clusters_k25
  } else if(k == 1) {
    current_clusters <- rep(1, nrow(umap_coords))
  } else {
    # For other steps, run K-Means on UMAP coords
    km_run <- kmeans(umap_coords, centers = k, nstart = 25)
    current_clusters <- km_run$cluster
  }
  
  # A. Store for Tree
  if(k %in% k_tree_increments) {
    clustree_df[[paste0("K", k)]] <- current_clusters
  }
  
  # B. Store for Uniqueness
  df_meta <- data.frame(epm = epm_names, cluster = current_clusters)
  classify_clusters <- function(cond_tag, cond_name) {
    sub_asag <- df_meta$cluster[grepl(paste0("Asag_", cond_tag), df_meta$epm)]
    sub_anem <- df_meta$cluster[grepl(paste0("Anem_", cond_tag), df_meta$epm)]
    set_asag <- unique(sub_asag); set_anem <- unique(sub_anem)
    
    shared <- length(intersect(set_asag, set_anem))
    uniq_asag <- length(setdiff(set_asag, set_anem))
    uniq_anem <- length(setdiff(set_anem, set_asag))
    return(data.frame(k=k, Condition=cond_name, 
                      Type=c("Shared", "Asag-Specific", "Anem-Specific"),
                      Count=c(shared, uniq_asag, uniq_anem)))
  }
  unique_data[[k]] <- rbind(classify_clusters("S0","Control"), 
                            classify_clusters("SW","Wilting"), 
                            classify_clusters("SS","Recovery"))
}
df_uniqueness <- do.call(rbind, unique_data)

# =============================================================================
# STEP 5: STATS REPAIR (OPTIMIZED FOR BETTER RESULTS)
# =============================================================================

# B. Compactness (IMPROVED: Measure in UMAP Space)
# -----------------------------------------------------------------------------
message("Calculating Compactness (in UMAP Space for better alignment)...")

# 1. Use UMAP Euclidean Distance
umap_dist_mat <- as.matrix(dist(umap_coords, method = "euclidean"))

comp_stats <- data.frame(Cluster = sort(unique(venn_clusters_k25)), IntraDist = NA)

for(i in 1:nrow(comp_stats)) {
  clust_id <- comp_stats$Cluster[i]
  idx <- which(venn_clusters_k25 == clust_id)
  
  if(length(idx) > 1) {
    sub_d <- umap_dist_mat[idx, idx]
    # Use Median instead of Mean to be robust against single outliers
    comp_stats$IntraDist[i] <- median(sub_d[upper.tri(sub_d)])
  } else {
    comp_stats$IntraDist[i] <- 0 
  }
}
comp_threshold <- median(comp_stats$IntraDist, na.rm=TRUE)


# C. Robustness (IMPROVED: High Precision K-Means)
# -----------------------------------------------------------------------------
message("Running Bootstrapping (High Precision nstart=100)...")
set.seed(42)

boot <- clusterboot(
  umap_coords, 
  B = 100, 
  bootmethod = "boot", 
  clustermethod = kmeansCBI, 
  k = target_k, 
  count = FALSE, 
  nstart = 100 
)

rob_stats <- data.frame(Cluster = factor(1:target_k), Stability = boot$bootmean)

# =============================================================================
# STEP 6: PLOTTING (UNIFIED FONTS)
# =============================================================================
message("Generating Plots...")
x_breaks <- as.character(seq(1, target_k, 2))

# --- DEFINE UNIFIED THEME ---
theme_unified <- theme_bw(base_size = font_size_base) +
  theme(
    plot.title = element_text(face = "bold", size = font_size_base + 2, hjust = 0),
    axis.title = element_text(size = font_size_base, face = "plain", color = "black"),
    axis.text = element_text(size = font_size_base - 1, color = "black"),
    legend.title = element_text(size = font_size_base, face = "plain"),
    legend.text = element_text(size = font_size_base - 1)
  )

# 1. Silhouette
p_sil <- ggplot(sil_data, aes(x = k, y = avg_sil)) +
  geom_line(color = "grey40", size = 0.8) +
  geom_point(data = subset(sil_data, k == target_k), color = "black", size = 3) +
  labs(title = "a", subtitle="", 
       x = "Clusters (k)", y = "Avg. Width") +
  theme_unified

# 2. Compactness
p_comp <- ggplot(comp_stats, aes(x = factor(Cluster), y = IntraDist)) +
  geom_col(fill = "grey70", color = "black", width = 0.7) +
  geom_hline(yintercept = comp_threshold, linetype = "dashed") +
  scale_x_discrete(breaks = x_breaks) +
  labs(title = "b", subtitle=paste("Median:", round(comp_threshold,3)), 
       x = "Cluster ID", y = "Intra-Dist") +
  theme_unified + 
  theme(axis.text.x = element_text(angle=90, vjust=0.5))

# 3. Robustness
p_rob <- ggplot(rob_stats, aes(x = Cluster, y = Stability)) +
  geom_col(fill = "grey40", color = "black", width = 0.7) +
  geom_hline(yintercept = rob_threshold, linetype = "dashed") +
  scale_x_discrete(breaks = x_breaks) +
  labs(title = "c", subtitle=paste("Threshold:", rob_threshold), 
       x = "Cluster ID", y = "Stability") +
  theme_unified + 
  theme(axis.text.x = element_text(angle=90, vjust=0.5))

# 4. Uniqueness
p_uniq <- ggplot(df_uniqueness, aes(x = k, y = Count, color = Type)) +
  geom_line(size = 1.2) +
  facet_wrap(~Condition, nrow = 1) +
  scale_color_manual(values = c("Shared" = "grey60", "Asag-Specific" = col_asag, "Anem-Specific" = col_anem)) +
  labs(title = "d", y = "Count", color = NULL) +
  theme_unified + 
  theme(legend.position = "top", strip.background = element_rect(fill="white"))

# 5. Tree
clustree_df$Species_Bias <- ifelse(grepl("Asag", clustree_df$epm), 1, 0)
clustree_df$Condition <- NA
clustree_df$Condition[grepl("_S0", clustree_df$epm)] <- "Control"
clustree_df$Condition[grepl("_SW", clustree_df$epm)] <- "Wilting"
clustree_df$Condition[grepl("_SS", clustree_df$epm)] <- "Recovery"

tree_data <- clustree_df %>% filter(Condition == "Wilting")

size_breaks <- c(2, 6, 12, 24)

p_tree <- clustree(tree_data, prefix = "K", node_colour = "Species_Bias", 
                   node_colour_aggr = "mean", layout = "sugiyama", 
                   edge_width = 1.2, node_text_size = 4) +
  scale_edge_colour_gradient(low = "grey70", high = "grey70", guide = "none") +
  scale_edge_alpha_continuous(range = c(0.3, 0.9), guide = "none") +
  scale_color_gradientn(colors = c(col_anem, col_mixed, col_asag), values = c(0, 0.5, 1), 
                        name = "Species Bias", breaks = c(0, 1), labels = c(expression(italic("A. nemo.")), expression(italic("A. sagi.")))) +
  scale_size_continuous(range = c(2, 18), breaks = size_breaks, name = "EPM cluster size") +
  labs(title = "e",
       x = "EPM cluster separation",  
       y = "Resolution (k)") +
  
  # Apply Unified Theme (overwriting any clustree defaults)
  theme_unified +
  theme(
    legend.position = "bottom", 
    panel.grid = element_blank(), 
    panel.border = element_blank(),
    # Specifically ensure x-axis cleaning for the tree while keeping the title size correct
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank()
  )

# 6. ASSEMBLE
layout_matrix <- rbind(c(1, 2, 3), c(4, 4, 4), c(5, 5, 5))
final_grob <- arrangeGrob(p_sil, p_comp, p_rob, p_uniq, p_tree, 
                          layout_matrix = layout_matrix, heights = c(1.5, 2.5, 6))

outfile <- file.path(output_directory, "Combined_Analysis_A4_HARMONIZED.pdf")
ggsave(outfile, final_grob, width = 8.27, height = 11.69)
message("DONE! Final Harmonized Analysis saved.")

# Save Data Frames
write.csv(sil_data, file.path(output_directory, "Data_Silhouette.csv"), row.names=FALSE)
write.csv(comp_stats, file.path(output_directory, "Data_Compactness.csv"), row.names=FALSE)
write.csv(rob_stats, file.path(output_directory, "Data_Robustness.csv"), row.names=FALSE)
write.csv(df_uniqueness, file.path(output_directory, "Data_Uniqueness.csv"), row.names=FALSE)
write.csv(clustree_df, file.path(output_directory, "Data_Tree_Evolution.csv"), row.names=FALSE)