###--Differential expression analysis using DESeq2--#######################
###' Abdul Saboor Khan (abdul.suboor123@yahoo.com) -- first created on 16.Jan.2023 - 
###' Modified on 11th July 2025
###' Last modified on 20th February 2026.
###' This script runs for Arabis sagittata and nemorensis counts files to find differential expression
###' first, get the counts file from the mapping (sam) files by using htseq-counts in bash
#### here the name "survival" used has been renamed as "recovery" in the paper

#rm(list=ls()) #clean the global environment
#dev.off() # clean the plot window

#install.packages("BiocManager")
#BiocManager::install("DESeq2")

#####-- These libraries are required--#######
library(BiocManager)
library(KEGGREST)
library("org.At.tair.db")
library(Rgraphviz)
library(topGO)
library(biomaRt)
library(ggplot2)
library(AnnotationDbi)
library(clusterProfiler)
require(ggplot2)
library(scales)
library(pasilla)
library(grid)
library(ashr)
library(apeglm)
library(data.table)
library(dplyr)
library(GGally)
library(ggplot2)
library(stringr)
library(corrplot)
library("DESeq2")
library("BiocParallel")
library("IHW")
library("iSEE")
library("Glimma")
library("vsn")
library("pheatmap")
library("RColorBrewer")
library("EnhancedVolcano")
library(edgeR)
library(readr)

#######################################################################################################
###############                         ###############################################################
############### all sampples.           ###############################################################
###############                         ###############################################################
#######################################################################################################
setwd("/Users/Shared/Files From d.localized/PhD-Uni-Koeln_2021-2024/PhD_work/seeds data/2nd Experiment 2022/transcriptome/Arabis_RNA_raw_data/90-774066047/anaylsis_drought_mRNA_both_species_with_new_two_genomes/analysis_both_species_with_nem_genomes/")
directory <- "/Users/Shared/Files From d.localized/PhD-Uni-Koeln_2021-2024/PhD_work/seeds data/2nd Experiment 2022/transcriptome/Arabis_RNA_raw_data/90-774066047/anaylsis_drought_mRNA_both_species_with_new_two_genomes/analysis_both_species_with_nem_genomes/"

###-- define the pattern of files to be analysed, the file should end as ".txt" ######
sampleFiles <- grep("txt",list.files(directory),value=TRUE)
condition <- c("control", "control", "control", "control", "survival", "survival", "survival", "wilting", "wilting", "wilting", "wilting", "control", "control", "control", "control", "survival", "survival", "survival", "wilting", "wilting", "wilting", "wilting")
genotype <- c('nemorensis', 'nemorensis', 'nemorensis', 'nemorensis', 'nemorensis', 'nemorensis', 'nemorensis', 'nemorensis', 'nemorensis','nemorensis', 'nemorensis', 'sagitatta', 'sagitatta', 'sagitatta','sagitatta','sagitatta','sagitatta','sagitatta','sagitatta','sagitatta','sagitatta','sagitatta')
sampleTable <- data.frame(sampleName = sampleFiles, fileName = sampleFiles, condition = condition, genotype = genotype)
sampleTable$condition <- factor(sampleTable$condition)
sampleTable$genotype <- factor(sampleTable$genotype)

##--deseq_from_htseqcount######
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = directory, design = ~ genotype + condition + genotype:condition) #--(design= ~ genotype + condition + genotype:condition)--##

# Filter out genes with low counts across all samples
#--- get DEG ---#

count_threshold <- 100  # Minimum average count threshold
num_samples <- 22       # Total number of samples
ddsHTSeq <- ddsHTSeq[rowSums(counts(ddsHTSeq)) / num_samples > count_threshold, ]
#count_threshold <- 100  # Minimum total counts across samples
#ddsHTSeq <- ddsHTSeq[rowSums(counts(ddsHTSeq)) > count_threshold, ]
# Run the DESeq pipeline
#ddsHTSeq <- DESeq(ddsHTSeq)
dds <- DESeq(ddsHTSeq, fitType = "mean")

res_dds_all <- results(dds)
#write.csv(res_dds_all, "script_for_juliette/new_go_output_from_deseqdrought_all/all_samples_res_dds_9315.csv")

##Filtering genes with low counts_normalization ###
dds_lowcount <- estimateSizeFactors(dds)

sizeFactors(dds_lowcount)

normalized_counts <- counts(dds_lowcount, normalized = TRUE)

#write.csv(normalized_counts, "script_for_juliette/new_go_output_from_deseqdrought_all/all_samples_normalized_9315.csv")

norm_deg <- normalizeBetweenArrays(normalized_counts, method="scale")
k <- 4
kmeans_result <- kmeans(t(norm_deg), centers=k)
d <- pheatmap(norm_deg,
              scale = "row",
              clustering_distance_rows = "euclidean",
              clustering_method = "complete",
              add.clusters = kmeans_result$cluster,
              color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
              fontsize_row = 8,
              show_rownames = FALSE,
              show_colnames = TRUE)

d
#ggsave("script_for_juliette/new_go_output_from_deseqdrought_all/kmeans_all_samples_9315_hclust.png", plot = d, width = 7, height = 7, dpi = 300)
#ggsave("script_for_juliette/new_go_output_from_deseqdrought_all/kmeans_all_samples_9315_hclust.pdf", plot = d, width = 7, height = 7, dpi = 300)

### hsclust k-means for DEGs

# Further filter res_dds_filtered to keep only those genes with padj < 0.05
res_dds_significant <- subset(res_dds_all, padj < 0.05)
#write.csv(res_dds_significant, "script_for_juliette/new_go_output_from_deseqdrought_all/all_samples_res_dds_DEGs_3526.csv")
significant_genes_df <- as.data.frame(res_dds_significant)
# filter the normalized_counts to keep only the genes that are in significant_genes_df
significant_gene_ids <- rownames(significant_genes_df)  # Get the gene IDs that are significant
normalized_counts_df_significant <- normalized_counts[rownames(normalized_counts) %in% significant_gene_ids, ]
#write.csv(normalized_counts_df_significant, "script_for_juliette/new_go_output_from_deseqdrought_all/all_samples_DEGs_normalized_3526.csv")

norm_deg <- normalizeBetweenArrays(normalized_counts_df_significant, method="scale")
k <- 4
kmeans_result <- kmeans(t(norm_deg), centers=k)
d <- pheatmap(norm_deg,
              scale = "row",
              clustering_distance_rows = "euclidean",
              clustering_method = "complete",
              add.clusters = kmeans_result$cluster,
              color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
              fontsize_row = 8,
              show_rownames = FALSE,
              show_colnames = TRUE)

d
#ggsave("script_for_juliette/new_go_output_from_deseqdrought_all/kmeans_DEGs_all_samples_3526_hclust.png", plot = d, width = 7, height = 7, dpi = 300)
#ggsave("script_for_juliette/new_go_output_from_deseqdrought_all/kmeans_DEGs_all_samples_3526_hclust.pdf", plot = d, width = 7, height = 7, dpi = 300)

#####
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
resultsNames(dds)
dds$group <- factor(paste0(dds$genotype, dds$condition))
design(dds) <- ~ group
levels(dds$group)

#####repeat for all levels to have all pairwise foldchange estimates#######
####ref "sagitattawilting" ####
dds$group<- relevel(dds$group, ref = "sagitattawilting")
levels(dds$group)
dds_test <- DESeq(dds)
resultsNames(dds_test)

# Extract contrasts for control vs. wilting within each species
nemwilt_sagwilt <- results(dds_test, name = "group_nemorensiswilting_vs_sagitattawilting")
#write.csv(nemwilt_sagwilt, "script_for_juliette/new_go_output_from_deseqdrought_all/nemwilt_sagwilt.csv")

# Extract contrasts for recovery vs. wilting within each species
sagsurv_sagwilt <- results(dds_test, name = "group_sagitattasurvival_vs_sagitattawilting")
#write.csv(sagsurv_sagwilt, "script_for_juliette/new_go_output_from_deseqdrought_all/sagsurv_sagwilt.csv")

####ref "nemorensiswilting" ####
dds$group<- relevel(dds$group, ref = "nemorensiswilting")
levels(dds$group)
dds_test <- DESeq(dds)
resultsNames(dds_test)

sagwilt_nemwilt <- results(dds_test, name = "group_sagitattawilting_vs_nemorensiswilting")
#write.csv(sagwilt_nemwilt, "script_for_juliette/new_go_output_from_deseqdrought_all/sagwilt_nemwilt.csv")

nemsurv_nemwilt <- results(dds_test, name = "group_nemorensissurvival_vs_nemorensiswilting")
#write.csv(nemsurv_nemwilt, "script_for_juliette/new_go_output_from_deseqdrought_all/nemsurv_nemwilt.csv")


####ref "nemorensiscontrol" ####
dds$group<- relevel(dds$group, ref = "nemorensiscontrol")
levels(dds$group)
dds_test <- DESeq(dds)
resultsNames(dds_test)

sagctrl_nemctrl <- results(dds_test, name = "group_sagitattacontrol_vs_nemorensiscontrol")
#write.csv(sagctrl_nemctrl, "script_for_juliette/new_go_output_from_deseqdrought_all/sagctrl_nemctrl.csv")

####ref "sagitattacontrol" ####
dds$group<- relevel(dds$group, ref = "sagitattacontrol")
levels(dds$group)
dds_test <- DESeq(dds)
resultsNames(dds_test)

# Extract contrasts for control vs. wilting within each species
nemctrl_sagctrl <- results(dds_test, name = "group_nemorensiscontrol_vs_sagitattacontrol")
#write.csv(nemctrl_sagctrl, "script_for_juliette/new_go_output_from_deseqdrought_all/nemctrl_sagctrl.csv")

# Extract contrasts for control vs. wilting within each species
sagwilt_sagctrl <- results(dds_test, name = "group_sagitattawilting_vs_sagitattacontrol")
#write.csv(sagwilt_sagctrl, "script_for_juliette/new_go_output_from_deseqdrought_all/sagwilt_sagctrl.csv")

sagsurv_sagctrl <- results(dds_test, name = "group_sagitattasurvival_vs_sagitattacontrol")
#write.csv(sagsurv_sagctrl, "script_for_juliette/new_go_output_from_deseqdrought_all/sagsurv_sagctrl.csv")

####ref "sagitattasurvival" ####
dds$group<- relevel(dds$group, ref = "sagitattasurvival")
levels(dds$group)
dds_test <- DESeq(dds)
resultsNames(dds_test)

# Extract contrasts for control vs. wilting within each species
nemsurv_sagsurv <- results(dds_test, name = "group_nemorensissurvival_vs_sagitattasurvival")
#write.csv(nemsurv_sagsurv, "script_for_juliette/new_go_output_from_deseqdrought_all/nemsurv_sagsurv.csv")


####ref "nemorensissurvival" ####
dds$group<- relevel(dds$group, ref = "nemorensissurvival")
levels(dds$group)
dds_test <- DESeq(dds)
resultsNames(dds_test)

# Extract contrasts for control vs. wilting within each species
sagsurv_nemsurv <- results(dds_test, name = "group_sagitattasurvival_vs_nemorensissurvival")
#write.csv(sagsurv_nemsurv, "script_for_juliette/new_go_output_from_deseqdrought_all/sagsurv_nemsurv.csv")

####ref "nemorensiscontrol" ####
dds$group<- relevel(dds$group, ref = "nemorensiscontrol")
levels(dds$group)
dds_test <- DESeq(dds)
resultsNames(dds_test)

# Extract contrasts for control vs. wilting within each species
nemwilt_nemctrl <- results(dds_test, name = "group_nemorensiswilting_vs_nemorensiscontrol")
#write.csv(nemwilt_nemctrl, "script_for_juliette/new_go_output_from_deseqdrought_all/nemwilt_nemctrl.csv")

nemsurv_nemctrl <- results(dds_test, name = "group_nemorensissurvival_vs_nemorensiscontrol")
#write.csv(nemsurv_nemctrl, "script_for_juliette/new_go_output_from_deseqdrought_all/nemsurv_nemctrl.csv")



###################### volcano plots for A. sagittata in wilt and survival vs control ################

#### A. sagittata wilting vs. ctrl

# Prepare the volcano plot data from res_dds_filtered_sag
volcano_data_sag <- data.frame(
  gene = rownames(sagwilt_sagctrl),                         # Gene names
  log2FoldChange = sagwilt_sagctrl$log2FoldChange,          # Log2 fold change values
  log10pvalue = -log10(sagwilt_sagctrl$pvalue)              # -log10 of p-value
)

# Ensure the data doesn't contain infinite values due to log10 transformation
volcano_data_sag <- volcano_data_sag[is.finite(volcano_data_sag$log10pvalue), ]

# Define custom color breaks and labels
my_breaks <- c(-Inf, -10, -5, -1, 0, 1, 5, 10, Inf)
my_labels <- c("< -10", "-10 to -5", "-5 to -1", "-1 to 0", "0 to 1", "1 to 5", "5 to 10", "> 10")

# Define a custom color palette
my_colors <- c("#762a83", "#67001f", "#b2182b", "#d6604d", "#f4a582", "#92c5de", "#4393c3", "#2166ac")

# Create a categorical column for coloring based on log2 fold change
volcano_data_sag$L2FC <- cut(volcano_data_sag$log2FoldChange, breaks = my_breaks, labels = my_labels, include.lowest = TRUE)

# Define the midpoint for the color scale
mid_value <- median(volcano_data_sag$log2FoldChange, na.rm = TRUE)

# Plot the volcano plot using ggplot2
f_sag_wilt <- ggplot(volcano_data_sag, aes(x = log2FoldChange, y = log10pvalue, color = L2FC)) +
  geom_point(data = subset(volcano_data_sag, log10pvalue > -log10(0.05)), 
             shape = 20, 
             size = 1.5, 
             alpha = 0.7) +
  scale_color_manual(values = my_colors) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +    # Add dashed vertical lines for log2FC thresholds
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") + # Add dashed horizontal line for significance threshold
  labs(x = "Log2 Fold Change", y = "-Log10 P-value", fill = "Log2 Fold Change", title = expression(paste("Exp. regulation in wilting  - ", italic("A. sagittata")))) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 24),
        axis.text.x = element_text(size = 24),
        axis.title = element_text(size = 24),
        legend.text = element_text(size = 24),
        plot.title = element_text(size = 24),
        legend.title = element_text(size = 24))

# Add legend adjustments and set axis limits
g_sag_wilt <- f_sag_wilt + guides(color = guide_legend(override.aes = list(size = 5))) +
  ylim(0, 100) +  # Adjust y-axis limit for the plot
  xlim(-10, 10)   # Adjust x-axis limit for the plot

# Print the final plot for sag
print(g_sag_wilt)

#ggsave("script_for_juliette/new_go_output_from_deseqdrought_all/volcano_sag_wilt_ctrl_samples.png", plot = g_sag_wilt, width = 9, height = 7, dpi = 300)
#ggsave("script_for_juliette/new_go_output_from_deseqdrought_all/volcano_sag_wilt_ctrl_samples.pdf", plot = g_sag_wilt, width = 9, height = 7, dpi = 300)


##### A. sagittata survival vs. control


# Prepare the volcano plot data from res_dds_filtered_sag
volcano_data_sag <- data.frame(
  gene = rownames(sagsurv_sagctrl),                         # Gene names
  log2FoldChange = sagsurv_sagctrl$log2FoldChange,          # Log2 fold change values
  log10pvalue = -log10(sagsurv_sagctrl$pvalue)              # -log10 of p-value
)

# Ensure the data doesn't contain infinite values due to log10 transformation
volcano_data_sag <- volcano_data_sag[is.finite(volcano_data_sag$log10pvalue), ]

# Define custom color breaks and labels
my_breaks <- c(-Inf, -10, -5, -1, 0, 1, 5, 10, Inf)
my_labels <- c("< -10", "-10 to -5", "-5 to -1", "-1 to 0", "0 to 1", "1 to 5", "5 to 10", "> 10")

# Define a custom color palette
my_colors <- c("#762a83", "#67001f", "#b2182b", "#d6604d", "#f4a582", "#92c5de", "#4393c3", "#2166ac")

# Create a categorical column for coloring based on log2 fold change
volcano_data_sag$L2FC <- cut(volcano_data_sag$log2FoldChange, breaks = my_breaks, labels = my_labels, include.lowest = TRUE)

# Define the midpoint for the color scale
mid_value <- median(volcano_data_sag$log2FoldChange, na.rm = TRUE)

# Plot the volcano plot using ggplot2
f_sag_surv <- ggplot(volcano_data_sag, aes(x = log2FoldChange, y = log10pvalue, color = L2FC)) +
  geom_point(data = subset(volcano_data_sag, log10pvalue > -log10(0.05)), 
             shape = 20, 
             size = 1.5, 
             alpha = 0.7) +
  scale_color_manual(values = my_colors) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +    # Add dashed vertical lines for log2FC thresholds
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") + # Add dashed horizontal line for significance threshold
  labs(x = "Log2 Fold Change", y = "-Log10 P-value", fill = "Log2 Fold Change", title = expression(paste("Exp. regulation in survival  - ", italic("A. sagittata")))) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 24),
        axis.text.x = element_text(size = 24),
        axis.title = element_text(size = 24),
        legend.text = element_text(size = 24),
        plot.title = element_text(size = 24),
        legend.title = element_text(size = 24))

# Add legend adjustments and set axis limits
g_sag_surv <- f_sag_surv + guides(color = guide_legend(override.aes = list(size = 5))) +
  ylim(0, 100) +  # Adjust y-axis limit for the plot
  xlim(-10, 10)   # Adjust x-axis limit for the plot

# Print the final plot for sag
print(g_sag_surv)

#ggsave("script_for_juliette/new_go_output_from_deseqdrought_all/volcano_sag_surv_ctrl_samples.png", plot = g_sag_surv, width = 9, height = 7, dpi = 300)
#ggsave("script_for_juliette/new_go_output_from_deseqdrought_all/volcano_sag_surv_ctrl_samples.pdf", plot = g_sag_surv, width = 9, height = 7, dpi = 300)

############################################################################################

###################### volcano plots for A. nemorensis in wilt and survival vs control ################

#### A. nemorensis wilting vs. ctrl

# Prepare the volcano plot data from res_dds_filtered_sag
volcano_data_nem <- data.frame(
  gene = rownames(nemwilt_nemctrl),                         # Gene names
  log2FoldChange = nemwilt_nemctrl$log2FoldChange,          # Log2 fold change values
  log10pvalue = -log10(nemwilt_nemctrl$pvalue)              # -log10 of p-value
)

# Ensure the data doesn't contain infinite values due to log10 transformation
volcano_data_nem <- volcano_data_nem[is.finite(volcano_data_nem$log10pvalue), ]

# Define custom color breaks and labels
my_breaks <- c(-Inf, -10, -5, -1, 0, 1, 5, 10, Inf)
my_labels <- c("< -10", "-10 to -5", "-5 to -1", "-1 to 0", "0 to 1", "1 to 5", "5 to 10", "> 10")

# Define a custom color palette
my_colors <- c("#762a83", "#67001f", "#b2182b", "#d6604d", "#f4a582", "#92c5de", "#4393c3", "#2166ac")

# Create a categorical column for coloring based on log2 fold change
volcano_data_nem$L2FC <- cut(volcano_data_nem$log2FoldChange, breaks = my_breaks, labels = my_labels, include.lowest = TRUE)

# Define the midpoint for the color scale
mid_value <- median(volcano_data_nem$log2FoldChange, na.rm = TRUE)

# Plot the volcano plot using ggplot2
f_nem_wilt <- ggplot(volcano_data_nem, aes(x = log2FoldChange, y = log10pvalue, color = L2FC)) +
  geom_point(data = subset(volcano_data_nem, log10pvalue > -log10(0.05)), 
             shape = 20, 
             size = 1.5, 
             alpha = 0.7) +
  scale_color_manual(values = my_colors) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +    # Add dashed vertical lines for log2FC thresholds
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") + # Add dashed horizontal line for significance threshold
  labs(x = "Log2 Fold Change", y = "-Log10 P-value", fill = "Log2 Fold Change", title = expression(paste("Exp. regulation in wilting  - ", italic("A. nemorensis")))) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 24),
        axis.text.x = element_text(size = 24),
        axis.title = element_text(size = 24),
        legend.text = element_text(size = 24),
        plot.title = element_text(size = 24),
        legend.title = element_text(size = 24))

# Add legend adjustments and set axis limits
g_nem_wilt <- f_nem_wilt + guides(color = guide_legend(override.aes = list(size = 5))) +
  ylim(0, 100) +  # Adjust y-axis limit for the plot
  xlim(-10, 10)   # Adjust x-axis limit for the plot

# Print the final plot for nem
print(g_nem_wilt)

#ggsave("script_for_juliette/new_go_output_from_deseqdrought_all/volcano_nem_wilt_ctrl_samples.png", plot = g_nem_wilt, width = 9, height = 7, dpi = 300)
#ggsave("script_for_juliette/new_go_output_from_deseqdrought_all/volcano_nem_wilt_ctrl_samples.pdf", plot = g_nem_wilt, width = 9, height = 7, dpi = 300)


##### A. nemorensis survival vs. control


# Prepare the volcano plot data from res_dds_filtered_sag
volcano_data_nem <- data.frame(
  gene = rownames(nemsurv_nemctrl),                         # Gene names
  log2FoldChange = nemsurv_nemctrl$log2FoldChange,          # Log2 fold change values
  log10pvalue = -log10(nemsurv_nemctrl$pvalue)              # -log10 of p-value
)

# Ensure the data doesn't contain infinite values due to log10 transformation
volcano_data_nem <- volcano_data_nem[is.finite(volcano_data_nem$log10pvalue), ]

# Define custom color breaks and labels
my_breaks <- c(-Inf, -10, -5, -1, 0, 1, 5, 10, Inf)
my_labels <- c("< -10", "-10 to -5", "-5 to -1", "-1 to 0", "0 to 1", "1 to 5", "5 to 10", "> 10")

# Define a custom color palette
my_colors <- c("#762a83", "#67001f", "#b2182b", "#d6604d", "#f4a582", "#92c5de", "#4393c3", "#2166ac")

# Create a categorical column for coloring based on log2 fold change
volcano_data_nem$L2FC <- cut(volcano_data_nem$log2FoldChange, breaks = my_breaks, labels = my_labels, include.lowest = TRUE)

# Define the midpoint for the color scale
mid_value <- median(volcano_data_nem$log2FoldChange, na.rm = TRUE)

# Plot the volcano plot using ggplot2
f_nem_surv <- ggplot(volcano_data_nem, aes(x = log2FoldChange, y = log10pvalue, color = L2FC)) +
  geom_point(data = subset(volcano_data_nem, log10pvalue > -log10(0.05)), 
             shape = 20, 
             size = 1.5, 
             alpha = 0.7) +
  scale_color_manual(values = my_colors) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +    # Add dashed vertical lines for log2FC thresholds
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") + # Add dashed horizontal line for significance threshold
  labs(x = "Log2 Fold Change", y = "-Log10 P-value", fill = "Log2 Fold Change", title = expression(paste("Exp. regulation in survival  - ", italic("A. nemorensis")))) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 24),
        axis.text.x = element_text(size = 24),
        axis.title = element_text(size = 24),
        legend.text = element_text(size = 24),
        plot.title = element_text(size = 24),
        legend.title = element_text(size = 24))

# Add legend adjustments and set axis limits
g_nem_surv <- f_nem_surv + guides(color = guide_legend(override.aes = list(size = 5))) +
  ylim(0, 100) +  # Adjust y-axis limit for the plot
  xlim(-10, 10)   # Adjust x-axis limit for the plot

# Print the final plot for sag
print(g_nem_surv)

#ggsave("script_for_juliette/new_go_output_from_deseqdrought_all/volcano_nem_surv_ctrl_samples.png", plot = g_nem_surv, width = 9, height = 7, dpi = 300)
#ggsave("script_for_juliette/new_go_output_from_deseqdrought_all/volcano_nem_surv_ctrl_samples.pdf", plot = g_nem_surv, width = 9, height = 7, dpi = 300)

############################################################################################

######## PCA analysis of all samples

vsd <- vst(dds, blind=FALSE)

#####Principal component plot of the samples#####
plotPCA(vsd, intgroup=c("genotype", "condition"))

###It is also possible to customize the PCA plot using the ggplot function###
pcaData <- plotPCA(vsd, intgroup=c("genotype", "condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=genotype)) +
  geom_point(size=2) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

p <- ggplot(pcaData, aes(PC1, PC2, color = condition, shape = genotype)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() + theme(
    panel.background = element_rect(fill = "white", color = NA),  # Set the background color to white
    panel.grid.major = element_line(color = "gray", linewidth = 0.2), panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))

p


x <- ggplot(pcaData, aes(PC1, PC2, color = condition, shape = genotype)) +
  geom_point(size = 8, alpha = 0.5 ) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() + theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "gray", linewidth = 0.2),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.key.size = unit(2, "lines"),  # Adjust the size of the legend
    legend.text = element_text(size = 12),  # Adjust the font size of legend text
    legend.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 18), axis.title = element_text(size = 18)  # Adjust the font size and style of the legend title
  )


x 

#ggsave("pca_all_samples_transparent.png", plot = x, width = 9, height = 7, dpi = 300)
#ggsave("pca_all_samples_transparent.svg", plot = x, width = 9, height = 7, dpi = 300)
#ggsave("pca_all_samples_transparent.pdf", plot = x, width = 9, height = 7, dpi = 300)

#########################################################################
#########################################################################
###############                         #################################
############### GxE wilting and control #################################
###############                         #################################
#########################################################################
#########################################################################

### A nemorensis and A sagittata #####
setwd("/Users/Shared/Files From d.localized/PhD-Uni-Koeln_2021-2024/PhD_work/seeds data/2nd Experiment 2022/transcriptome/Arabis_RNA_raw_data/90-774066047/anaylsis_drought_mRNA_both_species_with_new_two_genomes/analysis_both_species_with_nem_genomes/GxE_ctrl_wilt/")
directory <- "/Users/Shared/Files From d.localized/PhD-Uni-Koeln_2021-2024/PhD_work/seeds data/2nd Experiment 2022/transcriptome/Arabis_RNA_raw_data/90-774066047/anaylsis_drought_mRNA_both_species_with_new_two_genomes/analysis_both_species_with_nem_genomes/GxE_ctrl_wilt/"
###-- define the pattern of files to be analysed, the file should end as ".txt" ######
sampleFiles <- grep("txt",list.files(directory),value=TRUE)
condition <- c("control", "control", "control", "control", "wilting", "wilting", "wilting", "wilting", "control", "control", "control", "control", "wilting", "wilting", "wilting", "wilting")
genotype <- c('nemorensis', 'nemorensis', 'nemorensis', 'nemorensis', 'nemorensis', 'nemorensis','nemorensis', 'nemorensis', 'sagitatta', 'sagitatta','sagitatta','sagitatta','sagitatta','sagitatta','sagitatta','sagitatta')
sampleTable <- data.frame(sampleName = sampleFiles, fileName = sampleFiles, condition = condition, genotype = genotype)
sampleTable$condition <- factor(sampleTable$condition)
sampleTable$genotype <- factor(sampleTable$genotype)

##--deseq_from_htseqcount######
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = directory, design = ~ genotype + condition + genotype:condition) #--(design= ~ genotype + condition + genotype:condition)--##

# Filter out genes with low counts across all samples
#--- get DEG ---#
#count_threshold <- 100  # Minimum total counts across samples
count_threshold <- 100  # Minimum average count threshold
num_samples <- 16       # Total number of samples
ddsHTSeq <- ddsHTSeq[rowSums(counts(ddsHTSeq)) / num_samples > count_threshold, ]
# Run the DESeq pipeline
dds <- DESeq(ddsHTSeq, fitType = "mean")
res_dds_gxe <- results(dds)

##Filtering genes with low counts_normalization ###
dds_lowcount <- estimateSizeFactors(dds)
sizeFactors(dds_lowcount)
normalized_counts <- counts(dds_lowcount, normalized = TRUE)

# Further filter res_dds_filtered to keep only those genes with padj < 0.05
res_dds_gxe_significant <- subset(res_dds_gxe, padj < 0.05)

#write.csv(res_dds_gxe_significant, "res_dds_gxe_wilt_ctrl_significant_3980.csv")

##########################################################################
##########################################################################
###############                          #################################
############### GxE survival and control #################################
###############                          #################################
##########################################################################
##########################################################################

### A nemorensis and A sagittata #####
setwd("/Users/Shared/Files From d.localized/PhD-Uni-Koeln_2021-2024/PhD_work/seeds data/2nd Experiment 2022/transcriptome/Arabis_RNA_raw_data/90-774066047/anaylsis_drought_mRNA_both_species_with_new_two_genomes/analysis_both_species_with_nem_genomes/GxE_ctrl_surv/")
directory <- "/Users/Shared/Files From d.localized/PhD-Uni-Koeln_2021-2024/PhD_work/seeds data/2nd Experiment 2022/transcriptome/Arabis_RNA_raw_data/90-774066047/anaylsis_drought_mRNA_both_species_with_new_two_genomes/analysis_both_species_with_nem_genomes/GxE_ctrl_surv/"
###-- define the pattern of files to be analysed, the file should end as ".txt" ######
sampleFiles <- grep("txt",list.files(directory),value=TRUE)
condition <- c("control", "control", "control", "control", "survival","survival","survival",  "control", "control", "control", "control", "survival","survival","survival")
genotype <- c('nemorensis', 'nemorensis', 'nemorensis', 'nemorensis', 'nemorensis','nemorensis', 'nemorensis', 'sagitatta', 'sagitatta','sagitatta','sagitatta','sagitatta','sagitatta','sagitatta')
sampleTable <- data.frame(sampleName = sampleFiles, fileName = sampleFiles, condition = condition, genotype = genotype)
sampleTable$condition <- factor(sampleTable$condition)
sampleTable$genotype <- factor(sampleTable$genotype)

##--deseq_from_htseqcount######
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = directory, design = ~ genotype + condition + genotype:condition) #--(design= ~ genotype + condition + genotype:condition)--##

# Filter out genes with low counts across all samples
#--- get DEG ---#
#count_threshold <- 100  # Minimum total counts across samples
count_threshold <- 100  # Minimum average count threshold
num_samples <- 14       # Total number of samples
ddsHTSeq <- ddsHTSeq[rowSums(counts(ddsHTSeq)) / num_samples > count_threshold, ]
# Run the DESeq pipeline
dds <- DESeq(ddsHTSeq, fitType = "mean")
res_dds_gxe_surv <- results(dds)

##Filtering genes with low counts_normalization ###
dds_lowcount <- estimateSizeFactors(dds)
sizeFactors(dds_lowcount)
normalized_counts <- counts(dds_lowcount, normalized = TRUE)

# Further filter res_dds_filtered to keep only those genes with padj < 0.05
res_dds_gxe_surv_significant <- subset(res_dds_gxe_surv, padj < 0.05)

#write.csv(res_dds_gxe_surv_significant, "res_dds_gxe_surv_ctrl_significant_1973.csv")




##########################################################################
##########################################################################
###############                          #################################
############### GxE wilting and Recovery #################################
###############                          #################################
##########################################################################
##########################################################################

### A nemorensis and A sagittata #####
setwd("/Users/Shared/Files From d.localized/PhD-Uni-Koeln_2021-2024/PhD_work/seeds data/2nd Experiment 2022/transcriptome/Arabis_RNA_raw_data/90-774066047/anaylsis_drought_mRNA_both_species_with_new_two_genomes/analysis_both_species_with_nem_genomes/GxE_surv_wilt/")
directory <- "/Users/Shared/Files From d.localized/PhD-Uni-Koeln_2021-2024/PhD_work/seeds data/2nd Experiment 2022/transcriptome/Arabis_RNA_raw_data/90-774066047/anaylsis_drought_mRNA_both_species_with_new_two_genomes/analysis_both_species_with_nem_genomes/GxE_surv_wilt/"
###-- define the pattern of files to be analysed, the file should end as ".txt" ######
sampleFiles <- grep("txt",list.files(directory),value=TRUE)
condition <- c("survival", "survival", "survival", "wilting", "wilting", "wilting", "wilting", "survival", "survival", "survival", "wilting", "wilting", "wilting", "wilting")
genotype <- c('nemorensis', 'nemorensis', 'nemorensis', 'nemorensis', 'nemorensis','nemorensis', 'nemorensis', 'sagitatta','sagitatta','sagitatta','sagitatta','sagitatta','sagitatta','sagitatta')
sampleTable <- data.frame(sampleName = sampleFiles, fileName = sampleFiles, condition = condition, genotype = genotype)
sampleTable$condition <- factor(sampleTable$condition)
sampleTable$genotype <- factor(sampleTable$genotype)

##--deseq_from_htseqcount######
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = directory, design = ~ genotype + condition + genotype:condition) #--(design= ~ genotype + condition + genotype:condition)--##

# Filter out genes with low counts across all samples
#--- get DEG ---#
#count_threshold <- 100  # Minimum total counts across samples
count_threshold <- 100  # Minimum average count threshold
num_samples <- 14       # Total number of samples
ddsHTSeq <- ddsHTSeq[rowSums(counts(ddsHTSeq)) / num_samples > count_threshold, ]
# Run the DESeq pipeline
dds <- DESeq(ddsHTSeq, fitType = "mean")
res_dds_gxe_surv_wilt <- results(dds)

##Filtering genes with low counts_normalization ###
dds_lowcount <- estimateSizeFactors(dds)
sizeFactors(dds_lowcount)
normalized_counts <- counts(dds_lowcount, normalized = TRUE)

# Further filter res_dds_filtered to keep only those genes with padj < 0.05
res_dds_gxe_surv_wilt_significant <- subset(res_dds_gxe_surv_wilt, padj < 0.05)

#write.csv(res_dds_gxe_significant, "res_dds_gxe_wilt_ctrl_significant_3980.csv")



##################################################################
####################### plot survival against wilt ################
# Ensure the gene IDs are aligned ####
####### corrected after thesis submission, 6th May 2025 , this is the final plot #######
# Ensure the gene IDs are aligned
common_genes_sag_nem_surv_wilt <- intersect(rownames(sagsurv_sagwilt), rownames(nemsurv_nemwilt))

# Subset to common genes
sag_surv_wilt_common <- sagsurv_sagwilt[common_genes_sag_nem_surv_wilt, ]
nem_surv_wilt_common <- nemsurv_nemwilt[common_genes_sag_nem_surv_wilt, ]


# Combine data from wilt and control contrasts
combined_species_df_sag_nem_surv_wilt <- data.frame(
  gene_id = common_genes_sag_nem_surv_wilt,
  log2FC_sag = sag_surv_wilt_common$log2FoldChange,
  padj_sag = sag_surv_wilt_common$padj,
  log2FC_nem = nem_surv_wilt_common$log2FoldChange,
  padj_nem = nem_surv_wilt_common$padj
)

# Add GxE information
gxe_common <- res_dds_gxe_surv_wilt[common_genes_sag_nem_surv_wilt, ]
combined_species_df_sag_nem_surv_wilt$padj_gxe <- gxe_common$padj

# Remove NA values from the combined_species_df_sag_nem
combined_species_df_sag_nem_surv_wilt <- na.omit(combined_species_df_sag_nem_surv_wilt)

# Filter out genes that do not differ in expression between species in at least one time point
#combined_species_df_sag_nem <- combined_species_df_sag_nem[
#  combined_species_df_sag_nem$padj_wilt < 0.05 | combined_species_df_sag_nem$padj_ctrl < 0.05, 
#]

# Redefine the color categories
combined_species_df_sag_nem_surv_wilt$color <- ifelse(
  combined_species_df_sag_nem_surv_wilt$padj_gxe < 0.05, "red", # Significant in GxE (padj < 0.001)
  ifelse(
    combined_species_df_sag_nem_surv_wilt$padj_sag < 0.05 & combined_species_df_sag_nem_surv_wilt$padj_nem < 0.05, "green", # Significant in both wilt and recovery
    "gray" # Non-significant in both
  )
)

#write.csv(combined_species_df_sag_nem, "combined_species_df_sag_nem_wilt_ctrl.csv")


# Plot with prioritized layers
p <- ggplot(combined_species_df_sag_nem_surv_wilt) +
  geom_point(data = subset(combined_species_df_sag_nem_surv_wilt, color == "gray"), aes(x = log2FC_nem, y = log2FC_sag, color = color), size = 0.7, alpha = 0.7) +
  geom_point(data = subset(combined_species_df_sag_nem_surv_wilt, color == "green"), aes(x = log2FC_nem, y = log2FC_sag, color = color), size = 0.7, alpha = 0.7) +
  geom_point(data = subset(combined_species_df_sag_nem_surv_wilt, color == "red"), aes(x = log2FC_nem, y = log2FC_sag, color = color), size = 0.7, alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  # Central vertical line
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Central horizontal line
  labs(
    y = expression(Log[2]~FC~italic("A. sagittata")~"(recov vs stress)"),
    x = expression(Log[2]~FC~italic("A. nemorensis")~"(recov vs stress)"),
    color = "Significance", title = expression(paste("Exp. diff. b/w ",italic("A. nemorensis"), " and", italic(" A. sagittata"), " at recovery"))
  ) +
  scale_color_manual(
    values = c("gray" = "gray", "green" = "darkgreen", "red" = "darkred"),
    labels = c("gray" = "NS", "green" = "Significant E ", "red" = "GxE (padj < 0.001)")
  ) +
  guides(
    color = guide_legend(
      override.aes = list(size = 4)  # Increase the size of legend dots
    )
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 24),
    axis.title = element_text(size = 24, face = "bold"),
    legend.text = element_text(size = 20),
    legend.title = element_blank(), plot.title = element_text(size= 26)
  ) +
  ylim(-10, 10) + xlim(-10, 10)

# Print the plot
print(p)

#ggsave("recov_vs_stress_species_comparison.png", plot = p, width = 11, height = 7, dpi = 300)
#ggsave("recov_vs_stress_species_comparison.pdf", plot = p, width = 11, height = 7, dpi = 300)
#ggsave("recov_vs_stress_species_comparison.svg", plot = p, width = 11, height = 7, dpi = 300)




##################################################################
####################### plot wilt against control ################
# Ensure the gene IDs are aligned ####
####### corrected after thesis submission, 6th May 2025 , this is the final plot #######
# Ensure the gene IDs are aligned
common_genes_sag_nem <- intersect(rownames(sagwilt_sagctrl), rownames(nemwilt_nemctrl))

# Subset to common genes
sag_common <- sagwilt_sagctrl[common_genes_sag_nem, ]
nem_common <- nemwilt_nemctrl[common_genes_sag_nem, ]


# Combine data from wilt and control contrasts
combined_species_df_sag_nem <- data.frame(
  gene_id = common_genes_sag_nem,
  log2FC_sag = sag_common$log2FoldChange,
  padj_sag = sag_common$padj,
  log2FC_nem = nem_common$log2FoldChange,
  padj_nem = nem_common$padj
)

# Add GxE information
gxe_common <- res_dds_gxe[common_genes_sag_nem, ]
combined_species_df_sag_nem$padj_gxe <- gxe_common$padj

# Remove NA values from the combined_species_df_sag_nem
combined_species_df_sag_nem <- na.omit(combined_species_df_sag_nem)

# Filter out genes that do not differ in expression between species in at least one time point
#combined_species_df_sag_nem <- combined_species_df_sag_nem[
#  combined_species_df_sag_nem$padj_wilt < 0.05 | combined_species_df_sag_nem$padj_ctrl < 0.05, 
#]

# Redefine the color categories
combined_species_df_sag_nem$color <- ifelse(
  combined_species_df_sag_nem$padj_gxe < 0.001, "red", # Significant in GxE (padj < 0.001)
  ifelse(
    combined_species_df_sag_nem$padj_sag < 0.05 & combined_species_df_sag_nem$padj_nem < 0.05, "green", # Significant in both wilt and control
    "gray" # Non-significant in both
  )
)

#write.csv(combined_species_df_sag_nem, "combined_species_df_sag_nem_wilt_ctrl.csv")


# Plot with prioritized layers
p <- ggplot(combined_species_df_sag_nem) +
  geom_point(data = subset(combined_species_df_sag_nem, color == "gray"), aes(x = log2FC_nem, y = log2FC_sag, color = color), size = 0.7, alpha = 0.7) +
  geom_point(data = subset(combined_species_df_sag_nem, color == "green"), aes(x = log2FC_nem, y = log2FC_sag, color = color), size = 0.7, alpha = 0.7) +
  geom_point(data = subset(combined_species_df_sag_nem, color == "red"), aes(x = log2FC_nem, y = log2FC_sag, color = color), size = 0.7, alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  # Central vertical line
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Central horizontal line
  labs(
    y = expression(Log[2]~FC~italic("A. sagittata")~"(stress vs control)"),
    x = expression(Log[2]~FC~italic("A. nemorensis")~"(stress vs control)"),
    color = "Significance", title = expression(paste("Exp. diff. b/w ",italic("A. nemorensis"), " and", italic(" A. sagittata"), " at stress"))
  ) +
  scale_color_manual(
    values = c("gray" = "gray", "green" = "darkgreen", "red" = "darkred"),
    labels = c("gray" = "NS", "green" = "Significant E ", "red" = "GxE (padj < 0.001)")
  ) +
  guides(
    color = guide_legend(
      override.aes = list(size = 4)  # Increase the size of legend dots
    )
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 24),
    axis.title = element_text(size = 24, face = "bold"),
    legend.text = element_text(size = 20),
    legend.title = element_blank(), plot.title = element_text(size= 26)
  ) +
  ylim(-10, 10) + xlim(-10, 10)

# Print the plot
print(p)

#ggsave("wilt_vs_ctrl_species_comparison.png", plot = p, width = 11, height = 7, dpi = 300)
#ggsave("wilt_vs_ctrl_species_comparison.pdf", plot = p, width = 11, height = 7, dpi = 300)
#ggsave("wilt_vs_ctrl_species_comparison.svg", plot = p, width = 11, height = 7, dpi = 300)




#######################################################################
####################### plot survival against control ################
# Ensure the gene IDs are aligned
common_genes_sag_nem_surv <- intersect(rownames(sagsurv_sagctrl), rownames(nemsurv_nemctrl))

# Subset to common genes
sag_common <- sagsurv_sagctrl[common_genes_sag_nem_surv, ]
nem_common <- nemsurv_nemctrl[common_genes_sag_nem_surv, ]


# Combine data from wilt and control contrasts
combined_species_df_sag_nem_surv <- data.frame(
  gene_id = common_genes_sag_nem_surv,
  log2FC_sag = sag_common$log2FoldChange,
  padj_sag = sag_common$padj,
  log2FC_nem = nem_common$log2FoldChange,
  padj_nem = nem_common$padj
)

# Add GxE information
gxe_common <- res_dds_gxe_surv[common_genes_sag_nem_surv, ]
combined_species_df_sag_nem_surv$padj_gxe <- gxe_common$padj

# Remove NA values from the combined_species_df_sag_nem_surv
combined_species_df_sag_nem_surv <- na.omit(combined_species_df_sag_nem_surv)


# Redefine the color categories
combined_species_df_sag_nem_surv$color <- ifelse(
  combined_species_df_sag_nem_surv$padj_gxe < 0.05, "red", # Significant in GxE (padj < 0.001)
  ifelse(
    combined_species_df_sag_nem_surv$padj_sag < 0.05 & combined_species_df_sag_nem_surv$padj_nem < 0.05, "green", # Significant in both wilt and control
    "gray" # Non-significant in both
  )
)

#write.csv(combined_species_df_sag_nem_surv, "combined_species_df_sag_nem_surv_ctrl.csv")

# Plot with prioritized layers
v <- ggplot(combined_species_df_sag_nem_surv) +
  geom_point(data = subset(combined_species_df_sag_nem_surv, color == "gray"), aes(x = log2FC_nem, y = log2FC_sag, color = color), size = 0.7, alpha = 0.7) +
  geom_point(data = subset(combined_species_df_sag_nem_surv, color == "green"), aes(x = log2FC_nem, y = log2FC_sag, color = color), size = 0.7, alpha = 0.7) +
  geom_point(data = subset(combined_species_df_sag_nem_surv, color == "red"), aes(x = log2FC_nem, y = log2FC_sag, color = color), size = 0.7, alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  # Central vertical line
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Central horizontal line
  labs(
    y = expression(Log[2]~FC~italic("A. sagittata")~"(recovery vs control)"),
    x = expression(Log[2]~FC~italic("A. nemorensis")~"(recovery vs control)"),
    color = "Significance", title = expression(paste("Exp. diff. b/w ",italic("A. nemorensis"), " and", italic(" A. sagittata"), " at recovery"))
  ) +
  scale_color_manual(
    values = c("gray" = "gray", "green" = "darkgreen", "red" = "darkred"),
    labels = c("gray" = "NS", "green" = "Significant E", "red" = "GxE (padj < 0.05)")
  ) +
  guides(
    color = guide_legend(
      override.aes = list(size = 4)  # Increase the size of legend dots
    )
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 24),
    axis.title = element_text(size = 24, face = "bold"),
    legend.text = element_text(size = 20),
    legend.title = element_blank(), plot.title = element_text(size= 26)
  ) +
  ylim(-10, 10) + xlim(-10, 10)

# Print the plot
print(v)

#setwd("/Users/Shared/Files From d.localized/PhD-Uni-Koeln_2021-2024/PhD_work/seeds data/2nd Experiment 2022/transcriptome/Arabis_RNA_raw_data/90-774066047/anaylsis_drought_mRNA_both_species_with_new_two_genomes/analysis_both_species_with_nem_genomes/script_for_juliette/new_go_output_from_deseqdrought_all/")
#ggsave("surv_vs_ctrl_species_comparison.png", plot = v, width = 11, height = 7, dpi = 300)
#ggsave("surv_vs_ctrl_species_comparison.pdf", plot = v, width = 11, height = 7, dpi = 300)
#ggsave("surv_vs_ctrl_species_comparison.svg", plot = v, width = 11, height = 7, dpi = 300)


######### SPLIT THE WILT_AND_CTRL OBJECT INTO QUADRANTS #########
# Add columns to classify points based on their positions relative to vline, hline, and diagonal
combined_species_df_sag_nem <- combined_species_df_sag_nem %>%
  mutate(
    above_diag = log2FC_sag > log2FC_nem,
    above_hline = log2FC_sag > 0,
    right_vline = log2FC_nem > 0,
    quadrant = case_when(
      above_diag & above_hline & right_vline ~ "Q1",
      above_diag & above_hline & !right_vline ~ "Q2",
      !above_diag & above_hline & !right_vline ~ "Q3",
      !above_diag & !above_hline & !right_vline ~ "Q4",
      !above_diag & !above_hline & right_vline ~ "Q5",
      above_diag & !above_hline & right_vline ~ "Q6",
      above_diag & !above_hline & !right_vline ~ "Q7",
      !above_diag & above_hline & right_vline ~ "Q8"
    )
  )

setwd("/Users/Shared/Files From d.localized/PhD-Uni-Koeln_2021-2024/PhD_work/seeds data/2nd Experiment 2022/transcriptome/Arabis_RNA_raw_data/90-774066047/anaylsis_drought_mRNA_both_species_with_new_two_genomes/analysis_both_species_with_nem_genomes/script_for_juliette/new_go_output_from_deseqdrought_all/")
#write.csv(combined_species_df_sag_nem, "combined_species_df_sag_nem_wilt_ctrl_quadrants.csv")
# Save each quadrant's data to a CSV file
#for (q in unique(combined_log2fc_df$quadrant)) {
#  quad_data <- combined_log2fc_df %>% filter(quadrant == q)
#  write.csv(quad_data, paste0("quadrant_", q, ".csv"), row.names = FALSE)
#}


######### SPLIT THE REC_AND_CTRL OBJECT INTO QUADRANTS #########
# Add columns to classify points based on their positions relative to vline, hline, and diagonal
combined_species_df_sag_nem_surv <- combined_species_df_sag_nem_surv %>%
  mutate(
    above_diag = log2FC_sag > log2FC_nem,
    above_hline = log2FC_sag > 0,
    right_vline = log2FC_nem > 0,
    quadrant = case_when(
      above_diag & above_hline & right_vline ~ "Q1",
      above_diag & above_hline & !right_vline ~ "Q2",
      !above_diag & above_hline & !right_vline ~ "Q3",
      !above_diag & !above_hline & !right_vline ~ "Q4",
      !above_diag & !above_hline & right_vline ~ "Q5",
      above_diag & !above_hline & right_vline ~ "Q6",
      above_diag & !above_hline & !right_vline ~ "Q7",
      !above_diag & above_hline & right_vline ~ "Q8"
    )
  )

#write.csv(combined_species_df_sag_nem_surv, "combined_species_df_sag_nem_surv_ctrl_quadrants.csv")



######### SPLIT THE REC_AND_CTRL OBJECT INTO QUADRANTS #########
# Add columns to classify points based on their positions relative to vline, hline, and diagonal
combined_species_df_sag_nem_surv_wilt <- combined_species_df_sag_nem_surv_wilt %>%
  mutate(
    above_diag = log2FC_sag > log2FC_nem,
    above_hline = log2FC_sag > 0,
    right_vline = log2FC_nem > 0,
    quadrant = case_when(
      above_diag & above_hline & right_vline ~ "Q1",
      above_diag & above_hline & !right_vline ~ "Q2",
      !above_diag & above_hline & !right_vline ~ "Q3",
      !above_diag & !above_hline & !right_vline ~ "Q4",
      !above_diag & !above_hline & right_vline ~ "Q5",
      above_diag & !above_hline & right_vline ~ "Q6",
      above_diag & !above_hline & !right_vline ~ "Q7",
      !above_diag & above_hline & right_vline ~ "Q8"
    )
  )

#write.csv(combined_species_df_sag_nem_surv_wilt, "New_GO_analysis_6thMay2025/combined_species_df_sag_nem_surv_wilt_quadrants.csv")




######## GO TEST analysis (wilting)
# Load necessary libraries
library(BiocManager)
library(KEGGREST)
library(org.At.tair.db)
library(Rgraphviz)
library(topGO)
library(biomaRt)
library(ggplot2)
library(AnnotationDbi)
library(clusterProfiler)
require(ggplot2)
library(scales)
library(topGO)
library(dplyr)
library(readr)

#########################################################-------------
########## GO wilting and control
#Import the orthologues CSV file
# ========================== STEP 1: Load Data ==========================
orthologues <- read.csv("/Users/Shared/Files From d.localized/PhD-Uni-Koeln_2021-2024/PhD_work/seeds data/2nd Experiment 2022/transcriptome/Arabis_RNA_raw_data/90-774066047/anaylsis_drought_mRNA_both_species_with_new_two_genomes/analysis_both_species_with_nem_genomes/orthologues_cleaned.csv", header = TRUE)
new_dataframe1 <- data.frame(orthologues)

# Merge expression data with orthologue mapping
combined_species_df_sag_nem <- combined_species_df_sag_nem %>%
  left_join(new_dataframe1, by = c("gene_id" = "arabis_cleaned"))

#write.csv(combined_species_df_sag_nem, "combined_species_df_sag_nem_with_ortho.csv")


# Remove genes without an orthologue
combined_species_df_sag_nem <- combined_species_df_sag_nem %>%
  filter(!is.na(At))


# Merge expression data with orthologue mapping
combined_species_df_sag_nem_surv <- combined_species_df_sag_nem_surv %>%
  left_join(new_dataframe1, by = c("gene_id" = "arabis_cleaned"))

#write.csv(combined_species_df_sag_nem_surv, "combined_species_df_sag_nem_surv_with_ortho.csv")


# Remove genes without an orthologue
combined_species_df_sag_nem_surv <- combined_species_df_sag_nem_surv %>%
  filter(!is.na(At))


# Fix column names like At.x, At.x.x, At.y etc.
#colnames(combined_species_df_sag_nem_surv) <- gsub("\\.x+|\\.y+", "", colnames(combined_species_df_sag_nem_surv))

#combined_species_df_sag_nem_surv <- combined_species_df_sag_nem_surv[ , !grepl("^At\\.\\d+$", colnames(combined_species_df_sag_nem_surv))]

#colnames(combined_species_df_sag_nem_surv) <- make.unique(colnames(combined_species_df_sag_nem_surv))

# ========================== STEP 2: Define Universes ==========================
### split on diagonal
##### A. sag and A. nem (stress)
# Genes ABOVE the diagonal
upper_universe_sag_nem <- combined_species_df_sag_nem %>%
  filter(log2FC_sag > log2FC_nem) 

# Genes BELOW the diagonal
lower_universe_sag_nem <- combined_species_df_sag_nem %>%
  filter(log2FC_sag < log2FC_nem)

##### A. sag and A. nem (recovery)
# Genes ABOVE the diagonal
upper_universe_sag_nem_surv <- combined_species_df_sag_nem_surv %>%
  filter(log2FC_sag > log2FC_nem) 

# Genes BELOW the diagonal
lower_universe_sag_nem_surv <- combined_species_df_sag_nem_surv %>%
  filter(log2FC_sag < log2FC_nem)

##### split on vertical line
##### A. sag and A. nem (stress)
# Genes Left the vline
left_universe_sag_nem <- combined_species_df_sag_nem %>%
  filter(log2FC_nem < 0) 

# Genes right the vline
right_universe_sag_nem <- combined_species_df_sag_nem %>%
  filter(log2FC_nem > 0)

##### A. sag and A. nem (recovery)
# Genes ABOVE the diagonal
left_universe_sag_nem_surv <- combined_species_df_sag_nem_surv %>%
  filter(log2FC_nem < 0) 

# Genes BELOW the diagonal
right_universe_sag_nem_surv <- combined_species_df_sag_nem_surv %>%
  filter(log2FC_nem > 0)

#write.csv(left_universe_sag_nem, "left_universe_sag_nem_wilt.csv")
#write.csv(right_universe_sag_nem, "right_universe_sag_nem_wilt.csv")


#write.csv(left_universe_sag_nem_surv, "left_universe_sag_nem_recovery.csv")
#write.csv(right_universe_sag_nem_surv, "right_universe_sag_nem_recovery.csv")


# ========================== STEP 3: Perform GO Enrichments ==========================
setwd("/Users/Shared/Files From d.localized/PhD-Uni-Koeln_2021-2024/PhD_work/seeds data/2nd Experiment 2022/transcriptome/Arabis_RNA_raw_data/90-774066047/anaylsis_drought_mRNA_both_species_with_new_two_genomes/analysis_both_species_with_nem_genomes/script_for_juliette/new_go_output_from_deseqdrought_all/New_GO_analysis_6thMay2025//")
## First we checked shared response in both Arabis species to drought and recovery, 
# so first we will look for the similar response in both species or stress response mechanisms
# in both species together.
##### left universrse ##### A. nem and A. sag in drought
allGenes_numeric <- ifelse((left_universe_sag_nem$log2FC_nem < 0 & left_universe_sag_nem$log2FC_sag < 0 ) & left_universe_sag_nem$color == "green", 0, 1) ## check for common enrichment in green genes.

names(allGenes_numeric) <- left_universe_sag_nem$At

# Remove NA values if necessary
allGenes_numeric <- allGenes_numeric[!is.na(allGenes_numeric)]

head(allGenes_numeric)

# Create topGO data object
tGOdata_sag_nem_left_wilt_ctrl <- new("topGOdata",
                            description = "Enrichment Analysis for Q1",
                            ontology = "BP",  
                            allGenes = allGenes_numeric,
                            geneSel = function(x) x == 0,  # This marks genes of interest (green in Q1) as TRUE
                            nodeSize = 5,  # Minimum number of genes for a GO term to be considered
                            mapping = "org.At.tair.db",  # database for Arabidopsis thaliana
                            annot = annFUN.org)

# enrichment test KS (KS test targets specific and significant enrichment)
results.fisher_sag_nem_left_wilt_ctrl <- runTest(tGOdata_sag_nem_left_wilt_ctrl, algorithm="elim", statistic="fisher")

goEnrichmentQ_shared_sag_nem_left_wilt_ctrl <- GenTable(tGOdata_sag_nem_left_wilt_ctrl, fisher=results.fisher_sag_nem_left_wilt_ctrl, orderBy="fisher", topNodes=50)


# Convert p-values to numeric (handling "<1e-30" etc.)
goEnrichmentQ_shared_sag_nem_left_wilt_ctrl$fisher <- suppressWarnings(as.numeric(goEnrichmentQ_shared_sag_nem_left_wilt_ctrl$fisher))
goEnrichmentQ_shared_sag_nem_left_wilt_ctrl$fisher[is.na(goEnrichmentQ_shared_sag_nem_left_wilt_ctrl$fisher)] <- 1e-30  # Replace NAs with small value

# Add -log10 p-value for color gradient
goEnrichmentQ_shared_sag_nem_left_wilt_ctrl$log10_p <- -log10(goEnrichmentQ_shared_sag_nem_left_wilt_ctrl$fisher)

# Optional: reorder GO terms by significance
goEnrichmentQ_shared_sag_nem_left_wilt_ctrl$Term <- factor(goEnrichmentQ_shared_sag_nem_left_wilt_ctrl$Term, levels = rev(goEnrichmentQ_shared_sag_nem_left_wilt_ctrl$Term))


ggplot(goEnrichmentQ_shared_sag_nem_left_wilt_ctrl, aes(x = as.numeric(Significant), 
                             y = Term,
                             fill = log10_p)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "lightblue", high = "red", name = "-log10(p)") +
  labs(x = "Number of Significant Genes in table_shared_sag_nem_left_wilt_ctrl",
       y = "GO Term",
       title = "GO Enrichment (elimFisher)") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 14, face = "bold"))


#write.csv(go_table_NC_DOWN, "DEG_groups/go_table_NC_DOWN.csv")

# Get gene IDs for enriched GO terms
goEnrichmentQ_shared_sag_nem_left_wilt_ctrl_ID <- goEnrichmentQ_shared_sag_nem_left_wilt_ctrl$GO.ID

# Get significant genes per GO term
# Corwiltted function
getSigGenesFromGO <- function(GOterm, GOdata) {
  all_genes <- genesInTerm(GOdata)[[GOterm]]
  sig_genes <- all_genes[all_genes %in% sigGenes(GOdata)]
  return(sig_genes)
}

goEnrichmentQ_shared_sag_nem_left_wilt_ctrl_genes <- lapply(goEnrichmentQ_shared_sag_nem_left_wilt_ctrl_ID, getSigGenesFromGO, GOdata = tGOdata_sag_nem_left_wilt_ctrl)
names(goEnrichmentQ_shared_sag_nem_left_wilt_ctrl_genes) <- goEnrichmentQ_shared_sag_nem_left_wilt_ctrl_ID

# Flatten to a data frame
goEnrichmentQ_shared_sag_nem_left_wilt_ctrl_genes_df <- stack(goEnrichmentQ_shared_sag_nem_left_wilt_ctrl_genes)
colnames(goEnrichmentQ_shared_sag_nem_left_wilt_ctrl_genes_df) <- c("GeneID", "GO.ID")

# Save
#write.csv(goEnrichmentQ_shared_sag_nem_left_wilt_ctrl_genes_df, "goEnrichmentQ_shared_sag_nem_left_wilt_ctrl_genes.csv", row.names = FALSE)


# Collapse all gene IDs per GO term into a single string
goEnrichmentQ_shared_sag_nem_left_wilt_ctrl_genes_lists_by_go <- goEnrichmentQ_shared_sag_nem_left_wilt_ctrl_genes_df %>%
  group_by(GO.ID) %>%
  summarise(GeneList = paste(GeneID, collapse = ", "), .groups = "drop")


# Merge gene list into the main GO table by GO.ID
goEnrichmentQ_shared_sag_nem_left_wilt_ctrl_annotated <- left_join(goEnrichmentQ_shared_sag_nem_left_wilt_ctrl, goEnrichmentQ_shared_sag_nem_left_wilt_ctrl_genes_lists_by_go, by = "GO.ID")


#write.csv(goEnrichmentQ_shared_sag_nem_left_wilt_ctrl_annotated, "goEnrichmentQ_shared_sag_nem_left_wilt_ctrl_genes_annotated.csv", row.names = FALSE)

#write.csv(goEnrichmentQ_shared_sag_nem_left_wilt_ctrl, "goEnrichmentQ_shared_sag_nem_left_wilt_ctrl.csv")



####################
# -------------------------------
# 1. GO enrichment table (already created)
# -------------------------------

goEnrichmentQ_shared_sag_nem_left_wilt_ctrl <- GenTable(
  tGOdata_sag_nem_left_wilt_ctrl,
  fisher = results.fisher_sag_nem_left_wilt_ctrl,
  orderBy = "fisher",
  topNodes = 50
)

# Convert p-values to numeric
goEnrichmentQ_shared_sag_nem_left_wilt_ctrl$fisher <- suppressWarnings(
  as.numeric(goEnrichmentQ_shared_sag_nem_left_wilt_ctrl$fisher)
)
goEnrichmentQ_shared_sag_nem_left_wilt_ctrl$fisher[is.na(goEnrichmentQ_shared_sag_nem_left_wilt_ctrl$fisher)] <- 1e-30

# Add -log10 p-value
goEnrichmentQ_shared_sag_nem_left_wilt_ctrl$log10_p <-
  -log10(goEnrichmentQ_shared_sag_nem_left_wilt_ctrl$fisher)

# -------------------------------
# 2. Function to extract annotated & significant genes
# -------------------------------

getGenesFromGO <- function(GOterm, GOdata) {
  
  annotated_genes <- genesInTerm(GOdata)[[GOterm]]
  significant_genes <- intersect(annotated_genes, sigGenes(GOdata))
  
  tibble(
    GO.ID = GOterm,
    AnnotatedGenes = paste(annotated_genes, collapse = ", "),
    SignificantGenes = paste(significant_genes, collapse = ", "),
    n_annotated = length(annotated_genes),
    n_significant = length(significant_genes)
  )
}

# -------------------------------
# 3. Apply function to all enriched GO terms
# -------------------------------

go_ids <- goEnrichmentQ_shared_sag_nem_left_wilt_ctrl$GO.ID

genes_per_go <- bind_rows(
  lapply(go_ids, getGenesFromGO, GOdata = tGOdata_sag_nem_left_wilt_ctrl)
)

# -------------------------------
# 4. Merge gene lists back into GO table
# -------------------------------

goEnrichmentQ_shared_sag_nem_left_wilt_ctrl_annotated1 <- left_join(
  goEnrichmentQ_shared_sag_nem_left_wilt_ctrl,
  genes_per_go,
  by = "GO.ID"
)

# -------------------------------
# 5. Optional: long-format gene table (one gene per row)
# -------------------------------

#getGenesLong <- function(GOterm, GOdata) {
#  
#  annotated_genes <- genesInTerm(GOdata)[[GOterm]]
#  significant_genes <- intersect(annotated_genes, sigGenes(GOdata))
  
#  tibble(
#    GO.ID = GOterm,
#    GeneID = annotated_genes,
#    IsSignificant = annotated_genes %in% significant_genes
#  )
#}

#go_genes_long <- bind_rows(
#  lapply(go_ids, getGenesLong, GOdata = tGOdata_sag_nem_left_wilt_ctrl)
#)

# -------------------------------
# 6. Save results (optional)
# -------------------------------

# write.csv(goEnrichmentQ_shared_sag_nem_left_wilt_ctrl_annotated1,
#           "goEnrichmentQ_shared_sag_nem_left_wilt_ctrl_annotated_and_significant.csv",
#           row.names = FALSE)

# write.csv(go_genes_long,
#           "GO_annotated_genes_long_format.csv",
#           row.names = FALSE)

###################


######################################
##### right universrse ##### A. nem and A. sag in drought
allGenes_numeric <- ifelse((right_universe_sag_nem$log2FC_nem > 0 & right_universe_sag_nem$log2FC_sag > 0 ) & right_universe_sag_nem$color == "green", 0, 1) ## check for common enrichment in green genes.

names(allGenes_numeric) <- right_universe_sag_nem$At

# Remove NA values if necessary
allGenes_numeric <- allGenes_numeric[!is.na(allGenes_numeric)]

head(allGenes_numeric)

# Create topGO data object
tGOdata_sag_nem_right_wilt_ctrl <- new("topGOdata",
                                      description = "Enrichment Analysis for Q1",
                                      ontology = "BP",  
                                      allGenes = allGenes_numeric,
                                      geneSel = function(x) x == 0,  # This marks genes of interest (green in Q1) as TRUE
                                      nodeSize = 5,  # Minimum number of genes for a GO term to be considered
                                      mapping = "org.At.tair.db",  # database for Arabidopsis thaliana
                                      annot = annFUN.org)

# enrichment test KS (KS test targets specific and significant enrichment)
results.fisher_sag_nem_right_wilt_ctrl <- runTest(tGOdata_sag_nem_right_wilt_ctrl, algorithm="elim", statistic="fisher")

goEnrichmentQ_shared_sag_nem_right_wilt_ctrl <- GenTable(tGOdata_sag_nem_right_wilt_ctrl, fisher=results.fisher_sag_nem_right_wilt_ctrl, orderBy="fisher", topNodes=20)

#write.csv(goEnrichmentQ_shared_sag_nem_right_wilt_ctrl, "goEnrichmentQ_shared_sag_nem_right_wilt_ctrl.csv")

# Convert p-values to numeric (handling "<1e-30" etc.)
goEnrichmentQ_shared_sag_nem_right_wilt_ctrl$fisher <- suppressWarnings(as.numeric(goEnrichmentQ_shared_sag_nem_right_wilt_ctrl$fisher))
goEnrichmentQ_shared_sag_nem_right_wilt_ctrl$fisher[is.na(goEnrichmentQ_shared_sag_nem_right_wilt_ctrl$fisher)] <- 1e-30  # Replace NAs with small value

# Add -log10 p-value for color gradient
goEnrichmentQ_shared_sag_nem_right_wilt_ctrl$log10_p <- -log10(goEnrichmentQ_shared_sag_nem_right_wilt_ctrl$fisher)

# Optional: reorder GO terms by significance
goEnrichmentQ_shared_sag_nem_right_wilt_ctrl$Term <- factor(goEnrichmentQ_shared_sag_nem_right_wilt_ctrl$Term, levels = rev(goEnrichmentQ_shared_sag_nem_right_wilt_ctrl$Term))


ggplot(goEnrichmentQ_shared_sag_nem_right_wilt_ctrl, aes(x = as.numeric(Significant), 
                             y = Term,
                             fill = log10_p)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "lightblue", high = "red", name = "-log10(p)") +
  labs(x = "Number of Significant Genes in table_shared_sag_nem_right_wilt_ctrl",
       y = "GO Term",
       title = "GO Enrichment (elimFisher)") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 14, face = "bold"))


#write.csv(go_table_NC_DOWN, "DEG_groups/go_table_NC_DOWN.csv")


# Get gene IDs for enriched GO terms
goEnrichmentQ_shared_sag_nem_right_wilt_ctrl_ID <- goEnrichmentQ_shared_sag_nem_right_wilt_ctrl$GO.ID

# Get significant genes per GO term
# Corwiltted function
getSigGenesFromGO <- function(GOterm, GOdata) {
  all_genes <- genesInTerm(GOdata)[[GOterm]]
  sig_genes <- all_genes[all_genes %in% sigGenes(GOdata)]
  return(sig_genes)
}

goEnrichmentQ_shared_sag_nem_right_wilt_ctrl_genes <- lapply(goEnrichmentQ_shared_sag_nem_right_wilt_ctrl_ID, getSigGenesFromGO, GOdata = tGOdata_sag_nem_right_wilt_ctrl)
names(goEnrichmentQ_shared_sag_nem_right_wilt_ctrl_genes) <- goEnrichmentQ_shared_sag_nem_right_wilt_ctrl_ID

# Flatten to a data frame
goEnrichmentQ_shared_sag_nem_right_wilt_ctrl_genes_df <- stack(goEnrichmentQ_shared_sag_nem_right_wilt_ctrl_genes)
colnames(goEnrichmentQ_shared_sag_nem_right_wilt_ctrl_genes_df) <- c("GeneID", "GO.ID")

# Save
#write.csv(goEnrichmentQ_shared_sag_nem_right_wilt_ctrl_genes_df, "goEnrichmentQ_shared_sag_nem_right_wilt_ctrl_genes.csv", row.names = FALSE)


# Collapse all gene IDs per GO term into a single string
goEnrichmentQ_shared_sag_nem_right_wilt_ctrl_genes_lists_by_go <- goEnrichmentQ_shared_sag_nem_right_wilt_ctrl_genes_df %>%
  group_by(GO.ID) %>%
  summarise(GeneList = paste(GeneID, collapse = ", "), .groups = "drop")


# Merge gene list into the main GO table by GO.ID
goEnrichmentQ_shared_sag_nem_right_wilt_ctrl_annotated <- left_join(goEnrichmentQ_shared_sag_nem_right_wilt_ctrl, goEnrichmentQ_shared_sag_nem_right_wilt_ctrl_genes_lists_by_go, by = "GO.ID")


#write.csv(goEnrichmentQ_shared_sag_nem_right_wilt_ctrl_annotated, "goEnrichmentQ_shared_sag_nem_right_wilt_ctrl_genes_annotated.csv", row.names = FALSE)


####################
# -------------------------------
# 1. GO enrichment table (already created)
# -------------------------------

goEnrichmentQ_shared_sag_nem_right_wilt_ctrl <- GenTable(
  tGOdata_sag_nem_right_wilt_ctrl,
  fisher = results.fisher_sag_nem_right_wilt_ctrl,
  orderBy = "fisher",
  topNodes = 50
)

# Convert p-values to numeric
goEnrichmentQ_shared_sag_nem_right_wilt_ctrl$fisher <- suppressWarnings(
  as.numeric(goEnrichmentQ_shared_sag_nem_right_wilt_ctrl$fisher)
)
goEnrichmentQ_shared_sag_nem_right_wilt_ctrl$fisher[is.na(goEnrichmentQ_shared_sag_nem_right_wilt_ctrl$fisher)] <- 1e-30

# Add -log10 p-value
goEnrichmentQ_shared_sag_nem_right_wilt_ctrl$log10_p <-
  -log10(goEnrichmentQ_shared_sag_nem_right_wilt_ctrl$fisher)

# -------------------------------
# 2. Function to extract annotated & significant genes
# -------------------------------

getGenesFromGO <- function(GOterm, GOdata) {
  
  annotated_genes <- genesInTerm(GOdata)[[GOterm]]
  significant_genes <- intersect(annotated_genes, sigGenes(GOdata))
  
  tibble(
    GO.ID = GOterm,
    AnnotatedGenes = paste(annotated_genes, collapse = ", "),
    SignificantGenes = paste(significant_genes, collapse = ", "),
    n_annotated = length(annotated_genes),
    n_significant = length(significant_genes)
  )
}

# -------------------------------
# 3. Apply function to all enriched GO terms
# -------------------------------

go_ids <- goEnrichmentQ_shared_sag_nem_right_wilt_ctrl$GO.ID

genes_per_go <- bind_rows(
  lapply(go_ids, getGenesFromGO, GOdata = tGOdata_sag_nem_right_wilt_ctrl)
)

# -------------------------------
# 4. Merge gene lists back into GO table
# -------------------------------

goEnrichmentQ_shared_sag_nem_right_wilt_ctrl_annotated1 <- left_join(
  goEnrichmentQ_shared_sag_nem_right_wilt_ctrl,
  genes_per_go,
  by = "GO.ID"
)

# -------------------------------
# 5. Optional: long-format gene table (one gene per row)
# -------------------------------

#getGenesLong <- function(GOterm, GOdata) {
#  
#  annotated_genes <- genesInTerm(GOdata)[[GOterm]]
#  significant_genes <- intersect(annotated_genes, sigGenes(GOdata))

#  tibble(
#    GO.ID = GOterm,
#    GeneID = annotated_genes,
#    IsSignificant = annotated_genes %in% significant_genes
#  )
#}

#go_genes_long <- bind_rows(
#  lapply(go_ids, getGenesLong, GOdata = tGOdata_sag_nem_left_wilt_ctrl)
#)

# -------------------------------
# 6. Save results (optional)
# -------------------------------

# write.csv(goEnrichmentQ_shared_sag_nem_right_wilt_ctrl_annotated1,
#            "goEnrichmentQ_shared_sag_nem_right_wilt_ctrl_annotated_and_significant.csv",
#            row.names = FALSE)

# write.csv(go_genes_long,
#           "GO_annotated_genes_long_format.csv",
#           row.names = FALSE)

###################

######################################
##### A. sag and A. nem (stress)
##### left universrse
# for down regulated in sag use Q7, for more down in sag use Q4, in nem use Q4, more down in nem use Q7)
allGenes_numeric <- ifelse((left_universe_sag_nem$quadrant == "Q2") & left_universe_sag_nem$color == "red", 0, 1)
names(allGenes_numeric) <- left_universe_sag_nem$At

# Remove NA values if necessary
allGenes_numeric <- allGenes_numeric[!is.na(allGenes_numeric)]

head(allGenes_numeric)

# Create topGO2 data object
tGOdata2 <- new("topGOdata",
               description = "Enrichment Analysis for Q2",
               ontology = "BP",  
               allGenes = allGenes_numeric,
               geneSel = function(x) x == 0,  # This marks genes of interest (red in Q1) as TRUE
               nodeSize = 10,  # Minimum number of genes for a GO term to be considered
               mapping = "org.At.tair.db",  # database for Arabidopsis thaliana
               annot = annFUN.org)

# enrichment test KS (KS test targets specific and significant enrichment)
results.fisher2 <- runTest(tGOdata2, algorithm="elim", statistic="fisher")
goEnrichmentQ_sag_up_against_nem_wiltQ2 <- GenTable(tGOdata2, KS=results.fisher2, orderBy="KS", topNodes=50)


##################################
# for down regulated in sag use Q7, for more down in sag use Q4, in nem use Q4, more down in nem use Q7)
allGenes_numeric <- ifelse((left_universe_sag_nem$quadrant == "Q4") & left_universe_sag_nem$color == "red", 0, 1)
names(allGenes_numeric) <- left_universe_sag_nem$At

# Remove NA values if necessary
allGenes_numeric <- allGenes_numeric[!is.na(allGenes_numeric)]

head(allGenes_numeric)

# Create topGO4 data object
tGOdata4 <- new("topGOdata",
                description = "Enrichment Analysis for Q4",
                ontology = "BP",  
                allGenes = allGenes_numeric,
                geneSel = function(x) x == 0,  # This marks genes of interest (red in Q1) as TRUE
                nodeSize = 10,  # Minimum number of genes for a GO term to be considered
                mapping = "org.At.tair.db",  # database for Arabidopsis thaliana
                annot = annFUN.org)

# enrichment test KS (KS test targets specific and significant enrichment)
results.fisher4 <- runTest(tGOdata4, algorithm="elim", statistic="fisher")
goEnrichmentQ_sag_down_more_against_nem_wiltQ4 <- GenTable(tGOdata4, KS=results.fisher4, orderBy="KS", topNodes=50)


######################################
# for down regulated in sag use Q7, for more down in sag use Q4, in nem use Q4, more down in nem use Q7)
allGenes_numeric <- ifelse((left_universe_sag_nem$quadrant == "Q7") & left_universe_sag_nem$color == "red", 0, 1)
names(allGenes_numeric) <- left_universe_sag_nem$At

# Remove NA values if necessary
allGenes_numeric <- allGenes_numeric[!is.na(allGenes_numeric)]

head(allGenes_numeric)

# Create topGO7 data object
tGOdata7 <- new("topGOdata",
                description = "Enrichment Analysis for Q7",
                ontology = "BP",  
                allGenes = allGenes_numeric,
                geneSel = function(x) x == 0,  # This marks genes of interest (red in Q1) as TRUE
                nodeSize = 10,  # Minimum number of genes for a GO term to be considered
                mapping = "org.At.tair.db",  # database for Arabidopsis thaliana
                annot = annFUN.org)

# enrichment test KS (KS test targets specific and significant enrichment)
results.fisher7 <- runTest(tGOdata7, algorithm="elim", statistic="fisher")



###### right universe
# for up regulated in nem use Q5, for more up in sag use Q1, more up in nem use Q8)
allGenes_numeric <- ifelse((right_universe_sag_nem$quadrant == "Q1") & right_universe_sag_nem$color == "red", 0, 1)
names(allGenes_numeric) <- right_universe_sag_nem$At

# Remove NA values if necessary
allGenes_numeric <- allGenes_numeric[!is.na(allGenes_numeric)]

head(allGenes_numeric)

# Create topGO2 data object
tGOdata1 <- new("topGOdata",
                description = "Enrichment Analysis for Q1",
                ontology = "BP",  
                allGenes = allGenes_numeric,
                geneSel = function(x) x == 0,  # This marks genes of interest (red in Q1) as TRUE
                nodeSize = 10,  # Minimum number of genes for a GO term to be considered
                mapping = "org.At.tair.db",  # database for Arabidopsis thaliana
                annot = annFUN.org)

# enrichment test KS (KS test targets specific and significant enrichment)
results.fisher1 <- runTest(tGOdata1, algorithm="elim", statistic="fisher")

goEnrichmentQ_sag_up_more_against_nem_wiltQ1 <- GenTable(tGOdata1, fisher=results.fisher1, orderBy="fisher", topNodes=50)


# for up regulated in nem use Q5, for more up in sag use Q1, more up in nem use Q8)
allGenes_numeric <- ifelse((right_universe_sag_nem$quadrant == "Q5") & right_universe_sag_nem$color == "red", 0, 1)
names(allGenes_numeric) <- right_universe_sag_nem$At

# Remove NA values if necessary
allGenes_numeric <- allGenes_numeric[!is.na(allGenes_numeric)]

head(allGenes_numeric)

# Create topGO4 data object
tGOdata5 <- new("topGOdata",
                description = "Enrichment Analysis for Q5",
                ontology = "BP",  
                allGenes = allGenes_numeric,
                geneSel = function(x) x == 0,  # This marks genes of interest (red in Q1) as TRUE
                nodeSize = 10,  # Minimum number of genes for a GO term to be considered
                mapping = "org.At.tair.db",  # database for Arabidopsis thaliana
                annot = annFUN.org)

# enrichment test KS (KS test targets specific and significant enrichment)
results.fisher5 <- runTest(tGOdata5, algorithm="elim", statistic="fisher")



# for up regulated in nem use Q5, for more up in sag use Q1, more up in nem use Q8)
allGenes_numeric <- ifelse((right_universe_sag_nem$quadrant == "Q8") & right_universe_sag_nem$color == "red", 0, 1)
names(allGenes_numeric) <- right_universe_sag_nem$At

# Remove NA values if necessary
allGenes_numeric <- allGenes_numeric[!is.na(allGenes_numeric)]

head(allGenes_numeric)

# Create topGO7 data object
tGOdata8 <- new("topGOdata",
                description = "Enrichment Analysis for Q8",
                ontology = "BP",  
                allGenes = allGenes_numeric,
                geneSel = function(x) x == 0,  # This marks genes of interest (red in Q1) as TRUE
                nodeSize = 10,  # Minimum number of genes for a GO term to be considered
                mapping = "org.At.tair.db",  # database for Arabidopsis thaliana
                annot = annFUN.org)

# enrichment test KS (KS / fisher test targets specific and significant enrichment)
results.fisher8 <- runTest(tGOdata8, algorithm="elim", statistic="fisher")

goEnrichmentQ_sag_up_more_against_nem_wiltQ1 <- GenTable(tGOdata1, fisher=results.fisher1, orderBy="fisher", topNodes=50)
goEnrichmentQ_sag_down_more_against_nem_wiltQ4 <- GenTable(tGOdata4, fisher=results.fisher4, orderBy="fisher", topNodes=50)
goEnrichmentQ_sag_up_and_nem_down_wiltQ2 <- GenTable(tGOdata2, fisher=results.fisher2, orderBy="fisher", topNodes=50)
goEnrichmentQ_nem_up_more_against_sag_wiltQ8 <- GenTable(tGOdata8, fisher=results.fisher8, orderBy="fisher", topNodes=50)
goEnrichmentQ_nem_up_and_sag_down_wiltQ5 <- GenTable(tGOdata5, fisher=results.fisher5, orderBy="fisher", topNodes=50)
goEnrichmentQ_nem_down_more_against_sag_wiltQ7 <- GenTable(tGOdata7, fisher=results.fisher7, orderBy="fisher", topNodes=50)


#write.csv(goEnrichmentQ_sag_up_more_against_nem_wiltQ1, "goEnrichmentQ_sag_up_more_against_nem_wiltQ1.csv")
#write.csv(goEnrichmentQ_sag_down_more_against_nem_wiltQ4, "goEnrichmentQ_sag_down_more_against_nem_wiltQ4.csv")
#write.csv(goEnrichmentQ_sag_up_and_nem_down_wiltQ2, "goEnrichmentQ_sag_up_and_nem_down_wiltQ2.csv")
#write.csv(goEnrichmentQ_nem_up_more_against_sag_wiltQ8, "goEnrichmentQ_nem_up_more_against_sag_wiltQ8.csv")
#write.csv(goEnrichmentQ_nem_up_and_sag_down_wiltQ5, "goEnrichmentQ_nem_up_and_sag_down_wiltQ5.csv")
#write.csv(goEnrichmentQ_nem_down_more_against_sag_wiltQ7, "goEnrichmentQ_nem_down_more_against_sag_wiltQ7.csv")

###################################################
###################################################
##### A. sag and A. nem (Recovery)

## First we checked shared response in both Arabis species to and recovery, 
# so first we will look for the similar response in both species or recover response mechanisms
# in both species together.
##### left universrse ##### A. nem and A. sag in rcovery

allGenes_numeric <- ifelse((left_universe_sag_nem_surv$log2FC_nem < 0 & left_universe_sag_nem_surv$log2FC_sag < 0 ) & left_universe_sag_nem_surv$color == "green", 0, 1) ## check for common enrichment in green genes.

names(allGenes_numeric) <- left_universe_sag_nem_surv$At

# Remove NA values if necessary
allGenes_numeric <- allGenes_numeric[!is.na(allGenes_numeric)]

head(allGenes_numeric)

# Create topGO data object
tGOdata_sag_nem_left_surv_ctrl <- new("topGOdata",
                                      description = "Enrichment Analysis for Q1",
                                      ontology = "BP",  
                                      allGenes = allGenes_numeric,
                                      geneSel = function(x) x == 0,  # This marks genes of interest (green in Q1) as TRUE
                                      nodeSize = 5,  # Minimum number of genes for a GO term to be considered
                                      mapping = "org.At.tair.db",  # database for Arabidopsis thaliana
                                      annot = annFUN.org)

# enrichment test KS (KS test targets specific and significant enrichment)
results.fisher_sag_nem_left_surv_ctrl <- runTest(tGOdata_sag_nem_left_surv_ctrl, algorithm="elim", statistic="fisher")

goEnrichmentQ_shared_sag_nem_left_surv_ctrl <- GenTable(tGOdata_sag_nem_left_surv_ctrl, fisher=results.fisher_sag_nem_left_surv_ctrl, orderBy="fisher", topNodes=50)
#write.csv(goEnrichmentQ_shared_sag_nem_left_surv_ctrl, "goEnrichmentQ_shared_sag_nem_left_surv_ctrl.csv")



##### right universrse ##### A. nem and A. sag in rcovery

allGenes_numeric <- ifelse((right_universe_sag_nem_surv$log2FC_nem > 0 & right_universe_sag_nem_surv$log2FC_sag > 0 ) & right_universe_sag_nem_surv$color == "green", 0, 1) ## check for common enrichment in green genes.

names(allGenes_numeric) <- right_universe_sag_nem_surv$At

# Remove NA values if necessary
allGenes_numeric <- allGenes_numeric[!is.na(allGenes_numeric)]

head(allGenes_numeric)

# Create topGO data object
tGOdata_sag_nem_right_surv_ctrl <- new("topGOdata",
                                       description = "Enrichment Analysis for Q1",
                                       ontology = "BP",  
                                       allGenes = allGenes_numeric,
                                       geneSel = function(x) x == 0,  # This marks genes of interest (green in Q1) as TRUE
                                       nodeSize = 5,  # Minimum number of genes for a GO term to be considered
                                       mapping = "org.At.tair.db",  # database for Arabidopsis thaliana
                                       annot = annFUN.org)

# enrichment test KS (KS test targets specific and significant enrichment)
results.fisher_sag_nem_right_surv_ctrl <- runTest(tGOdata_sag_nem_right_surv_ctrl, algorithm="elim", statistic="fisher")

goEnrichmentQ_shared_sag_nem_right_surv_ctrl <- GenTable(tGOdata_sag_nem_right_surv_ctrl, fisher=results.fisher_sag_nem_right_surv_ctrl, orderBy="fisher", topNodes=50)
#write.csv(goEnrichmentQ_shared_sag_nem_right_surv_ctrl, "goEnrichmentQ_shared_sag_nem_right_surv_ctrl.csv")


##### left universrse
# for down regulated in sag use Q7, for more down in sag use Q4, in nem use Q4, more down in nem use Q7)
allGenes_numeric <- ifelse((left_universe_sag_nem_surv$quadrant == "Q7") & left_universe_sag_nem_surv$color == "red", 0, 1)
names(allGenes_numeric) <- left_universe_sag_nem_surv$At

# Remove NA values if necessary
allGenes_numeric <- allGenes_numeric[!is.na(allGenes_numeric)]

head(allGenes_numeric)

# Create topGO data object
tGOdata_rec7 <- new("topGOdata",
               description = "Enrichment Analysis for Q1",
               ontology = "BP",  
               allGenes = allGenes_numeric,
               geneSel = function(x) x == 0,  # This marks genes of interest (red in Q1) as TRUE
               nodeSize = 15,  # Minimum number of genes for a GO term to be considered
               mapping = "org.At.tair.db",  # database for Arabidopsis thaliana
               annot = annFUN.org)

# enrichment test KS (KS test targets specific and significant enrichment)
results.fisher7 <- runTest(tGOdata_rec7, algorithm="elim", statistic="fisher")




# for down regulated in sag use Q7, for more down in sag use Q4, in nem use Q4, more down in nem use Q7)
allGenes_numeric <- ifelse((left_universe_sag_nem_surv$quadrant == "Q4") & left_universe_sag_nem_surv$color == "red", 0, 1)
names(allGenes_numeric) <- left_universe_sag_nem_surv$At

# Remove NA values if necessary
allGenes_numeric <- allGenes_numeric[!is.na(allGenes_numeric)]

head(allGenes_numeric)

# Create topGO data object
tGOdata_rec4 <- new("topGOdata",
                description = "Enrichment Analysis for Q1",
                ontology = "BP",  
                allGenes = allGenes_numeric,
                geneSel = function(x) x == 0,  # This marks genes of interest (red in Q1) as TRUE
                nodeSize = 15,  # Minimum number of genes for a GO term to be considered
                mapping = "org.At.tair.db",  # database for Arabidopsis thaliana
                annot = annFUN.org)

# enrichment test KS (KS test targets specific and significant enrichment)
results.fisher4 <- runTest(tGOdata_rec4, algorithm="elim", statistic="fisher")




# for down regulated in sag use Q7, for more down in sag use Q4, in nem use Q4, more down in nem use Q7)
allGenes_numeric <- ifelse((left_universe_sag_nem_surv$quadrant == "Q2") & left_universe_sag_nem_surv$color == "red", 0, 1)
names(allGenes_numeric) <- left_universe_sag_nem_surv$At

# Remove NA values if necessary
allGenes_numeric <- allGenes_numeric[!is.na(allGenes_numeric)]

head(allGenes_numeric)

# Create topGO data object
tGOdata_rec2 <- new("topGOdata",
                description = "Enrichment Analysis for Q1",
                ontology = "BP",  
                allGenes = allGenes_numeric,
                geneSel = function(x) x == 0,  # This marks genes of interest (red in Q1) as TRUE
                nodeSize = 15,  # Minimum number of genes for a GO term to be considered
                mapping = "org.At.tair.db",  # database for Arabidopsis thaliana
                annot = annFUN.org)

# enrichment test KS (KS test targets specific and significant enrichment)
results.fisher2 <- runTest(tGOdata_rec2, algorithm="elim", statistic="fisher")


##### right universrse

# for down regulated in sag use Q7, for more down in sag use Q4, in nem use Q4, more down in nem use Q7)
allGenes_numeric <- ifelse((right_universe_sag_nem_surv$quadrant == "Q1") & right_universe_sag_nem_surv$color == "red", 0, 1)
names(allGenes_numeric) <- right_universe_sag_nem_surv$At

# Remove NA values if necessary 
allGenes_numeric <- allGenes_numeric[!is.na(allGenes_numeric)]

head(allGenes_numeric)

# Create topGO data object
tGOdata_rec1 <- new("topGOdata",
                description = "Enrichment Analysis for Q1",
                ontology = "BP",  
                allGenes = allGenes_numeric,
                geneSel = function(x) x == 0,  # This marks genes of interest (red in Q1) as TRUE
                nodeSize = 15,  # Minimum number of genes for a GO term to be considered
                mapping = "org.At.tair.db",  # database for Arabidopsis thaliana
                annot = annFUN.org)

# enrichment test KS (KS test targets specific and significant enrichment)
results.fisher1 <- runTest(tGOdata_rec1, algorithm="elim", statistic="fisher")

# for down regulated in sag use Q7, for more down in sag use Q4, in nem use Q4, more down in nem use Q7)
allGenes_numeric <- ifelse((right_universe_sag_nem_surv$quadrant == "Q5") & right_universe_sag_nem_surv$color == "red", 0, 1)
names(allGenes_numeric) <- right_universe_sag_nem_surv$At

# Remove NA values if necessary
allGenes_numeric <- allGenes_numeric[!is.na(allGenes_numeric)]

head(allGenes_numeric)

# Create topGO data object
tGOdata_rec5 <- new("topGOdata",
                description = "Enrichment Analysis for Q1",
                ontology = "BP",  
                allGenes = allGenes_numeric,
                geneSel = function(x) x == 0,  # This marks genes of interest (red in Q1) as TRUE
                nodeSize = 15,  # Minimum number of genes for a GO term to be considered
                mapping = "org.At.tair.db",  # database for Arabidopsis thaliana
                annot = annFUN.org)

# enrichment test KS (KS test targets specific and significant enrichment)
results.fisher5 <- runTest(tGOdata_rec5, algorithm="elim", statistic="fisher")

# for down regulated in sag use Q7, for more down in sag use Q4, in nem use Q4, more down in nem use Q7)
allGenes_numeric <- ifelse((right_universe_sag_nem_surv$quadrant == "Q8") & right_universe_sag_nem_surv$color == "red", 0, 1)
names(allGenes_numeric) <- right_universe_sag_nem_surv$At

# Remove NA values if necessary
allGenes_numeric <- allGenes_numeric[!is.na(allGenes_numeric)]

head(allGenes_numeric)

# Create topGO data object
tGOdata_rec8 <- new("topGOdata",
                description = "Enrichment Analysis for Q1",
                ontology = "BP",  
                allGenes = allGenes_numeric,
                geneSel = function(x) x == 0,  # This marks genes of interest (red in Q1) as TRUE
                nodeSize = 15,  # Minimum number of genes for a GO term to be considered
                mapping = "org.At.tair.db",  # database for Arabidopsis thaliana
                annot = annFUN.org)

# enrichment test KS (KS test targets specific and significant enrichment)
results.fisher8 <- runTest(tGOdata_rec8, algorithm="elim", statistic="fisher")

# Generate table of enriched GO terms
goEnrichmentQ_sag_up_more_against_nem_rec_Q1 <- GenTable(tGOdata_rec1, fisher=results.fisher1, orderBy="fisher", topNodes=50)
goEnrichmentQ_sag_up_and_nem_down_rec_Q2 <- GenTable(tGOdata_rec2, fisher=results.fisher2, orderBy="fisher", topNodes=50)
goEnrichmentQ_sag_down_more_against_nem_rec_Q4 <- GenTable(tGOdata_rec4, fisher=results.fisher4, orderBy="fisher", topNodes=50)
goEnrichmentQ_nem_up_and_sag_down_rec_Q5 <- GenTable(tGOdata_rec5, fisher=results.fisher5, orderBy="fisher", topNodes=50)
goEnrichmentQ_nem_down_more_against_sag_rec_Q7 <- GenTable(tGOdata_rec7, fisher=results.fisher7, orderBy="fisher", topNodes=50)
goEnrichmentQ_nem_up_more_against_sag_rec_Q8 <- GenTable(tGOdata_rec8, fisher=results.fisher8, orderBy="fisher", topNodes=50)


#write.csv(goEnrichmentQ_sag_up_more_against_nem_rec_Q1, "goEnrichmentQ_sag_up_more_against_nem_recQ1.csv")
#write.csv(goEnrichmentQ_sag_down_more_against_nem_rec_Q4, "goEnrichmentQ_sag_down_more_against_nem_recQ4.csv")
#write.csv(goEnrichmentQ_sag_up_and_nem_down_rec_Q2, "goEnrichmentQ_sag_up_and_nem_down_rec_Q2.csv")
#write.csv(goEnrichmentQ_nem_up_and_sag_down_rec_Q5, "goEnrichmentQ_nem_up_and_sag_down_rec_Q5.csv")
#write.csv(goEnrichmentQ_nem_down_more_against_sag_rec_Q7, "goEnrichmentQ_nem_down_more_against_sag_rec_Q7.csv")
#write.csv(goEnrichmentQ_nem_up_more_against_sag_rec_Q8, "goEnrichmentQ_nem_up_more_against_sag_rec_Q8.csv")



###############################################################
###############################################################
#########################################################-------------
########## GO Recovery and Stress
#Import the orthologues CSV file
# ========================== STEP 1: Load Data ==========================
orthologues <- read.csv("/Users/Shared/Files From d.localized/PhD-Uni-Koeln_2021-2024/PhD_work/seeds data/2nd Experiment 2022/transcriptome/Arabis_RNA_raw_data/90-774066047/anaylsis_drought_mRNA_both_species_with_new_two_genomes/analysis_both_species_with_nem_genomes/orthologues_cleaned.csv", header = TRUE)
new_dataframe1 <- data.frame(orthologues)

# Merge expression data with orthologue mapping
combined_species_df_sag_nem_surv_wilt <- combined_species_df_sag_nem_surv_wilt %>%
  left_join(new_dataframe1, by = c("gene_id" = "arabis_cleaned"))

write.csv(combined_species_df_sag_nem_surv_wilt, "combined_species_df_sag_nem_surv_wilt_with_ortho.csv")


# Remove genes without an orthologue
combined_species_df_sag_nem_surv_wilt <- combined_species_df_sag_nem_surv_wilt %>%
  filter(!is.na(At))


# Merge expression data with orthologue mapping
#combined_species_df_sag_nem_surv_wilt <- combined_species_df_sag_nem_surv_wilt %>%
#  left_join(new_dataframe1, by = c("gene_id" = "arabis_cleaned"))

# Fix column names like At.x, At.x.x, At.y etc.
colnames(combined_species_df_sag_nem_surv_wilt) <- gsub("\\.x+|\\.y+", "", colnames(combined_species_df_sag_nem_surv_wilt))

combined_species_df_sag_nem_surv_wilt <- combined_species_df_sag_nem_surv_wilt[ , !grepl("^At\\.\\d+$", colnames(combined_species_df_sag_nem_surv_wilt))]

colnames(combined_species_df_sag_nem_surv_wilt) <- make.unique(colnames(combined_species_df_sag_nem_surv_wilt))

# Remove genes without an orthologue
#combined_species_df_sag_nem_surv_wilt <- combined_species_df_sag_nem_surv_wilt %>%
#  filter(!is.na(At))

# ========================== STEP 2: Define Universes ==========================
### split on diagonal
##### A. sag and A. nem (stress)
# Genes ABOVE the diagonal
upper_universe_sag_nem_surv_wilt <- combined_species_df_sag_nem_surv_wilt %>%
  filter(log2FC_sag > log2FC_nem) 

# Genes BELOW the diagonal
lower_universe_sag_nem_surv_wilt <- combined_species_df_sag_nem_surv_wilt %>%
  filter(log2FC_sag < log2FC_nem)

##### A. sag and A. nem (recovery)
# Genes ABOVE the diagonal
upper_universe_sag_nem_surv_wilt <- combined_species_df_sag_nem_surv_wilt %>%
  filter(log2FC_sag > log2FC_nem) 

# Genes BELOW the diagonal
lower_universe_sag_nem_surv_wilt <- combined_species_df_sag_nem_surv_wilt %>%
  filter(log2FC_sag < log2FC_nem)

##### split on vertical line
##### A. sag and A. nem (stress)
# Genes Left the vline
left_universe_sag_nem_surv_wilt <- combined_species_df_sag_nem_surv_wilt %>%
  filter(log2FC_nem < 0) 

# Genes right the vline
right_universe_sag_nem_surv_wilt <- combined_species_df_sag_nem_surv_wilt %>%
  filter(log2FC_nem > 0)

##### A. sag and A. nem (recovery)
# Genes ABOVE the diagonal
left_universe_sag_nem_surv_wilt <- combined_species_df_sag_nem_surv_wilt %>%
  filter(log2FC_nem < 0) 

# Genes BELOW the diagonal
right_universe_sag_nem_surv_wilt <- combined_species_df_sag_nem_surv_wilt %>%
  filter(log2FC_nem > 0)


write.csv(left_universe_sag_nem_surv_wilt, "left_universe_sag_nem_recovery_wilt.csv")

write.csv(right_universe_sag_nem_surv_wilt, "right_universe_sag_nem_recovery_wilt.csv")

# ========================== STEP 3: Perform GO Enrichments ==========================
setwd("/Users/Shared/Files From d.localized/PhD-Uni-Koeln_2021-2024/PhD_work/seeds data/2nd Experiment 2022/transcriptome/Arabis_RNA_raw_data/90-774066047/anaylsis_drought_mRNA_both_species_with_new_two_genomes/analysis_both_species_with_nem_genomes/script_for_juliette/new_go_output_from_deseqdrought_all/New_GO_analysis_6thMay2025//")

## First we checked shared response in both Arabis species to drought and recovery, 
# so first we will look for the similar response in both species or stress response mechanisms
# in both species together.
##### left universrse ##### A. nem and A. sag in drought
allGenes_numeric <- ifelse((left_universe_sag_nem_surv_wilt$log2FC_nem < 0 & left_universe_sag_nem_surv_wilt$log2FC_sag < 0 ) & left_universe_sag_nem_surv_wilt$color == "green", 0, 1) ## check for common enrichment in green genes.

names(allGenes_numeric) <- left_universe_sag_nem_surv_wilt$At

# Remove NA values if necessary
allGenes_numeric <- allGenes_numeric[!is.na(allGenes_numeric)]

head(allGenes_numeric)

# Create topGO data object
tGOdata_sag_nem_left_surv_wilt <- new("topGOdata",
                                      description = "Enrichment Analysis for Q1",
                                      ontology = "BP",  
                                      allGenes = allGenes_numeric,
                                      geneSel = function(x) x == 0,  # This marks genes of interest (green in Q1) as TRUE
                                      nodeSize = 5,  # Minimum number of genes for a GO term to be considered
                                      mapping = "org.At.tair.db",  # database for Arabidopsis thaliana
                                      annot = annFUN.org)

# enrichment test KS (KS test targets specific and significant enrichment)
results.fisher_sag_nem_left_surv_wilt <- runTest(tGOdata_sag_nem_left_surv_wilt, algorithm="elim", statistic="fisher")

goEnrichmentQ_shared_sag_nem_left_surv_wilt <- GenTable(tGOdata_sag_nem_left_surv_wilt, fisher=results.fisher_sag_nem_left_surv_wilt, orderBy="fisher", topNodes=50)

#write.csv(goEnrichmentQ_shared_sag_nem_left_wilt_ctrl, "goEnrichmentQ_shared_sag_nem_left_wilt_ctrl.csv")



##### right universrse ##### A. nem and A. sag in drought

allGenes_numeric <- ifelse((right_universe_sag_nem_surv_wilt$log2FC_nem > 0 & right_universe_sag_nem_surv_wilt$log2FC_sag > 0 ) & right_universe_sag_nem$color == "green", 0, 1) ## check for common enrichment in green genes.

names(allGenes_numeric) <- right_universe_sag_nem_surv_wilt$At

# Remove NA values if necessary
allGenes_numeric <- allGenes_numeric[!is.na(allGenes_numeric)]

head(allGenes_numeric)

# Create topGO data object
tGOdata_sag_nem_right_wilt_ctrl <- new("topGOdata",
                                       description = "Enrichment Analysis for Q1",
                                       ontology = "BP",  
                                       allGenes = allGenes_numeric,
                                       geneSel = function(x) x == 0,  # This marks genes of interest (green in Q1) as TRUE
                                       nodeSize = 5,  # Minimum number of genes for a GO term to be considered
                                       mapping = "org.At.tair.db",  # database for Arabidopsis thaliana
                                       annot = annFUN.org)

# enrichment test KS (KS test targets specific and significant enrichment)
results.fisher_sag_nem_right_surv_wilt <- runTest(tGOdata_sag_nem_right_surv_wilt, algorithm="elim", statistic="fisher")

goEnrichmentQ_shared_sag_nem_right_surv_wilt <- GenTable(tGOdata_sag_nem_right_surv_wilt, fisher=results.fisher_sag_nem_right_wilt_ctrl, orderBy="fisher", topNodes=50)
#write.csv(goEnrichmentQ_shared_sag_nem_right_wilt_ctrl, "goEnrichmentQ_shared_sag_nem_right_wilt_ctrl.csv")



##### A. sag and A. nem (Recov and stress)
##### left universrse
# for down regulated in sag use Q7, for more down in sag use Q4, in nem use Q4, more down in nem use Q7)
allGenes_numeric <- ifelse((left_universe_sag_nem_surv_wilt$quadrant == "Q2") & left_universe_sag_nem_surv_wilt$color == "red", 0, 1)
names(allGenes_numeric) <- left_universe_sag_nem_surv_wilt$At

# Remove NA values if necessary
allGenes_numeric <- allGenes_numeric[!is.na(allGenes_numeric)]

head(allGenes_numeric)

# Create topGO2 data object
tGOdata2 <- new("topGOdata",
                description = "Enrichment Analysis for Q2",
                ontology = "BP",  
                allGenes = allGenes_numeric,
                geneSel = function(x) x == 0,  # This marks genes of interest (red in Q1) as TRUE
                nodeSize = 10,  # Minimum number of genes for a GO term to be considered
                mapping = "org.At.tair.db",  # database for Arabidopsis thaliana
                annot = annFUN.org)

# enrichment test KS (KS test targets specific and significant enrichment)
results.fisher2 <- runTest(tGOdata2, algorithm="elim", statistic="fisher")
goEnrichmentQ_sag_up_against_nem_surv_wiltQ2 <- GenTable(tGOdata2, KS=results.fisher2, orderBy="KS", topNodes=50)


# for down regulated in sag use Q7, for more down in sag use Q4, in nem use Q4, more down in nem use Q7)
allGenes_numeric <- ifelse((left_universe_sag_nem_surv_wilt$quadrant == "Q4") & left_universe_sag_nem_surv_wilt$color == "red", 0, 1)
names(allGenes_numeric) <- left_universe_sag_nem_surv_wilt$At

# Remove NA values if necessary
allGenes_numeric <- allGenes_numeric[!is.na(allGenes_numeric)]

head(allGenes_numeric)

# Create topGO4 data object
tGOdata4 <- new("topGOdata",
                description = "Enrichment Analysis for Q4",
                ontology = "BP",  
                allGenes = allGenes_numeric,
                geneSel = function(x) x == 0,  # This marks genes of interest (red in Q1) as TRUE
                nodeSize = 10,  # Minimum number of genes for a GO term to be considered
                mapping = "org.At.tair.db",  # database for Arabidopsis thaliana
                annot = annFUN.org)

# enrichment test KS (KS test targets specific and significant enrichment)
results.fisher4 <- runTest(tGOdata4, algorithm="elim", statistic="fisher")
goEnrichmentQ_sag_down_more_against_nem_surv_wiltQ4 <- GenTable(tGOdata4, KS=results.fisher4, orderBy="KS", topNodes=50)


# for down regulated in sag use Q7, for more down in sag use Q4, in nem use Q4, more down in nem use Q7)
allGenes_numeric <- ifelse((left_universe_sag_nem_surv_wilt$quadrant == "Q7") & left_universe_sag_nem_surv_wilt$color == "red", 0, 1)
names(allGenes_numeric) <- left_universe_sag_nem_surv_wilt$At

# Remove NA values if necessary
allGenes_numeric <- allGenes_numeric[!is.na(allGenes_numeric)]

head(allGenes_numeric)

# Create topGO7 data object
tGOdata7 <- new("topGOdata",
                description = "Enrichment Analysis for Q7",
                ontology = "BP",  
                allGenes = allGenes_numeric,
                geneSel = function(x) x == 0,  # This marks genes of interest (red in Q1) as TRUE
                nodeSize = 10,  # Minimum number of genes for a GO term to be considered
                mapping = "org.At.tair.db",  # database for Arabidopsis thaliana
                annot = annFUN.org)

# enrichment test KS (KS test targets specific and significant enrichment)
results.fisher7 <- runTest(tGOdata7, algorithm="elim", statistic="fisher")

goEnrichmentQ_sag_down_more_against_nem_surv_wiltQ7 <- GenTable(tGOdata7, KS=results.fisher7, orderBy="KS", topNodes=50)

###### right universe

# for up regulated in nem use Q5, for more up in sag use Q1, more up in nem use Q8)
allGenes_numeric <- ifelse((right_universe_sag_nem_surv_wilt$quadrant == "Q1") & right_universe_sag_nem_surv_wilt$color == "red", 0, 1)
names(allGenes_numeric) <- right_universe_sag_nem_surv_wilt$At

# Remove NA values if necessary
allGenes_numeric <- allGenes_numeric[!is.na(allGenes_numeric)]

head(allGenes_numeric)

# Create topGO2 data object
tGOdata1 <- new("topGOdata",
                description = "Enrichment Analysis for Q1",
                ontology = "BP",  
                allGenes = allGenes_numeric,
                geneSel = function(x) x == 0,  # This marks genes of interest (red in Q1) as TRUE
                nodeSize = 10,  # Minimum number of genes for a GO term to be considered
                mapping = "org.At.tair.db",  # database for Arabidopsis thaliana
                annot = annFUN.org)

# enrichment test KS (KS test targets specific and significant enrichment)
results.fisher1 <- runTest(tGOdata1, algorithm="elim", statistic="fisher")

goEnrichmentQ_sag_up_more_against_nem_surv_wiltQ1 <- GenTable(tGOdata1, fisher=results.fisher1, orderBy="fisher", topNodes=50)


# for up regulated in nem use Q5, for more up in sag use Q1, more up in nem use Q8)
allGenes_numeric <- ifelse((right_universe_sag_nem_surv_wilt$quadrant == "Q5") & right_universe_sag_nem_surv_wilt$color == "red", 0, 1)
names(allGenes_numeric) <- right_universe_sag_nem_surv_wilt$At

# Remove NA values if necessary
allGenes_numeric <- allGenes_numeric[!is.na(allGenes_numeric)]

head(allGenes_numeric)

# Create topGO4 data object
tGOdata5 <- new("topGOdata",
                description = "Enrichment Analysis for Q5",
                ontology = "BP",  
                allGenes = allGenes_numeric,
                geneSel = function(x) x == 0,  # This marks genes of interest (red in Q1) as TRUE
                nodeSize = 10,  # Minimum number of genes for a GO term to be considered
                mapping = "org.At.tair.db",  # database for Arabidopsis thaliana
                annot = annFUN.org)

# enrichment test KS (KS test targets specific and significant enrichment)
results.fisher5 <- runTest(tGOdata5, algorithm="elim", statistic="fisher")

goEnrichmentQ_sag_up_more_against_nem_surv_wiltQ5 <- GenTable(tGOdata5, fisher=results.fisher5, orderBy="fisher", topNodes=50)



# for up regulated in nem use Q5, for more up in sag use Q1, more up in nem use Q8)
allGenes_numeric <- ifelse((right_universe_sag_nem_surv_wilt$quadrant == "Q8") & right_universe_sag_nem_surv_wilt$color == "red", 0, 1)
names(allGenes_numeric) <- right_universe_sag_nem_surv_wilt$At

# Remove NA values if necessary
allGenes_numeric <- allGenes_numeric[!is.na(allGenes_numeric)]

head(allGenes_numeric)

# Create topGO7 data object
tGOdata8 <- new("topGOdata",
                description = "Enrichment Analysis for Q8",
                ontology = "BP",  
                allGenes = allGenes_numeric,
                geneSel = function(x) x == 0,  # This marks genes of interest (red in Q1) as TRUE
                nodeSize = 10,  # Minimum number of genes for a GO term to be considered
                mapping = "org.At.tair.db",  # database for Arabidopsis thaliana
                annot = annFUN.org)

# enrichment test KS (KS / fisher test targets specific and significant enrichment)
results.fisher8 <- runTest(tGOdata8, algorithm="elim", statistic="fisher")
goEnrichmentQ_sag_up_more_against_nem_surv_wiltQ8 <- GenTable(tGOdata8, fisher=results.fisher8, orderBy="fisher", topNodes=50)




goEnrichmentQ_sag_up_more_against_nem_surv_wiltQ1 <- GenTable(tGOdata1, fisher=results.fisher1, orderBy="fisher", topNodes=50)
goEnrichmentQ_sag_down_more_against_nem_surv_wiltQ4 <- GenTable(tGOdata4, fisher=results.fisher4, orderBy="fisher", topNodes=50)
goEnrichmentQ_sag_up_and_nem_down_surv_wiltQ2 <- GenTable(tGOdata2, fisher=results.fisher2, orderBy="fisher", topNodes=50)
goEnrichmentQ_nem_up_more_against_surv_sag_wiltQ8 <- GenTable(tGOdata8, fisher=results.fisher8, orderBy="fisher", topNodes=50)
goEnrichmentQ_nem_up_and_sag_down_surv_wiltQ5 <- GenTable(tGOdata5, fisher=results.fisher5, orderBy="fisher", topNodes=50)
goEnrichmentQ_nem_down_more_against_sag_surv_wiltQ7 <- GenTable(tGOdata7, fisher=results.fisher7, orderBy="fisher", topNodes=50)


#write.csv(goEnrichmentQ_sag_up_more_against_nem_surv_wiltQ1, "goEnrichmentQ_sag_up_more_against_nem_surv_wiltQ1.csv")
#write.csv(goEnrichmentQ_sag_down_more_against_nem_surv_wiltQ4, "goEnrichmentQ_sag_down_more_against_nem_surv_wiltQ4.csv")
#write.csv(goEnrichmentQ_sag_up_and_nem_down_surv_wiltQ2, "goEnrichmentQ_sag_up_and_nem_down_surv_wiltQ2.csv")
#write.csv(goEnrichmentQ_nem_up_more_against_surv_sag_wiltQ8, "goEnrichmentQ_nem_up_more_against_surv_sag_wiltQ8.csv")
#write.csv(goEnrichmentQ_nem_up_and_sag_down_surv_wiltQ5, "goEnrichmentQ_nem_up_and_sag_down_surv_wiltQ5.csv")
#write.csv(goEnrichmentQ_nem_down_more_against_sag_surv_wiltQ7, "goEnrichmentQ_nem_down_more_against_sag_surv_wiltQ7.csv")


################################################################
################################################################
################################################################
############################  gene shared by GO terms #########

library(igraph)
library(AnnotationDbi)
library(org.At.tair.db)
library(GO.db)
library(GO.db)
library(dplyr)
library(purrr)
library(tibble)
######################################

# Example GO enrichment result Q1 / stress
sig_go_wiltQ1_sag_up_more <- goEnrichmentQ_sag_up_more_against_nem_wiltQ1$GO.ID

# Get significant genes per GO term
# Corrected function
getSigGenesFromGO <- function(GOterm, GOdata) {
  all_genes <- genesInTerm(GOdata)[[GOterm]]
  sig_genes <- all_genes[all_genes %in% sigGenes(GOdata)]
  return(sig_genes)
}

goEnrichmentQ_sag_up_more_against_nem_wiltQ1_genes <- lapply(sig_go_wiltQ1_sag_up_more, getSigGenesFromGO, GOdata = tGOdata1)
names(goEnrichmentQ_sag_up_more_against_nem_wiltQ1_genes) <- sig_go_wiltQ1_sag_up_more

# Flatten to a data frame
goEnrichmentQ_sag_up_more_against_nem_wiltQ1_genes_df <- stack(goEnrichmentQ_sag_up_more_against_nem_wiltQ1_genes)
colnames(goEnrichmentQ_sag_up_more_against_nem_wiltQ1_genes_df) <- c("GeneID", "GO.ID")

# Save
write.csv(goEnrichmentQ_sag_up_more_against_nem_wiltQ1_genes_df, "goEnrichmentQ_sag_up_more_against_nem_wiltQ1_genes.csv", row.names = FALSE)


# Example GO enrichment result Q2 / stress
sig_go_wiltQ2 <- goEnrichmentQ_sag_up_and_nem_down_wiltQ2$GO.ID

# Get significant genes per GO term
# Corwiltted function
getSigGenesFromGO <- function(GOterm, GOdata) {
  all_genes <- genesInTerm(GOdata)[[GOterm]]
  sig_genes <- all_genes[all_genes %in% sigGenes(GOdata)]
  return(sig_genes)
}

goEnrichmentQ_sag_up_and_nem_down_wiltQ2_genes <- lapply(sig_go_wiltQ2, getSigGenesFromGO, GOdata = tGOdata2)
names(goEnrichmentQ_sag_up_and_nem_down_wiltQ2_genes) <- sig_go_wiltQ2

# Flatten to a data frame
goEnrichmentQ_sag_up_and_nem_down_wiltQ2_genes_df <- stack(goEnrichmentQ_sag_up_and_nem_down_wiltQ2_genes)
colnames(goEnrichmentQ_sag_up_and_nem_down_wiltQ2_genes_df) <- c("GeneID", "GO.ID")

# Save
write.csv(goEnrichmentQ_sag_up_and_nem_down_wiltQ2_genes_df, "goEnrichmentQ_sag_up_and_nem_down_wiltQ2_genes.csv", row.names = FALSE)


# Example GO enrichment result Q4 / stress
sig_go_wiltQ4 <- goEnrichmentQ_sag_down_more_against_nem_wiltQ4$GO.ID

# Get significant genes per GO term
# Corwiltted function
getSigGenesFromGO <- function(GOterm, GOdata) {
  all_genes <- genesInTerm(GOdata)[[GOterm]]
  sig_genes <- all_genes[all_genes %in% sigGenes(GOdata)]
  return(sig_genes)
}

goEnrichmentQ_sag_down_more_against_nem_wiltQ4_genes <- lapply(sig_go_wiltQ4, getSigGenesFromGO, GOdata = tGOdata4)
names(goEnrichmentQ_sag_down_more_against_nem_wiltQ4_genes) <- sig_go_wiltQ4

# Flatten to a data frame
goEnrichmentQ_sag_down_more_against_nem_wiltQ4_genes_df <- stack(goEnrichmentQ_sag_down_more_against_nem_wiltQ4_genes)
colnames(goEnrichmentQ_sag_down_more_against_nem_wiltQ4_genes_df) <- c("GeneID", "GO.ID")

# Save
write.csv(goEnrichmentQ_sag_down_more_against_nem_wiltQ4_genes_df, "goEnrichmentQ_sag_down_more_against_nem_wiltQ4_genes.csv", row.names = FALSE)

# Example GO enrichment result Q5 / stress
sig_go_wiltQ5 <- goEnrichmentQ_nem_up_and_sag_down_wiltQ5$GO.ID

# Get significant genes per GO term
# Corwiltted function
getSigGenesFromGO <- function(GOterm, GOdata) {
  all_genes <- genesInTerm(GOdata)[[GOterm]]
  sig_genes <- all_genes[all_genes %in% sigGenes(GOdata)]
  return(sig_genes)
}

goEnrichmentQ_nem_up_and_sag_down_wiltQ5_genes <- lapply(sig_go_wiltQ5, getSigGenesFromGO, GOdata = tGOdata5)
names(goEnrichmentQ_nem_up_and_sag_down_wiltQ5_genes) <- sig_go_wiltQ5

# Flatten to a data frame
goEnrichmentQ_nem_up_and_sag_down_wiltQ5_genes_df <- stack(goEnrichmentQ_nem_up_and_sag_down_wiltQ5_genes)
colnames(goEnrichmentQ_nem_up_and_sag_down_wiltQ5_genes_df) <- c("GeneID", "GO.ID")

# Save
write.csv(goEnrichmentQ_nem_up_and_sag_down_wiltQ5_genes_df, "goEnrichmentQ_nem_up_and_sag_down_wiltQ5_genes.csv", row.names = FALSE)


# Example GO enrichment result Q7 / stress
sig_go_wiltQ7 <- goEnrichmentQ_nem_down_more_against_sag_wiltQ7$GO.ID

# Get significant genes per GO term
# Corwiltted function
getSigGenesFromGO <- function(GOterm, GOdata) {
  all_genes <- genesInTerm(GOdata)[[GOterm]]
  sig_genes <- all_genes[all_genes %in% sigGenes(GOdata)]
  return(sig_genes)
}

goEnrichmentQ_nem_down_more_against_sag_wiltQ7_genes <- lapply(sig_go_wiltQ7, getSigGenesFromGO, GOdata = tGOdata7)
names(goEnrichmentQ_nem_down_more_against_sag_wiltQ7_genes) <- sig_go_wiltQ7

# Flatten to a data frame
goEnrichmentQ_nem_down_more_against_sag_wiltQ7_genes_df <- stack(goEnrichmentQ_nem_down_more_against_sag_wiltQ7_genes)
colnames(goEnrichmentQ_nem_down_more_against_sag_wiltQ7_genes_df) <- c("GeneID", "GO.ID")

# Save
write.csv(goEnrichmentQ_nem_down_more_against_sag_wiltQ7_genes_df, "goEnrichmentQ_nem_down_more_against_sag_wiltQ7_genes.csv", row.names = FALSE)




# Example GO enrichment result Q8 / stress
sig_go_wiltQ8 <- goEnrichmentQ_nem_up_more_against_sag_wiltQ8$GO.ID

# Get significant genes per GO term
# Corwiltted function
getSigGenesFromGO <- function(GOterm, GOdata) {
  all_genes <- genesInTerm(GOdata)[[GOterm]]
  sig_genes <- all_genes[all_genes %in% sigGenes(GOdata)]
  return(sig_genes)
}

goEnrichmentQ_nem_up_more_against_sag_wiltQ8_genes <- lapply(sig_go_wiltQ8, getSigGenesFromGO, GOdata = tGOdata8)
names(goEnrichmentQ_nem_up_more_against_sag_wiltQ8_genes) <- sig_go_wiltQ8

# Flatten to a data frame
goEnrichmentQ_nem_up_more_against_sag_wiltQ8_genes_df <- stack(goEnrichmentQ_nem_up_more_against_sag_wiltQ8_genes)
colnames(goEnrichmentQ_nem_up_more_against_sag_wiltQ8_genes_df) <- c("GeneID", "GO.ID")

# Save
write.csv(goEnrichmentQ_nem_up_more_against_sag_wiltQ8_genes_df, "goEnrichmentQ_nem_up_more_against_sag_wiltQ8_genes.csv", row.names = FALSE)




##################
######### recovery


# Example GO enrichment result Q1 / recover
sig_go_recQ1 <- goEnrichmentQ_sag_up_more_against_nem_rec_Q1$GO.ID

# Get significant genes per GO term
# Corrected function
getSigGenesFromGO <- function(GOterm, GOdata) {
  all_genes <- genesInTerm(GOdata)[[GOterm]]
  sig_genes <- all_genes[all_genes %in% sigGenes(GOdata)]
  return(sig_genes)
}

goEnrichmentQ_sag_up_more_against_nem_rec_Q1_genes <- lapply(sig_go_recQ1, getSigGenesFromGO, GOdata = tGOdata_rec1)
names(goEnrichmentQ_sag_up_more_against_nem_rec_Q1_genes) <- sig_go_recQ1

# Flatten to a data frame
goEnrichmentQ_sag_up_more_against_nem_rec_Q1_genes_df <- stack(goEnrichmentQ_sag_up_more_against_nem_rec_Q1_genes)
colnames(goEnrichmentQ_sag_up_more_against_nem_rec_Q1_genes_df) <- c("GeneID", "GO.ID")

# Save
write.csv(goEnrichmentQ_sag_up_more_against_nem_rec_Q1_genes_df, "goEnrichmentQ_sag_up_more_against_nem_recQ1_genes.csv", row.names = FALSE)



# Example GO enrichment result Q2 / recover
sig_go_recQ2 <- goEnrichmentQ_sag_up_and_nem_down_rec_Q2$GO.ID

# Get significant genes per GO term
# Corrected function
getSigGenesFromGO <- function(GOterm, GOdata) {
  all_genes <- genesInTerm(GOdata)[[GOterm]]
  sig_genes <- all_genes[all_genes %in% sigGenes(GOdata)]
  return(sig_genes)
}

goEnrichmentQ_sag_up_and_nem_down_rec_Q2_genes <- lapply(sig_go_recQ2, getSigGenesFromGO, GOdata = tGOdata_rec2)
names(goEnrichmentQ_sag_up_and_nem_down_rec_Q2_genes) <- sig_go_recQ2

# Flatten to a data frame
goEnrichmentQ_sag_up_and_nem_down_rec_Q2_genes_df <- stack(goEnrichmentQ_sag_up_and_nem_down_rec_Q2_genes)
colnames(goEnrichmentQ_sag_up_and_nem_down_rec_Q2_genes_df) <- c("GeneID", "GO.ID")

# Save
write.csv(goEnrichmentQ_sag_up_and_nem_down_rec_Q2_genes_df, "goEnrichmentQ_sag_up_and_nem_down_recQ2_genes.csv", row.names = FALSE)


# Example GO enrichment result Q4 / recover
sig_go_recQ4 <- goEnrichmentQ_sag_down_more_against_nem_rec_Q4$GO.ID

# Get significant genes per GO term
# Corrected function
getSigGenesFromGO <- function(GOterm, GOdata) {
  all_genes <- genesInTerm(GOdata)[[GOterm]]
  sig_genes <- all_genes[all_genes %in% sigGenes(GOdata)]
  return(sig_genes)
}

goEnrichmentQ_sag_down_more_against_nem_recQ4_genes <- lapply(sig_go_recQ4, getSigGenesFromGO, GOdata = tGOdata_rec4)
names(goEnrichmentQ_sag_down_more_against_nem_recQ4_genes) <- sig_go_recQ4

# Flatten to a data frame
goEnrichmentQ_sag_down_more_against_nem_recQ4_genes_df <- stack(goEnrichmentQ_sag_down_more_against_nem_recQ4_genes)
colnames(goEnrichmentQ_sag_down_more_against_nem_recQ4_genes_df) <- c("GeneID", "GO.ID")

# Save
write.csv(goEnrichmentQ_sag_down_more_against_nem_recQ4_genes_df, "goEnrichmentQ_sag_down_more_against_nem_recQ4_genes.csv", row.names = FALSE)



# Example GO enrichment result Q5 / recover
sig_go_recQ5 <- goEnrichmentQ_nem_up_and_sag_down_rec_Q5$GO.ID

# Get significant genes per GO term
# Corrected function
getSigGenesFromGO <- function(GOterm, GOdata) {
  all_genes <- genesInTerm(GOdata)[[GOterm]]
  sig_genes <- all_genes[all_genes %in% sigGenes(GOdata)]
  return(sig_genes)
}

goEnrichmentQ_nem_up_and_sag_down_recQ5_genes <- lapply(sig_go_recQ5, getSigGenesFromGO, GOdata = tGOdata_rec5)
names(goEnrichmentQ_nem_up_and_sag_down_recQ5_genes) <- sig_go_recQ5

# Flatten to a data frame
goEnrichmentQ_nem_up_and_sag_down_recQ5_genes_df <- stack(goEnrichmentQ_nem_up_and_sag_down_recQ5_genes)
colnames(goEnrichmentQ_nem_up_and_sag_down_recQ5_genes_df) <- c("GeneID", "GO.ID")

# Save
write.csv(goEnrichmentQ_nem_up_and_sag_down_recQ5_genes_df, "goEnrichmentQ_nem_up_and_sag_down_recQ5_genes.csv", row.names = FALSE)


# Example GO enrichment result Q7 / recover
sig_go_recQ7 <- goEnrichmentQ_nem_down_more_against_sag_rec_Q7$GO.ID

# Get significant genes per GO term
# Corrected function
getSigGenesFromGO <- function(GOterm, GOdata) {
  all_genes <- genesInTerm(GOdata)[[GOterm]]
  sig_genes <- all_genes[all_genes %in% sigGenes(GOdata)]
  return(sig_genes)
}

goEnrichmentQ_nem_down_more_against_sag_recQ7_genes <- lapply(sig_go_recQ7, getSigGenesFromGO, GOdata = tGOdata_rec7)
names(goEnrichmentQ_nem_down_more_against_sag_recQ7_genes) <- sig_go_recQ7

# Flatten to a data frame
goEnrichmentQ_nem_down_more_against_sag_recQ7_genes_df <- stack(goEnrichmentQ_nem_down_more_against_sag_recQ7_genes)
colnames(goEnrichmentQ_nem_down_more_against_sag_recQ7_genes_df) <- c("GeneID", "GO.ID")

# Save
write.csv(goEnrichmentQ_nem_down_more_against_sag_recQ7_genes_df, "goEnrichmentQ_nem_down_more_against_sag_recQ7_genes.csv", row.names = FALSE)


# Example GO enrichment result Q8 / recover
sig_go_recQ8 <- goEnrichmentQ_nem_up_more_against_sag_rec_Q8$GO.ID

# Get significant genes per GO term
# Corrected function
getSigGenesFromGO <- function(GOterm, GOdata) {
  all_genes <- genesInTerm(GOdata)[[GOterm]]
  sig_genes <- all_genes[all_genes %in% sigGenes(GOdata)]
  return(sig_genes)
}

goEnrichmentQ_nem_up_more_against_sag_recQ8_genes <- lapply(sig_go_recQ8, getSigGenesFromGO, GOdata = tGOdata_rec8)
names(goEnrichmentQ_nem_up_more_against_sag_recQ8_genes) <- sig_go_recQ8

# Flatten to a data frame
goEnrichmentQ_nem_up_more_against_sag_recQ8_genes_df <- stack(goEnrichmentQ_nem_up_more_against_sag_recQ8_genes)
colnames(goEnrichmentQ_nem_up_more_against_sag_recQ8_genes_df) <- c("GeneID", "GO.ID")

# Save
write.csv(goEnrichmentQ_nem_up_more_against_sag_recQ8_genes_df, "goEnrichmentQ_nem_up_more_against_sag_recQ8_genes.csv", row.names = FALSE)



########################################## End here on 11th July 2025 ################################################
########################################## End here on 11th July 2025 ################################################
########################################## End here on 11th July 2025 ################################################
########################################## End here on 11th July 2025 ################################################
########################################## End here on 11th July 2025 ################################################


##### starting 6th January 2026, here is the new analysis #######
##### A. nem and A. sag in drought vs control --- shared response in up-up regulation
# 1) Keep only genes with valid Arabidopsis IDs for org.At.tair.db
#combined_species_df_sag_nem_df <- combined_species_df_sag_nem[!is.na(combined_species_df_sag_nem$At), ]

# 2) Define genes of interest: green AND quadrant in Q1 or Q8
genes_of_interest_sag_nem_df <- combined_species_df_sag_nem_df$color == "green" & combined_species_df_sag_nem_df$quadrant %in% c("Q1", "Q8")

# 3) Build geneList for topGO (universe = all genes in df)
allGenes_numeric <- as.integer(genes_of_interest_sag_nem_df)
names(allGenes_numeric) <- df$At

# (optional but good) remove duplicated TAIR IDs if they exist
allGenes_numeric <- allGenes_numeric[!duplicated(names(allGenes_numeric))]

# 4) Create topGO object
tGOdata_sag_nem_df <- new("topGOdata",
                          description = "Green genes in Q1 or Q8",
                          ontology = "BP",
                          allGenes = allGenes_numeric,
                          geneSel = function(x) x == 1,
                          nodeSize = 5,
                          mapping = "org.At.tair.db",
                          annot = annFUN.org)

# 5) Run enrichment (Fisher + elim)
result_fisher_sag_nem_df <- runTest(tGOdata_sag_nem_df, algorithm = "elim", statistic = "fisher")

# 6) Output table
goTable_sag_nem_upup_df <- GenTable(tGOdata_sag_nem_df,
                                    fisher = result_fisher_sag_nem_df,
                                    orderBy = "fisher",
                                    topNodes = 20)

goTable_sag_nem_upup_df

write.csv(goTable_sag_nem_upup_df, "goTable_sag_nem_upup_df.csv")


##### A. nem and A. sag in drought vs control --- shared response in down-down regulation
# 1) Keep only genes with valid Arabidopsis IDs for org.At.tair.db
#combined_species_df_sag_nem_df <- combined_species_df_sag_nem[!is.na(combined_species_df_sag_nem$At), ]

# 2) Define genes of interest: green AND quadrant in Q1 or Q8
genes_of_interest_sag_nem_df <- combined_species_df_sag_nem_df$color == "green" & combined_species_df_sag_nem_df$quadrant %in% c("Q7", "Q4")

# 3) Build geneList for topGO (universe = all genes in df)
allGenes_numeric <- as.integer(genes_of_interest_sag_nem_df)
names(allGenes_numeric) <- df$At

# (optional but good) remove duplicated TAIR IDs if they exist
allGenes_numeric <- allGenes_numeric[!duplicated(names(allGenes_numeric))]

# 4) Create topGO object
tGOdata_sag_nem_df <- new("topGOdata",
                          description = "Green genes in Q7 or Q4",
                          ontology = "BP",
                          allGenes = allGenes_numeric,
                          geneSel = function(x) x == 1,
                          nodeSize = 5,
                          mapping = "org.At.tair.db",
                          annot = annFUN.org)

# 5) Run enrichment (Fisher + elim)
result_fisher_sag_nem_df <- runTest(tGOdata_sag_nem_df, algorithm = "elim", statistic = "fisher")

# 6) Output table
goTable_sag_nem_downdown_df <- GenTable(tGOdata_sag_nem_df,
                                        fisher = result_fisher_sag_nem_df,
                                        orderBy = "fisher",
                                        topNodes = 20)

goTable_sag_nem_downdown_df

write.csv(goTable_sag_nem_downdown_df, "goTable_sag_nem_downdown_df.csv")


##### A. nem and A. sag in drought vs control --- response more up regulation in A. sagittata
# 1) Keep only genes with valid Arabidopsis IDs for org.At.tair.db
#combined_species_df_sag_nem_df <- combined_species_df_sag_nem[!is.na(combined_species_df_sag_nem$At), ]

# 2) Define genes of interest: green AND quadrant in Q1 or Q8
genes_of_interest_sag_nem_df <- combined_species_df_sag_nem_df$color == "red" & combined_species_df_sag_nem_df$quadrant %in% c("Q1")

# 3) Build geneList for topGO (universe = all genes in df)
allGenes_numeric <- as.integer(genes_of_interest_sag_nem_df)
names(allGenes_numeric) <- df$At

# (optional but good) remove duplicated TAIR IDs if they exist
allGenes_numeric <- allGenes_numeric[!duplicated(names(allGenes_numeric))]

# 4) Create topGO object
tGOdata_sag_nem_df <- new("topGOdata",
                          description = "Green genes in Q1",
                          ontology = "BP",
                          allGenes = allGenes_numeric,
                          geneSel = function(x) x == 1,
                          nodeSize = 5,
                          mapping = "org.At.tair.db",
                          annot = annFUN.org)

# 5) Run enrichment (Fisher + elim)
result_fisher_sag_nem_df <- runTest(tGOdata_sag_nem_df, algorithm = "elim", statistic = "fisher")

# 6) Output table
goTable_stress_more_up_sag_df <- GenTable(tGOdata_sag_nem_df,
                                          fisher = result_fisher_sag_nem_df,
                                          orderBy = "fisher",
                                          topNodes = 20)

goTable_stress_more_up_sag_df


##### A. nem and A. sag in drought vs control --- response more up regulation in A. nemorensis
# 1) Keep only genes with valid Arabidopsis IDs for org.At.tair.db
#combined_species_df_sag_nem_df <- combined_species_df_sag_nem[!is.na(combined_species_df_sag_nem$At), ]

# 2) Define genes of interest: green AND quadrant in Q1 or Q8
genes_of_interest_sag_nem_df <- combined_species_df_sag_nem_df$color == "red" & combined_species_df_sag_nem_df$quadrant %in% c("Q8")

# 3) Build geneList for topGO (universe = all genes in df)
allGenes_numeric <- as.integer(genes_of_interest_sag_nem_df)
names(allGenes_numeric) <- df$At

# (optional but good) remove duplicated TAIR IDs if they exist
allGenes_numeric <- allGenes_numeric[!duplicated(names(allGenes_numeric))]

# 4) Create topGO object
tGOdata_sag_nem_df <- new("topGOdata",
                          description = "Green genes in Q8",
                          ontology = "BP",
                          allGenes = allGenes_numeric,
                          geneSel = function(x) x == 1,
                          nodeSize = 5,
                          mapping = "org.At.tair.db",
                          annot = annFUN.org)

# 5) Run enrichment (Fisher + elim)
result_fisher_sag_nem_df <- runTest(tGOdata_sag_nem_df, algorithm = "elim", statistic = "fisher")

# 6) Output table
goTable_stress_more_up_nem_df <- GenTable(tGOdata_sag_nem_df,
                                          fisher = result_fisher_sag_nem_df,
                                          orderBy = "fisher",
                                          topNodes = 20)

goTable_stress_more_up_nem_df


##### A. nem and A. sag in drought vs control --- response up regulation in A. sagittata, and down in A. nemorensis
# 1) Keep only genes with valid Arabidopsis IDs for org.At.tair.db
#combined_species_df_sag_nem_df <- combined_species_df_sag_nem[!is.na(combined_species_df_sag_nem$At), ]

# 2) Define genes of interest: green AND quadrant in Q1 or Q8
genes_of_interest_sag_nem_df <- combined_species_df_sag_nem_df$color == "red" & combined_species_df_sag_nem_df$quadrant %in% c("Q2")

# 3) Build geneList for topGO (universe = all genes in df)
allGenes_numeric <- as.integer(genes_of_interest_sag_nem_df)
names(allGenes_numeric) <- df$At

# (optional but good) remove duplicated TAIR IDs if they exist
allGenes_numeric <- allGenes_numeric[!duplicated(names(allGenes_numeric))]

# 4) Create topGO object
tGOdata_sag_nem_df <- new("topGOdata",
                          description = "Green genes in Q2",
                          ontology = "BP",
                          allGenes = allGenes_numeric,
                          geneSel = function(x) x == 1,
                          nodeSize = 5,
                          mapping = "org.At.tair.db",
                          annot = annFUN.org)

# 5) Run enrichment (Fisher + elim)
result_fisher_sag_nem_df <- runTest(tGOdata_sag_nem_df, algorithm = "elim", statistic = "fisher")

# 6) Output table
goTable_stress_up_sag_down_nem_df <- GenTable(tGOdata_sag_nem_df,
                                              fisher = result_fisher_sag_nem_df,
                                              orderBy = "fisher",
                                              topNodes = 20)

goTable_stress_up_sag_down_nem_df



##### A. nem and A. sag in drought vs control --- response up regulation in A. nemorensis, and down in A. sagittata
# 1) Keep only genes with valid Arabidopsis IDs for org.At.tair.db
#combined_species_df_sag_nem_df <- combined_species_df_sag_nem[!is.na(combined_species_df_sag_nem$At), ]

# 2) Define genes of interest: green AND quadrant in Q1 or Q8
genes_of_interest_sag_nem_df <- combined_species_df_sag_nem_df$color == "red" & combined_species_df_sag_nem_df$quadrant %in% c("Q5")

# 3) Build geneList for topGO (universe = all genes in df)
allGenes_numeric <- as.integer(genes_of_interest_sag_nem_df)
names(allGenes_numeric) <- df$At

# (optional but good) remove duplicated TAIR IDs if they exist
allGenes_numeric <- allGenes_numeric[!duplicated(names(allGenes_numeric))]

# 4) Create topGO object
tGOdata_sag_nem_df <- new("topGOdata",
                          description = "Green genes in Q2",
                          ontology = "BP",
                          allGenes = allGenes_numeric,
                          geneSel = function(x) x == 1,
                          nodeSize = 5,
                          mapping = "org.At.tair.db",
                          annot = annFUN.org)

# 5) Run enrichment (Fisher + elim)
result_fisher_sag_nem_df <- runTest(tGOdata_sag_nem_df, algorithm = "elim", statistic = "fisher")

# 6) Output table
goTable_stress_up_nem_down_sag_df <- GenTable(tGOdata_sag_nem_df,
                                              fisher = result_fisher_sag_nem_df,
                                              orderBy = "fisher",
                                              topNodes = 20)

goTable_stress_up_nem_down_sag_df



##### A. nem and A. sag in drought vs control --- response more down regulation in A. nemorensis, and less down in A. sagittata
# 1) Keep only genes with valid Arabidopsis IDs for org.At.tair.db
combined_species_df_sag_nem_df <- combined_species_df_sag_nem[!is.na(combined_species_df_sag_nem$At), ]

# 2) Define genes of interest: green AND quadrant in Q1 or Q8
genes_of_interest_sag_nem_df <- combined_species_df_sag_nem_df$color == "red" & combined_species_df_sag_nem_df$quadrant %in% c("Q7")

# 3) Build geneList for topGO (universe = all genes in df)
allGenes_numeric <- as.integer(genes_of_interest_sag_nem_df)
names(allGenes_numeric) <- df$At

# (optional but good) remove duplicated TAIR IDs if they exist
allGenes_numeric <- allGenes_numeric[!duplicated(names(allGenes_numeric))]

# 4) Create topGO object
tGOdata_sag_nem_df <- new("topGOdata",
                          description = "Green genes in Q2",
                          ontology = "BP",
                          allGenes = allGenes_numeric,
                          geneSel = function(x) x == 1,
                          nodeSize = 5,
                          mapping = "org.At.tair.db",
                          annot = annFUN.org)

# 5) Run enrichment (Fisher + elim)
result_fisher_sag_nem_df <- runTest(tGOdata_sag_nem_df, algorithm = "elim", statistic = "fisher")

# 6) Output table
goTable_stress_more_down_nem_df <- GenTable(tGOdata_sag_nem_df,
                                            fisher = result_fisher_sag_nem_df,
                                            orderBy = "fisher",
                                            topNodes = 20)

goTable_stress_more_down_nem_df


##### A. nem and A. sag in drought vs control --- response more down regulation in A. sagittata, and less down in A. nemorensis
# 1) Keep only genes with valid Arabidopsis IDs for org.At.tair.db
#combined_species_df_sag_nem_df <- combined_species_df_sag_nem[!is.na(combined_species_df_sag_nem$At), ]

# 2) Define genes of interest: green AND quadrant in Q1 or Q8
genes_of_interest_sag_nem_df <- combined_species_df_sag_nem_df$color == "red" & combined_species_df_sag_nem_df$quadrant %in% c("Q4")

# 3) Build geneList for topGO (universe = all genes in df)
allGenes_numeric <- as.integer(genes_of_interest_sag_nem_df)
names(allGenes_numeric) <- df$At

# (optional but good) remove duplicated TAIR IDs if they exist
allGenes_numeric <- allGenes_numeric[!duplicated(names(allGenes_numeric))]

# 4) Create topGO object
tGOdata_sag_nem_df <- new("topGOdata",
                          description = "Green genes in Q2",
                          ontology = "BP",
                          allGenes = allGenes_numeric,
                          geneSel = function(x) x == 1,
                          nodeSize = 5,
                          mapping = "org.At.tair.db",
                          annot = annFUN.org)

# 5) Run enrichment (Fisher + elim)
result_fisher_sag_nem_df <- runTest(tGOdata_sag_nem_df, algorithm = "elim", statistic = "fisher")

# 6) Output table
goTable_stress_more_down_sag_df <- GenTable(tGOdata_sag_nem_df,
                                            fisher = result_fisher_sag_nem_df,
                                            orderBy = "fisher",
                                            topNodes = 20)

goTable_stress_more_down_sag_df




write.csv(goTable_stress_up_nem_down_sag_df, "goTable_stress_up_nem_down_sag_df.csv")
write.csv(goTable_stress_up_sag_down_nem_df,  "goTable_stress_up_sag_down_nem_df.csv")
write.csv(goTable_stress_more_up_sag_df, "goTable_stress_more_up_sag_df.csv")
write.csv(goTable_stress_more_up_nem_df, "goTable_stress_more_up_nem_df.csv")
write.csv(goTable_stress_more_down_sag_df, "goTable_stress_more_down_sag_df.csv")
write.csv(goTable_stress_more_down_sag_df, "goTable_stress_more_down_sag_df.csv")





#######################################################
##### A. nem and A. sag in recovery vs control --- shared response in up-up regulation
# 1) Keep only genes with valid Arabidopsis IDs for org.At.tair.db
#combined_species_df_sag_nem_surv_df <- combined_species_df_sag_nem_surv[!is.na(combined_species_df_sag_nem_surv$At), ]

# 2) Define genes of interest: green AND quadrant in Q1 or Q8
genes_of_interest_sag_nem_surv_df <- combined_species_df_sag_nem_surv_df$color == "green" & combined_species_df_sag_nem_surv_df$quadrant %in% c("Q1", "Q8")

# 3) Build geneList for topGO (universe = all genes in df)
allGenes_numeric <- as.integer(genes_of_interest_sag_nem_surv_df)
names(allGenes_numeric) <- df$At

# (optional but good) remove duplicated TAIR IDs if they exist
allGenes_numeric <- allGenes_numeric[!duplicated(names(allGenes_numeric))]

# 4) Create topGO object
tGOdata_sag_nem_surv_df <- new("topGOdata",
                               description = "Green genes in Q1 or Q8",
                               ontology = "BP",
                               allGenes = allGenes_numeric,
                               geneSel = function(x) x == 1,
                               nodeSize = 5,
                               mapping = "org.At.tair.db",
                               annot = annFUN.org)

# 5) Run enrichment (Fisher + elim)
result_fisher_sag_nem_surv_df <- runTest(tGOdata_sag_nem_surv_df, algorithm = "elim", statistic = "fisher")

# 6) Output table
goTable_recovery_sag_nem_upup_df <- GenTable(tGOdata_sag_nem_surv_df,
                                             fisher = result_fisher_sag_nem_surv_df,
                                             orderBy = "fisher",
                                             topNodes = 20)

goTable_recovery_sag_nem_upup_df


##### A. nem and A. sag in recover vs control --- shared response in down-down regulation
# 1) Keep only genes with valid Arabidopsis IDs for org.At.tair.db
#combined_species_df_sag_nem_surv_df <- combined_species_df_sag_nem_surv[!is.na(combined_species_df_sag_nem_surv$At), ]

# 2) Define genes of interest: green AND quadrant in Q1 or Q8
genes_of_interest_sag_nem_surv_df <- combined_species_df_sag_nem_surv_df$color == "green" & combined_species_df_sag_nem_surv_df$quadrant %in% c("Q7", "Q4")

# 3) Build geneList for topGO (universe = all genes in df)
allGenes_numeric <- as.integer(genes_of_interest_sag_nem_surv_df)
names(allGenes_numeric) <- df$At

# (optional but good) remove duplicated TAIR IDs if they exist
allGenes_numeric <- allGenes_numeric[!duplicated(names(allGenes_numeric))]

# 4) Create topGO object
tGOdata_sag_nem_surv_df <- new("topGOdata",
                               description = "Green genes in Q7 or Q4",
                               ontology = "BP",
                               allGenes = allGenes_numeric,
                               geneSel = function(x) x == 1,
                               nodeSize = 5,
                               mapping = "org.At.tair.db",
                               annot = annFUN.org)

# 5) Run enrichment (Fisher + elim)
result_fisher_sag_nem_surv_df <- runTest(tGOdata_sag_nem_surv_df, algorithm = "elim", statistic = "fisher")

# 6) Output table
goTable_recovery_sag_nem_downdown_df <- GenTable(tGOdata_sag_nem_surv_df,
                                                 fisher = result_fisher_sag_nem_surv_df,
                                                 orderBy = "fisher",
                                                 topNodes = 20)

goTable_recovery_sag_nem_downdown_df



##### A. nem and A. sag in recovery vs control --- response more up regulation in A. sagittata
# 1) Keep only genes with valid Arabidopsis IDs for org.At.tair.db
#combined_species_df_sag_nem_surv_df <- combined_species_df_sag_nem_surv[!is.na(combined_species_df_sag_nem_surv$At), ]

# 2) Define genes of interest: green AND quadrant in Q1 or Q8
genes_of_interest_sag_nem_surv_df <- combined_species_df_sag_nem_surv_df$color == "red" & combined_species_df_sag_nem_surv_df$quadrant %in% c("Q1")

# 3) Build geneList for topGO (universe = all genes in df)
allGenes_numeric <- as.integer(genes_of_interest_sag_nem_surv_df)
names(allGenes_numeric) <- df$At

# (optional but good) remove duplicated TAIR IDs if they exist
allGenes_numeric <- allGenes_numeric[!duplicated(names(allGenes_numeric))]

# 4) Create topGO object
tGOdata_sag_nem_surv_df <- new("topGOdata",
                               description = "Green genes in Q1",
                               ontology = "BP",
                               allGenes = allGenes_numeric,
                               geneSel = function(x) x == 1,
                               nodeSize = 5,
                               mapping = "org.At.tair.db",
                               annot = annFUN.org)

# 5) Run enrichment (Fisher + elim)
result_fisher_sag_nem_surv_df <- runTest(tGOdata_sag_nem_surv_df, algorithm = "elim", statistic = "fisher")

# 6) Output table
goTable_recovery_more_up_sag_df <- GenTable(tGOdata_sag_nem_surv_df,
                                            fisher = result_fisher_sag_nem_surv_df,
                                            orderBy = "fisher",
                                            topNodes = 20)

goTable_recovery_more_up_sag_df


##### A. nem and A. sag in recovery vs control --- response more up regulation in A. nemorensis
# 1) Keep only genes with valid Arabidopsis IDs for org.At.tair.db
#combined_species_df_sag_nem_surv_df <- combined_species_df_sag_nem_surv[!is.na(combined_species_df_sag_nem_surv$At), ]

# 2) Define genes of interest: green AND quadrant in Q1 or Q8
genes_of_interest_sag_nem_surv_df <- combined_species_df_sag_nem_surv_df$color == "red" & combined_species_df_sag_nem_surv_df$quadrant %in% c("Q8")

# 3) Build geneList for topGO (universe = all genes in df)
allGenes_numeric <- as.integer(genes_of_interest_sag_nem_surv_df)
names(allGenes_numeric) <- df$At

# (optional but good) remove duplicated TAIR IDs if they exist
allGenes_numeric <- allGenes_numeric[!duplicated(names(allGenes_numeric))]

# 4) Create topGO object
tGOdata_sag_nem_surv_df <- new("topGOdata",
                               description = "Green genes in Q8",
                               ontology = "BP",
                               allGenes = allGenes_numeric,
                               geneSel = function(x) x == 1,
                               nodeSize = 5,
                               mapping = "org.At.tair.db",
                               annot = annFUN.org)

# 5) Run enrichment (Fisher + elim)
result_fisher_sag_nem_surv_df <- runTest(tGOdata_sag_nem_surv_df, algorithm = "elim", statistic = "fisher")

# 6) Output table
goTable_recovery_more_up_nem_df <- GenTable(tGOdata_sag_nem_surv_df,
                                            fisher = result_fisher_sag_nem_surv_df,
                                            orderBy = "fisher",
                                            topNodes = 20)

goTable_recovery_more_up_nem_df


##### A. nem and A. sag in recovery vs control --- response up regulation in A. sagittata, and down in A. nemorensis
# 1) Keep only genes with valid Arabidopsis IDs for org.At.tair.db
#combined_species_df_sag_nem_surv_df <- combined_species_df_sag_nem_surv[!is.na(combined_species_df_sag_nem_surv$At), ]

# 2) Define genes of interest: green AND quadrant in Q1 or Q8
genes_of_interest_sag_nem_surv_df <- combined_species_df_sag_nem_surv_df$color == "red" & combined_species_df_sag_nem_surv_df$quadrant %in% c("Q2")

# 3) Build geneList for topGO (universe = all genes in df)
allGenes_numeric <- as.integer(genes_of_interest_sag_nem_surv_df)
names(allGenes_numeric) <- df$At

# (optional but good) remove duplicated TAIR IDs if they exist
allGenes_numeric <- allGenes_numeric[!duplicated(names(allGenes_numeric))]

# 4) Create topGO object
tGOdata_sag_nem_surv_df <- new("topGOdata",
                               description = "Green genes in Q2",
                               ontology = "BP",
                               allGenes = allGenes_numeric,
                               geneSel = function(x) x == 1,
                               nodeSize = 5,
                               mapping = "org.At.tair.db",
                               annot = annFUN.org)

# 5) Run enrichment (Fisher + elim)
result_fisher_sag_nem_surv_df <- runTest(tGOdata_sag_nem_surv_df, algorithm = "elim", statistic = "fisher")

# 6) Output table
goTable_recovery_up_sag_down_nem_df <- GenTable(tGOdata_sag_nem_surv_df,
                                                fisher = result_fisher_sag_nem_surv_df,
                                                orderBy = "fisher",
                                                topNodes = 20)

goTable_recovery_up_sag_down_nem_df



##### A. nem and A. sag in recovery vs control --- response up regulation in A. nemorensis, and down in A. sagittata
# 1) Keep only genes with valid Arabidopsis IDs for org.At.tair.db
combined_species_df_sag_nem_surv_df <- combined_species_df_sag_nem_surv[!is.na(combined_species_df_sag_nem_surv$At), ]

# 2) Define genes of interest: green AND quadrant in Q1 or Q8
genes_of_interest_sag_nem_surv_df <- combined_species_df_sag_nem_surv_df$color == "red" & combined_species_df_sag_nem_surv_df$quadrant %in% c("Q5")

# 3) Build geneList for topGO (universe = all genes in df)
allGenes_numeric <- as.integer(genes_of_interest_sag_nem_surv_df)
names(allGenes_numeric) <- df$At

# (optional but good) remove duplicated TAIR IDs if they exist
allGenes_numeric <- allGenes_numeric[!duplicated(names(allGenes_numeric))]

# 4) Create topGO object
tGOdata_sag_nem_surv_df <- new("topGOdata",
                               description = "Green genes in Q2",
                               ontology = "BP",
                               allGenes = allGenes_numeric,
                               geneSel = function(x) x == 1,
                               nodeSize = 5,
                               mapping = "org.At.tair.db",
                               annot = annFUN.org)

# 5) Run enrichment (Fisher + elim)
result_fisher_sag_nem_surv_df <- runTest(tGOdata_sag_nem_surv_df, algorithm = "elim", statistic = "fisher")

# 6) Output table
goTable_recovery_up_nem_down_sag_df <- GenTable(tGOdata_sag_nem_surv_df,
                                                fisher = result_fisher_sag_nem_surv_df,
                                                orderBy = "fisher",
                                                topNodes = 20)

goTable_recovery_up_nem_down_sag_df



##### A. nem and A. sag in recovery vs control --- response more down regulation in A. nemorensis, and less down in A. sagittata
# 1) Keep only genes with valid Arabidopsis IDs for org.At.tair.db
combined_species_df_sag_nem_surv_df <- combined_species_df_sag_nem_surv[!is.na(combined_species_df_sag_nem_surv$At), ]

# 2) Define genes of interest: green AND quadrant in Q1 or Q8
genes_of_interest_sag_nem_surv_df <- combined_species_df_sag_nem_surv_df$color == "red" & combined_species_df_sag_nem_surv_df$quadrant %in% c("Q7")

# 3) Build geneList for topGO (universe = all genes in df)
allGenes_numeric <- as.integer(genes_of_interest_sag_nem_surv_df)
names(allGenes_numeric) <- df$At

# (optional but good) remove duplicated TAIR IDs if they exist
allGenes_numeric <- allGenes_numeric[!duplicated(names(allGenes_numeric))]

# 4) Create topGO object
tGOdata_sag_nem_surv_df <- new("topGOdata",
                               description = "Green genes in Q2",
                               ontology = "BP",
                               allGenes = allGenes_numeric,
                               geneSel = function(x) x == 1,
                               nodeSize = 5,
                               mapping = "org.At.tair.db",
                               annot = annFUN.org)

# 5) Run enrichment (Fisher + elim)
result_fisher_sag_nem_surv_df <- runTest(tGOdata_sag_nem_surv_df, algorithm = "elim", statistic = "fisher")

# 6) Output table
goTable_recovery_more_down_nem_df <- GenTable(tGOdata_sag_nem_surv_df,
                                              fisher = result_fisher_sag_nem_surv_df,
                                              orderBy = "fisher",
                                              topNodes = 20)

goTable_recovery_more_down_nem_df



##### A. nem and A. sag in recovery vs control --- response more down regulation in A. sagittata, and less down in A. nemorensis
# 1) Keep only genes with valid Arabidopsis IDs for org.At.tair.db
combined_species_df_sag_nem_surv_df <- combined_species_df_sag_nem_surv[!is.na(combined_species_df_sag_nem_surv$At), ]

# 2) Define genes of interest: green AND quadrant in Q1 or Q8
genes_of_interest_sag_nem_surv_df <- combined_species_df_sag_nem_surv_df$color == "red" & combined_species_df_sag_nem_surv_df$quadrant %in% c("Q4")

# 3) Build geneList for topGO (universe = all genes in df)
allGenes_numeric <- as.integer(genes_of_interest_sag_nem_surv_df)
names(allGenes_numeric) <- df$At

# (optional but good) remove duplicated TAIR IDs if they exist
allGenes_numeric <- allGenes_numeric[!duplicated(names(allGenes_numeric))]

# 4) Create topGO object
tGOdata_sag_nem_surv_df <- new("topGOdata",
                               description = "Green genes in Q2",
                               ontology = "BP",
                               allGenes = allGenes_numeric,
                               geneSel = function(x) x == 1,
                               nodeSize = 5,
                               mapping = "org.At.tair.db",
                               annot = annFUN.org)

# 5) Run enrichment (Fisher + elim)
result_fisher_sag_nem_surv_df <- runTest(tGOdata_sag_nem_surv_df, algorithm = "elim", statistic = "fisher")

# 6) Output table
goTable_recovery_more_down_sag_df <- GenTable(tGOdata_sag_nem_surv_df,
                                              fisher = result_fisher_sag_nem_surv_df,
                                              orderBy = "fisher",
                                              topNodes = 20)

goTable_recovery_more_down_sag_df


write.csv(goTable_recovery_up_nem_down_sag_df, "goTable_recovery_up_nem_down_sag_df.csv")
write.csv(goTable_recovery_up_sag_down_nem_df,  "goTable_recovery_up_sag_down_nem_df.csv")
write.csv(goTable_recovery_more_up_sag_df, "goTable_recovery_more_up_sag_df.csv")
write.csv(goTable_recovery_more_up_nem_df, "goTable_recovery_more_up_nem_df.csv")
write.csv(goTable_recovery_more_down_sag_df, "goTable_recovery_more_down_sag_df.csv")
write.csv(goTable_recovery_more_down_sag_df, "goTable_recovery_more_down_sag_df.csv")





#..........................................below is pending...

#######################################################
##### A. nem and A. sag in recovery vs stress --- shared response in up-up regulation
# 1) Keep only genes with valid Arabidopsis IDs for org.At.tair.db

# Merge expression data with orthologue mapping
#combined_species_df_sag_nem_surv_wilt <- combined_species_df_sag_nem_surv_wilt %>%
#  left_join(new_dataframe1, by = c("gene_id" = "arabis_cleaned"))


combined_species_df_sag_nem_surv_wilt_df <- combined_species_df_sag_nem_surv_wilt[!is.na(combined_species_df_sag_nem_surv_wilt$At), ]

# 2) Define genes of interest: green AND quadrant in Q1 or Q8
genes_of_interest_sag_nem_surv_wilt_df <- combined_species_df_sag_nem_surv_wilt_df$color == "green" & combined_species_df_sag_nem_surv_wilt_df$quadrant %in% c("Q1", "Q8")

# 3) Build geneList for topGO (universe = all genes in df)
allGenes_numeric <- as.integer(genes_of_interest_sag_nem_surv_wilt_df)
names(allGenes_numeric) <- df$At

# (optional but good) remove duplicated TAIR IDs if they exist
allGenes_numeric <- allGenes_numeric[!duplicated(names(allGenes_numeric))]

# 4) Create topGO object
tGOdata_sag_nem_surv_wilt_df <- new("topGOdata",
                                    description = "Green genes in Q1 or Q8",
                                    ontology = "BP",
                                    allGenes = allGenes_numeric,
                                    geneSel = function(x) x == 1,
                                    nodeSize = 5,
                                    mapping = "org.At.tair.db",
                                    annot = annFUN.org)

# 5) Run enrichment (Fisher + elim)
result_fisher_sag_nem_surv_wilt_df <- runTest(tGOdata_sag_nem_surv_wilt_df, algorithm = "elim", statistic = "fisher")

# 6) Output table
goTable_recovery_stress_sag_nem_upup_df <- GenTable(tGOdata_sag_nem_surv_wilt_df,
                                                    fisher = result_fisher_sag_nem_surv_wilt_df,
                                                    orderBy = "fisher",
                                                    topNodes = 20)

goTable_recovery_stress_sag_nem_upup_df


##### A. nem and A. sag in recover vs control --- shared response in down-down regulation
# 1) Keep only genes with valid Arabidopsis IDs for org.At.tair.db
combined_species_df_sag_nem_surv_wilt_df <- combined_species_df_sag_nem_surv_wilt[!is.na(combined_species_df_sag_nem_surv_wilt$At), ]

# 2) Define genes of interest: green AND quadrant in Q1 or Q8
genes_of_interest_sag_nem_surv_wilt_df <- combined_species_df_sag_nem_surv_wilt_df$color == "green" & combined_species_df_sag_nem_surv_wilt_df$quadrant %in% c("Q7", "Q4")

# 3) Build geneList for topGO (universe = all genes in df)
allGenes_numeric <- as.integer(genes_of_interest_sag_nem_surv_wilt_df)
names(allGenes_numeric) <- df$At

# (optional but good) remove duplicated TAIR IDs if they exist
allGenes_numeric <- allGenes_numeric[!duplicated(names(allGenes_numeric))]

# 4) Create topGO object
tGOdata_sag_nem_surv_wilt_df <- new("topGOdata",
                                    description = "Green genes in Q7 or Q4",
                                    ontology = "BP",
                                    allGenes = allGenes_numeric,
                                    geneSel = function(x) x == 1,
                                    nodeSize = 5,
                                    mapping = "org.At.tair.db",
                                    annot = annFUN.org)

# 5) Run enrichment (Fisher + elim)
result_fisher_sag_nem_surv_wilt_df <- runTest(tGOdata_sag_nem_surv_wilt_df, algorithm = "elim", statistic = "fisher")

# 6) Output table
goTable_recovery_stress_sag_nem_downdown_df <- GenTable(tGOdata_sag_nem_surv_wilt_df,
                                                        fisher = result_fisher_sag_nem_surv_wilt_df,
                                                        orderBy = "fisher",
                                                        topNodes = 20)

goTable_recovery_stress_sag_nem_downdown_df



##### A. nem and A. sag in recovery vs stress --- response more up regulation in A. sagittata
# 1) Keep only genes with valid Arabidopsis IDs for org.At.tair.db
combined_species_df_sag_nem_surv_wilt_df <- combined_species_df_sag_nem_surv_wilt[!is.na(combined_species_df_sag_nem_surv_wilt$At), ]

# 2) Define genes of interest: green AND quadrant in Q1 or Q8
genes_of_interest_sag_nem_surv_wilt_df <- combined_species_df_sag_nem_surv_wilt_df$color == "red" & combined_species_df_sag_nem_surv_wilt_df$quadrant %in% c("Q1")

# 3) Build geneList for topGO (universe = all genes in df)
allGenes_numeric <- as.integer(genes_of_interest_sag_nem_surv_wilt_df)
names(allGenes_numeric) <- df$At

# (optional but good) remove duplicated TAIR IDs if they exist
allGenes_numeric <- allGenes_numeric[!duplicated(names(allGenes_numeric))]

# 4) Create topGO object
tGOdata_sag_nem_surv_wilt_df <- new("topGOdata",
                                    description = "Green genes in Q1",
                                    ontology = "BP",
                                    allGenes = allGenes_numeric,
                                    geneSel = function(x) x == 1,
                                    nodeSize = 5,
                                    mapping = "org.At.tair.db",
                                    annot = annFUN.org)

# 5) Run enrichment (Fisher + elim)
result_fisher_sag_nem_surv_wilt_df <- runTest(tGOdata_sag_nem_surv_wilt_df, algorithm = "elim", statistic = "fisher")

# 6) Output table
goTable_recovery_stress_more_up_sag_df <- GenTable(tGOdata_sag_nem_surv_wilt_df,
                                                   fisher = result_fisher_sag_nem_surv_wilt_df,
                                                   orderBy = "fisher",
                                                   topNodes = 20)

goTable_recovery_stress_more_up_sag_df


##### A. nem and A. sag in recovery vs stress --- response more up regulation in A. nemorensis
# 1) Keep only genes with valid Arabidopsis IDs for org.At.tair.db
combined_species_df_sag_nem_surv_wilt_df <- combined_species_df_sag_nem_surv_wilt[!is.na(combined_species_df_sag_nem_surv_wilt$At), ]

# 2) Define genes of interest: green AND quadrant in Q1 or Q8
genes_of_interest_sag_nem_surv_wilt_df <- combined_species_df_sag_nem_surv_wilt_df$color == "red" & combined_species_df_sag_nem_surv_wilt_df$quadrant %in% c("Q8")

# 3) Build geneList for topGO (universe = all genes in df)
allGenes_numeric <- as.integer(genes_of_interest_sag_nem_surv_wilt_df)
names(allGenes_numeric) <- df$At

# (optional but good) remove duplicated TAIR IDs if they exist
allGenes_numeric <- allGenes_numeric[!duplicated(names(allGenes_numeric))]

# 4) Create topGO object
tGOdata_sag_nem_surv_wilt_df <- new("topGOdata",
                                    description = "Green genes in Q8",
                                    ontology = "BP",
                                    allGenes = allGenes_numeric,
                                    geneSel = function(x) x == 1,
                                    nodeSize = 5,
                                    mapping = "org.At.tair.db",
                                    annot = annFUN.org)

# 5) Run enrichment (Fisher + elim)
result_fisher_sag_nem_surv_wilt_df <- runTest(tGOdata_sag_nem_surv_wilt_df, algorithm = "elim", statistic = "fisher")

# 6) Output table
goTable_recovery_stress_more_up_nem_df <- GenTable(tGOdata_sag_nem_surv_wilt_df,
                                                   fisher = result_fisher_sag_nem_surv_wilt_df,
                                                   orderBy = "fisher",
                                                   topNodes = 20)

goTable_recovery_stress_more_up_nem_df


##### A. nem and A. sag in recovery vs stress --- response up regulation in A. sagittata, and down in A. nemorensis
# 1) Keep only genes with valid Arabidopsis IDs for org.At.tair.db
combined_species_df_sag_nem_surv_wilt_df <- combined_species_df_sag_nem_surv_wilt[!is.na(combined_species_df_sag_nem_surv_wilt$At), ]

# 2) Define genes of interest: green AND quadrant in Q1 or Q8
genes_of_interest_sag_nem_surv_wilt_df <- combined_species_df_sag_nem_surv_wilt_df$color == "red" & combined_species_df_sag_nem_surv_wilt_df$quadrant %in% c("Q2")

# 3) Build geneList for topGO (universe = all genes in df)
allGenes_numeric <- as.integer(genes_of_interest_sag_nem_surv_wilt_df)
names(allGenes_numeric) <- df$At

# (optional but good) remove duplicated TAIR IDs if they exist
allGenes_numeric <- allGenes_numeric[!duplicated(names(allGenes_numeric))]

# 4) Create topGO object
tGOdata_sag_nem_surv_wilt_df <- new("topGOdata",
                                    description = "Green genes in Q2",
                                    ontology = "BP",
                                    allGenes = allGenes_numeric,
                                    geneSel = function(x) x == 1,
                                    nodeSize = 5,
                                    mapping = "org.At.tair.db",
                                    annot = annFUN.org)

# 5) Run enrichment (Fisher + elim)
result_fisher_sag_nem_surv_wilt_df <- runTest(tGOdata_sag_nem_surv_wilt_df, algorithm = "elim", statistic = "fisher")

# 6) Output table
goTable_recovery_stress_up_sag_down_nem_df <- GenTable(tGOdata_sag_nem_surv_wilt_df,
                                                       fisher = result_fisher_sag_nem_surv_wilt_df,
                                                       orderBy = "fisher",
                                                       topNodes = 20)

goTable_recovery_stress_up_sag_down_nem_df



##### A. nem and A. sag in recovery vs stress --- response up regulation in A. nemorensis, and down in A. sagittata
# 1) Keep only genes with valid Arabidopsis IDs for org.At.tair.db
combined_species_df_sag_nem_surv_wilt_df <- combined_species_df_sag_nem_surv_wilt[!is.na(combined_species_df_sag_nem_surv_wilt$At), ]

# 2) Define genes of interest: green AND quadrant in Q1 or Q8
genes_of_interest_sag_nem_surv_wilt_df <- combined_species_df_sag_nem_surv_wilt_df$color == "red" & combined_species_df_sag_nem_surv_wilt_df$quadrant %in% c("Q5")

# 3) Build geneList for topGO (universe = all genes in df)
allGenes_numeric <- as.integer(genes_of_interest_sag_nem_surv_wilt_df)
names(allGenes_numeric) <- df$At

# (optional but good) remove duplicated TAIR IDs if they exist
allGenes_numeric <- allGenes_numeric[!duplicated(names(allGenes_numeric))]

# 4) Create topGO object
tGOdata_sag_nem_surv_wilt_df <- new("topGOdata",
                                    description = "Green genes in Q2",
                                    ontology = "BP",
                                    allGenes = allGenes_numeric,
                                    geneSel = function(x) x == 1,
                                    nodeSize = 5,
                                    mapping = "org.At.tair.db",
                                    annot = annFUN.org)

# 5) Run enrichment (Fisher + elim)
result_fisher_sag_nem_surv_wilt_df <- runTest(tGOdata_sag_nem_surv_wilt_df, algorithm = "elim", statistic = "fisher")

# 6) Output table
goTable_recovery_stress_up_nem_down_sag_df <- GenTable(tGOdata_sag_nem_surv_wilt_df,
                                                       fisher = result_fisher_sag_nem_surv_wilt_df,
                                                       orderBy = "fisher",
                                                       topNodes = 20)

goTable_recovery_stress_up_nem_down_sag_df



##### A. nem and A. sag in recovery vs stress --- response more down regulation in A. nemorensis, and less down in A. sagittata
# 1) Keep only genes with valid Arabidopsis IDs for org.At.tair.db
combined_species_df_sag_nem_surv_wilt_df <- combined_species_df_sag_nem_surv_wilt[!is.na(combined_species_df_sag_nem_surv_wilt$At), ]

# 2) Define genes of interest: green AND quadrant in Q1 or Q8
genes_of_interest_sag_nem_surv_wilt_df <- combined_species_df_sag_nem_surv_wilt_df$color == "red" & combined_species_df_sag_nem_surv_wilt_df$quadrant %in% c("Q7")

# 3) Build geneList for topGO (universe = all genes in df)
allGenes_numeric <- as.integer(genes_of_interest_sag_nem_surv_wilt_df)
names(allGenes_numeric) <- df$At

# (optional but good) remove duplicated TAIR IDs if they exist
allGenes_numeric <- allGenes_numeric[!duplicated(names(allGenes_numeric))]

# 4) Create topGO object
tGOdata_sag_nem_surv_wilt_df <- new("topGOdata",
                                    description = "Green genes in Q2",
                                    ontology = "BP",
                                    allGenes = allGenes_numeric,
                                    geneSel = function(x) x == 1,
                                    nodeSize = 5,
                                    mapping = "org.At.tair.db",
                                    annot = annFUN.org)

# 5) Run enrichment (Fisher + elim)
result_fisher_sag_nem_surv_wilt_df <- runTest(tGOdata_sag_nem_surv_wilt_df, algorithm = "elim", statistic = "fisher")

# 6) Output table
goTable_recovery_stress_more_down_nem_df <- GenTable(tGOdata_sag_nem_surv_wilt_df,
                                                     fisher = result_fisher_sag_nem_surv_wilt_df,
                                                     orderBy = "fisher",
                                                     topNodes = 20)

goTable_recovery_stress_more_down_nem_df



##### A. nem and A. sag in recovery vs stress --- response more down regulation in A. sagittata, and less down in A. nemorensis
# 1) Keep only genes with valid Arabidopsis IDs for org.At.tair.db
combined_species_df_sag_nem_surv_wilt_df <- combined_species_df_sag_nem_surv_wilt[!is.na(combined_species_df_sag_nem_surv_wilt$At), ]

# 2) Define genes of interest: green AND quadrant in Q1 or Q8
genes_of_interest_sag_nem_surv_wilt_df <- combined_species_df_sag_nem_surv_wilt_df$color == "red" & combined_species_df_sag_nem_surv_wilt_df$quadrant %in% c("Q4")

# 3) Build geneList for topGO (universe = all genes in df)
allGenes_numeric <- as.integer(genes_of_interest_sag_nem_surv_wilt_df)
names(allGenes_numeric) <- df$At

# (optional but good) remove duplicated TAIR IDs if they exist
allGenes_numeric <- allGenes_numeric[!duplicated(names(allGenes_numeric))]

# 4) Create topGO object
tGOdata_sag_nem_surv_wilt_df <- new("topGOdata",
                                    description = "Green genes in Q2",
                                    ontology = "BP",
                                    allGenes = allGenes_numeric,
                                    geneSel = function(x) x == 1,
                                    nodeSize = 5,
                                    mapping = "org.At.tair.db",
                                    annot = annFUN.org)

# 5) Run enrichment (Fisher + elim)
result_fisher_sag_nem_surv_wilt_df <- runTest(tGOdata_sag_nem_surv_wilt_df, algorithm = "elim", statistic = "fisher")

# 6) Output table
goTable_recovery_stress_more_down_sag_df <- GenTable(tGOdata_sag_nem_surv_wilt_df,
                                                     fisher = result_fisher_sag_nem_surv_wilt_df,
                                                     orderBy = "fisher",
                                                     topNodes = 20)

goTable_recovery_stress_more_down_sag_df




#write.csv(goTable_recovery_stress_more_down_sag_df, "goTable_recovery_stress_more_down_sag_df.csv")
#write.csv(goTable_recovery_stress_more_down_nem_df,  "goTable_recovery_stress_more_down_nem_df.csv")
#write.csv(goTable_recovery_stress_more_up_sag_df, "goTable_recovery_stress_more_up_sag_df.csv")
#write.csv(goTable_recovery_stress_more_up_nem_df, "goTable_recovery_stress_more_up_nem_df.csv")
#write.csv(goTable_recovery_stress_more_down_sag_df, "goTable_recovery_stress_more_down_sag_df.csv")
#write.csv(goTable_recovery_stress_more_down_sag_df, "goTable_recovery_stress_more_down_sag_df.csv")




##### Ending on 6th January 2026, here is the end of the new analysis #############
#########################################################################################################

####### 16th February 2026, Task: Change in genes in each quadrant from Stress to Recover #######

library(dplyr)
library(ggplot2)

# ----------------------------
# 0) Pick quadrant levels
# ----------------------------
# Full (8) quadrants (Q3/Q6 will just be 0 if absent)
q_levels_8 <- paste0("Q", 1:8)

# If you truly want only the 6 quadrants you see (remove Q3 & Q6):
q_levels_6 <- c("Q1","Q2","Q4","Q5","Q7","Q8")

use_levels <- q_levels_6   # <-- change to q_levels_6 if you want 6x6 matrix


# ----------------------------
# 1) Join stress + recovery by gene
# ----------------------------
# Use gene_id (or At). gene_id is safest if unique across your data.
trans_df <- combined_species_df_sag_nem %>%
  select(gene_id, quadrant_stress = quadrant,
         log2FC_sag_stress = log2FC_sag, log2FC_nem_stress = log2FC_nem,
         padj_sag_stress = padj_sag, padj_nem_stress = padj_nem, padj_gxe_stress = padj_gxe) %>%
  inner_join(
    combined_species_df_sag_nem_surv %>%
      select(gene_id, quadrant_recovery = quadrant,
             log2FC_sag_rec = log2FC_sag, log2FC_nem_rec = log2FC_nem,
             padj_sag_rec = padj_sag, padj_nem_rec = padj_nem, padj_gxe_rec = padj_gxe),
    by = "gene_id"
  ) %>%
  mutate(
    quadrant_stress   = factor(quadrant_stress, levels = use_levels),
    quadrant_recovery = factor(quadrant_recovery, levels = use_levels)
  )

# quick check
table(trans_df$quadrant_stress, useNA = "ifany")
table(trans_df$quadrant_recovery, useNA = "ifany")


# ----------------------------
# 2) Function to make transition matrices + heatmaps
# ----------------------------
make_transition <- function(df, title_prefix = "All genes") {
  
  # counts
  mat_counts <- table(df$quadrant_stress, df$quadrant_recovery)
  
  # row-wise probabilities: P(recovery quadrant | stress quadrant)
  mat_prop <- prop.table(mat_counts, margin = 1)
  
  # long format for ggplot
  counts_long <- as.data.frame(mat_counts) %>%
    rename(Q_stress = Var1, Q_recovery = Var2, n = Freq)
  
  prop_long <- as.data.frame(mat_prop) %>%
    rename(Q_stress = Var1, Q_recovery = Var2, p = Freq)
  
  # Heatmap (counts)
  p_counts <- ggplot(counts_long, aes(Q_recovery, Q_stress, fill = n)) +
    geom_tile(color = "white") +
    geom_text(aes(label = n), size = 3) +
    labs(
      title = paste0(title_prefix, " — Stress quadrant → Recovery quadrant (counts)"),
      x = "Recovery vs Control quadrant",
      y = "Stress vs Control quadrant"
    ) +
    theme_minimal()
  
  # Heatmap (row probabilities)
  p_prop <- ggplot(prop_long, aes(Q_recovery, Q_stress, fill = p)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", p)), size = 3) +
    scale_fill_gradient(low = "white", high = "saddlebrown") +
    labs(
      title = paste0(title_prefix, " - Stress → Recovery"),
      x = "Recovery vs Control",
      y = "Stress vs Control", fill = "Proportion"
    ) +
    theme_minimal()
  
  list(
    counts_matrix = mat_counts,
    prop_matrix   = mat_prop,
    plot_counts   = p_counts,
    plot_prop     = p_prop
  )
}


# ----------------------------
# 3) A) Matrix for ALL genes (universe)
# ----------------------------
res_all <- make_transition(trans_df, "All genes")
res_all$counts_matrix
res_all$prop_matrix
res_all$plot_counts
res_all$plot_prop

ggsave(
  filename = "All_genes_transition_plot.pdf",
  plot = res_all$plot_prop,
  width = 7,
  height = 6
)




# ----------------------------
# 4) B) “Independently” for sagittata: subset genes significant in sagittata
# ----------------------------
# choose which contrast’s significance defines “sagittata genes”
# Here I use stress OR recovery significant in SAG (you can tighten this if needed).
trans_df_sag <- trans_df %>%
  filter(padj_sag_stress < 0.05 | padj_sag_rec < 0.05)

res_sag <- make_transition(trans_df_sag, "A_sagittata-significant genes")
res_sag$counts_matrix
res_sag$plot_prop

ggsave(
  filename = "A_sagittata_transition_plot.pdf",
  plot = res_sag$plot_prop,
  width = 7,
  height = 6
)


# ----------------------------
# 5) C) “Independently” for nemorensis: subset genes significant in nemorensis
# ----------------------------
trans_df_nem <- trans_df %>%
  filter(padj_nem_stress < 0.05 | padj_nem_rec < 0.05)

res_nem <- make_transition(trans_df_nem, "A_nemorensis-significant genes")
res_nem$counts_matrix
res_nem$plot_prop

ggsave(
  filename = "A_nemorensis_transition_plot.pdf",
  plot = res_nem$plot_prop,
  width = 7,
  height = 6
)

# ----------------------------
# 6) Optional: focus only on GxE genes (interaction)
# ----------------------------
trans_df_gxe <- trans_df %>%
  filter(padj_gxe_stress < 0.001 | padj_gxe_rec < 0.001)

res_gxe <- make_transition(trans_df_gxe, "GxE genes")
res_gxe$plot_prop

######## End on 16th February 2026 ########################

########### Start on 20th February  ---- Task: Slightly changed, added Quadrant Q3 and Q6 (for genes with no significant interaction)  #######
############Change in genes in each quadrant from Stress to Recover #######

# ----------------------------
# 0) Quadrant levels (8)
# ----------------------------
q_levels_8 <- paste0("Q", 1:8)
use_levels <- q_levels_8


# ----------------------------
# Helper: build "quadrant_group" based on color rules
# ----------------------------
# Rules:
# - Q3: color green AND log2FC_sag > 0
# - Q6: color green AND log2FC_sag < 0
# - Else: keep original quadrant ONLY for red genes
# - Everything else becomes NA (excluded later)
make_quadrant_group <- function(df, quad_col = "quadrant", color_col = "color", fc_col = "log2FC_sag") {
  
  df %>%
    mutate(
      quadrant_group = case_when(
        .data[[color_col]] == "green" & .data[[fc_col]] > 0 ~ "Q3",
        .data[[color_col]] == "green" & .data[[fc_col]] < 0 ~ "Q6",
        .data[[color_col]] == "red" ~ as.character(.data[[quad_col]]),
        TRUE ~ NA_character_
      ),
      quadrant_group = factor(quadrant_group, levels = use_levels)
    )
}


# ----------------------------
# 1) Prepare stress and recovery data with quadrant_group
# ----------------------------
stress_df <- combined_species_df_sag_nem %>%
  make_quadrant_group(quad_col = "quadrant", color_col = "color", fc_col = "log2FC_sag") %>%
  select(
    gene_id,
    quadrant_stress = quadrant_group,
    log2FC_sag_stress = log2FC_sag, log2FC_nem_stress = log2FC_nem,
    padj_sag_stress = padj_sag, padj_nem_stress = padj_nem, padj_gxe_stress = padj_gxe
  )

recovery_df <- combined_species_df_sag_nem_surv %>%
  make_quadrant_group(quad_col = "quadrant", color_col = "color", fc_col = "log2FC_sag") %>%
  select(
    gene_id,
    quadrant_recovery = quadrant_group,
    log2FC_sag_rec = log2FC_sag, log2FC_nem_rec = log2FC_nem,
    padj_sag_rec = padj_sag, padj_nem_rec = padj_nem, padj_gxe_rec = padj_gxe
  )

# ----------------------------
# 2) Join stress + recovery and drop NA quadrant groups
# ----------------------------
trans_df <- stress_df %>%
  inner_join(recovery_df, by = "gene_id") %>%
  filter(!is.na(quadrant_stress), !is.na(quadrant_recovery))

# quick check
table(trans_df$quadrant_stress, useNA = "ifany")
table(trans_df$quadrant_recovery, useNA = "ifany")


# ----------------------------
# 3) Function to make transition matrices + heatmaps
# ----------------------------
make_transition <- function(df, title_prefix = "All genes") {
  
  mat_counts <- table(df$quadrant_stress, df$quadrant_recovery)
  mat_prop   <- prop.table(mat_counts, margin = 1)
  
  counts_long <- as.data.frame(mat_counts) %>%
    rename(Q_stress = Var1, Q_recovery = Var2, n = Freq)
  
  prop_long <- as.data.frame(mat_prop) %>%
    rename(Q_stress = Var1, Q_recovery = Var2, p = Freq)
  
  p_counts <- ggplot(counts_long, aes(Q_recovery, Q_stress, fill = n)) +
    geom_tile(color = "white") +
    geom_text(aes(label = n), size = 3) +
    labs(
      title = paste0(title_prefix, " — Stress quadrant → Recovery quadrant (counts)"),
      x = "Recovery vs Control quadrant",
      y = "Stress vs Control quadrant"
    ) +
    theme_minimal()
  
  p_prop <- ggplot(prop_long, aes(Q_recovery, Q_stress, fill = p)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", p)), size = 3) +
    scale_fill_gradient(low = "white", high = "saddlebrown") +
    labs(
      title = paste0(title_prefix, " - Stress → Recovery"),
      x = "Recovery vs Control",
      y = "Stress vs Control",
      fill = "Proportion"
    ) +
    theme_minimal()
  
  list(
    counts_matrix = mat_counts,
    prop_matrix   = mat_prop,
    plot_counts   = p_counts,
    plot_prop     = p_prop
  )
}

# ----------------------------
# 4) Run transitions
# ----------------------------
res_all <- make_transition(trans_df, "Quadrant groups (red + green-split Q3/Q6)")
res_all$counts_matrix
res_all$plot_prop




# ----------------------------
# 4) B) “Independently” for sagittata: subset genes significant in sagittata
# ----------------------------
# choose which contrast’s significance defines “sagittata genes”
# Here I use stress OR recovery significant in SAG (you can tighten this if needed).
trans_df_sag <- trans_df %>%
  filter(padj_sag_stress < 0.05 | padj_sag_rec < 0.05)

res_sag <- make_transition(trans_df_sag, "A_sagittata-significant genes")
res_sag$counts_matrix
res_sag$plot_prop

ggsave(
  filename = "A_sagittata_transition_plot.pdf",
  plot = res_sag$plot_prop,
  width = 7,
  height = 6
)


# ----------------------------
# 5) C) “Independently” for nemorensis: subset genes significant in nemorensis
# ----------------------------
trans_df_nem <- trans_df %>%
  filter(padj_nem_stress < 0.05 | padj_nem_rec < 0.05)

res_nem <- make_transition(trans_df_nem, "A_nemorensis-significant genes")
res_nem$counts_matrix
res_nem$plot_prop

ggsave(
  filename = "A_nemorensis_transition_plot.pdf",
  plot = res_nem$plot_prop,
  width = 7,
  height = 6
)


################ End on 20th Febuary 2026 ###########################



