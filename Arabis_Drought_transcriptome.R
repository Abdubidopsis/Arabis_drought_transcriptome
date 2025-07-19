###--Differential expression analysis of HTSeqcounts using DESeq2--#######################
###' Abdul Saboor Khan (abdul.suboor123@yahoo.com) -- 16.Jan.2023
###' This script runs for Arabis sagittata and nemorensis counts files to find differential expression
###' first, get the counts file from the mapping (sam) files by htseq-counts in bash

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

############################################################################################################################################################################################################
#################################################################################################################################################################################################################
#################################################################################################################################################################################################################
###############                         #########################################################################################################################################################################
############### all sampples.           #########################################################################################################################################################################
###############                         #########################################################################################################################################################################
#################################################################################################################################################################################################################
#################################################################################################################################################################################################################
#############################################################################################################################################################################################################################################
#############################################################################################################################################################################################################################################
#############################################################################################################################################################################################################################################

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


####ref "nemorensiswilting" ####
dds$group<- relevel(dds$group, ref = "nemorensiswilting")
levels(dds$group)
dds_test <- DESeq(dds)
resultsNames(dds_test)

sagwilt_nemwilt <- results(dds_test, name = "group_sagitattawilting_vs_nemorensiswilting")
#write.csv(sagwilt_nemwilt, "script_for_juliette/new_go_output_from_deseqdrought_all/sagwilt_nemwilt.csv")

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
  geom_point(size = 5) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() + theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "gray", linewidth = 0.2),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.key.size = unit(2, "lines"),  # Adjust the size of the legend
    legend.text = element_text(size = 12),  # Adjust the font size of legend text
    legend.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 14), axis.title = element_text(size = 14)  # Adjust the font size and style of the legend title
  )


x 

#ggsave("script_for_juliette/new_go_output_from_deseqdrought_all/pca_all_samples.png", plot = x, width = 9, height = 7, dpi = 300)
#ggsave("script_for_juliette/new_go_output_from_deseqdrought_all/pca_all_samples.pdf", plot = x, width = 9, height = 7, dpi = 300)



#### Random forest ----- Get expression matrix#############
#vsd <- vst(dds, blind=FALSE)
#expr <- t(assay(vsd))  # samples x genes

# Make sure it's all numeric and well-formatted
#expr_df <- as.data.frame(expr)
#expr_df$genotype <- as.factor(colData(vsd)$genotype)  # label column

# Check for NAs or non-numeric columns (just in case)
#stopifnot(all(sapply(expr_df[, -ncol(expr_df)], is.numeric)))

# Run Random Forest
#rf_model <- randomForest(species ~ ., data = expr_df, importance = TRUE, ntree = 500)

# 3. View importance of genes
#varImpPlot(rf_model)
#important_genes <- importance(rf_model)


library(xgboost)
library(tidyverse)
library(e1071)              # Load it in your script or session
library(nnet)

# Format data
#expr <- t(normalized_counts)  # Samples as rows
#labels <- colData(dds)$genotype  # Factor label e.g., "sag", "nem"
#labels_numeric <- as.numeric(as.factor(labels)) - 1  # Convert to 0/1

#dtrain <- xgb.DMatrix(data = expr, label = labels_numeric)

# Train model
#xgb_model <- xgboost(data = dtrain, nrounds = 50, objective = "binary:logistic", verbose = 0)

# Feature importance
#importance <- xgb.importance(model = xgb_model)
#xgb.plot.importance(importance)


# Data prep
df <- as.data.frame(t(normalized_counts))
df$label <- colData(dds)$genotype

# Train/test split
set.seed(123)
train_idx <- sample(1:nrow(df), size = 0.7 * nrow(df))
train <- df[train_idx, ]
test <- df[-train_idx, ]

# SVM training
svm_model <- svm(label ~ ., data = train, kernel = "linear")

# Predict & accuracy
pred <- predict(svm_model, test)
table(pred, test$label)

# Simulate multi-condition timepoint label
multi_label <- paste(colData(dds)$genotype, colData(dds)$timepoint, sep = "_")
df <- as.data.frame(t(normalized_counts))
df$label <- as.factor(multi_label)

# Fit model
multinom_model <- multinom(label ~ ., data = df)
summary(multinom_model)


# 1. Calculate variance of each gene
var_genes <- apply(normalized_counts, 1, var)

# 2. Select top 100 most variable genes
top_genes <- names(sort(var_genes, decreasing = TRUE))[1:100]

# 3. Subset your data
df <- as.data.frame(t(normalized_counts[top_genes, ]))
df$label <- as.factor(paste(colData(dds)$genotype, colData(dds)$timepoint, sep = "_"))
multinom_model <- multinom(label ~ ., data = df)
summary(multinom_model)



# Boruta requires data.frame input
df <- as.data.frame(t(normalized_counts))
df$label <- as.factor(colData(dds)$genotype)

# Run Boruta
set.seed(123)
boruta_result <- Boruta(label ~ ., data = df, doTrace = 2)

# Final selection
final <- TentativeRoughFix(boruta_result)
plot(final)
getSelectedAttributes(final)


library(survival)
library(randomForestSRC)

# Survival data: assuming you have 'days_survived' and 'status'
surv_data <- data.frame(t(normalized_counts))
surv_data$days <- survival
surv_data$status <- your_status_vector  # 1 = event (death), 0 = censored

# Fit Cox model
cox_model <- coxph(Surv(days, status) ~ ., data = surv_data)
summary(cox_model)



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
#write.csv(combined_species_df_sag_nem_wilt_ctrl, "combined_species_df_sag_nem_wilt_ctrl_quadrants.csv")
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

# Remove genes without an orthologue
combined_species_df_sag_nem <- combined_species_df_sag_nem %>%
  filter(!is.na(At))


# Merge expression data with orthologue mapping
combined_species_df_sag_nem_surv <- combined_species_df_sag_nem_surv %>%
  left_join(new_dataframe1, by = c("gene_id" = "arabis_cleaned"))

# Remove genes without an orthologue
combined_species_df_sag_nem_surv <- combined_species_df_sag_nem_surv %>%
  filter(!is.na(At))

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
write.csv(goEnrichmentQ_shared_sag_nem_left_wilt_ctrl, "goEnrichmentQ_shared_sag_nem_left_wilt_ctrl.csv")



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

goEnrichmentQ_shared_sag_nem_right_wilt_ctrl <- GenTable(tGOdata_sag_nem_right_wilt_ctrl, fisher=results.fisher_sag_nem_right_wilt_ctrl, orderBy="fisher", topNodes=50)
write.csv(goEnrichmentQ_shared_sag_nem_right_wilt_ctrl, "goEnrichmentQ_shared_sag_nem_right_wilt_ctrl.csv")



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


write.csv(goEnrichmentQ_sag_up_more_against_nem_wiltQ1, "goEnrichmentQ_sag_up_more_against_nem_wiltQ1.csv")
write.csv(goEnrichmentQ_sag_down_more_against_nem_wiltQ4, "goEnrichmentQ_sag_down_more_against_nem_wiltQ4.csv")
write.csv(goEnrichmentQ_sag_up_and_nem_down_wiltQ2, "goEnrichmentQ_sag_up_and_nem_down_wiltQ2.csv")
write.csv(goEnrichmentQ_nem_up_more_against_sag_wiltQ8, "goEnrichmentQ_nem_up_more_against_sag_wiltQ8.csv")
write.csv(goEnrichmentQ_nem_up_and_sag_down_wiltQ5, "goEnrichmentQ_nem_up_and_sag_down_wiltQ5.csv")
write.csv(goEnrichmentQ_nem_down_more_against_sag_wiltQ7, "goEnrichmentQ_nem_down_more_against_sag_wiltQ7.csv")

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
write.csv(goEnrichmentQ_shared_sag_nem_left_surv_ctrl, "goEnrichmentQ_shared_sag_nem_left_surv_ctrl.csv")



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
write.csv(goEnrichmentQ_shared_sag_nem_right_surv_ctrl, "goEnrichmentQ_shared_sag_nem_right_surv_ctrl.csv")


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


write.csv(goEnrichmentQ_sag_up_more_against_nem_rec_Q1, "goEnrichmentQ_sag_up_more_against_nem_recQ1.csv")
write.csv(goEnrichmentQ_sag_down_more_against_nem_rec_Q4, "goEnrichmentQ_sag_down_more_against_nem_recQ4.csv")
write.csv(goEnrichmentQ_sag_up_and_nem_down_rec_Q2, "goEnrichmentQ_sag_up_and_nem_down_rec_Q2.csv")
write.csv(goEnrichmentQ_nem_up_and_sag_down_rec_Q5, "goEnrichmentQ_nem_up_and_sag_down_rec_Q5.csv")
write.csv(goEnrichmentQ_nem_down_more_against_sag_rec_Q7, "goEnrichmentQ_nem_down_more_against_sag_rec_Q7.csv")
write.csv(goEnrichmentQ_nem_up_more_against_sag_rec_Q8, "goEnrichmentQ_nem_up_more_against_sag_rec_Q8.csv")

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
