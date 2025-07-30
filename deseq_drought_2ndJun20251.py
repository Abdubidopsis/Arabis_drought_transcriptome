#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Abdul Saboor Khan, PhD - University of Cologne in April 2023 in R and then Coverted to Python on 2nd June 2025.
Script for DESeq2-based transcriptome analysis of Arabis species, and Gene ontology analysis and GO term enrichment visualization.
Conversion of the provided R DESeq2-based transcriptome analysis pipeline into Python,
leveraging rpy2 for DESeq2 functions and Python libraries (pandas, numpy, sklearn, seaborn, matplotlib)
for downstream analyses (filtering, k-means clustering, heatmaps, volcano plots, PCA, and ML models).
"""

import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from xgboost import XGBClassifier, plot_importance
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from boruta import BorutaPy
import lifelines
from lifelines import CoxPHFitter

# Import rpy2 to access R and DESeq2
import rpy2.robjects as ro
from rpy2.robjects import default_converter, Formula
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri

# Load R libraries (make sure these packages are installed in your R environment)
base     = importr('base')
utils    = importr('utils')
stats    = importr('stats')
DESeq2   = importr('DESeq2')
limma    = importr('limma')
pheatmap = importr('pheatmap')

# Note: Instead of pandas2ri.activate(), we will wrap any pandas<->R conversion
#       in a `with localconverter(default_converter + pandas2ri.converter):` block.
#       For example:
#
#       with localconverter(default_converter + pandas2ri.converter):
#           r_df = ro.conversion.py2rpy(pandas_df)
#
#       and vice versa:
#
#       with localconverter(default_converter + pandas2ri.converter):
#           pandas_df = ro.conversion.rpy2py(r_df)

# The following R packages are needed for machine learning and survival in R context, 
# but we'll do most modeling in Python. We still import them in case we want R-based transformations.
# importr('randomForest')
# importr('xgboost')
# importr('e1071')
# importr('nnet')
# importr('Boruta')
# importr('survival')
# importr('randomForestSRC')


# =============================================================================
# PART I: ALL SAMPLES ANALYSIS (22 samples: 11 nemorensis + 11 sagitatta)
# =============================================================================

# 1. Set working directory in Python
wd = "/Users/Shared/Files From d.localized/PhD-Uni-Koeln_2021-2024/PhD_work/seeds data/" \
     "2nd Experiment 2022/transcriptome/Arabis_RNA_raw_data/90-774066047/" \
     "anaylsis_drought_mRNA_both_species_with_new_two_genomes/" \
     "analysis_both_species_with_nem_genomes/script_for_juliette/"
os.chdir(wd)

# 2. List all HTSeq-count ".txt" files in directory
sample_files = sorted([os.path.basename(f) for f in glob.glob(os.path.join(wd, "*.txt"))])

# 3. Define sample metadata: condition and genotype vectors
condition = [
    "control", "control", "control", "control",
    "survival", "survival", "survival",
    "wilting", "wilting", "wilting", "wilting",
    "control", "control", "control", "control",
    "survival", "survival", "survival",
    "wilting", "wilting", "wilting", "wilting"
]
genotype = [
    "nemorensis", "nemorensis", "nemorensis", "nemorensis",
    "nemorensis", "nemorensis", "nemorensis",
    "nemorensis", "nemorensis", "nemorensis", "nemorensis",
    "sagitatta", "sagitatta", "sagitatta", "sagitatta",
    "sagitatta", "sagitatta", "sagitatta",
    "sagitatta", "sagitatta", "sagitatta", "sagitatta"
]

if len(sample_files) != len(condition) or len(sample_files) != len(genotype):
    raise ValueError("Length of sample_files, condition, and genotype must match (22).")

# 4. Build a pandas DataFrame for sampleTable (R's data.frame)
sample_table = pd.DataFrame({
    'sampleName': sample_files,
    'fileName': sample_files,
    'condition': pd.Categorical(condition, categories=["control", "survival", "wilting"]),
    'genotype': pd.Categorical(genotype, categories=["nemorensis", "sagitatta"])
})

# Ensure all columns are of standard types
# (Make sure your pandas columns are plain Python strings or pandas.Categorical,
#  not some exotic dtype that rpy2 can’t handle directly.)

sample_table['sampleName'] = sample_table['sampleName'].astype(str)
sample_table['fileName']   = sample_table['fileName'].astype(str)
sample_table['condition']  = sample_table['condition'].astype('category')
sample_table['genotype']   = sample_table['genotype'].astype('category')

# Now convert to R data.frame
with localconverter(default_converter + pandas2ri.converter):
    r_sample_table = ro.conversion.py2rpy(sample_table)


# 6. Create a DESeqDataSetFromHTSeqCount via rpy2
#    design = ~ genotype + condition + genotype:condition
design_formula = Formula('~ genotype + condition + genotype:condition')
ddsHTSeq = DESeq2.DESeqDataSetFromHTSeqCount(
    sampleTable=r_sample_table,
    directory=wd,
    design=design_formula
)


# Call the R generic `counts(ddsHTSeq)` via ro.r('counts'), then convert it to a pandas DataFrame
count_matrix_R = ro.r('counts')(ddsHTSeq)

# 1. Call the R generic 'counts' to get an R matrix of raw counts:
count_matrix_R = ro.r('counts')(ddsHTSeq)

# 2. Convert the R matrix to a NumPy array:
count_arr = np.array(count_matrix_R)

# 3. Extract row and column names from the R matrix:
rownames_R = ro.r('rownames')(count_matrix_R)
colnames_R = ro.r('colnames')(count_matrix_R)

#    Convert those to Python lists of strings:
row_index = list(map(str, list(rownames_R)))
col_columns = list(map(str, list(colnames_R)))

# 4. Build a pandas DataFrame from the NumPy array, using those names:
counts_df = pd.DataFrame(count_arr, index=row_index, columns=col_columns)

# 5. If you truly only need a bare NumPy array, you can now do:
count_matrix = counts_df.values
print(counts_df.head(10))


# But in most workflows, it's easier to keep 'counts_df' as a DataFrame,
# since it still has genes as the index and sample names as the columns.



# (Optionally, you can keep counts_df and use its index/columns rather than a bare array.)

# 7. Filter out genes with low counts (average over samples > 100)
#    In R: ddsHTSeq[rowSums(counts(ddsHTSeq)) / num_samples > count_threshold, ]
num_samples = 22
count_threshold = 100

# Retrieve count matrix (R object) and convert to numpy
#count_matrix = np.array(DESeq2.counts(ddsHTSeq))
# count_matrix shape: genes x samples
# Compute row means
row_means = count_matrix.sum(axis=1) / float(num_samples)
keep_idx = np.where(row_means > count_threshold)[0]

# Subset the R DESeqDataSet object: there is no direct R subsetting via rpy2 indexing,
# so we call the R function 'subset' or use base indexing. E.g. ddsHTSeq[keep, ]
# In rpy2, one can do:
r_keep = ro.BoolVector([i in keep_idx for i in range(count_matrix.shape[0])])
#ddsHTSeq_filtered = ddsHTSeq.rx[r_keep, True]  # keep rows for which r_keep is TRUE; True means keep all columns

# 1) Put the DESeqDataSet into R’s global environment under a name, say “dds”
ro.globalenv['dds'] = ddsHTSeq

# 2) Also pass your Boolean “keep” vector into R:
ro.globalenv['keep_vec'] = r_keep

# 3) Instruct R to subset row-wise: dds_filtered <- dds[keep_vec, ]
ro.r('dds_filtered <- dds[keep_vec, ]')

# 4) Pull the filtered DESeqDataSet back into Python
ddsHTSeq_filtered = ro.globalenv['dds_filtered']


# 8. Run DESeq: (design already set)
dds = DESeq2.DESeq(ddsHTSeq_filtered, fitType="mean")

# 9. Extract results for all coefficients (no contrast specified)
res_dds_all = DESeq2.results(dds)

# 2) In R, coerce that S4 “DESeqResults” into a data.frame
# 1) Coerce DESeqResults into an R data.frame
res_dds_all_Rdf = ro.r('as.data.frame')(res_dds_all)

# 2) Extract the row names from that R data.frame
r_rownames = ro.r('rownames')(res_dds_all_Rdf)         # an R character vector
py_rownames = list(map(str, list(r_rownames)))         # convert to Python list of strings

# 3) Convert the R data.frame to a pandas DataFrame
from rpy2.robjects import default_converter
from rpy2.robjects.conversion import localconverter

with localconverter(default_converter + pandas2ri.converter):
    res_dds_all_df = ro.conversion.rpy2py(res_dds_all_Rdf)

# 4) Now assign the Python list of row names as the DataFrame index
res_dds_all_df.index = py_rownames

# Verify
print(res_dds_all_df.head())
# Count the number of rows in res_dds_all_df
num_rows = res_dds_all_df.shape[0]
num_rows
num_rows = res_dds_all_df.shape[0]
print(num_rows)
len(res_dds_all_df)
# Optionally save the results DataFrame to CSV
# res_dds_all_df.to_csv("all_samples_res_dds_9315.csv")

# =============================================================================
# Filtering genes with low counts (size factor estimation & normalization)
# =============================================================================

# In R: dds_lowcount <- estimateSizeFactors(dds)
# 1) Suppose you have already run:
#    dds = DESeq2.DESeq(ddsHTSeq_filtered, fitType="mean")

# 2) Call the R function estimateSizeFactors on that dds object:
dds_lowcount = ro.r('estimateSizeFactors')(dds)

# 3) If you need to fetch the size factors themselves (just for inspection):
size_factors_R = ro.r('sizeFactors')(dds_lowcount)
with localconverter(default_converter + pandas2ri.converter):
    size_factors = ro.conversion.rpy2py(size_factors_R)
print(size_factors)

# Get normalized counts: counts(dds_lowcount, normalized=TRUE)
#norm_counts_R = DESeq2.counts(dds_lowcount, normalized=True)  # returns R matrix
# 1) After you have dds_lowcount (the DESeqDataSet with size factors)
#    call the R generic 'counts' with normalized=TRUE:
norm_counts_R = ro.r('counts')(dds_lowcount, True)

# 2) Convert the R matrix to a pandas DataFrame:
count_arr   = np.array(norm_counts_R)
rownames_R  = ro.r('rownames')(norm_counts_R)
colnames_R  = ro.r('colnames')(norm_counts_R)
row_index   = list(map(str, list(rownames_R)))
col_columns = list(map(str, list(colnames_R)))

normalized_counts_df = pd.DataFrame(count_arr, index=row_index, columns=col_columns)

# 3) You can now inspect the first few rows:
print(normalized_counts_df.head())



# Convert to pandas DataFrame: rows=genes, cols=samples
norm_counts = pandas2ri.rpy2py(norm_counts_R)
# Set row and column names
gene_ids = np.array(base.rownames(norm_counts_R))
sample_ids = np.array(base.colnames(norm_counts_R))
norm_counts = pd.DataFrame(norm_counts, index=gene_ids, columns=sample_ids)

# Optionally save normalized counts
# norm_counts.to_csv("all_samples_normalized_9315.csv")

# =============================================================================
# K-means clustering + Heatmap (all genes passing count filter)
# =============================================================================

# In R: norm_deg <- normalizeBetweenArrays(normalized_counts, method="scale")
# We replicate row-wise scaling (center to mean=0, unit-variance) in Python:
#    norm_counts: genes x samples => scale rows
scaler = StandardScaler(with_mean=True, with_std=True)
# StandardScaler operates on features (columns) by default. To normalize rows, we transpose, scale, then transpose back.
norm_deg_array = scaler.fit_transform(norm_counts.values.T).T
norm_deg = pd.DataFrame(norm_deg_array, index=norm_counts.index, columns=norm_counts.columns)

# K-means clustering on transposed data (samples x genes)
k = 4
kmeans = KMeans(n_clusters=k, random_state=42)
# Fit on samples: rows = samples, columns = genes
# So we take norm_deg.T
kmeans.fit(norm_deg.values.T)
cluster_labels = kmeans.labels_  # length = number of samples (22)

# Create a clustermap (heatmap with hierarchical clustering) using seaborn
# We annotate columns (samples) by their cluster label.
col_colors = pd.Series(cluster_labels, index=norm_deg.columns).map(
    dict(zip(range(k), sns.color_palette("tab10", k)))
)

g = sns.clustermap(
    norm_deg,
    z_score=0,                      # z-score along rows (genes)
    method='complete',
    metric='euclidean',
    cmap="vlag",
    figsize=(8, 8),
    col_colors=col_colors,
    yticklabels=False
)
plt.setp(g.ax_heatmap.get_yticklabels(), fontsize=0)  # hide row labels
g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha='right')
# ← Add this line:
plt.show()
# g.savefig("kmeans_all_samples_9315_hclust.png", dpi=300)
# g.savefig("kmeans_all_samples_9315_hclust.pdf", dpi=300)

# =============================================================================
# Filter significant DEGs (padj < 0.05), then K-means + Heatmap on DEGs
# =============================================================================

# Convert res_dds_all to pandas DataFrame if not already
res_dds_all_df = res_dds_all_df.rename_axis('gene').reset_index().set_index('gene')

# Subset to padj < 0.05
res_dds_significant_df = res_dds_all_df[res_dds_all_df['padj'] < 0.05]

sig_gene_ids = res_dds_significant_df.index.values

# Subset normalized counts to significant genes
norm_counts_sig = norm_counts.loc[sig_gene_ids]

# Row-wise scaling of significant-gene matrix
norm_deg_sig_array = scaler.fit_transform(norm_counts_sig.values.T).T
norm_deg_sig = pd.DataFrame(norm_deg_sig_array, index=norm_counts_sig.index, columns=norm_counts_sig.columns)

# K-means (k=4) on significant genes: clustering samples by their expression of DEGs
k_sig = 4
kmeans_sig = KMeans(n_clusters=k_sig, random_state=42)
kmeans_sig.fit(norm_deg_sig.values.T)
cluster_labels_sig = kmeans_sig.labels_

col_colors_sig = pd.Series(cluster_labels_sig, index=norm_deg_sig.columns).map(
    dict(zip(range(k_sig), sns.color_palette("tab10", k_sig)))
)

g_sig = sns.clustermap(
    norm_deg_sig,
    z_score=0,
    method='complete',
    metric='euclidean',
    cmap="vlag",
    figsize=(8, 8),
    col_colors=col_colors_sig,
    yticklabels=False
)
plt.setp(g_sig.ax_heatmap.get_yticklabels(), fontsize=0)
g_sig.ax_heatmap.set_xticklabels(g_sig.ax_heatmap.get_xticklabels(), rotation=90, ha='right')
plt.show()
# g_sig.savefig("kmeans_DEGs_all_samples_3526_hclust.png", dpi=300)
# g_sig.savefig("kmeans_DEGs_all_samples_3526_hclust.pdf", dpi=300)

# =============================================================================
# Plot cooks distances (boxplot) and set up group factor for contrasts
# =============================================================================

# In R: boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
# Extract cooks from dds
cooks_R = DESeq2.assays(dds).rx2("cooks")
cooks = pandas2ri.rpy2py(cooks_R)
# Log10 transform (add small epsilon to avoid log(0))
eps = 1e-8
log10_cooks = np.log10(cooks + eps)

plt.figure(figsize=(10, 4))
plt.boxplot(log10_cooks, whis=0)  # range=0 means whiskers at min/max
plt.xticks(rotation=90)
plt.ylabel("log10(Cooks distance)")
plt.tight_layout()
# plt.savefig("cooks_boxplot.png", dpi=300)

# Create a new 'group' factor: genotype + condition
# In R: dds$group <- factor(paste0(dds$genotype, dds$condition))
# Get colData from dds
col_data_R = DESeq2.colData(vsd if 'vsd' in locals() else dds)  # but at this point, no vsd; use dds
col_data = pandas2ri.rpy2py(col_data_R)
col_data = pd.DataFrame(col_data, index=sample_ids)
col_data['genotype'] = col_data['genotype'].astype(str)
col_data['condition'] = col_data['condition'].astype(str)
col_data['group'] = col_data['genotype'] + col_data['condition']

# Redefine design(dss) = ~ group
# In rpy2: design(dds) <- ~ group
dds = DESeq2.design(dds, Formula('~ group'))

# =============================================================================
# Generate all pairwise contrasts by releveling the group reference in R
# =============================================================================

# Convenience function to run DESeq and extract a named contrast result from dds
def get_contrast(dds_obj, ref_level, contrast_name):
    """
    Take an existing DESeqDataSet with a 'group' factor, relevel the reference category to ref_level,
    rerun DESeq, and extract results for contrast_name.
    Returns a pandas DataFrame of results (indexed by gene).
    """
    # Relevel the factor in R: dds$group <- relevel(dds$group, ref=ref_level)
    ro.r('')  # no-op to separate commands
    relevel_cmd = f'dds$group <- relevel(dds$group, ref = "{ref_level}")'
    ro.r(relevel_cmd)

    # Rerun DESeq on dds (in place)
    dds_test = DESeq2.DESeq(dds_obj)

    # Extract results for contrast_name
    res_R = DESeq2.results(dds_test, name=contrast_name)
    res_df = pandas2ri.rpy2py(res_R).set_index('rownames')
    return res_df

# List contrasts based on R script:
# 1) ref "sagitattawilting": contrast "group_nemorensiswilting_vs_sagitattawilting"
nemwilt_sagwilt_df = get_contrast(dds, "sagitattawilting", "group_nemorensiswilting_vs_sagitattawilting")
# nemwilt_sagwilt_df.to_csv("nemwilt_sagwilt.csv")

# 2) ref "nemorensiswilting": contrast "group_sagitattawilting_vs_nemorensiswilting"
sagwilt_nemwilt_df = get_contrast(dds, "nemorensiswilting", "group_sagitattawilting_vs_nemorensiswilting")
# sagwilt_nemwilt_df.to_csv("sagwilt_nemwilt.csv")

# 3) ref "nemorensiscontrol": contrast "group_sagitattacontrol_vs_nemorensiscontrol"
sagctrl_nemctrl_df = get_contrast(dds, "nemorensiscontrol", "group_sagitattacontrol_vs_nemorensiscontrol")
# sagctrl_nemctrl_df.to_csv("sagctrl_nemctrl.csv")

# 4) ref "sagitattacontrol": contrasts
#    (a) "group_nemorensiscontrol_vs_sagitattacontrol"
nemctrl_sagctrl_df = get_contrast(dds, "sagitattacontrol", "group_nemorensiscontrol_vs_sagitattacontrol")
#    (b) "group_sagitattawilting_vs_sagitattacontrol"
sagwilt_sagctrl_df = get_contrast(dds, "sagitattacontrol", "group_sagitattawilting_vs_sagitattacontrol")
#    (c) "group_sagitattasurvival_vs_sagitattacontrol"
sagsurv_sagctrl_df = get_contrast(dds, "sagitattacontrol", "group_sagitattasurvival_vs_sagitattacontrol")
# Save as desired:
# nemctrl_sagctrl_df.to_csv("nemctrl_sagctrl.csv")
# sagwilt_sagctrl_df.to_csv("sagwilt_sagctrl.csv")
# sagsurv_sagctrl_df.to_csv("sagsurv_sagctrl.csv")

# 5) ref "sagitattasurvival": "group_nemorensissurvival_vs_sagitattasurvival"
nemsurv_sagsurv_df = get_contrast(dds, "sagitattasurvival", "group_nemorensissurvival_vs_sagitattasurvival")
# nemsurv_sagsurv_df.to_csv("nemsurv_sagsurv.csv")

# 6) ref "nemorensissurvival": "group_sagitattasurvival_vs_nemorensissurvival"
sagsurv_nemsurv_df = get_contrast(dds, "nemorensissurvival", "group_sagitattasurvival_vs_nemorensissurvival")
# sagsurv_nemsurv_df.to_csv("sagsurv_nemsurv.csv")

# 7) ref "nemorensiscontrol" (again): contrasts
#    (a) "group_nemorensiswilting_vs_nemorensiscontrol"
nemwilt_nemctrl_df = get_contrast(dds, "nemorensiscontrol", "group_nemorensiswilting_vs_nemorensiscontrol")
#    (b) "group_nemorensissurvival_vs_nemorensiscontrol"
nemsurv_nemctrl_df = get_contrast(dds, "nemorensiscontrol", "group_nemorensissurvival_vs_nemorensiscontrol")
# nemwilt_nemctrl_df.to_csv("nemwilt_nemctrl.csv")
# nemsurv_nemctrl_df.to_csv("nemsurv_nemctrl.csv")

# =============================================================================
# Volcano plots for A. sagittata (wilting vs. control, survival vs. control)
# =============================================================================

def plot_volcano(res_df, title, output_prefix=None):
    """
    Creates a volcano plot given a DataFrame (indexed by gene) with columns:
    'log2FoldChange' and 'pvalue'. Saves figure with output_prefix if provided.
    """
    # Prepare data
    df_plot = pd.DataFrame({
        'gene': res_df.index.values,
        'log2FoldChange': res_df['log2FoldChange'],
        'log10pvalue': -np.log10(res_df['pvalue'] + 1e-300)  # avoid log10(0)
    }).dropna(subset=['log10pvalue'])

    # Define color bins
    my_breaks = np.array([-np.inf, -10, -5, -1, 0, 1, 5, 10, np.inf])
    my_labels = ["< -10", "-10 to -5", "-5 to -1", "-1 to 0",
                 "0 to 1", "1 to 5", "5 to 10", "> 10"]
    cmap = [
        "#762a83", "#67001f", "#b2182b", "#d6604d",
        "#f4a582", "#92c5de", "#4393c3", "#2166ac"
    ]

    # Categorize by log2FoldChange
    bin_indices = np.digitize(df_plot['log2FoldChange'], my_breaks) - 1
    df_plot['L2FC'] = [my_labels[i] for i in bin_indices]

    # Plot
    plt.figure(figsize=(9, 7))
    for lbl, color in zip(my_labels, cmap):
        sub = df_plot[df_plot['L2FC'] == lbl]
        plt.scatter(
            sub['log2FoldChange'], sub['log10pvalue'],
            c=color, label=lbl, s=15, alpha=0.7, edgecolors='none'
        )

    # Draw threshold lines
    plt.axvline(-1, linestyle='--', color='black')
    plt.axvline(1, linestyle='--', color='black')
    plt.axhline(-np.log10(0.05), linestyle='--', color='black')

    plt.xlabel("Log2 Fold Change", fontsize=14)
    plt.ylabel("-Log10 P-value", fontsize=14)
    plt.title(title, fontsize=16)
    plt.xlim(-10, 10)
    plt.ylim(0, 100)
    plt.legend(title="L2FC bins", bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()

    if output_prefix:
        plt.savefig(f"{output_prefix}.png", dpi=300)
        plt.savefig(f"{output_prefix}.pdf", dpi=300)
    plt.show()


# A. sagittata wilting vs. control
plot_volcano(
    sagwilt_sagctrl_df,
    title="Exp. regulation in wilting  - A. sagittata",
    output_prefix="volcano_sag_wilt_ctrl"
)

# A. sagittata survival vs. control
plot_volcano(
    sagsurv_sagctrl_df,
    title="Exp. regulation in survival  - A. sagittata",
    output_prefix="volcano_sag_surv_ctrl"
)

# =============================================================================
# Volcano plots for A. nemorensis (wilting vs. control, survival vs. control)
# =============================================================================

# A. nemorensis wilting vs. control
plot_volcano(
    nemwilt_nemctrl_df,
    title="Exp. regulation in wilting  - A. nemorensis",
    output_prefix="volcano_nem_wilt_ctrl"
)

# A. nemorensis survival vs. control
plot_volcano(
    nemsurv_nemctrl_df,
    title="Exp. regulation in survival  - A. nemorensis",
    output_prefix="volcano_nem_surv_ctrl"
)

# =============================================================================
# PCA analysis of all samples
# =============================================================================

# In R: vsd <- vst(dds, blind=FALSE)
vst = importr('DESeq2').vst  # vsd function in DESeq2

vsd_R = vst(dds, blind=False)
# Extract transformed counts: assay(vsd)
assay_R = DESeq2.assay(vsd_R)
vsd_matrix = pandas2ri.rpy2py(assay_R)  # genes x samples
vsd_df = pd.DataFrame(vsd_matrix, index=np.array(base.rownames(assay_R)),
                      columns=np.array(base.colnames(assay_R)))

# Transpose for PCA: samples x genes
vsd_t = vsd_df.T

# Run PCA in Python
pca = PCA(n_components=2)
pca_result = pca.fit_transform(vsd_t.values)
percent_var = np.round(100 * pca.explained_variance_ratio_)

pca_df = pd.DataFrame({
    'PC1': pca_result[:, 0],
    'PC2': pca_result[:, 1],
    'genotype': col_data['genotype'].values,
    'condition': col_data['condition'].values
}, index=vsd_t.index)

# Basic PCA scatterplot
plt.figure(figsize=(7, 6))
sns.scatterplot(
    data=pca_df, x='PC1', y='PC2',
    hue='condition', style='genotype', s=100, palette='Set1'
)
plt.xlabel(f"PC1: {percent_var[0]}% variance", fontsize=14)
plt.ylabel(f"PC2: {percent_var[1]}% variance", fontsize=14)
plt.gca().set_aspect('equal', 'box')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
# plt.savefig("pca_all_samples.png", dpi=300)
# plt.savefig("pca_all_samples.pdf", dpi=300)
plt.show()

# =============================================================================
# Machine Learning Models on normalized_counts
# =============================================================================

# Prepare expression matrix: samples x genes
expr = norm_counts.T  # shape: 22 samples x ~genes

# Add genotype labels (0/1) for Random Forest and XGBoost
labels = col_data['genotype'].values
labels_binary = (labels == "sagitatta").astype(int)  # 1 if sagitatta, 0 if nemorensis

# 1) Random Forest
rf = RandomForestClassifier(n_estimators=500, random_state=42)
rf.fit(expr.values, labels_binary)
# Plot variable importance: use feature_importances_
importances = rf.feature_importances_
indices = np.argsort(importances)[::-1][:20]  # top 20 genes
top_genes = expr.columns[indices]
plt.figure(figsize=(8, 6))
sns.barplot(x=importances[indices], y=top_genes, palette='viridis')
plt.xlabel("Feature Importance (Mean Decrease in Gini)")
plt.ylabel("Gene")
plt.title("Top 20 Important Genes (Random Forest)")
plt.tight_layout()
plt.show()

# 2) XGBoost (binary classification)
xgb = XGBClassifier(n_estimators=50, use_label_encoder=False, eval_metric='logloss', random_state=42)
xgb.fit(expr.values, labels_binary)
plt.figure(figsize=(8, 6))
plot_importance(xgb, max_num_features=20)
plt.title("Top 20 Important Genes (XGBoost)")
plt.tight_layout()
plt.show()

# 3) SVM (linear)
svm_model = SVC(kernel='linear', probability=True, random_state=42)
# Simple 70/30 train-test split
from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(expr.values, labels_binary, test_size=0.3, random_state=123)
svm_model.fit(X_train, y_train)
y_pred = svm_model.predict(X_test)
from sklearn.metrics import confusion_matrix, accuracy_score
print("SVM Accuracy:", accuracy_score(y_test, y_pred))
print("SVM Confusion Matrix:\n", confusion_matrix(y_test, y_pred))

# 4) Multinomial logistic regression (for multi-label if we had multi-condition/timepoint)
#    Here, simulate multi-label: genotype_condition (e.g., "nemorensis_control", etc.)
multi_label = col_data['genotype'] + "_" + col_data.get('timepoint', pd.Series([""]*len(col_data))).astype(str)
# Note: In original R code, colData(dds)$timepoint was used, but not defined. If available, include.
df_ml = pd.DataFrame(expr.values, index=expr.index, columns=expr.columns)
df_ml['label'] = pd.Categorical(multi_label)
# Fit multinomial logistic regression on all genes (may be slow). We'll do it on top 100 most variable genes instead:
# 5) Select top 100 most variable genes
gene_variances = expr.var(axis=0)
top100 = gene_variances.sort_values(ascending=False).index[:100]
df_top100 = df_ml[top100].copy()
df_top100['label'] = df_ml['label']
# Fit model
multinom = LogisticRegression(multi_class='multinomial', solver='saga', max_iter=1000, random_state=42)
multinom.fit(df_top100[top100], df_top100['label'])
print("Multinomial model classes:", multinom.classes_)

# 6) Boruta feature selection (binary genotype)
boruta_selector = BorutaPy(
    estimator=RandomForestClassifier(n_estimators=100, random_state=42),
    n_estimators='auto',
    verbose=2,
    random_state=42
)
boruta_selector.fit(expr.values, labels_binary)
# Selected features (genes)
selected_genes = expr.columns[boruta_selector.support_].tolist()
print("Boruta selected genes:", selected_genes)

# 7) Survival analysis (Cox Proportional Hazards) using lifelines
#    Assume we have 'days_survived' and 'status' vectors for each sample
#    Here we create placeholders; user should replace with actual data.
days_survived = np.random.randint(5, 100, size=expr.shape[0])  # placeholder
status = np.random.binomial(1, 0.5, size=expr.shape[0])        # placeholder

surv_df = pd.DataFrame(expr.values, index=expr.index, columns=expr.columns)
surv_df['days'] = days_survived
surv_df['status'] = status

# Fit Cox model: use only top 20 variable genes for speed
top20 = gene_variances.sort_values(ascending=False).index[:20]
cox_df = surv_df[top20 + ['days', 'status']]
cph = CoxPHFitter()
cph.fit(cox_df, duration_col='days', event_col='status')
cph.print_summary()

# =============================================================================
# PART II: GxE WILTING VS CONTROL (16 samples)
# =============================================================================

# 1. Set working directory to GxE_ctrl_wilt directory
wd_gxe_wilt = wd + "GxE_ctrl_wilt/"
os.chdir(wd_gxe_wilt)

# 2. List HTSeq-count ".txt" files (16 samples)
sample_files_gxe_wilt = sorted([os.path.basename(f) for f in glob.glob(os.path.join(wd_gxe_wilt, "*.txt"))])

# 3. Define metadata for GxE wilt vs. control
condition_wilt = [
    "control", "control", "control", "control",
    "wilting", "wilting", "wilting", "wilting",
    "control", "control", "control", "control",
    "wilting", "wilting", "wilting", "wilting"
]
genotype_wilt = [
    "nemorensis", "nemorensis", "nemorensis", "nemorensis",
    "nemorensis", "nemorensis", "nemorensis", "nemorensis",
    "sagitatta", "sagitatta", "sagitatta", "sagitatta",
    "sagitatta", "sagitatta", "sagitatta", "sagitatta"
]

if len(sample_files_gxe_wilt) != len(condition_wilt):
    raise ValueError("Length mismatch for GxE wilt vs control metadata.")

sample_table_wilt = pd.DataFrame({
    'sampleName': sample_files_gxe_wilt,
    'fileName': sample_files_gxe_wilt,
    'condition': pd.Categorical(condition_wilt, categories=["control", "wilting"]),
    'genotype': pd.Categorical(genotype_wilt, categories=["nemorensis", "sagitatta"])
})
r_sample_table_wilt = pandas2ri.py2rpy(sample_table_wilt)

ddsHTSeq_gxe_wilt = DESeq2.DESeqDataSetFromHTSeqCount(
    sampleTable=r_sample_table_wilt,
    directory=wd_gxe_wilt,
    design=Formula('~ genotype + condition + genotype:condition')
)

# Filter low-count genes: average count > 100 over 16 samples
num_samples_wilt = 16
count_matrix_wilt = np.array(DESeq2.counts(ddsHTSeq_gxe_wilt))
row_means_wilt = count_matrix_wilt.sum(axis=1) / float(num_samples_wilt)
keep_idx_wilt = np.where(row_means_wilt > count_threshold)[0]
r_keep_wilt = ro.BoolVector([i in keep_idx_wilt for i in range(count_matrix_wilt.shape[0])])
ddsHTSeq_gxe_wilt_filtered = ddsHTSeq_gxe_wilt.rx[r_keep_wilt, True]

# DESeq analysis
dds_gxe_wilt = DESeq2.DESeq(ddsHTSeq_gxe_wilt_filtered, fitType="mean")
res_dds_gxe_wilt = DESeq2.results(dds_gxe_wilt)
res_dds_gxe_wilt_df = pandas2ri.rpy2py(res_dds_gxe_wilt).set_index('rownames')

# Normalization
dds_lowcount_gxe_wilt = DESeq2.estimateSizeFactors(dds_gxe_wilt)
norm_counts_gxe_wilt = pandas2ri.rpy2py(DESeq2.counts(dds_lowcount_gxe_wilt, normalized=True))
norm_counts_gxe_wilt = pd.DataFrame(
    norm_counts_gxe_wilt,
    index=np.array(base.rownames(DESeq2.counts(dds_lowcount_gxe_wilt, normalized=True))),
    columns=np.array(base.colnames(DESeq2.counts(dds_lowcount_gxe_wilt, normalized=True)))
)

# Filter significant DEGs (padj < 0.05)
res_dds_gxe_wilt_signif = res_dds_gxe_wilt_df[res_dds_gxe_wilt_df['padj'] < 0.05]
# res_dds_gxe_wilt_signif.to_csv("res_dds_gxe_wilt_ctrl_significant_3980.csv")

# =============================================================================
# PART III: GxE SURVIVAL VS CONTROL (14 samples)
# =============================================================================

wd_gxe_surv = wd + "GxE_ctrl_surv/"
os.chdir(wd_gxe_surv)

sample_files_gxe_surv = sorted([os.path.basename(f) for f in glob.glob(os.path.join(wd_gxe_surv, "*.txt"))])
condition_surv = [
    "control", "control", "control", "control",
    "survival", "survival", "survival",
    "control", "control", "control", "control",
    "survival", "survival", "survival"
]
genotype_surv = [
    "nemorensis", "nemorensis", "nemorensis", "nemorensis",
    "nemorensis", "nemorensis", "nemorensis",
    "sagitatta", "sagitatta", "sagitatta", "sagitatta",
    "sagitatta", "sagitatta", "sagitatta"
]

if len(sample_files_gxe_surv) != len(condition_surv):
    raise ValueError("Length mismatch for GxE survival vs control metadata.")

sample_table_surv = pd.DataFrame({
    'sampleName': sample_files_gxe_surv,
    'fileName': sample_files_gxe_surv,
    'condition': pd.Categorical(condition_surv, categories=["control", "survival"]),
    'genotype': pd.Categorical(genotype_surv, categories=["nemorensis", "sagitatta"])
})
r_sample_table_surv = pandas2ri.py2rpy(sample_table_surv)

ddsHTSeq_gxe_surv = DESeq2.DESeqDataSetFromHTSeqCount(
    sampleTable=r_sample_table_surv,
    directory=wd_gxe_surv,
    design=Formula('~ genotype + condition + genotype:condition')
)

# Filter low-count genes: average count > 100 over 14 samples
num_samples_surv = 14
count_matrix_surv = np.array(DESeq2.counts(ddsHTSeq_gxe_surv))
row_means_surv = count_matrix_surv.sum(axis=1) / float(num_samples_surv)
keep_idx_surv = np.where(row_means_surv > count_threshold)[0]
r_keep_surv = ro.BoolVector([i in keep_idx_surv for i in range(count_matrix_surv.shape[0])])
ddsHTSeq_gxe_surv_filtered = ddsHTSeq_gxe_surv.rx[r_keep_surv, True]

# DESeq analysis
dds_gxe_surv = DESeq2.DESeq(ddsHTSeq_gxe_surv_filtered, fitType="mean")
res_dds_gxe_surv = DESeq2.results(dds_gxe_surv)
res_dds_gxe_surv_df = pandas2ri.rpy2py(res_dds_gxe_surv).set_index('rownames')

# Normalization
dds_lowcount_gxe_surv = DESeq2.estimateSizeFactors(dds_gxe_surv)
norm_counts_gxe_surv = pandas2ri.rpy2py(DESeq2.counts(dds_lowcount_gxe_surv, normalized=True))
norm_counts_gxe_surv = pd.DataFrame(
    norm_counts_gxe_surv,
    index=np.array(base.rownames(DESeq2.counts(dds_lowcount_gxe_surv, normalized=True))),
    columns=np.array(base.colnames(DESeq2.counts(dds_lowcount_gxe_surv, normalized=True)))
)

# Filter significant DEGs (padj < 0.05)
res_dds_gxe_surv_signif = res_dds_gxe_surv_df[res_dds_gxe_surv_df['padj'] < 0.05]
# res_dds_gxe_surv_signif.to_csv("res_dds_gxe_surv_ctrl_significant_1973.csv")

# =============================================================================
# PART IV: Plot wilt vs. control across species
# =============================================================================

# Extract common genes between nemwilt_sagwilt and nemctrl_sagctrl
common_genes_wilt_ctrl = sig_gene_ids = set(nemwilt_sagwilt_df.index).intersection(set(nemctrl_sagctrl_df.index))

# Subset DataFrames
wilt_common = nemwilt_sagwilt_df.loc[common_genes_wilt_ctrl]
ctrl_common = nemctrl_sagctrl_df.loc[common_genes_wilt_ctrl]

# Build combined DataFrame
combined_df = pd.DataFrame({
    'gene_id': list(common_genes_wilt_ctrl),
    'log2FC_wilt': wilt_common['log2FoldChange'],
    'padj_wilt': wilt_common['padj'],
    'log2FC_ctrl': ctrl_common['log2FoldChange'],
    'padj_ctrl': ctrl_common['padj']
}).set_index('gene_id')

# Define color: red if both padj < 0.05, green if one < 0.05, gray otherwise
def color_mapper(row):
    if (row['padj_wilt'] < 0.05) and (row['padj_ctrl'] < 0.05):
        return "red"
    elif (row['padj_wilt'] < 0.05) or (row['padj_ctrl'] < 0.05):
        return "green"
    else:
        return "gray"

combined_df['color'] = combined_df.apply(color_mapper, axis=1)

# Scatter plot
plt.figure(figsize=(8, 6))
plt.scatter(
    combined_df['log2FC_wilt'],
    combined_df['log2FC_ctrl'],
    c=combined_df['color'],
    s=15,
    alpha=0.7,
    edgecolors='none'
)
plt.xlabel("Log2FC: A. nemorensis wilting vs control", fontsize=12)
plt.ylabel("Log2FC: A. sagitatta wilting vs control", fontsize=12)
plt.title("Wilt vs Control: A. nemorensis vs A. sagitatta", fontsize=14)
plt.axhline(0, color='black', linestyle='--')
plt.axvline(0, color='black', linestyle='--')
plt.tight_layout()
# plt.savefig("wilt_ctrl_scatter.png", dpi=300)
plt.show()


# -------------------------------------------------------------------------
# PART V: Species comparison plots (wilting vs. control and survival vs. control)
# -------------------------------------------------------------------------

import matplotlib.pyplot as plt
import seaborn as sns

# Assumes the following pandas DataFrames already exist from previous steps:
#   sagwilt_sagctrl_df:   results for A. sagittata wilting vs. control (indexed by gene)
#   nemwilt_nemctrl_df:   results for A. nemorensis wilting vs. control (indexed by gene)
#   res_dds_gxe:          GxE results (wilting) for all genes (indexed by gene)
#   sagsurv_sagctrl_df:   results for A. sagittata survival vs. control (indexed by gene)
#   nemsurv_nemctrl_df:   results for A. nemorensis survival vs. control (indexed by gene)
#   res_dds_gxe_surv:     GxE results (survival) for all genes (indexed by gene)

# --- 1) Plot wilting vs. control species comparison with GxE overlay ---

# 1.1 Find common genes between sagittata-wilt-ctrl and nemorensis-wilt-ctrl
common_genes_sag_nem = sagwilt_sagctrl_df.index.intersection(nemwilt_nemctrl_df.index)

# 1.2 Subset the two result DataFrames to those common genes
sag_common = sagwilt_sagctrl_df.loc[common_genes_sag_nem]
nem_common = nemwilt_nemctrl_df.loc[common_genes_sag_nem]

# 1.3 Build combined DataFrame
combined_species_df_sag_nem = pd.DataFrame({
    'gene_id':       common_genes_sag_nem,
    'log2FC_sag':    sag_common['log2FoldChange'],
    'padj_sag':      sag_common['padj'],
    'log2FC_nem':    nem_common['log2FoldChange'],
    'padj_nem':      nem_common['padj']
}).set_index('gene_id')

# 1.4 Add GxE padj values (wilting) for those common genes
#     (res_dds_gxe has GxE results from Part II, indexed by gene)
gxe_common = res_dds_gxe.loc[common_genes_sag_nem]
combined_species_df_sag_nem['padj_gxe'] = gxe_common['padj']

# 1.5 Drop any rows with NA (if any)
combined_species_df_sag_nem = combined_species_df_sag_nem.dropna(subset=['padj_gxe', 'padj_sag', 'padj_nem'])

# 1.6 Assign color categories:
#       - "red"   if padj_gxe < 0.001
#       - "green" if padj_sag < 0.05 AND padj_nem < 0.05
#       - "gray"  otherwise
def assign_color_wilt(row):
    if row['padj_gxe'] < 0.001:
        return "red"
    elif (row['padj_sag'] < 0.05) and (row['padj_nem'] < 0.05):
        return "green"
    else:
        return "gray"

combined_species_df_sag_nem['color'] = combined_species_df_sag_nem.apply(assign_color_wilt, axis=1)

# 1.7 Create scatter plot with layered points (gray first, then green, then red)
plt.figure(figsize=(11, 7))
# Plot gray points
gray_df = combined_species_df_sag_nem[combined_species_df_sag_nem['color'] == 'gray']
plt.scatter(
    gray_df['log2FC_nem'], gray_df['log2FC_sag'],
    c='gray', s=0.7, alpha=0.7, label='NS'
)
# Plot green points
green_df = combined_species_df_sag_nem[combined_species_df_sag_nem['color'] == 'green']
plt.scatter(
    green_df['log2FC_nem'], green_df['log2FC_sag'],
    c='darkgreen', s=0.7, alpha=0.7, label='Significant E'
)
# Plot red points
red_df = combined_species_df_sag_nem[combined_species_df_sag_nem['color'] == 'red']
plt.scatter(
    red_df['log2FC_nem'], red_df['log2FC_sag'],
    c='darkred', s=0.7, alpha=0.7, label='GxE (padj < 0.001)'
)

# Add dashed center lines
plt.axvline(0, linestyle='--', color='black')
plt.axhline(0, linestyle='--', color='black')

# Labels, title, limits, and legend
plt.xlabel("Log₂FC A. nemorensis (stress vs control)", fontsize=24, fontweight='bold')
plt.ylabel("Log₂FC A. sagittata (stress vs control)", fontsize=24, fontweight='bold')
plt.title("Exp. diff. b/w A. nemorensis and A. sagittata at stress", fontsize=26)
plt.xlim(-10, 10)
plt.ylim(-10, 10)

# Custom legend
legend_handles = [
    plt.Line2D([0], [0], marker='o', color='w', label='NS',
               markerfacecolor='gray', markersize=8),
    plt.Line2D([0], [0], marker='o', color='w', label='Significant E',
               markerfacecolor='darkgreen', markersize=8),
    plt.Line2D([0], [0], marker='o', color='w', label='GxE (padj < 0.001)',
               markerfacecolor='darkred', markersize=8)
]
plt.legend(handles=legend_handles, fontsize=20, frameon=False, title="Significance")

plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.tight_layout()

plt.savefig("wilt_vs_ctrl_species_comparison.png", dpi=300)
plt.savefig("wilt_vs_ctrl_species_comparison.pdf", dpi=300)
plt.savefig("wilt_vs_ctrl_species_comparison.svg", dpi=300)
plt.close()


# --- 2) Plot survival vs. control species comparison with GxE overlay ---

# 2.1 Find common genes between sagittata-surv-ctrl and nemorensis-surv-ctrl
common_genes_sag_nem_surv = sagsurv_sagctrl_df.index.intersection(nemsurv_nemctrl_df.index)

# 2.2 Subset the two DataFrames
sag_common_surv = sagsurv_sagctrl_df.loc[common_genes_sag_nem_surv]
nem_common_surv = nemsurv_nemctrl_df.loc[common_genes_sag_nem_surv]

# 2.3 Build combined DataFrame for survival
combined_species_df_sag_nem_surv = pd.DataFrame({
    'gene_id':       common_genes_sag_nem_surv,
    'log2FC_sag':    sag_common_surv['log2FoldChange'],
    'padj_sag':      sag_common_surv['padj'],
    'log2FC_nem':    nem_common_surv['log2FoldChange'],
    'padj_nem':      nem_common_surv['padj']
}).set_index('gene_id')

# 2.4 Add GxE padj values (survival) for those common genes
gxe_common_surv = res_dds_gxe_surv.loc[common_genes_sag_nem_surv]
combined_species_df_sag_nem_surv['padj_gxe'] = gxe_common_surv['padj']

# 2.5 Drop any rows with NA
combined_species_df_sag_nem_surv = combined_species_df_sag_nem_surv.dropna(subset=['padj_gxe', 'padj_sag', 'padj_nem'])

# 2.6 Assign color categories:
#       - "red"   if padj_gxe < 0.05
#       - "green" if padj_sag < 0.05 AND padj_nem < 0.05
#       - "gray"  otherwise
def assign_color_surv(row):
    if row['padj_gxe'] < 0.05:
        return "red"
    elif (row['padj_sag'] < 0.05) and (row['padj_nem'] < 0.05):
        return "green"
    else:
        return "gray"

combined_species_df_sag_nem_surv['color'] = combined_species_df_sag_nem_surv.apply(assign_color_surv, axis=1)

# 2.7 Create scatter plot with layered points (gray, then green, then red)
plt.figure(figsize=(11, 7))
# Plot gray points
gray_df_surv = combined_species_df_sag_nem_surv[combined_species_df_sag_nem_surv['color'] == 'gray']
plt.scatter(
    gray_df_surv['log2FC_nem'], gray_df_surv['log2FC_sag'],
    c='gray', s=0.7, alpha=0.7, label='NS'
)
# Plot green points
green_df_surv = combined_species_df_sag_nem_surv[combined_species_df_sag_nem_surv['color'] == 'green']
plt.scatter(
    green_df_surv['log2FC_nem'], green_df_surv['log2FC_sag'],
    c='darkgreen', s=0.7, alpha=0.7, label='Significant E'
)
# Plot red points
red_df_surv = combined_species_df_sag_nem_surv[combined_species_df_sag_nem_surv['color'] == 'red']
plt.scatter(
    red_df_surv['log2FC_nem'], red_df_surv['log2FC_sag'],
    c='darkred', s=0.7, alpha=0.7, label='GxE (padj < 0.05)'
)

# Add dashed center lines
plt.axvline(0, linestyle='--', color='black')
plt.axhline(0, linestyle='--', color='black')

# Labels, title, limits, and legend
plt.xlabel("Log₂FC A. nemorensis (recovery vs control)", fontsize=24, fontweight='bold')
plt.ylabel("Log₂FC A. sagittata (recovery vs control)", fontsize=24, fontweight='bold')
plt.title("Exp. diff. b/w A. nemorensis and A. sagittata at recovery", fontsize=26)
plt.xlim(-10, 10)
plt.ylim(-10, 10)

# Custom legend
legend_handles_surv = [
    plt.Line2D([0], [0], marker='o', color='w', label='NS',
               markerfacecolor='gray', markersize=8),
    plt.Line2D([0], [0], marker='o', color='w', label='Significant E',
               markerfacecolor='darkgreen', markersize=8),
    plt.Line2D([0], [0], marker='o', color='w', label='GxE (padj < 0.05)',
               markerfacecolor='darkred', markersize=8)
]
plt.legend(handles=legend_handles_surv, fontsize=20, frameon=False, title="Significance")

plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.tight_layout()

plt.savefig("surv_vs_ctrl_species_comparison.png", dpi=300)
plt.savefig("surv_vs_ctrl_species_comparison.pdf", dpi=300)
plt.savefig("surv_vs_ctrl_species_comparison.svg", dpi=300)
plt.close()


# -------------------------------------------------------------------------
# PART V (continued): Quadrant classification, orthologue merging, universe splitting,
# and GO enrichment analyses (wilting vs. control and survival vs. control)
# -------------------------------------------------------------------------

import os
import pandas as pd
import numpy as np

# For GO enrichment, we use goatools
from goatools.obo_parser import GODag
from goatools.go_enrichment import GOEnrichmentStudy

# -----------------------------------------------------------------------------
# 1) CLASSIFY POINTS INTO QUADRANTS
# -----------------------------------------------------------------------------

# Assume the following pandas DataFrames already exist from earlier steps:
#   combined_species_df_sag_nem       (stress vs. control comparison)
#   combined_species_df_sag_nem_surv  (recovery vs. control comparison)
#
# Each has columns: 'log2FC_sag', 'log2FC_nem', 'padj_sag', 'padj_nem', 'padj_gxe', 'color'

# Create a helper function to assign quadrants
def assign_quadrant(df):
    df = df.copy()
    df['above_diag'] = df['log2FC_sag'] > df['log2FC_nem']
    df['above_hline'] = df['log2FC_sag'] > 0
    df['right_vline'] = df['log2FC_nem'] > 0

    def quadrant_label(row):
        if     row['above_diag']  and row['above_hline']  and row['right_vline']:
            return "Q1"
        elif   row['above_diag']  and row['above_hline']  and not row['right_vline']:
            return "Q2"
        elif  not row['above_diag'] and row['above_hline']  and not row['right_vline']:
            return "Q3"
        elif  not row['above_diag'] and not row['above_hline'] and not row['right_vline']:
            return "Q4"
        elif  not row['above_diag'] and not row['above_hline'] and row['right_vline']:
            return "Q5"
        elif   row['above_diag']  and not row['above_hline'] and row['right_vline']:
            return "Q6"
        elif   row['above_diag']  and not row['above_hline'] and not row['right_vline']:
            return "Q7"
        elif  not row['above_diag'] and row['above_hline']  and row['right_vline']:
            return "Q8"
        else:
            return np.nan

    df['quadrant'] = df.apply(quadrant_label, axis=1)
    return df

# Apply quadrant classification to both stress and recovery DataFrames
combined_species_df_sag_nem = assign_quadrant(combined_species_df_sag_nem)
combined_species_df_sag_nem_surv = assign_quadrant(combined_species_df_sag_nem_surv)

# -----------------------------------------------------------------------------
# 2) MERGE WITH ORTHOLOGUES AND FILTER
# -----------------------------------------------------------------------------

# Path to the orthologues CSV
orthologues_path = (
    "/Users/Shared/Files From d.localized/PhD-Uni-Koeln_2021-2024/PhD_work/"
    "seeds data/2nd Experiment 2022/transcriptome/Arabis_RNA_raw_data/"
    "90-774066047/anaylsis_drought_mRNA_both_species_with_new_two_genomes/"
    "analysis_both_species_with_nem_genomes/orthologues_cleaned.csv"
)

orthologues = pd.read_csv(orthologues_path, header=0)
# orthologues must contain columns: 'arabis_cleaned' (gene_id) and 'At' (Arabidopsis ID)

# Merge for stress DataFrame
combined_species_df_sag_nem = (
    combined_species_df_sag_nem
    .merge(
        orthologues,
        left_on='gene_id',
        right_on='arabis_cleaned',
        how='left'
    )
    .dropna(subset=['At'])
)

# Merge for recovery DataFrame
combined_species_df_sag_nem_surv = (
    combined_species_df_sag_nem_surv
    .merge(
        orthologues,
        left_on='gene_id',
        right_on='arabis_cleaned',
        how='left'
    )
    .dropna(subset=['At'])
)

# -----------------------------------------------------------------------------
# 3) DEFINE UNIVERSES BASED ON STRESS vs. CONTROL
# -----------------------------------------------------------------------------

# Stress DataFrame: combined_species_df_sag_nem
upper_universe_sag_nem = combined_species_df_sag_nem[
    combined_species_df_sag_nem['log2FC_sag'] > combined_species_df_sag_nem['log2FC_nem']
].copy()

lower_universe_sag_nem = combined_species_df_sag_nem[
    combined_species_df_sag_nem['log2FC_sag'] < combined_species_df_sag_nem['log2FC_nem']
].copy()

left_universe_sag_nem = combined_species_df_sag_nem[
    combined_species_df_sag_nem['log2FC_nem'] < 0
].copy()

right_universe_sag_nem = combined_species_df_sag_nem[
    combined_species_df_sag_nem['log2FC_nem'] > 0
].copy()

# -----------------------------------------------------------------------------
# 4) DEFINE UNIVERSES BASED ON RECOVERY vs. CONTROL
# -----------------------------------------------------------------------------

# Recovery DataFrame: combined_species_df_sag_nem_surv
upper_universe_sag_nem_surv = combined_species_df_sag_nem_surv[
    combined_species_df_sag_nem_surv['log2FC_sag'] > combined_species_df_sag_nem_surv['log2FC_nem']
].copy()

lower_universe_sag_nem_surv = combined_species_df_sag_nem_surv[
    combined_species_df_sag_nem_surv['log2FC_sag'] < combined_species_df_sag_nem_surv['log2FC_nem']
].copy()

left_universe_sag_nem_surv = combined_species_df_sag_nem_surv[
    combined_species_df_sag_nem_surv['log2FC_nem'] < 0
].copy()

right_universe_sag_nem_surv = combined_species_df_sag_nem_surv[
    combined_species_df_sag_nem_surv['log2FC_nem'] > 0
].copy()

# -----------------------------------------------------------------------------
# 5) GO ENRICHMENT SETUP
# -----------------------------------------------------------------------------

# Path to GO OBO file and gene2go annotation file (user must supply correct paths)
go_obo_path = "/path/to/go-basic.obo"
gene2go_path = "/path/to/gene2go"  # can be a pickled dict or a tab-delimited file

# Load GO DAG
obodag = GODag(go_obo_path)

# Load gene-to-GO associations.
# For example, if 'gene2go_path' is a tab-delimited file with "gene_id\tGO:xxxxx" lines:
def load_gene2go(filepath):
    gene2go = {}
    with open(filepath) as f:
        for line in f:
            gene, go_term = line.strip().split('\t')
            gene2go.setdefault(gene, set()).add(go_term)
    return gene2go

# If gene2go_path is already a pickled Python dict, load with:
# import pickle
# geneid2gos = pickle.load(open(gene2go_path, "rb"))

# Otherwise, read from a file:
geneid2gos = load_gene2go(gene2go_path)

# -----------------------------------------------------------------------------
# 6) FUNCTION TO RUN GO ENRICHMENT FOR A GIVEN UNIVERSE AND QUADRANT
# -----------------------------------------------------------------------------

def run_go_enrichment(universe_df, quadrant_label, out_fname):
    """
    universe_df: pandas DataFrame for a particular universe (left or right).
                 Must have columns 'quadrant', 'color', and 'At'.
    quadrant_label: e.g. "Q1", "Q2", etc.
    out_fname: string path to save CSV results (top 50 GO terms).
    """
    # Determine study genes (At IDs where quadrant == quadrant_label and color == 'red')
    study_genes = universe_df[
        (universe_df['quadrant'] == quadrant_label) & 
        (universe_df['color'] == 'red')
    ]['At'].astype(str).tolist()

    # Population (background) genes: all At IDs in this universe
    population_genes = universe_df['At'].astype(str).tolist()

    # If no study genes, skip
    if len(study_genes) == 0:
        print(f"No red genes in {quadrant_label} for this universe. Skipping GO enrichment.")
        return None

    # Initialize GOEnrichmentStudy
    goe = GOEnrichmentStudy(
        population_genes,
        geneid2gos,
        obodag,
        propagate_counts=False,  # topGO's "elim" is roughly analogous to not propagating counts
        alpha=0.05,
        methods=['fisher']
    )

    # Run the enrichment
    results = goe.run_study(study_genes)

    # Convert results to pandas DataFrame and select top 50 by p_fdr_bh
    records = []
    for r in results:
        if r.p_fdr_bh is not None:
            records.append({
                'GO': r.GO,
                'name': r.name,
                'NS': r.NS,
                'n_Pop': r.n_pop,
                'n_S': r.n_s,
                'p_uncorrected': r.p_uncorrected,
                'p_fdr_bh': r.p_fdr_bh
            })
    if not records:
        print(f"No enriched GO terms found for {quadrant_label}.")
        return None

    res_df = pd.DataFrame.from_records(records)
    res_df = res_df.sort_values('p_fdr_bh').head(50)
    res_df.to_csv(out_fname, index=False)
    print(f"GO enrichment for {quadrant_label} saved to {out_fname}")
    return res_df

# -----------------------------------------------------------------------------
# 7) EXECUTE GO ENRICHMENTS FOR WILTING vs. CONTROL (STRESS UNIVERSES)
# -----------------------------------------------------------------------------

# Stress left universe: Q2, Q4, Q7
stress_left_tasks = [
    ("Q2", "goEnrichmentQ_sag_up_and_nem_down_wiltQ2.csv"),
    ("Q4", "goEnrichmentQ_sag_down_more_against_nem_wiltQ4.csv"),
    ("Q7", "goEnrichmentQ_nem_down_more_against_sag_wiltQ7.csv"),
]
for quad, fname in stress_left_tasks:
    run_go_enrichment(left_universe_sag_nem, quad, fname)

# Stress right universe: Q1, Q5, Q8
stress_right_tasks = [
    ("Q1", "goEnrichmentQ_sag_up_more_against_nem_wiltQ1.csv"),
    ("Q5", "goEnrichmentQ_nem_up_and_sag_down_wiltQ5.csv"),
    ("Q8", "goEnrichmentQ_nem_up_more_against_sag_wiltQ8.csv"),
]
for quad, fname in stress_right_tasks:
    run_go_enrichment(right_universe_sag_nem, quad, fname)

# -----------------------------------------------------------------------------
# 8) EXECUTE GO ENRICHMENTS FOR SURVIVAL vs. CONTROL (RECOVERY UNIVERSES)
# -----------------------------------------------------------------------------

# Recovery left universe: Q7, Q4, Q2
rec_left_tasks = [
    ("Q7", "goEnrichmentQ_sag_down_more_against_nem_recQ7.csv"),
    ("Q4", "goEnrichmentQ_sag_down_more_against_nem_recQ4.csv"),
    ("Q2", "goEnrichmentQ_sag_up_and_nem_down_recQ2.csv"),
]
for quad, fname in rec_left_tasks:
    run_go_enrichment(left_universe_sag_nem_surv, quad, fname)

# Recovery right universe: Q1, Q5, Q8
rec_right_tasks = [
    ("Q1", "goEnrichmentQ_sag_up_more_against_nem_recQ1.csv"),
    ("Q5", "goEnrichmentQ_nem_up_and_sag_down_recQ5.csv"),
    ("Q8", "goEnrichmentQ_nem_up_more_against_sag_recQ8.csv"),
]
for quad, fname in rec_right_tasks:
    run_go_enrichment(right_universe_sag_nem_surv, quad, fname)

# -------------------------------------------------------------------------
# PART VI: Extract significant genes per GO term from GO enrichment results
# -------------------------------------------------------------------------

import pandas as pd

# We assume that, for each quadrant (e.g., Q1, Q2, …), you have already run GOEnrichmentStudy
# and have:
#   1) A pandas DataFrame `goEnrichDF_Qx` containing the top-enriched GO terms, with a column 'GO' (GO ID).
#   2) A list `go_results_Qx` of GOEnrichmentRecord objects (the full result of GOEnrichmentStudy.run_study).
#
# For example:
#   goEnrichDF_Q1 = pd.read_csv("goEnrichmentQ_sag_up_more_against_nem_wiltQ1.csv")
#   go_results_Q1 = goe_Q1_results  # the list returned by GOEnrichmentStudy.run_study for Q1
#
# If you followed the earlier `run_go_enrichment` function, you can modify it to return both the DataFrame
# and the full results list. Here, we illustrate how to extract the “genes in each significant GO” for each Q.

# 1) Helper function to flatten GO → study_items into a DataFrame
def extract_sig_genes_per_go(go_enrich_df: pd.DataFrame, go_results: list, out_csv: str):
    """
    Given:
      - go_enrich_df: DataFrame with at least a column 'GO' listing significant GO IDs.
      - go_results:   list of GOEnrichmentRecord objects (from goatools) for that quadrant.
      - out_csv:      path to save the flattened gene‐GO mapping CSV.
    This function finds, for each GO ID in go_enrich_df['GO'], the set of study genes
    annotated to that GO (record.study_items), then writes a CSV with columns [GeneID, GO].
    """
    # 1.1 Build a mapping GO_ID → set_of_genes (study_items)
    go2genes = {}
    sig_go_list = go_enrich_df['GO'].astype(str).tolist()

    # Create a dict for quick lookup: GO_ID → GOEnrichmentRecord
    results_dict = {r.GO: r for r in go_results}

    for go_id in sig_go_list:
        if go_id in results_dict:
            record = results_dict[go_id]
            # .study_items is a set of gene IDs in the study that map to this GO
            go2genes[go_id] = record.study_items.copy()
        else:
            # If the GO ID is not in the results list, skip
            go2genes[go_id] = set()

    # 1.2 Flatten into a list of (GeneID, GO_ID) rows
    rows = []
    for go_id, geneset in go2genes.items():
        for gene in geneset:
            rows.append({'GeneID': gene, 'GO.ID': go_id})

    # 1.3 Convert to DataFrame
    go2genes_df = pd.DataFrame(rows)

    # 1.4 Save to CSV
    go2genes_df.to_csv(out_csv, index=False)
    print(f"Saved GO‐gene mappings to {out_csv}")

    return go2genes_df


# -----------------------------------------------------------------------------
# 2) EXAMPLE USAGE for “stress” quadrants (wilting vs. control)
# -----------------------------------------------------------------------------

# 2.1 Load the GO enrichment results DataFrames that you previously saved (top 50 GO terms):
goEnrichDF_Q1 = pd.read_csv("goEnrichmentQ_sag_up_more_against_nem_wiltQ1.csv")
goEnrichDF_Q2 = pd.read_csv("goEnrichmentQ_sag_up_and_nem_down_wiltQ2.csv")
goEnrichDF_Q4 = pd.read_csv("goEnrichmentQ_sag_down_more_against_nem_wiltQ4.csv")
goEnrichDF_Q5 = pd.read_csv("goEnrichmentQ_nem_up_and_sag_down_wiltQ5.csv")
goEnrichDF_Q7 = pd.read_csv("goEnrichmentQ_nem_down_more_against_sag_wiltQ7.csv")
goEnrichDF_Q8 = pd.read_csv("goEnrichmentQ_nem_up_more_against_sag_wiltQ8.csv")

# 2.2 Assume that when you ran GOEnrichmentStudy, you captured the full results list:
# For example:
#   results_Q1 = goe_Q1.run_study(study_genes_Q1)
#   results_Q2 = goe_Q2.run_study(study_genes_Q2)
#   ...
# (Each is a list of GOEnrichmentRecord objects.)

# Here we illustrate with placeholders; replace these with your actual variables:
results_Q1 = go_results_Q1  # list of GOEnrichmentRecord for quadrant Q1
results_Q2 = go_results_Q2  # list of GOEnrichmentRecord for quadrant Q2
results_Q4 = go_results_Q4  # list of GOEnrichmentRecord for quadrant Q4
results_Q5 = go_results_Q5  # list of GOEnrichmentRecord for quadrant Q5
results_Q7 = go_results_Q7  # list of GOEnrichmentRecord for quadrant Q7
results_Q8 = go_results_Q8  # list of GOEnrichmentRecord for quadrant Q8

# 2.3 Extract and save the gene lists for each quadrant
go2genes_Q1_df = extract_sig_genes_per_go(
    go_enrich_df=goEnrichDF_Q1,
    go_results=results_Q1,
    out_csv="goEnrichmentQ_sag_up_more_against_nem_wiltQ1_genes.csv"
)

go2genes_Q2_df = extract_sig_genes_per_go(
    go_enrich_df=goEnrichDF_Q2,
    go_results=results_Q2,
    out_csv="goEnrichmentQ_sag_up_and_nem_down_wiltQ2_genes.csv"
)

go2genes_Q4_df = extract_sig_genes_per_go(
    go_enrich_df=goEnrichDF_Q4,
    go_results=results_Q4,
    out_csv="goEnrichmentQ_sag_down_more_against_nem_wiltQ4_genes.csv"
)

go2genes_Q5_df = extract_sig_genes_per_go(
    go_enrich_df=goEnrichDF_Q5,
    go_results=results_Q5,
    out_csv="goEnrichmentQ_nem_up_and_sag_down_wiltQ5_genes.csv"
)

go2genes_Q7_df = extract_sig_genes_per_go(
    go_enrich_df=goEnrichDF_Q7,
    go_results=results_Q7,
    out_csv="goEnrichmentQ_nem_down_more_against_sag_wiltQ7_genes.csv"
)

go2genes_Q8_df = extract_sig_genes_per_go(
    go_enrich_df=goEnrichDF_Q8,
    go_results=results_Q8,
    out_csv="goEnrichmentQ_nem_up_more_against_sag_wiltQ8_genes.csv"
)

# -----------------------------------------------------------------------------
# 3) EXAMPLE USAGE for “recovery” quadrants (survival vs. control)
# -----------------------------------------------------------------------------

goEnrichDF_recQ1 = pd.read_csv("goEnrichmentQ_sag_up_more_against_nem_recQ1.csv")
goEnrichDF_recQ2 = pd.read_csv("goEnrichmentQ_sag_up_and_nem_down_recQ2.csv")
goEnrichDF_recQ4 = pd.read_csv("goEnrichmentQ_sag_down_more_against_nem_recQ4.csv")
goEnrichDF_recQ5 = pd.read_csv("goEnrichmentQ_nem_up_and_sag_down_recQ5.csv")
goEnrichDF_recQ7 = pd.read_csv("goEnrichmentQ_nem_down_more_against_sag_recQ7.csv")
goEnrichDF_recQ8 = pd.read_csv("goEnrichmentQ_nem_up_more_against_sag_recQ8.csv")

# Similarly, replace these with your actual GOEnrichmentRecord lists:
results_recQ1 = go_results_recQ1
results_recQ2 = go_results_recQ2
results_recQ4 = go_results_recQ4
results_recQ5 = go_results_recQ5
results_recQ7 = go_results_recQ7
results_recQ8 = go_results_recQ8

# Extract and save gene lists:
go2genes_recQ1_df = extract_sig_genes_per_go(
    go_enrich_df=goEnrichDF_recQ1,
    go_results=results_recQ1,
    out_csv="goEnrichmentQ_sag_up_more_against_nem_recQ1_genes.csv"
)

go2genes_recQ2_df = extract_sig_genes_per_go(
    go_enrich_df=goEnrichDF_recQ2,
    go_results=results_recQ2,
    out_csv="goEnrichmentQ_sag_up_and_nem_down_recQ2_genes.csv"
)

go2genes_recQ4_df = extract_sig_genes_per_go(
    go_enrich_df=goEnrichDF_recQ4,
    go_results=results_recQ4,
    out_csv="goEnrichmentQ_sag_down_more_against_nem_recQ4_genes.csv"
)

go2genes_recQ5_df = extract_sig_genes_per_go(
    go_enrich_df=goEnrichDF_recQ5,
    go_results=results_recQ5,
    out_csv="goEnrichmentQ_nem_up_and_sag_down_recQ5_genes.csv"
)

go2genes_recQ7_df = extract_sig_genes_per_go(
    go_enrich_df=goEnrichDF_recQ7,
    go_results=results_recQ7,
    out_csv="goEnrichmentQ_nem_down_more_against_sag_recQ7_genes.csv"
)

go2genes_recQ8_df = extract_sig_genes_per_go(
    go_enrich_df=goEnrichDF_recQ8,
    go_results=results_recQ8,
    out_csv="goEnrichmentQ_nem_up_more_against_sag_recQ8_genes.csv"
)

# =============================================================================
# END OF PIPELINE
# =============================================================================
