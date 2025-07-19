# Arabis Drought Transcriptome

This repository contains supplemental data and scripts associated with the analysis of transcriptomic responses to drought stress in Arabis species.

## 🧪 Project Description

This study explores the differential gene expression and regulatory patterns in **Arabis sagittata** and **Arabis nemorensis** under drought stress. It includes high-resolution transcriptome analysis, functional annotation, and exploration of miRNA-target interactions.

## 📁 Repository Contents

| File Name                                      | Description |
|-----------------------------------------------|-------------|
| `Arabis_Drought_transcriptome.R`              | Main R script used for transcriptome analysis and visualization. |
| `Supp-stat-1.xlsx`                             | Summary statistics of RNA-seq reads. |
| `Supp-stat-2.xlsx`                             | DEGs under drought in A. sagittata. |
| `Supp-stat-3_annotation.xlsx`                 | Gene annotations of DEGs. |
| `Supp-stat-4_orthologues.xlsx`                | Orthologue mapping between Arabis and Arabidopsis. |
| `Supp-stat-5_all_stress_targets_Asagittata.xlsx` | miRNA target genes under stress in A. sagittata. |
| `Supp-stat-6_EPM_filtered_occurence.xlsx',`Supp-stat-7.ods`, `Supp-stat-8.ods`, `Supp-stat-9.ods`Supp-stat-10.ods`,`Supp-stat-11.ods`, `Supp-stat-12.xlsx`, `Supp-stat-13.ods` | Expression Predictive Motifs related analysis |

## 📊 Methods

- Transcript quantification using DESeq2
- Functional enrichment using topGO
- Machine Learning based Expression Predictvie Motif clusters
- miRNA prediction and target pairing
- Orthologue identification with BLAST and OrthoFinder

## 💻 How to Use

To reproduce the results:
1. Go through the research paper, get the raw data and process it until you have the count data.
2. Install R and required packages listed in `Arabis_Drought_transcriptome.R`.
3. Open and run the script.
4. Reference supplemental files for interpretation and exploration.

## 📜 License

This project is distributed under (See license in the researcg paper) — feel free to use and adapt with proper citation.

## 📬 Contact

For questions, please contact: corresponding author in the research paper

---
