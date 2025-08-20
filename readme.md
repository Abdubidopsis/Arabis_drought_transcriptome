# Arabis Drought Transcriptome

This repository contains supplemental data and scripts associated with the analysis of transcriptomic responses to drought stress in Arabis species.

## 🧪 Project Description

This study explores the differential gene expression and regulatory patterns in **Arabis sagittata** and **Arabis nemorensis** under drought stress. It includes high-resolution transcriptome analysis, Expression prediction and enrichments, and exploration of miRNA-target interactions.

## 📁 Repository Contents

| File Name                                      | Description |
|-----------------------------------------------|-------------|
| `Arabis_Drought_transcriptome.R`               | Main R script used for transcriptome analysis and visualization. |
| `Table-S1.xlsx`                                | Phenotypes of Arabis nemonresis and Arabis sagittata measured during the dry down experiment. | 
| `Table-S3.xlsx`                                | Analysis of differential gene expression in Arabis nemonresis and Arabis sagittata in stress and recovery. |
| `Table-S7.ods`                                 | Details of extracted EPMs from control, wilting and recovery treatment of DeepCRE-like models including cluster assignment, positional preferences and significant similarities to JASPAR TF database. |
| `Table-S9.ods`                                 | Filtered annotations of EPM occurrence within the genomes of A. nemorensis and A. sagitatta |
| `Table-S10.ods`                                | Analyzes gene expression data linking differentially expressed genes to Gene Ontology terms with EPMs for A. nemorensis and A. sagitatta. |
| `Table-S12.ods`                                | EPM occurrence across Arabis species. |
| `Table-S13.ods`                                | Occurrence of EPMs and characterized TFBS in miRNA408 of A. nemorensis and A. sagitatta. |
| `Table-S14.ods`                                | Occurrence of EPMs and characterized TFBS in potential miRNA408 targets of A. nemorensis and A. sagitatta. |
| `File-S1.fas`                                  | Sequence alignment of miRNA408 with EPMs for A. thaliana, A. nemorensis and A. sagitatta. |
| `File-S2.fas`                                  | Sequence alignment of potential miRNA408 target with EPMs for A. thaliana, A. nemorensis and A. sagitatta. |

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
