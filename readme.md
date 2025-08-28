# 🌱 Arabis Drought Transcriptome

This repository contains supplemental data and scripts associated with the analysis of transcriptomic responses to drought stress in *Arabis* species (*A. sagittata* and *A. nemorensis*).

---

## 🧪 Project Description
This study explores how closely related *Arabis* species respond to drought stress at the molecular level.  
It includes:
- High-resolution transcriptome analysis
- Differential gene expression and functional enrichment
- miRNA–target interactions (with a focus on miR408)
- Expression Predictive Motifs (EPMs) discovered using machine learning
- Comparative analysis of regulatory elements across species

---

## 📁 Repository Contents

### 🔬 Analysis Scripts
R
	•	Arabis_Drought_transcriptome.R — primary workflow for drought transcriptome analysis and visualization
	•	EPM_DGEenrichment_analyses.v3.2.R — functional enrichment of EPM-associated differentially expressed genes
	•	EPM_Quadrant_enrichment_2025F.v2.1.R — quadrant-based enrichment of motif–gene associations
	•	EPM_cluster_venn_analysis.sz.v1.k25.R — cluster overlap and Venn analysis across treatments/species
	•	GO-EPMs_species_comp.sz.v1.2.cluster.R — comparative Gene Ontology enrichment linked to EPMs
	•	add_logMedTPM_no0_balanced_classification.2.R — balanced classification pipeline for expression prediction
	•	FDR_calc_multi.v1.0.R — multiple-testing FDR calculations
	•	corr_prediction_multi.v1.0.R — correlation-based expression prediction across datasets

Python
	•	epm_to_reference_alingment.0.py — maps EPM sequences to the reference (alignment)
	•	extract_fasta_seq.0.py — extracts sequences from FASTA by ID
	•	generate_deg_go_analyses_with_plots.2.py — runs DEG + GO analyses and produces plots
	•	generate_saliency_map_results_plots.up-do.py — visualizes saliency-map outputs
	•	generate_training_results_plots.0.py — generates model training results plots

### 📊 Supplemental Data
- **Table-S1.xlsx** – Phenotypic data from dry-down experiment (*A. nemorensis* & *A. sagittata*)  
- **Table-S3.xlsx** – Differential expression analysis across stress and recovery phases  
- **Table-S7.ods** – Extracted EPMs from treatments with cluster assignment & TF database matches  
- **Table-S9.ods** – Filtered genome-wide annotations of EPM occurrence  
- **Table-S10.ods** – DEGs linked to Gene Ontology terms and associated EPMs
- **Table-S12.ods** – EPM occurrence across Arabis species   
- **Table-S13.ods** – EPM and TFBS occurrence in *miR408* locus  
- **Table-S14.ods** – EPM and TFBS occurrence in predicted *miR408* targets  

### 🧬 Sequence Alignments
- **File-S1.fas** – Alignment of *miR408* with EPMs across *A. thaliana*, *A. nemorensis*, *A. sagittata*  
---

## 📊 Methods
- Transcript quantification with **DESeq2**  
- Functional enrichment using **topGO**  
- Machine learning (deepCRE-like models) for **Expression Predictive Motifs (EPMs)**  
- miRNA prediction & target pairing  
- Orthologue identification using **BLAST** and **OrthoFinder**

---

## 💻 Usage
To reproduce results:
1. Retrieve raw RNA-seq data (see paper).  
2. Process reads → count tables.  
3. Install R and required packages (see `Arabis_Drought_transcriptome.R`).  
4. Run analysis scripts as needed.  
5. Reference supplemental tables for detailed results.  

---

## 📜 License
This repository is distributed under the license specified in the related research paper. Please cite appropriately when reusing data or code.  

---

## 📬 Contact
For questions, please contact the **corresponding author** listed in the associated publication.  
