# Overview
This repository contains the code for TCGA and ICB cohorts' transcriptome data analysis used in the paper titled “Identification of a Tertiary Lymphoid Structure Signature for Predicting Tumor Outcomes through Transcriptomics Analysis”.

# Software Dependencies

Core Languages:  
R (≥ 3.6.0) with Bioconductor packages  

Key R Packages:  
library(tidyverse)      # Data wrangling and visualization  
library(data.table)     # Efficient data handling  
library(GSVA)           # Gene set variation analysis  
library(clusterProfiler) # Functional enrichment  
library(survival)       # Survival models  
library(survminer)      # Survival visualization  
library(ggplot2)        # Publication-quality plots  
library(ComplexHeatmap) # Advanced heatmaps  
library(patchwork)      # Plot arrangement  
library(meta)           # Meta analysis  
library(metafor)  
library(forestplot)  
library(dmetar)  


# Directory Organization
code/  
├── data_processing/            # Module 1: Data preprocessing and quality control  
├── denovo_Single_Gene/         # Module 2: Single-gene exploratory analysis  
├── enrichment_analysis/        # Module 3: Functional enrichment analysis  
├── get_infiltration/          # Module 4: Immune cell infiltration estimation  
├── meta_analysis/             # Module 5: Multi-cohort meta-analysis  
├── signatures_ICB/            # Module 6: Multi-signature analysis for ICB  
├── signatures_ICB_univariate/ # Module 7: Univariate ICB signature analysis  
├── signatures_TCGA/           # Module 8: Multi-signature analysis for TCGA  
├── signatures_TCGA_univariate/# Module 9: Univariate TCGA signature analysis  
└── summary_figures/           # Module 10: Result visualization and reporting  

# License:
This project is licensed under the MIT License - see the LICENSE file for details.

# Citation:

If you use this code in your research, please cite:  

Zhou et.al. Identification of a Tertiary Lymphoid Structure Signature for Predicting Tumor Outcomes through Transcriptomics Analysis, 2026, Submmited.
