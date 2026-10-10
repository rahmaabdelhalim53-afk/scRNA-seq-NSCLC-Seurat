# scRNA-seq-NSCLC-Seurat
# Single-Cell RNA Sequencing Analysis (NSCLC)

This repository contains an end-to-end single-cell RNA sequencing (scRNA-seq) analysis pipeline for Non-Small Cell Lung Cancer (NSCLC) datasets using **R** and **Seurat v5**.

## 📊 Workflow & Steps
1. **Quality Control (QC):** Filtering cells based on `nFeature_RNA`, `nCount_RNA`, and mitochondrial gene percentage (`percent.mt`).
2. **Normalization & Scaling:** Standard log-normalization and identification of highly variable features (HVGs).
3. **Dimensionality Reduction & Clustering:** PCA, Elbow Plot analysis, and graph-based clustering using UMAP.
4. **Cell Type Annotation:** Identifying major cell populations (Tumor cells, T/NK cells, Monocytes/Macrophages, Fibroblasts, Plasma cells, pDCs).

## 📁 Repository Structure
- `results/plots/`: Generated quality control, PCA, and UMAP visualization plots (PDF/PNG).
- `sessionInfo.txt`: R session details and package versions for reproducibility.

## 💻 Software & Environment
- **Language:** R (version 4.6.1)
- **Core Package:** Seurat (v5.5.1)
- See `sessionInfo.txt` for the complete dependency list.
