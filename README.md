# Single-Cell RNA-Seq Quality Control and Clustering Pipeline

A reproducible R pipeline using **Seurat** for processing, quality filtering, dimensional reduction (PCA/UMAP), graph-based clustering, and marker gene identification on single-cell RNA sequencing data.

## 🔬 Pipeline Workflow Steps:
1. **Quality Control (QC):** Filters cells based on library complexity (`nFeature_RNA`, `nCount_RNA`) and mitochondrial gene expression percentages (`percent.mt`) to remove dying cells and doublets.
2. **Normalization & Feature Selection:** Log-normalization and identification of highly variable genes (HVGs) to capture biological heterogeneity.
3. **Dimensionality Reduction:** Linear dimensionality reduction via PCA followed by non-linear reduction via UMAP.
4. **Graph-Based Clustering:** Modularity optimization (Louvain algorithm) to group cells into distinct functional clusters.
5. **Biomarker Discovery:** Differential expression analysis (`FindAllMarkers`) to identify cluster-specific marker genes.

## 🛠️ Requirements & Dependencies:
* R (>= 4.0)
* Seurat
* ggplot2
* dplyr

## 📁 Repository Structure:
* `scRNA_analysis.R`: The complete R script containing the full analytical workflow.
