# =====================================================================
# Project: Single-Cell RNA-Seq Quality Control and Clustering Pipeline
# Author: Nikoleta Raptopoulou
# Description: A clean, reproducible Seurat workflow for scRNA-seq analysis.
# =====================================================================

# 1. Load Required Libraries
library(Seurat)
library(ggplot2)
library(dplyr)

# 2. Create Output Directory for Results
if (!dir.exists("output")) {
  dir.create("output")
}

# 3. Load Data & Initialize Seurat Object
# (Note: Ensure your working directory points to the folder containing the 10X matrices)
pbmc.data <- Read10X(data.dir = "../data/pbmc3k/filtered_gene_bc_matrices/hg19/")

pbmc <- CreateSeuratObject(counts = pbmc.data, 
                           project = "pbmc3k", 
                           min.cells = 3, 
                           min.features = 200)

# 4. Quality Control (QC) & Filtering
# Calculate mitochondrial gene percentage (high % indicates dying/damaged cells)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Save QC Violin Plots
pdf("output/qc_violin_plots.pdf", width = 8, height = 5)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# Filter cells based on gene counts and mitochondrial thresholds
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# 5. Normalization and Feature Selection
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Scale data while regressing out mitochondrial variation
pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")

# 6. Dimensionality Reduction & Clustering
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Run UMAP for 2D visualization
pbmc <- RunUMAP(pbmc, dims = 1:10)

# Save UMAP Cluster Plot
umap_plot <- DimPlot(pbmc, reduction = "umap", label = TRUE) + NoLegend()
ggsave("output/umap_clusters.pdf", plot = umap_plot, width = 6, height = 5)

# 7. Biomarker Discovery (Find Marker Genes per Cluster)
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

# Export top markers to a CSV file
write.csv(top10, "output/top_cluster_markers.csv", row.names = FALSE)

# Generate and save Heatmap for top markers
pdf("output/marker_heatmap.pdf", width = 10, height = 8)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()
dev.off()

cat("Single-cell analysis pipeline executed successfully! Results saved in 'output/'
