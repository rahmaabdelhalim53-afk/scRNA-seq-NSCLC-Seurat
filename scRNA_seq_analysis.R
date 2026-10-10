## =====================================================================
## NSCLC scRNA-seq analysis (10x 20k NSCLC DTC 3', donor 1) - Seurat v5
## Run from the project folder that contains the .h5 file.
## Outputs: results/plots, results/tables, results/*.rds, sessionInfo.txt
## =====================================================================

library(Seurat)
library(dplyr)
library(ggplot2)

## ---- 0. Setup ------------------------------------------------------
dir.create("results/plots",  recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

save_plot <- function(p, name, w = 8, h = 6)
  ggsave(file.path("results/plots", paste0(name, ".png")), p,
         width = w, height = h, dpi = 300)

qc_feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt")

## ---- 1. Load data and create object ---------------------------------
h5   <- Read10X_h5("20k_NSCLC_DTC_3p_nextgem_donor_1_count_sample_feature_bc_matrix.h5")
cts  <- h5[["Gene Expression"]]          # drop antibody / multiplexing modalities
obj  <- CreateSeuratObject(cts, project = "NSCLC", min.cells = 3, min.features = 200)
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

## ---- 2. QC: before filtering ----------------------------------------
save_plot(VlnPlot(obj, qc_feats, ncol = 3, pt.size = 0),
          "01_QC_violin_before_filtering", 12, 5)
save_plot(FeatureScatter(obj, "nCount_RNA", "nFeature_RNA") + geom_smooth(method = "lm"),
          "02_QC_scatter_before_filtering")

## ---- 3. Filtering ---------------------------------------------------
obj <- subset(obj, subset = nFeature_RNA > 500 & nFeature_RNA < 7500 &
                nCount_RNA < 50000 & percent.mt < 10)

## ---- 4. QC: after filtering (saved now, while identity is still "NSCLC")
save_plot(VlnPlot(obj, qc_feats, ncol = 3, pt.size = 0),
          "03_QC_violin_after_filtering", 12, 5)
save_plot(FeatureScatter(obj, "nCount_RNA", "nFeature_RNA") + geom_smooth(method = "lm"),
          "04_QC_scatter_after_filtering")

## ---- 5. Normalize, HVGs, scale, PCA ---------------------------------
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
save_plot(LabelPoints(VariableFeaturePlot(obj),
                      points = head(VariableFeatures(obj), 10), repel = TRUE),
          "05_variable_features")

obj <- ScaleData(obj, features = rownames(obj))
obj <- RunPCA(obj, features = VariableFeatures(obj))
save_plot(DimHeatmap(obj, dims = 1, cells = 500, balanced = TRUE, fast = FALSE),
          "06_PCA_heatmap_PC1", 10, 8)
save_plot(ElbowPlot(obj), "07_elbow_plot")

## ---- 6. Clustering and UMAP -----------------------------------------
obj <- FindNeighbors(obj, dims = 1:15)
obj <- FindClusters(obj, resolution = c(0.1, 0.3, 0.5, 0.7, 1))
Idents(obj) <- "RNA_snn_res.0.3"         # chosen resolution
obj <- RunUMAP(obj, dims = 1:15)
save_plot(DimPlot(obj, reduction = "umap", label = TRUE), "08_UMAP_clusters_res0.3")

## ---- 7. Markers -----------------------------------------------------
markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top5    <- markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 5)

save_plot(FeaturePlot(obj, c("EPCAM", "PTPRC", "CD3E", "CD68", "COL1A1", "PECAM1")),
          "09_featureplot_markers", 10, 8)

## ---- 8. Annotation (order = cluster 0..6 at res 0.3) -----------------
new.names <- c("T/NK cells", "Tumor 1", "Tumor 2", "Monocytes/Macrophages",
               "Plasma cells", "Fibroblasts", "pDC")
names(new.names) <- levels(obj)
obj <- RenameIdents(obj, new.names)
obj$cell_type <- Idents(obj)

save_plot(DimPlot(obj, reduction = "umap", label = TRUE, repel = TRUE) + NoLegend(),
          "10_UMAP_annotated")
save_plot(DotPlot(obj, features = c("EPCAM", "CD3E", "NKG7", "CD68", "FCN1",
                                    "MZB1", "COL1A1", "DCN", "CLEC4C")) + RotatedAxis(),
          "11_dotplot_markers", 10, 6)

## ---- 9. Save tables, object, session info ---------------------------
write.csv(markers, "results/tables/all_markers.csv", row.names = FALSE)
write.csv(top5,    "results/tables/top5_markers_per_cluster.csv", row.names = FALSE)
write.csv(as.data.frame(table(obj$cell_type)),
          "results/tables/cell_type_counts.csv", row.names = FALSE)

saveRDS(obj, "results/nsclc_seurat_obj.rds")
writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")