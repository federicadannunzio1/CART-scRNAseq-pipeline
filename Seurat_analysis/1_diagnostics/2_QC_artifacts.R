rm(list = ls())

library(Seurat)
library(dplyr)
library(xlsx)
library(patchwork)

seurat_obj1 <- seurat_obj$Me_samples_bone$AB

# drivers of the first 2 PCA
VizDimLoadings(seurat_obj1, dims = 1:2, reduction = "pca", ncol = 2)
#################################
# mt and rb genes percentage, nfeatures, ncounts and malat1 violin plots

seurat_obj1[["percent.mt"]] <- PercentageFeatureSet(seurat_obj1, pattern = "^MT-")
seurat_obj1[["percent.rb"]] <- PercentageFeatureSet(seurat_obj1, pattern = "^RP[SL]")

VlnPlot(seurat_obj1, 
        features = c("nCount_RNA","nFeature_RNA", "percent.mt", "percent.rb", "MALAT1"), 
        group.by = "predicted.celltype.l2_bonemarrow", # check it
        pt.size = 0, 
        ncol = 3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  plot_annotation(title = "QC Metrics per Cluster")

#################################
# tables

qc_stats <- seurat_obj1@meta.data %>%
  group_by(predicted.celltype.l3_pbmc) %>%
  summarise(
    cell_count = n(),
    median_nFeature = median(nFeature_RNA),
    median_nCount = median(nCount_RNA),
    mean_MALAT1 = mean(seurat_obj1@assays$RNA@features["MALAT1", ])
  ) %>%
  arrange(median_nFeature)

print(qc_stats)

qc_stats <- as.data.frame(qc_stats)

write.xlsx(qc_stats, "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/Code/Seurat_analysis/1_make_markers_featureplots/qc_stats_table.xlsx" )
