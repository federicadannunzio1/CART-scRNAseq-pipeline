library(Seurat)
library(xlsx)
library(dplyr)

seurat_list <- readRDS("~/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/caruana-project/Data/seurat_obj_list/seurat_samples_sctype_azimuth_pbmc_bonemarrow_clonalvdj_CAR.rds")

sample <- seurat_list$Bo_samples_blood$I
sample_name <- "Bo_I"
############################
counts <- table(sample@meta.data$predicted.celltype.l2_pbmc)

bp <- barplot(counts, 
              las = 2, 
              ylim = c(0, max(counts) * 1.1), 
              main = "Cellule per Tipo",
              col = "lightblue")


text(x = bp, y = counts, label = counts, pos = 3, cex = 0.8, col = "black")
###########################
DimPlot(sample, group.by = "predicted.celltype.l3_pbmc")

Idents(sample) <- "predicted.celltype.l3_pbmc"

markers_for_cluster <- FindAllMarkers(
  sample,
  assay = "RNA",        
  only.pos = F,
  min.pct = 0.25,
  logfc.threshold = 0.1
)

markers_for_cluster <- markers_for_cluster[markers_for_cluster$p_val_adj <= 0.05,]
markers_for_cluster <- markers_for_cluster[order(abs(markers_for_cluster$avg_log2FC), decreasing = TRUE),]
markers_for_cluster$direction <- ifelse(markers_for_cluster$avg_log2FC > 0, "UP", "DOWN")

markers_for_cluster <- markers_for_cluster %>% 
  select(gene, cluster, avg_log2FC, p_val, p_val_adj, everything())

write.xlsx(markers_for_cluster, paste0("~/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/caruana-project/Code/Seurat_analysis/2_find_markers/", sample_name , "_markers_per_cluster.xlsx"), row.names = F)
