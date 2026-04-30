library(Azimuth)
library(ggplot2)
library(Seurat)
library(xlsx)
library(dplyr)
library(openxlsx)
library(purrr)

# 1. Caricamento del dato pulito (Precedentemente salvato con res 0.7)
Me_I_clean <- readRDS("~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/Code/Seurat_analysis/3_data_cleaning/Me_I/Me_I_clean_pre_annotation.rds")

# Assicuriamoci di essere sui cluster a res 0.7 (quelli che hanno prodotto 8 cluster)
Idents(Me_I_clean) <- "seurat_clusters"

# Calcolo dei marcatori per ogni cluster (essenziale per l'annotazione manuale)
all_markers_me <- FindAllMarkers(Me_I_clean, 
                                 only.pos = FALSE, 
                                 min.pct = 0.25, 
                                 logfc.threshold = 0.25)

# Salvataggio file Excel completo
write.xlsx(all_markers_me, "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/Code/Seurat_analysis/3_data_cleaning/Me_I/Markers_Res07_Me_I.xlsx")

# Filtraggio e creazione file Excel multi-foglio (Top 30 UP e DOWN)
significant_markers_me <- all_markers_me %>% filter(p_val_adj < 0.05)
cell_types_me <- unique(significant_markers_me$cluster)
list_to_save_me <- list()

for (cell in cell_types_me) {
  cluster_data <- significant_markers_me %>% filter(cluster == cell)
  
  top_up <- cluster_data %>% filter(avg_log2FC > 0) %>% slice_max(order_by = avg_log2FC, n = 30)
  top_down <- cluster_data %>% filter(avg_log2FC < 0) %>% slice_min(order_by = avg_log2FC, n = 30)
  
  list_to_save_me[[as.character(cell)]] <- rbind(top_up, top_down)
}

write.xlsx(list_to_save_me, 
           file = "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/Code/Seurat_analysis/3_data_cleaning/Me_I/Top_Markers_Per_CellType_Me_I.xlsx")

# 1. Definizione delle etichette per i 5 cluster
new_labels_me <- c(
  "0" = "CD8_Cytotoxic_Effector",
  "1" = "Proliferating_T_S_Phase",
  "2" = "CD4_Naive_Memory_Treg",
  "3" = "Proliferating_T_G2M_Phase",
  "4" = "Proliferating_CD8_Cytotoxic"
)

# 2. Applicazione all'oggetto
Me_I_clean <- RenameIdents(Me_I_clean, new_labels_me)
Me_I_clean$manual_annotation <- Idents(Me_I_clean)

# 3. Verifica visiva finale
DimPlot(Me_I_clean, reduction = "umap", label = TRUE, repel = TRUE) + 
  ggtitle("Me_I")

# 4. Salvataggio oggetto finale
saveRDS(Me_I_clean, "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/Code/Seurat_analysis/3_data_cleaning/Me_I/Me_I_clean_annotated_final.rds")
