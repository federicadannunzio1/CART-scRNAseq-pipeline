library(Azimuth)
library(ggplot2)
library(Seurat)
library(xlsx)
library(dplyr)
library(openxlsx)
library(purrr)

# 1. Caricamento del dato pulito (Precedentemente salvato con res 1.0)
Bo_I_clean <- readRDS("~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/Code/Seurat_analysis/3_data_cleaning/Bo_I/Bo_I_clean.rds")

# Assicuriamoci di essere sui cluster a res 1.0
Idents(Bo_I_clean) <- "seurat_clusters"

# Calcolo dei marcatori per ogni cluster (necessario per l'annotazione di precisione)
all_markers_bo <- FindAllMarkers(Bo_I_clean, 
                                 only.pos = FALSE, 
                                 min.pct = 0.25, 
                                 logfc.threshold = 0.25)

# Salvataggio file Excel completo
write.xlsx(all_markers_bo, "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/Code/Seurat_analysis/3_data_cleaning/Bo_I/Markers_Res1_Bo_I.xlsx")

# Filtraggio e creazione file Excel multi-foglio (Top 30 UP e DOWN)
significant_markers_bo <- all_markers_bo %>% filter(p_val_adj < 0.05)
cell_types_bo <- unique(significant_markers_bo$cluster)
list_to_save_bo <- list()

for (cell in cell_types_bo) {
  cluster_data <- significant_markers_bo %>% filter(cluster == cell)
  
  top_up <- cluster_data %>% filter(avg_log2FC > 0) %>% slice_max(order_by = avg_log2FC, n = 30)
  top_down <- cluster_data %>% filter(avg_log2FC < 0) %>% slice_min(order_by = avg_log2FC, n = 30)
  
  list_to_save_bo[[as.character(cell)]] <- rbind(top_up, top_down)
}

write.xlsx(list_to_save_bo, 
           file = "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/Code/Seurat_analysis/3_data_cleaning/Bo_I/Top_Markers_Per_CellType_Bo_I.xlsx")

###################################
# 2. Annotazione Manuale Consolidata
###################################
# In Bo_I, i cluster 0-3 sono quasi identici (CD4 in fortissima proliferazione)
# 1. Definizione dei nuovi nomi basata sull'evidenza della tabella
new_labels_bo <- c(
  "0" = "Proliferating_CD4_T",
  "1" = "Proliferating_CD4_T",
  "2" = "Proliferating_CD4_T",
  "3" = "Proliferating_CD4_T",
  "4" = "Effector_CD8_T",
  "5" = "Activated_CD4_T",
  "6" = "Monocyte_Contam."
)

# 2. Applicazione all'oggetto Seurat
Bo_I_clean <- RenameIdents(Bo_I_clean, new_labels_bo)
Bo_I_clean$manual_annotation <- Idents(Bo_I_clean)

# 3. Verifica Visiva
p_final_bo <- DimPlot(Bo_I_clean, reduction = "umap", label = TRUE, repel = TRUE) + 
  ggtitle("Bo_I")
print(p_final_bo)

# 4. Salvataggio Oggetto Annotato
saveRDS(Bo_I_clean, "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/Code/Seurat_analysis/3_data_cleaning/Bo_I/Bo_I_clean_with_annotation.rds")

