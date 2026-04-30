library(Azimuth)
library(ggplot2)
library(Seurat)
library(xlsx)
library(dplyr)
library(openxlsx)
library(purrr)

# 1. Caricamento del dato pulito (precedentemente salvato)
Ca_I_clean <- readRDS("~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/Code/Seurat_analysis/3_data_cleaning/Ca_I_clean.rds")

# Calcolo dei marcatori per ogni cluster
all_markers_ca <- FindAllMarkers(Ca_I_clean, 
                                 only.pos = FALSE, 
                                 min.pct = 0.25, 
                                 logfc.threshold = 0.25)

# Salvataggio file Excel completo
write.xlsx(all_markers_ca, "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/Code/Seurat_analysis/3_data_cleaning/Ca_I/Markers_Azimuth_Ca_I.xlsx")

# Filtraggio e creazione file Excel multi-foglio (Top 30 UP e DOWN)
significant_markers_ca <- all_markers_ca %>% filter(p_val_adj < 0.05)
cell_types_ca <- unique(significant_markers_ca$cluster)
list_to_save_ca <- list()

for (cell in cell_types_ca) {
  cluster_data <- significant_markers_ca %>% filter(cluster == cell)
  
  top_up <- cluster_data %>% filter(avg_log2FC > 0) %>% slice_max(order_by = avg_log2FC, n = 30)
  top_down <- cluster_data %>% filter(avg_log2FC < 0) %>% slice_min(order_by = avg_log2FC, n = 30)
  
  list_to_save_ca[[as.character(cell)]] <- rbind(top_up, top_down)
}

write.xlsx(list_to_save_ca, 
           file = "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/Code/Seurat_analysis/3_data_cleaning/Ca_I/Top_Markers_Per_CellType_Ca_I.xlsx")

# Definizione dei nomi basata sui marker della tabella
new_labels <- c(
  "0" = "Proliferating_T",
  "1" = "CD4_Naive",
  "2" = "CD4_Memory",
  "3" = "CD8_Cytotoxic",
  "4" = "CD4_Effector",
  "5" = "Regulatory_T"
)

# Applicazione all'oggetto Seurat
Ca_I_clean <- RenameIdents(Ca_I_clean, new_labels)
Ca_I_clean$manual_annotation <- Idents(Ca_I_clean)

# Verifica visiva
DimPlot(Ca_I_clean, reduction = "umap", label = TRUE) + 
  ggtitle("Ca_I")

saveRDS(Ca_I_clean, "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/Code/Seurat_analysis/3_data_cleaning/Ca_I/Ca_I_clean_with_annotation.rds")

####################################
##################################
# 2. Annotazione con Azimuth (PBMC Reference) NOT USED
###################################

# Lanciamo l'annotazione sul reference delle PBMC (pbmcref)
Ca_I_clean <- RunAzimuth(Ca_I_clean, reference = "pbmcref")

# Visualizziamo le nuove etichette (livello 2) sul nuovo UMAP
DimPlot(Ca_I_clean, group.by = "predicted.celltype.l2", label = TRUE) + 
  ggtitle("Ca_I - predicted.celltype.l2")

# Controlliamo quante CAR-T sono rimaste e in che tipi cellulari ricadono
print(table(Ca_I_clean$predicted.celltype.l2, Ca_I_clean$IS_CAR_ALLIN_scREP))

# Visualizza il punteggio di predizione per ogni cluster
VlnPlot(Ca_I_clean, features = "predicted.celltype.l2.score", group.by = "predicted.celltype.l2") + 
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Prediction_score - Azimuth (Ca_I)")

# Visualizza il livello 1 per vedere se è più pulito
DimPlot(Ca_I_clean, group.by = "predicted.celltype.l1", label = TRUE) + 
  ggtitle("Azimuth Livello 1: Macro-popolazioni")
###################################
# 3. Analisi dei Marcatori (Differential Expression)
###################################

# Impostiamo l'identità sui cluster di Azimuth per il calcolo dei markers
#Idents(Ca_I_clean) <- "predicted.celltype.l2" #no
