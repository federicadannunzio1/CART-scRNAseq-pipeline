library(Azimuth)
library(ggplot2)
library(Seurat)
library(xlsx)
library(dplyr)
library(openxlsx)

seurat_list <- readRDS("~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/Data/seurat_obj_list/seurat_samples_sctype_azimuth_pbmc_bonemarrow_clonalvdj_CAR.rds")

Bo_I_raw <- seurat_list$Bo_samples_blood$I
Bo_I_clean <- readRDS("~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/Code/Seurat_analysis/3_data_cleaning/Bo_I_clean.rds")

###################################
# Lanciamo l'annotazione sul reference delle PBMC (pbmcref)
Bo_I_clean <- RunAzimuth(Bo_I_clean, reference = "pbmcref")

# Visualizziamo le nuove etichette (livello 2) sul nuovo UMAP
DimPlot(Bo_I_clean, group.by = "predicted.celltype.l2", label = TRUE) + ggtitle("Bo_I - predicted.celltype.l2")

# Controlliamo quante CAR-T sono rimaste e in che tipi cellulari ricadono
table(Bo_I_clean$predicted.celltype.l2, Bo_I_clean$IS_CAR_ALLIN_scREP)
#########################################

# Visualizza il punteggio di predizione per ogni cluster
VlnPlot(Bo_I_clean, features = "predicted.celltype.l2.score", group.by = "predicted.celltype.l2") + 
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Prediction_score - Azimuth")

# Impostiamo l'identità sui cluster di Azimuth
Idents(Bo_I_clean) <- "predicted.celltype.l2"

# Calcolo dei marcatori
# Nota: logfc.threshold a 0.25 è lo standard, min.pct assicura che il gene sia presente
all_markers <- FindAllMarkers(Bo_I_clean, 
                              only.pos = FALSE, 
                              min.pct = 0.25, 
                              logfc.threshold = 0.25)

# Salviamo i risultati per aprirli in Excel/CSV
write.xlsx(all_markers, "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/Code/Seurat_analysis/3_data_cleaning/Markers_Azimuth_Bo_I.xlsx")

# 1. Filtriamo per p-value adjusted significativo (es. < 0.05)
significant_markers <- all_markers %>%
  filter(p_val_adj < 0.05)

# 2. Creiamo una lista di dataframe, uno per ogni tipo cellulare
cell_types <- unique(significant_markers$cluster)
list_to_save <- list()

for (cell in cell_types) {
  # Filtriamo per il cluster attuale
  cluster_data <- significant_markers %>% filter(cluster == cell)
  
  # Prendiamo i top 30 UP (ordinati per avg_log2FC decrescente)
  top_up <- cluster_data %>%
    filter(avg_log2FC > 0) %>%
    slice_max(order_by = avg_log2FC, n = 30)
  
  # Prendiamo i top 30 DOWN (ordinati per avg_log2FC crescente)
  top_down <- cluster_data %>%
    filter(avg_log2FC < 0) %>%
    slice_min(order_by = avg_log2FC, n = 30)
  
  # Uniamo i due set per questo specifico foglio
  list_to_save[[as.character(cell)]] <- rbind(top_up, top_down)
}

# 3. Salviamo il file Excel con un foglio per ogni tipo cellulare
write.xlsx(list_to_save, 
           file = "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/Code/Seurat_analysis/3_data_cleaning/Top_Markers_Per_CellType_Bo_I.xlsx")

####################################
# Vediamo i primi 3 geni UP e i primi 3 geni DOWN per ogni cluster
top_genes <- all_markers %>%
  group_by(cluster) %>%
  filter(p_val_adj < 0.05) %>%
  # Prendiamo i 3 col log2FC più alto (UP) e i 3 col log2FC più basso (DOWN)
  group_by(cluster) %>%
  slice(c(1:3, (n()-2):n())) %>% 
  arrange(cluster, desc(avg_log2FC))

DoHeatmap(Bo_I_clean, features = top_genes$gene, size = 3) + 
  scale_fill_gradientn(colors = c("blue", "white", "red")) +
  theme(axis.text.y = element_text(size = 7)) +
  ggtitle("Marcatori Up & Down per ogni Cluster Azimuth")

plot_features <- top_genes %>%
  group_by(cluster) %>%
  slice(c(1:2, (n()-1):n())) %>%
  pull(gene) %>%
  unique()

DotPlot(Bo_I_clean, features = plot_features, cols = c("blue", "red")) + 
  RotatedAxis() + 
  theme(axis.text.x = element_text(size = 8)) +
  ggtitle("Markers per cluster")
# 1. Definizione dei marcatori fondamentali per lineage
# Li dividiamo per categorie così è più facile leggere il grafico
canonical_markers <- list(
  'T-cell' = c("CD3D", "CD3E", "TRAC"),
  'CD4-T'  = c("CD4", "IL7R", "CCR7"),
  'CD8-T'  = c("CD8A", "CD8B", "GZMK", "GZMB"),
  'NK'     = c("NCAM1", "KLRB1", "NKG7", "GNLY"),
  'B-cell' = c("CD19", "MS4A1", "CD79A", "IGKC"),
  'Mono-CD14' = c("CD14", "LYZ", "S100A8"),
  'Mono-CD16' = c("FCGR3A", "MS4A7"),
  'DC'     = c("LILRA4", "CD1C", "CLEC4C"),
  'Plat/Mega' = c("PPBP", "PF4"),
  'Prolif' = c("MKI67", "TOP2A", "UHRF2")
)

# 2. Creazione del DotPlot
# Usiamo 'unlist' per passare la lista piatta al comando
DotPlot(Bo_I_clean, 
        features = canonical_markers, 
        group.by = "predicted.celltype.l2", 
        dot.scale = 6) + 
  RotatedAxis() + 
  scale_colour_gradient2(low = "blue", mid = "white", high = "red") +
  theme(axis.text.x = element_text(size = 9, face = "italic")) +
  ggtitle("Validation with known markers")

#cerco le cd8
FeaturePlot(Bo_I_clean, 
            features = c("CD4", "CD8A"), 
            ncol = 2, 
            min.cutoff = "q9", 
            cols = c("lightgrey", "red"),
            label = TRUE) # Aggiungiamo le label per orientarci tra i cluster

# Se vuoi vedere se le CAR-T sono finite in cluster specifici
FeaturePlot(Bo_I_clean, features = "IS_CAR_ALLIN_scREP", label = TRUE)
# Vediamo la distribuzione delle CAR-T tra i cluster definiti da Azimuth
ggplot(Bo_I_clean@meta.data, aes(x = predicted.celltype.l2, fill = IS_CAR_ALLIN_scREP)) +
  geom_bar(position = "fill") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Proporzione") +
  ggtitle("Distribuzione CAR-T nei cluster Azimuth")
