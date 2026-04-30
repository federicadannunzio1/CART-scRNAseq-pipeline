rm(list = ls())

library(Seurat)
library(dplyr)
library(ggplot2)

seurat_list <- readRDS("~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/Data/seurat_obj_list/seurat_samples_sctype_azimuth_pbmc_bonemarrow_clonalvdj_CAR.rds")

# 1. Filtering 
# Alzoo le soglie per isolare solo le cellule con buona complessità genica
Bo_I_raw <- seurat_list$Bo_samples_blood$I
Bo_I_filtered <- subset(Bo_I_raw, subset = nFeature_RNA > 800 & percent.mt < 7)

# 2. Rimozione dell'artefatto tecnico dnT (se presente nei metadati)
# Le dnT in questo campione sono detriti
Bo_I_clean <- subset(Bo_I_filtered, subset = azimuth_class != "dnT")

# 3. Normalizzazione e Selezione dei Variable Features
Bo_I_clean <- NormalizeData(Bo_I_clean)
Bo_I_clean <- FindVariableFeatures(Bo_I_clean, selection.method = "vst", nfeatures = 2000)

# 4. Esclusione Geni di Stress dalla PCA
# Identifichiamo i geni che guidano il bias tecnico
stress_genes <- c("MALAT1", "NEAT1", grep("^MT-", rownames(Bo_I_clean), value = TRUE))
VariableFeatures(Bo_I_clean) <- setdiff(VariableFeatures(Bo_I_clean), stress_genes)

# 5. Scaling e Dimensionality Reduction
# Regrediamo anche la percentuale mitocondriale per pulire ulteriormente il dato
Bo_I_clean <- ScaleData(Bo_I_clean, vars.to.regress = "percent.mt")

Bo_I_clean <- RunPCA(Bo_I_clean, npcs = 30, verbose = FALSE)
VizDimLoadings(Bo_I_clean, dims = 1:2, reduction = "pca")
DimHeatmap(Bo_I_clean, dims = 1, cells = 500, balanced = TRUE)

Bo_I_clean <- RunUMAP(Bo_I_clean, dims = 1:20)

# 6. Clustering
Bo_I_clean <- FindNeighbors(Bo_I_clean, dims = 1:20)
Bo_I_clean <- FindClusters(Bo_I_clean, resolution = 0.7)

# 7. Verifica Finale
p1 <- DimPlot(Bo_I_clean, reduction = "umap", label = TRUE) + ggtitle("S429_I processed")
p2 <- DimPlot(Bo_I_raw, reduction = "umap", label = TRUE) + ggtitle("S429_I_not_processed")
p1+p2

saveRDS(Bo_I_clean, "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/Code/Seurat_analysis/3_data_cleaning/Bo_I_clean.rds")

##############################
table(Bo_I_clean$CAR) # Controllo quante CAR T sono rimaste [cite: 336]

car_raw <- WhichCells(Bo_I_raw, expression = IS_CAR_ALLIN_scREP == "YES")
car_clean <- WhichCells(Bo_I_clean, expression = IS_CAR_ALLIN_scREP == "YES")
car_lost <- setdiff(car_raw, car_clean)
# 1. Creiamo una nuova colonna nei metadati inizializzata a "Other"
Bo_I_raw$QC_Comparison <- "Other"

# 2. Assegniamo le etichette basandoci sui barcode che hai estratto
Bo_I_raw@meta.data[car_lost, "QC_Comparison"] <- "Lost"
Bo_I_raw@meta.data[car_clean, "QC_Comparison"] <- "Kept"

# 3. Filtriamo l'oggetto solo per vedere le CAR-T (escludendo "Other")
subset_CAR <- subset(Bo_I_raw, subset = QC_Comparison %in% c("Lost", "Kept"))

# 4. Generiamo il Violin Plot
VlnPlot(subset_CAR, 
        features = c("nFeature_RNA", "percent.mt"), 
        group.by = "QC_Comparison", 
        cols = c("Lost" = "red", "Kept" = "green")) 

