rm(list = ls())

library(Seurat)
library(ggplot2)
# --- ANALISI PER Ca_I ---
seurat_list <- readRDS("~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/Data/seurat_obj_list/seurat_samples_sctype_azimuth_pbmc_bonemarrow_clonalvdj_CAR.rds")

# 1. Filtering 
# Utilizziamo le stesse soglie: nFeature > 800 e percent.mt < 7
Ca_I_raw <- seurat_list$Ca_samples_blood$I
Ca_I_filtered <- subset(Ca_I_raw, subset = nFeature_RNA > 800 & percent.mt < 7)

# 2. Rimozione dell'artefatto tecnico dnT
# Escludiamo i detriti identificati da Azimuth
Ca_I_clean <- subset(Ca_I_filtered, subset = azimuth_class != "dnT")

# 3. Normalizzazione e Selezione dei Variable Features
Ca_I_clean <- NormalizeData(Ca_I_clean)
Ca_I_clean <- FindVariableFeatures(Ca_I_clean, selection.method = "vst", nfeatures = 2000)

# 4. Esclusione Geni di Stress dalla PCA
stress_genes <- c("MALAT1", "NEAT1", grep("^MT-", rownames(Ca_I_clean), value = TRUE))
VariableFeatures(Ca_I_clean) <- setdiff(VariableFeatures(Ca_I_clean), stress_genes)

# 5. Scaling e Dimensionality Reduction
# Regrediamo percent.mt per mitigare l'influenza della sofferenza cellulare sulle distanze
Ca_I_clean <- ScaleData(Ca_I_clean, vars.to.regress = "percent.mt")

Ca_I_clean <- RunPCA(Ca_I_clean, npcs = 30, verbose = FALSE)
ElbowPlot(Ca_I_clean)
VizDimLoadings(Ca_I_clean, dims = 1:2, reduction = "pca")

Ca_I_clean <- RunUMAP(Ca_I_clean, dims = 1:30)

# 6. Clustering
Ca_I_clean <- FindNeighbors(Ca_I_clean, dims = 1:30)
Ca_I_clean <- FindClusters(Ca_I_clean, resolution = 0.7)
n_cells <- table(Ca_I_clean$seurat_clusters)

bp <- barplot(
  n_cells,
  main = "Ca_I - Number of cells in each cluster",
  ylim = c(0, max(n_cells) * 1.25),
  col = colorRampPalette(c(
    "salmon",      # 0
    "goldenrod",   # 1
    "limegreen",   # 2
    "turquoise3",  # 3
    "cornflowerblue", # 4
    "hotpink"      # 5
  ))(length(n_cells))
)

text(
  x = bp,
  y = n_cells,
  labels = n_cells,
  pos = 3,      # sopra la barra
  cex = 0.8
)
# 7. Verifica Finale e Salvataggio
p1_ca <- DimPlot(Ca_I_clean, reduction = "umap", label = TRUE) + ggtitle("Ca_I processed")
p2_ca <- DimPlot(Ca_I_raw, reduction = "umap", label = TRUE) + ggtitle("Ca_I_not_processed")
p1_ca + p2_ca

saveRDS(Ca_I_clean, "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/Code/Seurat_analysis/3_data_cleaning/Ca_I/Ca_I_clean.rds")

###########################################
###########################################
# ANALISI QC SULLE CAR-T PERSE (Ca_I)
###########################################

# Controllo quante CAR T sono rimaste
print(table(Ca_I_clean$CAR))

# Identificazione barcode CAR+
car_raw_ca <- WhichCells(Ca_I_raw, expression = IS_CAR_ALLIN_scREP == "YES")
car_clean_ca <- WhichCells(Ca_I_clean, expression = IS_CAR_ALLIN_scREP == "YES")
car_lost_ca <- setdiff(car_raw_ca, car_clean_ca)

# Creazione colonna QC nei metadati dell'oggetto raw
Ca_I_raw$QC_Comparison <- "Other"
Ca_I_raw@meta.data[car_lost_ca, "QC_Comparison"] <- "Lost"
Ca_I_raw@meta.data[car_clean_ca, "QC_Comparison"] <- "Kept"

# Subset per visualizzazione
subset_CAR_ca <- subset(Ca_I_raw, subset = QC_Comparison %in% c("Lost", "Kept"))

# Generazione Violin Plot
VlnPlot(subset_CAR_ca, 
        features = c("nFeature_RNA", "percent.mt"), 
        group.by = "QC_Comparison", 
        cols = c("Lost" = "red", "Kept" = "green")) +
  patchwork::plot_annotation(title = "QC Metrics: CAR-T cells Kept vs Lost (Ca_I)")

