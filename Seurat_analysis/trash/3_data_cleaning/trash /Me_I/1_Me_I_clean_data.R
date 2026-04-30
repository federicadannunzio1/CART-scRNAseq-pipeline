rm(list = ls())

library(Seurat)
library(ggplot2)
library(dplyr)
library(xlsx)

# --- ANALISI PER Me_I ---
seurat_list <- readRDS("~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/Data/seurat_obj_list/seurat_samples_sctype_azimuth_pbmc_bonemarrow_clonalvdj_CAR.rds")

# 1. Filtering (Soglie standard per i campioni Initial Product)
Me_I_raw <- seurat_list$Me_samples_bone$I
Me_I_clean <- subset(Me_I_raw, subset = nFeature_RNA > 800 & percent.mt < 7)

# 2. Rimozione dell'artefatto tecnico dnT
Me_I_clean <- subset(Me_I_clean, subset = azimuth_class != "dnT")

# 3. Normalizzazione e Selezione Variable Features
Me_I_clean <- NormalizeData(Me_I_clean)
Me_I_clean <- FindVariableFeatures(Me_I_clean, selection.method = "vst", nfeatures = 2000)

# 4. Esclusione Geni di Stress dalla PCA
stress_genes <- c("MALAT1", "NEAT1", grep("^MT-", rownames(Me_I_clean), value = TRUE))
VariableFeatures(Me_I_clean) <- setdiff(VariableFeatures(Me_I_clean), stress_genes)

# 5. Scaling (Regressione mitocondriale) e PCA a 30 PC
Me_I_clean <- ScaleData(Me_I_clean, vars.to.regress = "percent.mt")
Me_I_clean <- RunPCA(Me_I_clean, npcs = 30, verbose = FALSE)

# Check rapido (Opzionale)
ElbowPlot(Me_I_clean)

# 6. UMAP e Clustering (Resolution 1.0 per gestire l'omogeneità dei CD4)
Me_I_clean <- RunUMAP(Me_I_clean, dims = 1:30)
Me_I_clean <- FindNeighbors(Me_I_clean, dims = 1:30)
Me_I_clean <- FindClusters(Me_I_clean, resolution = 0.5)

# Visualizzazione distribuzione numerica
n_cells <- table(Me_I_clean$seurat_clusters)

bp <- barplot(
  n_cells,
  main = "Me_I - Number of cells in each cluster",
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
p1 <- DimPlot(Me_I_clean, reduction = "umap", label = TRUE) + ggtitle("Me_I processed")
p2 <- DimPlot(Me_I_raw, reduction = "umap", label = TRUE) + ggtitle("Me_I_not_processed")
p1 + p2


# 7. Salvataggio intermedio e Generazione Marcatori
saveRDS(Me_I_clean, "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/Code/Seurat_analysis/3_data_cleaning/Me_I/Me_I_clean_pre_annotation.rds")
