rm(list = ls())

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(openxlsx)

# ==========================================
# 1. SETUP E CARICAMENTO DATI
# ==========================================
base_dir <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/CART/Code/Seurat_analysis/3_data_cleaning/"
# SUGGERIMENTO: se il caricamento da Drive dà problemi, usa un percorso locale come "~/Desktop/seurat_list_raw.rds"
seurat_list_path <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/CART/Data/seurat_obj_list/seurat_samples_sctype_azimuth_pbmc_bonemarrow_clonalvdj_CAR.rds"

cat("Caricamento dataset globale...\n")
seurat_list <- readRDS(seurat_list_path)

# Estrazione dei 3 campioni raw
raw_samples <- list(
  Bo_I = seurat_list$Bo_samples_blood$I,
  Ca_I = seurat_list$Ca_samples_blood$I,
  Me_I = seurat_list$Me_samples_bone$I
)

# Definizione dei parametri specifici per ogni campione
params <- list(
  Bo_I = list(dims = 1:20, res = 0.7),
  Ca_I = list(dims = 1:30, res = 0.7),
  Me_I = list(dims = 1:30, res = 0.5)
)

# ==========================================
# 2. FUNZIONE DI PROCESSING AUTOMATIZZATA
# ==========================================
process_and_find_markers <- function(sample_name, raw_obj, p) {
  cat(paste0("\n-----------------------------------\n"))
  cat(paste0("Inizio processamento per: ", sample_name, "\n"))
  cat(paste0("Parametri -> Dims: 1:", max(p$dims), " | Resolution: ", p$res, "\n"))
  
  # 1-2. Filtering & Rimozione dnT
  clean_obj <- subset(raw_obj, subset = nFeature_RNA > 800 & percent.mt < 7 & azimuth_class != "dnT")
  
  # 3. Normalizzazione e Variable Features
  clean_obj <- NormalizeData(clean_obj, verbose = FALSE)
  clean_obj <- FindVariableFeatures(clean_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  
  # 4. Esclusione Geni di Stress
  stress_genes <- c("MALAT1", "NEAT1", grep("^MT-", rownames(clean_obj), value = TRUE))
  VariableFeatures(clean_obj) <- setdiff(VariableFeatures(clean_obj), stress_genes)
  
  # 5. Scaling e PCA
  clean_obj <- ScaleData(clean_obj, vars.to.regress = "percent.mt", verbose = FALSE)
  clean_obj <- RunPCA(clean_obj, npcs = 30, verbose = FALSE)
  
  # 6. UMAP e Clustering (con parametri specifici)
  clean_obj <- RunUMAP(clean_obj, dims = p$dims, verbose = FALSE)
  clean_obj <- FindNeighbors(clean_obj, dims = p$dims, verbose = FALSE)
  clean_obj <- FindClusters(clean_obj, resolution = p$res, verbose = FALSE)
  
  # 7. Plot di Controllo (UMAP Pre vs Post)
  p1 <- DimPlot(clean_obj, reduction = "umap", label = TRUE) + ggtitle(paste0(sample_name, " processed"))
  p2 <- DimPlot(raw_obj, reduction = "umap", label = TRUE) + ggtitle(paste0(sample_name, " not_processed"))
  print(p1 + p2)
  
  # 8. Creazione Directory per i file Excel
  out_dir <- paste0(base_dir, sample_name, "/")
  if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # ==========================================
  # 9. Calcolo e Salvataggio Marcatori
  # ==========================================
  cat("Calcolo dei marcatori (FindAllMarkers)...\n")
  Idents(clean_obj) <- "seurat_clusters"
  markers <- FindAllMarkers(clean_obj, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)
  
  # Salvataggio file Excel completo
  write.xlsx(markers, paste0(out_dir, "Markers_Res", gsub("\\.", "", p$res), "_", sample_name, ".xlsx"))
  
  # Filtraggio Top 30 e file Excel multi-foglio
  sig_markers <- markers %>% filter(p_val_adj < 0.05)
  cell_types <- unique(sig_markers$cluster)
  list_to_save <- list()
  
  for (cell in cell_types) {
    cluster_data <- sig_markers %>% filter(cluster == cell)
    top_up <- cluster_data %>% filter(avg_log2FC > 0) %>% slice_max(order_by = avg_log2FC, n = 30)
    top_down <- cluster_data %>% filter(avg_log2FC < 0) %>% slice_min(order_by = avg_log2FC, n = 30)
    list_to_save[[as.character(cell)]] <- rbind(top_up, top_down)
  }
  
  write.xlsx(list_to_save, paste0(out_dir, "Top_Markers_Per_CellType_", sample_name, ".xlsx"))
  cat("Estrazione marcatori completata.\n")
  
  # Restituiamo solo l'oggetto pulito da inserire nella lista
  return(clean_obj)
}

# ==========================================
# 3. ESECUZIONE DEL LOOP SU TUTTI I CAMPIONI
# ==========================================
processed_samples <- list()

for (s_name in names(raw_samples)) {
  processed_samples[[s_name]] <- process_and_find_markers(
    sample_name = s_name, 
    raw_obj = raw_samples[[s_name]], 
    p = params[[s_name]]
  )
}

# ==========================================
# 4. SALVATAGGIO DELLA LISTA UNICA
# ==========================================
cat("\nPipeline completata! Salvataggio della lista unica in corso...\n")
list_save_path <- paste0(base_dir, "all_samples_clean_pre_annotation.rds")

saveRDS(processed_samples, list_save_path)
cat("Lista salvata con successo in:", list_save_path, "\n")