library(Seurat)

# Funzione per pulire un singolo oggetto Seurat
clean_seurat_object <- function(obj) {
  
  # 1. Consolidamento dei layer (Fix tipico di Seurat v5)
  # Rimuove la frammentazione mantenendo solo 'counts', 'data' e 'scale.data'
  if (length(grep("^counts\\.", Layers(obj))) > 0) {
    obj <- JoinLayers(obj)
  }
  
  # 2. Pulizia dei metadati
  # Definiamo solo le colonne essenziali usate per plotting e downstream analysis
  cols_to_keep <- c(
    "orig.ident", 
    "nCount_RNA", 
    "nFeature_RNA", 
    "percent.mt",
    "sample", 
    "group", 
    "seurat_clusters", 
    "cell_type",
    "CAR", 
    "derived_CAR", 
    "clonal", 
    "barcode_screpertoire",
    "IS_CONSERVED_scRepertoire", 
    "IS_CAR_ALLIN_scREP"
  )
  
  # Filtriamo il dataframe dei metadati
  obj@meta.data <- obj@meta.data[, intersect(colnames(obj@meta.data), cols_to_keep)]
  
  # 3. Rimozione dei grafi SNN e NN
  # Liberiamo una quantità enorme di memoria
  obj@graphs <- list()
  
  # 4. Rimozione dello storico dei comandi
  obj@commands <- list()
  
  # 5. Rimozione di scale.data (OPZIONALE MA RACCOMANDATO)
  # Decommenta la riga sotto se NON devi fare DoHeatmap() in futuro. 
  # scale.data pesa moltissimo. FeaturePlot, VlnPlot e DotPlot usano il layer 'data'.
  # obj[["RNA"]]$scale.data <- NULL
  
  return(obj)
}

# Applichiamo la funzione a tutta la tua lista di oggetti (I e AB)
# Sostituisci 'all_samples' con il nome reale della tua lista caricata dall'ambiente
cleaned_seurat_list <- lapply(all_samples, clean_seurat_object)

# Controllo post-pulizia per verificare il risultato
cat("\nStruttura post-pulizia:\n")
for (nm in names(cleaned_seurat_list)) {
  cat(sprintf("  %-20s | Metadati: %2d colonne | Memoria ridotta\n",
              nm, ncol(cleaned_seurat_list[[nm]]@meta.data)))
}

saveRDS(cleaned_seurat_list, "~/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/caruana-project/CART/Code/Seurat_analysis/2_annotation/cleaned_annotated_seurat.rds")
