library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

# ==========================================
# 1. CARICAMENTO DATI
# ==========================================
base_dir <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/CART/Code/Seurat_analysis/3_data_cleaning/"
list_path <- paste0(base_dir, "all_samples_clean_pre_annotation.rds")

cat("Caricamento lista oggetti Seurat pre-annotazione...\n")
seurat_list <- readRDS(list_path)

Bo_I <- seurat_list$Bo_I
Ca_I <- seurat_list$Ca_I
Me_I <- seurat_list$Me_I

# ==========================================
# 2. ANNOTAZIONE ARMONIZZATA E SPECIFICA
# ==========================================
cat("Applicazione della nomenclatura coerente ai campioni...\n")

# ----- Bo_I -----
Idents(Bo_I) <- "seurat_clusters"
labels_bo <- c(
  "0" = "CD4_T_Cells",            # Corretto: Non stanno proliferando
  "1" = "CD4_T_Cells",            # Corretto: Non stanno proliferando
  "2" = "CD8_Cytotoxic_Effector", # Corretto in base al nuovo clustering
  "3" = "CD4_T_Cells",            # Corretto: Non stanno proliferando
  "4" = "Dendritic_Cells"         # Corretto in base al nuovo clustering (HLA-DRA+, LYZ+)
)
Bo_I <- RenameIdents(Bo_I, labels_bo)
Bo_I$annotazione_armonizzata <- Idents(Bo_I)

# ----- Ca_I -----
Idents(Ca_I) <- "seurat_clusters"
labels_ca <- c(
  "0" = "Proliferating_CD8_T",    # Corretto e specifico: CD4 negativo, NKG7+ e GZMA+
  "1" = "CD4_Naive",
  "2" = "CD4_Memory",
  "3" = "CD8_Cytotoxic_Effector", 
  "4" = "CD4_Effector",
  "5" = "CD4_Regulatory"          
)
Ca_I <- RenameIdents(Ca_I, labels_ca)
Ca_I$annotazione_armonizzata <- Idents(Ca_I)

# ----- Me_I -----
Idents(Me_I) <- "seurat_clusters"
labels_me <- c(
  "0" = "CD8_Cytotoxic_Effector",
  "1" = "Proliferating_CD4_T",    
  "2" = "CD4_Naive_Memory_Regulatory", # Etichetta descrittiva precisa senza subclustering
  "3" = "Proliferating_CD4_T",    
  "4" = "Proliferating_CD8_T"     
)
Me_I <- RenameIdents(Me_I, labels_me)
Me_I$annotazione_armonizzata <- Idents(Me_I)

# ==========================================
# 3. VERIFICA UMAP E SALVATAGGIO
# ==========================================
cat("Generazione degli UMAP annotati e salvataggio...\n")

# Disegno gli UMAP
p_bo_umap <- DimPlot(Bo_I, reduction = "umap", label = TRUE, repel = TRUE, label.size = 3.5) + ggtitle("Bo_I") + NoLegend()
p_ca_umap <- DimPlot(Ca_I, reduction = "umap", label = TRUE, repel = TRUE, label.size = 3.5) + ggtitle("Ca_I") + NoLegend()
p_me_umap <- DimPlot(Me_I, reduction = "umap", label = TRUE, repel = TRUE, label.size = 3.5) + ggtitle("Me_I") + NoLegend()

print(p_bo_umap | p_ca_umap | p_me_umap)

# Salvataggio nella lista finale
seurat_list_annotated <- list(Bo_I = Bo_I, Ca_I = Ca_I, Me_I = Me_I)
save_path <- paste0(base_dir, "all_samples_annotated_final.rds")
saveRDS(seurat_list_annotated, save_path)
cat(paste0("Annotazione completata in modo accurato. Dati salvati in: ", save_path, "\n"))

# ==========================================
# 4. CONTROLLO DI SICUREZZA 
# ==========================================
# Stampa un DotPlot per confermare le etichette assegnate 
# e mostrare l'espressione mista nel cluster 2 di Me_I
check_markers <- c("CD4", "CD8A", "MKI67", "TOP2A", "CCR7", "IL7R", "FOXP3", "GZMB", "HLA-DRA")
p_check <- DotPlot(Bo_I, features = check_markers, group.by = "annotazione_armonizzata") + 
  RotatedAxis() + 
  ggtitle("Verifica Marcatori Bo")
print(p_check)

p_check1 <- DotPlot(Me_I, features = check_markers, group.by = "annotazione_armonizzata") + 
  RotatedAxis() + 
  ggtitle("Verifica Marcatori  Me")
print(p_check1)

p_check2 <- DotPlot(Ca_I, features = check_markers, group.by = "annotazione_armonizzata") + 
  RotatedAxis() + 
  ggtitle("Verifica Marcatori Ca")
print(p_check2)



