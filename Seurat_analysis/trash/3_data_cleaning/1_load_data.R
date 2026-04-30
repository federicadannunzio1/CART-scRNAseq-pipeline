# ============================================================
#  STEP 1 – Caricamento e preparazione oggetti Seurat
#  Input:  all_samples_clean_pre_annotation.rds  (da 1_process_data.R)
#  Output: Bo_I, Ca_I, Me_I in environment (layers unificati)
# ============================================================

library(Seurat)

# ── UNICO PUNTO DA MODIFICARE ────────────────────────────────
base_dir <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/CART/Code/Seurat_analysis/3_data_cleaning/"
# ─────────────────────────────────────────────────────────────

rds_in <- paste0(base_dir, "all_samples_clean_pre_annotation.rds")

cat("Caricamento:", rds_in, "\n")
seurat_list <- readRDS(rds_in)

# ── Estrazione oggetti ────────────────────────────────────────
Bo_I <- seurat_list$Bo_I
Ca_I <- seurat_list$Ca_I
Me_I <- seurat_list$Me_I

# ── Verifica seurat_clusters (creata da FindClusters in 1_process_data.R) ──
check_clusters <- function(obj, name) {
  if (!"seurat_clusters" %in% colnames(obj@meta.data))
    stop(paste0("[ERRORE] ", name, ": 'seurat_clusters' non trovata. ",
                "Verifica che 1_process_data.R sia completato."))
  Idents(obj) <- "seurat_clusters"
  tbl <- table(obj$seurat_clusters)
  cat(paste0("[OK] ", name, " | ", ncol(obj), " cellule | ",
             length(tbl), " cluster: ", paste(names(tbl), collapse = ", "), "\n"))
  return(obj)
}

Bo_I <- check_clusters(Bo_I, "Bo_I")
Ca_I <- check_clusters(Ca_I, "Ca_I")
Me_I <- check_clusters(Me_I, "Me_I")

# ── JoinLayers (Seurat v5) ────────────────────────────────────
# 1_process_data.R produce layer "counts.S429_I" etc. (splittati per campione).
# JoinLayers li unifica in un unico layer "counts", necessario per
# AddModuleScore e per garantire compatibilità con tutte le funzioni Seurat.
join_layers <- function(obj, name) {
  split_counts <- grep("^counts\\.", Layers(obj), value = TRUE)
  if (length(split_counts) > 0) {
    cat(paste0("[", name, "] JoinLayers(): ", paste(split_counts, collapse = ", "),
               " → counts\n"))
    obj <- JoinLayers(obj)
  } else {
    cat(paste0("[", name, "] Layer già unificati.\n"))
  }
  return(obj)
}

Bo_I <- join_layers(Bo_I, "Bo_I")
Ca_I <- join_layers(Ca_I, "Ca_I")
Me_I <- join_layers(Me_I, "Me_I")

cat(paste0(
  "\n", strrep("=", 55), "\n",
  "  STEP 1 COMPLETATO\n",
  "  Oggetti in environment: Bo_I | Ca_I | Me_I\n",
  "  Prossimo step: esegui STEP_2_resolve_doubts.R\n",
  strrep("=", 55), "\n"
))

