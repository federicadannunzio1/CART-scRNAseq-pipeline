# ====== Parametri ======
base_dir <- "D:/Progetti/Ignazio_2/Riallineamento/cellranger_output/gex_1"
out_dir  <- "D:/Progetti/Ignazio_2/Data/Post_reallineamento"

car_gene <- "CAR"          # nome del transgene
car_umi_threshold <- 1L    # soglia UMI per CAR+

# ====== Librerie ======
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

# ====== Scansione campioni ======
samples <- list.dirs(base_dir, full.names = FALSE, recursive = FALSE)
samples <- samples[nzchar(samples)]
has_matrix <- function(s) dir.exists(file.path(base_dir, s, "count", "sample_filtered_feature_bc_matrix"))
samples <- samples[vapply(samples, has_matrix, logical(1))]
if (length(samples) == 0) {
  stop("Nessun campione trovato in ", base_dir,
       " con subpath 'count/sample_filtered_feature_bc_matrix'.")
}

# ====== Output containers ======
seurat_list <- list()
labels_list <- list()  # <- LISTA di data frame, uno per esperimento

message("Campioni trovati: ", paste(samples, collapse = ", "))

for (s in samples) {
  message("==> Elaboro: ", s)
  mtx_dir <- file.path(base_dir, s, "count", "sample_filtered_feature_bc_matrix")
  
  # 1) Leggi counts
  counts <- Read10X(data.dir = mtx_dir)
  if (is.list(counts)) {
    counts <- if ("Gene Expression" %in% names(counts)) counts[["Gene Expression"]] else counts[[1]]
  }
  
  # 2) Controllo gene CAR
  if (!(car_gene %in% rownames(counts))) {
    stop("Il gene '", car_gene, "' non è presente nel campione '", s, "'.")
  }
  
  # 3) Crea Seurat object
  seu <- CreateSeuratObject(counts = counts, project = s, min.cells = 0, min.features = 0)
  
  # 4) Rinomina cellule come SAMPLE_BARCODE
  colnames(seu) <- paste0(s, "_", colnames(seu))
  
  # 5) Calcola UMI CAR e metadato
  car_umi <- Matrix::colSums(GetAssayData(seu, slot = "counts")[car_gene, , drop = FALSE])
  seu$CAR_UMI <- car_umi
  seu$IS_CAR  <- ifelse(car_umi >= car_umi_threshold, "YES", "NO")
  
  # 6) Popola la LISTA di data frame (uno per sample)
  labels_list[[s]] <- data.frame(
    cell  = colnames(seu),
    IS_CAR = seu$IS_CAR,
    stringsAsFactors = FALSE
  )
  
  # 7) Aggiungi alla lista Seurat
  seurat_list[[s]] <- seu
  
  # 8) Pulizia RAM
  rm(counts, seu); gc()
}

# ====== Salvataggi ======
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# a) Lista di Seurat
saveRDS(seurat_list, file = file.path(out_dir, "seurat_list_with_IS_CAR.rds"))

# b) Lista di data frame (uno per esperimento)
saveRDS(labels_list, file = file.path(out_dir, "IS_CAR_labels_list.rds"))

# c) (comodo) CSV per singolo esperimento
csv_dir <- file.path(out_dir, "IS_CAR_labels_per_sample_csv")
if (!dir.exists(csv_dir)) dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)
for (s in names(labels_list)) {
  write.csv(labels_list[[s]],
            file = file.path(csv_dir, paste0("IS_CAR_labels_", s, ".csv")),
            row.names = FALSE)
}

message("Fatto. Salvati:")
message(" - Lista Seurat: ", file.path(out_dir, "seurat_list_with_IS_CAR.rds"))
message(" - Lista data frame: ", file.path(out_dir, "IS_CAR_labels_list.rds"))
message(" - CSV per sample in: ", csv_dir)
