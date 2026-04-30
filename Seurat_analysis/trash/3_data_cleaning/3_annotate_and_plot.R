# ============================================================
#  STEP 3 – Annotazione manuale e UMAP
#
#  Prerequisito: STEP_2_resolve_doubts.R già eseguito
#  Input:  Bo_I, Ca_I, Me_I + final_annotations in environment
#  Output: annotated_list in environment
#          all_samples_annotated_final.rds su disco
#          PNG UMAP in base_dir/Annotation_UMAP/
# ============================================================

library(Seurat)
library(ggplot2)
library(patchwork)

# ── UNICO PUNTO DA MODIFICARE ────────────────────────────────
base_dir <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/CART/Code/Seurat_analysis/3_data_cleaning/"
# ─────────────────────────────────────────────────────────────

out <- paste0(base_dir, "Annotation_UMAP/")
dir.create(out, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 1. DIZIONARI DI ANNOTAZIONE
#    Usa final_annotations da STEP 2 se disponibile,
#    altrimenti usa i dizionari di fallback.
# ============================================================

if (exists("final_annotations")) {
  cat(">> Uso 'final_annotations' da STEP_2_resolve_doubts.R\n\n")
  annotation_Bo <- final_annotations$Bo_I
  annotation_Ca <- final_annotations$Ca_I
  annotation_Me <- final_annotations$Me_I
} else {
  cat(">> ATTENZIONE: 'final_annotations' non trovato.\n")
  cat("   Uso dizionari di fallback. Per risultati ottimali\n")
  cat("   esegui prima STEP_2_resolve_doubts.R\n\n")

  annotation_Bo <- c(
    "0" = "CD4+ T cells",
    "1" = "CD8+ T cells",
    "2" = "Cytotoxic CD8+ T cells",
    "3" = "Memory T cells",
    "4" = "Dendritic Cells"
  )
  annotation_Ca <- c(
    "0" = "Proliferating T cells",
    "1" = "Naive CD4+ T cells",
    "2" = "Effector CD4+ T cells",
    "3" = "Cytotoxic CD8+ T cells",
    "4" = "Tregs",
    "5" = "NK cells"
  )
  annotation_Me <- c(
    "0" = "Cytotoxic CD8+ T cells",
    "1" = "Proliferating T cells",
    "2" = "Tregs",
    "3" = "Proliferating T cells",
    "4" = "Proliferating CD8+ T cells"
  )
}

# ============================================================
# 2. PALETTE COLORI CONDIVISA
#    Copre tutti i possibili label (inclusi quelli con fase
#    del ciclo cellulare e sottotipi DC da STEP 2).
# ============================================================

base_palette <- c(
  "CD4+ T cells"                     = "#4E9AF1",
  "Naive CD4+ T cells"               = "#A8D1F5",
  "Effector CD4+ T cells"            = "#1565C0",
  "CD8+ T cells"                     = "#E87D3E",
  "Cytotoxic CD8+ T cells"           = "#C0392B",
  "Memory T cells"                   = "#27AE60",
  "Proliferating T cells"            = "#8E44AD",
  "Proliferating T cells (S)"        = "#9B59B6",
  "Proliferating T cells (G2M)"      = "#6C3483",
  "Proliferating CD8+ T cells"       = "#4A235A",
  "Tregs"                            = "#F39C12",
  "NK cells"                         = "#16A085",
  "Dendritic Cells"                  = "#E74C3C",
  "cDC1 (Dendritic Cells)"           = "#E74C3C",
  "cDC2 (Dendritic Cells)"           = "#EC407A",
  "Monocyte-derived Dendritic Cells" = "#D81B60",
  "Plasmacytoid Dendritic Cells"     = "#AD1457"
)

get_colors <- function(obj) {
  types <- sort(unique(obj$cell_type))
  cols  <- base_palette[types]
  cols[is.na(cols)] <- "#999999"   # grigio per tipi non in palette
  return(cols)
}

# ============================================================
# 3. FUNZIONE DI ANNOTAZIONE
# ============================================================

annotate_sample <- function(obj, annotation_map, sample_name) {
  Idents(obj) <- "seurat_clusters"
  cluster_ids <- as.character(Idents(obj))

  # Controlla che tutti i cluster abbiano un'annotazione
  missing <- setdiff(unique(cluster_ids), names(annotation_map))
  if (length(missing) > 0)
    stop(paste0("[ERRORE] ", sample_name, ": cluster senza annotazione: ",
                paste(missing, collapse = ", ")))

  # Costruisce il vettore nominato per barcode.
  # unname() rimuove i cluster ID dai nomi del vettore intermedio;
  # riassegniamo i barcode come nomi cosi Seurat trova il match
  # corretto e non solleva "No cell overlap".
  cell_labels <- unname(annotation_map[cluster_ids])
  names(cell_labels) <- colnames(obj)

  obj <- AddMetaData(obj, metadata = cell_labels, col.name = "cell_type")
  Idents(obj) <- "cell_type"

  cat(paste0("\n[", sample_name, "] Annotazione completata:\n"))
  print(table(obj$cell_type))
  return(obj)
}

# ============================================================
# 4. APPLICA L'ANNOTAZIONE
# ============================================================

Bo_I_ann <- annotate_sample(Bo_I, annotation_Bo, "Bo_I")
Ca_I_ann <- annotate_sample(Ca_I, annotation_Ca, "Ca_I")
Me_I_ann <- annotate_sample(Me_I, annotation_Me, "Me_I")

# ============================================================
# 5. UMAP ANNOTATE – con label e senza (per pubblicazione)
# ============================================================

plot_umap <- function(obj, sample_name) {
  cols <- get_colors(obj)
  Idents(obj) <- "cell_type"

  p_label <- DimPlot(obj, reduction = "umap", label = TRUE,
                     label.size = 3.5, repel = TRUE,
                     cols = cols, pt.size = 0.5) +
    ggtitle(paste0(sample_name, " – Con label")) +
    theme_classic(base_size = 12) +
    theme(plot.title      = element_text(hjust = 0.5, face = "bold", size = 13),
          legend.text     = element_text(size = 9),
          legend.key.size = unit(0.4, "cm")) +
    guides(color = guide_legend(override.aes = list(size = 3)))

  p_clean <- DimPlot(obj, reduction = "umap", label = FALSE,
                     cols = cols, pt.size = 0.5) +
    ggtitle(paste0(sample_name, " – Senza label")) +
    theme_classic(base_size = 12) +
    theme(plot.title      = element_text(hjust = 0.5, size = 12),
          legend.text     = element_text(size = 9),
          legend.key.size = unit(0.4, "cm")) +
    guides(color = guide_legend(override.aes = list(size = 3)))

  combined <- p_label | p_clean
  path <- paste0(out, sample_name, "_UMAP_annotated.png")
  ggsave(path, plot = combined, width = 16, height = 7, dpi = 300, bg = "white")
  cat(paste0("[", sample_name, "] UMAP salvata → ", path, "\n"))

  return(p_label)   # restituisce solo la versione con label per il pannello
}

umap_Bo <- plot_umap(Bo_I_ann, "Bo_I")
umap_Ca <- plot_umap(Ca_I_ann, "Ca_I")
umap_Me <- plot_umap(Me_I_ann, "Me_I")

# Pannello combinato (tutti e tre i campioni in verticale)
combined_all <- umap_Bo / umap_Ca / umap_Me
path_all <- paste0(out, "ALL_samples_UMAP_annotated.png")
ggsave(path_all, plot = combined_all, width = 10, height = 20,
       dpi = 300, bg = "white")
cat(paste0("[ALL] Pannello combinato → ", path_all, "\n"))

# ============================================================
# 6. SALVATAGGIO SU DISCO
# ============================================================

annotated_list <- list(Bo_I = Bo_I_ann, Ca_I = Ca_I_ann, Me_I = Me_I_ann)

rds_out <- paste0(base_dir, "all_samples_annotated_final.rds")
saveRDS(annotated_list, rds_out)
cat(paste0("\nRDS salvato → ", rds_out, "\n"))

cat(paste0(
  "\n", strrep("=", 55), "\n",
  "  STEP 3 COMPLETATO\n",
  "  'annotated_list' disponibile in environment.\n",
  "  Accesso: annotated_list$Bo_I | $Ca_I | $Me_I\n",
  "  Metadata cell type: obj$cell_type\n",
  "  RDS: ", rds_out, "\n",
  "  UMAP: ", out, "\n",
  "  Prossimo step: esegui STEP_4_DC_analysis.R\n",
  strrep("=", 55), "\n"
))
