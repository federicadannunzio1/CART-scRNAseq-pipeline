# ============================================================
#  STEP 4 – Analisi marcatori Dendritic Cells
#
#  Prerequisito: STEP_3_annotate_and_plot.R già eseguito
#  Input:  all_samples_annotated_final.rds  (da STEP 3)
#          oppure annotated_list in environment
#  Output: PNG in base_dir/DC_analysis/
#          Stima % cellule DC-positive per campione (console)
#
#  Logica: i marker DC derivano dal cluster confermato di Bo_I
#  (C4: HLA-DRA+, LYZ+, S100A8+, VCAN+, AQP9+) + letteratura.
#  Verifichiamo la presenza attesa di DC anche in Ca_I e Me_I.
# ============================================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(scales)

# ── UNICO PUNTO DA MODIFICARE ────────────────────────────────
base_dir <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/CART/Code/Seurat_analysis/3_data_cleaning/"
# ─────────────────────────────────────────────────────────────

out <- paste0(base_dir, "DC_analysis/")
dir.create(out, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 1. CARICAMENTO
#    Usa annotated_list in environment se disponibile,
#    altrimenti carica da disco.
# ============================================================

if (exists("annotated_list")) {
  cat(">> Uso 'annotated_list' da environment (STEP 3).\n\n")
  obj_Bo <- annotated_list$Bo_I
  obj_Ca <- annotated_list$Ca_I
  obj_Me <- annotated_list$Me_I
} else {
  rds_in <- paste0(base_dir, "all_samples_annotated_final.rds")
  cat(">> Carico da disco:", rds_in, "\n\n")
  ann_list <- readRDS(rds_in)
  obj_Bo <- ann_list$Bo_I
  obj_Ca <- ann_list$Ca_I
  obj_Me <- ann_list$Me_I
}

# Imposta identità
for (nm in c("obj_Bo", "obj_Ca", "obj_Me")) {
  obj <- get(nm)
  Idents(obj) <- if ("cell_type" %in% colnames(obj@meta.data)) "cell_type" else "seurat_clusters"
  assign(nm, obj)
}

# ============================================================
# 2. MARCATORI DC
#    Combinazione di: marker emersi da Bo_I C4 + pannello
#    canonico della letteratura (HLA-DRA, CD1C, CLEC9A, etc.)
# ============================================================

dc_markers_BoI <- c(
  "AQP9",    # moDC / neutrofili infiammatori
  "VCAN",    # DC immaturi / mieloidi
  "S100A12", # DC infiammatori / mieloidi
  "LILRB2",  # recettore inibitore mieloide
  "MS4A6A",  # marker mieloide
  "FPR1",    # recettore formil-peptidi
  "HCK"      # chinasi mieloide
)

dc_markers_canonical <- c(
  "HLA-DRA", "HLA-DPB1", "HLA-DQA1",  # MHC II (tutti i DC)
  "CD1C", "FCER1A", "ITGAX",           # cDC2
  "CLEC9A", "XCR1", "BATF3",           # cDC1
  "LILRA4",                             # pDC
  "LYZ", "S100A8", "S100A9"            # mieloidi generici
)

all_dc <- unique(c(dc_markers_BoI, dc_markers_canonical))
cat("Marcatori DC usati:", paste(all_dc, collapse = ", "), "\n\n")

# ============================================================
# 3. DC MODULE SCORE (AddModuleScore)
# ============================================================

add_dc_score <- function(obj, sample_name) {
  present <- all_dc[all_dc %in% rownames(obj)]
  absent  <- setdiff(all_dc, rownames(obj))
  if (length(absent) > 0)
    cat(paste0("[", sample_name, "] Geni assenti (saltati): ",
               paste(absent, collapse = ", "), "\n"))

  obj <- AddModuleScore(obj, features = list(present), name = "DC_score")
  obj$DC_score  <- obj$DC_score1
  obj$DC_score1 <- NULL
  return(obj)
}

obj_Bo <- add_dc_score(obj_Bo, "Bo_I")
obj_Ca <- add_dc_score(obj_Ca, "Ca_I")
obj_Me <- add_dc_score(obj_Me, "Me_I")

# ============================================================
# 4. FEATUREPLOTS – espressione geni DC sui UMAP
# ============================================================

plot_dc_features <- function(obj, sample_name) {
  features_ok <- all_dc[all_dc %in% rownames(obj)]
  if (length(features_ok) == 0) return(invisible(NULL))

  p <- FeaturePlot(obj, features = features_ok, ncol = 4, pt.size = 0.3,
                   min.cutoff = "q05", max.cutoff = "q95",
                   cols = c("lightgrey", "#C0392B")) &
    theme_classic(base_size = 9) &
    theme(plot.title = element_text(size = 9, face = "bold"))

  path <- paste0(out, sample_name, "_DC_FeaturePlot.png")
  ggsave(path, plot = p,
         width  = 16,
         height = ceiling(length(features_ok) / 4) * 4,
         dpi = 300, bg = "white")
  cat(paste0("[", sample_name, "] FeaturePlot DC → ", path, "\n"))
  return(invisible(p))
}

plot_dc_features(obj_Bo, "Bo_I")
plot_dc_features(obj_Ca, "Ca_I")
plot_dc_features(obj_Me, "Me_I")

# ============================================================
# 5. DC SCORE – UMAP + VlnPlot per cell type
# ============================================================

plot_dc_score <- function(obj, sample_name) {
  p_feat <- FeaturePlot(obj, features = "DC_score", pt.size = 0.5,
                         min.cutoff = "q05", max.cutoff = "q95",
                         cols = c("lightgrey", "#C0392B")) +
    ggtitle(paste0(sample_name, " – DC module score")) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  p_vln <- VlnPlot(obj, features = "DC_score", pt.size = 0.1,
                    cols = hue_pal()(length(unique(Idents(obj))))) +
    ggtitle(paste0(sample_name, " – DC score per cluster")) +
    xlab("") + ylab("DC module score") +
    theme_classic(base_size = 10) +
    theme(plot.title      = element_text(hjust = 0.5, face = "bold"),
          axis.text.x     = element_text(angle = 35, hjust = 1, size = 9),
          legend.position = "none")

  path <- paste0(out, sample_name, "_DC_score.png")
  ggsave(path, plot = p_feat | p_vln, width = 14, height = 6,
         dpi = 300, bg = "white")
  cat(paste0("[", sample_name, "] DC score plot → ", path, "\n"))
  return(invisible(NULL))
}

plot_dc_score(obj_Bo, "Bo_I")
plot_dc_score(obj_Ca, "Ca_I")
plot_dc_score(obj_Me, "Me_I")

# ============================================================
# 6. DOTPLOT – marcatori DC per cell type
# ============================================================

plot_dc_dotplot <- function(obj, sample_name) {
  features_ok <- all_dc[all_dc %in% rownames(obj)]
  if (length(features_ok) < 3) return(invisible(NULL))

  p <- DotPlot(obj, features = features_ok) +
    RotatedAxis() +
    ggtitle(paste0(sample_name, " – DC markers per cell type")) +
    theme_classic(base_size = 10) +
    theme(plot.title  = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_text(size = 8)) +
    scale_color_gradient2(low = "white", mid = "#9B59B6",
                          high = "#C0392B", midpoint = 0)

  path <- paste0(out, sample_name, "_DC_DotPlot.png")
  ggsave(path, plot = p,
         width  = 12,
         height = max(4, length(unique(Idents(obj))) * 0.6 + 2),
         dpi = 300, bg = "white")
  cat(paste0("[", sample_name, "] DotPlot DC → ", path, "\n"))
  return(invisible(p))
}

plot_dc_dotplot(obj_Bo, "Bo_I")
plot_dc_dotplot(obj_Ca, "Ca_I")
plot_dc_dotplot(obj_Me, "Me_I")

# ============================================================
# 7. STIMA QUANTITATIVA % CELLULE DC-POSITIVE
#    Soglia conservativa: DC_score > mean + 2*sd
# ============================================================

cat(paste0("\n", strrep("=", 55), "\n",
           "  STIMA PRESENZA DC (score > mean + 2sd)\n",
           strrep("=", 55), "\n"))

for (nm in c("obj_Bo", "obj_Ca", "obj_Me")) {
  obj   <- get(nm)
  sname <- sub("obj_", "", nm)
  scores    <- obj$DC_score
  threshold <- mean(scores) + 2 * sd(scores)
  n_high    <- sum(scores > threshold)
  pct       <- round(100 * n_high / ncol(obj), 2)

  cat(paste0("\n[", sname, "]\n",
             "  Threshold (mean+2sd): ", round(threshold, 3), "\n",
             "  Cellule DC-high:      ", n_high, " / ", ncol(obj),
             " (", pct, "%)\n"))

  group_var <- if ("cell_type" %in% colnames(obj@meta.data)) "cell_type" else "seurat_clusters"

  summary_group <- data.frame(group    = obj@meta.data[[group_var]],
                               dc_score = obj$DC_score) %>%
    group_by(group) %>%
    summarise(mean_score = round(mean(dc_score), 3),
              n_high     = sum(dc_score > threshold),
              pct_high   = round(100 * sum(dc_score > threshold) / n(), 2),
              .groups    = "drop") %>%
    arrange(desc(mean_score))

  cat("  Distribuzione per cell type:\n")
  print(as.data.frame(summary_group))
}

cat(paste0(
  "\n", strrep("=", 55), "\n",
  "  STEP 4 COMPLETATO\n",
  "  Output in: ", out, "\n",
  strrep("=", 55), "\n"
))
