# ============================================================
#  STEP 3b – Aggiornamento annotazioni e palette finale
#
#  Prerequisiti: STEP_1, STEP_2, STEP_3 già eseguiti.
#  Input:  annotated_list  (da STEP 3, in environment o su disco)
#          Me_I originale  (da STEP 1, in environment)
#          Bo_I originale  (da STEP 1, in environment)
#
#  Modifica 4 punti basati sull'ispezione visiva degli output:
#
#  P1 – Bo_I C0 vs C1: entrambi "CD4+ T cells" → distingui
#       "Naive CD4+ T cells" vs "Effector CD4+ T cells"
#       usando AverageExpression di CCR7/SELL vs CD44/GZMK
#
#  P2 – Palette vivace e consistente tra tutti e 3 i campioni
#       (stessa cellula = stesso colore in tutti gli UMAP)
#
#  P3 – DC rinominate "Dendritic Cells" (senza sottotipo)
#
#  P4 – Me_I: tentativo di isolare cluster DC con risoluzione
#       più alta (res 0.5 → 1.0). Ca_I esclusa: espressione
#       DC troppo diffusa per definire un cluster.
#
#  Output: annotated_list_v2 in environment
#          all_samples_annotated_v2.rds su disco
#          PNG UMAP in base_dir/Annotation_UMAP_v2/
# ============================================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

# ── UNICO PUNTO DA MODIFICARE ────────────────────────────────
base_dir <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/CART/Code/Seurat_analysis/3_data_cleaning/"
# ─────────────────────────────────────────────────────────────

out <- paste0(base_dir, "Annotation_UMAP_v2/")
dir.create(out, showWarnings = FALSE, recursive = TRUE)

section <- function(title) {
  cat(paste0("\n", strrep("=", 65), "\n  ", title, "\n", strrep("=", 65), "\n"))
}

# ── Caricamento ───────────────────────────────────────────────
if (exists("annotated_list")) {
  cat(">> Uso annotated_list da environment (STEP 3).\n")
} else {
  cat(">> Carico all_samples_annotated_final.rds da disco.\n")
  annotated_list <- readRDS(paste0(base_dir, "all_samples_annotated_final.rds"))
}

if (!exists("Bo_I") || !exists("Me_I")) {
  cat(">> Bo_I / Me_I non trovati: ricarico all_samples_clean_pre_annotation.rds\n")
  sl <- readRDS(paste0(base_dir, "all_samples_clean_pre_annotation.rds"))
  if (!exists("Bo_I")) Bo_I <- sl$Bo_I
  if (!exists("Me_I")) Me_I <- sl$Me_I
  rm(sl)
}

# ============================================================
# HELPER: applica annotazione a un oggetto Seurat
# ============================================================

apply_annotation <- function(obj, annotation_map, col_name = "cell_type") {
  Idents(obj) <- "seurat_clusters"
  cluster_ids  <- as.character(Idents(obj))
  missing      <- setdiff(unique(cluster_ids), names(annotation_map))
  if (length(missing) > 0)
    stop("Cluster senza annotazione: ", paste(missing, collapse = ", "))
  labels <- unname(annotation_map[cluster_ids])
  names(labels) <- colnames(obj)
  obj <- AddMetaData(obj, metadata = labels, col.name = col_name)
  Idents(obj) <- col_name
  return(obj)
}

# ============================================================
# P1 – Bo_I: Naive vs Effector CD4+
#
#  Geni naive:    CCR7, SELL, IL7R, TCF7, LEF1
#  Geni effector: CD44, GZMK, S100A4, LGALS1
#
#  Calcola score medio per C0 e C1; il cluster con score
#  naive > effector riceve "Naive CD4+ T cells", l'altro
#  "Effector CD4+ T cells". Fallback: cluster più grande = Naive.
# ============================================================
section("P1 | Bo_I – C0 vs C1: Naive vs Effector CD4+")

Idents(Bo_I) <- "seurat_clusters"

naive_genes    <- c("CCR7", "SELL", "IL7R", "TCF7", "LEF1")
effector_genes <- c("CD44", "GZMK", "S100A4", "LGALS1")
all_genes      <- unique(c(naive_genes, effector_genes))
all_genes      <- all_genes[all_genes %in% rownames(Bo_I)]

avg <- AverageExpression(
  Bo_I, features = all_genes,
  group.by = "seurat_clusters", assay = "RNA", slot = "data"
)$RNA

# Fix nomi colonna Seurat v5
colnames(avg) <- gsub("^g", "", colnames(avg))
colnames(avg) <- gsub("^RNA_snn_res\\.[0-9.]+_", "", colnames(avg))

cat("\nAverage expression C0 e C1:\n")
print(round(avg[, intersect(c("0","1"), colnames(avg)), drop = FALSE], 3))

get_mean <- function(mat, genes, cl) {
  g <- genes[genes %in% rownames(mat)]
  if (length(g) == 0 || !cl %in% colnames(mat)) return(0)
  mean(mat[g, cl])
}

s_naive_c0    <- get_mean(avg, naive_genes,    "0")
s_effector_c0 <- get_mean(avg, effector_genes, "0")
s_naive_c1    <- get_mean(avg, naive_genes,    "1")
s_effector_c1 <- get_mean(avg, effector_genes, "1")

cat(paste0(
  "\nC0 – Naive score: ",    round(s_naive_c0, 3),
  "  |  Effector score: ",   round(s_effector_c0, 3), "\n",
  "C1 – Naive score: ",      round(s_naive_c1, 3),
  "  |  Effector score: ",   round(s_effector_c1, 3), "\n"
))

label_c0 <- if (s_naive_c0 >= s_effector_c0) "Naive CD4+ T cells" else "Effector CD4+ T cells"
label_c1 <- if (s_naive_c1 >= s_effector_c1) "Naive CD4+ T cells" else "Effector CD4+ T cells"

# Fallback: se lo score assegna lo stesso label a entrambi,
# il cluster più grande è Naive (biologicamente più abbondanti)
if (label_c0 == label_c1) {
  n0 <- sum(Bo_I$seurat_clusters == "0")
  n1 <- sum(Bo_I$seurat_clusters == "1")
  label_c0 <- if (n0 >= n1) "Naive CD4+ T cells" else "Effector CD4+ T cells"
  label_c1 <- if (n0 >= n1) "Effector CD4+ T cells" else "Naive CD4+ T cells"
  cat("[FALLBACK] Score simili → cluster più grande = Naive\n")
}

cat(paste0("\n[DECISIONE P1]  C0 → ", label_c0, "  |  C1 → ", label_c1, "\n"))

# ============================================================
# P4 – Me_I: tenta di isolare cluster DC con res = 1.0
#
#  Me_I originale aveva 5 cluster a res=0.5 (dims=1:30).
#  I FeaturePlot mostrano VCAN/S100A12/LYZ/S100A8/S100A9
#  focalizzati nello stesso angolo dell'UMAP → possibile
#  piccolo cluster DC assorbito. Proviamo res=1.0.
#
#  Strategia:
#   1. FindClusters res=1.0 su Me_I (neighbors già calcolati)
#   2. Per ogni nuovo cluster, calcola DC module score
#   3. Se un cluster piccolo ha score >> media → DC cluster
#   4. Mappa i nuovi cluster a cell types
# ============================================================
section("P4 | Me_I – Reclustering a risoluzione maggiore per DC")

# Geni DC per il module score
dc_genes_all <- c(
  "AQP9", "VCAN", "S100A12", "LILRB2", "MS4A6A", "FPR1", "HCK",
  "HLA-DRA", "HLA-DPB1", "HLA-DQA1", "CD1C", "CLEC9A", "LILRA4",
  "FCER1A", "ITGAX", "XCR1", "BATF3", "LYZ", "S100A8", "S100A9"
)
dc_genes_ok <- dc_genes_all[dc_genes_all %in% rownames(Me_I)]

# Re-run FindNeighbors (stessi parametri di 1_process_data.R)
Me_I_hires <- FindNeighbors(Me_I, dims = 1:30, verbose = FALSE)
Me_I_hires <- FindClusters(Me_I_hires, resolution = 1.0, verbose = FALSE)

cat(paste0("\nMe_I – nuovi cluster a res=1.0: ",
           length(unique(Me_I_hires$seurat_clusters)), "\n"))
print(table(Me_I_hires$seurat_clusters))

# DC module score su tutti i nuovi cluster
Me_I_hires <- AddModuleScore(Me_I_hires,
                              features = list(dc_genes_ok),
                              name     = "DC_score")
Me_I_hires$DC_score <- Me_I_hires$DC_score1
Me_I_hires$DC_score1 <- NULL

dc_by_cluster <- Me_I_hires@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(
    n_cells     = n(),
    mean_DC     = round(mean(DC_score), 4),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_DC))

cat("\nDC score per nuovo cluster (ordinato):\n")
print(as.data.frame(dc_by_cluster))

# Soglia automatica: un cluster è DC se:
#  (a) ha DC score > media + 2*sd di tutti i cluster
#  (b) ha < 10% delle cellule totali (cluster piccolo)
global_mean_dc <- mean(dc_by_cluster$mean_DC)
global_sd_dc   <- sd(dc_by_cluster$mean_DC)
dc_threshold   <- global_mean_dc + 2 * global_sd_dc
total_cells    <- ncol(Me_I_hires)
max_pct_dc     <- 0.10

dc_candidates <- dc_by_cluster %>%
  filter(mean_DC > dc_threshold,
         n_cells / total_cells < max_pct_dc)

cat(paste0("\nSoglia DC score (mean+2sd): ", round(dc_threshold, 4), "\n"))

if (nrow(dc_candidates) > 0) {
  dc_clusters_me <- as.character(dc_candidates$seurat_clusters)
  cat(paste0("[TROVATO] Cluster DC in Me_I a res=1.0: ",
             paste(dc_clusters_me, collapse = ", "),
             " (", paste(dc_candidates$n_cells, collapse="+"), " cellule)\n"))
  use_hires_me <- TRUE
} else {
  cat("[NON TROVATO] Nessun cluster DC definibile in Me_I anche a res=1.0.\n")
  cat("  Il segnale DC è diffuso e non sufficiente per un cluster separato.\n")
  use_hires_me <- FALSE
}

# ============================================================
# COSTRUZIONE DIZIONARI FINALI
# ============================================================
section("Costruzione dizionari annotazione aggiornati")

# ── Bo_I ─────────────────────────────────────────────────────
# P1: C0/C1 distinti come Naive/Effector
# P3: C4 → "Dendritic Cells" (era "Monocyte-derived Dendritic Cells")
annotation_Bo_v2 <- c(
  "0" = label_c0,
  "1" = label_c1,
  "2" = "Cytotoxic CD8+ T cells",
  "3" = "Memory T cells",
  "4" = "Dendritic Cells"
)
cat("\n--- Bo_I v2 ---\n"); print(annotation_Bo_v2)

# ── Ca_I ─────────────────────────────────────────────────────
# Nessuna modifica ai cluster. C5 rimane "Tregs".
# P3: nessuna DC presente → nessuna modifica.
annotation_Ca_v2 <- c(
  "0" = "Proliferating T cells",
  "1" = "Naive CD4+ T cells",
  "2" = "Effector CD4+ T cells",
  "3" = "Cytotoxic CD8+ T cells",
  "4" = "Tregs",
  "5" = "Tregs"
)
cat("\n--- Ca_I v2 ---\n"); print(annotation_Ca_v2)

# ── Me_I ─────────────────────────────────────────────────────
# Se il reclustering ha trovato un cluster DC, usa i nuovi cluster.
# Altrimenti mantieni i 5 cluster originali.

if (use_hires_me) {
  # Identifica il cluster T proliferante a res=1.0 (più numerosi)
  # Mappa automatica: cluster NON-DC → cell type per similarità con
  # annotazione originale (usa il cell_type di maggioranza da STEP 3)
  Me_I_hires$cell_type_orig <- annotated_list$Me_I$cell_type[colnames(Me_I_hires)]

  cluster_map_me <- Me_I_hires@meta.data %>%
    group_by(seurat_clusters) %>%
    count(cell_type_orig) %>%
    slice_max(n, n = 1) %>%
    ungroup() %>%
    select(seurat_clusters, majority_type = cell_type_orig)

  cat("\nMapping cluster hires → cell type (maggioranza):\n")
  print(as.data.frame(cluster_map_me))

  annotation_Me_v2 <- setNames(
    cluster_map_me$majority_type,
    as.character(cluster_map_me$seurat_clusters)
  )
  # Sovrascrivi cluster DC
  annotation_Me_v2[dc_clusters_me] <- "Dendritic Cells"

  cat("\n--- Me_I v2 (hires, con DC) ---\n")
  print(annotation_Me_v2)

} else {
  # Mantieni annotazione originale con 5 cluster
  annotation_Me_v2 <- c(
    "0" = "Cytotoxic CD8+ T cells",
    "1" = "Proliferating T cells",
    "2" = "Tregs",
    "3" = "Proliferating T cells",
    "4" = "Proliferating CD8+ T cells"
  )
  cat("\n--- Me_I v2 (cluster originali, nessun DC) ---\n")
  print(annotation_Me_v2)
}

# ============================================================
# APPLICA LE ANNOTAZIONI
# ============================================================
section("Applicazione annotazioni v2")

Bo_I_ann_v2 <- apply_annotation(annotated_list$Bo_I, annotation_Bo_v2)
Ca_I_ann_v2 <- apply_annotation(annotated_list$Ca_I, annotation_Ca_v2)

if (use_hires_me) {
  Me_I_ann_v2 <- apply_annotation(Me_I_hires, annotation_Me_v2)
} else {
  Me_I_ann_v2 <- apply_annotation(annotated_list$Me_I, annotation_Me_v2)
}

for (nm in c("Bo_I_ann_v2", "Ca_I_ann_v2", "Me_I_ann_v2")) {
  obj <- get(nm)
  cat(paste0("\n[", nm, "]:\n"))
  print(table(obj$cell_type))
}

# ============================================================
# P2 – PALETTE VIVACE E CONSISTENTE
#
#  Gerarchia cromatica:
#   Blu     → CD4 (Naive = blu scuro, Effector = blu medio)
#   Rosso   → CD8 citotossici
#   Verde   → Memory T
#   Viola   → Proliferating (T generico = viola chiaro,
#              CD8 proliferanti = viola scuro)
#   Arancio → Tregs
#   Fucsia  → Dendritic Cells
# ============================================================
section("P2 | Palette vivace e consistente")

vivid_palette <- c(
  "Naive CD4+ T cells"         = "#1565C0",   # blu intenso
  "Effector CD4+ T cells"      = "#42A5F5",   # blu cielo
  "Cytotoxic CD8+ T cells"     = "#C62828",   # rosso intenso
  "Memory T cells"             = "#2E7D32",   # verde scuro
  "Proliferating T cells"      = "#AB47BC",   # viola medio
  "Proliferating CD8+ T cells" = "#4A148C",   # viola scurissimo
  "Tregs"                      = "#E65100",   # arancio bruciato
  "Dendritic Cells"            = "#E91E63"    # fucsia acceso
)

get_colors <- function(obj) {
  types <- sort(unique(as.character(obj$cell_type)))
  cols  <- vivid_palette[types]
  # Grigio per tipi non in palette (non dovrebbe succedere)
  cols[is.na(cols)] <- "#9E9E9E"
  return(cols)
}

# ============================================================
# PLOT UMAP FINALI
# ============================================================
section("UMAP finali con annotazioni v2")

plot_umap_final <- function(obj, sample_name, palette) {
  Idents(obj) <- "cell_type"
  cols <- get_colors(obj)

  p_label <- DimPlot(
    obj, reduction = "umap", label = TRUE,
    label.size = 3.5, repel = TRUE,
    cols = cols, pt.size = 0.6
  ) +
    ggtitle(paste0(sample_name, " – Con label")) +
    theme_classic(base_size = 12) +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 13),
      legend.text     = element_text(size = 9),
      legend.key.size = unit(0.45, "cm")
    ) +
    guides(color = guide_legend(override.aes = list(size = 3.5)))

  p_clean <- DimPlot(
    obj, reduction = "umap", label = FALSE,
    cols = cols, pt.size = 0.6
  ) +
    ggtitle(paste0(sample_name, " – Senza label")) +
    theme_classic(base_size = 12) +
    theme(
      plot.title      = element_text(hjust = 0.5, size = 12),
      legend.text     = element_text(size = 9),
      legend.key.size = unit(0.45, "cm")
    ) +
    guides(color = guide_legend(override.aes = list(size = 3.5)))

  combined <- p_label | p_clean
  path <- paste0(out, sample_name, "_UMAP_v2.png")
  ggsave(path, plot = combined, width = 16, height = 7, dpi = 300, bg = "white")
  cat(paste0("[", sample_name, "] UMAP v2 salvata → ", path, "\n"))
  return(p_label)
}

p_Bo <- plot_umap_final(Bo_I_ann_v2, "Bo_I", vivid_palette)
p_Ca <- plot_umap_final(Ca_I_ann_v2, "Ca_I", vivid_palette)
p_Me <- plot_umap_final(Me_I_ann_v2, "Me_I", vivid_palette)

# Pannello combinato verticale
panel_all <- p_Bo / p_Ca / p_Me
path_all  <- paste0(out, "ALL_samples_UMAP_v2.png")
ggsave(path_all, plot = panel_all, width = 10, height = 21,
       dpi = 300, bg = "white")
cat(paste0("[ALL] Pannello combinato → ", path_all, "\n"))

# ============================================================
# SALVATAGGIO
# ============================================================
section("Salvataggio")

annotated_list_v2 <- list(
  Bo_I = Bo_I_ann_v2,
  Ca_I = Ca_I_ann_v2,
  Me_I = Me_I_ann_v2
)

rds_out <- paste0(base_dir, "all_samples_annotated_v2.rds")
saveRDS(annotated_list_v2, rds_out)

cat(paste0(
  "\n", strrep("=", 65), "\n",
  "  STEP 3b COMPLETATO\n\n",
  "  Modifiche applicate:\n",
  "  P1 – Bo_I C0 → ", label_c0, "\n",
  "       Bo_I C1 → ", label_c1, "\n",
  "  P2 – Palette vivace e consistente tra i 3 campioni\n",
  "  P3 – DC rinominate 'Dendritic Cells'\n",
  "  P4 – Me_I reclustering res=1.0: DC ",
  ifelse(use_hires_me, "TROVATE", "non trovate"), "\n\n",
  "  Oggetto: annotated_list_v2  ($Bo_I | $Ca_I | $Me_I)\n",
  "  RDS:     ", rds_out, "\n",
  "  UMAP:    ", out, "\n",
  strrep("=", 65), "\n"
))
