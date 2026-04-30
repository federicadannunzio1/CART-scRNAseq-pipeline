# ============================================================
#  STEP 5 – Classificazione cellule DC-like in Ca_I e Me_I
#
#  Prerequisito: STEP_3b già eseguito (annotated_list_v2)
#
#  Motivazione:
#    I FeaturePlot di STEP 4 mostrano che in Ca_I e Me_I i
#    marker DC sono presenti ma diffusi su tutto l'UMAP, senza
#    un cluster separato. Questo pattern indica la presenza di
#    rare cellule mieloidi/DC disperse tra i linfociti T, non
#    sufficienti numericamente per formare un cluster distinto.
#
#    Approccio alternativo: classificazione a livello di
#    singola cellula tramite DC module score. Le cellule con
#    score > soglia stringente vengono etichettate "DC-like"
#    e visualizzate sull'UMAP esistente con cells.highlight.
#    Questo approccio è riportato in letteratura per identificare
#    "contaminating myeloid cells" in dataset di TIL.
#
#  Output per Ca_I e Me_I:
#    A) UMAP highlight: DC-like cells evidenziate in fucsia
#       su sfondo grigio (base annotation preservata)
#    B) DotPlot: marker DC per cell type annotato
#    C) Tabella % DC-like per campione e per cell type
#    D) UMAP con DC_score continuo (gradiente rosso)
#    E) VlnPlot DC_score per cell type
#
#  Output in: base_dir/DC_classification/
# ============================================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

# ── UNICO PUNTO DA MODIFICARE ────────────────────────────────
base_dir <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/CART/Code/Seurat_analysis/3_data_cleaning/"
# ─────────────────────────────────────────────────────────────

out <- paste0(base_dir, "DC_classification/")
dir.create(out, showWarnings = FALSE, recursive = TRUE)

section <- function(title) {
  cat(paste0("\n", strrep("=", 65), "\n  ", title, "\n", strrep("=", 65), "\n"))
}

# ── Caricamento ───────────────────────────────────────────────
if (exists("annotated_list_v2")) {
  cat(">> Uso annotated_list_v2 da environment.\n")
} else {
  path_rds <- paste0(base_dir, "all_samples_annotated_v2.rds")
  cat(">> Carico:", path_rds, "\n")
  annotated_list_v2 <- readRDS(path_rds)
}

obj_Ca <- annotated_list_v2$Ca_I
obj_Me <- annotated_list_v2$Me_I
# Bo_I ha già il cluster DC definito → non serve qui

Idents(obj_Ca) <- "cell_type"
Idents(obj_Me) <- "cell_type"

# ============================================================
# 1. MARCATORI DC
#    Pannello derivato da Bo_I C4 (cluster DC confermato)
#    + marcatori canonici della letteratura
# ============================================================

dc_genes <- c(
  # Marker emersi da Bo_I C4 (moDC-like)
  "AQP9", "VCAN", "S100A12", "LILRB2", "MS4A6A", "FPR1", "HCK",
  # MHC II – tutti i DC
  "HLA-DRA", "HLA-DPB1", "HLA-DQA1",
  # cDC2
  "CD1C", "FCER1A", "ITGAX",
  # cDC1
  "CLEC9A", "XCR1", "BATF3",
  # pDC
  "LILRA4",
  # Mieloidi generici
  "LYZ", "S100A8", "S100A9"
)

# ============================================================
# 2. DC MODULE SCORE
#    AddModuleScore calcola uno score per ogni cellula basato
#    sull'espressione dei geni DC rispetto a geni di controllo
#    random (stesso numero, matchati per espressione media).
# ============================================================

add_dc_score <- function(obj, sample_name) {
  genes_ok <- dc_genes[dc_genes %in% rownames(obj)]
  genes_absent <- setdiff(dc_genes, rownames(obj))
  if (length(genes_absent) > 0)
    cat(paste0("[", sample_name, "] Geni assenti: ",
               paste(genes_absent, collapse=", "), "\n"))
  cat(paste0("[", sample_name, "] Geni usati per DC score: ",
             length(genes_ok), "/", length(dc_genes), "\n"))

  obj <- AddModuleScore(obj, features = list(genes_ok), name = "DC_score")
  obj$DC_score  <- obj$DC_score1
  obj$DC_score1 <- NULL
  return(obj)
}

obj_Ca <- add_dc_score(obj_Ca, "Ca_I")
obj_Me <- add_dc_score(obj_Me, "Me_I")

# ============================================================
# 3. DEFINIZIONE SOGLIA E CLASSIFICAZIONE DC-LIKE
#
#  Soglia: percentile 97.5 del DC_score (top 2.5% delle cellule).
#  Più interpretabile rispetto a mean+2sd perché:
#   - non dipende dalla distribuzione (spesso asimmetrica)
#   - seleziona sempre una % fissa di cellule per confronto
#   - corrisponde a una scelta conservativa ma visivamente utile
#
#  Una cellula è classificata "DC-like" se ha DC_score ≥ soglia.
# ============================================================
section("Classificazione DC-like (soglia = percentile 97.5)")

classify_dc <- function(obj, sample_name, percentile = 0.975) {
  scores    <- obj$DC_score
  threshold <- quantile(scores, percentile)
  dc_cells  <- colnames(obj)[scores >= threshold]
  n_dc      <- length(dc_cells)
  pct_dc    <- round(100 * n_dc / ncol(obj), 2)

  obj$dc_classification <- ifelse(scores >= threshold, "DC-like", "Other")

  cat(paste0(
    "\n[", sample_name, "]\n",
    "  Soglia DC score (p", percentile*100, "): ", round(threshold, 4), "\n",
    "  Cellule DC-like: ", n_dc, " / ", ncol(obj), " (", pct_dc, "%)\n"
  ))

  # Distribuzione DC-like per cell type annotato
  df <- data.frame(
    cell_type = as.character(obj$cell_type),
    dc_class  = obj$dc_classification
  )
  tbl <- df %>%
    group_by(cell_type) %>%
    summarise(
      n_total  = n(),
      n_dclike = sum(dc_class == "DC-like"),
      pct_dc   = round(100 * sum(dc_class == "DC-like") / n(), 2),
      .groups  = "drop"
    ) %>%
    arrange(desc(pct_dc))

  cat("  DC-like per cell type:\n")
  print(as.data.frame(tbl))

  return(list(obj = obj, dc_cells = dc_cells, threshold = threshold, table = tbl))
}

res_Ca <- classify_dc(obj_Ca, "Ca_I")
res_Me <- classify_dc(obj_Me, "Me_I")

obj_Ca    <- res_Ca$obj
obj_Me    <- res_Me$obj
dc_Ca     <- res_Ca$dc_cells
dc_Me     <- res_Me$dc_cells

# ============================================================
# 4. VISUALIZZAZIONI
# ============================================================

section("Generazione plot")

# ── Palette rainbow (stessa di STEP 3c per consistenza) ──────
rainbow_palette <- c(
  "Naive CD4+ T cells"         = "#E63946",
  "Effector CD4+ T cells"      = "#F4A261",
  "Memory T cells"             = "#2A9D8F",
  "Cytotoxic CD8+ T cells"     = "#264653",
  "Proliferating T cells"      = "#457B9D",
  "Proliferating CD8+ T cells" = "#6A0572",
  "Tregs"                      = "#E9C46A",
  "Dendritic Cells"            = "#E76F51"
)

get_colors <- function(obj) {
  types <- sort(unique(as.character(obj$cell_type)))
  cols  <- rainbow_palette[types]
  cols[is.na(cols)] <- "#9E9E9E"
  return(cols)
}

# ── A) UMAP highlight DC-like cells ──────────────────────────
# Mostra l'annotazione esistente in grigio chiaro + DC-like in fucsia
# cells.highlight mette in evidenza cellule specifiche sul UMAP

plot_dc_highlight <- function(obj, dc_cells, sample_name) {

  # Plot 1: highlight DC-like su sfondo grigio
  p_highlight <- DimPlot(
    obj,
    reduction       = "umap",
    cells.highlight = dc_cells,
    cols.highlight  = "#E91E63",   # fucsia per DC-like
    cols            = "#D3D3D3",   # grigio per le altre
    pt.size         = 0.7,
    sizes.highlight = 1.5
  ) +
    ggtitle(paste0(sample_name, " – DC-like cells (top 2.5% DC score)")) +
    theme_classic(base_size = 12) +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 12),
      legend.position = "right"
    ) +
    scale_color_manual(
      values = c("Highlighted" = "#E91E63", "Unhighlighted" = "#D3D3D3"),
      labels = c(paste0("DC-like (n=", length(dc_cells), ")"), "Other cells"),
      name   = ""
    )

  # Plot 2: annotazione rainbow + DC-like in nero sopra
  # Usiamo due layer: prima il DimPlot annotato, poi i punti DC sopra
  p_annot <- DimPlot(
    obj, reduction = "umap", cols = get_colors(obj),
    pt.size = 0.6, label = FALSE
  ) +
    ggtitle(paste0(sample_name, " – Annotazione + DC-like evidenziate")) +
    theme_classic(base_size = 12) +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 12),
      legend.text     = element_text(size = 8),
      legend.key.size = unit(0.4, "cm")
    ) +
    guides(color = guide_legend(override.aes = list(size = 3)))

  # Sovrapponi i punti DC in fucsia (geom_point sulle coordinate UMAP)
  umap_coords <- as.data.frame(Embeddings(obj, reduction = "umap"))
  dc_coords   <- umap_coords[dc_cells, , drop = FALSE]
  colnames(dc_coords) <- c("UMAP_1", "UMAP_2")

  p_annot <- p_annot +
    geom_point(
      data = dc_coords,
      aes(x = UMAP_1, y = UMAP_2),
      color  = "#E91E63",
      size   = 2,
      shape  = 16,
      inherit.aes = FALSE
    )

  combined <- p_highlight | p_annot
  path     <- paste0(out, sample_name, "_DC_highlight.png")
  ggsave(path, plot = combined, width = 16, height = 7,
         dpi = 300, bg = "white")
  cat(paste0("[", sample_name, "] UMAP highlight → ", path, "\n"))
  return(invisible(combined))
}

plot_dc_highlight(obj_Ca, dc_Ca, "Ca_I")
plot_dc_highlight(obj_Me, dc_Me, "Me_I")

# ── B) DC_score continuo sull'UMAP ───────────────────────────
plot_dc_score_umap <- function(obj, sample_name) {
  p_score <- FeaturePlot(
    obj, features = "DC_score", pt.size = 0.5,
    min.cutoff = "q05", max.cutoff = "q95",
    cols = c("lightgrey", "#C0392B")
  ) +
    ggtitle(paste0(sample_name, " – DC module score (gradiente)")) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  path <- paste0(out, sample_name, "_DC_score_umap.png")
  ggsave(path, plot = p_score, width = 7, height = 7,
         dpi = 300, bg = "white")
  cat(paste0("[", sample_name, "] DC score UMAP → ", path, "\n"))
  return(invisible(p_score))
}

plot_dc_score_umap(obj_Ca, "Ca_I")
plot_dc_score_umap(obj_Me, "Me_I")

# ── C) VlnPlot DC_score per cell type ────────────────────────
plot_dc_vln <- function(obj, sample_name) {
  cols <- get_colors(obj)

  p <- VlnPlot(
    obj, features = "DC_score",
    cols    = cols,
    pt.size = 0.1
  ) +
    ggtitle(paste0(sample_name, " – DC score per cell type")) +
    xlab("") + ylab("DC module score") +
    theme_classic(base_size = 10) +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold"),
      axis.text.x     = element_text(angle = 35, hjust = 1, size = 9),
      legend.position = "none"
    ) +
    geom_hline(
      yintercept = quantile(obj$DC_score, 0.975),
      linetype   = "dashed",
      color      = "#E91E63",
      linewidth  = 0.8
    ) +
    annotate("text", x = 0.6,
             y    = quantile(obj$DC_score, 0.975) + 0.02,
             label = "soglia p97.5",
             color = "#E91E63", size = 3)

  path <- paste0(out, sample_name, "_DC_score_vln.png")
  ggsave(path, plot = p, width = 8, height = 5, dpi = 300, bg = "white")
  cat(paste0("[", sample_name, "] VlnPlot DC score → ", path, "\n"))
  return(invisible(p))
}

plot_dc_vln(obj_Ca, "Ca_I")
plot_dc_vln(obj_Me, "Me_I")

# ── D) DotPlot marker DC per cell type ───────────────────────
# Mostra quali cell types esprimono marker DC
# Utile per identificare cluster con "contaminazione" mieloide

plot_dc_dotplot <- function(obj, sample_name) {
  genes_ok <- dc_genes[dc_genes %in% rownames(obj)]

  # Raggruppa i geni per categoria
  gene_order <- c(
    genes_ok[genes_ok %in% c("AQP9","VCAN","S100A12","LILRB2","MS4A6A","FPR1","HCK")],
    genes_ok[genes_ok %in% c("HLA-DRA","HLA-DPB1","HLA-DQA1")],
    genes_ok[genes_ok %in% c("CD1C","FCER1A","ITGAX")],
    genes_ok[genes_ok %in% c("CLEC9A","XCR1","BATF3")],
    genes_ok[genes_ok %in% c("LILRA4")],
    genes_ok[genes_ok %in% c("LYZ","S100A8","S100A9")]
  )
  gene_order <- unique(gene_order)

  p <- DotPlot(obj, features = gene_order) +
    RotatedAxis() +
    ggtitle(paste0(sample_name, " – Marker DC per cell type annotato")) +
    theme_classic(base_size = 10) +
    theme(
      plot.title  = element_text(hjust = 0.5, face = "bold", size = 11),
      axis.text.x = element_text(size = 8)
    ) +
    scale_color_gradient2(
      low      = "white",
      mid      = "#9B59B6",
      high     = "#C0392B",
      midpoint = 0,
      name     = "Avg\nExpression"
    )

  path <- paste0(out, sample_name, "_DC_DotPlot.png")
  ggsave(path, plot = p,
         width  = max(10, length(gene_order) * 0.45 + 3),
         height = max(4,  length(unique(Idents(obj))) * 0.55 + 2),
         dpi    = 300, bg = "white")
  cat(paste0("[", sample_name, "] DotPlot DC → ", path, "\n"))
  return(invisible(p))
}

plot_dc_dotplot(obj_Ca, "Ca_I")
plot_dc_dotplot(obj_Me, "Me_I")

# ── E) Pannello comparativo Ca_I vs Me_I ─────────────────────
# Un singolo plot che mette a confronto i DC score dei due campioni

df_ca <- data.frame(
  DC_score  = obj_Ca$DC_score,
  cell_type = as.character(obj_Ca$cell_type),
  sample    = "Ca_I",
  dc_class  = obj_Ca$dc_classification
)
df_me <- data.frame(
  DC_score  = obj_Me$DC_score,
  cell_type = as.character(obj_Me$cell_type),
  sample    = "Me_I",
  dc_class  = obj_Me$dc_classification
)
df_all <- rbind(df_ca, df_me)

p_comp <- ggplot(df_all, aes(x = cell_type, y = DC_score, fill = sample)) +
  geom_violin(scale = "width", alpha = 0.7, trim = TRUE) +
  geom_hline(
    data = data.frame(
      sample    = c("Ca_I", "Me_I"),
      threshold = c(res_Ca$threshold, res_Me$threshold)
    ),
    aes(yintercept = threshold, color = sample),
    linetype = "dashed", linewidth = 0.8, show.legend = FALSE
  ) +
  scale_fill_manual(values  = c("Ca_I" = "#457B9D", "Me_I" = "#E63946")) +
  scale_color_manual(values = c("Ca_I" = "#264653", "Me_I" = "#C62828")) +
  facet_wrap(~ sample, scales = "free_x", ncol = 2) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x  = element_text(angle = 40, hjust = 1, size = 8),
    strip.text   = element_text(face = "bold", size = 11),
    legend.position = "none",
    plot.title   = element_text(hjust = 0.5, face = "bold", size = 12)
  ) +
  labs(
    title = "DC module score per cell type (Ca_I vs Me_I)",
    x     = "",
    y     = "DC module score"
  )

path_comp <- paste0(out, "Ca_Me_DC_score_comparison.png")
ggsave(path_comp, plot = p_comp, width = 14, height = 6,
       dpi = 300, bg = "white")
cat(paste0("[COMPARATIVO] Ca_I vs Me_I → ", path_comp, "\n"))

# ============================================================
# 5. RIEPILOGO QUANTITATIVO FINALE
# ============================================================
section("RIEPILOGO QUANTITATIVO")

for (nm in c("Ca_I", "Me_I")) {
  res <- if (nm == "Ca_I") res_Ca else res_Me
  cat(paste0("\n[", nm, "]\n"))
  cat(paste0("  Soglia DC score (p97.5): ", round(res$threshold, 4), "\n"))
  cat(paste0("  Cellule DC-like: ",
             length(res$dc_cells), " (", 
             round(100*length(res$dc_cells)/
                     if(nm=="Ca_I") ncol(obj_Ca) else ncol(obj_Me), 2),
             "%)\n"))
  cat("  Distribuzione per cell type:\n")
  print(as.data.frame(res$table))
}

cat(paste0(
  "\n", strrep("=", 65), "\n",
  "  STEP 5 COMPLETATO\n\n",
  "  Output in: ", out, "\n\n",
  "  File generati:\n",
  "  Ca_I_DC_highlight.png    – UMAP con DC-like evidenziate\n",
  "  Me_I_DC_highlight.png    – UMAP con DC-like evidenziate\n",
  "  Ca_I_DC_score_umap.png   – DC score continuo\n",
  "  Me_I_DC_score_umap.png   – DC score continuo\n",
  "  Ca_I_DC_score_vln.png    – VlnPlot per cell type\n",
  "  Me_I_DC_score_vln.png    – VlnPlot per cell type\n",
  "  Ca_I_DC_DotPlot.png      – Marker DC per cell type\n",
  "  Me_I_DC_DotPlot.png      – Marker DC per cell type\n",
  "  Ca_Me_DC_score_comparison.png – Comparativo\n",
  strrep("=", 65), "\n"
))
