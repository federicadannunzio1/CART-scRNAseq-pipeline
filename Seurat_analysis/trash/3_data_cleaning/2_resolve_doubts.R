# ============================================================
#  STEP 2 – Risoluzione dubbi di annotazione
#
#  Prerequisito: STEP_1_load_data.R già eseguito
#  Input:  Bo_I, Ca_I, Me_I in environment
#  Output: oggetto `final_annotations` in environment
#          PNG diagnostici in base_dir/Doubt_resolution/
#
#  Risolve 4 dubbi:
#    D1 – Bo_I C0 vs C1: CD4 o CD8?
#    D2 – Ca_I C5: NK cells o secondo subset Treg?
#    D3 – Me_I C1 vs C3: stessa fase o fasi distinte del ciclo?
#    D4 – Bo_I C4: sottotipo DC (cDC1 / cDC2 / moDC / pDC)?
# ============================================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(scales)

# ── UNICO PUNTO DA MODIFICARE ────────────────────────────────
base_dir <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/CART/Code/Seurat_analysis/3_data_cleaning/"
# ─────────────────────────────────────────────────────────────

out <- paste0(base_dir, "Doubt_resolution/")
dir.create(out, showWarnings = FALSE, recursive = TRUE)

# ── Helper functions ──────────────────────────────────────────
save_plot <- function(p, filename, w = 14, h = 7) {
  path <- paste0(out, filename)
  ggsave(path, plot = p, width = w, height = h, dpi = 300, bg = "white")
  cat("  → Salvato:", path, "\n")
}

section <- function(title) {
  cat(paste0("\n", strrep("=", 65), "\n  ", title, "\n", strrep("=", 65), "\n"))
}

# ============================================================
# D1 – Bo_I C0 vs C1: CD4 o CD8?
# ============================================================
section("D1 | Bo_I – C0 vs C1: CD4 o CD8?")

Idents(Bo_I) <- "seurat_clusters"

# FindMarkers bidirezionale con soglia abbassata per catturare segnale debole
markers_01 <- FindMarkers(
  Bo_I, ident.1 = "0", ident.2 = "1",
  logfc.threshold = 0.1, min.pct = 0.1, only.pos = FALSE
)
markers_01$gene <- rownames(markers_01)

cd_genes <- c("CD4", "CD8A", "CD8B", "CD3D", "CD3E", "IL7R",
               "CCR7", "TCF7", "GZMK", "GZMB", "PRF1", "NKG7",
               "FOXP3", "CD44", "CD69", "SELL", "HLA-DRA")

cat("\nGeni T canonici differenziali tra C0 e C1",
    "(+ = più alto in C0, - = più alto in C1):\n")
markers_focus <- markers_01 %>%
  filter(gene %in% cd_genes) %>%
  arrange(avg_log2FC)
print(markers_focus[, c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")])

# Average expression per cluster
avg_expr <- AverageExpression(
  Bo_I,
  features = cd_genes[cd_genes %in% rownames(Bo_I)],
  group.by = "seurat_clusters", assay = "RNA", slot = "data"
)$RNA

# ── FIX Seurat v5: AverageExpression può restituire nomi colonna con
#    prefisso (es. "g0","g1") invece di "0","1".
#    Normalizziamo sempre i nomi colonna ai numeri puri del cluster.
colnames(avg_expr) <- gsub("^g", "", colnames(avg_expr))   # rimuove "g" iniziale
colnames(avg_expr) <- gsub("^RNA_snn_res\\.[0-9.]+_", "",  # rimuove prefisso res
                            colnames(avg_expr))

cat("\nNomi colonna avg_expr dopo normalizzazione:", paste(colnames(avg_expr), collapse=", "), "\n")
cat("\nMedia espressione per cluster:\n")
print(round(avg_expr, 3))

# Helper sicuro per estrarre un valore da avg_expr
get_avg <- function(mat, gene, cluster) {
  if (!gene    %in% rownames(mat)) return(0)
  if (!cluster %in% colnames(mat)) return(0)
  mat[gene, cluster]
}

# VlnPlot + FeaturePlot
key_d1 <- c("CD4", "CD8A", "CD8B", "IL7R", "GZMK", "GZMB", "CCR7", "NKG7")
key_d1 <- key_d1[key_d1 %in% rownames(Bo_I)]

p_vln <- VlnPlot(Bo_I, features = key_d1, idents = c("0", "1"),
                 ncol = 4, pt.size = 0, cols = c("#4E9AF1", "#E87D3E")) &
  theme_classic(base_size = 10) &
  theme(axis.title.x = element_blank(),
        plot.title   = element_text(face = "bold", size = 10))

p_feat <- FeaturePlot(Bo_I, features = key_d1, ncol = 4, pt.size = 0.3,
                      min.cutoff = "q05", cols = c("lightgrey", "#C0392B")) &
  theme_classic(base_size = 9)

save_plot(
  (p_vln / p_feat) +
    plot_annotation(
      title    = "D1 | Bo_I – C0 vs C1: CD4 vs CD8",
      subtitle = "VlnPlot (sopra): C0=blu, C1=arancio | FeaturePlot (sotto): tutti i cluster",
      theme    = theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5))
    ),
  "D1_BoI_C0vsC1_CD4vsCD8.png", w = 16, h = 14
)

# Decisione automatica basata su average expression
avg_c0_cd4 <- get_avg(avg_expr, "CD4",  "0")
avg_c1_cd4 <- get_avg(avg_expr, "CD4",  "1")
avg_c0_cd8 <- get_avg(avg_expr, "CD8A", "0")
avg_c1_cd8 <- get_avg(avg_expr, "CD8A", "1")

d1_label_c0 <- if (avg_c0_cd4 >= avg_c0_cd8) "CD4+ T cells" else "CD8+ T cells"
d1_label_c1 <- if (avg_c1_cd4 >= avg_c1_cd8) "CD4+ T cells" else "CD8+ T cells"
# Fallback: se entrambi i valori sono quasi zero, assegna per esclusione
if (max(avg_c0_cd4, avg_c0_cd8) < 0.05) d1_label_c0 <- "CD4+ T cells"
if (max(avg_c1_cd4, avg_c1_cd8) < 0.05) d1_label_c1 <- "CD8+ T cells"

cat(paste0("\n[DECISIONE D1]\n",
           "  C0 → ", d1_label_c0, "  (CD4=", round(avg_c0_cd4,3),
           " | CD8A=", round(avg_c0_cd8,3), ")\n",
           "  C1 → ", d1_label_c1, "  (CD4=", round(avg_c1_cd4,3),
           " | CD8A=", round(avg_c1_cd8,3), ")\n"))

# ============================================================
# D2 – Ca_I C5: NK cells o secondo subset Treg?
# ============================================================
section("D2 | Ca_I – C5: NK cells o Treg?")

Idents(Ca_I) <- "seurat_clusters"

nk_genes   <- c("NKG7", "NCAM1", "NCR1", "KLRD1", "KLRB1", "GNLY",
                 "FCGR3A", "TYROBP", "CD247", "FCER1G", "S1PR5",
                 "KLRC1", "KLRC2", "KLRF1", "CX3CR1", "GZMB", "PRF1")
treg_genes <- c("FOXP3", "IL2RA", "CTLA4", "TIGIT", "IKZF2",
                "TNFRSF18", "TNFRSF4", "ENTPD1", "BATF", "RTKN2",
                "CCR8", "LAYN", "DUSP4")

nk_present   <- nk_genes[nk_genes %in% rownames(Ca_I)]
treg_present <- treg_genes[treg_genes %in% rownames(Ca_I)]

# Module scores
Ca_I <- AddModuleScore(Ca_I, features = list(nk_present),   name = "NK_score")
Ca_I <- AddModuleScore(Ca_I, features = list(treg_present), name = "Treg_score")
Ca_I$NK_score   <- Ca_I$NK_score1;   Ca_I$NK_score1   <- NULL
Ca_I$Treg_score <- Ca_I$Treg_score1; Ca_I$Treg_score1 <- NULL

score_summary <- Ca_I@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(mean_NK   = round(mean(NK_score), 3),
            mean_Treg = round(mean(Treg_score), 3),
            n = n())
cat("\nNK score vs Treg score per cluster:\n")
print(as.data.frame(score_summary))

# Plots
p_vln_score <- VlnPlot(Ca_I, features = c("NK_score", "Treg_score"),
                        idents = c("4", "5"), ncol = 2, pt.size = 0.2,
                        cols = c("#F39C12", "#16A085")) &
  theme_classic(base_size = 10) &
  theme(axis.title.x = element_blank(),
        plot.title   = element_text(face = "bold", size = 10))

nk_feat   <- c("NKG7","NCAM1","GNLY","KLRD1","KLRB1","NCR1","S1PR5","CX3CR1")
treg_feat <- c("FOXP3","IL2RA","CTLA4","TIGIT","IKZF2","TNFRSF18","CCR8","RTKN2")

p_feat_nk <- FeaturePlot(Ca_I, features = nk_feat[nk_feat %in% rownames(Ca_I)],
                          ncol = 4, pt.size = 0.3, min.cutoff = "q05",
                          cols = c("lightgrey", "#16A085")) &
  theme_classic(base_size = 9)

p_feat_treg <- FeaturePlot(Ca_I, features = treg_feat[treg_feat %in% rownames(Ca_I)],
                            ncol = 4, pt.size = 0.3, min.cutoff = "q05",
                            cols = c("lightgrey", "#F39C12")) &
  theme_classic(base_size = 9)

save_plot(
  (p_vln_score / p_feat_nk / p_feat_treg) +
    plot_annotation(
      title    = "D2 | Ca_I – C5: NK cells vs Treg?",
      subtitle = "Score C4 vs C5 (sopra) | geni NK (centro) | geni Treg (sotto)",
      theme    = theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5))
    ),
  "D2_CaI_C5_NKvsTreg.png", w = 16, h = 20
)

# Decisione automatica
c5_nk   <- score_summary$mean_NK[score_summary$seurat_clusters == "5"]
c5_treg <- score_summary$mean_Treg[score_summary$seurat_clusters == "5"]
d2_label <- if (c5_nk > c5_treg) "NK cells" else "Tregs"
cat(paste0("\n[DECISIONE D2]\n",
           "  C5 NK score:   ", c5_nk, "\n",
           "  C5 Treg score: ", c5_treg, "\n",
           "  C5 → ", d2_label, "\n"))

# ============================================================
# D3 – Me_I C1 vs C3: fasi del ciclo cellulare distinte?
# ============================================================
section("D3 | Me_I – C1 vs C3: stessa fase del ciclo cellulare?")

Idents(Me_I) <- "seurat_clusters"

s_ok   <- cc.genes$s.genes[cc.genes$s.genes %in% rownames(Me_I)]
g2m_ok <- cc.genes$g2m.genes[cc.genes$g2m.genes %in% rownames(Me_I)]
cat(paste0("Geni S presenti: ", length(s_ok), " | G2M presenti: ", length(g2m_ok), "\n"))

Me_I <- CellCycleScoring(Me_I, s.features = s_ok, g2m.features = g2m_ok,
                          set.ident = FALSE)

phase_dist <- Me_I@meta.data %>%
  filter(seurat_clusters %in% c("1", "3")) %>%
  group_by(seurat_clusters, Phase) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(seurat_clusters) %>%
  mutate(pct = round(100 * n / sum(n), 1))

cat("\nDistribuzione fasi del ciclo in C1 e C3:\n")
print(as.data.frame(phase_dist))

cat("\nScore medio per tutti i cluster:\n")
print(as.data.frame(
  Me_I@meta.data %>%
    group_by(seurat_clusters) %>%
    summarise(mean_S   = round(mean(S.Score), 3),
              mean_G2M = round(mean(G2M.Score), 3),
              pct_S    = round(100 * mean(Phase == "S"), 1),
              pct_G2M  = round(100 * mean(Phase == "G2M"), 1),
              pct_G1   = round(100 * mean(Phase == "G1"), 1))
))

# Plots
p_umap_phase <- DimPlot(Me_I, group.by = "Phase", pt.size = 0.5,
                         cols = c("G1"="#BDC3C7","S"="#3498DB","G2M"="#E74C3C")) +
  ggtitle("Me_I – Cell cycle phase") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

p_scatter <- ggplot(
  Me_I@meta.data %>% filter(seurat_clusters %in% c("1", "3")),
  aes(x = S.Score, y = G2M.Score, color = Phase, shape = seurat_clusters)
) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("G1"="#BDC3C7","S"="#3498DB","G2M"="#E74C3C")) +
  theme_classic(base_size = 11) +
  labs(title = "C1 vs C3 – S.Score vs G2M.Score",
       shape = "Cluster", color = "Phase") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

p_bar <- phase_dist %>%
  ggplot(aes(x = seurat_clusters, y = pct, fill = Phase)) +
  geom_col(position = "fill", width = 0.6) +
  scale_fill_manual(values = c("G1"="#BDC3C7","S"="#3498DB","G2M"="#E74C3C")) +
  scale_y_continuous(labels = percent_format()) +
  theme_classic(base_size = 11) +
  labs(title = "Proporzione fasi C1 vs C3", x = "Cluster", y = "% cellule") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

cc_diag <- c("MKI67","TOP2A","STMN1","PCNA","MCM5","CCNB1","CDK1","CCNA2","BUB1B","TYMS")
cc_diag <- cc_diag[cc_diag %in% rownames(Me_I)]
p_vln_cc <- VlnPlot(Me_I, features = cc_diag, idents = c("1","3"),
                     ncol = 5, pt.size = 0,
                     cols = c("#3498DB","#E74C3C")) &
  theme_classic(base_size = 9) &
  theme(axis.title.x = element_blank(),
        plot.title   = element_text(face = "bold", size = 9))

save_plot(
  (p_umap_phase | p_scatter | p_bar) / p_vln_cc +
    plot_annotation(
      title    = "D3 | Me_I – C1 vs C3: fasi del ciclo cellulare",
      subtitle = "Blu = C1, Rosso = C3 nel VlnPlot",
      theme    = theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5))
    ) +
    plot_layout(heights = c(1.2, 1)),
  "D3_MeI_C1vsC3_CellCycle.png", w = 18, h = 13
)

# Decisione automatica: label distinto solo se fase dominante è diversa
c1_dom <- phase_dist %>% filter(seurat_clusters == "1") %>%
  slice_max(pct, n = 1) %>% pull(Phase)
c3_dom <- phase_dist %>% filter(seurat_clusters == "3") %>%
  slice_max(pct, n = 1) %>% pull(Phase)

if (c1_dom != c3_dom) {
  d3_label_c1 <- paste0("Proliferating T cells (", c1_dom, ")")
  d3_label_c3 <- paste0("Proliferating T cells (", c3_dom, ")")
} else {
  d3_label_c1 <- "Proliferating T cells"
  d3_label_c3 <- "Proliferating T cells"
}
cat(paste0("\n[DECISIONE D3]\n",
           "  C1 fase dominante: ", c1_dom, " → ", d3_label_c1, "\n",
           "  C3 fase dominante: ", c3_dom, " → ", d3_label_c3, "\n"))

# ============================================================
# D4 – Bo_I C4: sottotipo DC
# ============================================================
section("D4 | Bo_I – C4: sottotipo DC (cDC1 / cDC2 / moDC / pDC)?")

Idents(Bo_I) <- "seurat_clusters"

dc_subtypes <- list(
  cDC1 = c("CLEC9A","XCR1","CADM1","BATF3","IDO1","WDFY4","SNX22","CPNE3","CLNK"),
  cDC2 = c("CD1C","FCER1A","CLEC10A","FCGR2B","CX3CR1","ITGAX","CD14","S100A4","ANXA2"),
  moDC = c("VCAN","S100A12","S100A8","S100A9","FCN1","CXCL8","THBS1",
            "AQP9","PLBD1","FPR1","LILRB2","MS4A6A","HCK"),
  pDC  = c("LILRA4","CLEC4C","PTGDS","JCHAIN","TCF4","IL3RA","GZMB","ITM2C","MZB1")
)

all_dc_genes <- unique(unlist(dc_subtypes))
all_dc_genes <- all_dc_genes[all_dc_genes %in% rownames(Bo_I)]

avg_bo <- AverageExpression(Bo_I, features = all_dc_genes,
                             group.by = "seurat_clusters",
                             assay = "RNA", slot = "data")$RNA
# Normalizza nomi colonna (stesso fix di D1)
colnames(avg_bo) <- gsub("^g", "", colnames(avg_bo))
colnames(avg_bo) <- gsub("^RNA_snn_res\\.[0-9.]+_", "", colnames(avg_bo))
cat("\nAverage expression geni DC per cluster:\n")
print(round(avg_bo, 4))

# Module score per sottotipo
score_cols_computed <- c()
for (st in names(dc_subtypes)) {
  genes_ok <- dc_subtypes[[st]][dc_subtypes[[st]] %in% rownames(Bo_I)]
  if (length(genes_ok) < 2) {
    cat(paste0("  [WARN] ", st, ": solo ", length(genes_ok), " geni presenti. Salto.\n"))
    next
  }
  Bo_I <- AddModuleScore(Bo_I, features = list(genes_ok),
                          name = paste0(st, "_score"))
  # AddModuleScore aggiunge "1" al nome → rinomina
  Bo_I[[paste0(st, "_score")]] <- Bo_I[[paste0(st, "_score1")]]
  Bo_I[[paste0(st, "_score1")]] <- NULL
  score_cols_computed <- c(score_cols_computed, paste0(st, "_score"))
}

score_dc_summary <- Bo_I@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(across(all_of(score_cols_computed), ~ round(mean(.x), 4)), n = n())
cat("\nScore DC per sottotipo (media per cluster):\n")
print(as.data.frame(score_dc_summary))

# Plots
p_vln_dc <- VlnPlot(Bo_I, features = score_cols_computed, idents = "4",
                     ncol = 2, pt.size = 0.5, cols = "#E74C3C") &
  theme_classic(base_size = 10) &
  theme(axis.title.x = element_blank(),
        plot.title   = element_text(face = "bold", size = 10),
        legend.position = "none")

top_genes_dc <- c("CLEC9A","XCR1","CD1C","FCER1A","VCAN","S100A12","LILRA4","CLEC4C")
top_genes_dc <- top_genes_dc[top_genes_dc %in% rownames(Bo_I)]

p_feat_dc <- FeaturePlot(Bo_I, features = top_genes_dc, ncol = 4, pt.size = 0.5,
                          min.cutoff = "q05", cols = c("lightgrey","#E74C3C")) &
  theme_classic(base_size = 9)

p_dot_dc <- DotPlot(Bo_I, features = rev(all_dc_genes)) +
  RotatedAxis() +
  ggtitle("Bo_I – Geni DC per sottotipo (tutti i cluster)") +
  theme_classic(base_size = 9) +
  theme(plot.title  = element_text(face = "bold", hjust = 0.5, size = 11),
        axis.text.x = element_text(size = 7)) +
  scale_color_gradient2(low = "white", mid = "#9B59B6", high = "#C0392B", midpoint = 0)

save_plot(
  (p_vln_dc | p_feat_dc) / p_dot_dc +
    plot_annotation(
      title    = "D4 | Bo_I – C4: sottotipo DC",
      subtitle = "VlnPlot score (sinistra) | FeaturePlot geni chiave (destra) | DotPlot tutti cluster (sotto)",
      theme    = theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5))
    ) +
    plot_layout(heights = c(1, 1.4)),
  "D4_BoI_C4_DC_subtype.png", w = 18, h = 18
)

# Decisione automatica: sottotipo con score più alto in C4
c4_scores <- score_dc_summary %>%
  filter(seurat_clusters == "4") %>%
  select(all_of(score_cols_computed))

dc_label_map <- c(
  cDC1_score = "cDC1 (Dendritic Cells)",
  cDC2_score = "cDC2 (Dendritic Cells)",
  moDC_score = "Monocyte-derived Dendritic Cells",
  pDC_score  = "Plasmacytoid Dendritic Cells"
)

if (ncol(c4_scores) > 0 && nrow(c4_scores) > 0) {
  best_col  <- names(which.max(c4_scores[1, ]))
  d4_label  <- dc_label_map[best_col]
  cat(paste0("\n[DECISIONE D4]\n",
             "  Score massimo in C4: ", best_col, " = ",
             round(c4_scores[[best_col]], 4), "\n",
             "  C4 → ", d4_label, "\n"))
} else {
  d4_label <- "Dendritic Cells"
  cat("[DECISIONE D4] Score non calcolabili → label generico: Dendritic Cells\n")
}

# ============================================================
# RIEPILOGO – Costruzione di final_annotations
# ============================================================
section("RIEPILOGO DECISIONI E DIZIONARI FINALI")

# unname() è necessario su d4_label e d2_label perché provengono da
# subsetting di named vector (dc_label_map[best_subtype] e
# if/else su score_summary): senza unname() portano il nome originale
# della chiave (es. "moDC_score") che sovrascrive il nome del cluster
# nella lista, rendendo i cluster 4/5 non trovabili da annotate_sample.
annotation_Bo <- c(
  "0" = unname(d1_label_c0),
  "1" = unname(d1_label_c1),
  "2" = "Cytotoxic CD8+ T cells",
  "3" = "Memory T cells",
  "4" = unname(d4_label)
)

annotation_Ca <- c(
  "0" = "Proliferating T cells",
  "1" = "Naive CD4+ T cells",
  "2" = "Effector CD4+ T cells",
  "3" = "Cytotoxic CD8+ T cells",
  "4" = "Tregs",
  "5" = unname(d2_label)
)

annotation_Me <- c(
  "0" = "Cytotoxic CD8+ T cells",
  "1" = unname(d3_label_c1),
  "2" = "Tregs",
  "3" = unname(d3_label_c3),
  "4" = "Proliferating CD8+ T cells"
)

cat("\n--- Bo_I ---\n"); print(annotation_Bo)
cat("  Chiavi:", paste(names(annotation_Bo), collapse=", "), "\n")
cat("\n--- Ca_I ---\n"); print(annotation_Ca)
cat("  Chiavi:", paste(names(annotation_Ca), collapse=", "), "\n")
cat("\n--- Me_I ---\n"); print(annotation_Me)
cat("  Chiavi:", paste(names(annotation_Me), collapse=", "), "\n")

final_annotations <- list(
  Bo_I = annotation_Bo,
  Ca_I = annotation_Ca,
  Me_I = annotation_Me
)

cat(paste0(
  "\n", strrep("=", 55), "\n",
  "  STEP 2 COMPLETATO\n",
  "  Oggetto 'final_annotations' disponibile in environment.\n",
  "  Plot diagnostici salvati in: ", out, "\n",
  "  Prossimo step: esegui STEP_3_annotate_and_plot.R\n",
  strrep("=", 55), "\n"
))
