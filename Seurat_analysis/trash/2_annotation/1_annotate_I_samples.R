
# ============================================================
#  PIPELINE 1 – Campioni I (infusione)
#  Copre: 0_process_data + 1_load_data + 2_resolve_doubts +
#         3_annotate_and_plot + 3b_update_annotations +
#         ricerca cluster mieloidi in TUTTI i campioni I
#
#  Sostituisce: script 0-3b + script 4 + script 5
#  Modifiche rispetto alla pipeline originale:
#    - D4 (sottotipo DC) → ricerca mieloidi (non DC)
#    - Firma mieloide ampliata (LYZ/S100A8/VCAN/CSF1R/CD68)
#    - Label "Dendritic Cells" → "Myeloid cells" ovunque
#    - Th helper subtypes (Th1/Th2/Th17/Tfh) nel ramo CD4
#    - Legenda consistente: tutti i tipi in ordine biologico
#      fisso; tipi assenti mostrati in grigio con "—"
#    - Ricerca mieloidi estesa a Ca_I e Me_I
#
#  Output (tutti in base_dir/):
#    Pipeline_I/
#      Doubt_resolution/   – plot diagnostici D1-D3 + mieloidi
#      Annotation_UMAP/    – UMAP finali con legenda consistente
#      Myeloid_search/     – FeaturePlot + VlnPlot mieloidi
#    all_I_samples_annotated.rds
# ============================================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(openxlsx)
library(scales)

# ── UNICO PUNTO DA MODIFICARE ────────────────────────────────
base_dir <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/CART/Code/Seurat_analysis/2_annotation/"
seurat_list_path <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/CART/Data/seurat_obj_list/seurat_samples_sctype_azimuth_pbmc_bonemarrow_clonalvdj_CAR.rds"
# ─────────────────────────────────────────────────────────────

out_doubt  <- paste0(base_dir, "Pipeline_I/Doubt_resolution/")
out_umap   <- paste0(base_dir, "Pipeline_I/Annotation_UMAP/")
out_mye    <- paste0(base_dir, "Pipeline_I/Myeloid_search/")
for (d in c(out_doubt, out_umap, out_mye))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

section <- function(title)
  cat(paste0("\n", strrep("=", 65), "\n  ", title, "\n",
             strrep("=", 65), "\n"))

save_plot <- function(p, path, w = 16, h = 7) {
  ggsave(path, plot = p, width = w, height = h,
         dpi = 300, bg = "white")
  cat("  → Salvato:", path, "\n")
}

# ============================================================
# STEP 0 – PROCESSING RAW I SAMPLES
# ============================================================
section("STEP 0 | Processing campioni I")

if (!exists("seurat_list")) {
  cat("Caricamento seurat_list...\n")
  seurat_list <- readRDS(seurat_list_path)
}

raw_I <- list(
  Bo_I = seurat_list$Bo_samples_blood$I,
  Ca_I = seurat_list$Ca_samples_blood$I,
  Me_I = seurat_list$Me_samples_bone$I
)

params_I <- list(
  Bo_I = list(dims = 1:20, res = 0.7),
  Ca_I = list(dims = 1:30, res = 0.7),
  Me_I = list(dims = 1:30, res = 0.5)
)

# Controlla se i campioni processati esistono già su disco
rds_preproc <- paste0(base_dir, "all_samples_clean_pre_annotation.rds")

if (file.exists(rds_preproc)) {
  cat(">> Trovato all_samples_clean_pre_annotation.rds – salto il processing.\n")
  processed_I <- readRDS(rds_preproc)
} else {
  cat(">> File non trovato – processing da zero.\n")

  process_sample <- function(s_name, raw_obj, p) {
    cat(paste0("\n  Processing: ", s_name,
               " | dims 1:", max(p$dims), " | res ", p$res, "\n"))

    obj <- subset(raw_obj,
                  subset = nFeature_RNA > 800 &
                           percent.mt < 7 &
                           azimuth_class != "dnT")
    obj <- NormalizeData(obj, verbose = FALSE)
    obj <- FindVariableFeatures(obj, selection.method = "vst",
                                nfeatures = 2000, verbose = FALSE)
    stress <- c("MALAT1","NEAT1",
                grep("^MT-", rownames(obj), value = TRUE))
    VariableFeatures(obj) <- setdiff(VariableFeatures(obj), stress)
    obj <- ScaleData(obj, vars.to.regress = "percent.mt",
                     verbose = FALSE)
    obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
    obj <- RunUMAP(obj, dims = p$dims, verbose = FALSE)
    obj <- FindNeighbors(obj, dims = p$dims, verbose = FALSE)
    obj <- FindClusters(obj, resolution = p$res, verbose = FALSE)

    out_xls <- paste0(base_dir, s_name, "/")
    dir.create(out_xls, showWarnings = FALSE, recursive = TRUE)
    Idents(obj) <- "seurat_clusters"
    markers <- FindAllMarkers(obj, only.pos = FALSE,
                              min.pct = 0.25,
                              logfc.threshold = 0.25,
                              verbose = FALSE)
    write.xlsx(markers,
               paste0(out_xls, "Markers_", s_name, ".xlsx"))
    return(obj)
  }

  processed_I <- list()
  for (s in names(raw_I))
    processed_I[[s]] <- process_sample(s, raw_I[[s]], params_I[[s]])

  saveRDS(processed_I, rds_preproc)
  cat(">> Salvato:", rds_preproc, "\n")
}

# ── Estrazione e JoinLayers ───────────────────────────────────
join_safe <- function(obj, nm) {
  if (length(grep("^counts\\.", Layers(obj), value = TRUE)) > 0)
    obj <- JoinLayers(obj)
  Idents(obj) <- "seurat_clusters"
  cat(sprintf("[%s] %d cellule | %d cluster\n",
              nm, ncol(obj),
              length(unique(obj$seurat_clusters))))
  obj
}

Bo_I <- join_safe(processed_I$Bo_I, "Bo_I")
Ca_I <- join_safe(processed_I$Ca_I, "Ca_I")
Me_I <- join_safe(processed_I$Me_I, "Me_I")

# ============================================================
# FIRME GENICHE CONDIVISE  (usate in D1-D3, mieloidi, Th)
# ============================================================

SIG <- list(
  cd3      = c("CD3D","CD3E","CD3G"),
  cd4      = c("CD4"),
  cd8      = c("CD8A","CD8B"),
  cd19     = c("CD19","MS4A1"),
  cd56     = c("NCAM1"),
  cd14     = c("CD14"),
  cd16     = c("FCGR3A"),

  naive    = c("CCR7","SELL","IL7R","TCF7","LEF1","KLF2"),
  effector = c("GZMK","CD44","S100A4","LGALS1","DUSP2"),
  cytotox  = c("GZMB","NKG7","PRF1","GNLY","GZMA","FGFBP2"),
  treg     = c("FOXP3","IL2RA","CTLA4","TIGIT","IKZF2",
               "TNFRSF18","TNFRSF4","ENTPD1","CCR8","RTKN2"),
  prolif   = c("MKI67","TOP2A","STMN1","PCNA","CCNB1"),
  nk       = c("NCAM1","NCR1","KLRD1","GNLY","TYROBP",
               "S1PR5","CX3CR1","XCL1","FCGR3A"),

  # T helper subtypes
  # Th1: TBX21 (T-bet) + CXCR3 i più stabili a riposo
  th1  = c("TBX21","CXCR3","CCR5","IL12RB2","STAT4",
           "IFNG","HAVCR2","PHLPP1","TNFSF10"),
  # Th2: GATA3 presente anche NK/ILC → protetto da cd4 > 0.3
  th2  = c("GATA3","CCR4","MAF","PTGDR2","IL4R",
           "IL4","IL13","HPGDS"),
  # Th17: RORC + CCR6 affidabili; IL17A spesso bassa a riposo
  th17 = c("RORC","CCR6","IL23R","RORA","FURIN",
           "IL17A","IL17F","TMEM176A","TMEM176B","STAT3"),
  # Tfh: CXCR5 + BCL6 più specifici in scRNA-seq
  tfh  = c("CXCR5","BCL6","ICOS","PDCD1","SH2D1A",
           "IL21","CXCL13","TOX2","MAF","TIGIT"),

  # Myeloid cells – NON dendritic cells.
  # Bo_I C4 esprimeva: LYZ, VCAN, S100A8/9, AQP9, FPR1, HCK,
  # LILRB2, MS4A6A – tutti marker mieloidi generici/moDC.
  # I marker DC classici (CLEC9A, XCR1, CD1C, LILRA4) erano
  # assenti o minimi: il cluster è mieloide, non DC.
  myeloid = c(
    "LYZ","VCAN","S100A8","S100A9","AQP9","FPR1",
    "LILRB2","MS4A6A","HCK","CSF1R","CD68","ITGAM",
    "FCN1","CXCL8","THBS1","PLBD1","TYROBP","SPI1",
    "MRC1","C1QA","C1QB","C1QC",
    "HLA-DRA","HLA-DPB1","HLA-DQA1",
    "CD1C","FCER1A","CLEC9A","XCR1","LILRA4",
    "ITGAX","BATF3","CLEC10A","SIGLEC6"
  ),

  # NK
  nkt  = c("NCAM1","CD3D","KLRB1","ZBTB16","NKG7","GNLY","GZMB")
)

# ── helper per average expression ────────────────────────────
get_avg <- function(mat, genes, cluster) {
  g <- genes[genes %in% rownames(mat)]
  if (length(g) == 0 || !cluster %in% colnames(mat)) return(0)
  mean(mat[g, cluster])
}

fix_colnames <- function(mat) {
  colnames(mat) <- gsub("^g", "", colnames(mat))
  colnames(mat) <- gsub("^RNA_snn_res\\.[0-9.]+_", "", colnames(mat))
  mat
}

# ============================================================
# STEP 1 – PALETTE E LEGENDA CONSISTENTE
# ============================================================
section("STEP 1 | Palette e legenda")

# Colori fissi: invariati tra tutti i campioni I e AB
PALETTE <- c(
  "Naive CD4+ T cells"         = "#E63946",
  "Th1 cells"                  = "#C1121F",
  "Th2 cells"                  = "#FF99C8",
  "Th17 cells"                 = "#FB5607",
  "Tfh cells"                  = "#800F2F",
  "Effector CD4+ T cells"      = "#F4A261",
  "Memory T cells"             = "#2A9D8F",
  "Tregs"                      = "#E9C46A",
  "Cytotoxic CD8+ T cells"     = "#264653",
  "Naive CD8+ T cells"         = "#577590",
  "Proliferating T cells"      = "#457B9D",
  "Proliferating CD4+ T cells" = "#023E8A",
  "Proliferating CD8+ T cells" = "#6A0572",
  "NK cells"                   = "#43AA8B",
  "NKT cells"                  = "#277DA1",
  "gamma-delta T cells"        = "#4D908E",
  "MAIT cells"                 = "#F3722C",
  "ILC"                        = "#F9C74F",
  "B cells"                    = "#90BE6D",
  "Memory B cells"             = "#52B788",
  "Plasma cells"               = "#C77DFF",
  "CD14 Monocytes"             = "#9D0208",
  "CD16 Monocytes"             = "#DC2F02",
  "Myeloid cells"              = "#E76F51",
  "Basophils"                  = "#6A4C93",
  "HSPC"                       = "#B5838D",
  "Erythroid cells"            = "#FFAFCC",
  "Platelets"                  = "#CDB4DB"
)

# Ordine biologico canonico – invariato in tutti i pannelli
CANONICAL_ORDER <- c(
  "Naive CD4+ T cells","Th1 cells","Th2 cells","Th17 cells","Tfh cells",
  "Effector CD4+ T cells","Memory T cells","Tregs",
  "Cytotoxic CD8+ T cells","Naive CD8+ T cells",
  "Proliferating T cells","Proliferating CD4+ T cells","Proliferating CD8+ T cells",
  "NKT cells","gamma-delta T cells","MAIT cells",
  "NK cells","ILC",
  "B cells","Memory B cells","Plasma cells",
  "CD14 Monocytes","CD16 Monocytes","Myeloid cells","Basophils",
  "HSPC","Erythroid cells","Platelets"
)

get_colors <- function(present_types) {
  cols    <- PALETTE[present_types]
  missing <- present_types[is.na(cols)]
  if (length(missing) > 0) {
    extra <- setNames(hue_pal()(length(missing)), missing)
    cols[missing] <- extra
    PALETTE[missing] <<- extra
    cat("[INFO] Colori auto-generati:", paste(missing, collapse=", "), "\n")
  }
  cols
}

# Pannello legenda con tipi assenti segnalati con "—"
make_full_legend <- function(present_types) {
  all_t   <- c(CANONICAL_ORDER,
               setdiff(names(PALETTE), CANONICAL_ORDER))
  all_t   <- all_t[all_t %in% names(PALETTE)]
  n       <- length(all_t)
  present <- all_t %in% present_types

  df <- data.frame(
    y       = rev(seq_len(n)),
    label   = ifelse(present, all_t, paste0(all_t, " \u2014")),
    pt_col  = ifelse(present, PALETTE[all_t], "#CCCCCC"),
    txt_col = ifelse(present, "black", "gray60"),
    alpha   = ifelse(present, 1.0, 0.45),
    stringsAsFactors = FALSE
  )

  ggplot(df) +
    geom_point(aes(x = 0, y = y), color = df$pt_col,
               alpha = df$alpha, size = 3, shape = 16) +
    geom_text(aes(x = 0.2, y = y, label = label),
              hjust = 0, size = 2.6, color = df$txt_col) +
    xlim(-0.2, 5) +
    theme_void() +
    theme(plot.background = element_rect(fill = "white",
                                         color = NA),
          plot.margin = margin(8, 4, 8, 4))
}

# UMAP con legenda esterna consistente
plot_umap_consistent <- function(obj, sample_name, out_path) {
  Idents(obj) <- "cell_type"
  present <- sort(unique(as.character(obj$cell_type)))
  cols    <- get_colors(present)

  make_dp <- function(label) {
    DimPlot(obj, reduction = "umap", group.by = "cell_type",
            label = label, label.size = 3.2, repel = TRUE,
            cols = cols, pt.size = 0.6) +
      ggtitle(paste0(sample_name,
                     ifelse(label, " – Con label", " – Senza label"))) +
      theme_classic(base_size = 12) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold",
                                      size = 13),
            legend.position = "none")
  }

  combined <- (make_dp(TRUE) | make_dp(FALSE) | make_full_legend(present)) +
    plot_layout(widths = c(5, 5, 2.5)) +
    plot_annotation(title = sample_name,
                    theme = theme(plot.title =
                                    element_text(face = "bold",
                                                 hjust = 0.5,
                                                 size = 14)))
  ggsave(out_path, plot = combined,
         width = 20, height = 7, dpi = 300, bg = "white")
  cat("  → UMAP:", out_path, "\n")
  return(invisible(combined))
}

# ============================================================
# STEP 2 – RISOLUZIONE DUBBI D1-D3
# (D4 rimosso: era analisi DC → sostituita dalla ricerca
#  mieloidi in STEP 3)
# ============================================================
section("STEP 2 | Risoluzione dubbi D1-D3")

# ── D1 – Bo_I C0 vs C1: Naive vs Effector CD4+ ──────────────
cat("\n[D1] Bo_I C0 vs C1 – Naive vs Effector CD4+\n")
Idents(Bo_I) <- "seurat_clusters"

naive_eff_genes <- unique(c(SIG$naive, SIG$effector))
naive_eff_genes <- naive_eff_genes[naive_eff_genes %in% rownames(Bo_I)]

avg_d1 <- fix_colnames(
  AverageExpression(Bo_I, features = naive_eff_genes,
                    group.by = "seurat_clusters",
                    assay = "RNA", slot = "data")$RNA)

s_naive_c0    <- get_avg(avg_d1, SIG$naive,    "0")
s_effector_c0 <- get_avg(avg_d1, SIG$effector, "0")
s_naive_c1    <- get_avg(avg_d1, SIG$naive,    "1")
s_effector_c1 <- get_avg(avg_d1, SIG$effector, "1")

label_Bo_c0 <- if (s_naive_c0 >= s_effector_c0) "Naive CD4+ T cells" else "Effector CD4+ T cells"
label_Bo_c1 <- if (s_naive_c1 >= s_effector_c1) "Naive CD4+ T cells" else "Effector CD4+ T cells"

# Fallback: score quasi uguali → cluster più grande = Naive
if (label_Bo_c0 == label_Bo_c1) {
  n0 <- sum(Bo_I$seurat_clusters == "0")
  n1 <- sum(Bo_I$seurat_clusters == "1")
  label_Bo_c0 <- if (n0 >= n1) "Naive CD4+ T cells" else "Effector CD4+ T cells"
  label_Bo_c1 <- if (n0 >= n1) "Effector CD4+ T cells" else "Naive CD4+ T cells"
  cat("  [FALLBACK] Score simili – cluster più grande = Naive\n")
}
cat(sprintf("  C0 → %s (naive=%.3f, eff=%.3f)\n",
            label_Bo_c0, s_naive_c0, s_effector_c0))
cat(sprintf("  C1 → %s (naive=%.3f, eff=%.3f)\n",
            label_Bo_c1, s_naive_c1, s_effector_c1))

# Plot diagnostico D1
key_d1 <- c("CD4","CD8A","CCR7","SELL","IL7R","GZMK","CD44","S100A4")
key_d1 <- key_d1[key_d1 %in% rownames(Bo_I)]
p_d1 <- VlnPlot(Bo_I, features = key_d1, idents = c("0","1"),
                 ncol = 4, pt.size = 0,
                 cols = c("#1565C0","#C62828")) &
  theme_classic(base_size = 10) &
  theme(axis.title.x = element_blank(),
        plot.title = element_text(face = "bold", size = 10))
save_plot(p_d1 + plot_annotation(
  title = "D1 | Bo_I – C0 (blu) vs C1 (rosso): Naive vs Effector CD4+"),
  paste0(out_doubt, "D1_BoI_NaiveVsEffector.png"), w = 16, h = 6)

# ── D2 – Ca_I C5: NK o Tregs? ────────────────────────────────
cat("\n[D2] Ca_I C5 – NK cells o Tregs?\n")
Idents(Ca_I) <- "seurat_clusters"

nk_genes   <- c("NKG7","NCAM1","NCR1","KLRD1","KLRB1","GNLY",
                "FCGR3A","TYROBP","S1PR5","CX3CR1","GZMB","PRF1")
treg_genes <- c("FOXP3","IL2RA","CTLA4","TIGIT","IKZF2",
                "TNFRSF18","TNFRSF4","ENTPD1","CCR8","RTKN2")

Ca_I <- AddModuleScore(Ca_I, features = list(nk_genes[nk_genes %in% rownames(Ca_I)]),
                       name = "NK_score")
Ca_I <- AddModuleScore(Ca_I, features = list(treg_genes[treg_genes %in% rownames(Ca_I)]),
                       name = "Treg_score")
Ca_I$NK_score   <- Ca_I$NK_score1;   Ca_I$NK_score1   <- NULL
Ca_I$Treg_score <- Ca_I$Treg_score1; Ca_I$Treg_score1 <- NULL

score_d2 <- Ca_I@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(NK   = round(mean(NK_score), 3),
            Treg = round(mean(Treg_score), 3),
            n    = n(), .groups = "drop")
cat("\n  NK vs Treg score per cluster:\n")
print(as.data.frame(score_d2))

c5_nk   <- score_d2$NK[score_d2$seurat_clusters   == "5"]
c5_treg <- score_d2$Treg[score_d2$seurat_clusters == "5"]
label_Ca_c5 <- if (length(c5_nk) > 0 && c5_nk > c5_treg) "NK cells" else "Tregs"
cat(sprintf("  C5 → %s (NK=%.3f, Treg=%.3f)\n",
            label_Ca_c5, c5_nk, c5_treg))

p_d2 <- VlnPlot(Ca_I, features = c("NK_score","Treg_score"),
                 idents = c("4","5"), ncol = 2, pt.size = 0.2,
                 cols = c("#F39C12","#16A085")) &
  theme_classic(base_size = 10) &
  theme(axis.title.x = element_blank())
save_plot(p_d2 + plot_annotation(
  title = "D2 | Ca_I – C5: NK cells vs Tregs?"),
  paste0(out_doubt, "D2_CaI_C5_NKvsTreg.png"), w = 10, h = 5)

# Caratterizzazione CD4 vs CD8 di Ca_I C0 (cluster proliferante).
# La funzione characterize_prolif() viene definita nella sezione D3
# → qui viene chiamata dopo D3; Ca_I C0 è aggiornato lì sotto.
# (label_Ca_c0_prolif viene assegnata subito dopo D3)

# ── D3 – Me_I C1 vs C3: fasi del ciclo cellulare ────────────
cat("\n[D3] Me_I C1 vs C3 – fasi ciclo cellulare\n")
Idents(Me_I) <- "seurat_clusters"

s_ok   <- cc.genes$s.genes[cc.genes$s.genes %in% rownames(Me_I)]
g2m_ok <- cc.genes$g2m.genes[cc.genes$g2m.genes %in% rownames(Me_I)]
Me_I   <- CellCycleScoring(Me_I, s.features = s_ok,
                           g2m.features = g2m_ok, set.ident = FALSE)

phase_d3 <- Me_I@meta.data %>%
  filter(seurat_clusters %in% c("1","3")) %>%
  group_by(seurat_clusters, Phase) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(seurat_clusters) %>%
  mutate(pct = round(100 * n / sum(n), 1))
print(as.data.frame(phase_d3))

c1_dom <- phase_d3 %>% filter(seurat_clusters == "1") %>%
  slice_max(pct, n = 1) %>% pull(Phase)
c3_dom <- phase_d3 %>% filter(seurat_clusters == "3") %>%
  slice_max(pct, n = 1) %>% pull(Phase)

# ── Funzione di caratterizzazione CD4 vs CD8 per cluster
#    proliferanti. Usata per Me_I C1/C3 e Ca_I C0.
#    Logica:
#      CD4 > CD8 → "Proliferating CD4+ T cells"
#      CD8 > CD4 → "Proliferating CD8+ T cells"
#      Ambiguo (differenza < soglia) → "Proliferating T cells"
#    Soglia: 20% di differenza relativa per essere sicuri;
#    al di sotto del 20% il segnale è troppo simile per concludere.
characterize_prolif <- function(obj, cluster_id) {
  Idents(obj) <- "seurat_clusters"

  lineage_genes <- c("CD4","CD8A","CD8B")
  lineage_genes <- lineage_genes[lineage_genes %in% rownames(obj)]

  avg <- fix_colnames(
    AverageExpression(obj, features = lineage_genes,
                      group.by = "seurat_clusters",
                      assay = "RNA", slot = "data")$RNA)

  cl <- as.character(cluster_id)
  if (!cl %in% colnames(avg)) {
    cat(sprintf("  [WARN] Cluster %s non trovato in avg_expr.\n", cl))
    return("Proliferating T cells")
  }

  avg_cd4 <- if ("CD4"  %in% rownames(avg)) avg["CD4",  cl] else 0
  avg_cd8 <- max(
    if ("CD8A" %in% rownames(avg)) avg["CD8A", cl] else 0,
    if ("CD8B" %in% rownames(avg)) avg["CD8B", cl] else 0
  )

  cat(sprintf("  C%s | CD4=%.4f | CD8=%.4f", cl, avg_cd4, avg_cd8))

  denom <- max(avg_cd4, avg_cd8, 1e-6)
  delta <- abs(avg_cd4 - avg_cd8) / denom   # differenza relativa

  label <- if (delta < 0.20) {
    cat(" → AMBIGUO (delta<20%) → Proliferating T cells\n")
    "Proliferating T cells"
  } else if (avg_cd4 > avg_cd8) {
    cat(" → CD4 dominante → Proliferating CD4+ T cells\n")
    "Proliferating CD4+ T cells"
  } else {
    cat(" → CD8 dominante → Proliferating CD8+ T cells\n")
    "Proliferating CD8+ T cells"
  }
  label
}

cat("\n[D3] Caratterizzazione CD4 vs CD8 cluster proliferanti in Me_I:\n")
label_Me_c1 <- characterize_prolif(Me_I, "1")
label_Me_c3 <- characterize_prolif(Me_I, "3")

# Se le fasi dominanti del ciclo differiscono, aggiunge il suffisso
# (S)/(G2M) al label per distinguere ulteriormente i due cluster
if (!identical(c1_dom, c3_dom)) {
  label_Me_c1 <- sub("T cells$",
                     paste0("T cells (", c1_dom, ")"), label_Me_c1)
  label_Me_c3 <- sub("T cells$",
                     paste0("T cells (", c3_dom, ")"), label_Me_c3)
}
cat(sprintf("  C1 → %s | C3 → %s\n", label_Me_c1, label_Me_c3))

# Caratterizzazione Ca_I C0 (cluster proliferante)
cat("\n[D2+] Caratterizzazione CD4 vs CD8 cluster proliferante in Ca_I:\n")
label_Ca_c0 <- characterize_prolif(Ca_I, "0")

lineage_ok_ca <- c("CD4","CD8A","CD8B")[c("CD4","CD8A","CD8B") %in% rownames(Ca_I)]
p_prolif_ca <- VlnPlot(Ca_I, features = lineage_ok_ca,
                        idents = "0", ncol = 3,
                        pt.size = 0.3,
                        cols = "#457B9D") &
  theme_classic(base_size = 10) &
  theme(axis.title.x = element_blank(),
        plot.title = element_text(face = "bold", size = 10))
save_plot(
  p_prolif_ca + plot_annotation(
    title    = "Ca_I – C0: CD4 vs CD8 (cluster proliferante)",
    theme    = theme(plot.title = element_text(face = "bold",
                                               size = 12, hjust = 0.5))),
  paste0(out_doubt, "CaI_C0_prolif_CD4vsCD8.png"), w = 9, h = 5)

# Plot diagnostico: VlnPlot CD4/CD8A/CD8B sui cluster proliferanti
lineage_ok <- c("CD4","CD8A","CD8B")[c("CD4","CD8A","CD8B") %in% rownames(Me_I)]
p_prolif_me <- VlnPlot(Me_I, features = lineage_ok,
                        idents = c("1","3"), ncol = 3,
                        pt.size = 0.3,
                        cols = c("#3498DB","#E74C3C")) &
  theme_classic(base_size = 10) &
  theme(axis.title.x = element_blank(),
        plot.title = element_text(face = "bold", size = 10))
save_plot(
  p_prolif_me + plot_annotation(
    title    = "D3 | Me_I – C1 (blu) vs C3 (rosso): CD4 vs CD8",
    subtitle = "Caratterizzazione lineage cluster proliferanti",
    theme    = theme(plot.title = element_text(face = "bold",
                                               size = 12, hjust = 0.5))),
  paste0(out_doubt, "D3_MeI_prolif_CD4vsCD8.png"), w = 9, h = 5)

# ============================================================
# STEP 3 – RICERCA CLUSTER MIELOIDI IN TUTTI I CAMPIONI I
#
# Approccio:
#  a) Calcola myeloid_score su tutti i cluster con AddModuleScore
#  b) Stampa tabella ordinata per myeloid_score decrescente
#  c) Criteri di assegnazione:
#       CONFERMATO: myeloid_score > 0.35 E cd3_mean < 0.3
#       PROPOSTA:   myeloid_score > 0.20 E cd3_mean < 0.5
#                   (non riassegnato in automatico, segnalato)
#  d) Genera FeaturePlot + VlnPlot + MyeloidScore vs Cluster
#     in Pipeline_I/Myeloid_search/
#
# Per Me_I: se nessun cluster è CONFERMATO a res=0.5, tenta
# un reclustering a res=1.0 prima di dichiarare assenza.
# ============================================================
section("STEP 3 | Ricerca cluster mieloidi in tutti i campioni I")

# Marcatori diagnostici da visualizzare
MYELOID_FEAT <- c("LYZ","S100A8","S100A9","VCAN","CD14",
                  "FCGR3A","HLA-DRA","CSF1R","CD68","ITGAM",
                  "CLEC9A","CD1C","FCER1A")

search_myeloid <- function(obj, sample_name, out_dir,
                           try_hires = FALSE) {
  cat(paste0("\n", strrep("-", 55), "\n",
             "  Ricerca mieloidi: ", sample_name, "\n",
             strrep("-", 55), "\n"))
  cat("Annotazione corrente:\n")
  print(table(obj$seurat_clusters))

  Idents(obj) <- "seurat_clusters"

  # ── Module scores ─────────────────────────────────────────
  score_sigs <- list(myeloid = SIG$myeloid,
                     mono14  = SIG$mono14 %||%
                       c("CD14","LYZ","S100A8","S100A9","FCN1",
                         "VCAN","CXCL8","THBS1","PLBD1","CLEC7A"),
                     mono16  = SIG$mono16 %||%
                       c("FCGR3A","MS4A7","LILRB2","CDKN1C",
                         "CX3CR1","CSF1R","HMOX1"))

  for (sn in names(score_sigs)) {
    g <- score_sigs[[sn]][score_sigs[[sn]] %in% rownames(obj)]
    if (length(g) >= 3) {
      obj <- AddModuleScore(obj, features = list(g),
                            name = paste0(sn, "_Ms"))
      obj[[paste0(sn, "_Ms")]] <- obj[[paste0(sn, "_Ms1")]]
      obj[[paste0(sn, "_Ms1")]] <- NULL
    }
  }

  # CD3 medio da AverageExpression
  avg_cd3 <- fix_colnames(
    AverageExpression(obj,
                      features = SIG$cd3[SIG$cd3 %in% rownames(obj)],
                      group.by = "seurat_clusters",
                      assay = "RNA", slot = "data")$RNA)

  score_cols <- paste0(names(score_sigs), "_Ms")
  score_cols <- score_cols[score_cols %in% colnames(obj@meta.data)]

  tbl <- obj@meta.data %>%
    group_by(seurat_clusters) %>%
    summarise(across(all_of(score_cols),
                     ~ round(mean(.x, na.rm = TRUE), 4)),
              n_cells = n(), .groups = "drop") %>%
    mutate(seurat_clusters = as.character(as.integer(as.character(seurat_clusters))),
           cd3_mean = round(sapply(seurat_clusters, function(cl) {
             v <- avg_cd3[, colnames(avg_cd3) == cl, drop = FALSE]
             if (ncol(v) == 0) NA_real_ else mean(v, na.rm = TRUE)
           }), 4))

  if ("myeloid_Ms" %in% colnames(tbl))
    tbl <- tbl %>% arrange(desc(myeloid_Ms))

  cat("\nScore per cluster (myeloid_score decrescente):\n")
  print(as.data.frame(tbl), row.names = FALSE)

  # ── Classificazione ───────────────────────────────────────
  confirmed <- character(0)
  proposed  <- character(0)

  if ("myeloid_Ms" %in% colnames(tbl)) {
    cd3_v <- if ("cd3_mean" %in% colnames(tbl))
               replace(tbl$cd3_mean, is.na(tbl$cd3_mean), 0) else
               rep(0, nrow(tbl))

    for (i in seq_len(nrow(tbl))) {
      ms  <- tbl$myeloid_Ms[i]
      cd3 <- cd3_v[i]
      cl  <- tbl$seurat_clusters[i]

      if (ms > 0.35 && cd3 < 0.3) {
        confirmed <- c(confirmed, cl)
        cat(sprintf("  [CONFERMATO] C%s [%d celle] myeloid=%.3f cd3=%.3f\n",
                    cl, tbl$n_cells[i], ms, cd3))
      } else if (ms > 0.20 && cd3 < 0.5) {
        proposed <- c(proposed, cl)
        cat(sprintf("  [PROPOSTA  ] C%s [%d celle] myeloid=%.3f cd3=%.3f\n",
                    cl, tbl$n_cells[i], ms, cd3))
      }
    }
  }

  # ── Reclustering a res=1.0 se nessun CONFERMATO ──────────
  # (solo per campioni specificati con try_hires=TRUE)
  if (length(confirmed) == 0 && try_hires) {
    cat("\n  Nessun cluster CONFERMATO a res corrente.\n")
    cat("  Provo reclustering a res=1.0...\n")

    obj_hi <- FindNeighbors(obj, dims = 1:30, verbose = FALSE)
    obj_hi <- FindClusters(obj_hi, resolution = 1.0, verbose = FALSE)
    cat(paste0("  Nuovi cluster a res=1.0: ",
               length(unique(obj_hi$seurat_clusters)), "\n"))

    g_mye <- SIG$myeloid[SIG$myeloid %in% rownames(obj_hi)]
    if (length(g_mye) >= 3) {
      obj_hi <- AddModuleScore(obj_hi, features = list(g_mye),
                               name = "myeloid_Ms")
      obj_hi$myeloid_Ms <- obj_hi$myeloid_Ms1
      obj_hi$myeloid_Ms1 <- NULL
    }

    avg_cd3_hi <- fix_colnames(
      AverageExpression(obj_hi,
                        features = SIG$cd3[SIG$cd3 %in% rownames(obj_hi)],
                        group.by = "seurat_clusters",
                        assay = "RNA", slot = "data")$RNA)

    tbl_hi <- obj_hi@meta.data %>%
      group_by(seurat_clusters) %>%
      summarise(myeloid_Ms = round(mean(myeloid_Ms, na.rm = TRUE), 4),
                n_cells    = n(), .groups = "drop") %>%
      mutate(seurat_clusters = as.character(as.integer(as.character(seurat_clusters))),
             cd3_mean = round(sapply(seurat_clusters, function(cl) {
               v <- avg_cd3_hi[, colnames(avg_cd3_hi) == cl, drop = FALSE]
               if (ncol(v) == 0) NA_real_ else mean(v, na.rm = TRUE)
             }), 4)) %>%
      arrange(desc(myeloid_Ms))

    cat("\nScore a res=1.0:\n")
    print(as.data.frame(tbl_hi), row.names = FALSE)

    cd3_hi <- replace(tbl_hi$cd3_mean, is.na(tbl_hi$cd3_mean), 0)
    for (i in seq_len(nrow(tbl_hi))) {
      ms  <- tbl_hi$myeloid_Ms[i]
      cd3 <- cd3_hi[i]
      cl  <- tbl_hi$seurat_clusters[i]
      # Soglia leggermente più lasca al reclustering
      if (ms > 0.30 && cd3 < 0.35 &&
          tbl_hi$n_cells[i] / ncol(obj_hi) < 0.10) {
        confirmed <- c(confirmed, cl)
        cat(sprintf("  [CONFERMATO hires] C%s [%d celle] myeloid=%.3f\n",
                    cl, tbl_hi$n_cells[i], ms))
      }
    }

    if (length(confirmed) > 0) {
      # Trasferisci annotazione originale e rinomina cluster confermati
      obj_orig_ct <- obj$seurat_clusters
      obj_hi$cell_type_orig <- names(obj_orig_ct)[
        match(colnames(obj_hi), colnames(obj))]

      cl_map_hi <- obj_hi@meta.data %>%
        group_by(seurat_clusters) %>%
        count(cell_type_orig) %>%
        slice_max(n, n = 1) %>% ungroup() %>%
        mutate(label = if_else(
          as.character(seurat_clusters) %in% confirmed,
          "Myeloid cells",
          as.character(cell_type_orig)))

      cat("\n  Mappa cluster hires → label:\n")
      print(as.data.frame(cl_map_hi[, c("seurat_clusters","label")]))

      ann_map <- setNames(cl_map_hi$label,
                          as.character(cl_map_hi$seurat_clusters))
      new_ct  <- unname(ann_map[as.character(obj_hi$seurat_clusters)])
      names(new_ct) <- colnames(obj_hi)
      obj_hi <- AddMetaData(obj_hi, new_ct, "cell_type")
      # Restituisce l'oggetto hires annotato
      obj <- obj_hi
    } else {
      cat("  -> Nessun cluster mieloide trovato neanche a res=1.0.\n")
      cat("     Segnale mieloide diffuso o assente in", sample_name, "\n")
    }
  }

  if (length(confirmed) == 0 && length(proposed) == 0) {
    cat("  -> Nessun segnale mieloide rilevante in", sample_name, "\n")
    cat("     (campione arricchito in T cells; mieloidi assenti o <1%)\n")
  }

  if (length(confirmed) > 0 && !"cell_type" %in% colnames(obj@meta.data)) {
    # Rinomina su clustering originale
    curr  <- if ("cell_type" %in% colnames(obj@meta.data))
               as.character(obj$cell_type) else
               as.character(obj$seurat_clusters)
    cl_id <- as.character(obj$seurat_clusters)
    new_ct <- ifelse(cl_id %in% confirmed, "Myeloid cells", curr)
    names(new_ct) <- colnames(obj)
    obj <- AddMetaData(obj, new_ct, "cell_type")
    cat("  Annotazione cell_type aggiornata con Myeloid cells.\n")
  }

  # ── Plot diagnostici ──────────────────────────────────────
  feat_ok <- MYELOID_FEAT[MYELOID_FEAT %in% rownames(obj)]

  if (length(feat_ok) >= 4) {
    p_fp <- FeaturePlot(obj, features = feat_ok, reduction = "umap",
                        ncol = 4, min.cutoff = "q05",
                        max.cutoff = "q95", order = TRUE,
                        pt.size = 0.4) &
      theme_classic(base_size = 9) &
      theme(axis.text = element_blank(),
            axis.ticks = element_blank(),
            plot.title = element_text(size = 9))
    save_plot(p_fp,
              paste0(out_dir, sample_name, "_myeloid_FeaturePlot.png"),
              w = 16, h = ceiling(length(feat_ok) / 4) * 4)
  }

  if (length(feat_ok) >= 2) {
    Idents(obj) <- "seurat_clusters"
    p_vln <- VlnPlot(obj, features = feat_ok, ncol = 4,
                     pt.size = 0, fill.by = "ident") &
      theme(axis.text.x = element_text(size = 7),
            axis.title   = element_blank(),
            plot.title   = element_text(size = 9))
    save_plot(p_vln,
              paste0(out_dir, sample_name, "_myeloid_VlnPlot.png"),
              w = 16, h = ceiling(length(feat_ok) / 4) * 3)
    Idents(obj) <- if ("cell_type" %in% colnames(obj@meta.data))
                     "cell_type" else "seurat_clusters"
  }

  if ("myeloid_Ms" %in% colnames(obj@meta.data)) {
    Idents(obj) <- "seurat_clusters"
    p_sc <- FeaturePlot(obj, features = "myeloid_Ms",
                        reduction = "umap", pt.size = 0.5,
                        order = TRUE) +
      scale_color_gradientn(
        colors = c("lightgrey","#FFF176","#FB8C00","#B71C1C"),
        name = "Myeloid\nscore") +
      ggtitle(paste0(sample_name, " – Myeloid score")) +
      theme_classic(base_size = 11) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))

    p_cl <- DimPlot(obj, reduction = "umap", label = TRUE,
                    label.size = 3.5, repel = TRUE, pt.size = 0.5) +
      ggtitle(paste0(sample_name, " – Clusters")) +
      theme_classic(base_size = 11) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            legend.position = "none")

    save_plot(p_sc | p_cl,
              paste0(out_dir, sample_name, "_myeloid_score_vs_clusters.png"),
              w = 14, h = 6)
    Idents(obj) <- if ("cell_type" %in% colnames(obj@meta.data))
                     "cell_type" else "seurat_clusters"
  }

  if (length(proposed) > 0) {
    cat(paste0(
      "\n  [ATTENZIONE] Cluster PROPOSTA: ",
      paste(proposed, collapse = ", "), "\n",
      "  Segnale mieloide borderline – verificare i plot in:\n",
      "  ", out_dir, "\n",
      "  Non riassegnati automaticamente.\n"))
  }

  return(list(obj = obj, confirmed = confirmed, proposed = proposed))
}

# Operatore %||% (null-coalesce) per i sig opzionali
`%||%` <- function(a, b) if (!is.null(a)) a else b

res_mye_Bo <- search_myeloid(Bo_I, "Bo_I", out_mye, try_hires = FALSE)
res_mye_Ca <- search_myeloid(Ca_I, "Ca_I", out_mye, try_hires = TRUE)
res_mye_Me <- search_myeloid(Me_I, "Me_I", out_mye, try_hires = TRUE)

# Recupera oggetti eventualmente aggiornati con cluster mieloidi
Bo_I <- res_mye_Bo$obj
Ca_I <- res_mye_Ca$obj
Me_I <- res_mye_Me$obj

# ============================================================
# STEP 4 – ANNOTAZIONE FINALE
# ============================================================
section("STEP 4 | Annotazione finale campioni I")

# Funzione che applica l'annotazione a un oggetto Seurat.
# Preserva eventuali Myeloid cells già assegnate da search_myeloid.
apply_annotation <- function(obj, annotation_map, sample_name) {
  Idents(obj) <- "seurat_clusters"
  cl_ids  <- as.character(Idents(obj))
  missing <- setdiff(unique(cl_ids), names(annotation_map))
  if (length(missing) > 0) {
    cat(sprintf("[WARN] %s: cluster senza mappa (verifica): %s\n",
                sample_name, paste(missing, collapse = ", ")))
    # Assegna fallback ai cluster mancanti
    for (m in missing) annotation_map[[m]] <- "Unknown"
  }

  # Se search_myeloid ha già assegnato cell_type per alcuni cluster
  # (quelli confermati come Myeloid cells), preserva quella assegnazione
  existing_ct <- if ("cell_type" %in% colnames(obj@meta.data))
                   as.character(obj$cell_type) else rep(NA, ncol(obj))
  myeloid_barcodes <- colnames(obj)[!is.na(existing_ct) &
                                      existing_ct == "Myeloid cells"]

  labels <- unname(annotation_map[cl_ids])
  names(labels) <- colnames(obj)

  # Ripristina Myeloid cells dove erano già state assegnate
  if (length(myeloid_barcodes) > 0)
    labels[myeloid_barcodes] <- "Myeloid cells"

  obj <- AddMetaData(obj, metadata = labels, col.name = "cell_type")
  Idents(obj) <- "cell_type"
  cat(paste0("[", sample_name, "] Annotazione:\n"))
  print(table(obj$cell_type))
  return(obj)
}

# ── Dizionari annotazione ────────────────────────────────────
# Bo_I: C4 → Myeloid cells (ex DC, confermato in D4 originale)
annotation_Bo <- c(
  "0" = label_Bo_c0,
  "1" = label_Bo_c1,
  "2" = "Cytotoxic CD8+ T cells",
  "3" = "Memory T cells",
  "4" = "Myeloid cells"
)

# Ca_I: C0 = cluster proliferante caratterizzato (CD4 o CD8)
#        C5 da D2 automatico (NK o Tregs)
annotation_Ca <- c(
  "0" = label_Ca_c0,
  "1" = "Naive CD4+ T cells",
  "2" = "Effector CD4+ T cells",
  "3" = "Cytotoxic CD8+ T cells",
  "4" = "Tregs",
  "5" = label_Ca_c5
)

# Me_I: C4 = secondo cluster proliferante → caratterizzato come C1/C3
# (viene calcolato qui perché characterize_prolif è già disponibile)
cat("\n[D3+] Caratterizzazione CD4 vs CD8 cluster proliferante in Me_I C4:\n")
label_Me_c4 <- characterize_prolif(Me_I, "4")

# Me_I: 5 cluster base (se search_myeloid ha trovato cluster
# mieloidi a res=1.0 l'oggetto Me_I è già stato aggiornato)
annotation_Me <- c(
  "0" = "Cytotoxic CD8+ T cells",
  "1" = label_Me_c1,
  "2" = "Tregs",
  "3" = label_Me_c3,
  "4" = label_Me_c4
)

# Se Me_I è già stato riclusterizzato a res=1.0, i cluster
# sono già annotati in cell_type: salta l'applicazione manuale
if (length(res_mye_Me$confirmed) > 0 &&
    "cell_type" %in% colnames(Me_I@meta.data) &&
    length(unique(Me_I$seurat_clusters)) > 5) {
  cat("[Me_I] Usa annotazione da reclustering hires (già in cell_type).\n")
  Me_I_ann <- Me_I
  Idents(Me_I_ann) <- "cell_type"
} else {
  Me_I_ann <- apply_annotation(Me_I, annotation_Me, "Me_I")
}

Bo_I_ann <- apply_annotation(Bo_I, annotation_Bo, "Bo_I")
Ca_I_ann <- apply_annotation(Ca_I, annotation_Ca, "Ca_I")

# ── Th subtypes: non applicabili nei campioni I (pochi cluster)
# Nei campioni I i cluster CD4 sono già distinti come Naive/Effector.
# I Th subtypes verranno cercati negli AB dove ci sono più cluster.

# ============================================================
# STEP 5 – UMAP CON LEGENDA CONSISTENTE
# ============================================================
section("STEP 5 | UMAP finali campioni I")

plot_umap_consistent(Bo_I_ann, "Bo_I",
                     paste0(out_umap, "Bo_I_UMAP_annotated.png"))
plot_umap_consistent(Ca_I_ann, "Ca_I",
                     paste0(out_umap, "Ca_I_UMAP_annotated.png"))
plot_umap_consistent(Me_I_ann, "Me_I",
                     paste0(out_umap, "Me_I_UMAP_annotated.png"))

# Pannello combinato I (tutti e tre in verticale, legenda globale)
all_present_I <- unique(c(
  as.character(Bo_I_ann$cell_type),
  as.character(Ca_I_ann$cell_type),
  as.character(Me_I_ann$cell_type)
))

make_dp_mini <- function(obj, nm) {
  cols <- get_colors(sort(unique(as.character(obj$cell_type))))
  DimPlot(obj, reduction = "umap", group.by = "cell_type",
          label = TRUE, label.size = 3, repel = TRUE,
          cols = cols, pt.size = 0.5) +
    ggtitle(nm) + theme_classic(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "none")
}

p_I_all <- (make_dp_mini(Bo_I_ann,"Bo_I") /
            make_dp_mini(Ca_I_ann,"Ca_I") /
            make_dp_mini(Me_I_ann,"Me_I") |
            make_full_legend(all_present_I)) +
  plot_layout(widths = c(10, 2.5)) +
  plot_annotation(title = "Campioni I – panoramica",
                  theme = theme(plot.title =
                                  element_text(face = "bold",
                                               hjust = 0.5,
                                               size = 14)))

save_plot(p_I_all,
          paste0(out_umap, "ALL_I_samples_UMAP.png"),
          w = 14, h = 21)

# ============================================================
# STEP 6 – SALVATAGGIO
# ============================================================
section("STEP 6 | Salvataggio")

annotated_I <- list(
  Bo_bone_I = Bo_I_ann,
  Ca_bone_I = Ca_I_ann,
  Me_bone_I = Me_I_ann
)

rds_out <- paste0(base_dir, "all_I_samples_annotated.rds")
saveRDS(annotated_I, rds_out)
cat(paste0("\n>> RDS salvato: ", rds_out, "\n"))

cat(paste0(
  "\n", strrep("=", 65), "\n",
  "  PIPELINE 1 COMPLETATA\n\n",
  "  Campioni I annotati: Bo_bone_I | Ca_bone_I | Me_bone_I\n",
  "  RDS: all_I_samples_annotated.rds\n",
  "  Dubbi risolti:\n",
  "    D1 – Bo_I C0 → ", label_Bo_c0, "\n",
  "         Bo_I C1 → ", label_Bo_c1, "\n",
  "    D2 – Ca_I C5 → ", label_Ca_c5, "\n",
  "    D3 – Me_I C1 → ", label_Me_c1, "\n",
  "         Me_I C3 → ", label_Me_c3, "\n",
  "  Mieloidi Bo_I C4 → Myeloid cells (confermato)\n",
  "  Ricerca mieloidi Ca_I e Me_I:\n",
  "    Ca_I confermati: ", paste(res_mye_Ca$confirmed, collapse=", "),
  if (length(res_mye_Ca$confirmed)==0) "nessuno" else "", "\n",
  "    Me_I confermati: ", paste(res_mye_Me$confirmed, collapse=", "),
  if (length(res_mye_Me$confirmed)==0) "nessuno" else "", "\n",
  "  Prossimo step: PIPELINE_2_annotate_AB_samples.R\n",
  strrep("=", 65), "\n"
))
