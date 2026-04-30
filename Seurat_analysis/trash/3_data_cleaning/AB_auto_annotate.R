# ============================================================
#  NEW_AB – Annotazione automatica campioni AB
#
#  Approccio: replica della logica di 2_resolve_doubts.R
#  (AverageExpression + AddModuleScore) applicata a tutti i
#  cluster di tutti i campioni AB.
#
#  POPOLAZIONI CONSIDERATE:
#  ── Già viste nei campioni I (colori FISSI) ───────────────
#    Naive CD4+ T cells, Effector CD4+ T cells, Memory T cells,
#    Cytotoxic CD8+ T cells, Proliferating T cells,
#    Proliferating CD8+ T cells, Tregs, NK cells,
#    Dendritic Cells
#  ── Nuove (possibili negli AB, colori dinamici) ───────────
#    Naive CD8+ T cells, NKT cells, gamma-delta T cells,
#    MAIT cells, ILC,
#    B cells, Memory B cells, Plasma cells,
#    CD14 Monocytes, CD16 Monocytes, Basophils,
#    HSPC, Erythroid cells, Platelets
#
#  GERARCHIA DECISIONALE:
#   1.  Platelets       (PPBP, PF4 – CD3/CD19 assenti)
#   2.  Erythroid       (HBB, HBA1 – CD3/CD19 assenti)
#   3.  HSPC            (CD34 – CD3/CD19/CD56 assenti)
#   4.  Plasma cells    (JCHAIN, MZB1 – CD3 assente)
#   5.  B cells         (CD19, MS4A1 – CD3 assente)
#   6.  Basophils       (CPA3, TPSAB1 – CD3/CD19 assenti)
#   7.  CD14 Monocytes  (CD14 alto, DC score basso)
#   8.  CD16 Monocytes  (FCGR3A alto, CD14 basso, DC score basso)
#   9.  Dendritic Cells (DC score > 0.35, CD3 < 0.3, effector < 0.4)
#   10. NKT cells       (CD3+ e CD56+ insieme)
#   11. gamma-delta T   (TRDC/TRGC – CD3+, CD4-/CD8-)
#   12. MAIT cells      (SLC4A10, KLRB1 – CD3+)
#   13. ILC             (GATA3/RORC – CD3- e CD56-)
#   14. NK cells        (NCAM1/KLRD1 – CD3 basso, CD4/CD8 bassi)
#   15. Tregs           (FOXP3/IL2RA – CD4+)
#   16. Proliferating   (MKI67/TOP2A alti)
#   17. CD8+ T cells    (CD8A > CD4)
#   18. CD4+ T cells    (CD4 alto)
#   19. Fallback        (score T cell massimo)
#
#  Output:
#    base_dir/AB_annotation/
#      <sample>_annotation_decisions.xlsx
#      <sample>_UMAP_annotated.png
#    ALL_AB_samples_UMAP.png
#    all_AB_samples_annotated.rds
#    all_samples_annotated_COMPLETE.rds  (I + AB unificati)
# ============================================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(openxlsx)
library(scales)

# ── UNICO PUNTO DA MODIFICARE ────────────────────────────────
base_dir <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/CART/Code/Seurat_analysis/3_data_cleaning/"
# ─────────────────────────────────────────────────────────────

out_dir <- paste0(base_dir, "AB_annotation/")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

section <- function(title) {
  cat(paste0("\n", strrep("=", 65), "\n  ", title, "\n", strrep("=", 65), "\n"))
}

# ── Caricamento seurat_list ───────────────────────────────────
if (!exists("seurat_list")) {
  seurat_list_path <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/CART/Data/seurat_obj_list/seurat_samples_sctype_azimuth_pbmc_bonemarrow_clonalvdj_CAR.rds"
  cat("Caricamento seurat_list...\n")
  seurat_list <- readRDS(seurat_list_path)
}

ab_samples <- list(
  Ca_blood_AB = seurat_list$Ca_samples_blood$AB,
  Ca_bone_AB  = seurat_list$Ca_samples_bone$AB,
  Bo_blood_AB = seurat_list$Bo_samples_blood$AB,
  Bo_bone_AB  = seurat_list$Bo_samples_bone$AB,
  Me_bone_AB  = seurat_list$Me_samples_bone$AB
)

cat("\nCampioni AB da annotare:\n")
for (nm in names(ab_samples)) {
  cat(sprintf("  %-15s %d cellule | %d cluster\n",
              nm, ncol(ab_samples[[nm]]),
              length(unique(ab_samples[[nm]]$seurat_clusters))))
}

# ============================================================
# HELPERS
# ============================================================

join_layers_safe <- function(obj, name) {
  split_counts <- grep("^counts\\.", Layers(obj), value = TRUE)
  if (length(split_counts) > 0) {
    obj <- JoinLayers(obj)
    cat(paste0("[", name, "] JoinLayers applicato\n"))
  }
  Idents(obj) <- "seurat_clusters"
  return(obj)
}

mean_genes <- function(mat, genes, cluster) {
  genes_ok <- genes[genes %in% rownames(mat)]
  if (length(genes_ok) == 0 || !cluster %in% colnames(mat)) return(0)
  mean(mat[genes_ok, cluster])
}

get_score <- function(score_summary, cluster_id, score_name) {
  col <- paste0(score_name, "_score")
  if (!col %in% colnames(score_summary)) return(-99)
  val <- score_summary[[col]][
    as.character(score_summary$seurat_clusters) == as.character(cluster_id)]
  if (length(val) == 0) return(-99)
  val[1]
}

# ============================================================
# FIRME GENICHE COMPLETE
# ============================================================

SIG <- list(

  # Lineage markers (usati direttamente da AverageExpression)
  cd3    = c("CD3D","CD3E","CD3G"),
  cd4    = c("CD4"),
  cd8    = c("CD8A","CD8B"),
  cd19   = c("CD19","MS4A1"),
  cd56   = c("NCAM1"),
  cd14   = c("CD14"),
  cd16   = c("FCGR3A"),
  cd34   = c("CD34"),

  # T cell states (stessi di 2_resolve_doubts.R)
  naive    = c("CCR7","SELL","IL7R","TCF7","LEF1","KLF2"),
  effector = c("GZMK","CD44","S100A4","LGALS1","DUSP2"),
  cytotox  = c("GZMB","NKG7","PRF1","GNLY","GZMA","FGFBP2"),
  treg     = c("FOXP3","IL2RA","CTLA4","TIGIT","IKZF2",
               "TNFRSF18","TNFRSF4","ENTPD1","CCR8","RTKN2"),
  prolif   = c("MKI67","TOP2A","STMN1","PCNA","CCNB1","BIRC5","UBE2C"),

  # NK / innate lymphoid
  nk       = c("NCAM1","NCR1","KLRD1","KLRB1","GNLY",
               "FCGR3A","TYROBP","S1PR5","CX3CR1","XCL1","XCL2"),
  ilc      = c("GATA3","RORC","ICOS","IL1R1","KIT",
               "AREG","IL13","HPGDS"),

  # NKT cells
  nkt      = c("NCAM1","CD3D","KLRB1","ZBTB16","NKG7","GNLY","GZMB"),

  # gamma-delta T cells
  gdt      = c("TRDC","TRGC1","TRGC2","TRDV1","TRDV2","TRDV3",
               "TRGV9","KLRC1","KLRC2"),

  # MAIT cells
  mait     = c("SLC4A10","KLRB1","NCR3","RORC","IL18RAP","CXCR6"),

  # B cells
  bcell    = c("MS4A1","CD19","CD79A","CD79B","PAX5",
               "BANK1","IGHM","IGHD","TCL1A"),

  # Memory B cells
  bmem     = c("CD27","AIM2","TNFRSF13B","FCRL4","FCRL5","CD80"),

  # Plasma cells
  plasma   = c("JCHAIN","MZB1","IGHG1","IGHG2","IGKC","IGLC2",
               "SDC1","CD38","PRDM1","XBP1","DERL3"),

  # Monocytes CD14+
  mono14   = c("CD14","LYZ","S100A8","S100A9","FCN1","VCAN",
               "CXCL8","THBS1","PLBD1","CLEC7A","CTSS"),

  # Monocytes CD16+
  mono16   = c("FCGR3A","MS4A7","LILRB2","CDKN1C","CX3CR1",
               "LYPD2","CSF1R","HMOX1"),

  # Dendritic cells (stessa lista D4 in 2_resolve_doubts.R)
  dc       = c("HLA-DRA","HLA-DPB1","HLA-DQA1",
               "LYZ","VCAN","S100A8","S100A9",
               "CLEC9A","XCR1","CD1C","FCER1A","LILRA4",
               "AQP9","FPR1","LILRB2","MS4A6A","ITGAX",
               "BATF3","CLEC10A","SIGLEC6"),

  # Basophils / Mast cells
  baso     = c("CPA3","TPSAB1","TPSB2","MS4A2","GATA2",
               "HDC","SLC18A2","KIT","FCER1A"),

  # HSPC
  hspc     = c("CD34","SPINK2","AVP","CRHBP","HOPX",
               "MLLT3","MECOM","PROM1"),

  # Erythroid
  erythro  = c("HBB","HBA1","HBA2","GYPA","GYPB",
               "ALAS2","SLC4A1","CA1","AHSP","KLF1"),

  # Platelets / Megakaryocytes
  platelet = c("PPBP","PF4","GP1BA","ITGA2B","ITGB3",
               "GNG11","TUBB1","TREML1","CMTM5")
)

# ============================================================
# PALETTE COLORI
#
# palette_fixed: colori IDENTICI a STEP_3c – non modificare mai
# palette_new:   colori pre-definiti per tipi nuovi negli AB
# full_palette:  unione, usato come <<- per aggiornamenti globali
# ============================================================

palette_fixed <- c(
  "Naive CD4+ T cells"         = "#E63946",
  "Effector CD4+ T cells"      = "#F4A261",
  "Memory T cells"             = "#2A9D8F",
  "Cytotoxic CD8+ T cells"     = "#264653",
  "Proliferating T cells"      = "#457B9D",
  "Proliferating CD8+ T cells" = "#6A0572",
  "Tregs"                      = "#E9C46A",
  "Dendritic Cells"            = "#E76F51",
  "NK cells"                   = "#43AA8B"
)

palette_new <- c(
  "Naive CD8+ T cells"  = "#577590",
  "NKT cells"           = "#277DA1",
  "gamma-delta T cells" = "#4D908E",
  "MAIT cells"          = "#F3722C",
  "ILC"                 = "#F9C74F",
  "B cells"             = "#90BE6D",
  "Memory B cells"      = "#52B788",
  "Plasma cells"        = "#C77DFF",
  "CD14 Monocytes"      = "#9D0208",
  "CD16 Monocytes"      = "#DC2F02",
  "Basophils"           = "#6A4C93",
  "HSPC"                = "#B5838D",
  "Erythroid cells"     = "#FFAFCC",
  "Platelets"           = "#CDB4DB"
)

full_palette <- c(palette_fixed, palette_new)

# Restituisce i colori per l'oggetto; se trova tipi non in palette
# genera colori automatici via hue_pal() e li aggiunge alla palette
# globale (<<-) per consistenza tra campioni
get_colors <- function(obj) {
  types   <- sort(unique(as.character(obj$cell_type)))
  cols    <- full_palette[types]
  missing <- types[is.na(cols)]
  if (length(missing) > 0) {
    cat("[INFO] Tipi non in palette, colori auto-generati:",
        paste(missing, collapse = ", "), "\n")
    extra <- setNames(hue_pal()(length(missing)), missing)
    cols[missing]  <- extra
    full_palette   <<- c(full_palette, extra)
  }
  return(cols)
}

# ============================================================
# FUNZIONE DI ANNOTAZIONE AUTOMATICA
# ============================================================

auto_annotate <- function(obj, sample_name) {

  section(paste0("AUTO-ANNOTAZIONE: ", sample_name))

  Idents(obj) <- "seurat_clusters"

  # Fix Seurat v5: as.character() prima di as.integer()
  # evita il bug dove as.numeric() su factor perde il cluster 0
  clusters_chr <- as.character(
    sort(as.integer(as.character(unique(obj$seurat_clusters)))))

  # ── 1. AverageExpression ───────────────────────────────────
  all_sig_genes <- unique(unlist(SIG))
  genes_present <- all_sig_genes[all_sig_genes %in% rownames(obj)]

  avg <- AverageExpression(
    obj, features = genes_present,
    group.by = "seurat_clusters", assay = "RNA", slot = "data"
  )$RNA

  colnames(avg) <- gsub("^g", "", colnames(avg))
  colnames(avg) <- gsub("^RNA_snn_res\\.[0-9.]+_", "", colnames(avg))

  # ── 2. AddModuleScore per tutte le firme ───────────────────
  score_names <- c("naive","effector","cytotox","treg","prolif",
                   "nk","ilc","nkt","gdt","mait",
                   "bcell","bmem","plasma",
                   "mono14","mono16","dc","baso","hspc","erythro","platelet")

  for (sig_name in score_names) {
    genes_ok <- SIG[[sig_name]][SIG[[sig_name]] %in% rownames(obj)]
    if (length(genes_ok) >= 3) {
      score_col <- paste0(sig_name, "_score")
      obj <- AddModuleScore(obj, features = list(genes_ok), name = score_col)
      obj[[score_col]]               <- obj[[paste0(score_col, "1")]]
      obj[[paste0(score_col, "1")]]  <- NULL
    }
  }

  score_cols <- paste0(score_names, "_score")
  score_cols <- score_cols[score_cols %in% colnames(obj@meta.data)]

  score_summary <- obj@meta.data %>%
    group_by(seurat_clusters) %>%
    summarise(across(all_of(score_cols),
                     ~ round(mean(.x, na.rm = TRUE), 4)),
              n_cells = n(), .groups = "drop") %>%
    mutate(seurat_clusters = as.character(seurat_clusters)) %>%
    arrange(as.integer(seurat_clusters))

  cat("\nScore medio per cluster:\n")
  print(as.data.frame(score_summary))

  gs <- function(cl, sname) get_score(score_summary, cl, sname)

  # ── 3. DECISIONE PER CLUSTER ───────────────────────────────
  decisions <- data.frame(
    cluster    = clusters_chr,
    label      = NA_character_,
    confidence = NA_character_,
    rationale  = NA_character_,
    stringsAsFactors = FALSE
  )

  cat(paste0("\n[", sample_name, "] Decisioni per cluster:\n"))

  for (cl in clusters_chr) {

    s_cd3  <- mean_genes(avg, SIG$cd3,  cl)
    s_cd4  <- mean_genes(avg, SIG$cd4,  cl)
    s_cd8  <- mean_genes(avg, SIG$cd8,  cl)
    s_cd19 <- mean_genes(avg, SIG$cd19, cl)
    s_cd56 <- mean_genes(avg, SIG$cd56, cl)
    s_cd14 <- mean_genes(avg, SIG$cd14, cl)
    s_cd16 <- mean_genes(avg, SIG$cd16, cl)
    s_cd34 <- mean_genes(avg, SIG$cd34, cl)
    s_mk67 <- mean_genes(avg, c("MKI67","TOP2A"), cl)

    s_naive    <- gs(cl, "naive")
    s_effector <- gs(cl, "effector")
    s_cytotox  <- gs(cl, "cytotox")
    s_treg     <- gs(cl, "treg")
    s_prolif   <- gs(cl, "prolif")
    s_nk       <- gs(cl, "nk")
    s_ilc      <- gs(cl, "ilc")
    s_nkt      <- gs(cl, "nkt")
    s_gdt      <- gs(cl, "gdt")
    s_mait     <- gs(cl, "mait")
    s_bcell    <- gs(cl, "bcell")
    s_bmem     <- gs(cl, "bmem")
    s_plasma   <- gs(cl, "plasma")
    s_mono14   <- gs(cl, "mono14")
    s_mono16   <- gs(cl, "mono16")
    s_dc       <- gs(cl, "dc")
    s_baso     <- gs(cl, "baso")
    s_hspc     <- gs(cl, "hspc")
    s_erythro  <- gs(cl, "erythro")
    s_platelet <- gs(cl, "platelet")

    label      <- NA_character_
    confidence <- "ALTA"
    rationale  <- ""

    # ─── GERARCHIA DECISIONALE ────────────────────────────────

    # 1. PLATELETS
    if (s_platelet > 0.30 && s_cd3 < 0.3 && s_cd19 < 0.3) {
      label     <- "Platelets"
      rationale <- sprintf("platelet_score=%.3f, CD3=%.3f", s_platelet, s_cd3)
      if (s_platelet < 0.50) confidence <- "MEDIA"

    # 2. ERYTHROID
    } else if (s_erythro > 0.30 && s_cd3 < 0.3 && s_cd19 < 0.3) {
      label     <- "Erythroid cells"
      rationale <- sprintf("erythro_score=%.3f, CD3=%.3f", s_erythro, s_cd3)
      if (s_erythro < 0.50) confidence <- "MEDIA"

    # 3. HSPC
    } else if (s_hspc > 0.15 && s_cd34 > 0.3 &&
               s_cd3 < 0.3 && s_cd19 < 0.3) {
      label     <- "HSPC"
      rationale <- sprintf("hspc_score=%.3f, CD34=%.3f", s_hspc, s_cd34)
      if (s_hspc < 0.25) confidence <- "MEDIA"

    # 4. PLASMA CELLS
    } else if (s_plasma > 0.20 && s_cd3 < 0.3) {
      label     <- "Plasma cells"
      rationale <- sprintf("plasma_score=%.3f, CD19=%.3f, CD3=%.3f",
                           s_plasma, s_cd19, s_cd3)
      if (s_plasma < 0.35) confidence <- "MEDIA"

    # 5. B CELLS
    } else if (s_cd19 > 0.5 && s_cd3 < 0.5) {
      if (s_bmem > s_bcell && s_bmem > 0.05) {
        label     <- "Memory B cells"
        rationale <- sprintf("CD19=%.3f, bmem_score=%.3f > bcell_score=%.3f",
                             s_cd19, s_bmem, s_bcell)
      } else {
        label     <- "B cells"
        rationale <- sprintf("CD19=%.3f, bcell_score=%.3f, CD3=%.3f",
                             s_cd19, s_bcell, s_cd3)
      }
      if (s_cd19 < 0.8) confidence <- "MEDIA"

    # 6. BASOPHILS
    } else if (s_baso > 0.15 && s_cd3 < 0.3 && s_cd19 < 0.3) {
      label     <- "Basophils"
      rationale <- sprintf("baso_score=%.3f, CD3=%.3f", s_baso, s_cd3)
      if (s_baso < 0.30) confidence <- "MEDIA"

    # 7. CD14 MONOCYTES
    } else if (s_mono14 > 0.15 && s_cd14 > 0.5 &&
               s_cd3 < 0.3 && s_cd56 < 0.3 && s_dc < 0.35) {
      label     <- "CD14 Monocytes"
      rationale <- sprintf("mono14_score=%.3f, CD14=%.3f, DC_score=%.3f",
                           s_mono14, s_cd14, s_dc)
      if (s_mono14 < 0.25) confidence <- "MEDIA"

    # 8. CD16 MONOCYTES
    } else if (s_mono16 > 0.10 && s_cd16 > 0.3 &&
               s_cd3 < 0.3 && s_cd14 < 0.5 && s_dc < 0.35) {
      label     <- "CD16 Monocytes"
      rationale <- sprintf("mono16_score=%.3f, CD16=%.3f, CD14=%.3f",
                           s_mono16, s_cd16, s_cd14)
      if (s_mono16 < 0.20) confidence <- "MEDIA"

    # 9. DENDRITIC CELLS
    # Criteri stringenti: DC_score > 0.35, CD3 < 0.3, effector < 0.4
    # (evita falsi positivi con T cells attivate che esprimono HLA-DRA)
    } else if (s_dc > 0.35 && s_cd3 < 0.3 && s_effector < 0.4) {
      label     <- "Dendritic Cells"
      rationale <- sprintf("DC_score=%.3f, CD3=%.3f, effector=%.3f",
                           s_dc, s_cd3, s_effector)
      if (s_dc < 0.50) confidence <- "MEDIA"

    # 10. NKT CELLS (CD3+ e CD56+ insieme)
    } else if (s_cd3 > 0.5 && s_cd56 > 0.3 && s_nkt > s_nk) {
      label     <- "NKT cells"
      rationale <- sprintf("CD3=%.3f, CD56=%.3f, nkt_score=%.3f",
                           s_cd3, s_cd56, s_nkt)
      if (s_nkt < 0.05) confidence <- "MEDIA"

    # 11. GAMMA-DELTA T CELLS (TRDC alto, CD3+, CD4-/CD8-)
    } else if (s_gdt > 0.05 && s_cd3 > 0.5 &&
               s_cd4 < 0.5 && s_cd8 < 0.5) {
      label     <- "gamma-delta T cells"
      rationale <- sprintf("gdt_score=%.3f, CD3=%.3f, CD4=%.3f, CD8A=%.3f",
                           s_gdt, s_cd3, s_cd4, s_cd8)
      if (s_gdt < 0.10) confidence <- "MEDIA"

    # 12. MAIT CELLS (SLC4A10/KLRB1, CD3+)
    } else if (s_mait > 0.05 && s_cd3 > 0.5 &&
               s_mait > s_naive && s_mait > s_effector) {
      label     <- "MAIT cells"
      rationale <- sprintf("mait_score=%.3f, CD3=%.3f, CD8A=%.3f",
                           s_mait, s_cd3, s_cd8)
      if (s_mait < 0.10) confidence <- "MEDIA"

    # 13. ILC (GATA3/RORC, CD3- e CD56-)
    } else if (s_ilc > 0.05 && s_cd3 < 0.5 &&
               s_cd56 < 0.3 && s_ilc > s_nk) {
      label     <- "ILC"
      rationale <- sprintf("ilc_score=%.3f, CD3=%.3f, NK_score=%.3f",
                           s_ilc, s_cd3, s_nk)
      if (s_ilc < 0.10) confidence <- "MEDIA"

    # 14. NK CELLS
    } else if (s_nk > 0.10 && s_cd3 < 0.8 &&
               s_cd4 < 0.5 && s_cd8 < 0.5) {
      label     <- "NK cells"
      rationale <- sprintf("NK_score=%.3f, CD3=%.3f, CD4=%.3f, CD8A=%.3f",
                           s_nk, s_cd3, s_cd4, s_cd8)
      if (s_nk < 0.20) confidence <- "MEDIA"

    # 15. TREGS (FOXP3/IL2RA, CD4+)
    } else if (s_treg > 0.05 && s_treg > s_nk && s_cd4 > 0.3) {
      label     <- "Tregs"
      rationale <- sprintf("Treg_score=%.3f > NK_score=%.3f, CD4=%.3f",
                           s_treg, s_nk, s_cd4)
      if (s_treg < 0.10) confidence <- "MEDIA"

    # 16. PROLIFERATING
    } else if (s_mk67 > 0.5 || s_prolif > 0.05) {
      if (s_cd8 > s_cd4 && s_cd8 > 0.3) {
        label     <- "Proliferating CD8+ T cells"
        rationale <- sprintf("MKI67=%.3f, prolif=%.3f, CD8A=%.3f",
                             s_mk67, s_prolif, s_cd8)
      } else {
        label     <- "Proliferating T cells"
        rationale <- sprintf("MKI67=%.3f, prolif=%.3f, CD4=%.3f, CD8A=%.3f",
                             s_mk67, s_prolif, s_cd4, s_cd8)
      }
      if (s_mk67 < 1.0 && s_prolif < 0.10) confidence <- "MEDIA"

    # 17. CD8+ T CELLS
    } else if (s_cd8 > s_cd4 && s_cd8 > 0.3) {
      if (s_cytotox > s_naive && s_cytotox > 0.0) {
        label     <- "Cytotoxic CD8+ T cells"
        rationale <- sprintf("CD8A=%.3f, cytotox=%.3f > naive=%.3f",
                             s_cd8, s_cytotox, s_naive)
      } else if (s_naive > 0.1) {
        label      <- "Naive CD8+ T cells"
        rationale  <- sprintf("CD8A=%.3f, naive=%.3f, cytotox=%.3f",
                              s_cd8, s_naive, s_cytotox)
        confidence <- "MEDIA"
      } else {
        label      <- "Memory T cells"
        rationale  <- sprintf("CD8A=%.3f, naive=%.3f (low), cytotox=%.3f (low)",
                              s_cd8, s_naive, s_cytotox)
        confidence <- "MEDIA"
      }

    # 18. CD4+ T CELLS
    } else if (s_cd4 > 0.3) {
      if (s_naive > s_effector) {
        label     <- "Naive CD4+ T cells"
        rationale <- sprintf("CD4=%.3f, naive=%.3f > effector=%.3f",
                             s_cd4, s_naive, s_effector)
      } else {
        label     <- "Effector CD4+ T cells"
        rationale <- sprintf("CD4=%.3f, effector=%.3f > naive=%.3f",
                             s_cd4, s_effector, s_naive)
      }
      if (abs(s_naive - s_effector) < 0.02) confidence <- "BASSA"

    # 19. FALLBACK
    } else {
      scores_T <- c(
        "Naive CD4+ T cells"     = s_naive,
        "Effector CD4+ T cells"  = s_effector,
        "Cytotoxic CD8+ T cells" = s_cytotox,
        "Tregs"                  = s_treg,
        "NK cells"               = s_nk
      )
      label      <- names(which.max(scores_T))
      rationale  <- sprintf(
        "AMBIGUO: CD4=%.3f, CD8A=%.3f, CD3=%.3f. Fallback su score max (%s=%.3f)",
        s_cd4, s_cd8, s_cd3, label, max(scores_T))
      confidence <- "BASSA - VERIFICA MANUALE"
    }

    n_cells <- score_summary$n_cells[
      as.character(score_summary$seurat_clusters) == as.character(cl)]
    n_cells <- if (length(n_cells) == 0) "?" else n_cells[1]

    cat(sprintf("  C%s [%s cellule] -> %-35s [%s]\n    %s\n",
                cl, n_cells, label, confidence, rationale))

    decisions$label[decisions$cluster == cl]      <- label
    decisions$confidence[decisions$cluster == cl] <- confidence
    decisions$rationale[decisions$cluster == cl]  <- rationale
  }

  # ── 4. APPLICA ANNOTAZIONE ─────────────────────────────────
  Idents(obj) <- "seurat_clusters"
  cluster_ids    <- as.character(Idents(obj))
  annotation_map <- setNames(decisions$label, decisions$cluster)

  cell_labels <- unname(annotation_map[cluster_ids])
  names(cell_labels) <- colnames(obj)
  obj <- AddMetaData(obj, metadata = cell_labels, col.name = "cell_type")
  Idents(obj) <- "cell_type"

  cat(paste0("\n[", sample_name, "] Tabella annotazione finale:\n"))
  print(table(obj$cell_type))

  # ── 5. SALVA EXCEL DECISIONI ───────────────────────────────
  decisions_full <- merge(
    decisions,
    score_summary %>%
      mutate(cluster = as.character(seurat_clusters)) %>%
      select(-seurat_clusters),
    by = "cluster"
  )
  write.xlsx(decisions_full,
             paste0(out_dir, sample_name, "_annotation_decisions.xlsx"))
  cat(paste0("  -> ", sample_name, "_annotation_decisions.xlsx salvato\n"))

  return(list(obj = obj, decisions = decisions))
}

# ============================================================
# FUNZIONE UMAP
# ============================================================

plot_umap_ab <- function(obj, sample_name) {
  Idents(obj) <- "cell_type"
  cols <- get_colors(obj)

  p_label <- DimPlot(
    obj, reduction = "umap", label = TRUE,
    label.size = 3.5, repel = TRUE, cols = cols, pt.size = 0.6
  ) +
    ggtitle(paste0(sample_name, " - Con label")) +
    theme_classic(base_size = 12) +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 13),
      legend.text     = element_text(size = 9),
      legend.key.size = unit(0.45, "cm")
    ) +
    guides(color = guide_legend(override.aes = list(size = 4)))

  p_clean <- DimPlot(
    obj, reduction = "umap", label = FALSE,
    cols = cols, pt.size = 0.6
  ) +
    ggtitle(paste0(sample_name, " - Senza label")) +
    theme_classic(base_size = 12) +
    theme(
      plot.title      = element_text(hjust = 0.5, size = 12),
      legend.text     = element_text(size = 9),
      legend.key.size = unit(0.45, "cm")
    ) +
    guides(color = guide_legend(override.aes = list(size = 4)))

  path <- paste0(out_dir, sample_name, "_UMAP_annotated.png")
  ggsave(path, plot = p_label | p_clean,
         width = 16, height = 7, dpi = 300, bg = "white")
  cat(paste0("[", sample_name, "] UMAP -> ", path, "\n"))
  return(p_label)
}

# ============================================================
# LOOP PRINCIPALE
# ============================================================

annotated_AB <- list()
umap_list    <- list()

for (nm in names(ab_samples)) {
  obj    <- join_layers_safe(ab_samples[[nm]], nm)
  result <- auto_annotate(obj, nm)
  annotated_AB[[nm]] <- result$obj
  umap_list[[nm]]    <- plot_umap_ab(result$obj, nm)
}

# ============================================================
# PANNELLI COMPARATIVI
# ============================================================

section("Pannelli comparativi")

if (all(c("Ca_blood_AB","Ca_bone_AB") %in% names(umap_list))) {
  p_ca <- umap_list$Ca_blood_AB / umap_list$Ca_bone_AB +
    plot_annotation(title = "Ca - Blood AB vs Bone AB",
                    theme = theme(plot.title = element_text(
                      face="bold", size=14, hjust=0.5)))
  ggsave(paste0(out_dir, "Ca_AB_blood_vs_bone.png"),
         plot = p_ca, width = 9, height = 14, dpi = 300, bg = "white")
}

if (all(c("Bo_blood_AB","Bo_bone_AB") %in% names(umap_list))) {
  p_bo <- umap_list$Bo_blood_AB / umap_list$Bo_bone_AB +
    plot_annotation(title = "Bo - Blood AB vs Bone AB",
                    theme = theme(plot.title = element_text(
                      face="bold", size=14, hjust=0.5)))
  ggsave(paste0(out_dir, "Bo_AB_blood_vs_bone.png"),
         plot = p_bo, width = 9, height = 14, dpi = 300, bg = "white")
}

p_all <- patchwork::wrap_plots(umap_list, ncol = 1)
ggsave(paste0(out_dir, "ALL_AB_samples_UMAP.png"),
       plot = p_all, width = 10,
       height = length(umap_list) * 7, dpi = 300, bg = "white")
cat("[ALL] Pannello completo AB -> ALL_AB_samples_UMAP.png\n")

# ============================================================
# SALVATAGGIO AB
# ============================================================
section("Salvataggio")

rds_out <- paste0(base_dir, "all_AB_samples_annotated.rds")
saveRDS(annotated_AB, rds_out)
cat(paste0(">> AB salvato: ", rds_out, "\n"))

# ============================================================
# OGGETTO UNIFICATO I + AB
# ============================================================
section("Costruzione oggetto unificato I + AB")

rds_I <- NULL
for (candidate in c("all_samples_annotated_v3.rds",
                     "all_samples_annotated_v2.rds",
                     "all_samples_annotated_final.rds")) {
  path_candidate <- paste0(base_dir, candidate)
  if (file.exists(path_candidate)) {
    rds_I <- path_candidate
    cat(paste0(">> Campioni I trovati in: ", candidate, "\n"))
    break
  }
}

if (!is.null(rds_I)) {
  annotated_I <- readRDS(rds_I)

  name_map_I <- c("Bo_I" = "Bo_bone_I",
                  "Ca_I" = "Ca_bone_I",
                  "Me_I" = "Me_bone_I")
  names(annotated_I) <- ifelse(
    names(annotated_I) %in% names(name_map_I),
    name_map_I[names(annotated_I)],
    names(annotated_I)
  )
  cat("Campioni I rinominati:", paste(names(annotated_I), collapse=", "), "\n")

  all_samples_annotated <- c(annotated_I, annotated_AB)

  cat("\nOggetto unificato - riepilogo:\n")
  for (nm in names(all_samples_annotated)) {
    obj   <- all_samples_annotated[[nm]]
    types <- sort(unique(as.character(obj$cell_type)))
    cat(sprintf("  %-20s %5d cellule | %2d cell types: %s\n",
                nm, ncol(obj), length(types), paste(types, collapse=", ")))
  }

  rds_all <- paste0(base_dir, "all_samples_annotated_COMPLETE.rds")
  saveRDS(all_samples_annotated, rds_all)
  cat(paste0("\n>> Oggetto completo: all_samples_annotated_COMPLETE.rds\n"))

} else {
  cat("[WARN] File campioni I non trovato. Controlla che uno tra questi esista:\n")
  cat("  all_samples_annotated_v3.rds / v2.rds / final.rds\n")
  all_samples_annotated <- annotated_AB
}

# ============================================================
# RIEPILOGO CLUSTER A BASSA CONFIDENZA
# ============================================================
section("CLUSTER A BASSA CONFIDENZA - VERIFICA CONSIGLIATA")

for (nm in names(annotated_AB)) {
  cat(paste0("\n[", nm, "]\n"))
  tbl <- annotated_AB[[nm]]@meta.data %>%
    distinct(seurat_clusters, cell_type) %>%
    arrange(as.integer(as.character(seurat_clusters)))
  print(as.data.frame(tbl))
}

# ============================================================
# RIEPILOGO PALETTE FINALE
# ============================================================
section("PALETTE COLORI FINALE")

cat("Tipi FISSI (campioni I - colori invariati rispetto a STEP_3c):\n")
for (nm in names(palette_fixed))
  cat(sprintf("  %-35s %s\n", nm, palette_fixed[nm]))

new_found <- setdiff(names(full_palette), names(palette_fixed))
if (length(new_found) > 0) {
  cat("\nNuovi tipi trovati negli AB (colori aggiuntivi):\n")
  for (nm in new_found)
    cat(sprintf("  %-35s %s\n", nm, full_palette[nm]))
} else {
  cat("\nNessun nuovo tipo cellulare trovato negli AB.\n")
  cat("Tutti i cluster corrispondono a popolazioni gia presenti nei campioni I.\n")
}

cat(paste0(
  "\n", strrep("=", 65), "\n",
  "  NEW_AB COMPLETATO\n\n",
  "  Output: ", out_dir, "\n",
  "  RDS AB-only:  all_AB_samples_annotated.rds\n",
  "  RDS completo: all_samples_annotated_COMPLETE.rds\n",
  "  UMAP:         ALL_AB_samples_UMAP.png\n\n",
  "  Controlla i cluster BASSA nell'xlsx per verifica manuale.\n",
  strrep("=", 65), "\n"
))
