# ============================================================
#  PIPELINE 2 – Campioni AB (dopo infusione)
#  Copre: AB_auto_annotate.R con tutte le modifiche richieste
#
#  Prerequisito: PIPELINE_1 eseguita (all_I_samples_annotated.rds)
#
#  Modifiche rispetto alla versione precedente:
#    - "Dendritic Cells" → "Myeloid cells" + firma ampliata
#    - Th helper subtypes: Th1, Th2, Th17, Tfh nel ramo CD4
#    - Legenda consistente: CANONICAL_ORDER fisso,
#      tipi assenti marcati con "—"
#    - Palette condivisa con Pipeline 1 (stessi colori)
#
#  Decisione automatica per cluster:
#    ALTA:   segnale netto, soglie superate con ampio margine
#    MEDIA:  score borderline o Th subtypes
#            (citochine IL4/IL17A/IFNG basse a riposo)
#    BASSA:  ambiguo → verifica manuale richiesta
#
#  Output (tutti in base_dir/AB_annotation/):
#    <sample>_annotation_decisions.xlsx
#    <sample>_UMAP_annotated.png   (label | clean | legenda)
#    Ca_AB_blood_vs_bone.png
#    Bo_AB_blood_vs_bone.png
#    ALL_AB_samples_UMAP.png
#    all_AB_samples_annotated.rds
#    all_samples_annotated_COMPLETE.rds  (I + AB)
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

out_dir <- paste0(base_dir, "AB_annotation/")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

section <- function(title)
  cat(paste0("\n", strrep("=", 65), "\n  ", title,
             "\n", strrep("=", 65), "\n"))

# ============================================================
# CARICAMENTO DATI
# ============================================================
section("Caricamento dati")

if (!exists("seurat_list")) {
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

cat("\nCampioni AB:\n")
for (nm in names(ab_samples))
  cat(sprintf("  %-15s %d cellule | %d cluster\n",
              nm, ncol(ab_samples[[nm]]),
              length(unique(ab_samples[[nm]]$seurat_clusters))))

# JoinLayers Seurat v5
join_safe <- function(obj, nm) {
  if (length(grep("^counts\\.", Layers(obj), value = TRUE)) > 0) {
    obj <- JoinLayers(obj)
    cat(paste0("[", nm, "] JoinLayers applicato\n"))
  }
  Idents(obj) <- "seurat_clusters"
  obj
}

# ============================================================
# FIRME GENICHE
# ============================================================

SIG <- list(

  # Lineage markers diretti
  cd3    = c("CD3D","CD3E","CD3G"),
  cd4    = c("CD4"),
  cd8    = c("CD8A","CD8B"),
  cd19   = c("CD19","MS4A1"),
  cd56   = c("NCAM1"),
  cd14   = c("CD14"),
  cd16   = c("FCGR3A"),
  cd34   = c("CD34"),

  # T cell states
  naive    = c("CCR7","SELL","IL7R","TCF7","LEF1","KLF2"),
  effector = c("GZMK","CD44","S100A4","LGALS1","DUSP2"),
  cytotox  = c("GZMB","NKG7","PRF1","GNLY","GZMA","FGFBP2"),
  treg     = c("FOXP3","IL2RA","CTLA4","TIGIT","IKZF2",
               "TNFRSF18","TNFRSF4","ENTPD1","CCR8","RTKN2"),
  prolif   = c("MKI67","TOP2A","STMN1","PCNA","CCNB1","BIRC5","UBE2C"),

  # ── T helper subtypes ────────────────────────────────────
  # Th1: TBX21 (T-bet) + CXCR3 più stabili a riposo in scRNA-seq
  th1  = c("TBX21","CXCR3","CCR5","IL12RB2","STAT4",
           "IFNG","HAVCR2","PHLPP1","TNFSF10"),

  # Th2: GATA3 espresso anche da NK/ILC → protetto da CD4>0.3
  th2  = c("GATA3","CCR4","MAF","PTGDR2","IL4R",
           "IL4","IL13","HPGDS"),

  # Th17: RORC+CCR6 affidabili; IL17A/F spesso basse a riposo
  th17 = c("RORC","CCR6","IL23R","RORA","FURIN",
           "IL17A","IL17F","TMEM176A","TMEM176B","STAT3"),

  # Tfh: CXCR5+BCL6 i più specifici in scRNA-seq
  tfh  = c("CXCR5","BCL6","ICOS","PDCD1","SH2D1A",
           "IL21","CXCL13","TOX2","MAF","TIGIT"),

  # NK / innate lymphoid
  nk   = c("NCAM1","NCR1","KLRD1","KLRB1","GNLY",
            "FCGR3A","TYROBP","S1PR5","CX3CR1","XCL1","XCL2"),
  ilc  = c("GATA3","RORC","ICOS","IL1R1","KIT",
            "AREG","IL13","HPGDS"),
  nkt  = c("NCAM1","CD3D","KLRB1","ZBTB16","NKG7","GNLY","GZMB"),
  gdt  = c("TRDC","TRGC1","TRGC2","TRDV1","TRDV2",
            "TRGV9","KLRC1","KLRC2"),
  mait = c("SLC4A10","KLRB1","NCR3","RORC","IL18RAP","CXCR6"),

  # B lineage
  bcell  = c("MS4A1","CD19","CD79A","CD79B","PAX5",
              "BANK1","IGHM","IGHD","TCL1A"),
  bmem   = c("CD27","AIM2","TNFRSF13B","FCRL4","FCRL5","CD80"),
  plasma = c("JCHAIN","MZB1","IGHG1","IGHG2","IGKC","IGLC2",
              "SDC1","CD38","PRDM1","XBP1","DERL3"),

  # Monocytes
  mono14 = c("CD14","LYZ","S100A8","S100A9","FCN1","VCAN",
              "CXCL8","THBS1","PLBD1","CLEC7A","CTSS"),
  mono16 = c("FCGR3A","MS4A7","LILRB2","CDKN1C","CX3CR1",
              "LYPD2","CSF1R","HMOX1"),

  # ── Myeloid cells (ex "Dendritic Cells") ─────────────────
  # Bo_I C4 esprimeva LYZ/VCAN/S100A8-9/AQP9/FPR1/HCK/LILRB2
  # → marcatori mieloidi generici/moDC, non DC classici.
  # CLEC9A/XCR1/CD1C/LILRA4 assenti o minimi: cluster mieloide.
  # Nota: NO check effector<0.4 (S100A4 è fisiologico nei mieloidi)
  myeloid = c(
    "LYZ","VCAN","S100A8","S100A9","AQP9","FPR1",
    "LILRB2","MS4A6A","HCK","CSF1R","CD68","ITGAM",
    "FCN1","CXCL8","THBS1","PLBD1","TYROBP","SPI1",
    "MRC1","C1QA","C1QB","C1QC",
    "HLA-DRA","HLA-DPB1","HLA-DQA1",
    "CD1C","FCER1A","CLEC9A","XCR1","LILRA4",
    "ITGAX","BATF3","CLEC10A","SIGLEC6"
  ),

  # Other
  baso     = c("CPA3","TPSAB1","TPSB2","MS4A2","GATA2",
               "HDC","SLC18A2","KIT","FCER1A"),
  hspc     = c("CD34","SPINK2","AVP","CRHBP","HOPX",
               "MLLT3","MECOM","PROM1"),
  erythro  = c("HBB","HBA1","HBA2","GYPA","GYPB",
               "ALAS2","SLC4A1","CA1","AHSP","KLF1"),
  platelet = c("PPBP","PF4","GP1BA","ITGA2B","ITGB3",
               "GNG11","TUBB1","TREML1","CMTM5")
)

# ============================================================
# PALETTE E LEGENDA CONSISTENTE
# ============================================================

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
    cat("[INFO] Colori auto:", paste(missing, collapse=", "), "\n")
  }
  cols
}

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
    theme(plot.background = element_rect(fill = "white", color = NA),
          plot.margin = margin(8, 4, 8, 4))
}

# ============================================================
# HELPERS
# ============================================================

mean_genes <- function(mat, genes, cluster) {
  g <- genes[genes %in% rownames(mat)]
  if (length(g) == 0 || !cluster %in% colnames(mat)) return(0)
  mean(mat[g, cluster])
}

get_score <- function(score_summary, cluster_id, score_name) {
  col <- paste0(score_name, "_score")
  if (!col %in% colnames(score_summary)) return(-99)
  val <- score_summary[[col]][
    as.character(score_summary$seurat_clusters) == as.character(cluster_id)]
  if (length(val) == 0) return(-99)
  val[1]
}

fix_colnames <- function(mat) {
  colnames(mat) <- gsub("^g", "", colnames(mat))
  colnames(mat) <- gsub("^RNA_snn_res\\.[0-9.]+_", "", colnames(mat))
  mat
}

# ============================================================
# FUNZIONE DI ANNOTAZIONE AUTOMATICA
# Gerarchia a 19 livelli (in ordine di priorità decrescente):
#  1. Platelets     2. Erythroid    3. HSPC
#  4. Plasma cells  5. B cells      6. Basophils
#  7. CD14 Mono     8. CD16 Mono    9. Myeloid cells
# 10. NKT          11. gamma-delta 12. MAIT
# 13. ILC          14. NK          15. Tregs
# 16. Proliferating 17. CD8+ T     18. CD4+ T (Th subtypes)
# 19. Fallback
# ============================================================

auto_annotate <- function(obj, sample_name) {

  section(paste0("AUTO-ANNOTAZIONE: ", sample_name))

  Idents(obj) <- "seurat_clusters"
  # Fix Seurat v5: as.integer(as.character()) evita il bug su cluster 0
  clusters_chr <- as.character(
    sort(as.integer(as.character(unique(obj$seurat_clusters)))))

  # ── AverageExpression ───────────────────────────────────────
  all_sig_genes <- unique(unlist(SIG))
  genes_present <- all_sig_genes[all_sig_genes %in% rownames(obj)]

  avg <- fix_colnames(
    AverageExpression(obj, features = genes_present,
                      group.by = "seurat_clusters",
                      assay = "RNA", slot = "data")$RNA)

  # ── AddModuleScore per tutte le firme ──────────────────────
  score_names <- c("naive","effector","cytotox","treg","prolif",
                   "th1","th2","th17","tfh",
                   "nk","ilc","nkt","gdt","mait",
                   "bcell","bmem","plasma",
                   "mono14","mono16","myeloid",
                   "baso","hspc","erythro","platelet")

  for (sn in score_names) {
    genes_ok <- SIG[[sn]][SIG[[sn]] %in% rownames(obj)]
    if (length(genes_ok) >= 3) {
      score_col <- paste0(sn, "_score")
      obj <- AddModuleScore(obj, features = list(genes_ok),
                            name = score_col)
      obj[[score_col]]              <- obj[[paste0(score_col, "1")]]
      obj[[paste0(score_col, "1")]] <- NULL
    }
  }

  score_cols   <- paste0(score_names, "_score")
  score_cols   <- score_cols[score_cols %in% colnames(obj@meta.data)]

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

  # ── Decisione per cluster ───────────────────────────────────
  decisions <- data.frame(
    cluster    = clusters_chr,
    label      = NA_character_,
    confidence = NA_character_,
    rationale  = NA_character_,
    stringsAsFactors = FALSE
  )

  cat(paste0("\n[", sample_name, "] Decisioni:\n"))

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

    s_naive    <- gs(cl,"naive");    s_effector <- gs(cl,"effector")
    s_cytotox  <- gs(cl,"cytotox");  s_treg     <- gs(cl,"treg")
    s_prolif   <- gs(cl,"prolif")
    s_th1      <- gs(cl,"th1");      s_th2      <- gs(cl,"th2")
    s_th17     <- gs(cl,"th17");     s_tfh      <- gs(cl,"tfh")
    s_nk       <- gs(cl,"nk");       s_ilc      <- gs(cl,"ilc")
    s_nkt      <- gs(cl,"nkt");      s_gdt      <- gs(cl,"gdt")
    s_mait     <- gs(cl,"mait")
    s_bcell    <- gs(cl,"bcell");    s_bmem     <- gs(cl,"bmem")
    s_plasma   <- gs(cl,"plasma")
    s_mono14   <- gs(cl,"mono14");   s_mono16   <- gs(cl,"mono16")
    s_myeloid  <- gs(cl,"myeloid")
    s_baso     <- gs(cl,"baso");     s_hspc     <- gs(cl,"hspc")
    s_erythro  <- gs(cl,"erythro"); s_platelet <- gs(cl,"platelet")

    label      <- NA_character_
    confidence <- "ALTA"
    rationale  <- ""

    # ── 1. PLATELETS ─────────────────────────────────────────
    if (s_platelet > 0.30 && s_cd3 < 0.3 && s_cd19 < 0.3) {
      label     <- "Platelets"
      rationale <- sprintf("platelet=%.3f, CD3=%.3f", s_platelet, s_cd3)
      if (s_platelet < 0.50) confidence <- "MEDIA"

    # ── 2. ERYTHROID ─────────────────────────────────────────
    } else if (s_erythro > 0.30 && s_cd3 < 0.3 && s_cd19 < 0.3) {
      label     <- "Erythroid cells"
      rationale <- sprintf("erythro=%.3f, CD3=%.3f", s_erythro, s_cd3)
      if (s_erythro < 0.50) confidence <- "MEDIA"

    # ── 3. HSPC ──────────────────────────────────────────────
    } else if (s_hspc > 0.15 && s_cd34 > 0.3 &&
               s_cd3 < 0.3 && s_cd19 < 0.3) {
      label     <- "HSPC"
      rationale <- sprintf("hspc=%.3f, CD34=%.3f", s_hspc, s_cd34)
      if (s_hspc < 0.25) confidence <- "MEDIA"

    # ── 4. PLASMA CELLS ──────────────────────────────────────
    } else if (s_plasma > 0.20 && s_cd3 < 0.3) {
      label     <- "Plasma cells"
      rationale <- sprintf("plasma=%.3f, CD3=%.3f", s_plasma, s_cd3)
      if (s_plasma < 0.35) confidence <- "MEDIA"

    # ── 5. B CELLS / MEMORY B ────────────────────────────────
    } else if (s_cd19 > 0.5 && s_cd3 < 0.5) {
      if (s_bmem > s_bcell && s_bmem > 0.05) {
        label     <- "Memory B cells"
        rationale <- sprintf("CD19=%.3f, bmem=%.3f > bcell=%.3f",
                             s_cd19, s_bmem, s_bcell)
      } else {
        label     <- "B cells"
        rationale <- sprintf("CD19=%.3f, bcell=%.3f", s_cd19, s_bcell)
      }
      if (s_cd19 < 0.8) confidence <- "MEDIA"

    # ── 6. BASOPHILS ─────────────────────────────────────────
    } else if (s_baso > 0.15 && s_cd3 < 0.3 && s_cd19 < 0.3) {
      label     <- "Basophils"
      rationale <- sprintf("baso=%.3f, CD3=%.3f", s_baso, s_cd3)
      if (s_baso < 0.30) confidence <- "MEDIA"

    # ── 7. CD14 MONOCYTES ────────────────────────────────────
    } else if (s_cd14 > 0.5 && s_mono14 > 0.15 &&
               s_cd3 < 0.3 && s_cd56 < 0.3) {
      label     <- "CD14 Monocytes"
      rationale <- sprintf("CD14=%.3f, mono14=%.3f", s_cd14, s_mono14)
      if (s_mono14 < 0.25) confidence <- "MEDIA"

    # ── 8. CD16 MONOCYTES ────────────────────────────────────
    } else if (s_cd16 > 0.3 && s_mono16 > 0.10 &&
               s_cd3 < 0.3 && s_cd14 < 0.5) {
      label     <- "CD16 Monocytes"
      rationale <- sprintf("CD16=%.3f, mono16=%.3f", s_cd16, s_mono16)
      if (s_mono16 < 0.20) confidence <- "MEDIA"

    # ── 9. MYELOID CELLS ─────────────────────────────────────
    # Ex "Dendritic Cells". Nessun check su effector_score:
    # S100A4 (firma effector) è fisiologico nei mieloidi.
    } else if (s_myeloid > 0.35 &&
               s_cd3 < 0.3 && s_cd19 < 0.3 && s_cd56 < 0.3) {
      label     <- "Myeloid cells"
      rationale <- sprintf("myeloid=%.3f, CD3=%.3f, CD19=%.3f",
                           s_myeloid, s_cd3, s_cd19)
      if (s_myeloid < 0.55) confidence <- "MEDIA"

    # ── 10. NKT ──────────────────────────────────────────────
    } else if (s_cd3 > 0.5 && s_cd56 > 0.3 && s_nkt > s_nk) {
      label     <- "NKT cells"
      rationale <- sprintf("CD3=%.3f, CD56=%.3f, nkt=%.3f",
                           s_cd3, s_cd56, s_nkt)
      if (s_nkt < 0.05) confidence <- "MEDIA"

    # ── 11. GAMMA-DELTA T ────────────────────────────────────
    } else if (s_gdt > 0.05 && s_cd3 > 0.5 &&
               s_cd4 < 0.5 && s_cd8 < 0.5) {
      label     <- "gamma-delta T cells"
      rationale <- sprintf("gdt=%.3f, CD3=%.3f", s_gdt, s_cd3)
      if (s_gdt < 0.10) confidence <- "MEDIA"

    # ── 12. MAIT ─────────────────────────────────────────────
    } else if (s_mait > 0.05 && s_cd3 > 0.5 &&
               s_mait > s_naive && s_mait > s_effector) {
      label     <- "MAIT cells"
      rationale <- sprintf("mait=%.3f, CD3=%.3f", s_mait, s_cd3)
      if (s_mait < 0.10) confidence <- "MEDIA"

    # ── 13. ILC ──────────────────────────────────────────────
    } else if (s_ilc > 0.05 && s_cd3 < 0.5 &&
               s_cd56 < 0.3 && s_ilc > s_nk) {
      label     <- "ILC"
      rationale <- sprintf("ilc=%.3f, CD3=%.3f", s_ilc, s_cd3)
      if (s_ilc < 0.10) confidence <- "MEDIA"

    # ── 14. NK CELLS ─────────────────────────────────────────
    } else if (s_nk > 0.10 && s_cd3 < 0.8 &&
               s_cd4 < 0.5 && s_cd8 < 0.5) {
      label     <- "NK cells"
      rationale <- sprintf("NK=%.3f, CD3=%.3f", s_nk, s_cd3)
      if (s_nk < 0.20) confidence <- "MEDIA"

    # ── 15. TREGS ────────────────────────────────────────────
    } else if (s_treg > 0.05 && s_treg > s_nk && s_cd4 > 0.3) {
      label     <- "Tregs"
      rationale <- sprintf("Treg=%.3f, NK=%.3f, CD4=%.3f",
                           s_treg, s_nk, s_cd4)
      if (s_treg < 0.10) confidence <- "MEDIA"

    # ── 16. PROLIFERATING ────────────────────────────────────
    # Logica CD4/CD8:
    #   - Calcola differenza relativa tra CD4 e CD8
    #   - Se delta > 20% → assegna al lineage dominante
    #   - Se delta ≤ 20% (ambiguo) → "Proliferating T cells"
    # Nota: MKI67/TOP2A bassa espressione media per cluster è
    # normale in scRNA-seq (solo le cellule in fase S/G2M sono
    # alte, ma vengono diluiti nella media del cluster).
    } else if (s_mk67 > 0.5 || s_prolif > 0.05) {
      prolif_denom <- max(s_cd4, s_cd8, 1e-6)
      prolif_delta <- abs(s_cd4 - s_cd8) / prolif_denom

      if (prolif_delta >= 0.20 && s_cd8 > s_cd4 && s_cd8 > 0.2) {
        label     <- "Proliferating CD8+ T cells"
        rationale <- sprintf("MKI67=%.3f, prolif=%.3f, CD8=%.3f>CD4=%.3f (delta=%.0f%%)",
                             s_mk67, s_prolif, s_cd8, s_cd4,
                             prolif_delta * 100)
      } else if (prolif_delta >= 0.20 && s_cd4 > s_cd8 && s_cd4 > 0.2) {
        label     <- "Proliferating CD4+ T cells"
        rationale <- sprintf("MKI67=%.3f, prolif=%.3f, CD4=%.3f>CD8=%.3f (delta=%.0f%%)",
                             s_mk67, s_prolif, s_cd4, s_cd8,
                             prolif_delta * 100)
      } else {
        label     <- "Proliferating T cells"
        rationale <- sprintf("MKI67=%.3f, prolif=%.3f, CD4=%.3f≈CD8=%.3f (delta=%.0f%% <20%%)",
                             s_mk67, s_prolif, s_cd4, s_cd8,
                             prolif_delta * 100)
      }
      if (s_mk67 < 1.0 && s_prolif < 0.10) confidence <- "MEDIA"

    # ── 17. CD8+ T CELLS ─────────────────────────────────────
    } else if (s_cd8 > s_cd4 && s_cd8 > 0.3) {
      if (s_cytotox > s_naive && s_cytotox > 0.0) {
        label     <- "Cytotoxic CD8+ T cells"
        rationale <- sprintf("CD8=%.3f, cytotox=%.3f > naive=%.3f",
                             s_cd8, s_cytotox, s_naive)
      } else if (s_naive > 0.1) {
        label      <- "Naive CD8+ T cells"
        confidence <- "MEDIA"
        rationale  <- sprintf("CD8=%.3f, naive=%.3f", s_cd8, s_naive)
      } else {
        label      <- "Memory T cells"
        confidence <- "MEDIA"
        rationale  <- sprintf("CD8=%.3f, low naive+cytotox", s_cd8)
      }

    # ── 18. CD4+ T CELLS con sub-gerarchia T helper ──────────
    #
    # Ordine: Naive (CCR7/SELL) → Tfh (CXCR5/BCL6) →
    #         Th17 (RORC/CCR6) → Th1 (TBX21/CXCR3) →
    #         Th2 (GATA3/CCR4) → Effector CD4+ (fallback)
    #
    # Tutti i Th subtypes → confidence MEDIA perché citochine
    # effettrici (IFNG, IL17A, IL4) basse a riposo in scRNA-seq.
    } else if (s_cd4 > 0.3) {
      base_conf <- if (s_cd4 < 0.5) "MEDIA" else "ALTA"

      if (s_naive > s_effector && s_naive > 0.05) {
        label      <- "Naive CD4+ T cells"
        confidence <- base_conf
        rationale  <- sprintf("CD4=%.3f, naive=%.3f > eff=%.3f",
                               s_cd4, s_naive, s_effector)

      } else if (s_tfh > 0.08 && s_tfh > s_naive) {
        label      <- "Tfh cells"
        confidence <- "MEDIA"
        rationale  <- sprintf("CD4=%.3f, tfh=%.3f > naive=%.3f",
                               s_cd4, s_tfh, s_naive)

      } else if (s_th17 > 0.06 && s_th17 > s_naive &&
                 s_th17 > s_th1 && s_th17 > s_th2) {
        label      <- "Th17 cells"
        confidence <- "MEDIA"
        rationale  <- sprintf("CD4=%.3f, th17=%.3f>th1=%.3f,th2=%.3f",
                               s_cd4, s_th17, s_th1, s_th2)

      } else if (s_th1 > 0.06 && s_th1 > s_naive && s_th1 > s_th2) {
        label      <- "Th1 cells"
        confidence <- "MEDIA"
        rationale  <- sprintf("CD4=%.3f, th1=%.3f>th2=%.3f,naive=%.3f",
                               s_cd4, s_th1, s_th2, s_naive)

      } else if (s_th2 > 0.06 && s_th2 > s_naive) {
        label      <- "Th2 cells"
        confidence <- "MEDIA"
        rationale  <- sprintf("CD4=%.3f, th2=%.3f>naive=%.3f",
                               s_cd4, s_th2, s_naive)

      } else {
        label      <- "Effector CD4+ T cells"
        confidence <- if (abs(s_naive - s_effector) < 0.02)
                        "BASSA" else base_conf
        rationale  <- sprintf("CD4=%.3f, eff=%.3f≥naive=%.3f (no Th)",
                               s_cd4, s_effector, s_naive)
      }

    # ── 19. FALLBACK ─────────────────────────────────────────
    } else {
      sc_T <- c("Naive CD4+ T cells"     = s_naive,
                "Effector CD4+ T cells"  = s_effector,
                "Cytotoxic CD8+ T cells" = s_cytotox,
                "Tregs"                  = s_treg,
                "NK cells"               = s_nk)
      label      <- names(which.max(sc_T))
      confidence <- "BASSA - VERIFICA MANUALE"
      rationale  <- sprintf("AMBIGUO: CD4=%.3f,CD8=%.3f,CD3=%.3f. Max=%s",
                             s_cd4, s_cd8, s_cd3, label)
    }

    n_cells <- score_summary$n_cells[
      as.character(score_summary$seurat_clusters) == as.character(cl)]
    n_cells <- if (length(n_cells) == 0) "?" else n_cells[1]

    cat(sprintf("  C%s [%s celle] → %-35s [%s]\n    %s\n",
                cl, n_cells, label, confidence, rationale))

    decisions$label[decisions$cluster == cl]      <- label
    decisions$confidence[decisions$cluster == cl] <- confidence
    decisions$rationale[decisions$cluster == cl]  <- rationale
  }

  # ── Applica annotazione ────────────────────────────────────
  Idents(obj) <- "seurat_clusters"
  cl_ids      <- as.character(Idents(obj))
  ann_map     <- setNames(decisions$label, decisions$cluster)
  cell_labels <- unname(ann_map[cl_ids])
  names(cell_labels) <- colnames(obj)
  obj <- AddMetaData(obj, cell_labels, "cell_type")
  Idents(obj) <- "cell_type"

  cat(paste0("\n[", sample_name, "] Tabella finale:\n"))
  print(table(obj$cell_type))

  # ── Salva Excel decisioni ──────────────────────────────────
  decisions_full <- merge(
    decisions,
    score_summary %>%
      mutate(cluster = as.character(seurat_clusters)) %>%
      select(-seurat_clusters),
    by = "cluster"
  )
  write.xlsx(decisions_full,
             paste0(out_dir, sample_name, "_annotation_decisions.xlsx"))
  cat(paste0("  → ", sample_name, "_annotation_decisions.xlsx\n"))

  return(list(obj = obj, decisions = decisions))
}

# ============================================================
# FUNZIONE UMAP CON LEGENDA CONSISTENTE
# ============================================================

plot_umap_ab <- function(obj, sample_name) {
  present <- sort(unique(as.character(obj$cell_type)))
  cols    <- get_colors(present)

  make_dp <- function(show_label) {
    DimPlot(obj, reduction = "umap", group.by = "cell_type",
            label = show_label, label.size = 3.2, repel = TRUE,
            cols = cols, pt.size = 0.6) +
      ggtitle(paste0(sample_name,
                     if (show_label) " – Con label" else " – Senza label")) +
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
                                                 hjust = 0.5, size = 14)))
  path <- paste0(out_dir, sample_name, "_UMAP_annotated.png")
  ggsave(path, plot = combined, width = 20, height = 7,
         dpi = 300, bg = "white")
  cat(paste0("[", sample_name, "] UMAP → ", path, "\n"))
  return(invisible(combined))
}

# ============================================================
# LOOP PRINCIPALE
# ============================================================
section("Annotazione campioni AB")

annotated_AB <- list()

for (nm in names(ab_samples)) {
  obj                <- join_safe(ab_samples[[nm]], nm)
  result             <- auto_annotate(obj, nm)
  annotated_AB[[nm]] <- result$obj
  plot_umap_ab(result$obj, nm)
}

# ============================================================
# PANNELLI COMPARATIVI
# ============================================================
section("Pannelli comparativi")

make_mini_dp <- function(obj, title_lbl) {
  cols <- get_colors(sort(unique(as.character(obj$cell_type))))
  DimPlot(obj, reduction = "umap", group.by = "cell_type",
          label = TRUE, label.size = 3, repel = TRUE,
          cols = cols, pt.size = 0.5) +
    ggtitle(title_lbl) + theme_classic(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "none")
}

make_comparison <- function(nm1, nm2, title, filename) {
  if (!all(c(nm1, nm2) %in% names(annotated_AB))) return(invisible(NULL))
  all_pres <- union(unique(as.character(annotated_AB[[nm1]]$cell_type)),
                    unique(as.character(annotated_AB[[nm2]]$cell_type)))
  combined <- (make_mini_dp(annotated_AB[[nm1]], nm1) |
               make_mini_dp(annotated_AB[[nm2]], nm2) |
               make_full_legend(all_pres)) +
    plot_layout(widths = c(5, 5, 2.5)) +
    plot_annotation(title = title,
                    theme = theme(plot.title =
                                    element_text(face = "bold",
                                                 hjust = 0.5, size = 14)))
  ggsave(paste0(out_dir, filename), plot = combined,
         width = 20, height = 7, dpi = 300, bg = "white")
  cat(paste0("[", title, "] → ", filename, "\n"))
}

make_comparison("Ca_blood_AB","Ca_bone_AB",
                "Ca – Blood AB vs Bone AB","Ca_AB_blood_vs_bone.png")
make_comparison("Bo_blood_AB","Bo_bone_AB",
                "Bo – Blood AB vs Bone AB","Bo_AB_blood_vs_bone.png")

# Pannello completo tutti gli AB
all_pres_global <- unique(unlist(lapply(annotated_AB, function(o)
  unique(as.character(o$cell_type)))))

p_list <- lapply(names(annotated_AB), function(nm)
  make_mini_dp(annotated_AB[[nm]], nm))

p_all <- (patchwork::wrap_plots(p_list, ncol = 1) |
          make_full_legend(all_pres_global)) +
  plot_layout(widths = c(10, 2.5))

ggsave(paste0(out_dir, "ALL_AB_samples_UMAP.png"),
       plot = p_all,
       width = 14,
       height = length(annotated_AB) * 7,
       dpi = 300, bg = "white")
cat("[ALL] → ALL_AB_samples_UMAP.png\n")

# ============================================================
# SALVATAGGIO
# ============================================================
section("Salvataggio")

saveRDS(annotated_AB,
        paste0(base_dir, "all_AB_samples_annotated.rds"))
cat(paste0(">> AB: all_AB_samples_annotated.rds\n"))

############################################################
# Unifica con campioni I
rds_I_path <- paste0(base_dir, "all_I_samples_annotated.rds")
if (file.exists(rds_I_path)) {
  annotated_I <- readRDS(rds_I_path)
  all_samples <- c(annotated_I, annotated_AB)
  saveRDS(all_samples,
          paste0(base_dir, "all_samples_annotated_COMPLETE.rds"))
  cat(">> COMPLETE: all_samples_annotated_COMPLETE.rds\n")

  cat("\nRiepilogo oggetto unificato:\n")
  for (nm in names(all_samples)) {
    types <- sort(unique(as.character(all_samples[[nm]]$cell_type)))
    cat(sprintf("  %-20s %5d celle | %d tipi: %s\n",
                nm, ncol(all_samples[[nm]]), length(types),
                paste(types, collapse=", ")))
  }
} else {
  cat("[WARN] all_I_samples_annotated.rds non trovato.\n")
  cat("       Esegui PIPELINE_1 prima di questa pipeline.\n")
}

# Riepilogo cluster a bassa confidenza
section("CLUSTER A BASSA CONFIDENZA – VERIFICA CONSIGLIATA")
for (nm in names(annotated_AB)) {
  bassa <- annotated_AB[[nm]]@meta.data %>%
    distinct(seurat_clusters, cell_type) %>%
    arrange(as.integer(as.character(seurat_clusters)))
  cat(paste0("\n[", nm, "]:\n"))
  print(as.data.frame(bassa))
}

cat(paste0(
  "\n", strrep("=", 65), "\n",
  "  PIPELINE 2 COMPLETATA\n\n",
  "  Output: ", out_dir, "\n",
  "  RDS AB:       all_AB_samples_annotated.rds\n",
  "  RDS COMPLETE: all_samples_annotated_COMPLETE.rds\n\n",
  "  Modifiche applicate:\n",
  "  - Dendritic Cells → Myeloid cells\n",
  "  - T helper subtypes: Th1/Th2/Th17/Tfh nel ramo CD4\n",
  "  - Legenda consistente con tipi assenti marcati (—)\n\n",
  "  Prossimo step: PIPELINE_3_compare_I_vs_AB.R\n",
  strrep("=", 65), "\n"
))
