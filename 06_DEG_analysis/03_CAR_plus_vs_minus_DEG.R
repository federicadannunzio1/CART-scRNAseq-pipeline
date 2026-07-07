# CAR+ vs CAR- Differential Expression — Top Genes
# Converted from P3_top_genes_CAR_plus.Rmd

# ======================================================================
# LOAD LIBRARIES
# ======================================================================
suppressMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(readr)
  library(DT)
  library(knitr)
  library(ggrepel)
  library(RColorBrewer)
  library(Matrix)
})

# Enrichment tools
if (!requireNamespace("clusterProfiler", quietly = TRUE))
  BiocManager::install("clusterProfiler", ask = FALSE)
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
  BiocManager::install("org.Hs.eg.db", ask = FALSE)
if (!requireNamespace("AnnotationDbi", quietly = TRUE))
  BiocManager::install("AnnotationDbi", ask = FALSE)

suppressMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

if (!requireNamespace("UpSetR", quietly = TRUE))
  install.packages("UpSetR", repos = "https://cran.rstudio.com/")
if (!requireNamespace("writexl", quietly = TRUE))
  install.packages("writexl", repos = "https://cran.rstudio.com/")

suppressMessages({
  library(UpSetR)
  library(writexl)
})

# ======================================================================
# DEFINE PATHS
# ======================================================================
# ── Percorsi dati ────────────────────────────────────────────
RDS_I  <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/CART/Code/Seurat_analysis/2_annotation/all_I_samples_annotated.rds"
RDS_AB <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/CART/Code/Seurat_analysis/2_annotation/all_AB_samples_annotated.rds"

# ── T cell types: CAR- è definito come sole cellule T ────────
# Il gruppo CAR+ è per definizione interamente composto da T cells.
# Restringiamo il gruppo CAR- alle sole T cells per evitare il confounding
# dovuto alla composizione cellulare eterogenea (monociti, NK, B cells).
T_CELL_TYPES <- c(
  "Cytotoxic CD8+ T cells", "Effector CD4+ T cells", "Memory T cells",
  "Naive CD4+ T cells", "Naive CD8+ T cells", "Proliferating CD4+ T cells",
  "Proliferating CD8+ T cells", "Tregs", "NKT cells", "gamma-delta T cells"
)

# ── Geni da escludere ────────────────────────────────────────
# Non rimossi dalla matrice (preserva la normalizzazione),
# filtrati dai risultati DE.
NOISE_PATTERNS <- c("^MT-", "^RPS", "^RPL")
NOISE_EXACT    <- c("MALAT1", "NEAT1")

filter_noise <- function(genes) {
  is_noise <- grepl(paste(NOISE_PATTERNS, collapse = "|"), genes) |
              genes %in% NOISE_EXACT
  genes[!is_noise]
}

# ── Struttura output ─────────────────────────────────────────
OUT_DIR <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/CART/Code/Seurat_analysis/11_integrated_pipeline/P3_results"

# ── Parametri analisi ───────────────────────────────────────
MIN_CAR_CELLS <- 10     # minimo CAR+ per includere il campione
DE_LOG2FC     <- 0.5    # soglia log2FC per DE markers
DE_PADJ       <- 0.05   # soglia FDR per DE markers
ENRICH_PVAL   <- 0.05   # soglia p.adjust per arricchimento
MIN_GENE_SET  <- 10
MAX_GENE_SET  <- 500

# ── Colori ──────────────────────────────────────────────────
PATIENT_COLORS <- c(Bo = "#E64B35", Ca = "#4DBBD5", Me = "#00A087")
TIMEPOINT_COLORS <- c(
  "I (Infusion)"   = "#7B68EE",
  "Blood AB"       = "#E64B35",
  "Bone marrow AB" = "#00A087"
)

# ======================================================================
# HELPERS
# ======================================================================
get_car_status <- function(meta) {
  for (col in c("IS_CAR_ALLIN_scREP","IS_CAR","CAR")) {
    if (col %in% colnames(meta)) {
      vals <- as.character(meta[[col]])
      return(ifelse(grepl("^(YES|TRUE|yes|true|1)$", vals), "CAR+", "CAR-"))
    }
  }
  rep("CAR-", nrow(meta))
}

symbols_to_entrez <- function(gene_symbols) {
  bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID",
       OrgDb = org.Hs.eg.db, drop = TRUE)
}

run_go_enrichment <- function(gene_symbols, label = "") {
  eg <- symbols_to_entrez(gene_symbols)
  if (nrow(eg) < 5) { message(label, ": troppo pochi geni mappabili"); return(NULL) }
  res <- enrichGO(
    gene = eg$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP",
    pAdjustMethod = "BH", pvalueCutoff = ENRICH_PVAL, qvalueCutoff = 0.2,
    readable = TRUE, minGSSize = MIN_GENE_SET, maxGSSize = MAX_GENE_SET
  )
  if (is.null(res) || nrow(as.data.frame(res)) == 0) {
    message(label, ": nessun termine GO arricchito"); return(NULL)
  }
  res
}

run_kegg_enrichment <- function(gene_symbols, label = "") {
  eg <- symbols_to_entrez(gene_symbols)
  if (nrow(eg) < 5) { message(label, ": troppo pochi geni mappabili"); return(NULL) }
  res <- enrichKEGG(
    gene = eg$ENTREZID, organism = "hsa",
    pAdjustMethod = "BH", pvalueCutoff = ENRICH_PVAL, qvalueCutoff = 0.2,
    minGSSize = MIN_GENE_SET, maxGSSize = MAX_GENE_SET
  )
  if (is.null(res) || nrow(as.data.frame(res)) == 0) {
    message(label, ": nessun pathway KEGG arricchito"); return(NULL)
  }
  tryCatch(
    setReadable(res, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"),
    error = function(e) res
  )
}

# ======================================================================
# LOAD I
# ======================================================================
cat("Caricamento campioni I...\n")
I_samples_all <- readRDS(RDS_I)

for (sname in names(I_samples_all)) {
  meta <- I_samples_all[[sname]]@meta.data
  I_samples_all[[sname]]$car_status <- get_car_status(meta)
  I_samples_all[[sname]]$patient    <- sub("_I$|_bone_I$","", sname)
}

summary_I <- bind_rows(lapply(names(I_samples_all), function(s) {
  meta  <- I_samples_all[[s]]@meta.data
  n_car <- sum(meta$car_status == "CAR+")
  data.frame(Campione = s, Paziente = unique(meta$patient),
             Totale = nrow(meta), CAR_pos = n_car,
             Perc_CAR = round(100 * n_car / nrow(meta), 1),
             Incluso = ifelse(n_car >= MIN_CAR_CELLS, "SI", "NO"))
}))

I_valid <- I_samples_all[sapply(I_samples_all, function(s)
  sum(s@meta.data$car_status == "CAR+") >= MIN_CAR_CELLS)]

cat("Campioni I validi:", paste(names(I_valid), collapse = ", "), "\n")

# ======================================================================
# SUMMARY I
# ======================================================================
kable(summary_I, caption = "Campioni I — overview CAR+", align = "c")

# ======================================================================
# LOAD BLOOD AB
# ======================================================================
cat("Caricamento campioni blood AB...\n")
AB_raw <- readRDS(RDS_AB)

blood_names <- grep("blood_AB", names(AB_raw), value = TRUE)
AB_blood    <- AB_raw[blood_names]

for (sname in names(AB_blood)) {
  meta <- AB_blood[[sname]]@meta.data
  AB_blood[[sname]]$car_status <- get_car_status(meta)
  AB_blood[[sname]]$patient    <- sub("_blood_AB$","", sname)
}

summary_blood <- bind_rows(lapply(names(AB_blood), function(s) {
  meta  <- AB_blood[[s]]@meta.data
  n_car <- sum(meta$car_status == "CAR+")
  data.frame(Campione = s, Paziente = unique(meta$patient),
             Totale = nrow(meta), CAR_pos = n_car,
             Perc_CAR = round(100 * n_car / nrow(meta), 1),
             Incluso = ifelse(n_car >= MIN_CAR_CELLS, "SI", "NO"))
}))

AB_blood_valid <- AB_blood[sapply(AB_blood, function(s)
  sum(s@meta.data$car_status == "CAR+") >= MIN_CAR_CELLS)]

cat("Campioni blood validi:", paste(names(AB_blood_valid), collapse = ", "), "\n")

# ======================================================================
# SUMMARY BLOOD
# ======================================================================
kable(summary_blood, caption = "Campioni blood AB — overview CAR+", align = "c")

# ======================================================================
# LOAD BONE AB
# ======================================================================
cat("Caricamento campioni bone marrow AB...\n")
bone_names <- grep("bone_AB", names(AB_raw), value = TRUE)
AB_bone    <- AB_raw[bone_names]
rm(AB_raw); invisible(gc())

for (sname in names(AB_bone)) {
  meta <- AB_bone[[sname]]@meta.data
  AB_bone[[sname]]$car_status <- get_car_status(meta)
  AB_bone[[sname]]$patient    <- sub("_bone_AB$","", sname)
}

summary_bone <- bind_rows(lapply(names(AB_bone), function(s) {
  meta  <- AB_bone[[s]]@meta.data
  n_car <- sum(meta$car_status == "CAR+")
  data.frame(Campione = s, Paziente = unique(meta$patient),
             Totale = nrow(meta), CAR_pos = n_car,
             Perc_CAR = round(100 * n_car / nrow(meta), 1),
             Incluso = ifelse(n_car >= MIN_CAR_CELLS, "SI", "NO"))
}))

AB_bone_valid <- AB_bone[sapply(AB_bone, function(s)
  sum(s@meta.data$car_status == "CAR+") >= MIN_CAR_CELLS)]

cat("Campioni bone validi:", paste(names(AB_bone_valid), collapse = ", "), "\n")

# ======================================================================
# SUMMARY BONE
# ======================================================================
kable(summary_bone, caption = "Campioni bone marrow AB — overview CAR+", align = "c")

# ======================================================================
# CAR MINUS COMPOSITION
# ======================================================================
# Funzione: report composizione cell_type nel gruppo CAR-
report_car_minus_composition <- function(seurat_list, timepoint_label) {
  bind_rows(lapply(names(seurat_list), function(sname) {
    meta <- seurat_list[[sname]]@meta.data
    if (!"cell_type" %in% colnames(meta)) return(NULL)
    meta_minus <- meta[meta$car_status == "CAR-", ]
    if (nrow(meta_minus) == 0) return(NULL)
    ct_counts <- sort(table(meta_minus$cell_type), decreasing = TRUE)
    data.frame(
      Campione       = sname,
      Paziente       = unique(meta$patient),
      Timepoint      = timepoint_label,
      Tipo_cellulare = names(ct_counts),
      N_cellule      = as.integer(ct_counts),
      Perc           = round(100 * as.integer(ct_counts) / nrow(meta_minus), 1),
      Is_T_cell      = names(ct_counts) %in% T_CELL_TYPES
    )
  }))
}

comp_I     <- report_car_minus_composition(I_valid,        "I (Infusion)")
comp_blood <- report_car_minus_composition(AB_blood_valid, "Blood AB")
comp_bone  <- report_car_minus_composition(AB_bone_valid,  "Bone marrow AB")

comp_all <- bind_rows(comp_I, comp_blood, comp_bone)

# Salva tabella
DIR_PP_COMP_TABLES <- file.path(OUT_DIR, "05_per_patient", "comparison", "tables")
dir.create(DIR_PP_COMP_TABLES, recursive = TRUE, showWarnings = FALSE)
write_xlsx(comp_all,
           file.path(DIR_PP_COMP_TABLES, "CAR_minus_composition_before_filter.xlsx"))

cat("Composizione CAR- (prima del filtraggio ai T cells):\n")
cat("Totale righe nel report:", nrow(comp_all), "\n")

# ======================================================================
# CAR MINUS TABLE
# ======================================================================
for (tp in c("I (Infusion)", "Blood AB", "Bone marrow AB")) {
  df_tp <- comp_all[comp_all$Timepoint == tp, ]
  if (nrow(df_tp) == 0) next
  cat("\n### ", tp, "\n\n")
  # Totale T cells vs non-T cells per campione
  summary_tp <- df_tp %>%
    group_by(Campione, Paziente) %>%
    summarise(
      N_tot       = sum(N_cellule),
      N_T_cells   = sum(N_cellule[Is_T_cell]),
      N_non_T     = sum(N_cellule[!Is_T_cell]),
      Perc_T      = round(100 * sum(N_cellule[Is_T_cell]) / sum(N_cellule), 1),
      .groups = "drop"
    )
  print(kable(summary_tp,
              caption = paste("Composizione CAR- —", tp),
              align = "c"))
  cat("\n\n")
  print(
    DT::datatable(
      df_tp[, c("Campione","Paziente","Tipo_cellulare","N_cellule","Perc","Is_T_cell")],
      caption  = paste("Dettaglio tipi cellulari CAR- —", tp),
      filter   = "top",
      options  = list(pageLength = 15, scrollX = TRUE),
      rownames = FALSE
    )
  )
  cat("\n\n")
}

# ======================================================================
# SETUP PER PATIENT DIRS
# ======================================================================
DIR_PP <- file.path(OUT_DIR, "05_per_patient")

for (pt in c("Bo","Ca","Me")) {
  dir.create(file.path(DIR_PP, pt, "figures"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(DIR_PP, pt, "tables"),  recursive = TRUE, showWarnings = FALSE)
}
dir.create(file.path(DIR_PP, "comparison", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DIR_PP, "comparison", "tables"),  recursive = TRUE, showWarnings = FALSE)

# ======================================================================
# PER PATIENT DE
# ======================================================================
# Mappa paziente → lista di oggetti Seurat per timepoint
# Un oggetto Seurat per campione (non merged)
get_pt_obj <- function(seurat_list, patient) {
  matches <- Filter(Negate(is.null), lapply(seurat_list, function(obj) {
    if (patient %in% unique(obj@meta.data$patient)) obj else NULL
  }))
  if (length(matches) == 0) return(NULL)
  matches[[1]]
}

run_de_single <- function(obj, patient, timepoint_label) {
  cat("  DE", patient, "@", timepoint_label, "... ")
  obj <- tryCatch(JoinLayers(obj), error = function(e) obj)

  # Restringe CAR- ai soli T cells
  # CAR+ sono già per definizione T cells; CAR- viene filtrato per evitare confounding
  if ("cell_type" %in% colnames(obj@meta.data)) {
    n_minus_total <- sum(obj@meta.data$car_status == "CAR-")
    keep_cells <- rownames(obj@meta.data)[
      obj@meta.data$car_status == "CAR+" |
      (obj@meta.data$car_status == "CAR-" & obj@meta.data$cell_type %in% T_CELL_TYPES)
    ]
    obj <- subset(obj, cells = keep_cells)
    obj <- tryCatch(JoinLayers(obj), error = function(e) obj)
    n_minus_t <- sum(obj@meta.data$car_status == "CAR-")
    cat(sprintf("(CAR- filtrato: %d T cells / %d totale) ", n_minus_t, n_minus_total))
  } else {
    cat("(cell_type non disponibile — CAR- non filtrato) ")
  }

  Idents(obj) <- obj@meta.data$car_status
  n_plus  <- sum(obj@meta.data$car_status == "CAR+")
  n_minus <- sum(obj@meta.data$car_status == "CAR-")
  cat("CAR+:", n_plus, "CAR-:", n_minus, "\n")
  if (n_plus < 3 || n_minus < 3) { cat("  → skipped\n"); return(NULL) }

  markers <- tryCatch(
    FindMarkers(obj, ident.1 = "CAR+", ident.2 = "CAR-",
                test.use = "wilcox", logfc.threshold = DE_LOG2FC,
                min.pct = 0.1, only.pos = FALSE),
    error = function(e) { cat("  → errore:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(markers)) return(NULL)
  markers$gene      <- rownames(markers)
  markers$patient   <- patient
  markers$timepoint <- timepoint_label
  markers <- markers[markers$p_val_adj < DE_PADJ, ]
  markers <- markers[markers$gene %in% filter_noise(markers$gene), ]
  markers$direction <- ifelse(markers$avg_log2FC > 0, "up", "down")
  markers[order(markers$avg_log2FC, decreasing = TRUE), ]
}

# Mappa campioni disponibili per paziente × timepoint
PATIENT_TP <- list(
  Bo = list(I = get_pt_obj(I_valid, "Bo"),
            Blood = get_pt_obj(AB_blood_valid, "Bo"),
            Bone  = get_pt_obj(AB_bone_valid,  "Bo")),
  Ca = list(I = get_pt_obj(I_valid, "Ca"),
            Blood = get_pt_obj(AB_blood_valid, "Ca"),
            Bone  = NULL),
  Me = list(I = get_pt_obj(I_valid, "Me"),
            Blood = NULL,
            Bone  = get_pt_obj(AB_bone_valid, "Me"))
)

TP_LABELS <- c(I = "I (Infusion)", Blood = "Blood AB", Bone = "Bone marrow AB")

# Esegui DE per ogni combinazione paziente × timepoint
de_pp <- lapply(names(PATIENT_TP), function(pt) {
  cat("\n=== Paziente:", pt, "===\n")
  setNames(lapply(names(TP_LABELS), function(tp) {
    obj <- PATIENT_TP[[pt]][[tp]]
    if (is.null(obj)) return(NULL)
    run_de_single(obj, pt, TP_LABELS[[tp]])
  }), names(TP_LABELS))
})
names(de_pp) <- names(PATIENT_TP)

# Salva tabelle DE per paziente
for (pt in names(de_pp)) {
  for (tp in names(de_pp[[pt]])) {
    df <- de_pp[[pt]][[tp]]
    if (!is.null(df) && nrow(df) > 0) {
      fname <- paste0("DE_", tp, ".xlsx")
      write_xlsx(df, file.path(DIR_PP, pt, "tables", fname))
    }
  }
}
cat("\nTabelle DE per paziente salvate in:", DIR_PP, "\n")

# ======================================================================
# PER PATIENT ENRICHMENT
# ======================================================================
# GO e KEGG per ogni paziente × timepoint
enrich_pp <- lapply(names(de_pp), function(pt) {
  cat("\n=== Enrichment paziente:", pt, "===\n")
  setNames(lapply(names(de_pp[[pt]]), function(tp) {
    df <- de_pp[[pt]][[tp]]
    genes <- if (!is.null(df) && nrow(df) >= 10) df$gene else NULL
    if (is.null(genes)) { cat(" ", tp, ": skip (DE < 10)\n"); return(list(go=NULL,kegg=NULL)) }
    cat(" ", tp, ": GO... ")
    go   <- run_go_enrichment(genes, label = paste(pt, tp))
    cat(if (!is.null(go)) nrow(as.data.frame(go)) else 0, "termini | KEGG... ")
    kegg <- run_kegg_enrichment(genes, label = paste(pt, tp))
    cat(if (!is.null(kegg)) nrow(as.data.frame(kegg)) else 0, "pathways\n")

    # Salva tabelle
    if (!is.null(go))
      write_xlsx(as.data.frame(go),
                 file.path(DIR_PP, pt, "tables", paste0("GO_BP_", tp, ".xlsx")))
    if (!is.null(kegg))
      write_xlsx(as.data.frame(kegg),
                 file.path(DIR_PP, pt, "tables", paste0("KEGG_", tp, ".xlsx")))

    list(go = go, kegg = kegg)
  }), names(de_pp[[pt]]))
})
names(enrich_pp) <- names(de_pp)

# ======================================================================
# VOLCANO PER PATIENT
# ======================================================================
PATIENT_COLORS_FULL <- c(Bo = "#E64B35", Ca = "#4DBBD5", Me = "#00A087")

for (pt in names(de_pp)) {
  plots <- lapply(names(TP_LABELS), function(tp) {
    df    <- de_pp[[pt]][[tp]]
    color <- PATIENT_COLORS_FULL[pt]
    label <- TP_LABELS[[tp]]
    if (is.null(df) || nrow(df) == 0)
      return(ggplot() +
               labs(title = paste(pt, "—", label, "— n.d.")) +
               theme_void(base_size = 10))
    df$neg_log10 <- -log10(df$p_val_adj + 1e-300)
    df$direction <- ifelse(df$avg_log2FC > 0, "up", "down")
    n_up   <- sum(df$direction == "up")
    n_down <- sum(df$direction == "down")
    # label top 10 up + top 10 down by abs(log2FC)
    top_up   <- head(df[df$direction == "up",   ], 10)
    top_down <- head(df[df$direction == "down",  ], 10)
    top_lab  <- rbind(top_up, top_down)
    dir_colors <- c(up = color, down = "steelblue")
    ggplot(df, aes(avg_log2FC, neg_log10, color = direction)) +
      geom_point(alpha = 0.5, size = 1.2) +
      scale_color_manual(values = dir_colors, guide = "none") +
      geom_text_repel(data = top_lab, aes(label = gene), size = 2.5,
                      max.overlaps = 20, color = "black") +
      geom_vline(xintercept = c(-DE_LOG2FC, DE_LOG2FC),
                 linetype = "dashed", color = "grey60") +
      geom_hline(yintercept = -log10(DE_PADJ),
                 linetype = "dashed", color = "grey60") +
      geom_vline(xintercept = 0, color = "grey40", linewidth = 0.3) +
      labs(title    = paste(pt, "—", label),
           subtitle = paste0(n_up, " up (CAR+)  |  ", n_down, " down (CAR−)"),
           x = "avg log2FC  (positive = higher in CAR+)", y = "-log10(adj p-val)") +
      theme_classic(base_size = 9)
  })
  p_vol <- wrap_plots(plots, nrow = 1)
  ggsave(file.path(DIR_PP, pt, "figures", "volcano_timepoints.png"),
         plot = p_vol, width = 13, height = 5, dpi = 150)
  print(p_vol)
  cat("\n")
}

# ======================================================================
# GO PER PATIENT
# ======================================================================
for (pt in names(enrich_pp)) {
  for (tp in names(enrich_pp[[pt]])) {
    go_res <- enrich_pp[[pt]][[tp]]$go
    if (is.null(go_res)) next
    p <- dotplot(go_res, showCategory = 15,
                 title = paste("GO BP —", pt, "—", TP_LABELS[[tp]])) +
      theme(axis.text.y = element_text(size = 8))
    ggsave(file.path(DIR_PP, pt, "figures", paste0("dotplot_GO_", tp, ".png")),
           plot = p, width = 10, height = 7, dpi = 150)
    print(p)
    cat("\n")
  }
}

# ======================================================================
# KEGG PER PATIENT
# ======================================================================
for (pt in names(enrich_pp)) {
  for (tp in names(enrich_pp[[pt]])) {
    kegg_res <- enrich_pp[[pt]][[tp]]$kegg
    if (is.null(kegg_res)) next
    p <- dotplot(kegg_res, showCategory = 12,
                 title = paste("KEGG —", pt, "—", TP_LABELS[[tp]])) +
      theme(axis.text.y = element_text(size = 8))
    ggsave(file.path(DIR_PP, pt, "figures", paste0("dotplot_KEGG_", tp, ".png")),
           plot = p, width = 10, height = 6, dpi = 150)
    print(p)
    cat("\n")
  }
}

# ======================================================================
# DT DE PER PATIENT
# ======================================================================
de_pp_all <- bind_rows(lapply(names(de_pp), function(pt) {
  bind_rows(lapply(de_pp[[pt]], function(df) df))
}))

if (nrow(de_pp_all) > 0) {
  DT::datatable(
    de_pp_all[, c("patient","timepoint","gene","avg_log2FC","pct.1","pct.2","p_val_adj")],
    caption  = "DE markers CAR+ vs CAR- — analisi per paziente",
    filter   = "top",
    options  = list(pageLength = 20, scrollX = TRUE),
    rownames = FALSE
  ) %>%
    DT::formatRound(c("avg_log2FC","pct.1","pct.2","p_val_adj"), digits = 4)
}

# ======================================================================
# UPSET PATIENTS BY TP
# ======================================================================
for (tp in names(TP_LABELS)) {
  gene_sets_pt <- Filter(function(x) length(x) > 0,
    setNames(lapply(names(de_pp), function(pt) {
      df <- de_pp[[pt]][[tp]]
      if (!is.null(df) && nrow(df) > 0) df$gene else character(0)
    }), names(de_pp))
  )
  if (length(gene_sets_pt) < 2) next

  fname <- paste0("upset_patients_", tp, ".png")
  png(file.path(DIR_PP, "comparison", "figures", fname),
      width = 9, height = 5, units = "in", res = 150)
  upset(fromList(gene_sets_pt), order.by = "freq",
        sets = names(gene_sets_pt),
        sets.bar.color = unname(PATIENT_COLORS_FULL[names(gene_sets_pt)]),
        main.bar.color = "#555555", text.scale = 1.3,
        mainbar.y.label = "N. geni DE condivisi",
        sets.x.label   = paste("DE in", TP_LABELS[[tp]]))
  dev.off()

  cat("\n###", TP_LABELS[[tp]], "\n")
  upset(fromList(gene_sets_pt), order.by = "freq",
        sets = names(gene_sets_pt),
        sets.bar.color = unname(PATIENT_COLORS_FULL[names(gene_sets_pt)]),
        main.bar.color = "#555555", text.scale = 1.3,
        mainbar.y.label = "N. geni DE condivisi",
        sets.x.label   = paste("DE in", TP_LABELS[[tp]]))
}

# ======================================================================
# GO HEATMAP PATIENTS
# ======================================================================
go_all_pp <- bind_rows(lapply(names(enrich_pp), function(pt) {
  bind_rows(lapply(names(enrich_pp[[pt]]), function(tp) {
    res <- enrich_pp[[pt]][[tp]]$go
    if (is.null(res)) return(NULL)
    df <- as.data.frame(res)
    if (nrow(df) == 0) return(NULL)
    df %>% arrange(p.adjust) %>% head(10) %>%
      mutate(patient = pt, tp_label = TP_LABELS[[tp]],
             group = paste(pt, TP_LABELS[[tp]], sep = " — "),
             neg_log10 = -log10(p.adjust))
  }))
}))

if (!is.null(go_all_pp) && nrow(go_all_pp) > 0) {
  top_terms_pp <- go_all_pp %>% count(Description) %>%
    arrange(desc(n)) %>% head(30) %>% pull(Description)

  go_cmp_pp <- go_all_pp %>%
    filter(Description %in% top_terms_pp) %>%
    mutate(patient   = factor(patient, levels = c("Bo","Ca","Me")),
           tp_label  = factor(tp_label,
                        levels = c("I (Infusion)","Blood AB","Bone marrow AB")))

  p_go_hm <- ggplot(go_cmp_pp,
                    aes(x = tp_label, y = Description,
                        size = Count, color = neg_log10)) +
    geom_point(alpha = 0.85) +
    facet_grid(. ~ patient, scales = "free_x", space = "free_x") +
    scale_color_gradient(low = "#4DBBD5", high = "#E64B35",
                         name = "-log10\n(adj p-val)") +
    scale_size_continuous(range = c(1, 8), name = "Count") +
    labs(title    = "GO BP — confronto pazienti × timepoint",
         subtitle = "Top 30 termini arricchiti in almeno una combinazione",
         x = NULL, y = NULL) +
    theme_classic(base_size = 9) +
    theme(axis.text.x  = element_text(angle = 30, hjust = 1, size = 8),
          axis.text.y  = element_text(size = 7),
          strip.text   = element_text(face = "bold", size = 10),
          strip.background = element_rect(fill = "#F0F0F0"))

  ggsave(file.path(DIR_PP, "comparison", "figures", "GO_heatmap_patients_timepoints.png"),
         plot = p_go_hm, width = 12, height = 9, dpi = 150)
  p_go_hm
}

# ======================================================================
# KEGG HEATMAP PATIENTS
# ======================================================================
kegg_all_pp <- bind_rows(lapply(names(enrich_pp), function(pt) {
  bind_rows(lapply(names(enrich_pp[[pt]]), function(tp) {
    res <- enrich_pp[[pt]][[tp]]$kegg
    if (is.null(res)) return(NULL)
    df <- as.data.frame(res)
    if (nrow(df) == 0) return(NULL)
    df %>% arrange(p.adjust) %>% head(10) %>%
      mutate(patient = pt, tp_label = TP_LABELS[[tp]],
             neg_log10 = -log10(p.adjust))
  }))
}))

if (!is.null(kegg_all_pp) && nrow(kegg_all_pp) > 0) {
  top_kegg_pp <- kegg_all_pp %>% count(Description) %>%
    arrange(desc(n)) %>% head(25) %>% pull(Description)

  kegg_cmp_pp <- kegg_all_pp %>%
    filter(Description %in% top_kegg_pp) %>%
    mutate(patient  = factor(patient, levels = c("Bo","Ca","Me")),
           tp_label = factor(tp_label,
                       levels = c("I (Infusion)","Blood AB","Bone marrow AB")))

  p_kegg_hm <- ggplot(kegg_cmp_pp,
                      aes(x = tp_label, y = Description,
                          size = Count, color = neg_log10)) +
    geom_point(alpha = 0.85) +
    facet_grid(. ~ patient, scales = "free_x", space = "free_x") +
    scale_color_gradient(low = "#4DBBD5", high = "#E64B35",
                         name = "-log10\n(adj p-val)") +
    scale_size_continuous(range = c(1, 8), name = "Count") +
    labs(title    = "KEGG — confronto pazienti × timepoint",
         subtitle = "Top 25 pathway in almeno una combinazione",
         x = NULL, y = NULL) +
    theme_classic(base_size = 9) +
    theme(axis.text.x  = element_text(angle = 30, hjust = 1, size = 8),
          axis.text.y  = element_text(size = 7),
          strip.text   = element_text(face = "bold", size = 10),
          strip.background = element_rect(fill = "#F0F0F0"))

  ggsave(file.path(DIR_PP, "comparison", "figures", "KEGG_heatmap_patients_timepoints.png"),
         plot = p_kegg_hm, width = 12, height = 8, dpi = 150)
  p_kegg_hm
}

# ======================================================================
# BARPLOT N DE PER PATIENT
# ======================================================================
n_de_df <- bind_rows(lapply(names(de_pp), function(pt) {
  bind_rows(lapply(names(TP_LABELS), function(tp) {
    df <- de_pp[[pt]][[tp]]
    n_up   <- if (!is.null(df)) sum(df$avg_log2FC > 0) else 0
    n_down <- if (!is.null(df)) sum(df$avg_log2FC < 0) else 0
    bind_rows(
      data.frame(patient=pt, timepoint=TP_LABELS[[tp]], n_DE=n_up,   direction="up (CAR+)"),
      data.frame(patient=pt, timepoint=TP_LABELS[[tp]], n_DE=n_down, direction="down (CAR-)")
    )
  }))
})) %>%
  mutate(patient   = factor(patient, levels = c("Bo","Ca","Me")),
         timepoint = factor(timepoint,
                     levels = c("I (Infusion)","Blood AB","Bone marrow AB")),
         direction = factor(direction, levels = c("up (CAR+)", "down (CAR-)")))

p_n_de <- ggplot(n_de_df, aes(x = timepoint, y = n_DE, fill = patient, alpha = direction)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = ifelse(n_DE > 0, n_DE, "")),
            position = position_dodge(width = 0.9), vjust = -0.4, size = 3) +
  scale_fill_manual(values = PATIENT_COLORS_FULL) +
  scale_alpha_manual(values = c("up (CAR+)" = 0.9, "down (CAR-)" = 0.45)) +
  labs(title = "N. geni DE (CAR+ vs CAR- T cells) per paziente e timepoint",
       subtitle = paste0("Wilcoxon, |log2FC| > ", DE_LOG2FC, ", FDR < ", DE_PADJ,
                         " | Esclusi: MT, ribosomiali, MALAT1/NEAT1"),
       x = NULL, y = "N. geni DE", fill = "Paziente", alpha = "Direzione") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

ggsave(file.path(DIR_PP, "comparison", "figures", "n_DE_per_patient_timepoint.png"),
       plot = p_n_de, width = 7, height = 4, dpi = 150)
p_n_de

# ======================================================================
# DT ENRICH PER PATIENT
# ======================================================================
enrich_summary <- bind_rows(lapply(names(enrich_pp), function(pt) {
  bind_rows(lapply(names(enrich_pp[[pt]]), function(tp) {
    go_n   <- if (!is.null(enrich_pp[[pt]][[tp]]$go))
                nrow(as.data.frame(enrich_pp[[pt]][[tp]]$go)) else 0
    kegg_n <- if (!is.null(enrich_pp[[pt]][[tp]]$kegg))
                nrow(as.data.frame(enrich_pp[[pt]][[tp]]$kegg)) else 0
    de_n   <- if (!is.null(de_pp[[pt]][[tp]])) nrow(de_pp[[pt]][[tp]]) else 0
    data.frame(Paziente = pt, Timepoint = TP_LABELS[[tp]],
               DE_genes = de_n, GO_terms = go_n, KEGG_pathways = kegg_n)
  }))
}))

write_xlsx(enrich_summary,
           file.path(DIR_PP, "comparison", "tables", "enrichment_summary_per_patient.xlsx"))

kable(enrich_summary, align = "c",
      caption = "Riepilogo analisi per paziente × timepoint")

# ======================================================================
# SUMMARY SECTION
# ======================================================================
cat("=== RIEPILOGO ANALISI P3 ===\n\n")

cat("APPROCCIO: DE markers CAR+ vs CAR- (solo T cells), per paziente × timepoint\n")
cat("           Wilcoxon, log2FC >", DE_LOG2FC, ", FDR <", DE_PADJ, "\n")
cat("           Esclusi: MT, ribosomiali, MALAT1/NEAT1\n\n")

cat("CAMPIONI PER TIMEPOINT:\n")
cat("  I:             ", paste(sapply(I_valid,        function(o) unique(o@meta.data$patient)), collapse=", "), "\n")
cat("  Blood AB:      ", paste(sapply(AB_blood_valid, function(o) unique(o@meta.data$patient)), collapse=", "), "\n")
cat("  Bone marrow AB:", paste(sapply(AB_bone_valid,  function(o) unique(o@meta.data$patient)), collapse=", "), "\n\n")

cat("RISULTATI DE PER PAZIENTE:\n")
for (pt in names(de_pp)) {
  cat(" ", pt, ":\n")
  for (tp in names(TP_LABELS)) {
    df <- de_pp[[pt]][[tp]]
    n  <- if (!is.null(df)) nrow(df) else 0
    cat(sprintf("    %-20s  %d geni DE\n", TP_LABELS[[tp]], n))
  }
}

cat("\nFILE SALVATI IN:", file.path(OUT_DIR, "05_per_patient"), "\n")
cat("  comparison/tables/ ←", length(list.files(file.path(OUT_DIR,"05_per_patient","comparison","tables"))), "tabelle\n")
