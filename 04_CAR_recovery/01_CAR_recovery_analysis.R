# ============================================================
#  4_CAR_recovery_analysis.R
#
#  Analisi del bottleneck e recovery delle cellule CAR+
#  perse durante il QC dei campioni I.
#
#  Implementa la proposta del Prof.:
#    - Trova clonotipi CAR+ in A/B non presenti in I annotato
#    - Quantifica trasparentemente la perdita
#    - Non richiede riallineamento
#
#  Prerequisiti:
#    - seurat_list originale (raw, pre-filtro QC)
#    - all_samples_annotated_COMPLETE.rds (già disponibile)
#
#  Output in <out_dir>/CAR_recovery/:
#    QC audit:
#      XX_QC_KEPT_vs_LOST.png
#      QC_audit_CAR_in_removed_cells.xlsx
#    Cross-referencing:
#      XX_CAR_recovery_barplot.png
#      XX_UMAP_CAR_recovery.png
#      CAR_cross_reference_I_vs_AB.xlsx
# ============================================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(scales)
library(ggrepel)
library(openxlsx)

# ── UNICO PUNTO DA MODIFICARE ────────────────────────────────
seurat_list_path   <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/CART/Data/seurat_obj_list/seurat_samples_sctype_azimuth_pbmc_bonemarrow_clonalvdj_CAR.rds"
annotated_rds_path <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/CART/Code/Seurat_analysis/2_annotation/all_samples_annotated_COMPLETE.rds"
out_dir            <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/CART/Code/Seurat_analysis/4_CAR_recovery/"
# ─────────────────────────────────────────────────────────────

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

section <- function(title)
  cat(paste0("\n", strrep("=", 65), "\n  ", title,
             "\n", strrep("=", 65), "\n"))

# ── MAPPA PAZIENTE → CAMPIONI ─────────────────────────────────
# ADATTA i nomi in base a names(all_samples)
patient_map <- list(
  Ca = list(I = "Ca_bone_I",  AB = c("Ca_blood_AB", "Ca_bone_AB")),
  Bo = list(I = "Bo_bone_I",  AB = c("Bo_blood_AB", "Bo_bone_AB")),
  Me = list(I = "Me_bone_I",  AB = c("Me_bone_AB"))
)

# ============================================================
# CARICAMENTO
# ============================================================
section("Caricamento oggetti")

cat("Caricamento seurat_list (raw)...\n")
seurat_list <- readRDS(seurat_list_path)

cat("Caricamento all_samples annotato...\n")
all_samples <- readRDS(annotated_rds_path)

cat("\nCampioni disponibili in all_samples:\n")
for (nm in names(all_samples))
  cat(sprintf("  %-20s | %5d celle | cell_type: %s | IS_CAR: %s\n",
              nm, ncol(all_samples[[nm]]),
              if ("cell_type" %in% colnames(all_samples[[nm]]@meta.data)) "OK" else "MANCANTE",
              if ("IS_CAR_ALLIN_scREP" %in% colnames(all_samples[[nm]]@meta.data)) "OK" else "MANCANTE"))

# Oggetti I raw (pre-filtro QC)
raw_I_map <- list(
  Bo_bone_I = seurat_list$Bo_samples_blood$I,
  Ca_bone_I = seurat_list$Ca_samples_blood$I,
  Me_bone_I = seurat_list$Me_samples_bone$I
)

# ============================================================
# STEP 1-2: QC AUDIT
# Identifica cellule rimosse dal QC in I e verifica se
# avevano segnale CAR
# ============================================================
section("STEP 1-2 | QC Audit campioni I")

audit_I_sample <- function(nm, raw_obj, annotated_obj) {

  cat(paste0("\n── ", nm, " ", strrep("─", 40), "\n"))

  if (is.null(raw_obj)) {
    cat("[SKIP] Oggetto raw non trovato\n")
    return(NULL)
  }
  if (is.null(annotated_obj)) {
    cat("[SKIP] Oggetto annotato non trovato\n")
    return(NULL)
  }

  # Barcode delle due versioni
  bc_raw  <- Cells(raw_obj)
  bc_kept <- Cells(annotated_obj)
  bc_lost <- setdiff(bc_raw, bc_kept)

  cat(sprintf("  Raw totale: %d\n  Mantenuti:  %d\n  Rimossi:    %d (%.1f%%)\n",
              length(bc_raw), length(bc_kept), length(bc_lost),
              100 * length(bc_lost) / length(bc_raw)))

  if (length(bc_lost) == 0) {
    cat("  [INFO] Nessuna cellula rimossa. Niente da verificare.\n")
    return(data.frame(sample = nm, n_raw = length(bc_raw),
                      n_kept = length(bc_kept), n_lost = 0))
  }

  meta_raw <- raw_obj@meta.data

  # Colonne CAR presenti nel raw
  car_cols <- intersect(
    c("CAR", "derived_CAR", "IS_CAR_ALLIN_scREP",
      "IS_CONSERVED_scRepertoire", "clonal", "barcode_screpertoire"),
    colnames(meta_raw)
  )
  cat("  Colonne CAR nel raw:", paste(car_cols, collapse = ", "), "\n")

  # Aggiunge colonna stato
  meta_raw$qc_status <- ifelse(rownames(meta_raw) %in% bc_kept, "KEPT", "LOST")

  # ── Violinplot QC ──────────────────────────────────────────
  qc_vars <- intersect(c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                       colnames(meta_raw))

  violin_plots <- lapply(qc_vars, function(var) {
    med_kept <- median(meta_raw[[var]][meta_raw$qc_status == "KEPT"], na.rm = TRUE)
    med_lost <- median(meta_raw[[var]][meta_raw$qc_status == "LOST"], na.rm = TRUE)

    ggplot(meta_raw, aes(x = qc_status, y = .data[[var]], fill = qc_status)) +
      geom_violin(alpha = 0.7, trim = FALSE) +
      geom_boxplot(width = 0.12, fill = "white", outlier.size = 0.4,
                   outlier.alpha = 0.3) +
      scale_fill_manual(values = c(KEPT = "#58a6ff", LOST = "#f85149"),
                        guide = "none") +
      annotate("text", x = 1, y = Inf,
               label = paste0("med=", round(med_kept, 1)),
               vjust = 2, size = 3, color = "#58a6ff") +
      annotate("text", x = 2, y = Inf,
               label = paste0("med=", round(med_lost, 1)),
               vjust = 2, size = 3, color = "#f85149") +
      labs(title = var, x = NULL, y = NULL) +
      theme_classic(base_size = 11) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5))
  })

  # Aggiunge filtri come linee di riferimento
  if ("nFeature_RNA" %in% qc_vars) {
    violin_plots[[which(qc_vars == "nFeature_RNA")]] <-
      violin_plots[[which(qc_vars == "nFeature_RNA")]] +
      geom_hline(yintercept = 800, linetype = "dashed",
                 color = "orange", linewidth = 0.8) +
      annotate("text", x = 0.5, y = 820, label = "cutoff 800",
               size = 2.8, color = "orange", hjust = 0)
  }
  if ("percent.mt" %in% qc_vars) {
    violin_plots[[which(qc_vars == "percent.mt")]] <-
      violin_plots[[which(qc_vars == "percent.mt")]] +
      geom_hline(yintercept = 7, linetype = "dashed",
                 color = "orange", linewidth = 0.8) +
      annotate("text", x = 0.5, y = 7.3, label = "cutoff 7%",
               size = 2.8, color = "orange", hjust = 0)
  }

  p_violin <- wrap_plots(violin_plots, nrow = 1) +
    plot_annotation(
      title    = paste0(nm, " – QC: KEPT (blu) vs LOST (rosso)"),
      subtitle = paste0("Rimossi: ", length(bc_lost), " (",
                        round(100 * length(bc_lost) / length(bc_raw), 1),
                        "%) | Filtri: nFeature_RNA > 800 & percent.mt < 7"),
      theme = theme(
        plot.title    = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(color = "gray40", size = 10)))

  ggsave(paste0(out_dir, nm, "_QC_KEPT_vs_LOST.png"),
         p_violin, width = 13, height = 5, dpi = 300, bg = "white")
  cat(paste0("  → ", nm, "_QC_KEPT_vs_LOST.png\n"))

  # ── CAR nelle cellule rimosse ──────────────────────────────
  meta_lost <- meta_raw[rownames(meta_raw) %in% bc_lost, , drop = FALSE]
  car_audit <- data.frame(
    sample = nm,
    n_raw  = length(bc_raw),
    n_kept = length(bc_kept),
    n_lost = length(bc_lost),
    pct_lost = round(100 * length(bc_lost) / length(bc_raw), 1),
    stringsAsFactors = FALSE
  )

  for (col in car_cols) {
    vals  <- as.character(meta_lost[[col]])
    n_yes <- sum(vals %in% c("YES", "TRUE", "1", "yes", "true"), na.rm = TRUE)
    car_audit[[paste0("n_", col)]] <- n_yes
    if (n_yes > 0)
      cat(sprintf("  ⚠ CAR nelle RIMOSSE: col '%s' = %d cellule\n", col, n_yes))
  }

  # Controlla anche espressione genica del costrutto CAR
  if ("CAR" %in% rownames(raw_obj)) {
    tryCatch({
      car_counts <- GetAssayData(raw_obj[, bc_lost], slot = "counts")["CAR", ]
      n_expr     <- sum(car_counts > 0, na.rm = TRUE)
      car_audit$n_CAR_expr_counts_positive <- n_expr
      cat(sprintf("  ⚠ Cellule rimosse con counts CAR > 0: %d\n", n_expr))
    }, error = function(e)
      cat(paste0("  [WARN] Impossibile accedere ai counts CAR: ", e$message, "\n")))
  }

  car_audit
}

# Esegui audit per tutti i campioni I
audit_results <- list()
for (nm in names(raw_I_map)) {
  ann_obj <- all_samples[[nm]]
  if (is.null(ann_obj)) {
    cat(paste0("[SKIP] ", nm, " non trovato in all_samples. Verifica il nome.\n"))
    next
  }
  audit_results[[nm]] <- audit_I_sample(nm, raw_I_map[[nm]], ann_obj)
}

# Tabella riepilogativa QC audit
audit_tbl <- bind_rows(Filter(Negate(is.null), audit_results))
cat("\nRIEPILOGO QC AUDIT:\n")
print(audit_tbl, row.names = FALSE)

write.xlsx(audit_tbl,
           paste0(out_dir, "QC_audit_CAR_in_removed_cells.xlsx"),
           overwrite = TRUE)
cat(paste0("→ QC_audit_CAR_in_removed_cells.xlsx\n"))

# ============================================================
# STEP 3-4: CROSS-REFERENCING CLONOTIPI I ↔ A/B
# Proposta del Prof.: trova clonotipi CAR in A/B non in I
# ============================================================
section("STEP 3-4 | Cross-referencing clonotipi I ↔ A/B")

# Funzione per estrarre info clonotipo da un oggetto Seurat
get_clonotype_df <- function(obj, source_label) {
  meta   <- obj@meta.data

  # Colonna identificatore del clonotipo – priorità: clonal > barcode_screpertoire
  id_col <- if ("clonal" %in% colnames(meta)) "clonal" else
            if ("barcode_screpertoire" %in% colnames(meta)) "barcode_screpertoire" else NULL

  if (is.null(id_col)) {
    warning(paste0("Nessuna colonna clonotipo in ", source_label))
    return(NULL)
  }

  # Colonna CAR
  car_col <- if ("IS_CAR_ALLIN_scREP" %in% colnames(meta)) "IS_CAR_ALLIN_scREP" else
             if ("CAR" %in% colnames(meta)) "CAR" else NULL
  car_vals <- if (!is.null(car_col)) as.character(meta[[car_col]]) else
              rep(NA_character_, nrow(meta))

  data.frame(
    barcode       = rownames(meta),
    clonotype     = as.character(meta[[id_col]]),
    IS_CAR        = car_vals,
    cell_type     = if ("cell_type" %in% colnames(meta)) as.character(meta$cell_type) else NA_character_,
    source_sample = source_label,
    nFeature      = if ("nFeature_RNA" %in% colnames(meta)) meta$nFeature_RNA else NA_real_,
    percent_mt    = if ("percent.mt"   %in% colnames(meta)) meta$percent.mt   else NA_real_,
    stringsAsFactors = FALSE
  ) %>%
    filter(!is.na(clonotype), clonotype != "", clonotype != "NA",
           !is.na(IS_CAR))
}

cross_reference_clonotypes <- function(patient_id, pm, all_samples, out_dir) {

  cat(paste0("\n── Paziente: ", patient_id,
             " ──────────────────────────────────\n"))

  obj_I  <- all_samples[[pm$I]]
  nms_AB <- pm$AB[pm$AB %in% names(all_samples)]

  if (is.null(obj_I)) {
    cat("[SKIP] campione I non trovato. Verifica patient_map.\n")
    return(NULL)
  }
  if (length(nms_AB) == 0) {
    cat("[SKIP] nessun campione AB disponibile.\n")
    return(NULL)
  }

  df_I  <- get_clonotype_df(obj_I, pm$I)
  if (is.null(df_I) || nrow(df_I) == 0) {
    cat("[SKIP] nessuna info clonotipo in I.\n")
    return(NULL)
  }

  df_AB <- bind_rows(lapply(nms_AB, function(nm)
    get_clonotype_df(all_samples[[nm]], nm)))
  if (is.null(df_AB) || nrow(df_AB) == 0) {
    cat("[SKIP] nessuna info clonotipo in AB.\n")
    return(NULL)
  }

  # Clonotipi CAR+ nei due insiemi
  car_clono_I  <- unique(df_I$clonotype[df_I$IS_CAR == "YES"])
  car_clono_AB <- unique(df_AB$clonotype[df_AB$IS_CAR == "YES"])

  only_in_AB <- setdiff(car_clono_AB, car_clono_I)
  in_both    <- intersect(car_clono_AB, car_clono_I)
  only_in_I  <- setdiff(car_clono_I,  car_clono_AB)

  cat(sprintf("  Clonotipi CAR+:\n"))
  cat(sprintf("    In I annotato:              %d\n", length(car_clono_I)))
  cat(sprintf("    In A/B:                     %d\n", length(car_clono_AB)))
  cat(sprintf("    Confermati (in entrambi):   %d\n", length(in_both)))
  cat(sprintf("    Solo in A/B (PERSI da I):   %d  ← potenzialmente filtrati\n",
              length(only_in_AB)))
  cat(sprintf("    Solo in I (non trovati AB): %d\n", length(only_in_I)))

  # Tabella Venn per paziente
  venn_tbl <- data.frame(
    patient             = patient_id,
    CAR_clono_in_I      = length(car_clono_I),
    CAR_clono_in_AB     = length(car_clono_AB),
    confirmed_in_both   = length(in_both),
    only_in_AB_lost_I   = length(only_in_AB),
    only_in_I_notAB     = length(only_in_I),
    cells_recovered_AB  = sum(df_AB$IS_CAR == "YES" &
                                df_AB$clonotype %in% only_in_AB)
  )

  # ── Tag metadata negli oggetti AB ───────────────────────────
  recovery_lookup <- c(
    setNames(rep("CAR_confirmed_in_both",           length(in_both)),    in_both),
    setNames(rep("CAR_present_in_AB_lost_from_I",   length(only_in_AB)), only_in_AB)
  )

  updated_AB <- list()
  for (nm in nms_AB) {
    obj    <- all_samples[[nm]]
    meta   <- obj@meta.data
    id_col <- if ("clonal" %in% colnames(meta)) "clonal" else "barcode_screpertoire"
    clono  <- as.character(meta[[id_col]])
    obj$CAR_recovery_status <- ifelse(
      clono %in% names(recovery_lookup),
      recovery_lookup[clono],
      "not_CAR"
    )
    updated_AB[[nm]] <- obj
    n_rec  <- sum(obj$CAR_recovery_status == "CAR_present_in_AB_lost_from_I", na.rm = TRUE)
    n_conf <- sum(obj$CAR_recovery_status == "CAR_confirmed_in_both",          na.rm = TRUE)
    cat(sprintf("  %s: %d confermati | %d recuperati\n", nm, n_conf, n_rec))
  }

  # ── Barplot breakdown ────────────────────────────────────────
  df_bar <- bind_rows(lapply(nms_AB, function(nm) {
    obj <- updated_AB[[nm]]
    as.data.frame(table(CAR_status = obj$CAR_recovery_status)) %>%
      mutate(sample = nm,
             pct    = round(100 * Freq / ncol(obj), 2))
  }))

  status_colors <- c(
    CAR_confirmed_in_both          = "#3fb950",
    CAR_present_in_AB_lost_from_I  = "#d29922",
    not_CAR                        = "#30363d"
  )
  status_labels <- c(
    CAR_confirmed_in_both          = "CAR confermata (in I e AB)",
    CAR_present_in_AB_lost_from_I  = "CAR recuperata (AB, persa da I)",
    not_CAR                        = "Non-CAR"
  )

  p_bar <- ggplot(df_bar %>% filter(CAR_status != "not_CAR"),
                  aes(x = sample, y = Freq, fill = CAR_status)) +
    geom_col(position = "stack", width = 0.65) +
    geom_text(aes(label = Freq), position = position_stack(vjust = 0.5),
              color = "white", size = 3.5, fontface = "bold") +
    scale_fill_manual(values = status_colors, labels = status_labels,
                      name = NULL) +
    labs(title    = paste0(patient_id, " – CAR recovery status in A/B"),
         subtitle = paste0("Arancione = clonotipi CAR in A/B non trovati in I annotato (n=",
                           length(only_in_AB), " clonotipi, ",
                           sum(df_AB$IS_CAR == "YES" & df_AB$clonotype %in% only_in_AB),
                           " cellule)"),
         x = NULL, y = "n cellule") +
    theme_classic(base_size = 12) +
    theme(legend.position = "bottom",
          axis.text.x     = element_text(angle = 30, hjust = 1),
          plot.title       = element_text(face = "bold", size = 13),
          plot.subtitle    = element_text(size = 10, color = "gray40"))

  ggsave(paste0(out_dir, patient_id, "_CAR_recovery_barplot.png"),
         p_bar, width = 9, height = 6, dpi = 300, bg = "white")
  cat(paste0("  → ", patient_id, "_CAR_recovery_barplot.png\n"))

  # Dettaglio cellule recuperate
  df_detail <- bind_rows(
    df_AB %>% filter(IS_CAR == "YES", clonotype %in% only_in_AB) %>%
      mutate(status = "CAR_present_in_AB_lost_from_I"),
    df_AB %>% filter(IS_CAR == "YES", clonotype %in% in_both) %>%
      mutate(status = "CAR_confirmed_in_both")
  )

  list(venn = venn_tbl, detail = df_detail, updated_AB = updated_AB)
}

# Esegui cross-referencing per tutti i pazienti
cross_results <- list()
for (pid in names(patient_map)) {
  cross_results[[pid]] <- cross_reference_clonotypes(
    pid, patient_map[[pid]], all_samples, out_dir
  )
}

# Tabella riepilogativa
summary_venn <- bind_rows(lapply(cross_results, function(r)
  if (!is.null(r)) r$venn))

cat("\nRIEPILOGO CROSS-REFERENCING:\n")
print(summary_venn, row.names = FALSE)

# Salva Excel
wb <- createWorkbook()
addWorksheet(wb, "Riepilogo")
writeData(wb, "Riepilogo", summary_venn)
for (pid in names(cross_results)) {
  r <- cross_results[[pid]]
  if (!is.null(r) && !is.null(r$detail) && nrow(r$detail) > 0) {
    sheet_nm <- substr(pid, 1, 31)
    addWorksheet(wb, sheet_nm)
    writeData(wb, sheet_nm, r$detail)
  }
}
saveWorkbook(wb, paste0(out_dir, "CAR_cross_reference_I_vs_AB.xlsx"),
             overwrite = TRUE)
cat("→ CAR_cross_reference_I_vs_AB.xlsx\n")

# ============================================================
# STEP 5: UMAP OVERLAY CON CAR RECOVERY STATUS
# ============================================================
section("STEP 5 | UMAP recovery status")

plot_recovery_umap <- function(obj, sample_name, out_dir) {

  if (!"CAR_recovery_status" %in% colnames(obj@meta.data)) {
    cat(paste0("[SKIP] ", sample_name, ": CAR_recovery_status non trovata\n"))
    return(NULL)
  }

  # Trova UMAP
  umap_key <- NULL
  for (nm in c("umap","wnn.umap","umap.harmony","umap2","RNA.umap","integrated.umap"))
    if (nm %in% names(obj@reductions)) { umap_key <- nm; break }
  if (is.null(umap_key)) {
    cat(paste0("[SKIP] ", sample_name, ": nessuna riduzione UMAP\n"))
    return(NULL)
  }

  coords <- as.data.frame(Embeddings(obj, umap_key)[, 1:2])
  colnames(coords) <- c("UMAP1", "UMAP2")
  coords$status    <- as.character(obj$CAR_recovery_status)
  coords$cell_type <- as.character(obj$cell_type)

  # Ordine plot: not_CAR sotto, CAR sopra
  coords <- coords[order(match(coords$status,
    c("not_CAR","CAR_confirmed_in_both","CAR_present_in_AB_lost_from_I"))), ]

  status_colors <- c(
    not_CAR                       = "#e0e0e0",
    CAR_confirmed_in_both         = "#3fb950",
    CAR_present_in_AB_lost_from_I = "#d29922"
  )
  status_labels <- c(
    not_CAR                       = "Non-CAR",
    CAR_confirmed_in_both         = "CAR confermata (in I e AB)",
    CAR_present_in_AB_lost_from_I = "CAR recuperata (AB, persa da I)"
  )

  centroids <- coords %>%
    group_by(cell_type) %>%
    summarise(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2), .groups = "drop")

  n_rec  <- sum(coords$status == "CAR_present_in_AB_lost_from_I")
  n_conf <- sum(coords$status == "CAR_confirmed_in_both")

  p <- ggplot(coords, aes(x = UMAP1, y = UMAP2)) +
    # Non-CAR sotto
    geom_point(data = coords[coords$status == "not_CAR",],
               color = "#e0e0e0", size = 0.4, alpha = 0.4, shape = 16) +
    # CAR confermata
    geom_point(data = coords[coords$status == "CAR_confirmed_in_both",],
               color = "white", size = 2.8, alpha = 0.9, shape = 16) +
    geom_point(data = coords[coords$status == "CAR_confirmed_in_both",],
               color = "#3fb950", size = 1.8, alpha = 0.9, shape = 16) +
    # CAR recuperata (persa da I)
    geom_point(data = coords[coords$status == "CAR_present_in_AB_lost_from_I",],
               color = "white", size = 2.8, alpha = 0.9, shape = 16) +
    geom_point(data = coords[coords$status == "CAR_present_in_AB_lost_from_I",],
               color = "#d29922", size = 1.8, alpha = 0.9, shape = 16) +
    # Label popolazioni
    ggrepel::geom_label_repel(
      data = centroids,
      aes(x = UMAP1, y = UMAP2, label = cell_type),
      size = 3, fontface = "bold",
      fill = scales::alpha("white", 0.7),
      color = "black", label.size = 0.2,
      label.padding = unit(0.15, "lines"),
      max.overlaps = 25, seed = 42, force = 2
    ) +
    # Legenda manuale
    annotate("point", x = -Inf, y = Inf, color = "#3fb950", size = 3) +
    annotate("point", x = -Inf, y = Inf, color = "#d29922", size = 3) +
    scale_x_continuous(expand = expansion(mult = 0.05)) +
    scale_y_continuous(expand = expansion(mult = 0.05)) +
    ggtitle(
      paste0(sample_name, " – CAR recovery status"),
      subtitle = paste0("Verde: CAR confermata (n=", n_conf, ")  |  ",
                        "Arancione: CAR recuperata/persa da I (n=", n_rec, ")")
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title    = element_text(face = "bold", hjust = 0.5, size = 13),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray40"),
      axis.text = element_blank(), axis.ticks = element_blank(),
      axis.title = element_text(size = 9)
    )

  ggsave(paste0(out_dir, sample_name, "_UMAP_CAR_recovery.png"),
         p, width = 10, height = 9, dpi = 300, bg = "white")
  cat(paste0("  → ", sample_name, "_UMAP_CAR_recovery.png\n"))
}

# Applica a tutti i campioni AB aggiornati
for (pid in names(cross_results)) {
  r <- cross_results[[pid]]
  if (is.null(r)) next
  for (nm in names(r$updated_AB))
    plot_recovery_umap(r$updated_AB[[nm]], nm, out_dir)
}

# ============================================================
# RIEPILOGO FINALE
# ============================================================
section("Riepilogo finale")

cat("\n📋 QC AUDIT (cellule rimosse dai filtri di I):\n")
print(audit_tbl, row.names = FALSE)

cat("\n\n📋 CROSS-REFERENCING I ↔ A/B:\n")
print(summary_venn, row.names = FALSE)

cat(paste0(
  "\n", strrep("=", 65), "\n",
  "  ANALISI COMPLETATA\n\n",
  "  Output in: ", out_dir, "\n\n",
  "  QC Audit:\n",
  "    XX_QC_KEPT_vs_LOST.png         – violin QC KEPT vs LOST\n",
  "    QC_audit_CAR_in_removed_cells.xlsx\n\n",
  "  Cross-referencing:\n",
  "    XX_CAR_recovery_barplot.png    – CAR confermata vs recuperata\n",
  "    XX_UMAP_CAR_recovery.png       – UMAP con recovery status\n",
  "    CAR_cross_reference_I_vs_AB.xlsx\n",
  strrep("=", 65), "\n"
))

