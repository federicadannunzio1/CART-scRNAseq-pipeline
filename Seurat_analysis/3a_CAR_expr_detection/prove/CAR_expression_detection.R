# ============================================================
#  CAR-T DETECTION BY CO-EXPRESSION (STRATEGIA BOTH)
#
#  Identifica cellule CAR+ richiedendo che ENTRAMBI i geni
#  del costrutto siano sopra soglia contemporaneamente:
#    TNFRSF9  = 4-1BB (CD137) – dominio costimatatorio
#    CD247    = CD3zeta        – dominio di segnalazione
#
#  PERCHÉ SOLO CO-ESPRESSIONE:
#  ─────────────────────────────────────────────────────────
#  CD247 è espresso endogenamente da tutti i linfociti T
#  come parte del complesso TCR → soglia singola = molti
#  falsi positivi. TNFRSF9 viene upregolato anche sui T
#  effettori endogeni attivati. Richiedere che entrambi
#  siano sopra soglia riduce drasticamente i falsi positivi:
#  la probabilità che un T endogeno abbia overespressione
#  simultanea di entrambi è molto bassa.
#
#  NOTA: lo script richiede che ENTRAMBI i geni siano
#  presenti nella matrice. Se uno manca, probabilmente il
#  genoma di riferimento in CellRanger non includeva le
#  sequenze del costrutto CAR → l'analisi non è applicabile.
#
#  APPROCCIO:
#  1. Per ogni campione: estrae espressione TNFRSF9 e CD247
#  2. Calcola soglia al percentile configurabile
#     → calcolata sulle cellule esprimenti (> 0), non su tutte,
#       per non abbassare artificialmente la soglia verso zero
#     → calcolata sull'intero campione, non per cluster
#  3. CAR_BOTH = CAR+ solo se ENTRAMBI > soglia
#  4. VlnPlot per cluster (TNFRSF9 e CD247, con linea soglia)
#  5. Scatter plot TNFRSF9 vs CD247 con quadrante CAR+
#  6. UMAP con 4 categorie vs scREP:
#       Overlap | Solo scREP | Solo BOTH | CAR- entrambi
#  7. Concordanza quantitativa con IS_CAR_ALLIN_scREP
#  8. Excel riepilogativo
#
#  Input:  all_samples_annotated_COMPLETE.rds
#  Output: <out_dir>/CAR_expr_detection/
# ============================================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(openxlsx)
library(scales)
library(ggrepel)

# ── PARAMETRI CONFIGURABILI ──────────────────────────────────

rds_path <- "/Users/federicadannunzio/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/caruana-project/CART/Code/Seurat_analysis/2_annotation/all_samples_annotated_COMPLETE.rds"
out_dir  <- "/Users/federicadannunzio/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/caruana-project/CART/Code/Seurat_analysis/3a_CAR_expr_detection/res"

if (!dir.exists(out_dir)) {
  dir.create(out_dir)} else {
    print(paste(out_dir, "already exists"))
  }
# Geni del costrutto CAR (servono ENTRAMBI)
GENE_4BB  <- "TNFRSF9"   # 4-1BB / CD137
GENE_CD3Z <- "CD247"     # CD3zeta

# Percentile soglia applicato a entrambi i geni.
# Calcolato sulle cellule che esprimono il gene (> 0).
# 90 = top 10% delle cellule esprimenti (meno stringente)
# 95 = top  5% delle cellule esprimenti (più stringente)
PERCENTILE <- 70

# ─────────────────────────────────────────────────────────────

section <- function(title)
  cat(paste0("\n", strrep("=", 65), "\n  ", title,
             "\n", strrep("=", 65), "\n"))

# ============================================================
# PALETTE
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

get_color <- function(x) {
  col <- PALETTE[x]
  if (is.na(col)) "#888888" else col
}

# ============================================================
# HELPERS
# ============================================================

find_umap <- function(obj) {
  for (nm in c("umap","wnn.umap","umap.harmony","RNA.umap"))
    if (nm %in% names(obj@reductions)) return(nm)
  NULL
}

# Calcola soglia al percentile sulle cellule esprimenti (> 0).
# Escludere gli zeri evita che l'alta quota di dropout
# abbassi artificialmente la soglia verso zero.
calc_threshold <- function(vec, pct) {
  vec_nz <- vec[vec > 0]
  if (length(vec_nz) < 10) {
    cat("    [NOTA] < 10 cellule esprimenti: soglia su tutte.\n")
    return(quantile(vec, pct / 100))
  }
  quantile(vec_nz, pct / 100)
}

# ============================================================
# FUNZIONE PRINCIPALE PER CAMPIONE
# ============================================================

analyze_car_expression <- function(obj, sample_name, out_dir,
                                   pct = PERCENTILE) {

  cat(paste0("\n", strrep("-", 55), "\n",
             "  Campione: ", sample_name, "\n",
             strrep("-", 55), "\n"))

  DefaultAssay(obj) <- "RNA"

  # ── Controlla presenza di entrambi i geni ─────────────────
  missing <- setdiff(c(GENE_4BB, GENE_CD3Z), rownames(obj))
  if (length(missing) > 0) {
    cat(sprintf("  [WARN] Geni mancanti: %s\n",
                paste(missing, collapse = ", ")))
  }
  if (!all(c(GENE_4BB, GENE_CD3Z) %in% rownames(obj))) {
    cat("  [SKIP] Servono ENTRAMBI i geni per la strategia BOTH.\n")
    cat("  Verifica che il gtf usato in CellRanger includesse\n")
    cat("  le sequenze del costrutto CAR.\n")
    return(NULL)
  }

  # ── Estrai espressione ─────────────────────────────────────
  v_4bb  <- as.numeric(GetAssayData(obj, slot = "data")[GENE_4BB,  ])
  v_cd3z <- as.numeric(GetAssayData(obj, slot = "data")[GENE_CD3Z, ])

  # ── Statistiche di base ───────────────────────────────────
  cat(sprintf(
    "  %-8s | Cellule>0: %5d (%5.1f%%) | Med expr: %.3f | Max: %.2f\n",
    GENE_4BB,  sum(v_4bb  > 0), mean(v_4bb  > 0)*100,
    median(v_4bb[v_4bb > 0]),   max(v_4bb)))
  cat(sprintf(
    "  %-8s | Cellule>0: %5d (%5.1f%%) | Med expr: %.3f | Max: %.2f\n",
    GENE_CD3Z, sum(v_cd3z > 0), mean(v_cd3z > 0)*100,
    median(v_cd3z[v_cd3z > 0]), max(v_cd3z)))

  # ── Calcola soglie ────────────────────────────────────────
  thr_4bb  <- calc_threshold(v_4bb,  pct)
  thr_cd3z <- calc_threshold(v_cd3z, pct)
  cat(sprintf("  Soglia %-8s (%d° pct esprimenti): %.4f\n",
              GENE_4BB,  pct, thr_4bb))
  cat(sprintf("  Soglia %-8s (%d° pct esprimenti): %.4f\n",
              GENE_CD3Z, pct, thr_cd3z))

  # ── Classificazione CAR_BOTH ──────────────────────────────
  is_car_both <- v_4bb > thr_4bb & v_cd3z > thr_cd3z
  n_car_both  <- sum(is_car_both)
  cat(sprintf("  CAR+ BOTH (%d° pct): %d cellule (%.2f%%)\n",
              pct, n_car_both, n_car_both / ncol(obj) * 100))

  # ── Metadati ──────────────────────────────────────────────
  meta              <- obj@meta.data
  meta$cell_type    <- as.character(meta$cell_type)
  meta$expr_TNFRSF9 <- v_4bb
  meta$expr_CD247   <- v_cd3z
  meta$CAR_BOTH     <- ifelse(is_car_both, "CAR+", "CAR-")

  has_screp <- "IS_CAR_ALLIN_scREP" %in% colnames(meta)
  if (has_screp) {
    meta$CAR_scREP <- ifelse(meta$IS_CAR_ALLIN_scREP == "YES",
                             "CAR+", "CAR-")
    n_screp <- sum(meta$CAR_scREP == "CAR+")
    cat(sprintf("  CAR+ scREP:         %d cellule (%.2f%%)\n",
                n_screp, n_screp / ncol(obj) * 100))
  }

  pop_order <- sort(unique(meta$cell_type))
  pop_cols  <- setNames(sapply(pop_order, get_color), pop_order)
  Idents(obj) <- "cell_type"

  # ── 1. VlnPlot TNFRSF9 e CD247 per cluster ───────────────
  # La linea tratteggiata mostra la soglia BOTH per quel gene.
  # I cluster dove molte cellule la superano sono arricchiti
  # in CAR-T. Atteso: Cytotoxic CD8+ e Proliferating CD8+
  # dovrebbero avere la coda superiore più lunga.

  make_vln <- function(gene, thr, line_col) {
    VlnPlot(obj,
            features = gene,
            group.by = "cell_type",
            cols     = pop_cols,
            pt.size  = 0.3) +
      geom_hline(yintercept = thr,
                 linetype = "dashed", color = line_col,
                 linewidth = 0.9) +
      annotate("text",
               x     = length(pop_order) * 0.97,
               y     = thr * 1.1,
               label = paste0(pct, "° pct\n(soglia BOTH)"),
               color = line_col, size = 3,
               hjust = 1, fontface = "italic") +
      labs(
        title    = paste0(sample_name, " – ", gene),
        subtitle = paste0(
          "Soglia ", pct, "° percentile (calcolata su cellule esprimenti)\n",
          "Cellule sopra la linea contribuiscono alla call CAR+ BOTH"),
        y = "Espressione (log-norm)", x = NULL) +
      theme_classic(base_size = 10) +
      theme(
        axis.text.x   = element_text(angle = 45, hjust = 1,
                                      size = 8, face = "bold"),
        plot.title    = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, size = 8,
                                      color = "gray40"),
        legend.position = "none")
  }

  p_vln <- make_vln(GENE_4BB, thr_4bb, "#E63946") /
           make_vln(GENE_CD3Z, thr_cd3z, "#4361EE") +
    plot_annotation(
      title = paste0(sample_name,
                     " – Distribuzione geni CAR per cluster"),
      theme = theme(plot.title = element_text(face = "bold",
                                               hjust = 0.5,
                                               size = 13)))

  vln_w <- max(14, length(pop_order) * 0.7 + 3)
  ggsave(paste0(out_dir, sample_name, "_vln_CAR_genes.png"),
         plot = p_vln, width = vln_w, height = 11,
         dpi = 300, bg = "white")
  cat(paste0("  → ", sample_name, "_vln_CAR_genes.png\n"))

  # ── 2. Scatter TNFRSF9 vs CD247 ──────────────────────────
  # Il quadrante in alto a destra = CAR+ BOTH.
  # I punti sono colorati per categoria di concordanza con scREP.

  df_sc <- data.frame(
    TNFRSF9   = v_4bb,
    CD247     = v_cd3z,
    cell_type = meta$cell_type,
    CAR_BOTH  = meta$CAR_BOTH,
    stringsAsFactors = FALSE)

  if (has_screp) {
    df_sc$CAR_scREP <- meta$CAR_scREP
    df_sc$category  <- case_when(
      df_sc$CAR_BOTH == "CAR+" & df_sc$CAR_scREP == "CAR+" ~
        "Overlap (BOTH + scREP)",
      df_sc$CAR_BOTH == "CAR+" & df_sc$CAR_scREP == "CAR-" ~
        "Solo BOTH (nuovo)",
      df_sc$CAR_BOTH == "CAR-" & df_sc$CAR_scREP == "CAR+" ~
        "Solo scREP (non catturato)",
      TRUE ~ "CAR-")
  } else {
    df_sc$category <- ifelse(df_sc$CAR_BOTH == "CAR+",
                             "CAR+ BOTH", "CAR-")
  }

  cat_levels <- c("CAR-",
                  "Solo scREP (non catturato)",
                  "Solo BOTH (nuovo)",
                  "Overlap (BOTH + scREP)",
                  "CAR+ BOTH")
  df_sc$category <- factor(df_sc$category,
                            levels = cat_levels[
                              cat_levels %in% df_sc$category])
  df_sc <- df_sc[order(df_sc$category), ]

  cat_colors <- c(
    "CAR-"                        = "#DDDDDD",
    "Solo scREP (non catturato)"  = "#4361EE",
    "Solo BOTH (nuovo)"           = "#F4A261",
    "Overlap (BOTH + scREP)"      = "#B00020",
    "CAR+ BOTH"                   = "#B00020")
  cat_sizes <- c(
    "CAR-"                        = 0.25,
    "Solo scREP (non catturato)"  = 1.0,
    "Solo BOTH (nuovo)"           = 1.0,
    "Overlap (BOTH + scREP)"      = 1.2,
    "CAR+ BOTH"                   = 1.2)

  p_scatter <- ggplot(df_sc,
    aes(x = TNFRSF9, y = CD247,
        color = category, size = category)) +
    geom_point(alpha = 0.6) +
    geom_vline(xintercept = thr_4bb,
               color = "#E63946", linetype = "dashed",
               linewidth = 0.8) +
    geom_hline(yintercept = thr_cd3z,
               color = "#4361EE", linetype = "dashed",
               linewidth = 0.8) +
    annotate("rect",
             xmin = thr_4bb, xmax = Inf,
             ymin = thr_cd3z, ymax = Inf,
             fill = "#B00020", alpha = 0.05) +
    annotate("text",
             x = max(v_4bb) * 0.72, y = max(v_cd3z) * 0.92,
             label = paste0("CAR+ BOTH\nn = ", n_car_both),
             color = "#B00020", size = 4.5, fontface = "bold") +
    scale_color_manual(
      values = cat_colors[names(cat_colors) %in%
                            levels(df_sc$category)],
      name = NULL) +
    scale_size_manual(
      values = cat_sizes[names(cat_sizes) %in%
                           levels(df_sc$category)],
      guide = "none") +
    labs(
      title    = paste0(sample_name,
                        " – TNFRSF9 vs CD247 (co-espressione BOTH)"),
      subtitle = paste0(
        "Quadrante rosso = CAR+ BOTH (sopra entrambe le soglie)\n",
        "Linea rossa = soglia ", pct, "° pct TNFRSF9  |  ",
        "Linea blu = soglia ", pct, "° pct CD247"),
      x = paste0(GENE_4BB, " (log-norm)"),
      y = paste0(GENE_CD3Z, " (log-norm)")) +
    theme_classic(base_size = 11) +
    theme(
      plot.title    = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = 8.5,
                                    color = "gray40"),
      legend.position = "bottom",
      legend.text     = element_text(size = 9)) +
    guides(color = guide_legend(
      override.aes = list(size = 3, alpha = 1)))

  ggsave(paste0(out_dir, sample_name, "_scatter_BOTH.png"),
         plot = p_scatter, width = 7, height = 7,
         dpi = 300, bg = "white")
  cat(paste0("  → ", sample_name, "_scatter_BOTH.png\n"))

  # ── 3. UMAP – 4 categorie ────────────────────────────────
  umap_key <- find_umap(obj)

  if (!is.null(umap_key)) {
    coords <- as.data.frame(Embeddings(obj, umap_key)[, 1:2])
    colnames(coords) <- c("UMAP1","UMAP2")
    coords$category  <- as.character(df_sc$category)
    coords$category  <- factor(coords$category,
                                levels = levels(df_sc$category))

    centroids <- data.frame(
      cell_type = meta$cell_type,
      UMAP1 = coords$UMAP1, UMAP2 = coords$UMAP2) %>%
      group_by(cell_type) %>%
      summarise(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2),
                .groups = "drop")

    coords_sorted <- coords[order(coords$category), ]

    n_overlap    <- sum(coords$category == "Overlap (BOTH + scREP)")
    n_solo_screp <- sum(coords$category == "Solo scREP (non catturato)")
    n_solo_both  <- sum(coords$category == "Solo BOTH (nuovo)")
    n_car_only   <- sum(coords$category == "CAR+ BOTH")

    subtitle_umap <- if (has_screp) {
      paste0("Rosso = overlap (", n_overlap, ")  |  ",
             "Blu = solo scREP (", n_solo_screp,
             ") = non catturato da expr\n",
             "Arancio = solo BOTH (", n_solo_both,
             ") = potenziali nuovi CAR+  |  Grigio = CAR-")
    } else {
      paste0("Rosso = CAR+ BOTH (", n_car_only, ")  |  Grigio = CAR-")
    }

    umap_alpha <- c("CAR-"                        = 0.25,
                    "Solo scREP (non catturato)"  = 0.85,
                    "Solo BOTH (nuovo)"           = 0.85,
                    "Overlap (BOTH + scREP)"      = 1.0,
                    "CAR+ BOTH"                   = 1.0)
    umap_size  <- c("CAR-"                        = 0.3,
                    "Solo scREP (non catturato)"  = 1.2,
                    "Solo BOTH (nuovo)"           = 1.2,
                    "Overlap (BOTH + scREP)"      = 1.5,
                    "CAR+ BOTH"                   = 1.5)

    lev <- levels(df_sc$category)

    p_umap <- ggplot(coords_sorted,
      aes(x = UMAP1, y = UMAP2,
          color = category,
          size  = category,
          alpha = category)) +
      geom_point() +
      scale_color_manual(
        values = cat_colors[lev], name = NULL) +
      scale_size_manual(
        values = umap_size[lev], guide = "none") +
      scale_alpha_manual(
        values = umap_alpha[lev], guide = "none") +
      ggrepel::geom_label_repel(
        data = centroids,
        aes(x = UMAP1, y = UMAP2, label = cell_type),
        inherit.aes   = FALSE,
        size          = 2.8, fontface = "bold",
        fill          = alpha("white", 0.65), color = "black",
        label.size    = 0.12,
        label.padding = unit(0.1, "lines"),
        max.overlaps  = 25, seed = 42) +
      labs(title    = paste0(sample_name,
                             " – CAR+ BOTH su UMAP"),
           subtitle = subtitle_umap) +
      theme_classic(base_size = 11) +
      theme(
        plot.title    = element_text(face = "bold",
                                      hjust = 0.5, size = 12),
        plot.subtitle = element_text(hjust = 0.5, size = 8.5,
                                      color = "gray40"),
        axis.text     = element_blank(),
        axis.ticks    = element_blank(),
        legend.position = "bottom",
        legend.text     = element_text(size = 9)) +
      guides(color = guide_legend(
        override.aes = list(size = 3, alpha = 1)))

    ggsave(paste0(out_dir, sample_name,
                  "_UMAP_BOTH_vs_scREP.png"),
           plot = p_umap, width = 9, height = 8,
           dpi = 300, bg = "white")
    cat(paste0("  → ", sample_name,
               "_UMAP_BOTH_vs_scREP.png\n"))
  }

  # ── 4. Concordanza quantitativa ───────────────────────────
  concordance_row <- data.frame(
    campione         = sample_name,
    percentile       = pct,
    n_totale_cellule = ncol(obj),
    n_CAR_BOTH       = n_car_both,
    pct_CAR_BOTH     = round(n_car_both / ncol(obj) * 100, 2),
    stringsAsFactors = FALSE)

  if (has_screp) {
    is_screp_pos <- meta$CAR_scREP == "CAR+"
    n_screp_pos  <- sum(is_screp_pos)
    overlap      <- sum(is_car_both &  is_screp_pos)
    solo_screp   <- sum(!is_car_both & is_screp_pos)
    solo_both    <- sum(is_car_both  & !is_screp_pos)
    true_neg     <- sum(!is_car_both & !is_screp_pos)
    sensitivity  <- if (n_screp_pos > 0)
      round(overlap / n_screp_pos * 100, 1) else NA
    specificity  <- if ((ncol(obj) - n_screp_pos) > 0)
      round(true_neg / (ncol(obj) - n_screp_pos) * 100, 1) else NA

    concordance_row$n_CAR_scREP    <- n_screp_pos
    concordance_row$overlap         <- overlap
    concordance_row$solo_scREP      <- solo_screp
    concordance_row$solo_BOTH       <- solo_both
    concordance_row$sensitivity_pct <- sensitivity
    concordance_row$specificity_pct <- specificity

    cat(sprintf("  Overlap (BOTH ∩ scREP):      %d\n", overlap))
    cat(sprintf("  Solo scREP (non catturato):  %d\n", solo_screp))
    cat(sprintf("  Solo BOTH (nuovi candidati): %d\n", solo_both))
    cat(sprintf("  Sensibilità vs scREP: %.1f%%\n", sensitivity))
    cat(sprintf("  Specificità vs scREP: %.1f%%\n", specificity))
  }

  return(list(
    sample      = sample_name,
    meta        = meta,
    concordance = concordance_row,
    thr_4bb     = thr_4bb,
    thr_cd3z    = thr_cd3z
  ))
}

# ============================================================
# CARICAMENTO
# ============================================================
section("Caricamento")

cat("Caricamento:", rds_path, "\n")
all_samples <- readRDS(rds_path)
if (inherits(all_samples, "Seurat")) {
  nm <- unique(all_samples$orig.ident)
  nm <- if (length(nm) == 1) nm else "Sample"
  all_samples <- setNames(list(all_samples), nm)
}
cat(sprintf("Campioni trovati: %d\n", length(all_samples)))
cat(sprintf("Strategia: co-espressione BOTH | Percentile: %d°\n",
            PERCENTILE))

# ============================================================
# LOOP PRINCIPALE
# ============================================================
section(paste0("Analisi BOTH [",
               GENE_4BB, " AND ", GENE_CD3Z, "]"))

results <- list()
for (nm in names(all_samples)) {
  results[[nm]] <- analyze_car_expression(
    all_samples[[nm]], nm, out_dir, pct = PERCENTILE)
}
results <- Filter(Negate(is.null), results)

# ============================================================
# EXCEL RIEPILOGATIVO
# ============================================================
section("Excel riepilogativo")

wb <- createWorkbook()

all_concordance <- bind_rows(lapply(results, `[[`, "concordance"))
addWorksheet(wb, "Concordanza_Globale")
writeData(wb, "Concordanza_Globale", all_concordance)

thresh_df <- bind_rows(lapply(names(results), function(nm) {
  r <- results[[nm]]
  data.frame(campione       = nm,
             percentile     = PERCENTILE,
             soglia_TNFRSF9 = round(r$thr_4bb,  4),
             soglia_CD247   = round(r$thr_cd3z, 4),
             stringsAsFactors = FALSE)
}))
addWorksheet(wb, "Soglie")
writeData(wb, "Soglie", thresh_df)

for (nm in names(results)) {
  r        <- results[[nm]]
  sheet_nm <- substr(nm, 1, 31)
  addWorksheet(wb, sheet_nm)
  cols_keep <- c("cell_type", "expr_TNFRSF9", "expr_CD247",
                 "CAR_BOTH",
                 if ("CAR_scREP" %in% colnames(r$meta))
                   "CAR_scREP" else NULL)
  writeData(wb, sheet_nm,
            r$meta[, cols_keep[cols_keep %in% colnames(r$meta)]])
}

xlsx_path <- paste0(out_dir, "CAR_expression_detection.xlsx")
saveWorkbook(wb, xlsx_path, overwrite = TRUE)
cat(paste0("  → ", xlsx_path, "\n"))

# ============================================================
# RIEPILOGO FINALE
# ============================================================
section("Riepilogo finale")

cat("\nRisultati co-espressione BOTH per campione:\n\n")
has_screp_global <- "sensitivity_pct" %in%
                    colnames(all_concordance)

if (has_screp_global) {
  cat(sprintf("  %-18s | %7s | %7s | %9s | %9s | %6s | %6s\n",
              "Campione","n_BOTH","n_scREP",
              "Overlap","SoloSCREP","Sens%","Spec%"))
  cat(paste0("  ", strrep("-", 76), "\n"))
  for (i in seq_len(nrow(all_concordance))) {
    r <- all_concordance[i, ]
    cat(sprintf(
      "  %-18s | %7d | %7d | %9d | %9d | %5.1f%% | %5.1f%%\n",
      r$campione, r$n_CAR_BOTH, r$n_CAR_scREP,
      r$overlap, r$solo_scREP,
      r$sensitivity_pct, r$specificity_pct))
  }
} else {
  cat(sprintf("  %-18s | %8s | %6s\n",
              "Campione","n_BOTH","pct"))
  cat(paste0("  ", strrep("-", 36), "\n"))
  for (i in seq_len(nrow(all_concordance))) {
    r <- all_concordance[i, ]
    cat(sprintf("  %-18s | %8d | %5.2f%%\n",
                r$campione, r$n_CAR_BOTH, r$pct_CAR_BOTH))
  }
}

cat(paste0(
  "\n", strrep("=", 65), "\n",
  "  ANALISI COMPLETATA\n\n",
  "  Output: ", out_dir, "\n\n",
  "  Per campione:\n",
  "    <sample>_vln_CAR_genes.png       violinplot TNFRSF9 e CD247\n",
  "    <sample>_scatter_BOTH.png        scatter con quadrante CAR+\n",
  "    <sample>_UMAP_BOTH_vs_scREP.png  UMAP confronto scREP\n",
  "  Globale:\n",
  "    CAR_expression_detection.xlsx\n\n",
  "  COME INTERPRETARE:\n",
  "  ─ Overlap alto    → i due metodi concordano bene\n",
  "  ─ Solo scREP alto → dropout o reads CAR non mappati:\n",
  "                      l'expr non cattura tutte le CAR-T\n",
  "  ─ Solo BOTH alto  → potenziali CAR-T mancate da scREP;\n",
  "                      controlla se sono nei cluster\n",
  "                      Cytotoxic CD8+ o Proliferating\n",
  "  ─ Sensibilità < 50% → prova ad abbassare PERCENTILE a 85\n",
  strrep("=", 65), "\n"))
