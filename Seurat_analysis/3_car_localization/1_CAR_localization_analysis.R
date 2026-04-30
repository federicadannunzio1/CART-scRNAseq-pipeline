# ============================================================
#  CAR-T CELL LOCALIZATION ANALYSIS
#
#  Identifica in quali popolazioni cellulari ricadono le
#  cellule CAR+ (IS_CAR_ALLIN_scREP == "YES") in ogni campione.
#
#  Input:  all_samples_annotated_COMPLETE.rds
#          (lista di oggetti Seurat con cell_type e UMAP)
#
#  Output in <base_dir>/CAR_analysis/:
#    Per campione:
#      <sample>_UMAP_CAR_overlay.png   – UMAP cell_type + CAR overlay
#      <sample>_CAR_barplot.png        – % CAR per popolazione
#    Combinato:
#      ALL_samples_CAR_overview.png    – pannello multi-campione
#      CAR_distribution_all_samples.xlsx – tabelle per ogni campione
#
#  Nota: IS_CAR_ALLIN_scREP può essere assente in alcuni campioni
#  (es. quelli senza dati VDJ/scREP). Il codice lo gestisce
#  saltando il campione con un avviso.
# ============================================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(openxlsx)
library(scales)
library(ggrepel)

# ── UNICO PUNTO DA MODIFICARE ────────────────────────────────
rds_path <- "/Users/federicadannunzio/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/caruana-project/CART/Code/Seurat_analysis/2_annotation/all_samples_annotated_COMPLETE.rds"
out_dir  <- "/Users/federicadannunzio/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/caruana-project/CART/Code/Seurat_analysis/3_car_localization/res/"
# ─────────────────────────────────────────────────────────────

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

section <- function(title)
  cat(paste0("\n", strrep("=", 65), "\n  ", title,
             "\n", strrep("=", 65), "\n"))

# ============================================================
# PALETTE E ORDINE CANONICO
# (identici alle Pipeline 1-3 per coerenza visiva)
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
  "Platelets"                  = "#CDB4DB",
  "Unknown"                    = "#AAAAAA"
)

get_colors <- function(types) {
  cols    <- PALETTE[types]
  missing <- types[is.na(cols)]
  if (length(missing) > 0) {
    extra        <- setNames(hue_pal()(length(missing)), missing)
    cols[missing] <- extra
  }
  cols
}

# ============================================================
# HELPERS
# ============================================================

# Trova la riduzione UMAP presente nell'oggetto (Seurat v4/v5)
find_umap <- function(obj) {
  for (nm in c("umap","wnn.umap","umap.harmony","umap2",
               "RNA.umap","integrated.umap")) {
    if (nm %in% names(obj@reductions)) return(nm)
  }
  stop("Nessuna riduzione UMAP trovata. Riduzioni disponibili: ",
       paste(names(obj@reductions), collapse=", "))
}

# Estrae coordinate UMAP + metadati in un data.frame pulito
get_umap_df <- function(obj, umap_key) {
  coords <- as.data.frame(
    Embeddings(obj, reduction = umap_key)[, 1:2])
  colnames(coords) <- c("UMAP1","UMAP2")
  coords$cell_type <- as.character(obj$cell_type)
  coords$IS_CAR    <- as.character(obj$IS_CAR_ALLIN_scREP)
  # Normalizza: qualsiasi valore non "YES" → "NO"
  coords$IS_CAR    <- ifelse(coords$IS_CAR == "YES", "YES", "NO")
  coords$barcode   <- rownames(coords)
  coords
}

# ============================================================
# FUNZIONE PRINCIPALE: analisi CAR per singolo campione
# ============================================================

analyze_car_sample <- function(obj, sample_name, out_dir) {

  cat(paste0("\n", strrep("-", 55), "\n",
             "  Campione: ", sample_name, "\n",
             strrep("-", 55), "\n"))

  # ── Controllo colonna IS_CAR_ALLIN_scREP ──────────────────
  if (!"IS_CAR_ALLIN_scREP" %in% colnames(obj@meta.data)) {
    cat(paste0("  [SKIP] IS_CAR_ALLIN_scREP non trovata in ",
               sample_name, "\n"))
    cat("  Colonne disponibili:\n")
    cat(paste0("  ", paste(colnames(obj@meta.data), collapse=", "), "\n"))
    return(NULL)
  }

  # ── Controllo cell_type ────────────────────────────────────
  if (!"cell_type" %in% colnames(obj@meta.data)) {
    cat(paste0("  [WARN] cell_type non trovata in ", sample_name,
               " – uso seurat_clusters.\n"))
    obj$cell_type <- as.character(obj$seurat_clusters)
  }

  umap_key <- tryCatch(find_umap(obj),
    error = function(e) {
      cat(paste0("  [SKIP] ", conditionMessage(e), "\n"))
      return(NULL)
    })
  if (is.null(umap_key)) return(NULL)

  df <- get_umap_df(obj, umap_key)

  n_total <- nrow(df)
  n_yes   <- sum(df$IS_CAR == "YES")
  n_no    <- sum(df$IS_CAR == "NO")
  pct_car <- round(100 * n_yes / n_total, 1)

  cat(sprintf("  Cellule totali : %d\n", n_total))
  cat(sprintf("  CAR+ (YES)     : %d (%.1f%%)\n", n_yes, pct_car))
  cat(sprintf("  CAR- (NO)      : %d (%.1f%%)\n", n_no,
              round(100 * n_no / n_total, 1)))

  if (n_yes == 0) {
    cat("  [INFO] Nessuna cellula CAR+ in questo campione.\n")
    return(NULL)
  }

  # ── Tabella statistica per popolazione ────────────────────
  tbl <- df %>%
    group_by(cell_type, IS_CAR) %>%
    summarise(n = n(), .groups = "drop") %>%
    pivot_wider(names_from = IS_CAR,
                values_from = n, values_fill = 0) %>%
    rename_with(~ paste0("n_", .), -cell_type) %>%
    mutate(
      n_YES = if ("n_YES" %in% names(.)) n_YES else 0L,
      n_NO  = if ("n_NO"  %in% names(.)) n_NO  else 0L,
      n_tot = n_YES + n_NO,
      pct_CAR_in_pop   = round(100 * n_YES / n_tot, 1),
      pct_of_all_CARs  = round(100 * n_YES / max(n_yes, 1), 1)
    ) %>%
    arrange(desc(n_YES))

  # Avvisa se CAR+ si trovano in popolazioni non-T
  non_t_pop <- c("B cells","Memory B cells","Plasma cells",
                 "NK cells","NKT cells","ILC","MAIT cells",
                 "CD14 Monocytes","CD16 Monocytes","Myeloid cells",
                 "Basophils","HSPC","Erythroid cells","Platelets")
  car_in_nonT <- tbl %>% filter(cell_type %in% non_t_pop & n_YES > 0)
  if (nrow(car_in_nonT) > 0) {
    cat(paste0("  [ATTENZIONE] CAR+ rilevate in popolazioni non-T:\n"))
    for (i in seq_len(nrow(car_in_nonT)))
      cat(sprintf("    %s: %d cellule (%.1f%% della pop)\n",
                  car_in_nonT$cell_type[i],
                  car_in_nonT$n_YES[i],
                  car_in_nonT$pct_CAR_in_pop[i]))
  }

  cat("\nDistribuzione CAR+ per popolazione:\n")
  print(as.data.frame(tbl), row.names = FALSE)

  # ── Costruzione UMAP ──────────────────────────────────────
  #
  # Layout del plot (due pannelli affiancati):
  #   A) UMAP colorato per cell_type, CAR+ sovrapposto
  #      come punti neri con alone bianco (size maggiore)
  #   B) UMAP dicotomico: CAR+ (rosso) vs CAR- (grigio chiaro)
  #      con label delle popolazioni in grassetto
  #
  # Perché due pannelli: A mostra il contesto biologico,
  # B mostra dove si trovano le CAR in modo pulito.

  present_types <- sort(unique(df$cell_type))
  cols_ct       <- get_colors(present_types)
  names(cols_ct) <- present_types

  df_no  <- df[df$IS_CAR == "NO",  ]
  df_yes <- df[df$IS_CAR == "YES", ]

  # Calcola centroidi delle popolazioni per label
  centroids <- df %>%
    group_by(cell_type) %>%
    summarise(UMAP1 = median(UMAP1),
              UMAP2 = median(UMAP2), .groups = "drop")

  # ── Pannello A: cell_type + overlay CAR ──────────────────
  p_A <- ggplot() +
    # Layer 1: tutte le cellule colorate per tipo
    geom_point(data = df_no,
               aes(x = UMAP1, y = UMAP2, color = cell_type),
               size = 0.5, alpha = 0.55, shape = 16) +
    # Layer 2: cellule CAR- colorate (più scure, dimensione normale)
    # [già incluse sopra: le CAR- sono tutte quelle IS_CAR=="NO"]
    # Layer 3: cellule CAR+ → alone bianco + punto rosso scuro
    geom_point(data = df_yes,
               aes(x = UMAP1, y = UMAP2),
               color = "white", size = 2.5, shape = 16, alpha = 0.9) +
    geom_point(data = df_yes,
               aes(x = UMAP1, y = UMAP2),
               color = "#B00020", size = 1.8, shape = 16, alpha = 0.9) +
    # Label popolazioni
    ggrepel::geom_label_repel(
      data        = centroids,
      aes(x = UMAP1, y = UMAP2, label = cell_type),
      size        = 3.0, fontface = "bold",
      fill        = alpha("white", 0.70),
      color       = "black", label.size  = 0.2,
      label.padding = unit(0.15, "lines"),
      max.overlaps = 25, seed = 42, force = 2) +
    scale_color_manual(values = cols_ct, guide = "none") +
    ggtitle(paste0(sample_name,
                   "\nPopulazioni + CAR+ (rosso, n=", n_yes, ")")) +
    theme_classic(base_size = 11) +
    theme(plot.title   = element_text(face = "bold", hjust = 0.5,
                                       size = 11),
          axis.text    = element_blank(),
          axis.ticks   = element_blank(),
          axis.title   = element_text(size = 9))

  # ── Pannello B: CAR+ vs CAR- dicotomico ──────────────────
  # Ordine strati: CAR- prima (grigio, sotto), CAR+ sopra (rosso)
  df_plot_B <- df %>%
    arrange(IS_CAR)   # NO prima, YES dopo → YES sovrapposto

  p_B <- ggplot() +
    geom_point(data = df_plot_B[df_plot_B$IS_CAR == "NO", ],
               aes(x = UMAP1, y = UMAP2),
               color = "#DDDDDD", size = 0.5, alpha = 0.5,
               shape = 16) +
    geom_point(data = df_plot_B[df_plot_B$IS_CAR == "YES", ],
               aes(x = UMAP1, y = UMAP2),
               color = "#B00020", size = 1.2, alpha = 0.85,
               shape = 16) +
    ggrepel::geom_label_repel(
      data        = centroids,
      aes(x = UMAP1, y = UMAP2, label = cell_type),
      size        = 3.0, fontface = "bold",
      fill        = alpha("white", 0.70),
      color       = "black", label.size  = 0.2,
      label.padding = unit(0.15, "lines"),
      max.overlaps = 25, seed = 42, force = 2) +
    ggtitle(paste0(sample_name,
                   "\nCAR+ (rosso) vs CAR- (grigio)")) +
    theme_classic(base_size = 11) +
    theme(plot.title   = element_text(face = "bold", hjust = 0.5,
                                       size = 11),
          axis.text    = element_blank(),
          axis.ticks   = element_blank(),
          axis.title   = element_text(size = 9))

  # ── Pannello C: barplot CAR+ per popolazione ─────────────
  # Mostra solo le popolazioni con almeno 1 cellula CAR+.
  # Due barre affiancate: % CAR nella pop (quanto è CAR+)
  # e % delle CAR totali che cadono in quella pop (composizione).

  tbl_bar <- tbl %>%
    filter(n_YES > 0) %>%
    arrange(desc(pct_of_all_CARs)) %>%
    mutate(cell_type = factor(cell_type, levels = rev(cell_type)))

  bar_cols <- get_colors(as.character(tbl_bar$cell_type))

  p_C <- ggplot(tbl_bar,
                aes(y = cell_type, x = pct_of_all_CARs,
                    fill = cell_type)) +
    geom_col(width = 0.7, show.legend = FALSE) +
    geom_text(aes(label = paste0(n_YES, " celle\n(",
                                  pct_CAR_in_pop, "% della pop)")),
              hjust = -0.05, size = 2.8, color = "gray20") +
    scale_fill_manual(values = bar_cols) +
    scale_x_continuous(
      expand = expansion(mult = c(0, 0.35)),
      labels = function(x) paste0(x, "%")) +
    labs(
      title = paste0(sample_name, " – Distribuzione CAR+"),
      subtitle = paste0("n CAR+ = ", n_yes,
                        " (", pct_car, "% del campione)"),
      x = "% delle CAR+ totali del campione",
      y = NULL) +
    theme_classic(base_size = 11) +
    theme(
      plot.title    = element_text(face = "bold", hjust = 0,
                                    size = 12),
      plot.subtitle = element_text(color = "gray30", size = 10),
      axis.text.y   = element_text(face = "bold", size = 9),
      axis.text.x   = element_text(size = 9),
      panel.grid.major.x = element_line(color = "gray90",
                                         linewidth = 0.3)
    )

  # ── Salvataggio ───────────────────────────────────────────
  n_pops_car <- nrow(tbl_bar)
  bar_h      <- max(4, n_pops_car * 0.55 + 2)

  # UMAP overlay
  umap_combined <- (p_A | p_B) +
    plot_annotation(
      title = sample_name,
      theme = theme(plot.title = element_text(face = "bold",
                                               hjust = 0.5,
                                               size = 13)))
  ggsave(paste0(out_dir, sample_name, "_UMAP_CAR_overlay.png"),
         plot = umap_combined,
         width = 16, height = 7, dpi = 300, bg = "white")
  cat(paste0("  → ", sample_name, "_UMAP_CAR_overlay.png\n"))

  # Barplot
  ggsave(paste0(out_dir, sample_name, "_CAR_barplot.png"),
         plot = p_C,
         width = 10, height = bar_h, dpi = 300, bg = "white")
  cat(paste0("  → ", sample_name, "_CAR_barplot.png\n"))

  return(list(tbl = tbl, n_yes = n_yes, n_total = n_total,
              pct = pct_car, obj = obj,
              umap_key = umap_key, df = df))
}

# ============================================================
# CARICAMENTO E LOOP PRINCIPALE
# ============================================================
section("Caricamento oggetto annotato")

cat("Caricamento:", rds_path, "\n")
all_samples <- readRDS(rds_path)

# Supporto sia lista che singolo oggetto Seurat
if (inherits(all_samples, "Seurat")) {
  cat("[INFO] Oggetto singolo Seurat rilevato.\n")
  # Prova a recuperare il nome del campione da orig.ident
  nm <- unique(all_samples$orig.ident)
  nm <- if (length(nm) == 1) nm else "Sample"
  all_samples <- list(all_samples)
  names(all_samples) <- nm
}

cat(paste0("Campioni trovati: ", length(all_samples), "\n"))
for (nm in names(all_samples)) {
  obj <- all_samples[[nm]]
  has_car <- "IS_CAR_ALLIN_scREP" %in% colnames(obj@meta.data)
  has_ct  <- "cell_type" %in% colnames(obj@meta.data)
  n_car   <- if (has_car)
    sum(as.character(obj$IS_CAR_ALLIN_scREP) == "YES") else NA
  cat(sprintf("  %-18s | %5d celle | cell_type: %s | IS_CAR: %s (n_YES=%s)\n",
              nm, ncol(obj),
              if (has_ct) "OK" else "MANCANTE",
              if (has_car) "OK" else "MANCANTE",
              if (is.na(n_car)) "N/A" else n_car))
}

# ============================================================
section("Analisi CAR per campione")

results <- list()
for (nm in names(all_samples)) {
  results[[nm]] <- analyze_car_sample(all_samples[[nm]], nm, out_dir)
}

# Rimuovi campioni saltati
results <- Filter(Negate(is.null), results)

if (length(results) == 0) {
  stop("[ERRORE] Nessun campione elaborato. Verifica il metadato IS_CAR_ALLIN_scREP.")
}

# ============================================================
# EXCEL RIEPILOGATIVO
# ============================================================
section("Salvataggio Excel")

wb <- createWorkbook()

# Foglio 1: riepilogo per campione
summ_df <- bind_rows(lapply(names(results), function(nm) {
  r <- results[[nm]]
  data.frame(campione  = nm,
             n_totale  = r$n_total,
             n_CAR_pos = r$n_yes,
             pct_CAR   = r$pct,
             stringsAsFactors = FALSE)
}))
addWorksheet(wb, "Riepilogo")
writeData(wb, "Riepilogo", summ_df)

# Un foglio per campione con tabella dettagliata
for (nm in names(results)) {
  sheet_nm <- substr(nm, 1, 31)  # Excel max 31 caratteri per sheet
  addWorksheet(wb, sheet_nm)
  writeData(wb, sheet_nm, results[[nm]]$tbl)
}

xlsx_path <- paste0(out_dir, "CAR_distribution_all_samples.xlsx")
saveWorkbook(wb, xlsx_path, overwrite = TRUE)
cat(paste0("  → ", xlsx_path, "\n"))

# ============================================================
# PANNELLO MULTI-CAMPIONE
# ============================================================
section("Pannello multi-campione")

# Mini-UMAP per ogni campione (CAR+ rosso, CAR- grigio)
# con label grassetto – formato compatto per visione d'insieme

mini_plots <- lapply(names(results), function(nm) {
  r  <- results[[nm]]
  df <- r$df

  centroids <- df %>%
    group_by(cell_type) %>%
    summarise(UMAP1 = median(UMAP1),
              UMAP2 = median(UMAP2), .groups = "drop")

  n_yes <- r$n_yes

  ggplot() +
    geom_point(data = df[df$IS_CAR == "NO", ],
               aes(x = UMAP1, y = UMAP2),
               color = "#E0E0E0", size = 0.35, alpha = 0.45,
               shape = 16) +
    geom_point(data = df[df$IS_CAR == "YES", ],
               aes(x = UMAP1, y = UMAP2),
               color = "#B00020", size = 0.8, alpha = 0.85,
               shape = 16) +
    ggrepel::geom_label_repel(
      data        = centroids,
      aes(x = UMAP1, y = UMAP2, label = cell_type),
      size        = 2.3, fontface = "bold",
      fill        = alpha("white", 0.65),
      color       = "black", label.size  = 0.15,
      label.padding = unit(0.10, "lines"),
      max.overlaps = 20, seed = 42, force = 1.5) +
    ggtitle(paste0(nm, "  [CAR+=", n_yes, "]")) +
    theme_classic(base_size = 9) +
    theme(plot.title   = element_text(face = "bold", hjust = 0.5,
                                       size = 9),
          axis.text    = element_blank(),
          axis.ticks   = element_blank(),
          axis.title   = element_text(size = 8))
})

n_plots  <- length(mini_plots)
n_cols   <- min(3, n_plots)
n_rows   <- ceiling(n_plots / n_cols)

p_overview <- patchwork::wrap_plots(mini_plots, ncol = n_cols) +
  plot_annotation(
    title    = "CAR+ cells in tutte le popolazioni (rosso = CAR+)",
    subtitle = "Le label mostrano la popolazione cellulare",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 14,
                                    hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5,
                                    color = "gray40")))

ggsave(paste0(out_dir, "ALL_samples_CAR_overview.png"),
       plot   = p_overview,
       width  = n_cols * 6,
       height = n_rows * 5.5,
       dpi    = 300, bg = "white")
cat(paste0("  → ALL_samples_CAR_overview.png\n"))

# ============================================================
# RIEPILOGO FINALE A CONSOLE
# ============================================================
section("Riepilogo finale")

cat("\nDistribuzione CAR+ per campione e popolazione:\n\n")
for (nm in names(results)) {
  r   <- results[[nm]]
  top <- r$tbl %>% filter(n_YES > 0) %>%
    arrange(desc(n_YES)) %>%
    head(5)
  cat(sprintf("  %-18s | CAR+ totali: %d (%.1f%%)\n",
              nm, r$n_yes, r$pct))
  for (i in seq_len(nrow(top)))
    cat(sprintf("    %-32s %4d celle | %.1f%% della pop | %.1f%% delle CAR\n",
                top$cell_type[i], top$n_YES[i],
                top$pct_CAR_in_pop[i], top$pct_of_all_CARs[i]))
  cat("\n")
}

cat(paste0(strrep("=", 65), "\n",
           "  ANALISI CAR COMPLETATA\n\n",
           "  Output: ", out_dir, "\n",
           "  Per campione:\n",
           "    <sample>_UMAP_CAR_overlay.png\n",
           "    <sample>_CAR_barplot.png\n",
           "  Globale:\n",
           "    ALL_samples_CAR_overview.png\n",
           "    CAR_distribution_all_samples.xlsx\n",
           strrep("=", 65), "\n"))

