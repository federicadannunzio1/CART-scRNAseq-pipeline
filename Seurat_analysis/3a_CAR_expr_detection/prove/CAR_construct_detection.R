# ============================================================
#  CAR-T DETECTION VIA FEATURE DEL COSTRUTTO
#
#  Logica: quando CellRanger allinea i reads al genoma custom
#  con il costrutto CAR aggiunto al gtf, crea una feature
#  dedicata nella matrice di conteggio. Quella feature conta
#  SOLO i reads del transgene, non dell'endogeno.
#
#  Questo script:
#  1. Scansiona TUTTE le feature dell'oggetto Seurat e stampa
#     quelle non-standard → permette di identificare il nome
#     della feature CAR nel tuo specifico dataset
#  2. Tenta un rilevamento automatico della feature CAR con
#     keyword comuni (configurabile)
#  3. Se la feature è trovata: usa la sua espressione per
#     chiamare CAR+ le cellule sopra soglia (percentile
#     calcolato sulle cellule esprimenti)
#  4. Confronto con IS_CAR_ALLIN_scREP:
#     ─ Overlap       : entrambi concordano → vera positivo
#     ─ Solo scREP    : reads costrutto assenti ma VDJ+ → ok
#     ─ Solo expr     : reads costrutto presenti, VDJ dropout
#                       → CELLULE PROBABILMENTE PERSE DA scREP
#     ─ CAR- entrambi : CAR-
#  5. VlnPlot + density plot della feature CAR per cluster
#  6. UMAP con 4 categorie
#  7. Excel con metadati e concordanza
#
#  PERCHÉ QUESTO APPROCCIO È PIÙ AFFIDABILE DI TNFRSF9/CD247:
#  ─────────────────────────────────────────────────────────
#  TNFRSF9 e CD247 sono geni endogeni: i reads del costrutto
#  CAR (che usa le stesse sequenze) si allineano al locus
#  endogeno, non a una feature separata. Quindi l'espressione
#  di quei geni riflette principalmente il fondo endogeno.
#  La feature del costrutto nel gtf custom è l'unica
#  rappresentazione trascrizionalmente specifica del transgene.
#
#  PERCHÉ scREP PUÒ AVER PERSO CELLULE CAR:
#  ─────────────────────────────────────────────────────────
#  scREP/IS_CAR_ALLIN_scREP si basa sull'analisi VDJ. Il VDJ
#  in scRNA-seq ha un dropout del 30-50%: cellule T reali che
#  non hanno un clonotype VDJ recuperato non vengono chiamate
#  CAR+ da scREP, anche se hanno reads del costrutto.
#  Le cellule "Solo expr" (sotto) sono i candidati mancati.
#
#  Input:  all_samples_annotated_COMPLETE.rds
#  Output: <out_dir>/CAR_construct_detection/
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
out_dir  <- "/Users/federicadannunzio/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/caruana-project/CART/Code/Seurat_analysis/3a_CAR_expr_detection/res_2/"

# Keyword per la ricerca automatica della feature CAR.
# Lo script cerca feature il cui nome contiene almeno una
# di queste stringhe (case-insensitive).
# Aggiungi il nome esatto se lo conosci già.
CAR_KEYWORDS <- c("CAR", "car", "FMC63", "fmc63",
                  "CD19", "scFv", "scfv",
                  "construct", "transgene",
                  "lentiviral", "retroviral",
                  "CART", "cart")

# Percentile soglia per chiamare CAR+ dalla feature expr.
# Calcolato sulle cellule con espressione > 0.
# Abbassa a 85 se la sensibilità vs scREP è bassa.
PERCENTILE <- 90

# Numero massimo di feature da stampare nella scansione
# iniziale (per non intasare la console)
MAX_NONSTANDARD_PRINT <- 50

# ─────────────────────────────────────────────────────────────

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

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
  col <- PALETTE[x]; if (is.na(col)) "#888888" else col
}

# ============================================================
# HELPERS
# ============================================================

find_umap <- function(obj) {
  for (nm in c("umap","wnn.umap","umap.harmony","RNA.umap"))
    if (nm %in% names(obj@reductions)) return(nm)
  NULL
}

calc_threshold <- function(vec, pct) {
  vec_nz <- vec[vec > 0]
  if (length(vec_nz) < 5) {
    cat("    [NOTA] Meno di 5 cellule esprimenti: soglia su tutte.\n")
    return(quantile(vec, pct / 100))
  }
  quantile(vec_nz, pct / 100)
}

# Determina se una feature è un gene umano standard hg38.
# Gene symbols umani tipici: lettere maiuscole + numeri,
# talvolta con trattino. ENSG... = Ensembl ID.
# Feature non-standard = potenziale costrutto esogeno.
is_human_gene <- function(name) {
  grepl("^[A-Z][A-Z0-9]{1,7}([-\\.][A-Z0-9]+)?$", name) |
  grepl("^ENSG[0-9]+", name) |
  grepl("^MT-", name)
}

# ============================================================
# STEP 1: SCANSIONE FEATURE — eseguita su tutti i campioni
# ============================================================

section("Step 1 – Scansione feature non-standard")

cat(paste0(
  "\nQuesta sezione stampa tutte le feature presenti nella\n",
  "matrice di espressione che NON sembrano geni umani standard.\n",
  "Cerca il nome del tuo costrutto CAR in questo elenco.\n",
  "Se non lo trovi, o se l'elenco è vuoto, significa che il\n",
  "costrutto CAR NON è stato aggiunto al gtf di CellRanger\n",
  "come feature separata → l'identificazione tramite espressione\n",
  "non è applicabile e IS_CAR_ALLIN_scREP rimane il metodo\n",
  "più affidabile disponibile.\n\n"))

cat("Caricamento:", rds_path, "\n")
all_samples <- readRDS(rds_path)
if (inherits(all_samples, "Seurat")) {
  nm <- unique(all_samples$orig.ident)
  nm <- if (length(nm) == 1) nm else "Sample"
  all_samples <- setNames(list(all_samples), nm)
}
cat(sprintf("Campioni trovati: %d\n\n", length(all_samples)))

# Raccoglie tutte le feature da tutti i campioni
all_features <- unique(unlist(lapply(all_samples, rownames)))
cat(sprintf("Feature totali (unione campioni): %d\n\n",
            length(all_features)))

nonstandard <- all_features[!is_human_gene(all_features)]
cat(sprintf("Feature non-standard trovate: %d\n", length(nonstandard)))

if (length(nonstandard) == 0) {
  cat(paste0(
    "\n  [ATTENZIONE] Nessuna feature non-standard trovata.\n",
    "  Possibili cause:\n",
    "  1. Il costrutto CAR non è stato aggiunto al gtf custom\n",
    "     come feature separata. I reads del transgene vengono\n",
    "     catturati dai loci endogeni (TNFRSF9, CD247) e non\n",
    "     sono distinguibili dall'espressione endogena.\n",
    "  2. Il nome del costrutto segue la stessa convenzione\n",
    "     dei geni umani standard → non è rilevabile in modo\n",
    "     automatico. Cerca manualmente tra le feature note.\n\n",
    "  CONSEGUENZA: senza una feature specifica del costrutto,\n",
    "  IS_CAR_ALLIN_scREP (basato su VDJ) rimane l'unico metodo\n",
    "  disponibile. La perdita di cellule CAR per VDJ dropout\n",
    "  non è recuperabile tramite espressione genica.\n\n",
    "  Puoi impostare manualmente CAR_FEATURE_NAME (vedi sotto)\n",
    "  se conosci il nome esatto della feature nel tuo gtf.\n"))
} else {
  n_print <- min(length(nonstandard), MAX_NONSTANDARD_PRINT)
  cat(paste0("\n  Prime ", n_print, " feature non-standard:\n"))
  for (f in head(nonstandard, n_print))
    cat(paste0("    ", f, "\n"))
  if (length(nonstandard) > MAX_NONSTANDARD_PRINT)
    cat(sprintf("    ... e altre %d\n",
                length(nonstandard) - MAX_NONSTANDARD_PRINT))
}

# Ricerca automatica con keyword
car_candidates <- nonstandard[
  sapply(nonstandard, function(f)
    any(sapply(CAR_KEYWORDS, function(k)
      grepl(k, f, ignore.case = TRUE))))]

# Cerca anche tra le feature standard (nel caso il nome del
# costrutto segua la stessa convenzione dei geni umani)
car_candidates_all <- all_features[
  sapply(all_features, function(f)
    any(sapply(CAR_KEYWORDS, function(k)
      grepl(k, f, ignore.case = TRUE))))]

cat(paste0("\n  Feature che contengono keyword CAR ",
           "(in tutte le feature, incluse standard):\n"))
if (length(car_candidates_all) == 0) {
  cat("    Nessuna trovata con le keyword attuali.\n")
  cat("    → Modifica CAR_KEYWORDS con il nome esatto del tuo costruttto.\n")
} else {
  for (f in car_candidates_all)
    cat(paste0("    ", f, "\n"))
}

# ── Imposta manualmente qui se la ricerca automatica fallisce
# Esempi: "CD19CAR", "CAR_construct", "FMC63", "EGFRt", ecc.
# Lascia NULL per usare il primo candidato trovato in automatico.
CAR_FEATURE_NAME <- NULL

# Selezione finale della feature
if (is.null(CAR_FEATURE_NAME)) {
  if (length(car_candidates_all) == 1) {
    CAR_FEATURE_NAME <- car_candidates_all[1]
    cat(sprintf("\n  Feature CAR selezionata automaticamente: %s\n",
                CAR_FEATURE_NAME))
  } else if (length(car_candidates_all) > 1) {
    CAR_FEATURE_NAME <- car_candidates_all[1]
    cat(sprintf(paste0(
      "\n  [ATTENZIONE] %d candidati trovati. Usato il primo: %s\n",
      "  Modifica CAR_FEATURE_NAME se non è quello corretto.\n"),
      length(car_candidates_all), CAR_FEATURE_NAME))
  } else {
    cat(paste0(
      "\n  [STOP] Nessuna feature CAR identificata.\n",
      "  Imposta CAR_FEATURE_NAME manualmente nel codice.\n",
      "  L'analisi si ferma qui.\n\n"))
    # Salva report scansione e termina
    sink(paste0(out_dir, "feature_scan_report.txt"))
    cat("Feature non-standard:\n")
    cat(paste(nonstandard, collapse = "\n"))
    cat("\n\nFeature con keyword CAR:\n")
    cat(paste(car_candidates_all, collapse = "\n"))
    sink()
    cat(paste0("  Report salvato: ", out_dir,
               "feature_scan_report.txt\n"))
    stop("CAR_FEATURE_NAME non determinato. ",
         "Imposta il nome manualmente nel codice.")
  }
}

cat(sprintf(paste0(
  "\n  Feature CAR usata per l'analisi: [%s]\n",
  "  Percentile soglia: %d°\n"),
  CAR_FEATURE_NAME, PERCENTILE))

# ============================================================
# FUNZIONE PRINCIPALE PER CAMPIONE
# ============================================================

analyze_car_construct <- function(obj, sample_name, out_dir,
                                  car_feature, pct) {

  cat(paste0("\n", strrep("-", 55), "\n",
             "  Campione: ", sample_name, "\n",
             strrep("-", 55), "\n"))

  DefaultAssay(obj) <- "RNA"

  # ── Controlla presenza della feature ─────────────────────
  if (!car_feature %in% rownames(obj)) {
    cat(sprintf(
      "  [SKIP] Feature [%s] non presente in questo campione.\n",
      car_feature))
    cat("  Features disponibili che contengono 'CAR':\n")
    cands <- rownames(obj)[grepl("CAR|car", rownames(obj),
                                  ignore.case = TRUE)]
    if (length(cands) > 0)
      cat(paste0("    ", cands, "\n"))
    else
      cat("    Nessuna.\n")
    return(NULL)
  }

  # ── Estrai espressione del costrutto ─────────────────────
  v_car <- as.numeric(
    GetAssayData(obj, slot = "data")[car_feature, ])

  n_expressing <- sum(v_car > 0)
  pct_expr     <- n_expressing / ncol(obj) * 100
  cat(sprintf(
    "  Feature [%s]: %d cellule esprimenti (%.1f%%)\n",
    car_feature, n_expressing, pct_expr))
  cat(sprintf(
    "  Mediana (esprimenti): %.3f | Max: %.3f\n",
    if (n_expressing > 0) median(v_car[v_car > 0]) else 0,
    max(v_car)))

  # ── Soglia e classificazione ──────────────────────────────
  thr <- calc_threshold(v_car, pct)
  cat(sprintf(
    "  Soglia %d° pct (cellule esprimenti): %.4f\n",
    pct, thr))

  is_car_expr <- v_car > thr
  n_car_expr  <- sum(is_car_expr)
  cat(sprintf(
    "  CAR+ (expr > soglia): %d cellule (%.2f%%)\n",
    n_car_expr, n_car_expr / ncol(obj) * 100))

  if (n_car_expr == 0) {
    cat(paste0(
      "  [WARN] Nessuna cellula sopra soglia.\n",
      "  Possibili cause:\n",
      "  ─ Tutti i reads del costrutto hanno espressione 0\n",
      "    → il feature mapping non ha funzionato correttamente\n",
      "  ─ La feature è presente ma senza reads mappati\n",
      "  → Verifica l'allineamento CellRanger su questo campione.\n"))
    return(NULL)
  }

  # ── Metadati ──────────────────────────────────────────────
  meta              <- obj@meta.data
  meta$cell_type    <- as.character(meta$cell_type)
  meta$expr_CAR     <- v_car
  meta$CAR_expr     <- ifelse(is_car_expr, "CAR+", "CAR-")

  has_screp <- "IS_CAR_ALLIN_scREP" %in% colnames(meta)
  if (has_screp) {
    meta$CAR_scREP <- ifelse(
      meta$IS_CAR_ALLIN_scREP == "YES", "CAR+", "CAR-")
    n_screp <- sum(meta$CAR_scREP == "CAR+")
    cat(sprintf(
      "  CAR+ scREP: %d cellule (%.2f%%)\n",
      n_screp, n_screp / ncol(obj) * 100))
  }

  pop_order <- sort(unique(meta$cell_type))
  pop_cols  <- setNames(sapply(pop_order, get_color), pop_order)
  Idents(obj) <- "cell_type"

  # ── 1. VlnPlot feature CAR per cluster ───────────────────
  # Mostra quanti reads del costrutto ha ogni cluster.
  # Atteso nei campioni I (infusione): quasi tutti i cluster
  # T avranno qualche cellula positiva se il prodotto è ricco.
  # Nei campioni AB: solo i cluster con CAR-T persistenti.

  obj$expr_CAR <- v_car
  p_vln <- VlnPlot(obj,
                   features = "expr_CAR",
                   group.by = "cell_type",
                   cols     = pop_cols,
                   pt.size  = 0.4) +
    geom_hline(yintercept = thr,
               linetype = "dashed", color = "#B00020",
               linewidth = 0.9) +
    annotate("text",
             x     = length(pop_order) * 0.97,
             y     = thr * 1.12,
             label = paste0(pct, "° pct\n(soglia)"),
             color = "#B00020", size = 3,
             hjust = 1, fontface = "italic") +
    labs(
      title    = paste0(sample_name,
                        " – Espressione feature [", car_feature,
                        "] per cluster"),
      subtitle = paste0(
        "Reads del costrutto CAR per cellula per cluster\n",
        "Linea = soglia ", pct,
        "° percentile (calcolata su cellule esprimenti, n=",
        n_expressing, ")"),
      y = "Espressione (log-norm)", x = NULL) +
    theme_classic(base_size = 10) +
    theme(
      axis.text.x   = element_text(angle = 45, hjust = 1,
                                    size = 8, face = "bold"),
      plot.title    = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = 8,
                                    color = "gray40"),
      legend.position = "none")

  vln_w <- max(12, length(pop_order) * 0.7 + 3)
  ggsave(paste0(out_dir, sample_name, "_vln_CAR_construct.png"),
         plot = p_vln, width = vln_w, height = 6,
         dpi = 300, bg = "white")
  cat(paste0("  → ", sample_name, "_vln_CAR_construct.png\n"))

  # ── 2. Density plot ───────────────────────────────────────
  # Distribuzione dell'espressione del costrutto su tutte le
  # cellule. Una bimodalità chiara (picco a 0 + coda/secondo
  # picco a destra) indica buona separazione CAR vs non-CAR
  # e valida l'approccio. Una distribuzione unimodale indica
  # che la feature ha reads dispersi → soglia meno affidabile.

  df_dens <- data.frame(expr = v_car)
  p_dens <- ggplot(df_dens, aes(x = expr)) +
    # Mostra solo le cellule esprimenti per vedere meglio
    # la distribuzione reale (escludiamo la massa di zeri)
    geom_density(
      data = df_dens[df_dens$expr > 0, , drop = FALSE],
      fill = "#FFCCD5", alpha = 0.7,
      color = "#B00020", linewidth = 0.8) +
    geom_vline(xintercept = thr,
               color = "#B00020", linetype = "dashed",
               linewidth = 0.9) +
    annotate("text",
             x = thr, y = Inf, vjust = 1.5, hjust = -0.15,
             label = paste0(pct, "° pct\n(n CAR+ = ",
                            n_car_expr, ")"),
             color = "#B00020", size = 3.5) +
    labs(
      title    = paste0(sample_name,
                        " – Distribuzione [", car_feature, "]"),
      subtitle = paste0(
        "Solo cellule con expr > 0 (n = ", n_expressing, ")\n",
        "Bimodalità = buona separazione CAR vs non-CAR"),
      x = "Espressione del costrutto (log-norm)",
      y = "Densità") +
    theme_classic(base_size = 11) +
    theme(
      plot.title    = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = 8.5,
                                    color = "gray40"))

  ggsave(paste0(out_dir, sample_name,
                "_density_CAR_construct.png"),
         plot = p_dens, width = 6, height = 5,
         dpi = 300, bg = "white")
  cat(paste0("  → ", sample_name,
             "_density_CAR_construct.png\n"))

  # ── 3. UMAP con 4 categorie di confronto ─────────────────
  # Overlap      = entrambi concordano → CAR-T certe
  # Solo scREP   = VDJ+ ma 0 reads costrutto → plausibile
  #                se dropout di reads del transgene
  # Solo expr    = reads costrutto ma VDJ dropout → MANCATE
  #                da scREP; queste sono le cellule da recuperare
  # CAR- entrambi = non CAR

  umap_key <- find_umap(obj)

  if (!is.null(umap_key)) {
    coords <- as.data.frame(Embeddings(obj, umap_key)[, 1:2])
    colnames(coords) <- c("UMAP1","UMAP2")

    if (has_screp) {
      coords$category <- case_when(
        is_car_expr & meta$CAR_scREP == "CAR+" ~
          "Overlap (expr + scREP)",
        is_car_expr & meta$CAR_scREP == "CAR-" ~
          "Solo expr (perse da scREP?)",
        !is_car_expr & meta$CAR_scREP == "CAR+" ~
          "Solo scREP (0 reads costrutto)",
        TRUE ~ "CAR-")
    } else {
      coords$category <- ifelse(is_car_expr, "CAR+ expr", "CAR-")
    }

    cat_levels_umap <- c(
      "CAR-",
      "Solo scREP (0 reads costrutto)",
      "Solo expr (perse da scREP?)",
      "Overlap (expr + scREP)",
      "CAR+ expr")
    coords$category <- factor(
      coords$category,
      levels = cat_levels_umap[
        cat_levels_umap %in% coords$category])
    coords <- coords[order(coords$category), ]

    cat_colors_umap <- c(
      "CAR-"                          = "#DDDDDD",
      "Solo scREP (0 reads costrutto)"= "#4361EE",
      "Solo expr (perse da scREP?)"   = "#F4A261",
      "Overlap (expr + scREP)"        = "#B00020",
      "CAR+ expr"                     = "#B00020")
    cat_sizes_umap <- c(
      "CAR-"                          = 0.3,
      "Solo scREP (0 reads costrutto)"= 1.2,
      "Solo expr (perse da scREP?)"   = 1.3,
      "Overlap (expr + scREP)"        = 1.5,
      "CAR+ expr"                     = 1.5)
    cat_alpha_umap <- c(
      "CAR-"                          = 0.25,
      "Solo scREP (0 reads costrutto)"= 0.85,
      "Solo expr (perse da scREP?)"   = 0.90,
      "Overlap (expr + scREP)"        = 1.0,
      "CAR+ expr"                     = 1.0)

    centroids <- data.frame(
      cell_type = meta$cell_type,
      UMAP1 = coords$UMAP1,
      UMAP2 = coords$UMAP2) %>%
      group_by(cell_type) %>%
      summarise(UMAP1 = median(UMAP1),
                UMAP2 = median(UMAP2),
                .groups = "drop")

    n_overlap    <- sum(coords$category ==
                          "Overlap (expr + scREP)")
    n_solo_screp <- sum(coords$category ==
                          "Solo scREP (0 reads costrutto)")
    n_solo_expr  <- sum(coords$category ==
                          "Solo expr (perse da scREP?)")

    subtitle_u <- if (has_screp) paste0(
      "Rosso = overlap (", n_overlap,
      ")  |  Blu = solo scREP (", n_solo_screp,
      ") = 0 reads costrutto\n",
      "Arancio = solo expr (", n_solo_expr,
      ") = probabilmente perse da scREP per VDJ dropout") else
      paste0("CAR+ da espressione costrutto: ", n_car_expr)

    lev_umap <- levels(coords$category)

    p_umap <- ggplot(coords,
      aes(x = UMAP1, y = UMAP2,
          color = category,
          size  = category,
          alpha = category)) +
      geom_point() +
      scale_color_manual(
        values = cat_colors_umap[lev_umap], name = NULL) +
      scale_size_manual(
        values = cat_sizes_umap[lev_umap], guide = "none") +
      scale_alpha_manual(
        values = cat_alpha_umap[lev_umap], guide = "none") +
      ggrepel::geom_label_repel(
        data = centroids,
        aes(x = UMAP1, y = UMAP2, label = cell_type),
        inherit.aes   = FALSE,
        size          = 2.8, fontface = "bold",
        fill          = alpha("white", 0.65), color = "black",
        label.size    = 0.12,
        label.padding = unit(0.1, "lines"),
        max.overlaps  = 25, seed = 42) +
      labs(
        title    = paste0(sample_name,
                          " – Feature [", car_feature,
                          "] su UMAP"),
        subtitle = subtitle_u) +
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
                  "_UMAP_CAR_construct.png"),
           plot = p_umap, width = 9, height = 8,
           dpi = 300, bg = "white")
    cat(paste0("  → ", sample_name,
               "_UMAP_CAR_construct.png\n"))

    # ── Extra: UMAP barplot cellule perse da scREP per cluster
    # Mostra in quali cluster si concentrano le cellule
    # "Solo expr" = candidate CAR-T perse da scREP.
    if (has_screp && n_solo_expr > 0) {
      df_solo <- data.frame(
        cell_type = meta$cell_type[
          coords$category == "Solo expr (perse da scREP?)"])

      df_solo_counts <- df_solo %>%
        count(cell_type, name = "n") %>%
        arrange(desc(n))

      bar_cols <- setNames(
        sapply(df_solo_counts$cell_type, get_color),
        df_solo_counts$cell_type)

      p_bar_solo <- ggplot(df_solo_counts,
        aes(x = reorder(cell_type, n),
            y = n, fill = cell_type)) +
        geom_col(show.legend = FALSE) +
        geom_text(aes(label = n),
                  hjust = -0.2, size = 3.5,
                  fontface = "bold") +
        scale_fill_manual(values = bar_cols) +
        coord_flip() +
        expand_limits(y = max(df_solo_counts$n) * 1.15) +
        labs(
          title    = paste0(sample_name,
                            " – Cellule 'Solo expr'",
                            " per cluster"),
          subtitle = paste0(
            "Reads del costrutto presenti ma non in scREP\n",
            "(probabile VDJ dropout: n totale = ",
            n_solo_expr, ")"),
          x = NULL, y = "n cellule") +
        theme_classic(base_size = 11) +
        theme(
          plot.title    = element_text(face = "bold",
                                        hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, size = 8.5,
                                        color = "gray40"))

      ggsave(paste0(out_dir, sample_name,
                    "_barplot_solo_expr.png"),
             plot = p_bar_solo,
             width = 8,
             height = max(4, nrow(df_solo_counts) * 0.5 + 2),
             dpi = 300, bg = "white")
      cat(paste0("  → ", sample_name,
                 "_barplot_solo_expr.png\n"))
    }
  }

  # ── 4. Concordanza quantitativa ───────────────────────────
  concordance_row <- data.frame(
    campione           = sample_name,
    car_feature        = car_feature,
    percentile         = pct,
    n_totale_cellule   = ncol(obj),
    n_CAR_expr         = n_car_expr,
    pct_CAR_expr       = round(n_car_expr / ncol(obj) * 100, 2),
    stringsAsFactors   = FALSE)

  if (has_screp) {
    is_screp_pos <- meta$CAR_scREP == "CAR+"
    n_screp_pos  <- sum(is_screp_pos)
    overlap      <- sum(is_car_expr &  is_screp_pos)
    solo_screp   <- sum(!is_car_expr & is_screp_pos)
    solo_expr    <- sum(is_car_expr  & !is_screp_pos)
    true_neg     <- sum(!is_car_expr & !is_screp_pos)
    sensitivity  <- if (n_screp_pos > 0)
      round(overlap / n_screp_pos * 100, 1) else NA
    specificity  <- if ((ncol(obj) - n_screp_pos) > 0)
      round(true_neg / (ncol(obj) - n_screp_pos) * 100, 1) else NA

    concordance_row$n_CAR_scREP    <- n_screp_pos
    concordance_row$overlap         <- overlap
    concordance_row$solo_scREP      <- solo_screp
    concordance_row$solo_expr       <- solo_expr
    concordance_row$sensitivity_pct <- sensitivity
    concordance_row$specificity_pct <- specificity

    cat(sprintf("  Overlap (expr ∩ scREP):             %d\n",
                overlap))
    cat(sprintf("  Solo scREP (0 reads costrutto):     %d\n",
                solo_screp))
    cat(sprintf("  Solo expr (perse da scREP?):        %d\n",
                solo_expr))
    cat(sprintf("  Sensibilità vs scREP: %.1f%%\n", sensitivity))
    cat(sprintf("  Specificità vs scREP: %.1f%%\n", specificity))
  }

  return(list(
    sample      = sample_name,
    meta        = meta,
    concordance = concordance_row,
    thr         = thr,
    n_car_expr  = n_car_expr
  ))
}

# ============================================================
# LOOP PRINCIPALE
# ============================================================
section(paste0("Analisi feature [", CAR_FEATURE_NAME, "]"))

results <- list()
for (nm in names(all_samples)) {
  results[[nm]] <- analyze_car_construct(
    obj         = all_samples[[nm]],
    sample_name = nm,
    out_dir     = out_dir,
    car_feature = CAR_FEATURE_NAME,
    pct         = PERCENTILE)
}
results <- Filter(Negate(is.null), results)

if (length(results) == 0) {
  cat(paste0(
    "\n[STOP] Nessun campione ha prodotto risultati.\n",
    "Verifica:\n",
    "1. CAR_FEATURE_NAME è corretto?\n",
    "2. Il genoma custom includeva il costrutto come feature?\n",
    "3. I reads del costrutto sono stati allineati correttamente?\n"))
  stop("Nessun risultato prodotto.")
}

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
  data.frame(campione    = nm,
             car_feature = CAR_FEATURE_NAME,
             percentile  = PERCENTILE,
             soglia      = round(r$thr, 4),
             stringsAsFactors = FALSE)
}))
addWorksheet(wb, "Soglie")
writeData(wb, "Soglie", thresh_df)

for (nm in names(results)) {
  r        <- results[[nm]]
  sheet_nm <- substr(nm, 1, 31)
  addWorksheet(wb, sheet_nm)
  cols_keep <- c("cell_type", "expr_CAR", "CAR_expr",
                 if ("CAR_scREP" %in% colnames(r$meta))
                   "CAR_scREP" else NULL)
  writeData(wb, sheet_nm,
            r$meta[, cols_keep[cols_keep %in%
                                  colnames(r$meta)]])
}

xlsx_path <- paste0(out_dir,
                    "CAR_construct_detection.xlsx")
saveWorkbook(wb, xlsx_path, overwrite = TRUE)
cat(paste0("  → ", xlsx_path, "\n"))

# ============================================================
# RIEPILOGO FINALE
# ============================================================
section("Riepilogo finale")

cat("\nRisultati per campione:\n\n")
has_screp_global <- "sensitivity_pct" %in%
                    colnames(all_concordance)

if (has_screp_global) {
  cat(sprintf(
    "  %-18s | %7s | %7s | %9s | %9s | %6s | %6s\n",
    "Campione","n_expr","n_scREP",
    "Overlap","SoloSCREP","Sens%","Spec%"))
  cat(paste0("  ", strrep("-", 76), "\n"))
  for (i in seq_len(nrow(all_concordance))) {
    r <- all_concordance[i, ]
    cat(sprintf(
      "  %-18s | %7d | %7d | %9d | %9d | %5.1f%% | %5.1f%%\n",
      r$campione, r$n_CAR_expr, r$n_CAR_scREP,
      r$overlap, r$solo_scREP,
      r$sensitivity_pct, r$specificity_pct))
  }
} else {
  for (i in seq_len(nrow(all_concordance))) {
    r <- all_concordance[i, ]
    cat(sprintf("  %-18s : %d CAR+ (%.2f%%)\n",
                r$campione, r$n_CAR_expr, r$pct_CAR_expr))
  }
}

cat(paste0(
  "\n", strrep("=", 65), "\n",
  "  ANALISI COMPLETATA\n\n",
  "  Feature usata: [", CAR_FEATURE_NAME, "]\n",
  "  Output: ", out_dir, "\n\n",
  "  Per campione:\n",
  "    _vln_CAR_construct.png       reads costrutto per cluster\n",
  "    _density_CAR_construct.png   distribuzione + soglia\n",
  "    _UMAP_CAR_construct.png      UMAP con 4 categorie\n",
  "    _barplot_solo_expr.png       cluster con cellule perse\n",
  "  Globale:\n",
  "    CAR_construct_detection.xlsx\n\n",
  "  INTERPRETAZIONE:\n",
  "  ─ Overlap alto       → metodi concordano, ottimo\n",
  "  ─ Solo scREP alto    → quelle cellule CAR+ non hanno\n",
  "                         reads della feature costrutto:\n",
  "                         normale per dropout tecnico\n",
  "  ─ Solo expr alto     → CELLULE PROBABILMENTE MANCATE\n",
  "                         DA scREP per VDJ dropout:\n",
  "                         controlla _barplot_solo_expr.png\n",
  "                         per vedere in quali cluster sono\n",
  "  ─ Se overlap ~ 0     → la feature potrebbe non essere\n",
  "                         quella corretta, oppure i reads\n",
  "                         del costrutto non si allineano\n",
  "                         a questa feature nel tuo gtf\n",
  strrep("=", 65), "\n"))
