# ============================================================
#  Q2: Caratterizzazione trascrizionale delle cellule CAR+
#      nel prodotto di infusione (I)
#
#  Domanda biologica:
#    Quali caratteristiche funzionali hanno le CART cells
#    nel prodotto di infusione?
#    - Sono effettrici? memory? esauste? naïve?
#    - In quali sottotipi cellulari si trovano principalmente?
#    - Come differiscono le CAR+ dalle CAR- nello stesso prodotto?
#
#  Approccio:
#    A) Distribuzione dei tipi cellulari nelle CAR+ vs CAR- in I
#    B) Module scores di programmi funzionali (effector, memory,
#       exhaustion, stemness, proliferazione)
#    C) DEG: CAR+ vs CAR- per paziente (FindMarkers, Wilcoxon)
#       NOTA: analisi per paziente individuale, non pooling;
#       pooling richiederebbe pseudobulk con n=3 che ha
#       potenza statistica molto limitata
#    D) DotPlot top marker per tipo cellulare × stato CAR
#
#  Prerequisiti:
#    all_samples_annotated_COMPLETE_IS_CAR_REVISED.rds
#
#  Output in out_dir/Q2_CART_in_I/:
#    Q2_<paziente>_I_UMAP_CAR_overlay.png
#    Q2_<paziente>_I_module_scores.png
#    Q2_<paziente>_I_celltype_distribution.png
#    Q2_ALL_I_DEG_CARpos_vs_CARneg.xlsx
#    Q2_ALL_I_module_score_summary.png
# ============================================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(scales)
library(openxlsx)
library(ggrepel)

# ── UNICO PUNTO DA MODIFICARE ────────────────────────────────
rds_path <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/CART/Code/Seurat_analysis/2_annotation/all_samples_annotated_COMPLETE_IS_CAR_REVISED.rds"
out_dir  <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/CART/Code/Seurat_analysis/10_CART_functional_analysis/Q2_CART_in_I/"
# ─────────────────────────────────────────────────────────────

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

section <- function(title)
  cat(paste0("\n", strrep("=", 65), "\n  ", title, "\n",
             strrep("=", 65), "\n"))

# ============================================================
# FIRME GENICHE PER I PROGRAMMI FUNZIONALI
# ============================================================
# Riferimenti: Gattinoni et al. (Tscm), Wherry (exhaustion),
# Blackburn (exhaustion), Philip et al., Sade-Feldman et al.
# Tutte le firme includono solo geni comuni in scRNA-seq
# (si evitano citochine tipo IFNG/IL2/IL17A che sono basse
#  a riposo - dropout frequente in 10x Chromium).

SIGNATURES <- list(

  # Effettore/Citotossicità
  # Geni robusti anche in scRNA-seq a riposo
  Effector = c(
    "GZMB",   # Granzyme B – citotossicità diretta
    "PRF1",   # Perforina
    "NKG7",   # Marker effettore CD8
    "GNLY",   # Granulysin
    "GZMA",   # Granzyme A
    "GZMK",   # Granzyme K – effettore early/EM
    "FGFBP2", # Marker CD8 effettore terminale
    "CX3CR1"  # Traffico effettore periferico
  ),

  # Memoria / Stemness
  # TCF7 (TCF1) è il marker chiave delle Tscm/Tcm
  # IL7R importante per survival
  Memory_Stemness = c(
    "TCF7",  # TCF1 – master regulator stemness/memory
    "CCR7",  # Homing linfonodali – Tnaive e Tcm
    "SELL",  # CD62L – Tnaive e Tcm
    "IL7R",  # CD127 – survival e memory
    "LEF1",  # Wnt/β-catenin – stemness
    "KLF2",  # Mantiene fenotipo naive/Tcm
    "BCL2",  # Sopravvivenza cellule memory
    "FOXO1"  # Regolatore transcrizionale memory
  ),

  # Esaurimento (Exhaustion)
  # Geni inibitori co-espressi nelle cellule esauste
  Exhaustion = c(
    "PDCD1",   # PD-1 – checkpoint inibizione
    "LAG3",    # LAG3 – inibizione citotossicità
    "HAVCR2",  # TIM-3 – terminale esaurimento
    "TIGIT",   # Inibizione NK e T
    "TOX",     # Master regolatore esaurimento
    "TOX2",    # Co-regulatore
    "ENTPD1",  # CD39 – marker esaurimento
    "CTLA4",   # CTLA4 – inibizione precoce
    "BATF"     # TF esaurimento (down in effettori, up in esausti)
  ),

  # Attivazione precoce / effettore "early"
  # Distingue cellule attivate da quelle naïve
  Activation = c(
    "CD69",    # Marker attivazione precoce
    "CD44",    # Marker esperienza antigenica
    "TNFRSF9", # 4-1BB – attivazione T
    "IL2RA",   # CD25 – attivazione IL-2
    "ICOS",    # Costimolazione
    "CD38"     # Attivazione e differenziazione
  ),

  # Proliferazione
  Proliferation = c(
    "MKI67",  # Ki67 – marker proliferazione principale
    "TOP2A",  # Fase G2/M
    "PCNA",   # Fase S
    "CCNB1",  # Ciclina B1
    "STMN1",  # Marker proliferazione
    "UBE2C"   # Ubiquitina E2 – ciclo cellulare
  ),

  # Fenotipo naïve (per valutare immaturità)
  Naive = c(
    "CCR7", "SELL", "IL7R", "TCF7", "LEF1",
    "KLF2", "LDHB", "RCAN3"
  ),

  # Stem-like / progenitor exhausted (Tpex)
  # Importante in contesti immunoterapia
  Tpex_StemLike = c(
    "TCF7",    # Richiesto per Tpex
    "CXCR5",   # Marker Tpex follicolare
    "TOX",     # Co-espresso in Tpex
    "BCL6",    # Regolatore Tpex
    "SLAMF6",  # Marker Tpex
    "ID3"      # Stemness T
  ),

  # Effettore terminale / esausto terminale (Tex)
  Tex_Terminal = c(
    "HAVCR2",  # TIM-3 massimo in Tex terminale
    "TIGIT",
    "LAG3",
    "CD160",   # HVEM ligand – terminale
    "ENTPD1",  # CD39
    "PRDM1",   # BLIMP1 – differenziazione terminale
    "ZEB2"     # Tex terminale
  )
)

# ============================================================
# CARICAMENTO E ESTRAZIONE CAMPIONI I
# ============================================================
section("Caricamento dati")

all_samples <- readRDS(rds_path)

# Seleziona campioni I
I_sample_names <- grep("_I$", names(all_samples), value = TRUE)
cat("Campioni I trovati:", paste(I_sample_names, collapse = ", "), "\n")

if (length(I_sample_names) == 0)
  stop("Nessun campione I trovato. Verifica i nomi nel file RDS.")

# ── Helper: colonna CAR ─────────────────────────────────────
get_car_status <- function(obj, sample_name) {
  meta <- obj@meta.data
  for (col in c("IS_CAR_ALLIN_scREP", "IS_CAR", "CAR")) {
    if (col %in% colnames(meta)) {
      vals    <- as.character(meta[[col]])
      car_pos <- grepl("^(YES|TRUE|yes|true|1)$", vals)
      cat(sprintf("  %s: '%s' | CAR+ = %d / %d (%.1f%%)\n",
                  sample_name, col,
                  sum(car_pos), length(car_pos),
                  100 * mean(car_pos)))
      return(ifelse(car_pos, "CAR+", "CAR-"))
    }
  }
  cat(sprintf("[WARN] %s: nessuna colonna CAR\n", sample_name))
  rep("CAR-", ncol(obj))
}

# ── Helper: verifica UMAP disponibile ──────────────────────
get_umap_key <- function(obj) {
  for (nm in c("umap", "wnn.umap", "umap.harmony", "RNA.umap"))
    if (nm %in% names(obj@reductions)) return(nm)
  return(NULL)
}

# ── Helper: filtra geni presenti nell'oggetto ───────────────
filter_genes <- function(genes, obj) {
  genes[genes %in% rownames(obj)]
}

# ============================================================
# STEP 1: AGGIUNGI CAR STATUS E MODULE SCORES
# ============================================================
section("STEP 1 | CAR status + Module scores in I")

I_objects <- list()

for (nm in I_sample_names) {
  cat(paste0("\n── ", nm, " ──\n"))
  obj <- all_samples[[nm]]

  # CAR status
  obj$CAR_status <- get_car_status(obj, nm)

  # Verifica livello minimo di CAR+ per analisi
  n_car_pos <- sum(obj$CAR_status == "CAR+")
  cat(sprintf("  CAR+ totali: %d\n", n_car_pos))
  if (n_car_pos < 5)
    cat("  [ATTENZIONE] Pochissime cellule CAR+. Risultati DEG non affidabili.\n")

  # JoinLayers se necessario (Seurat v5)
  if (length(grep("^counts\\.", Layers(obj), value = TRUE)) > 0)
    obj <- JoinLayers(obj)

  # ── Module scores ────────────────────────────────────────
  for (sig_name in names(SIGNATURES)) {
    genes_ok <- filter_genes(SIGNATURES[[sig_name]], obj)
    if (length(genes_ok) < 3) {
      cat(sprintf("  [SKIP] Firma '%s': solo %d geni trovati (min 3)\n",
                  sig_name, length(genes_ok)))
      next
    }
    col_nm <- paste0("Score_", sig_name)
    obj    <- AddModuleScore(obj, features = list(genes_ok),
                             name = col_nm, seed = 42)
    obj[[col_nm]]              <- obj[[paste0(col_nm, "1")]]
    obj[[paste0(col_nm, "1")]] <- NULL
    cat(sprintf("  Score_%s calcolato (%d geni)\n",
                sig_name, length(genes_ok)))
  }

  I_objects[[nm]] <- obj
}

# ============================================================
# STEP 2: PLOT UMAP CON OVERLAY CAR + MODULE SCORES
# ============================================================
section("STEP 2 | UMAP overlay (CAR + module scores)")

# Colori per stato CAR
CAR_COLORS <- c("CAR+" = "#E63946", "CAR-" = "#ADB5BD")

for (nm in names(I_objects)) {
  cat(paste0("\n── UMAP plot: ", nm, " ──\n"))
  obj     <- I_objects[[nm]]
  umap_key <- get_umap_key(obj)

  if (is.null(umap_key)) {
    cat("  [SKIP] Nessuna riduzione UMAP trovata\n")
    next
  }

  Idents(obj) <- "cell_type"

  # A) UMAP colorato per tipo cellulare + CAR overlay
  # Ordina le cellule: CAR- sotto, CAR+ sopra
  coords <- as.data.frame(Embeddings(obj, umap_key)[, 1:2])
  colnames(coords) <- c("UMAP1", "UMAP2")
  coords$car       <- obj$CAR_status
  coords$cell_type <- as.character(obj$cell_type)
  coords$sort_key  <- ifelse(coords$car == "CAR+", 1, 0)
  coords <- coords[order(coords$sort_key), ]

  centroids <- coords %>%
    group_by(cell_type) %>%
    summarise(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2), .groups = "drop")

  n_car_pos <- sum(coords$car == "CAR+")
  n_car_neg <- sum(coords$car == "CAR-")

  p_umap_car <- ggplot(coords, aes(x = UMAP1, y = UMAP2)) +
    geom_point(data = coords[coords$car == "CAR-",],
               color = "#D3D3D3", size = 0.4, alpha = 0.5) +
    geom_point(data = coords[coords$car == "CAR+",],
               color = "white", size = 2.5, alpha = 0.9) +
    geom_point(data = coords[coords$car == "CAR+",],
               color = "#E63946", size = 1.6, alpha = 0.9) +
    geom_label_repel(
      data = centroids,
      aes(x = UMAP1, y = UMAP2, label = cell_type),
      size = 3, fontface = "bold",
      fill = scales::alpha("white", 0.7), color = "black",
      label.size = 0.2, max.overlaps = 20, seed = 42
    ) +
    ggtitle(
      paste0(nm, " – CAR+ overlay su UMAP"),
      subtitle = paste0("Rosso: CAR+ (n=", n_car_pos, ")  |  ",
                        "Grigio: CAR- (n=", n_car_neg, ")")
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title    = element_text(face = "bold", hjust = 0.5, size = 12),
      plot.subtitle = element_text(hjust = 0.5, color = "gray40", size = 9),
      axis.text     = element_blank(), axis.ticks = element_blank()
    )

  # B) FeaturePlot dei module scores principali
  score_cols <- grep("^Score_", colnames(obj@meta.data), value = TRUE)
  # Mostra i 4 più rilevanti (effettore, memory, exhaustion, proliferazione)
  key_scores <- intersect(
    c("Score_Effector", "Score_Memory_Stemness", "Score_Exhaustion",
      "Score_Proliferation", "Score_Naive", "Score_Tpex_StemLike"),
    score_cols
  )

  if (length(key_scores) >= 2) {
    fp_list <- lapply(key_scores, function(sc) {
      FeaturePlot(obj, features = sc, reduction = umap_key,
                  pt.size = 0.4, order = TRUE,
                  min.cutoff = "q05", max.cutoff = "q95") +
        scale_color_gradientn(
          colors = c("lightgrey", "#FFF176", "#FB8C00", "#B71C1C"),
          name   = gsub("Score_", "", sc)
        ) +
        ggtitle(gsub("Score_", "", sc)) +
        theme_classic(base_size = 9) +
        theme(
          plot.title  = element_text(face = "bold", hjust = 0.5, size = 9),
          axis.text   = element_blank(), axis.ticks = element_blank(),
          legend.text = element_text(size = 7)
        )
    })

    p_scores <- wrap_plots(fp_list, ncol = 3)
  } else {
    p_scores <- NULL
  }

  # Combina e salva
  combined <- if (!is.null(p_scores))
    (p_umap_car / p_scores) +
      plot_annotation(title = paste0(nm, " – Module Scores CAR+"),
                      theme = theme(plot.title = element_text(face = "bold")))
  else
    p_umap_car

  out_path <- paste0(out_dir, "Q2_", nm, "_UMAP_CAR_modules.png")
  ggsave(out_path, combined,
         width = 16, height = if (!is.null(p_scores)) 14 else 7,
         dpi = 300, bg = "white")
  cat(paste0("  → ", out_path, "\n"))
}

# ============================================================
# STEP 3: DISTRIBUZIONE CELLULARE CAR+ VS CAR- IN I
# ============================================================
section("STEP 3 | Distribuzione tipo cellulare: CAR+ vs CAR- in I")

PALETTE_CT <- c(
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
  "Myeloid cells"              = "#E76F51"
)

bar_plots_ct <- lapply(names(I_objects), function(nm) {
  obj <- I_objects[[nm]]
  df  <- data.frame(
    cell_type  = as.character(obj$cell_type),
    car        = obj$CAR_status,
    stringsAsFactors = FALSE
  ) %>%
    group_by(car, cell_type) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(car) %>%
    mutate(pct = 100 * n / sum(n)) %>%
    ungroup() %>%
    mutate(car = factor(car, levels = c("CAR+", "CAR-")))

  types_present <- unique(df$cell_type)
  cols          <- PALETTE_CT[types_present]
  cols[is.na(cols)] <- "grey60"

  ggplot(df, aes(x = car, y = pct, fill = cell_type)) +
    geom_col(position = "stack", width = 0.7, color = "white",
             linewidth = 0.2) +
    scale_fill_manual(values = cols, name = NULL) +
    scale_y_continuous(labels = function(x) paste0(x, "%"),
                       limits = c(0, 105), expand = c(0, 0)) +
    labs(title = nm,
         x = "Stato CAR", y = "% cellule") +
    theme_classic(base_size = 10) +
    theme(
      plot.title    = element_text(face = "bold", hjust = 0.5, size = 11),
      axis.text.x   = element_text(size = 10),
      legend.text   = element_text(size = 7),
      legend.key.size = unit(0.3, "cm")
    )
})

p_ct_combined <- wrap_plots(bar_plots_ct, nrow = 1) +
  plot_annotation(
    title    = "Distribuzione tipi cellulari: CAR+ vs CAR- nel prodotto di infusione (I)",
    subtitle = "Proporzione di ogni tipo cellulare dentro ciascun gruppo (CAR+/CAR-)",
    theme    = theme(plot.title = element_text(face = "bold",
                                               hjust = 0.5, size = 13))
  )

out_ct <- paste0(out_dir, "Q2_ALL_I_celltype_distribution_CAR.png")
ggsave(out_ct, p_ct_combined,
       width = 5 * length(I_objects) + 2, height = 7,
       dpi = 300, bg = "white")
cat(paste0("  → ", out_ct, "\n"))

# ============================================================
# STEP 4: VIOLIN PLOT MODULE SCORES: CAR+ vs CAR- in I
# ============================================================
section("STEP 4 | VlnPlot module scores: CAR+ vs CAR-")

score_violin_list <- lapply(names(I_objects), function(nm) {
  obj        <- I_objects[[nm]]
  score_cols <- grep("^Score_", colnames(obj@meta.data), value = TRUE)

  if (length(score_cols) == 0) return(NULL)

  # Dati per violin
  df_vln <- obj@meta.data %>%
    select(all_of(score_cols), CAR_status) %>%
    mutate(CAR_status = factor(CAR_status, levels = c("CAR+", "CAR-"))) %>%
    pivot_longer(cols = all_of(score_cols),
                 names_to = "Signature",
                 values_to = "Score") %>%
    mutate(Signature = gsub("^Score_", "", Signature))

  # Test Mann-Whitney per ogni firma
  sig_names <- unique(df_vln$Signature)
  pvals <- sapply(sig_names, function(s) {
    d_pos <- df_vln$Score[df_vln$Signature == s & df_vln$CAR_status == "CAR+"]
    d_neg <- df_vln$Score[df_vln$Signature == s & df_vln$CAR_status == "CAR-"]
    if (length(d_pos) < 3 || length(d_neg) < 3) return(NA)
    tryCatch(
      wilcox.test(d_pos, d_neg, exact = FALSE)$p.value,
      error = function(e) NA
    )
  })

  cat(paste0("\n  [", nm, "] Mann-Whitney p-values (CAR+ vs CAR-):\n"))
  pval_df <- data.frame(
    signature = sig_names,
    p_value   = round(pvals, 4),
    sig       = ifelse(pvals < 0.05, "*", "ns"),
    stringsAsFactors = FALSE
  )
  print(pval_df, row.names = FALSE)
  cat("  NOTA: test su singolo campione - solo valore descrittivo.\n")

  ggplot(df_vln, aes(x = CAR_status, y = Score, fill = CAR_status)) +
    geom_violin(trim = TRUE, alpha = 0.7, scale = "width") +
    geom_boxplot(width = 0.12, fill = "white", outlier.size = 0.3,
                 outlier.alpha = 0.3) +
    scale_fill_manual(values = c("CAR+" = "#E63946", "CAR-" = "#ADB5BD"),
                      guide = "none") +
    facet_wrap(~ Signature, scales = "free_y", nrow = 2) +
    labs(
      title    = paste0(nm, " – Module scores: CAR+ vs CAR-"),
      subtitle = "p-value (Mann-Whitney) indicato nei dati; N=1 campione, solo descrittivo",
      x = NULL, y = "Module score"
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title    = element_text(face = "bold", hjust = 0.5, size = 12),
      plot.subtitle = element_text(hjust = 0.5, color = "gray40", size = 8),
      strip.text    = element_text(face = "bold", size = 9),
      axis.text.x   = element_text(size = 10)
    )
})

score_violin_list <- Filter(Negate(is.null), score_violin_list)
if (length(score_violin_list) > 0) {
  p_vln_all <- wrap_plots(score_violin_list, ncol = 1) +
    plot_annotation(
      title = "Module scores CAR+ vs CAR- – tutti i campioni I",
      theme = theme(plot.title = element_text(face = "bold", hjust = 0.5))
    )
  out_vln <- paste0(out_dir, "Q2_ALL_I_module_scores_violin.png")
  ggsave(out_vln, p_vln_all,
         width = 14, height = 8 * length(score_violin_list),
         dpi = 300, bg = "white", limitsize = FALSE)
  cat(paste0("  → ", out_vln, "\n"))
}

# ============================================================
# STEP 5: DEG – CAR+ vs CAR- per paziente
# ============================================================
section("STEP 5 | DEG: CAR+ vs CAR- in I (per paziente)")

cat("\n[NOTA METODOLOGICA]\n")
cat("  I DEG sono calcolati con FindMarkers (Wilcoxon) PER PAZIENTE.\n")
cat("  Ogni cellula è trattata come un'osservazione indipendente.\n")
cat("  Questo è un approccio comune ma tecnicamente problematico\n")
cat("  perché viola l'assunzione di indipendenza (pseudo-replication).\n")
cat("  Per risultati robusti: cerca geni consistentemente significativi\n")
cat("  in ≥2 pazienti. Un gene trovato in un solo paziente ha valore limitato.\n\n")

deg_results <- list()
wb_deg <- createWorkbook()

for (nm in names(I_objects)) {
  cat(paste0("── DEG: ", nm, " ──\n"))
  obj <- I_objects[[nm]]

  n_pos <- sum(obj$CAR_status == "CAR+")
  n_neg <- sum(obj$CAR_status == "CAR-")
  cat(sprintf("  Cellule CAR+: %d | CAR-: %d\n", n_pos, n_neg))

  if (n_pos < 5 || n_neg < 5) {
    cat("  [SKIP] Troppo poche cellule in un gruppo (min 5)\n")
    next
  }

  Idents(obj) <- "CAR_status"

  # FindMarkers: CAR+ vs CAR-
  tryCatch({
    deg <- FindMarkers(
      obj,
      ident.1         = "CAR+",
      ident.2         = "CAR-",
      min.pct         = 0.10,   # Gene presente in ≥10% di almeno un gruppo
      logfc.threshold = 0.25,   # Log2FC ≥ 0.25 (soglia conservativa)
      test.use        = "wilcox",
      verbose         = FALSE
    )

    deg$gene   <- rownames(deg)
    deg$sample <- nm
    deg <- deg[order(deg$avg_log2FC, decreasing = TRUE), ]

    cat(sprintf("  Geni significativi (p_adj<0.05): %d\n",
                sum(deg$p_val_adj < 0.05, na.rm = TRUE)))
    cat(sprintf("  Up in CAR+ (logFC>0.5, p<0.05): %d\n",
                sum(deg$avg_log2FC > 0.5 & deg$p_val_adj < 0.05, na.rm = TRUE)))
    cat(sprintf("  Down in CAR+ (logFC<-0.5, p<0.05): %d\n",
                sum(deg$avg_log2FC < -0.5 & deg$p_val_adj < 0.05, na.rm = TRUE)))

    cat("  Top 10 UP in CAR+:\n")
    top_up <- head(deg[deg$avg_log2FC > 0, ], 10)
    print(top_up[, c("gene","avg_log2FC","p_val_adj","pct.1","pct.2")],
          row.names = FALSE)

    deg_results[[nm]] <- deg
    addWorksheet(wb_deg, substr(nm, 1, 31))
    writeData(wb_deg, substr(nm, 1, 31), deg)

  }, error = function(e) {
    cat(sprintf("  [ERRORE] %s: %s\n", nm, conditionMessage(e)))
  })
}

# Geni consistenti tra pazienti
if (length(deg_results) >= 2) {
  section("Geni consistenti tra pazienti")
  cat("\n(geni significativi in ≥2 pazienti, stessa direzione)\n\n")

  all_deg_df <- bind_rows(deg_results)
  genes_sig <- all_deg_df %>%
    filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.25) %>%
    group_by(gene) %>%
    summarise(
      n_patients_up   = sum(avg_log2FC > 0),
      n_patients_down = sum(avg_log2FC < 0),
      mean_logFC      = round(mean(avg_log2FC), 3),
      .groups = "drop"
    ) %>%
    filter(n_patients_up >= 2 | n_patients_down >= 2) %>%
    arrange(desc(abs(mean_logFC)))

  cat("Geni up in CAR+ in ≥2 pazienti:\n")
  print(genes_sig %>% filter(n_patients_up >= 2), n = 20)
  cat("\nGeni down in CAR+ in ≥2 pazienti:\n")
  print(genes_sig %>% filter(n_patients_down >= 2), n = 20)

  addWorksheet(wb_deg, "Consistenti_multi_paz")
  writeData(wb_deg, "Consistenti_multi_paz", genes_sig)
}

saveWorkbook(wb_deg, paste0(out_dir, "Q2_DEG_CAR_pos_vs_neg_in_I.xlsx"),
             overwrite = TRUE)
cat(paste0("\n  → Q2_DEG_CAR_pos_vs_neg_in_I.xlsx\n"))

# ============================================================
# STEP 6: HEATMAP DEI GENI FUNZIONALI CHIAVE
#         (media per tipo cellulare, CAR+ vs CAR-)
# ============================================================
section("STEP 6 | Heatmap geni funzionali chiave")

# Geni candidati per visualizzazione (subset curato)
KEY_GENES <- c(
  # Stemness/memory
  "TCF7", "CCR7", "IL7R", "SELL", "LEF1",
  # Effettore
  "GZMB", "PRF1", "NKG7", "GZMA", "CX3CR1",
  # Exhaustion
  "PDCD1", "LAG3", "HAVCR2", "TIGIT", "TOX", "ENTPD1",
  # Proliferazione
  "MKI67", "TOP2A",
  # Attivazione
  "CD69", "CD44", "CD38"
)

for (nm in names(I_objects)) {
  obj       <- I_objects[[nm]]
  genes_ok  <- filter_genes(KEY_GENES, obj)
  if (length(genes_ok) < 5) next

  Idents(obj) <- "CAR_status"
  avg <- AverageExpression(obj,
                           features = genes_ok,
                           group.by = "CAR_status",
                           assay = "RNA", slot = "data")$RNA

  # Normalizza per gene (z-score per riga)
  mat_z <- t(scale(t(as.matrix(avg))))
  mat_z[is.nan(mat_z)] <- 0

  df_heat <- as.data.frame(mat_z) %>%
    tibble::rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "group", values_to = "z_score") %>%
    mutate(
      gene  = factor(gene, levels = rev(genes_ok)),
      group = factor(group, levels = sort(unique(group)))
    )

  p_heat_gene <- ggplot(df_heat, aes(x = group, y = gene, fill = z_score)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_gradientn(
      colors = c("#2196F3", "white", "#E63946"),
      name   = "Z-score\nexpr",
      limits = c(-3, 3), oob = scales::squish
    ) +
    labs(
      title    = paste0(nm, " – Geni funzionali chiave: CAR+ vs CAR-"),
      subtitle = "Z-score dell'espressione media; gene × stato CAR",
      x = NULL, y = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title  = element_text(face = "bold", hjust = 0.5, size = 11),
      axis.text.x = element_text(size = 11, face = "bold"),
      axis.text.y = element_text(size = 9),
      panel.grid  = element_blank()
    )

  out_heat <- paste0(out_dir, "Q2_", nm, "_key_genes_heatmap.png")
  ggsave(out_heat, p_heat_gene, width = 6, height = 10,
         dpi = 300, bg = "white")
  cat(paste0("  → ", out_heat, "\n"))
}

cat(paste0(
  "\n", strrep("=", 65), "\n",
  "  Q2 COMPLETATA\n\n",
  "  INTERPRETAZIONE ATTESA:\n",
  "  - TCF7 alto + GZMB basso nelle CAR+ → fenotipo stemlike,\n",
  "    buon segno prognostico per persistenza in vivo.\n",
  "  - GZMB/PRF1 alti nelle CAR+ → già effettrici alla somministrazione;\n",
  "    rischio di esaurimento più rapido in vivo.\n",
  "  - PDCD1/LAG3/HAVCR2 già elevati in I → esaurimento pre-esistente\n",
  "    nel prodotto; può limitare l'efficacia.\n",
  "  - MKI67 elevato in CAR+ → buona frazione proliferante.\n",
  "  NOTA: Confronta i pattern tra pazienti per valutare\n",
  "        se le differenze nella qualità del prodotto correlano\n",
  "        con la successiva espansione/persistenza in vivo.\n",
  strrep("=", 65), "\n"
))
