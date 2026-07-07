# ============================================================
#  CAR-T DETECTION DA OGGETTO SEURAT
#  Metodi alternativi a IS_CAR_ALLIN_scREP (VDJ-based)
#
#  ADATTATO ALLA STRUTTURA REALE DELL'OGGETTO:
#  ─ Lista nominata di 8 oggetti Seurat (già splittata)
#  ─ Riduzioni disponibili: pca (50 dim) + umap (2 dim)
#  ─ Nessuna harmony
#  ─ Seurat v5 (Assay5): layers = counts.*, data, scale.data
#  ─ Metadato CAR: IS_CAR_ALLIN_scREP (YES / NO)
#  ─ Ca_bone_AB ha 0 CAR+ → viene skippato automaticamente
#
#  METODO A – Firma trascrittomica (DEG + module score)
#    DEG tra CAR+ e CAR- confermati → AddModuleScore su tutti
#    i T cells → soglia = 95° pct score dei CAR- (no data leakage)
#
#  METODO B – kNN neighborhood score
#    Per ogni T cell: frazione dei k=20 vicini più prossimi
#    in spazio PCA (30 dim) che sono IS_CAR_ALLIN_scREP = YES.
#    Non assume firma globale, sfrutta struttura locale del
#    manifold trascrittomica.
#
#  INTEGRAZIONE: A ∩ B = candidati ad alta confidenza
#
#  LIMITE INVALICABILE:
#    Entrambi i metodi trovano cellule che ASSOMIGLIANO alle
#    CAR-T confermate. Senza prova molecolare (costrutto o VDJ)
#    non si può escludere T cells endogeni con fenotipo simile.
#
#  Input:   all_samples_annotated_COMPLETE.rds
#  Output:  <out_dir>/
# ============================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(openxlsx)

# ── PARAMETRI — modifica solo questa sezione ─────────────────

rds_path <- "/Users/federicadannunzio/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/caruana-project/CART/Code/Seurat_analysis/2_annotation/all_samples_annotated_COMPLETE.rds"

out_dir <- "/Users/federicadannunzio/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/caruana-project/CART/Code/Seurat_analysis/3a_CAR_expr_detection/seurat_methods/"

# Colonna metadato gold standard (valori: "YES" / "NO")
META_CAR_COL <- "IS_CAR_ALLIN_scREP"
META_CAR_POS <- "YES"

# Colonna cell type
META_CELLTYPE_COL <- "cell_type"

# Tipi cellulari su cui eseguire l'analisi.
# Limitare ai T cells riduce drasticamente i falsi positivi.
T_CELL_TYPES <- c(
  "Cytotoxic CD8+ T cells",
  "Naive CD8+ T cells",
  "Naive CD4+ T cells",
  "Memory T cells",
  "Effector CD4+ T cells",
  "Proliferating T cells",
  "Proliferating CD4+ T cells",
  "Proliferating CD8+ T cells",
  "Th1 cells", "Th2 cells", "Th17 cells",
  "Tfh cells", "Tregs",
  "NKT cells", "MAIT cells",
  "gamma-delta T cells"
)

# ── Parametri Metodo A ────────────────────────────────────────
N_GENES_SIGNATURE  <- 50    # geni top DEG per la firma
FDR_THRESHOLD      <- 0.05
LFC_THRESHOLD      <- 0.4   # log2FC minimo
SCORE_PERCENTILE_A <- 95    # soglia: pct dei CAR- confermati

# ── Parametri Metodo B ────────────────────────────────────────
# Riduzione PCA: unica disponibile nell'oggetto (no harmony)
REDUCTION_NAME <- "pca"
REDUCTION_DIMS <- 30  # coerente con FindNeighbors del tuo workflow

# k vicini per il neighborhood score
KNN_K <- 20

# Soglia kNN score: frazione di vicini CAR+ sopra cui
# classificare la cellula come candidata CAR+
KNN_THRESHOLD_B <- 0.30

# ─────────────────────────────────────────────────────────────

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

section <- function(title)
  cat(paste0("\n", strrep("=", 65), "\n  ", title,
             "\n", strrep("=", 65), "\n"))

# ============================================================
# CARICAMENTO
# L'oggetto è già una lista nominata di 8 Seurat objects.
# Non è necessario nessun SplitObject.
# ============================================================

section("Caricamento oggetto Seurat")

cat("Caricamento:", rds_path, "\n")
all_samples <- readRDS(rds_path)

# Verifica che sia effettivamente una lista nominata
if (!is.list(all_samples) || inherits(all_samples, "Seurat")) {
  stop(paste0(
    "L'oggetto caricato non e' una lista di Seurat objects.\n",
    "Struttura trovata: ", class(all_samples), "\n",
    "Controlla rds_path."))
}

cat(sprintf("Campioni trovati: %d\n", length(all_samples)))
cat("Nomi campioni:\n")
for (nm in names(all_samples))
  cat(sprintf("  %s: %d cellule\n", nm, ncol(all_samples[[nm]])))

# ============================================================
# HELPER: kNN score
# Calcola per ogni cellula la frazione di k nearest neighbors
# che sono IS_CAR_ALLIN_scREP = YES.
# Usa le embeddings PCA (30 dim) già presenti nell'oggetto.
# ============================================================

compute_knn_score <- function(embeddings, is_car_pos, k) {

  if (!requireNamespace("FNN", quietly = TRUE)) {
    cat("  [INFO] Pacchetto FNN non trovato.\n")
    cat("         Installa con: install.packages('FNN')\n")
    cat("         Fallback: dist() base R (lento su dataset grandi).\n")
    D <- as.matrix(dist(embeddings))
    diag(D) <- Inf
    knn_idx <- t(apply(D, 1, function(x) order(x)[seq_len(k)]))
  } else {
    # k+1 perche' il primo vicino e' la cellula stessa (dist=0)
    knn_res <- FNN::get.knnx(embeddings, embeddings, k = k + 1)
    knn_idx <- knn_res$nn.index[, -1, drop = FALSE]
  }

  is_pos_num <- as.numeric(is_car_pos)
  knn_score  <- apply(knn_idx, 1,
                      function(idx) mean(is_pos_num[idx]))
  return(knn_score)
}

# ============================================================
# FUNZIONE PRINCIPALE PER CAMPIONE
# ============================================================

analyze_sample <- function(obj, sample_name, out_dir) {

  cat(paste0("\n", strrep("-", 60), "\n",
             "  Campione: ", sample_name, "\n",
             strrep("-", 60), "\n"))

  DefaultAssay(obj) <- "RNA"
  meta <- obj@meta.data

  # ── Verifica colonne necessarie ───────────────────────────
  missing_cols <- setdiff(c(META_CAR_COL, META_CELLTYPE_COL),
                          colnames(meta))
  if (length(missing_cols) > 0) {
    cat(sprintf("  [SKIP] Colonne mancanti nei metadati: %s\n",
                paste(missing_cols, collapse = ", ")))
    return(NULL)
  }

  # ── Gold standard VDJ ────────────────────────────────────
  is_car_screp <- meta[[META_CAR_COL]] == META_CAR_POS
  is_car_screp[is.na(is_car_screp)] <- FALSE
  n_car_screp  <- sum(is_car_screp)
  n_total      <- ncol(obj)

  cat(sprintf("  Cellule totali:             %d\n", n_total))
  cat(sprintf("  CAR+ (IS_CAR_ALLIN_scREP): %d (%.1f%%)\n",
              n_car_screp, n_car_screp / n_total * 100))

  if (n_car_screp == 0) {
    cat(paste0(
      "  [SKIP] Nessuna cellula CAR+ di riferimento.\n",
      "         Impossibile costruire firma o kNN score.\n",
      "         Per Ca_bone_AB: verificare se il file VDJ\n",
      "         e' stato processato correttamente.\n"))
    return(NULL)
  }

  # ── Subset T cells ────────────────────────────────────────
  cell_types_in_sample <- unique(as.character(
    meta[[META_CELLTYPE_COL]]))
  t_types_found <- intersect(T_CELL_TYPES, cell_types_in_sample)

  if (length(t_types_found) == 0) {
    cat("  [WARN] Nessun tipo T cell trovato. Uso tutte le cellule.\n")
    cat("  Tipi presenti:", paste(cell_types_in_sample,
                                   collapse = ", "), "\n")
    is_tcell <- rep(TRUE, n_total)
    t_cells  <- obj
  } else {
    is_tcell <- meta[[META_CELLTYPE_COL]] %in% t_types_found
    t_cells  <- subset(obj, cells = colnames(obj)[is_tcell])
    cat(sprintf("  T cells: %d / %d (%.1f%%)\n",
                sum(is_tcell), n_total,
                sum(is_tcell) / n_total * 100))
  }

  meta_t         <- t_cells@meta.data
  is_car_screp_t <- meta_t[[META_CAR_COL]] == META_CAR_POS
  is_car_screp_t[is.na(is_car_screp_t)] <- FALSE
  n_car_t <- sum(is_car_screp_t)
  n_tcell <- ncol(t_cells)

  cat(sprintf("  CAR+ nei T cells: %d / %d totali\n",
              n_car_t, n_car_screp))

  if (sum(is_car_screp & !is_tcell) > 0)
    cat(sprintf("  [NOTA] %d CAR+ fuori dai T cell types definiti.\n",
                sum(is_car_screp & !is_tcell)))

  if (n_car_t < 5)
    cat(sprintf(
      "  [WARN] Solo %d CAR+ nei T cells: DEG poco affidabili.\n",
      n_car_t))

  res <- list(
    sample      = sample_name,
    n_total     = n_total,
    n_car_screp = n_car_screp,
    n_tcell     = n_tcell,
    n_car_t     = n_car_t
  )

  # ===========================================================
  # SEURAT v5: JoinLayers prima di FindMarkers
  # Necessario quando i counts sono in layer separati (Assay5).
  # JoinLayers e' idempotente se i layer sono gia' uniti.
  # ===========================================================

  tryCatch({
    t_cells <- JoinLayers(t_cells)
    cat("  JoinLayers: OK\n")
  }, error = function(e) {
    # Seurat < 5: JoinLayers non esiste, nessun problema
  })

  method_A_ok <- FALSE
  method_B_ok <- FALSE

  # ===========================================================
  # METODO A: FIRMA TRASCRITTOMICA
  # ===========================================================

  cat("\n  ── METODO A: Firma trascrittomica ──\n")

  tryCatch({

    Idents(t_cells) <- ifelse(is_car_screp_t, "CAR_pos", "CAR_neg")

    deg <- FindMarkers(
      t_cells,
      ident.1         = "CAR_pos",
      ident.2         = "CAR_neg",
      min.pct         = 0.05,
      logfc.threshold = LFC_THRESHOLD,
      test.use        = "wilcox",
      verbose         = FALSE
    )
    deg$gene <- rownames(deg)

    sig_up <- deg %>%
      dplyr::filter(p_val_adj < FDR_THRESHOLD,
                    avg_log2FC > LFC_THRESHOLD) %>%
      dplyr::arrange(desc(avg_log2FC)) %>%
      head(N_GENES_SIGNATURE)

    cat(sprintf(
      "  DEG upregolati CAR+ (FDR<%.2f, lFC>%.1f): %d\n",
      FDR_THRESHOLD, LFC_THRESHOLD, nrow(sig_up)))

    if (nrow(sig_up) == 0) {
      cat(paste0(
        "  [WARN] Nessun DEG trovato.\n",
        "         Cause: poche cellule CAR+ o profilo CAR-T\n",
        "         indistinguibile da T endogeni in questo campione.\n",
        "  → Metodo A non applicabile. Uso solo Metodo B.\n"))

    } else {

      # AddModuleScore aggiunge automaticamente "CAR_sig_score1"
      t_cells <- AddModuleScore(
        t_cells,
        features = list(sig_up$gene),
        name     = "CAR_sig_score",
        ctrl     = 100,
        seed     = 42
      )
      score_col <- "CAR_sig_score1"

      # Soglia sui CAR- (no data leakage)
      scores_carneg <- t_cells@meta.data[[score_col]][!is_car_screp_t]
      thr_A <- quantile(scores_carneg,
                        SCORE_PERCENTILE_A / 100,
                        na.rm = TRUE)
      cat(sprintf("  Soglia A (%d° pct su CAR-): %.4f\n",
                  SCORE_PERCENTILE_A, thr_A))

      is_car_A  <- t_cells@meta.data[[score_col]] > thr_A
      n_car_A   <- sum(is_car_A)
      overlap_A <- sum(is_car_A & is_car_screp_t)
      new_A     <- sum(is_car_A & !is_car_screp_t)
      sens_A    <- round(overlap_A / n_car_t * 100, 1)

      cat(sprintf("  CAR+ da A:       %d (%.1f%% dei T cells)\n",
                  n_car_A, n_car_A / n_tcell * 100))
      cat(sprintf("  Overlap scREP:   %d / %d (sens. %.1f%%)\n",
                  overlap_A, n_car_t, sens_A))
      cat(sprintf("  Nuovi candidati: %d\n", new_A))

      t_cells$CAR_method_A <- ifelse(is_car_A, "CAR+", "CAR-")

      res$n_car_A   <- n_car_A
      res$new_A     <- new_A
      res$sens_A    <- sens_A
      res$thr_A     <- thr_A
      res$deg_genes <- sig_up$gene
      method_A_ok   <- TRUE

      write.csv(deg,
                paste0(out_dir, sample_name,
                       "_DEG_CARpos_vs_CARneg.csv"),
                row.names = TRUE)
      cat(sprintf("  → %s_DEG_CARpos_vs_CARneg.csv\n",
                  sample_name))

      # Plot score
      df_score <- data.frame(
        score = t_cells@meta.data[[score_col]],
        group = ifelse(is_car_screp_t,
                       "CAR+ (scREP)", "CAR- (scREP)"),
        stringsAsFactors = FALSE)

      p_score <- ggplot(df_score,
                        aes(x = score, fill = group)) +
        geom_density(alpha = 0.6, color = NA) +
        geom_vline(xintercept = thr_A,
                   linetype = "dashed", color = "#B00020",
                   linewidth = 0.9) +
        annotate("text",
                 x = thr_A, y = Inf, vjust = 1.5,
                 hjust = -0.1,
                 label = sprintf("%d° pct\n(soglia A)",
                                 SCORE_PERCENTILE_A),
                 color = "#B00020", size = 3.5) +
        scale_fill_manual(
          values = c("CAR+ (scREP)" = "#264653",
                     "CAR- (scREP)" = "#E9C46A")) +
        labs(title    = paste0(sample_name,
                               " - Module score CAR"),
             subtitle = paste0("Score su ",
                               nrow(sig_up), " DEG"),
             x = "Module score", y = "Densita'",
             fill = NULL) +
        theme_classic(base_size = 11) +
        theme(
          plot.title    = element_text(face = "bold",
                                        hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5,
                                        size = 9,
                                        color = "gray40"),
          legend.position = "top")

      ggsave(paste0(out_dir, sample_name,
                    "_score_density_A.png"),
             plot = p_score, width = 8, height = 5,
             dpi = 300, bg = "white")
      cat(sprintf("  -> %s_score_density_A.png\n",
                  sample_name))
    }

  }, error = function(e) {
    cat(sprintf("  [ERRORE Metodo A] %s\n", e$message))
  })

  # ===========================================================
  # METODO B: kNN NEIGHBORHOOD SCORE
  # Usa PCA (30 dim) gia' calcolata nell'oggetto.
  # ===========================================================

  cat("\n  ── METODO B: kNN neighborhood score ──\n")

  tryCatch({

    if (!REDUCTION_NAME %in% names(obj@reductions))
      stop(sprintf(
        "Riduzione '%s' non trovata. Disponibili: %s",
        REDUCTION_NAME,
        paste(names(obj@reductions), collapse = ", ")))

    avail_dims  <- ncol(
      obj@reductions[[REDUCTION_NAME]]@cell.embeddings)
    actual_dims <- min(REDUCTION_DIMS, avail_dims)
    cat(sprintf("  Riduzione: %s, %d dim\n",
                REDUCTION_NAME, actual_dims))

    # Embeddings PCA per i T cells
    emb_all <- obj@reductions[[REDUCTION_NAME]]@cell.embeddings
    emb_t   <- emb_all[colnames(t_cells),
                        seq_len(actual_dims),
                        drop = FALSE]

    cat(sprintf(
      "  Calcolo kNN score (k=%d) su %d T cells...\n",
      KNN_K, nrow(emb_t)))

    knn_score <- compute_knn_score(
      embeddings = emb_t,
      is_car_pos = is_car_screp_t,
      k          = KNN_K
    )
    cat("  Completato.\n")
    cat(sprintf("  Score medio CAR+ (scREP): %.3f\n",
                mean(knn_score[is_car_screp_t])))
    cat(sprintf("  Score medio CAR- (scREP): %.3f\n",
                mean(knn_score[!is_car_screp_t])))

    is_car_B  <- knn_score >= KNN_THRESHOLD_B
    n_car_B   <- sum(is_car_B)
    overlap_B <- sum(is_car_B & is_car_screp_t)
    new_B     <- sum(is_car_B & !is_car_screp_t)
    sens_B    <- round(overlap_B / n_car_t * 100, 1)

    cat(sprintf(
      "  CAR+ da B (soglia %.2f): %d (%.1f%% dei T cells)\n",
      KNN_THRESHOLD_B, n_car_B, n_car_B / n_tcell * 100))
    cat(sprintf("  Overlap scREP:           %d / %d (sens. %.1f%%)\n",
                overlap_B, n_car_t, sens_B))
    cat(sprintf("  Nuovi candidati:         %d\n", new_B))

    t_cells$CAR_knn_score <- knn_score
    t_cells$CAR_method_B  <- ifelse(is_car_B, "CAR+", "CAR-")

    res$knn_score <- knn_score
    res$n_car_B   <- n_car_B
    res$new_B     <- new_B
    res$sens_B    <- sens_B
    method_B_ok   <- TRUE

    # Plot kNN score
    df_knn <- data.frame(
      knn_score = knn_score,
      group     = ifelse(is_car_screp_t,
                         "CAR+ (scREP)", "CAR- (scREP)"),
      stringsAsFactors = FALSE)

    p_knn <- ggplot(df_knn,
                    aes(x = knn_score, fill = group)) +
      geom_density(alpha = 0.6, color = NA) +
      geom_vline(xintercept = KNN_THRESHOLD_B,
                 linetype = "dashed", color = "#B00020",
                 linewidth = 0.9) +
      annotate("text",
               x = KNN_THRESHOLD_B, y = Inf,
               vjust = 1.5, hjust = -0.1,
               label = sprintf("soglia\n%.2f",
                               KNN_THRESHOLD_B),
               color = "#B00020", size = 3.5) +
      scale_fill_manual(
        values = c("CAR+ (scREP)" = "#264653",
                   "CAR- (scREP)" = "#E9C46A")) +
      labs(title    = paste0(sample_name,
                             " - kNN CAR score"),
           subtitle = paste0("Frazione k=", KNN_K,
                             " vicini IS_CAR+"),
           x = "kNN score", y = "Densita'",
           fill = NULL) +
      theme_classic(base_size = 11) +
      theme(
        plot.title    = element_text(face = "bold",
                                      hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5,
                                      size = 9,
                                      color = "gray40"),
        legend.position = "top")

    ggsave(paste0(out_dir, sample_name,
                  "_knn_score_B.png"),
           plot = p_knn, width = 8, height = 5,
           dpi = 300, bg = "white")
    cat(sprintf("  -> %s_knn_score_B.png\n", sample_name))

  }, error = function(e) {
    cat(sprintf("  [ERRORE Metodo B] %s\n", e$message))
  })

  # ===========================================================
  # INTEGRAZIONE A e B
  # ===========================================================

  cat("\n  ── INTEGRAZIONE A e B ──\n")

  has_A <- method_A_ok &&
           "CAR_method_A" %in% colnames(t_cells@meta.data)
  has_B <- method_B_ok &&
           "CAR_method_B" %in% colnames(t_cells@meta.data)

  if (!has_A && !has_B) {
    cat("  Nessun metodo riuscito per questo campione.\n")
    return(res)
  }

  if (has_A && has_B) {
    is_A <- t_cells$CAR_method_A == "CAR+"
    is_B <- t_cells$CAR_method_B == "CAR+"
    cat(sprintf(
      "  Nuovi: solo A=%d | solo B=%d | A intersez. B=%d\n",
      sum(is_A & !is_B & !is_car_screp_t),
      sum(!is_A & is_B & !is_car_screp_t),
      sum(is_A & is_B & !is_car_screp_t)))
    t_cells$CAR_integrated <- dplyr::case_when(
      is_car_screp_t ~ "scREP_confirmed",
      is_A & is_B    ~ "new_A_and_B",
      is_A & !is_B   ~ "new_A_only",
      !is_A & is_B   ~ "new_B_only",
      TRUE           ~ "CAR_negative")
    res$n_new_integrated <-
      sum(is_A & is_B & !is_car_screp_t)
    res$n_only_A <- sum(is_A & !is_B & !is_car_screp_t)
    res$n_only_B <- sum(!is_A & is_B & !is_car_screp_t)

  } else if (has_B) {
    cat("  Solo Metodo B disponibile.\n")
    is_B <- t_cells$CAR_method_B == "CAR+"
    t_cells$CAR_integrated <- dplyr::case_when(
      is_car_screp_t ~ "scREP_confirmed",
      is_B           ~ "new_B_only",
      TRUE           ~ "CAR_negative")
    res$n_new_integrated <- sum(is_B & !is_car_screp_t)

  } else {
    cat("  Solo Metodo A disponibile.\n")
    is_A <- t_cells$CAR_method_A == "CAR+"
    t_cells$CAR_integrated <- dplyr::case_when(
      is_car_screp_t ~ "scREP_confirmed",
      is_A           ~ "new_A_only",
      TRUE           ~ "CAR_negative")
    res$n_new_integrated <- sum(is_A & !is_car_screp_t)
  }

  # ===========================================================
  # UMAP
  # IMPORTANTE: nell'oggetto i nomi delle colonne UMAP sono
  # "umap_1" e "umap_2" (Seurat v5, non UMAP_1/UMAP_2)
  # ===========================================================

  if ("umap" %in% names(obj@reductions)) {

    umap_emb <- obj@reductions[["umap"]]@cell.embeddings
    # Rinomina per sicurezza (nomi originali: umap_1, umap_2)
    colnames(umap_emb) <- c("UMAP_1", "UMAP_2")

    umap_t <- as.data.frame(
      umap_emb[colnames(t_cells), , drop = FALSE])

    df_umap <- data.frame(
      UMAP_1 = umap_t$UMAP_1,
      UMAP_2 = umap_t$UMAP_2,
      class  = t_cells$CAR_integrated,
      stringsAsFactors = FALSE)

    # Ordina: CAR- sotto, CAR+ sopra
    priority <- c("CAR_negative", "new_A_only",
                  "new_B_only", "new_A_and_B",
                  "scREP_confirmed")
    df_umap$class <- factor(df_umap$class, levels = priority)
    df_umap <- df_umap[order(df_umap$class), ]

    int_colors <- c(
      "scREP_confirmed" = "#264653",
      "new_A_and_B"     = "#E63946",
      "new_A_only"      = "#F4A261",
      "new_B_only"      = "#2A9D8F",
      "CAR_negative"    = "#CCCCCC"
    )
    int_labels <- c(
      "scREP_confirmed" = "CAR+ scREP (gold standard)",
      "new_A_and_B"     = "Nuovi CAR+ A+B (alta conf.)",
      "new_A_only"      = "Nuovi CAR+ solo firma",
      "new_B_only"      = "Nuovi CAR+ solo kNN",
      "CAR_negative"    = "CAR-"
    )

    n_new_int <- if (!is.null(res$n_new_integrated))
      res$n_new_integrated else 0

    p_umap <- ggplot(df_umap,
                     aes(x = UMAP_1, y = UMAP_2,
                         color = class)) +
      geom_point(
        data  = df_umap[df_umap$class == "CAR_negative", ],
        size  = 0.3, alpha = 0.3) +
      geom_point(
        data  = df_umap[df_umap$class != "CAR_negative", ],
        size  = 1.2, alpha = 0.95) +
      scale_color_manual(values = int_colors,
                         labels = int_labels,
                         drop   = FALSE) +
      labs(
        title    = paste0(sample_name,
                          " - Classificazione CAR"),
        subtitle = paste0(
          "T cells: ", n_tcell,
          "  |  scREP: ", n_car_t,
          "  |  Nuovi A+B: ", n_new_int),
        color = NULL) +
      guides(color = guide_legend(
        override.aes = list(size = 3, alpha = 1))) +
      theme_void(base_size = 11) +
      theme(
        plot.title      = element_text(face = "bold",
                                        hjust = 0.5),
        plot.subtitle   = element_text(hjust = 0.5,
                                        size = 9,
                                        color = "gray40"),
        legend.position = "right",
        legend.text     = element_text(size = 8))

    ggsave(paste0(out_dir, sample_name,
                  "_UMAP_CAR_integrated.png"),
           plot = p_umap, width = 10, height = 7,
           dpi = 300, bg = "white")
    cat(sprintf("  -> %s_UMAP_CAR_integrated.png\n",
                sample_name))
  }

  # ===========================================================
  # EXPORT METADATI
  # ===========================================================

  cols_keep <- intersect(
    c(META_CELLTYPE_COL, META_CAR_COL,
      "CAR_method_A", "CAR_sig_score1",
      "CAR_knn_score", "CAR_method_B",
      "CAR_integrated"),
    colnames(t_cells@meta.data))

  meta_export <- t_cells@meta.data[, cols_keep, drop = FALSE]
  meta_export$barcode     <- rownames(meta_export)
  meta_export$sample_name <- sample_name

  write.csv(meta_export,
            paste0(out_dir, sample_name,
                   "_CAR_classification.csv"),
            row.names = FALSE)
  cat(sprintf("  -> %s_CAR_classification.csv\n", sample_name))

  res$meta_export <- meta_export
  return(res)
}

# ============================================================
# LOOP PRINCIPALE
# ============================================================

section("Analisi per campione")

all_results <- list()
for (nm in names(all_samples)) {
  all_results[[nm]] <- analyze_sample(
    obj         = all_samples[[nm]],
    sample_name = nm,
    out_dir     = out_dir)
}
all_results <- Filter(Negate(is.null), all_results)

cat(sprintf(
  "\nCampioni processati con successo: %d / %d\n",
  length(all_results), length(all_samples)))

# ============================================================
# TABELLA RIEPILOGATIVA
# ============================================================

section("Riepilogo concordanza")

summary_rows <- lapply(names(all_results), function(nm) {
  r <- all_results[[nm]]
  data.frame(
    campione        = nm,
    n_totale        = r$n_total,
    n_tcell         = r$n_tcell,
    n_CAR_scREP     = r$n_car_screp,
    pct_CAR_scREP   = round(r$n_car_screp / r$n_total * 100, 2),
    n_CAR_A         = if (!is.null(r$n_car_A)) r$n_car_A else NA,
    sens_A_pct      = if (!is.null(r$sens_A)) r$sens_A else NA,
    new_CAR_A       = if (!is.null(r$new_A)) r$new_A else NA,
    n_CAR_B         = if (!is.null(r$n_car_B)) r$n_car_B else NA,
    sens_B_pct      = if (!is.null(r$sens_B)) r$sens_B else NA,
    new_CAR_B       = if (!is.null(r$new_B)) r$new_B else NA,
    new_CAR_A_and_B = if (!is.null(r$n_new_integrated))
                        r$n_new_integrated else NA,
    stringsAsFactors = FALSE)
})
summary_df <- dplyr::bind_rows(summary_rows)

cat("\n")
print(summary_df, row.names = FALSE)

# ============================================================
# EXCEL RIEPILOGATIVO
# ============================================================

section("Excel riepilogativo")

wb <- createWorkbook()
addWorksheet(wb, "Riepilogo")
writeData(wb, "Riepilogo", summary_df)

for (nm in names(all_results)) {
  r <- all_results[[nm]]
  if (!is.null(r$meta_export)) {
    sheet_nm <- substr(nm, 1, 31)
    addWorksheet(wb, sheet_nm)
    writeData(wb, sheet_nm, r$meta_export)
  }
}

xlsx_path <- paste0(out_dir,
                    "CAR_seurat_methods_concordance.xlsx")
saveWorkbook(wb, xlsx_path, overwrite = TRUE)
cat(paste0("  -> ", xlsx_path, "\n"))

# ============================================================
# RIEPILOGO FINALE
# ============================================================

section("Completato")

cat(paste0(
  "\n  Output per campione:\n",
  "    _DEG_CARpos_vs_CARneg.csv   DEG usati per la firma\n",
  "    _score_density_A.png        distribuzione module score\n",
  "    _knn_score_B.png            distribuzione kNN score\n",
  "    _UMAP_CAR_integrated.png    UMAP con classificazione\n",
  "    _CAR_classification.csv     metadati per cellula\n",
  "\n  Output globale:\n",
  "    CAR_seurat_methods_concordance.xlsx\n",
  "\n  INTERPRETAZIONE:\n",
  "  - sens_A/B > 70%  -> metodo riconosce i CAR+ gia' noti\n",
  "  - new_A_and_B     -> alta confidenza (due metodi indip.)\n",
  "  - new_A_only      -> controlla DEG: geni di attivazione\n",
  "                       generica = probabli falsi positivi\n",
  "  - new_B_only      -> controlla UMAP: se isolati dai\n",
  "                       cluster CAR+ = probabile noise\n",
  "  - Ca_bone_AB      -> SKIPPATO (0 CAR+ di riferimento)\n",
  "\n  LIMITE: i candidati nuovi sono cellule simili alle\n",
  "  CAR-T confermate. Senza prova molecolare non si puo'\n",
  "  escludere T cells endogeni con fenotipo analogo.\n",
  strrep("=", 65), "\n"))

