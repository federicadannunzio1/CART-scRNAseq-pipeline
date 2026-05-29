# ============================================================
#  Q1b: Composizione per STATI FUNZIONALI (approccio agnostico)
#
#  Motivazione:
#    L'annotazione CD4/CD8 da scRNA-seq non è affidabile in
#    questi dati (FACS mostra più CD4 CART ma scRNA-seq mostra
#    più CD8 → probabile dropout del trascritto CD4, artefatto
#    tecnico noto in 10x Chromium).
#
#  Approccio alternativo:
#    Raggruppa le cellule T in STATI FUNZIONALI trasversali
#    al lineage CD4/CD8:
#      Naive-like   → Naive CD4+ + Naive CD8+
#      Memory-like  → Memory T + Th1/Th2/Th17/Tfh
#      Effector     → Effector CD4+ + Cytotoxic CD8+
#      Regulatory   → Tregs
#      Proliferating→ Proliferating CD4+ + Proliferating CD8+
#
#    1. Proporzioni degli stati funzionali: I vs AB, per paziente
#       stratificate per stato CAR (CAR+ / CAR-)
#    2. Module scores funzionali: CAR+ vs CAR- in ciascuno
#       stato (I vs AB)
#    3. Sub-clustering T cells: risoluzione fine sulle sole
#       cellule T per identificare popolazioni funzionali
#       in base all'espressione genica (indipendente da CD4/CD8)
#    4. Trajectory / ordinamento cellule CAR+ per pseudotime
#       (opzionale, richiede monocle3)
#
#  Prerequisiti:
#    all_samples_annotated_COMPLETE_IS_CAR_REVISED.rds
#
#  Output in out_dir/:
#    Q1b_<paziente>_functional_states_I_vs_AB.png
#    Q1b_ALL_functional_state_heatmap.png
#    Q1b_<paziente>_module_scores_CARpos_vs_CARneg.png
#    Q1b_ALL_subcluster_UMAP.png
#    Q1b_functional_summary.xlsx
# ============================================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(scales)
library(openxlsx)

# ── UNICO PUNTO DA MODIFICARE ────────────────────────────────
rds_path <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/CART/Code/Seurat_analysis/2_annotation/all_samples_annotated_COMPLETE_IS_CAR_REVISED.rds"
out_dir  <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/CART/Code/Seurat_analysis/10_CART_functional_analysis/Q1b_functional_states/"
# ─────────────────────────────────────────────────────────────

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

section <- function(title)
  cat(paste0("\n", strrep("=", 65), "\n  ", title, "\n",
             strrep("=", 65), "\n"))

# ============================================================
# MAPPATURA STATI FUNZIONALI (agnostica per lineage)
# ============================================================
# Ogni tipo annotato viene assegnato a uno stato funzionale.
# L'assegnazione è basata sul COMPORTAMENTO FUNZIONALE atteso
# (quiescenza, homing, effettore, soppressione, ciclo cellulare)
# NON sulla molecola di superficie CD4 o CD8.

FUNCTIONAL_STATE_MAP <- list(
  "Naive-like" = c(
    "Naive CD4+ T cells",
    "Naive CD8+ T cells"
  ),
  "Memory-like" = c(
    "Memory T cells",  # Tcm/Tem con bassa espressione citotossica
    "Th1 cells",       # Effettori helper CD4 (inclusi per stato memoire)
    "Th2 cells",
    "Th17 cells",
    "Tfh cells"
  ),
  "Effector" = c(
    "Effector CD4+ T cells",
    "Cytotoxic CD8+ T cells"
  ),
  "Regulatory" = c(
    "Tregs"
  ),
  "Proliferating" = c(
    "Proliferating CD4+ T cells",
    "Proliferating CD8+ T cells"
  )
)

# Ordine biologico per plot (da più quiescente a più differenziato)
FUNCTIONAL_ORDER <- c(
  "Naive-like", "Memory-like", "Effector", "Regulatory", "Proliferating"
)

# Palette colori per gli stati funzionali
STATE_PALETTE <- c(
  "Naive-like"   = "#4DBBD5",
  "Memory-like"  = "#00A087",
  "Effector"     = "#E64B35",
  "Regulatory"   = "#F39B7F",
  "Proliferating"= "#7E6148"
)

# ============================================================
# FIRME GENICHE FUNZIONALI (identiche a Q2 per comparabilità)
# ============================================================
SIGNATURES <- list(

  Effector = c(
    "GZMB", "PRF1", "NKG7", "GNLY",
    "GZMA", "GZMK", "FGFBP2", "CX3CR1"
  ),

  Memory_Stemness = c(
    "TCF7", "CCR7", "SELL", "IL7R",
    "LEF1", "KLF2", "BCL2", "FOXO1"
  ),

  Exhaustion = c(
    "PDCD1", "LAG3", "HAVCR2", "TIGIT",
    "TOX", "TOX2", "ENTPD1", "CTLA4", "BATF"
  ),

  Activation = c(
    "CD69", "CD44", "TNFRSF9", "IL2RA", "ICOS", "CD38"
  ),

  Proliferation = c(
    "MKI67", "TOP2A", "PCNA", "CCNB1", "STMN1", "UBE2C"
  ),

  Tpex_StemLike = c(
    "TCF7", "CXCR5", "TOX", "BCL6", "SLAMF6", "ID3"
  ),

  Tex_Terminal = c(
    "HAVCR2", "TIGIT", "LAG3", "CD160", "ENTPD1", "PRDM1", "ZEB2"
  )
)

# ============================================================
# MAPPA CAMPIONI PER PAZIENTE
# ============================================================
PATIENT_MAP <- list(
  Bo = list(
    I  = c("Bo_bone_I"),
    AB = c("Bo_blood_AB", "Bo_bone_AB")
  ),
  Ca = list(
    I  = c("Ca_bone_I"),
    AB = c("Ca_blood_AB", "Ca_bone_AB")
  ),
  Me = list(
    I  = c("Me_bone_I"),
    AB = c("Me_bone_AB")
  )
)

# ============================================================
# HELPER FUNCTIONS
# ============================================================

# Recupera stato CAR dalla migliore colonna disponibile
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
  warning(sprintf("  %s: nessuna colonna CAR trovata. Tutte le cellule = CAR-", sample_name))
  return(rep("CAR-", ncol(obj)))
}

# Mappa il tipo cellulare annotato → stato funzionale
map_to_functional_state <- function(cell_types) {
  state <- rep("Other", length(cell_types))
  for (fs_name in names(FUNCTIONAL_STATE_MAP)) {
    hits <- cell_types %in% FUNCTIONAL_STATE_MAP[[fs_name]]
    state[hits] <- fs_name
  }
  state
}

# Calcola proporzioni per stato funzionale × CAR status
# all'interno di ciascun campione
extract_functional_props <- function(obj, sample_name, timepoint) {
  meta      <- obj@meta.data
  car_vec   <- get_car_status(obj, sample_name)
  cell_type <- as.character(meta$cell_type)

  fs_vec <- map_to_functional_state(cell_type)

  df <- data.frame(
    cell_type       = cell_type,
    functional_state= fs_vec,
    car_status      = car_vec,
    sample          = sample_name,
    timepoint       = timepoint,
    stringsAsFactors = FALSE
  )

  # Includi solo cellule T annotate (non "Other")
  df_t <- df[df$functional_state != "Other", ]

  if (nrow(df_t) == 0) {
    warning(sprintf("  %s: nessuna cellula T trovata con gli stati definiti.", sample_name))
    return(NULL)
  }

  # Proporzioni DENTRO ogni CAR status (CAR+ separato da CAR-)
  props <- df_t %>%
    group_by(car_status, functional_state) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(car_status) %>%
    mutate(
      total = sum(n),
      prop  = n / total
    ) %>%
    ungroup() %>%
    mutate(
      sample    = sample_name,
      timepoint = timepoint,
      functional_state = factor(functional_state, levels = FUNCTIONAL_ORDER)
    )

  props
}

# ============================================================
# CARICAMENTO DATI
# ============================================================
section("Caricamento dati")
all_samples <- readRDS(rds_path)
cat("Campioni disponibili:", paste(names(all_samples), collapse = ", "), "\n")

# ============================================================
# AGGIUNTA MODULE SCORES
# ============================================================
section("Calcolo module scores")

for (sname in names(all_samples)) {
  cat(sprintf("  Module scores: %s\n", sname))
  for (sig_name in names(SIGNATURES)) {
    genes_ok <- intersect(SIGNATURES[[sig_name]],
                          rownames(all_samples[[sname]]))
    if (length(genes_ok) < 2) {
      cat(sprintf("    [WARN] %s: solo %d gene/i per '%s', skip\n",
                  sname, length(genes_ok), sig_name))
      next
    }
    score_col <- paste0("score_", sig_name)
    all_samples[[sname]] <- AddModuleScore(
      all_samples[[sname]],
      features = list(genes_ok),
      name     = score_col,
      seed     = 42
    )
    # Rinomina: AddModuleScore aggiunge "1" in fondo
    old_col <- paste0(score_col, "1")
    all_samples[[sname]]@meta.data[[score_col]] <-
      all_samples[[sname]]@meta.data[[old_col]]
    all_samples[[sname]]@meta.data[[old_col]] <- NULL
  }
}

# ============================================================
# SEZIONE 1: PROPORZIONI STATI FUNZIONALI (I vs AB)
# ============================================================
section("Sezione 1: Proporzioni stati funzionali I vs AB")

all_props <- list()

for (patient in names(PATIENT_MAP)) {
  cat(sprintf("\nPaziente: %s\n", patient))
  pat_props <- list()

  for (tp in c("I", "AB")) {
    samples_tp <- PATIENT_MAP[[patient]][[tp]]
    for (sname in samples_tp) {
      if (!sname %in% names(all_samples)) {
        cat(sprintf("  [SKIP] %s non trovato\n", sname))
        next
      }
      p <- extract_functional_props(all_samples[[sname]], sname, tp)
      if (!is.null(p)) pat_props[[sname]] <- p
    }
  }

  if (length(pat_props) == 0) next

  pat_df <- bind_rows(pat_props)
  all_props[[patient]] <- pat_df

  # ── Barplot per paziente ─────────────────────────────────
  # Aggiunge etichetta campione per i campioni AB (blood vs bone)
  pat_df_plot <- pat_df %>%
    mutate(
      # Crea etichetta leggibile: timepoint + nome campione ridotto
      tp_label = case_when(
        timepoint == "I"  ~ "I",
        TRUE              ~ gsub(paste0(patient, "_"), "", sample)
      ),
      tp_label = factor(tp_label,
                        levels = unique(tp_label[order(timepoint, sample)]))
    )

  # Facet per CAR status
  p_bar <- ggplot(pat_df_plot,
                  aes(x = tp_label, y = prop,
                      fill = functional_state)) +
    geom_col(width = 0.75, color = "white", linewidth = 0.3) +
    facet_wrap(~ car_status, ncol = 2) +
    scale_fill_manual(values = STATE_PALETTE, drop = FALSE,
                      name = "Stato funzionale") +
    scale_y_continuous(labels = percent_format(), expand = c(0, 0)) +
    labs(
      title    = sprintf("%s — Composizione stati funzionali T cells", patient),
      subtitle = "Proporzioni calcolate DENTRO ogni stato CAR (CD4/CD8-agnostico)",
      x = NULL, y = "Proporzione"
    ) +
    theme_classic(base_size = 12) +
    theme(
      strip.background = element_rect(fill = "#F0F0F0"),
      strip.text       = element_text(face = "bold"),
      legend.position  = "right",
      axis.text.x      = element_text(angle = 30, hjust = 1)
    )

  ggsave(file.path(out_dir,
                   sprintf("Q1b_%s_functional_states_I_vs_AB.png", patient)),
         p_bar, width = 8, height = 5, dpi = 300)
  cat(sprintf("  Salvato: Q1b_%s_functional_states_I_vs_AB.png\n", patient))
}

# ── Heatmap inter-paziente ────────────────────────────────────
section("Heatmap inter-paziente")

heat_df <- bind_rows(all_props, .id = "patient")

# Solo CAR+ cellule in I
heat_car_pos_I <- heat_df %>%
  filter(car_status == "CAR+", timepoint == "I") %>%
  group_by(patient, functional_state) %>%
  summarise(prop = mean(prop), .groups = "drop") %>%
  mutate(functional_state = factor(functional_state, levels = FUNCTIONAL_ORDER))

if (nrow(heat_car_pos_I) > 0) {
  p_heat_I <- ggplot(heat_car_pos_I,
                     aes(x = patient, y = functional_state, fill = prop)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.1f%%", 100 * prop)),
              size = 3.5, color = "black") +
    scale_fill_gradient(low = "white", high = "#E64B35",
                        labels = percent_format(), name = "Proporzione") +
    labs(title = "CAR+ cells in I — Proporzione stati funzionali",
         subtitle = "Heatmap inter-paziente (agnostico CD4/CD8)",
         x = "Paziente", y = NULL) +
    theme_classic(base_size = 12)

  ggsave(file.path(out_dir, "Q1b_ALL_heatmap_CARpos_I.png"),
         p_heat_I, width = 6, height = 4, dpi = 300)
}

# Solo CAR+ cellule in AB (media tra campioni AB per paziente)
heat_car_pos_AB <- heat_df %>%
  filter(car_status == "CAR+", timepoint == "AB") %>%
  group_by(patient, functional_state) %>%
  summarise(prop = mean(prop), .groups = "drop") %>%
  mutate(functional_state = factor(functional_state, levels = FUNCTIONAL_ORDER))

if (nrow(heat_car_pos_AB) > 0) {
  p_heat_AB <- ggplot(heat_car_pos_AB,
                      aes(x = patient, y = functional_state, fill = prop)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.1f%%", 100 * prop)),
              size = 3.5, color = "black") +
    scale_fill_gradient(low = "white", high = "#3C5488",
                        labels = percent_format(), name = "Proporzione") +
    labs(title = "CAR+ cells in AB — Proporzione stati funzionali",
         subtitle = "Heatmap inter-paziente (agnostico CD4/CD8)",
         x = "Paziente", y = NULL) +
    theme_classic(base_size = 12)

  ggsave(file.path(out_dir, "Q1b_ALL_heatmap_CARpos_AB.png"),
         p_heat_AB, width = 6, height = 4, dpi = 300)
}

# ============================================================
# SEZIONE 2: MODULE SCORES CAR+ vs CAR- (agnostico per stato)
# ============================================================
section("Sezione 2: Module scores CAR+ vs CAR-")

score_cols <- paste0("score_", names(SIGNATURES))

for (patient in names(PATIENT_MAP)) {
  cat(sprintf("\nPaziente: %s\n", patient))

  # Usa campioni I (confronto baseline)
  objs_I <- lapply(PATIENT_MAP[[patient]]$I, function(s) {
    if (!s %in% names(all_samples)) return(NULL)
    obj <- all_samples[[s]]
    obj$car_status <- get_car_status(obj, s)
    obj$sample     <- s
    obj
  })
  objs_I <- Filter(Negate(is.null), objs_I)
  if (length(objs_I) == 0) next

  # Merge se più campioni (non succede mai per I, ma difensivo)
  merged_I <- if (length(objs_I) == 1) objs_I[[1]] else
    merge(objs_I[[1]], y = objs_I[-1])

  meta_I <- merged_I@meta.data
  avail_scores <- intersect(score_cols, colnames(meta_I))

  if (length(avail_scores) == 0) {
    cat(sprintf("  [WARN] Nessun module score trovato per %s\n", patient))
    next
  }

  # Violin plots: un pannello per firma, CAR+ vs CAR-
  # Tutti i tipi cellulari T insieme (agnostico)
  meta_I_t <- meta_I %>%
    mutate(
      fs = map_to_functional_state(as.character(cell_type))
    ) %>%
    filter(fs != "Other")

  vln_plots <- lapply(avail_scores, function(sc) {
    sig_name <- gsub("score_", "", sc)
    df_vln <- meta_I_t %>%
      select(car_status, score = !!sym(sc)) %>%
      filter(!is.na(score))

    n_pos <- sum(df_vln$car_status == "CAR+")
    n_neg <- sum(df_vln$car_status == "CAR-")

    # Mann-Whitney test
    mw <- tryCatch({
      wt <- wilcox.test(score ~ car_status, data = df_vln, exact = FALSE)
      sprintf("p=%.3g", wt$p.value)
    }, error = function(e) "p=NA")

    ggplot(df_vln, aes(x = car_status, y = score, fill = car_status)) +
      geom_violin(trim = TRUE, alpha = 0.8, color = "white") +
      geom_boxplot(width = 0.15, outlier.shape = NA,
                   fill = "white", alpha = 0.6) +
      scale_fill_manual(values = c("CAR+" = "#E64B35", "CAR-" = "#4DBBD5"),
                        guide = "none") +
      labs(title = sig_name,
           subtitle = sprintf("CAR+ n=%d | CAR- n=%d | %s", n_pos, n_neg, mw),
           x = NULL, y = "Module score") +
      theme_classic(base_size = 10) +
      theme(plot.title = element_text(size = 9, face = "bold"),
            plot.subtitle = element_text(size = 7))
  })

  p_vln <- wrap_plots(vln_plots, ncol = 4) +
    plot_annotation(
      title    = sprintf("%s (I) — Module scores CAR+ vs CAR-", patient),
      subtitle = "Tutte le cellule T insieme, indipendentemente da CD4/CD8",
      theme    = theme(plot.title    = element_text(face = "bold"),
                       plot.subtitle = element_text(size = 9, color = "gray40"))
    )

  ggsave(file.path(out_dir,
                   sprintf("Q1b_%s_I_module_scores_CARpos_vs_CARneg.png", patient)),
         p_vln, width = 14, height = 8, dpi = 300)
  cat(sprintf("  Salvato: Q1b_%s_I_module_scores_CARpos_vs_CARneg.png\n", patient))
}

# ============================================================
# SEZIONE 3: SUB-CLUSTERING DELLE SOLE CELLULE T
#            (risoluzione fine, indipendente da CD4/CD8)
# ============================================================
section("Sezione 3: Sub-clustering T cells (agnostico)")

# Raccogli tutte le cellule T da tutti i campioni I
all_T_objs <- list()

for (patient in names(PATIENT_MAP)) {
  for (sname in PATIENT_MAP[[patient]]$I) {
    if (!sname %in% names(all_samples)) next
    obj <- all_samples[[sname]]
    car_vec <- get_car_status(obj, sname)

    # Aggiungi metadati
    obj$car_status <- car_vec
    obj$sample     <- sname
    obj$patient    <- patient

    # Filtra solo cellule T annotate (esclude B, NK, mieloidi, ecc.)
    meta <- obj@meta.data
    t_mask <- map_to_functional_state(as.character(meta$cell_type)) != "Other"
    cat(sprintf("  %s: %d cellule T su %d totali\n",
                sname, sum(t_mask), length(t_mask)))

    if (sum(t_mask) < 30) {
      cat(sprintf("    [SKIP] Troppo poche cellule T (%d)\n", sum(t_mask)))
      next
    }
    all_T_objs[[sname]] <- subset(obj, cells = which(t_mask))
  }
}

if (length(all_T_objs) == 0) {
  cat("[WARN] Nessun campione I con cellule T sufficienti per sub-clustering.\n")
} else {

  cat(sprintf("\nMerge di %d oggetti T-cell per sub-clustering\n",
              length(all_T_objs)))

  merged_T <- if (length(all_T_objs) == 1) {
    all_T_objs[[1]]
  } else {
    merge(all_T_objs[[1]], y = all_T_objs[-1],
          add.cell.ids = names(all_T_objs))
  }

  cat(sprintf("  Totale cellule T: %d\n", ncol(merged_T)))

  # Normalizzazione + PCA + UMAP
  # Usa SCTransform se disponibile, altrimenti NormalizeData
  merged_T <- NormalizeData(merged_T, verbose = FALSE)
  merged_T <- FindVariableFeatures(merged_T, nfeatures = 2000, verbose = FALSE)

  # Rimuovi geni CAR-associati (transgene) se presenti nel feature space
  # Questi potrebbero creare un cluster artificiale
  var_genes <- VariableFeatures(merged_T)
  car_genes_to_remove <- grep("^CAR|^GD2|^FMC63|SCFV|transgene",
                               var_genes, ignore.case = TRUE, value = TRUE)
  if (length(car_genes_to_remove) > 0) {
    cat(sprintf("  Rimossi %d geni CAR-like dai variable features\n",
                length(car_genes_to_remove)))
    var_genes <- setdiff(var_genes, car_genes_to_remove)
    VariableFeatures(merged_T) <- var_genes
  }

  merged_T <- ScaleData(merged_T, verbose = FALSE)
  merged_T <- RunPCA(merged_T, npcs = 30, verbose = FALSE)

  # Determina n PC ottimali con gomito
  n_pcs <- 20  # conservativo per n cellule medio-basso

  merged_T <- RunUMAP(merged_T, dims = 1:n_pcs, verbose = FALSE,
                      min.dist = 0.3, spread = 1.0)
  merged_T <- FindNeighbors(merged_T, dims = 1:n_pcs, verbose = FALSE)

  # Clustering a risoluzione fine (0.3-0.5 per campioni piccoli)
  merged_T <- FindClusters(merged_T, resolution = 0.4, verbose = FALSE)

  cat(sprintf("  Sub-cluster identificati: %s\n",
              paste(levels(merged_T$seurat_clusters), collapse = ", ")))

  # ── UMAP plots ──────────────────────────────────────────────

  # 1. Colorato per cluster
  p_umap_clust <- DimPlot(merged_T, group.by = "seurat_clusters",
                           label = TRUE, label.size = 4,
                           pt.size = 0.8) +
    labs(title = "Sub-cluster T cells (prodotto infusione I)",
         subtitle = "Risoluzione 0.4, agnostico CD4/CD8") +
    NoLegend()

  # 2. Colorato per stato CAR
  p_umap_car <- DimPlot(merged_T,
                         cells.highlight = WhichCells(merged_T,
                                            expression = car_status == "CAR+"),
                         cols.highlight = "#E64B35",
                         cols = "lightgrey",
                         pt.size = 0.8) +
    labs(title = "CAR+ cells (rosso) nei sub-cluster T",
         subtitle = "Prodotto infusione I") +
    theme(legend.position = "none")

  # 3. Colorato per paziente
  p_umap_pat <- DimPlot(merged_T, group.by = "patient",
                         pt.size = 0.8,
                         cols = c(Bo = "#E64B35", Ca = "#4DBBD5",
                                  Me = "#00A087")) +
    labs(title = "Distribuzione per paziente",
         subtitle = "Prodotto infusione I")

  # 4. Colorato per tipo cellulare originale
  p_umap_type <- DimPlot(merged_T, group.by = "cell_type",
                          pt.size = 0.6, label = FALSE) +
    labs(title = "Tipo cellulare annotato",
         subtitle = "Sovrapposizione con sub-cluster")

  p_umap_combined <- (p_umap_clust | p_umap_car) /
                     (p_umap_pat   | p_umap_type)

  ggsave(file.path(out_dir, "Q1b_ALL_Tcell_subcluster_UMAP.png"),
         p_umap_combined, width = 14, height = 12, dpi = 300)
  cat("  Salvato: Q1b_ALL_Tcell_subcluster_UMAP.png\n")

  # ── Feature plots module scores sul UMAP ───────────────────
  # Aggiungi module scores all'oggetto merged_T
  for (sig_name in names(SIGNATURES)) {
    genes_ok <- intersect(SIGNATURES[[sig_name]], rownames(merged_T))
    if (length(genes_ok) < 2) next
    score_col <- paste0("score_", sig_name)
    merged_T <- AddModuleScore(merged_T,
                               features = list(genes_ok),
                               name     = score_col,
                               seed     = 42)
    old_col <- paste0(score_col, "1")
    merged_T@meta.data[[score_col]] <- merged_T@meta.data[[old_col]]
    merged_T@meta.data[[old_col]]   <- NULL
  }

  fp_plots <- lapply(names(SIGNATURES), function(sig_name) {
    sc <- paste0("score_", sig_name)
    if (!sc %in% colnames(merged_T@meta.data)) return(NULL)
    FeaturePlot(merged_T, features = sc, pt.size = 0.5,
                order = TRUE, min.cutoff = "q10") +
      scale_colour_gradientn(colours = c("lightgrey", "#E64B35")) +
      labs(title = sig_name) +
      theme(legend.key.size = unit(0.4, "cm"),
            plot.title = element_text(size = 9, face = "bold"))
  })
  fp_plots <- Filter(Negate(is.null), fp_plots)

  p_fp <- wrap_plots(fp_plots, ncol = 4) +
    plot_annotation(
      title = "Module scores sul UMAP sub-cluster T cells (I)",
      theme = theme(plot.title = element_text(face = "bold"))
    )

  ggsave(file.path(out_dir, "Q1b_ALL_Tcell_subcluster_module_scores.png"),
         p_fp, width = 16, height = 10, dpi = 300)
  cat("  Salvato: Q1b_ALL_Tcell_subcluster_module_scores.png\n")

  # ── Marker per sub-cluster ──────────────────────────────────
  section("Marker per sub-cluster (FindAllMarkers)")

  # Seurat v5: i layer devono essere uniti prima di FindAllMarkers
  merged_T <- JoinLayers(merged_T)

  Idents(merged_T) <- "seurat_clusters"
  cluster_markers <- FindAllMarkers(
    merged_T,
    only.pos  = TRUE,
    min.pct   = 0.25,
    logfc.threshold = 0.25,
    test.use  = "wilcox",
    verbose   = FALSE
  )

  # Gestione caso dataframe vuoto o senza la colonna attesa
  if (nrow(cluster_markers) == 0 || !"p_val_adj" %in% colnames(cluster_markers)) {
    cat("  [WARN] FindAllMarkers non ha restituito marker significativi.\n")
    cluster_markers <- data.frame()
    top_markers     <- data.frame()
  } else {
    top_markers <- cluster_markers %>%
      filter(p_val_adj < 0.05) %>%
      group_by(cluster) %>%
      slice_max(order_by = avg_log2FC, n = 10) %>%
      ungroup()
    cat(sprintf("  Marker significativi: %d\n",
                nrow(cluster_markers[cluster_markers$p_val_adj < 0.05, ])))
  }

  # DotPlot top 5 per cluster
  top5_genes <- if (nrow(top_markers) > 0 && "gene" %in% colnames(top_markers)) {
    top_markers %>%
      group_by(cluster) %>%
      slice_max(avg_log2FC, n = 5) %>%
      pull(gene) %>%
      unique()
  } else {
    character(0)
  }

  if (length(top5_genes) > 0) {
    p_dot <- DotPlot(merged_T, features = top5_genes,
                     group.by = "seurat_clusters") +
      RotatedAxis() +
      scale_color_gradientn(colours = c("lightgrey", "#E64B35")) +
      labs(title = "Top marker per sub-cluster T cells",
           subtitle = "Prodotto infusione I — agnostico CD4/CD8") +
      theme(axis.text.x = element_text(size = 8))

    ggsave(file.path(out_dir, "Q1b_ALL_subcluster_dotplot_markers.png"),
           p_dot, width = max(12, length(top5_genes) * 0.5), height = 6, dpi = 300)
    cat("  Salvato: Q1b_ALL_subcluster_dotplot_markers.png\n")
  }

  # ── Proporzione CAR+ per sub-cluster ───────────────────────
  clust_car_df <- merged_T@meta.data %>%
    group_by(seurat_clusters, car_status) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(seurat_clusters) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup()

  p_clust_car <- ggplot(
    clust_car_df %>% filter(car_status == "CAR+"),
    aes(x = seurat_clusters, y = prop, fill = seurat_clusters)
  ) +
    geom_col(show.legend = FALSE) +
    geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", 100 * prop, n)),
              vjust = -0.3, size = 3.5) +
    scale_y_continuous(labels = percent_format(),
                       expand = expansion(mult = c(0, 0.15))) +
    labs(title = "Proporzione CAR+ per sub-cluster",
         subtitle = "Prodotto infusione I (tutte le cellule T, agnostico CD4/CD8)",
         x = "Sub-cluster", y = "% CAR+") +
    theme_classic(base_size = 12)

  ggsave(file.path(out_dir, "Q1b_ALL_CARpos_proportion_per_subcluster.png"),
         p_clust_car, width = 8, height = 5, dpi = 300)
  cat("  Salvato: Q1b_ALL_CARpos_proportion_per_subcluster.png\n")

  # Salva markers in Excel
  wb <- createWorkbook()
  addWorksheet(wb, "SubclusterMarkers")
  writeData(wb, "SubclusterMarkers", cluster_markers)
  saveWorkbook(wb, file.path(out_dir, "Q1b_subcluster_markers.xlsx"),
               overwrite = TRUE)
  cat("  Salvato: Q1b_subcluster_markers.xlsx\n")
}

# ============================================================
# SEZIONE 4: PROPORZIONI STATI FUNZIONALI PER PAZIENTE
#            CONFRONTO DIRETTO I vs AB (solo CAR+)
# ============================================================
section("Sezione 4: Confronto I vs AB per stati funzionali — solo CAR+")

if (length(all_props) > 0) {
  combined_df <- bind_rows(all_props, .id = "patient") %>%
    filter(car_status == "CAR+") %>%
    group_by(patient, timepoint, functional_state) %>%
    summarise(prop = mean(prop), .groups = "drop") %>%
    mutate(
      functional_state = factor(functional_state, levels = FUNCTIONAL_ORDER),
      timepoint        = factor(timepoint, levels = c("I", "AB"))
    )

  p_all_patients <- ggplot(combined_df,
                           aes(x = timepoint, y = prop, fill = functional_state)) +
    geom_col(width = 0.7, color = "white", linewidth = 0.3) +
    facet_wrap(~ patient, ncol = 3) +
    scale_fill_manual(values = STATE_PALETTE, drop = FALSE,
                      name = "Stato funzionale") +
    scale_y_continuous(labels = percent_format(), expand = c(0, 0)) +
    labs(
      title    = "CAR+ cells — Composizione stati funzionali: I vs AB",
      subtitle = "Approccio agnostico per lineage CD4/CD8 | Proporzioni dentro CAR+",
      x = "Timepoint", y = "Proporzione"
    ) +
    theme_classic(base_size = 12) +
    theme(
      strip.background = element_rect(fill = "#F0F0F0"),
      strip.text       = element_text(face = "bold", size = 12),
      legend.position  = "right"
    )

  ggsave(file.path(out_dir, "Q1b_ALL_CARpos_I_vs_AB_functional_states.png"),
         p_all_patients, width = 12, height = 5, dpi = 300)
  cat("Salvato: Q1b_ALL_CARpos_I_vs_AB_functional_states.png\n")
}

# ============================================================
# EXPORT EXCEL SUMMARY
# ============================================================
section("Export Excel")

if (length(all_props) > 0) {
  wb <- createWorkbook()
  summary_df <- bind_rows(all_props, .id = "patient")
  addWorksheet(wb, "Functional_Props")
  writeData(wb, "Functional_Props", summary_df)
  saveWorkbook(wb, file.path(out_dir, "Q1b_functional_summary.xlsx"),
               overwrite = TRUE)
  cat("Salvato: Q1b_functional_summary.xlsx\n")
}

section("COMPLETATO")
cat(sprintf("Output in: %s\n", out_dir))
cat("\nFile prodotti:\n")
for (f in list.files(out_dir, full.names = FALSE)) {
  cat(sprintf("  %s\n", f))
}
