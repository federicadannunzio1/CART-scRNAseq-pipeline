# ============================================================
#  9b – Classificazione CAR-T a livello di singola cellula
#
#  Motivazione: lo script precedente (9_cd4_threshold_reannotation)
#  ha dimostrato che abbassare la soglia CD4 nella media di cluster
#  non cambia nulla. Il problema è che le cellule CAR-T CD4+ sono
#  distribuite in cluster a dominanza CD8+ dove la MEDIA di cluster
#  oscura il segnale CD4+. Serve operare a livello di singola cellula.
#
#  Strategia:
#    1. Per ogni cellula IS_CAR_ALLIN_scREP nei campioni AB:
#       a. Estrai espressione normalizzata di CD4, CD8A, CD8B via FetchData
#       b. Usa i module score PER CELLULA già in metadata
#          (naive, cytotox, treg, th1/2/17/tfh, prolif, effector)
#    2. Applica una gerarchia di classificazione per singola cellula:
#       1. Treg:           treg_score > 0.15
#       2. Proliferating:  prolif_score > 0.10  → CD4/CD8 via marker
#       3. Cytotoxic CD8+: cytotox_score > 0.15
#       4. Th subtypes:    th1/th2/th17/tfh elevati → CD4+
#       5. CD8 detected:   max(CD8A, CD8B) > soglia_expr
#       6. CD4 detected:   CD4 > soglia_expr
#       7. Memory T:       effector_score alto, cytotox basso
#       8. Fallback:       naive vs cytotox
#    3. Confronta annotazione cluster-level vs single-cell per CAR-T
#    4. Testa soglie di espressione genica: 0, 0.1, 0.5, 1.0
#    5. Salva tabella e grafici comparativi
#
#  NB: per i campioni I l'annotazione manuale è già a buona risoluzione,
#  quindi questa analisi è focalizzata sui campioni AB post-infusione.
# ============================================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(readxl)
library(writexl)
library(scales)

# ── PERCORSI ──────────────────────────────────────────────────
base_dir <- path.expand("~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/CART/Code/Seurat_analysis/")
rds_path <- file.path(base_dir, "2_annotation", "all_samples_annotated_COMPLETE.rds")
out_dir  <- file.path(base_dir, "9_CAR_final", "res")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Soglie di espressione genica (log-normalizzata) da testare per CD4/CD8
EXPR_THRESHOLDS <- c(0, 0.1, 0.5, 1.0)

section <- function(title)
  cat(paste0("\n", strrep("=", 65), "\n  ", title,
             "\n", strrep("=", 65), "\n"))

# ── Classificazione CD4/CD8 aggregata ─────────────────────────
group_label <- function(lbl) {
  cd4_types <- c("Naive CD4+ T cells","Th1 cells","Th2 cells","Th17 cells",
                 "Tfh cells","Effector CD4+ T cells","Tregs",
                 "Proliferating CD4+ T cells")
  cd8_types <- c("Cytotoxic CD8+ T cells","Naive CD8+ T cells",
                 "Proliferating CD8+ T cells")
  mem_types <- c("Memory T cells")
  if (is.na(lbl))           return("Other")
  if (lbl %in% cd4_types)  return("CD4+")
  if (lbl %in% cd8_types)  return("CD8+")
  if (lbl %in% mem_types)  return("Memory T")
  return("Other")
}

# ============================================================
# 1. CARICAMENTO DATI
# ============================================================
section("Caricamento dati")

cat(sprintf("rds_path esiste: %s\n", file.exists(rds_path)))
all_samples <- readRDS(rds_path)

ab_names <- c("Ca_blood_AB", "Ca_bone_AB", "Bo_blood_AB", "Bo_bone_AB", "Me_bone_AB")
i_names  <- c("Bo_bone_I", "Ca_bone_I", "Me_bone_I")

cat("\nCampioni AB disponibili:\n")
for (nm in ab_names) {
  obj <- all_samples[[nm]]
  if (is.null(obj)) { cat(sprintf("  %-15s [MANCANTE]\n", nm)); next }
  n_car <- sum(!is.na(obj$IS_CAR_ALLIN_scREP) & obj$IS_CAR_ALLIN_scREP == "YES")
  cat(sprintf("  %-15s | %d celle totali | %d CAR+\n", nm, ncol(obj), n_car))
}

# ============================================================
# 2. FUNZIONE DI CLASSIFICAZIONE PER SINGOLA CELLULA
# ============================================================

classify_single_cell <- function(cd4_expr, cd8a_expr, cd8b_expr,
                                  treg_s, prolif_s, cytotox_s, naive_s,
                                  effector_s, th1_s, th2_s, th17_s, tfh_s,
                                  expr_thr) {

  cd8_expr <- max(cd8a_expr, cd8b_expr, na.rm = TRUE)

  # 1. Treg: treg_score dominante
  if (!is.na(treg_s) && treg_s > 0.15) {
    return("Tregs")
  }

  # 2. Proliferating: forte segnale di ciclo cellulare
  if (!is.na(prolif_s) && prolif_s > 0.10) {
    # Subtipizza con CD4/CD8 espressione + cytotox
    if (!is.na(cd8_expr) && cd8_expr > cd4_expr && cd8_expr > expr_thr) {
      return("Proliferating CD8+ T cells")
    } else if (!is.na(cd4_expr) && cd4_expr > expr_thr) {
      return("Proliferating CD4+ T cells")
    } else if (!is.na(cytotox_s) && cytotox_s > 0.05) {
      return("Proliferating CD8+ T cells")
    } else {
      return("Proliferating CD4+ T cells")  # default per T cell prodotto CAR
    }
  }

  # 3. Citotossico: cytotox_score elevato → CD8+
  #    (GZMB, PRF1, NKG7 sono quasi esclusivi di CD8 effettori)
  if (!is.na(cytotox_s) && cytotox_s > 0.15) {
    return("Cytotoxic CD8+ T cells")
  }

  # 4. T helper subtypes → questi score sono CD4-specifici
  #    (TBX21, GATA3, RORC, BCL6 → non espressi in CD8 normali)
  th_max <- max(c(th1_s, th2_s, th17_s, tfh_s), na.rm = TRUE)
  if (!is.na(th_max) && th_max > 0.08) {
    if      (!is.na(tfh_s)  && tfh_s  == th_max && tfh_s  > 0.08) return("Tfh cells")
    else if (!is.na(th17_s) && th17_s == th_max && th17_s > 0.08) return("Th17 cells")
    else if (!is.na(th1_s)  && th1_s  == th_max && th1_s  > 0.08) return("Th1 cells")
    else if (!is.na(th2_s)  && th2_s  == th_max && th2_s  > 0.08) return("Th2 cells")
  }

  # 5. Espressione genica CD8: CD8A o CD8B rilevati → CD8+
  if (!is.na(cd8_expr) && cd8_expr > expr_thr) {
    # CD8 rilevato: naive vs memoria/effettore
    if (!is.na(naive_s) && naive_s > 0.05) {
      return("Naive CD8+ T cells")
    } else if (!is.na(cytotox_s) && cytotox_s > 0) {
      return("Cytotoxic CD8+ T cells")
    } else {
      return("Memory T cells")
    }
  }

  # 6. Espressione genica CD4: CD4 rilevato → CD4+
  if (!is.na(cd4_expr) && cd4_expr > expr_thr) {
    if (!is.na(naive_s) && !is.na(effector_s) &&
        naive_s > effector_s && naive_s > 0.02) {
      return("Naive CD4+ T cells")
    } else {
      return("Effector CD4+ T cells")
    }
  }

  # 7. Nessun CD4/CD8 rilevato: usa i module score come guida
  #    CD8 effettore → cytotox_s
  #    CD4 naive/memoria → naive_s, effector_s
  if (!is.na(cytotox_s) && !is.na(naive_s)) {
    if (cytotox_s > naive_s && cytotox_s > 0) {
      return("Cytotoxic CD8+ T cells")
    }
    if (naive_s > 0.05) {
      # Ambiguo CD4/CD8 naive: classifichiamo come Memory T (CD4 dropout)
      return("Memory T cells")
    }
  }

  # 8. Effector senza lineage chiaro → Memory T (molto comune in CAR-T attivati)
  if (!is.na(effector_s) && effector_s > 0.05) {
    return("Memory T cells")
  }

  return("Memory T cells")  # fallback per CAR-T con profilo ambiguo
}

# ============================================================
# 3. CLASSIFICAZIONE PER SOGLIA DI ESPRESSIONE
# ============================================================
section("Classificazione single-cell CAR-T")

results_sc   <- list()
cell_details <- list()

for (nm in ab_names) {
  cat(sprintf("\n--- Campione: %s ---\n", nm))

  obj <- all_samples[[nm]]
  if (is.null(obj)) { cat("  [SKIP]\n"); next }

  # JoinLayers Seurat v5
  if (length(grep("^counts\\.", Layers(obj), value = TRUE)) > 0) {
    obj <- JoinLayers(obj)
  }

  # Seleziona cellule CAR+
  car_barcodes <- rownames(obj@meta.data)[
    !is.na(obj$IS_CAR_ALLIN_scREP) & obj$IS_CAR_ALLIN_scREP == "YES"
  ]
  n_car <- length(car_barcodes)
  cat(sprintf("  Cellule CAR+: %d\n", n_car))
  if (n_car == 0) { cat("  [SKIP] Nessuna cellula CAR+\n"); next }

  # Estrai espressione genica (log-normalizzata) per CD4/CD8
  expr_genes <- c("CD4", "CD8A", "CD8B")
  expr_genes_ok <- expr_genes[expr_genes %in% rownames(obj)]
  cat(sprintf("  Geni disponibili: %s\n", paste(expr_genes_ok, collapse=", ")))

  expr_df <- FetchData(obj, vars = expr_genes_ok, cells = car_barcodes,
                       layer = "data")
  # Aggiungi colonne mancanti come 0
  for (g in setdiff(expr_genes, colnames(expr_df)))
    expr_df[[g]] <- 0

  # Metadata con module score
  meta_car <- obj@meta.data[car_barcodes, ]

  # Colonne dei module score
  score_cols <- c("naive_score","effector_score","cytotox_score","treg_score",
                  "prolif_score","th1_score","th2_score","th17_score","tfh_score")
  # Verifica che ci siano
  score_cols_ok <- score_cols[score_cols %in% colnames(meta_car)]
  cat(sprintf("  Module score disponibili: %d/%d\n",
              length(score_cols_ok), length(score_cols)))

  gs <- function(col) {
    if (col %in% colnames(meta_car)) meta_car[[col]] else rep(0, nrow(meta_car))
  }

  # Classifica a ogni soglia di espressione
  for (expr_thr in EXPR_THRESHOLDS) {

    labels <- mapply(
      classify_single_cell,
      cd4_expr  = expr_df[["CD4"]],
      cd8a_expr = expr_df[["CD8A"]],
      cd8b_expr = expr_df[["CD8B"]],
      treg_s    = gs("treg_score"),
      prolif_s  = gs("prolif_score"),
      cytotox_s = gs("cytotox_score"),
      naive_s   = gs("naive_score"),
      effector_s= gs("effector_score"),
      th1_s     = gs("th1_score"),
      th2_s     = gs("th2_score"),
      th17_s    = gs("th17_score"),
      tfh_s     = gs("tfh_score"),
      MoreArgs  = list(expr_thr = expr_thr),
      SIMPLIFY  = TRUE
    )

    tbl <- as.data.frame(table(labels), stringsAsFactors = FALSE)
    colnames(tbl) <- c("cell_type_sc", "n_CAR")
    tbl$group   <- vapply(tbl$cell_type_sc, group_label, character(1L))
    tbl$sample  <- nm
    tbl$expr_thr <- expr_thr

    results_sc[[paste0(nm, "_", expr_thr)]] <- tbl

    n_cd4 <- sum(tbl$n_CAR[tbl$group == "CD4+"])
    n_cd8 <- sum(tbl$n_CAR[tbl$group == "CD8+"])
    n_mem <- sum(tbl$n_CAR[tbl$group == "Memory T"])
    cat(sprintf("  expr_thr=%.1f | CD4+: %d | CD8+: %d | Memory: %d | Other: %d\n",
                expr_thr, n_cd4, n_cd8, n_mem,
                sum(tbl$n_CAR[tbl$group == "Other"])))
  }

  # Dettaglio singola cellula (per expr_thr=0.1, il più bilanciato)
  labels_detail <- mapply(
    classify_single_cell,
    cd4_expr  = expr_df[["CD4"]],
    cd8a_expr = expr_df[["CD8A"]],
    cd8b_expr = expr_df[["CD8B"]],
    treg_s    = gs("treg_score"),
    prolif_s  = gs("prolif_score"),
    cytotox_s = gs("cytotox_score"),
    naive_s   = gs("naive_score"),
    effector_s= gs("effector_score"),
    th1_s     = gs("th1_score"),
    th2_s     = gs("th2_score"),
    th17_s    = gs("th17_score"),
    tfh_s     = gs("tfh_score"),
    MoreArgs  = list(expr_thr = 0.1),
    SIMPLIFY  = TRUE
  )

  cell_details[[nm]] <- data.frame(
    barcode      = car_barcodes,
    sample       = nm,
    cell_type_cluster = meta_car$cell_type,
    cell_type_sc = labels_detail,
    group_cluster = vapply(meta_car$cell_type, group_label, character(1L)),
    group_sc      = vapply(labels_detail, group_label, character(1L)),
    cd4_expr     = expr_df[["CD4"]],
    cd8a_expr    = expr_df[["CD8A"]],
    cd8b_expr    = expr_df[["CD8B"]],
    cytotox_score = gs("cytotox_score"),
    naive_score   = gs("naive_score"),
    treg_score    = gs("treg_score"),
    prolif_score  = gs("prolif_score"),
    stringsAsFactors = FALSE
  )
}

# ============================================================
# 4. CAMPIONI I – RIFERIMENTO (ANNOTAZIONE MANUALE)
# ============================================================
section("Campioni I – riferimento")

i_results_sc <- list()
for (nm in i_names) {
  obj <- all_samples[[nm]]
  if (is.null(obj)) next
  meta_df <- obj@meta.data
  if (!"IS_CAR_ALLIN_scREP" %in% colnames(meta_df) ||
      !"cell_type" %in% colnames(meta_df)) next

  car_cells <- meta_df[!is.na(meta_df$IS_CAR_ALLIN_scREP) &
                         meta_df$IS_CAR_ALLIN_scREP == "YES", ]
  tbl <- as.data.frame(table(car_cells$cell_type), stringsAsFactors = FALSE)
  colnames(tbl) <- c("cell_type_sc", "n_CAR")
  tbl$group    <- vapply(tbl$cell_type_sc, group_label, character(1L))
  tbl$sample   <- nm
  tbl$expr_thr <- NA_real_
  i_results_sc[[nm]] <- tbl

  n_cd4 <- sum(tbl$n_CAR[tbl$group == "CD4+"])
  n_cd8 <- sum(tbl$n_CAR[tbl$group == "CD8+"])
  cat(sprintf("  %s | CD4+: %d | CD8+: %d\n", nm, n_cd4, n_cd8))
}

# ============================================================
# 5. AGGREGAZIONE E EXPORT
# ============================================================
section("Aggregazione e salvataggio")

all_sc <- bind_rows(results_sc)

# Sintesi per soglia (tutti AB)
summary_sc <- all_sc %>%
  group_by(expr_thr, group) %>%
  summarise(n_CAR = sum(n_CAR), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = n_CAR, values_fill = 0) %>%
  arrange(expr_thr) %>%
  mutate(
    Total_T  = rowSums(across(any_of(c("CD4+","CD8+","Memory T")))),
    pct_CD4  = round(100 * `CD4+` / pmax(Total_T, 1), 1),
    pct_CD8  = round(100 * `CD8+` / pmax(Total_T, 1), 1),
    pct_Mem  = round(100 * `Memory T` / pmax(Total_T, 1), 1)
  )

cat("\nSintesi single-cell (AB aggregati) per soglia espressione:\n")
print(as.data.frame(summary_sc))

# Sintesi per campione (expr_thr = 0.1)
summary_per_sample <- all_sc %>%
  filter(expr_thr == 0.1) %>%
  group_by(sample, group) %>%
  summarise(n_CAR = sum(n_CAR), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = n_CAR, values_fill = 0)
cat("\nPer campione (expr_thr=0.1):\n")
print(as.data.frame(summary_per_sample))

# Confronto cluster-level vs single-cell (expr_thr=0.1)
all_details <- bind_rows(cell_details)
if (nrow(all_details) > 0) {
  comparison <- all_details %>%
    group_by(sample, group_cluster, group_sc) %>%
    summarise(n = n(), .groups = "drop") %>%
    arrange(sample, group_cluster, desc(n))

  cat("\nConfronto cluster-level vs single-cell (expr_thr=0.1):\n")
  print(as.data.frame(comparison))
}

# I samples riferimento
i_summary_sc <- bind_rows(i_results_sc) %>%
  group_by(group) %>%
  summarise(n_CAR = sum(n_CAR), .groups = "drop")
cat("\nCampioni I (riferimento manuale):\n")
print(as.data.frame(i_summary_sc))

# Export Excel
out_excel <- file.path(out_dir, "CAR_singlecell_classification.xlsx")
sheets <- list(
  "Sintesi_per_soglia_expr"  = as.data.frame(summary_sc),
  "Per_campione_expr0.1"     = as.data.frame(summary_per_sample),
  "I_samples_riferimento"    = as.data.frame(bind_rows(i_results_sc))
)
if (nrow(all_details) > 0)
  sheets[["Confronto_cluster_vs_sc"]] <- as.data.frame(comparison)

writexl::write_xlsx(sheets, path = out_excel)
cat(sprintf("\nExcel → %s\n", out_excel))

# ============================================================
# 6. GRAFICI
# ============================================================
section("Grafici")

PALETTE_GROUP <- c(
  "CD4+"     = "#E63946",
  "CD8+"     = "#264653",
  "Memory T" = "#2A9D8F",
  "Other"    = "#AAAAAA"
)

# ── P5: Sintesi single-cell per soglia espressione ────────────
p5_data <- summary_sc %>%
  select(expr_thr, `CD4+`, `CD8+`, `Memory T`) %>%
  pivot_longer(cols = c(`CD4+`,`CD8+`,`Memory T`),
               names_to = "group", values_to = "n_CAR")

p5 <- ggplot(p5_data, aes(x = factor(expr_thr), y = n_CAR, fill = group)) +
  geom_col(position = "stack", width = 0.7) +
  scale_fill_manual(values = PALETTE_GROUP) +
  labs(
    title    = "CAR-T: classificazione single-cell per soglia di espressione CD4/CD8",
    subtitle = "Campioni AB – IS_CAR_ALLIN_scREP = YES",
    x        = "Soglia espressione genica (log-norm)",
    y        = "Numero cellule", fill = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, color = "gray40"))

ggsave(file.path(out_dir, "P5_singlecell_by_expr_threshold.png"),
       p5, width = 9, height = 6, dpi = 300, bg = "white")
cat("  P5 salvato\n")

# ── P6: % CD4 vs % CD8 al variare della soglia ──────────────
p6_data <- summary_sc %>%
  select(expr_thr, pct_CD4, pct_CD8, pct_Mem) %>%
  pivot_longer(cols = c(pct_CD4, pct_CD8, pct_Mem),
               names_to = "lineage", values_to = "pct") %>%
  mutate(lineage = recode(lineage,
    pct_CD4 = "CD4+", pct_CD8 = "CD8+", pct_Mem = "Memory T"))

p6 <- ggplot(p6_data, aes(x = expr_thr, y = pct,
                           color = lineage, group = lineage)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(values = PALETTE_GROUP) +
  labs(
    title  = "% CAR-T per lineage al variare della soglia di espressione",
    x      = "Soglia espressione genica (log-norm)",
    y      = "% cellule CAR", color = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave(file.path(out_dir, "P6_pct_lineage_vs_expr_threshold.png"),
       p6, width = 9, height = 6, dpi = 300, bg = "white")
cat("  P6 salvato\n")

# ── P7: Per campione (expr_thr=0.1) ──────────────────────────
p7_data <- all_sc %>%
  filter(expr_thr == 0.1) %>%
  group_by(sample, group) %>%
  summarise(n_CAR = sum(n_CAR), .groups = "drop") %>%
  filter(group %in% c("CD4+","CD8+","Memory T"))

p7 <- ggplot(p7_data, aes(x = group, y = n_CAR, fill = group)) +
  geom_col(width = 0.6) +
  facet_wrap(~ sample, nrow = 1) +
  scale_fill_manual(values = PALETTE_GROUP) +
  labs(
    title = "CAR-T single-cell classification per campione (expr_thr = 0.1)",
    x = NULL, y = "Numero cellule CAR", fill = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        strip.text = element_text(face = "bold"),
        legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(out_dir, "P7_singlecell_per_sample.png"),
       p7, width = 14, height = 6, dpi = 300, bg = "white")
cat("  P7 salvato\n")

# ── P8: Confronto cluster-level vs single-cell (tutti AB) ────
if (nrow(all_details) > 0) {
  p8_data <- all_details %>%
    select(sample, group_cluster, group_sc) %>%
    pivot_longer(cols = c(group_cluster, group_sc),
                 names_to = "method", values_to = "group") %>%
    group_by(method, group) %>%
    summarise(n = n(), .groups = "drop") %>%
    filter(group %in% c("CD4+","CD8+","Memory T")) %>%
    mutate(method = recode(method,
      group_cluster = "Cluster-level\n(metodo originale)",
      group_sc      = "Single-cell\n(questo script)"))

  p8 <- ggplot(p8_data, aes(x = group, y = n, fill = group)) +
    geom_col(width = 0.6) +
    facet_wrap(~ method) +
    scale_fill_manual(values = PALETTE_GROUP) +
    labs(
      title    = "Confronto: annotazione cluster-level vs single-cell per CAR-T",
      subtitle = "Campioni AB aggregati – expr_thr=0.1",
      x = NULL, y = "Numero cellule CAR", fill = NULL
    ) +
    theme_classic(base_size = 13) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
          strip.text = element_text(face = "bold"),
          legend.position = "none")

  ggsave(file.path(out_dir, "P8_cluster_vs_singlecell_comparison.png"),
         p8, width = 10, height = 6, dpi = 300, bg = "white")
  cat("  P8 salvato\n")
}

# ── P9: Distribuzione espressione CD4/CD8 nelle cellule CAR+ ─
if (nrow(all_details) > 0) {
  expr_long <- all_details %>%
    select(sample, group_sc, cd4_expr, cd8a_expr, cd8b_expr) %>%
    mutate(cd8_max = pmax(cd8a_expr, cd8b_expr)) %>%
    select(sample, group_sc, CD4 = cd4_expr, CD8_max = cd8_max) %>%
    pivot_longer(cols = c(CD4, CD8_max),
                 names_to = "gene", values_to = "expr") %>%
    filter(group_sc %in% c("CD4+","CD8+","Memory T"))

  p9 <- ggplot(expr_long, aes(x = group_sc, y = expr, fill = group_sc)) +
    geom_violin(alpha = 0.7, scale = "width") +
    geom_boxplot(width = 0.15, outlier.size = 0.5, fill = "white") +
    facet_wrap(~ gene, scales = "free_y") +
    scale_fill_manual(values = PALETTE_GROUP) +
    labs(
      title    = "Espressione CD4 e CD8 (log-norm) nelle cellule CAR+ per lineage assegnato",
      subtitle = "Tutti campioni AB – expr_thr=0.1",
      x = NULL, y = "Espressione (log-norm)", fill = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
          legend.position = "none")

  ggsave(file.path(out_dir, "P9_CD4_CD8_expr_by_lineage.png"),
         p9, width = 12, height = 6, dpi = 300, bg = "white")
  cat("  P9 salvato\n")
}

# ============================================================
# 7. RIEPILOGO FINALE
# ============================================================
section("RIEPILOGO FINALE")

cat("\n── Campioni AB: single-cell (expr_thr=0.1) ──\n")
sc_01 <- all_sc %>%
  filter(expr_thr == 0.1) %>%
  group_by(group) %>%
  summarise(n = sum(n_CAR), .groups = "drop")
print(as.data.frame(sc_01))

cat("\n── Campioni AB: cluster-level (originale) ──\n")
if (nrow(all_details) > 0) {
  cl_orig <- all_details %>%
    group_by(group_cluster) %>%
    summarise(n = n(), .groups = "drop")
  print(as.data.frame(cl_orig))
}

cat("\n── Campioni I: annotazione manuale (riferimento) ──\n")
print(as.data.frame(i_summary_sc))

cat(paste0(
  "\n", strrep("-",65), "\n",
  "INTERPRETAZIONE CHIAVE:\n",
  "  Cluster-level: usa la media del cluster → CD4 mRNA dropout\n",
  "                 diluce il segnale nei cluster misti\n",
  "  Single-cell:   ogni cellula classificata sui propri score\n",
  "                 → recupera cellule CD4+ in cluster CD8+-dominati\n",
  strrep("-",65), "\n"
))

cat(sprintf("\nOutput in: %s\n", out_dir))
for (f in c("CAR_singlecell_classification.xlsx",
            "P5_singlecell_by_expr_threshold.png",
            "P6_pct_lineage_vs_expr_threshold.png",
            "P7_singlecell_per_sample.png",
            "P8_cluster_vs_singlecell_comparison.png",
            "P9_CD4_CD8_expr_by_lineage.png"))
  cat(sprintf("  - %s\n", f))
