# ============================================================
#  Q3: Caratteristiche delle cellule CAR+ nel midollo osseo
#      post-infusione (AB) + confronto con prodotto di infusione (I)
#
#  Domanda biologica:
#    Come si comportano le CART cells nel midollo osseo
#    dopo l'infusione? Cosa cambia rispetto al prodotto I?
#    Ci sono segni di esaurimento, attivazione, perdita di stemness?
#
#  Approccio:
#    A) Module scores nelle cellule CAR+ del midollo AB
#       (stesso pannello di Q2 per confronto diretto)
#    B) Confronto diretto I vs bone_AB nelle CAR+:
#       - Variazione dei module scores
#       - DotPlot geni funzionali chiave
#    C) DEG: CAR+ in bone_AB vs CAR+ in I (per paziente)
#       Risponde a: cosa è cambiato trascrizionalmente nelle
#       CART cells tra il prodotto di infusione e il midollo?
#
#  Prerequisiti:
#    all_samples_annotated_COMPLETE_IS_CAR_REVISED.rds
#
#  Output in out_dir/Q3_CART_in_bone_AB/:
#    Q3_<paziente>_bone_AB_module_scores.png
#    Q3_<paziente>_I_vs_bone_AB_comparison.png
#    Q3_ALL_patients_CAR_scores_boxplot.png
#    Q3_DEG_boneAB_vs_I.xlsx
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
out_dir  <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/CART/Code/Seurat_analysis/10_CART_functional_analysis/Q3_CART_in_bone_AB/"
# ─────────────────────────────────────────────────────────────

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

section <- function(title)
  cat(paste0("\n", strrep("=", 65), "\n  ", title, "\n",
             strrep("=", 65), "\n"))

# ============================================================
# FIRME GENICHE (identiche a Q2 per confronto diretto)
# ============================================================
SIGNATURES <- list(
  Effector = c("GZMB","PRF1","NKG7","GNLY","GZMA","GZMK","FGFBP2","CX3CR1"),
  Memory_Stemness = c("TCF7","CCR7","SELL","IL7R","LEF1","KLF2","BCL2","FOXO1"),
  Exhaustion = c("PDCD1","LAG3","HAVCR2","TIGIT","TOX","TOX2","ENTPD1","CTLA4","BATF"),
  Activation = c("CD69","CD44","TNFRSF9","IL2RA","ICOS","CD38"),
  Proliferation = c("MKI67","TOP2A","PCNA","CCNB1","STMN1","UBE2C"),
  Naive = c("CCR7","SELL","IL7R","TCF7","LEF1","KLF2","LDHB","RCAN3"),
  Tpex_StemLike = c("TCF7","CXCR5","TOX","BCL6","SLAMF6","ID3"),
  Tex_Terminal = c("HAVCR2","TIGIT","LAG3","CD160","ENTPD1","PRDM1","ZEB2")
)

# Mappa paziente → campioni
BONE_AB_MAP <- list(
  Bo = list(I = "Bo_bone_I", bone_AB = "Bo_bone_AB"),
  Ca = list(I = "Ca_bone_I", bone_AB = "Ca_bone_AB"),
  Me = list(I = "Me_bone_I", bone_AB = "Me_bone_AB")
)

# ── Helper: colonna CAR ─────────────────────────────────────
get_car_status <- function(obj, sample_name) {
  meta <- obj@meta.data
  for (col in c("IS_CAR_ALLIN_scREP", "IS_CAR", "CAR")) {
    if (col %in% colnames(meta)) {
      vals <- as.character(meta[[col]])
      car_pos <- grepl("^(YES|TRUE|yes|true|1)$", vals)
      cat(sprintf("  %s | '%s': CAR+ = %d (%.1f%%)\n",
                  sample_name, col, sum(car_pos), 100*mean(car_pos)))
      return(ifelse(car_pos, "CAR+", "CAR-"))
    }
  }
  cat(sprintf("[WARN] %s: nessuna colonna CAR\n", sample_name))
  rep("CAR-", ncol(obj))
}

filter_genes <- function(genes, obj)
  genes[genes %in% rownames(obj)]

get_umap_key <- function(obj) {
  for (nm in c("umap","wnn.umap","umap.harmony","RNA.umap"))
    if (nm %in% names(obj@reductions)) return(nm)
  return(NULL)
}

# ============================================================
# CARICAMENTO
# ============================================================
section("Caricamento dati")

all_samples <- readRDS(rds_path)
cat("Campioni disponibili:", paste(names(all_samples), collapse = ", "), "\n")

# ============================================================
# STEP 1: AGGIUNTA CAR STATUS E MODULE SCORES
#         A TUTTI I CAMPIONI RILEVANTI (I + bone_AB)
# ============================================================
section("STEP 1 | CAR status + Module scores (I + bone_AB)")

add_scores <- function(obj, nm) {
  obj$CAR_status <- get_car_status(obj, nm)
  if (length(grep("^counts\\.", Layers(obj), value = TRUE)) > 0)
    obj <- JoinLayers(obj)

  for (sig_name in names(SIGNATURES)) {
    genes_ok <- filter_genes(SIGNATURES[[sig_name]], obj)
    if (length(genes_ok) < 3) next
    col_nm <- paste0("Score_", sig_name)
    obj <- AddModuleScore(obj, features = list(genes_ok),
                          name = col_nm, seed = 42)
    obj[[col_nm]]              <- obj[[paste0(col_nm, "1")]]
    obj[[paste0(col_nm, "1")]] <- NULL
  }
  obj
}

processed <- list()
for (pid in names(BONE_AB_MAP)) {
  pm <- BONE_AB_MAP[[pid]]
  for (role in c("I", "bone_AB")) {
    nm <- pm[[role]]
    if (!nm %in% names(all_samples)) {
      cat(sprintf("[WARN] %s (%s) non trovato\n", nm, role))
      next
    }
    key <- paste0(pid, "_", role)
    cat(paste0("\n── ", key, " (", nm, ") ──\n"))
    processed[[key]] <- add_scores(all_samples[[nm]], nm)
    processed[[key]]$patient  <- pid
    processed[[key]]$timepoint <- role  # "I" o "bone_AB"
    processed[[key]]$sample_nm <- nm
  }
}

# ============================================================
# STEP 2: CONFRONTO SCORE CAR+ : I vs bone_AB
# ============================================================
section("STEP 2 | Confronto module scores: I vs bone_AB nelle CAR+")

# Estrai metadata CAR+ da tutti i campioni processati
score_cols <- paste0("Score_", names(SIGNATURES))

meta_car_pos <- bind_rows(lapply(names(processed), function(key) {
  obj    <- processed[[key]]
  cols_ok <- intersect(score_cols, colnames(obj@meta.data))
  meta   <- obj@meta.data[obj$CAR_status == "CAR+", , drop = FALSE]
  if (nrow(meta) == 0) return(NULL)
  meta %>%
    select(any_of(c("patient", "timepoint", "sample_nm", "cell_type",
                    cols_ok))) %>%
    mutate(key = key)
}))

n_car_summary <- meta_car_pos %>%
  group_by(patient, timepoint) %>%
  summarise(n_CAR_pos = n(), .groups = "drop")
cat("\nCellule CAR+ per paziente e timepoint:\n")
print(n_car_summary)

if (nrow(meta_car_pos) == 0) {
  cat("[ERRORE] Nessuna cellula CAR+ trovata. Verifica la colonna CAR.\n")
  stop()
}

# ── Boxplot comparativo I vs bone_AB per paziente ────────────
score_long <- meta_car_pos %>%
  pivot_longer(
    cols      = starts_with("Score_"),
    names_to  = "Signature",
    values_to = "Score"
  ) %>%
  mutate(
    Signature = gsub("^Score_", "", Signature),
    timepoint = factor(timepoint, levels = c("I", "bone_AB")),
    patient   = factor(patient, levels = sort(unique(patient)))
  )

p_comparison <- ggplot(
  score_long,
  aes(x = timepoint, y = Score, fill = timepoint, color = timepoint)
) +
  geom_violin(alpha = 0.5, trim = TRUE, scale = "width") +
  geom_boxplot(width = 0.18, fill = "white", outlier.size = 0.3,
               outlier.alpha = 0.3, color = "black") +
  # Connetti medie per paziente (mostra tendenza)
  stat_summary(
    fun = mean,
    geom = "point", shape = 21, size = 2.5,
    fill = "white", color = "black", stroke = 0.8
  ) +
  scale_fill_manual(values = c("I" = "#457B9D", "bone_AB" = "#E63946"),
                    guide = "none") +
  scale_color_manual(values = c("I" = "#457B9D", "bone_AB" = "#E63946"),
                     guide = "none") +
  facet_grid(Signature ~ patient, scales = "free_y") +
  labs(
    title    = "Module scores nelle CAR+: prodotto I vs midollo osseo AB",
    subtitle = "Ogni colonna = un paziente | Blu = prodotto I | Rosso = midollo post-infusione",
    x = NULL, y = "Module score"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5, size = 12),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40", size = 9),
    strip.text.x  = element_text(face = "bold", size = 10),
    strip.text.y  = element_text(face = "bold", size = 8),
    axis.text.x   = element_text(angle = 30, hjust = 1, size = 9)
  )

n_sigs <- length(unique(score_long$Signature))
n_pats <- length(unique(score_long$patient))
out_cmp <- paste0(out_dir, "Q3_ALL_patients_I_vs_bone_AB_scores.png")
ggsave(out_cmp, p_comparison,
       width = max(10, n_pats * 4 + 2),
       height = max(12, n_sigs * 2 + 3),
       dpi = 300, bg = "white")
cat(paste0("  → ", out_cmp, "\n"))

# ── Test Mann-Whitney per ogni paziente × firma ───────────────
cat("\n[Variazioni module score I → bone_AB, solo cellule CAR+]\n")
cat("[Mann-Whitney per paziente – solo descrittivo con N=1/gruppo]\n\n")

comparison_tbl <- score_long %>%
  group_by(patient, Signature) %>%
  summarise(
    n_I       = sum(timepoint == "I"),
    n_bone_AB = sum(timepoint == "bone_AB"),
    mean_I    = round(mean(Score[timepoint == "I"], na.rm = TRUE), 4),
    mean_bAB  = round(mean(Score[timepoint == "bone_AB"], na.rm = TRUE), 4),
    delta     = round(mean_bAB - mean_I, 4),
    p_mw      = tryCatch({
      d_i   <- Score[timepoint == "I"]
      d_bab <- Score[timepoint == "bone_AB"]
      if (length(d_i) < 3 || length(d_bab) < 3) NA_real_
      else round(wilcox.test(d_bab, d_i, exact = FALSE)$p.value, 4)
    }, error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  mutate(
    direction = case_when(
      delta > 0.05  ~ "UP in bone_AB",
      delta < -0.05 ~ "DOWN in bone_AB",
      TRUE          ~ "stable"
    )
  ) %>%
  arrange(patient, Signature)

cat("Tabella delta module score (midollo AB vs I) nelle CAR+:\n")
print(as.data.frame(comparison_tbl), row.names = FALSE)

# ============================================================
# STEP 3: UMAP DEL MIDOLLO BONE_AB CON CAR OVERLAY
# ============================================================
section("STEP 3 | UMAP midollo AB con CAR overlay + module scores")

for (pid in names(BONE_AB_MAP)) {
  key <- paste0(pid, "_bone_AB")
  if (!key %in% names(processed)) next
  obj <- processed[[key]]

  umap_key <- get_umap_key(obj)
  if (is.null(umap_key)) {
    cat(sprintf("[SKIP] %s: nessuna UMAP\n", key))
    next
  }

  coords <- as.data.frame(Embeddings(obj, umap_key)[, 1:2])
  colnames(coords) <- c("UMAP1", "UMAP2")
  coords$car       <- obj$CAR_status
  coords$cell_type <- as.character(obj$cell_type)
  coords           <- coords[order(coords$car == "CAR+"), ]

  centroids <- coords %>%
    group_by(cell_type) %>%
    summarise(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2), .groups = "drop")

  n_pos <- sum(coords$car == "CAR+")
  n_neg <- sum(coords$car == "CAR-")

  p_umap <- ggplot(coords, aes(UMAP1, UMAP2)) +
    geom_point(data = coords[coords$car == "CAR-",],
               color = "#D3D3D3", size = 0.4, alpha = 0.4) +
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
      paste0(pid, " – midollo bone_AB – CAR+ overlay"),
      subtitle = paste0("Rosso: CAR+ (n=", n_pos, ")  Grigio: CAR- (n=", n_neg, ")")
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title    = element_text(face = "bold", hjust = 0.5, size = 12),
      plot.subtitle = element_text(hjust = 0.5, color = "gray40", size = 9),
      axis.text     = element_blank(), axis.ticks = element_blank()
    )

  # Module score overlay (Exhaustion e Memory)
  key_scores_to_plot <- intersect(
    c("Score_Exhaustion", "Score_Memory_Stemness",
      "Score_Effector", "Score_Proliferation"),
    colnames(obj@meta.data)
  )

  fp_list <- lapply(key_scores_to_plot, function(sc) {
    FeaturePlot(obj, features = sc, reduction = umap_key,
                pt.size = 0.4, order = TRUE,
                min.cutoff = "q05", max.cutoff = "q95") +
      scale_color_gradientn(
        colors = c("lightgrey", "#FFF176", "#FB8C00", "#B71C1C"),
        name   = gsub("Score_", "", sc)
      ) +
      ggtitle(gsub("Score_", "", sc)) +
      theme_classic(base_size = 9) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 9),
            axis.text  = element_blank(), axis.ticks = element_blank())
  })

  p_combined <- wrap_plots(c(list(p_umap), fp_list), ncol = 3) +
    plot_annotation(
      title = paste0(pid, " – midollo bone_AB – CAR+ e module scores"),
      theme = theme(plot.title = element_text(face = "bold", hjust = 0.5))
    )

  out_path <- paste0(out_dir, "Q3_", pid, "_bone_AB_UMAP_scores.png")
  ggsave(out_path, p_combined,
         width = 16, height = ceil((1 + length(fp_list)) / 3) * 6,
         dpi = 300, bg = "white")
  cat(paste0("  → ", out_path, "\n"))
}

ceil <- function(x) as.integer(x) + as.integer(x != round(x))

# ============================================================
# STEP 4: DEG – CAR+ in bone_AB vs CAR+ in I (per paziente)
# ============================================================
section("STEP 4 | DEG: CAR+ bone_AB vs CAR+ I")

cat("\n[INTERPRETAZIONE BIOLOGICA]\n")
cat("  I geni UP nelle CAR+ di bone_AB rispetto a I =\n")
cat("  geni indotti in vivo dopo l'infusione nel midollo.\n")
cat("  Esempi attesi:\n")
cat("  - UP: GZMB/PRF1 (attivazione effettrice in vivo)\n")
cat("  - DOWN: TCF7/SELL/CCR7 (perdita di stemness)\n")
cat("  - UP: PDCD1/LAG3 (esaurimento progressivo)\n\n")

wb_deg <- createWorkbook()
deg_bone_I <- list()

for (pid in names(BONE_AB_MAP)) {
  cat(paste0("── DEG paziente: ", pid, " ──\n"))
  key_I   <- paste0(pid, "_I")
  key_bAB <- paste0(pid, "_bone_AB")

  if (!key_I %in% names(processed) || !key_bAB %in% names(processed)) {
    cat(sprintf("  [SKIP] %s: campione I o bone_AB mancante\n", pid))
    next
  }

  obj_I   <- processed[[key_I]]
  obj_bAB <- processed[[key_bAB]]

  # Estrai solo cellule CAR+ da entrambi
  bc_I   <- Cells(obj_I)[obj_I$CAR_status == "CAR+"]
  bc_bAB <- Cells(obj_bAB)[obj_bAB$CAR_status == "CAR+"]

  cat(sprintf("  CAR+ in I: %d | CAR+ in bone_AB: %d\n",
              length(bc_I), length(bc_bAB)))

  if (length(bc_I) < 5 || length(bc_bAB) < 5) {
    cat("  [SKIP] Troppo poche cellule CAR+ in un gruppo (min 5)\n")
    next
  }

  # Subset
  obj_I_car   <- subset(obj_I, cells = bc_I)
  obj_bAB_car <- subset(obj_bAB, cells = bc_bAB)

  # Assicura stesso set di geni (intersezione)
  genes_common <- intersect(rownames(obj_I_car), rownames(obj_bAB_car))
  cat(sprintf("  Geni comuni: %d\n", length(genes_common)))

  # Merge temporaneo per FindMarkers
  obj_I_car$compare_group   <- paste0(pid, "_I_CAR+")
  obj_bAB_car$compare_group <- paste0(pid, "_boneAB_CAR+")

  tryCatch({
    merged <- merge(obj_I_car, obj_bAB_car,
                    add.cell.ids = c("I", "boneAB"))
    if (length(grep("^counts\\.", Layers(merged), value = TRUE)) > 0)
      merged <- JoinLayers(merged)

    Idents(merged) <- "compare_group"

    deg <- FindMarkers(
      merged,
      ident.1         = paste0(pid, "_boneAB_CAR+"),  # test: bone_AB
      ident.2         = paste0(pid, "_I_CAR+"),        # reference: I
      min.pct         = 0.10,
      logfc.threshold = 0.25,
      test.use        = "wilcox",
      verbose         = FALSE
    )
    deg$gene      <- rownames(deg)
    deg$patient   <- pid
    deg <- deg[order(deg$avg_log2FC, decreasing = TRUE), ]

    cat(sprintf("  Sig (p_adj<0.05): %d geni\n",
                sum(deg$p_val_adj < 0.05, na.rm = TRUE)))
    cat("  Top 5 UP in bone_AB (attivazione/esaurimento in vivo):\n")
    top_up <- head(deg[deg$avg_log2FC > 0 & deg$p_val_adj < 0.05, ], 5)
    if (nrow(top_up) > 0)
      print(top_up[, c("gene","avg_log2FC","p_val_adj","pct.1","pct.2")],
            row.names = FALSE)
    cat("  Top 5 DOWN in bone_AB (perdita vs I):\n")
    top_dn <- tail(deg[deg$p_val_adj < 0.05, ], 5)
    if (nrow(top_dn) > 0)
      print(top_dn[, c("gene","avg_log2FC","p_val_adj","pct.1","pct.2")],
            row.names = FALSE)

    deg_bone_I[[pid]] <- deg
    addWorksheet(wb_deg, substr(pid, 1, 31))
    writeData(wb_deg, substr(pid, 1, 31), deg)

    rm(merged)
  }, error = function(e) {
    cat(sprintf("  [ERRORE] %s\n", conditionMessage(e)))
  })
}

# Geni consistenti tra pazienti (bone_AB vs I)
if (length(deg_bone_I) >= 2) {
  all_deg <- bind_rows(deg_bone_I)
  consistent <- all_deg %>%
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

  cat("\nGeni consistentemente variati (bone_AB vs I, ≥2 pazienti):\n")
  print(as.data.frame(consistent), n = 30)

  addWorksheet(wb_deg, "Consistenti_multi_paz")
  writeData(wb_deg, "Consistenti_multi_paz", consistent)
  addWorksheet(wb_deg, "Comparison_table")
  writeData(wb_deg, "Comparison_table", comparison_tbl)
}

saveWorkbook(wb_deg, paste0(out_dir, "Q3_DEG_boneAB_vs_I_CARpos.xlsx"),
             overwrite = TRUE)
cat(paste0("\n  → Q3_DEG_boneAB_vs_I_CARpos.xlsx\n"))

# ============================================================
# STEP 5: DOTPLOT – GENI FUNZIONALI CHIAVE (I vs bone_AB)
# ============================================================
section("STEP 5 | DotPlot geni funzionali (I vs bone_AB, solo CAR+)")

KEY_GENES <- c(
  # Stemness/memory
  "TCF7", "CCR7", "IL7R", "SELL", "LEF1", "BCL2",
  # Effettore
  "GZMB", "PRF1", "NKG7", "GZMA", "CX3CR1", "GNLY",
  # Exhaustion
  "PDCD1", "LAG3", "HAVCR2", "TIGIT", "TOX", "ENTPD1",
  # Proliferazione
  "MKI67", "TOP2A",
  # Attivazione
  "CD69", "CD44", "CD38", "TNFRSF9"
)

for (pid in names(BONE_AB_MAP)) {
  key_I   <- paste0(pid, "_I")
  key_bAB <- paste0(pid, "_bone_AB")
  if (!all(c(key_I, key_bAB) %in% names(processed))) next

  obj_I   <- processed[[key_I]]
  obj_bAB <- processed[[key_bAB]]

  bc_I   <- Cells(obj_I)[obj_I$CAR_status   == "CAR+"]
  bc_bAB <- Cells(obj_bAB)[obj_bAB$CAR_status == "CAR+"]

  if (length(bc_I) < 5 || length(bc_bAB) < 5) next

  obj_I$compare_group   <- "I_CAR+"
  obj_bAB$compare_group <- "bone_AB_CAR+"

  tryCatch({
    merged <- merge(
      subset(obj_I, cells = bc_I),
      subset(obj_bAB, cells = bc_bAB),
      add.cell.ids = c("I", "boneAB")
    )
    if (length(grep("^counts\\.", Layers(merged), value = TRUE)) > 0)
      merged <- JoinLayers(merged)

    Idents(merged) <- "compare_group"
    genes_ok <- filter_genes(KEY_GENES, merged)

    if (length(genes_ok) >= 5) {
      p_dot <- DotPlot(merged, features = genes_ok,
                       cols = c("lightgrey", "#B71C1C")) +
        RotatedAxis() +
        ggtitle(paste0(pid, " – Geni funzionali: I_CAR+ vs bone_AB_CAR+")) +
        theme_classic(base_size = 10) +
        theme(
          plot.title  = element_text(face = "bold", hjust = 0.5, size = 11),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
        )

      out_dot <- paste0(out_dir, "Q3_", pid, "_dotplot_I_vs_bone_AB.png")
      ggsave(out_dot, p_dot,
             width = max(10, length(genes_ok) * 0.5 + 3),
             height = 5, dpi = 300, bg = "white")
      cat(paste0("  → ", out_dot, "\n"))
    }
    rm(merged)
  }, error = function(e) cat(sprintf("  [WARN] %s: %s\n", pid, conditionMessage(e))))
}

cat(paste0(
  "\n", strrep("=", 65), "\n",
  "  Q3 COMPLETATA\n\n",
  "  INTERPRETAZIONE ATTESA:\n",
  "  - Se TCF7/IL7R DOWN in bone_AB vs I:\n",
  "    perdita di stemness = le CART si sono differenziate\n",
  "    in effettrici in vivo (normale e desiderabile se non eccessivo).\n",
  "  - Se PDCD1/LAG3/HAVCR2 UP in bone_AB:\n",
  "    esaurimento progressivo in vivo = problema per persistenza.\n",
  "  - Se GZMB/PRF1 UP: attivazione effettrice in vivo.\n",
  "  - Se MKI67 UP: proliferazione post-infusione.\n\n",
  "  CONFRONTO CON Q2:\n",
  "  - Un prodotto I con TCF7 alto (Q2) che poi acquisisce\n",
  "    GZMB in bone_AB (Q3) = dinamica ottimale (stem→effector).\n",
  "  - TCF7 già basso in I + PDCD1 alto in bone_AB =\n",
  "    segnale sfavorevole (esaurimento precoce).\n",
  strrep("=", 65), "\n"
))
