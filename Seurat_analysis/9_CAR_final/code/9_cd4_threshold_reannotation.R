# ============================================================
#  9_CAR_FINAL – Analisi sensibilità soglia CD4 in auto_annotate
#
#  Obiettivo: diagnosticare perché IS_CAR_ALLIN_scREP è sbilanciata
#  verso CD8+ invece di CD4+, testando diverse soglie CD4 nella
#  funzione auto_annotate dei campioni AB.
#
#  Strategia (efficiente, NO re-run di AddModuleScore):
#    1. Carica all_samples_annotated_COMPLETE.rds
#       → per avere IS_CAR_ALLIN_scREP, seurat_clusters, sample per cella
#    2. Per ogni campione AB: estrai AverageExpression dei soli geni
#       lineage (CD4, CD8A/B, CD3D/E/G, NCAM1, CD19, CD14, CD16,
#       CD34, MKI67, TOP2A) → rapido, ~sec per campione
#    3. Carica gli score di modulo già calcolati dai file Excel in
#       2_annotation/AB_annotation/
#    4. Re-applica la gerarchia decisionale con CD4_THR variabile:
#       {0.30 (originale), 0.25, 0.20, 0.15, 0.10}
#    5. Mappa le nuove etichette alle singole cellule
#    6. Conta le cellule IS_CAR_ALLIN_scREP per cell_type / soglia / campione
#    7. Salva tabella Excel e grafici in 9_CAR_final/res/
#
#  NB: i campioni I hanno annotazione manuale (PIPELINE_1) → non
#  dipendono dalla soglia CD4 di auto_annotate. Vengono inclusi nel
#  riepilogo finale come controllo immutabile.
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
base_dir  <- path.expand("~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/CART/Code/Seurat_analysis/")
rds_path  <- file.path(base_dir, "2_annotation", "all_samples_annotated_COMPLETE.rds")
excel_dir <- file.path(base_dir, "2_annotation", "AB_annotation")
out_dir   <- file.path(base_dir, "9_CAR_final", "res")
# Verifica percorso
cat(sprintf("base_dir: %s\n", base_dir))
cat(sprintf("rds_path esiste: %s\n", file.exists(rds_path)))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Soglie da testare
CD4_THRESHOLDS <- c(0.30, 0.25, 0.20, 0.15, 0.10)

section <- function(title)
  cat(paste0("\n", strrep("=", 65), "\n  ", title,
             "\n", strrep("=", 65), "\n"))

# ============================================================
# 1. CARICAMENTO DATI
# ============================================================
section("Caricamento dati")

cat("Caricamento all_samples_annotated_COMPLETE.rds ...\n")
all_samples <- readRDS(rds_path)

# Nomi campioni AB e I
ab_names <- c("Ca_blood_AB", "Ca_bone_AB", "Bo_blood_AB", "Bo_bone_AB", "Me_bone_AB")
i_names  <- c("Bo_bone_I", "Ca_bone_I", "Me_bone_I")

cat("\nCampioni disponibili:\n")
for (nm in names(all_samples))
  cat(sprintf("  %-20s | %d celle | cols: %s\n",
              nm, ncol(all_samples[[nm]]),
              paste(colnames(all_samples[[nm]]@meta.data), collapse=", ")))

# ============================================================
# 2. HELPER: DECISION TREE con CD4_THR parametrizzabile
# ============================================================

# Geni lineage da calcolare con AverageExpression
LINEAGE_GENES <- c("CD4","CD8A","CD8B","CD3D","CD3E","CD3G",
                   "NCAM1","CD19","MS4A1","CD14","FCGR3A",
                   "CD34","MKI67","TOP2A")

# Estrai media per gene su cluster (equivalente a mean_genes in PIPELINE_2)
mean_gene_cluster <- function(avg_mat, genes, cluster_id) {
  g   <- genes[genes %in% rownames(avg_mat)]
  cid <- as.character(cluster_id)
  if (length(g) == 0 || !cid %in% colnames(avg_mat)) return(0)
  vals <- as.numeric(avg_mat[g, cid])
  if (length(vals) == 0 || all(is.na(vals))) return(0)
  mean(vals, na.rm = TRUE)
}

# Funzione di decisione per singolo cluster (versione parametrizzata)
decide_cluster <- function(cl, avg_mat, sc_row, CD4_THR) {
  # Valori lineage da AverageExpression
  s_cd3  <- mean_gene_cluster(avg_mat, c("CD3D","CD3E","CD3G"), cl)
  s_cd4  <- mean_gene_cluster(avg_mat, c("CD4"),                cl)
  s_cd8  <- mean_gene_cluster(avg_mat, c("CD8A","CD8B"),        cl)
  s_cd19 <- mean_gene_cluster(avg_mat, c("CD19","MS4A1"),       cl)
  s_cd56 <- mean_gene_cluster(avg_mat, c("NCAM1"),              cl)
  s_cd14 <- mean_gene_cluster(avg_mat, c("CD14"),               cl)
  s_cd16 <- mean_gene_cluster(avg_mat, c("FCGR3A"),             cl)
  s_cd34 <- mean_gene_cluster(avg_mat, c("CD34"),               cl)
  s_mk67 <- mean_gene_cluster(avg_mat, c("MKI67","TOP2A"),      cl)

  # Score di modulo dalla riga del dataframe Excel (sc_row)
  gs <- function(name) {
    col <- paste0(name, "_score")
    if (col %in% colnames(sc_row)) sc_row[[col]][1] else -99
  }
  s_naive   <- gs("naive");   s_effector <- gs("effector")
  s_cytotox <- gs("cytotox"); s_treg     <- gs("treg")
  s_prolif  <- gs("prolif")
  s_th1     <- gs("th1");     s_th2      <- gs("th2")
  s_th17    <- gs("th17");    s_tfh      <- gs("tfh")
  s_nk      <- gs("nk");      s_ilc      <- gs("ilc")
  s_nkt     <- gs("nkt");     s_gdt      <- gs("gdt")
  s_mait    <- gs("mait")
  s_bcell   <- gs("bcell");   s_bmem     <- gs("bmem")
  s_plasma  <- gs("plasma")
  s_mono14  <- gs("mono14");  s_mono16   <- gs("mono16")
  s_myeloid <- gs("myeloid")
  s_baso    <- gs("baso");    s_hspc     <- gs("hspc")
  s_erythro <- gs("erythro"); s_platelet <- gs("platelet")

  label <- NA_character_

  # ── 1. PLATELETS ─────────────────────────────────────────────
  if (s_platelet > 0.30 && s_cd3 < 0.3 && s_cd19 < 0.3) {
    label <- "Platelets"
  } else if (s_erythro > 0.30 && s_cd3 < 0.3 && s_cd19 < 0.3) {
    label <- "Erythroid cells"
  } else if (s_hspc > 0.15 && s_cd34 > 0.3 &&
             s_cd3 < 0.3 && s_cd19 < 0.3) {
    label <- "HSPC"
  } else if (s_plasma > 0.20 && s_cd3 < 0.3) {
    label <- "Plasma cells"
  } else if (s_cd19 > 0.5 && s_cd3 < 0.5) {
    label <- if (s_bmem > s_bcell && s_bmem > 0.05) "Memory B cells" else "B cells"
  } else if (s_baso > 0.15 && s_cd3 < 0.3 && s_cd19 < 0.3) {
    label <- "Basophils"
  } else if (s_cd14 > 0.5 && s_mono14 > 0.15 &&
             s_cd3 < 0.3 && s_cd56 < 0.3) {
    label <- "CD14 Monocytes"
  } else if (s_cd16 > 0.3 && s_mono16 > 0.10 &&
             s_cd3 < 0.3 && s_cd14 < 0.5) {
    label <- "CD16 Monocytes"
  } else if (s_myeloid > 0.35 &&
             s_cd3 < 0.3 && s_cd19 < 0.3 && s_cd56 < 0.3) {
    label <- "Myeloid cells"
  } else if (s_cd3 > 0.5 && s_cd56 > 0.3 && s_nkt > s_nk) {
    label <- "NKT cells"
  } else if (s_gdt > 0.05 && s_cd3 > 0.5 &&
             s_cd4 < 0.5 && s_cd8 < 0.5) {
    label <- "gamma-delta T cells"
  } else if (s_mait > 0.05 && s_cd3 > 0.5 &&
             s_mait > s_naive && s_mait > s_effector) {
    label <- "MAIT cells"
  } else if (s_ilc > 0.05 && s_cd3 < 0.5 &&
             s_cd56 < 0.3 && s_ilc > s_nk) {
    label <- "ILC"
  } else if (s_nk > 0.10 && s_cd3 < 0.8 &&
             s_cd4 < 0.5 && s_cd8 < 0.5) {
    label <- "NK cells"
  } else if (s_treg > 0.05 && s_treg > s_nk && s_cd4 > CD4_THR) {
    label <- "Tregs"
  } else if (s_mk67 > 0.5 || s_prolif > 0.05) {
    # Proliferating: logica CD4/CD8 invariata
    if (max(s_cd4, s_cd8) < 0.05) {
      label <- if (s_cytotox > 0.05) "Proliferating CD8+ T cells" else "Proliferating CD4+ T cells"
    } else if (s_cd8 > s_cd4) {
      label <- "Proliferating CD8+ T cells"
    } else if (s_cd4 > s_cd8) {
      label <- "Proliferating CD4+ T cells"
    } else {
      label <- if (s_cytotox > 0.05) "Proliferating CD8+ T cells" else "Proliferating CD4+ T cells"
    }
  # ── 17. CD8+ ────────────────────────────────────────────────
  } else if (s_cd8 > s_cd4 && s_cd8 > 0.3) {
    if (s_cytotox > s_naive && s_cytotox > 0.0) {
      label <- "Cytotoxic CD8+ T cells"
    } else if (s_naive > 0.1) {
      label <- "Naive CD8+ T cells"
    } else {
      label <- "Memory T cells"
    }
  # ── 18. CD4+ con soglia parametrizzata ──────────────────────
  } else if (s_cd4 > CD4_THR) {
    if (s_naive > s_effector && s_naive > 0.05) {
      label <- "Naive CD4+ T cells"
    } else if (s_tfh > 0.08 && s_tfh > s_naive) {
      label <- "Tfh cells"
    } else if (s_th17 > 0.06 && s_th17 > s_naive &&
               s_th17 > s_th1 && s_th17 > s_th2) {
      label <- "Th17 cells"
    } else if (s_th1 > 0.06 && s_th1 > s_naive && s_th1 > s_th2) {
      label <- "Th1 cells"
    } else if (s_th2 > 0.06 && s_th2 > s_naive) {
      label <- "Th2 cells"
    } else {
      label <- "Effector CD4+ T cells"
    }
  # ── 19. FALLBACK ────────────────────────────────────────────
  } else {
    sc_T <- c("Naive CD4+ T cells"    = s_naive,
              "Effector CD4+ T cells" = s_effector,
              "Cytotoxic CD8+ T cells"= s_cytotox,
              "Tregs"                 = s_treg,
              "NK cells"              = s_nk)
    label <- names(which.max(sc_T))
  }

  data.frame(
    cluster   = as.character(cl),
    label_new = label,
    s_cd4     = s_cd4,
    s_cd8     = s_cd8,
    s_cd3     = s_cd3,
    stringsAsFactors = FALSE
  )
}

# ============================================================
# 3. REANNOTAZIONE AB A DIVERSE SOGLIE
# ============================================================
section("Re-annotazione campioni AB a diverse soglie CD4")

# Raggruppa etichette CD4+ / CD8+ / altro per sintesi
group_label <- function(lbl) {
  cd4_types <- c("Naive CD4+ T cells","Th1 cells","Th2 cells","Th17 cells",
                 "Tfh cells","Effector CD4+ T cells","Tregs",
                 "Proliferating CD4+ T cells")
  cd8_types <- c("Cytotoxic CD8+ T cells","Naive CD8+ T cells",
                 "Proliferating CD8+ T cells")
  mem_types <- c("Memory T cells")
  if (lbl %in% cd4_types) return("CD4+")
  if (lbl %in% cd8_types) return("CD8+")
  if (lbl %in% mem_types) return("Memory T")
  return("Other")
}

# Lista per raccogliere risultati
results_list <- list()

for (nm in ab_names) {
  cat(sprintf("\n--- Campione: %s ---\n", nm))

  obj <- all_samples[[nm]]
  if (is.null(obj)) { cat("  [SKIP] Non trovato\n"); next }

  # JoinLayers se necessario (Seurat v5)
  if (length(grep("^counts\\.", Layers(obj), value = TRUE)) > 0) {
    obj <- JoinLayers(obj)
  }
  Idents(obj) <- "seurat_clusters"

  # ── Calcolo AverageExpression (geni lineage) ─────────────────
  genes_ok <- LINEAGE_GENES[LINEAGE_GENES %in% rownames(obj)]
  avg_raw  <- as.matrix(AverageExpression(obj, features = genes_ok,
                                          group.by = "seurat_clusters",
                                          assay = "RNA", slot = "data")$RNA)
  # Fix Seurat v5: rimuovi prefisso "g" dai nomi di colonna se presente
  colnames(avg_raw) <- gsub("^g", "", colnames(avg_raw))
  colnames(avg_raw) <- gsub("^RNA_snn_res\\.[0-9.]+_", "", colnames(avg_raw))
  storage.mode(avg_raw) <- "double"   # garantisce numerici
  cat(sprintf("  AverageExpression calcolata: %d geni x %d cluster\n",
              nrow(avg_raw), ncol(avg_raw)))

  # ── Carica score da Excel ─────────────────────────────────────
  excel_path <- file.path(excel_dir, paste0(nm, "_annotation_decisions.xlsx"))
  if (!file.exists(excel_path)) {
    cat(sprintf("  [SKIP] Excel non trovato: %s\n", excel_path)); next
  }
  sc_df <- as.data.frame(read_excel(excel_path))
  sc_df$cluster <- as.character(sc_df$cluster)
  cat(sprintf("  Score Excel caricati: %d cluster\n", nrow(sc_df)))

  clusters_chr <- sort(unique(as.character(obj$seurat_clusters)),
                       decreasing = FALSE)
  # Sort numerico
  clusters_chr <- as.character(sort(as.integer(clusters_chr)))

  # ── Per ogni soglia, re-annota ────────────────────────────────
  for (CD4_THR in CD4_THRESHOLDS) {

    new_labels <- data.frame()
    for (cl in clusters_chr) {
      sc_row <- sc_df[sc_df$cluster == cl, , drop = FALSE]
      if (nrow(sc_row) == 0) {
        cat(sprintf("  [WARN] Cluster %s non in Excel\n", cl))
        next
      }
      res <- decide_cluster(cl, avg_raw, sc_row, CD4_THR)
      new_labels <- rbind(new_labels, res)
    }

    # Mappa etichetta su ogni cellula
    ann_map   <- setNames(new_labels$label_new, new_labels$cluster)
    meta_df   <- obj@meta.data
    meta_df$cell_type_new <- unname(ann_map[as.character(meta_df$seurat_clusters)])

    # IS_CAR_ALLIN_scREP check
    if (!"IS_CAR_ALLIN_scREP" %in% colnames(meta_df)) {
      cat("  [WARN] IS_CAR_ALLIN_scREP non presente nei metadati\n")
      next
    }

    # Conta cellule CAR per cell_type_new
    car_cells <- meta_df[!is.na(meta_df$IS_CAR_ALLIN_scREP) &
                           meta_df$IS_CAR_ALLIN_scREP == "YES", ]

    tbl <- car_cells %>%
      group_by(cell_type_new) %>%
      summarise(n_CAR = n(), .groups = "drop") %>%
      mutate(
        sample      = nm,
        CD4_THR     = CD4_THR,
        group       = vapply(cell_type_new, group_label, character(1L)),
        n_total     = nrow(meta_df),
        n_CAR_total = nrow(car_cells)
      )

    results_list[[paste0(nm, "_", CD4_THR)]] <- tbl

    cat(sprintf("  CD4_THR=%.2f | CAR cells: %d → CD4+: %d | CD8+: %d | Memory: %d | Other: %d\n",
                CD4_THR,
                nrow(car_cells),
                sum(tbl$group == "CD4+"),
                sum(tbl$group == "CD8+"),
                sum(tbl$group == "Memory T"),
                sum(tbl$group == "Other")))
  }
}

# ============================================================
# 4. CAMPIONI I – annotazione manuale (riferimento immutabile)
# ============================================================
section("Campioni I – annotazione manuale (invariata)")

i_results <- list()
for (nm in i_names) {
  obj <- all_samples[[nm]]
  if (is.null(obj)) { cat(sprintf("  [SKIP] %s non trovato\n", nm)); next }

  meta_df <- obj@meta.data
  if (!"IS_CAR_ALLIN_scREP" %in% colnames(meta_df)) next
  if (!"cell_type" %in% colnames(meta_df)) next

  car_cells <- meta_df[!is.na(meta_df$IS_CAR_ALLIN_scREP) &
                         meta_df$IS_CAR_ALLIN_scREP == "YES", ]

  tbl <- car_cells %>%
    group_by(cell_type_new = cell_type) %>%
    summarise(n_CAR = n(), .groups = "drop") %>%
    mutate(
      sample   = nm,
      CD4_THR  = NA_real_,
      group    = vapply(cell_type_new, group_label, character(1L)),
      n_total  = nrow(meta_df),
      n_CAR_total = nrow(car_cells)
    )

  i_results[[nm]] <- tbl
  cat(sprintf("  %s | CAR total: %d → %s\n",
              nm, nrow(car_cells),
              paste(tbl$cell_type_new, tbl$n_CAR, sep=":", collapse=" | ")))
}

# ============================================================
# 5. AGGREGAZIONE RISULTATI E EXPORT
# ============================================================
section("Aggregazione e salvataggio risultati")

all_results <- bind_rows(results_list)

# ── Tabella sintesi per soglia (AB) ──────────────────────────
summary_by_thr <- all_results %>%
  group_by(CD4_THR, group) %>%
  summarise(n_CAR = sum(n_CAR), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = n_CAR, values_fill = 0) %>%
  arrange(CD4_THR) %>%
  mutate(
    Total_T   = rowSums(across(c(any_of(c("CD4+","CD8+","Memory T"))))),
    pct_CD4   = round(100 * `CD4+` / pmax(Total_T, 1), 1),
    pct_CD8   = round(100 * `CD8+` / pmax(Total_T, 1), 1)
  )

cat("\nSintesi AB per soglia CD4:\n")
print(as.data.frame(summary_by_thr))

# ── Tabella per campione e soglia ──────────────────────────────
summary_by_sample <- all_results %>%
  group_by(sample, CD4_THR, group) %>%
  summarise(n_CAR = sum(n_CAR), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = n_CAR, values_fill = 0)

# ── Tabella dettagliata per tipo cellulare ──────────────────────
detail_table <- all_results %>%
  select(sample, CD4_THR, cell_type_new, group, n_CAR, n_CAR_total) %>%
  arrange(sample, CD4_THR, desc(n_CAR))

# ── Campioni I (controllo) ──────────────────────────────────────
i_summary <- bind_rows(i_results) %>%
  group_by(group) %>%
  summarise(n_CAR = sum(n_CAR), .groups = "drop") %>%
  mutate(sample = "I_samples_all", CD4_THR = NA_real_)

# ── Export Excel multi-foglio ──────────────────────────────────
out_excel <- file.path(out_dir, "CAR_CD4_threshold_sensitivity.xlsx")
writexl::write_xlsx(
  list(
    "Sintesi_per_soglia"     = as.data.frame(summary_by_thr),
    "Per_campione_soglia"    = as.data.frame(summary_by_sample),
    "Dettaglio_celltype"     = as.data.frame(detail_table),
    "I_samples_riferimento"  = as.data.frame(bind_rows(i_results))
  ),
  path = out_excel
)
cat(sprintf("\nExcel salvato: %s\n", out_excel))

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

# ── Plot 1: CAR per soglia (tutti campioni AB aggregati) ────────
p1_data <- summary_by_thr %>%
  select(CD4_THR, `CD4+`, `CD8+`, `Memory T`) %>%
  pivot_longer(cols = c(`CD4+`, `CD8+`, `Memory T`),
               names_to = "group", values_to = "n_CAR") %>%
  mutate(CD4_THR_lbl = paste0("THR=", CD4_THR))

p1 <- ggplot(p1_data, aes(x = factor(CD4_THR), y = n_CAR, fill = group)) +
  geom_col(position = "stack", width = 0.7) +
  scale_fill_manual(values = PALETTE_GROUP) +
  labs(
    title    = "Distribuzione CAR-T per soglia CD4 (campioni AB aggregati)",
    subtitle = "IS_CAR_ALLIN_scREP = YES",
    x        = "Soglia CD4 (auto_annotate)",
    y        = "Numero cellule CAR",
    fill     = "Tipo cellulare"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40")
  )

ggsave(file.path(out_dir, "P1_CAR_by_CD4_threshold_stacked.png"),
       p1, width = 10, height = 6, dpi = 300, bg = "white")
cat("  P1 salvato\n")

# ── Plot 2: % CD4 vs % CD8 al variare della soglia ─────────────
p2_data <- summary_by_thr %>%
  select(CD4_THR, pct_CD4, pct_CD8) %>%
  pivot_longer(cols = c(pct_CD4, pct_CD8),
               names_to = "lineage", values_to = "pct") %>%
  mutate(lineage = recode(lineage, pct_CD4 = "CD4+", pct_CD8 = "CD8+"))

p2 <- ggplot(p2_data, aes(x = CD4_THR, y = pct, color = lineage, group = lineage)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0.30, linetype = "dashed", color = "gray50") +
  annotate("text", x = 0.30, y = max(p2_data$pct) * 1.05,
           label = "soglia\noriginale", size = 3.5, hjust = -0.1, color = "gray40") +
  scale_color_manual(values = c("CD4+" = "#E63946", "CD8+" = "#264653")) +
  scale_x_reverse(breaks = CD4_THRESHOLDS) +
  labs(
    title    = "% CAR-T CD4+ vs CD8+ al variare della soglia CD4",
    subtitle = "Campioni AB aggregati — IS_CAR_ALLIN_scREP = YES",
    x        = "Soglia CD4 (più bassa → più sensibile)",
    y        = "% su T cells CAR",
    color    = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40")
  )

ggsave(file.path(out_dir, "P2_pct_CD4_CD8_vs_threshold.png"),
       p2, width = 9, height = 6, dpi = 300, bg = "white")
cat("  P2 salvato\n")

# ── Plot 3: Per campione, soglia originale vs ottimale ──────────
thr_compare <- all_results %>%
  filter(CD4_THR %in% c(0.30, 0.15)) %>%
  group_by(sample, CD4_THR, group) %>%
  summarise(n_CAR = sum(n_CAR), .groups = "drop") %>%
  filter(group %in% c("CD4+","CD8+"))

p3 <- ggplot(thr_compare,
             aes(x = group, y = n_CAR, fill = group)) +
  geom_col(width = 0.6) +
  facet_grid(CD4_THR ~ sample, labeller = labeller(
    CD4_THR = function(x) paste0("THR=", x))) +
  scale_fill_manual(values = PALETTE_GROUP) +
  labs(
    title    = "CAR-T CD4+ vs CD8+ per campione: soglia 0.30 vs 0.15",
    x        = NULL, y = "Numero cellule CAR", fill = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title   = element_text(face = "bold", hjust = 0.5),
    strip.text   = element_text(face = "bold"),
    legend.position = "none"
  )

ggsave(file.path(out_dir, "P3_CAR_per_sample_thr030_vs_015.png"),
       p3, width = 14, height = 7, dpi = 300, bg = "white")
cat("  P3 salvato\n")

# ── Plot 4: Lollipop variazione CD4 per cluster (Bo_blood_AB) ──
# Mostra quanto cambia la classificazione cluster per cluster
sample_demo <- "Bo_blood_AB"
if (sample_demo %in% ab_names) {
  detail_demo <- detail_table %>%
    filter(sample == sample_demo, CD4_THR %in% c(0.30, 0.15)) %>%
    mutate(CD4_THR_lbl = paste0("THR=", CD4_THR))

  p4 <- ggplot(detail_demo,
               aes(x = reorder(cell_type_new, n_CAR), y = n_CAR,
                   fill = group)) +
    geom_col(width = 0.6) +
    facet_wrap(~ CD4_THR_lbl, ncol = 1) +
    scale_fill_manual(values = PALETTE_GROUP) +
    coord_flip() +
    labs(
      title = paste0(sample_demo, " – distribuzione CAR per cell_type"),
      x = NULL, y = "Numero cellule CAR", fill = NULL
    ) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))

  ggsave(file.path(out_dir, paste0("P4_", sample_demo, "_cluster_detail.png")),
         p4, width = 10, height = 8, dpi = 300, bg = "white")
  cat("  P4 salvato\n")
}

# ============================================================
# 7. RIEPILOGO FINALE
# ============================================================
section("RIEPILOGO FINALE")

cat("\n--- Campioni AB: cellule IS_CAR_ALLIN_scREP per soglia ---\n")
print(as.data.frame(summary_by_thr))

cat("\n--- Campioni I (annotazione manuale, invariata) ---\n")
i_ref <- bind_rows(i_results) %>%
  group_by(group) %>%
  summarise(n_CAR = sum(n_CAR), .groups = "drop")
print(as.data.frame(i_ref))

cat(paste0(
  "\n",
  strrep("-", 65), "\n",
  "INTERPRETAZIONE:\n",
  "  Soglia originale (0.30): baseline corrente\n",
  "  Soglie ridotte (≤0.20): cattura cluster CD4+ con mRNA dropout\n",
  "  Se il numero di CAR CD4+ aumenta abbassando la soglia,\n",
  "  questo conferma che la soglia 0.30 era troppo alta per\n",
  "  cellule CD4+ attivate (noto problema di mRNA dropout).\n",
  strrep("-", 65), "\n"
))

cat(paste0("\nOutput salvati in:\n  ", out_dir, "\n"))
cat("  - CAR_CD4_threshold_sensitivity.xlsx\n")
cat("  - P1_CAR_by_CD4_threshold_stacked.png\n")
cat("  - P2_pct_CD4_CD8_vs_threshold.png\n")
cat("  - P3_CAR_per_sample_thr030_vs_015.png\n")
cat("  - P4_*_cluster_detail.png\n")
