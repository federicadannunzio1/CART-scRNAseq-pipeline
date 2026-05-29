# ============================================================
#  Fix: riprende il sub-clustering da zero e completa
#  le sezioni mancanti dopo il crash (JoinLayers mancante)
# ============================================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(scales)
library(openxlsx)

rds_path <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/CART/Code/Seurat_analysis/2_annotation/all_samples_annotated_COMPLETE_IS_CAR_REVISED.rds"
out_dir  <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/CART/Code/Seurat_analysis/10_CART_functional_analysis/Q1b_functional_states/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

section <- function(title)
  cat(paste0("\n", strrep("=", 65), "\n  ", title, "\n",
             strrep("=", 65), "\n"))

PATIENT_MAP <- list(
  Bo = list(I = c("Bo_bone_I"), AB = c("Bo_blood_AB", "Bo_bone_AB")),
  Ca = list(I = c("Ca_bone_I"), AB = c("Ca_blood_AB", "Ca_bone_AB")),
  Me = list(I = c("Me_bone_I"), AB = c("Me_bone_AB"))
)

SIGNATURES <- list(
  Effector        = c("GZMB","PRF1","NKG7","GNLY","GZMA","GZMK","FGFBP2","CX3CR1"),
  Memory_Stemness = c("TCF7","CCR7","SELL","IL7R","LEF1","KLF2","BCL2","FOXO1"),
  Exhaustion      = c("PDCD1","LAG3","HAVCR2","TIGIT","TOX","TOX2","ENTPD1","CTLA4","BATF"),
  Activation      = c("CD69","CD44","TNFRSF9","IL2RA","ICOS","CD38"),
  Proliferation   = c("MKI67","TOP2A","PCNA","CCNB1","STMN1","UBE2C"),
  Tpex_StemLike   = c("TCF7","CXCR5","TOX","BCL6","SLAMF6","ID3"),
  Tex_Terminal    = c("HAVCR2","TIGIT","LAG3","CD160","ENTPD1","PRDM1","ZEB2")
)

FUNCTIONAL_STATE_MAP <- list(
  "Naive-like"   = c("Naive CD4+ T cells", "Naive CD8+ T cells"),
  "Memory-like"  = c("Memory T cells","Th1 cells","Th2 cells","Th17 cells","Tfh cells"),
  "Effector"     = c("Effector CD4+ T cells","Cytotoxic CD8+ T cells"),
  "Regulatory"   = c("Tregs"),
  "Proliferating"= c("Proliferating CD4+ T cells","Proliferating CD8+ T cells")
)
FUNCTIONAL_ORDER <- c("Naive-like","Memory-like","Effector","Regulatory","Proliferating")
STATE_PALETTE    <- c("Naive-like"="#4DBBD5","Memory-like"="#00A087","Effector"="#E64B35",
                      "Regulatory"="#F39B7F","Proliferating"="#7E6148")

get_car_status <- function(obj, sample_name) {
  meta <- obj@meta.data
  for (col in c("IS_CAR_ALLIN_scREP","IS_CAR","CAR")) {
    if (col %in% colnames(meta)) {
      vals    <- as.character(meta[[col]])
      car_pos <- grepl("^(YES|TRUE|yes|true|1)$", vals)
      cat(sprintf("  %s: '%s' | CAR+ = %d / %d (%.1f%%)\n",
                  sample_name, col, sum(car_pos), length(car_pos), 100*mean(car_pos)))
      return(ifelse(car_pos, "CAR+", "CAR-"))
    }
  }
  rep("CAR-", ncol(obj))
}

map_to_functional_state <- function(cell_types) {
  state <- rep("Other", length(cell_types))
  for (fs_name in names(FUNCTIONAL_STATE_MAP)) {
    state[cell_types %in% FUNCTIONAL_STATE_MAP[[fs_name]]] <- fs_name
  }
  state
}

# ── Carica dati ──────────────────────────────────────────────
section("Caricamento dati")
all_samples <- readRDS(rds_path)
cat("Campioni:", paste(names(all_samples), collapse = ", "), "\n")

# ── Sub-clustering ───────────────────────────────────────────
section("Sub-clustering T cells (con JoinLayers)")

all_T_objs <- list()
for (patient in names(PATIENT_MAP)) {
  for (sname in PATIENT_MAP[[patient]]$I) {
    if (!sname %in% names(all_samples)) next
    obj <- all_samples[[sname]]
    obj$car_status <- get_car_status(obj, sname)
    obj$sample     <- sname
    obj$patient    <- patient
    t_mask <- map_to_functional_state(as.character(obj@meta.data$cell_type)) != "Other"
    cat(sprintf("  %s: %d cellule T\n", sname, sum(t_mask)))
    if (sum(t_mask) < 30) next
    all_T_objs[[sname]] <- subset(obj, cells = which(t_mask))
  }
}

merged_T <- if (length(all_T_objs) == 1) all_T_objs[[1]] else
  merge(all_T_objs[[1]], y = all_T_objs[-1], add.cell.ids = names(all_T_objs))

cat(sprintf("  Totale cellule T: %d\n", ncol(merged_T)))

merged_T <- NormalizeData(merged_T, verbose = FALSE)
merged_T <- FindVariableFeatures(merged_T, nfeatures = 2000, verbose = FALSE)
var_genes <- VariableFeatures(merged_T)
car_genes_to_remove <- grep("^CAR|^GD2|^FMC63|SCFV|transgene",
                             var_genes, ignore.case = TRUE, value = TRUE)
if (length(car_genes_to_remove) > 0) VariableFeatures(merged_T) <-
  setdiff(var_genes, car_genes_to_remove)

merged_T <- ScaleData(merged_T, verbose = FALSE)
merged_T <- RunPCA(merged_T, npcs = 30, verbose = FALSE)
merged_T <- RunUMAP(merged_T, dims = 1:20, verbose = FALSE, min.dist = 0.3)
merged_T <- FindNeighbors(merged_T, dims = 1:20, verbose = FALSE)
merged_T <- FindClusters(merged_T, resolution = 0.4, verbose = FALSE)

cat(sprintf("  Sub-cluster: %s\n",
            paste(levels(merged_T$seurat_clusters), collapse = ", ")))

# CRITICO: JoinLayers prima di FindAllMarkers (Seurat v5)
merged_T <- JoinLayers(merged_T)

# ── Module scores sul merged_T ───────────────────────────────
for (sig_name in names(SIGNATURES)) {
  genes_ok <- intersect(SIGNATURES[[sig_name]], rownames(merged_T))
  if (length(genes_ok) < 2) next
  score_col <- paste0("score_", sig_name)
  merged_T <- AddModuleScore(merged_T, features = list(genes_ok),
                             name = score_col, seed = 42)
  old_col <- paste0(score_col, "1")
  merged_T@meta.data[[score_col]] <- merged_T@meta.data[[old_col]]
  merged_T@meta.data[[old_col]]   <- NULL
}

# ── UMAP plots ───────────────────────────────────────────────
p_umap_clust <- DimPlot(merged_T, group.by = "seurat_clusters",
                         label = TRUE, label.size = 4, pt.size = 0.8) +
  labs(title = "Sub-cluster T cells (I)", subtitle = "Risoluzione 0.4, agnostico CD4/CD8") +
  NoLegend()

p_umap_car <- DimPlot(merged_T,
                       cells.highlight = WhichCells(merged_T,
                                          expression = car_status == "CAR+"),
                       cols.highlight = "#E64B35", cols = "lightgrey", pt.size = 0.8) +
  labs(title = "CAR+ cells (rosso)") + theme(legend.position = "none")

p_umap_pat <- DimPlot(merged_T, group.by = "patient", pt.size = 0.8,
                       cols = c(Bo = "#E64B35", Ca = "#4DBBD5", Me = "#00A087")) +
  labs(title = "Per paziente")

p_umap_type <- DimPlot(merged_T, group.by = "cell_type", pt.size = 0.6) +
  labs(title = "Tipo cellulare annotato")

p_umap_combined <- (p_umap_clust | p_umap_car) / (p_umap_pat | p_umap_type)
ggsave(file.path(out_dir, "Q1b_ALL_Tcell_subcluster_UMAP.png"),
       p_umap_combined, width = 14, height = 12, dpi = 300)
cat("  Salvato: Q1b_ALL_Tcell_subcluster_UMAP.png\n")

# ── Feature plots ────────────────────────────────────────────
fp_plots <- lapply(names(SIGNATURES), function(sig_name) {
  sc <- paste0("score_", sig_name)
  if (!sc %in% colnames(merged_T@meta.data)) return(NULL)
  FeaturePlot(merged_T, features = sc, pt.size = 0.5,
              order = TRUE, min.cutoff = "q10") +
    scale_colour_gradientn(colours = c("lightgrey","#E64B35")) +
    labs(title = sig_name) +
    theme(legend.key.size = unit(0.4,"cm"),
          plot.title = element_text(size = 9, face = "bold"))
})
fp_plots <- Filter(Negate(is.null), fp_plots)
p_fp <- wrap_plots(fp_plots, ncol = 4) +
  plot_annotation(title = "Module scores sub-cluster T cells (I)",
                  theme = theme(plot.title = element_text(face = "bold")))
ggsave(file.path(out_dir, "Q1b_ALL_Tcell_subcluster_module_scores.png"),
       p_fp, width = 16, height = 10, dpi = 300)
cat("  Salvato: Q1b_ALL_Tcell_subcluster_module_scores.png\n")

# ── FindAllMarkers ───────────────────────────────────────────
section("Marker per sub-cluster (FindAllMarkers)")

Idents(merged_T) <- "seurat_clusters"
cluster_markers <- FindAllMarkers(
  merged_T,
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.25,
  test.use        = "wilcox",
  verbose         = FALSE
)

if (nrow(cluster_markers) == 0 || !"p_val_adj" %in% colnames(cluster_markers)) {
  cat("  [WARN] Nessun marker trovato.\n")
  top_markers  <- data.frame()
  top5_genes   <- character(0)
} else {
  top_markers <- cluster_markers %>%
    filter(p_val_adj < 0.05) %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = 10) %>%
    ungroup()
  cat(sprintf("  Marker significativi: %d\n",
              sum(cluster_markers$p_val_adj < 0.05)))

  top5_genes <- top_markers %>%
    group_by(cluster) %>%
    slice_max(avg_log2FC, n = 5) %>%
    pull(gene) %>%
    unique()
}

if (length(top5_genes) > 0) {
  p_dot <- DotPlot(merged_T, features = top5_genes,
                   group.by = "seurat_clusters") +
    RotatedAxis() +
    scale_color_gradientn(colours = c("lightgrey","#E64B35")) +
    labs(title = "Top marker per sub-cluster T cells",
         subtitle = "Prodotto infusione I — agnostico CD4/CD8") +
    theme(axis.text.x = element_text(size = 8))
  ggsave(file.path(out_dir, "Q1b_ALL_subcluster_dotplot_markers.png"),
         p_dot, width = max(12, length(top5_genes) * 0.5), height = 6, dpi = 300)
  cat("  Salvato: Q1b_ALL_subcluster_dotplot_markers.png\n")
}

# ── Proporzione CAR+ per sub-cluster ────────────────────────
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
  geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", 100*prop, n)),
            vjust = -0.3, size = 3.5) +
  scale_y_continuous(labels = percent_format(),
                     expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Proporzione CAR+ per sub-cluster",
       subtitle = "Prodotto infusione I (cellule T, agnostico CD4/CD8)",
       x = "Sub-cluster", y = "% CAR+") +
  theme_classic(base_size = 12)

ggsave(file.path(out_dir, "Q1b_ALL_CARpos_proportion_per_subcluster.png"),
       p_clust_car, width = 8, height = 5, dpi = 300)
cat("  Salvato: Q1b_ALL_CARpos_proportion_per_subcluster.png\n")

# ── Excel markers ────────────────────────────────────────────
if (nrow(cluster_markers) > 0) {
  wb <- createWorkbook()
  addWorksheet(wb, "SubclusterMarkers")
  writeData(wb, "SubclusterMarkers", cluster_markers)
  saveWorkbook(wb, file.path(out_dir, "Q1b_subcluster_markers.xlsx"),
               overwrite = TRUE)
  cat("  Salvato: Q1b_subcluster_markers.xlsx\n")
}

# ── Sezione 4: I vs AB CAR+ ──────────────────────────────────
section("Sezione 4: I vs AB CAR+ — stati funzionali")

extract_functional_props <- function(obj, sample_name, timepoint) {
  meta    <- obj@meta.data
  car_vec <- get_car_status(obj, sample_name)
  fs_vec  <- map_to_functional_state(as.character(meta$cell_type))
  df <- data.frame(functional_state = fs_vec, car_status = car_vec,
                   sample = sample_name, timepoint = timepoint,
                   stringsAsFactors = FALSE)
  df_t <- df[df$functional_state != "Other", ]
  if (nrow(df_t) == 0) return(NULL)
  df_t %>%
    group_by(car_status, functional_state) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(car_status) %>%
    mutate(total = sum(n), prop = n / total) %>%
    ungroup() %>%
    mutate(sample = sample_name, timepoint = timepoint,
           functional_state = factor(functional_state, levels = FUNCTIONAL_ORDER))
}

all_props <- list()
for (patient in names(PATIENT_MAP)) {
  pat_props <- list()
  for (tp in c("I","AB")) {
    for (sname in PATIENT_MAP[[patient]][[tp]]) {
      if (!sname %in% names(all_samples)) next
      p <- extract_functional_props(all_samples[[sname]], sname, tp)
      if (!is.null(p)) pat_props[[sname]] <- p
    }
  }
  if (length(pat_props) > 0) all_props[[patient]] <- bind_rows(pat_props)
}

combined_df <- bind_rows(all_props, .id = "patient") %>%
  filter(car_status == "CAR+") %>%
  group_by(patient, timepoint, functional_state) %>%
  summarise(prop = mean(prop), .groups = "drop") %>%
  mutate(functional_state = factor(functional_state, levels = FUNCTIONAL_ORDER),
         timepoint = factor(timepoint, levels = c("I","AB")))

p_all_patients <- ggplot(combined_df,
                          aes(x = timepoint, y = prop, fill = functional_state)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.3) +
  facet_wrap(~ patient, ncol = 3) +
  scale_fill_manual(values = STATE_PALETTE, drop = FALSE, name = "Stato funzionale") +
  scale_y_continuous(labels = percent_format(), expand = c(0,0)) +
  labs(title    = "CAR+ cells — Stati funzionali: I vs AB",
       subtitle = "Agnostico CD4/CD8 | Proporzioni dentro CAR+",
       x = "Timepoint", y = "Proporzione") +
  theme_classic(base_size = 12) +
  theme(strip.background = element_rect(fill = "#F0F0F0"),
        strip.text = element_text(face = "bold", size = 12))

ggsave(file.path(out_dir, "Q1b_ALL_CARpos_I_vs_AB_functional_states.png"),
       p_all_patients, width = 12, height = 5, dpi = 300)
cat("Salvato: Q1b_ALL_CARpos_I_vs_AB_functional_states.png\n")

# ── Excel summary ────────────────────────────────────────────
wb2 <- createWorkbook()
addWorksheet(wb2, "Functional_Props")
writeData(wb2, "Functional_Props", bind_rows(all_props, .id = "patient"))
saveWorkbook(wb2, file.path(out_dir, "Q1b_functional_summary.xlsx"),
             overwrite = TRUE)
cat("Salvato: Q1b_functional_summary.xlsx\n")

section("COMPLETATO")
cat("File prodotti:\n")
for (f in list.files(out_dir, full.names = FALSE))
  cat(sprintf("  %s\n", f))
