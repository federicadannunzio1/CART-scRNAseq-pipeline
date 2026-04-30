library(Seurat)
library(dplyr)
library(ggplot2)
library(openxlsx)

# ==========================================
# 1. PARAMETRI GENERALI
# ==========================================

min_pct <- 0.25
logfc_threshold <- 0.25
top_n_markers <- 30

canonical_markers <- c(
  # T cells
  "CD3D","CD3E","CD4","CD8A","CD8B",
  # Naive / memory
  "CCR7","IL7R","TCF7","GZMK",
  # Cytotoxic
  "GZMB","PRF1","NKG7","IFNG",
  # Proliferation
  "MKI67","TOP2A","STMN1",
  # Treg
  "FOXP3","IL2RA","CTLA4","TIGIT",
  # APC / DC
  "HLA-DRA","HLA-DPB1","HLA-DQA1",
  "CD1C","CLEC9A","LILRA4","FCER1A","ITGAX","BATF3","XCR1",
  # Monocytes
  "LYZ","S100A8","S100A9","FCGR3A"
)

# ==========================================
# 2. FUNZIONE PER ANALISI COMPLETA
# ==========================================

analyze_sample <- function(seurat_obj, sample_name, output_dir) {
  
  Idents(seurat_obj) <- "seurat_clusters"
  
  cat(paste0("\nAnalisi sample: ", sample_name, "\n"))
  
  # ----------------------------------------
  # A. MARKER DIFFERENZIALI
  # ----------------------------------------
  markers <- FindAllMarkers(
    seurat_obj,
    only.pos = TRUE,
    min.pct = min_pct,
    logfc.threshold = logfc_threshold
  )
  
  markers <- markers %>%
    arrange(cluster, desc(avg_log2FC))
  
  # Top marker per cluster
  top_markers <- markers %>%
    group_by(cluster) %>%
    slice_max(n = top_n_markers, order_by = avg_log2FC)
  
  # ----------------------------------------
  # B. PERCENTUALE CELLULE PER CLUSTER
  # ----------------------------------------
  cluster_freq <- as.data.frame(table(Idents(seurat_obj)))
  colnames(cluster_freq) <- c("cluster", "n_cells")
  cluster_freq$percentage <- round(100 * cluster_freq$n_cells / sum(cluster_freq$n_cells), 2)
  
  # ----------------------------------------
  # C. DOTPLOT MARKER CANONICI
  # ----------------------------------------
  p_dot <- DotPlot(
    seurat_obj,
    features = canonical_markers
  ) + RotatedAxis() + ggtitle(paste0(sample_name, " - Canonical markers"))
  
  ggsave(
    filename = paste0(output_dir, "/", sample_name, "_DotPlot_canonical.png"),
    plot = p_dot,
    width = 14,
    height = 8
  )
  
  # ----------------------------------------
  # D. HEATMAP TOP MARKER
  # ----------------------------------------
  p_heat <- DoHeatmap(
    seurat_obj,
    features = unique(top_markers$gene),
    size = 3
  ) + NoLegend()
  
  ggsave(
    filename = paste0(output_dir, "/", sample_name, "_Heatmap_topMarkers.png"),
    plot = p_heat,
    width = 12,
    height = 10
  )
  
  # ----------------------------------------
  # E. ESPORTAZIONE EXCEL
  # ----------------------------------------
  wb <- createWorkbook()
  
  addWorksheet(wb, "All_markers")
  writeData(wb, "All_markers", markers)
  
  addWorksheet(wb, "Top_markers")
  writeData(wb, "Top_markers", top_markers)
  
  addWorksheet(wb, "Cluster_frequency")
  writeData(wb, "Cluster_frequency", cluster_freq)
  
  saveWorkbook(
    wb,
    file = paste0(output_dir, "/", sample_name, "_Annotation_Output.xlsx"),
    overwrite = TRUE
  )
  
  return(list(
    markers = markers,
    top_markers = top_markers,
    cluster_freq = cluster_freq
  ))
}

# ==========================================
# 3. ESECUZIONE PER TUTTI I CAMPIONI
# ==========================================

output_dir <- "Manual_annotation_output"
dir.create(output_dir, showWarnings = FALSE)

results_Bo <- analyze_sample(Bo_I, "Bo_I", output_dir)
results_Ca <- analyze_sample(Ca_I, "Ca_I", output_dir)
results_Me <- analyze_sample(Me_I, "Me_I", output_dir)

cat("\nOutput completo generato.\n")
