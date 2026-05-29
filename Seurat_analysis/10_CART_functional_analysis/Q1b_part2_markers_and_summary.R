# ============================================================
#  Q1b — Parte 2: Sub-clustering (con JoinLayers) + I vs AB
#
#  Usa file RDS separati per gestire la memoria (8 GB RAM):
#    BLOCCO A: all_I_samples_annotated.rds (244 MB)
#              → sub-clustering, FindAllMarkers, CAR+ per cluster
#    BLOCCO B: all_AB_samples_annotated.rds (807 MB)
#              → proporzioni stati funzionali I vs AB
#
#  Prerequisiti: output della Sezione 1 e 2 già salvati
#                da Q1b_functional_states_agnostic.R
# ============================================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(scales)
library(openxlsx)

out_dir <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/CART/Code/Seurat_analysis/10_CART_functional_analysis/Q1b_functional_states/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

rds_I  <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/CART/Code/Seurat_analysis/2_annotation/all_I_samples_annotated.rds"
rds_AB <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/CART/Code/Seurat_analysis/2_annotation/all_AB_samples_annotated.rds"

section <- function(title)
  cat(paste0("\n", strrep("=", 65), "\n  ", title, "\n",
             strrep("=", 65), "\n"))

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

SIGNATURES <- list(
  Effector        = c("GZMB","PRF1","NKG7","GNLY","GZMA","GZMK","FGFBP2","CX3CR1"),
  Memory_Stemness = c("TCF7","CCR7","SELL","IL7R","LEF1","KLF2","BCL2","FOXO1"),
  Exhaustion      = c("PDCD1","LAG3","HAVCR2","TIGIT","TOX","TOX2","ENTPD1","CTLA4","BATF"),
  Activation      = c("CD69","CD44","TNFRSF9","IL2RA","ICOS","CD38"),
  Proliferation   = c("MKI67","TOP2A","PCNA","CCNB1","STMN1","UBE2C"),
  Tpex_StemLike   = c("TCF7","CXCR5","TOX","BCL6","SLAMF6","ID3"),
  Tex_Terminal    = c("HAVCR2","TIGIT","LAG3","CD160","ENTPD1","PRDM1","ZEB2")
)

PATIENT_MAP_I  <- list(Bo="Bo_bone_I", Ca="Ca_bone_I", Me="Me_bone_I")
PATIENT_MAP_AB <- list(
  Bo = c("Bo_blood_AB","Bo_bone_AB"),
  Ca = c("Ca_blood_AB","Ca_bone_AB"),
  Me = c("Me_bone_AB")
)

get_car_status <- function(obj, sample_name) {
  meta <- obj@meta.data
  for (col in c("IS_CAR_ALLIN_scREP","IS_CAR","CAR")) {
    if (col %in% colnames(meta)) {
      vals <- as.character(meta[[col]])
      car_pos <- grepl("^(YES|TRUE|yes|true|1)$", vals)
      cat(sprintf("  %s: CAR+ = %d / %d (%.1f%%)\n",
                  sample_name, sum(car_pos), length(car_pos), 100*mean(car_pos)))
      return(ifelse(car_pos, "CAR+", "CAR-"))
    }
  }
  rep("CAR-", ncol(obj))
}

map_fs <- function(cell_types) {
  state <- rep("Other", length(cell_types))
  for (nm in names(FUNCTIONAL_STATE_MAP))
    state[cell_types %in% FUNCTIONAL_STATE_MAP[[nm]]] <- nm
  state
}

# ============================================================
# BLOCCO A: I samples → sub-clustering + FindAllMarkers
# ============================================================
section("BLOCCO A: Caricamento campioni I (244 MB)")

I_samples <- readRDS(rds_I)
cat("Campioni I:", paste(names(I_samples), collapse=", "), "\n")

# Filtra T cells, aggiungi metadati
all_T_objs <- list()
for (sname in names(I_samples)) {
  obj <- I_samples[[sname]]
  patient <- sub("_bone_I$","", sname)
  obj$car_status <- get_car_status(obj, sname)
  obj$sample     <- sname
  obj$patient    <- patient
  t_mask <- map_fs(as.character(obj@meta.data$cell_type)) != "Other"
  cat(sprintf("  %s: %d T cells\n", sname, sum(t_mask)))
  if (sum(t_mask) < 30) next
  all_T_objs[[sname]] <- subset(obj, cells = which(t_mask))
}

# Libera memoria campioni I originali
rm(I_samples); invisible(gc())

section("Sub-clustering T cells")

merged_T <- if (length(all_T_objs) == 1) all_T_objs[[1]] else
  merge(all_T_objs[[1]], y = all_T_objs[-1], add.cell.ids = names(all_T_objs))
rm(all_T_objs); invisible(gc())

cat(sprintf("  Totale T cells: %d\n", ncol(merged_T)))

merged_T <- NormalizeData(merged_T, verbose=FALSE)
merged_T <- FindVariableFeatures(merged_T, nfeatures=2000, verbose=FALSE)
var_genes <- VariableFeatures(merged_T)
car_rm <- grep("^CAR|^GD2|^FMC63|SCFV|transgene", var_genes, ignore.case=TRUE, value=TRUE)
if (length(car_rm)>0) VariableFeatures(merged_T) <- setdiff(var_genes, car_rm)
merged_T <- ScaleData(merged_T, verbose=FALSE)
merged_T <- RunPCA(merged_T, npcs=30, verbose=FALSE)
merged_T <- RunUMAP(merged_T, dims=1:20, verbose=FALSE, min.dist=0.3)
merged_T <- FindNeighbors(merged_T, dims=1:20, verbose=FALSE)
merged_T <- FindClusters(merged_T, resolution=0.4, verbose=FALSE)
cat(sprintf("  Sub-cluster: %s\n", paste(levels(merged_T$seurat_clusters), collapse=", ")))

# Module scores
for (sig_name in names(SIGNATURES)) {
  genes_ok <- intersect(SIGNATURES[[sig_name]], rownames(merged_T))
  if (length(genes_ok) < 2) next
  sc <- paste0("score_", sig_name)
  merged_T <- AddModuleScore(merged_T, features=list(genes_ok), name=sc, seed=42)
  merged_T@meta.data[[sc]] <- merged_T@meta.data[[paste0(sc,"1")]]
  merged_T@meta.data[[paste0(sc,"1")]] <- NULL
}

# ── UMAP plots ───────────────────────────────────────────────
p1 <- DimPlot(merged_T, group.by="seurat_clusters", label=TRUE, label.size=4, pt.size=0.8) +
  labs(title="Sub-cluster T cells (I)", subtitle="Risoluzione 0.4, agnostico CD4/CD8") + NoLegend()

p2 <- DimPlot(merged_T,
              cells.highlight=WhichCells(merged_T, expression=car_status=="CAR+"),
              cols.highlight="#E64B35", cols="lightgrey", pt.size=0.8) +
  labs(title="CAR+ cells (rosso)") + theme(legend.position="none")

p3 <- DimPlot(merged_T, group.by="patient", pt.size=0.8,
              cols=c(Bo="#E64B35",Ca="#4DBBD5",Me="#00A087")) +
  labs(title="Per paziente")

p4 <- DimPlot(merged_T, group.by="cell_type", pt.size=0.6) +
  labs(title="Tipo cellulare annotato")

ggsave(file.path(out_dir,"Q1b_ALL_Tcell_subcluster_UMAP.png"),
       (p1|p2)/(p3|p4), width=14, height=12, dpi=300)
cat("  Salvato: Q1b_ALL_Tcell_subcluster_UMAP.png\n")

# ── Feature plots module scores ───────────────────────────────
fp_plots <- lapply(names(SIGNATURES), function(sig_name) {
  sc <- paste0("score_", sig_name)
  if (!sc %in% colnames(merged_T@meta.data)) return(NULL)
  FeaturePlot(merged_T, features=sc, pt.size=0.5, order=TRUE, min.cutoff="q10") +
    scale_colour_gradientn(colours=c("lightgrey","#E64B35")) +
    labs(title=sig_name) +
    theme(legend.key.size=unit(0.4,"cm"),
          plot.title=element_text(size=9,face="bold"))
})
fp_plots <- Filter(Negate(is.null), fp_plots)
ggsave(file.path(out_dir,"Q1b_ALL_Tcell_subcluster_module_scores.png"),
       wrap_plots(fp_plots,ncol=4) +
         plot_annotation(title="Module scores sub-cluster T cells (I)",
                         theme=theme(plot.title=element_text(face="bold"))),
       width=16, height=10, dpi=300)
cat("  Salvato: Q1b_ALL_Tcell_subcluster_module_scores.png\n")

# ── JoinLayers + FindAllMarkers ───────────────────────────────
section("FindAllMarkers (con JoinLayers)")
merged_T <- JoinLayers(merged_T)
Idents(merged_T) <- "seurat_clusters"

cluster_markers <- FindAllMarkers(
  merged_T, only.pos=TRUE, min.pct=0.25,
  logfc.threshold=0.25, test.use="wilcox", verbose=FALSE
)

if (nrow(cluster_markers)==0 || !"p_val_adj" %in% colnames(cluster_markers)) {
  cat("  [WARN] Nessun marker trovato.\n")
  top_markers <- data.frame(); top5_genes <- character(0)
} else {
  top_markers <- cluster_markers %>%
    filter(p_val_adj < 0.05) %>%
    group_by(cluster) %>%
    slice_max(avg_log2FC, n=10) %>% ungroup()
  cat(sprintf("  Marker significativi: %d\n", sum(cluster_markers$p_val_adj < 0.05)))

  top5_genes <- top_markers %>%
    group_by(cluster) %>% slice_max(avg_log2FC, n=5) %>%
    pull(gene) %>% unique()
}

if (length(top5_genes) > 0) {
  p_dot <- DotPlot(merged_T, features=top5_genes, group.by="seurat_clusters") +
    RotatedAxis() +
    scale_color_gradientn(colours=c("lightgrey","#E64B35")) +
    labs(title="Top marker per sub-cluster T cells",
         subtitle="Prodotto infusione I — agnostico CD4/CD8") +
    theme(axis.text.x=element_text(size=8))
  ggsave(file.path(out_dir,"Q1b_ALL_subcluster_dotplot_markers.png"),
         p_dot, width=max(12,length(top5_genes)*0.5), height=6, dpi=300)
  cat("  Salvato: Q1b_ALL_subcluster_dotplot_markers.png\n")
}

# ── Proporzione CAR+ per cluster ─────────────────────────────
clust_car_df <- merged_T@meta.data %>%
  group_by(seurat_clusters, car_status) %>%
  summarise(n=n(), .groups="drop") %>%
  group_by(seurat_clusters) %>%
  mutate(prop=n/sum(n)) %>% ungroup()

p_clust_car <- ggplot(
  clust_car_df %>% filter(car_status=="CAR+"),
  aes(x=seurat_clusters, y=prop, fill=seurat_clusters)
) +
  geom_col(show.legend=FALSE) +
  geom_text(aes(label=sprintf("%.1f%%\n(n=%d)", 100*prop, n)),
            vjust=-0.3, size=3.5) +
  scale_y_continuous(labels=percent_format(),
                     expand=expansion(mult=c(0,0.15))) +
  labs(title="Proporzione CAR+ per sub-cluster",
       subtitle="Prodotto infusione I — agnostico CD4/CD8",
       x="Sub-cluster", y="% CAR+") +
  theme_classic(base_size=12)

ggsave(file.path(out_dir,"Q1b_ALL_CARpos_proportion_per_subcluster.png"),
       p_clust_car, width=8, height=5, dpi=300)
cat("  Salvato: Q1b_ALL_CARpos_proportion_per_subcluster.png\n")

# Salva metadata del clustering per uso futuro
merged_T_meta <- merged_T@meta.data
if (nrow(cluster_markers) > 0) {
  wb <- createWorkbook()
  addWorksheet(wb,"SubclusterMarkers"); writeData(wb,"SubclusterMarkers",cluster_markers)
  addWorksheet(wb,"ClusterMeta"); writeData(wb,"ClusterMeta",merged_T_meta)
  saveWorkbook(wb, file.path(out_dir,"Q1b_subcluster_markers.xlsx"), overwrite=TRUE)
  cat("  Salvato: Q1b_subcluster_markers.xlsx\n")
}

# Libera merged_T prima di caricare AB
rm(merged_T, merged_T_meta, cluster_markers, top_markers, all_T_objs)
invisible(gc()); invisible(gc())

# ============================================================
# BLOCCO B: AB samples → proporzioni I vs AB (CAR+)
# ============================================================
section("BLOCCO B: Caricamento campioni AB (807 MB)")

# Calcola proporzioni per campioni I (ricarica solo il file piccolo)
cat("Ricaricamento campioni I...\n")
I_samples <- readRDS(rds_I)

extract_props <- function(obj, sname, tp) {
  meta <- obj@meta.data
  for (col in c("IS_CAR_ALLIN_scREP","IS_CAR","CAR")) {
    if (col %in% colnames(meta)) {
      vals <- as.character(meta[[col]])
      car_vec <- ifelse(grepl("^(YES|TRUE|yes|true|1)$",vals),"CAR+","CAR-")
      fs_vec  <- map_fs(as.character(meta$cell_type))
      df <- data.frame(functional_state=fs_vec, car_status=car_vec,
                       sample=sname, timepoint=tp, stringsAsFactors=FALSE)
      df_t <- df[df$functional_state!="Other",]
      if (nrow(df_t)==0) return(NULL)
      return(df_t %>%
        group_by(car_status, functional_state) %>%
        summarise(n=n(), .groups="drop") %>%
        group_by(car_status) %>%
        mutate(total=sum(n), prop=n/total) %>% ungroup() %>%
        mutate(sample=sname, timepoint=tp,
               functional_state=factor(functional_state,levels=FUNCTIONAL_ORDER)))
    }
  }
  NULL
}

props_I <- list()
for (sname in names(I_samples)) {
  patient <- sub("_bone_I$","",sname)
  p <- extract_props(I_samples[[sname]], sname, "I")
  if (!is.null(p)) props_I[[sname]] <- mutate(p, patient=patient)
}
rm(I_samples); invisible(gc())

cat("Caricamento campioni AB...\n")
AB_samples <- readRDS(rds_AB)
cat("Campioni AB:", paste(names(AB_samples), collapse=", "), "\n")

props_AB <- list()
for (patient in names(PATIENT_MAP_AB)) {
  for (sname in PATIENT_MAP_AB[[patient]]) {
    if (!sname %in% names(AB_samples)) next
    p <- extract_props(AB_samples[[sname]], sname, "AB")
    if (!is.null(p)) props_AB[[sname]] <- mutate(p, patient=patient)
  }
}
rm(AB_samples); invisible(gc())

all_props <- bind_rows(c(props_I, props_AB))

# Barplot per paziente (I vs AB, solo CAR+)
combined_df <- all_props %>%
  filter(car_status=="CAR+") %>%
  group_by(patient, timepoint, functional_state) %>%
  summarise(prop=mean(prop), .groups="drop") %>%
  mutate(functional_state=factor(functional_state, levels=FUNCTIONAL_ORDER),
         timepoint=factor(timepoint, levels=c("I","AB")))

p_all <- ggplot(combined_df,
                aes(x=timepoint, y=prop, fill=functional_state)) +
  geom_col(width=0.7, color="white", linewidth=0.3) +
  facet_wrap(~patient, ncol=3) +
  scale_fill_manual(values=STATE_PALETTE, drop=FALSE, name="Stato funzionale") +
  scale_y_continuous(labels=percent_format(), expand=c(0,0)) +
  labs(title="CAR+ cells — Stati funzionali: I vs AB",
       subtitle="Agnostico CD4/CD8 | Proporzioni dentro CAR+",
       x="Timepoint", y="Proporzione") +
  theme_classic(base_size=12) +
  theme(strip.background=element_rect(fill="#F0F0F0"),
        strip.text=element_text(face="bold", size=12))

ggsave(file.path(out_dir,"Q1b_ALL_CARpos_I_vs_AB_functional_states.png"),
       p_all, width=12, height=5, dpi=300)
cat("Salvato: Q1b_ALL_CARpos_I_vs_AB_functional_states.png\n")

# Heatmap I (CAR+)
heat_I <- all_props %>%
  filter(car_status=="CAR+", timepoint=="I") %>%
  group_by(patient, functional_state) %>%
  summarise(prop=mean(prop), .groups="drop") %>%
  mutate(functional_state=factor(functional_state,levels=FUNCTIONAL_ORDER))

p_hI <- ggplot(heat_I, aes(x=patient,y=functional_state,fill=prop)) +
  geom_tile(color="white",linewidth=0.5) +
  geom_text(aes(label=sprintf("%.1f%%",100*prop)),size=3.5) +
  scale_fill_gradient(low="white",high="#E64B35",labels=percent_format(),name="Prop") +
  labs(title="CAR+ in I — stati funzionali",x="Paziente",y=NULL) +
  theme_classic(base_size=12)
ggsave(file.path(out_dir,"Q1b_ALL_heatmap_CARpos_I.png"),
       p_hI, width=6, height=4, dpi=300)
cat("Salvato: Q1b_ALL_heatmap_CARpos_I.png\n")

# Heatmap AB (CAR+)
heat_AB <- all_props %>%
  filter(car_status=="CAR+", timepoint=="AB") %>%
  group_by(patient, functional_state) %>%
  summarise(prop=mean(prop), .groups="drop") %>%
  mutate(functional_state=factor(functional_state,levels=FUNCTIONAL_ORDER))

if (nrow(heat_AB) > 0) {
  p_hAB <- ggplot(heat_AB, aes(x=patient,y=functional_state,fill=prop)) +
    geom_tile(color="white",linewidth=0.5) +
    geom_text(aes(label=sprintf("%.1f%%",100*prop)),size=3.5) +
    scale_fill_gradient(low="white",high="#3C5488",labels=percent_format(),name="Prop") +
    labs(title="CAR+ in AB — stati funzionali",x="Paziente",y=NULL) +
    theme_classic(base_size=12)
  ggsave(file.path(out_dir,"Q1b_ALL_heatmap_CARpos_AB.png"),
         p_hAB, width=6, height=4, dpi=300)
  cat("Salvato: Q1b_ALL_heatmap_CARpos_AB.png\n")
}

# Excel summary
wb2 <- createWorkbook()
addWorksheet(wb2,"Functional_Props"); writeData(wb2,"Functional_Props",all_props)
saveWorkbook(wb2, file.path(out_dir,"Q1b_functional_summary.xlsx"), overwrite=TRUE)
cat("Salvato: Q1b_functional_summary.xlsx\n")

section("COMPLETATO")
cat("Output in:", out_dir, "\n")
cat("File prodotti:\n")
for (f in list.files(out_dir, full.names=FALSE))
  cat(sprintf("  %s\n", f))
