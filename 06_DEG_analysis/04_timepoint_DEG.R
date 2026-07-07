# Timepoint Comparison DEG Analysis
# Converted from P4_timepoint_comparison.Rmd

# ======================================================================
# LOAD LIBRARIES
# ======================================================================
suppressMessages({
  library(Seurat)
  library(DESeq2)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(DT)
  library(knitr)
  library(ggrepel)
  library(Matrix)
})

suppressMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

if (!requireNamespace("writexl", quietly = TRUE))
  install.packages("writexl", repos = "https://cran.rstudio.com/")
suppressMessages(library(writexl))

# ======================================================================
# DEFINE PATHS
# ======================================================================
RDS_I  <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/CART/Code/Seurat_analysis/2_annotation/all_I_samples_annotated.rds"
RDS_AB <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/CART/Code/Seurat_analysis/2_annotation/all_AB_samples_annotated.rds"

OUT_DIR <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/CART/Code/Seurat_analysis/11_integrated_pipeline/P4_results"

for (cmp in c("I_vs_Blood", "I_vs_Bone")) {
  for (sub in c("figures", "tables"))
    dir.create(file.path(OUT_DIR, cmp, sub), recursive = TRUE, showWarnings = FALSE)
}

MIN_CAR_CELLS  <- 10
PADJ_CUTOFF    <- 0.05
LFC_CUTOFF     <- 0.5     # log2FC minimo per filtraggio biologico (non usato nel test DESeq2)
ENRICH_PVAL    <- 0.05
MIN_GENE_SET   <- 10
MAX_GENE_SET   <- 500

# Geni da escludere dai risultati
NOISE_PATTERNS <- c("^MT-", "^RPS", "^RPL", "^TRAV", "^TRB")
NOISE_EXACT    <- c("MALAT1", "NEAT1")

filter_noise <- function(genes) {
  genes[!grepl(paste(NOISE_PATTERNS, collapse="|"), genes) & !genes %in% NOISE_EXACT]
}

PATIENT_COLORS <- c(Bo = "#E64B35", Ca = "#4DBBD5", Me = "#00A087")

# ======================================================================
# HELPERS
# ======================================================================
get_car_status <- function(meta) {
  for (col in c("IS_CAR_ALLIN_scREP","IS_CAR","CAR")) {
    if (col %in% colnames(meta)) {
      vals <- as.character(meta[[col]])
      return(ifelse(grepl("^(YES|TRUE|yes|true|1)$", vals), "CAR+", "CAR-"))
    }
  }
  rep("CAR-", nrow(meta))
}

# Pseudobulk: somma raw counts delle CAR+ cells
make_pseudobulk <- function(obj, sample_id) {
  cells_car <- rownames(obj@meta.data)[obj@meta.data$car_status == "CAR+"]
  if (length(cells_car) < MIN_CAR_CELLS) return(NULL)
  mat <- GetAssayData(obj, assay = "RNA", layer = "counts")[, cells_car, drop = FALSE]
  pb  <- Matrix::rowSums(mat)
  cat("  ", sample_id, ": CAR+ cells =", length(cells_car),
      "| geni con counts > 0 =", sum(pb > 0), "\n")
  list(counts = pb, n_cells = length(cells_car), sample_id = sample_id)
}

run_go <- function(genes, label = "") {
  if (length(genes) < 5) return(NULL)
  eg <- tryCatch(bitr(genes, fromType="SYMBOL", toType="ENTREZID",
                      OrgDb=org.Hs.eg.db, drop=TRUE), error=function(e) NULL)
  if (is.null(eg) || nrow(eg) < 5) return(NULL)
  res <- enrichGO(gene=eg$ENTREZID, OrgDb=org.Hs.eg.db, ont="BP",
                  pAdjustMethod="BH", pvalueCutoff=ENRICH_PVAL, qvalueCutoff=0.2,
                  readable=TRUE, minGSSize=MIN_GENE_SET, maxGSSize=MAX_GENE_SET)
  if (is.null(res) || nrow(as.data.frame(res))==0) return(NULL)
  res
}

run_kegg <- function(genes, label = "") {
  if (length(genes) < 5) return(NULL)
  eg <- tryCatch(bitr(genes, fromType="SYMBOL", toType="ENTREZID",
                      OrgDb=org.Hs.eg.db, drop=TRUE), error=function(e) NULL)
  if (is.null(eg) || nrow(eg) < 5) return(NULL)
  res <- tryCatch(
    enrichKEGG(gene=eg$ENTREZID, organism="hsa",
               pAdjustMethod="BH", pvalueCutoff=ENRICH_PVAL, qvalueCutoff=0.2,
               minGSSize=MIN_GENE_SET, maxGSSize=MAX_GENE_SET),
    error = function(e) { cat(label, ": KEGG non raggiungibile (", conditionMessage(e), ")\n"); NULL }
  )
  if (is.null(res) || nrow(as.data.frame(res))==0) return(NULL)
  tryCatch(setReadable(res, OrgDb=org.Hs.eg.db, keyType="ENTREZID"), error=function(e) res)
}

run_deseq2 <- function(pb_list, coldata, design_formula, contrast,
                        comparison_label) {
  # Costruisci matrice counts: geni x campioni
  all_genes <- Reduce(intersect, lapply(pb_list, function(x) names(x$counts)))
  cat("Geni in comune tra campioni:", length(all_genes), "\n")

  count_mat <- do.call(cbind, lapply(pb_list, function(x) x$counts[all_genes]))
  colnames(count_mat) <- sapply(pb_list, function(x) x$sample_id)

  # Rimuovi geni a bassa espressione: espressi in almeno 2 campioni con counts >= 5
  keep <- rowSums(count_mat >= 5) >= 2
  count_mat <- count_mat[keep, ]
  cat("Geni dopo filtro di espressione minima:", nrow(count_mat), "\n")

  count_mat <- round(count_mat)
  storage.mode(count_mat) <- "integer"

  dds <- DESeqDataSetFromMatrix(
    countData = count_mat,
    colData   = coldata,
    design    = design_formula
  )

  dds <- tryCatch(
    DESeq(dds, fitType = "local", quiet = TRUE),
    error = function(e) {
      cat("DESeq2 errore con fitType=local, provo mean:\n", conditionMessage(e), "\n")
      DESeq(dds, fitType = "mean", quiet = TRUE)
    }
  )

  res <- results(dds, contrast = contrast, alpha = PADJ_CUTOFF)
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  res_df <- res_df[!is.na(res_df$padj), ]

  # Filtro significatività
  sig <- res_df[res_df$padj < PADJ_CUTOFF, ]
  cat("Geni DE significativi (padj <", PADJ_CUTOFF, "):", nrow(sig), "\n")

  # Filtro geni rumore / TCR
  sig <- sig[sig$gene %in% filter_noise(sig$gene), ]
  cat("Dopo filtro rumore/TCR:", nrow(sig), "\n")

  sig[order(sig$log2FoldChange, decreasing = TRUE), ]
}

# ======================================================================
# LOAD DATA
# ======================================================================
cat("Caricamento campioni I...\n")
I_samples_all <- readRDS(RDS_I)
for (sname in names(I_samples_all)) {
  meta <- I_samples_all[[sname]]@meta.data
  I_samples_all[[sname]]$car_status <- get_car_status(meta)
  I_samples_all[[sname]]$patient    <- sub("_I$|_bone_I$", "", sname)
}

cat("Caricamento campioni AB...\n")
AB_raw <- readRDS(RDS_AB)

AB_blood <- AB_raw[grep("blood_AB", names(AB_raw), value=TRUE)]
for (sname in names(AB_blood)) {
  meta <- AB_blood[[sname]]@meta.data
  AB_blood[[sname]]$car_status <- get_car_status(meta)
  AB_blood[[sname]]$patient    <- sub("_blood_AB$", "", sname)
}

AB_bone <- AB_raw[grep("bone_AB", names(AB_raw), value=TRUE)]
rm(AB_raw); invisible(gc())
for (sname in names(AB_bone)) {
  meta <- AB_bone[[sname]]@meta.data
  AB_bone[[sname]]$car_status <- get_car_status(meta)
  AB_bone[[sname]]$patient    <- sub("_bone_AB$", "", sname)
}

get_pt_obj <- function(seurat_list, patient) {
  m <- Filter(Negate(is.null), lapply(seurat_list, function(obj) {
    if (patient %in% unique(obj@meta.data$patient)) obj else NULL
  }))
  if (length(m) == 0) return(NULL)
  m[[1]]
}

# ======================================================================
# CAR COUNTS
# ======================================================================
count_df <- bind_rows(lapply(c("Bo","Ca","Me"), function(pt) {
  get_n <- function(lst) {
    obj <- get_pt_obj(lst, pt)
    if (is.null(obj)) return(NA)
    sum(obj@meta.data$car_status == "CAR+")
  }
  data.frame(Paziente=pt,
             CAR_I     = get_n(I_samples_all),
             CAR_Blood = get_n(AB_blood),
             CAR_Bone  = get_n(AB_bone))
}))
kable(count_df, caption="N. cellule CAR+ per paziente × timepoint", align="c")

# ======================================================================
# PSEUDOBULK BLOOD
# ======================================================================
cat("=== Pseudobulk: I vs Blood AB ===\n\n")
cat("Aggregazione CAR+ cells:\n")

pb_Bo_I     <- make_pseudobulk(get_pt_obj(I_samples_all, "Bo"), "Bo_I")
pb_Bo_blood <- make_pseudobulk(get_pt_obj(AB_blood,      "Bo"), "Bo_Blood")
pb_Ca_I     <- make_pseudobulk(get_pt_obj(I_samples_all, "Ca"), "Ca_I")
pb_Ca_blood <- make_pseudobulk(get_pt_obj(AB_blood,      "Ca"), "Ca_Blood")

pb_blood_list <- Filter(Negate(is.null),
                        list(pb_Bo_I, pb_Bo_blood, pb_Ca_I, pb_Ca_blood))

coldata_blood <- data.frame(
  sample    = c("Bo_I", "Bo_Blood", "Ca_I", "Ca_Blood"),
  patient   = factor(c("Bo","Bo","Ca","Ca")),
  timepoint = factor(c("I","Blood","I","Blood"), levels=c("I","Blood"))
)
rownames(coldata_blood) <- coldata_blood$sample

cat("\nColdata:\n")
print(coldata_blood[, c("patient","timepoint")])

cat("\nDESeq2 ~ patient + timepoint:\n")
de_blood <- run_deseq2(
  pb_list          = pb_blood_list,
  coldata          = coldata_blood,
  design_formula   = ~ patient + timepoint,
  contrast         = c("timepoint", "I", "Blood"),
  comparison_label = "I_vs_Blood"
)

if (!is.null(de_blood) && nrow(de_blood) > 0) {
  write_xlsx(de_blood, file.path(OUT_DIR, "I_vs_Blood", "tables", "DE_I_vs_Blood_DESeq2.xlsx"))
  cat("\nTabella salvata.\n")
}

# ======================================================================
# SUMMARY BLOOD
# ======================================================================
if (!is.null(de_blood) && nrow(de_blood) > 0) {
  n_up <- sum(de_blood$log2FoldChange > 0)
  n_dn <- sum(de_blood$log2FoldChange < 0)
  cat("Geni up at I (log2FC > 0):", n_up, "\n")
  cat("Geni up at Blood (log2FC < 0):", n_dn, "\n\n")
  cat("Top 10 geni up at I:\n")
  print(head(de_blood[de_blood$log2FoldChange > 0,
                      c("gene","log2FoldChange","padj")], 10))
  cat("\nTop 10 geni up at Blood:\n")
  tmp <- de_blood[de_blood$log2FoldChange < 0, ]
  print(head(tmp[order(tmp$log2FoldChange), c("gene","log2FoldChange","padj")], 10))
} else {
  cat("Nessun gene DE significativo trovato.\n")
}

# ======================================================================
# MA BLOOD
# ======================================================================
if (!is.null(de_blood) && nrow(de_blood) > 0) {
  de_blood$sig <- de_blood$padj < PADJ_CUTOFF
  top_up <- head(de_blood[de_blood$log2FoldChange > 0, ], 15)
  top_dn <- head(de_blood[order(de_blood$log2FoldChange), ], 15)
  top_lab <- bind_rows(top_up, top_dn)

  p <- ggplot(de_blood, aes(x = baseMean, y = log2FoldChange,
                             color = log2FoldChange > 0)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_text_repel(data = top_lab, aes(label = gene),
                    size = 2.8, max.overlaps = 20, show.legend = FALSE) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    scale_x_log10() +
    scale_color_manual(values = c("TRUE"="#E64B35","FALSE"="#4DBBD5"),
                       labels = c("TRUE"=paste0("Up at I (n=",sum(de_blood$log2FoldChange>0),")"),
                                  "FALSE"=paste0("Up at Blood (n=",sum(de_blood$log2FoldChange<0),")")),
                       name = NULL) +
    labs(title    = "MA plot — I vs Blood AB",
         subtitle = "DESeq2 pseudobulk · design ~ patient + timepoint",
         x = "Mean expression (log10)", y = "log₂ Fold Change (I / Blood)") +
    theme_classic(base_size = 11) +
    theme(legend.position = "bottom")

  ggsave(file.path(OUT_DIR, "I_vs_Blood", "figures", "MA_plot.png"),
         plot = p, width = 8, height = 5, dpi = 150)
  p
} else {
  cat("Nessun dato da plottare.\n")
}

# ======================================================================
# ENRICH BLOOD
# ======================================================================
if (!is.null(de_blood) && nrow(de_blood) > 0) {
  genes_up_I     <- de_blood$gene[de_blood$log2FoldChange > 0]
  genes_up_blood <- de_blood$gene[de_blood$log2FoldChange < 0]

  cat("Geni up at I:", length(genes_up_I),
      "| Geni up at Blood:", length(genes_up_blood), "\n")

  go_up_I     <- run_go(genes_up_I,     "up_I")
  go_up_blood <- run_go(genes_up_blood, "up_Blood")
  kk_up_I     <- run_kegg(genes_up_I,     "up_I")
  kk_up_blood <- run_kegg(genes_up_blood, "up_Blood")

  d <- file.path(OUT_DIR, "I_vs_Blood", "tables")
  if (!is.null(go_up_I))
    write_xlsx(as.data.frame(go_up_I),     file.path(d, "GO_up_at_I.xlsx"))
  if (!is.null(go_up_blood))
    write_xlsx(as.data.frame(go_up_blood), file.path(d, "GO_up_at_Blood.xlsx"))
  if (!is.null(kk_up_I))
    write_xlsx(as.data.frame(kk_up_I),     file.path(d, "KEGG_up_at_I.xlsx"))
  if (!is.null(kk_up_blood))
    write_xlsx(as.data.frame(kk_up_blood), file.path(d, "KEGG_up_at_Blood.xlsx"))

  enrich_list <- list(
    "GO — up at I"     = go_up_I,
    "GO — up at Blood" = go_up_blood,
    "KEGG — up at I"   = kk_up_I,
    "KEGG — up at Blood" = kk_up_blood
  )
  for (nm in names(enrich_list)) {
    res <- enrich_list[[nm]]
    if (is.null(res)) { cat(nm, ": nessun termine\n"); next }
    p <- dotplot(res, showCategory=12, title=paste("I vs Blood —", nm)) +
      theme(axis.text.y=element_text(size=8))
    fname <- paste0("dotplot_", gsub("[^A-Za-z0-9]","_", nm), ".png")
    ggsave(file.path(OUT_DIR, "I_vs_Blood", "figures", fname),
           plot=p, width=10, height=7, dpi=150)
    print(p); cat("\n")
  }
} else {
  cat("Nessun gene DE — enrichment saltato.\n")
}

# ======================================================================
# DT BLOOD
# ======================================================================
if (!is.null(de_blood) && nrow(de_blood) > 0) {
  de_blood$direction <- ifelse(de_blood$log2FoldChange > 0, "up_at_I", "up_at_Blood")
  DT::datatable(
    de_blood[, c("gene","direction","log2FoldChange","baseMean","padj")],
    caption  = "DE genes — I vs Blood AB (DESeq2 pseudobulk)",
    filter   = "top",
    options  = list(pageLength=20, scrollX=TRUE),
    rownames = FALSE
  ) %>% DT::formatRound(c("log2FoldChange","baseMean","padj"), digits=4)
}

# ======================================================================
# PSEUDOBULK BONE
# ======================================================================
cat("=== Pseudobulk: I vs Bone AB ===\n\n")
cat("Aggregazione CAR+ cells:\n")

pb_Me_I    <- make_pseudobulk(get_pt_obj(I_samples_all, "Me"), "Me_I")
pb_Bo_bone <- make_pseudobulk(get_pt_obj(AB_bone,       "Bo"), "Bo_Bone")
pb_Me_bone <- make_pseudobulk(get_pt_obj(AB_bone,       "Me"), "Me_Bone")

pb_bone_list <- Filter(Negate(is.null),
                       list(pb_Bo_I, pb_Bo_bone, pb_Me_I, pb_Me_bone))

coldata_bone <- data.frame(
  sample    = c("Bo_I", "Bo_Bone", "Me_I", "Me_Bone"),
  patient   = factor(c("Bo","Bo","Me","Me")),
  timepoint = factor(c("I","Bone","I","Bone"), levels=c("I","Bone"))
)
rownames(coldata_bone) <- coldata_bone$sample

cat("\nColdata:\n")
print(coldata_bone[, c("patient","timepoint")])

cat("\nDESeq2 ~ patient + timepoint:\n")
de_bone <- run_deseq2(
  pb_list          = pb_bone_list,
  coldata          = coldata_bone,
  design_formula   = ~ patient + timepoint,
  contrast         = c("timepoint", "I", "Bone"),
  comparison_label = "I_vs_Bone"
)

if (!is.null(de_bone) && nrow(de_bone) > 0) {
  write_xlsx(de_bone, file.path(OUT_DIR, "I_vs_Bone", "tables", "DE_I_vs_Bone_DESeq2.xlsx"))
  cat("\nTabella salvata.\n")
}

# ======================================================================
# SUMMARY BONE
# ======================================================================
if (!is.null(de_bone) && nrow(de_bone) > 0) {
  n_up <- sum(de_bone$log2FoldChange > 0)
  n_dn <- sum(de_bone$log2FoldChange < 0)
  cat("Geni up at I (log2FC > 0):", n_up, "\n")
  cat("Geni up at Bone (log2FC < 0):", n_dn, "\n\n")
  cat("Top 10 geni up at I:\n")
  print(head(de_bone[de_bone$log2FoldChange > 0,
                     c("gene","log2FoldChange","padj")], 10))
  cat("\nTop 10 geni up at Bone:\n")
  tmp <- de_bone[de_bone$log2FoldChange < 0, ]
  print(head(tmp[order(tmp$log2FoldChange), c("gene","log2FoldChange","padj")], 10))
} else {
  cat("Nessun gene DE significativo trovato.\n")
}

# ======================================================================
# MA BONE
# ======================================================================
if (!is.null(de_bone) && nrow(de_bone) > 0) {
  top_up  <- head(de_bone[de_bone$log2FoldChange > 0, ], 15)
  top_dn  <- head(de_bone[order(de_bone$log2FoldChange), ], 15)
  top_lab <- bind_rows(top_up, top_dn)

  p <- ggplot(de_bone, aes(x=baseMean, y=log2FoldChange,
                            color=log2FoldChange > 0)) +
    geom_point(alpha=0.6, size=1.5) +
    geom_text_repel(data=top_lab, aes(label=gene),
                    size=2.8, max.overlaps=20, show.legend=FALSE) +
    geom_hline(yintercept=0, linetype="dashed", color="grey40") +
    scale_x_log10() +
    scale_color_manual(values=c("TRUE"="#E64B35","FALSE"="#00A087"),
                       labels=c("TRUE"=paste0("Up at I (n=",sum(de_bone$log2FoldChange>0),")"),
                                "FALSE"=paste0("Up at Bone (n=",sum(de_bone$log2FoldChange<0),")")),
                       name=NULL) +
    labs(title    = "MA plot — I vs Bone marrow AB",
         subtitle = "DESeq2 pseudobulk · design ~ patient + timepoint",
         x="Mean expression (log10)", y="log₂ Fold Change (I / Bone)") +
    theme_classic(base_size=11) +
    theme(legend.position="bottom")

  ggsave(file.path(OUT_DIR, "I_vs_Bone", "figures", "MA_plot.png"),
         plot=p, width=8, height=5, dpi=150)
  p
} else {
  cat("Nessun dato da plottare.\n")
}

# ======================================================================
# ENRICH BONE
# ======================================================================
if (!is.null(de_bone) && nrow(de_bone) > 0) {
  genes_up_I    <- de_bone$gene[de_bone$log2FoldChange > 0]
  genes_up_bone <- de_bone$gene[de_bone$log2FoldChange < 0]

  cat("Geni up at I:", length(genes_up_I),
      "| Geni up at Bone:", length(genes_up_bone), "\n")

  go_up_I    <- run_go(genes_up_I,    "up_I")
  go_up_bone <- run_go(genes_up_bone, "up_Bone")
  kk_up_I    <- run_kegg(genes_up_I,    "up_I")
  kk_up_bone <- run_kegg(genes_up_bone, "up_Bone")

  d <- file.path(OUT_DIR, "I_vs_Bone", "tables")
  if (!is.null(go_up_I))
    write_xlsx(as.data.frame(go_up_I),    file.path(d, "GO_up_at_I.xlsx"))
  if (!is.null(go_up_bone))
    write_xlsx(as.data.frame(go_up_bone), file.path(d, "GO_up_at_Bone.xlsx"))
  if (!is.null(kk_up_I))
    write_xlsx(as.data.frame(kk_up_I),    file.path(d, "KEGG_up_at_I.xlsx"))
  if (!is.null(kk_up_bone))
    write_xlsx(as.data.frame(kk_up_bone), file.path(d, "KEGG_up_at_Bone.xlsx"))

  enrich_list <- list(
    "GO — up at I"    = go_up_I,
    "GO — up at Bone" = go_up_bone,
    "KEGG — up at I"  = kk_up_I,
    "KEGG — up at Bone" = kk_up_bone
  )
  for (nm in names(enrich_list)) {
    res <- enrich_list[[nm]]
    if (is.null(res)) { cat(nm, ": nessun termine\n"); next }
    p <- dotplot(res, showCategory=12, title=paste("I vs Bone —", nm)) +
      theme(axis.text.y=element_text(size=8))
    fname <- paste0("dotplot_", gsub("[^A-Za-z0-9]","_", nm), ".png")
    ggsave(file.path(OUT_DIR, "I_vs_Bone", "figures", fname),
           plot=p, width=10, height=7, dpi=150)
    print(p); cat("\n")
  }
} else {
  cat("Nessun gene DE — enrichment saltato.\n")
}

# ======================================================================
# DT BONE
# ======================================================================
if (!is.null(de_bone) && nrow(de_bone) > 0) {
  de_bone$direction <- ifelse(de_bone$log2FoldChange > 0, "up_at_I", "up_at_Bone")
  DT::datatable(
    de_bone[, c("gene","direction","log2FoldChange","baseMean","padj")],
    caption  = "DE genes — I vs Bone AB (DESeq2 pseudobulk)",
    filter   = "top",
    options  = list(pageLength=20, scrollX=TRUE),
    rownames = FALSE
  ) %>% DT::formatRound(c("log2FoldChange","baseMean","padj"), digits=4)
}

# ======================================================================
# OVERLAP
# ======================================================================
if (!is.null(de_blood) && !is.null(de_bone) &&
    nrow(de_blood) > 0 && nrow(de_bone) > 0) {

  genes_blood <- de_blood$gene
  genes_bone  <- de_bone$gene
  shared      <- intersect(genes_blood, genes_bone)

  cat("DEG in I vs Blood:", length(genes_blood), "\n")
  cat("DEG in I vs Bone: ", length(genes_bone),  "\n")
  cat("Geni condivisi:   ", length(shared), "\n\n")

  if (length(shared) > 0) {
    shared_df <- merge(
      de_blood[de_blood$gene %in% shared,
               c("gene","log2FoldChange","padj")],
      de_bone[de_bone$gene  %in% shared,
              c("gene","log2FoldChange","padj")],
      by="gene", suffixes=c("_blood","_bone")
    )
    shared_df <- shared_df[sign(shared_df$log2FoldChange_blood) ==
                           sign(shared_df$log2FoldChange_bone), ]
    cat("Geni concordanti (stessa direzione in entrambe):", nrow(shared_df), "\n")
    cat("\nTop geni concordanti up at I:\n")
    tmp <- shared_df[shared_df$log2FoldChange_blood > 0, ]
    print(head(tmp[order(-tmp$log2FoldChange_blood),
                   c("gene","log2FoldChange_blood","log2FoldChange_bone",
                     "padj_blood","padj_bone")], 10))
    cat("\nTop geni concordanti up at AB:\n")
    tmp2 <- shared_df[shared_df$log2FoldChange_blood < 0, ]
    print(head(tmp2[order(tmp2$log2FoldChange_blood),
                    c("gene","log2FoldChange_blood","log2FoldChange_bone",
                      "padj_blood","padj_bone")], 10))

    write_xlsx(shared_df,
               file.path(OUT_DIR, "shared_DEG_concordant.xlsx"))
  }
}

# ======================================================================
# VENN SETUP
# ======================================================================
suppressMessages(library(ggVennDiagram))

if (!is.null(de_blood) && !is.null(de_bone) &&
    nrow(de_blood) > 0 && nrow(de_bone) > 0) {

  up_I_blood  <- de_blood$gene[de_blood$log2FoldChange > 0]
  up_I_bone   <- de_bone$gene [de_bone$log2FoldChange  > 0]
  up_AB_blood <- de_blood$gene[de_blood$log2FoldChange < 0]
  up_AB_bone  <- de_bone$gene [de_bone$log2FoldChange  < 0]

  cat("Up at I  — Blood:", length(up_I_blood),
      "| Bone:", length(up_I_bone),
      "| Overlap:", length(intersect(up_I_blood, up_I_bone)), "\n")
  cat("Up at AB — Blood:", length(up_AB_blood),
      "| Bone:", length(up_AB_bone),
      "| Overlap:", length(intersect(up_AB_blood, up_AB_bone)), "\n")
}

# ======================================================================
# VENN UP I
# ======================================================================
if (!is.null(de_blood) && !is.null(de_bone) &&
    nrow(de_blood) > 0 && nrow(de_bone) > 0) {

  venn_sets_I <- list(
    "Up at I\n(I vs Blood)" = up_I_blood,
    "Up at I\n(I vs Bone)"  = up_I_bone
  )

  p_venn_I <- ggVennDiagram(venn_sets_I, label_alpha = 0,
                              set_color = c("#7B68EE","#7B68EE")) +
    scale_fill_gradient(low = "#EEF0FF", high = "#7B68EE") +
    labs(title    = "Geni up a I — overlap tra le due comparazioni",
         subtitle = "I vs Blood AB  ∩  I vs Bone AB") +
    theme(legend.position = "none",
          plot.title    = element_text(size=13, face="bold"),
          plot.subtitle = element_text(size=11))

  ggsave(file.path(OUT_DIR, "venn_up_at_I.png"),
         plot = p_venn_I, width = 7, height = 5, dpi = 150)
  p_venn_I
}

# ======================================================================
# VENN UP AB
# ======================================================================
if (!is.null(de_blood) && !is.null(de_bone) &&
    nrow(de_blood) > 0 && nrow(de_bone) > 0) {

  venn_sets_AB <- list(
    "Up at Blood AB\n(I vs Blood)" = up_AB_blood,
    "Up at Bone AB\n(I vs Bone)"   = up_AB_bone
  )

  p_venn_AB <- ggVennDiagram(venn_sets_AB, label_alpha = 0,
                               set_color = c("#E64B35","#00A087")) +
    scale_fill_gradient(low = "#FFF5F0", high = "#E64B35") +
    labs(title    = "Geni up ad AB — overlap tra le due comparazioni",
         subtitle = "up Blood AB  ∩  up Bone AB") +
    theme(legend.position = "none",
          plot.title    = element_text(size=13, face="bold"),
          plot.subtitle = element_text(size=11))

  ggsave(file.path(OUT_DIR, "venn_up_at_AB.png"),
         plot = p_venn_AB, width = 7, height = 5, dpi = 150)
  p_venn_AB
}

# ======================================================================
# TABLE SHARED UP I
# ======================================================================
if (!is.null(de_blood) && !is.null(de_bone) &&
    nrow(de_blood) > 0 && nrow(de_bone) > 0) {

  shared_up_I <- intersect(up_I_blood, up_I_bone)

  if (length(shared_up_I) > 0) {
    df_shared_I <- merge(
      de_blood[de_blood$gene %in% shared_up_I,
               c("gene","log2FoldChange","padj")],
      de_bone [de_bone$gene  %in% shared_up_I,
               c("gene","log2FoldChange","padj")],
      by = "gene", suffixes = c("_blood","_bone")
    )
    df_shared_I$lfc_mean <- rowMeans(
      df_shared_I[, c("log2FoldChange_blood","log2FoldChange_bone")])
    df_shared_I <- df_shared_I[order(-df_shared_I$lfc_mean), ]

    write_xlsx(df_shared_I,
               file.path(OUT_DIR, "shared_up_at_I.xlsx"))

    DT::datatable(
      df_shared_I[, c("gene","log2FoldChange_blood","padj_blood",
                      "log2FoldChange_bone","padj_bone","lfc_mean")],
      caption  = paste0("Geni condivisi up a I — n=", nrow(df_shared_I),
                        " (presenti in entrambe le comparazioni)"),
      filter   = "top",
      options  = list(pageLength=20, scrollX=TRUE),
      rownames = FALSE
    ) %>% DT::formatRound(c("log2FoldChange_blood","padj_blood",
                             "log2FoldChange_bone","padj_bone","lfc_mean"),
                           digits=3)
  } else {
    cat("Nessun gene condiviso up a I.\n")
  }
}

# ======================================================================
# TABLE SHARED UP AB
# ======================================================================
if (!is.null(de_blood) && !is.null(de_bone) &&
    nrow(de_blood) > 0 && nrow(de_bone) > 0) {

  shared_up_AB <- intersect(up_AB_blood, up_AB_bone)

  if (length(shared_up_AB) > 0) {
    df_shared_AB <- merge(
      de_blood[de_blood$gene %in% shared_up_AB,
               c("gene","log2FoldChange","padj")],
      de_bone [de_bone$gene  %in% shared_up_AB,
               c("gene","log2FoldChange","padj")],
      by = "gene", suffixes = c("_blood","_bone")
    )
    df_shared_AB$lfc_mean <- rowMeans(
      df_shared_AB[, c("log2FoldChange_blood","log2FoldChange_bone")])
    df_shared_AB <- df_shared_AB[order(df_shared_AB$lfc_mean), ]

    write_xlsx(df_shared_AB,
               file.path(OUT_DIR, "shared_up_at_AB.xlsx"))

    DT::datatable(
      df_shared_AB[, c("gene","log2FoldChange_blood","padj_blood",
                       "log2FoldChange_bone","padj_bone","lfc_mean")],
      caption  = paste0("Geni condivisi up ad AB — n=", nrow(df_shared_AB),
                        " (presenti in entrambe le comparazioni)"),
      filter   = "top",
      options  = list(pageLength=20, scrollX=TRUE),
      rownames = FALSE
    ) %>% DT::formatRound(c("log2FoldChange_blood","padj_blood",
                             "log2FoldChange_bone","padj_bone","lfc_mean"),
                           digits=3)
  } else {
    cat("Nessun gene condiviso up ad AB.\n")
  }
}

# ======================================================================
# TABLE UNIQUE GENES
# ======================================================================
make_unique_table <- function(genes_unique, de_df, comparison_label) {
  if (length(genes_unique) == 0) {
    cat("Nessun gene unico per", comparison_label, "\n")
    return(invisible(NULL))
  }
  df <- de_df[de_df$gene %in% genes_unique,
              c("gene", "log2FoldChange", "padj", "baseMean")]
  df <- df[order(-abs(df$log2FoldChange)), ]
  df
}

if (!is.null(de_blood) && !is.null(de_bone) &&
    nrow(de_blood) > 0 && nrow(de_bone) > 0) {

  # — geni up a I unici a I vs Blood ——————————————————————————————————————————
  unique_up_I_blood  <- setdiff(up_I_blood,  up_I_bone)
  unique_up_I_bone   <- setdiff(up_I_bone,   up_I_blood)
  unique_up_AB_blood <- setdiff(up_AB_blood, up_AB_bone)
  unique_up_AB_bone  <- setdiff(up_AB_bone,  up_AB_blood)

  cat("=== Geni up a I ===\n")
  cat("  Condivisi (I vs Blood AND I vs Bone):", length(intersect(up_I_blood, up_I_bone)), "\n")
  cat("  Unici a I vs Blood:                  ", length(unique_up_I_blood), "\n")
  cat("  Unici a I vs Bone:                   ", length(unique_up_I_bone),  "\n\n")

  cat("=== Geni up ad AB ===\n")
  cat("  Condivisi (I vs Blood AND I vs Bone):", length(intersect(up_AB_blood, up_AB_bone)), "\n")
  cat("  Unici a I vs Blood:                  ", length(unique_up_AB_blood), "\n")
  cat("  Unici a I vs Bone:                   ", length(unique_up_AB_bone),  "\n\n")

  # Costruisce e salva tutte e 4 le tabelle
  df_unique_up_I_blood  <- make_unique_table(unique_up_I_blood,  de_blood, "up-I · I vs Blood")
  df_unique_up_I_bone   <- make_unique_table(unique_up_I_bone,   de_bone,  "up-I · I vs Bone")
  df_unique_up_AB_blood <- make_unique_table(unique_up_AB_blood, de_blood, "up-AB · I vs Blood")
  df_unique_up_AB_bone  <- make_unique_table(unique_up_AB_bone,  de_bone,  "up-AB · I vs Bone")

  if (!is.null(df_unique_up_I_blood))
    write_xlsx(df_unique_up_I_blood,
               file.path(OUT_DIR, "unique_up_at_I_only_Blood.xlsx"))
  if (!is.null(df_unique_up_I_bone))
    write_xlsx(df_unique_up_I_bone,
               file.path(OUT_DIR, "unique_up_at_I_only_Bone.xlsx"))
  if (!is.null(df_unique_up_AB_blood))
    write_xlsx(df_unique_up_AB_blood,
               file.path(OUT_DIR, "unique_up_at_AB_only_Blood.xlsx"))
  if (!is.null(df_unique_up_AB_bone))
    write_xlsx(df_unique_up_AB_bone,
               file.path(OUT_DIR, "unique_up_at_AB_only_Bone.xlsx"))
}

# ======================================================================
# DT UNIQUE UP I BLOOD
# ======================================================================
print(DT::datatable(df_unique_up_I_blood,
  caption  = paste0("Geni up a I SOLO in I vs Blood — n=", nrow(df_unique_up_I_blood)),
  filter="top", options=list(pageLength=20, scrollX=TRUE), rownames=FALSE) %>%
  DT::formatRound(c("log2FoldChange","padj","baseMean"), digits=3))

# ======================================================================
# DT UNIQUE UP I BONE
# ======================================================================
print(DT::datatable(df_unique_up_I_bone,
  caption  = paste0("Geni up a I SOLO in I vs Bone — n=", nrow(df_unique_up_I_bone)),
  filter="top", options=list(pageLength=20, scrollX=TRUE), rownames=FALSE) %>%
  DT::formatRound(c("log2FoldChange","padj","baseMean"), digits=3))

# ======================================================================
# DT UNIQUE UP AB BLOOD
# ======================================================================
print(DT::datatable(df_unique_up_AB_blood,
  caption  = paste0("Geni up ad AB SOLO in I vs Blood — n=", nrow(df_unique_up_AB_blood)),
  filter="top", options=list(pageLength=20, scrollX=TRUE), rownames=FALSE) %>%
  DT::formatRound(c("log2FoldChange","padj","baseMean"), digits=3))

# ======================================================================
# DT UNIQUE UP AB BONE
# ======================================================================
print(DT::datatable(df_unique_up_AB_bone,
  caption  = paste0("Geni up ad AB SOLO in I vs Bone — n=", nrow(df_unique_up_AB_bone)),
  filter="top", options=list(pageLength=20, scrollX=TRUE), rownames=FALSE) %>%
  DT::formatRound(c("log2FoldChange","padj","baseMean"), digits=3))

# ======================================================================
# FINAL SUMMARY
# ======================================================================
cat("=== RIEPILOGO P4 — PSEUDOBULK DESeq2 ===\n\n")
cat("Design: ~ patient + timepoint\n")
cat("padj <", PADJ_CUTOFF, "| Filtro: MT/RPS/RPL/MALAT1/NEAT1/TRAV/TRB\n\n")

for (nm in c("I vs Blood AB","I vs Bone AB")) {
  df <- if (nm=="I vs Blood AB") de_blood else de_bone
  if (is.null(df) || nrow(df)==0) {
    cat(nm, ": nessun gene DE\n")
  } else {
    cat(nm, ":", nrow(df), "DEG totali |",
        sum(df$log2FoldChange>0), "up at I |",
        sum(df$log2FoldChange<0), "up at AB\n")
  }
}
