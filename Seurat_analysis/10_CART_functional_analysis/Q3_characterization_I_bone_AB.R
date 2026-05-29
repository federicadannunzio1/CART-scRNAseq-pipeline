# ============================================================
#  Q3: Caratterizzazione CAR+ in I e nel midollo AB
#      + Analisi di selezione (cosa predice la persistenza?)
#
#  Domande biologiche:
#    A. Che caratteristiche hanno le CAR+ nel prodotto di
#       infusione (I)? Differiscono tra pazienti?
#    B. Che caratteristiche hanno le CAR+ nel midollo osseo
#       post-infusione (AB)?
#    C. Esiste qualcosa nel prodotto I che ha indotto
#       selezione/influenzato negativamente la persistenza?
#
#  Logica:
#    - N=3 → NO test statistico formale inter-paziente.
#      Tutti i confronti sono descrittivi + coerenza.
#    - Per DEG intra-paziente (CAR+ vs CAR-) si usa Wilcoxon
#      ma i risultati vanno interpretati con cautela (n piccolo).
#    - La firma di "selezione" è costruita confrontando
#      Bo_CAR+_I (best responder) vs Ca+Me_CAR+_I (poor/non
#      responder) — approccio esplorativo.
#    - Approccio agnostico CD4/CD8 (come in Q1b).
#
#  MEMORIA: usa file RDS separati per rispettare il limite
#  di 8 GB RAM (244 MB per I, 807 MB per AB).
#
#  Output in out_dir/:
#    Q3_A_interpatient_module_scores_I.png
#    Q3_A_Bo_vs_CaMepool_DEG.png / .xlsx
#    Q3_B_bone_AB_module_scores.png
#    Q3_B_Bo_vs_Me_bone_AB_DEG.xlsx
#    Q3_B_CAR_longitudinal_I_vs_bone_AB.png
#    Q3_C_selection_summary.png
#    Q3_C_persistence_signature_heatmap.png
# ============================================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(scales)
library(openxlsx)
library(ggrepel)

# ── Percorsi ─────────────────────────────────────────────────
rds_I  <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/CART/Code/Seurat_analysis/2_annotation/all_I_samples_annotated.rds"
rds_AB <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/CART/Code/Seurat_analysis/2_annotation/all_AB_samples_annotated.rds"
out_dir <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/CART/Code/Seurat_analysis/10_CART_functional_analysis/Q3_characterization/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

section <- function(title)
  cat(paste0("\n", strrep("=", 65), "\n  ", title, "\n",
             strrep("=", 65), "\n"))

# ── Definizioni ───────────────────────────────────────────────

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
  "Naive-like"   = c("Naive CD4+ T cells","Naive CD8+ T cells"),
  "Memory-like"  = c("Memory T cells","Th1 cells","Th2 cells","Th17 cells","Tfh cells"),
  "Effector"     = c("Effector CD4+ T cells","Cytotoxic CD8+ T cells"),
  "Regulatory"   = c("Tregs"),
  "Proliferating"= c("Proliferating CD4+ T cells","Proliferating CD8+ T cells")
)

# Esiti clinici per contesto interpretativo
OUTCOME <- c(Bo = "Expansion (1.9%→22%)", Ca = "Failure (19.1%→0%)", Me = "Partial (9.9%→5%)")
OUTCOME_COLOR <- c(Bo = "#E64B35", Ca = "#4DBBD5", Me = "#00A087")

get_car_status <- function(obj) {
  meta <- obj@meta.data
  for (col in c("IS_CAR_ALLIN_scREP","IS_CAR","CAR")) {
    if (col %in% colnames(meta)) {
      vals <- as.character(meta[[col]])
      return(ifelse(grepl("^(YES|TRUE|yes|true|1)$", vals), "CAR+", "CAR-"))
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

add_scores <- function(obj) {
  for (sig in names(SIGNATURES)) {
    genes_ok <- intersect(SIGNATURES[[sig]], rownames(obj))
    if (length(genes_ok) < 2) next
    sc <- paste0("score_", sig)
    obj <- AddModuleScore(obj, features = list(genes_ok), name = sc, seed = 42)
    obj@meta.data[[sc]] <- obj@meta.data[[paste0(sc,"1")]]
    obj@meta.data[[paste0(sc,"1")]] <- NULL
  }
  obj
}

# ============================================================
# BLOCCO A — CAMPIONI I
# ============================================================
section("Caricamento campioni I")
I_samples <- readRDS(rds_I)
cat("Campioni:", paste(names(I_samples), collapse=", "), "\n")

# Aggiungi metadati e module scores per ogni campione
for (sname in names(I_samples)) {
  patient <- sub("_bone_I$","", sname)
  I_samples[[sname]]$car_status <- get_car_status(I_samples[[sname]])
  I_samples[[sname]]$patient    <- patient
  I_samples[[sname]]$fs         <- map_fs(as.character(I_samples[[sname]]@meta.data$cell_type))
  cat(sprintf("  %s: CAR+ = %d, CAR- = %d\n",
              sname,
              sum(I_samples[[sname]]$car_status == "CAR+"),
              sum(I_samples[[sname]]$car_status == "CAR-")))
  I_samples[[sname]] <- add_scores(I_samples[[sname]])
}

# ── SEZIONE A1: Module scores inter-paziente (CAR+ in I) ─────
section("A1: Module scores inter-paziente — CAR+ in I")

score_cols <- paste0("score_", names(SIGNATURES))

# Costruisci dataframe CAR+ di tutti i pazienti in I
meta_I_carpos <- lapply(names(I_samples), function(sname) {
  meta <- I_samples[[sname]]@meta.data
  meta_car <- meta[meta$car_status == "CAR+", ]
  meta_car$sname <- sname
  meta_car
}) %>% bind_rows()

# Heatmap media module scores × paziente
heat_scores <- meta_I_carpos %>%
  group_by(patient) %>%
  summarise(across(all_of(intersect(score_cols, colnames(meta_I_carpos))), mean, na.rm = TRUE),
            .groups = "drop") %>%
  pivot_longer(-patient, names_to = "signature", values_to = "score") %>%
  mutate(signature = gsub("score_","", signature))

p_heat_scores <- ggplot(heat_scores,
                        aes(x = patient, y = signature, fill = score)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.3f", score)), size = 3.5) +
  scale_fill_gradient2(low = "#4DBBD5", mid = "white", high = "#E64B35",
                       midpoint = 0, name = "Score medio") +
  labs(
    title    = "Module scores CAR+ nel prodotto di infusione (I)",
    subtitle = "Confronto inter-paziente — agnostico CD4/CD8\nBo=espansione | Ca=fallimento | Me=parziale",
    x = "Paziente", y = "Firma genica"
  ) +
  theme_classic(base_size = 12)

ggsave(file.path(out_dir, "Q3_A1_interpatient_module_scores_I.png"),
       p_heat_scores, width = 7, height = 6, dpi = 300)
cat("Salvato: Q3_A1_interpatient_module_scores_I.png\n")

# Violin plot per ogni firma: 3 pazienti affiancati (solo CAR+)
vln_list <- lapply(intersect(score_cols, colnames(meta_I_carpos)), function(sc) {
  sig_name <- gsub("score_","", sc)
  df <- meta_I_carpos %>% select(patient, score = !!sym(sc)) %>% filter(!is.na(score))
  ggplot(df, aes(x = patient, y = score, fill = patient)) +
    geom_violin(trim = TRUE, alpha = 0.85, color = "white") +
    geom_boxplot(width = 0.12, fill = "white", alpha = 0.7, outlier.shape = NA) +
    scale_fill_manual(values = OUTCOME_COLOR, guide = "none") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.4) +
    labs(title = sig_name, x = NULL, y = "Score") +
    theme_classic(base_size = 9) +
    theme(plot.title = element_text(face = "bold", size = 9))
})

p_vln_A <- wrap_plots(vln_list, ncol = 4) +
  plot_annotation(
    title    = "CAR+ cells in I — Module scores per paziente",
    subtitle = "Rosso=Bo(espansione) | Blu=Ca(fallimento) | Verde=Me(parziale)\nN=3 pazienti → solo confronto descrittivo",
    theme    = theme(plot.title    = element_text(face = "bold"),
                     plot.subtitle = element_text(size = 9, color = "gray40"))
  )

ggsave(file.path(out_dir, "Q3_A1_violin_module_scores_I_bypatient.png"),
       p_vln_A, width = 14, height = 8, dpi = 300)
cat("Salvato: Q3_A1_violin_module_scores_I_bypatient.png\n")

# ── SEZIONE A2: DEG inter-paziente CAR+ in I ─────────────────
# Bo (best responder) vs Ca+Me pool (poor responders)
# NOTA METODOLOGICA: N=3, questo è esplorativo.
# Non è un test statistico formale ma una ricerca di pattern.
section("A2: DEG Bo_CAR+_I vs Ca+Me_CAR+_I (firma di persistenza)")

# Estrai cellule CAR+
Bo_I_carpos  <- subset(I_samples$Bo_bone_I,  cells = WhichCells(I_samples$Bo_bone_I,  expression = car_status == "CAR+"))
Ca_I_carpos  <- subset(I_samples$Ca_bone_I,  cells = WhichCells(I_samples$Ca_bone_I,  expression = car_status == "CAR+"))
Me_I_carpos  <- subset(I_samples$Me_bone_I,  cells = WhichCells(I_samples$Me_bone_I,  expression = car_status == "CAR+"))

cat(sprintf("  Bo CAR+: %d | Ca CAR+: %d | Me CAR+: %d\n",
            ncol(Bo_I_carpos), ncol(Ca_I_carpos), ncol(Me_I_carpos)))

# Merge poor responders
pool_CaMe <- merge(Ca_I_carpos, y = Me_I_carpos,
                   add.cell.ids = c("Ca","Me"))
pool_CaMe <- JoinLayers(pool_CaMe)

# Merge tutto per FindMarkers
all_I_carpos <- merge(Bo_I_carpos, y = pool_CaMe,
                      add.cell.ids = c("Bo","CaMe"))
all_I_carpos <- JoinLayers(all_I_carpos)
all_I_carpos$responder_group <- ifelse(
  grepl("^Bo", colnames(all_I_carpos)), "Bo_expander", "CaMe_poor"
)
Idents(all_I_carpos) <- "responder_group"

# Normalizzazione
all_I_carpos <- NormalizeData(all_I_carpos, verbose = FALSE)

deg_persistence <- FindMarkers(
  all_I_carpos,
  ident.1  = "Bo_expander",
  ident.2  = "CaMe_poor",
  test.use = "wilcox",
  min.pct  = 0.1,
  logfc.threshold = 0.2,
  verbose  = FALSE
)
deg_persistence$gene <- rownames(deg_persistence)
deg_persistence <- deg_persistence %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC))

cat(sprintf("  DEG sign. (Bo vs Ca+Me): %d geni\n", nrow(deg_persistence)))
cat("  Top 10 UP in Bo:\n")
print(head(deg_persistence[deg_persistence$avg_log2FC > 0, c("gene","avg_log2FC","p_val_adj")], 10))
cat("  Top 10 DOWN in Bo (= UP in Ca+Me):\n")
print(head(deg_persistence[deg_persistence$avg_log2FC < 0, c("gene","avg_log2FC","p_val_adj")], 10))

# Volcano plot
deg_plot <- deg_persistence %>%
  mutate(
    direction = case_when(
      avg_log2FC >  0.5 & p_val_adj < 0.01 ~ "Bo_higher",
      avg_log2FC < -0.5 & p_val_adj < 0.01 ~ "CaMe_higher",
      TRUE ~ "ns"
    ),
    label = ifelse(abs(avg_log2FC) > 1 & p_val_adj < 1e-5, gene, NA)
  )

p_volcano <- ggplot(deg_plot, aes(x = avg_log2FC, y = -log10(p_val_adj),
                                   color = direction, label = label)) +
  geom_point(size = 1.5, alpha = 0.7) +
  geom_text_repel(size = 3, max.overlaps = 15, na.rm = TRUE) +
  scale_color_manual(values = c(Bo_higher = "#E64B35", CaMe_higher = "#4DBBD5", ns = "gray70"),
                     labels = c("Bo_higher"   = "Higher in Bo (expander)",
                                "CaMe_higher" = "Higher in Ca+Me (poor)",
                                "ns"          = "n.s."),
                     name = NULL) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "gray60") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "gray60") +
  labs(
    title    = "Firma di persistenza: Bo_CAR+_I vs Ca+Me_CAR+_I",
    subtitle = "Geni più alti in Bo = potenziali predittori di espansione in vivo\nN.B.: n=3 pazienti → analisi esplorativa",
    x        = "avg log2FC (Bo vs Ca+Me)",
    y        = "-log10(p_val_adj)"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "Q3_A2_Bo_vs_CaMe_volcano.png"),
       p_volcano, width = 9, height = 7, dpi = 300)
cat("Salvato: Q3_A2_Bo_vs_CaMe_volcano.png\n")

# Salva DEG in Excel
wb <- createWorkbook()
addWorksheet(wb, "Persistence_DEG")
writeData(wb, "Persistence_DEG", deg_persistence)
saveWorkbook(wb, file.path(out_dir, "Q3_A2_persistence_signature_DEG.xlsx"),
             overwrite = TRUE)
cat("Salvato: Q3_A2_persistence_signature_DEG.xlsx\n")

# Heatmap geni chiave (top 20 per direzione)
key_genes_Bo   <- head(deg_persistence$gene[deg_persistence$avg_log2FC > 0], 20)
key_genes_CaMe <- head(deg_persistence$gene[deg_persistence$avg_log2FC < 0], 20)
key_genes_all  <- c(key_genes_Bo, key_genes_CaMe)

if (length(key_genes_all) >= 4) {
  p_dot_deg <- DotPlot(all_I_carpos,
                       features  = key_genes_all,
                       group.by  = "responder_group") +
    RotatedAxis() +
    scale_color_gradientn(colours = c("#4DBBD5","white","#E64B35")) +
    labs(title    = "Top geni: Bo_expander vs Ca+Me_poor (CAR+ in I)",
         subtitle = "Sinistra = geni più alti in Bo | Destra = più alti in Ca+Me") +
    theme(axis.text.x = element_text(size = 8))

  ggsave(file.path(out_dir, "Q3_A2_Bo_vs_CaMe_dotplot.png"),
         p_dot_deg,
         width  = max(14, length(key_genes_all) * 0.6),
         height = 5, dpi = 300)
  cat("Salvato: Q3_A2_Bo_vs_CaMe_dotplot.png\n")
}

# Libera memoria prima di caricare AB
rm(I_samples, Bo_I_carpos, Ca_I_carpos, Me_I_carpos,
   pool_CaMe, all_I_carpos, meta_I_carpos, deg_persistence)
invisible(gc()); invisible(gc())

# ============================================================
# BLOCCO B — MIDOLLO OSSEO AB
# ============================================================
section("Caricamento campioni AB (solo bone marrow)")

AB_samples <- readRDS(rds_AB)
cat("Campioni AB:", paste(names(AB_samples), collapse=", "), "\n")

# Filtra solo bone marrow (non blood AB)
bone_AB_names <- grep("bone_AB$", names(AB_samples), value = TRUE)
cat("Bone marrow AB:", paste(bone_AB_names, collapse=", "), "\n")

for (sname in bone_AB_names) {
  patient <- sub("_bone_AB$","", sname)
  AB_samples[[sname]]$car_status <- get_car_status(AB_samples[[sname]])
  AB_samples[[sname]]$patient    <- patient
  AB_samples[[sname]]$fs         <- map_fs(as.character(AB_samples[[sname]]@meta.data$cell_type))
  n_car <- sum(AB_samples[[sname]]$car_status == "CAR+")
  cat(sprintf("  %s: CAR+ = %d\n", sname, n_car))
  if (n_car < 10) {
    cat(sprintf("    [SKIP] Troppo poche cellule CAR+ (%d) per analisi\n", n_car))
    next
  }
  AB_samples[[sname]] <- add_scores(AB_samples[[sname]])
}

# ── SEZIONE B1: Module scores CAR+ nel midollo AB ────────────
section("B1: Module scores CAR+ nel midollo AB")

# Campioni con CAR+ sufficienti: Bo_bone_AB e Me_bone_AB
bone_AB_usable <- names(which(sapply(bone_AB_names, function(s)
  sum(AB_samples[[s]]$car_status == "CAR+") >= 10)))

meta_bone_AB_carpos <- lapply(bone_AB_usable, function(sname) {
  meta <- AB_samples[[sname]]@meta.data
  meta[meta$car_status == "CAR+", ]
}) %>% bind_rows()

if (nrow(meta_bone_AB_carpos) > 0) {
  # Heatmap scores × paziente in bone AB
  heat_AB <- meta_bone_AB_carpos %>%
    group_by(patient) %>%
    summarise(across(all_of(intersect(paste0("score_", names(SIGNATURES)),
                                      colnames(meta_bone_AB_carpos))),
                     mean, na.rm = TRUE),
              .groups = "drop") %>%
    pivot_longer(-patient, names_to = "signature", values_to = "score") %>%
    mutate(signature = gsub("score_","", signature))

  p_heat_AB <- ggplot(heat_AB, aes(x = patient, y = signature, fill = score)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.3f", score)), size = 3.5) +
    scale_fill_gradient2(low = "#4DBBD5", mid = "white", high = "#E64B35",
                         midpoint = 0, name = "Score medio") +
    labs(title    = "Module scores CAR+ nel midollo osseo (AB)",
         subtitle = "Solo pazienti con ≥10 cellule CAR+ in bone AB\n(Ca_bone_AB = 0 cellule CAR+, escluso)",
         x = "Paziente", y = "Firma genica") +
    theme_classic(base_size = 12)

  ggsave(file.path(out_dir, "Q3_B1_module_scores_bone_AB.png"),
         p_heat_AB, width = 6, height = 6, dpi = 300)
  cat("Salvato: Q3_B1_module_scores_bone_AB.png\n")

  # Violin: Bo vs Me in bone AB (CAR+)
  vln_AB <- lapply(intersect(paste0("score_",names(SIGNATURES)),
                              colnames(meta_bone_AB_carpos)), function(sc) {
    sig_name <- gsub("score_","", sc)
    df <- meta_bone_AB_carpos %>% select(patient, score = !!sym(sc)) %>% filter(!is.na(score))
    ggplot(df, aes(x = patient, y = score, fill = patient)) +
      geom_violin(trim = TRUE, alpha = 0.85, color = "white") +
      geom_boxplot(width = 0.12, fill = "white", alpha = 0.7, outlier.shape = NA) +
      scale_fill_manual(values = OUTCOME_COLOR, guide = "none") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.4) +
      labs(title = sig_name, x = NULL, y = "Score") +
      theme_classic(base_size = 9) +
      theme(plot.title = element_text(face = "bold", size = 9))
  })

  p_vln_AB <- wrap_plots(vln_AB, ncol = 4) +
    plot_annotation(
      title    = "CAR+ cells nel midollo AB — Module scores per paziente",
      subtitle = "Bo (espansione) vs Me (parziale) | Ca escluso (0 cellule)",
      theme    = theme(plot.title = element_text(face = "bold"),
                       plot.subtitle = element_text(size = 9, color = "gray40"))
    )

  ggsave(file.path(out_dir, "Q3_B1_violin_module_scores_bone_AB.png"),
         p_vln_AB, width = 14, height = 8, dpi = 300)
  cat("Salvato: Q3_B1_violin_module_scores_bone_AB.png\n")
}

# ── SEZIONE B2: DEG Bo_bone_AB vs Me_bone_AB (CAR+ only) ─────
section("B2: DEG Bo_bone_AB vs Me_bone_AB (CAR+)")

bo_bone_AB_carpos  <- subset(AB_samples$Bo_bone_AB,
                              cells = WhichCells(AB_samples$Bo_bone_AB,
                                                 expression = car_status == "CAR+"))
me_bone_AB_carpos  <- NULL
if ("Me_bone_AB" %in% names(AB_samples)) {
  me_n <- sum(AB_samples$Me_bone_AB$car_status == "CAR+")
  if (me_n >= 10)
    me_bone_AB_carpos <- subset(AB_samples$Me_bone_AB,
                                 cells = WhichCells(AB_samples$Me_bone_AB,
                                                    expression = car_status == "CAR+"))
}

if (!is.null(me_bone_AB_carpos)) {
  cat(sprintf("  Bo_bone_AB CAR+: %d | Me_bone_AB CAR+: %d\n",
              ncol(bo_bone_AB_carpos), ncol(me_bone_AB_carpos)))

  merged_bone_AB <- merge(bo_bone_AB_carpos, y = me_bone_AB_carpos,
                          add.cell.ids = c("Bo","Me"))
  merged_bone_AB <- JoinLayers(merged_bone_AB)
  merged_bone_AB <- NormalizeData(merged_bone_AB, verbose = FALSE)
  Idents(merged_bone_AB) <- "patient"

  deg_bone_AB <- FindMarkers(
    merged_bone_AB,
    ident.1  = "Bo",
    ident.2  = "Me",
    test.use = "wilcox",
    min.pct  = 0.1,
    logfc.threshold = 0.2,
    verbose  = FALSE
  )
  deg_bone_AB$gene <- rownames(deg_bone_AB)
  deg_bone_AB <- deg_bone_AB %>% filter(p_val_adj < 0.05) %>%
    arrange(desc(avg_log2FC))

  cat(sprintf("  DEG sign. (Bo vs Me bone AB): %d geni\n", nrow(deg_bone_AB)))

  wb_b <- createWorkbook()
  addWorksheet(wb_b, "BoneAB_DEG")
  writeData(wb_b, "BoneAB_DEG", deg_bone_AB)
  saveWorkbook(wb_b, file.path(out_dir, "Q3_B2_BoneAB_Bo_vs_Me_DEG.xlsx"),
               overwrite = TRUE)
  cat("Salvato: Q3_B2_BoneAB_Bo_vs_Me_DEG.xlsx\n")

  # Dotplot top geni
  top_bo  <- head(deg_bone_AB$gene[deg_bone_AB$avg_log2FC > 0], 15)
  top_me  <- head(deg_bone_AB$gene[deg_bone_AB$avg_log2FC < 0], 15)
  top_all <- c(top_bo, top_me)

  if (length(top_all) >= 4) {
    p_dot_B2 <- DotPlot(merged_bone_AB, features = top_all, group.by = "patient") +
      RotatedAxis() +
      scale_color_gradientn(colours = c("#4DBBD5","white","#E64B35")) +
      labs(title    = "Top DEG: Bo_bone_AB vs Me_bone_AB (CAR+ only)",
           subtitle = "Sinistra = più alto in Bo | Destra = più alto in Me") +
      theme(axis.text.x = element_text(size = 8))

    ggsave(file.path(out_dir, "Q3_B2_BoneAB_dotplot.png"),
           p_dot_B2,
           width  = max(12, length(top_all) * 0.6),
           height = 5, dpi = 300)
    cat("Salvato: Q3_B2_BoneAB_dotplot.png\n")
  }

  rm(merged_bone_AB, deg_bone_AB)
}

rm(bo_bone_AB_carpos, me_bone_AB_carpos)

# ── SEZIONE B3: Confronto longitudinale CAR+ I vs bone AB ────
section("B3: Confronto longitudinale CAR+ (I vs bone AB)")

# Ricarica I per confronto longitudinale
rm(AB_samples); invisible(gc()); invisible(gc())
cat("Ricaricamento I per confronto longitudinale...\n")
I_samples <- readRDS(rds_I)

for (sname in names(I_samples)) {
  I_samples[[sname]]$car_status <- get_car_status(I_samples[[sname]])
  I_samples[[sname]]$patient    <- sub("_bone_I$","", sname)
  I_samples[[sname]] <- add_scores(I_samples[[sname]])
}

# Ricarica AB
AB_samples <- readRDS(rds_AB)
bone_AB_names <- grep("bone_AB$", names(AB_samples), value = TRUE)
for (sname in bone_AB_names) {
  AB_samples[[sname]]$car_status <- get_car_status(AB_samples[[sname]])
  AB_samples[[sname]]$patient    <- sub("_bone_AB$","", sname)
  AB_samples[[sname]] <- add_scores(AB_samples[[sname]])
}

# Costruisci dataframe longitudinale (CAR+ only)
long_df <- bind_rows(
  lapply(names(I_samples), function(s) {
    m <- I_samples[[s]]@meta.data
    m[m$car_status == "CAR+", ] %>% mutate(timepoint = "I")
  }),
  lapply(bone_AB_names, function(s) {
    m <- AB_samples[[s]]@meta.data
    cp <- m[m$car_status == "CAR+", ]
    if (nrow(cp) < 10) return(NULL)
    cp %>% mutate(timepoint = "bone_AB")
  })
) %>% filter(!is.null(.))

sc_cols <- intersect(paste0("score_",names(SIGNATURES)), colnames(long_df))

long_summary <- long_df %>%
  group_by(patient, timepoint) %>%
  summarise(across(all_of(sc_cols), mean, na.rm = TRUE), .groups = "drop") %>%
  pivot_longer(all_of(sc_cols), names_to = "signature", values_to = "score") %>%
  mutate(signature = gsub("score_","", signature),
         timepoint = factor(timepoint, levels = c("I","bone_AB")))

p_long <- ggplot(long_summary,
                 aes(x = timepoint, y = score,
                     color = patient, group = patient)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(values = OUTCOME_COLOR) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  facet_wrap(~ signature, ncol = 4, scales = "free_y") +
  labs(
    title    = "Evoluzione module scores CAR+ : I → midollo AB",
    subtitle = "Solo pazienti con CAR+ rilevabili in bone AB (Ca escluso)\nRosso=Bo | Verde=Me",
    x = "Timepoint", y = "Score medio CAR+"
  ) +
  theme_classic(base_size = 11) +
  theme(strip.background = element_rect(fill = "#F0F0F0"),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")

ggsave(file.path(out_dir, "Q3_B3_longitudinal_I_to_bone_AB.png"),
       p_long, width = 14, height = 10, dpi = 300)
cat("Salvato: Q3_B3_longitudinal_I_to_bone_AB.png\n")

# ============================================================
# SEZIONE C: ANALISI DI SELEZIONE
# ============================================================
section("C: Analisi di selezione — predittori di persistenza")

# C1: Proporzione Treg e Proliferating nelle CAR+ in I
# (le proporzioni sono già in Q1b_functional_summary.xlsx
#  ma le ricalcoliamo qui per completezza)

prop_I_carpos <- lapply(names(I_samples), function(sname) {
  meta <- I_samples[[sname]]@meta.data
  patient <- unique(meta$patient)
  meta_car <- meta[meta$car_status == "CAR+", ]
  if (nrow(meta_car) == 0) return(NULL)
  fs <- map_fs(as.character(meta_car$cell_type))
  data.frame(
    patient   = patient,
    Regulatory   = mean(fs == "Regulatory"),
    Proliferating= mean(fs == "Proliferating"),
    Effector     = mean(fs == "Effector"),
    Naive_like   = mean(fs == "Naive-like"),
    Memory_like  = mean(fs == "Memory-like"),
    n_carpos     = nrow(meta_car),
    # Esito clinico (% CAR+ in AB)
    AB_expansion = c(Bo=22, Ca=0, Me=4.9)[patient]
  )
}) %>% bind_rows()

cat("\n=== Proporzioni stati funzionali nelle CAR+ in I ===\n")
print(prop_I_carpos)

# Plot correlazione: % Treg in I vs espansione AB
prop_long <- prop_I_carpos %>%
  pivot_longer(c(Regulatory, Proliferating, Effector, Naive_like, Memory_like),
               names_to = "state", values_to = "prop_in_I")

p_selection <- ggplot(prop_long,
                      aes(x = 100 * prop_in_I, y = AB_expansion,
                          color = patient, label = patient)) +
  geom_point(size = 4) +
  geom_text_repel(size = 4, fontface = "bold") +
  scale_color_manual(values = OUTCOME_COLOR, guide = "none") +
  facet_wrap(~ state, ncol = 3, scales = "free_x") +
  labs(
    title    = "Predittori di persistenza/espansione: proporzioni CAR+ in I",
    subtitle = "Y = % CAR+ rilevata nel sangue/midollo AB\n(N=3 pazienti — SOLO descrittivo, no test statistico)",
    x = "% cellule nello stato nel prodotto I (CAR+ only)",
    y = "% CAR+ in AB (sangue)"
  ) +
  theme_classic(base_size = 11) +
  theme(strip.background = element_rect(fill = "#F0F0F0"),
        strip.text = element_text(face = "bold"))

ggsave(file.path(out_dir, "Q3_C_selection_predictors.png"),
       p_selection, width = 12, height = 8, dpi = 300)
cat("Salvato: Q3_C_selection_predictors.png\n")

# C2: Module score Effector in I vs espansione
eff_scores <- long_df %>%
  filter(timepoint == "I") %>%
  group_by(patient) %>%
  summarise(mean_Effector = mean(score_Effector, na.rm=TRUE),
            mean_Exhaustion = mean(score_Exhaustion, na.rm=TRUE),
            mean_Proliferation = mean(score_Proliferation, na.rm=TRUE),
            mean_Memory_Stemness = mean(score_Memory_Stemness, na.rm=TRUE),
            .groups = "drop") %>%
  mutate(AB_expansion = c(Bo=22, Ca=0, Me=4.9)[patient])

cat("\n=== Module scores CAR+ in I vs esito AB ===\n")
print(eff_scores)

# Salva tabella predittori
wb_c <- createWorkbook()
addWorksheet(wb_c, "FunctionalState_Props")
writeData(wb_c, "FunctionalState_Props", prop_I_carpos)
addWorksheet(wb_c, "ModuleScores_vs_Outcome")
writeData(wb_c, "ModuleScores_vs_Outcome", eff_scores)
saveWorkbook(wb_c, file.path(out_dir, "Q3_C_selection_summary.xlsx"),
             overwrite = TRUE)
cat("Salvato: Q3_C_selection_summary.xlsx\n")

section("COMPLETATO")
cat("Output in:", out_dir, "\n\nFile prodotti:\n")
for (f in list.files(out_dir, full.names = FALSE))
  cat(sprintf("  %s\n", f))
