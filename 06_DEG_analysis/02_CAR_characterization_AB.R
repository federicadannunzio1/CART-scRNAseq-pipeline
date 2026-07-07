# CAR+ Cell Characterization — AB (Bone Marrow)
# Converted from P2_CAR_in_bone_AB.Rmd

# ======================================================================
# OUTCOMES TABLE
# ======================================================================
library(knitr)
outcomes <- data.frame(
  Paziente        = c("Bo", "Ca", "Me"),
  CAR_pos_in_I    = c("142 (1.9%)", "215 (19.1%)", "78 (9.9%)"),
  CAR_pos_in_AB   = c("958 in midollo", "0 (escluso)", "260 in midollo"),
  Expansion_index = c("22% (sangue+midollo)", "~0% (fallimento)", "4.9% (midollo)"),
  Incluso_P2      = c("SI", "NO (0 CAR+)", "SI")
)
kable(outcomes,
      caption = "Clinical outcomes: Ca excluded from P2 — no CAR+ in bone marrow AB",
      align = "c")

# ======================================================================
# LOAD LIBRARIES
# ======================================================================
suppressMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(readr)
  library(DT)
  library(knitr)
  library(ggrepel)
  library(RColorBrewer)
})

# ======================================================================
# DEFINE PATHS
# ======================================================================
# ── Percorsi ────────────────────────────────────────────────
RDS_AB  <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/CART/Code/Seurat_analysis/2_annotation/all_AB_samples_annotated.rds"
RDS_I   <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/CART/Code/Seurat_analysis/2_annotation/all_I_samples_annotated.rds"

VDJ_BASE <- "/Users/federicadannunzio/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/caruana-project/CART/Data/output_allineamento_original_no_CAR_si_VDJ"

# ── Definizioni biologiche ───────────────────────────────────
PATIENT_COLORS <- c(Bo = "#E64B35", Ca = "#4DBBD5", Me = "#00A087")
OUTCOME_LABEL  <- c(Bo = "Bo (expansion)", Ca = "Ca (failure)", Me = "Me (partial)")

# T cell types: CAR- is defined as T cells without CAR (excludes monocytes, NK, plasma, erythroid, etc.)
T_CELL_TYPES <- c(
  "Cytotoxic CD8+ T cells", "Effector CD4+ T cells", "Memory T cells",
  "Naive CD4+ T cells", "Naive CD8+ T cells", "Proliferating CD4+ T cells",
  "Proliferating CD8+ T cells", "Tregs", "NKT cells", "gamma-delta T cells"
)

# Solo pazienti analizzabili in P2 (bone AB)
PATIENTS_AB <- c("Bo", "Me")

# Retained for TCR × functional state integration
FUNCTIONAL_STATE_MAP <- list(
  "Naive-like"  = c("Naive CD4+ T cells","Naive CD8+ T cells"),
  "Memory-like" = c("Memory T cells","Th1 cells","Th2 cells","Th17 cells","Tfh cells"),
  "Effector"    = c("Effector CD4+ T cells","Cytotoxic CD8+ T cells"),
  "Regulatory"  = c("Tregs")
)
FUNCTIONAL_ORDER <- c("Naive-like","Memory-like","Effector","Regulatory")
STATE_PALETTE    <- c("Naive-like"="#4DBBD5","Memory-like"="#00A087","Effector"="#E64B35",
                      "Regulatory"="#F39B7F")

# ── Firme geniche ────────────────────────────────────────────
# Exhaustion = unione di Exhaustion + Tex_Terminal
# Aggiunte Naive_like e Regulatory
# Proliferazione calcolata separatamente per tipo cellulare
SIGNATURES <- list(
  Effector        = c("GZMB","PRF1","NKG7","GNLY","GZMA","GZMK","FGFBP2","CX3CR1"),
  Memory_Stemness = c("TCF7","CCR7","SELL","IL7R","LEF1","KLF2","BCL2","FOXO1"),
  Naive_like      = c("CCR7","SELL","IL7R","TCF7","LEF1","KLF2","S1PR1"),
  Exhaustion      = c("PDCD1","LAG3","HAVCR2","TIGIT","TOX","TOX2","ENTPD1",
                      "CTLA4","BATF","CD160","PRDM1","ZEB2"),
  Activation      = c("CD69","CD44","TNFRSF9","IL2RA","ICOS","CD38"),
  Tpex_StemLike   = c("TCF7","CXCR5","TOX","BCL6","SLAMF6","ID3"),
  Regulatory      = c("FOXP3","IL2RA","CTLA4","IKZF2","TNFRSF18","TIGIT","LAYN")
)

PROLIFERATION_GENES <- c("MKI67","TOP2A","PCNA","CCNB1","STMN1","UBE2C")

# Mappa campioni VDJ I-stage (per confronto longitudinale)
# NOTA: Bo (S429_I) ha il file con trattino: filtered_contig_annotations-.csv
VDJ_I_MAP <- list(
  Bo = list(folder = "2/S429_I", patient = "Bo"),
  Ca = list(folder = "1/S151_I", patient = "Ca"),
  Me = list(folder = "1/S345_I", patient = "Me")
)

# Mappa campioni VDJ B-stage (midollo AB)
VDJ_B_MAP <- list(
  Bo = list(folder = "2/S435_B", patient = "Bo"),
  Ca = list(folder = "1/S393_B", patient = "Ca"),
  Me = list(folder = "2/S431_B", patient = "Me")
)

# ======================================================================
# HELPERS
# ======================================================================
map_fs <- function(cell_types) {
  state <- rep("Other", length(cell_types))
  for (nm in names(FUNCTIONAL_STATE_MAP))
    state[cell_types %in% FUNCTIONAL_STATE_MAP[[nm]]] <- nm
  state
}

get_car_status <- function(meta) {
  for (col in c("IS_CAR_ALLIN_scREP","IS_CAR","CAR")) {
    if (col %in% colnames(meta)) {
      vals <- as.character(meta[[col]])
      return(ifelse(grepl("^(YES|TRUE|yes|true|1)$", vals), "CAR+", "CAR-"))
    }
  }
  rep("CAR-", nrow(meta))
}

# Legge file VDJ — gestisce nome con trattino (Bo I-stage) e separatore ; o ,
read_vdj <- function(folder) {
  f1 <- file.path(VDJ_BASE, folder, "vdj_t", "filtered_contig_annotations.csv")
  f2 <- file.path(VDJ_BASE, folder, "vdj_t", "filtered_contig_annotations-.csv")
  f  <- if (file.exists(f1)) f1 else if (file.exists(f2)) f2 else NULL
  if (is.null(f)) { warning("VDJ file non trovato: ", folder); return(NULL) }

  first_line <- readLines(f, n = 1)
  delim <- if (grepl(";", first_line)) ";" else ","
  df <- read_delim(f, delim = delim, show_col_types = FALSE)
  colnames(df) <- make.names(colnames(df))

  df %>%
    filter(
      is_cell        %in% c("true","True","TRUE",TRUE),
      high_confidence %in% c("true","True","TRUE",TRUE),
      tolower(as.character(productive)) == "true",
      full_length    %in% c("true","True","TRUE",TRUE),
      chain          %in% c("TRA","TRB"),
      !is.na(cdr3), nchar(as.character(cdr3)) > 0, cdr3 != "None"
    )
}

# Best contig per barcode per catena (max UMI)
best_contig <- function(d) {
  d %>%
    group_by(barcode) %>%
    arrange(desc(umis), desc(reads)) %>%
    slice(1) %>%
    ungroup()
}

# Pairing TRA+TRB e pulizia barcode
make_paired_vdj <- function(raw_df, patient_id) {
  tra <- raw_df %>% filter(chain == "TRA") %>% best_contig() %>%
    select(barcode, TRA_v = v_gene, TRA_cdr3 = cdr3)
  trb <- raw_df %>% filter(chain == "TRB") %>% best_contig() %>%
    select(barcode, TRB_v = v_gene, TRB_cdr3 = cdr3)

  inner_join(tra, trb, by = "barcode") %>%
    mutate(
      patient        = patient_id,
      Clone_ID_CDR3  = paste(TRA_cdr3, TRB_cdr3, sep = "_"),
      clean_barcode  = gsub("-[0-9]+$","", barcode)
    )
}

# Calcola metriche diversità da vettore di clonotipi
diversity_metrics <- function(clone_vec) {
  clone_vec <- na.omit(clone_vec[clone_vec != "" & clone_vec != "NA"])
  n_total <- length(clone_vec)
  if (n_total == 0)
    return(list(n_cells=0, n_clones=0, shannon=NA, simpson=NA,
                clonality=NA, pct_top1=NA))
  tbl  <- table(clone_vec)
  freq <- as.numeric(tbl) / n_total
  H     <- -sum(freq * log(freq + 1e-12))
  H_max <- log(length(tbl))
  Si    <- 1 - sum(freq^2)
  Cl    <- if (H_max > 0) 1 - H / H_max else NA
  top1  <- max(tbl) / n_total
  list(n_cells=n_total, n_clones=length(tbl),
       shannon=round(H,3), simpson=round(Si,3),
       clonality=round(Cl,3), pct_top1=round(top1,3))
}

# ======================================================================
# LOAD SEURAT AB
# ======================================================================
cat("Caricamento Seurat AB (solo midollo)...\n")
AB_samples_all <- readRDS(RDS_AB)

# Filtra solo campioni midollo osseo
bone_AB_names <- grep("bone_AB", names(AB_samples_all), value = TRUE)
cat("Campioni midollo AB:", paste(bone_AB_names, collapse = ", "), "\n")

AB_bone <- AB_samples_all[bone_AB_names]
rm(AB_samples_all); invisible(gc())

# Aggiungi metadati e conta CAR+
for (sname in names(AB_bone)) {
  meta <- AB_bone[[sname]]@meta.data
  AB_bone[[sname]]$car_status <- get_car_status(meta)
  pt <- sub("_bone_AB$","", sname)
  AB_bone[[sname]]$patient    <- pt
  AB_bone[[sname]]$fs         <- map_fs(as.character(meta$cell_type))
}

# Riepilogo
summary_tbl <- bind_rows(lapply(names(AB_bone), function(s) {
  meta <- AB_bone[[s]]@meta.data
  n_car <- sum(meta$car_status == "CAR+")
  data.frame(
    Campione  = s,
    Paziente  = unique(meta$patient),
    Totale    = nrow(meta),
    CAR_pos   = n_car,
    CAR_neg   = sum(meta$car_status == "CAR-"),
    Perc_CAR  = round(100 * n_car / nrow(meta), 1),
    Incluso   = ifelse(n_car >= 10, "SI", "NO (escluso)")
  )
}))

# ======================================================================
# OVERVIEW TABLE
# ======================================================================
kable(summary_tbl,
      caption = "Overview of bone marrow AB samples: CAR+ and CAR- cells",
      align = "c")

# ======================================================================
# FILTER AB SAMPLES
# ======================================================================
# Mantieni solo campioni con CAR+ sufficienti
AB_valid <- AB_bone[sapply(AB_bone, function(s) {
  sum(s@meta.data$car_status == "CAR+") >= 10
})]

cat("Campioni validi con CAR+ >= 10:", paste(names(AB_valid), collapse=", "), "\n")

# ======================================================================
# COMPUTE MODULE SCORES AB
# ======================================================================
for (sname in names(AB_valid)) {
  for (sig in names(SIGNATURES)) {
    genes_ok <- intersect(SIGNATURES[[sig]], rownames(AB_valid[[sname]]))
    if (length(genes_ok) < 2) next
    sc <- paste0("score_", sig)
    AB_valid[[sname]] <- AddModuleScore(
      AB_valid[[sname]], features = list(genes_ok), name = sc, seed = 42
    )
    AB_valid[[sname]]@meta.data[[sc]] <-
      AB_valid[[sname]]@meta.data[[paste0(sc,"1")]]
    AB_valid[[sname]]@meta.data[[paste0(sc,"1")]] <- NULL
  }
}

# ======================================================================
# MODULE SCORES BY CELLTYPE AB
# ======================================================================
score_cols <- paste0("score_", names(SIGNATURES))

meta_all <- bind_rows(lapply(names(AB_valid), function(s) {
  AB_valid[[s]]@meta.data
}))

sc_avail <- intersect(score_cols, colnames(meta_all))

# Aggrega per paziente × tipo cellulare × CAR status (solo T cells)
heat_ct <- meta_all %>%
  filter(!is.na(cell_type), cell_type %in% T_CELL_TYPES) %>%
  group_by(patient, car_status, cell_type) %>%
  summarise(across(all_of(sc_avail), \(x) mean(x, na.rm = TRUE)),
            n_cells = n(), .groups = "drop") %>%
  filter(n_cells >= 5) %>%
  pivot_longer(all_of(sc_avail), names_to = "signature", values_to = "score") %>%
  mutate(signature = gsub("score_","", signature),
         signature = factor(signature, levels = names(SIGNATURES)))

ggplot(heat_ct, aes(x = patient, y = cell_type, fill = score)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_gradient2(low = "#4DBBD5", mid = "white", high = "#E64B35",
                       midpoint = 0, name = "Mean\nscore") +
  facet_grid(car_status ~ signature) +
  labs(title = "Module scores per cell type in bone marrow AB",
       subtitle = "CD4/CD8 agnostic | Separated by CAR status",
       x = NULL, y = NULL) +
  theme_classic(base_size = 9) +
  theme(axis.text.x  = element_text(angle = 25, hjust = 1),
        axis.text.y  = element_text(size = 7),
        strip.background = element_rect(fill = "#F0F0F0"),
        strip.text   = element_text(face = "bold", size = 8))

# ======================================================================
# MODULE SCORES HEATMAP AB
# ======================================================================
meta_all_carpos <- bind_rows(lapply(names(AB_valid), function(s) {
  m <- AB_valid[[s]]@meta.data
  m[m$car_status == "CAR+", ]
}))

heat_df <- meta_all_carpos %>%
  group_by(patient) %>%
  summarise(across(all_of(intersect(score_cols, colnames(meta_all_carpos))),
                   \(x) mean(x, na.rm = TRUE)),
            .groups = "drop") %>%
  pivot_longer(-patient, names_to = "signature", values_to = "score") %>%
  mutate(signature = gsub("score_","", signature),
         signature = factor(signature, levels = names(SIGNATURES)))

ggplot(heat_df, aes(x = patient, y = signature, fill = score)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.3f", score)), size = 3.5) +
  scale_fill_gradient2(low = "#4DBBD5", mid = "white", high = "#E64B35",
                       midpoint = 0, name = "Mean\nscore") +
  scale_x_discrete(labels = c(Bo = "Bo (expansion)", Me = "Me (partial)")) +
  labs(title = "Module scores CAR+ in bone marrow AB",
       subtitle = "Bo vs Me comparison (Ca excluded — 0 CAR+)",
       x = NULL, y = NULL) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

# ======================================================================
# MODULE SCORES VIOLIN AB
# ======================================================================
sc_avail_ab <- intersect(score_cols, colnames(meta_all))

vln_plots <- lapply(sc_avail_ab, function(sc) {
  sig_name <- gsub("score_","", sc)
  df <- meta_all %>%
    filter(fs != "Other") %>%
    select(patient, car_status, score = !!rlang::sym(sc)) %>%
    filter(!is.na(score))

  ggplot(df, aes(x = car_status, y = score, fill = car_status)) +
    geom_violin(trim = TRUE, alpha = 0.85, color = "white") +
    geom_boxplot(width = 0.12, fill = "white", alpha = 0.7, outlier.shape = NA) +
    scale_fill_manual(values = c("CAR+" = "#E64B35","CAR-" = "#8FBCDB"), guide = "none") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.4) +
    facet_wrap(~ patient, ncol = 3,
               labeller = labeller(patient = c(Bo="Bo",Me="Me"))) +
    labs(title = sig_name, x = NULL, y = "Score") +
    theme_classic(base_size = 9) +
    theme(plot.title = element_text(face="bold", size=9),
          strip.background = element_rect(fill = "#F0F0F0"))
})

wrap_plots(vln_plots, ncol = 2) +
  plot_annotation(
    title    = "Module scores: CAR+ vs CAR- in bone marrow AB",
    subtitle = "All T cells, CD4/CD8 agnostic",
    theme    = theme(plot.title = element_text(face="bold"))
  )

# ======================================================================
# PROLIFERATION BY CELLTYPE AB
# ======================================================================
# Calcola score di proliferazione (su tutte le cellule)
for (sname in names(AB_valid)) {
  genes_ok <- intersect(PROLIFERATION_GENES, rownames(AB_valid[[sname]]))
  if (length(genes_ok) < 2) next
  AB_valid[[sname]] <- AddModuleScore(
    AB_valid[[sname]], features = list(genes_ok), name = "score_Prolif", seed = 42
  )
  AB_valid[[sname]]@meta.data[["score_Prolif"]] <-
    AB_valid[[sname]]@meta.data[["score_Prolif1"]]
  AB_valid[[sname]]@meta.data[["score_Prolif1"]] <- NULL
}

meta_prolif_ab <- bind_rows(lapply(names(AB_valid), function(s) {
  AB_valid[[s]]@meta.data
}))

prolif_by_ct_ab <- meta_prolif_ab %>%
  filter(!is.na(cell_type), !is.na(score_Prolif), cell_type %in% T_CELL_TYPES) %>%
  group_by(patient, cell_type, car_status) %>%
  summarise(mean_prolif = mean(score_Prolif, na.rm = TRUE),
            n = n(), .groups = "drop") %>%
  filter(n >= 5)

# Ordina i tipi cellulari per proliferazione media complessiva
ct_order_ab <- prolif_by_ct_ab %>%
  group_by(cell_type) %>%
  summarise(overall = mean(mean_prolif), .groups = "drop") %>%
  arrange(overall) %>%
  pull(cell_type)

prolif_by_ct_ab <- prolif_by_ct_ab %>%
  mutate(cell_type = factor(cell_type, levels = ct_order_ab))

ggplot(prolif_by_ct_ab, aes(x = cell_type, y = mean_prolif, fill = patient)) +
  geom_col(position = position_dodge(0.8), width = 0.75,
           color = "white", linewidth = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  scale_fill_manual(values = PATIENT_COLORS[PATIENTS_AB],
                    labels = OUTCOME_LABEL[PATIENTS_AB], name = NULL) +
  coord_flip() +
  facet_wrap(~ car_status, ncol = 2) +
  labs(title = "Proliferation score per cell type — Bone marrow AB",
       subtitle = paste0("Genes: ", paste(PROLIFERATION_GENES, collapse = ", ")),
       x = NULL, y = "Mean proliferation score") +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "#F0F0F0"),
        strip.text = element_text(face = "bold"))

# ======================================================================
# PROLIF DOTPLOT AB
# ======================================================================
prolif_dot_ab <- meta_prolif_ab %>%
  filter(!is.na(cell_type), !is.na(score_Prolif), car_status == "CAR+", cell_type %in% T_CELL_TYPES) %>%
  group_by(patient, cell_type) %>%
  summarise(mean_prolif = mean(score_Prolif, na.rm = TRUE),
            n = n(), .groups = "drop") %>%
  filter(n >= 3) %>%
  mutate(cell_type = factor(cell_type, levels = ct_order_ab))

ggplot(prolif_dot_ab, aes(x = patient, y = cell_type,
                           size = n, color = mean_prolif)) +
  geom_point(alpha = 0.85) +
  scale_color_gradient2(low = "#4DBBD5", mid = "white", high = "#E64B35",
                        midpoint = 0, name = "Proliferation\nscore") +
  scale_size_continuous(name = "N cells", range = c(2, 10)) +
  scale_x_discrete(labels = OUTCOME_LABEL[PATIENTS_AB]) +
  labs(title = "Proliferation per cell type — CAR+ in bone marrow AB",
       x = NULL, y = NULL) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

# ======================================================================
# LOAD VDJ B
# ======================================================================
vdj_B_all <- list()

for (pt in names(VDJ_B_MAP)) {
  info <- VDJ_B_MAP[[pt]]
  df   <- read_vdj(info$folder)
  if (is.null(df) || nrow(df) == 0) {
    cat(pt, ": VDJ B-stage non caricato\n")
    next
  }

  paired <- make_paired_vdj(df, pt)
  vdj_B_all[[pt]] <- paired
  cat(sprintf("%s (B-stage): %d cellule paired TRA+TRB\n", pt, nrow(paired)))
}

vdj_B_df <- bind_rows(vdj_B_all)

# ======================================================================
# DIVERSITY METRICS B
# ======================================================================
div_results <- list()
for (pt in names(vdj_B_all)) {
  sname <- paste0(pt, "_bone_AB")
  if (!sname %in% names(AB_valid)) next

  meta  <- AB_valid[[sname]]@meta.data
  meta$barcode_clean <- gsub("-[0-9]+$","", sub("^([^_]+_){2}","", rownames(meta)))

  vdj_pt <- vdj_B_all[[pt]] %>%
    left_join(meta %>% select(barcode_clean, car_status),
              by = c("clean_barcode" = "barcode_clean"))

  for (cst in c("CAR+","CAR-")) {
    clones_sub <- vdj_pt$Clone_ID_CDR3[!is.na(vdj_pt$car_status) &
                                          vdj_pt$car_status == cst]
    dm <- diversity_metrics(clones_sub)
    div_results[[paste0(pt,"_",cst)]] <- data.frame(
      patient    = pt,
      car_status = cst,
      n_cells    = dm$n_cells,
      n_clones   = dm$n_clones,
      shannon    = dm$shannon,
      clonality  = dm$clonality,
      simpson    = dm$simpson,
      pct_top1   = round(100 * dm$pct_top1, 1)
    )
  }
}

div_df <- bind_rows(div_results)

p_div <- ggplot(div_df, aes(x = patient, y = clonality, fill = car_status)) +
  geom_col(position = "dodge", width = 0.6, color = "white") +
  scale_fill_manual(values = c("CAR+" = "#E64B35","CAR-" = "#8FBCDB"),
                    name = "CAR status") +
  scale_x_discrete(labels = c(Bo="Bo (expansion)", Me="Me (partial)")) +
  labs(title = "Clonality index in bone marrow AB",
       subtitle = "Clonality = 1 - H_normalized | 0 = max diversity, 1 = clonal monopoly",
       x = NULL, y = "Clonality") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

p_n <- ggplot(div_df, aes(x = patient, y = n_clones, fill = car_status)) +
  geom_col(position = "dodge", width = 0.6, color = "white") +
  geom_text(aes(label = n_clones), position = position_dodge(0.6),
            vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = c("CAR+" = "#E64B35","CAR-" = "#8FBCDB"), guide = "none") +
  labs(title = "Unique clones in bone marrow AB",
       x = NULL, y = "No. clones (TRA+TRB paired)") +
  scale_x_discrete(labels = c(Bo="Bo (expansion)", Me="Me (partial)")) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

p_div | p_n

# ======================================================================
# DIVERSITY TABLE B
# ======================================================================
kable(div_df,
      caption = "TCR diversity metrics in bone marrow AB",
      digits = 3, align = "c")

# ======================================================================
# RANK FREQ PLOT B
# ======================================================================
rank_df <- bind_rows(lapply(names(vdj_B_all), function(pt) {
  sname <- paste0(pt, "_bone_AB")
  if (!sname %in% names(AB_valid)) return(NULL)

  meta  <- AB_valid[[sname]]@meta.data
  meta$bc_clean <- gsub("-[0-9]+$","", sub("^([^_]+_){2}","", rownames(meta)))

  vdj_pt <- vdj_B_all[[pt]] %>%
    left_join(meta %>% select(bc_clean, car_status),
              by = c("clean_barcode" = "bc_clean")) %>%
    filter(!is.na(car_status), car_status == "CAR+")

  if (nrow(vdj_pt) == 0) return(NULL)

  tbl <- sort(table(vdj_pt$Clone_ID_CDR3), decreasing = TRUE)
  data.frame(
    patient  = pt,
    rank     = seq_along(tbl),
    n_cells  = as.numeric(tbl),
    freq     = as.numeric(tbl) / sum(tbl),
    clone_id = names(tbl)
  )
}))

if (!is.null(rank_df) && nrow(rank_df) > 0) {
  ggplot(rank_df, aes(x = rank, y = freq, color = patient)) +
    geom_line(linewidth = 1) +
    geom_point(data = rank_df %>% filter(rank <= 5), size = 2) +
    scale_color_manual(values = PATIENT_COLORS[PATIENTS_AB],
                       labels = OUTCOME_LABEL[PATIENTS_AB], name = NULL) +
    scale_y_log10(labels = percent_format()) +
    scale_x_log10() +
    labs(title = "Clonal frequency distribution — CAR+ in bone marrow AB",
         subtitle = "Steep curve → oligoclonal | flat curve → polyclonal",
         x = "Clone rank (1 = most frequent)", y = "Relative frequency") +
    theme_classic(base_size = 12) +
    theme(legend.position = "bottom")
}

# ======================================================================
# TOP CLONES TABLE B
# ======================================================================
if (!is.null(rank_df) && nrow(rank_df) > 0) {
  top_clones <- rank_df %>%
    filter(rank <= 10) %>%
    mutate(
      CDR3_TRA = gsub("_.*","", clone_id),
      CDR3_TRB = gsub(".*_","", clone_id),
      freq_pct = sprintf("%.1f%%", 100 * freq)
    ) %>%
    select(Paziente=patient, Rank=rank, CDR3_TRA, CDR3_TRB,
           N_cellule=n_cells, Frequenza=freq_pct)

  DT::datatable(top_clones,
                caption = "Top 10 clones (CAR+) per patient in bone marrow AB",
                filter = "top",
                options = list(pageLength = 15, scrollX = TRUE))
}

# ======================================================================
# LOAD VDJ I FOR COMPARISON
# ======================================================================
vdj_I_all <- list()

for (pt in PATIENTS_AB) {  # Solo Bo e Me (Ca non ha CAR+ in AB)
  info <- VDJ_I_MAP[[pt]]
  df   <- read_vdj(info$folder)
  if (is.null(df) || nrow(df) == 0) {
    cat(pt, ": VDJ I-stage non caricato\n")
    next
  }
  paired <- make_paired_vdj(df, pt)
  vdj_I_all[[pt]] <- paired
  cat(sprintf("%s (I-stage): %d cellule paired TRA+TRB\n", pt, nrow(paired)))
}

# ======================================================================
# CLONAL FATE ANALYSIS
# ======================================================================
fate_results <- list()

for (pt in PATIENTS_AB) {
  if (!pt %in% names(vdj_I_all) || !pt %in% names(vdj_B_all)) next

  sname_AB <- paste0(pt, "_bone_AB")

  clones_I  <- unique(vdj_I_all[[pt]]$Clone_ID_CDR3)

  # Cloni CAR+ in AB
  vdj_ab_pt <- vdj_B_all[[pt]]
  if (sname_AB %in% names(AB_valid)) {
    meta_ab <- AB_valid[[sname_AB]]@meta.data
    meta_ab$bc_clean <- gsub("-[0-9]+$","", sub("^([^_]+_){2}","", rownames(meta_ab)))
    vdj_ab_pt <- vdj_ab_pt %>%
      left_join(meta_ab %>% select(bc_clean, car_status),
                by = c("clean_barcode" = "bc_clean"))
    clones_AB_carpos <- unique(vdj_ab_pt$Clone_ID_CDR3[
      !is.na(vdj_ab_pt$car_status) & vdj_ab_pt$car_status == "CAR+"])
    clones_AB_all    <- unique(vdj_ab_pt$Clone_ID_CDR3)
  } else {
    clones_AB_carpos <- character(0)
    clones_AB_all    <- unique(vdj_ab_pt$Clone_ID_CDR3)
  }

  # Fate classification
  n_survived     <- length(intersect(clones_I, clones_AB_carpos))
  n_denovo_carpos<- length(setdiff(clones_AB_carpos, clones_I))
  n_lost         <- length(setdiff(clones_I, clones_AB_all))

  fate_results[[pt]] <- data.frame(
    patient            = pt,
    n_clones_I_all     = length(clones_I),
    n_clones_AB_carpos = length(clones_AB_carpos),
    n_survived         = n_survived,
    n_denovo_carpos    = n_denovo_carpos,
    n_lost_from_I      = n_lost,
    pct_survived       = round(100 * n_survived / max(1, length(clones_AB_carpos)), 1),
    pct_denovo         = round(100 * n_denovo_carpos / max(1, length(clones_AB_carpos)), 1)
  )

  cat(sprintf(
    "\n%s: Cloni I=%d, CAR+ in AB=%d\n  Survived I→AB: %d (%.1f%%)\n  De novo in AB: %d (%.1f%%)\n",
    pt,
    length(clones_I), length(clones_AB_carpos),
    n_survived, 100 * n_survived / max(1, length(clones_AB_carpos)),
    n_denovo_carpos, 100 * n_denovo_carpos / max(1, length(clones_AB_carpos))
  ))
}

fate_df <- bind_rows(fate_results)
kable(fate_df,
      caption = "Fate clonale: proporzione di cloni I-stage ritrovati nel midollo AB",
      align = "c")

# ======================================================================
# FATE PLOT
# ======================================================================
if (nrow(fate_df) > 0) {
  fate_long <- fate_df %>%
    select(patient, Survived = n_survived, De_novo = n_denovo_carpos) %>%
    pivot_longer(-patient, names_to = "fate", values_to = "n_clones") %>%
    mutate(fate = factor(fate, levels = c("Survived","De_novo"),
                         labels = c("Survived from I","De novo in AB")))

  ggplot(fate_long, aes(x = patient, y = n_clones, fill = fate)) +
    geom_col(position = "dodge", width = 0.6, color = "white") +
    geom_text(aes(label = n_clones),
              position = position_dodge(0.6), vjust = -0.3, size = 4) +
    scale_fill_manual(
      values = c("Survived from I" = "#E64B35","De novo in AB" = "#8FBCDB"),
      name = NULL
    ) +
    scale_x_discrete(labels = c(Bo="Bo (expansion)", Me="Me (partial)")) +
    labs(title = "Clonal fate: origin of CAR+ clones in bone marrow AB",
         subtitle = "CDR3-matched clones from I-stage vs clones not detected in I",
         x = NULL, y = "No. clones") +
    theme_classic(base_size = 12) +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 15, hjust = 1))
}

# ======================================================================
# SURVIVED RANK PLOT
# ======================================================================
if (!is.null(rank_df) && nrow(rank_df) > 0 && length(vdj_I_all) > 0) {

  rank_with_fate <- bind_rows(lapply(PATIENTS_AB, function(pt) {
    if (!pt %in% names(vdj_I_all) || !pt %in% names(vdj_B_all)) return(NULL)

    clones_I_all <- unique(vdj_I_all[[pt]]$Clone_ID_CDR3)

    rank_df %>%
      filter(patient == pt) %>%
      mutate(fate = ifelse(clone_id %in% clones_I_all,
                           "Survived from I", "De novo in AB"))
  }))

  if (!is.null(rank_with_fate) && nrow(rank_with_fate) > 0) {
    ggplot(rank_with_fate, aes(x = rank, y = freq, color = fate)) +
      geom_point(alpha = 0.7, size = 1.5) +
      geom_point(data = rank_with_fate %>% filter(rank <= 3), size = 3) +
      scale_color_manual(
        values = c("Survived from I" = "#E64B35", "De novo in AB" = "#8FBCDB"),
        name = NULL
      ) +
      scale_y_log10(labels = percent_format()) +
      scale_x_log10() +
      facet_wrap(~ patient, ncol = 2,
                 labeller = labeller(patient = OUTCOME_LABEL)) +
      labs(title = "Rank-frequency: survived vs de novo clones (CAR+ in bone marrow AB)",
           subtitle = "Red = clone also detected in infusion product I | Blue = not detected in I",
           x = "Rank (1 = most frequent in AB)", y = "Relative frequency") +
      theme_classic(base_size = 12) +
      theme(legend.position = "bottom",
            strip.background = element_rect(fill = "#F0F0F0"),
            strip.text = element_text(face = "bold"))
  }
}

# ======================================================================
# SURVIVED CLONE SEQUENCES
# ======================================================================
survived_detail <- list()

for (pt in PATIENTS_AB) {
  if (!pt %in% names(vdj_I_all) || !pt %in% names(vdj_B_all)) next

  sname_AB <- paste0(pt, "_bone_AB")
  if (!sname_AB %in% names(AB_valid)) next

  meta_ab <- AB_valid[[sname_AB]]@meta.data
  meta_ab$bc_clean <- gsub("-[0-9]+$","", sub("^([^_]+_){2}","", rownames(meta_ab)))

  # Cloni CAR+ in AB con frequenza
  vdj_ab_car <- vdj_B_all[[pt]] %>%
    left_join(meta_ab %>% select(bc_clean, car_status),
              by = c("clean_barcode" = "bc_clean")) %>%
    filter(!is.na(car_status), car_status == "CAR+")

  freq_ab <- vdj_ab_car %>%
    count(Clone_ID_CDR3, TRA_cdr3, TRB_cdr3, TRA_v, TRB_v, name = "n_AB") %>%
    mutate(freq_AB = round(100 * n_AB / sum(n_AB), 2)) %>%
    arrange(desc(n_AB)) %>%
    mutate(rank_AB = row_number())

  # Cloni presenti in I (con frequenza in I)
  freq_I <- vdj_I_all[[pt]] %>%
    count(Clone_ID_CDR3, name = "n_I") %>%
    mutate(freq_I = round(100 * n_I / sum(n_I), 2))

  # Sopravvissuti = in entrambi
  survived_pt <- freq_ab %>%
    inner_join(freq_I, by = "Clone_ID_CDR3") %>%
    mutate(patient = pt) %>%
    select(patient, rank_AB, Clone_ID_CDR3,
           TRA_V = TRA_v, CDR3_alpha = TRA_cdr3,
           TRB_V = TRB_v, CDR3_beta  = TRB_cdr3,
           N_in_I = n_I, Freq_I_pct = freq_I,
           N_in_AB = n_AB, Freq_AB_pct = freq_AB)

  survived_detail[[pt]] <- survived_pt
  cat(sprintf("%s: %d cloni sopravvissuti da I a AB (su %d CAR+ in AB)\n",
              pt, nrow(survived_pt), nrow(freq_ab)))
}

survived_all <- bind_rows(survived_detail)

if (nrow(survived_all) > 0) {
  DT::datatable(
    survived_all,
    caption  = "Clonotypes survived from I to bone AB — full TRA+TRB sequences",
    filter   = "top",
    rownames = FALSE,
    options  = list(pageLength = 20, scrollX = TRUE)
  ) %>%
    DT::formatStyle("Freq_AB_pct",
                    background = DT::styleColorBar(range(survived_all$Freq_AB_pct), "#E64B3540"),
                    backgroundSize = "100% 80%",
                    backgroundRepeat = "no-repeat",
                    backgroundPosition = "center")
}

# ======================================================================
# SURVIVED FREQ PLOT
# ======================================================================
if (nrow(survived_all) > 0) {
  top_survived <- survived_all %>%
    group_by(patient) %>%
    slice_max(Freq_AB_pct, n = 15) %>%
    ungroup() %>%
    mutate(clone_label = paste0(substr(CDR3_alpha, 1, 10), "…/",
                                 substr(CDR3_beta,  1, 10), "…"),
           clone_label = reorder(clone_label, Freq_AB_pct))

  ggplot(top_survived, aes(x = clone_label, y = Freq_AB_pct, fill = patient)) +
    geom_col(color = "white") +
    scale_fill_manual(values = PATIENT_COLORS[PATIENTS_AB], guide = "none") +
    coord_flip() +
    facet_wrap(~ patient, scales = "free_y", ncol = 2,
               labeller = labeller(patient = OUTCOME_LABEL)) +
    labs(title = "Top survived clones (I → bone AB): frequency in bone marrow AB",
         subtitle = "Only clones with CDR3 TRA+TRB identified at both timepoints",
         x = NULL, y = "Frequency in bone marrow AB (%)") +
    theme_classic(base_size = 11) +
    theme(strip.background = element_rect(fill = "#F0F0F0"),
          strip.text = element_text(face = "bold"))
}

# ======================================================================
# CROSS PATIENT SHARING
# ======================================================================
if (length(survived_detail) >= 2 &&
    "Bo" %in% names(survived_detail) &&
    "Me" %in% names(survived_detail)) {

  bo_clones <- survived_detail[["Bo"]]$Clone_ID_CDR3
  me_clones <- survived_detail[["Me"]]$Clone_ID_CDR3

  shared_bo_me <- intersect(bo_clones, me_clones)

  cat(sprintf("Cloni sopravvissuti in Bo: %d\n", length(bo_clones)))
  cat(sprintf("Cloni sopravvissuti in Me: %d\n", length(me_clones)))
  cat(sprintf("Cloni condivisi Bo ∩ Me (sopravvissuti in entrambi): %d\n",
              length(shared_bo_me)))

  if (length(shared_bo_me) > 0) {
    shared_detail <- survived_all %>%
      filter(Clone_ID_CDR3 %in% shared_bo_me) %>%
      select(patient, CDR3_alpha, CDR3_beta, TRA_V, TRB_V,
             N_in_I, Freq_I_pct, N_in_AB, Freq_AB_pct)
    kable(shared_detail,
          caption = "Public clonotypes: same CDR3α+CDR3β survived in both Bo and Me",
          align = "c")
  } else {
    cat("Nessun clonotipo identico (CDR3α+CDR3β) condiviso tra Bo e Me.\n\n")
    cat("Interpretazione: le espansioni di Bo e Me sono indipendenti —\n")
    cat("i cloni che si sono espansi nei due pazienti hanno specificità TCR diverse.\n")
  }

  # Controllo: CDR3β solo (più permissivo — stessa catena β, α diversa)
  bo_trb <- survived_detail[["Bo"]]$CDR3_beta
  me_trb <- survived_detail[["Me"]]$CDR3_beta
  shared_trb_only <- intersect(bo_trb, me_trb)
  cat(sprintf("\nCDR3β condivisi (catena β sola, α può differire): %d\n",
              length(shared_trb_only)))
  if (length(shared_trb_only) > 0) {
    cat("CDR3β condivisi:", paste(head(shared_trb_only, 5), collapse=", "), "\n")
  }
}

# ======================================================================
# TCR FUNCTIONAL AB
# ======================================================================
integrated_ab_df <- bind_rows(lapply(names(AB_valid), function(sname) {
  pt <- sub("_bone_AB$","", sname)
  if (!pt %in% names(vdj_B_all)) return(NULL)

  meta  <- AB_valid[[sname]]@meta.data
  meta$bc_clean <- gsub("-[0-9]+$","", sub("^([^_]+_){2}","", rownames(meta)))

  vdj_pt <- vdj_B_all[[pt]] %>%
    left_join(meta %>% select(bc_clean, car_status, fs,
                               starts_with("score_")),
              by = c("clean_barcode" = "bc_clean"))

  # Rank clonale per CAR+ cells
  tbl_carpos <- sort(
    table(vdj_pt$Clone_ID_CDR3[!is.na(vdj_pt$car_status) &
                                  vdj_pt$car_status == "CAR+"],
          useNA = "no"),
    decreasing = TRUE
  )
  rank_map <- data.frame(
    Clone_ID_CDR3 = names(tbl_carpos),
    clone_rank    = seq_along(tbl_carpos),
    clone_freq    = as.numeric(tbl_carpos) / sum(tbl_carpos)
  )

  vdj_pt %>%
    filter(!is.na(car_status), car_status == "CAR+") %>%
    left_join(rank_map, by = "Clone_ID_CDR3") %>%
    mutate(patient = pt,
           rank_cat = case_when(
             clone_rank == 1    ~ "Clone #1 (dominant)",
             clone_rank <= 5    ~ "Clone #2-5",
             clone_rank <= 20   ~ "Clone #6-20",
             !is.na(clone_rank) ~ "Clone #21+",
             TRUE               ~ "No TCR"
           ))
}))

if (!is.null(integrated_ab_df) && nrow(integrated_ab_df) > 0) {
  plot_df <- integrated_ab_df %>%
    filter(!is.na(fs), fs != "Other") %>%
    group_by(patient, rank_cat, fs) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(patient, rank_cat) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup() %>%
    mutate(
      fs = factor(fs, levels = FUNCTIONAL_ORDER),
      rank_cat = factor(rank_cat,
                        levels = c("Clone #1 (dominant)","Clone #2-5",
                                   "Clone #6-20","Clone #21+","No TCR"))
    )

  ggplot(plot_df, aes(x = rank_cat, y = prop, fill = fs)) +
    geom_col(width = 0.75, color = "white", linewidth = 0.3) +
    facet_wrap(~ patient, ncol = 2,
               labeller = labeller(patient = OUTCOME_LABEL)) +
    scale_fill_manual(values = STATE_PALETTE, drop = FALSE,
                      name = "Functional state") +
    scale_y_continuous(labels = percent_format(), expand = c(0,0)) +
    labs(
      title    = "Functional state of CAR+ by clonal rank in bone marrow AB",
      subtitle = "What functional state are the most frequent clones (rank #1) in?",
      x        = "Clonal rank", y = "Proportion"
    ) +
    theme_classic(base_size = 11) +
    theme(axis.text.x  = element_text(angle = 30, hjust = 1),
          strip.background = element_rect(fill = "#F0F0F0"),
          strip.text = element_text(face = "bold"))
}

# ======================================================================
# CLONE MODULE SCORES
# ======================================================================
if (!is.null(integrated_ab_df) && nrow(integrated_ab_df) > 0) {
  sc_avail_int <- intersect(paste0("score_", names(SIGNATURES)),
                             colnames(integrated_ab_df))

  if (length(sc_avail_int) > 0) {
    score_by_rank <- integrated_ab_df %>%
      filter(!is.na(clone_rank)) %>%
      mutate(rank_cat = factor(rank_cat,
                               levels = c("Clone #1 (dominant)","Clone #2-5",
                                          "Clone #6-20","Clone #21+"))) %>%
      group_by(patient, rank_cat) %>%
      summarise(across(all_of(sc_avail_int), \(x) mean(x, na.rm=TRUE)),
                n = n(), .groups = "drop") %>%
      pivot_longer(all_of(sc_avail_int), names_to = "signature", values_to = "score") %>%
      mutate(signature = gsub("score_","", signature))

    ggplot(score_by_rank, aes(x = rank_cat, y = score, color = patient, group = patient)) +
      geom_line(linewidth = 1.2) +
      geom_point(size = 3) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
      scale_color_manual(values = PATIENT_COLORS[PATIENTS_AB],
                         labels = OUTCOME_LABEL[PATIENTS_AB], name = NULL) +
      facet_wrap(~ signature, scales = "free_y", ncol = 3) +
      labs(title = "Mean module scores by clonal rank — CAR+ bone marrow AB",
           subtitle = "Do dominant clones show different signatures from rare clones?",
           x = "Clonal rank", y = "Mean score") +
      theme_classic(base_size = 10) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1),
            strip.background = element_rect(fill = "#F0F0F0"),
            strip.text = element_text(face = "bold"),
            legend.position = "bottom")
  }
}

# ======================================================================
# LOAD BLOOD AB
# ======================================================================
# Pre-calcola AB_means per la sezione longitudinale (prima di liberare AB_valid)
.score_cols_tmp <- paste0("score_", names(SIGNATURES))
AB_means_precomputed <- bind_rows(lapply(names(AB_valid), function(s) {
  m  <- AB_valid[[s]]@meta.data
  pt <- unique(m$patient)
  m  <- m[m$car_status == "CAR+", ]
  if (nrow(m) == 0) return(NULL)
  sc_avail_l <- intersect(.score_cols_tmp, colnames(m))
  means <- colMeans(m[, sc_avail_l, drop = FALSE], na.rm = TRUE)
  data.frame(patient = pt, timepoint = "Bone marrow AB", t(means))
}))
rm(.score_cols_tmp)

# Pre-calcola n_carpos_ab per summary-panel (prima del rilascio memoria)
n_carpos_ab_precomputed <- bind_rows(lapply(names(AB_valid), function(s) {
  data.frame(patient  = sub("_bone_AB$","", s),
             n_carpos = sum(AB_valid[[s]]@meta.data$car_status == "CAR+"))
}))

# Libera memoria dalla sezione bone AB prima di ricaricare RDS_AB
for (.obj in c("AB_bone","AB_valid","meta_all","meta_prolif_ab","meta_all_carpos")) {
  if (exists(.obj)) { rm(list = .obj); invisible(gc()) }
}

cat("Caricamento Seurat blood AB...\n")
AB_samples_blood_raw <- readRDS(RDS_AB)

blood_AB_names <- grep("blood_AB", names(AB_samples_blood_raw), value = TRUE)
cat("Campioni blood AB trovati:", paste(blood_AB_names, collapse = ", "), "\n")

AB_blood <- AB_samples_blood_raw[blood_AB_names]
rm(AB_samples_blood_raw); invisible(gc())

# Aggiungi metadati
for (sname in names(AB_blood)) {
  meta <- AB_blood[[sname]]@meta.data
  AB_blood[[sname]]$car_status <- get_car_status(meta)
  pt <- sub("_blood_AB$","", sname)
  AB_blood[[sname]]$patient    <- pt
  AB_blood[[sname]]$fs         <- map_fs(as.character(meta$cell_type))
}

# Riepilogo
summary_blood <- bind_rows(lapply(names(AB_blood), function(s) {
  meta <- AB_blood[[s]]@meta.data
  n_car <- sum(meta$car_status == "CAR+")
  data.frame(
    Campione  = s,
    Paziente  = unique(meta$patient),
    Totale    = nrow(meta),
    CAR_pos   = n_car,
    CAR_neg   = sum(meta$car_status == "CAR-"),
    Perc_CAR  = round(100 * n_car / nrow(meta), 1),
    Incluso   = ifelse(n_car >= 10, "SI", "NO (escluso)")
  )
}))

kable(summary_blood,
      caption = "Overview of blood AB samples: CAR+ and CAR- cells",
      align = "c")

# ======================================================================
# FILTER BLOOD SAMPLES
# ======================================================================
AB_blood_valid <- AB_blood[sapply(AB_blood, function(s) {
  sum(s@meta.data$car_status == "CAR+") >= 10
})]

cat("Campioni blood validi con CAR+ >= 10:", paste(names(AB_blood_valid), collapse=", "), "\n")

# ======================================================================
# COMPUTE MODULE SCORES BLOOD
# ======================================================================
for (sname in names(AB_blood_valid)) {
  for (sig in names(SIGNATURES)) {
    genes_ok <- intersect(SIGNATURES[[sig]], rownames(AB_blood_valid[[sname]]))
    if (length(genes_ok) < 2) next
    sc <- paste0("score_", sig)
    AB_blood_valid[[sname]] <- AddModuleScore(
      AB_blood_valid[[sname]], features = list(genes_ok), name = sc, seed = 42
    )
    AB_blood_valid[[sname]]@meta.data[[sc]] <-
      AB_blood_valid[[sname]]@meta.data[[paste0(sc,"1")]]
    AB_blood_valid[[sname]]@meta.data[[paste0(sc,"1")]] <- NULL
  }
}

# ======================================================================
# MODULE SCORES BY CELLTYPE BLOOD
# ======================================================================
score_cols_blood <- paste0("score_", names(SIGNATURES))

meta_blood_all <- bind_rows(lapply(names(AB_blood_valid), function(s) {
  AB_blood_valid[[s]]@meta.data
}))

sc_avail_blood <- intersect(score_cols_blood, colnames(meta_blood_all))

# Aggrega per paziente × tipo cellulare × CAR status (solo T cells)
heat_ct_blood <- meta_blood_all %>%
  filter(!is.na(cell_type), cell_type %in% T_CELL_TYPES) %>%
  group_by(patient, car_status, cell_type) %>%
  summarise(across(all_of(sc_avail_blood), \(x) mean(x, na.rm = TRUE)),
            n_cells = n(), .groups = "drop") %>%
  filter(n_cells >= 5) %>%
  pivot_longer(all_of(sc_avail_blood), names_to = "signature", values_to = "score") %>%
  mutate(signature = gsub("score_","", signature),
         signature = factor(signature, levels = names(SIGNATURES)))

ggplot(heat_ct_blood, aes(x = patient, y = cell_type, fill = score)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_gradient2(low = "#4DBBD5", mid = "white", high = "#E64B35",
                       midpoint = 0, name = "Mean\nscore") +
  facet_grid(car_status ~ signature) +
  labs(title = "Module scores per cell type — Blood AB",
       subtitle = "CD4/CD8 agnostic | Separated by CAR status",
       x = NULL, y = NULL) +
  theme_classic(base_size = 9) +
  theme(axis.text.x  = element_text(angle = 25, hjust = 1),
        axis.text.y  = element_text(size = 7),
        strip.background = element_rect(fill = "#F0F0F0"),
        strip.text   = element_text(face = "bold", size = 8))

# ======================================================================
# MODULE SCORES HEATMAP BLOOD
# ======================================================================
meta_blood_carpos <- meta_blood_all[meta_blood_all$car_status == "CAR+", ]

if (nrow(meta_blood_carpos) > 0) {
  heat_blood <- meta_blood_carpos %>%
    group_by(patient) %>%
    summarise(across(all_of(intersect(score_cols_blood, colnames(meta_blood_carpos))),
                     \(x) mean(x, na.rm = TRUE)),
              .groups = "drop") %>%
    pivot_longer(-patient, names_to = "signature", values_to = "score") %>%
    mutate(signature = gsub("score_","", signature),
           signature = factor(signature, levels = names(SIGNATURES)))

  ggplot(heat_blood, aes(x = patient, y = signature, fill = score)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(aes(label = sprintf("%.3f", score)), size = 3.5) +
    scale_fill_gradient2(low = "#4DBBD5", mid = "white", high = "#E64B35",
                         midpoint = 0, name = "Mean\nscore") +
    scale_x_discrete(labels = OUTCOME_LABEL) +
    labs(title = "Module scores CAR+ in blood AB",
         subtitle = "Inter-patient comparison (descriptive)",
         x = NULL, y = NULL) +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))
}

# ======================================================================
# MODULE SCORES VIOLIN BLOOD
# ======================================================================
vln_plots_blood <- lapply(sc_avail_blood, function(sc) {
  sig_name <- gsub("score_","", sc)
  df <- meta_blood_all %>%
    filter(fs != "Other") %>%
    select(patient, car_status, score = !!rlang::sym(sc)) %>%
    filter(!is.na(score))

  if (nrow(df) == 0) return(NULL)

  ggplot(df, aes(x = car_status, y = score, fill = car_status)) +
    geom_violin(trim = TRUE, alpha = 0.85, color = "white") +
    geom_boxplot(width = 0.12, fill = "white", alpha = 0.7, outlier.shape = NA) +
    scale_fill_manual(values = c("CAR+" = "#E64B35","CAR-" = "#8FBCDB"), guide = "none") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.4) +
    facet_wrap(~ patient, ncol = 3) +
    labs(title = sig_name, x = NULL, y = "Score") +
    theme_classic(base_size = 9) +
    theme(plot.title = element_text(face="bold", size=9),
          strip.background = element_rect(fill = "#F0F0F0"))
})

vln_plots_blood <- Filter(Negate(is.null), vln_plots_blood)

if (length(vln_plots_blood) > 0) {
  wrap_plots(vln_plots_blood, ncol = 2) +
    plot_annotation(
      title    = "Module scores: CAR+ vs CAR- in blood AB",
      subtitle = "All T cells, CD4/CD8 agnostic",
      theme    = theme(plot.title = element_text(face="bold"))
    )
}

# ======================================================================
# PROLIFERATION BLOOD
# ======================================================================
for (sname in names(AB_blood_valid)) {
  genes_ok <- intersect(PROLIFERATION_GENES, rownames(AB_blood_valid[[sname]]))
  if (length(genes_ok) < 2) next
  AB_blood_valid[[sname]] <- AddModuleScore(
    AB_blood_valid[[sname]], features = list(genes_ok), name = "score_Prolif", seed = 42
  )
  AB_blood_valid[[sname]]@meta.data[["score_Prolif"]] <-
    AB_blood_valid[[sname]]@meta.data[["score_Prolif1"]]
  AB_blood_valid[[sname]]@meta.data[["score_Prolif1"]] <- NULL
}

meta_prolif_blood <- bind_rows(lapply(names(AB_blood_valid), function(s) {
  AB_blood_valid[[s]]@meta.data
}))

prolif_by_ct_blood <- meta_prolif_blood %>%
  filter(!is.na(cell_type), !is.na(score_Prolif), cell_type %in% T_CELL_TYPES) %>%
  group_by(patient, cell_type, car_status) %>%
  summarise(mean_prolif = mean(score_Prolif, na.rm = TRUE),
            n = n(), .groups = "drop") %>%
  filter(n >= 5)

ct_order_blood <- prolif_by_ct_blood %>%
  group_by(cell_type) %>%
  summarise(overall = mean(mean_prolif), .groups = "drop") %>%
  arrange(overall) %>%
  pull(cell_type)

prolif_by_ct_blood <- prolif_by_ct_blood %>%
  mutate(cell_type = factor(cell_type, levels = ct_order_blood))

ggplot(prolif_by_ct_blood, aes(x = cell_type, y = mean_prolif, fill = patient)) +
  geom_col(position = position_dodge(0.8), width = 0.75,
           color = "white", linewidth = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  scale_fill_manual(values = PATIENT_COLORS, labels = OUTCOME_LABEL, name = NULL) +
  coord_flip() +
  facet_wrap(~ car_status, ncol = 2) +
  labs(title = "Proliferation score per cell type — Blood AB",
       subtitle = paste0("Genes: ", paste(PROLIFERATION_GENES, collapse = ", ")),
       x = NULL, y = "Mean proliferation score") +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "#F0F0F0"),
        strip.text = element_text(face = "bold"))

# ======================================================================
# PROLIF DOTPLOT BLOOD
# ======================================================================
prolif_dot_blood <- meta_prolif_blood %>%
  filter(!is.na(cell_type), !is.na(score_Prolif), car_status == "CAR+", cell_type %in% T_CELL_TYPES) %>%
  group_by(patient, cell_type) %>%
  summarise(mean_prolif = mean(score_Prolif, na.rm = TRUE),
            n = n(), .groups = "drop") %>%
  filter(n >= 3) %>%
  mutate(cell_type = factor(cell_type, levels = ct_order_blood))

if (nrow(prolif_dot_blood) > 0) {
  ggplot(prolif_dot_blood, aes(x = patient, y = cell_type,
                                size = n, color = mean_prolif)) +
    geom_point(alpha = 0.85) +
    scale_color_gradient2(low = "#4DBBD5", mid = "white", high = "#E64B35",
                          midpoint = 0, name = "Proliferation\nscore") +
    scale_size_continuous(name = "N cells", range = c(2, 10)) +
    scale_x_discrete(labels = OUTCOME_LABEL) +
    labs(title = "Proliferation per cell type — CAR+ in blood AB",
         x = NULL, y = NULL) +
    theme_classic(base_size = 11) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))
}

# ======================================================================
# LOAD I FOR LONGITUDINAL
# ======================================================================
cat("Caricamento Seurat I per confronto longitudinale...\n")
I_samples_all <- readRDS(RDS_I)

# Carica tutti i pazienti per I — Ca ha dati I validi anche se non ha midollo AB
all_patients <- c("Bo", "Me", "Ca")
I_valid <- I_samples_all[paste0(all_patients, "_bone_I")]
I_valid <- Filter(Negate(is.null), I_valid)
rm(I_samples_all); invisible(gc())

# Aggiungi metadati
for (sname in names(I_valid)) {
  meta <- I_valid[[sname]]@meta.data
  I_valid[[sname]]$car_status <- get_car_status(meta)
  I_valid[[sname]]$patient    <- sub("_bone_I$","", sname)
  I_valid[[sname]]$fs         <- map_fs(as.character(meta$cell_type))
}

# Calcola module scores per I-stage
for (sname in names(I_valid)) {
  for (sig in names(SIGNATURES)) {
    genes_ok <- intersect(SIGNATURES[[sig]], rownames(I_valid[[sname]]))
    if (length(genes_ok) < 2) next
    sc <- paste0("score_", sig)
    I_valid[[sname]] <- AddModuleScore(
      I_valid[[sname]], features = list(genes_ok), name = sc, seed = 42
    )
    I_valid[[sname]]@meta.data[[sc]] <-
      I_valid[[sname]]@meta.data[[paste0(sc,"1")]]
    I_valid[[sname]]@meta.data[[paste0(sc,"1")]] <- NULL
  }
}

# ======================================================================
# LONGITUDINAL HEATMAP
# ======================================================================
score_cols_long <- paste0("score_", names(SIGNATURES))

# Estrai medie CAR+ per I-stage
I_means <- bind_rows(lapply(names(I_valid), function(s) {
  m  <- I_valid[[s]]@meta.data
  pt <- unique(m$patient)
  m  <- m[m$car_status == "CAR+", ]
  if (nrow(m) == 0) return(NULL)
  sc_avail_l <- intersect(score_cols_long, colnames(m))
  means <- colMeans(m[, sc_avail_l, drop=FALSE], na.rm=TRUE)
  data.frame(patient=pt, timepoint="Infusion I", t(means))
}))

# Estrai medie CAR+ per sangue AB (solo pazienti con campione valido)
blood_means <- bind_rows(lapply(names(AB_blood_valid), function(s) {
  m  <- AB_blood_valid[[s]]@meta.data
  pt <- unique(m$patient)
  m  <- m[m$car_status == "CAR+", ]
  if (nrow(m) == 0) return(NULL)
  sc_avail_l <- intersect(score_cols_long, colnames(m))
  means <- colMeans(m[, sc_avail_l, drop=FALSE], na.rm=TRUE)
  data.frame(patient=pt, timepoint="Blood AB", t(means))
}))

# Medie CAR+ midollo AB (pre-calcolate prima del rilascio memoria)
AB_means <- AB_means_precomputed

long_means <- bind_rows(I_means, blood_means, AB_means) %>%
  pivot_longer(-c(patient, timepoint), names_to = "signature", values_to = "score") %>%
  mutate(signature = gsub("score_","", signature),
         signature = factor(signature, levels = names(SIGNATURES)),
         timepoint = factor(timepoint, levels = c("Infusion I","Blood AB","Bone marrow AB")))

ggplot(long_means, aes(x = timepoint, y = signature, fill = score)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.3f", score)), size = 3) +
  scale_fill_gradient2(low = "#4DBBD5", mid = "white", high = "#E64B35",
                       midpoint = 0, name = "Mean\nscore") +
  facet_wrap(~ patient, ncol = 3,
             labeller = labeller(patient = OUTCOME_LABEL)) +
  labs(title = "Module scores CAR+: Infusion I → Blood AB → Bone marrow AB",
       subtitle = "Longitudinal comparison | Blood AB disponibile solo per Bo e Ca",
       x = NULL, y = NULL) +
  theme_classic(base_size = 11) +
  theme(strip.background = element_rect(fill = "#F0F0F0"),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 20, hjust = 1))

# ======================================================================
# DELTA SCORE PLOT
# ======================================================================
# Delta I → Blood (solo Bo e Ca)
delta_I_blood <- long_means %>%
  filter(timepoint %in% c("Infusion I","Blood AB")) %>%
  pivot_wider(names_from = timepoint, values_from = score) %>%
  filter(!is.na(`Blood AB`)) %>%
  mutate(delta = `Blood AB` - `Infusion I`,
         tratto = "I → Blood AB")

# Delta I → Bone marrow (tutti i pazienti)
delta_I_bone <- long_means %>%
  filter(timepoint %in% c("Infusion I","Bone marrow AB")) %>%
  pivot_wider(names_from = timepoint, values_from = score) %>%
  filter(!is.na(`Bone marrow AB`)) %>%
  mutate(delta = `Bone marrow AB` - `Infusion I`,
         tratto = "I → Bone marrow AB")

delta_all <- bind_rows(delta_I_blood, delta_I_bone) %>%
  mutate(tratto = factor(tratto, levels = c("I → Blood AB",
                                             "I → Bone marrow AB")))

ggplot(delta_all, aes(x = signature, y = delta, fill = delta > 0)) +
  geom_col(color = "white") +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  scale_y_continuous(expand = expansion(mult = c(0.20, 0.15))) +
  geom_text(aes(label = sprintf("%+.3f", delta),
                vjust = ifelse(delta >= 0, -0.5, 1.6)),
            size = 2.8) +
  scale_fill_manual(values = c("TRUE" = "#E64B35","FALSE" = "#4DBBD5"), guide = "none") +
  facet_grid(patient ~ tratto,
             labeller = labeller(patient = c(Bo="Bo", Ca="Ca", Me="Me"))) +
  labs(title = "Longitudinal delta-score of CAR+ functional signatures",
       subtitle = "Red = increased signature | Blue = decreased signature",
       x = NULL, y = "Δ score") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        strip.background = element_rect(fill = "#F0F0F0"),
        strip.text = element_text(face = "bold", size = 9))

# ======================================================================
# CLEANUP I
# ======================================================================
rm(I_valid); invisible(gc())

# ======================================================================
# SUMMARY PANEL
# ======================================================================
# Pannello riassuntivo: metriche chiave a confronto

# 1. N° CAR+ cells nel midollo AB (pre-calcolato prima del rilascio memoria)
n_carpos_ab <- n_carpos_ab_precomputed

# 2. Clonality index CAR+ in AB
clon_ab <- div_df %>%
  filter(car_status == "CAR+") %>%
  select(patient, clonality)

# 3. Fate clonale: % sopravvissuto
if (nrow(fate_df) > 0) {
  fate_summ <- fate_df %>% select(patient, pct_survived)
} else {
  fate_summ <- data.frame(patient = PATIENTS_AB, pct_survived = NA)
}

summary_panel <- n_carpos_ab %>%
  left_join(clon_ab, by = "patient") %>%
  left_join(fate_summ, by = "patient")

kable(summary_panel,
      col.names = c("Patient","N CAR+ in bone AB","Clonality CAR+","% clones survived from I"),
      caption = "Summary of key metrics in bone marrow AB",
      align = "c",
      digits = 3)

# ======================================================================
# SUMMARY BARPLOT
# ======================================================================
summary_long <- summary_panel %>%
  pivot_longer(-patient, names_to = "metrica", values_to = "valore") %>%
  mutate(metrica = factor(metrica,
                          levels = c("n_carpos","clonality","pct_survived"),
                          labels = c("No. CAR+ in bone AB",
                                     "Clonality TCR (CAR+)",
                                     "% clones survived from I")))

ggplot(summary_long, aes(x = patient, y = valore, fill = patient)) +
  geom_col(width = 0.6, color = "white") +
  geom_text(aes(label = round(valore, 2)), vjust = -0.3, size = 4) +
  scale_fill_manual(values = PATIENT_COLORS[PATIENTS_AB],
                    labels = OUTCOME_LABEL[PATIENTS_AB], guide = "none") +
  scale_x_discrete(labels = OUTCOME_LABEL[PATIENTS_AB]) +
  facet_wrap(~ metrica, scales = "free_y", ncol = 3) +
  labs(title = "Bo vs Me comparison: key metrics in bone marrow AB",
       x = NULL, y = NULL) +
  theme_classic(base_size = 12) +
  theme(strip.background = element_rect(fill = "#F0F0F0"),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 15, hjust = 1))

# ======================================================================
# CLEANUP FINAL
# ======================================================================
rm(AB_blood_valid, AB_blood); invisible(gc())
