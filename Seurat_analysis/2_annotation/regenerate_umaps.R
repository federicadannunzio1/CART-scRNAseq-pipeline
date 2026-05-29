#!/usr/bin/env Rscript
# Regenerate UMAP PNGs without Italian panel subtitles ("Con label"/"Senza label").
# Loads pre-computed annotated Seurat objects; only reruns plotting + ggsave.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

base_dir <- file.path(
  path.expand("~"),
  "federica.dannunzio@uniroma1.it - Google Drive",
  "Drive condivisi", "caruana-project", "CART",
  "Code", "Seurat_analysis", "2_annotation"
)

out_I  <- file.path(base_dir, "Pipeline_I",  "Annotation_UMAP")
out_AB <- file.path(base_dir, "AB_annotation")
dir.create(out_I,  showWarnings = FALSE, recursive = TRUE)
dir.create(out_AB, showWarnings = FALSE, recursive = TRUE)

# ── Palette and canonical order (identical to original scripts) ───────────────
PALETTE <- c(
  "Naive CD4+ T cells"         = "#E63946",
  "Th1 cells"                  = "#C1121F",
  "Th2 cells"                  = "#FF99C8",
  "Th17 cells"                 = "#FB5607",
  "Tfh cells"                  = "#800F2F",
  "Effector CD4+ T cells"      = "#F4A261",
  "Memory T cells"             = "#2A9D8F",
  "Tregs"                      = "#E9C46A",
  "Cytotoxic CD8+ T cells"     = "#264653",
  "Naive CD8+ T cells"         = "#577590",
  "Proliferating T cells"      = "#457B9D",
  "Proliferating CD4+ T cells" = "#023E8A",
  "Proliferating CD8+ T cells" = "#6A0572",
  "NK cells"                   = "#43AA8B",
  "NKT cells"                  = "#277DA1",
  "gamma-delta T cells"        = "#4D908E",
  "MAIT cells"                 = "#F3722C",
  "ILC"                        = "#F9C74F",
  "B cells"                    = "#90BE6D",
  "Memory B cells"             = "#52B788",
  "Plasma cells"               = "#C77DFF",
  "CD14 Monocytes"             = "#9D0208",
  "CD16 Monocytes"             = "#DC2F02",
  "Myeloid cells"              = "#E76F51",
  "Basophils"                  = "#6A4C93",
  "HSPC"                       = "#B5838D",
  "Erythroid cells"            = "#FFAFCC",
  "Platelets"                  = "#CDB4DB"
)

CANONICAL_ORDER <- c(
  "Naive CD4+ T cells","Th1 cells","Th2 cells","Th17 cells","Tfh cells",
  "Effector CD4+ T cells","Memory T cells","Tregs",
  "Cytotoxic CD8+ T cells","Naive CD8+ T cells",
  "Proliferating T cells","Proliferating CD4+ T cells","Proliferating CD8+ T cells",
  "NKT cells","gamma-delta T cells","MAIT cells",
  "NK cells","ILC",
  "B cells","Memory B cells","Plasma cells",
  "CD14 Monocytes","CD16 Monocytes","Myeloid cells","Basophils",
  "HSPC","Erythroid cells","Platelets"
)

get_colors <- function(present_types) {
  cols    <- PALETTE[present_types]
  missing <- present_types[is.na(cols)]
  if (length(missing) > 0) {
    extra <- setNames(hue_pal()(length(missing)), missing)
    cols[missing] <- extra
    PALETTE[missing] <<- extra
  }
  cols
}

make_full_legend <- function(present_types) {
  all_t   <- c(CANONICAL_ORDER, setdiff(names(PALETTE), CANONICAL_ORDER))
  all_t   <- all_t[all_t %in% names(PALETTE)]
  n       <- length(all_t)
  present <- all_t %in% present_types
  df <- data.frame(
    y       = rev(seq_len(n)),
    label   = ifelse(present, all_t, paste0(all_t, " \u2014")),
    pt_col  = ifelse(present, PALETTE[all_t], "#CCCCCC"),
    txt_col = ifelse(present, "black", "gray60"),
    alpha   = ifelse(present, 1.0, 0.45),
    stringsAsFactors = FALSE
  )
  ggplot(df) +
    geom_point(aes(x = 0, y = y), color = df$pt_col,
               alpha = df$alpha, size = 3, shape = 16) +
    geom_text(aes(x = 0.2, y = y, label = label),
              hjust = 0, size = 2.6, color = df$txt_col) +
    xlim(-0.2, 5) +
    theme_void() +
    theme(plot.background = element_rect(fill = "white", color = NA),
          plot.margin = margin(8, 4, 8, 4))
}

add_bold_labels <- function(p, obj, id = "cell_type", size = 4.5) {
  LabelClusters(
    plot     = p,
    id       = id,
    repel    = TRUE,
    fontface = "bold",
    size     = size,
    box      = TRUE,
    fill     = scales::alpha("white", 0.65),
    color    = "black",
    label.padding = ggplot2::unit(0.15, "lines")
  )
}

# ── Core plot function ────────────────────────────────────────────────────────
# panel_title: displayed as ggtitle on the labeled UMAP (e.g. "Bo - I")
plot_umap_clean <- function(obj, panel_title, out_path, pt.size = 0.6) {
  Idents(obj) <- "cell_type"
  present <- sort(unique(as.character(obj$cell_type)))
  cols    <- get_colors(present)

  make_dp <- function(show_label) {
    p <- DimPlot(obj, reduction = "umap", group.by = "cell_type",
                 label = FALSE, cols = cols, pt.size = pt.size) +
      theme_classic(base_size = 12) +
      theme(legend.position = "none")
    if (show_label) {
      p <- p +
        ggtitle(panel_title) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))
      p <- add_bold_labels(p, obj)
    } else {
      p <- p + theme(plot.title = element_blank())
    }
    p
  }

  combined <- (make_dp(TRUE) | make_dp(FALSE) | make_full_legend(present)) +
    plot_layout(widths = c(5, 5, 2.5))

  ggsave(out_path, plot = combined, width = 20, height = 7,
         dpi = 300, bg = "white")
  cat("  saved:", out_path, "\n")
}

# ── 1. Infusion-product samples ───────────────────────────────────────────────
cat("\n=== Infusion product samples (I) ===\n")
I_list <- readRDS(file.path(base_dir, "all_I_samples_annotated.rds"))
cat("Keys in RDS:", paste(names(I_list), collapse = ", "), "\n")

# RDS key → (panel title, output filename)
I_map <- list(
  Bo_bone_I = list(title = "Bo - I", file = "Bo_I_UMAP_annotated.png"),
  Ca_bone_I = list(title = "Ca - I", file = "Ca_I_UMAP_annotated.png"),
  Me_bone_I = list(title = "Me - I", file = "Me_I_UMAP_annotated.png")
)

for (key in names(I_map)) {
  if (!key %in% names(I_list)) {
    cat("  [SKIP] key not found:", key, "\n"); next
  }
  m <- I_map[[key]]
  plot_umap_clean(I_list[[key]], m$title, file.path(out_I, m$file))
}

# ── 2. Post-infusion samples (AB) ─────────────────────────────────────────────
cat("\n=== Post-infusion samples (AB) ===\n")
AB_list <- readRDS(file.path(base_dir, "all_AB_samples_annotated.rds"))
cat("Keys in RDS:", paste(names(AB_list), collapse = ", "), "\n")

# RDS key → panel title
AB_titles <- list(
  Bo_blood_AB = "Bo - A",
  Bo_bone_AB  = "Bo - B",
  Ca_blood_AB = "Ca - A",
  Ca_bone_AB  = "Ca - B",
  Me_bone_AB  = "Me - B"
)

for (nm in names(AB_list)) {
  title    <- if (nm %in% names(AB_titles)) AB_titles[[nm]] else nm
  out_path <- file.path(out_AB, paste0(nm, "_UMAP_annotated.png"))
  plot_umap_clean(AB_list[[nm]], title, out_path)
}

cat("\nDone.\n")
