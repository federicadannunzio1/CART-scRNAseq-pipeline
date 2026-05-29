#!/usr/bin/env Rscript
# Generate one composite PNG per patient:
#   (I_labeled | A_labeled | B_labeled | shared_legend)
# No "Senza label" panel. Titles sized for Word display at full page width.

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

out_dir <- file.path(base_dir, "Patient_composites")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ── Palette / legend (identical to annotation scripts) ────────────────────────
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
              hjust = 0, size = 3.2, color = df$txt_col) +
    xlim(-0.2, 5) +
    theme_void() +
    theme(plot.background = element_rect(fill = "white", color = NA),
          plot.margin = margin(8, 4, 8, 4))
}

add_bold_labels <- function(p, obj, id = "cell_type", size = 4.5) {
  LabelClusters(
    plot = p, id = id, repel = TRUE, fontface = "bold", size = size,
    box = TRUE, fill = scales::alpha("white", 0.65), color = "black",
    label.padding = ggplot2::unit(0.15, "lines")
  )
}

# ── Single labeled panel for one sample ───────────────────────────────────────
make_panel <- function(obj, panel_title, pt.size = 0.6, label_size = 4.5) {
  Idents(obj) <- "cell_type"
  present <- sort(unique(as.character(obj$cell_type)))
  cols    <- get_colors(present)

  p <- DimPlot(obj, reduction = "umap", group.by = "cell_type",
               label = FALSE, cols = cols, pt.size = pt.size) +
    ggtitle(panel_title) +
    theme_classic(base_size = 12) +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 22),
      legend.position = "none"
    )
  add_bold_labels(p, obj, size = label_size)
}

# ── Composite figure: panels + shared legend ───────────────────────────────────
# panels_list: named list, e.g. list("Infusion product" = obj_I, "Blood" = obj_A, ...)
# patient_title: overall title (e.g. "Patient Bo")
# out_path: where to save
make_composite <- function(panels_list, patient_title, out_path,
                           fig_w = 20, fig_h = 7) {
  all_present <- unique(unlist(lapply(panels_list, function(o)
    as.character(o$cell_type))))

  panel_plots <- lapply(names(panels_list), function(nm)
    make_panel(panels_list[[nm]], nm))

  n <- length(panel_plots)
  combined <- Reduce(`|`, panel_plots) | make_full_legend(all_present)

  widths <- c(rep(5, n), 2.5)
  combined <- combined +
    plot_layout(widths = widths) +
    plot_annotation(
      title = patient_title,
      theme = theme(plot.title = element_text(
        face = "bold", hjust = 0.5, size = 32))
    )

  ggsave(out_path, plot = combined, width = fig_w, height = fig_h,
         dpi = 300, bg = "white")
  cat("  saved:", out_path, "\n")
}

# ── Load annotated objects ────────────────────────────────────────────────────
cat("\nLoading annotated objects …\n")
I_list  <- readRDS(file.path(base_dir, "all_I_samples_annotated.rds"))
AB_list <- readRDS(file.path(base_dir, "all_AB_samples_annotated.rds"))
cat("I keys: ",  paste(names(I_list),  collapse = ", "), "\n")
cat("AB keys:", paste(names(AB_list), collapse = ", "), "\n")

# ── Patient Bo ────────────────────────────────────────────────────────────────
cat("\n[Bo]\n")
make_composite(
  panels_list = list(
    "Infusion product"      = I_list[["Bo_bone_I"]],
    "Blood (day +30)"       = AB_list[["Bo_blood_AB"]],
    "Bone marrow (day +30)" = AB_list[["Bo_bone_AB"]]
  ),
  patient_title = "Patient Bo",
  out_path = file.path(out_dir, "Bo_composite.png"),
  fig_w = 20, fig_h = 7
)

# ── Patient Ca ────────────────────────────────────────────────────────────────
cat("\n[Ca]\n")
make_composite(
  panels_list = list(
    "Infusion product"       = I_list[["Ca_bone_I"]],
    "Blood (day +100)"       = AB_list[["Ca_blood_AB"]],
    "Bone marrow (day +100)" = AB_list[["Ca_bone_AB"]]
  ),
  patient_title = "Patient Ca",
  out_path = file.path(out_dir, "Ca_composite.png"),
  fig_w = 20, fig_h = 7
)

# ── Patient Me (no blood sample) ──────────────────────────────────────────────
cat("\n[Me]\n")
make_composite(
  panels_list = list(
    "Infusion product"        = I_list[["Me_bone_I"]],
    "Bone marrow (day +200)"  = AB_list[["Me_bone_AB"]]
  ),
  patient_title = "Patient Me",
  out_path = file.path(out_dir, "Me_composite.png"),
  fig_w = 14, fig_h = 7
)

cat("\nDone.\n")
