# ============================================================
#  PIPELINE 3 – Confronto I vs AB per paziente
#
#  Prerequisiti:
#    PIPELINE_1 → all_I_samples_annotated.rds
#    PIPELINE_2 → all_AB_samples_annotated.rds
#             oppure all_samples_annotated_COMPLETE.rds
#
#  Genera per ogni paziente (Ca, Bo, Me):
#    A) Pannello UMAP: campioni I + AB affiancati
#       con legenda unificata (stessi tipi = stesso colore)
#    B) Barplot: proporzione tipi cellulari per campione
#    C) DotPlot: marker chiave per confronto tra timepoint
#
#  Output in base_dir/Pipeline_I_vs_AB/
# ============================================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(scales)

# ── UNICO PUNTO DA MODIFICARE ────────────────────────────────
base_dir <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/CART/Code/Seurat_analysis/3_data_cleaning/"
# ─────────────────────────────────────────────────────────────

out_dir <- paste0(base_dir, "Pipeline_I_vs_AB/")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

section <- function(title)
  cat(paste0("\n", strrep("=", 65), "\n  ", title,
             "\n", strrep("=", 65), "\n"))

# ============================================================
# CARICAMENTO
# ============================================================
section("Caricamento oggetti annotati")

# Prova prima il file COMPLETE, poi i singoli RDS
complete_rds <- paste0(base_dir, "all_samples_annotated_COMPLETE.rds")
I_rds        <- paste0(base_dir, "all_I_samples_annotated.rds")
AB_rds       <- paste0(base_dir, "all_AB_samples_annotated.rds")

if (file.exists(complete_rds)) {
  cat(">> Carico all_samples_annotated_COMPLETE.rds\n")
  all_samples <- readRDS(complete_rds)
  nms         <- names(all_samples)
  annotated_I  <- all_samples[grep("_I$", nms)]
  annotated_AB <- all_samples[grep("_AB$", nms)]
} else if (file.exists(I_rds) && file.exists(AB_rds)) {
  cat(">> Carico I e AB separatamente\n")
  annotated_I  <- readRDS(I_rds)
  annotated_AB <- readRDS(AB_rds)
  all_samples  <- c(annotated_I, annotated_AB)
} else {
  stop(paste0(
    "File non trovati. Assicurati che PIPELINE_1 e PIPELINE_2\n",
    "siano state eseguite prima di questa pipeline.\n",
    "Attesi:\n  ", complete_rds, "\n  oppure\n  ",
    I_rds, "\n  e\n  ", AB_rds))
}

cat("\nCampioni caricati:\n")
for (nm in names(all_samples))
  cat(sprintf("  %-20s %5d celle | tipi: %s\n",
              nm, ncol(all_samples[[nm]]),
              paste(sort(unique(as.character(
                all_samples[[nm]]$cell_type))), collapse=", ")))

# ============================================================
# PALETTE E LEGENDA CONSISTENTE
# ============================================================

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
  all_t   <- c(CANONICAL_ORDER,
               setdiff(names(PALETTE), CANONICAL_ORDER))
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

make_mini_dp <- function(obj, title_lbl, shared_types) {
  cols <- get_colors(sort(unique(as.character(obj$cell_type))))
  DimPlot(obj, reduction = "umap", group.by = "cell_type",
          label = TRUE, label.size = 3, repel = TRUE,
          cols = cols, pt.size = 0.5) +
    ggtitle(title_lbl) +
    theme_classic(base_size = 11) +
    theme(plot.title      = element_text(hjust = 0.5, face = "bold"),
          legend.position = "none")
}

# ============================================================
# MAPPA PAZIENTE → CAMPIONI
# ============================================================

patient_map <- list(
  Ca = list(
    I  = "Ca_bone_I",
    AB = c("Ca_blood_AB","Ca_bone_AB")
  ),
  Bo = list(
    I  = "Bo_bone_I",
    AB = c("Bo_blood_AB","Bo_bone_AB")
  ),
  Me = list(
    I  = "Me_bone_I",
    AB = c("Me_bone_AB")
  )
)

# ============================================================
# FUNZIONE: PANNELLO COMPLETO PER PAZIENTE
# ============================================================

# Marcatori di riferimento per DotPlot comparativo
DOT_MARKERS <- c(
  # T cells
  "CD3D","CD4","CD8A","CCR7","SELL",  # Naive
  "GZMB","PRF1","NKG7",               # Cytotoxic
  "FOXP3","IL2RA",                    # Treg
  "MKI67",                            # Prolif
  # Th subtypes
  "TBX21","CXCR3",                    # Th1
  "GATA3","CCR4",                     # Th2
  "RORC","CCR6",                      # Th17
  "CXCR5","BCL6",                     # Tfh
  # NK
  "NCAM1","KLRD1",
  # Myeloid
  "LYZ","S100A8","S100A9","VCAN",
  "CD14","FCGR3A","CSF1R"
)

make_patient_panel <- function(patient_id, pm, all_samples,
                               out_dir) {
  cat(paste0("\n── Paziente ", patient_id,
             " ──────────────────────────────────\n"))

  nm_I   <- pm$I
  nm_ABs <- pm$AB[pm$AB %in% names(all_samples)]

  if (!nm_I %in% names(all_samples)) {
    cat(paste0("[SKIP] ", nm_I, " non trovato negli oggetti caricati.\n"))
    return(invisible(NULL))
  }
  if (length(nm_ABs) == 0) {
    cat(paste0("[SKIP] Nessun campione AB per ", patient_id, "\n"))
    return(invisible(NULL))
  }

  all_objs    <- c(list(all_samples[[nm_I]]),
                   all_samples[nm_ABs])
  all_present <- unique(unlist(lapply(all_objs, function(o)
    unique(as.character(o$cell_type)))))

  # ── A) Pannello UMAP ────────────────────────────────────
  panels <- c(
    list(make_mini_dp(all_samples[[nm_I]], nm_I, all_present)),
    lapply(nm_ABs, function(nm)
      make_mini_dp(all_samples[[nm]], nm, all_present))
  )
  n_panels <- length(panels)
  widths   <- c(rep(5, n_panels), 2.5)

  p_umap <- patchwork::wrap_plots(
    c(panels, list(make_full_legend(all_present))),
    nrow = 1, widths = widths) +
    plot_annotation(
      title = paste0(patient_id, " – I vs AB (legenda unificata)"),
      theme = theme(plot.title =
                      element_text(face = "bold", size = 14,
                                   hjust = 0.5)))

  umap_w <- (n_panels * 5 + 2.5) * 1.6
  ggsave(paste0(out_dir, patient_id, "_I_vs_AB_UMAP.png"),
         plot = p_umap, width = umap_w, height = 7,
         dpi = 300, bg = "white")
  cat(paste0("  → ", patient_id, "_I_vs_AB_UMAP.png\n"))

  # ── B) Barplot proporzione tipi per campione ─────────────
  # Costruisce dataframe con % per tipo per campione
  all_nms <- c(nm_I, nm_ABs)
  df_bar  <- bind_rows(lapply(all_nms, function(nm) {
    ct <- as.character(all_samples[[nm]]$cell_type)
    data.frame(sample    = nm,
               cell_type = ct,
               stringsAsFactors = FALSE)
  })) %>%
    group_by(sample, cell_type) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(sample) %>%
    mutate(pct = round(100 * n / sum(n), 2)) %>%
    ungroup() %>%
    mutate(
      sample    = factor(sample, levels = all_nms),
      cell_type = factor(cell_type,
                         levels = rev(intersect(CANONICAL_ORDER,
                                                unique(cell_type))))
    )

  # Colori solo per i tipi presenti
  bar_types  <- levels(df_bar$cell_type)
  bar_colors <- get_colors(as.character(bar_types))

  p_bar <- ggplot(df_bar,
                  aes(x = sample, y = pct, fill = cell_type)) +
    geom_col(position = "stack", width = 0.75) +
    scale_fill_manual(values = bar_colors, name = "Cell type") +
    scale_y_continuous(labels = percent_format(scale = 1),
                       expand = c(0, 0)) +
    theme_classic(base_size = 11) +
    theme(
      axis.text.x  = element_text(angle = 40, hjust = 1, size = 10),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 10),
      legend.text  = element_text(size = 8),
      legend.key.size = unit(0.35, "cm"),
      plot.title   = element_text(face = "bold", hjust = 0.5, size = 12)
    ) +
    labs(title = paste0(patient_id, " – Composizione cellulare"),
         y = "% cellule")

  ggsave(paste0(out_dir, patient_id, "_I_vs_AB_barplot.png"),
         plot = p_bar, width = max(6, n_panels * 2 + 3),
         height = 7, dpi = 300, bg = "white")
  cat(paste0("  → ", patient_id, "_I_vs_AB_barplot.png\n"))

  # ── C) DotPlot comparativo (merge temporaneo) ───────────
  # Aggiunge colonna "sample_id" e unisce i campioni del paziente
  # in un unico oggetto temporaneo per il DotPlot
  objs_tmp <- lapply(all_nms, function(nm) {
    o <- all_samples[[nm]]
    o$sample_id <- nm
    Idents(o) <- "sample_id"
    o
  })

  # merge() Seurat: unisce i layer mantenendo cell_type
  tryCatch({
    merged <- merge(objs_tmp[[1]],
                    y = if (length(objs_tmp) > 1)
                          objs_tmp[-1] else NULL,
                    add.cell.ids = all_nms)

    # Cell_type label: "sample | tipo" per il DotPlot
    merged$plot_group <- paste0(
      merged$sample_id, "\n",
      as.character(merged$cell_type))
    Idents(merged) <- "plot_group"

    dot_genes <- DOT_MARKERS[DOT_MARKERS %in% rownames(merged)]
    if (length(dot_genes) > 0) {
      p_dot <- DotPlot(merged, features = dot_genes) +
        RotatedAxis() +
        ggtitle(paste0(patient_id,
                       " – Marker chiave (I vs AB)")) +
        theme_classic(base_size = 9) +
        theme(
          plot.title  = element_text(face = "bold", hjust = 0.5,
                                     size = 11),
          axis.text.x = element_text(size = 7)
        ) +
        scale_color_gradient2(low = "lightgrey", mid = "#9B59B6",
                              high = "#C0392B", midpoint = 0,
                              name = "Avg\nExpr")

      dot_h <- max(5, length(unique(Idents(merged))) * 0.4 + 2)
      dot_w <- max(10, length(dot_genes) * 0.45 + 3)
      ggsave(paste0(out_dir, patient_id, "_I_vs_AB_dotplot.png"),
             plot = p_dot, width = dot_w, height = dot_h,
             dpi = 300, bg = "white")
      cat(paste0("  → ", patient_id, "_I_vs_AB_dotplot.png\n"))
    }
    rm(merged, objs_tmp)
  }, error = function(e) {
    cat(paste0("  [WARN] DotPlot comparativo non generato: ",
               conditionMessage(e), "\n"))
  })

  # ── Tabella numerica di sintesi ──────────────────────────
  cat(paste0("\n  Composizione cellulare (%) per campione:\n"))
  df_wide <- df_bar %>%
    select(sample, cell_type, pct) %>%
    pivot_wider(names_from = sample, values_from = pct,
                values_fill = 0) %>%
    arrange(match(as.character(cell_type),
                  intersect(CANONICAL_ORDER,
                            as.character(cell_type))))
  print(as.data.frame(df_wide), row.names = FALSE)

  return(invisible(NULL))
}

# ============================================================
# ESECUZIONE PER TUTTI I PAZIENTI
# ============================================================
section("Pannelli per paziente")

for (patient_id in names(patient_map)) {
  make_patient_panel(patient_id, patient_map[[patient_id]],
                     all_samples, out_dir)
}

# ============================================================
# PANNELLO GLOBALE – tutti i campioni in colonna
# ============================================================
section("Pannello globale")

all_present_global <- unique(unlist(lapply(all_samples, function(o)
  unique(as.character(o$cell_type)))))

p_global_list <- lapply(names(all_samples), function(nm) {
  cols <- get_colors(sort(unique(as.character(
    all_samples[[nm]]$cell_type))))
  DimPlot(all_samples[[nm]], reduction = "umap",
          group.by = "cell_type",
          label = TRUE, label.size = 2.8, repel = TRUE,
          cols = cols, pt.size = 0.4) +
    ggtitle(nm) + theme_classic(base_size = 10) +
    theme(plot.title      = element_text(hjust = 0.5, face = "bold"),
          legend.position = "none")
})

p_global <- (patchwork::wrap_plots(p_global_list, ncol = 1) |
             make_full_legend(all_present_global)) +
  plot_layout(widths = c(10, 2.5)) +
  plot_annotation(
    title = "Tutti i campioni – I e AB",
    theme = theme(plot.title = element_text(face = "bold",
                                            hjust = 0.5, size = 15)))

ggsave(paste0(out_dir, "ALL_samples_I_and_AB.png"),
       plot = p_global,
       width = 14,
       height = length(all_samples) * 7,
       dpi = 300, bg = "white")
cat(paste0("[ALL] → ALL_samples_I_and_AB.png\n"))

# ============================================================
# TABELLA RIASSUNTIVA GLOBALE
# ============================================================
section("Tabella riassuntiva globale")

df_global <- bind_rows(lapply(names(all_samples), function(nm) {
  o <- all_samples[[nm]]
  ct <- as.character(o$cell_type)
  data.frame(sample    = nm,
             cell_type = ct,
             n_total   = ncol(o),
             stringsAsFactors = FALSE)
})) %>%
  group_by(sample, cell_type, n_total) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  mutate(pct = round(100 * n_cells / n_total, 2)) %>%
  select(sample, cell_type, n_cells, pct) %>%
  arrange(sample, match(cell_type, CANONICAL_ORDER))

# Salva tabella riassuntiva
library(openxlsx)
write.xlsx(df_global,
           paste0(out_dir, "cell_type_composition_all_samples.xlsx"))
cat("  → cell_type_composition_all_samples.xlsx\n")

cat(paste0(
  "\n", strrep("=", 65), "\n",
  "  PIPELINE 3 COMPLETATA\n\n",
  "  Output in: ", out_dir, "\n\n",
  "  Per paziente:\n",
  "    Ca_I_vs_AB_UMAP.png      – UMAP affiancati\n",
  "    Ca_I_vs_AB_barplot.png   – proporzione tipi cellulari\n",
  "    Ca_I_vs_AB_dotplot.png   – marker chiave\n",
  "    (idem per Bo e Me)\n\n",
  "  Globale:\n",
  "    ALL_samples_I_and_AB.png\n",
  "    cell_type_composition_all_samples.xlsx\n",
  strrep("=", 65), "\n"
))
