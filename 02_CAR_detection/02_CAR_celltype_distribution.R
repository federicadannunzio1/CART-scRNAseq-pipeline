# ============================================================
#  DISTRIBUZIONE PER TIPO CELLULARE — NUOVE CAR-T IDENTIFICATE
#
#  Questo script legge i file _CAR_classification.csv prodotti
#  dall'analisi precedente e risponde alla domanda:
#
#    "In che tipo cellulare cadono le nuove CAR-T identificate
#     dai Metodi A e B, rispetto al gold standard scREP?"
#
#  Output:
#  ─ barplot stacked per campione (proporzioni per cell_type)
#  ─ heatmap celltype × categoria CAR (frequenze normalizzate)
#  ─ tabella riassuntiva counts + percentuali
#  ─ Excel con tutti i dati
#
#  Input:  cartella con i file *_CAR_classification.csv
#          (prodotti da CAR_detection_from_seurat.R)
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(openxlsx)

# ── PARAMETRI ─────────────────────────────────────────────────

# Cartella dove si trovano i _CAR_classification.csv
csv_dir <- "/Users/federicadannunzio/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/caruana-project/CART/Code/Seurat_analysis/3a_CAR_expr_detection/seurat_methods/"

# Cartella output (può essere la stessa)
out_dir <- csv_dir

# Categorie da analizzare (ordine per il plot)
CAR_LEVELS <- c(
  "scREP_confirmed",
  "new_A_and_B",
  "new_A_only",
  "new_B_only",
  "CAR_negative"
)

# Colori per categoria CAR
CAR_COLORS <- c(
  "scREP_confirmed" = "#264653",
  "new_A_and_B"     = "#E63946",
  "new_A_only"      = "#F4A261",
  "new_B_only"      = "#2A9D8F",
  "CAR_negative"    = "#CCCCCC"
)

# Categorie "nuove" da mettere in evidenza
NEW_CAR_CATS <- c("new_A_and_B", "new_A_only", "new_B_only")

# ─────────────────────────────────────────────────────────────

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 1. CARICAMENTO E UNIONE DI TUTTI I CSV
# ============================================================

cat("Lettura file CSV...\n")

csv_files <- list.files(
  csv_dir,
  pattern = "_CAR_classification\\.csv$",
  full.names = TRUE
)

if (length(csv_files) == 0)
  stop(paste0(
    "Nessun file *_CAR_classification.csv trovato in:\n",
    csv_dir, "\n",
    "Assicurati di aver eseguito prima CAR_detection_from_seurat.R"))

cat(sprintf("  Trovati %d file:\n", length(csv_files)))
for (f in csv_files) cat(sprintf("    %s\n", basename(f)))

all_data <- dplyr::bind_rows(lapply(csv_files, read.csv,
                                     stringsAsFactors = FALSE))

cat(sprintf("\nCellule totali caricate: %d\n", nrow(all_data)))

# ── Verifica colonne ──────────────────────────────────────────
required_cols <- c("cell_type", "CAR_integrated",
                   "sample_name", "IS_CAR_ALLIN_scREP")
missing <- setdiff(required_cols, colnames(all_data))
if (length(missing) > 0)
  stop(paste0("Colonne mancanti: ", paste(missing, collapse=", ")))

# ── Pulizia e fattori ─────────────────────────────────────────
all_data <- all_data %>%
  dplyr::filter(!is.na(CAR_integrated), !is.na(cell_type)) %>%
  dplyr::mutate(
    CAR_integrated = factor(CAR_integrated, levels = CAR_LEVELS),
    # Accorcia i nomi dei cell types molto lunghi per i plot
    cell_type_short = dplyr::recode(cell_type,
      "Cytotoxic CD8+ T cells"   = "CD8+ Cytotox.",
      "Naive CD8+ T cells"       = "CD8+ Naive",
      "Naive CD4+ T cells"       = "CD4+ Naive",
      "Memory T cells"           = "Memory T",
      "Effector CD4+ T cells"    = "CD4+ Effector",
      "Proliferating T cells"    = "Prolif. T",
      "Proliferating CD8+ T cells" = "CD8+ Prolif.",
      "Proliferating CD4+ T cells" = "CD4+ Prolif.",
      "NKT cells"                = "NKT",
      "MAIT cells"               = "MAIT",
      "gamma-delta T cells"      = "γδ T",
      "Tregs"                    = "Tregs",
      "Th1 cells"                = "Th1",
      "Th2 cells"                = "Th2",
      "Th17 cells"               = "Th17",
      "Tfh cells"                = "Tfh",
      .default = cell_type
    )
  )

cat(sprintf("Cellule dopo pulizia:   %d\n", nrow(all_data)))
cat(sprintf("Campioni:               %d\n",
            length(unique(all_data$sample_name))))
cat(sprintf("Tipi cellulari:         %d\n",
            length(unique(all_data$cell_type))))
cat("\nDistribuzione categorie CAR:\n")
print(table(all_data$CAR_integrated))

# ============================================================
# 2. TABELLA RIEPILOGATIVA: celltype × categoria CAR
# ============================================================

cat("\n\n── Tabella celltype × categoria CAR ──\n")

# Counts assoluti
counts_tbl <- all_data %>%
  dplyr::count(cell_type_short, CAR_integrated,
               name = "n_cells") %>%
  tidyr::pivot_wider(
    names_from  = CAR_integrated,
    values_from = n_cells,
    values_fill = 0
  )

# Aggiungi totale riga
counts_tbl <- counts_tbl %>%
  dplyr::mutate(TOTAL = rowSums(dplyr::select(., -cell_type_short))) %>%
  dplyr::arrange(desc(TOTAL))

cat("\nConteggi per tipo cellulare:\n")
print(as.data.frame(counts_tbl), row.names = FALSE)

# Percentuali: per ogni categoria CAR, % per cell_type
pct_tbl <- all_data %>%
  dplyr::count(cell_type_short, CAR_integrated,
               name = "n_cells") %>%
  dplyr::group_by(CAR_integrated) %>%
  dplyr::mutate(
    pct_within_class = round(n_cells / sum(n_cells) * 100, 1)
  ) %>%
  dplyr::ungroup()

# ============================================================
# 3. BARPLOT STACKED — proporzione celltype per categoria CAR
#    Domanda: "Dentro ogni categoria, che cellule ci sono?"
# ============================================================

cat("\n── Plot 1: distribuzione cell_type dentro ogni categoria ──\n")

# Ordina i cell types per frequenza totale
ct_order <- all_data %>%
  dplyr::count(cell_type_short, sort = TRUE) %>%
  dplyr::pull(cell_type_short)

pct_tbl$cell_type_short <- factor(pct_tbl$cell_type_short,
                                   levels = rev(ct_order))

# Palette cell types (distinguibile fino a 16 categorie)
n_ct <- length(ct_order)
ct_palette <- setNames(
  colorRampPalette(c(
    "#264653","#2A9D8F","#52B788","#80B918",
    "#F4A261","#E76F51","#E9C46A","#8338EC",
    "#3A86FF","#FB5607","#FF006E","#8AC926",
    "#1982C4","#6A4C93","#FFBE0B","#B5838D"
  ))(n_ct),
  ct_order
)

# Solo le categorie con cellule (esclude CAR_negative per focus)
plot_data_1 <- pct_tbl %>%
  dplyr::filter(CAR_integrated != "CAR_negative")

p1 <- ggplot(plot_data_1,
             aes(x = CAR_integrated,
                 y = pct_within_class,
                 fill = cell_type_short)) +
  geom_bar(stat = "identity", color = "white",
           linewidth = 0.3, width = 0.75) +
  scale_fill_manual(values = ct_palette, name = "Tipo cellulare") +
  scale_y_continuous(expand = c(0, 0),
                     labels = function(x) paste0(x, "%")) +
  scale_x_discrete(labels = c(
    "scREP_confirmed" = "scREP\nconfirmato\n(gold)",
    "new_A_and_B"     = "Nuovi\nA ∩ B\n(alta conf.)",
    "new_A_only"      = "Nuovi\nsolo A\n(firma)",
    "new_B_only"      = "Nuovi\nsolo B\n(kNN)"
  )) +
  labs(
    title    = "Distribuzione per tipo cellulare nelle categorie CAR",
    subtitle = "% di cellule per tipo all'interno di ogni categoria — tutti i campioni aggregati",
    x        = NULL,
    y        = "% cellule nella categoria"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold", size = 12,
                                    hjust = 0.5),
    plot.subtitle   = element_text(size = 9, hjust = 0.5,
                                    color = "gray40"),
    axis.text.x     = element_text(size = 9, color = "black",
                                    lineheight = 1.2),
    legend.position = "right",
    legend.text     = element_text(size = 8),
    legend.key.size = unit(0.4, "cm"),
    panel.grid.major.y = element_line(color = "gray90",
                                       linewidth = 0.4)
  ) +
  guides(fill = guide_legend(ncol = 1, reverse = TRUE))

ggsave(paste0(out_dir, "celltype_distribution_per_CARclass.png"),
       plot = p1, width = 12, height = 7, dpi = 300, bg = "white")
cat("  → celltype_distribution_per_CARclass.png\n")

# ============================================================
# 4. DOTPLOT / HEATMAP — enrichment di ogni celltype per classe
#    Domanda: "Quale celltype è arricchito nelle nuove CAR-T
#               rispetto al background CAR-?"
# ============================================================

cat("\n── Plot 2: enrichment per cell_type vs CAR_negative ──\n")

# Calcola proporzione di ogni cell_type in ogni classe CAR
# e normalizza rispetto alla proporzione nel CAR_negative
# (log2 fold enrichment)

bg_pct <- all_data %>%
  dplyr::filter(CAR_integrated == "CAR_negative") %>%
  dplyr::count(cell_type_short, name = "n_bg") %>%
  dplyr::mutate(pct_bg = n_bg / sum(n_bg) * 100)

enrich_data <- pct_tbl %>%
  dplyr::filter(CAR_integrated != "CAR_negative") %>%
  dplyr::left_join(bg_pct, by = "cell_type_short") %>%
  dplyr::mutate(
    pct_bg   = dplyr::coalesce(pct_bg, 0.001),  # evita log(0)
    log2_enrich = log2((pct_within_class + 0.1) /
                       (pct_bg + 0.1))
  )

# Ordina i cell types per enrichment medio nelle nuove CAR
ct_enrich_order <- enrich_data %>%
  dplyr::filter(CAR_integrated == "new_A_and_B") %>%
  dplyr::arrange(desc(log2_enrich)) %>%
  dplyr::pull(cell_type_short)

# Aggiungi quelli mancanti in coda
ct_enrich_order <- c(
  ct_enrich_order,
  setdiff(unique(enrich_data$cell_type_short), ct_enrich_order)
)

enrich_data$cell_type_short <- factor(
  enrich_data$cell_type_short,
  levels = ct_enrich_order
)
enrich_data$CAR_integrated <- factor(
  enrich_data$CAR_integrated,
  levels = c("scREP_confirmed","new_A_and_B",
             "new_A_only","new_B_only")
)

p2 <- ggplot(enrich_data,
             aes(x = CAR_integrated,
                 y = cell_type_short,
                 fill = log2_enrich,
                 size = pct_within_class)) +
  geom_point(shape = 21, color = "white", stroke = 0.5) +
  scale_fill_gradient2(
    low      = "#2166AC",
    mid      = "#F7F7F7",
    high     = "#B2182B",
    midpoint = 0,
    name     = "log2 enrichment\nvs CAR-",
    limits   = c(-3, 3),
    oob      = scales::squish
  ) +
  scale_size_continuous(
    name   = "% nella classe",
    range  = c(1, 10),
    breaks = c(1, 5, 15, 30, 50)
  ) +
  scale_x_discrete(labels = c(
    "scREP_confirmed" = "gold\nstandard",
    "new_A_and_B"     = "A ∩ B\nalta conf.",
    "new_A_only"      = "solo A\nfirma",
    "new_B_only"      = "solo B\nkNN"
  )) +
  labs(
    title    = "Enrichment dei tipi cellulari nelle nuove CAR-T",
    subtitle = "Rosso = arricchito rispetto al background CAR-  |  Blu = impoverito",
    x        = NULL,
    y        = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold", size = 12,
                                    hjust = 0.5),
    plot.subtitle   = element_text(size = 9, hjust = 0.5,
                                    color = "gray40"),
    axis.text.x     = element_text(size = 9, lineheight = 1.2),
    axis.text.y     = element_text(size = 9),
    panel.grid      = element_line(color = "gray90",
                                    linewidth = 0.3),
    legend.position = "right"
  )

ggsave(paste0(out_dir, "celltype_enrichment_dotplot.png"),
       plot = p2, width = 10, height = 7, dpi = 300, bg = "white")
cat("  → celltype_enrichment_dotplot.png\n")

# ============================================================
# 5. BARPLOT PER CAMPIONE — focus sulle nuove CAR-T
#    Domanda: "Campione per campione, in quale celltype
#               cadono i nuovi candidati A∩B?"
# ============================================================

cat("\n── Plot 3: per campione — nuove CAR-T per cell_type ──\n")

new_car_data <- all_data %>%
  dplyr::filter(CAR_integrated %in% NEW_CAR_CATS) %>%
  dplyr::count(sample_name, cell_type_short,
               CAR_integrated, name = "n_cells") %>%
  dplyr::group_by(sample_name, CAR_integrated) %>%
  dplyr::mutate(pct = n_cells / sum(n_cells) * 100) %>%
  dplyr::ungroup()

# Ordina campioni
sample_order <- c("Bo_bone_I","Ca_bone_I","Me_bone_I",
                  "Ca_blood_AB","Bo_blood_AB",
                  "Bo_bone_AB","Me_bone_AB")
new_car_data$sample_name <- factor(new_car_data$sample_name,
                                    levels = sample_order)
new_car_data$cell_type_short <- factor(new_car_data$cell_type_short,
                                        levels = ct_order)
new_car_data$CAR_integrated <- factor(new_car_data$CAR_integrated,
  levels = c("new_A_and_B","new_A_only","new_B_only"),
  labels = c("A ∩ B (alta conf.)","Solo A (firma)","Solo B (kNN)")
)

p3 <- ggplot(new_car_data,
             aes(x = sample_name,
                 y = pct,
                 fill = cell_type_short)) +
  geom_bar(stat = "identity", color = "white",
           linewidth = 0.2, width = 0.8) +
  facet_wrap(~ CAR_integrated, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = ct_palette, name = NULL) +
  scale_y_continuous(expand = c(0, 0),
                     labels = function(x) paste0(x, "%")) +
  scale_x_discrete(guide = guide_axis(angle = 35)) +
  labs(
    title    = "Cell type delle nuove CAR-T per campione",
    subtitle = "Ogni pannello = una categoria di confidenza",
    x        = NULL,
    y        = "% cellule nel campione × categoria"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title       = element_text(face = "bold", size = 12,
                                     hjust = 0.5),
    plot.subtitle    = element_text(size = 9, hjust = 0.5,
                                     color = "gray40"),
    strip.background = element_rect(fill = "#F0F4F8",
                                     color = "gray80"),
    strip.text       = element_text(face = "bold", size = 9),
    legend.position  = "right",
    legend.text      = element_text(size = 8),
    legend.key.size  = unit(0.4, "cm"),
    panel.grid.major.y = element_line(color = "gray90",
                                       linewidth = 0.3)
  ) +
  guides(fill = guide_legend(ncol = 1, reverse = TRUE))

ggsave(paste0(out_dir, "celltype_per_sample_new_CAR.png"),
       plot = p3, width = 11, height = 11, dpi = 300, bg = "white")
cat("  → celltype_per_sample_new_CAR.png\n")

# ============================================================
# 6. CONFRONTO DIRETTO gold vs new_A_and_B per cell_type
#    Domanda: "I nuovi A∩B hanno la stessa composizione
#               del gold standard, o sono diversi?"
# ============================================================

cat("\n── Plot 4: gold standard vs new A∩B — confronto composizione ──\n")

compare_data <- all_data %>%
  dplyr::filter(CAR_integrated %in%
                  c("scREP_confirmed", "new_A_and_B")) %>%
  dplyr::count(CAR_integrated, cell_type_short,
               name = "n_cells") %>%
  dplyr::group_by(CAR_integrated) %>%
  dplyr::mutate(pct = n_cells / sum(n_cells) * 100) %>%
  dplyr::ungroup()

compare_data$cell_type_short <- factor(
  compare_data$cell_type_short, levels = ct_order)
compare_data$CAR_integrated <- factor(
  compare_data$CAR_integrated,
  levels = c("scREP_confirmed","new_A_and_B"),
  labels = c("gold standard\n(scREP confirmed)","nuovi candidati\n(A ∩ B)")
)

p4 <- ggplot(compare_data,
             aes(x = cell_type_short,
                 y = pct,
                 fill = CAR_integrated)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.75),
           width = 0.65) +
  scale_fill_manual(
    values = c("gold standard\n(scREP confirmed)" = "#264653",
               "nuovi candidati\n(A ∩ B)"         = "#E63946"),
    name = NULL
  ) +
  scale_y_continuous(expand = c(0, 0),
                     labels = function(x) paste0(x, "%")) +
  labs(
    title    = "Composizione per cell_type: gold standard vs nuovi A∩B",
    subtitle = "Distribuzioni simili = i nuovi candidati assomigliano biologicamente al gold standard",
    x        = NULL,
    y        = "% del gruppo"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title      = element_text(face = "bold", size = 12,
                                    hjust = 0.5),
    plot.subtitle   = element_text(size = 9, hjust = 0.5,
                                    color = "gray40"),
    axis.text.x     = element_text(angle = 40, hjust = 1,
                                    size = 9),
    legend.position = "top",
    panel.grid.major.y = element_line(color = "gray90",
                                       linewidth = 0.4)
  )

ggsave(paste0(out_dir, "celltype_gold_vs_new_AandB.png"),
       plot = p4, width = 12, height = 6, dpi = 300, bg = "white")
cat("  → celltype_gold_vs_new_AandB.png\n")

# ============================================================
# 7. EXPORT EXCEL
# ============================================================

cat("\n── Export Excel ──\n")

wb <- createWorkbook()

# Sheet 1: counts per celltype × categoria CAR
addWorksheet(wb, "Counts_celltype_CARclass")
writeData(wb, "Counts_celltype_CARclass", counts_tbl)

# Sheet 2: percentuali per campione
pct_per_sample <- all_data %>%
  dplyr::count(sample_name, cell_type, CAR_integrated,
               name = "n_cells") %>%
  dplyr::group_by(sample_name, CAR_integrated) %>%
  dplyr::mutate(pct_in_class = round(n_cells / sum(n_cells) * 100, 2)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(sample_name, CAR_integrated, desc(n_cells))

addWorksheet(wb, "Pct_per_sample")
writeData(wb, "Pct_per_sample", pct_per_sample)

# Sheet 3: enrichment
enrich_export <- enrich_data %>%
  dplyr::select(cell_type_short, CAR_integrated,
                n_cells, pct_within_class,
                pct_bg, log2_enrich) %>%
  dplyr::arrange(CAR_integrated, desc(log2_enrich))

addWorksheet(wb, "Enrichment_vs_CARneg")
writeData(wb, "Enrichment_vs_CARneg", enrich_export)

# Sheet 4: dati grezzi nuove CAR-T
new_raw <- all_data %>%
  dplyr::filter(CAR_integrated %in% NEW_CAR_CATS) %>%
  dplyr::select(barcode, sample_name, cell_type,
                IS_CAR_ALLIN_scREP, CAR_integrated,
                dplyr::any_of(c("CAR_method_A",
                                "CAR_sig_score1",
                                "CAR_knn_score",
                                "CAR_method_B"))) %>%
  dplyr::arrange(sample_name, CAR_integrated)

addWorksheet(wb, "New_CAR_raw")
writeData(wb, "New_CAR_raw", new_raw)

xlsx_path <- paste0(out_dir,
                    "CAR_celltype_distribution.xlsx")
saveWorkbook(wb, xlsx_path, overwrite = TRUE)
cat(paste0("  → ", xlsx_path, "\n"))

# ============================================================
# 8. RIEPILOGO TESTUALE
# ============================================================

cat(paste0(
  "\n", strrep("=", 60), "\n",
  "  RIEPILOGO\n",
  strrep("=", 60), "\n\n"
))

cat("Top cell types nelle NUOVE CAR-T (A∩B, alta confidenza):\n")
top_new <- pct_tbl %>%
  dplyr::filter(CAR_integrated == "new_A_and_B") %>%
  dplyr::arrange(desc(pct_within_class)) %>%
  head(5)
for (i in seq_len(nrow(top_new))) {
  cat(sprintf("  %d. %-30s %5.1f%%  (n=%d)\n",
              i,
              top_new$cell_type_short[i],
              top_new$pct_within_class[i],
              top_new$n_cells[i]))
}

cat("\nTop cell types nel GOLD STANDARD (scREP confirmed):\n")
top_gold <- pct_tbl %>%
  dplyr::filter(CAR_integrated == "scREP_confirmed") %>%
  dplyr::arrange(desc(pct_within_class)) %>%
  head(5)
for (i in seq_len(nrow(top_gold))) {
  cat(sprintf("  %d. %-30s %5.1f%%  (n=%d)\n",
              i,
              top_gold$cell_type_short[i],
              top_gold$pct_within_class[i],
              top_gold$n_cells[i]))
}

cat(paste0(
  "\n  Output prodotti:\n",
  "    celltype_distribution_per_CARclass.png\n",
  "      → % celltype dentro ogni categoria CAR\n",
  "    celltype_enrichment_dotplot.png\n",
  "      → quali celltypes sono arricchiti vs background CAR-\n",
  "    celltype_per_sample_new_CAR.png\n",
  "      → nuovi candidati per campione e categoria\n",
  "    celltype_gold_vs_new_AandB.png\n",
  "      → confronto composizione gold vs nuovi A∩B\n",
  "    CAR_celltype_distribution.xlsx\n",
  "      → tutti i dati in tabella\n",
  "\n  COME INTERPRETARE:\n",
  "  Se la distribuzione dei nuovi A∩B assomiglia al gold\n",
  "  standard (dominata da CD8+ citotossici / effettori)\n",
  "  è un buon segnale di specificità biologica.\n",
  "  Se invece new_A_only è dominata da molti tipi diversi\n",
  "  (CD4+, NK-T, γδ...) probabilmente sta catturando\n",
  "  attivazione T generica, non il costrutto CAR.\n",
  strrep("=", 60), "\n"
))

