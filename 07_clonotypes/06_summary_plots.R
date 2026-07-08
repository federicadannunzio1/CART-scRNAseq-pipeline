# ==============================================================================
# 06_summary_plots.R
# Grafici di sintesi per tutti i risultati clonotipi CAR-T
# Dipende da: 01_build_clonotypes.R + 03_expansion_dynamics.R + 04_conserved_families.R
# ==============================================================================

suppressMessages({
  library(dplyr); library(tidyr); library(ggplot2)
  library(patchwork); library(readxl); library(stringr); library(scales)
  library(forcats)
})

OUT <- "/Users/federicadannunzio/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/caruana-project/CART/Code/07_clonotypes/results/figures"
TAB <- "/Users/federicadannunzio/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/caruana-project/CART/Code/07_clonotypes/results/tables"

PAT_COL  <- c(Bo = "#E64B35", Ca = "#4DBBD5", Me = "#00A087")
PAT_LABEL <- c(Bo = "Bo (expansion)", Ca = "Ca (failure)", Me = "Me (partial)")

# ==============================================================================
# FIGURA 1 — Panoramica cellule CAR+ per paziente × stage
# ==============================================================================
cell_counts <- tibble::tribble(
  ~patient, ~stage, ~n_car, ~n_paired,
  "Bo",     "I",     374,    374,
  "Bo",     "A",    1563,   1563,
  "Bo",     "B",     951,    951,
  "Ca",     "A",      10,     10,
  "Ca",     "I",     545,    545,
  "Me",     "I",     106,    106,
  "Me",     "B",     259,    259
)

paired <- tibble::tribble(
  ~patient, ~n_complete,
  "Bo", 2888,
  "Ca",  555,
  "Me",  365
)

p_cells <- ggplot(cell_counts,
                  aes(x = stage, y = n_car, fill = patient)) +
  geom_col(position = "dodge", width = 0.6, color = "white") +
  geom_text(aes(label = n_car), position = position_dodge(0.6),
            vjust = -0.4, size = 3.5, fontface = "bold") +
  facet_wrap(~patient, labeller = as_labeller(PAT_LABEL), scales = "free_x") +
  scale_fill_manual(values = PAT_COL, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold", size = 11),
        panel.grid.major.x = element_blank()) +
  labs(title = "CAR+ cells per patient and stage",
       subtitle = "Cells with complete TRA+TRB pairing available for clonotype analysis",
       x = "Stage", y = "N cells")

p_paired <- ggplot(paired, aes(x = patient, y = n_complete, fill = patient)) +
  geom_col(width = 0.5, color = "white") +
  geom_text(aes(label = n_complete), vjust = -0.4, size = 4, fontface = "bold") +
  scale_fill_manual(values = PAT_COL, labels = PAT_LABEL, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.x = element_blank()) +
  labs(title = "Complete clonotypes (TRA+TRB paired)",
       x = NULL, y = "N complete clonotypes")

fig1 <- p_cells / p_paired +
  plot_annotation(tag_levels = "A",
                  title = "Figure 1 — CAR+ cell overview")

ggsave(file.path(OUT, "Fig1_CAR_cell_overview.png"),
       fig1, width = 10, height = 9, dpi = 300, bg = "white")
message("Fig1 saved")

# ==============================================================================
# FIGURA 2 — Contaminanti inter-paziente (6 cloni condivisi)
# ==============================================================================
contam <- read_xlsx(file.path(TAB, "contamination_report.xlsx"))

contam_plot <- contam %>%
  mutate(
    clone_label = paste0(str_trunc(TRA_cdr3, 14), " / ", str_trunc(TRB_cdr3, 20)),
    clone_label = factor(clone_label, levels = rev(clone_label))
  )

p_contam <- ggplot(contam_plot,
                   aes(x = dominant_n, y = clone_label, fill = dominant_patient)) +
  geom_col(width = 0.6, color = "white") +
  geom_text(aes(label = paste0(pazienti, "  (dominant: ", dominant_patient, ")")),
            x = 5, hjust = 0, size = 3.2, color = "grey20") +
  scale_fill_manual(values = PAT_COL, name = "Dominant\npatient") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.y = element_blank(),
        legend.position = "right") +
  labs(title = "Clones shared between patients (contamination filter)",
       subtitle = "CDR3_nt identical between patients — removed from non-dominant patient",
       x = "N cells (dominant patient)", y = NULL)

ggsave(file.path(OUT, "Fig2_contamination_report.png"),
       p_contam, width = 11, height = 5, dpi = 300, bg = "white")
message("Fig2 saved")

# ==============================================================================
# FIGURA 3 — Espansione: categorie per paziente
# ==============================================================================
exp_data <- read_xlsx(file.path(TAB, "RISULTATI_expansion_dynamics.xlsx"),
                      sheet = "01_Dinamica_completa")

cat_counts <- exp_data %>%
  filter(categoria != "Assente") %>%
  count(patient, categoria) %>%
  mutate(
    patient  = factor(patient, levels = c("Bo", "Ca", "Me")),
    categoria = factor(categoria,
                       levels = c("Espanso (non rilevato in I)", "Espanso (FC>=2)",
                                  "Stabile", "Contratto"))
  )

cat_colors <- c("Espanso (non rilevato in I)" = "#E41A1C",
                "Espanso (FC>=2)"             = "#FF7F00",
                "Stabile"                     = "#4DAF4A",
                "Contratto"                   = "#377EB8")

p_cat <- ggplot(cat_counts, aes(x = patient, y = n, fill = categoria)) +
  geom_col(position = "fill", width = 0.6, color = "white") +
  geom_text(aes(label = ifelse(n > 5, n, "")),
            position = position_fill(vjust = 0.5),
            size = 3.5, color = "white", fontface = "bold") +
  scale_fill_manual(values = cat_colors, name = "Category") +
  scale_x_discrete(labels = PAT_LABEL) +
  scale_y_continuous(labels = percent) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 20, hjust = 1)) +
  labs(title = "Clonal expansion categories (I -> B)",
       subtitle = "After contamination filter (251 cells removed)",
       x = NULL, y = "Proportion of clonotypes")

p_abs <- ggplot(cat_counts, aes(x = patient, y = n, fill = categoria)) +
  geom_col(position = "stack", width = 0.6, color = "white") +
  geom_text(aes(label = ifelse(n > 5, n, "")),
            position = position_stack(vjust = 0.5),
            size = 3.5, color = "white", fontface = "bold") +
  scale_fill_manual(values = cat_colors, name = "Category") +
  scale_x_discrete(labels = PAT_LABEL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 20, hjust = 1)) +
  labs(title = "Clonal expansion (absolute)",
       x = NULL, y = "N clonotypes")

fig3 <- (p_cat | p_abs) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A",
                  title = "Figure 3 — Clonal expansion dynamics I -> B")

ggsave(file.path(OUT, "Fig3_expansion_categories.png"),
       fig3, width = 12, height = 6, dpi = 300, bg = "white")
message("Fig3 saved")

# ==============================================================================
# FIGURA 4 — Top cloni espansi in B per paziente
# ==============================================================================
espansi <- read_xlsx(file.path(TAB, "RISULTATI_expansion_dynamics.xlsx"),
                     sheet = "02_Cloni_espansi_in_B")

top_exp <- espansi %>%
  group_by(patient) %>%
  slice_max(n_cells_B, n = 10, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    clone_label = paste0(str_trunc(Gene_Label, 22), "\n",
                         str_trunc(TRB_cdr3, 18)),
    patient = factor(patient, levels = c("Bo", "Ca", "Me")),
    FC_label = ifelse(categoria == "Espanso (non rilevato in I)", "undetected in I",
                      paste0("FC=", round(as.numeric(FC_I_to_B), 1)))
  )

p_top <- ggplot(top_exp,
                aes(x = n_cells_B,
                    y = reorder(clone_label, n_cells_B),
                    fill = categoria)) +
  geom_col(width = 0.7, color = "white") +
  geom_text(aes(label = FC_label), hjust = -0.1, size = 3) +
  facet_wrap(~patient, scales = "free", labeller = as_labeller(PAT_LABEL), ncol = 1) +
  scale_fill_manual(values = cat_colors, name = "Category") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.25))) +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold"),
        panel.grid.major.y = element_blank()) +
  labs(title = "Top 10 expanded clones in stage B (per patient)",
       subtitle = "Fold-change relative to stage I shown; 'undetected in I' = below detection limit",
       x = "N cells in stage B", y = NULL)

ggsave(file.path(OUT, "Fig4_top_expanded_clones.png"),
       p_top, width = 11, height = 14, dpi = 300, bg = "white")
message("Fig4 saved")

# ==============================================================================
# FIGURA 5 — VDJdb: specificita' antigenica per paziente
# ==============================================================================
vdjdb <- read_xlsx(file.path(TAB, "VDJdb_search_results.xlsx"))

ag_counts <- vdjdb %>%
  filter(!is.na(antigen.species)) %>%
  distinct(patient, TRB_cdr3, antigen.species) %>%
  count(patient, antigen.species) %>%
  group_by(patient) %>%
  mutate(pct = n / sum(n)) %>%
  ungroup() %>%
  mutate(
    patient = factor(patient, levels = c("Bo", "Ca", "Me")),
    antigen.species = fct_reorder(antigen.species, n, sum)
  )

ag_colors <- c(
  "CMV"        = "#E64B35",
  "SARS-CoV-2" = "#F39B7F",
  "InfluenzaA" = "#4DBBD5",
  "EBV"        = "#00A087",
  "YFV"        = "#3C5488",
  "HCV"        = "#B09C85"
)

p_ag <- ggplot(ag_counts,
               aes(x = patient, y = n, fill = antigen.species)) +
  geom_col(position = "stack", width = 0.55, color = "white") +
  geom_text(aes(label = ifelse(n >= 1, n, "")),
            position = position_stack(vjust = 0.5),
            size = 3.5, color = "white", fontface = "bold") +
  scale_fill_manual(values = ag_colors, name = "Antigen") +
  scale_x_discrete(labels = PAT_LABEL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 20, hjust = 1)) +
  labs(title = "VDJdb antigen specificity matches",
       subtitle = "N distinct CDR3b sequences with VDJdb match per patient",
       x = NULL, y = "N matched CDR3b sequences")

# Tabella: dettaglio Bo (clone dominante CMV)
bo_top <- vdjdb %>%
  filter(patient == "Bo") %>%
  distinct(TRB_cdr3, antigen.species, antigen.gene, antigen.epitope, mhc.a) %>%
  arrange(antigen.species, antigen.epitope)

message("Bo VDJdb matches:\n")
print(bo_top)

ggsave(file.path(OUT, "Fig5_VDJdb_antigen_specificity.png"),
       p_ag, width = 9, height = 6, dpi = 300, bg = "white")
message("Fig5 saved")

# ==============================================================================
# FIGURA 6 — Panel riassuntivo: outcomes × espansione × VDJdb
# ==============================================================================
outcomes <- tibble::tribble(
  ~patient, ~car_in_I_pct, ~car_in_B_pct, ~n_expanded, ~outcome,
  "Bo",      1.9,           22.0,           27,          "Expansion",
  "Ca",     19.1,            0.0,            0,          "Failure",
  "Me",      9.9,            4.9,            1,          "Partial"
)

out_col <- c(Expansion = "#E64B35", Failure = "#4DBBD5", Partial = "#00A087")

p_o1 <- ggplot(outcomes, aes(x = patient, fill = outcome)) +
  geom_col(aes(y = car_in_I_pct), width = 0.4, color = "white") +
  geom_col(aes(y = car_in_B_pct), width = 0.4, alpha = 0.5,
           position = position_nudge(x = 0.42), color = "white") +
  scale_fill_manual(values = out_col, name = "Outcome") +
  annotate("text", x = 0.6, y = 23, label = "I", size = 3.5, color = "grey40") +
  annotate("text", x = 1.0, y = 23, label = "B", size = 3.5, color = "grey40") +
  annotate("text", x = 1.6, y = 23, label = "I", size = 3.5, color = "grey40") +
  annotate("text", x = 2.0, y = 23, label = "B", size = 3.5, color = "grey40") +
  annotate("text", x = 2.6, y = 23, label = "I", size = 3.5, color = "grey40") +
  annotate("text", x = 3.0, y = 23, label = "B", size = 3.5, color = "grey40") +
  scale_x_discrete(labels = PAT_LABEL) +
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0, 0.12))) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 20, hjust = 1)) +
  labs(title = "CAR+ frequency: I vs B", x = NULL, y = "% CAR+ cells")

p_o2 <- ggplot(outcomes, aes(x = patient, y = n_expanded, fill = outcome)) +
  geom_col(width = 0.5, color = "white") +
  geom_text(aes(label = n_expanded), vjust = -0.4, size = 4, fontface = "bold") +
  scale_fill_manual(values = out_col, guide = "none") +
  scale_x_discrete(labels = PAT_LABEL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 20, hjust = 1)) +
  labs(title = "Expanded clones in B (>=5 cells, FC>=2)", x = NULL, y = "N clones")

fig6 <- (p_o1 | p_o2 | p_ag) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "Figure 6 — Summary: CAR-T outcomes, clonal expansion, and antigen specificity",
    subtitle = "3 patients (Bo=expansion, Ca=failure, Me=partial); timepoints I, A, B",
    tag_levels = "A"
  )

ggsave(file.path(OUT, "Fig6_summary_panel.png"),
       fig6, width = 16, height = 6, dpi = 300, bg = "white")
message("Fig6 saved")

message("\n=== TUTTI I GRAFICI SALVATI IN: ", OUT, " ===")
