# ==============================================================================
# 11_shared_vgene_not_expanded.R
#
# Messaggio: le macrofamiglie V-gene dei cloni espansi in Bo sono parzialmente
# presenti anche in Ca/Me, ma NON si espandono in quegli stessi pazienti.
# La condivisione è a livello di V-gene (repertorio convergente), non di clone.
#
# Dipende da: final_clone_sequences.xlsx + RISULTATI_expansion_dynamics.xlsx
# ==============================================================================

suppressMessages({
  library(dplyr); library(tidyr); library(ggplot2)
  library(readxl); library(stringr); library(scales)
})

TAB <- "/Users/federicadannunzio/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/caruana-project/CART/Code/07_clonotypes/results/tables"
FIG <- "/Users/federicadannunzio/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/caruana-project/CART/Code/07_clonotypes/results/figures"

macro <- function(x) str_remove(x, "-[0-9]+$")

# ── Carica dati ───────────────────────────────────────────────────────────────
esp <- read_xlsx(file.path(TAB, "RISULTATI_expansion_dynamics.xlsx"),
                 sheet = "02_Cloni_espansi_in_B")
fin <- read_xlsx(file.path(TAB, "final_clone_sequences.xlsx"))

# Coppie macrofamiglia V-gene per i cloni espansi in Bo
bo_exp <- esp %>%
  mutate(macro_TRA = macro(TRA_v_gene),
         macro_TRB = macro(TRB_v_gene),
         pair      = paste0(macro_TRA, " + ", macro_TRB)) %>%
  group_by(pair) %>%
  summarise(n_expanded_Bo = n(),
            n_cells_Bo_B  = sum(n_cells_B),
            .groups = "drop")

# Tieni solo le coppie che appaiono in Ca o Me — separato per stage
# Ca: ha solo stage I (nessun stage B) → mostriamo stage I
# Me: ha stage I e stage B → mostriamo entrambi separatamente
came_vgene <- fin %>%
  filter(patient %in% c("Ca", "Me")) %>%
  mutate(macro_TRA = macro(TRA_v_gene),
         macro_TRB = macro(TRB_v_gene),
         pair      = paste0(macro_TRA, " + ", macro_TRB)) %>%
  group_by(patient, stage, pair) %>%
  summarise(n_cells = sum(n_cells), .groups = "drop")

# Coppie degli espansi Bo presenti in Ca/Me
shared_pairs <- bo_exp %>%
  semi_join(came_vgene, by = "pair") %>%
  arrange(desc(n_cells_Bo_B)) %>%
  pull(pair)

# Costruisce dataset per bubble plot: Bo (espansi) + Ca/Me (presenti)
# Tutti i pazienti e stage da final_clone_sequences (post-decontaminazione)
all_vgene_stage <- fin %>%
  mutate(macro_TRA = macro(TRA_v_gene),
         macro_TRB = macro(TRB_v_gene),
         pair      = paste0(macro_TRA, " + ", macro_TRB)) %>%
  filter(pair %in% shared_pairs) %>%
  group_by(patient, stage, pair) %>%
  summarise(n_cells = sum(n_cells), .groups = "drop")

# Flag "expanded": Bo stage B con ≥5 cellule
expanded_flag <- bo_exp %>%
  filter(pair %in% shared_pairs) %>%
  select(pair) %>%
  mutate(is_expanded = TRUE)

plot_data <- all_vgene_stage %>%
  left_join(expanded_flag, by = "pair") %>%
  mutate(
    is_expanded = ifelse(!is.na(is_expanded) & patient == "Bo" & stage == "B",
                         TRUE, FALSE),
    # Colonna asse x: paziente × stage
    pat_stage = paste0(patient, "\n", stage),
    pat_stage = factor(pat_stage, levels = c(
      "Bo\nI", "Bo\nA", "Bo\nB",
      "Ca\nI", "Ca\nA",
      "Me\nI", "Me\nB"
    )),
    pair = factor(pair, levels = rev(shared_pairs))
  ) %>%
  filter(!is.na(pat_stage))

# ── Bubble plot ───────────────────────────────────────────────────────────────
# Separatore verticale tra pazienti
vlines <- c(3.5, 5.5)  # dopo Bo-B e dopo Ca-A

p_bubble <- ggplot(plot_data,
                   aes(x = pat_stage, y = pair,
                       size  = n_cells,
                       color = is_expanded)) +
  geom_vline(xintercept = vlines, color = "grey70", linetype = "dashed", linewidth = 0.5) +
  geom_point(alpha = 0.85) +
  scale_size_area(max_size = 14, name = "N cellule") +
  scale_color_manual(
    values = c("TRUE" = "#D73027", "FALSE" = "#4575B4"),
    labels = c("TRUE" = "Expanded in B", "FALSE" = "Present (not expanded)"),
    name   = NULL
  ) +
  annotate("text", x = 2,   y = length(shared_pairs) + 0.9,
           label = "Bo (expansion)",  fontface = "bold", size = 3.8, color = "#E64B35") +
  annotate("text", x = 4,   y = length(shared_pairs) + 0.9,
           label = "Ca (failure)",    fontface = "bold", size = 3.8, color = "#4DBBD5") +
  annotate("text", x = 5.5, y = length(shared_pairs) + 0.9,
           label = "Me (partial)",    fontface = "bold", size = 3.8, color = "#00A087") +
  scale_x_discrete(labels = c(
    "Bo\nI" = "Bo\nI", "Bo\nA" = "Bo\nA", "Bo\nB" = "Bo\nB",
    "Ca\nI" = "Ca\nI", "Ca\nA" = "Ca\nA",
    "Me\nI" = "Me\nI", "Me\nB" = "Me\nB"
  )) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major   = element_line(color = "grey92"),
    axis.text.x        = element_text(size = 10),
    axis.text.y        = element_text(size = 9),
    legend.position    = "bottom",
    plot.margin        = margin(t = 30, r = 60, b = 10, l = 10)
  ) +
  labs(
    title    = "V-gene macrofamily pairs: shared but not co-expanded",
    subtitle = paste0(
      "V\u03b1+V\u03b2 pairs of Bo's expanded clones (stage B) also present in Ca/Me\n",
      "Same V-gene family, different CDR3 \u2014 not expanded in Ca or Me\n",
      "Post-contamination filter (6 cross-patient clones removed)"
    ),
    x = NULL, y = "V\u03b1 + V\u03b2 macrofamily pair"
  )

ggsave(file.path(FIG, "Fig11_shared_vgene_not_expanded.png"),
       p_bubble, width = 14, height = 11, dpi = 300, bg = "white")
message("Salvata: Fig11_shared_vgene_not_expanded.png")

# ── Tabella riassuntiva ────────────────────────────────────────────────────────
summary_out <- bo_exp %>%
  filter(pair %in% shared_pairs) %>%
  left_join(
    came_vgene %>%
      group_by(pair) %>%
      summarise(patients_CaMe   = paste(sort(unique(patient)), collapse = "+"),
                n_cells_CaMe    = sum(n_cells),
                .groups = "drop"),
    by = "pair"
  ) %>%
  arrange(desc(n_cells_Bo_B))

writexl::write_xlsx(summary_out,
                    file.path(TAB, "11_shared_vgene_not_expanded.xlsx"))
message("Salvata: 11_shared_vgene_not_expanded.xlsx")

message("\nRiepilogo:")
message("  Coppie V-gene espansi Bo totali: ", nrow(bo_exp))
message("  Di cui presenti in Ca/Me:        ", length(shared_pairs))
message("  Esclusive di Bo:                 ", nrow(bo_exp) - length(shared_pairs))
