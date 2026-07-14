# ==============================================================================
# 17_trbv_expansion.R
#
# Analisi: ruolo di TRBV7 e TRBV6 nell'espansione CAR-T
#
# Motivazione: TRBV7 è la catena beta della coppia più espansa in Bo
#   (TRAV20+TRBV7, FC=25x) e in Me (TRAV38+TRBV7, FC=15x) — alpha diverse.
#   TRBV6 compare anch'essa in entrambi (TRAV12+TRBV6 Bo; TRAV26+TRBV6 Me).
#   Si indaga se questo riflette una selezione a livello di catena beta,
#   indipendente dalla catena alpha.
#
# Due analisi:
#
#   A. Frequenza e FC della singola catena TRBV (non di coppia):
#      - Quanto è abbondante TRBV7 (e TRBV6) in stage I vs B in ogni paziente?
#      - Il FC I→B di TRBV7 è più alto rispetto agli altri TRBV?
#      - È un pattern comune a Bo e Me ma assente in Ca?
#
#   B. Diversità alpha di TRBV7 (e TRBV6) per stage:
#      - Con quante TRAV diverse si accoppia TRBV7 in I vs B?
#      - L'espansione coinvolge poche alpha (oligoclonale) o molte (beta-driven)?
#      - Stessa domanda per Ca (che non si espande): confronto con Bo e Me.
#
# Dipende da: final_clone_sequences.xlsx (post-decontaminazione)
#
# Output: Fig17a_trbv_frequency_stages.png
#         Fig17b_trbv_fc_comparison.png
#         Fig17c_trbv7_alpha_diversity.png
#         Fig17d_trbv7_alpha_partners.png
#         17_trbv_expansion.xlsx
# ==============================================================================

suppressMessages({
  library(dplyr); library(tidyr); library(ggplot2)
  library(readxl); library(writexl); library(stringr)
  library(scales); library(patchwork); library(forcats)
})

TAB <- "/Users/federicadannunzio/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/caruana-project/CART/Code/07_clonotypes/results/tables"
FIG <- "/Users/federicadannunzio/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/caruana-project/CART/Code/07_clonotypes/results/figures"

PAT_COL   <- c(Bo = "#E64B35", Ca = "#4DBBD5", Me = "#00A087")
PAT_LABEL <- c(Bo = "Bo (expansion)", Ca = "Ca (failure)", Me = "Me (partial)")

TRBV_FOCUS <- c("TRBV7", "TRBV6")   # catene di interesse

macro <- function(x) str_remove(x, "-[0-9]+$")

# ── STEP 1: Carica dati ──────────────────────────────────────────────────────
message("\n--- STEP 1: Caricamento final_clone_sequences.xlsx ---")

fin <- read_xlsx(file.path(TAB, "final_clone_sequences.xlsx")) %>%
  mutate(macro_TRA = macro(TRA_v_gene),
         macro_TRB = macro(TRB_v_gene))

message(sprintf("  Record: %d | Pazienti: %s",
                nrow(fin), paste(sort(unique(fin$patient)), collapse=", ")))

# ── STEP 2: Frequenza TRBV macrofamiglia per paziente × stage ────────────────
message("\n--- STEP 2: Frequenza TRBV per paziente × stage ---")

trbv_freq <- fin %>%
  group_by(patient, stage, macro_TRB) %>%
  summarise(n_cells = sum(n_cells), .groups = "drop") %>%
  group_by(patient, stage) %>%
  mutate(n_tot   = sum(n_cells),
         freq    = n_cells / n_tot) %>%
  ungroup() %>%
  mutate(stage   = factor(stage, levels = c("I","A","B")),
         patient = factor(patient, levels = c("Bo","Ca","Me")),
         focus   = macro_TRB %in% TRBV_FOCUS)

# FC I→B per TRBV macrofamiglia (solo pazienti con stage B: Bo e Me)
trbv_wide <- trbv_freq %>%
  select(patient, stage, macro_TRB, freq) %>%
  pivot_wider(names_from = stage, values_from = freq, values_fill = 0) %>%
  mutate(
    FC_I_to_B = case_when(
      I == 0 & B == 0 ~ NA_real_,
      I == 0 & B >  0 ~ NA_real_,   # non rilevato in I = artefatto campionamento
      TRUE            ~ B / I
    ),
    log2FC_I_to_B = ifelse(is.na(FC_I_to_B), NA_real_, log2(FC_I_to_B)),
    focus = macro_TRB %in% TRBV_FOCUS
  ) %>%
  filter(!is.na(B))   # tieni solo pazienti con stage B

message("  TRBV con FC>=2 in B per paziente:")
print(trbv_wide %>% filter(!is.na(FC_I_to_B), FC_I_to_B >= 2) %>%
        select(patient, macro_TRB, I, B, FC_I_to_B, log2FC_I_to_B, focus) %>%
        arrange(patient, desc(FC_I_to_B)))

message("\n  TRBV7 e TRBV6 — frequenze per paziente × stage:")
print(trbv_freq %>% filter(macro_TRB %in% TRBV_FOCUS) %>%
        select(patient, stage, macro_TRB, n_cells, n_tot, freq) %>%
        arrange(macro_TRB, patient, stage))

# ── STEP 3: Diversità catena alpha per TRBV7 e TRBV6 ────────────────────────
message("\n--- STEP 3: Diversità alpha (TRAV partners) di TRBV7 e TRBV6 ---")

alpha_div <- fin %>%
  filter(macro_TRB %in% TRBV_FOCUS) %>%
  group_by(patient, stage, macro_TRB) %>%
  summarise(
    n_cells_tot       = sum(n_cells),
    n_trav_partners   = n_distinct(macro_TRA),
    n_cloni_cdr3      = n_distinct(paste(TRA_cdr3, TRB_cdr3)),
    trav_list         = paste(sort(unique(macro_TRA)), collapse="; "),
    top_trav          = {
      tb <- fin %>%
        filter(macro_TRB == macro_TRB[1], patient == patient[1], stage == stage[1]) %>%
        group_by(macro_TRA) %>% summarise(n=sum(n_cells)) %>%
        slice_max(n, n=1, with_ties=FALSE)
      paste0(tb$macro_TRA, " (", tb$n, " cells)")
    },
    .groups = "drop"
  ) %>%
  mutate(stage   = factor(stage, levels = c("I","A","B")),
         patient = factor(patient, levels = c("Bo","Ca","Me")))

message("  Diversità alpha per TRBV di focus:")
print(alpha_div %>%
        select(patient, stage, macro_TRB, n_cells_tot, n_trav_partners, n_cloni_cdr3, top_trav))

# Dettaglio: TRAV partner × paziente × stage per TRBV7 e TRBV6
alpha_detail <- fin %>%
  filter(macro_TRB %in% TRBV_FOCUS) %>%
  group_by(patient, stage, macro_TRB, macro_TRA) %>%
  summarise(n_cells       = sum(n_cells),
            n_cloni_cdr3  = n_distinct(paste(TRA_cdr3, TRB_cdr3)),
            .groups = "drop") %>%
  group_by(patient, stage, macro_TRB) %>%
  mutate(freq_within = n_cells / sum(n_cells)) %>%
  ungroup() %>%
  mutate(stage   = factor(stage, levels = c("I","A","B")),
         patient = factor(patient, levels = c("Bo","Ca","Me")))

# ── STEP 4: Figure ───────────────────────────────────────────────────────────
message("\n--- STEP 4: Figure ---")

# ── Figura A: frequenza TRBV7 e TRBV6 attraverso gli stage ──────────────────
focus_freq <- trbv_freq %>% filter(focus)

p_freq <- ggplot(focus_freq,
                 aes(x = stage, y = freq * 100,
                     color = patient, group = patient)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 4) +
  geom_text(aes(label = paste0(round(freq*100, 1), "%")),
            vjust = -0.9, size = 3.5, fontface = "bold") +
  facet_wrap(~ macro_TRB, ncol = 2) +
  scale_color_manual(values = PAT_COL, labels = PAT_LABEL, name = NULL) +
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0, 0.2))) +
  theme_minimal(base_size = 12) +
  theme(strip.text       = element_text(face = "bold", size = 13),
        legend.position  = "bottom",
        panel.grid.minor = element_blank()) +
  labs(title    = "TRBV7 and TRBV6 frequency across stages",
       subtitle = "% of all CAR+ cells using that TRBV macrofamily | TRBV7 is the top-expanded beta in both Bo and Me",
       x = "Stage", y = "Frequency (%)")

ggsave(file.path(FIG, "Fig17a_trbv_frequency_stages.png"),
       p_freq, width = 10, height = 6, dpi = 300, bg = "white")
message("Salvata: Fig17a_trbv_frequency_stages.png")

# ── Figura B: FC I→B tutti i TRBV, evidenziando TRBV7 e TRBV6 ──────────────
fc_plot_data <- trbv_wide %>%
  filter(!is.na(FC_I_to_B), I > 0) %>%
  mutate(
    label   = ifelse(focus, macro_TRB, NA_character_),
    alpha_v = ifelse(focus, 1, 0.5),
    size_v  = ifelse(focus, 3.5, 2)
  )

p_fc <- ggplot(fc_plot_data,
               aes(x = log2FC_I_to_B,
                   y = fct_reorder(macro_TRB, log2FC_I_to_B, .fun = median),
                   color = patient, size = focus)) +
  geom_vline(xintercept = 0, color = "grey70") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_point(alpha = 0.8) +
  geom_text(aes(label = label), hjust = -0.3, size = 3.5,
            fontface = "bold", show.legend = FALSE) +
  facet_wrap(~ patient, ncol = 1, labeller = as_labeller(PAT_LABEL),
             scales = "free_y") +
  scale_color_manual(values = PAT_COL, guide = "none") +
  scale_size_manual(values = c("FALSE" = 2, "TRUE" = 4),
                    labels = c("FALSE" = "Other TRBV", "TRUE" = "TRBV7 / TRBV6"),
                    name = NULL) +
  theme_minimal(base_size = 11) +
  theme(strip.text      = element_text(face = "bold"),
        legend.position = "bottom",
        panel.grid.major.y = element_blank()) +
  labs(title    = "TRBV macrofamily expansion: log₂(FC I→B)",
       subtitle = "Large dots = TRBV7 and TRBV6 | dashed = FC=2 threshold",
       x = "log₂(FC I→B)", y = "TRBV macrofamily")

ggsave(file.path(FIG, "Fig17b_trbv_fc_comparison.png"),
       p_fc, width = 9, height = 12, dpi = 300, bg = "white")
message("Salvata: Fig17b_trbv_fc_comparison.png")

# ── Figura C: diversità alpha (n TRAV partners) per TRBV7 e TRBV6 ───────────
p_adiv <- ggplot(alpha_div,
                 aes(x = stage, y = n_trav_partners,
                     color = patient, group = patient)) +
  geom_line(linewidth = 1.3) +
  geom_point(size = 4) +
  geom_text(aes(label = n_trav_partners), vjust = -0.9, size = 3.5, fontface = "bold") +
  facet_wrap(~ macro_TRB, ncol = 2) +
  scale_color_manual(values = PAT_COL, labels = PAT_LABEL, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2)), limits = c(0, NA)) +
  theme_minimal(base_size = 12) +
  theme(strip.text       = element_text(face = "bold", size = 13),
        legend.position  = "bottom",
        panel.grid.minor = element_blank()) +
  labs(
    title    = "Alpha chain diversity of TRBV7 and TRBV6 across stages",
    subtitle = paste0(
      "N unique TRAV macrofamily partners\n",
      "High diversity in B = expansion involves many alpha chains (beta-driven selection)\n",
      "Low diversity in B = oligoclonal expansion (few alpha partners)"
    ),
    x = "Stage", y = "N unique TRAV partners"
  )

ggsave(file.path(FIG, "Fig17c_trbv7_alpha_diversity.png"),
       p_adiv, width = 10, height = 6, dpi = 300, bg = "white")
message("Salvata: Fig17c_trbv7_alpha_diversity.png")

# ── Figura D: heatmap TRAV partners per TRBV7 × paziente × stage ─────────────
# Mostra quali catene alpha si accoppiano con TRBV7 in ciascun stage e paziente
d_heat <- alpha_detail %>%
  filter(macro_TRB == "TRBV7") %>%
  mutate(pat_stage = paste0(patient, "\n", stage),
         pat_stage = factor(pat_stage,
                            levels = c("Bo\nI","Bo\nA","Bo\nB",
                                       "Ca\nI","Ca\nA","Me\nI","Me\nB"))) %>%
  filter(!is.na(pat_stage))

# Ordina TRAV per totale cellule
trav_order <- d_heat %>%
  group_by(macro_TRA) %>%
  summarise(tot = sum(n_cells)) %>%
  arrange(desc(tot)) %>%
  pull(macro_TRA)

d_heat <- d_heat %>%
  mutate(macro_TRA = factor(macro_TRA, levels = rev(trav_order)))

p_heat <- ggplot(d_heat,
                 aes(x = pat_stage, y = macro_TRA, fill = freq_within)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(n_cells >= 2,
                               paste0(n_cells, "\n(", round(freq_within*100), "%)"),
                               "")),
            size = 2.5, color = "grey20", lineheight = 0.9) +
  geom_vline(xintercept = c(3.5, 5.5), color = "grey50",
             linetype = "dashed", linewidth = 0.5) +
  scale_fill_gradientn(
    colors = c("white","#FEE8C8","#FC8D59","#D73027"),
    name   = "Fraction\nwithin TRBV7",
    labels = percent_format()
  ) +
  scale_x_discrete(position = "bottom") +
  theme_minimal(base_size = 10) +
  theme(panel.grid    = element_blank(),
        axis.text.x   = element_text(size = 9),
        axis.text.y   = element_text(size = 9),
        legend.position = "right") +
  labs(
    title    = "TRAV partners of TRBV7 across patients and stages",
    subtitle = "Each cell = N cells (% of TRBV7 cells) | Vertical dashed lines separate patients",
    x = NULL, y = "TRAV macrofamily"
  )

ggsave(file.path(FIG, "Fig17d_trbv7_alpha_partners.png"),
       p_heat,
       width  = 11,
       height = max(6, length(trav_order) * 0.45 + 3),
       dpi = 300, bg = "white")
message("Salvata: Fig17d_trbv7_alpha_partners.png")

# ── STEP 5: Salva tabelle ─────────────────────────────────────────────────────
message("\n--- STEP 5: Salvataggio tabelle ---")

write_xlsx(list(
  "01_TRBV_freq_per_stage"     = trbv_freq %>%
    select(patient, stage, macro_TRB, n_cells, n_tot, freq, focus),
  "02_TRBV_FC_I_to_B"         = trbv_wide,
  "03_Alpha_diversity_TRBV_focus" = alpha_div,
  "04_Alpha_detail_TRBV7"     = alpha_detail %>% filter(macro_TRB == "TRBV7"),
  "05_Alpha_detail_TRBV6"     = alpha_detail %>% filter(macro_TRB == "TRBV6")
), file.path(TAB, "17_trbv_expansion.xlsx"))

message("Salvata: 17_trbv_expansion.xlsx")

# ── STEP 6: Riepilogo ────────────────────────────────────────────────────────
message("\n=== RIEPILOGO ===")
message("Frequenza TRBV7 per paziente:")
print(trbv_freq %>% filter(macro_TRB == "TRBV7") %>%
        select(patient, stage, freq) %>%
        mutate(freq = round(freq*100, 2)) %>%
        arrange(patient, stage))

message("\nFrequenza TRBV6 per paziente:")
print(trbv_freq %>% filter(macro_TRB == "TRBV6") %>%
        select(patient, stage, freq) %>%
        mutate(freq = round(freq*100, 2)) %>%
        arrange(patient, stage))

message("\nFC TRBV7 in Bo e Me:")
print(trbv_wide %>% filter(macro_TRB == "TRBV7") %>%
        select(patient, I, B, FC_I_to_B, log2FC_I_to_B))

message("\nFC TRBV6 in Bo e Me:")
print(trbv_wide %>% filter(macro_TRB == "TRBV6") %>%
        select(patient, I, B, FC_I_to_B, log2FC_I_to_B))

message("\nDiversità alpha TRBV7 in stage B vs I:")
print(alpha_div %>% filter(macro_TRB == "TRBV7") %>%
        select(patient, stage, n_trav_partners, n_cloni_cdr3, top_trav))
