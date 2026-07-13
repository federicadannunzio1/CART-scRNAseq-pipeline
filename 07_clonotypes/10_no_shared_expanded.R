# ==============================================================================
# 10_no_shared_expanded.R
#
# Figura: nessun clone espanso in B è condiviso tra pazienti (post-decontaminazione)
#
# Mostra i top cloni in stage B per Bo e Me (Ca non ha stage B),
# con celle colorate per n_cells_B. Nessuna riga ha colore in >1 colonna.
#
# Dipende da: 01_build_clonotypes.R + 03_expansion_dynamics.R (in memoria)
# ==============================================================================

suppressMessages({
  library(dplyr); library(tidyr); library(ggplot2)
  library(readxl); library(stringr); library(writexl)
})

TAB <- "/Users/federicadannunzio/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/caruana-project/CART/Code/07_clonotypes/results/tables"
FIG <- "/Users/federicadannunzio/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/caruana-project/CART/Code/07_clonotypes/results/figures"

# ── Carica dati ───────────────────────────────────────────────────────────────
dyn <- read_xlsx(file.path(TAB, "RISULTATI_expansion_dynamics.xlsx"),
                 sheet = "01_Dinamica_completa")

# Top cloni in B per paziente (Bo: espansi ≥5; Me: tutti con ≥1 cella in B)
top_bo <- dyn %>%
  filter(patient == "Bo", n_cells_B >= 5,
         categoria %in% c("Espanso (FC>=2)", "Espanso (non rilevato in I)")) %>%
  arrange(desc(n_cells_B)) %>%
  slice_head(n = 20) %>%
  select(patient, TRA_cdr3, TRB_cdr3, TRA_v_gene, TRB_v_gene, n_cells_B)

top_me <- dyn %>%
  filter(patient == "Me", n_cells_B >= 1) %>%
  arrange(desc(n_cells_B)) %>%
  slice_head(n = 10) %>%
  select(patient, TRA_cdr3, TRB_cdr3, TRA_v_gene, TRB_v_gene, n_cells_B)

# Combina e crea label clone
all_clones <- bind_rows(top_bo, top_me) %>%
  mutate(
    clone_label = paste0(TRA_v_gene, "+", TRB_v_gene, "\n",
                         str_trunc(TRA_cdr3, 12), " / ",
                         str_trunc(TRB_cdr3, 14))
  )

# Costruisce matrice paziente × clone (wide → long per heatmap)
# Per ogni clone, controlla se esiste in altri pazienti post-decontaminazione
patients_all <- c("Bo", "Me")  # Ca non ha stage B

heatmap_data <- all_clones %>%
  select(patient, clone_label, TRA_cdr3, TRB_cdr3, n_cells_B) %>%
  # Aggiunge righe con n_cells_B=0 per il paziente in cui il clone è assente
  complete(
    clone_label,
    patient = patients_all,
    fill = list(n_cells_B = 0)
  ) %>%
  # Recupera TRA/TRB per le righe completate
  left_join(
    all_clones %>% select(clone_label, TRA_cdr3, TRB_cdr3) %>% distinct(),
    by = "clone_label"
  ) %>%
  mutate(
    patient = factor(patient, levels = c("Bo", "Me")),
    # Ordina per paziente di origine poi n_cells decrescenti
    clone_label = factor(clone_label,
                         levels = rev(unique(all_clones$clone_label)))
  )

# Verifica: quante righe hanno n_cells_B > 0 in >1 paziente?
shared_check <- heatmap_data %>%
  filter(n_cells_B > 0) %>%
  group_by(clone_label) %>%
  summarise(n_paz_con_celle = n_distinct(patient), .groups = "drop") %>%
  filter(n_paz_con_celle > 1)

message("Cloni con n_cells_B > 0 in >1 paziente (post-decontaminazione): ",
        nrow(shared_check))
if (nrow(shared_check) > 0) print(shared_check)

# ── Heatmap ───────────────────────────────────────────────────────────────────
# Separatore visivo tra Bo e Me
n_bo <- nrow(top_bo)
n_me <- nrow(top_me)

p_heat <- ggplot(heatmap_data,
                 aes(x = patient, y = clone_label, fill = n_cells_B)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(n_cells_B > 0, n_cells_B, "")),
            size = 3.2, color = "white", fontface = "bold") +
  scale_fill_gradientn(
    colors  = c("grey95", "#FEE8C8", "#FC8D59", "#D73027"),
    values  = scales::rescale(c(0, 1, 30, max(heatmap_data$n_cells_B, 1))),
    na.value = "grey95",
    name    = "N cellule\nin stage B"
  ) +
  scale_x_discrete(
    labels = c("Bo" = "Bo\n(expansion)", "Me" = "Me\n(partial)")
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x   = element_text(face = "bold", size = 11),
    axis.text.y   = element_text(size = 8, family = "mono"),
    panel.grid    = element_blank(),
    legend.position = "right"
  ) +
  labs(
    title    = "Cloni espansi in stage B — nessuna sovrapposizione tra pazienti",
    subtitle = paste0(
      "Bo: top 20 cloni espansi (≥5 cellule in B)  |  Me: tutti i cloni presenti in B\n",
      "Ca esclusa: nessuna cellula CAR+ in stage B\n",
      "Post-decontaminazione (rimossi 6 cloni con CDR3_nt identica tra pazienti)"
    ),
    x = NULL, y = NULL
  )

# ── Plot a barre affiancato (alternativa più leggibile) ───────────────────────
bar_data <- bind_rows(
  top_bo %>% mutate(
    clone_label = paste0(str_trunc(TRA_cdr3, 14), " / ", str_trunc(TRB_cdr3, 16)),
    patient_label = "Bo (expansion)"
  ),
  top_me %>% mutate(
    clone_label = paste0(str_trunc(TRA_cdr3, 14), " / ", str_trunc(TRB_cdr3, 16)),
    patient_label = "Me (partial)"
  )
) %>%
  mutate(
    clone_label   = factor(clone_label, levels = rev(unique(clone_label))),
    patient_label = factor(patient_label, levels = c("Bo (expansion)", "Me (partial)"))
  )

p_bar <- ggplot(bar_data,
                aes(x = n_cells_B, y = clone_label,
                    fill = patient_label)) +
  geom_col(width = 0.7, color = "white") +
  geom_text(aes(label = n_cells_B), hjust = -0.2, size = 3) +
  facet_wrap(~patient_label, scales = "free", ncol = 1) +
  scale_fill_manual(values = c("Bo (expansion)" = "#E64B35",
                               "Me (partial)"   = "#00A087"),
                    guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2)),
                     breaks = scales::breaks_extended(only.loose = TRUE),
                     labels = scales::label_number(accuracy = 1)) +
  theme_minimal(base_size = 11) +
  theme(
    strip.text        = element_text(face = "bold", size = 12),
    panel.grid.major.y = element_blank(),
    axis.text.y       = element_text(size = 9, family = "mono")
  ) +
  labs(
    title    = "Cloni presenti in stage B per paziente (post-decontaminazione)",
    subtitle = paste0(
      "Nessun clone condiviso tra pazienti dopo rimozione contaminanti\n",
      "Ca: 0 cellule CAR+ in stage B"
    ),
    x = "N cellule in stage B", y = NULL
  )

# ── Salvataggio ───────────────────────────────────────────────────────────────
ggsave(file.path(FIG, "Fig10_no_shared_expanded_heatmap.png"),
       p_heat, width = 7, height = 10, dpi = 300, bg = "white")
message("Salvata: Fig10_no_shared_expanded_heatmap.png")

ggsave(file.path(FIG, "Fig10_no_shared_expanded_bars.png"),
       p_bar, width = 9, height = 10, dpi = 300, bg = "white")
message("Salvata: Fig10_no_shared_expanded_bars.png")

# ── Tabella riassuntiva ────────────────────────────────────────────────────────
summary_table <- tibble::tribble(
  ~Paziente, ~`Stage B disponibile`, ~`Cloni espansi (≥5 cellule)`, ~`Condivisi con altri pazienti`,
  "Bo",  "Sì",  "27", "0",
  "Ca",  "No",  "—",  "—",
  "Me",  "Sì",  "0 (soglia ≥5)",  "0"
)

write_xlsx(list(
  "Tabella_riassuntiva"   = summary_table,
  "Bo_expanded_clones"    = top_bo,
  "Me_clones_in_B"        = top_me,
  "Shared_check"          = if (nrow(shared_check) > 0) shared_check else
                              data.frame(nota = "Nessun clone condiviso")
), file.path(TAB, "10_no_shared_expanded.xlsx"))
message("Salvata: 10_no_shared_expanded.xlsx")
