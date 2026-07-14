# ==============================================================================
# 13_shared_clones_AB.R
#
# Analisi: coppie TRAV+TRBV macrofamiglia condivise tra tutti i pazienti
# negli stage di trattamento (A e/o B)
#
# Unità di analisi: coppia TRAV_macrofamiglia + TRBV_macrofamiglia
#   macro(TRAV12-1) = TRAV12,  macro(TRBV7-9) = TRBV7
#   pair = "TRAV12 + TRBV7"
#
# Stage considerati per ciascun paziente:
#   Bo → stage A e B
#   Ca → stage A (Ca non ha cellule CAR+ in stage B)
#   Me → stage B (Me non ha stage A)
#
# Per ogni coppia trovata in tutti e 3:
#   - Quante cellule? Quanti CDR3 distinti (cloni diversi con stessa V-gene pair)?
#   - La coppia si espande in Bo ma non in Ca/Me? (collegamento a Fig11)
#   - Caratterizzazione funzionale
#
# Dipende da: final_clone_sequences.xlsx (post-decontaminazione)
#
# Output: Fig13_shared_vgene_pairs_AB.png
#         13_shared_clones_AB.xlsx
# ==============================================================================

suppressMessages({
  library(dplyr); library(tidyr); library(ggplot2)
  library(readxl); library(writexl); library(stringr); library(scales); library(patchwork)
})

TAB <- "/Users/federicadannunzio/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/caruana-project/CART/Code/07_clonotypes/results/tables"
FIG <- "/Users/federicadannunzio/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/caruana-project/CART/Code/07_clonotypes/results/figures"

PAT_COL   <- c(Bo = "#E64B35", Ca = "#4DBBD5", Me = "#00A087")
PAT_LABEL <- c(Bo = "Bo (expansion)", Ca = "Ca (failure)", Me = "Me (partial)")

macro <- function(x) str_remove(x, "-[0-9]+$")

# ── STEP 1: Carica dati ─────────────────────────────────────────────────────
message("\n--- STEP 1: Caricamento final_clone_sequences.xlsx ---")

fin <- read_xlsx(file.path(TAB, "final_clone_sequences.xlsx")) %>%
  mutate(macro_TRA  = macro(TRA_v_gene),
         macro_TRB  = macro(TRB_v_gene),
         pair_vgene = paste0(macro_TRA, " + ", macro_TRB))

message(sprintf("  Record totali: %d", nrow(fin)))
message("  Coppie V-gene uniche per paziente:")
print(fin %>% group_by(patient) %>% summarise(n_pairs=n_distinct(pair_vgene), .groups="drop"))

# ── STEP 2: Filtra stage di trattamento ─────────────────────────────────────
message("\n--- STEP 2: Selezione stage di trattamento (A o B) ---")

ab_data <- fin %>%
  filter(
    (patient == "Bo" & stage %in% c("A","B")) |
    (patient == "Ca" & stage == "A") |
    (patient == "Me" & stage == "B")
  )

message("  Record selezionati (stage A/B):")
print(ab_data %>% count(patient, stage))

# ── STEP 3: Coppie V-gene per paziente ──────────────────────────────────────
message("\n--- STEP 3: Aggregazione per coppia V-gene × paziente ---")

# Aggrega: per ogni coppia V-gene, n_cells e n_cloni_distinti (CDR3 diversi)
pair_per_paz <- ab_data %>%
  group_by(patient, stage, pair_vgene, macro_TRA, macro_TRB) %>%
  summarise(
    n_cells         = sum(n_cells),
    n_cloni_distinti = n_distinct(paste(TRA_cdr3, TRB_cdr3)),
    .groups = "drop"
  )

# Pivot: una riga per coppia, colonne per paziente
pair_wide <- pair_per_paz %>%
  group_by(pair_vgene) %>%
  summarise(
    n_pazienti      = n_distinct(patient),
    pazienti        = paste(sort(unique(patient)), collapse=" + "),
    n_cells_Bo      = sum(n_cells[patient=="Bo"], na.rm=TRUE),
    n_cells_Ca      = sum(n_cells[patient=="Ca"], na.rm=TRUE),
    n_cells_Me      = sum(n_cells[patient=="Me"], na.rm=TRUE),
    n_cloni_Bo      = sum(n_cloni_distinti[patient=="Bo"], na.rm=TRUE),
    n_cloni_Ca      = sum(n_cloni_distinti[patient=="Ca"], na.rm=TRUE),
    n_cloni_Me      = sum(n_cloni_distinti[patient=="Me"], na.rm=TRUE),
    stages          = paste(sort(unique(paste(patient, stage, sep="-"))), collapse="; "),
    .groups = "drop"
  ) %>%
  arrange(desc(n_pazienti), desc(n_cells_Bo + n_cells_Ca + n_cells_Me))

shared_all3  <- pair_wide %>% filter(n_pazienti == 3)
shared_2paz  <- pair_wide %>% filter(n_pazienti == 2)
private_1paz <- pair_wide %>% filter(n_pazienti == 1)

message(sprintf("\n  Coppie V-gene in tutti e 3 i pazienti (A/B): %d", nrow(shared_all3)))
message(sprintf("  Coppie V-gene in esattamente 2 pazienti:    %d", nrow(shared_2paz)))
message(sprintf("  Coppie V-gene private (1 paziente):         %d", nrow(private_1paz)))

if (nrow(shared_all3) > 0) {
  message("\n  Coppie condivise tra tutti e 3 (A/B):")
  print(shared_all3 %>% select(pair_vgene, pazienti,
                                n_cells_Bo, n_cells_Ca, n_cells_Me,
                                n_cloni_Bo, n_cloni_Ca, n_cloni_Me))
}

# ── STEP 4: Collegamento a espansione Bo ─────────────────────────────────────
message("\n--- STEP 4: Le coppie condivise sono espanse in Bo? ---")

espansi <- read_xlsx(file.path(TAB, "RISULTATI_expansion_dynamics.xlsx"),
                     sheet = "02_Cloni_espansi_in_B")

bo_exp_pairs <- espansi %>%
  filter(patient == "Bo") %>%
  mutate(pair_vgene = paste0(macro(TRA_v_gene), " + ", macro(TRB_v_gene))) %>%
  group_by(pair_vgene) %>%
  summarise(n_cloni_espansi_Bo = n(),
            n_cells_Bo_B = sum(n_cells_B), .groups="drop")

shared_all3 <- shared_all3 %>%
  left_join(bo_exp_pairs, by="pair_vgene") %>%
  mutate(
    espansa_in_Bo = !is.na(n_cloni_espansi_Bo),
    n_cloni_espansi_Bo = replace_na(n_cloni_espansi_Bo, 0),
    n_cells_Bo_B       = replace_na(n_cells_Bo_B, 0)
  )

message("  Coppie condivise in tutti e 3 — espanse in Bo:")
print(shared_all3 %>% count(espansa_in_Bo))

# ── STEP 5: Figura ──────────────────────────────────────────────────────────
message("\n--- STEP 5: Figure ---")

# Panel A: panoramica sharing
sharing_summ <- tibble(
  condivisione = c("Tutti e 3 i paz.", "Esattamente 2 paz.", "Private (1 paz.)"),
  n_coppie    = c(nrow(shared_all3), nrow(shared_2paz), nrow(private_1paz))
) %>%
  mutate(condivisione = factor(condivisione,
                                levels=c("Tutti e 3 i paz.","Esattamente 2 paz.","Private (1 paz.)")))

p_panoramica <- ggplot(sharing_summ, aes(x="Stage A/B", y=n_coppie, fill=condivisione)) +
  geom_col(width=0.5, color="white") +
  geom_text(aes(label=n_coppie), position=position_stack(vjust=0.5),
            size=4.5, color="white", fontface="bold") +
  scale_fill_manual(
    values=c("Tutti e 3 i paz."="#D73027","Esattamente 2 paz."="#FC8D59","Private (1 paz.)"="#CCCCCC"),
    name="Condivisione"
  ) +
  scale_y_continuous(expand=expansion(mult=c(0,0.08))) +
  theme_minimal(base_size=12) +
  theme(panel.grid.major.x=element_blank(), legend.position="right") +
  labs(title="TRAV+TRBV pair sharing in treatment stages (A/B)",
       x=NULL, y="N coppie V-gene")

# Panel B: bubble plot coppie condivise tra 3 pazienti — n cellule per paziente
if (nrow(shared_all3) > 0) {
  bubble_data <- shared_all3 %>%
    select(pair_vgene, n_cells_Bo, n_cells_Ca, n_cells_Me, espansa_in_Bo) %>%
    pivot_longer(cols=c(n_cells_Bo, n_cells_Ca, n_cells_Me),
                 names_to="patient", values_to="n_cells") %>%
    mutate(
      patient    = str_remove(patient, "n_cells_"),
      patient    = factor(patient, levels=c("Bo","Ca","Me")),
      pair_vgene = factor(pair_vgene,
                          levels=rev(shared_all3$pair_vgene[order(shared_all3$n_cells_Bo)]))
    )

  p_bubble <- ggplot(bubble_data,
                     aes(x=patient, y=pair_vgene, size=n_cells,
                         color=interaction(patient, espansa_in_Bo))) +
    geom_point(alpha=0.85) +
    scale_size_area(max_size=14, name="N cellule") +
    scale_color_manual(
      values=c("Bo.TRUE"="#D73027","Bo.FALSE"="#FC8D59",
               "Ca.TRUE"="#4DBBD5","Ca.FALSE"="#4DBBD5",
               "Me.TRUE"="#00A087","Me.FALSE"="#00A087"),
      labels=c("Bo.TRUE"="Bo espanso","Bo.FALSE"="Bo presente",
               "Ca.TRUE"="Ca","Ca.FALSE"="Ca",
               "Me.TRUE"="Me","Me.FALSE"="Me"),
      name=NULL, guide="none"
    ) +
    scale_x_discrete(labels=PAT_LABEL) +
    theme_minimal(base_size=11) +
    theme(panel.grid.major=element_line(color="grey92"),
          axis.text.x=element_text(angle=20, hjust=1),
          axis.text.y=element_text(size=9),
          legend.position="bottom") +
    labs(
      title    = "Cells per TRAV+TRBV pair shared across all 3 patients",
      subtitle = "Stage A/B only | Red = pairs expanded in Bo stage B",
      x=NULL, y="TRAV + TRBV macrofamily pair"
    )

  fig13 <- p_panoramica + p_bubble +
    plot_layout(widths=c(1,2)) +
    plot_annotation(
      title    = "Figure 13 — TRAV+TRBV pair sharing across patients in treatment stages",
      subtitle = "Bo: A+B | Ca: A | Me: B | Post-contamination filter",
      tag_levels = "A"
    )
} else {
  fig13 <- p_panoramica +
    plot_annotation(title="Figure 13 — TRAV+TRBV pair sharing in treatment stages (A/B)")
}

ggsave(file.path(FIG, "Fig13_shared_vgene_pairs_AB.png"),
       fig13, width=14, height=8, dpi=300, bg="white")
message("Salvata: Fig13_shared_vgene_pairs_AB.png")

# ── STEP 6: Salva tabelle ────────────────────────────────────────────────────
write_xlsx(list(
  "01_Coppie_tutti_3_pazienti"  = if(nrow(shared_all3)>0)  shared_all3  else data.frame(nota="nessuna"),
  "02_Coppie_2_pazienti"        = if(nrow(shared_2paz)>0)  shared_2paz  else data.frame(nota="nessuna"),
  "03_Coppie_private"           = private_1paz,
  "04_Tutte_coppie_AB"          = pair_wide,
  "05_Dati_per_paz_stage"       = pair_per_paz
), file.path(TAB, "13_shared_clones_AB.xlsx"))

message("Salvata: 13_shared_clones_AB.xlsx")
message("\nRiepilogo:")
message("  Coppie TRAV+TRBV in tutti e 3 i pazienti (stage A/B): ", nrow(shared_all3))
message("  Di cui espanse in Bo stage B: ", sum(shared_all3$espansa_in_Bo, na.rm=TRUE))
message("  Di cui NON espanse in Bo (presenti ma non dominanti): ",
        sum(!shared_all3$espansa_in_Bo, na.rm=TRUE))
