# ==============================================================================
# 14_shared_chains.R
#
# Analisi: famiglie TRAV e TRBV condivise tra i tre pazienti (tutti gli stage)
#
# Unità di analisi: macrofamiglia V-gene
#   TRAV singolo  (es. TRAV12)
#   TRBV singolo  (es. TRBV7)
#   Coppia TRAV + TRBV (es. TRAV12 + TRBV7)
#
# Tre livelli:
#   1. Quali TRBV macrofamilies appaiono in tutti e 3 i pazienti?
#   2. Quali TRAV macrofamilies appaiono in tutti e 3 i pazienti?
#   3. Quali coppie TRAV+TRBV appaiono in tutti e 3 i pazienti?
#
# Per ciascun livello: quante cellule, quanti CDR3 distinti, in quali stage
#
# Dipende da: final_clone_sequences.xlsx (post-decontaminazione)
#
# Output: Fig14a_shared_vgenes_overview.png
#         Fig14b_shared_trbv_usage.png
#         Fig14c_shared_pairs_bubble.png
#         14_shared_chains.xlsx
# ==============================================================================

suppressMessages({
  library(dplyr); library(tidyr); library(ggplot2)
  library(readxl); library(writexl); library(stringr); library(scales); library(patchwork)
  library(forcats)
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

message(sprintf("  Record: %d | Pazienti: %s",
                nrow(fin), paste(sort(unique(fin$patient)), collapse=", ")))

# Helper: classifica sharing
sharing_label <- function(n) {
  case_when(n == 3 ~ "Tutti e 3 i paz.", n == 2 ~ "Esattamente 2 paz.",
            TRUE   ~ "Privata (1 paz.)")
}

SHARE_COLORS <- c("Tutti e 3 i paz." = "#D73027",
                  "Esattamente 2 paz." = "#FC8D59",
                  "Privata (1 paz.)" = "#CCCCCC")

# ── STEP 2: TRBV macrofamiglia ───────────────────────────────────────────────
message("\n--- STEP 2: TRBV macrofamiglia condivisa ---")

trbv_all <- fin %>%
  group_by(macro_TRB) %>%
  summarise(
    n_pazienti       = n_distinct(patient),
    pazienti         = paste(sort(unique(patient)), collapse=" + "),
    stages           = paste(sort(unique(paste(patient, stage, sep="-"))), collapse="; "),
    n_cells_tot      = sum(n_cells),
    n_cells_Bo       = sum(n_cells[patient=="Bo"], na.rm=TRUE),
    n_cells_Ca       = sum(n_cells[patient=="Ca"], na.rm=TRUE),
    n_cells_Me       = sum(n_cells[patient=="Me"], na.rm=TRUE),
    n_cloni_distinti = n_distinct(paste(TRA_cdr3, TRB_cdr3)),
    .groups = "drop"
  ) %>%
  mutate(sharing = sharing_label(n_pazienti)) %>%
  arrange(desc(n_pazienti), desc(n_cells_tot))

message(sprintf("  TRBV in tutti e 3: %d | in 2: %d | private: %d",
                sum(trbv_all$n_pazienti==3), sum(trbv_all$n_pazienti==2),
                sum(trbv_all$n_pazienti==1)))
if (any(trbv_all$n_pazienti==3)) {
  message("  TRBV in tutti e 3:")
  print(trbv_all %>% filter(n_pazienti==3) %>%
          select(macro_TRB, pazienti, n_cells_Bo, n_cells_Ca, n_cells_Me, n_cloni_distinti))
}

# ── STEP 3: TRAV macrofamiglia ───────────────────────────────────────────────
message("\n--- STEP 3: TRAV macrofamiglia condivisa ---")

trav_all <- fin %>%
  group_by(macro_TRA) %>%
  summarise(
    n_pazienti       = n_distinct(patient),
    pazienti         = paste(sort(unique(patient)), collapse=" + "),
    stages           = paste(sort(unique(paste(patient, stage, sep="-"))), collapse="; "),
    n_cells_tot      = sum(n_cells),
    n_cells_Bo       = sum(n_cells[patient=="Bo"], na.rm=TRUE),
    n_cells_Ca       = sum(n_cells[patient=="Ca"], na.rm=TRUE),
    n_cells_Me       = sum(n_cells[patient=="Me"], na.rm=TRUE),
    n_cloni_distinti = n_distinct(paste(TRA_cdr3, TRB_cdr3)),
    .groups = "drop"
  ) %>%
  mutate(sharing = sharing_label(n_pazienti)) %>%
  arrange(desc(n_pazienti), desc(n_cells_tot))

message(sprintf("  TRAV in tutti e 3: %d | in 2: %d | private: %d",
                sum(trav_all$n_pazienti==3), sum(trav_all$n_pazienti==2),
                sum(trav_all$n_pazienti==1)))

# ── STEP 4: Coppie TRAV+TRBV ────────────────────────────────────────────────
message("\n--- STEP 4: Coppie TRAV+TRBV condivise ---")

pair_all <- fin %>%
  group_by(pair_vgene, macro_TRA, macro_TRB) %>%
  summarise(
    n_pazienti       = n_distinct(patient),
    pazienti         = paste(sort(unique(patient)), collapse=" + "),
    stages           = paste(sort(unique(paste(patient, stage, sep="-"))), collapse="; "),
    n_cells_tot      = sum(n_cells),
    n_cells_Bo       = sum(n_cells[patient=="Bo"], na.rm=TRUE),
    n_cells_Ca       = sum(n_cells[patient=="Ca"], na.rm=TRUE),
    n_cells_Me       = sum(n_cells[patient=="Me"], na.rm=TRUE),
    n_cloni_distinti = n_distinct(paste(TRA_cdr3, TRB_cdr3)),
    .groups = "drop"
  ) %>%
  mutate(sharing = sharing_label(n_pazienti)) %>%
  arrange(desc(n_pazienti), desc(n_cells_tot))

message(sprintf("  Coppie in tutti e 3: %d | in 2: %d | private: %d",
                sum(pair_all$n_pazienti==3), sum(pair_all$n_pazienti==2),
                sum(pair_all$n_pazienti==1)))

# ── STEP 5: V-gene usage per paziente (normalizzato) ────────────────────────
message("\n--- STEP 5: V-gene usage per paziente ---")

# Frequenza relativa di ogni TRBV macrofamiglia per paziente
trbv_usage <- fin %>%
  group_by(patient, macro_TRB) %>%
  summarise(n_cells = sum(n_cells), .groups="drop") %>%
  group_by(patient) %>%
  mutate(freq = n_cells / sum(n_cells),
         patient = factor(patient, levels=c("Bo","Ca","Me"))) %>%
  ungroup()

trav_usage <- fin %>%
  group_by(patient, macro_TRA) %>%
  summarise(n_cells = sum(n_cells), .groups="drop") %>%
  group_by(patient) %>%
  mutate(freq = n_cells / sum(n_cells),
         patient = factor(patient, levels=c("Bo","Ca","Me"))) %>%
  ungroup()

# ── STEP 6: Figure ──────────────────────────────────────────────────────────
message("\n--- STEP 6: Figure ---")

# Figura A — panoramica sharing per livello (TRAV / TRBV / coppia)
summ_sharing <- bind_rows(
  trbv_all %>% count(sharing) %>% mutate(catena="TRBV (β)"),
  trav_all %>% count(sharing) %>% mutate(catena="TRAV (α)"),
  pair_all %>% count(sharing) %>% mutate(catena="Coppia α+β")
) %>%
  mutate(
    sharing = factor(sharing, levels=c("Tutti e 3 i paz.","Esattamente 2 paz.","Privata (1 paz.)")),
    catena  = factor(catena,  levels=c("TRBV (β)","TRAV (α)","Coppia α+β"))
  )

p_overview <- ggplot(summ_sharing, aes(x=catena, y=n, fill=sharing)) +
  geom_col(width=0.6, color="white") +
  geom_text(aes(label=n), position=position_stack(vjust=0.5),
            size=4, color="white", fontface="bold") +
  scale_fill_manual(values=SHARE_COLORS, name="Condivisione") +
  scale_y_continuous(expand=expansion(mult=c(0,0.05))) +
  theme_minimal(base_size=12) +
  theme(panel.grid.major.x=element_blank(), legend.position="right") +
  labs(title="V-gene family sharing between patients (all stages)",
       subtitle="Post-contamination filter | CDR3 sequences are unique per patient",
       x=NULL, y="N famiglie V-gene (macrofamiglia)")

ggsave(file.path(FIG, "Fig14a_shared_vgenes_overview.png"),
       p_overview, width=10, height=6, dpi=300, bg="white")
message("Salvata: Fig14a_shared_vgenes_overview.png")

# Figura B — uso TRBV per paziente (top 20) colorato per sharing
top_trbv <- trbv_all %>%
  slice_max(n_cells_tot, n=20) %>%
  pull(macro_TRB)

p_trbv_usage <- trbv_usage %>%
  filter(macro_TRB %in% top_trbv) %>%
  left_join(trbv_all %>% select(macro_TRB, sharing), by="macro_TRB") %>%
  mutate(macro_TRB = fct_reorder(macro_TRB, freq, sum)) %>%
  ggplot(aes(x=freq, y=macro_TRB, fill=patient, alpha=sharing)) +
  geom_col(position="dodge", width=0.7, color="white") +
  scale_fill_manual(values=PAT_COL, labels=PAT_LABEL, name="Patient") +
  scale_alpha_manual(values=c("Tutti e 3 i paz."=1,"Esattamente 2 paz."=0.7,"Privata (1 paz.)"=0.4),
                     name="Sharing") +
  scale_x_continuous(labels=percent_format(accuracy=0.1)) +
  facet_wrap(~patient, labeller=as_labeller(PAT_LABEL), ncol=1, scales="free_x") +
  theme_minimal(base_size=11) +
  theme(strip.text=element_text(face="bold"),
        panel.grid.major.y=element_blank(),
        legend.position="bottom") +
  labs(title="TRBV macrofamily usage per patient (top 20 by total cells)",
       subtitle="Transparency = sharing level | Full opacity = in all 3 patients",
       x="Relative frequency", y="TRBV macrofamily")

ggsave(file.path(FIG, "Fig14b_trbv_usage_by_patient.png"),
       p_trbv_usage, width=10, height=12, dpi=300, bg="white")
message("Salvata: Fig14b_trbv_usage_by_patient.png")

# Figura C — bubble plot coppie condivise in tutti e 3
pairs_3paz <- pair_all %>% filter(n_pazienti==3)
if (nrow(pairs_3paz) > 0) {
  bubble_data <- pairs_3paz %>%
    select(pair_vgene, n_cells_Bo, n_cells_Ca, n_cells_Me) %>%
    pivot_longer(cols=starts_with("n_cells_"),
                 names_to="patient", values_to="n_cells") %>%
    mutate(
      patient    = str_remove(patient, "n_cells_"),
      patient    = factor(patient, levels=c("Bo","Ca","Me")),
      pair_vgene = factor(pair_vgene,
                          levels=rev(pairs_3paz$pair_vgene[order(pairs_3paz$n_cells_Bo)]))
    )

  p_bubble3 <- ggplot(bubble_data,
                      aes(x=patient, y=pair_vgene, size=n_cells, color=patient)) +
    geom_point(alpha=0.85) +
    scale_size_area(max_size=14, name="N cellule") +
    scale_color_manual(values=PAT_COL, labels=PAT_LABEL, guide="none") +
    scale_x_discrete(labels=PAT_LABEL) +
    theme_minimal(base_size=11) +
    theme(panel.grid.major=element_line(color="grey92"),
          axis.text.x=element_text(angle=20, hjust=1),
          axis.text.y=element_text(size=9),
          legend.position="bottom") +
    labs(
      title    = "TRAV+TRBV pairs shared among all 3 patients (any stage)",
      subtitle = paste0(
        "Each dot = cells in that patient using that V-gene pair\n",
        "Note: CDR3 sequences are DIFFERENT between patients — same V-gene, different clones"
      ),
      x=NULL, y="TRAV + TRBV macrofamily pair"
    )

  ggsave(file.path(FIG, "Fig14c_shared_pairs_all3_bubble.png"),
         p_bubble3, width=9, height=max(6, nrow(pairs_3paz)*0.6 + 3),
         dpi=300, bg="white")
  message("Salvata: Fig14c_shared_pairs_all3_bubble.png")
}

# ── STEP 7: Salva tabelle ────────────────────────────────────────────────────
write_xlsx(list(
  "01_TRBV_tutti_3_pazienti"  = trbv_all %>% filter(n_pazienti==3),
  "02_TRBV_2_pazienti"        = trbv_all %>% filter(n_pazienti==2),
  "03_TRAV_tutti_3_pazienti"  = trav_all %>% filter(n_pazienti==3),
  "04_TRAV_2_pazienti"        = trav_all %>% filter(n_pazienti==2),
  "05_Coppie_tutti_3_pazienti"= pairs_3paz,
  "06_Coppie_2_pazienti"      = pair_all %>% filter(n_pazienti==2),
  "07_TRBV_completo"          = trbv_all,
  "08_TRAV_completo"          = trav_all,
  "09_Coppie_complete"        = pair_all
), file.path(TAB, "14_shared_chains.xlsx"))

message("Salvata: 14_shared_chains.xlsx")
message("\nRiepilogo:")
message("  TRBV macrofamiglie in tutti e 3: ", sum(trbv_all$n_pazienti==3))
message("  TRAV macrofamiglie in tutti e 3: ", sum(trav_all$n_pazienti==3))
message("  Coppie TRAV+TRBV in tutti e 3:  ", sum(pair_all$n_pazienti==3))
message("  Coppie in 2 pazienti:            ", sum(pair_all$n_pazienti==2))
