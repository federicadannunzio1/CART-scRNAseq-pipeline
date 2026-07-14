# ==============================================================================
# 14_shared_chains.R
#
# Analisi: famiglie TRAV e TRBV condivise tra pazienti — due categorie biologiche
#
# CATEGORIA 2 — Condivise in A/B (selezione in vivo)
#   Coppie TRAV+TRBV presenti in ≥2 pazienti negli stage di trattamento (A e/o B).
#   Bo → A+B | Ca → A | Me → B (ma solo 7 cellule genuine post-decontaminazione)
#   Interpretazione: potenzialmente selezionate in vivo dalla terapia o
#   dall'ambiente tumorale.
#
# CATEGORIA 3 — Pre-esistenti in I e mantenute/espanse in A/B
#   Coppie TRAV+TRBV presenti in stage I di ≥2 pazienti E osservate in A/B
#   per almeno uno di quegli stessi pazienti.
#   Interpretazione: repertorio pre-infusione condiviso che viene favorito
#   dall'espansione in vivo.
#
# Unità di analisi: macrofamiglia V-gene (macro(TRAV12-1) = TRAV12)
#
# Dipende da: final_clone_sequences.xlsx (post-decontaminazione)
#
# Output: Fig14a_cat2_shared_AB.png
#         Fig14b_cat3_preexisting.png
#         Fig14c_overlap_cat2_cat3.png
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
message("  Distribuzione per paziente × stage:")
print(fin %>% group_by(patient, stage) %>% summarise(n_cloni=n(), n_cells=sum(n_cells), .groups="drop"))

# Dati per stage
infusion_data  <- fin %>% filter(stage == "I")
treatment_data <- fin %>% filter(
  (patient == "Bo" & stage %in% c("A","B")) |
  (patient == "Ca" & stage == "A") |
  (patient == "Me" & stage == "B")
)

message("\n  Cellule per stage:")
message(sprintf("    Infusion (I):    %d cellule, %d pazienti",
                sum(infusion_data$n_cells), n_distinct(infusion_data$patient)))
message(sprintf("    Treatment (A/B): %d cellule, %d pazienti",
                sum(treatment_data$n_cells), n_distinct(treatment_data$patient)))
message("  NOTA: Me-B ha solo 7 cellule genuine post-decontaminazione.")

# Helper: aggrega per pair_vgene in un dataset, restituisce condivisione
summarise_pairs <- function(df) {
  df %>%
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
    arrange(desc(n_pazienti), desc(n_cells_tot))
}

# ── STEP 2: CATEGORIA 2 — Condivise in A/B (selezione in vivo) ──────────────
message("\n--- STEP 2: CATEGORIA 2 — Coppie condivise in stage A/B ---")

cat2_pairs <- summarise_pairs(treatment_data)

cat2_shared <- cat2_pairs %>% filter(n_pazienti >= 2)
cat2_priv   <- cat2_pairs %>% filter(n_pazienti == 1)

message(sprintf("  Coppie in ≥2 pazienti (A/B): %d", nrow(cat2_shared)))
message(sprintf("  Coppie private (1 paz.):      %d", nrow(cat2_priv)))

if (nrow(cat2_shared) > 0) {
  message("  Coppie condivise (A/B):")
  print(cat2_shared %>% select(pair_vgene, pazienti, n_cells_Bo, n_cells_Ca, n_cells_Me,
                                n_cloni_distinti, stages))
}

# ── STEP 3: CATEGORIA 3 — Pre-esistenti in I, mantenute/espanse in A/B ───────
message("\n--- STEP 3: CATEGORIA 3 — Pre-esistenti in I, presenti in A/B ---")

# Per ogni paziente: quali pair_vgene sono in I?
pairs_in_I <- infusion_data %>%
  group_by(patient, pair_vgene) %>%
  summarise(n_cells_I = sum(n_cells), n_cloni_I = n_distinct(paste(TRA_cdr3, TRB_cdr3)),
            .groups="drop")

# Per ogni paziente: quali pair_vgene sono in A/B?
pairs_in_AB <- treatment_data %>%
  group_by(patient, pair_vgene) %>%
  summarise(n_cells_AB = sum(n_cells), n_cloni_AB = n_distinct(paste(TRA_cdr3, TRB_cdr3)),
            .groups="drop")

# Coppie in I per ≥2 pazienti
I_shared <- pairs_in_I %>%
  group_by(pair_vgene) %>%
  filter(n_distinct(patient) >= 2) %>%
  summarise(
    pazienti_in_I    = paste(sort(unique(patient)), collapse=" + "),
    n_paz_I          = n_distinct(patient),
    n_cells_I_Bo     = sum(n_cells_I[patient=="Bo"], na.rm=TRUE),
    n_cells_I_Ca     = sum(n_cells_I[patient=="Ca"], na.rm=TRUE),
    n_cells_I_Me     = sum(n_cells_I[patient=="Me"], na.rm=TRUE),
    .groups = "drop"
  )

# Join con A/B: queste coppie compaiono poi in A/B?
cat3_pairs <- I_shared %>%
  left_join(
    pairs_in_AB %>%
      group_by(pair_vgene) %>%
      summarise(
        pazienti_in_AB   = paste(sort(unique(patient)), collapse=" + "),
        n_paz_AB         = n_distinct(patient),
        n_cells_AB_Bo    = sum(n_cells_AB[patient=="Bo"], na.rm=TRUE),
        n_cells_AB_Ca    = sum(n_cells_AB[patient=="Ca"], na.rm=TRUE),
        n_cells_AB_Me    = sum(n_cells_AB[patient=="Me"], na.rm=TRUE),
        .groups="drop"
      ),
    by = "pair_vgene"
  ) %>%
  mutate(
    presente_in_AB = !is.na(pazienti_in_AB),
    # FC approssimativo per Bo (unico paziente con dati AB affidabili)
    FC_Bo = case_when(
      n_cells_I_Bo == 0 ~ NA_real_,
      TRUE              ~ n_cells_AB_Bo / n_cells_I_Bo
    )
  ) %>%
  arrange(desc(n_paz_I), desc(presente_in_AB), desc(FC_Bo))

cat3_maintained <- cat3_pairs %>% filter(presente_in_AB)
cat3_lost       <- cat3_pairs %>% filter(!presente_in_AB)

message(sprintf("  Coppie in I di ≥2 paz.: %d", nrow(cat3_pairs)))
message(sprintf("    → presenti anche in A/B: %d (mantenute/espanse)", nrow(cat3_maintained)))
message(sprintf("    → assenti in A/B:        %d (perse post-infusione)", nrow(cat3_lost)))

if (nrow(cat3_maintained) > 0) {
  message("  Coppie pre-esistenti e mantenute (Cat. 3):")
  print(cat3_maintained %>%
          select(pair_vgene, pazienti_in_I, pazienti_in_AB,
                 n_cells_I_Bo, n_cells_AB_Bo, FC_Bo,
                 n_cells_I_Ca, n_cells_AB_Ca,
                 n_cells_I_Me, n_cells_AB_Me))
}

# ── STEP 4: Overlap Cat2 ∩ Cat3 ─────────────────────────────────────────────
message("\n--- STEP 4: Overlap categorie 2 e 3 ---")

overlap <- intersect(cat2_shared$pair_vgene, cat3_maintained$pair_vgene)
message(sprintf("  Coppie in entrambe Cat2 e Cat3: %d", length(overlap)))
if (length(overlap) > 0) {
  message("  Overlap:")
  cat(paste(" ", overlap, collapse="\n"), "\n")
}

only_cat2 <- setdiff(cat2_shared$pair_vgene, cat3_maintained$pair_vgene)
only_cat3 <- setdiff(cat3_maintained$pair_vgene, cat2_shared$pair_vgene)
message(sprintf("  Solo Cat2 (selezionate in vivo, non pre-esistenti condivise): %d", length(only_cat2)))
message(sprintf("  Solo Cat3 (pre-esistenti, mantenute ma non shared in A/B):    %d", length(only_cat3)))

# ── STEP 5: Figure ──────────────────────────────────────────────────────────
message("\n--- STEP 5: Figure ---")

CAT_COL <- c("Cat2: shared A/B" = "#D73027",
             "Cat3: pre-existing I→AB" = "#4575B4",
             "Cat2 ∩ Cat3" = "#762A83",
             "Private" = "#CCCCCC")

# Figura A — Categoria 2: bubble plot coppie shared in A/B
if (nrow(cat2_shared) > 0) {
  bub2 <- cat2_shared %>%
    select(pair_vgene, n_cells_Bo, n_cells_Ca, n_cells_Me, pazienti) %>%
    pivot_longer(cols=starts_with("n_cells_"), names_to="patient", values_to="n_cells") %>%
    mutate(
      patient    = str_remove(patient, "n_cells_"),
      patient    = factor(patient, levels=c("Bo","Ca","Me")),
      pair_vgene = fct_reorder(pair_vgene, n_cells, sum)
    )

  p_cat2 <- ggplot(bub2, aes(x=patient, y=pair_vgene, size=n_cells, color=patient)) +
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
      title    = "Categoria 2 — TRAV+TRBV shared in treatment stages (A/B)",
      subtitle = "Coppie V-gene in ≥2 pazienti in stage A/B | CDR3 diverse tra pazienti",
      x=NULL, y="TRAV + TRBV macrofamily"
    )

  ggsave(file.path(FIG, "Fig14a_cat2_shared_AB.png"),
         p_cat2, width=9, height=max(5, nrow(cat2_shared)*0.5 + 3),
         dpi=300, bg="white")
  message("Salvata: Fig14a_cat2_shared_AB.png")
} else {
  message("  [Cat2] Nessuna coppia condivisa in A/B — figura non generata")
}

# Figura B — Categoria 3: coppie pre-esistenti in I, mantenute in A/B
if (nrow(cat3_maintained) > 0) {
  # Mostra: n_cells in I vs AB per Bo (paziente con dati affidabili)
  bub3 <- cat3_maintained %>%
    select(pair_vgene, n_cells_I_Bo, n_cells_AB_Bo, n_cells_I_Ca, n_cells_AB_Ca,
           n_cells_I_Me, n_cells_AB_Me, FC_Bo) %>%
    pivot_longer(
      cols = -c(pair_vgene, FC_Bo),
      names_to = c("stage_type", "patient"),
      names_pattern = "n_cells_(I|AB)_(Bo|Ca|Me)"
    ) %>%
    mutate(
      stage_type = factor(stage_type, levels=c("I","AB")),
      patient    = factor(patient, levels=c("Bo","Ca","Me")),
      pair_vgene = fct_reorder(pair_vgene, value, sum)
    )

  p_cat3 <- ggplot(bub3, aes(x=stage_type, y=pair_vgene, size=value, color=patient)) +
    geom_point(alpha=0.8, position=position_dodge(width=0.5)) +
    scale_size_area(max_size=12, name="N cellule") +
    scale_color_manual(values=PAT_COL, labels=PAT_LABEL, name="Patient") +
    facet_wrap(~patient, labeller=as_labeller(PAT_LABEL), ncol=3) +
    theme_minimal(base_size=11) +
    theme(panel.grid.major=element_line(color="grey92"),
          axis.text.y=element_text(size=9),
          strip.text=element_text(face="bold"),
          legend.position="none") +
    labs(
      title    = "Categoria 3 — TRAV+TRBV pre-esistenti in I, presenti anche in A/B",
      subtitle = "Coppie in I di ≥2 pazienti e osservate in A/B | I=infusion, AB=treatment",
      x="Stage", y="TRAV + TRBV macrofamily"
    )

  ggsave(file.path(FIG, "Fig14b_cat3_preexisting.png"),
         p_cat3, width=12, height=max(5, nrow(cat3_maintained)*0.45 + 3),
         dpi=300, bg="white")
  message("Salvata: Fig14b_cat3_preexisting.png")
} else {
  message("  [Cat3] Nessuna coppia mantenuta I→AB — figura non generata")
}

# Figura C — Venn/upset: overlap Cat2 e Cat3
venn_data <- tibble(
  categoria = c("Solo Cat2\n(selected in vivo)", "Cat2 ∩ Cat3\n(pre-existing+selected)",
                "Solo Cat3\n(pre-existing, maintained)"),
  n = c(length(only_cat2), length(overlap), length(only_cat3)),
  fill = c("#D73027","#762A83","#4575B4")
) %>% filter(n > 0)

if (nrow(venn_data) > 0) {
  p_venn <- ggplot(venn_data, aes(x=categoria, y=n, fill=fill)) +
    geom_col(width=0.5, color="white") +
    geom_text(aes(label=n), vjust=-0.4, size=5, fontface="bold") +
    scale_fill_identity() +
    scale_y_continuous(expand=expansion(mult=c(0,0.15))) +
    theme_minimal(base_size=12) +
    theme(panel.grid.major.x=element_blank()) +
    labs(
      title    = "Overlap tra Categoria 2 e Categoria 3",
      subtitle = "Cat2 = shared in A/B | Cat3 = pre-existing in I, maintained in A/B",
      x=NULL, y="N coppie TRAV+TRBV"
    )

  ggsave(file.path(FIG, "Fig14c_overlap_cat2_cat3.png"),
         p_venn, width=8, height=5, dpi=300, bg="white")
  message("Salvata: Fig14c_overlap_cat2_cat3.png")
}

# ── STEP 6: TRBV frequenze relative in A/B ──────────────────────────────────
message("\n--- STEP 6: TRBV frequenza relativa negli stage A/B ---")

# Frequenza relativa per TRBV per paziente in A/B
trbv_freq_ab <- treatment_data %>%
  group_by(patient, macro_TRB) %>%
  summarise(n_cells = sum(n_cells), .groups="drop") %>%
  group_by(patient) %>%
  mutate(n_tot = sum(n_cells), freq = n_cells / n_tot) %>%
  ungroup()

# Totali per trasparenza
message("  Totale cellule per paziente in A/B:")
print(trbv_freq_ab %>% distinct(patient, n_tot) %>% arrange(patient))

# Pivot wide per confronto cross-paziente
trbv_freq_wide <- trbv_freq_ab %>%
  select(patient, macro_TRB, freq) %>%
  pivot_wider(names_from=patient, values_from=freq, values_fill=0) %>%
  mutate(
    n_paz_present = (Bo > 0) + (ifelse(exists("Ca", where=cur_data()), Ca, 0) > 0) +
                    (ifelse(exists("Me", where=cur_data()), Me, 0) > 0)
  )

# Assicura colonne Ca e Me esistano
if (!"Ca" %in% names(trbv_freq_wide)) trbv_freq_wide$Ca <- 0
if (!"Me" %in% names(trbv_freq_wide)) trbv_freq_wide$Me <- 0

trbv_freq_wide <- trbv_freq_wide %>%
  mutate(n_paz_present = (Bo > 0) + (Ca > 0) + (Me > 0)) %>%
  arrange(desc(n_paz_present), desc(Bo + Ca + Me))

trbv_freq_shared <- trbv_freq_wide %>% filter(n_paz_present >= 2)

message(sprintf("  TRBV in >=2 paz. in A/B (freq. relativa): %d", nrow(trbv_freq_shared)))
message("  Top TRBV shared per freq. relativa:")
print(trbv_freq_shared %>%
        mutate(Bo=round(Bo*100,1), Ca=round(Ca*100,1), Me=round(Me*100,1)) %>%
        select(macro_TRB, n_paz_present, Bo, Ca, Me))

# Nota su incertezza statistica
message("\n  NOTA: Ca-A tot=", unique(trbv_freq_ab$n_tot[trbv_freq_ab$patient=="Ca"]),
        " celle, Me-B tot=", unique(trbv_freq_ab$n_tot[trbv_freq_ab$patient=="Me"]),
        " celle genuine.")
message("  Ogni cellula in Me vale ", round(100/7, 1),
        "% del repertorio — IC ampi per Me e Ca.")

# Figura D — frequenze relative TRBV in A/B per paziente (dot plot)
trbv_freq_plot <- trbv_freq_ab %>%
  filter(macro_TRB %in% trbv_freq_shared$macro_TRB) %>%
  mutate(
    patient    = factor(patient, levels=c("Bo","Ca","Me")),
    macro_TRB  = fct_reorder(macro_TRB, freq, max),
    n_cells_label = paste0("n=", n_cells, "/", n_tot)
  )

p_freq <- ggplot(trbv_freq_plot,
                 aes(x=freq, y=macro_TRB, color=patient, size=n_cells)) +
  geom_point(alpha=0.85) +
  geom_text(aes(label=n_cells_label), hjust=-0.15, size=3, color="grey30") +
  scale_x_continuous(labels=percent_format(accuracy=0.1),
                     expand=expansion(mult=c(0.02, 0.25))) +
  scale_color_manual(values=PAT_COL, labels=PAT_LABEL, name="Patient") +
  scale_size_area(max_size=10, name="N cellule (abs.)") +
  theme_minimal(base_size=12) +
  theme(panel.grid.major.y=element_line(color="grey92"),
        panel.grid.minor=element_blank(),
        legend.position="bottom") +
  labs(
    title    = "TRBV relative frequency in treatment stages (A/B) — shared across patients",
    subtitle = paste0("Solo TRBV presenti in ≥2 pazienti | Annotazione: n cellule / totale\n",
                      "Attenzione: Ca-A=8 celle, Me-B=7 celle genuine — IC statistici molto ampi"),
    x="Frequenza relativa", y="TRBV macrofamiglia"
  )

ggsave(file.path(FIG, "Fig14d_trbv_freq_AB_shared.png"),
       p_freq, width=10, height=max(5, nrow(trbv_freq_shared)*0.55 + 3),
       dpi=300, bg="white")
message("Salvata: Fig14d_trbv_freq_AB_shared.png")

# ── STEP 7: Salva tabelle ────────────────────────────────────────────────────
message("\n--- STEP 7: Salvataggio tabelle ---")

write_xlsx(list(
  "01_Cat2_shared_AB"        = cat2_shared,
  "02_Cat2_private_AB"       = cat2_priv,
  "03_Cat3_maintained_I_AB"  = cat3_maintained,
  "04_Cat3_lost_postI"       = cat3_lost,
  "05_Cat3_tutte_I_shared"   = cat3_pairs,
  "06_Overlap_Cat2_Cat3"     = cat3_maintained %>%
                                 filter(pair_vgene %in% overlap) %>%
                                 left_join(cat2_shared %>%
                                             select(pair_vgene, pazienti_AB=pazienti,
                                                    n_cells_AB_tot=n_cells_tot),
                                           by="pair_vgene"),
  "07_TRBV_freq_AB"          = trbv_freq_wide %>%
                                 mutate(Bo_pct=round(Bo*100,2),
                                        Ca_pct=round(Ca*100,2),
                                        Me_pct=round(Me*100,2)) %>%
                                 select(macro_TRB, n_paz_present, Bo_pct, Ca_pct, Me_pct)
), file.path(TAB, "14_shared_chains.xlsx"))

message("Salvata: 14_shared_chains.xlsx")

message("\n=== RIEPILOGO ===")
message(sprintf("  CATEGORIA 2 — Shared in A/B (≥2 paz.):          %d coppie", nrow(cat2_shared)))
message(sprintf("  CATEGORIA 3 — Pre-esistenti I→AB (≥2 paz.):      %d coppie", nrow(cat3_maintained)))
message(sprintf("  Overlap Cat2 ∩ Cat3:                              %d coppie", length(overlap)))
message(sprintf("  TRBV shared in A/B per freq. relativa (≥2 paz.): %d", nrow(trbv_freq_shared)))
