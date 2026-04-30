library(dplyr)
library(ggplot2)
library(readxl)
library(writexl)
library(stringr)

# ==============================================================================
# 03_private_clones.R
#
# COSA FA:
#   1. Identifica i cloni "condivisi" tra pazienti (identità nucleotidica 100%)
#      → questi sono probabilmente dallo stesso lotto CAR-T, non convergenza
#   2. Produce un report dettagliato di contaminazione
#   3. Filtra i cloni "privati" (non condivisi) per analizzare
#      la vera espansione biologica intra-paziente
#   4. Plotta le famiglie TCR di interesse (target_families)
#
# DIPENDE DA: 01_build_clonotypes.R (full_data in memoria)
#
# CORREZIONI rispetto a 2_fixed_plot_unique_tcr.R:
#   - Rimosso library(xlsx) → sostituito con readxl + writexl
#   - Non legge più da file xlsx: usa full_data direttamente dalla memoria
#   - Colonne aggiornate: TRA_cdr3/TRB_cdr3/TRA_v_gene/TRB_v_gene (minuscolo)
#   - Aggiunto controllo preventivo colonne
# ==============================================================================

# ── Configurazione ─────────────────────────────────────────────────────────────
output_dir <- "/Users/federicadannunzio/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/caruana-project/CART/Code/Seurat_analysis/4_clonotypes_expansion_analysis/res/res_clone_sharing"

# Famiglie beta di interesse biologico — modifica secondo le tue ipotesi
target_families <- c("TRBV7-8", "TRBV2", "TRBV5-1", "TRBV28")

# ── Controllo colonne ──────────────────────────────────────────────────────────
if (!exists("full_data")) {
  stop("full_data non trovato in memoria.\nEsegui prima 01_build_clonotypes.R")
}

cols_attese <- c("TRA_cdr3","TRB_cdr3","TRA_cdr3_nt","TRB_cdr3_nt",
                 "TRA_v_gene","TRB_v_gene","patient","stage","Clone_Quality")
cols_mancanti <- setdiff(cols_attese, colnames(full_data))
if (length(cols_mancanti) > 0) {
  stop("Colonne mancanti in full_data: ", paste(cols_mancanti, collapse=", "),
       "\nRiesegui 01_build_clonotypes.R (versione corretta)")
}

message("✓ Colonne verificate")
message("  Cellule totali in full_data: ", nrow(full_data))

# ── STEP 1: Identificazione cloni condivisi (identità nucleotidica) ────────────
message("\n--- STEP 1: Identificazione cloni condivisi tra pazienti ---")

# Un clone è "condiviso" se la stessa coppia CDR3 aa appare in >1 paziente
# con CDR3 nucleotidica identica → stesso riarrangiamento → stesso lotto CAR-T
shared_clones <- full_data %>%
  filter(Clone_Quality == "Complete") %>%
  group_by(TRA_cdr3, TRB_cdr3, TRA_cdr3_nt, TRB_cdr3_nt) %>%
  filter(n_distinct(patient) > 1) %>%
  summarise(
    pazienti_coinvolti = paste(sort(unique(patient)), collapse=" & "),
    cellule_totali     = n(),
    gene_label         = first(Gene_Label),
    stages_Bo          = paste(sort(unique(stage[patient=="Bo"])), collapse=","),
    stages_Ca          = paste(sort(unique(stage[patient=="Ca"])), collapse=","),
    stages_Me          = paste(sort(unique(stage[patient=="Me"])), collapse=","),
    .groups            = "drop"
  ) %>%
  arrange(desc(cellule_totali))

message("Cloni condivisi trovati: ", nrow(shared_clones))
if (nrow(shared_clones) > 0) {
  message("Riepilogo:")
  print(shared_clones %>%
          select(gene_label, TRA_cdr3, TRB_cdr3,
                 pazienti_coinvolti, cellule_totali))
}
# ── STEP 2: Filtraggio cloni privati ──────────────────────────────────────────
message("\n--- STEP 2: Filtraggio cloni privati ---")

shared_keys <- paste0(shared_clones$TRA_cdr3, "_", shared_clones$TRB_cdr3)

private_data <- full_data %>%
  filter(Clone_Quality == "Complete") %>%
  mutate(clone_key = paste0(TRA_cdr3, "_", TRB_cdr3)) %>%
  filter(!(clone_key %in% shared_keys)) %>%
  select(-clone_key)

message("Cellule nei cloni condivisi (escluse): ",
        nrow(full_data %>% filter(Clone_Quality=="Complete")) - nrow(private_data))
message("Cellule nei cloni privati (mantenute): ", nrow(private_data))
message("% cellule private: ",
        round(nrow(private_data) /
              nrow(full_data %>% filter(Clone_Quality=="Complete")) * 100, 1), "%")

# ── STEP 3: Top cloni privati per paziente ─────────────────────────────────────
message("\n--- STEP 3: Top cloni privati ---")

top_private <- private_data %>%
  group_by(patient, stage, TRA_v_gene, TRB_v_gene,
           Gene_Label, TRA_cdr3, TRB_cdr3) %>%
  summarise(n_cells = n(), .groups="drop") %>%
  arrange(patient, stage, desc(n_cells)) %>%
  group_by(patient) %>%
  slice_head(n=10)

message("Top cloni privati:")
print(top_private %>% select(patient, Gene_Label, TRA_cdr3, TRB_cdr3, n_cells))

# ── STEP 4: Plot famiglie target ───────────────────────────────────────────────
message("\n--- STEP 4: Plot famiglie TCR target ---")

families_available <- unique(private_data$TRB_v_gene)
families_found <- intersect(target_families, families_available)
families_missing <- setdiff(target_families, families_available)

if (length(families_missing) > 0) {
  message("⚠ Famiglie non trovate nei dati privati: ",
          paste(families_missing, collapse=", "))
}
if (length(families_found) == 0) {
  message("⚠ Nessuna delle famiglie target trovata. Plot non generato.")
  plot_data <- data.frame()
} else {
  message("Famiglie trovate: ", paste(families_found, collapse=", "))

  plot_data <- private_data %>%
    filter(TRB_v_gene %in% families_found) %>%
    group_by(patient, TRB_v_gene, TRA_v_gene, TRA_cdr3, TRB_cdr3) %>%
    summarise(espansione = n(), .groups="drop") %>%
    arrange(desc(espansione)) %>%
    group_by(patient) %>%
    slice_head(n=5) %>%
    ungroup() %>%
    mutate(
      clone_label = paste0(TRA_v_gene, " + ", TRB_v_gene, "\n",
                           str_trunc(TRB_cdr3, 18))
    )
}

if (nrow(plot_data) > 0) {
  p <- ggplot(plot_data,
              aes(x = reorder(clone_label, espansione),
                  y = espansione,
                  fill = TRB_v_gene)) +
    geom_bar(stat="identity", color="black", width=0.7) +
    facet_wrap(~patient, scales="free_y", ncol=1) +
    coord_flip() +
    scale_fill_brewer(palette="Set1") +
    theme_minimal() +
    theme(
      strip.text      = element_text(face="bold", size=12),
      legend.position = "bottom",
      axis.text.y     = element_text(size=8)
    ) +
    labs(
      title    = "Famiglie TCR target — Espansione cloni privati",
      subtitle = paste0("Esclusi cloni condivisi tra pazienti (n=",
                        nrow(shared_clones), ")\n",
                        "Famiglie: ", paste(families_found, collapse=", ")),
      x    = "Clone (V alpha + V beta | CDR3 beta)",
      y    = "Numero cellule",
      fill = "Famiglia V beta"
    )

  print(p)
  ggsave(file.path(output_dir, "Grafico_Famiglie_Target_Privato.png"),
         p, width=10, height=10, dpi=300, bg="white")
  message("Grafico salvato.")
} else {
  message("Nessun dato per il grafico con le famiglie selezionate.")
}

# ── STEP 5: Salvataggio ────────────────────────────────────────────────────────
message("\n--- STEP 5: Salvataggio ---")

write_xlsx(list(
  "Cloni_condivisi"      = if(nrow(shared_clones)>0) shared_clones else
                             data.frame(nota="nessuno"),
  "Dati_privati"         = private_data,
  "Top10_cloni_privati"  = top_private,
  "Famiglie_target"      = if(nrow(plot_data)>0) plot_data else
                             data.frame(nota="nessun dato per famiglie selezionate")
), file.path(output_dir, "REPORT_CONTAMINAZIONE_E_PRIVATI.xlsx"))

message("Salvato: REPORT_CONTAMINAZIONE_E_PRIVATI.xlsx in ", output_dir)
