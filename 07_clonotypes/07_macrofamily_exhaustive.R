# ==============================================================================
# 07_macrofamily_exhaustive.R
#
# Analisi esaustiva: ci sono macrofamiglie (V+J, alpha+beta) condivise tra
# i cloni espansi in Bo e i cloni presenti in Ca e Me?
#
# 3 approcci complementari:
#   A) V+J macrofamily pairs (alpha+beta)
#   B) Edit distance CDR3 beta <= 2 (near-identical TCR)
#   C) VDJdb: specificità antigenica condivisa tra Bo-espansi e Ca/Me
#
# Input: results/tables/ (xlsx files da 03, 04, 05)
# ==============================================================================

suppressMessages({
  library(dplyr); library(tidyr); library(readxl); library(writexl)
  library(stringr); library(ggplot2); library(patchwork)
})

BASE    <- "/Users/federicadannunzio/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/caruana-project/CART/Code/07_clonotypes"
TAB     <- file.path(BASE, "results", "tables")
FIG     <- file.path(BASE, "results", "figures")

# ── Carica dati ───────────────────────────────────────────────────────────────
esp <- read_xlsx(file.path(TAB, "RISULTATI_expansion_dynamics.xlsx"),
                 sheet = "02_Cloni_espansi_in_B")   # cloni espansi Bo (e altri)
fin <- read_xlsx(file.path(TAB, "final_clone_sequences.xlsx"))       # tutti cloni con V+J
vdb <- read_xlsx(file.path(TAB, "VDJdb_search_results.xlsx"))        # VDJdb hits

# Funzione per estrarre macrofamiglia (es. TRBV12-4 → TRBV12)
macro <- function(x) str_remove(x, "-[0-9]+$")

# ── Enricchisce esp con J-gene da fin ─────────────────────────────────────────
# esp non contiene j_gene; fin sì (aggregato per paziente+CDR3)
j_info <- fin %>%
  select(patient, TRA_cdr3, TRB_cdr3, TRA_j_gene, TRB_j_gene) %>%
  distinct()

esp_full <- esp %>%
  left_join(j_info, by = c("patient","TRA_cdr3","TRB_cdr3")) %>%
  mutate(
    macro_TRA_V = macro(TRA_v_gene),
    macro_TRB_V = macro(TRB_v_gene),
    macro_TRA_J = macro(TRA_j_gene),
    macro_TRB_J = macro(TRB_j_gene),
    macro_VJ_pair = paste0(macro_TRA_V,"+",macro_TRB_V,
                           " / ",macro_TRA_J,"+",macro_TRB_J)
  )

# Tutti i cloni (non solo espansi) con macrofamiglie
fin_macro <- fin %>%
  mutate(
    macro_TRA_V = macro(TRA_v_gene),
    macro_TRB_V = macro(TRB_v_gene),
    macro_TRA_J = macro(TRA_j_gene),
    macro_TRB_J = macro(TRB_j_gene),
    macro_VJ_pair = paste0(macro_TRA_V,"+",macro_TRB_V,
                           " / ",macro_TRA_J,"+",macro_TRB_J)
  )

# ==============================================================================
# ANALISI A — V+J macrofamily pairs condivise tra Bo-espansi e Ca/Me
# ==============================================================================
message("\n=== ANALISI A: V+J macrofamily pairs ===")

bo_exp <- esp_full %>%
  filter(patient == "Bo", categoria %in% c("Espanso (FC>=2)","Espanso (non rilevato in I)"))

cat_me  <- fin_macro %>% filter(patient %in% c("Ca","Me"))

bo_vj_pairs <- bo_exp %>%
  count(macro_VJ_pair, macro_TRA_V, macro_TRB_V, macro_TRA_J, macro_TRB_J,
        name = "n_cloni_bo") %>%
  filter(!is.na(macro_VJ_pair), !str_detect(macro_VJ_pair,"NA"))

came_vj_pairs <- cat_me %>%
  count(patient, macro_VJ_pair, name = "n_cloni") %>%
  filter(!is.na(macro_VJ_pair), !str_detect(macro_VJ_pair,"NA"))

shared_vj <- bo_vj_pairs %>%
  inner_join(
    came_vj_pairs %>%
      group_by(macro_VJ_pair) %>%
      summarise(patients_CaMe = paste(sort(unique(patient)), collapse="+"),
                n_cloni_CaMe  = sum(n_cloni), .groups="drop"),
    by = "macro_VJ_pair"
  ) %>%
  arrange(desc(n_cloni_bo))

message("V+J pairs condivise tra Bo-espansi e Ca/Me: ", nrow(shared_vj))
print(shared_vj)

# ==============================================================================
# ANALISI B — Edit distance CDR3 beta <= 2 tra Bo-espansi e Ca/Me
# ==============================================================================
message("\n=== ANALISI B: Edit distance CDR3 beta ≤ 2 ===")

bo_cdr3b <- bo_exp$TRB_cdr3
came_df   <- fin %>%
  filter(patient %in% c("Ca","Me")) %>%
  select(patient, stage, TRB_cdr3, TRA_cdr3, TRB_v_gene, n_cells) %>%
  distinct()

# Calcola distanza per ogni coppia Bo vs Ca/Me
results_edit <- data.frame()
for (b in unique(bo_cdr3b)) {
  if (is.na(b) || b == "") next
  for (i in seq_len(nrow(came_df))) {
    q <- came_df$TRB_cdr3[i]
    if (is.na(q) || q == "") next
    d <- adist(b, q)[1,1]
    if (d <= 2) {
      results_edit <- rbind(results_edit, data.frame(
        bo_TRB_cdr3   = b,
        came_TRB_cdr3 = q,
        came_patient  = came_df$patient[i],
        came_stage    = came_df$stage[i],
        came_TRB_Vgene= came_df$TRB_v_gene[i],
        edit_dist     = d,
        stringsAsFactors = FALSE
      ))
    }
  }
}

if (nrow(results_edit) > 0) {
  results_edit <- results_edit %>% distinct() %>% arrange(edit_dist, bo_TRB_cdr3)
} else {
  results_edit <- data.frame(bo_TRB_cdr3=character(), came_TRB_cdr3=character(),
                             came_patient=character(), came_stage=character(),
                             came_TRB_Vgene=character(), edit_dist=integer())
}

message("Coppie CDR3 beta con edit distance ≤2: ", nrow(results_edit))
if (nrow(results_edit) > 0) {
  message("Dettaglio:")
  print(results_edit)
} else {
  message("  → Nessuna CDR3 beta near-identical (edit dist ≤2) tra Bo-espansi e Ca/Me")
}

# ==============================================================================
# ANALISI C — VDJdb: Bo-espansi hanno specificità antigenica presente anche in Ca/Me?
# ==============================================================================
message("\n=== ANALISI C: VDJdb specificità antigenica ===")

# VDJdb hits per ogni paziente
vdb_bo_exp <- vdb %>%
  filter(patient == "Bo") %>%
  semi_join(bo_exp, by = c("TRB_cdr3")) %>%
  select(patient, TRB_cdr3, antigen.species, antigen.gene, antigen.epitope)

vdb_came <- vdb %>%
  filter(patient %in% c("Ca","Me")) %>%
  select(patient, TRB_cdr3, antigen.species, antigen.gene, antigen.epitope)

message("VDJdb hits nei Bo-espansi: ", n_distinct(vdb_bo_exp$TRB_cdr3), " CDR3b distinti")
message("VDJdb hits in Ca/Me: ", n_distinct(vdb_came$TRB_cdr3), " CDR3b distinti")

# Specificità antigeniche di Bo-espansi presenti anche in Ca o Me
ag_bo <- vdb_bo_exp %>%
  distinct(antigen.species, antigen.gene, antigen.epitope)

ag_shared <- vdb_came %>%
  inner_join(ag_bo, by = c("antigen.species","antigen.gene","antigen.epitope")) %>%
  distinct(patient, antigen.species, antigen.gene, antigen.epitope, TRB_cdr3) %>%
  arrange(antigen.species)

message("Specificità antigeniche di Bo-espansi trovate anche in Ca/Me:")
if (nrow(ag_shared) > 0) {
  print(ag_shared)
} else {
  message("  → Nessuna specificità antigenica condivisa")
}

# ==============================================================================
# RIEPILOGO
# ==============================================================================
message("\n========== RIEPILOGO ==========")
message("A) V+J macrofamily pairs condivise Bo-espansi ↔ Ca/Me:  ", nrow(shared_vj))
message("B) CDR3 beta near-identical (edit≤2) Bo-espansi ↔ Ca/Me: ", nrow(results_edit %>% distinct(bo_TRB_cdr3, came_TRB_cdr3)))
message("C) Specificità antigeniche VDJdb condivise:              ", nrow(ag_shared))

# ==============================================================================
# SALVATAGGIO
# ==============================================================================
write_xlsx(list(
  "A_VJ_pairs_shared"    = shared_vj,
  "B_edit_dist_CDR3b"    = results_edit,
  "C_VDJdb_ag_shared"    = ag_shared,
  "C_VDJdb_bo_expanded"  = vdb_bo_exp
), file.path(TAB, "07_macrofamily_exhaustive.xlsx"))

message("\nSalvato: results/tables/07_macrofamily_exhaustive.xlsx")
