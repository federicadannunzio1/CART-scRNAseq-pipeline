library(dplyr)
library(tidyr)
library(ggplot2)
library(writexl)
library(stringr)

# ==============================================================================
# 04_expansion_dynamics.R
#
# COSA FA:
#   Risponde alla domanda biologica principale:
#   "Quali popolazioni T si espandono in risposta al CAR-T?"
#
#   Per ogni clone calcola n cellule per stage, frequenza relativa,
#   fold-change I→B, e classifica ogni clone in:
#     De novo in B   → assente in I, compare dopo trattamento
#     Espanso FC≥2   → raddoppia o più
#     Stabile        → presente ma non cambia
#     Contratto      → si riduce
#
#   Cerca inoltre cloni espansi indipendentemente in più pazienti.
#
# DIPENDE DA: 01_build_clonotypes.R (full_data in memoria)
# ==============================================================================

# ── Configurazione ─────────────────────────────────────────────────────────────
output_dir <- "/Users/federicadannunzio/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/caruana-project/CART/Code/Seurat_analysis/4_clonotypes_expansion_analysis/res/res_expansion_dynamics"

soglia_cellule <- 5  # n minimo cellule in stage B per considerare un clone

# ── Controllo colonne ──────────────────────────────────────────────────────────
if (!exists("full_data")) {
  stop("full_data non trovato in memoria.\nEsegui prima 01_build_clonotypes.R")
}

cols_attese <- c("TRA_v_gene","TRB_v_gene","TRA_cdr3","TRB_cdr3",
                 "stage","patient","Clone_Quality","Clone_ID_CDR3","Gene_Label")
cols_mancanti <- setdiff(cols_attese, colnames(full_data))
if (length(cols_mancanti) > 0) {
  stop("Colonne mancanti: ", paste(cols_mancanti, collapse=", "),
       "\nRiesegui 01_build_clonotypes.R (versione corretta)")
}

message("✓ Colonne verificate")
message("  Pazienti: ", paste(sort(unique(full_data$patient)), collapse=", "))
message("  Stage:    ", paste(sort(unique(full_data$stage)),   collapse=", "))

# ── STEP 1: Conteggi per clone × paziente × stage ─────────────────────────────
message("\n--- STEP 1: Conteggi per stage ---")

clone_dynamics <- full_data %>%
  filter(Clone_Quality == "Complete") %>%
  group_by(patient, Clone_ID_CDR3, Gene_Label,
           TRA_cdr3, TRB_cdr3, TRA_v_gene, TRB_v_gene, stage) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  # Completa con 0 per stage in cui il clone non è presente
  complete(
    nesting(patient, Clone_ID_CDR3, Gene_Label,
            TRA_cdr3, TRB_cdr3, TRA_v_gene, TRB_v_gene),
    stage = c("I", "A", "B"),
    fill  = list(n_cells = 0)
  ) %>%
  filter(!(patient == "Me" & stage == "A")) %>%
  # Frequenza relativa per stage: normalizza la profondità di sequenziamento
  group_by(patient, stage) %>%
  mutate(
    n_tot_stage = sum(n_cells),
    freq        = if_else(n_tot_stage > 0, n_cells / n_tot_stage, 0)
  ) %>%
  ungroup()

# ── STEP 2: Pivot wide + fold-change ──────────────────────────────────────────
message("\n--- STEP 2: Calcolo fold-change I → B ---")

clone_wide <- clone_dynamics %>%
  select(patient, Clone_ID_CDR3, Gene_Label,
         TRA_cdr3, TRB_cdr3, TRA_v_gene, TRB_v_gene,
         stage, n_cells, freq) %>%
  pivot_wider(names_from  = stage,
              values_from = c(n_cells, freq),
              values_fill = 0)

# Assicura che tutte le colonne esistano anche se uno stage manca nel dataset
for (s in c("I","A","B")) {
  if (!paste0("n_cells_",s) %in% colnames(clone_wide))
    clone_wide[[paste0("n_cells_",s)]] <- 0L
  if (!paste0("freq_",s) %in% colnames(clone_wide))
    clone_wide[[paste0("freq_",s)]] <- 0
}

clone_wide <- clone_wide %>%
  mutate(
    FC_I_to_B = case_when(
      freq_I == 0 & freq_B  > 0 ~ Inf,       # De novo
      freq_I == 0 & freq_B == 0 ~ NA_real_,   # Assente in entrambi
      TRUE                      ~ freq_B / freq_I
    ),
    FC_I_to_A = case_when(
      freq_I == 0 & freq_A  > 0 ~ Inf,
      freq_I == 0 & freq_A == 0 ~ NA_real_,
      TRUE                      ~ freq_A / freq_I
    ),
    categoria = case_when(
      n_cells_I == 0 & n_cells_B  > 0 ~ "De novo in B",
      n_cells_I == 0 & n_cells_B == 0 ~ "Assente",
      FC_I_to_B >= 2                  ~ "Espanso (FC>=2)",
      FC_I_to_B >= 1                  ~ "Stabile",
      FC_I_to_B <  1                  ~ "Contratto",
      TRUE                            ~ "Altro"
    )
  ) %>%
  arrange(patient, desc(n_cells_B))

message("Categorie per paziente:")
print(table(clone_wide$categoria, clone_wide$patient))

# ── STEP 3: Cloni espansi in B ────────────────────────────────────────────────
message("\n--- STEP 3: Cloni espansi in B (soglia: ", soglia_cellule, " cellule) ---")

cloni_espansi <- clone_wide %>%
  filter(n_cells_B >= soglia_cellule,
         categoria %in% c("Espanso (FC>=2)", "De novo in B")) %>%
  arrange(patient, desc(n_cells_B))

message("Cloni espansi per paziente:")
print(table(cloni_espansi$patient))

# Cloni espansi in più pazienti indipendentemente
# (esclusi i cloni noti per essere dallo stesso lotto CAR-T)
cloni_multi_paz <- cloni_espansi %>%
  group_by(TRA_cdr3, TRB_cdr3) %>%
  filter(n_distinct(patient) > 1) %>%
  ungroup() %>%
  arrange(TRA_cdr3, patient)

message("\nCloni espansi in B in >1 paziente: ",
        n_distinct(paste(cloni_multi_paz$TRA_cdr3, cloni_multi_paz$TRB_cdr3)))

# ── STEP 4: Plot andamento I→A→B (line plot) ──────────────────────────────────
message("\n--- STEP 4: Plot espansione ---")

top_B_per_paziente <- clone_wide %>%
  group_by(patient) %>%
  slice_max(n_cells_B, n=10, with_ties=FALSE) %>%
  ungroup()

plot_dyn <- clone_dynamics %>%
  semi_join(top_B_per_paziente, by=c("patient","Clone_ID_CDR3")) %>%
  left_join(top_B_per_paziente %>%
              select(patient, Clone_ID_CDR3, n_cells_B, FC_I_to_B, categoria),
            by=c("patient","Clone_ID_CDR3"))

p_lines <- ggplot(plot_dyn,
                  aes(x     = stage,
                      y     = n_cells,
                      group = Clone_ID_CDR3,
                      color = categoria)) +
  geom_line(linewidth=1.2, alpha=0.8) +
  geom_point(size=3) +
  facet_wrap(~patient, scales="free_y", ncol=1) +
  scale_color_manual(values=c(
    "De novo in B"   = "#E41A1C",
    "Espanso (FC>=2)"= "#FF7F00",
    "Stabile"        = "#4DAF4A",
    "Contratto"      = "#377EB8",
    "Assente"        = "grey80",
    "Altro"          = "grey60"
  )) +
  scale_x_discrete(limits=c("I","A","B")) +
  theme_minimal() +
  theme(strip.text      = element_text(face="bold", size=12),
        legend.position = "bottom") +
  labs(
    title    = "Dinamica espansione clonale: stage I → A → B",
    subtitle = paste0("Top 10 cloni per stage B — solo cellule CAR+\n",
                      "Rosso=De novo, Arancio=Espanso, Verde=Stabile, Blu=Contratto"),
    x="Stage", y="N. cellule", color="Categoria"
  )

# ── STEP 5: Heatmap frequenza relativa ────────────────────────────────────────

top20_overall <- clone_wide %>%
  group_by(TRA_cdr3, TRB_cdr3) %>%
  summarise(n_tot_B = sum(n_cells_B), .groups="drop") %>%
  slice_max(n_tot_B, n=20, with_ties=FALSE)

heat_data <- clone_dynamics %>%
  semi_join(top20_overall, by=c("TRA_cdr3","TRB_cdr3")) %>%
  left_join(clone_wide %>% select(patient, Clone_ID_CDR3, categoria),
            by=c("patient","Clone_ID_CDR3")) %>%
  mutate(
    clone_label    = paste0(str_trunc(Gene_Label,20),
                            " | ", str_trunc(TRB_cdr3,16)),
    paziente_stage = factor(paste(patient, stage, sep="-"),
                            levels=c("Bo-I","Bo-A","Bo-B",
                                     "Ca-I","Ca-A","Ca-B",
                                     "Me-I","Me-B"))
  )

p_heat <- ggplot(heat_data,
                 aes(x    = paziente_stage,
                     y    = reorder(clone_label, freq),
                     fill = freq)) +
  geom_tile(color="white", linewidth=0.3) +
  geom_text(aes(label=ifelse(n_cells>0, n_cells, "")),
            size=3, color="white", fontface="bold") +
  scale_fill_gradientn(
    colors = c("white","#FEE8C8","#FC8D59","#D73027"),
    name   = "Freq.\nrelativa"
  ) +
  theme_minimal() +
  theme(axis.text.x  = element_text(angle=45, hjust=1, size=10),
        axis.text.y  = element_text(size=9),
        panel.grid   = element_blank()) +
  labs(
    title    = "Frequenza clonale relativa — Top 20 cloni in stage B",
    subtitle = "Numero nelle celle = cellule assolute",
    x="", y=""
  )

print(p_lines)
print(p_heat)

# ── STEP 6: Tabella per pubblicazione ─────────────────────────────────────────

tabella_pub <- clone_wide %>%
  filter(n_cells_B > 0) %>%
  transmute(
    Paziente          = patient,
    `V alpha`         = TRA_v_gene,
    `V beta`          = TRB_v_gene,
    `CDR3 alpha`      = TRA_cdr3,
    `CDR3 beta`       = TRB_cdr3,
    `Cellule stage I` = n_cells_I,
    `Cellule stage B` = n_cells_B,
    `Freq rel I`      = round(freq_I, 4),
    `Freq rel B`      = round(freq_B, 4),
    `Fold-change I→B` = round(FC_I_to_B, 2),
    `Categoria`       = categoria
  ) %>%
  arrange(Paziente, desc(`Cellule stage B`))

# ── STEP 7: Salvataggio ────────────────────────────────────────────────────────
message("\n--- STEP 7: Salvataggio ---")

ggsave(file.path(output_dir, "Expansion_dynamics_lineplot.png"),
       p_lines, width=12, height=14, dpi=300, bg="white")
ggsave(file.path(output_dir, "Expansion_heatmap.png"),
       p_heat,  width=14, height=10, dpi=300, bg="white")

write_xlsx(list(
  "01_Dinamica_completa"    = clone_wide,
  "02_Cloni_espansi_in_B"   = cloni_espansi,
  "03_Espansi_multi_paz"    = if(nrow(cloni_multi_paz)>0) cloni_multi_paz else
                                data.frame(nota="nessuno"),
  "04_Tabella_pubblicazione" = tabella_pub
), file.path(output_dir, "RISULTATI_expansion_dynamics.xlsx"))

message("Salvato in: ", output_dir)
message("\nFile prodotti:")
message("  Expansion_dynamics_lineplot.png")
message("  Expansion_heatmap.png")
message("  RISULTATI_expansion_dynamics.xlsx")

