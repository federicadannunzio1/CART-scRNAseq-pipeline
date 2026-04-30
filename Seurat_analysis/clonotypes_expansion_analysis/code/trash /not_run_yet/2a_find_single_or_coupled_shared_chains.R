library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(writexl)
# ==============================================================================
# SCRIPT 2 - VERSIONE CORRETTA
# FIX PRINCIPALI:
# 1. Distingue tra "condivisione V-gene" e "condivisione clone vero (CDR3)"
# 2. Non separa più le righe con catene multiple (evita conta duplicata)
# 3. Analizza sia convergenza di V-gene (pubblici) che veri public clones
# 4. Gestisce meglio le cellule con dati mancanti
# ==============================================================================

message("\n--- 🧬 ANALISI CONDIVISIONE: V-GENE vs VERI CLONI ---\n")
output_dir <- "/Users/federicadannunzio/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/caruana-project/CART/Code/Seurat_analysis/4_clonotypes_expansion_analysis/res"

# ==============================================================================
# 0. PREPARAZIONE DATI
# ==============================================================================

# Creiamo due dataset separati per le due analisi:

# DATASET 1: Per analisi V-gene (permette catene multiple)
# Qui separiamo le catene multiple perché ci interessa vedere se lo stesso
# V-gene appare in pazienti diversi, anche se parte di cloni diversi
data_for_Vgene <- full_data %>%
  filter(Gene_Label != "? + ?") %>%
  # Separiamo le catene multiple (es. "TRAV1/TRAV2" diventa 2 righe)
  separate_rows(TRA_V, sep = "/") %>%
  separate_rows(TRB_V, sep = "/") %>%
  # Rimuoviamo valori vuoti
  filter(TRA_V != "" | TRB_V != "")

# DATASET 2: Per analisi veri cloni (NO separazione, cellule intere)
# Qui teniamo le cellule intatte perché ci interessa vedere se lo STESSO
# clone (stesso CDR3) è condiviso tra pazienti
data_for_Clones <- full_data %>%
  filter(Clone_Quality == "Complete") # Solo cloni con dati completi
  # NON separiamo le righe - ogni riga = 1 cellula = 1 clone

# ==============================================================================
# 1. ANALISI CONDIVISIONE V-GENE (CONVERGENZA PUBBLICA)
# ==============================================================================
message("📊 PARTE 1: CONVERGENZA V-GENE (public TCR segments)")
message("Questa analisi mostra V-genes usati da pazienti diversi")
message("⚠️ NOTA: Stesso V-gene ≠ stesso clone!\n")

# Funzione per trovare elementi condivisi da N pazienti
get_shared_elements <- function(data, col_name, min_patients = 3) {
  data %>%
    filter(!!sym(col_name) != "?" & !is.na(!!sym(col_name)) & !!sym(col_name) != "") %>%
    distinct(patient, !!sym(col_name)) %>%
    group_by(!!sym(col_name)) %>%
    summarise(
      n_patients = n_distinct(patient),
      patients_list = paste(sort(unique(patient)), collapse = ", "),
      .groups = 'drop'
    ) %>%
    filter(n_patients >= min_patients) %>%
    arrange(desc(n_patients), !!sym(col_name))
}

# Trova V-genes condivisi
shared_TRA_V <- get_shared_elements(data_for_Vgene, "TRA_V", min_patients = 3)
shared_TRB_V <- get_shared_elements(data_for_Vgene, "TRB_V", min_patients = 3)

message(paste("✅ Trovati", nrow(shared_TRA_V), "geni TRAV condivisi da 3 pazienti"))
message(paste("✅ Trovati", nrow(shared_TRB_V), "geni TRBV condivisi da 3 pazienti"))

# ==============================================================================
# 2. ANALISI VERI PUBLIC CLONES (CONDIVISIONE CDR3)
# ==============================================================================
message("\n📊 PARTE 2: VERI PUBLIC CLONES (stesso CDR3)")
message("Questi sono cloni IDENTICI condivisi tra pazienti\n")

# Trova cloni completi identici condivisi
shared_Clones <- data_for_Clones %>%
  group_by(Clone_ID_CDR3, Gene_Label, TRA_CDR3, TRB_CDR3) %>%
  summarise(
    n_patients = n_distinct(patient),
    patients_list = paste(sort(unique(patient)), collapse = ", "),
    total_cells = n(),
    .groups = 'drop'
  ) %>%
  filter(n_patients >= 2) %>%  # Almeno 2 pazienti (3 pazienti è rarissimo)
  arrange(desc(n_patients), desc(total_cells))

message(paste("✅ Trovati", nrow(shared_Clones %>% filter(n_patients == 3)), 
              "cloni identici condivisi da 3 pazienti"))
message(paste("✅ Trovati", nrow(shared_Clones %>% filter(n_patients == 2)), 
              "cloni identici condivisi da 2 pazienti"))

if(nrow(shared_Clones %>% filter(n_patients == 3)) > 0) {
  message("\n🎯 CLONI PUBBLICI (3 pazienti):")
  print(shared_Clones %>% filter(n_patients == 3) %>% 
          select(Gene_Label, TRB_CDR3, patients_list, total_cells))
}

# ==============================================================================
# 3. PLOT 1: CATENE TRB-V CONDIVISE (CONVERGENZA PUBBLICA)
# ==============================================================================
message("\n--- 📊 Generazione Plot 1: TRB-V condivisi ---")

if(nrow(shared_TRB_V) > 0) {
  
  # Prepara dati per il plot: conta cellule per ogni TRB-V condiviso
  plot_data_TRBV <- data_for_Vgene %>%
    filter(TRB_V %in% shared_TRB_V[[1]]) %>%  # Prima colonna contiene i nomi TRB_V
    group_by(patient, stage, TRB_V) %>%
    summarise(total_cells = n(), .groups = 'drop') %>%
    # Completa con zeri per stage mancanti
    complete(TRB_V, nesting(patient, stage), fill = list(total_cells = 0)) %>%
    mutate(stage = factor(stage, levels = c("I", "A", "B")))
  
  p1 <- ggplot(plot_data_TRBV, aes(x = patient, y = total_cells, fill = stage)) +
    geom_bar(stat = "identity", position = "dodge", color="black", size=0.2) +
    facet_wrap(~TRB_V, scales = "free_y", ncol = 4) +
    scale_fill_manual(values = c("I"="#619CFF", "A"="#F8766D", "B"="#00BA38")) +
    theme_bw() +
    theme(
      strip.text = element_text(face = "bold", size = 9, color = "white"),
      strip.background = element_rect(fill = "#0073C2"),
      axis.text.x = element_text(angle = 0),
      legend.position = "bottom"
    ) +
    labs(
      title = "Catene Beta (TRBV) Condivise - CONVERGENZA PUBBLICA",
      subtitle = "⚠️ Stesso V-gene ≠ Stesso clone! Questi sono FAMIGLIE di cloni convergenti",
      y = "N. Cellule (somma di tutti i cloni con questo V-gene)", 
      x = "Paziente",
      fill = "Stage"
    )
  
  print(p1)
  
} else {
  message("⚠️ Nessuna catena Beta condivisa da tutti e 3 i pazienti.")
}

# ==============================================================================
# 4. PLOT 2: CATENE TRA-V CONDIVISE
# ==============================================================================
message("\n--- 📊 Generazione Plot 2: TRA-V condivisi ---")

if(nrow(shared_TRA_V) > 0) {
  
  plot_data_TRAV <- data_for_Vgene %>%
    filter(TRA_V %in% shared_TRA_V[[1]]) %>%
    group_by(patient, stage, TRA_V) %>%
    summarise(total_cells = n(), .groups = 'drop') %>%
    complete(TRA_V, nesting(patient, stage), fill = list(total_cells = 0)) %>%
    mutate(stage = factor(stage, levels = c("I", "A", "B")))
  
  p2 <- ggplot(plot_data_TRAV, aes(x = patient, y = total_cells, fill = stage)) +
    geom_bar(stat = "identity", position = "dodge", color="black", size=0.2) +
    facet_wrap(~TRA_V, scales = "free_y", ncol = 4) +
    scale_fill_manual(values = c("I"="#619CFF", "A"="#F8766D", "B"="#00BA38")) +
    theme_bw() +
    theme(
      strip.text = element_text(face = "bold", size = 9),
      axis.text.x = element_text(angle = 0),
      legend.position = "bottom"
    ) +
    labs(
      title = "Catene Alpha (TRAV) Condivise",
      subtitle = "V-genes Alpha usati da tutti e 3 i pazienti",
      y = "N. Cellule", 
      x = "Paziente"
    )
  
  print(p2)
  
} else {
  message("⚠️ Nessuna catena Alpha condivisa da tutti e 3 i pazienti.")
}

# ==============================================================================
# 5. PLOT 3: VERI PUBLIC CLONES (COPPIE CDR3 IDENTICHE)
# ==============================================================================
message("\n--- 📊 Generazione Plot 3: Veri Public Clones (CDR3 identici) ---")

if(nrow(shared_Clones) > 0) {
  
  # Prepara dati: conta cellule per ogni clone condiviso
  plot_data_clones <- data_for_Clones %>%
    inner_join(
      shared_Clones %>% select(Clone_ID_CDR3, n_patients),
      by = "Clone_ID_CDR3"
    ) %>%
    group_by(patient, stage, Clone_ID_CDR3, Gene_Label, TRB_CDR3, n_patients) %>%
    summarise(total_cells = n(), .groups = 'drop') %>%
    complete(
      nesting(Clone_ID_CDR3, Gene_Label, TRB_CDR3, n_patients), 
      nesting(patient, stage), 
      fill = list(total_cells = 0)
    ) %>%
    mutate(
      stage = factor(stage, levels = c("I", "A", "B")),
      # Label descrittiva per il plot
      Clone_Label = paste0(Gene_Label, "\n", substr(TRB_CDR3, 1, 10), "..."),
      # Ordina per numero di pazienti e poi per abbondanza
      Clone_Label = factor(Clone_Label)
    )
  
  # Ordina i cloni per rilevanza (prima quelli condivisi da 3, poi 2 pazienti)
  clone_order <- plot_data_clones %>%
    group_by(Clone_Label, n_patients) %>%
    summarise(total = sum(total_cells), .groups = 'drop') %>%
    arrange(desc(n_patients), desc(total)) %>%
    pull(Clone_Label)
  
  plot_data_clones <- plot_data_clones %>%
    mutate(Clone_Label = factor(Clone_Label, levels = clone_order))
  
  p3 <- ggplot(plot_data_clones, aes(x = patient, y = total_cells, fill = stage)) +
    geom_bar(stat = "identity", position = "dodge", color="black", size=0.2) +
    facet_wrap(~Clone_Label, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = c("I"="#619CFF", "A"="#F8766D", "B"="#00BA38")) +
    theme_bw() +
    theme(
      strip.text = element_text(face = "bold", size = 8, color = "white"),
      strip.background = element_rect(fill = "#D55E00"),  # Arancione per evidenziare
      axis.text.x = element_text(angle = 0),
      legend.position = "bottom",
      panel.spacing = unit(1, "lines")
    ) +
    labs(
      title = "VERI PUBLIC CLONES - CDR3 IDENTICI tra Pazienti ✅",
      subtitle = "Questi sono cloni REALMENTE identici (stesso V+J+CDR3). Rarissimi ma biologicamente significativi!",
      y = "N. Cellule", 
      x = "Paziente",
      fill = "Stage"
    )
  
  print(p3)
  
  # Salva anche la lista dei public clones
  write.csv(shared_Clones, "RISULTATI_Public_Clones_CDR3_identici.csv", row.names = FALSE)
  
} else {
  message("⚠️ Nessun clone identico (CDR3) condiviso da 2+ pazienti.")
  message("   Questo è normale: i veri public clones sono rarissimi!")
}

# ==============================================================================
# 6. ANALISI STATISTICA: QUANTI CLONI UNICI vs CONDIVISI?
# ==============================================================================
message("\n--- 📊 STATISTICHE FINALI ---\n")

# Per ogni paziente, quanti dei loro top cloni usano V-genes pubblici?
stats_vgene <- data_for_Vgene %>%
  group_by(patient) %>%
  summarise(
    Total_Cells = n(),
    Cells_with_Public_TRBV = sum(TRB_V %in% shared_TRB_V[[1]]),
    Percent_Public_TRBV = round(100 * Cells_with_Public_TRBV / Total_Cells, 1),
    .groups = 'drop'
  )

message("📊 Percentuale cellule che usano TRBV pubblici (condivisi da 3 pazienti):")
print(stats_vgene)

# Cloni unici vs condivisi
stats_clones <- data_for_Clones %>%
  group_by(patient) %>%
  summarise(
    Total_Unique_Clones = n_distinct(Clone_ID_CDR3),
    Total_Cells = n(),
    .groups = 'drop'
  )

message("\n📊 Numero di cloni unici per paziente:")
print(stats_clones)

# ==============================================================================
# 7. SALVATAGGIO RISULTATI
# ==============================================================================
write_xlsx(shared_TRB_V, paste0(output_dir, "/RISULTATI_TRBV_condivisi_3pazienti.xlsx"))

write_xlsx(shared_TRA_V, paste0(output_dir, "/RISULTATI_TRAV_condivisi_3pazienti.xlsx"))

message("\n✅ Script completato!")
message("\n🎯 INTERPRETAZIONE:")
message("   - V-genes condivisi = CONVERGENZA (pazienti usano stessi segmenti genetici)")
message("   - CDR3 identici = VERI PUBLIC CLONES (cloni esattamente identici, rarissimi)")
message("\n⚠️ Per analisi biologiche, distingui sempre tra i due livelli!")
