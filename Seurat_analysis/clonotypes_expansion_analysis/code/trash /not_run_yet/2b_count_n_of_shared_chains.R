library(dplyr)
library(tidyr)
library(ggplot2)

# ==============================================================================
# SCRIPT 2B - VERSIONE CORRETTA
# FIX PRINCIPALI:
# 1. Distingue chiaramente tra condivisione V-gene e condivisione CDR3
# 2. Aggiunge metriche di diversità clonale
# 3. Calcola sovrapposizione reale tra repertori
# 4. Visualizzazioni multiple per interpretazione completa
# ==============================================================================

message("\n--- 🧬 ANALISI QUANTITATIVA CONDIVISIONE CLONALE ---\n")

# ==============================================================================
# 1. PREPARAZIONE DATI (come nello script 2)
# ==============================================================================

# Dataset per analisi V-gene
data_for_Vgene <- full_data %>%
  filter(Gene_Label != "? + ?") %>%
  separate_rows(TRA_V, sep = "/") %>%
  separate_rows(TRB_V, sep = "/") %>%
  filter(TRA_V != "" | TRB_V != "")

# Dataset per analisi cloni veri
data_for_Clones <- full_data %>%
  filter(Clone_Quality == "Complete")

# ==============================================================================
# 2. FUNZIONE HELPER MIGLIORATA
# ==============================================================================

analyze_sharing_advanced <- function(data, column_name, label, is_clone_level = FALSE) {
  
  message(paste0("\n[", label, "] Analisi in corso..."))
  
  # Prepara i dati unici per paziente
  unique_features <- data %>%
    rename(Feature = !!sym(column_name)) %>%
    filter(Feature != "?" & Feature != "" & !is.na(Feature)) %>%
    distinct(patient, Feature)
  
  # Conta in quanti pazienti appare ogni feature
  shared_counts <- unique_features %>%
    group_by(Feature) %>%
    summarise(
      n_Patients = n_distinct(patient),
      Patients_List = paste(sort(unique(patient)), collapse = ", "),
      .groups = 'drop'
    ) %>%
    arrange(desc(n_Patients), Feature)
  
  # Conta anche l'abbondanza (numero di cellule)
  abundance <- data %>%
    rename(Feature = !!sym(column_name)) %>%
    filter(Feature != "?" & Feature != "" & !is.na(Feature)) %>%
    group_by(Feature, patient) %>%
    summarise(n_cells = n(), .groups = 'drop')
  
  shared_with_abundance <- shared_counts %>%
    left_join(
      abundance %>% 
        group_by(Feature) %>% 
        summarise(Total_Cells = sum(n_cells), .groups = 'drop'),
      by = "Feature"
    )
  
  # Statistiche per livello di condivisione
  sharing_stats <- shared_counts %>%
    group_by(n_Patients) %>%
    summarise(
      Count = n(),
      Label = paste0(n_Patients, " paziente", ifelse(n_Patients > 1, "i", "")),
      .groups = 'drop'
    )
  
  message(paste0("  • Totale elementi unici: ", nrow(shared_counts)))
  message(paste0("  • Condivisi da 3 pazienti: ", nrow(shared_counts %>% filter(n_Patients == 3))))
  message(paste0("  • Condivisi da 2 pazienti: ", nrow(shared_counts %>% filter(n_Patients == 2))))
  message(paste0("  • Privati (1 paziente): ", nrow(shared_counts %>% filter(n_Patients == 1))))
  
  # Mostra esempi dei più abbondanti condivisi da 3
  if(nrow(shared_counts %>% filter(n_Patients == 3)) > 0) {
    message("\n  🎯 Top 5 elementi condivisi da 3 pazienti (per abbondanza):")
    top_shared <- shared_with_abundance %>% 
      filter(n_Patients == 3) %>%
      arrange(desc(Total_Cells)) %>%
      head(5)
    print(top_shared)
  }
  
  return(list(
    shared_counts = shared_counts,
    sharing_stats = sharing_stats,
    abundance = shared_with_abundance
  ))
}

# ==============================================================================
# 3. ANALISI COPPIA ESATTA (V-gene level)
# ==============================================================================
message("\n🔹 ANALISI 1: COPPIE V-GENE (TRA_V + TRB_V)")
message("Queste sono combinazioni di V-genes, NON cloni identici\n")

# Ricrea Gene_Label dai V-genes separati
data_vgene_pairs <- data_for_Vgene %>%
  mutate(VGene_Pair = paste(
    ifelse(TRA_V == "" | is.na(TRA_V), "?", TRA_V),
    "+",
    ifelse(TRB_V == "" | is.na(TRB_V), "?", TRB_V)
  ))

results_vgene_pairs <- analyze_sharing_advanced(
  data_vgene_pairs, 
  "VGene_Pair", 
  "Coppie V-Gene"
)

# ==============================================================================
# 4. ANALISI SOLO CATENA ALPHA (V-gene)
# ==============================================================================
message("\n🔹 ANALISI 2: SOLO CATENA ALPHA (TRA_V)")

results_trav <- analyze_sharing_advanced(
  data_for_Vgene, 
  "TRA_V", 
  "Solo TRA-V"
)

# ==============================================================================
# 5. ANALISI SOLO CATENA BETA (V-gene)
# ==============================================================================
message("\n🔹 ANALISI 3: SOLO CATENA BETA (TRB_V)")

results_trbv <- analyze_sharing_advanced(
  data_for_Vgene, 
  "TRB_V", 
  "Solo TRB-V"
)

# ==============================================================================
# 6. ANALISI VERI CLONI (CDR3-based) ⭐ LA PIÙ IMPORTANTE
# ==============================================================================
message("\n🔹 ANALISI 4: VERI CLONI (basato su Clone_ID_CDR3)")
message("Questi sono cloni REALMENTE identici\n")

results_clones <- analyze_sharing_advanced(
  data_for_Clones, 
  "Clone_ID_CDR3", 
  "Cloni CDR3",
  is_clone_level = TRUE
)

# ==============================================================================
# 7. PLOT 1: OVERVIEW CONDIVISIONE (tutti i livelli)
# ==============================================================================
message("\n--- 📊 Generazione Plot 1: Overview Condivisione ---")

plot_sharing_all <- bind_rows(
  results_vgene_pairs$sharing_stats %>% mutate(Type = "Coppia V-gene (A+B)"),
  results_trav$sharing_stats %>% mutate(Type = "Solo TRAV"),
  results_trbv$sharing_stats %>% mutate(Type = "Solo TRBV"),
  results_clones$sharing_stats %>% mutate(Type = "⭐ Veri Cloni (CDR3)")
) %>%
  mutate(
    Type = factor(Type, levels = c(
      "⭐ Veri Cloni (CDR3)", 
      "Coppia V-gene (A+B)", 
      "Solo TRBV", 
      "Solo TRAV"
    )),
    n_Patients = factor(n_Patients, levels = c(3, 2, 1)),
    Label = factor(Label)
  )

p1 <- ggplot(plot_sharing_all, aes(x = n_Patients, y = Count, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", size = 0.3) +
  geom_text(aes(label = Count), position = position_dodge(width = 0.9), 
            vjust = -0.5, size = 3) +
  scale_fill_manual(values = c(
    "⭐ Veri Cloni (CDR3)" = "#D55E00",  # Arancione acceso
    "Coppia V-gene (A+B)" = "#E7B800",
    "Solo TRBV" = "#00AFBB",
    "Solo TRAV" = "#FC4E07"
  )) +
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 11, face = "bold")
  ) +
  labs(
    title = "Livelli di Condivisione Clonale tra Pazienti",
    subtitle = "Confronto tra condivisione V-gene (convergenza) e CDR3 (veri public clones)",
    x = "Numero di Pazienti che condividono l'elemento",
    y = "Numero di Elementi Unici",
    fill = "Tipo di Analisi"
  ) +
  scale_y_log10()  # Log scale perché i numeri variano molto

print(p1)

# ==============================================================================
# 8. PLOT 2: FOCUS SU ELEMENTI CONDIVISI DA 3 PAZIENTI
# ==============================================================================
message("\n--- 📊 Generazione Plot 2: Focus Condivisione Completa (3 pazienti) ---")

# Conta solo gli elementi condivisi da tutti e 3
shared_by_3 <- data.frame(
  Type = c("Coppia V-gene", "TRAV", "TRBV", "⭐ Veri Cloni CDR3"),
  Count = c(
    nrow(results_vgene_pairs$shared_counts %>% filter(n_Patients == 3)),
    nrow(results_trav$shared_counts %>% filter(n_Patients == 3)),
    nrow(results_trbv$shared_counts %>% filter(n_Patients == 3)),
    nrow(results_clones$shared_counts %>% filter(n_Patients == 3))
  )
) %>%
  mutate(Type = factor(Type, levels = rev(c("Coppia V-gene", "TRAV", "TRBV", "⭐ Veri Cloni CDR3"))))

p2 <- ggplot(shared_by_3, aes(x = Type, y = Count, fill = Type)) +
  geom_bar(stat = "identity", color = "black", size = 0.3) +
  geom_text(aes(label = Count), hjust = -0.3, size = 5, fontface = "bold") +
  coord_flip() +
  scale_fill_manual(values = c(
    "⭐ Veri Cloni CDR3" = "#D55E00",
    "TRBV" = "#00AFBB",
    "TRAV" = "#FC4E07",
    "Coppia V-gene" = "#E7B800"
  )) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 12, face = "bold")
  ) +
  labs(
    title = "Convergenza Pubblica: Elementi Condivisi da TUTTI e 3 i Pazienti",
    subtitle = "Più alto il numero, più forte è la convergenza a quel livello",
    x = "",
    y = "Numero di Elementi Condivisi"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

print(p2)

# ==============================================================================
# 9. ANALISI DIVERSITÀ CLONALE (Shannon Diversity, Clonality)
# ==============================================================================
message("\n--- 📊 CALCOLO DIVERSITÀ CLONALE ---")

# Funzione per calcolare diversità Shannon
calc_shannon <- function(freqs) {
  freqs <- freqs[freqs > 0]
  p <- freqs / sum(freqs)
  -sum(p * log(p))
}

# Funzione per calcolare clonalità (1 - Shannon normalizzato)
calc_clonality <- function(freqs) {
  shannon <- calc_shannon(freqs)
  shannon_max <- log(length(freqs))
  if(shannon_max == 0) return(1)
  1 - (shannon / shannon_max)
}

diversity_stats <- data_for_Clones %>%
  group_by(patient, stage) %>%
  summarise(
    n_unique_clones = n_distinct(Clone_ID_CDR3),
    total_cells = n(),
    shannon_diversity = calc_shannon(table(Clone_ID_CDR3)),
    clonality = calc_clonality(table(Clone_ID_CDR3)),
    .groups = 'drop'
  ) %>%
  arrange(patient, stage)

message("\n📊 Diversità Clonale per Paziente e Stage:")
print(diversity_stats)

# ==============================================================================
# 10. PLOT 3: DIVERSITÀ CLONALE
# ==============================================================================
message("\n--- 📊 Generazione Plot 3: Diversità Clonale ---")

p3 <- ggplot(diversity_stats, aes(x = stage, y = shannon_diversity, fill = stage)) +
  geom_bar(stat = "identity", color = "black", size = 0.3) +
  geom_text(aes(label = round(shannon_diversity, 2)), vjust = -0.5, size = 3) +
  facet_wrap(~patient, scales = "free_x") +
  scale_fill_manual(values = c("I"="#619CFF", "A"="#F8766D", "B"="#00BA38")) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 11)
  ) +
  labs(
    title = "Diversità Clonale (Shannon Diversity Index)",
    subtitle = "Più alto = repertorio più diversificato. Più basso = espansione oligoclonale.",
    x = "Stage",
    y = "Shannon Diversity",
    fill = "Stage"
  )

print(p3)

# ==============================================================================
# 11. PLOT 4: CLONALITÀ (inverso della diversità)
# ==============================================================================
p4 <- ggplot(diversity_stats, aes(x = stage, y = clonality, fill = stage)) +
  geom_bar(stat = "identity", color = "black", size = 0.3) +
  geom_text(aes(label = round(clonality, 3)), vjust = -0.5, size = 3) +
  facet_wrap(~patient, scales = "free_x") +
  scale_fill_manual(values = c("I"="#619CFF", "A"="#F8766D", "B"="#00BA38")) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 11)
  ) +
  labs(
    title = "Clonalità (Espansione Oligoclonale)",
    subtitle = "Più alto = pochi cloni dominano. Valore ideale per CAR-T = 0.5-0.8",
    x = "Stage",
    y = "Clonality Index (0-1)",
    fill = "Stage"
  ) +
  scale_y_continuous(limits = c(0, 1))

print(p4)

# ==============================================================================
# 12. ANALISI OVERLAP INDEX (Jaccard) TRA PAZIENTI
# ==============================================================================
message("\n--- 📊 CALCOLO OVERLAP TRA REPERTORI ---")

# Funzione per calcolare Jaccard Index
calc_jaccard <- function(set1, set2) {
  intersection <- length(intersect(set1, set2))
  union <- length(union(set1, set2))
  if(union == 0) return(0)
  intersection / union
}

# Repertori per paziente (solo cloni completi)
repertoires <- data_for_Clones %>%
  group_by(patient) %>%
  summarise(clones = list(unique(Clone_ID_CDR3)), .groups = 'drop')

# Calcola Jaccard per tutte le coppie
patients <- repertoires$patient
overlap_matrix <- matrix(0, nrow = length(patients), ncol = length(patients))
rownames(overlap_matrix) <- patients
colnames(overlap_matrix) <- patients

for(i in 1:length(patients)) {
  for(j in 1:length(patients)) {
    overlap_matrix[i, j] <- calc_jaccard(
      repertoires$clones[[i]], 
      repertoires$clones[[j]]
    )
  }
}

message("\n📊 Jaccard Index tra Repertori (0=nessuna sovrapposizione, 1=identici):")
print(round(overlap_matrix, 4))

# Converti in formato long per plot
overlap_long <- as.data.frame(overlap_matrix) %>%
  mutate(Patient1 = rownames(.)) %>%
  pivot_longer(cols = -Patient1, names_to = "Patient2", values_to = "Jaccard") %>%
  mutate(
    Patient1 = factor(Patient1, levels = patients),
    Patient2 = factor(Patient2, levels = patients)
  )

# ==============================================================================
# 13. PLOT 5: HEATMAP OVERLAP
# ==============================================================================
p5 <- ggplot(overlap_long, aes(x = Patient1, y = Patient2, fill = Jaccard)) +
  geom_tile(color = "white", size = 1) +
  geom_text(aes(label = round(Jaccard, 3)), size = 5, fontface = "bold") +
  scale_fill_gradient2(
    low = "white", 
    mid = "#FFF4B3", 
    high = "#D55E00",
    midpoint = 0.05,
    limits = c(0, max(overlap_long$Jaccard))
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 12, face = "bold"),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Overlap tra Repertori TCR (Jaccard Index)",
    subtitle = "Basato su CDR3 identici. Valori bassi = repertori privati, alti = condivisione",
    x = "",
    y = "",
    fill = "Jaccard\nIndex"
  ) +
  coord_fixed()

print(p5)

# ==============================================================================
# 14. SUMMARY REPORT
# ==============================================================================
message("\n" , paste(rep("=", 70), collapse = ""))
message("📊 REPORT FINALE - CONDIVISIONE CLONALE")
message(paste(rep("=", 70), collapse = ""))

# V-gene sharing
message("\n🔸 LIVELLO V-GENE (Convergenza):")
message(sprintf("  • TRBV condivisi da 3 pazienti: %d", 
                nrow(results_trbv$shared_counts %>% filter(n_Patients == 3))))
message(sprintf("  • TRAV condivisi da 3 pazienti: %d", 
                nrow(results_trav$shared_counts %>% filter(n_Patients == 3))))

# Clone sharing
n_public_clones_3 <- nrow(results_clones$shared_counts %>% filter(n_Patients == 3))
n_public_clones_2 <- nrow(results_clones$shared_counts %>% filter(n_Patients == 2))

message("\n⭐ LIVELLO CLONE (Veri Public Clones - CDR3):")
message(sprintf("  • Cloni identici condivisi da 3 pazienti: %d", n_public_clones_3))
message(sprintf("  • Cloni identici condivisi da 2 pazienti: %d", n_public_clones_2))

# Diversity
message("\n🔸 DIVERSITÀ CLONALE:")
for(p in unique(diversity_stats$patient)) {
  pat_stats <- diversity_stats %>% filter(patient == p)
  message(sprintf("  • %s: Shannon range = %.2f-%.2f, Clonality range = %.3f-%.3f",
                  p,
                  min(pat_stats$shannon_diversity),
                  max(pat_stats$shannon_diversity),
                  min(pat_stats$clonality),
                  max(pat_stats$clonality)))
}

# Overlap
avg_overlap <- mean(overlap_matrix[upper.tri(overlap_matrix)])
message(sprintf("\n🔸 OVERLAP REPERTORI (Jaccard medio tra pazienti): %.4f", avg_overlap))

message("\n" , paste(rep("=", 70), collapse = ""))

# ==============================================================================
# 15. SALVATAGGIO TUTTI I RISULTATI
# ==============================================================================
message("\n--- 💾 Salvataggio risultati... ---")

write.csv(results_vgene_pairs$shared_counts, "RISULTATI_Sharing_VgenePairs.csv", row.names = FALSE)
write.csv(results_trbv$shared_counts, "RISULTATI_Sharing_TRBV.csv", row.names = FALSE)
write.csv(results_trav$shared_counts, "RISULTATI_Sharing_TRAV.csv", row.names = FALSE)
write.csv(results_clones$shared_counts, "RISULTATI_Sharing_Clones_CDR3.csv", row.names = FALSE)
write.csv(diversity_stats, "RISULTATI_Diversity_Analysis.csv", row.names = FALSE)
write.csv(overlap_long, "RISULTATI_Jaccard_Overlap.csv", row.names = FALSE)

message("✅ Script completato! Tutti i file salvati.")
message("\n💡 INTERPRETAZIONE CHIAVE:")
message("   1. V-gene condivisi = Convergenza verso segmenti genici pubblici")
message("   2. CDR3 identici = Veri public clones (rarissimi, molto significativi)")
message("   3. Shannon alto = repertorio diverso, Clonality alto = espansione oligoclonale")
message("   4. Jaccard basso = repertori privati (normale), alto = forte condivisione")
