library(dplyr)
library(ggplot2)
library(readr)
library(writexl)
library(stringr)
library(xlsx)

# 1. CARICAMENTO DATI
# Assicurati che i file siano nella tua directory di lavoro
full_data <- read.xlsx("RISULTATI_Cloni_Dati_Completi_con_CDR3.xlsx", sheetIndex = 1)

# 2. IDENTIFICAZIONE CONTAMINAZIONE (Stesso DNA tra pazienti)
# Definiamo i cloni che hanno mostrato identità nucleotidica al 100% tra pazienti diversi
contaminated_clones <- full_data %>%
  filter(Clone_Quality == "Complete") %>%
  group_by(TRA_cdr3, TRB_cdr3, TRA_cdr3_nt, TRB_cdr3_nt) %>%
  filter(n_distinct(patient) > 1) %>%
  summarise(
    pazienti_coinvolti = paste(unique(patient), collapse = " & "),
    cellule_totali = n(),
    gene_label = first(Gene_Label),
    .groups = "drop"
  )

# Salvataggio Report Contaminazione
write_xlsx(contaminated_clones, "REPORT_CONTAMINAZIONE_DETTAGLIATO.xlsx")

# 3. FILTRAGGIO PER "FAMIGLIE D'ORO" (Solo cloni privati)
# Escludiamo i cloni contaminati per vedere la vera fitness biologica
shared_keys <- paste0(contaminated_clones$TRA_cdr3, "_", contaminated_clones$TRB_cdr3)

private_data <- full_data %>%
  filter(Clone_Quality == "Complete") %>%
  mutate(temp_key = paste0(TRA_cdr3, "_", TRB_cdr3)) %>%
  filter(!(temp_key %in% shared_keys))

# Selezioniamo le famiglie target (TRBV7-8, TRBV2, TRBV5-1)
target_families <- c("TRBV7-8", "TRBV2", "TRBV5-1", "TRBV28")

plot_data <- private_data %>%
  filter(TRB_v_gene %in% target_families) %>%
  group_by(patient, TRB_v_gene, TRA_v_gene, TRA_cdr3, TRB_cdr3) %>%
  summarise(espansione = n(), .groups = "drop") %>%
  arrange(desc(espansione)) %>%
  group_by(patient) %>%
  slice_head(n = 5) # Prendiamo i top 5 privati per paziente in queste famiglie

# 4. GENERAZIONE GRAFICO
p <- ggplot(plot_data, aes(x = reorder(paste(TRA_v_gene, TRB_v_gene, sep=" + "), espansione), 
                           y = espansione, fill = TRB_v_gene)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  facet_wrap(~patient, scales = "free_y", ncol = 1) +
  coord_flip() +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() +
  labs(
    title = "Famiglie TCR 'Gold' - Espansione Reale (Privata)",
    subtitle = "Esclusi i cloni con identità nucleotidica tra pazienti (contaminazione)",
    x = "Coppia Alpha + Beta",
    y = "Numero di Cellule (Espansione)",
    fill = "Famiglia Beta"
  ) +
  theme(strip.text = element_text(face = "bold", size = 12))

# Salvataggio Grafico
ggsave("Grafico_Famiglie_Gold_Privato.png", p, width = 10, height = 8, dpi = 300)

message("Analisi completata. Generati: 'REPORT_CONTAMINAZIONE_DETTAGLIATO.xlsx' e 'Grafico_Famiglie_Gold_Privato.png'")
