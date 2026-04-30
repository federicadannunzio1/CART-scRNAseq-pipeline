library(dplyr)
library(readr)
library(readxl)
library(utils) # Per decomprimere lo zip

# ==============================================================================
# 1. CONFIGURAZIONE
# ==============================================================================
# Il tuo file con i cloni (VERIFICA CHE IL PERCORSO SIA CORRETTO)
file_path <- "/Users/federicadannunzio/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/caruana-project/CART/Code/Seurat_analysis/4_clonotypes_expansion_analysis/res/RISULTATI_Top10_Cloni_CDR3.xlsx"

if(!file.exists(file_path)) stop("❌ File Excel dei cloni non trovato! Controlla il percorso.")

# ==============================================================================
# 2. SCARICAMENTO E LETTURA VDJdb (AUTOMATICO)
# ==============================================================================
message("--- 1. Scaricamento Database VDJdb (Versione 2024-06-13) ---")

url_zip <- "https://github.com/antigenomics/vdjdb-db/releases/download/2024-06-13/vdjdb-2024-06-13.zip"
temp_zip <- tempfile(fileext = ".zip")
temp_dir <- tempdir()

# Scarica
download.file(url_zip, temp_zip, mode = "wb")

# Estrae
message("--- 2. Estrazione e Lettura ---")
unzip(temp_zip, exdir = temp_dir)
target_file <- list.files(temp_dir, pattern = "vdjdb.*\\.txt|vdjdb.*\\.tsv", full.names = TRUE, recursive = TRUE)[1]

# Legge il file grezzo
vdjdb_raw <- read_tsv(target_file, show_col_types = FALSE)

# ==============================================================================
# 3. PREPARAZIONE DATI (Separazione Alpha/Beta per evitare errori)
# ==============================================================================
message("--- 3. Ristrutturazione Database (Alpha vs Beta) ---")

# Creiamo il DB di riferimento per la BETA (usando cdr3.beta)
vdjdb_beta <- vdjdb_raw %>%
  filter(species == "HomoSapiens") %>%
  select(cdr3.beta, antigen.species, antigen.gene, antigen.epitope) %>%
  rename(
    CDR3_Sequence = cdr3.beta,
    Virus = antigen.species,
    Target_Protein = antigen.gene,
    Epitope_Seq = antigen.epitope
  ) %>%
  filter(!is.na(CDR3_Sequence) & CDR3_Sequence != "") %>%
  distinct()

# Creiamo il DB di riferimento per la ALPHA (usando cdr3.alpha)
vdjdb_alpha <- vdjdb_raw %>%
  filter(species == "HomoSapiens") %>%
  select(cdr3.alpha, antigen.species, antigen.gene, antigen.epitope) %>%
  rename(
    CDR3_Sequence = cdr3.alpha,
    Virus = antigen.species,
    Target_Protein = antigen.gene,
    Epitope_Seq = antigen.epitope
  ) %>%
  filter(!is.na(CDR3_Sequence) & CDR3_Sequence != "") %>%
  distinct()

# ==============================================================================
# 4. CARICAMENTO TUOI CLONI E MATCH
# ==============================================================================
message("--- 4. Ricerca Match con i tuoi cloni ---")

my_clones <- read_excel(file_path)

# Match sulla Beta
matches_beta <- my_clones %>%
  inner_join(vdjdb_beta, by = c("TRB_CDR3" = "CDR3_Sequence"))

# Match sulla Alpha
matches_alpha <- my_clones %>%
  inner_join(vdjdb_alpha, by = c("TRA_CDR3" = "CDR3_Sequence"))

# ==============================================================================
# 5. RISULTATI E OUTPUT
# ==============================================================================
message("\n==================================================")
message("               RISULTATI RICERCA                  ")
message("==================================================")

found_something <- FALSE

# Risultati BETA
if(nrow(matches_beta) > 0) {
  found_something <- TRUE
  message("✅ TROVATO MATCH PER LA CATENA BETA!")
  print(matches_beta %>% select(patient, TRB_CDR3, Virus, Target_Protein, Epitope_Seq))
  write.csv(matches_beta, file.path(dirname(file_path), "VDJdb_Match_Beta.csv"))
} else {
  message("❌ Nessun match esatto per la Beta.")
}

# Risultati ALPHA
if(nrow(matches_alpha) > 0) {
  found_something <- TRUE
  message("✅ TROVATO MATCH PER LA CATENA ALPHA!")
  print(matches_alpha %>% select(patient, TRA_CDR3, Virus, Target_Protein, Epitope_Seq))
  write.csv(matches_alpha, file.path(dirname(file_path), "VDJdb_Match_Alpha.csv"))
} else {
  message("❌ Nessun match esatto per la Alpha.")
}

if(!found_something) {
  message("\n--------------------------------------------------")
  message("CONCLUSIONE: I tuoi cloni sono 'Orfani'.")
  message("Non sono presenti nel database pubblico VDJdb (versione Giugno 2024).")
  message("Questo indica che sono cloni privati o specifici per l'HLA dei pazienti.")
}

# Pulizia
unlink(temp_zip)