library(dplyr)
library(readr)
library(readxl)
library(writexl)
library(utils)

# ==============================================================================
# 02_vdjdb_search.R
#
# COSA FA:
#   Scarica VDJdb (release 2024-06-13) e cerca le CDR3 dei tuoi cloni.
#   Cerca sia su catena alpha che beta.
#   Output: specificità antigenica nota (virus, proteine tumorali, epitopi).
#
# DIPENDE DA: 01_build_clonotypes.R (full_data e top_clones_CDR3 in memoria)
#
# CORREZIONI rispetto a 2_find_cd3_in_database.R:
#   - Colonne aggiornate: TRB_cdr3 / TRA_cdr3 (minuscolo, nuovo script)
#   - Aggiunto output xlsx riassuntivo
#   - Aggiunto fallback: legge da xlsx se full_data non è in memoria
# ==============================================================================

# ── Configurazione ─────────────────────────────────────────────────────────────
output_dir <- "/Users/federicadannunzio/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/caruana-project/CART/Code/Seurat_analysis/4_clonotypes_expansion_analysis/res/res_fixed"

xlsx_fallback <- file.path(output_dir, "RISULTATI_Top10_Cloni_CDR3.xlsx")

# ── Controllo: usa full_data in memoria oppure carica da xlsx ──────────────────
if (!exists("top_clones_CDR3")) {
  message("top_clones_CDR3 non trovato in memoria — carico da xlsx")
  if (!file.exists(xlsx_fallback)) {
    stop("File xlsx non trovato: ", xlsx_fallback,
         "\nEsegui prima 01_build_clonotypes.R")
  }
  top_clones_CDR3 <- read_xlsx(xlsx_fallback)
}

# Verifica colonne (gestisce sia vecchio che nuovo formato)
has_new_cols <- all(c("TRA_cdr3","TRB_cdr3") %in% colnames(top_clones_CDR3))
has_old_cols <- all(c("TRA_CDR3","TRB_CDR3") %in% colnames(top_clones_CDR3))

if (has_new_cols) {
  col_alpha <- "TRA_cdr3"
  col_beta  <- "TRB_cdr3"
  message("✓ Usando colonne nuovo formato (TRA_cdr3, TRB_cdr3)")
} else if (has_old_cols) {
  col_alpha <- "TRA_CDR3"
  col_beta  <- "TRB_CDR3"
  message("⚠ Usando colonne vecchio formato (TRA_CDR3, TRB_CDR3)")
  message("  Considera di rieseguire 01_build_clonotypes.R per il formato aggiornato")
} else {
  stop("Colonne CDR3 non trovate. Colonne disponibili: ",
       paste(colnames(top_clones_CDR3), collapse=", "))
}

message("Cloni da cercare: ", nrow(top_clones_CDR3))

# ── Download VDJdb ─────────────────────────────────────────────────────────────
message("\n--- Scaricamento VDJdb (2024-06-13) ---")
url_zip  <- "https://github.com/antigenomics/vdjdb-db/releases/download/2024-06-13/vdjdb-2024-06-13.zip"
temp_zip <- tempfile(fileext=".zip")
temp_dir <- tempdir()

tryCatch({
  download.file(url_zip, temp_zip, mode="wb", quiet=FALSE)
}, error=function(e) {
  stop("Download fallito. Controlla la connessione internet.\n", e$message)
})

message("--- Estrazione ---")
unzip(temp_zip, exdir=temp_dir)
target_file <- list.files(temp_dir,
                           pattern="vdjdb.*\\.txt|vdjdb.*\\.tsv",
                           full.names=TRUE, recursive=TRUE)[1]

if (is.na(target_file)) stop("File VDJdb non trovato nello zip estratto")
message("File VDJdb: ", target_file)

vdjdb_raw <- read_tsv(target_file, show_col_types=FALSE)
message("Righe VDJdb totali: ", nrow(vdjdb_raw))
message("Colonne disponibili: ", paste(colnames(vdjdb_raw), collapse=", "))

# ── Preparazione database ──────────────────────────────────────────────────────
# VDJdb ha una riga per coppia TCR — contiene cdr3.alpha e cdr3.beta
message("\n--- Preparazione database ---")

# DB per beta
vdjdb_beta <- vdjdb_raw %>%
  filter(species == "HomoSapiens") %>%
  select(cdr3.beta, antigen.species, antigen.gene, antigen.epitope,
         any_of(c("v.beta","j.beta","mhc.a","mhc.b","score"))) %>%
  rename(CDR3_beta = cdr3.beta,
         Antigen_species = antigen.species,
         Antigen_gene    = antigen.gene,
         Epitope         = antigen.epitope) %>%
  filter(!is.na(CDR3_beta), CDR3_beta != "") %>%
  distinct()

# DB per alpha
vdjdb_alpha <- vdjdb_raw %>%
  filter(species == "HomoSapiens") %>%
  select(cdr3.alpha, antigen.species, antigen.gene, antigen.epitope,
         any_of(c("v.alpha","j.alpha","mhc.a","mhc.b","score"))) %>%
  rename(CDR3_alpha      = cdr3.alpha,
         Antigen_species = antigen.species,
         Antigen_gene    = antigen.gene,
         Epitope         = antigen.epitope) %>%
  filter(!is.na(CDR3_alpha), CDR3_alpha != "") %>%
  distinct()

message("Sequenze beta nel DB:  ", n_distinct(vdjdb_beta$CDR3_beta))
message("Sequenze alpha nel DB: ", n_distinct(vdjdb_alpha$CDR3_alpha))

# ── Ricerca match ──────────────────────────────────────────────────────────────
message("\n--- Ricerca match ---")

my_beta  <- top_clones_CDR3 %>% pull(!!sym(col_beta))  %>% unique() %>% na.omit()
my_alpha <- top_clones_CDR3 %>% pull(!!sym(col_alpha)) %>% unique() %>% na.omit()

matches_beta <- top_clones_CDR3 %>%
  rename(CDR3_beta = !!sym(col_beta)) %>%
  inner_join(vdjdb_beta, by="CDR3_beta")

matches_alpha <- top_clones_CDR3 %>%
  rename(CDR3_alpha = !!sym(col_alpha)) %>%
  inner_join(vdjdb_alpha, by="CDR3_alpha")

# ── Risultati ──────────────────────────────────────────────────────────────────
message("\n==================================================")
message("               RISULTATI RICERCA                 ")
message("==================================================")
message("CDR3 beta cercate:  ", length(my_beta))
message("CDR3 alpha cercate: ", length(my_alpha))

found <- FALSE

if (nrow(matches_beta) > 0) {
  found <- TRUE
  message("\n✅ MATCH SU CATENA BETA (", nrow(matches_beta), " righe):")
  print(matches_beta %>%
          select(patient, any_of(c("Gene_Label","total_n")),
                 CDR3_beta, Antigen_species, Antigen_gene, Epitope))
} else {
  message("\n❌ Nessun match esatto sulla beta")
}

if (nrow(matches_alpha) > 0) {
  found <- TRUE
  message("\n✅ MATCH SU CATENA ALPHA (", nrow(matches_alpha), " righe):")
  print(matches_alpha %>%
          select(patient, any_of(c("Gene_Label","total_n")),
                 CDR3_alpha, Antigen_species, Antigen_gene, Epitope))
} else {
  message("\n❌ Nessun match esatto sulla alpha")
}

if (!found) {
  message("\n--------------------------------------------------")
  message("CONCLUSIONE: I tuoi cloni sono 'Privati'.")
  message("Non sono nel database VDJdb (giugno 2024).")
  message("Possibili spiegazioni:")
  message("  1. Cloni specifici per l'HLA dei pazienti (molto comune)")
  message("  2. Cloni specifici per antigeni tumorali non ancora in DB")
  message("  3. Cloni privati CAR-T non presenti in database pubblici")
}

# ── Salvataggio ────────────────────────────────────────────────────────────────
write_xlsx(list(
  "Match_Beta"  = if(nrow(matches_beta)>0)  matches_beta  else data.frame(nota="nessun match"),
  "Match_Alpha" = if(nrow(matches_alpha)>0) matches_alpha else data.frame(nota="nessun match"),
  "CDR3_cercate" = data.frame(
    catena = c(rep("beta",  length(my_beta)),
               rep("alpha", length(my_alpha))),
    cdr3   = c(my_beta, my_alpha)
  )
), file.path(output_dir, "VDJdb_Summary.xlsx"))

message("\nSalvato: VDJdb_Summary.xlsx in ", output_dir)

# Pulizia file temporanei
unlink(temp_zip)
