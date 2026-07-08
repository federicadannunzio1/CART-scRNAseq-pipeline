# Conserved Clonotype Families
# Converted from 05_conserved_families.Rmd

# ======================================================================
# LIBRARIES
# ======================================================================
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(writexl)
library(readr)
library(purrr)
library(DT)
library(ggrepel)
library(patchwork)

# ======================================================================
# PATHS
# ======================================================================
BASE    <- "/Users/federicadannunzio/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/caruana-project/CART/Code/07_clonotypes"
OUT_DIR <- file.path(BASE, "results")

dir.create(file.path(OUT_DIR, "tables"),  recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(OUT_DIR, "figures"), recursive=TRUE, showWarnings=FALSE)

# ======================================================================
# BUILD CLONOTYPES
# ======================================================================
source(file.path(BASE, "01_build_clonotypes.R"))
cat("\nfull_data caricato:", nrow(full_data), "cellule CAR+\n")
cat("Distribuzione paziente × stage:\n")
print(table(full_data$patient, full_data$stage))

# ======================================================================
# CONTAMINATION FILTER
# ======================================================================
# ── Passo 1: conta cellule per clone × campione (patient × stage) ─────────────
clone_counts <- full_data %>%
  filter(Clone_Quality == "Complete",
         !is.na(TRA_cdr3_nt), !is.na(TRB_cdr3_nt),
         TRA_cdr3_nt != "", TRB_cdr3_nt != "") %>%
  mutate(sample_id = paste(patient, stage, sep="_")) %>%
  group_by(TRA_cdr3_nt, TRB_cdr3_nt, TRA_cdr3, TRB_cdr3, sample_id, patient, stage) %>%
  summarise(n_cells = n(), .groups = "drop")

# ── Passo 2: identifica cloni condivisi tra pazienti (non tra stage dello stesso pz) ─
shared_nt <- clone_counts %>%
  group_by(TRA_cdr3_nt, TRB_cdr3_nt) %>%
  filter(n_distinct(patient) > 1) %>%
  ungroup()

cat("Cloni (CDR3_nt) condivisi tra ≥2 pazienti:", n_distinct(paste(shared_nt$TRA_cdr3_nt, shared_nt$TRB_cdr3_nt)), "\n")
cat("Cellule coinvolte:", nrow(full_data %>% filter(Clone_Quality=="Complete") %>%
  semi_join(shared_nt %>% distinct(TRA_cdr3_nt, TRB_cdr3_nt),
            by=c("TRA_cdr3_nt","TRB_cdr3_nt"))), "\n\n")

# ── Passo 3: per ogni clone condiviso, trova il campione con più cellule ────────
# Aggregato per paziente (somma su tutti i stage di quel paziente)
dominant_patient <- shared_nt %>%
  group_by(TRA_cdr3_nt, TRB_cdr3_nt, patient) %>%
  summarise(n_cells_tot = sum(n_cells), .groups = "drop") %>%
  group_by(TRA_cdr3_nt, TRB_cdr3_nt) %>%
  slice_max(n_cells_tot, n=1, with_ties=FALSE) %>%   # un solo dominante
  ungroup() %>%
  rename(dominant_patient = patient, dominant_n = n_cells_tot)

# Tabella di report
contamination_report <- shared_nt %>%
  group_by(TRA_cdr3_nt, TRB_cdr3_nt, TRA_cdr3, TRB_cdr3) %>%
  summarise(
    pazienti = paste(sort(unique(patient)), collapse=" & "),
    dettaglio = paste(paste(patient, stage, n_cells, sep="×"), collapse=" | "),
    .groups = "drop"
  ) %>%
  left_join(dominant_patient, by=c("TRA_cdr3_nt","TRB_cdr3_nt")) %>%
  arrange(desc(dominant_n))

cat("Report cloni condivisi:\n")
print(contamination_report %>%
        select(TRA_cdr3, TRB_cdr3, pazienti, dominant_patient, dominant_n, dettaglio))

write_xlsx(contamination_report,
           file.path(OUT_DIR, "tables", "contamination_report.xlsx"))

# ======================================================================
# APPLY FILTER
# ======================================================================
# ── Passo 4: escludi i cloni condivisi da TUTTI i pazienti ──────────────────────
# Motivazione: CDR3_nt identica tra pazienti con terapia autologa non condivisa
# indica cross-contaminazione durante la manifattura → il clone non è attribuibile
# con certezza a nessun paziente, quindi viene rimosso da tutti.
exclude_clones <- shared_nt %>%
  distinct(TRA_cdr3_nt, TRB_cdr3_nt)   # tutte le coppie nt condivise

n_before <- nrow(full_data %>% filter(Clone_Quality=="Complete"))

clean_data <- full_data %>%
  filter(Clone_Quality == "Complete",
         !is.na(TRA_cdr3_nt), !is.na(TRB_cdr3_nt)) %>%
  anti_join(exclude_clones, by = c("TRA_cdr3_nt","TRB_cdr3_nt"))

n_after  <- nrow(clean_data)
n_removed <- n_before - n_after

cat(sprintf("Cellule prima del filtro:   %d\n", n_before))
cat(sprintf("Cellule rimosse (contamin.): %d (%.1f%%)\n",
            n_removed, 100*n_removed/n_before))
cat(sprintf("Cellule mantenute:           %d\n\n", n_after))
cat("Distribuzione post-filtro (paziente × stage):\n")
print(table(clean_data$patient, clean_data$stage))

# ======================================================================
# STRICT SHARING
# ======================================================================
strict_shared <- clean_data %>%
  group_by(TRA_cdr3, TRB_cdr3) %>%
  filter(n_distinct(patient) > 1) %>%
  summarise(
    pazienti    = paste(sort(unique(patient)), collapse=" & "),
    n_pazienti  = n_distinct(patient),
    n_cellule   = n(),
    TRA_v_gene  = first(TRA_v_gene),
    TRB_v_gene  = first(TRB_v_gene),
    stages      = paste(sort(unique(paste(patient,stage,sep="-"))), collapse="; "),
    TRA_cdr3_nt = first(TRA_cdr3_nt),
    TRB_cdr3_nt = first(TRB_cdr3_nt),
    .groups = "drop"
  ) %>%
  arrange(desc(n_pazienti), desc(n_cellule))

cat("Clonotipi (TRA+TRB CDR3 aa) condivisi tra pazienti dopo decontaminazione:\n")
cat("N =", nrow(strict_shared), "\n\n")

if (nrow(strict_shared) > 0) {
  print(strict_shared %>% select(TRA_v_gene, TRB_v_gene, TRA_cdr3, TRB_cdr3,
                                  pazienti, n_cellule, stages))
  write_xlsx(strict_shared,
             file.path(OUT_DIR, "tables", "conserved_strict_TRA_TRB.xlsx"))
} else {
  cat("Nessun clone con CDR3 aa identica condiviso. Vedi sharing parziale.\n")
}

# ======================================================================
# BETA SHARING
# ======================================================================
beta_shared <- clean_data %>%
  group_by(TRB_cdr3) %>%
  filter(n_distinct(patient) > 1) %>%
  summarise(
    pazienti      = paste(sort(unique(patient)), collapse=" & "),
    n_pazienti    = n_distinct(patient),
    n_cellule     = n(),
    TRB_v_gene    = first(TRB_v_gene),
    TRB_j_gene    = first(TRB_j_gene),
    TRA_cdr3_list = paste(unique(TRA_cdr3), collapse=" / "),
    TRA_v_list    = paste(unique(TRA_v_gene), collapse=" / "),
    TRB_cdr3_nt   = first(TRB_cdr3_nt),
    .groups = "drop"
  ) %>%
  arrange(desc(n_pazienti), desc(n_cellule))

cat("Clonotipi con stessa CDR3 beta condivisi tra pazienti:\n")
cat("N =", nrow(beta_shared), "\n\n")

if (nrow(beta_shared) > 0) {
  print(beta_shared %>% select(TRB_v_gene, TRB_j_gene, TRB_cdr3,
                                n_pazienti, n_cellule, pazienti,
                                TRA_cdr3_list))
  write_xlsx(beta_shared,
             file.path(OUT_DIR, "tables", "conserved_beta_CDR3.xlsx"))
}

# ======================================================================
# ALPHA SHARING
# ======================================================================
alpha_shared <- clean_data %>%
  group_by(TRA_cdr3) %>%
  filter(n_distinct(patient) > 1) %>%
  summarise(
    pazienti      = paste(sort(unique(patient)), collapse=" & "),
    n_pazienti    = n_distinct(patient),
    n_cellule     = n(),
    TRA_v_gene    = first(TRA_v_gene),
    TRA_j_gene    = first(TRA_j_gene),
    TRB_cdr3_list = paste(unique(TRB_cdr3), collapse=" / "),
    TRB_v_list    = paste(unique(TRB_v_gene), collapse=" / "),
    TRA_cdr3_nt   = first(TRA_cdr3_nt),
    .groups = "drop"
  ) %>%
  arrange(desc(n_pazienti), desc(n_cellule))

cat("Clonotipi con stessa CDR3 alpha condivisi tra pazienti:\n")
cat("N =", nrow(alpha_shared), "\n\n")
if (nrow(alpha_shared) > 0) {
  print(alpha_shared %>% select(TRA_v_gene, TRA_j_gene, TRA_cdr3,
                                 n_pazienti, n_cellule, pazienti))
  write_xlsx(alpha_shared,
             file.path(OUT_DIR, "tables", "conserved_alpha_CDR3.xlsx"))
}

# ======================================================================
# MASTER TABLE
# ======================================================================
# Candidati: clonotipi condivisi per beta (più inclusivo)
# Se nessuno per beta, usa tutti i privati top espansi
candidates <- if (nrow(beta_shared) > 0) {
  clean_data %>%
    semi_join(beta_shared, by="TRB_cdr3")
} else {
  clean_data %>%
    group_by(patient, TRA_cdr3, TRB_cdr3) %>%
    mutate(n_cells=n()) %>% ungroup() %>%
    filter(n_cells >= 5)
}

master_table <- candidates %>%
  group_by(patient, stage, TRA_v_gene, TRA_j_gene,
           TRB_v_gene, TRB_d_gene, TRB_j_gene,
           TRA_cdr3, TRB_cdr3, TRA_cdr3_nt, TRB_cdr3_nt) %>%
  summarise(
    n_cells    = n(),
    shared_TRB = TRB_cdr3[1] %in% beta_shared$TRB_cdr3,
    shared_TRA = TRA_cdr3[1] %in% alpha_shared$TRA_cdr3,
    .groups    = "drop"
  ) %>%
  arrange(desc(shared_TRB & shared_TRA), desc(n_cells))

write_xlsx(master_table,
           file.path(OUT_DIR, "tables", "master_conserved_sequences.xlsx"))

datatable(
  master_table %>% select(patient, stage, TRA_v_gene, TRB_v_gene,
                           TRA_cdr3, TRB_cdr3, n_cells,
                           shared_TRB, shared_TRA),
  caption = "Sequenze dei clonotipi conservati (beta condivisa tra pazienti)",
  filter  = "top",
  options = list(pageLength=20, scrollX=TRUE),
  rownames= FALSE
)

# ======================================================================
# VDJDB DOWNLOAD
# ======================================================================
# Scarica VDJdb 2024-06-13
url_zip  <- "https://github.com/antigenomics/vdjdb-db/releases/download/2024-06-13/vdjdb-2024-06-13.zip"
temp_zip <- tempfile(fileext=".zip")
temp_dir <- tempdir()

dl_ok <- tryCatch({
  download.file(url_zip, temp_zip, mode="wb", quiet=TRUE)
  # Verify the file is a valid zip (not an error page)
  file.size(temp_zip) > 1e6
}, error=function(e) { message("⚠ Download fallito: ", e$message); FALSE })

if (isTRUE(dl_ok)) {
  unzip(temp_zip, exdir=temp_dir, overwrite=TRUE)
  vdjdb_file <- list.files(temp_dir,
                           pattern="vdjdb.*\\.txt|vdjdb.*\\.tsv",
                           full.names=TRUE, recursive=TRUE)[1]
  vdjdb_raw  <- read_tsv(vdjdb_file, show_col_types=FALSE)
  cat("VDJdb scaricato:", nrow(vdjdb_raw), "righe\n")
  unlink(temp_zip)
  vdjdb_ok <- TRUE
} else {
  message("⚠ VDJdb non disponibile — proseguo senza annotazione VDJdb (flag vdjdb_hit=FALSE)")
  vdjdb_raw <- data.frame(species=character(), cdr3.beta=character(), v.beta=character(),
                          j.beta=character(), cdr3.alpha=character(), v.alpha=character(),
                          j.alpha=character(), antigen.species=character(),
                          antigen.gene=character(), antigen.epitope=character(),
                          stringsAsFactors=FALSE)
  vdjdb_ok <- FALSE
}

# ======================================================================
# VDJDB PREP
# ======================================================================
# Filtra su Homo sapiens e prepara DB beta e alpha
vdjdb_hs <- vdjdb_raw %>% filter(species == "HomoSapiens")

vdjdb_beta <- vdjdb_hs %>%
  filter(!is.na(cdr3.beta), cdr3.beta != "") %>%
  select(cdr3.beta, v.beta, j.beta,
         antigen.species, antigen.gene, antigen.epitope,
         any_of(c("score","mhc.a","reference.id"))) %>%
  distinct() %>%
  rename(TRB_cdr3 = cdr3.beta)

vdjdb_alpha <- vdjdb_hs %>%
  filter(!is.na(cdr3.alpha), cdr3.alpha != "") %>%
  select(cdr3.alpha, v.alpha, j.alpha,
         antigen.species, antigen.gene, antigen.epitope,
         any_of(c("score","mhc.a","reference.id"))) %>%
  distinct() %>%
  rename(TRA_cdr3 = cdr3.alpha)

cat("Sequenze beta nel DB:  ", n_distinct(vdjdb_beta$TRB_cdr3), "\n")
cat("Sequenze alpha nel DB: ", n_distinct(vdjdb_alpha$TRA_cdr3), "\n")

# ======================================================================
# VDJDB SEARCH
# ======================================================================
# Cerca tutte le CDR3 dei cloni puliti (non solo top 10)
all_beta  <- unique(na.omit(clean_data$TRB_cdr3))
all_alpha <- unique(na.omit(clean_data$TRA_cdr3))

cat("CDR3 beta da cercare:  ", length(all_beta), "\n")
cat("CDR3 alpha da cercare: ", length(all_alpha), "\n\n")

# Join diretto
matches_beta <- clean_data %>%
  distinct(patient, TRA_v_gene, TRB_v_gene,
           TRA_cdr3, TRB_cdr3, TRA_cdr3_nt, TRB_cdr3_nt) %>%
  inner_join(vdjdb_beta, by="TRB_cdr3") %>%
  select(patient, TRB_v_gene, v.beta, TRA_v_gene,
         TRB_cdr3, TRA_cdr3, TRB_cdr3_nt, TRA_cdr3_nt,
         antigen.species, antigen.gene, antigen.epitope,
         any_of(c("score","mhc.a")))

matches_alpha <- clean_data %>%
  distinct(patient, TRA_v_gene, TRB_v_gene,
           TRA_cdr3, TRB_cdr3, TRA_cdr3_nt, TRB_cdr3_nt) %>%
  inner_join(vdjdb_alpha, by="TRA_cdr3") %>%
  select(patient, TRA_v_gene, v.alpha, TRB_v_gene,
         TRA_cdr3, TRB_cdr3, TRA_cdr3_nt, TRB_cdr3_nt,
         antigen.species, antigen.gene, antigen.epitope,
         any_of(c("score","mhc.a")))

cat("=== MATCH CATENA BETA ===\n")
if (nrow(matches_beta) > 0) {
  cat("Trovati:", nrow(matches_beta), "match su", n_distinct(matches_beta$TRB_cdr3), "CDR3 distinte\n")
  print(matches_beta)
} else {
  cat("Nessun match esatto su catena beta\n")
}

cat("\n=== MATCH CATENA ALPHA ===\n")
if (nrow(matches_alpha) > 0) {
  cat("Trovati:", nrow(matches_alpha), "match su", n_distinct(matches_alpha$TRA_cdr3), "CDR3 distinte\n")
  print(matches_alpha)
} else {
  cat("Nessun match esatto su catena alpha\n")
}

# ======================================================================
# VDJDB LEUKEMIA
# ======================================================================
# Filtra specificamente per antigeni leucemici / tumorali noti
leukemia_keywords <- c("leukemia","AML","CLL","CML","ALL","B-cell","myeloid",
                        "WT1","PRAME","survivin","NY-ESO","MAGE","AFP","RHAMM",
                        "CD19","CD33","CD123","FLT3","NPM1","IDH","BCMA","CD38",
                        "tumor","cancer","lymphoma","acute")

leuk_pattern <- paste(leukemia_keywords, collapse="|")

leukemia_beta <- vdjdb_beta %>%
  filter(str_detect(tolower(antigen.species), leuk_pattern) |
         str_detect(tolower(antigen.gene),    leuk_pattern) |
         str_detect(tolower(antigen.epitope), leuk_pattern))

leukemia_alpha <- vdjdb_alpha %>%
  filter(str_detect(tolower(antigen.species), leuk_pattern) |
         str_detect(tolower(antigen.gene),    leuk_pattern) |
         str_detect(tolower(antigen.epitope), leuk_pattern))

cat("\nCDR3 beta nel DB associate ad antigeni leucemici/tumorali:",
    n_distinct(leukemia_beta$TRB_cdr3), "\n")
cat("CDR3 alpha nel DB associate ad antigeni leucemici/tumorali:",
    n_distinct(leukemia_alpha$TRA_cdr3), "\n\n")

# Cerca se i nostri cloni matchano sequenze leucemia-associate
leuk_matches_beta <- clean_data %>%
  distinct(patient, TRA_v_gene, TRB_v_gene, TRA_cdr3, TRB_cdr3,
           TRA_cdr3_nt, TRB_cdr3_nt) %>%
  inner_join(leukemia_beta, by="TRB_cdr3")

leuk_matches_alpha <- clean_data %>%
  distinct(patient, TRA_v_gene, TRB_v_gene, TRA_cdr3, TRB_cdr3,
           TRA_cdr3_nt, TRB_cdr3_nt) %>%
  inner_join(leukemia_alpha, by="TRA_cdr3")

if (nrow(leuk_matches_beta) > 0) {
  cat("✅ MATCH LEUCEMIA — CATENA BETA:\n")
  print(leuk_matches_beta %>% select(patient, TRB_v_gene, TRB_cdr3, TRA_cdr3,
                                      antigen.species, antigen.gene, antigen.epitope))
} else {
  cat("Nessun match esatto su beta contro antigeni leucemici in VDJdb\n")
}

if (nrow(leuk_matches_alpha) > 0) {
  cat("\n✅ MATCH LEUCEMIA — CATENA ALPHA:\n")
  print(leuk_matches_alpha %>% select(patient, TRA_v_gene, TRA_cdr3, TRB_cdr3,
                                       antigen.species, antigen.gene, antigen.epitope))
} else {
  cat("Nessun match esatto su alpha contro antigeni leucemici in VDJdb\n")
}

# ======================================================================
# VDJDB SAVE
# ======================================================================
write_xlsx(list(
  "Matches_beta"        = if(nrow(matches_beta)>0)       matches_beta       else data.frame(nota="nessun match"),
  "Matches_alpha"       = if(nrow(matches_alpha)>0)      matches_alpha      else data.frame(nota="nessun match"),
  "Leukemia_beta"       = if(nrow(leuk_matches_beta)>0)  leuk_matches_beta  else data.frame(nota="nessun match"),
  "Leukemia_alpha"      = if(nrow(leuk_matches_alpha)>0) leuk_matches_alpha else data.frame(nota="nessun match"),
  "VDJdb_leuk_beta_ref" = leukemia_beta,
  "VDJdb_leuk_alpha_ref"= leukemia_alpha
), file.path(OUT_DIR, "tables", "VDJdb_search_results.xlsx"))

# ======================================================================
# FINAL SEQUENCES
# ======================================================================
# Top cloni per paziente (post-decontaminazione)
top_clones_clean <- clean_data %>%
  group_by(patient, stage, TRA_v_gene, TRA_j_gene,
           TRB_v_gene, TRB_d_gene, TRB_j_gene,
           TRA_cdr3, TRB_cdr3, TRA_cdr3_nt, TRB_cdr3_nt) %>%
  summarise(n_cells = n(), .groups="drop") %>%
  mutate(
    shared_beta  = TRB_cdr3 %in% beta_shared$TRB_cdr3,
    shared_alpha = TRA_cdr3 %in% alpha_shared$TRA_cdr3,
    shared_both  = shared_beta & shared_alpha,
    vdjdb_hit    = TRB_cdr3 %in% matches_beta$TRB_cdr3 |
                   TRA_cdr3 %in% matches_alpha$TRA_cdr3,
    leuk_hit     = TRB_cdr3 %in% leuk_matches_beta$TRB_cdr3 |
                   TRA_cdr3 %in% leuk_matches_alpha$TRA_cdr3
  ) %>%
  arrange(desc(leuk_hit), desc(vdjdb_hit),
          desc(shared_both), desc(shared_beta), desc(n_cells))

write_xlsx(top_clones_clean,
           file.path(OUT_DIR, "tables", "final_clone_sequences.xlsx"))

datatable(
  top_clones_clean,
  caption  = "Tutti i clonotipi CAR+ post-decontaminazione con annotazioni",
  filter   = "top",
  options  = list(pageLength=25, scrollX=TRUE),
  rownames = FALSE
) %>%
  DT::formatStyle("leuk_hit",
                  backgroundColor=DT::styleEqual(TRUE, "#FFE4E1"),
                  fontWeight=DT::styleEqual(TRUE, "bold")) %>%
  DT::formatStyle("shared_beta",
                  backgroundColor=DT::styleEqual(TRUE, "#E8F5E9"))

# ======================================================================
# SUMMARY BLOCK
# ======================================================================
cat("=== RIEPILOGO P5 ===\n\n")
cat("Cellule CAR+ totali (complete TRA+TRB):  ", nrow(full_data %>% filter(Clone_Quality=="Complete")), "\n")
cat("Cellule post-decontaminazione:           ", nrow(clean_data), "\n")
cat("Clonotipi unici post-filtro:             ", n_distinct(paste(clean_data$TRA_cdr3, clean_data$TRB_cdr3)), "\n\n")

cat("Famiglie conservate (CDR3 beta condivisa tra pazienti):  ", nrow(beta_shared), "\n")
cat("Famiglie conservate (CDR3 alpha condivisa tra pazienti): ", nrow(alpha_shared), "\n")
cat("Clonotipi con TRA+TRB identici tra pazienti:            ", nrow(strict_shared), "\n\n")

cat("Match VDJdb (tutte specificità):    beta=", nrow(matches_beta),
    "  alpha=", nrow(matches_alpha), "\n")
cat("Match VDJdb (leucemia/tumore):      beta=", nrow(leuk_matches_beta),
    "  alpha=", nrow(leuk_matches_alpha), "\n\n")

cat("File prodotti in:", OUT_DIR, "\n")
cat("  tables/contamination_report.xlsx\n")
cat("  tables/conserved_strict_TRA_TRB.xlsx\n")
cat("  tables/conserved_beta_CDR3.xlsx\n")
cat("  tables/conserved_alpha_CDR3.xlsx\n")
cat("  tables/master_conserved_sequences.xlsx\n")
cat("  tables/VDJdb_search_results.xlsx  (match + ref leucemia)\n")
cat("  tables/final_clone_sequences.xlsx  (tabella completa con flag)\n")
