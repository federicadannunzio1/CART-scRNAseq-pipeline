library(purrr); library(dplyr); library(tidyr); library(stringr)
library(readr); library(ggplot2); library(writexl)

# ==============================================================================
# SCRIPT: 2_verify_shared_clonotypes_CORRETTO.R
#
# CORREZIONI RISPETTO ALLO SCRIPT ORIGINALE:
#
#  BUG #1 — paste(..., collapse="/") sulle catene
#    Concatenava CDR3 e V gene di catene doppie della stessa cellula.
#    Questo produceva TRBV7-9/TRBV6-4 e Clone_ID ibridi.
#    FIX: si seleziona UNA catena per barcode (max UMI → max reads).
#
#  BUG #2 — Filtri incompleti
#    Filtrava solo per productive. Mancavano is_cell, high_confidence,
#    full_length — ogni filtro mancante introduce contigs spurii.
#    FIX: tutti i filtri con tracciamento righe rimosse.
#
#  NUOVO — CDR1 e CDR2 aa e nt
#    Auto-rilevate dal CSV e confrontate tra pazienti.
#
#  NUOVO — CDR3 nucleotidica
#    Test principale: nt identica → contaminazione; nt diversa → convergenza.
#
#  NUOVO — Barcode check
#    Stesso barcode in pazienti diversi = contaminazione tecnica certa.
# ==============================================================================


# ==============================================================================
# 1. PERCORSI (identici al tuo script originale)
# ==============================================================================
seurat_path <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/CART/Data/seurat_obj_list/seurat_samples_sctype_azimuth_pbmc_bonemarrow_clonalvdj_CAR.rds"

base_path <- "/Users/federicadannunzio/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/caruana-project/CART/Data/output_allineamento_original_no_CAR_si_VDJ"

output_dir  <- "/Users/federicadannunzio/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/caruana-project/CART/Code/07_clonotypes/results"
figures_dir <- file.path(output_dir, "figures")
tables_dir  <- file.path(output_dir, "tables")

dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir,  recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 2. STEP 1: CENSIMENTO CELLULE CAR+ (invariato)
# ==============================================================================
message("\n--- STEP 1: CENSIMENTO CELLULE CAR+ ---")
if (!file.exists(seurat_path)) stop("ERRORE: File Seurat non trovato!")
seurat_list <- readRDS(seurat_path)
flat_list   <- unlist(seurat_list, recursive = FALSE)

car_cells_map <- map_df(names(flat_list), function(obj_name) {
  meta  <- flat_list[[obj_name]]@meta.data
  is_car <- rep(FALSE, nrow(meta))
  if ("IS_CAR_ALLIN_scREP" %in% colnames(meta)) {
    is_car <- grepl("YES|TRUE", as.character(meta$IS_CAR_ALLIN_scREP), ignore.case=TRUE)
  } else if ("CAR" %in% colnames(meta)) {
    is_car <- grepl("YES|TRUE", as.character(meta$CAR), ignore.case=TRUE)
  }
  if (sum(is_car) == 0) return(NULL)
  meta %>% filter(is_car) %>%
    mutate(
      obj_name_r    = obj_name,
      full_barcode  = rownames(.),
      folder_name   = str_extract(full_barcode, "^.*(?=_[ACGT]+-[0-9]+)"),
      clean_barcode = str_extract(full_barcode, "[ACGT]+-[0-9]+")
    ) %>%
    select(obj_name_r, folder_name, clean_barcode)
})
message("Cellule CAR+ trovate per cartella:")
print(car_cells_map %>% count(folder_name))


# ==============================================================================
# 3. STEP 2: LETTURA VDJ — VERSIONE CORRETTA
#
# DIFFERENZA CHIAVE dal tuo script:
#   PRIMA: paste(unique(cdr3[chain=="TRA"]), collapse="/")
#     → concatena le CDR3 di eventuali 2 TRA/TRB con "/"
#     → produce TRBV7-9/TRBV6-4 e Clone_ID ibridi
#
#   ORA: per ogni barcode+chain prende UN solo contig (max UMI)
#     → barcode = 1 TRA + 1 TRB, puliti
# ==============================================================================
message("\n--- STEP 2: LETTURA VDJ CORRETTA ---")

read_vdj_corretto <- function(folder, base_dir) {
  f1     <- file.path(base_dir, "1", folder, "vdj_t", "filtered_contig_annotations.csv")
  f2     <- file.path(base_dir, "2", folder, "vdj_t", "filtered_contig_annotations.csv")
  # Fallback: S429_I ha il file con un trattino nel nome (filtered_contig_annotations-.csv)
  f1alt  <- file.path(base_dir, "1", folder, "vdj_t", "filtered_contig_annotations-.csv")
  f2alt  <- file.path(base_dir, "2", folder, "vdj_t", "filtered_contig_annotations-.csv")
  target <- if (file.exists(f1)) f1 else if (file.exists(f2)) f2 else
            if (file.exists(f1alt)) f1alt else if (file.exists(f2alt)) f2alt else NULL
  if (is.null(target)) { message("  MANCANTE: ", folder); return(NULL) }
  if (grepl("-\\.csv$", target)) message("  NOTA: usando file con nome alternativo: ", basename(target))
  message("  TROVATO: ", folder)
  
  tryCatch({
    l1  <- readLines(target, n=1)
    sep <- if (grepl(";", l1)) ";" else ","
    raw <- read_delim(target, delim=sep, show_col_types=FALSE)
    if (ncol(raw) <= 1) raw <- read_csv(target, show_col_types=FALSE)
    
    # Stampa colonne CDR disponibili (utile per debug)
    cdr_cols <- grep("^cdr", colnames(raw), value=TRUE, ignore.case=TRUE)
    message("    CDR disponibili: ", paste(cdr_cols, collapse=", "))
    
    # ── Filtri con tracciamento (BUG #2 FIX) ─────────────────────────────
    n0 <- nrow(raw)
    message("    Righe raw: ", n0)
    
    # F1: is_cell — MANCANTE nel vecchio script
    df <- raw %>% filter(is_cell %in% c(TRUE,"True","true","TRUE"))
    message("    Dopo is_cell:         ", nrow(df), " (rimossi:", n0-nrow(df), ")")
    
    # F2: high_confidence — MANCANTE nel vecchio script
    prev <- nrow(df)
    df <- df %>% filter(high_confidence %in% c(TRUE,"True","true","TRUE"))
    message("    Dopo high_confidence: ", nrow(df), " (rimossi:", prev-nrow(df), ")")
    
    # F3: productive (era presente nel vecchio script)
    prev <- nrow(df)
    df <- df %>% filter(tolower(as.character(productive)) == "true")
    message("    Dopo productive:      ", nrow(df), " (rimossi:", prev-nrow(df), ")")
    
    # F4: full_length — MANCANTE nel vecchio script
    prev <- nrow(df)
    df <- df %>% filter(full_length %in% c(TRUE,"True","true","TRUE"))
    message("    Dopo full_length:     ", nrow(df), " (rimossi:", prev-nrow(df), ")")
    
    # F5: solo TRA e TRB
    prev <- nrow(df)
    df <- df %>% filter(chain %in% c("TRA","TRB"))
    message("    Dopo chain TRA/TRB:   ", nrow(df), " (rimossi:", prev-nrow(df), ")")
    
    # F6: CDR3 valida (aa e nt)
    prev <- nrow(df)
    df <- df %>% filter(
      !is.na(cdr3),    nchar(as.character(cdr3))>0,    cdr3    != "None",
      !is.na(cdr3_nt), nchar(as.character(cdr3_nt))>0, cdr3_nt != "None"
    )
    message("    Dopo CDR3 valida:     ", nrow(df), " (rimossi:", prev-nrow(df), ")")
    
    if (nrow(df) == 0) { message("    NESSUN CONTIG VALIDO — skip"); return(NULL) }
    
    # Avvisi catene doppie (erano la fonte del bug nel vecchio script)
    n_dbl_TRA <- df %>% filter(chain=="TRA") %>% count(barcode) %>% filter(n>1) %>% nrow()
    n_dbl_TRB <- df %>% filter(chain=="TRB") %>% count(barcode) %>% filter(n>1) %>% nrow()
    message("    Cellule con >1 TRA: ", n_dbl_TRA,
            " (nel vecchio script: CDR3 concatenate con '/')")
    message("    Cellule con >1 TRB: ", n_dbl_TRB,
            " (nel vecchio script: produceva TRBV7-9/TRBV6-4)")
    
    # ── BUG #1 FIX: seleziona UNA catena per barcode ─────────────────────
    # Strategia: max UMI → max reads → primo contig (deterministico)
    best_chain <- function(d) {
      d %>% group_by(barcode) %>%
        arrange(desc(umis), desc(reads)) %>%
        slice(1) %>% ungroup()
    }
    
    # Auto-rileva tutte le colonne CDR presenti (cdr1, cdr1_nt, cdr2, cdr2_nt, ecc.)
    all_cdr_cols <- grep("^cdr", colnames(df), value=TRUE, ignore.case=TRUE)
    
    # TRA: miglior contig per cellula
    tra_sel <- intersect(
      c("barcode","v_gene","j_gene", all_cdr_cols, "umis","reads"),
      colnames(df)
    )
    tra <- df %>% filter(chain=="TRA") %>% best_chain() %>%
      select(all_of(tra_sel)) %>%
      rename_with(~paste0("TRA_",.), -barcode)
    
    # TRB: miglior contig per cellula
    trb_sel <- intersect(
      c("barcode","v_gene","d_gene","j_gene", all_cdr_cols, "umis","reads"),
      colnames(df)
    )
    trb <- df %>% filter(chain=="TRB") %>% best_chain() %>%
      select(all_of(trb_sel)) %>%
      rename_with(~paste0("TRB_",.), -barcode)
    
    # Pairing: solo cellule con entrambe le catene
    paired <- inner_join(tra, trb, by="barcode") %>%
      mutate(
        folder_name   = folder,
        # Gene_Label PULITO: non conterrà più "/" da doppia catena
        Gene_Label    = paste(
          if_else(is.na(TRA_v_gene)|TRA_v_gene=="","?",TRA_v_gene), "+",
          if_else(is.na(TRB_v_gene)|TRB_v_gene=="","?",TRB_v_gene)
        ),
        Clone_ID_CDR3 = paste(TRA_cdr3, TRB_cdr3, sep="_"),
        Clone_ID_Full = paste0(TRA_v_gene,":",TRA_j_gene,":",TRA_cdr3,"_",
                               TRB_v_gene,":",TRB_j_gene,":",TRB_cdr3),
        Has_Complete_TRA = (!is.na(TRA_v_gene) & TRA_v_gene!="" &
                              !is.na(TRA_cdr3)   & TRA_cdr3  !=""),
        Has_Complete_TRB = (!is.na(TRB_v_gene) & TRB_v_gene!="" &
                              !is.na(TRB_cdr3)   & TRB_cdr3  !=""),
        Clone_Quality = case_when(
          Has_Complete_TRA & Has_Complete_TRB  ~ "Complete",
          Has_Complete_TRB & !Has_Complete_TRA ~ "TRB_only",
          Has_Complete_TRA & !Has_Complete_TRB ~ "TRA_only",
          TRUE                                 ~ "Incomplete"
        )
      )
    
    message("    Cellule paired TRA+TRB: ", nrow(paired))
    paired
    
  }, error=function(e) { message("  ERRORE: ", folder, " — ", e$message); NULL })
}

unique_folders <- unique(na.omit(car_cells_map$folder_name))
vdj_database   <- map_df(unique_folders, ~read_vdj_corretto(.x, base_path))


# ==============================================================================
# 4. STEP 3: UNIONE SEURAT + VDJ (invariato)
# ==============================================================================
message("\n--- STEP 3: UNIONE DATI ---")

full_data <- car_cells_map %>%
  inner_join(vdj_database, by=c("folder_name","clean_barcode"="barcode")) %>%
  mutate(
    patient = case_when(
      grepl("S151|S393", folder_name) ~ "Ca",
      grepl("S429|S435", folder_name) ~ "Bo",
      grepl("S345|S431", folder_name) ~ "Me",
      TRUE ~ "Unknown"
    ),
    stage = case_when(
      grepl("_I", folder_name) ~ "I",
      grepl("_A", folder_name) ~ "A",
      grepl("_B", folder_name) ~ "B",
      TRUE ~ "Unknown"
    )
  )

message("Distribuzione Pazienti x Stage:")
print(table(full_data$patient, full_data$stage))
message("Qualita cloni:")
print(table(full_data$Clone_Quality, full_data$patient))

# ==============================================================================
# 5. STEP 4: TOP 10 CLONI (invariato nella logica)
# ==============================================================================
message("\n--- STEP 4: TOP 10 CLONI ---")

# Filtro contaminanti inter-paziente (CDR3_nt identica tra pazienti diversi)
# Identico al filtro usato in 03 e 04 — rimuove i cloni da tutti i pazienti.
# Motivazione: Bo e Me non hanno ricevuto lo stesso lotto CAR-T, quindi
# cloni con CDR3_nt identica tra pazienti sono contaminazione di laboratorio.
exclude_nt_plot <- full_data %>%
  filter(Clone_Quality=="Complete",
         !is.na(TRA_cdr3_nt), !is.na(TRB_cdr3_nt),
         TRA_cdr3_nt != "", TRB_cdr3_nt != "") %>%
  group_by(TRA_cdr3_nt, TRB_cdr3_nt, patient) %>%
  summarise(n=n(), .groups="drop") %>%
  group_by(TRA_cdr3_nt, TRB_cdr3_nt) %>%
  filter(n_distinct(patient) > 1) %>%
  ungroup() %>%
  distinct(TRA_cdr3_nt, TRB_cdr3_nt)

clean_data_plot <- full_data %>%
  filter(Clone_Quality=="Complete",
         !is.na(TRA_cdr3_nt), !is.na(TRB_cdr3_nt)) %>%
  anti_join(exclude_nt_plot, by=c("TRA_cdr3_nt","TRB_cdr3_nt"))

# Top cloni per stage I e stage B separatamente, poi unione.
# Motivo: i cloni dominanti in B (espansi) hanno centinaia di cellule e
# mascherano completamente i cloni presenti nel prodotto di infusione (I),
# che hanno poche cellule. L'unione permette di vedere il turnover clonale.
clone_counts_per_stage <- clean_data_plot %>%
  group_by(patient, stage, Clone_ID_CDR3, Gene_Label, TRA_cdr3, TRB_cdr3) %>%
  summarise(n_cells=n(), .groups="drop")

top_I_clones <- clone_counts_per_stage %>%
  filter(stage=="I") %>%
  group_by(patient) %>% slice_max(n_cells, n=5, with_ties=FALSE) %>%
  ungroup() %>% select(patient, Clone_ID_CDR3)

top_B_clones <- clone_counts_per_stage %>%
  filter(stage=="B") %>%
  group_by(patient) %>% slice_max(n_cells, n=5, with_ties=FALSE) %>%
  ungroup() %>% select(patient, Clone_ID_CDR3)

top_IB_union <- bind_rows(top_I_clones, top_B_clones) %>%
  distinct(patient, Clone_ID_CDR3)

top_clones_CDR3 <- clean_data_plot %>%
  semi_join(top_IB_union, by=c("patient","Clone_ID_CDR3")) %>%
  group_by(patient, Clone_ID_CDR3, Gene_Label, TRA_cdr3, TRB_cdr3) %>%
  summarise(total_n=n(), .groups="drop") %>%
  arrange(patient, desc(total_n)) %>%
  group_by(patient) %>% mutate(Rank=row_number()) %>% ungroup()

top_clones_Vgene <- full_data %>%
  group_by(patient, Gene_Label) %>%
  summarise(total_n=n(), .groups="drop") %>%
  arrange(patient, desc(total_n)) %>%
  group_by(patient) %>% slice_head(n=10) %>%
  mutate(Rank_ID=paste0("Clone ", row_number()))


# ==============================================================================
# 6. STEP 5: ANALISI CLONOTIPI CONDIVISI — VERSIONE CORRETTA
# ==============================================================================
message("\n--- STEP 5: ANALISI CLONOTIPI CONDIVISI TRA PAZIENTI ---")

tra_cdr_cols <- grep("^TRA_cdr", colnames(full_data), value=TRUE)
trb_cdr_cols <- grep("^TRB_cdr", colnames(full_data), value=TRUE)
gene_cols    <- grep("^TR[AB]_(v|d|j)_gene$", colnames(full_data), value=TRUE)
all_cdr_cols <- c(tra_cdr_cols, trb_cdr_cols)

# Funzione: valore più frequente (consensus clonotipo)
top1 <- function(x) {
  t <- table(x[!is.na(x) & x!=""])
  if (length(t)==0) return(NA_character_)
  names(sort(t, decreasing=TRUE))[1]
}

# Aggrega in clonotipi
clonotipi_summary <- full_data %>%
  filter(Clone_Quality=="Complete") %>%
  group_by(patient, stage, TRA_v_gene, TRA_j_gene,
           TRB_v_gene, TRB_d_gene, TRB_j_gene, TRA_cdr3, TRB_cdr3) %>%
  summarise(
    n_cells = n(),
    TRA_cdr3_nt_top      = top1(TRA_cdr3_nt),
    TRB_cdr3_nt_top      = top1(TRB_cdr3_nt),
    across(any_of(grep("^TR[AB]_cdr[12]",colnames(full_data),value=TRUE)),
           ~top1(.x), .names="{.col}__top"),
    .groups = "drop"
  ) %>%
  group_by(patient) %>%
  mutate(n_tot_paz=sum(n_cells), freq_rel=round(n_cells/n_tot_paz,4)) %>%
  ungroup()

# Identifica condivisi
shared_strict <- clonotipi_summary %>%
  group_by(TRA_cdr3, TRB_cdr3) %>%
  filter(n_distinct(patient)>1) %>% ungroup()

shared_trb_only <- clonotipi_summary %>%
  group_by(TRB_cdr3) %>% filter(n_distinct(patient)>1) %>% ungroup()

shared_tra_only <- clonotipi_summary %>%
  group_by(TRA_cdr3) %>% filter(n_distinct(patient)>1) %>% ungroup()

# Inizializza oggetti per evitare errori nel salvataggio se non ci sono condivisi
confronto_CDR  <- data.frame(nota="nessun clonotipo condiviso")
confronto_geni <- data.frame(nota="nessun clonotipo condiviso")
barcode_check  <- data.frame(nota="nessun barcode condiviso")

if (nrow(shared_strict) > 0) {
  
  shared_cells <- full_data %>%
    filter(Clone_Quality=="Complete") %>%
    semi_join(shared_strict, by=c("TRA_cdr3","TRB_cdr3"))
  
  # FIX PIVOT: Escludiamo le colonne ID dal pivot per non farle sparire
  cols_to_pivot <- setdiff(intersect(all_cdr_cols, colnames(shared_cells)), c("TRA_cdr3", "TRB_cdr3"))
  
  confronto_CDR <- shared_cells %>%
    select(patient, stage, clean_barcode, folder_name, TRA_cdr3, TRB_cdr3, all_of(cols_to_pivot)) %>%
    pivot_longer(cols = all_of(cols_to_pivot), names_to="regione", values_to="sequenza") %>%
    mutate(
      catena = str_extract(regione,"^TR[AB]"),
      tipo = case_when(
        str_detect(regione,"cdr1") & !str_detect(regione,"_nt") ~ "CDR1_aa",
        str_detect(regione,"cdr1") &  str_detect(regione,"_nt") ~ "CDR1_nt",
        str_detect(regione,"cdr2") & !str_detect(regione,"_nt") ~ "CDR2_aa",
        str_detect(regione,"cdr2") &  str_detect(regione,"_nt") ~ "CDR2_nt",
        str_detect(regione,"cdr3") &  str_detect(regione,"_nt") ~ "CDR3_nt",
        TRUE ~ regione
      )
    ) %>%
    group_by(TRA_cdr3, TRB_cdr3, catena, tipo, regione) %>%
    summarise(
      n_pazienti = n_distinct(patient),
      valori_per_paziente = paste(unique(paste(patient,stage,sequenza,sep=":")), collapse=" || "),
      IDENTICA_tra_paz = n_distinct(sequenza[!is.na(sequenza)])==1,
      .groups="drop"
    )
  
  confronto_geni <- shared_cells %>%
    select(patient, stage, TRA_cdr3, TRB_cdr3, all_of(intersect(gene_cols, colnames(.)))) %>%
    pivot_longer(cols = all_of(intersect(gene_cols, colnames(.))), names_to="tipo_gene", values_to="gene") %>%
    group_by(TRA_cdr3, TRB_cdr3, tipo_gene) %>%
    summarise(
      valori_per_paziente = paste(unique(paste(patient,stage,gene,sep=":")), collapse=" vs "),
      GENI_IDENTICI = n_distinct(gene[!is.na(gene)])==1,
      .groups="drop"
    )
  
  barcode_check <- shared_cells %>%
    group_by(clean_barcode) %>%
    filter(n_distinct(patient)>1) %>% ungroup()
}

# ==============================================================================
# 7. STEP 6: DATI PLOT (invariati)
# ==============================================================================
message("\n--- STEP 6: DATI PLOT ---")

plot_ready_Vgene <- full_data %>%
  inner_join(top_clones_Vgene %>% select(patient,Gene_Label,Rank_ID), by=c("patient","Gene_Label")) %>%
  group_by(patient,Rank_ID,Gene_Label,stage) %>%
  summarise(n_cells=n(), .groups="drop") %>%
  complete(nesting(patient,Rank_ID,Gene_Label), stage=c("I","A","B"), fill=list(n_cells=0)) %>%
  filter(!(patient=="Me" & stage=="A")) %>%
  mutate(stage=factor(stage,levels=c("I","A","B")),
         Rank_ID=factor(Rank_ID,levels=paste0("Clone ",10:1)))

# Indica l'origine del clone (top-I, top-B, o entrambi) per annotare il plot
clone_origin <- bind_rows(
  top_I_clones %>% mutate(in_topI=TRUE),
  top_B_clones %>% mutate(in_topB=TRUE)
) %>%
  group_by(patient, Clone_ID_CDR3) %>%
  summarise(in_topI=any(!is.na(in_topI)), in_topB=any(!is.na(in_topB)), .groups="drop") %>%
  mutate(origin = case_when(
    in_topI & in_topB ~ "[I+B]",
    in_topI           ~ "[I]",
    TRUE              ~ "[B]"
  ))

n_cloni_per_paz <- top_clones_CDR3 %>% group_by(patient) %>% summarise(n=n())

plot_ready_CDR3 <- clean_data_plot %>%
  inner_join(top_clones_CDR3 %>% select(patient,Clone_ID_CDR3,Rank), by=c("patient","Clone_ID_CDR3")) %>%
  left_join(clone_origin %>% select(patient,Clone_ID_CDR3,origin), by=c("patient","Clone_ID_CDR3")) %>%
  group_by(patient,Rank,Gene_Label,stage,TRB_cdr3,origin) %>%
  summarise(n_cells=n(), .groups="drop") %>%
  left_join(n_cloni_per_paz, by="patient") %>%
  complete(nesting(patient,Rank,Gene_Label,TRB_cdr3,origin), stage=c("I","A","B"), fill=list(n_cells=0)) %>%
  filter(!(patient=="Me" & stage=="A")) %>%
  mutate(stage=factor(stage,levels=c("I","A","B")),
         # Label: "Clone 1 [B]  CASXXX" — origin indica se il clone è top in I, B, o entrambi
         Rank_Label=factor(
           paste0("Clone ",Rank," ",origin,"  |  ",str_trunc(TRB_cdr3,14)),
           levels=rev(sort(unique(paste0("Clone ",Rank," ",origin,"  |  ",str_trunc(TRB_cdr3,14)))))
         )
  )


# ==============================================================================
# 8. STEP 7: PLOT
# ==============================================================================
message("\n--- STEP 7: PLOT ---")

p1 <- ggplot(plot_ready_Vgene, aes(x=Rank_ID, y=n_cells, fill=stage)) +
  geom_bar(stat="identity", position=position_dodge(preserve="single"), width=0.8, color="black", linewidth=0.2) +
  facet_wrap(~patient, scales="free", ncol=1) +
  coord_flip() +
  scale_fill_manual(values=c("I"="#619CFF","A"="#F8766D","B"="#00BA38")) +
  theme_minimal() + labs(title="Top 10 Cloni (V-Gene)", y="Cellule", x="")

p2 <- ggplot(plot_ready_CDR3, aes(x=Rank_Label, y=n_cells, fill=stage)) +
  geom_bar(stat="identity", position=position_dodge(preserve="single"), width=0.8, color="black", linewidth=0.2) +
  facet_wrap(~patient, scales="free", ncol=1) +
  coord_flip() +
  scale_fill_manual(values=c("I"="#619CFF","A"="#F8766D","B"="#00BA38")) +
  theme_minimal() + labs(title="Top 10 Cloni (CDR3)", y="Cellule", x="")


# ==============================================================================
# 9. STEP 8: SALVATAGGIO
# ==============================================================================
# ==============================================================================
# 8. STEP 7: PLOT
# ==============================================================================
message("\n--- STEP 7: PLOT ---")

# Plot 1: V-Gene
p1 <- ggplot(plot_ready_Vgene, aes(x=Rank_ID, y=n_cells, fill=stage)) +
  geom_bar(stat="identity", position=position_dodge(preserve="single"),
           width=0.8, color="black", linewidth=0.2) +
  geom_text(
    data=plot_ready_Vgene %>% group_by(patient,Rank_ID) %>% slice(1) %>%
      group_by(patient) %>% mutate(y_pos=-max(n_cells, na.rm=TRUE)*0.08),
    aes(x=Rank_ID, y=y_pos, label=Gene_Label),
    hjust=1, size=3.5, color="grey20", inherit.aes=FALSE
  ) +
  facet_wrap(~patient, scales="free", ncol=1) +
  coord_flip(clip="off", ylim=c(0,NA)) +
  scale_fill_manual(values=c("I"="#619CFF","A"="#F8766D","B"="#00BA38")) +
  theme_minimal() +
  theme(plot.margin=margin(10,20,10,200), axis.text.y=element_blank(),
        axis.ticks.y=element_blank(), legend.position="bottom",
        panel.grid.minor=element_blank(), panel.grid.major.y=element_blank(),
        strip.text=element_text(face="bold",size=12)) +
  labs(title="Top 10 Cloni CAR T (Metodo V-Gene)",
       y="Numero di Cellule", x="", fill="Stage")

# Plot 2: CDR3 (Il metodo più preciso)
# Rank_Label include già il CDR3 beta troncato e il tag [I]/[B]/[I+B]
# che indica in quale stage il clone è dominante
p2 <- ggplot(plot_ready_CDR3, aes(x=Rank_Label, y=n_cells, fill=stage)) +
  geom_bar(stat="identity", position=position_dodge(preserve="single"),
           width=0.8, color="black", linewidth=0.2) +
  facet_wrap(~patient, scales="free", ncol=1) +
  coord_flip(clip="off", ylim=c(0,NA)) +
  scale_fill_manual(values=c("I"="#619CFF","A"="#F8766D","B"="#00BA38")) +
  theme_minimal() +
  theme(plot.margin=margin(10,20,10,20), axis.text.y=element_text(size=9),
        legend.position="bottom",
        panel.grid.minor=element_blank(), panel.grid.major.y=element_blank(),
        strip.text=element_text(face="bold",size=12)) +
  labs(title="Top Cloni CAR T (CDR3) — Top 5 stage I + Top 5 stage B per paziente",
       subtitle="[I] = dominante in stage I  |  [B] = dominante in stage B  |  [I+B] = presente in entrambi",
       y="Numero di Cellule", x="", fill="Stage")

# I plot vengono salvati con ggsave — print() è omesso per evitare Rplots.pdf
if (interactive()) { print(p1); print(p2) }


# ==============================================================================
# 9. STEP 8: SALVATAGGIO COMPLETO
# ==============================================================================
message("\n--- STEP 8: SALVATAGGIO ---")

# 1. Salvataggio Grafici PNG
ggsave(file.path(figures_dir,"Top10_Cloni_Vgene_vertical_CORRETTO.png"),
       p1, width=12, height=12, dpi=300, bg="white")
ggsave(file.path(figures_dir,"Top10_Cloni_CDR3_vertical_CORRETTO.png"),
       p2, width=12, height=12, dpi=300, bg="white")

# 2. Salvataggio Excel Multi-foglio (Verifica CDR e Convergenza)
write_xlsx(list(
  "01_Dati_completi"          = full_data,
  "02_Top10_CDR3"             = top_clones_CDR3,
  "03_Top10_Vgene"            = top_clones_Vgene,
  "04_Clonotipi_summary"      = clonotipi_summary,
  "05_Condivisi_strict"       = if(nrow(shared_strict)>0) shared_strict else data.frame(nota="nessuno"),
  "06_Condivisi_solo_beta"    = if(nrow(shared_trb_only)>0) shared_trb_only else data.frame(nota="nessuno"),
  "07_Condivisi_solo_alpha"   = if(nrow(shared_tra_only)>0) shared_tra_only else data.frame(nota="nessuno"),
  "08_Confronto_CDR_completo" = confronto_CDR,
  "09_Confronto_VDJ_geni"     = confronto_geni,
  "10_Barcode_check"          = barcode_check
), file.path(tables_dir,"RISULTATI_verifica_CDR_completo_CORRETTO.xlsx"))

# 3. Salvataggio file Excel singoli (Compatibilità con analisi precedenti)
write_xlsx(full_data,        file.path(tables_dir,"RISULTATI_Cloni_Dati_Completi_con_CDR3.xlsx"))
write_xlsx(top_clones_CDR3,  file.path(tables_dir,"RISULTATI_Top10_Cloni_CDR3.xlsx"))
write_xlsx(top_clones_Vgene, file.path(tables_dir,"RISULTATI_Top10_Cloni_Vgene.xlsx"))
write_xlsx(plot_ready_CDR3,  file.path(tables_dir,"RISULTATI_Plot_CDR3.xlsx"))

message("Script completato! File in: ", tables_dir, " e ", figures_dir)
