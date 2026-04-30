################################################################################
#VERSIONE COMMENTATA E RESA PIù LEGGIBILE DA GPT

# PROJECT: CAR-T / VDJ Clonotype Integration Pipeline
# DESCRIPTION:
#   Pipeline modulare per combinare dati VDJ (scRepertoire) con dati singola cellula (Seurat),
#   identificare cloni espansi, cloni condivisi, cellule CAR derivate per clonotipo
#   e barcodes provenienti da scRepertoire.
#
# NOTE:
#   - Segue al 100% la logica del tuo script originale.
#   - Tutto è organizzato in funzioni, modulare e condivisibile.
#   - Commenti altamente esplicativi.
################################################################################


##############################################
# 📦 0. Setup ambiente e pacchetti
##############################################

setwd("D:/Progetti/Ignazio_2/")

library(scRepertoire)
library(data.table)
library(digest)
library(Seurat)
library(SeuratObject)
library(patchwork)
library(dplyr)


##############################################
# 📂 1. FUNZIONI UTILI
##############################################

# ------------------------------------------------------------------------------
# Funzione 1:
#   Genera codici hash univoci per clonotipo
#   (utile per evitare conflitti e standardizzare i nomi dei clonotipi)
# ------------------------------------------------------------------------------
generate_clonotype_codes <- function(clonotypes, algo = "sha1", n = 10, prefix = "CLN_") {
  clonotypes <- trimws(clonotypes)
  clonotypes[is.na(clonotypes)] <- ""
  
  sapply(clonotypes, function(s) {
    h <- digest(s, algo = algo, serialize = FALSE)
    paste0(prefix, substr(h, 1, n))
  }, USE.NAMES = FALSE)
}


##############################################
# 📂 2. CARICAMENTO VDJ + DEFINIZIONE GRUPPI
##############################################

# ---------------------------
# Carico file VDJ dei campioni
# ---------------------------
S151_I <- fread("Data/Preprocessing/1/S151_I/vdj_t/filtered_contig_annotations.csv")
S345_I <- fread("Data/Preprocessing/1/S345_I/vdj_t/filtered_contig_annotations.csv")
S393_A <- fread("Data/Preprocessing/1/S393_A/vdj_t/filtered_contig_annotations.csv")
S393_B <- fread("Data/Preprocessing/1/S393_B/vdj_t/filtered_contig_annotations.csv")
S429_I <- fread("Data/Preprocessing/2/S429_I/vdj_t/filtered_contig_annotations.csv")
S431_B <- fread("Data/Preprocessing/2/S431_B/vdj_t/filtered_contig_annotations.csv")
S435_A <- fread("Data/Preprocessing/2/S435_A/vdj_t/filtered_contig_annotations.csv")
S435_B <- fread("Data/Preprocessing/2/S435_B/vdj_t/filtered_contig_annotations.csv")

# Gruppi di campioni
Ca_samples <- c("S151_I", "S393_A", "S393_B")
Me_samples <- c("S345_I", "S431_B")
Bo_samples <- c("S429_I", "S435_A", "S435_B")

All_samples <- c(Ca_samples, Me_samples, Bo_samples)
samples <- All_samples

# Coppie per confronto tra sangue ↔ midollo
coppie_samples <- list(
  Ca_samples_blood = c("S151_I", "S393_A"),
  Ca_samples_bone  = c("S151_I", "S393_B"),
  Bo_samples_blood = c("S429_I", "S435_A"),
  Bo_samples_bone  = c("S429_I", "S435_B"),
  Me_samples_bone = c("S345_I", "S431_B")
)

# Creo lista VDJ
contig.list <- list(
  S151_I = S151_I,
  S345_I = S345_I,
  S393_A = S393_A,
  S393_B = S393_B,
  S429_I = S429_I,
  S431_B = S431_B,
  S435_A = S435_A,
  S435_B = S435_B
)


##############################################
# 📂 3. CARICO OGGETTI SEURAT CON METADATI VDJ
##############################################

sc_obj_list <- readRDS("D:/Progetti/Ignazio_2/Results/seurat_samples_sctype_azimuth_pbmc_bonemarrow_clonalvdj_immunos.rds")

sc_obj_list_Ca <- sc_obj_list[Ca_samples]
sc_obj_list_Me <- sc_obj_list[Me_samples]
sc_obj_list_Bo <- sc_obj_list[Bo_samples]


##############################################
# 📂 4. COMBINAZIONE VDJ CON scRepertoire
##############################################

combined.TCR_all <- combineTCR(
  contig.list,
  samples     = All_samples,
  removeNA    = TRUE,
  removeMulti = FALSE,
  filterMulti = FALSE
)

# Normalizzo i clonotipi (mantiene la logica originale)
combined.TCR_all <- lapply(combined.TCR_all, function(x) {
  x$CTstrict <- x$CTstrict
  x
})


##############################################
# 📂 5. ANALISI CLONALE: FREQUENZA & CONSERVAZIONE
##############################################

# ------------------------------------------------------------------------------
# 5.1 – Frequenza dei clonotipi in ogni campione
# ------------------------------------------------------------------------------
cloni_all_freq <- lapply(samples, function(smp) {
  
  clonotype <- na.omit(combined.TCR_all[[smp]]$CTstrict)
  
  tab <- table(clonotype)
  
  df <- data.frame(
    clonotype = names(tab),
    count     = as.integer(tab)
  )
  
  df
})
names(cloni_all_freq) <- samples


# ------------------------------------------------------------------------------
# 5.2 – Clonotipi condivisi tra campioni "accoppiati"
# ------------------------------------------------------------------------------
cloni_conservati_all <- lapply(coppie_samples, function(pair) {
  
  s1 <- pair[1]
  s2 <- pair[2]
  
  cln1 <- unique(combined.TCR_all[[s1]]$CTstrict)
  cln2 <- unique(combined.TCR_all[[s2]]$CTstrict)
  
  intersect(cln1, cln2)
})


# ------------------------------------------------------------------------------
# 5.3 – Frequenze dei clonotipi conservati
# ------------------------------------------------------------------------------
cloni_conservati_all_freq <- lapply(coppie_samples, function(pair) {
  
  s1 <- pair[1]
  s2 <- pair[2]
  
  df1 <- cloni_all_freq[[s1]]
  df2 <- cloni_all_freq[[s2]]
  
  comuni <- intersect(df1$clonotype, df2$clonotype)
  
  list(
    I  = df1[df1$clonotype %in% comuni, ],
    AB = df2[df2$clonotype %in% comuni, ]
  )
})


##############################################
# 📂 6. INTEGRAZIONE CON DATI SINGLE CELL
##############################################

# Carico oggetti Seurat "puri"
seu <- readRDS("D:/Progetti/Ignazio_2/Results/seurat_samples_sctype_azimuth_pbmc_bonemarrow_clonalvdj.rds")

# Carico liste CAR (da alignments gex_1/gex_2)
gex_1_CAR_list <- readRDS("Data/Post_reallineamento_2/gex_1/IS_CAR_labels_list.rds")
gex_2_CAR_list <- readRDS("Data/Post_reallineamento_2/gex_2/IS_CAR_labels_list.rds")
CAR_list <- c(gex_1_CAR_list, gex_2_CAR_list)

samples <- names(seu)


# ------------------------------------------------------------------------------
# 6.1 – Re-inserisco metadato CAR "diretto"
# ------------------------------------------------------------------------------
seu <- lapply(samples, function(smp) {
  CAR <- CAR_list[[smp]]
  obj <- seu[[smp]]
  
  obj$CAR <- "NO"
  obj$CAR[CAR$cell[CAR$IS_CAR == "YES"]] <- "YES"
  
  obj
})
names(seu) <- samples


# ------------------------------------------------------------------------------
# 6.2 – Identifico cellule CAR “derivate per clonotipo”
#       (cellule non etichettate come CAR, ma che condividono clonotipo con CAR)
# ------------------------------------------------------------------------------
seu <- lapply(samples, function(smp) {
  
  obj <- seu[[smp]]
  
  CAR_clonotypes <- na.omit(unique(obj$clonal[obj$CAR == "YES"]))
  
  obj$derived_CAR <- "NO"
  obj$derived_CAR[obj$clonal %in% CAR_clonotypes] <- "YES"
  
  obj
})
names(seu) <- samples


# ------------------------------------------------------------------------------
# 6.3 – Inserisco barcode scRepertoire → Seurat
# ------------------------------------------------------------------------------
seu <- lapply(samples, function(smp) {
  
  obj <- seu[[smp]]
  
  ct <- na.omit(combined.TCR_all[[smp]]$CTstrict)
  names(ct) <- combined.TCR_all[[smp]]$barcode
  
  obj$barcode_screpertoire <- NA
  obj$barcode_screpertoire[combined.TCR_all[[smp]]$barcode] <- ct[combined.TCR_all[[smp]]$barcode]
  
  obj
})
names(seu) <- samples


##############################################
# 📂 7. COSTRUZIONE OGGETTI PER COPPIE
##############################################

lista_seu_coppie <- lapply(names(coppie_samples), function(nm) {
  
  #nm <- names(coppie_samples)[[5]]
  
  s1 <- coppie_samples[[nm]][1]
  s2 <- coppie_samples[[nm]][2]
  
  cl1 <- cloni_conservati_all_freq[[nm]]$I$clonotype
  cl2 <- cloni_conservati_all_freq[[nm]]$AB$clonotype
  
  
  # Etichetta “clonotipo conservato”
  seu[[s1]]$IS_CONSERVED_scRepertoire <- "NO"
  seu[[s2]]$IS_CONSERVED_scRepertoire <- "NO"
  
  seu[[s1]]$IS_CONSERVED_scRepertoire[seu[[s1]]$barcode_screpertoire %in% cl1] <- "YES"
  seu[[s2]]$IS_CONSERVED_scRepertoire[seu[[s2]]$barcode_screpertoire %in% cl2] <- "YES"
  
  
  # Etichetta finale: CAR (derivato + conservato)
  seu[[s1]]$IS_CAR_ALLIN_scREP <- "NO"
  seu[[s2]]$IS_CAR_ALLIN_scREP <- "NO"
  
  seu[[s1]]$IS_CAR_ALLIN_scREP[
    seu[[s1]]$IS_CONSERVED_scRepertoire == "YES" |
      seu[[s1]]$derived_CAR == "YES"
  ] <- "YES"
  
  seu[[s2]]$IS_CAR_ALLIN_scREP[
    seu[[s2]]$IS_CONSERVED_scRepertoire == "YES" |
      seu[[s2]]$derived_CAR == "YES"
  ] <- "YES"
  
  list(
    I  = seu[[s1]],
    AB = seu[[s2]]
  )
})

names(lista_seu_coppie) <- names(coppie_samples)

saveRDS(object =lista_seu_coppie, file ="D:/Progetti/Ignazio_2/Results/seurat_samples_sctype_azimuth_pbmc_bonemarrow_clonalvdj_CAR_2.rds")

