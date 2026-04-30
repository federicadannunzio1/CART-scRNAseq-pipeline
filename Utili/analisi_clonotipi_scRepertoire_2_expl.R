################################################################################
#VERSIONE COMMENTATA E RESA PIù LEGGIBILE DA GPT

# PIPELINE: Analisi clonotipi (VDJ + Seurat) e visualizzazione
#
# OBIETTIVO
#   1) Leggere i contig VDJ e combinare i dati con scRepertoire
#   2) Assegnare un clonotipo (CTstrict) ad ogni barcode
#   3) Mappare i clonotipi sugli oggetti Seurat (sc_obj_list_out)
#   4) Trovare clonotipi comuni tra coppie di campioni (I vs A/B)
#   5) Visualizzare:
#        - clonalCompare (alluvial plot)
#        - UMAP colorate per clonotipo
#        - UMAP che evidenziano solo i clonotipi conservati
################################################################################

##############################################
# 0. SETUP: working directory e pacchetti
##############################################

setwd("D:/Progetti/Ignazio_2/")

# Installazioni (lasciate commentate per riferimento)
# BiocManager::install("immApex")
# BiocManager::install("scRepertoire")
# remotes::install_github(c("BorchLab/immApex", "BorchLab/scRepertoire"))

library(scRepertoire)
library(data.table)
library(digest)
library(Seurat)
library(SeuratObject)

##############################################
# 1. FUNZIONE UTILE: codifica dei clonotipi
##############################################

# generate_clonotype_codes:
#   - prende un vettore di clonotipi (stringhe)
#   - rimuove spazi iniziali/finali
#   - trasforma NA in stringa vuota
#   - per ciascun clonotipo calcola un hash (es. SHA1)
#   - restituisce un codice tipo "CLN_xxxxxxxxxx"
generate_clonotype_codes <- function(clonotypes,
                                     algo   = "sha1",
                                     n      = 10,
                                     prefix = "CLN_") {
  clonotypes <- trimws(clonotypes)
  clonotypes[is.na(clonotypes)] <- ""
  
  sapply(clonotypes, function(s) {
    h <- digest(s, algo = algo, serialize = FALSE)
    paste0(prefix, substr(h, 1, n))
  }, USE.NAMES = FALSE)
}

##############################################
# 2. CARICAMENTO CONTIG VDJ (10x)
##############################################

# Ogni oggetto è il filtered_contig_annotations.csv di 10x
S151_I <- fread("Data/Preprocessing/1/S151_I/vdj_t/filtered_contig_annotations.csv")
S345_I <- fread("Data/Preprocessing/1/S345_I/vdj_t/filtered_contig_annotations.csv")
S393_A <- fread("Data/Preprocessing/1/S393_A/vdj_t/filtered_contig_annotations.csv")
S393_B <- fread("Data/Preprocessing/1/S393_B/vdj_t/filtered_contig_annotations.csv")
S429_I <- fread("Data/Preprocessing/2/S429_I/vdj_t/filtered_contig_annotations.csv")
S431_B <- fread("Data/Preprocessing/2/S431_B/vdj_t/filtered_contig_annotations.csv")
S435_A <- fread("Data/Preprocessing/2/S435_A/vdj_t/filtered_contig_annotations.csv")
S435_B <- fread("Data/Preprocessing/2/S435_B/vdj_t/filtered_contig_annotations.csv")

# Gruppi di sample per paziente
Ca_samples <- c("S151_I", "S393_A", "S393_B")
Me_samples <- c("S345_I", "S431_B")
Bo_samples <- c("S429_I", "S435_A", "S435_B")

# Tutti i sample
All_samples <- c("S151_I", "S345_I", "S393_A", "S393_B",
                 "S429_I", "S431_B", "S435_A", "S435_B")

##############################################
# 3. CARICO GLI OGGETTI SEURAT (con VDJ immunos)
##############################################

in_rds <- "D:/Progetti/Ignazio_2/Results/seurat_samples_sctype_azimuth_pbmc_bonemarrow_clonalvdj_immunos.rds"
sc_obj_list <- readRDS(in_rds)

# Subset per paziente (se servono)
sc_obj_list_Ca <- sc_obj_list[Ca_samples]
sc_obj_list_Me <- sc_obj_list[Me_samples]
sc_obj_list_Bo <- sc_obj_list[Bo_samples]

##############################################
# 4. LISTA VDJ PER scRepertoire (contig.list)
##############################################

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
# 5. COMBINAZIONE VDJ CON scRepertoire (combineTCR)
##############################################

# combined.TCR_all:
#   - lista, un elemento per sample
#   - contiene almeno: barcode, CTstrict, CTaa, ecc.
combined.TCR_all <- combineTCR(
  contig.list,
  samples     = All_samples,
  removeNA    = TRUE,   # rimuove contig senza clonotipo valido
  removeMulti = FALSE,
  filterMulti = FALSE
)

# Opzionale: ricodifica dei clonotipi CTstrict con hash
# (mantiene relazioni ma li rende più compatti/uniformi)
combined.TCR_all <- lapply(combined.TCR_all, function(x) {
  x$CTstrict <- generate_clonotype_codes(x$CTstrict)
  x
})

# Subset per paziente (per comodità nei plot)
combined.TCR_Ca <- combined.TCR_all[Ca_samples]
combined.TCR_Me <- combined.TCR_all[Me_samples]
combined.TCR_Bo <- combined.TCR_all[Bo_samples]

##############################################
# 6. MAPPATURA: barcode → clonotipo (all_barcode_clonotype)
##############################################

# all_barcode_clonotype:
#   - lista, un data.frame per sample
#   - colonne: barcode, clonotype
#   - rownames = barcode (per merge rapido con Seurat)
all_barcode_clonotype <- lapply(combined.TCR_all, function(x) {
  clono_col <- "CTstrict"  # colonna di clonotipo da usare
  
  df <- data.frame(
    barcode   = x$barcode,
    clonotype = x[[clono_col]],
    stringsAsFactors = FALSE
  )
  
  df <- unique(df)
  rownames(df) <- df$barcode
  df
})

# Allinea i nomi agli stessi dei combined.TCR_all
names(all_barcode_clonotype) <- names(combined.TCR_all)

##############################################
# 7. AGGIUNTA CLONOTIPO AGLI OGGETTI SEURAT
##############################################

# Risultato: sc_obj_list_out
#   - lista di oggetti Seurat
#   - aggiunta colonna meta.data 'clonotype' (CTstrict) per i barcodes presenti
sc_obj_list_out <- lapply(All_samples, function(smp) {
  
  sc_obj   <- sc_obj_list[[smp]]
  bc_clono <- all_barcode_clonotype[[smp]]
  
  # inizializza metadato clonotype
  sc_obj$clonotype <- NA_character_
  
  # trova barcodes in comune tra Seurat e VDJ
  common <- intersect(colnames(sc_obj), rownames(bc_clono))
  
  if (length(common) > 0) {
    sc_obj@meta.data[common, "clonotype"] <- bc_clono[common, "clonotype"]
  }
  
  sc_obj
})
names(sc_obj_list_out) <- All_samples

################################################################################
# 8. DEFINIZIONE DELLE COPPIE DI CAMPIONI (I vs A/B)
################################################################################

Ca_samples_blood <- c("S151_I", "S393_A")
Ca_samples_bone  <- c("S151_I", "S393_B")
Me_samples_bone  <- c("S345_I", "S431_B")
Bo_samples_blood <- c("S429_I", "S435_A")
Bo_samples_bone  <- c("S429_I", "S435_B")

# Lista con tutte le coppie che vogliamo confrontare
Coppie_samples <- list(
  Ca_samples_blood = Ca_samples_blood,
  Ca_samples_bone  = Ca_samples_bone,
  Me_samples_bone  = Me_samples_bone,
  Bo_samples_blood = Bo_samples_blood,
  Bo_samples_bone  = Bo_samples_bone
)

################################################################################
# 9. CONTA DEI CLONOTIPI COMUNI (conta_clono)
################################################################################

# conta_clono:
#   - lista (stessa lunghezza di Coppie_samples)
#   - ogni elemento è una lista con:
#       [[ sample_I  ]] = data.frame(count clonotipo in I)
#       [[ sample_AB ]] = data.frame(count clonotipo in A/B)
#   - Var1 = clonotipo, Freq = conteggio
conta_clono <- lapply(Coppie_samples, function(x) {
  # x è un vettore c(sample_I, sample_AB)
  nomi <- x
  
  comparable <- all_barcode_clonotype[x]
  sample_1   <- comparable[[1]]  # baseline / I
  sample_2   <- comparable[[2]]  # timepoint / A-B
  
  # clonotipi presenti in entrambi i campioni
  clonotipi_comuni <- intersect(sample_1$clonotype, sample_2$clonotype)
  
  # filtra le righe per clonotipi comuni
  sample_I  <- sample_1[sample_1$clonotype %in% clonotipi_comuni, ]
  sample_AB <- sample_2[sample_2$clonotype %in% clonotipi_comuni, ]
  
  # tabelle di frequenza per clonotipo
  freq_I  <- as.data.frame(sort(table(sample_I$clonotype),  decreasing = TRUE))
  freq_AB <- as.data.frame(sort(table(sample_AB$clonotype), decreasing = TRUE))
  
  res <- list(freq_I, freq_AB)
  names(res) <- nomi
  res
})

################################################################################
# 10. PLOT scRepertoire: clonalCompare (alluvial)
################################################################################

# Usiamo le stesse coppie definite sopra

p_BO_blood <- clonalCompare(
  combined.TCR_Bo,
  top.clones = 10,
  samples    = Bo_samples_blood,
  cloneCall  = "strict",
  graph      = "alluvial"
) + ggtitle("Bo_samples_blood")

p_BO_bone <- clonalCompare(
  combined.TCR_Bo,
  top.clones = 10,
  samples    = Bo_samples_bone,
  cloneCall  = "strict",
  graph      = "alluvial"
) + ggtitle("Bo_samples_bone")

p_ME_bone <- clonalCompare(
  combined.TCR_Me,
  top.clones = 10,
  samples    = Me_samples_bone,
  cloneCall  = "strict",
  graph      = "alluvial"
) + ggtitle("Me_samples_bone")

p_CA_blood <- clonalCompare(
  combined.TCR_Ca,
  top.clones = 10,
  samples    = Ca_samples_blood,
  cloneCall  = "strict",
  graph      = "alluvial"
) + ggtitle("Ca_samples_blood")

p_CA_bone <- clonalCompare(
  combined.TCR_Ca,
  top.clones = 10,
  samples    = Ca_samples_bone,
  cloneCall  = "strict",
  graph      = "alluvial"
) + ggtitle("Ca_samples_bone")

# Visualizzazione (in un ambiente interattivo)
p_CA_blood
p_CA_bone
p_ME_bone
p_BO_blood
p_BO_bone

################################################################################
# 11. UMAP COLORATI PER CLONOTIPO (tutte le cellule)
################################################################################

# plot_umap_clono:
#   - lista di plot UMAP (uno per sample)
#   - colore = clonotype (CTstrict)
plot_umap_clono <- lapply(names(sc_obj_list_out), function(smp) {
  obj <- sc_obj_list_out[[smp]]
  DimPlot(obj, reduction = "umap", group.by = "clonotype") +
    ggtitle(smp) +
    theme_minimal()
})
names(plot_umap_clono) <- names(sc_obj_list_out)

# esempio
plot_umap_clono[[1]]

# Versione senza legenda (utile per figure multiple affiancate)
plot_umap_clono_no_legend <- lapply(names(sc_obj_list_out), function(smp) {
  obj <- sc_obj_list_out[[smp]]
  DimPlot(obj, reduction = "umap", group.by = "clonotype") +
    ggtitle(smp) +
    theme_minimal() +
    theme(legend.position = "none")
})
names(plot_umap_clono_no_legend) <- names(sc_obj_list_out)

# esempi
plot_umap_clono_no_legend[[1]]
plot_umap_clono_no_legend[[2]]

################################################################################
# 12. UMAP DEI SOLI CLONOTIPI COMUNI (conserved clones)
################################################################################

# Obiettivo:
#   - Per ogni coppia (I vs A/B) mostrare solo i clonotipi conservati
#   - Tutte le altre cellule vengono etichettate come "Other"
#   - Si ottengono due UMAP per coppia: una per I, una per A/B

plot_umap_conserved_clones <- lapply(conta_clono, function(clonos) {
  # clonos è la lista per una coppia: lista[[sample_I]], lista[[sample_AB]]
  nomi <- names(clonos)  # c("sample_I", "sample_AB")
  
  # Estrai i clonotipi da plottare (chiave = Var1 nella tabella)
  clonotipi_da_plottare_I  <- sort(as.character(clonos[[1]]$Var1))
  clonotipi_da_plottare_AB <- sort(as.character(clonos[[2]]$Var1))
  
  # Oggetti Seurat per I e A/B
  obj_I  <- sc_obj_list_out[[nomi[1]]]
  obj_AB <- sc_obj_list_out[[nomi[2]]]
  
  # Crea una nuova colonna "clono_highlight":
  #   - clonotipo specifico se appartiene ai clonotipi da plottare
  #   - "Other" altrimenti
  obj_I$clono_highlight <- ifelse(
    !is.na(obj_I$clonotype) & obj_I$clonotype %in% clonotipi_da_plottare_I,
    obj_I$clonotype,
    "Other"
  )
  
  obj_AB$clono_highlight <- ifelse(
    !is.na(obj_AB$clonotype) & obj_AB$clonotype %in% clonotipi_da_plottare_AB,
    obj_AB$clonotype,
    "Other"
  )
  
  # Ordina i livelli del fattore (clonotipi specifici + "Other" in coda)
  obj_I$clono_highlight <- factor(
    obj_I$clono_highlight,
    levels = c(sort(unique(clonotipi_da_plottare_I)), "Other")
  )
  
  obj_AB$clono_highlight <- factor(
    obj_AB$clono_highlight,
    levels = c(sort(unique(clonotipi_da_plottare_AB)), "Other")
  )
  
  # UMAP:
  #   - tutte le cellule mostrate
  #   - clonotipi conservati colorati distintamente
  #   - "Other" come categoria residuale
  plot_umap_I <- DimPlot(
    obj_I,
    reduction = "umap",
    group.by  = "clono_highlight"
  ) +
    ggtitle(nomi[1]) +
    theme_minimal()
  
  plot_umap_AB <- DimPlot(
    obj_AB,
    reduction = "umap",
    group.by  = "clono_highlight"
  ) +
    ggtitle(nomi[2]) +
    theme_minimal()
  
  res <- list(plot_umap_I, plot_umap_AB)
  names(res) <- nomi
  res
})

# Esempi di visualizzazione (per tutte le coppie)
plot_umap_conserved_clones[[1]][[1]]
plot_umap_conserved_clones[[1]][[2]]
plot_umap_conserved_clones[[2]][[1]]
plot_umap_conserved_clones[[2]][[2]]
plot_umap_conserved_clones[[3]][[1]]
plot_umap_conserved_clones[[3]][[2]]
plot_umap_conserved_clones[[4]][[1]]
plot_umap_conserved_clones[[4]][[2]]
plot_umap_conserved_clones[[5]][[1]]
plot_umap_conserved_clones[[5]][[2]]

################################################################################
# FINE SCRIPT
################################################################################
