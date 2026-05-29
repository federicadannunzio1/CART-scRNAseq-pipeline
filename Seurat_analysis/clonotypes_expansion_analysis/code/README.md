# Pipeline Analisi Clonotipi CAR-T

**Progetto:** caruana-project / CART\
**Autore:** Federica D'Annunzio

------------------------------------------------------------------------

## Schema generale

```         
INPUT
├── seurat_samples_...CAR.rds          ← cellule CAR+ con metadati
└── output_allineamento_.../
    └── {1,2}/{campione}/vdj_t/
        └── filtered_contig_annotations.csv

PIPELINE
         ┌─────────────────────────────┐
         │   01_build_clonotypes.R     │  ~5-10 min
         │   Costruisce full_data      │
         └──────────────┬──────────────┘
                        │ full_data (in memoria R)
          ┌─────────────┼─────────────┬──────────────┐
          ▼             ▼             ▼              ▼
   02_vdjdb_     03_private_   04_expansion_   (futuri script)
   search.R      clones.R      dynamics.R
   Cerca CDR3    Cloni privati  Fold-change
   in database   vs condivisi   I → A → B
   pubblico
```

------------------------------------------------------------------------

## Script da tenere (4 script attivi)

| \# | File | Funzione | Dipende da |
|----|----|----|----|
| 1 | `01_build_clonotypes.R` | Costruisce full_data con TCR puliti | — |
| 2 | `02_vdjdb_search.R` | Cerca CDR3 in VDJdb | 01 |
| 3 | `03_private_clones.R` | Cloni privati vs condivisi + plot famiglie | 01 |
| 4 | `04_expansion_dynamics.R` | Espansione I→A→B, fold-change, heatmap | 01 |

## Script da archiviare (NON usare)

| File | Motivo |
|----|----|
| `1_find_chains_colonotypes.R` | **BUG CRITICO**: `paste(collapse="/")` concatena CDR3 di catene doppie producendo cloni ibridi tipo `TRBV7-9/TRBV6-4`. Filtri incompleti (mancano `is_cell`, `high_confidence`, `full_length`). |
| `2_fixed_plot_unique_tcr.R` | Usa `library(xlsx)` deprecato. Legge da file invece che da memoria. Colonne in maiuscolo incompatibili con `01_build_clonotypes.R`. |
| `2_find_cd3_in_database.R` | Colonne in maiuscolo (`TRB_CDR3`) incompatibili. Nessun output xlsx riassuntivo. |

------------------------------------------------------------------------

## Ordine di esecuzione

``` r
# In RStudio, esegui in ordine:
source("01_build_clonotypes.R")   # produce full_data — eseguire SEMPRE per primo
source("02_vdjdb_search.R")       # richiede connessione internet
source("03_private_clones.R")
source("04_expansion_dynamics.R")
```

> **Nota:** `full_data` prodotto da `01_` rimane in memoria R per tutti gli script successivi. Se chiudi RStudio, devi rieseguire `01_` prima degli altri.

------------------------------------------------------------------------

## Dettaglio di ogni script

### 01_build_clonotypes.R

**Cosa fa:**

Legge i due input principali (Seurat + VDJ) e costruisce `full_data`, la tabella centrale di tutta l'analisi. Ogni riga è una cellula CAR+ con il suo TCR completo.

**Filtri applicati sul VDJ (con tracciamento righe):**

| Filtro                    | Scopo                                  |
|---------------------------|----------------------------------------|
| `is_cell == TRUE`         | Rimuove droplet vuoti                  |
| `high_confidence == TRUE` | Rimuove mapping ambigui                |
| `productive == TRUE`      | Rimuove riarrangiamenti non funzionali |
| `full_length == TRUE`     | Assicura sequenza V-J completa         |
| `chain TRA/TRB`           | Esclude catene non TCR                 |
| `CDR3 valida`             | Rimuove CDR3 mancanti o "None"         |

**Fix bug catene doppie:**\
Cellule con \>1 TRA o \>1 TRB → seleziona il contig con **max UMI** (determinismo garantito). Elimina la concatenazione con `/` del vecchio script.

**Oggetti R prodotti (in memoria):** - `full_data` — tabella principale, una riga per cellula CAR+ - `top_clones_CDR3` — top 10 per paziente - `top_clones_Vgene` — top 10 per paziente (metodo V-gene) - `shared_strict` — cloni con CDR3 identica in \>1 paziente

**File salvati:**

```         
res/
├── Top10_Cloni_CDR3_vertical_CORRETTO.png
├── Top10_Cloni_Vgene_vertical_CORRETTO.png
├── RISULTATI_Cloni_Dati_Completi_con_CDR3.xlsx   ← full_data completo
├── RISULTATI_Top10_Cloni_CDR3.xlsx
├── RISULTATI_Top10_Cloni_Vgene.xlsx
└── RISULTATI_verifica_CDR_completo_CORRETTO.xlsx
    ├── 01_Dati_completi
    ├── 05_Condivisi_strict
    ├── 08_Confronto_CDR_completo   ← IDENTICA_tra_paz = TRUE/FALSE
    ├── 09_Confronto_VDJ_geni
    └── 10_Barcode_check
```

------------------------------------------------------------------------

### 02_vdjdb_search.R

**Cosa fa:**

Scarica automaticamente VDJdb (database pubblico di TCR con specificità antigenica nota) e cerca se le CDR3 dei tuoi cloni vi compaiono.

**Logica di ricerca:** - Scarica `vdjdb-2024-06-13.zip` da GitHub - Filtra per `species == "HomoSapiens"` - Inner join tra le tue CDR3 e il database, separatamente per alpha e beta - Se non trova nulla → cloni "privati" (specifici per HLA del paziente o per antigeni non ancora catalogati)

**File salvati:**

```         
res/
└── VDJdb_Summary.xlsx
    ├── Match_Beta    ← CDR3 beta trovate nel DB con antigene associato
    ├── Match_Alpha   ← CDR3 alpha trovate nel DB
    └── CDR3_cercate  ← lista completa delle CDR3 cercate
```

------------------------------------------------------------------------

### 03_private_clones.R

**Cosa fa:**

Separa il repertorio in cloni "condivisi" (stesso riarrangiamento nucleotidico in \>1 paziente, probabilmente dallo stesso lotto CAR-T) e cloni "privati" (unici per paziente, vera espansione biologica individuale).

**Classificazione:** - **Condiviso**: stessa CDR3 aa **e** stessa CDR3 nt in \>1 paziente → stesso clone - **Privato**: tutto il resto → espansione biologica reale

Poi plotta le famiglie TCR beta di interesse (`target_families`) nei soli cloni privati.

**File salvati:**

```         
res/
├── Grafico_Famiglie_Target_Privato.png
└── REPORT_CONTAMINAZIONE_E_PRIVATI.xlsx
    ├── Cloni_condivisi     ← cloni da stesso lotto CAR-T
    ├── Dati_privati        ← full_data filtrato
    ├── Top10_cloni_privati
    └── Famiglie_target
```

------------------------------------------------------------------------

### 04_expansion_dynamics.R

**Cosa fa:**

Risponde alla domanda biologica principale: *quali cloni T si espandono in risposta al CAR-T?*

**Calcoli per ogni clone:** - `n_cells` per stage I, A, B - `freq` = frequenza relativa (normalizzata per profondità di sequenziamento) - `FC_I_to_B` = `freq_B / freq_I` (fold-change)

**Classificazione:**

| Categoria | Criterio | Interpretazione |
|----|----|----|
| De novo in B | Assente in I, presente in B | Espansione in vivo dopo trattamento |
| Espanso FC≥2 | `freq_B / freq_I ≥ 2` | Selezionato positivamente dal trattamento |
| Stabile | FC tra 1 e 2 | Presente ma non risponde |
| Contratto | FC \< 1 | Si riduce → possibile esaurimento |

**File salvati:**

```         
res/
├── Expansion_dynamics_lineplot.png   ← andamento I→A→B top 10 cloni
├── Expansion_heatmap.png             ← frequenza relativa per clone × campione
└── RISULTATI_expansion_dynamics.xlsx
    ├── 01_Dinamica_completa          ← tutti i cloni con FC e categoria
    ├── 02_Cloni_espansi_in_B
    ├── 03_Espansi_multi_paz          ← espansi in >1 paziente
    └── 04_Tabella_pubblicazione      ← pronta per paper
```

------------------------------------------------------------------------

## Struttura colonne di full_data

**Metadati cellula:**

| Colonna         | Descrizione                                 |
|-----------------|---------------------------------------------|
| `clean_barcode` | Barcode 10x                                 |
| `folder_name`   | Cartella CellRanger                         |
| `patient`       | Bo / Ca / Me                                |
| `stage`         | I / A / B                                   |
| `Clone_Quality` | Complete / TRA_only / TRB_only / Incomplete |

**Catena Alpha (TRA):**

| Colonna                    | Descrizione                         |
|----------------------------|-------------------------------------|
| `TRA_v_gene`               | es. TRAV17                          |
| `TRA_j_gene`               | es. TRAJ29                          |
| `TRA_cdr1` / `TRA_cdr1_nt` | CDR1 aa e nucleotidica              |
| `TRA_cdr2` / `TRA_cdr2_nt` | CDR2 aa e nucleotidica              |
| `TRA_cdr3` / `TRA_cdr3_nt` | CDR3 aa e nt — chiave del clonotipo |
| `TRA_umis`                 | UMI del contig selezionato          |

**Catena Beta (TRB):**\
Stesse colonne con prefisso `TRB_`, in più `TRB_d_gene`.

**Identificatori clone:**

| Colonna         | Esempio                                            |
|-----------------|----------------------------------------------------|
| `Gene_Label`    | `TRAV17 + TRBV7-9`                                 |
| `Clone_ID_CDR3` | `CATPVRGNTPLVF_CASSSTGWDSPYNYGYTF`                 |
| `Clone_ID_Full` | `TRAV17:TRAJ29:CATPVR..._TRBV7-9:TRBJ1-2:CASSS...` |

------------------------------------------------------------------------

## Interpretazione risultato chiave (foglio 08)

Il foglio `08_Confronto_CDR_completo` nel file `RISULTATI_verifica_CDR_completo_CORRETTO.xlsx` ha una colonna `IDENTICA_tra_paz`:

| Valore | Interpretazione |
|----|----|
| `TRUE` | Stessa sequenza nucleotidica in tutti i pazienti → stesso lotto CAR-T |
| `FALSE` | Sequenze diverse → convergenza immunologica reale |

Nei tuoi dati attuali tutti i cloni condivisi tra Bo e Me hanno `IDENTICA_tra_paz = TRUE` su CDR1, CDR2, CDR3 nt e geni V/J → provengono dallo stesso prodotto CAR-T, non da risposte immuni indipendenti.
