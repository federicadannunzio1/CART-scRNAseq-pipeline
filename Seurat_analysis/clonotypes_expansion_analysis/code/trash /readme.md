# GUIDA COMPLETA AGLI SCRIPT CORRETTI
## Analisi Repertorio TCR delle CAR-T Cells

---

## 📋 INDICE
1. [Problemi degli Script Originali](#problemi-script-originali)
2. [Modifiche Apportate](#modifiche-apportate)
3. [Come Eseguire gli Script](#come-eseguire)
4. [Interpretazione dei Risultati](#interpretazione-risultati)
5. [Concetti Chiave da Ricordare](#concetti-chiave)
6. [Prossimi Passi](#prossimi-passi)

---

## 🔴 PROBLEMI DEGLI SCRIPT ORIGINALI

### PROBLEMA PRINCIPALE: Mancanza di CDR3

**Cosa facevano gli script originali:**
- Definivano i "cloni" usando SOLO i geni V (es. TRBV2, TRAV21)
- Un "clone" era identificato come: `TRAV21 + TRBV2`

**Perché è sbagliato:**
```
Stesso V-gene ≠ Stesso clone!

Esempio:
Cellula 1: TRBV2 + CDR3="CASSLAPGTQYF"  -> Clone A
Cellula 2: TRBV2 + CDR3="CASSLQGAYEQYF" -> Clone B (DIVERSO!)
Cellula 3: TRBV2 + CDR3="CASSLAPGTQYF"  -> Clone A (UGUALE a 1)

Gli script vecchi le avrebbero contate TUTTE come "stesso clone TRBV2"
```

**Il CDR3 (Complementarity-Determining Region 3) è:**
- La regione ipervariabile del TCR
- Quella che determina la specificità antigenica
- L'UNICA vera firma molecolare di un clone
- Lunga ~10-15 amminoacidi, con miliardi di combinazioni possibili

**Conseguenze dell'errore:**
1. ❌ I tuoi "cloni condivisi" erano in realtà FAMIGLIE di cloni diversi
2. ❌ Non potevi distinguere espansioni di un clone vs espansioni policlonali
3. ❌ I numeri di cellule erano sovrastimati (raggruppavi cloni diversi)
4. ❌ I "public clones" non erano veri cloni pubblici

---

## ✅ MODIFICHE APPORTATE

### SCRIPT 1: `1_find_chains_colonotypes_FIXED.R`

#### Modifica 1: Estrazione CDR3 dai file VDJ
**Prima:**
```r
summarise(
  TRA = paste(unique(na.omit(v_gene[chain == "TRA"])), collapse = "/"),
  TRB = paste(unique(na.omit(v_gene[chain == "TRB"])), collapse = "/")
)
```

**Dopo:**
```r
summarise(
  # V-genes
  TRA_V = paste(unique(na.omit(v_gene[chain == "TRA"])), collapse = "/"),
  TRB_V = paste(unique(na.omit(v_gene[chain == "TRB"])), collapse = "/"),
  
  # J-genes (NEW)
  TRA_J = paste(unique(na.omit(j_gene[chain == "TRA"])), collapse = "/"),
  TRB_J = paste(unique(na.omit(j_gene[chain == "TRB"])), collapse = "/"),
  
  # CDR3 (NEW - ESSENZIALE!)
  TRA_CDR3 = paste(unique(na.omit(cdr3[chain == "TRA"])), collapse = "/"),
  TRB_CDR3 = paste(unique(na.omit(cdr3[chain == "TRB"])), collapse = "/")
)
```

#### Modifica 2: Creazione ID Clone Completo
```r
# ID clone completo (V + J + CDR3)
Clone_ID_Full = paste0(
  TRA_V, ":", TRA_J, ":", TRA_CDR3,
  "_",
  TRB_V, ":", TRB_J, ":", TRB_CDR3
)

# ID clone solo CDR3 (più semplice, ma ugualmente accurato)
Clone_ID_CDR3 = paste(TRA_CDR3, TRB_CDR3, sep = "_")
```

#### Modifica 3: Flag Qualità Dati
```r
Clone_Quality = case_when(
  Has_Complete_TRA & Has_Complete_TRB ~ "Complete",      # OK per analisi
  Has_Complete_TRB & !Has_Complete_TRA ~ "TRB_only",    # Parziale
  Has_Complete_TRA & !Has_Complete_TRB ~ "TRA_only",    # Parziale
  TRUE ~ "Incomplete"                                    # Scarta
)
```

Questo ti permette di:
- Filtrare solo cloni di alta qualità (`Clone_Quality == "Complete"`)
- Capire se hai problemi tecnici nel sequenziamento TRA
- Quantificare la qualità dei tuoi dati

#### Modifica 4: Due Analisi Parallele
Lo script ora genera DUE analisi:

**Analisi 1 - V-gene based (vecchio metodo):**
- Top 10 cloni per V-gene
- Utile per confronto con analisi precedenti
- Mostra "famiglie" di cloni con stesso V-gene

**Analisi 2 - CDR3 based (metodo corretto):**
- Top 10 VERI cloni per CDR3
- Questa è l'analisi da usare per conclusioni biologiche
- Ogni clone è univoco

#### Modifica 5: Fix geom_text Duplicato
**Prima:** Avevi DUE `geom_text()` che creavano confusione
**Dopo:** UN SOLO `geom_text()` ben posizionato con label comprensibili

---

### SCRIPT 2: `2_find_single_or_coupled_shared_chains_FIXED.R`

#### Modifica 1: Distinzione Chiara tra Livelli di Condivisione

**Livello 1 - Convergenza V-gene (Pubblica):**
```r
# Usa data_for_Vgene (con catene separate)
# Risponde: "Quali V-genes sono usati da più pazienti?"
# Interpretazione: Convergenza verso segmenti genici comuni
```

**Livello 2 - Veri Public Clones (CDR3 identici):**
```r
# Usa data_for_Clones (cloni completi, no separazione)
# Risponde: "Quali cloni IDENTICI sono condivisi?"
# Interpretazione: Veri cloni pubblici, biologicamente significativi
```

#### Modifica 2: Gestione Corretta delle Catene Multiple

**Prima:**
```r
# Separava SEMPRE le catene multiple
separate_rows(TRA, sep = "/") %>%
separate_rows(TRB, sep = "/")
# Problema: conta 1 cellula come 2+ cellule se ha catene multiple
```

**Dopo:**
```r
# PER V-GENE: separa (vogliamo vedere tutti i V-genes usati)
data_for_Vgene <- full_data %>%
  separate_rows(TRA_V, sep = "/") %>%
  separate_rows(TRB_V, sep = "/")

# PER CLONI: NON separa (1 cellula = 1 clone, anche se ha 2 catene)
data_for_Clones <- full_data %>%
  filter(Clone_Quality == "Complete")
  # Nessuna separazione!
```

#### Modifica 3: Plot Separati con Messaggi Chiari

Ora generi 3 plot distinti:

1. **Plot TRBV condivisi** - Convergenza catene Beta
2. **Plot TRAV condivisi** - Convergenza catene Alpha  
3. **Plot Cloni CDR3 identici** - Veri public clones (⭐ il più importante!)

Ogni plot ha titolo e sottotitolo che spiega cosa stai guardando.

---

### SCRIPT 3: `2b_count_n_of_shared_chains_FIXED.R`

#### Modifiche Principali:

1. **Analisi a 4 Livelli:**
   - Coppie V-gene
   - Solo TRAV
   - Solo TRBV
   - ⭐ Veri cloni CDR3

2. **Metriche di Diversità Clonale:**
   ```r
   Shannon Diversity = -Σ(pi × log(pi))
   - Più alto = repertorio diversificato
   - Più basso = pochi cloni dominano
   
   Clonality = 1 - (Shannon / Shannon_max)
   - Valore 0 = massima diversità
   - Valore 1 = monoespansione
   - Valore ideale CAR-T = 0.5-0.8
   ```

3. **Jaccard Index (Overlap tra Repertori):**
   ```r
   Jaccard = |A ∩ B| / |A ∪ B|
   
   Confronta quanto due repertori si sovrappongono:
   - 0 = nessun clone condiviso
   - 1 = repertori identici
   - 0.01-0.05 = tipico per pazienti diversi
   ```

4. **Visualizzazioni Multiple:**
   - Plot 1: Overview tutti i livelli
   - Plot 2: Focus condivisione 3 pazienti
   - Plot 3: Shannon diversity per stage
   - Plot 4: Clonality index per stage
   - Plot 5: Heatmap Jaccard overlap

---

## 🚀 COME ESEGUIRE GLI SCRIPT

### Ordine di Esecuzione:

```r
# 1. Script principale (genera dati base)
source("1_find_chains_colonotypes_FIXED.R")

# A questo punto hai l'oggetto 'full_data' in memoria

# 2. Analisi condivisione qualitativa
source("2_find_single_or_coupled_shared_chains_FIXED.R")

# 3. Analisi condivisione quantitativa
source("2b_count_n_of_shared_chains_FIXED.R")
```

### File Generati:

#### Dallo Script 1:
- `RISULTATI_Cloni_Dati_Completi_con_CDR3.csv` - Tutti i dati
- `RISULTATI_Top10_Cloni_CDR3.csv` - Top cloni (metodo corretto)
- `RISULTATI_Top10_Cloni_Vgene.csv` - Top cloni (vecchio metodo, confronto)
- `RISULTATI_Plot_CDR3.csv` - Dati per il plot principale

#### Dallo Script 2:
- `RISULTATI_TRBV_condivisi_3pazienti.csv`
- `RISULTATI_TRAV_condivisi_3pazienti.csv`
- `RISULTATI_Public_Clones_CDR3_identici.csv` ⭐ (il più importante!)

#### Dallo Script 3:
- `RISULTATI_Sharing_VgenePairs.csv`
- `RISULTATI_Sharing_TRBV.csv`
- `RISULTATI_Sharing_TRAV.csv`
- `RISULTATI_Sharing_Clones_CDR3.csv` ⭐
- `RISULTATI_Diversity_Analysis.csv`
- `RISULTATI_Jaccard_Overlap.csv`

---

## 📊 INTERPRETAZIONE DEI RISULTATI

### Cosa Cercare nei Tuoi Dati:

#### 1. QUALITÀ DATI (Console Output Script 1)
```
Quality:
  Complete = 450    ✅ Ottimo! Usa questi
  TRB_only = 120    ⚠️ Verifica perché TRA manca
  TRA_only = 15     ⚠️ Insolito
  Incomplete = 5    ❌ Scarta questi
```

**Se hai >70% Complete:** Dati OK  
**Se hai <50% Complete:** Problema tecnico, contatta il sequencing core

#### 2. TOP CLONI (Plot Script 1)

**Plot "Metodo V-Gene":**
- Mostra espansioni di famiglie di cloni
- Utile per vedere quali V-genes dominano
- ⚠️ NON usare per contare cloni!

**Plot "Metodo CDR3" (⭐ USARE QUESTO!):**
- Ogni barra = UN clone vero
- Altezza barra = quante cellule di QUEL clone esatto
- Confronto tra stage mostra dinamica clonale

**Cosa dedurre:**
```
Clone espanso solo in Stage I -> Espansione precoce, poi contrazione
Clone espanso in Stage A+B -> Espansione persistente (buon segno!)
Clone presente in tutti gli stage con numeri simili -> Clone stabile
Clone che appare solo in B -> Espansione tardiva (investigate!)
```

#### 3. CONVERGENZA V-GENE (Script 2, Plot 1-2)

**Esempio dai tuoi dati:**
```
TRBV2 condiviso da Bo, Ca, Me in tutti gli stage
```

**Interpretazione:**
- ✅ Forte selezione per TRBV2
- ✅ Probabilmente riconosce antigene comune (CAR target o epitopo tumorale)
- ⚠️ PERÒ: Diversi cloni con TRBV2 possono riconoscere antigeni DIVERSI
- 👉 Vai a vedere i CDR3 di questi TRBV2 per capire se sono lo stesso clone

#### 4. PUBLIC CLONES VERI (Script 2, Plot 3 - ⭐⭐⭐)

**Se trovi cloni CDR3-identici condivisi tra pazienti:**
```
Es: Clone "TRAV21+TRBV2 / CASSLAPGTQYF"
    Presente in Bo (20 celle), Ca (15 celle), Me (8 celle)
```

**Questo è GOLD!** Significa:
- ✅ Clone IDENTICO (stesso recettore)
- ✅ Riconosce lo STESSO epitopo
- ✅ Probabilmente epitopo pubblico del CAR target
- ✅ Potenziale biomarcatore di risposta
- 💡 Candidato per TCR engineering o terapia off-the-shelf

**Se NON trovi cloni condivisi:**
- Normale! I veri public clones sono rari (~1-5% dei casi)
- Repertori privati sono la norma
- Convergenza V-gene è già significativa

#### 5. DIVERSITÀ CLONALE (Script 3, Plot 3-4)

**Shannon Diversity:**
```
Alto (>3):  Repertorio diversificato, molti cloni
Medio (1-3): Espansione oligoclonale bilanciata ✅ (ideale CAR-T)
Basso (<1):  Pochi cloni dominano
```

**Clonality Index:**
```
0.0-0.3: Troppo diversificato (possibile espansione insufficiente)
0.4-0.7: Range ottimale per CAR-T ✅
0.8-1.0: Monoespansione (possibile esaurimento a lungo termine)
```

**Pattern Dinamici:**
```
I -> A -> B

Shannon diminuisce: Selezione clonale (normale)
Shannon aumenta: Ri-diversificazione (possibile risposta secondaria)
Shannon stabile: Equilibrio raggiunto
```

#### 6. OVERLAP REPERTORI (Script 3, Plot 5 - Heatmap Jaccard)

**Valori Tipici:**
```
Stesso paziente, stesso stage:   1.00 (ovvio)
Stesso paziente, stage diversi:  0.20-0.60 (overlap parziale)
Pazienti diversi:                0.00-0.05 (repertori privati)
```

**Se Jaccard tra pazienti >0.10:**
- ⚠️ Alto! Hai public clones significativi
- Controlla `RISULTATI_Public_Clones_CDR3_identici.csv`
- Questi cloni meritano analisi funzionale approfondita

---

## 🎯 CONCETTI CHIAVE DA RICORDARE

### 1. GERARCHIA DELLA SPECIFICITÀ

```
V-gene (famiglia)
  ↓
J-gene (sottogruppo)
  ↓
CDR3 (clone specifico) ⭐ <- USA QUESTO!
```

**Esempio:**
- **TRBV2** = Famiglia con ~10,000 CDR3 possibili
- **TRBV2 + TRBJ1-1** = Sottogruppo con ~1,000 CDR3
- **TRBV2 + TRBJ1-1 + CASSLAPGTQYF** = Clone UNICO

### 2. CONVERGENZA vs PUBLIC CLONES

| Tipo | Definizione | Frequenza | Significato |
|------|-------------|-----------|-------------|
| **Convergenza V-gene** | Stesso V-gene, CDR3 diversi | Comune (30-50% V-genes) | Selezione verso famiglie geniche |
| **Convergenza CDR3** | Stesso CDR3, qualche aa diverso | Rara (1-5%) | Convergenza funzionale |
| **Public Clones** | CDR3 100% identico | Rarissima (<1%) | Stesso clone, forte selezione |

### 3. QUANDO QUALCOSA È "CONDIVISO"?

❌ **SBAGLIATO:** "Ho TRBV2 in entrambi i pazienti, quindi stesso clone"
✅ **CORRETTO:** "Ho TRBV2 in entrambi. Vediamo i CDR3..."

```
Paziente 1 TRBV2:
  - Clone A: CASSLAPGTQYF (100 celle)
  - Clone B: CASSLQGAYEQYF (50 celle)
  
Paziente 2 TRBV2:
  - Clone C: CASSLAPGTQYF (80 celle)  <- STESSO CLONE di A! ✅
  - Clone D: CASSLRTGAYEQYF (30 celle)
  
Public Clone = Solo "CASSLAPGTQYF" (A = C)
```

---

## 🔬 PROSSIMI PASSI SUGGERITI

### Analisi Immediate (Puoi fare con questi script):

1. **Identifica i Top Public Clones**
   ```r
   # Dal file RISULTATI_Public_Clones_CDR3_identici.csv
   public <- read.csv("RISULTATI_Public_Clones_CDR3_identici.csv")
   top_public <- public %>% 
     filter(n_patients == 3) %>%
     arrange(desc(total_cells))
   ```

2. **Confronta Diversità tra Stage**
   ```r
   # Dal file RISULTATI_Diversity_Analysis.csv
   div <- read.csv("RISULTATI_Diversity_Analysis.csv")
   
   # Quale stage ha maggiore clonality?
   div %>% arrange(desc(clonality))
   ```

3. **Trova i TRBV più espansi per paziente**
   ```r
   full_data %>%
     filter(Clone_Quality == "Complete") %>%
     group_by(patient, TRB_V) %>%
     summarise(n_cells = n()) %>%
     arrange(patient, desc(n_cells))
   ```

### Analisi Avanzate (Richiedono ulteriori tool):

1. **Analisi Motif CDR3**
   - Cerca pattern ricorrenti nei CDR3 (es: "LAP", "QYF")
   - Usa GLIPH2 o TCRdist per clustering

2. **Predizione Specificità**
   - Usa TCRex, NetTCR, o ERGO per predire antigene riconosciuto
   - Confronta con database VDJdb

3. **Network Analysis**
   - Costruisci network di similarità tra cloni (Hamming distance <3)
   - Identifica cluster di cloni correlati

4. **Analisi Temporale**
   - Se I < A < B è ordine temporale, plot traiettorie clonali
   - Identifica cloni espandenti vs contraenti

5. **Confronto con Controlli**
   - Se hai campioni pre-CAR, confronta repertori
   - Identifica cloni CAR-induced vs pre-esistenti

---

## ⚠️ DISCLAIMER E LIMITI

### Limiti Tecnici:

1. **5' vs 3' sequencing:**
   - Se hai usato 5' chemistry, dovresti avere anche il V(D)J
   - Se hai usato 3' chemistry + VDJ enrichment, tutto OK
   - Verifica quale hai usato

2. **Paired vs Unpaired:**
   - Gli script assumono che hai paired TRA+TRB
   - Se hai solo TRB, modifica per usare solo `Clone_ID = TRB_CDR3`

3. **Cellule Dual-TCR:**
   - Cellule con 2 TRA o 2 TRB sono rare ma esistono
   - Gli script le gestiscono con "/" separator
   - Per analisi rigorosa, potresti volerle escludere

4. **CAR vs TCR endogeno:**
   - Se il CAR ha sequenza TCR, potrebbe apparire nei dati VDJ
   - Filtra i CDR3 che matchano la sequenza CAR

### Limiti Biologici:

1. **Espansione non significa funzionalità:**
   - Un clone espanso potrebbe essere exhausted
   - Serve analisi fenotipica (CD39, PD-1, ecc.)

2. **Public clones rari sono normali:**
   - Non trovare public clones NON significa fallimento
   - Convergenza V-gene è già biologicamente rilevante

3. **Stage A definizione:**
   - Verifica che "A" sia davvero Aferesi e non altro
   - L'ordine temporale cambia l'interpretazione

---

## 📚 RISORSE UTILI

### Database TCR:
- **VDJdb**: Database di TCR con specificità nota
- **McPAS-TCR**: TCR associati a patologie
- **IEDB**: Epitopi immunogenici

### Tool di Analisi:
- **immunarch**: R package per analisi repertori
- **GLIPH2**: Clustering TCR per motif
- **TCRdist**: Calcola distanze tra TCR
- **Scirpy**: Python package per scTCR-seq

### Letture Consigliate:
- Dash et al. 2017 "Quantifiable predictive features define epitope-specific T cell receptor repertoires"
- Glanville et al. 2017 "Identifying specificity groups in the T cell receptor repertoire"
- Pogorelyy et al. 2019 "Detecting T cell receptors involved in immune responses"

---

## ✅ CHECKLIST FINALE

Prima di procedere con le conclusioni, verifica:

- [ ] Ho eseguito tutti e 3 gli script FIXED
- [ ] Ho almeno 50% cloni "Complete" nei dati
- [ ] Ho controllato la tabella `patient x stage` (tutti presenti?)
- [ ] Ho verificato che gli stage siano nell'ordine temporale corretto
- [ ] Ho guardato il file `RISULTATI_Public_Clones_CDR3_identici.csv`
- [ ] Ho capito la differenza tra V-gene convergence e public clones
- [ ] Ho calcolato Shannon diversity e clonality per ogni campione
- [ ] Ho interpretato il Jaccard index correttamente
- [ ] Ho identificato i top 3 cloni più espansi per paziente
- [ ] Ho verificato se i public clones (se presenti) hanno pattern biologicamente sensati

---

## 🎓 CONCLUSIONE

**Messaggi da Portare a Casa:**

1. **CDR3 è essenziale** - Non fare mai analisi TCR senza CDR3
2. **V-gene ≠ Clone** - Convergenza V-gene è interessante ma non è clonalità
3. **Public clones sono rari** - Non aspettarti di trovarne molti (se ci sono, è GOLD!)
4. **Qualità > Quantità** - Meglio 100 cloni completi che 1000 parziali
5. **Context matters** - Shannon, Clonality, Jaccard vanno interpretati insieme

**Prossimo Step:**
Esegui gli script sui tuoi dati e guarda i risultati. Se hai domande specifiche sui pattern che trovi, sono qui per aiutarti a interpretarli!

---

*Documento generato insieme agli script corretti*  
*Versione 1.0 - Data: [Data esecuzione]*
