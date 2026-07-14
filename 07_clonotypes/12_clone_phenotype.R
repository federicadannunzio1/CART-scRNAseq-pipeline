# ==============================================================================
# 12_clone_phenotype.R
#
# Analisi: coppie TRAV+TRBV macrofamiglia → fenotipo/stato funzionale
#
# Unità di analisi: coppia di macrofamiglie V-gene
#   macro(TRAV12-1) = TRAV12,  macro(TRBV7-9) = TRBV7
#   pair = "TRAV12 + TRBV7"
#
# Risponde a: "Le cellule che usano le stesse coppie TRAV+TRBV degli espansi di Bo
# — che fenotipo/stato funzionale hanno in ciascun paziente e stage?"
#
# Dipende da:
#   RISULTATI_Cloni_Dati_Completi_con_CDR3.xlsx  (cellule individuali + barcodes)
#   RISULTATI_expansion_dynamics.xlsx            (cloni espansi Bo)
#   all_I_samples_annotated.rds                  (cell_type stage I)
#   all_AB_samples_annotated.rds                 (cell_type stage A/B)
#
# Output: Fig12a_vgene_pair_phenotype_bar.png
#         Fig12b_vgene_pair_phenotype_heatmap.png
#         12_clone_phenotype.xlsx
# ==============================================================================

suppressMessages({
  library(Seurat)
  library(dplyr); library(tidyr); library(ggplot2)
  library(readxl); library(writexl); library(stringr)
  library(scales); library(patchwork)
})

TAB <- "/Users/federicadannunzio/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/caruana-project/CART/Code/07_clonotypes/results/tables"
FIG <- "/Users/federicadannunzio/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/caruana-project/CART/Code/07_clonotypes/results/figures"

RDS_I  <- "/Users/federicadannunzio/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/caruana-project/CART/Code/01_seurat_annotation/results/all_I_samples_annotated.rds"
RDS_AB <- "/Users/federicadannunzio/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/caruana-project/CART/Code/01_seurat_annotation/results/all_AB_samples_annotated.rds"

# ── Costanti ────────────────────────────────────────────────────────────────
PAT_COL   <- c(Bo = "#E64B35", Ca = "#4DBBD5", Me = "#00A087")
PAT_LABEL <- c(Bo = "Bo (expansion)", Ca = "Ca (failure)", Me = "Me (partial)")

FUNCTIONAL_STATE_MAP <- list(
  "Naive-like"  = c("Naive CD4+ T cells", "Naive CD8+ T cells"),
  "Memory-like" = c("Memory T cells", "Th1 cells", "Th2 cells", "Th17 cells", "Tfh cells"),
  "Effector"    = c("Effector CD4+ T cells", "Cytotoxic CD8+ T cells"),
  "Regulatory"  = c("Tregs")
)
FUNCTIONAL_ORDER <- c("Naive-like", "Memory-like", "Effector", "Regulatory", "Other")
STATE_PALETTE    <- c("Naive-like"  = "#4DBBD5", "Memory-like" = "#00A087",
                      "Effector"    = "#E64B35", "Regulatory"  = "#F39B7F",
                      "Other"       = "grey80")

# Macrofamiglia: rimuove l'allele (es. TRAV12-1 → TRAV12)
macro <- function(x) str_remove(x, "-[0-9]+$")

map_fs <- function(cell_types) {
  state <- rep("Other", length(cell_types))
  for (nm in names(FUNCTIONAL_STATE_MAP))
    state[cell_types %in% FUNCTIONAL_STATE_MAP[[nm]]] <- nm
  state
}

get_car_status <- function(meta) {
  for (col in c("IS_CAR_ALLIN_scREP", "IS_CAR", "CAR")) {
    if (col %in% colnames(meta)) {
      vals <- as.character(meta[[col]])
      return(ifelse(grepl("^(YES|TRUE|yes|true|1)$", vals), "CAR+", "CAR-"))
    }
  }
  rep("CAR-", nrow(meta))
}

# ── STEP 1: Estrai metadati Seurat ───────────────────────────────────────────
message("\n--- STEP 1: Caricamento metadati Seurat ---")

extract_meta <- function(seurat_list_obj) {
  bind_rows(lapply(names(seurat_list_obj), function(sname) {
    obj  <- seurat_list_obj[[sname]]
    meta <- obj@meta.data
    bc   <- rownames(meta)
    bc_for_join <- gsub("-[0-9]+$", "", sub("^([^_]+_){2}", "", bc))
    folder  <- str_extract(bc, "^[^_]+_[^_]+")
    stage_bc <- case_when(
      grepl("_I", folder) ~ "I",
      grepl("_A", folder) ~ "A",
      grepl("_B", folder) ~ "B",
      TRUE                ~ NA_character_
    )
    pat <- if ("patient" %in% colnames(meta)) as.character(meta$patient) else
      rep(str_extract(sname, "^[A-Z][a-z]+"), nrow(meta))
    data.frame(bc_for_join = bc_for_join, patient = pat, stage = stage_bc,
               cell_type   = if ("cell_type" %in% colnames(meta)) as.character(meta$cell_type) else NA_character_,
               car_status  = get_car_status(meta), stringsAsFactors = FALSE)
  }))
}

message("  Caricamento all_I_samples_annotated.rds ...")
I_samples  <- readRDS(RDS_I)
message("  Caricamento all_AB_samples_annotated.rds ...")
AB_samples <- readRDS(RDS_AB)

seurat_meta <- bind_rows(extract_meta(I_samples), extract_meta(AB_samples)) %>%
  mutate(functional_state = map_fs(cell_type),
         functional_state = factor(functional_state, levels = FUNCTIONAL_ORDER))

message(sprintf("  Cellule totali nel meta: %d  |  CAR+: %d",
                nrow(seurat_meta), sum(seurat_meta$car_status == "CAR+", na.rm=TRUE)))

# ── STEP 2: Carica clonotipi + calcola coppia V-gene ────────────────────────
message("\n--- STEP 2: Caricamento dati clonotipi ---")

full_data_raw <- read_xlsx(file.path(TAB, "RISULTATI_Cloni_Dati_Completi_con_CDR3.xlsx"))

exclude_nt <- full_data_raw %>%
  filter(Clone_Quality == "Complete",
         !is.na(TRA_cdr3_nt), !is.na(TRB_cdr3_nt),
         TRA_cdr3_nt != "", TRB_cdr3_nt != "") %>%
  group_by(TRA_cdr3_nt, TRB_cdr3_nt, patient) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(TRA_cdr3_nt, TRB_cdr3_nt) %>%
  filter(n_distinct(patient) > 1) %>%
  ungroup() %>%
  distinct(TRA_cdr3_nt, TRB_cdr3_nt)

clean_full <- full_data_raw %>%
  filter(Clone_Quality == "Complete",
         !is.na(TRA_cdr3_nt), !is.na(TRB_cdr3_nt)) %>%
  anti_join(exclude_nt, by = c("TRA_cdr3_nt", "TRB_cdr3_nt")) %>%
  mutate(
    macro_TRA  = macro(TRA_v_gene),
    macro_TRB  = macro(TRB_v_gene),
    pair_vgene = paste0(macro_TRA, " + ", macro_TRB),
    bc_for_join = gsub("-[0-9]+$", "", clean_barcode)
  )

message(sprintf("  Cellule post-decontaminazione: %d", nrow(clean_full)))
message("  Coppie V-gene uniche per paziente:")
print(clean_full %>% count(patient, pair_vgene) %>% count(patient, name="n_pairs"))

# ── STEP 3: Coppie V-gene dei cloni espansi in Bo ───────────────────────────
message("\n--- STEP 3: Coppie V-gene dei cloni espansi in Bo ---")

espansi <- read_xlsx(file.path(TAB, "RISULTATI_expansion_dynamics.xlsx"),
                     sheet = "02_Cloni_espansi_in_B")

bo_exp_pairs <- espansi %>%
  filter(patient == "Bo") %>%
  mutate(macro_TRA  = macro(TRA_v_gene),
         macro_TRB  = macro(TRB_v_gene),
         pair_vgene = paste0(macro_TRA, " + ", macro_TRB)) %>%
  group_by(pair_vgene) %>%
  summarise(n_cloni_espansi = n(),
            n_cells_Bo_B    = sum(n_cells_B), .groups="drop") %>%
  arrange(desc(n_cells_Bo_B))

message("  Coppie V-gene degli espansi Bo: ", nrow(bo_exp_pairs))

# ── STEP 4: Join VDJ ↔ cell_type ────────────────────────────────────────────
message("\n--- STEP 4: Join barcode → cell_type ---")

clono_pheno <- clean_full %>%
  left_join(seurat_meta %>%
              select(patient, bc_for_join, cell_type, functional_state, car_status),
            by = c("patient", "bc_for_join"))

n_matched <- sum(!is.na(clono_pheno$cell_type))
message(sprintf("  Celle matchate: %d / %d (%.1f%%)",
                n_matched, nrow(clono_pheno), 100*n_matched/nrow(clono_pheno)))

if (n_matched == 0)
  stop("Nessuna cellula matchata. Controlla il formato dei barcodes.")

# ── STEP 5: Fenotipo per coppia V-gene × paziente × stage ───────────────────
message("\n--- STEP 5: Aggregazione fenotipo per coppia V-gene ---")

# Filtra alle coppie V-gene degli espansi Bo (presenti in tutti i pazienti)
pheno_vgene <- clono_pheno %>%
  filter(pair_vgene %in% bo_exp_pairs$pair_vgene,
         !is.na(functional_state)) %>%
  group_by(patient, stage, pair_vgene, functional_state) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  group_by(patient, stage, pair_vgene) %>%
  mutate(prop = n_cells / sum(n_cells),
         n_tot = sum(n_cells)) %>%
  ungroup() %>%
  mutate(
    stage   = factor(stage, levels = c("I","A","B")),
    patient = factor(patient, levels = c("Bo","Ca","Me"))
  ) %>%
  left_join(bo_exp_pairs %>% select(pair_vgene, n_cells_Bo_B), by="pair_vgene") %>%
  mutate(pair_vgene = factor(pair_vgene,
                              levels = bo_exp_pairs$pair_vgene))

# Top 10 coppie per n_cells_Bo_B
top10_pairs <- bo_exp_pairs %>% slice_head(n=10) %>% pull(pair_vgene)

# ── STEP 6: Figura A — heatmap fenotipo per coppia V-gene × paziente ─────────
message("\n--- STEP 6: Figure ---")

heat_data <- pheno_vgene %>%
  filter(pair_vgene %in% top10_pairs) %>%
  mutate(pat_stage = paste0(patient, "\n", stage),
         pat_stage = factor(pat_stage,
                             levels = c("Bo\nI","Bo\nA","Bo\nB",
                                        "Ca\nI","Ca\nA","Me\nI","Me\nB"))) %>%
  filter(!is.na(pat_stage))

p_heat <- ggplot(heat_data,
                 aes(x = pat_stage, y = pair_vgene, fill = prop)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(n_tot >= 3, percent(prop, accuracy=1), "")),
            size = 2.8, color = "grey20") +
  facet_wrap(~ functional_state, nrow = 1) +
  scale_fill_gradientn(
    colors = c("white","#FEE8C8","#FC8D59","#D73027"),
    name   = "Proportion", labels = percent_format()
  ) +
  geom_vline(xintercept = c(3.5, 5.5), color="grey50", linetype="dashed", linewidth=0.4) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x   = element_text(size=8),
        axis.text.y   = element_text(size=9),
        panel.grid    = element_blank(),
        strip.text    = element_text(face="bold"),
        legend.position = "right") +
  labs(
    title    = "Functional state of cells by TRAV+TRBV pair across patients and stages",
    subtitle = paste0(
      "Top 10 V-gene pairs of Bo's expanded clones — all patients shown\n",
      "Cells using the same V-gene pair as Bo's expanded clones, split by functional state"
    ),
    x = NULL, y = "TRAV + TRBV macrofamily pair"
  )

ggsave(file.path(FIG, "Fig12a_vgene_pair_phenotype_heatmap.png"),
       p_heat, width=16, height=8, dpi=300, bg="white")
message("Salvata: Fig12a_vgene_pair_phenotype_heatmap.png")

# Figura B — stacked bar: distribuzione stato funzionale per coppia V-gene × paziente
bar_data <- pheno_vgene %>%
  filter(pair_vgene %in% top10_pairs) %>%
  group_by(patient, pair_vgene, functional_state) %>%
  summarise(n_cells = sum(n_cells), .groups="drop") %>%
  group_by(patient, pair_vgene) %>%
  mutate(prop = n_cells / sum(n_cells)) %>%
  ungroup() %>%
  mutate(patient = factor(patient, levels=c("Bo","Ca","Me")))

p_bar <- ggplot(bar_data,
                aes(x = patient, y = prop, fill = functional_state)) +
  geom_col(width=0.7, color="white", linewidth=0.3) +
  geom_text(aes(label=ifelse(prop>=0.1, percent(prop,accuracy=1),"")),
            position=position_stack(vjust=0.5),
            size=2.8, color="white", fontface="bold") +
  facet_wrap(~ pair_vgene, ncol=5) +
  scale_fill_manual(values=STATE_PALETTE, name="Functional state", drop=FALSE) +
  scale_x_discrete(labels=PAT_LABEL) +
  scale_y_continuous(labels=percent_format(), expand=c(0,0)) +
  theme_minimal(base_size=10) +
  theme(strip.text = element_text(face="bold", size=8),
        axis.text.x = element_text(angle=20, hjust=1, size=8),
        panel.grid.major.x = element_blank(),
        legend.position = "bottom") +
  labs(
    title    = "Functional state by TRAV+TRBV pair — all stages combined",
    subtitle = "Top 10 V-gene pairs of Bo's expanded clones | All 3 patients",
    x=NULL, y="Proportion"
  )

ggsave(file.path(FIG, "Fig12b_vgene_pair_phenotype_bar.png"),
       p_bar, width=14, height=8, dpi=300, bg="white")
message("Salvata: Fig12b_vgene_pair_phenotype_bar.png")

# ── STEP 7: Salva tabelle ────────────────────────────────────────────────────
write_xlsx(list(
  "Fenotipo_per_pair_paz_stage" = pheno_vgene,
  "Fenotipo_aggregato_per_pair"  = bar_data,
  "Coppie_espansi_Bo"            = bo_exp_pairs
), file.path(TAB, "12_clone_phenotype.xlsx"))

message("Salvata: 12_clone_phenotype.xlsx")
message("\nRiepilogo:")
message("  Coppie V-gene degli espansi Bo: ", nrow(bo_exp_pairs))
message("  Di cui con celle matched in altri pazienti: ",
        n_distinct(pheno_vgene$pair_vgene))
