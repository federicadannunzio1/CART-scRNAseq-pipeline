# ==============================================================================
# 15_relative_clonality.R
#
# Analisi: clonalità relativa ed espansione per coppia TRAV+TRBV macrofamiglia
#
# Unità di analisi principale: coppia TRAV+TRBV macrofamiglia
#   (aggiuntivo: metriche a livello clone CDR3 per confronto)
#
# METRICHE USATE:
#
#   Clonalità V-gene (Simpson's index per coppie V-gene):
#     C_vgene = Σ p²ᵢ
#     dove pᵢ = freq relativa della coppia V-gene i nel paziente × stage
#     Rango: 0 (tutte le coppie V-gene ugualmente usate) → 1 (una coppia domina)
#
#   Espansione relativa di una coppia V-gene (fold-change):
#     FC = (n_cellule_pair_B / n_cellule_tot_B) / (n_cellule_pair_I / n_cellule_tot_I)
#        = freq_B / freq_I
#     Se freq_I = 0 → coppia "De novo in B" (non rilevata nel prodotto di infusione)
#     Log2(FC) usato per visualizzazione
#
#   Dominanza V-gene: freq relativa della coppia V-gene più abbondante
#
#   (Aggiuntivo) Clonalità CDR3 normalizzata per confronto:
#     C_cdr3 = 1 − H/log(N)  dove H=Shannon entropy, N=N cloni unici
#
# Dipende da: RISULTATI_expansion_dynamics.xlsx
#             final_clone_sequences.xlsx
#
# Output: Fig15a_vgene_clonality_stages.png
#         Fig15b_vgene_expansion_fc.png
#         Fig15c_vgene_rank_frequency.png
#         Fig15d_cdr3_clonality_comparison.png
#         15_relative_clonality.xlsx
# ==============================================================================

suppressMessages({
  library(dplyr); library(tidyr); library(ggplot2)
  library(readxl); library(writexl); library(stringr); library(scales); library(patchwork)
})

TAB <- "/Users/federicadannunzio/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/caruana-project/CART/Code/07_clonotypes/results/tables"
FIG <- "/Users/federicadannunzio/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/caruana-project/CART/Code/07_clonotypes/results/figures"

PAT_COL   <- c(Bo = "#E64B35", Ca = "#4DBBD5", Me = "#00A087")
PAT_LABEL <- c(Bo = "Bo (expansion)", Ca = "Ca (failure)", Me = "Me (partial)")

macro <- function(x) str_remove(x, "-[0-9]+$")

# ── STEP 1: Carica dati ─────────────────────────────────────────────────────
message("\n--- STEP 1: Caricamento dati ---")

fin <- read_xlsx(file.path(TAB, "final_clone_sequences.xlsx")) %>%
  mutate(macro_TRA  = macro(TRA_v_gene),
         macro_TRB  = macro(TRB_v_gene),
         pair_vgene = paste0(macro_TRA, " + ", macro_TRB))

clone_wide <- read_xlsx(file.path(TAB, "RISULTATI_expansion_dynamics.xlsx"),
                        sheet = "01_Dinamica_completa")

message(sprintf("  final_clone_sequences: %d record", nrow(fin)))
message(sprintf("  expansion_dynamics:    %d record", nrow(clone_wide)))

# ── STEP 2: Frequenza coppia V-gene per paziente × stage ─────────────────────
message("\n--- STEP 2: Frequenza coppia V-gene per paziente × stage ---")

vgene_counts <- fin %>%
  group_by(patient, stage, pair_vgene) %>%
  summarise(n_cells       = sum(n_cells),
            n_cloni_cdr3  = n_distinct(paste(TRA_cdr3, TRB_cdr3)),
            .groups="drop")

# Frequenza relativa per paziente × stage
vgene_freq <- vgene_counts %>%
  group_by(patient, stage) %>%
  mutate(n_tot_stage = sum(n_cells),
         freq        = n_cells / n_tot_stage) %>%
  ungroup()

# Completa con 0 per stage mancanti (no Me-A)
vgene_freq_complete <- vgene_freq %>%
  complete(nesting(patient, pair_vgene), stage=c("I","A","B"),
           fill=list(n_cells=0, freq=0, n_cloni_cdr3=0)) %>%
  filter(!(patient=="Me" & stage=="A"))

# Pivot wide per calcolo FC
vgene_wide <- vgene_freq_complete %>%
  select(patient, pair_vgene, stage, n_cells, freq) %>%
  pivot_wider(names_from=stage, values_from=c(n_cells,freq), values_fill=0) %>%
  mutate(
    FC_I_to_B  = case_when(
      freq_I==0 & freq_B==0 ~ NA_real_,
      freq_I==0 & freq_B>0  ~ NA_real_,   # non osservato in I = artefatto campionamento
      TRUE                  ~ freq_B / freq_I
    ),
    FC_I_to_A  = case_when(
      freq_I==0 & freq_A==0 ~ NA_real_,
      freq_I==0 & freq_A>0  ~ NA_real_,
      TRUE                  ~ freq_A / freq_I
    ),
    log2FC_I_to_B = case_when(
      is.na(FC_I_to_B) ~ NA_real_,
      TRUE             ~ log2(FC_I_to_B)
    ),
    # "Non rilevato in I" = la coppia V-gene non è stata catturata nel campionamento di I,
    # ma TRAV e TRBV singoli sono quasi sempre già presenti in I → artefatto di profondità
    non_rilevato_I = freq_I == 0 & freq_B > 0,
    categoria = case_when(
      freq_I==0 & freq_B>0  ~ "Non rilevato in I",
      is.na(FC_I_to_B)      ~ "Assente",
      FC_I_to_B >= 2        ~ "Espanso (FC>=2)",
      FC_I_to_B >= 1        ~ "Stabile",
      FC_I_to_B <  1        ~ "Contratto",
      TRUE                  ~ "Altro"
    )
  )

message("  Categorie di espansione V-gene per paziente:")
print(vgene_wide %>% filter(!is.na(categoria), categoria!="Assente") %>%
        count(patient, categoria) %>% arrange(patient, categoria))

# ── STEP 3: Indici di clonalità V-gene per paziente × stage ──────────────────
message("\n--- STEP 3: Indici di clonalità (V-gene pairs) ---")

calc_div <- function(freqs) {
  p <- freqs[freqs > 0]
  if (length(p)==0) return(list(n=0, simpson=NA, clonality_norm=NA, dominance=NA, shannon=NA))
  H   <- -sum(p * log(p + 1e-12))
  Hm  <- log(length(p))
  list(n             = length(p),
       simpson       = round(sum(p^2), 4),
       clonality_norm= round(if(Hm>0) 1-H/Hm else NA, 4),
       dominance     = round(max(p), 4),
       shannon       = round(H, 4))
}

div_vgene <- bind_rows(lapply(c("I","A","B"), function(st) {
  col_f <- paste0("freq_",st); col_n <- paste0("n_cells_",st)
  if (!col_f %in% colnames(vgene_wide)) return(NULL)
  lapply(unique(vgene_wide$patient), function(pt) {
    sub <- vgene_wide %>% filter(patient==pt, .data[[col_n]]>0) %>% pull(col_f)
    if (length(sub)==0 || all(sub==0)) return(NULL)
    dm <- calc_div(sub)
    data.frame(patient=pt, stage=st, n_pairs=dm$n, simpson=dm$simpson,
               clonality_norm=dm$clonality_norm, dominance=dm$dominance,
               shannon=dm$shannon)
  }) %>% bind_rows()
})) %>%
  mutate(stage=factor(stage, levels=c("I","A","B")),
         patient=factor(patient, levels=c("Bo","Ca","Me")))

message("  Indici clonalità V-gene:")
print(div_vgene %>% select(patient, stage, n_pairs, clonality_norm, simpson, dominance))

# Indici CDR3 (per confronto)
div_cdr3 <- bind_rows(lapply(c("I","A","B"), function(st) {
  col_f <- paste0("freq_",st); col_n <- paste0("n_cells_",st)
  if (!col_f %in% colnames(clone_wide)) return(NULL)
  lapply(unique(clone_wide$patient), function(pt) {
    sub <- clone_wide %>% filter(patient==pt, .data[[col_n]]>0) %>% pull(col_f)
    if (length(sub)==0) return(NULL)
    dm <- calc_div(sub[!is.na(sub)])
    data.frame(patient=pt, stage=st, n_clones=dm$n, simpson_cdr3=dm$simpson,
               clonality_norm_cdr3=dm$clonality_norm)
  }) %>% bind_rows()
})) %>% mutate(stage=factor(stage, levels=c("I","A","B")),
               patient=factor(patient, levels=c("Bo","Ca","Me")))

# ── STEP 4: Figure ──────────────────────────────────────────────────────────
message("\n--- STEP 4: Figure ---")

# Figura A — clonalità V-gene attraverso gli stage
p_simp_vg <- ggplot(div_vgene,
                    aes(x=stage, y=simpson, color=patient, group=patient)) +
  geom_line(linewidth=1.5) + geom_point(size=4) +
  geom_text(aes(label=round(simpson,3)), vjust=-0.8, size=3.5, fontface="bold") +
  scale_color_manual(values=PAT_COL, labels=PAT_LABEL, name=NULL) +
  scale_y_continuous(limits=c(0,NA)) +
  theme_minimal(base_size=12) +
  theme(legend.position="bottom", panel.grid.minor=element_blank()) +
  labs(title="V-gene pair clonality (Simpson's Σp²) across stages",
       subtitle="p = relative frequency of each TRAV+TRBV pair | 0=diverse → 1=monopoly",
       x="Stage", y="Simpson's clonality (Σp²)")

p_dom_vg <- ggplot(div_vgene,
                   aes(x=stage, y=dominance*100, color=patient, group=patient)) +
  geom_line(linewidth=1.5) + geom_point(size=4) +
  geom_text(aes(label=paste0(round(dominance*100,1),"%")), vjust=-0.8, size=3.5, fontface="bold") +
  scale_color_manual(values=PAT_COL, labels=PAT_LABEL, name=NULL) +
  scale_y_continuous(labels=function(x) paste0(x,"%"), expand=expansion(mult=c(0,0.15))) +
  theme_minimal(base_size=12) +
  theme(legend.position="bottom", panel.grid.minor=element_blank()) +
  labs(title="Dominant V-gene pair frequency across stages",
       subtitle="% cells using the most frequent TRAV+TRBV pair",
       x="Stage", y="Top V-gene pair frequency (%)")

p_npairs <- ggplot(div_vgene,
                   aes(x=stage, y=n_pairs, color=patient, group=patient)) +
  geom_line(linewidth=1.5) + geom_point(size=4) +
  geom_text(aes(label=n_pairs), vjust=-0.8, size=3.5, fontface="bold") +
  scale_color_manual(values=PAT_COL, labels=PAT_LABEL, name=NULL) +
  theme_minimal(base_size=12) +
  theme(legend.position="bottom", panel.grid.minor=element_blank()) +
  labs(title="Number of unique TRAV+TRBV pairs per stage",
       x="Stage", y="N unique V-gene pairs")

fig15a <- (p_simp_vg | p_dom_vg | p_npairs) +
  plot_layout(guides="collect") & theme(legend.position="bottom")
fig15a <- fig15a + plot_annotation(
  title    = "Figure 15A — V-gene pair clonality across stages I → A → B",
  tag_levels = "a"
)

ggsave(file.path(FIG, "Fig15a_vgene_clonality_stages.png"),
       fig15a, width=15, height=6, dpi=300, bg="white")
message("Salvata: Fig15a_vgene_clonality_stages.png")

# Figura B — FC espansione V-gene (log2FC, top coppie)
top_pairs_B <- vgene_wide %>%
  filter(n_cells_B > 0) %>%
  group_by(patient) %>%
  slice_max(n_cells_B, n=10, with_ties=FALSE) %>%
  ungroup()

p_fc_bar <- ggplot(top_pairs_B %>% filter(!non_rilevato_I, !is.na(log2FC_I_to_B)),
                   aes(x=log2FC_I_to_B, y=reorder(pair_vgene, log2FC_I_to_B),
                       fill=patient)) +
  geom_col(width=0.7, color="white") +
  geom_vline(xintercept=1, linetype="dashed", color="grey40") +
  geom_vline(xintercept=0, linetype="solid",  color="grey70") +
  facet_wrap(~patient, scales="free_y", ncol=1,
             labeller=as_labeller(PAT_LABEL)) +
  scale_fill_manual(values=PAT_COL, guide="none") +
  theme_minimal(base_size=11) +
  theme(strip.text=element_text(face="bold"),
        panel.grid.major.y=element_blank()) +
  labs(title="V-gene pair expansion: log₂(FC I→B) — top 10 in B per patient",
       subtitle="FC = freq_B / freq_I | Solo coppie osservate in I | dashed = FC=2 | solid = FC=1",
       x="log₂(FC I→B)", y="TRAV + TRBV macrofamily pair")

ggsave(file.path(FIG, "Fig15b_vgene_expansion_fc.png"),
       p_fc_bar, width=10, height=12, dpi=300, bg="white")
message("Salvata: Fig15b_vgene_expansion_fc.png")

# Figura C — rank-frequency V-gene pair per stage
rank_vg <- vgene_freq %>%
  filter(freq>0) %>%
  group_by(patient, stage) %>%
  arrange(desc(freq)) %>%
  mutate(rank=row_number(),
         stage=factor(stage, levels=c("I","A","B"))) %>%
  ungroup()

p_rank <- ggplot(rank_vg, aes(x=rank, y=freq, color=stage, group=stage)) +
  geom_line(linewidth=0.9, alpha=0.85) +
  scale_color_manual(values=c("I"="#619CFF","A"="#F8766D","B"="#00BA38"), name="Stage") +
  scale_y_log10(labels=percent_format(accuracy=0.1)) +
  facet_wrap(~patient, scales="free_x", labeller=as_labeller(PAT_LABEL)) +
  theme_minimal(base_size=12) +
  theme(strip.text=element_text(face="bold"), legend.position="bottom") +
  labs(title="V-gene pair rank-frequency distribution per stage",
       subtitle="Steeper slope in B = repertoire focused on fewer V-gene pairs",
       x="V-gene pair rank (1 = most used)", y="Relative frequency")

ggsave(file.path(FIG, "Fig15c_vgene_rank_frequency.png"),
       p_rank, width=12, height=5, dpi=300, bg="white")
message("Salvata: Fig15c_vgene_rank_frequency.png")

# Figura D — confronto clonalità CDR3 vs V-gene
div_merged <- div_vgene %>%
  select(patient, stage, simpson_vgene=simpson) %>%
  left_join(div_cdr3 %>% select(patient, stage, simpson_cdr3), by=c("patient","stage")) %>%
  pivot_longer(cols=c(simpson_vgene, simpson_cdr3),
               names_to="livello", values_to="simpson") %>%
  mutate(livello=recode(livello,
                        "simpson_vgene"="V-gene pair",
                        "simpson_cdr3"="CDR3 (clone)"))

p_compare <- ggplot(div_merged %>% filter(!is.na(simpson)),
                    aes(x=stage, y=simpson, color=patient,
                        group=interaction(patient, livello),
                        linetype=livello)) +
  geom_line(linewidth=1.2) + geom_point(size=3) +
  scale_color_manual(values=PAT_COL, labels=PAT_LABEL, name="Patient") +
  scale_linetype_manual(values=c("V-gene pair"="solid","CDR3 (clone)"="dashed"),
                        name="Level") +
  theme_minimal(base_size=12) +
  theme(legend.position="bottom") +
  labs(title="Clonality comparison: V-gene pair level vs CDR3 clone level",
       subtitle="Solid = V-gene pair Simpson | Dashed = CDR3 clone Simpson",
       x="Stage", y="Simpson's clonality (Σp²)")

ggsave(file.path(FIG, "Fig15d_clonality_comparison.png"),
       p_compare, width=10, height=6, dpi=300, bg="white")
message("Salvata: Fig15d_clonality_comparison.png")

# ── STEP 5: Coppie dominanti per stage ──────────────────────────────────────
dominant_pairs <- bind_rows(lapply(c("I","A","B"), function(st) {
  col_f <- paste0("freq_",st); col_n <- paste0("n_cells_",st)
  if (!col_f %in% colnames(vgene_wide)) return(NULL)
  vgene_wide %>%
    filter(.data[[col_n]]>0) %>%
    group_by(patient) %>%
    slice_max(.data[[col_f]], n=3, with_ties=FALSE) %>%
    mutate(stage=st, rank=row_number(), freq=.data[[col_f]], n_cells=.data[[col_n]]) %>%
    ungroup() %>%
    select(patient, stage, rank, pair_vgene, n_cells, freq)
}))

# ── STEP 6: Salva tabelle ───────────────────────────────────────────────────
write_xlsx(list(
  "01_Clonalita_Vgene"          = div_vgene,
  "02_Clonalita_CDR3"           = div_cdr3,
  "03_FC_Vgene_pairs"           = vgene_wide %>% filter(!is.na(categoria), categoria!="Assente"),
  "04_Coppie_dominanti_per_stage"= dominant_pairs,
  "05_Frequenze_Vgene"          = vgene_freq
), file.path(TAB, "15_relative_clonality.xlsx"))

message("Salvata: 15_relative_clonality.xlsx")
message("\nRiepilogo:")
message("  Livello di analisi: coppie TRAV+TRBV macrofamiglia")
message("  Metrica clonalità V-gene: Simpson's = Σ(freq_pair²)")
message("  Metrica espansione: FC = freq_B_pair / freq_I_pair")
message("")
message("  Coppie V-gene non rilevate in I (artefatto campionamento, FC non calcolabile):")
print(vgene_wide %>% filter(non_rilevato_I) %>% count(patient))
message("  Nota: TRAV e TRBV di queste coppie sono quasi sempre presenti in I in altre combinazioni.")
message("  V-gene pairs con FC>=2 in B (osservate in I, genuinamente espanse):")
print(vgene_wide %>% filter(categoria=="Espanso (FC>=2)") %>% count(patient))
