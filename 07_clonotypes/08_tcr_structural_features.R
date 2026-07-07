# ==============================================================================
# 08_tcr_structural_features.R
#
# Domanda: quali caratteristiche strutturali del TCR distinguono i cloni
# che si espandono da quelli che non si espandono?
#
# Confronto:
#   expanded     = 27 cloni di Bo con >=5 cellule in stage B
#   non_exp_Bo   = cloni di Bo in stage B con <5 cellule
#   non_exp_CaMe = tutti i cloni di Ca (stage I) e Me (stage I+B)
#
# Analisi:
#   1. Lunghezza CDR3 alpha e beta
#   2. V-gene usage (alpha e beta)
#   3. J-gene usage (alpha e beta)
#   4. D-gene usage (beta)
#   5. Composizione aminoacidica CDR3 beta (posizione per posizione)
#   6. Idrofobicità media CDR3 beta (scala Kyte-Doolittle)
#   7. Carica netta CDR3 beta
# ==============================================================================

suppressMessages({
  library(dplyr); library(tidyr); library(ggplot2)
  library(readxl); library(writexl); library(stringr)
  library(patchwork); library(forcats)
})

BASE <- "/Users/federicadannunzio/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/caruana-project/CART/Code/07_clonotypes"
TAB  <- file.path(BASE, "results", "tables")
FIG  <- file.path(BASE, "results", "figures")

# ── Carica dati ───────────────────────────────────────────────────────────────
fin <- read_xlsx(file.path(TAB, "final_clone_sequences.xlsx"))
esp <- read_xlsx(file.path(TAB, "RISULTATI_expansion_dynamics.xlsx"),
                 sheet = "02_Cloni_espansi_in_B")

bo_expanded_ids <- esp %>%
  filter(patient == "Bo") %>%
  select(TRA_cdr3, TRB_cdr3)

# ── Definisci gruppi ──────────────────────────────────────────────────────────
clones <- fin %>%
  mutate(
    group = case_when(
      patient == "Bo" & stage == "B" &
        paste(TRA_cdr3, TRB_cdr3) %in%
        paste(bo_expanded_ids$TRA_cdr3, bo_expanded_ids$TRB_cdr3) ~ "Expanded (Bo≥5)",
      patient == "Bo" & stage == "B"                               ~ "Non-exp Bo (B<5)",
      patient == "Ca"                                              ~ "Non-exp Ca",
      patient == "Me"                                              ~ "Non-exp Me",
      TRUE                                                         ~ NA_character_
    )
  ) %>%
  filter(!is.na(group)) %>%
  mutate(
    len_alpha = nchar(TRA_cdr3),
    len_beta  = nchar(TRB_cdr3),
    macro_TRA_V = str_remove(TRA_v_gene, "-[0-9]+$"),
    macro_TRB_V = str_remove(TRB_v_gene, "-[0-9]+$"),
    macro_TRA_J = str_remove(TRA_j_gene, "-[0-9]+$"),
    macro_TRB_J = str_remove(TRB_j_gene, "-[0-9]+$")
  )

group_colors <- c(
  "Expanded (Bo≥5)"   = "#E64B35",
  "Non-exp Bo (B<5)"  = "#F39B7F",
  "Non-exp Ca"        = "#4DBBD5",
  "Non-exp Me"        = "#00A087"
)

message("Cloni per gruppo:")
print(table(clones$group))

# ── 1. LUNGHEZZA CDR3 ─────────────────────────────────────────────────────────
message("\n=== 1. Lunghezza CDR3 ===")

len_stats <- clones %>%
  group_by(group) %>%
  summarise(
    median_alpha = median(len_alpha, na.rm=TRUE),
    median_beta  = median(len_beta,  na.rm=TRUE),
    mean_alpha   = round(mean(len_alpha, na.rm=TRUE),1),
    mean_beta    = round(mean(len_beta,  na.rm=TRUE),1),
    .groups="drop"
  )
print(len_stats)

p_len_a <- ggplot(clones, aes(x=len_alpha, fill=group)) +
  geom_histogram(binwidth=1, position="dodge", alpha=0.85) +
  scale_fill_manual(values=group_colors, name="Group") +
  theme_minimal(base_size=12) +
  labs(title="CDR3 alpha length", x="AA length", y="N clonotypes")

p_len_b <- ggplot(clones, aes(x=len_beta, fill=group)) +
  geom_histogram(binwidth=1, position="dodge", alpha=0.85) +
  scale_fill_manual(values=group_colors, name="Group") +
  theme_minimal(base_size=12) +
  labs(title="CDR3 beta length", x="AA length", y="N clonotypes")

p_len_box_a <- ggplot(clones, aes(x=group, y=len_alpha, fill=group)) +
  geom_boxplot(alpha=0.8, outlier.size=1) +
  geom_jitter(width=0.15, alpha=0.3, size=1.5) +
  scale_fill_manual(values=group_colors, guide="none") +
  theme_minimal(base_size=12) +
  theme(axis.text.x=element_text(angle=25,hjust=1)) +
  labs(title="CDR3 alpha length", x=NULL, y="AA length")

p_len_box_b <- ggplot(clones, aes(x=group, y=len_beta, fill=group)) +
  geom_boxplot(alpha=0.8, outlier.size=1) +
  geom_jitter(width=0.15, alpha=0.3, size=1.5) +
  scale_fill_manual(values=group_colors, guide="none") +
  theme_minimal(base_size=12) +
  theme(axis.text.x=element_text(angle=25,hjust=1)) +
  labs(title="CDR3 beta length", x=NULL, y="AA length")

fig_len <- (p_len_box_a | p_len_box_b) / (p_len_a | p_len_b) +
  plot_layout(guides="collect") +
  plot_annotation(
    title="Figure 8A — CDR3 length: expanded vs non-expanded clones",
    tag_levels="a"
  )

ggsave(file.path(FIG,"Fig8A_CDR3_length.png"),
       fig_len, width=13, height=10, dpi=300, bg="white")
message("Fig8A saved")

# ── 2. V-GENE USAGE ───────────────────────────────────────────────────────────
message("\n=== 2. V-gene usage ===")

plot_vgene <- function(df, chain="TRB", top_n=10) {
  col_v <- if(chain=="TRB") "macro_TRB_V" else "macro_TRA_V"
  df %>%
    count(group, !!sym(col_v)) %>%
    group_by(group) %>%
    mutate(freq=n/sum(n)) %>%
    ungroup() %>%
    group_by(!!sym(col_v)) %>%
    mutate(max_freq=max(freq)) %>%
    ungroup() %>%
    filter(dense_rank(desc(max_freq)) <= top_n) %>%
    ggplot(aes(x=reorder(!!sym(col_v), freq),
               y=freq, fill=group)) +
    geom_col(position="dodge", alpha=0.85) +
    scale_fill_manual(values=group_colors, name="Group") +
    scale_y_continuous(labels=scales::percent) +
    coord_flip() +
    theme_minimal(base_size=11) +
    labs(title=paste0(chain," V-gene usage (top ",top_n,")"),
         x=NULL, y="Frequency")
}

p_trb_v <- plot_vgene(clones, "TRB", top_n=12)
p_tra_v <- plot_vgene(clones, "TRA", top_n=12)

fig_vgene <- (p_tra_v | p_trb_v) +
  plot_layout(guides="collect") +
  plot_annotation(
    title="Figure 8B — V-gene usage: expanded vs non-expanded",
    tag_levels="a"
  )

ggsave(file.path(FIG,"Fig8B_Vgene_usage.png"),
       fig_vgene, width=15, height=8, dpi=300, bg="white")
message("Fig8B saved")

# ── 3. J-GENE USAGE ───────────────────────────────────────────────────────────
message("\n=== 3. J-gene usage ===")

plot_jgene <- function(df, chain="TRB") {
  col_j <- if(chain=="TRB") "macro_TRB_J" else "macro_TRA_J"
  df %>%
    filter(!is.na(!!sym(col_j)), !!sym(col_j)!="") %>%
    count(group, !!sym(col_j)) %>%
    group_by(group) %>%
    mutate(freq=n/sum(n)) %>%
    ungroup() %>%
    ggplot(aes(x=reorder(!!sym(col_j), freq),
               y=freq, fill=group)) +
    geom_col(position="dodge", alpha=0.85) +
    scale_fill_manual(values=group_colors, name="Group") +
    scale_y_continuous(labels=scales::percent) +
    coord_flip() +
    theme_minimal(base_size=11) +
    labs(title=paste0(chain," J-gene usage"),
         x=NULL, y="Frequency")
}

p_trb_j <- plot_jgene(clones, "TRB")
p_tra_j <- plot_jgene(clones, "TRA")

fig_jgene <- (p_tra_j | p_trb_j) +
  plot_layout(guides="collect") +
  plot_annotation(
    title="Figure 8C — J-gene usage: expanded vs non-expanded",
    tag_levels="a"
  )

ggsave(file.path(FIG,"Fig8C_Jgene_usage.png"),
       fig_jgene, width=15, height=8, dpi=300, bg="white")
message("Fig8C saved")

# ── 4. D-GENE USAGE (beta only) ───────────────────────────────────────────────
message("\n=== 4. D-gene usage ===")

clones_d <- clones %>%
  filter(!is.na(TRB_d_gene), TRB_d_gene != "") %>%
  mutate(macro_TRB_D = str_remove(TRB_d_gene, "-[0-9]+$"))

p_dg <- clones_d %>%
  count(group, macro_TRB_D) %>%
  group_by(group) %>%
  mutate(freq=n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=macro_TRB_D, y=freq, fill=group)) +
  geom_col(position="dodge", alpha=0.85) +
  scale_fill_manual(values=group_colors, name="Group") +
  scale_y_continuous(labels=scales::percent) +
  theme_minimal(base_size=12) +
  labs(title="Figure 8D — TRB D-gene usage: expanded vs non-expanded",
       x="D gene", y="Frequency")

ggsave(file.path(FIG,"Fig8D_Dgene_usage.png"),
       p_dg, width=10, height=5, dpi=300, bg="white")
message("Fig8D saved")

# ── 5. IDROFOBICITÀ CDR3 beta (Kyte-Doolittle) ────────────────────────────────
message("\n=== 5. Idrofobicità CDR3 beta ===")

kd <- c(I=4.5, V=4.2, L=3.8, F=2.8, C=2.5, M=1.9, A=1.8,
        G=-0.4, T=-0.7, W=-0.9, S=-0.8, Y=-1.3, P=-1.6,
        H=-3.2, E=-3.5, Q=-3.5, D=-3.5, N=-3.5, K=-3.9, R=-4.5)

hydrophobicity <- function(seq) {
  aas <- strsplit(seq, "")[[1]]
  scores <- kd[aas]
  scores <- scores[!is.na(scores)]
  if (length(scores) == 0) return(NA_real_)
  mean(scores)
}

charge_aa <- c(K=1, R=1, H=0.1, D=-1, E=-1)
net_charge <- function(seq) {
  aas <- strsplit(seq, "")[[1]]
  scores <- charge_aa[aas]
  scores <- scores[!is.na(scores)]
  sum(scores, na.rm=TRUE)
}

clones <- clones %>%
  rowwise() %>%
  mutate(
    hydro_beta  = hydrophobicity(TRB_cdr3),
    hydro_alpha = hydrophobicity(TRA_cdr3),
    charge_beta = net_charge(TRB_cdr3),
    charge_alpha= net_charge(TRA_cdr3)
  ) %>%
  ungroup()

hydro_stats <- clones %>%
  group_by(group) %>%
  summarise(
    median_hydro_beta  = round(median(hydro_beta, na.rm=TRUE), 3),
    median_hydro_alpha = round(median(hydro_alpha, na.rm=TRUE), 3),
    median_charge_beta = round(median(charge_beta, na.rm=TRUE), 2),
    .groups="drop"
  )
message("Idrofobicità e carica mediana:")
print(hydro_stats)

p_hydro <- ggplot(clones, aes(x=group, y=hydro_beta, fill=group)) +
  geom_boxplot(alpha=0.8, outlier.size=1) +
  geom_jitter(width=0.15, alpha=0.3, size=1.5) +
  scale_fill_manual(values=group_colors, guide="none") +
  theme_minimal(base_size=12) +
  theme(axis.text.x=element_text(angle=25,hjust=1)) +
  labs(title="CDR3 beta — mean hydrophobicity (Kyte-Doolittle)",
       x=NULL, y="Mean hydrophobicity score")

p_charge <- ggplot(clones, aes(x=group, y=charge_beta, fill=group)) +
  geom_boxplot(alpha=0.8, outlier.size=1) +
  geom_jitter(width=0.15, alpha=0.3, size=1.5) +
  scale_fill_manual(values=group_colors, guide="none") +
  theme_minimal(base_size=12) +
  theme(axis.text.x=element_text(angle=25,hjust=1)) +
  labs(title="CDR3 beta — net charge",
       x=NULL, y="Net charge (K/R=+1, D/E=−1)")

p_hydro_a <- ggplot(clones, aes(x=group, y=hydro_alpha, fill=group)) +
  geom_boxplot(alpha=0.8, outlier.size=1) +
  geom_jitter(width=0.15, alpha=0.3, size=1.5) +
  scale_fill_manual(values=group_colors, guide="none") +
  theme_minimal(base_size=12) +
  theme(axis.text.x=element_text(angle=25,hjust=1)) +
  labs(title="CDR3 alpha — mean hydrophobicity",
       x=NULL, y="Mean hydrophobicity score")

fig_phys <- (p_hydro | p_hydro_a | p_charge) +
  plot_layout(guides="collect") +
  plot_annotation(
    title="Figure 8E — CDR3 physicochemical properties: expanded vs non-expanded",
    tag_levels="a"
  )

ggsave(file.path(FIG,"Fig8E_CDR3_physicochemical.png"),
       fig_phys, width=15, height=6, dpi=300, bg="white")
message("Fig8E saved")

# ── 6. COMPOSIZIONE AA CDR3 BETA (posizione per posizione) ────────────────────
message("\n=== 6. Composizione AA CDR3 beta per posizione ===")

# Allinea tutte le CDR3 beta partendo dal C-term (posizione conservata C105 e F/W finale)
# Usa posizione dall'inizio e dalla fine
max_len <- max(nchar(clones$TRB_cdr3), na.rm=TRUE)

aa_by_pos <- clones %>%
  filter(!is.na(TRB_cdr3)) %>%
  mutate(split_aa = strsplit(TRB_cdr3, "")) %>%
  select(group, split_aa, len_beta) %>%
  unnest(split_aa) %>%
  group_by(group) %>%
  mutate(pos = row_number()) %>%
  ungroup()

# Più utile: confronta AA alla stessa posizione relativa (allineamento da N-term)
aa_pos_data <- clones %>%
  filter(!is.na(TRB_cdr3), group %in% c("Expanded (Bo≥5)", "Non-exp Ca", "Non-exp Me")) %>%
  rowwise() %>%
  mutate(aa_list = list(data.frame(
    pos = seq_along(strsplit(TRB_cdr3,"")[[1]]),
    aa  = strsplit(TRB_cdr3,"")[[1]]
  ))) %>%
  ungroup() %>%
  select(group, aa_list) %>%
  unnest(aa_list) %>%
  filter(pos >= 4, pos <= 14) %>%   # CDR3 core (escludi C iniziale e F/W finale fissi)
  count(group, pos, aa) %>%
  group_by(group, pos) %>%
  mutate(freq=n/sum(n)) %>%
  ungroup()

p_aa <- ggplot(aa_pos_data, aes(x=factor(pos), y=aa, fill=freq)) +
  geom_tile(color="white", linewidth=0.3) +
  scale_fill_gradientn(colors=c("white","#FEE8C8","#FC8D59","#D73027"),
                       name="Freq") +
  facet_wrap(~group, ncol=1) +
  theme_minimal(base_size=11) +
  theme(panel.grid=element_blank(),
        axis.text.y=element_text(size=9)) +
  labs(title="Figure 8F — CDR3 beta AA composition by position (core positions 4-14)",
       subtitle="Expanded vs non-expanded clones",
       x="Position in CDR3 beta", y="Amino acid")

ggsave(file.path(FIG,"Fig8F_CDR3_AA_composition.png"),
       p_aa, width=13, height=10, dpi=300, bg="white")
message("Fig8F saved")

# ── 7. TEST STATISTICI ────────────────────────────────────────────────────────
message("\n=== 7. Test statistici (Wilcoxon) ===")

exp_vals  <- clones %>% filter(group=="Expanded (Bo≥5)")
non_exp   <- clones %>% filter(group!="Expanded (Bo≥5)")

tests <- data.frame(
  feature = c("len_beta","len_alpha","hydro_beta","hydro_alpha","charge_beta"),
  p_value = sapply(c("len_beta","len_alpha","hydro_beta","hydro_alpha","charge_beta"), function(f) {
    x <- exp_vals[[f]];  y <- non_exp[[f]]
    x <- x[!is.na(x)];  y <- y[!is.na(y)]
    if(length(x)<3 | length(y)<3) return(NA)
    wilcox.test(x, y)$p.value
  })
) %>%
  mutate(
    median_expanded    = sapply(feature, function(f) round(median(exp_vals[[f]],na.rm=TRUE),2)),
    median_non_expanded= sapply(feature, function(f) round(median(non_exp[[f]],na.rm=TRUE),2)),
    significant        = p_value < 0.05
  )

message("Wilcoxon tests (expanded vs all non-expanded):")
print(tests)

# ── SALVATAGGIO TABELLE ───────────────────────────────────────────────────────
write_xlsx(list(
  "clones_with_features" = clones %>%
    select(group, patient, stage, TRA_v_gene, TRB_v_gene,
           TRA_j_gene, TRB_j_gene, TRB_d_gene,
           TRA_cdr3, TRB_cdr3, len_alpha, len_beta,
           hydro_beta, hydro_alpha, charge_beta, n_cells),
  "stats_summary"  = tests,
  "hydro_by_group" = hydro_stats,
  "length_by_group"= len_stats
), file.path(TAB, "08_tcr_structural_features.xlsx"))

message("\nFigure salvate in: ", FIG)
message("Tabelle salvate in: ", TAB, "/08_tcr_structural_features.xlsx")
message("\nFig8A — CDR3 length (alpha + beta)")
message("Fig8B — V-gene usage")
message("Fig8C — J-gene usage")
message("Fig8D — D-gene usage (beta)")
message("Fig8E — Hydrophobicity + charge")
message("Fig8F — AA composition by position")
