# ==============================================================================
# 09_within_expanded_features.R
#
# Domanda: cosa hanno in comune tra loro i cloni che SI espandono?
#
# Cloni analizzati:
#   - Bo expanded: 27 cloni con >=5 cellule in stage B
#   - Me De novo B: 7 cloni De novo in B (anche se <5 cellule, per confronto)
#
# Analisi:
#   1. J-gene bias (TRBJ2 vs TRBJ1)
#   2. V-gene usage nei soli espansi
#   3. V+J combination (TRBV+TRBJ)
#   4. CDR3 beta: prefissi condivisi / struttura
#   5. Sequence logo CDR3 beta (heatmap AA per posizione)
#   6. Clustering CDR3 beta (edit distance matrix)
#   7. Confronto Bo expanded vs Me De novo
# ==============================================================================

suppressMessages({
  library(dplyr); library(tidyr); library(ggplot2)
  library(readxl); library(writexl); library(stringr)
  library(patchwork); library(forcats); library(scales)
})

BASE <- "/Users/federicadannunzio/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/caruana-project/CART/Code/07_clonotypes"
TAB  <- file.path(BASE, "results", "tables")
FIG  <- file.path(BASE, "results", "figures")

# ── Carica dati ───────────────────────────────────────────────────────────────
esp <- read_xlsx(file.path(TAB, "RISULTATI_expansion_dynamics.xlsx"),
                 sheet = "02_Cloni_espansi_in_B")
fin <- read_xlsx(file.path(TAB, "final_clone_sequences.xlsx"))

j_info <- fin %>% select(TRA_cdr3, TRB_cdr3, TRA_j_gene, TRB_j_gene,
                          TRB_d_gene) %>% distinct()

# Bo expanded con info J-gene
bo_exp <- esp %>%
  filter(patient == "Bo") %>%
  left_join(j_info, by = c("TRA_cdr3","TRB_cdr3")) %>%
  mutate(
    group        = "Bo expanded (≥5)",
    macro_TRB_V  = str_remove(TRB_v_gene, "-[0-9]+$"),
    macro_TRA_V  = str_remove(TRA_v_gene, "-[0-9]+$"),
    macro_TRB_J  = str_remove(TRB_j_gene, "-[0-9]+$"),
    macro_TRA_J  = str_remove(TRA_j_gene, "-[0-9]+$"),
    len_beta     = nchar(TRB_cdr3),
    len_alpha    = nchar(TRA_cdr3)
  )

# Me De novo in B (anche se <5 cellule)
me_b <- fin %>%
  filter(patient == "Me", stage == "B") %>%
  mutate(
    group       = "Me De novo (B)",
    macro_TRB_V = str_remove(TRB_v_gene, "-[0-9]+$"),
    macro_TRA_V = str_remove(TRA_v_gene, "-[0-9]+$"),
    macro_TRB_J = str_remove(TRB_j_gene, "-[0-9]+$"),
    macro_TRA_J = str_remove(TRA_j_gene, "-[0-9]+$"),
    len_beta    = nchar(TRB_cdr3),
    len_alpha   = nchar(TRA_cdr3)
  ) %>%
  rename(TRA_v_gene=TRA_v_gene, TRB_v_gene=TRB_v_gene)

# Tutti i non-expanded (Ca + Me stage I) come background
bg <- fin %>%
  filter(!(patient == "Bo"), !(patient == "Me" & stage == "B")) %>%
  mutate(
    group       = "Non-expanded (Ca+Me I)",
    macro_TRB_V = str_remove(TRB_v_gene, "-[0-9]+$"),
    macro_TRA_V = str_remove(TRA_v_gene, "-[0-9]+$"),
    macro_TRB_J = str_remove(TRB_j_gene, "-[0-9]+$"),
    macro_TRA_J = str_remove(TRA_j_gene, "-[0-9]+$"),
    len_beta    = nchar(TRB_cdr3),
    len_alpha   = nchar(TRA_cdr3)
  )

cols_needed <- c("group","TRA_cdr3","TRB_cdr3","TRA_v_gene","TRB_v_gene",
                 "TRA_j_gene","TRB_j_gene","TRB_d_gene",
                 "macro_TRB_V","macro_TRA_V","macro_TRB_J","macro_TRA_J",
                 "len_beta","len_alpha")

combined <- bind_rows(
  bo_exp[, intersect(colnames(bo_exp), cols_needed)],
  me_b[,  intersect(colnames(me_b),   cols_needed)],
  bg[,    intersect(colnames(bg),     cols_needed)]
) %>%
  mutate(group = factor(group, levels = c("Bo expanded (≥5)",
                                          "Me De novo (B)",
                                          "Non-expanded (Ca+Me I)")))

group_colors <- c(
  "Bo expanded (≥5)"       = "#E64B35",
  "Me De novo (B)"         = "#00A087",
  "Non-expanded (Ca+Me I)" = "#4DBBD5"
)

message("Cloni per gruppo:")
print(table(combined$group))

# ══════════════════════════════════════════════════════════════════════════════
# 1. J-GENE BIAS (TRBJ1 vs TRBJ2)
# ══════════════════════════════════════════════════════════════════════════════
message("\n=== 1. J-gene bias ===")

j_family <- combined %>%
  filter(!is.na(macro_TRB_J)) %>%
  count(group, macro_TRB_J) %>%
  group_by(group) %>%
  mutate(freq = n/sum(n)) %>%
  ungroup()

message("Freq TRBJ per gruppo:")
print(j_family %>% arrange(group, macro_TRB_J))

p_j <- ggplot(j_family, aes(x=macro_TRB_J, y=freq, fill=group)) +
  geom_col(position="dodge", alpha=0.9, width=0.6) +
  geom_text(aes(label=paste0(round(freq*100),"%")),
            position=position_dodge(0.6), vjust=-0.3, size=3.5) +
  scale_fill_manual(values=group_colors, name=NULL) +
  scale_y_continuous(labels=percent, expand=expansion(mult=c(0,0.12))) +
  theme_minimal(base_size=13) +
  theme(panel.grid.major.x=element_blank(),
        legend.position="top") +
  labs(title="TRBJ family usage: expanded vs non-expanded",
       subtitle="TRBJ2 is dominant in expanded clones",
       x="TRBJ family", y="% clonotypes")

# Fisher test TRBJ2 vs TRBJ1: Bo expanded vs background
j_bo  <- bo_exp %>% filter(!is.na(macro_TRB_J)) %>%
  mutate(fam=macro_TRB_J) %>% count(fam)
j_bg_ <- bg %>% filter(!is.na(macro_TRB_J)) %>%
  mutate(fam=macro_TRB_J) %>% count(fam)

make_2x2 <- function(a, b, fam="TRBJ2") {
  a2 <- sum(a$n[a$fam==fam], na.rm=TRUE)
  a1 <- sum(a$n[a$fam!=fam], na.rm=TRUE)
  b2 <- sum(b$n[b$fam==fam], na.rm=TRUE)
  b1 <- sum(b$n[b$fam!=fam], na.rm=TRUE)
  matrix(c(a2,a1,b2,b1), nrow=2,
         dimnames=list(c("TRBJ2","other"),c("expanded","non-exp")))
}
m <- make_2x2(j_bo, j_bg_)
ft <- fisher.test(m)
message(sprintf("Fisher test TRBJ2 in Bo-expanded vs non-expanded: p=%.4f, OR=%.2f",
                ft$p.value, ft$estimate))

# ══════════════════════════════════════════════════════════════════════════════
# 2. V-GENE USAGE NEGLI ESPANSI
# ══════════════════════════════════════════════════════════════════════════════
message("\n=== 2. V-gene usage negli espansi ===")

# Solo Bo expanded: quali V-gene dominano?
bo_vb_counts <- bo_exp %>%
  count(macro_TRB_V, sort=TRUE) %>%
  mutate(freq=n/sum(n), gene=fct_reorder(macro_TRB_V, freq))

message("TRB V-gene in Bo expanded:")
print(bo_vb_counts)

p_vb_exp <- ggplot(bo_vb_counts, aes(x=gene, y=freq)) +
  geom_col(fill="#E64B35", alpha=0.85) +
  geom_text(aes(label=paste0(n, " (", round(freq*100), "%)")),
            hjust=-0.1, size=3.5) +
  scale_y_continuous(labels=percent, expand=expansion(mult=c(0,0.25))) +
  coord_flip() +
  theme_minimal(base_size=12) +
  labs(title="TRB V-gene in Bo expanded clones (27)",
       x=NULL, y="% clonotypes")

# Confronto Bo vs background
vb_compare <- combined %>%
  filter(!is.na(macro_TRB_V)) %>%
  count(group, macro_TRB_V) %>%
  group_by(group) %>%
  mutate(freq=n/sum(n)) %>%
  ungroup() %>%
  filter(macro_TRB_V %in% c("TRBV5","TRBV7","TRBV2","TRBV12","TRBV19","TRBV28","TRBV6"))

p_vb_cmp <- ggplot(vb_compare, aes(x=macro_TRB_V, y=freq, fill=group)) +
  geom_col(position="dodge", alpha=0.85) +
  scale_fill_manual(values=group_colors, name=NULL) +
  scale_y_continuous(labels=percent) +
  theme_minimal(base_size=12) +
  theme(legend.position="top") +
  labs(title="TRBV usage: top genes in expanded vs background",
       x=NULL, y="% clonotypes")

fig_vj <- (p_vb_exp | p_vb_cmp) / p_j +
  plot_annotation(
    title="Figure 9A — V-gene and J-gene usage in expanded clones",
    tag_levels="a"
  )

ggsave(file.path(FIG,"Fig9A_VJ_expanded.png"),
       fig_vj, width=15, height=12, dpi=300, bg="white")
message("Fig9A saved")

# ══════════════════════════════════════════════════════════════════════════════
# 3. V+J COMBINATION (TRBV + TRBJ) NEGLI ESPANSI
# ══════════════════════════════════════════════════════════════════════════════
message("\n=== 3. V+J combination ===")

vj_combo <- bo_exp %>%
  filter(!is.na(macro_TRB_V), !is.na(macro_TRB_J)) %>%
  count(macro_TRB_V, macro_TRB_J, sort=TRUE) %>%
  mutate(combo=paste0(macro_TRB_V, " / ", macro_TRB_J),
         freq=n/sum(n))

message("V+J combinations in Bo expanded:")
print(vj_combo)

p_vj_heat <- ggplot(bo_exp %>%
                      filter(!is.na(macro_TRB_V), !is.na(macro_TRB_J)) %>%
                      count(macro_TRB_V, macro_TRB_J),
                    aes(x=macro_TRB_J, y=macro_TRB_V, fill=n)) +
  geom_tile(color="white") +
  geom_text(aes(label=n), size=5, fontface="bold", color="white") +
  scale_fill_gradientn(colors=c("#FEE8C8","#FC8D59","#D73027"), name="N clones") +
  theme_minimal(base_size=12) +
  theme(panel.grid=element_blank()) +
  labs(title="Figure 9B — TRBV × TRBJ combinations in Bo expanded clones",
       x="TRBJ family", y="TRBV family")

ggsave(file.path(FIG,"Fig9B_VJ_heatmap_expanded.png"),
       p_vj_heat, width=8, height=7, dpi=300, bg="white")
message("Fig9B saved")

# ══════════════════════════════════════════════════════════════════════════════
# 4. CDR3 BETA: PREFIX CONDIVISI
# ══════════════════════════════════════════════════════════════════════════════
message("\n=== 4. Prefissi CDR3 beta (posizioni 1-5) ===")

bo_exp_prefix <- bo_exp %>%
  mutate(prefix5 = substr(TRB_cdr3, 1, 5)) %>%
  count(prefix5, sort=TRUE) %>%
  mutate(freq=n/sum(n))
message("Top prefissi CDR3 beta (5 aa):")
print(bo_exp_prefix %>% filter(n>1))

# ══════════════════════════════════════════════════════════════════════════════
# 5. SEQUENCE LOGO (heatmap AA frequency by position)
# ══════════════════════════════════════════════════════════════════════════════
message("\n=== 5. CDR3 beta AA composition ===")

aa_matrix <- function(seqs, group_lbl, min_pos=1, max_pos=18) {
  seqs <- seqs[!is.na(seqs) & nchar(seqs)>0]
  data.frame(seq=seqs, group=group_lbl) %>%
    rowwise() %>%
    mutate(aa_list=list(data.frame(
      pos=seq_len(nchar(seq)),
      aa=strsplit(seq,"")[[1]]
    ))) %>%
    ungroup() %>%
    select(group, aa_list) %>%
    unnest(aa_list) %>%
    filter(pos >= min_pos, pos <= max_pos)
}

aa_bo  <- aa_matrix(bo_exp$TRB_cdr3,  "Bo expanded (≥5)")
aa_me  <- aa_matrix(me_b$TRB_cdr3,    "Me De novo (B)")
aa_bg_ <- aa_matrix(bg$TRB_cdr3,      "Non-expanded (Ca+Me I)")

aa_all <- bind_rows(aa_bo, aa_me, aa_bg_) %>%
  mutate(group=factor(group, levels=c("Bo expanded (≥5)",
                                      "Me De novo (B)",
                                      "Non-expanded (Ca+Me I)"))) %>%
  count(group, pos, aa) %>%
  group_by(group, pos) %>%
  mutate(freq=n/sum(n)) %>%
  ungroup()

p_logo <- ggplot(aa_all %>% filter(pos >= 3, pos <= 16),
                 aes(x=factor(pos), y=aa, fill=freq)) +
  geom_tile(color="white", linewidth=0.2) +
  geom_text(aes(label=ifelse(freq>0.25, aa, "")),
            size=3, color="white", fontface="bold") +
  scale_fill_gradientn(
    colors=c("white","#FEE8C8","#FC8D59","#D73027"),
    name="AA\nfreq", limits=c(0,1), na.value="white"
  ) +
  facet_wrap(~group, ncol=1) +
  theme_minimal(base_size=11) +
  theme(panel.grid=element_blank(),
        strip.text=element_text(face="bold", size=11)) +
  labs(title="Figure 9C — CDR3 beta amino acid composition by position",
       subtitle="Positions 3-16 | Labelled if frequency >25%",
       x="Position in CDR3 beta", y="Amino acid")

ggsave(file.path(FIG,"Fig9C_CDR3_logo_heatmap.png"),
       p_logo, width=14, height=11, dpi=300, bg="white")
message("Fig9C saved")

# ══════════════════════════════════════════════════════════════════════════════
# 6. CLUSTERING CDR3 BETA (edit distance)
# ══════════════════════════════════════════════════════════════════════════════
message("\n=== 6. Clustering CDR3 beta (edit distance matrix) ===")

seqs_bo <- unique(bo_exp$TRB_cdr3[!is.na(bo_exp$TRB_cdr3)])
d_mat   <- adist(seqs_bo)
rownames(d_mat) <- colnames(d_mat) <- str_trunc(seqs_bo, 18, "right")

# Converti distanza in data frame per heatmap
d_long <- as.data.frame(d_mat) %>%
  mutate(seq1=rownames(d_mat)) %>%
  pivot_longer(-seq1, names_to="seq2", values_to="dist")

p_dist <- ggplot(d_long, aes(x=seq1, y=seq2, fill=dist)) +
  geom_tile(color="white", linewidth=0.2) +
  geom_text(aes(label=ifelse(dist>0 & dist<=3, dist, "")),
            size=2.5, color="white") +
  scale_fill_gradientn(
    colors=c("#D73027","#FC8D59","#FEE8C8","white","white"),
    values=scales::rescale(c(0,1,2,4,max(d_mat))),
    name="Edit\ndist"
  ) +
  theme_minimal(base_size=8) +
  theme(axis.text.x=element_text(angle=90, hjust=1, size=7),
        axis.text.y=element_text(size=7),
        panel.grid=element_blank()) +
  labs(title="Figure 9D — Edit distance between CDR3 beta of Bo expanded clones",
       subtitle="Red = near-identical (dist ≤ 2), white = unrelated",
       x=NULL, y=NULL)

ggsave(file.path(FIG,"Fig9D_CDR3_editdist_matrix.png"),
       p_dist, width=14, height=13, dpi=300, bg="white")
message("Fig9D saved")

# ══════════════════════════════════════════════════════════════════════════════
# 7. RIEPILOGO NUMERICO
# ══════════════════════════════════════════════════════════════════════════════
message("\n=== RIEPILOGO ===")
message("Bo expanded (27 cloni):")
message(sprintf("  TRBJ2: %d/%d (%.0f%%)",
                sum(bo_exp$macro_TRB_J=="TRBJ2", na.rm=TRUE), nrow(bo_exp),
                100*mean(bo_exp$macro_TRB_J=="TRBJ2", na.rm=TRUE)))
message(sprintf("  TRBV5: %d (%.0f%%)  TRBV7: %d (%.0f%%)",
                sum(bo_exp$macro_TRB_V=="TRBV5", na.rm=TRUE),
                100*mean(bo_exp$macro_TRB_V=="TRBV5", na.rm=TRUE),
                sum(bo_exp$macro_TRB_V=="TRBV7", na.rm=TRUE),
                100*mean(bo_exp$macro_TRB_V=="TRBV7", na.rm=TRUE)))
message(sprintf("  CASSP- prefix: %d cloni",
                sum(startsWith(bo_exp$TRB_cdr3,"CASSP"))))
message(sprintf("  CASSSD- prefix: %d cloni",
                sum(startsWith(bo_exp$TRB_cdr3,"CASSSD"))))
message(sprintf("  Fisher test TRBJ2: p=%.4f", ft$p.value))

message("\nMe De novo in B (7 cloni):")
message(sprintf("  TRBJ2: %d/%d (%.0f%%)",
                sum(me_b$macro_TRB_J=="TRBJ2", na.rm=TRUE), nrow(me_b),
                100*mean(me_b$macro_TRB_J=="TRBJ2", na.rm=TRUE)))

# Salva tabella riassuntiva
write_xlsx(list(
  "bo_expanded_genes"  = bo_exp %>%
    select(TRA_cdr3, TRB_cdr3, TRA_v_gene, TRB_v_gene,
           TRA_j_gene, TRB_j_gene, TRB_d_gene,
           macro_TRA_V, macro_TRB_V, macro_TRA_J, macro_TRB_J,
           len_alpha, len_beta, n_cells_B, categoria),
  "vj_combinations"    = vj_combo,
  "j_family_freq"      = j_family,
  "v_freq_bo_expanded" = bo_vb_counts
), file.path(TAB, "09_within_expanded_features.xlsx"))

message("\nSalvato: 09_within_expanded_features.xlsx")
message("Figure: Fig9A (V+J usage), Fig9B (VxJ heatmap), Fig9C (AA logo), Fig9D (edit dist matrix)")
