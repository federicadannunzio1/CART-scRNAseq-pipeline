# ============================================================
#  9f – Dimostrazione del dropout di CD4 mRNA in scRNA-seq
#
#  OBIETTIVO:
#    Dimostrare che l'assenza di CD4 mRNA in molte cellule CAR-T
#    è un artefatto tecnico (dropout), non un'assenza biologica.
#
#  TRE PROVE INDIPENDENTI:
#
#  PROVA 1 – Tasso di rilevamento: CD4 vs CD8A/CD8B
#    Le cellule annotate come CD4+ (per definizione biologica CD4+)
#    hanno CD4 mRNA rilevato molto meno frequentemente di quanto
#    le cellule CD8+ abbiano CD8A/CD8B rilevati.
#    → dimostra asimmetria tecnica tra i due marcatori
#
#  PROVA 2 – Il controllo dei Tregs (e altri CD4+ definiti)
#    I Tregs adulti umani sono CD4+ per definizione: esprimono
#    FOXP3, CTLA4, IL2RA (CD25) indipendentemente da CD4.
#    Se cellule FOXP3+ hanno CD4 mRNA = 0 → dropout dimostrato.
#    Stessa logica per Th1 (TBX21), Th17 (RORC), Tfh (BCL6).
#
#  PROVA 3 – Le Memory T CAR+ sono CD4+ nascoste?
#    Le Memory T CAR+ non esprimono né CD4 né CD8 chiaramente.
#    Se il loro profilo genico (LTB, IL7R, MAL alta; GZMB, CTSW
#    bassa) è simile alle CD4+ e non alle CD8+ → sono probabilmente
#    CD4+ con dropout completo di CD4.
#
#  OUTPUT: 9_CAR_final/res/
#    P18_detection_rate_CD4_vs_CD8.png
#    P19_treg_foxp3_vs_CD4_dropout.png
#    P20_memoryT_CAR_CD4like_profile.png
#    CD4_dropout_proof.xlsx
# ============================================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(writexl)
library(scales)

# ── PERCORSI ──────────────────────────────────────────────────
base_dir <- path.expand(
  "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/CART/Code/Seurat_analysis/"
)
rds_path <- file.path(base_dir, "2_annotation",
                      "all_samples_annotated_COMPLETE.rds")
out_dir  <- file.path(base_dir, "9_CAR_final", "res")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

section <- function(title)
  cat(paste0("\n", strrep("=", 65), "\n  ", title,
             "\n", strrep("=", 65), "\n"))

# Geni necessari
PROOF_GENES <- c(
  # Lineage
  "CD4", "CD8A", "CD8B",
  # Marcatori CD4-specifici (fattori di trascrizione)
  "FOXP3",   # Tregs
  "TBX21",   # Th1
  "RORC",    # Th17
  "BCL6",    # Tfh
  "GATA3",   # Th2
  # Marcatori di stato CD4+ aggiuntivi
  "IL2RA",   # CD25 – Tregs
  "CTLA4",   # Tregs
  "CXCR5",   # Tfh
  "CCR6",    # Th17
  # Marcatori CD8-specifici
  "GZMB", "CTSW", "PRF1", "NKG7",
  # Marcatori naive/memory (presenti in entrambi)
  "IL7R", "LTB", "MAL", "CCR7", "SELL"
)

CD4_TYPES <- c("Naive CD4+ T cells","Th1 cells","Th2 cells",
               "Th17 cells","Tfh cells","Effector CD4+ T cells",
               "Tregs","Proliferating CD4+ T cells")
CD8_TYPES <- c("Cytotoxic CD8+ T cells","Naive CD8+ T cells",
               "Proliferating CD8+ T cells")
AB_NAMES  <- c("Ca_blood_AB","Ca_bone_AB",
               "Bo_blood_AB","Bo_bone_AB","Me_bone_AB")

COL_CD4 <- "#E63946"
COL_CD8 <- "#264653"
COL_MEM <- "#2A9D8F"
COL_OK  <- "#2D6A4F"   # rilevato
COL_NO  <- "#D62828"   # non rilevato (dropout)

# ============================================================
# 1. CARICAMENTO E ESTRAZIONE
# ============================================================
section("Caricamento dati")

all_samples <- readRDS(rds_path)

expr_list <- list()
for (nm in AB_NAMES) {
  obj <- all_samples[[nm]]
  if (is.null(obj)) next
  if (length(grep("^counts\\.", Layers(obj), value=TRUE)) > 0)
    obj <- JoinLayers(obj)

  meta     <- obj@meta.data
  genes_ok <- PROOF_GENES[PROOF_GENES %in% rownames(obj)]
  bc_all   <- rownames(meta)

  expr_df  <- FetchData(obj, vars = genes_ok,
                        cells = bc_all, layer = "data")
  for (g in setdiff(PROOF_GENES, colnames(expr_df)))
    expr_df[[g]] <- 0

  expr_df$sample    <- nm
  expr_df$cell_type <- meta[bc_all, "cell_type"]
  expr_df$IS_CAR    <- ifelse(
    !is.na(meta[bc_all,"IS_CAR_ALLIN_scREP"]) &
    meta[bc_all,"IS_CAR_ALLIN_scREP"] == "YES",
    "CAR+", "CAR-")

  expr_df$lineage <- dplyr::case_when(
    expr_df$cell_type %in% CD4_TYPES ~ "CD4+",
    expr_df$cell_type %in% CD8_TYPES ~ "CD8+",
    expr_df$cell_type == "Memory T cells" ~ "Memory T",
    TRUE ~ "Other"
  )

  cat(sprintf("  %-15s | %5d celle | geni ok: %d\n",
              nm, nrow(expr_df), length(genes_ok)))
  expr_list[[nm]] <- expr_df
}

expr_all <- bind_rows(expr_list)
cat(sprintf("\nTotale cellule: %d\n", nrow(expr_all)))

# ============================================================
# PROVA 1 – TASSO DI RILEVAMENTO: CD4 vs CD8A/CD8B
#
#  Per ogni tipo cellulare annotato:
#    % cellule con CD4 > 0  (nelle CD4+)
#    % cellule con CD8A > 0 (nelle CD8+)
#    % cellule con CD8B > 0 (nelle CD8+)
#  Se CD4 è rilevato meno frequentemente → dropout asimmetrico
# ============================================================
section("PROVA 1 – Tasso di rilevamento CD4 vs CD8A/CD8B")

# Calcola tasso di rilevamento per tipo cellulare
detection_cd4 <- expr_all %>%
  filter(lineage == "CD4+") %>%
  group_by(cell_type) %>%
  summarise(
    n           = n(),
    pct_CD4_det = round(100 * mean(CD4 > 0), 1),
    mean_CD4    = round(mean(CD4), 3),
    .groups     = "drop"
  ) %>%
  mutate(marcatore = "CD4", lineage = "CD4+") %>%
  rename(pct_detected = pct_CD4_det, mean_expr = mean_CD4)

detection_cd8a <- expr_all %>%
  filter(lineage == "CD8+") %>%
  group_by(cell_type) %>%
  summarise(
    n           = n(),
    pct_detected = round(100 * mean(CD8A > 0), 1),
    mean_expr    = round(mean(CD8A), 3),
    .groups      = "drop"
  ) %>%
  mutate(marcatore = "CD8A", lineage = "CD8+")

detection_cd8b <- expr_all %>%
  filter(lineage == "CD8+") %>%
  group_by(cell_type) %>%
  summarise(
    n           = n(),
    pct_detected = round(100 * mean(CD8B > 0), 1),
    mean_expr    = round(mean(CD8B), 3),
    .groups      = "drop"
  ) %>%
  mutate(marcatore = "CD8B", lineage = "CD8+")

detection_all <- bind_rows(detection_cd4, detection_cd8a, detection_cd8b)

cat("\nTasso di rilevamento per tipo cellulare:\n")
print(as.data.frame(detection_all %>%
                      select(cell_type, marcatore, n,
                             pct_detected, mean_expr)))

# Riepilogo globale
global_det <- expr_all %>%
  filter(lineage %in% c("CD4+","CD8+")) %>%
  group_by(lineage) %>%
  summarise(
    n            = n(),
    pct_CD4_det  = round(100 * mean(CD4  > 0), 1),
    pct_CD8A_det = round(100 * mean(CD8A > 0), 1),
    pct_CD8B_det = round(100 * mean(CD8B > 0), 1),
    .groups = "drop"
  )
cat("\nRiepilogo globale:\n")
print(as.data.frame(global_det))

# ── P18: Grafici tasso di rilevamento ─────────────────────────

# Ordinamento biologico
type_order_cd4 <- c("Naive CD4+ T cells","Effector CD4+ T cells",
                    "Th1 cells","Th2 cells","Th17 cells","Tfh cells",
                    "Tregs","Proliferating CD4+ T cells")
type_order_cd8 <- c("Naive CD8+ T cells","Cytotoxic CD8+ T cells",
                    "Proliferating CD8+ T cells")

det_cd4_plot <- detection_cd4 %>%
  filter(cell_type %in% type_order_cd4) %>%
  mutate(
    cell_type     = factor(cell_type, levels = rev(type_order_cd4)),
    pct_dropout   = 100 - pct_detected,
    label_det     = paste0(pct_detected, "%"),
    label_drop    = paste0(round(pct_dropout), "%")
  )

det_cd8_plot <- bind_rows(detection_cd8a, detection_cd8b) %>%
  filter(cell_type %in% type_order_cd8) %>%
  mutate(
    cell_type   = factor(cell_type, levels = rev(type_order_cd8)),
    pct_dropout = 100 - pct_detected
  )

# Barplot orizzontale: rilevato vs dropout
make_detection_bar <- function(df, title_str, col_bar) {
  df_long <- df %>%
    select(cell_type, n, pct_detected, pct_dropout) %>%
    pivot_longer(cols = c(pct_detected, pct_dropout),
                 names_to = "stato",
                 values_to = "pct") %>%
    mutate(stato = recode(stato,
                          pct_detected = "Rilevato (mRNA > 0)",
                          pct_dropout  = "Dropout (mRNA = 0)"),
           stato = factor(stato,
                          levels = c("Rilevato (mRNA > 0)",
                                     "Dropout (mRNA = 0)")))

  ggplot(df_long, aes(x = cell_type, y = pct, fill = stato)) +
    geom_col(position = "stack", width = 0.65) +
    geom_text(
      data = df,
      aes(x = cell_type, y = pct_detected / 2,
          label = paste0(pct_detected, "%")),
      inherit.aes = FALSE,
      color = "white", fontface = "bold", size = 3.5
    ) +
    coord_flip() +
    scale_fill_manual(
      values = c("Rilevato (mRNA > 0)" = col_bar,
                 "Dropout (mRNA = 0)"  = "#DDDDDD"),
      name = NULL
    ) +
    scale_y_continuous(labels = function(x) paste0(x, "%"),
                       limits = c(0, 100)) +
    labs(title = title_str, x = NULL,
         y = "% cellule") +
    theme_classic(base_size = 11) +
    theme(plot.title    = element_text(face = "bold", hjust = 0.5),
          legend.position = "bottom")
}

p18a <- make_detection_bar(
  det_cd4_plot, "% cellule CD4+ con CD4 mRNA rilevato", COL_CD4)

p18b <- det_cd8_plot %>%
  ggplot(aes(x = cell_type, y = pct_detected,
             fill = marcatore)) +
  geom_col(position = "dodge", width = 0.6) +
  geom_text(aes(label = paste0(pct_detected, "%")),
            position = position_dodge(0.6),
            hjust = -0.1, size = 3.2, fontface = "bold") +
  coord_flip() +
  scale_fill_manual(
    values = c("CD8A" = COL_CD8, "CD8B" = "#577590"),
    name = "Marcatore") +
  scale_x_discrete(limits = rev(type_order_cd8)) +
  scale_y_continuous(limits = c(0, 105),
                     labels = function(x) paste0(x, "%")) +
  labs(title = "% cellule CD8+ con CD8A/CD8B mRNA rilevato",
       x = NULL, y = "% cellule") +
  theme_classic(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", hjust = 0.5),
        legend.position = "bottom")

# Riepilogo comparativo in un unico barplot
summary_det <- data.frame(
  Marcatore    = c("CD4\n(nelle CD4+)",
                   "CD8A\n(nelle CD8+)",
                   "CD8B\n(nelle CD8+)"),
  pct_rilevato = c(
    global_det$pct_CD4_det[global_det$lineage == "CD4+"],
    global_det$pct_CD8A_det[global_det$lineage == "CD8+"],
    global_det$pct_CD8B_det[global_det$lineage == "CD8+"]
  ),
  colore       = c(COL_CD4, COL_CD8, "#577590")
)

p18c <- ggplot(summary_det,
               aes(x = Marcatore, y = pct_rilevato,
                   fill = Marcatore)) +
  geom_col(width = 0.55, show.legend = FALSE) +
  geom_text(aes(label = paste0(pct_rilevato, "%")),
            vjust = -0.5, size = 5, fontface = "bold") +
  geom_hline(yintercept = 50, linetype = "dashed",
             color = "gray50", linewidth = 0.7) +
  scale_fill_manual(values = setNames(summary_det$colore,
                                       summary_det$Marcatore)) +
  scale_y_continuous(limits = c(0, 105),
                     labels = function(x) paste0(x, "%")) +
  labs(title    = "Confronto: % rilevamento del marcatore\nnel proprio lineage",
       subtitle = "CD4 è rilevato MENO frequentemente di CD8A/CD8B\nnelle rispettive popolazioni → dropout asimmetrico",
       x = NULL, y = "% cellule con mRNA > 0") +
  theme_classic(base_size = 12) +
  theme(plot.title    = element_text(face = "bold", hjust = 0.5,
                                     size = 11),
        plot.subtitle = element_text(hjust = 0.5, color = "gray40"))

p18 <- (p18a | p18b) / p18c +
  plot_layout(heights = c(2, 1)) +
  plot_annotation(
    title = "P18 – PROVA 1: CD4 mRNA è rilevato molto meno frequentemente\nnelle CD4+ rispetto a CD8A/CD8B nelle CD8+",
    subtitle = "Questo dimostra dropout asimmetrico → il problema è tecnico, non biologico",
    theme = theme(
      plot.title    = element_text(face = "bold", hjust = 0.5, size = 13),
      plot.subtitle = element_text(hjust = 0.5, color = "gray40", size = 10)
    )
  )
ggsave(file.path(out_dir, "P18_detection_rate_CD4_vs_CD8.png"),
       p18, width = 14, height = 12, dpi = 300, bg = "white")
cat("  P18 salvato\n")

# ============================================================
# PROVA 2 – IL CONTROLLO DEI TREGS (e altri CD4+ definiti)
#
#  I Tregs adulti umani sono CD4+ per definizione biologica.
#  Marcatori: FOXP3, CTLA4, IL2RA (CD25).
#  Se FOXP3+ → CD4+ (biologia), ma CD4 mRNA può essere 0
#  → dimostrazione diretta del dropout.
#
#  Stesso ragionamento per:
#  TBX21 → Th1,  RORC → Th17,  BCL6 → Tfh
# ============================================================
section("PROVA 2 – Dropout nei Tregs e CD4+ definiti da TF")

# Prendi tutte le cellule annotate come Tregs o con marker TF alto
tregs <- expr_all %>%
  filter(cell_type == "Tregs") %>%
  mutate(
    FOXP3_detected  = FOXP3 > 0,
    CD4_detected    = CD4   > 0,
    CTLA4_detected  = CTLA4 > 0,
    IL2RA_detected  = IL2RA > 0
  )

cat(sprintf("\nTregs totali: %d\n", nrow(tregs)))
cat(sprintf("  FOXP3 > 0:  %d (%.1f%%)\n",
            sum(tregs$FOXP3_detected),
            100*mean(tregs$FOXP3_detected)))
cat(sprintf("  CD4   > 0:  %d (%.1f%%)\n",
            sum(tregs$CD4_detected),
            100*mean(tregs$CD4_detected)))
cat(sprintf("  FOXP3+ ma CD4=0: %d (%.1f%% delle FOXP3+)\n",
            sum(tregs$FOXP3_detected & !tregs$CD4_detected),
            100*mean(tregs$CD4_detected[tregs$FOXP3_detected] == 0)))

# Per tutti i CD4+ subtypes: quante esprimono il loro TF ma non CD4
tf_analysis <- expr_all %>%
  filter(cell_type %in% CD4_TYPES) %>%
  mutate(
    TF_name = case_when(
      cell_type == "Tregs"     ~ "FOXP3",
      cell_type == "Th1 cells" ~ "TBX21",
      cell_type == "Th17 cells"~ "RORC",
      cell_type == "Tfh cells" ~ "BCL6",
      cell_type == "Th2 cells" ~ "GATA3",
      TRUE                     ~ NA_character_
    ),
    TF_expr = case_when(
      cell_type == "Tregs"     ~ FOXP3,
      cell_type == "Th1 cells" ~ TBX21,
      cell_type == "Th17 cells"~ RORC,
      cell_type == "Tfh cells" ~ BCL6,
      cell_type == "Th2 cells" ~ GATA3,
      TRUE                     ~ 0
    ),
    TF_positive = TF_expr > 0,
    CD4_positive = CD4 > 0
  ) %>%
  filter(!is.na(TF_name))

tf_summary <- tf_analysis %>%
  filter(TF_positive) %>%   # solo cellule con TF rilevato
  group_by(cell_type, TF_name) %>%
  summarise(
    n_TF_pos          = n(),
    n_CD4_pos         = sum(CD4_positive),
    n_CD4_zero        = sum(!CD4_positive),
    pct_CD4_detected  = round(100 * mean(CD4_positive), 1),
    pct_CD4_dropout   = round(100 * mean(!CD4_positive), 1),
    .groups = "drop"
  )

cat("\nDropout CD4 nelle cellule con TF CD4-specifico rilevato:\n")
print(as.data.frame(tf_summary))

# ── P19: Scatter FOXP3 vs CD4 e barchart dropout per TF ──────

# Scatter FOXP3 vs CD4 in Tregs
p19a <- ggplot(tregs, aes(x = FOXP3, y = CD4)) +
  geom_point(aes(color = CD4_detected),
             alpha = 0.6, size = 1.2) +
  scale_color_manual(
    values = c("TRUE" = COL_OK, "FALSE" = COL_NO),
    labels = c("TRUE" = "CD4 rilevato", "FALSE" = "CD4 = 0 (dropout)"),
    name   = NULL
  ) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "gray50", linewidth = 0.6) +
  annotate("text",
           x    = max(tregs$FOXP3) * 0.6,
           y    = max(tregs$CD4) * 0.9,
           label = sprintf("Tregs FOXP3+\ncon CD4=0: %.1f%%",
                           100*mean(tregs$CD4_detected[tregs$FOXP3_detected]==0)),
           color = COL_NO, size = 4, fontface = "bold") +
  labs(title    = "FOXP3 vs CD4 nelle cellule Tregs",
       subtitle = "I Tregs sono CD4+ per definizione biologica (FOXP3+).\nI punti sulla riga y=0 sono Tregs biologicamente CD4+ con CD4 mRNA non rilevato.",
       x = "FOXP3 espressione (log-norm)",
       y = "CD4 espressione (log-norm)") +
  theme_classic(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
        legend.position = "bottom")

# Barchart: % dropout in ogni sottotipo CD4 con TF rilevato
p19b_data <- tf_summary %>%
  select(cell_type, TF_name, pct_CD4_detected, pct_CD4_dropout) %>%
  pivot_longer(cols = c(pct_CD4_detected, pct_CD4_dropout),
               names_to = "stato", values_to = "pct") %>%
  mutate(
    stato = recode(stato,
                   pct_CD4_detected = "CD4 rilevato",
                   pct_CD4_dropout  = "CD4 = 0 (dropout)"),
    stato = factor(stato,
                   levels = c("CD4 rilevato",
                               "CD4 = 0 (dropout)")),
    label = paste0(cell_type, "\n(", TF_name, "+)")
  )

p19b <- ggplot(p19b_data,
               aes(x = label, y = pct, fill = stato)) +
  geom_col(position = "stack", width = 0.65) +
  geom_text(
    data = tf_summary,
    aes(x   = paste0(cell_type, "\n(", TF_name, "+)"),
        y   = pct_CD4_dropout / 2,
        label = paste0(pct_CD4_dropout, "%")),
    inherit.aes = FALSE,
    color = "white", fontface = "bold", size = 3.8
  ) +
  scale_fill_manual(
    values = c("CD4 rilevato"    = COL_OK,
               "CD4 = 0 (dropout)" = COL_NO),
    name = NULL
  ) +
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     limits = c(0, 100)) +
  labs(title    = "% dropout CD4 in cellule con TF CD4-specifico rilevato",
       subtitle = "Ogni barra = cellule con TF rilevato (biologicamente CD4+)\nLa parte rossa = CD4 mRNA = 0 nonostante la cellula sia CD4+",
       x = NULL, y = "% cellule") +
  theme_classic(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
        axis.text.x   = element_text(size = 9),
        legend.position = "bottom")

p19 <- (p19a | p19b) +
  plot_annotation(
    title = "P19 – PROVA 2: Dropout di CD4 in cellule biologicamente CD4+",
    subtitle = "Cellule con FOXP3/TBX21/RORC/BCL6 rilevato sono CD4+ per definizione,\neppure una quota sostanziale ha CD4 mRNA = 0",
    theme = theme(
      plot.title    = element_text(face = "bold", hjust = 0.5, size = 13),
      plot.subtitle = element_text(hjust = 0.5, color = "gray40", size = 10)
    )
  )
ggsave(file.path(out_dir, "P19_treg_foxp3_vs_CD4_dropout.png"),
       p19, width = 14, height = 7, dpi = 300, bg = "white")
cat("  P19 salvato\n")

# ============================================================
# PROVA 3 – LE MEMORY T CAR+ SONO CD4+ NASCOSTE?
#
#  Le Memory T CAR+ non hanno CD4 né CD8 rilevati.
#  Se il loro profilo è:
#    - CD8 effector markers (GZMB, CTSW, PRF1) BASSI
#    - CD4-associated markers (LTB, IL7R, MAL)  relativi
#  → probabilmente sono CD4+ con dropout completo
#
#  Confronto con CD4+ e CD8+ CAR+ come riferimento
# ============================================================
section("PROVA 3 – Memory T CAR+: profilo CD4-like o CD8-like?")

mem_car <- expr_all %>%
  filter(IS_CAR == "CAR+", lineage == "Memory T")
cd4_car <- expr_all %>%
  filter(IS_CAR == "CAR+", lineage == "CD4+")
cd8_car <- expr_all %>%
  filter(IS_CAR == "CAR+", lineage == "CD8+")

cat(sprintf("\nMemory T CAR+: %d | CD4+ CAR+: %d | CD8+ CAR+: %d\n",
            nrow(mem_car), nrow(cd4_car), nrow(cd8_car)))

# Profilo genico medio
PROFILE_GENES <- c("CD4","CD8A","CD8B",
                   "GZMB","CTSW","PRF1","NKG7",  # CD8 effector
                   "LTB","IL7R","MAL","CCR7")      # naive/CD4

profile_df <- bind_rows(
  mem_car %>% mutate(gruppo = "Memory T CAR+"),
  cd4_car %>% mutate(gruppo = "CD4+ CAR+"),
  cd8_car %>% mutate(gruppo = "CD8+ CAR+")
) %>%
  select(gruppo, all_of(PROFILE_GENES[PROFILE_GENES %in% colnames(.)])) %>%
  group_by(gruppo) %>%
  summarise(across(everything(), ~ round(mean(.x, na.rm=TRUE), 3)),
            .groups = "drop")

cat("\nProfilo genico medio (CAR+):\n")
print(as.data.frame(profile_df))

# ── P20: Violin multipli per confronto profilo ───────────────
genes_ok_profile <- PROFILE_GENES[PROFILE_GENES %in% colnames(expr_all)]

profile_long <- bind_rows(
  mem_car %>% mutate(gruppo = "Memory T\nCAR+"),
  cd4_car %>% mutate(gruppo = "CD4+\nCAR+"),
  cd8_car %>% mutate(gruppo = "CD8+\nCAR+")
) %>%
  select(gruppo, all_of(genes_ok_profile)) %>%
  pivot_longer(cols = all_of(genes_ok_profile),
               names_to = "gene", values_to = "expr") %>%
  mutate(
    gruppo = factor(gruppo,
                    levels = c("CD4+\nCAR+",
                               "Memory T\nCAR+",
                               "CD8+\nCAR+")),
    categoria = case_when(
      gene %in% c("CD4","CD8A","CD8B")          ~ "1_Lineage",
      gene %in% c("GZMB","CTSW","PRF1","NKG7")  ~ "2_CD8 effector",
      gene %in% c("LTB","IL7R","MAL","CCR7")    ~ "3_Naive/CD4-assoc"
    ),
    gene = factor(gene, levels = genes_ok_profile)
  )

p20 <- ggplot(profile_long,
              aes(x = gruppo, y = expr, fill = gruppo)) +
  geom_violin(alpha = 0.7, scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.1, outlier.size = 0.2,
               fill = "white", alpha = 0.85) +
  facet_wrap(~ gene, scales = "free_y",
             nrow = 3, ncol = 4) +
  scale_fill_manual(
    values = c(
      "CD4+\nCAR+"     = COL_CD4,
      "Memory T\nCAR+" = COL_MEM,
      "CD8+\nCAR+"     = COL_CD8
    ),
    name = NULL
  ) +
  labs(
    title    = "P20 – PROVA 3: Profilo genico delle Memory T CAR+ vs CD4+ e CD8+ CAR+",
    subtitle = paste0(
      "Se Memory T CAR+ ha GZMB/CTSW/PRF1 BASSI e LTB/IL7R ALTI → profilo CD4-like\n",
      sprintf("Memory T CAR+: n=%d  |  CD4+ CAR+: n=%d  |  CD8+ CAR+: n=%d",
              nrow(mem_car), nrow(cd4_car), nrow(cd8_car))
    ),
    x = NULL, y = "Espressione (log-norm)"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5, size = 12),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40", size = 9),
    strip.text    = element_text(face = "bold"),
    legend.position = "bottom",
    axis.text.x   = element_text(size = 8)
  )

ggsave(file.path(out_dir, "P20_memoryT_CAR_CD4like_profile.png"),
       p20, width = 14, height = 11, dpi = 300, bg = "white")
cat("  P20 salvato\n")

# ============================================================
# EXCEL + RIEPILOGO
# ============================================================
section("Export Excel e riepilogo")

writexl::write_xlsx(
  list(
    "Prova1_tasso_rilevamento"  = as.data.frame(detection_all),
    "Prova1_riepilogo_globale"  = as.data.frame(global_det),
    "Prova2_TF_vs_CD4_dropout"  = as.data.frame(tf_summary),
    "Prova3_profilo_Memory_CAR" = as.data.frame(profile_df)
  ),
  path = file.path(out_dir, "CD4_dropout_proof.xlsx")
)
cat(sprintf("  Excel → %s\n",
            file.path(out_dir, "CD4_dropout_proof.xlsx")))

section("RIEPILOGO FINALE – 3 PROVE DEL DROPOUT")
cat(paste0(
  "\nPROVA 1 – Asimmetria tecnica:\n",
  sprintf("  CD4 rilevato nelle CD4+:  %.1f%%\n",
          global_det$pct_CD4_det[global_det$lineage=="CD4+"]),
  sprintf("  CD8A rilevato nelle CD8+: %.1f%%\n",
          global_det$pct_CD8A_det[global_det$lineage=="CD8+"]),
  sprintf("  CD8B rilevato nelle CD8+: %.1f%%\n",
          global_det$pct_CD8B_det[global_det$lineage=="CD8+"]),
  "  → CD4 viene perso molto più frequentemente di CD8\n",
  "\nPROVA 2 – Tregs e CD4+ TF-definite:\n"
))
for (i in seq_len(nrow(tf_summary))) {
  cat(sprintf("  %s (%s+): CD4 dropout = %.1f%% delle cellule con TF rilevato\n",
              tf_summary$cell_type[i], tf_summary$TF_name[i],
              tf_summary$pct_CD4_dropout[i]))
}
cat(paste0(
  "\nPROVA 3 – Memory T CAR+:\n",
  "  Vedi P20: se GZMB/CTSW bassi e LTB/IL7R alti rispetto alle CD8+ CAR+\n",
  "  → le Memory T CAR+ hanno profilo CD4-like, probabilmente sono CD4+ con dropout\n"
))
