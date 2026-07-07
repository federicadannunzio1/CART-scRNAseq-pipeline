# ============================================================
#  9e – Identità delle cellule CD8αα-like + ricalcolo rapporto CAR CD4/CD8
#
#  OBIETTIVI:
#    1. Capire cosa sono le CD8αα-like (CD8A+, CD8B=0):
#       MAIT cells, NKT cells, γδ T cells, o altro?
#       → marcatori: SLC4A10, KLRB1, ZBTB16, TRDC, TRGC1, NCR3, DPP4
#
#    2. Ricalcolare il rapporto CAR-T CD4/CD8 in tre scenari:
#       BASELINE : annotazione attuale (CD8αα-like conteggiate come CD8+)
#       SCENARIO A: escludi le CD8αα-like dall'analisi
#       SCENARIO B: riannotale come MAIT/NKT/γδ (non più CD8+)
#       → quanto cambia il rapporto CD4/CD8?
#
#  OUTPUT: 9_CAR_final/res/
#    P16_CD8alike_identity.png
#    P17_ratio_CD4_CD8_scenarios.png
#    CD8alike_identity_and_ratio.xlsx
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

# ── Categorie ─────────────────────────────────────────────────
CD4_TYPES <- c("Naive CD4+ T cells","Th1 cells","Th2 cells",
               "Th17 cells","Tfh cells","Effector CD4+ T cells",
               "Tregs","Proliferating CD4+ T cells")
CD8_TYPES <- c("Cytotoxic CD8+ T cells","Naive CD8+ T cells",
               "Proliferating CD8+ T cells")

AB_NAMES  <- c("Ca_blood_AB","Ca_bone_AB",
               "Bo_blood_AB","Bo_bone_AB","Me_bone_AB")

# Geni per identificare le CD8αα-like
# Gerarchia di specificità:
#   TRDC/TRGC1  → γδ T (quasi esclusivi)
#   SLC4A10     → MAIT (quasi esclusivo negli adulti)
#   KLRB1+ZBTB16→ NKT (PLZF = master TF delle iNKT)
#   KLRB1 solo  → MAIT o NKT (CD161+, ambiguo)
#   NCR3        → MAIT + NK
#   DPP4        → iNKT marker aggiuntivo
IDENTITY_GENES <- c("SLC4A10","KLRB1","ZBTB16","TRDC","TRGC1",
                    "NCR3","DPP4","IL18R1","RORC")
LINEA_GENES    <- c("CD4","CD8A","CD8B")
ALL_GENES      <- unique(c(LINEA_GENES, IDENTITY_GENES))

# ============================================================
# 1. CARICAMENTO + ESTRAZIONE ESPRESSIONE
# ============================================================
section("Caricamento e estrazione dati")

all_samples <- readRDS(rds_path)

expr_list <- list()

for (nm in AB_NAMES) {
  obj <- all_samples[[nm]]
  if (is.null(obj)) next
  if (length(grep("^counts\\.", Layers(obj), value=TRUE)) > 0)
    obj <- JoinLayers(obj)

  meta     <- obj@meta.data
  genes_ok <- ALL_GENES[ALL_GENES %in% rownames(obj)]

  # Tutte le cellule (per il rapporto globale)
  bc_all <- rownames(meta)
  expr_df <- FetchData(obj, vars = genes_ok,
                       cells = bc_all, layer = "data")
  for (g in setdiff(ALL_GENES, colnames(expr_df))) expr_df[[g]] <- 0

  expr_df$sample    <- nm
  expr_df$cell_type <- meta[bc_all, "cell_type"]
  expr_df$IS_CAR    <- ifelse(
    !is.na(meta[bc_all,"IS_CAR_ALLIN_scREP"]) &
    meta[bc_all,"IS_CAR_ALLIN_scREP"] == "YES", "CAR+", "CAR-")
  expr_df$barcode   <- rownames(expr_df)

  expr_df$lineage <- dplyr::case_when(
    expr_df$cell_type %in% CD4_TYPES ~ "CD4+",
    expr_df$cell_type %in% CD8_TYPES ~ "CD8+",
    TRUE                              ~ "Other"
  )

  cat(sprintf("  %-15s | celle: %5d | CAR+: %4d | geni: %d\n",
              nm, nrow(expr_df),
              sum(expr_df$IS_CAR == "CAR+"),
              length(genes_ok)))
  expr_list[[nm]] <- expr_df
}

expr_all <- bind_rows(expr_list)

# ── Classifica CD8αβ vs CD8αα-like ───────────────────────────
expr_all <- expr_all %>%
  mutate(
    cd8_class = case_when(
      lineage == "CD8+" & CD8A > 0 & CD8B > 0  ~ "CD8ab_true",
      lineage == "CD8+" & CD8A > 0 & CD8B == 0 ~ "CD8aa_like",
      lineage == "CD8+"                         ~ "CD8_dropout",
      TRUE                                      ~ "non_CD8"
    )
  )

# Sottoinsieme: solo CD8αα-like
expr_aa <- expr_all %>%
  filter(cd8_class == "CD8aa_like")

cat(sprintf("\nCD8αα-like totali (AB): %d\n", nrow(expr_aa)))
cat(sprintf("di cui CAR+: %d\n", sum(expr_aa$IS_CAR == "CAR+")))

# ============================================================
# 2. CLASSIFICAZIONE DELL'IDENTITÀ DELLE CD8αα-LIKE
#
#  Gerarchia (dall'identificatore più specifico al meno):
#    1. TRDC > 0 OR TRGC1 > 0  → γδ T cell
#    2. SLC4A10 > 0             → MAIT cell (quasi esclusivo)
#    3. ZBTB16 > 0              → iNKT cell (PLZF)
#    4. KLRB1 > 0               → CD161+ (MAIT o NKT, ambiguo)
#    5. Nessuno                 → CD8αα innate-like NEC
#       (potrebbe essere CD8αβ con dropout di CD8B)
# ============================================================
section("Classificazione identità CD8αα-like")

classify_aa <- function(trdc, trgc1, slc4a10, zbtb16, klrb1) {
  dplyr::case_when(
    trdc > 0 | trgc1 > 0   ~ "γδ T cell",
    slc4a10 > 0             ~ "MAIT cell",
    zbtb16  > 0             ~ "iNKT cell",
    klrb1   > 0             ~ "CD161+ (MAIT/NKT ambiguo)",
    TRUE                    ~ "Innate-like NEC\n(CD8B dropout?)"
  )
}

expr_aa <- expr_aa %>%
  mutate(
    identity = classify_aa(TRDC, TRGC1, SLC4A10, ZBTB16, KLRB1)
  )

# Distribuzione identità
id_summary_all <- expr_aa %>%
  count(identity) %>%
  mutate(pct = round(100 * n / sum(n), 1),
         gruppo = "CD8αα-like totali")

id_summary_car <- expr_aa %>%
  filter(IS_CAR == "CAR+") %>%
  count(identity) %>%
  mutate(pct = round(100 * n / sum(n), 1),
         gruppo = "CD8αα-like CAR+")

cat("\nIdentità CD8αα-like (tutte):\n")
print(as.data.frame(id_summary_all))
cat("\nIdentità CD8αα-like (solo CAR+):\n")
print(as.data.frame(id_summary_car))

# Per campione + CAR status
id_per_sample <- expr_aa %>%
  group_by(sample, IS_CAR, identity) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(sample, IS_CAR) %>%
  mutate(pct = round(100 * n / sum(n), 1))

# ============================================================
# 3. RICALCOLO RAPPORTO CD4/CD8 CAR-T IN 3 SCENARI
# ============================================================
section("Ricalcolo rapporto CD4/CD8 CAR-T")

# Lavora solo sulle cellule CAR+ nei campioni AB
car_cells <- expr_all %>%
  filter(IS_CAR == "CAR+")

cat(sprintf("\nTotale CAR+ (AB): %d\n", nrow(car_cells)))
cat("\nDistribuzione per lineage (annotazione attuale):\n")
print(table(car_cells$lineage))
cat("\nDistribuzione per cell_type (annotazione attuale):\n")
print(sort(table(car_cells$cell_type), decreasing = TRUE))

# ── Funzione: calcola rapporto CD4/CD8 dato un dataframe CAR+ ─
calc_ratio <- function(df, label) {
  n_cd4   <- sum(df$lineage == "CD4+")
  n_cd8   <- sum(df$lineage == "CD8+")
  n_other <- sum(df$lineage == "Other")
  n_tot   <- nrow(df)

  # Solo cellule con lineage definito per il rapporto
  n_defined <- n_cd4 + n_cd8
  pct_cd4   <- if (n_defined > 0) round(100*n_cd4/n_defined, 1) else NA
  pct_cd8   <- if (n_defined > 0) round(100*n_cd8/n_defined, 1) else NA
  ratio_cd4_cd8 <- if (n_cd8 > 0) round(n_cd4/n_cd8, 3) else NA

  cat(sprintf(
    "\n[%s]\n  CAR+ totali:  %d\n  CD4+: %d (%.1f%%)  CD8+: %d (%.1f%%)  Other: %d\n  Rapporto CD4/CD8: %.3f\n",
    label, n_tot, n_cd4, pct_cd4, n_cd8, pct_cd8, n_other, ratio_cd4_cd8))

  data.frame(
    scenario       = label,
    n_CAR_totali   = n_tot,
    n_CD4          = n_cd4,
    n_CD8          = n_cd8,
    n_Other        = n_other,
    n_defined      = n_defined,
    pct_CD4        = pct_cd4,
    pct_CD8        = pct_cd8,
    ratio_CD4_CD8  = ratio_cd4_cd8,
    stringsAsFactors = FALSE
  )
}

# ── BASELINE: annotazione attuale ────────────────────────────
baseline <- calc_ratio(car_cells, "BASELINE\n(annotazione attuale)")

# ── SCENARIO A: escludi le CD8αα-like dall'analisi ───────────
# Le CD8αα-like vengono tolte dal dataset → non contano né come CD8 né come CD4
car_excl <- car_cells %>%
  filter(cd8_class != "CD8aa_like")
scenA <- calc_ratio(car_excl, "SCENARIO A\n(escludi CD8αα-like)")

# ── SCENARIO B: riannotale come MAIT/NKT/γδ ──────────────────
# Le CD8αα-like rimangono nel dataset ma il loro lineage diventa "Other"
# (non più CD8+, non CD4+) — quindi escono dal numeratore e denominatore
car_rean <- car_cells %>%
  mutate(
    lineage = if_else(cd8_class == "CD8aa_like", "Other", lineage)
  )
scenB <- calc_ratio(car_rean, "SCENARIO B\n(riannotate come MAIT/NKT/γδ)")

# ── Tabella comparativa ───────────────────────────────────────
ratio_table <- bind_rows(baseline, scenA, scenB)

cat("\n── TABELLA COMPARATIVA ──────────────────────────────────\n")
print(as.data.frame(ratio_table %>% select(scenario, n_CD4, n_CD8,
                                            pct_CD4, pct_CD8,
                                            ratio_CD4_CD8)))

# Delta rispetto al baseline
delta_A <- ratio_table$pct_CD4[2] - ratio_table$pct_CD4[1]
delta_B <- ratio_table$pct_CD4[3] - ratio_table$pct_CD4[1]
cat(sprintf(
  "\nDelta %% CD4 (Scenario A vs Baseline): %+.1f punti percentuali\n",
  delta_A))
cat(sprintf(
  "Delta %% CD4 (Scenario B vs Baseline): %+.1f punti percentuali\n",
  delta_B))

# ── Per campione (scenario B, il più informativo) ────────────
ratio_per_sample_B <- car_rean %>%
  group_by(sample) %>%
  summarise(
    n_CD4          = sum(lineage == "CD4+"),
    n_CD8          = sum(lineage == "CD8+"),
    n_Other        = sum(lineage == "Other"),
    .groups = "drop"
  ) %>%
  mutate(
    n_defined     = n_CD4 + n_CD8,
    pct_CD4       = round(100*n_CD4/pmax(n_defined,1), 1),
    pct_CD8       = round(100*n_CD8/pmax(n_defined,1), 1),
    ratio_CD4_CD8 = round(n_CD4/pmax(n_CD8,1), 3)
  )

ratio_per_sample_BL <- car_cells %>%
  group_by(sample) %>%
  summarise(
    n_CD4  = sum(lineage == "CD4+"),
    n_CD8  = sum(lineage == "CD8+"),
    n_Other= sum(lineage == "Other"),
    .groups = "drop"
  ) %>%
  mutate(
    n_defined     = n_CD4 + n_CD8,
    pct_CD4       = round(100*n_CD4/pmax(n_defined,1), 1),
    pct_CD8       = round(100*n_CD8/pmax(n_defined,1), 1),
    ratio_CD4_CD8 = round(n_CD4/pmax(n_CD8,1), 3)
  )

cat("\nRapporto per campione – BASELINE:\n")
print(as.data.frame(ratio_per_sample_BL))
cat("\nRapporto per campione – SCENARIO B:\n")
print(as.data.frame(ratio_per_sample_B))

# ============================================================
# 4. GRAFICI
# ============================================================
section("Grafici")

COL_MAIT  <- "#F3722C"
COL_NKT   <- "#277DA1"
COL_GDT   <- "#4D908E"
COL_CD161 <- "#F9C74F"
COL_NEC   <- "#AAAAAA"
COL_CD4   <- "#E63946"
COL_CD8   <- "#264653"
COL_OTH   <- "#2A9D8F"

PALETTE_ID <- c(
  "MAIT cell"                    = COL_MAIT,
  "iNKT cell"                    = COL_NKT,
  "γδ T cell"                    = COL_GDT,
  "CD161+ (MAIT/NKT ambiguo)"    = COL_CD161,
  "Innate-like NEC\n(CD8B dropout?)" = COL_NEC
)

# ── P16a: Identità CD8αα-like (barre aggregate) ──────────────
id_combined <- bind_rows(id_summary_all, id_summary_car) %>%
  mutate(identity = factor(identity, levels = names(PALETTE_ID)))

p16a <- ggplot(id_combined,
               aes(x = identity, y = pct, fill = identity)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_text(aes(label = paste0(pct, "%\n(n=", n, ")")),
            vjust = -0.3, size = 3.2) +
  facet_wrap(~ gruppo, nrow = 1) +
  scale_fill_manual(values = PALETTE_ID) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title    = "Identità delle cellule CD8αα-like",
       subtitle = "Sinistra: tutte le CD8αα-like  |  Destra: solo quelle CAR+",
       x = NULL, y = "% cellule CD8αα-like") +
  theme_classic(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
        strip.text    = element_text(face = "bold", size = 11),
        axis.text.x   = element_text(angle = 30, hjust = 1, size = 9))

# ── P16b: Per campione e CAR status ──────────────────────────
id_per_sample <- id_per_sample %>%
  mutate(identity = factor(identity, levels = names(PALETTE_ID)))

p16b <- ggplot(id_per_sample,
               aes(x = sample, y = pct, fill = identity)) +
  geom_col(position = "stack", width = 0.65) +
  geom_text(aes(label = ifelse(pct > 8, paste0(round(pct), "%"), "")),
            position = position_stack(vjust = 0.5),
            size = 2.8, color = "white", fontface = "bold") +
  facet_wrap(~ IS_CAR, nrow = 1) +
  scale_fill_manual(values = PALETTE_ID, name = NULL) +
  labs(title    = "Identità CD8αα-like per campione",
       subtitle = "CAR- = background biologico  |  CAR+ = quelle rilevanti",
       x = NULL, y = "% tra le CD8αα-like") +
  theme_classic(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
        strip.text    = element_text(face = "bold"),
        axis.text.x   = element_text(angle = 30, hjust = 1),
        legend.position = "bottom")

# ── P16c: Espressione marcatori nelle CD8αα-like per identità ─
marker_long <- expr_aa %>%
  select(identity, all_of(IDENTITY_GENES)) %>%
  pivot_longer(cols = all_of(IDENTITY_GENES),
               names_to = "gene", values_to = "expr") %>%
  filter(!is.na(expr)) %>%
  mutate(identity = factor(identity, levels = names(PALETTE_ID)),
         gene = factor(gene, levels = IDENTITY_GENES))

p16c <- ggplot(marker_long,
               aes(x = identity, y = expr, fill = identity)) +
  geom_violin(alpha = 0.75, scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.1, outlier.size = 0.2,
               fill = "white", alpha = 0.85) +
  facet_wrap(~ gene, nrow = 2, scales = "free_y") +
  scale_fill_manual(values = PALETTE_ID, guide = "none") +
  labs(title    = "Espressione marcatori identificativi per classe CD8αα-like",
       subtitle = "Conferma che la classificazione è consistente con i dati",
       x = NULL, y = "Espressione (log-norm)") +
  theme_classic(base_size = 10) +
  theme(plot.title    = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
        strip.text    = element_text(face = "bold"),
        axis.text.x   = element_text(angle = 35, hjust = 1, size = 8),
        legend.position = "none")

p16 <- (p16a / p16b / p16c) +
  plot_layout(heights = c(1.2, 1.2, 2)) +
  plot_annotation(
    title = "P16 – Identità delle cellule CD8αα-like: MAIT, NKT, γδ o altro?",
    theme = theme(plot.title = element_text(face = "bold",
                                            hjust = 0.5, size = 13))
  )
ggsave(file.path(out_dir, "P16_CD8alike_identity.png"),
       p16, width = 14, height = 18, dpi = 300, bg = "white")
cat("  P16 salvato\n")

# ── P17: Rapporto CD4/CD8 nei 3 scenari ──────────────────────
PALETTE_LIN <- c("CD4+" = COL_CD4, "CD8+" = COL_CD8, "Other" = COL_OTH)

# Dati in formato long per il barplot
ratio_long <- ratio_table %>%
  mutate(scenario = factor(scenario,
                           levels = c("BASELINE\n(annotazione attuale)",
                                      "SCENARIO A\n(escludi CD8αα-like)",
                                      "SCENARIO B\n(riannotate come MAIT/NKT/γδ)"))) %>%
  select(scenario, n_CD4, n_CD8, n_Other) %>%
  pivot_longer(cols = c(n_CD4, n_CD8, n_Other),
               names_to = "lineage",
               values_to = "n") %>%
  mutate(
    lineage = recode(lineage,
                     n_CD4 = "CD4+", n_CD8 = "CD8+", n_Other = "Other"),
    lineage = factor(lineage, levels = c("CD4+","CD8+","Other"))
  ) %>%
  group_by(scenario) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  ungroup()

# Tabella annotazioni: %CD4, %CD8, ratio
annot_df <- ratio_table %>%
  mutate(
    scenario = factor(scenario,
                      levels = c("BASELINE\n(annotazione attuale)",
                                 "SCENARIO A\n(escludi CD8αα-like)",
                                 "SCENARIO B\n(riannotate come MAIT/NKT/γδ)")),
    label = sprintf("CD4: %.1f%%\nCD8: %.1f%%\nratio: %.2f",
                    pct_CD4, pct_CD8, ratio_CD4_CD8)
  )

p17a <- ggplot(ratio_long,
               aes(x = scenario, y = n, fill = lineage)) +
  geom_col(position = "stack", width = 0.6) +
  geom_text(aes(label = ifelse(pct > 3,
                               paste0(pct, "%"), "")),
            position = position_stack(vjust = 0.5),
            size = 3.5, color = "white", fontface = "bold") +
  scale_fill_manual(values = PALETTE_LIN, name = NULL) +
  labs(title    = "Numero CAR-T per lineage nei 3 scenari",
       x = NULL, y = "Numero cellule CAR+") +
  theme_classic(base_size = 12) +
  theme(plot.title    = element_text(face = "bold", hjust = 0.5),
        axis.text.x   = element_text(size = 10),
        legend.position = "bottom")

# Grafico % CD4 vs % CD8
pct_long <- ratio_table %>%
  mutate(scenario = factor(scenario,
                           levels = c("BASELINE\n(annotazione attuale)",
                                      "SCENARIO A\n(escludi CD8αα-like)",
                                      "SCENARIO B\n(riannotate come MAIT/NKT/γδ)"))) %>%
  select(scenario, pct_CD4, pct_CD8) %>%
  pivot_longer(cols = c(pct_CD4, pct_CD8),
               names_to = "lineage", values_to = "pct") %>%
  mutate(lineage = recode(lineage, pct_CD4 = "CD4+", pct_CD8 = "CD8+"))

p17b <- ggplot(pct_long,
               aes(x = scenario, y = pct,
                   color = lineage, group = lineage)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 4) +
  geom_text(aes(label = paste0(pct, "%")),
            vjust = -1, size = 3.5, fontface = "bold") +
  scale_color_manual(values = c("CD4+" = COL_CD4, "CD8+" = COL_CD8),
                     name = NULL) +
  scale_y_continuous(limits = c(0, 100),
                     labels = function(x) paste0(x, "%")) +
  labs(title    = "% CD4+ e CD8+ tra le CAR-T nei 3 scenari",
       subtitle = sprintf(
         "Delta %% CD4: Scenario A = %+.1f pp  |  Scenario B = %+.1f pp",
         delta_A, delta_B),
       x = NULL, y = "% cellule CAR+") +
  theme_classic(base_size = 12) +
  theme(plot.title    = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
        legend.position = "bottom")

# Per campione: baseline vs scenario B
per_sample_compare <- bind_rows(
  ratio_per_sample_BL %>% mutate(scenario = "BASELINE"),
  ratio_per_sample_B  %>% mutate(scenario = "SCENARIO B")
) %>%
  select(scenario, sample, pct_CD4, pct_CD8, ratio_CD4_CD8) %>%
  pivot_longer(cols = c(pct_CD4, pct_CD8),
               names_to = "lineage", values_to = "pct") %>%
  mutate(lineage = recode(lineage, pct_CD4 = "CD4+", pct_CD8 = "CD8+"))

p17c <- ggplot(per_sample_compare,
               aes(x = sample, y = pct,
                   fill = lineage, alpha = scenario)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_text(aes(label = paste0(round(pct), "%")),
            position = position_dodge(0.8),
            vjust = -0.4, size = 2.8, fontface = "bold") +
  scale_fill_manual(values = c("CD4+" = COL_CD4, "CD8+" = COL_CD8),
                    name = NULL) +
  scale_alpha_manual(values = c("BASELINE" = 0.45, "SCENARIO B" = 1.0),
                     name = NULL) +
  facet_wrap(~ lineage, nrow = 1) +
  labs(title    = "% CD4+ e CD8+ CAR-T per campione: BASELINE vs SCENARIO B",
       subtitle = "Opaco = Scenario B (riannotate)  |  Trasparente = Baseline",
       x = NULL, y = "% CAR-T") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_classic(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
        strip.text    = element_text(face = "bold"),
        axis.text.x   = element_text(angle = 30, hjust = 1),
        legend.position = "bottom")

p17 <- (p17a | p17b) / p17c +
  plot_annotation(
    title = "P17 – Rapporto CAR-T CD4/CD8: Baseline vs Scenario A vs Scenario B",
    theme = theme(plot.title = element_text(face = "bold",
                                            hjust = 0.5, size = 13))
  )
ggsave(file.path(out_dir, "P17_ratio_CD4_CD8_scenarios.png"),
       p17, width = 16, height = 14, dpi = 300, bg = "white")
cat("  P17 salvato\n")

# ============================================================
# 5. EXCEL
# ============================================================
section("Export Excel")

writexl::write_xlsx(
  list(
    "Identita_CD8alike_totali"  = as.data.frame(id_summary_all),
    "Identita_CD8alike_CAR"     = as.data.frame(id_summary_car),
    "Identita_per_campione"     = as.data.frame(id_per_sample),
    "Rapporto_3_scenari"        = as.data.frame(ratio_table %>%
                                    select(scenario, n_CD4, n_CD8,
                                           n_Other, pct_CD4, pct_CD8,
                                           ratio_CD4_CD8)),
    "Rapporto_per_campione_BL"  = as.data.frame(ratio_per_sample_BL),
    "Rapporto_per_campione_ScB" = as.data.frame(ratio_per_sample_B)
  ),
  path = file.path(out_dir, "CD8alike_identity_and_ratio.xlsx")
)
cat(sprintf("  Excel → %s\n",
            file.path(out_dir, "CD8alike_identity_and_ratio.xlsx")))

# ============================================================
# 6. RIEPILOGO FINALE
# ============================================================
section("RIEPILOGO FINALE")

cat(sprintf("
─── Identità CD8αα-like (CAR+) ──────────────────────────────
"))
print(as.data.frame(id_summary_car %>% select(identity, n, pct, gruppo)))

cat(sprintf("
─── Rapporto CD4/CD8 CAR-T ──────────────────────────────────
  BASELINE   : CD4= %.1f%%  CD8= %.1f%%  ratio= %.3f
  SCENARIO A : CD4= %.1f%%  CD8= %.1f%%  ratio= %.3f  (Δ CD4: %+.1f pp)
  SCENARIO B : CD4= %.1f%%  CD8= %.1f%%  ratio= %.3f  (Δ CD4: %+.1f pp)
─────────────────────────────────────────────────────────────
",
  ratio_table$pct_CD4[1], ratio_table$pct_CD8[1], ratio_table$ratio_CD4_CD8[1],
  ratio_table$pct_CD4[2], ratio_table$pct_CD8[2], ratio_table$ratio_CD4_CD8[2], delta_A,
  ratio_table$pct_CD4[3], ratio_table$pct_CD8[3], ratio_table$ratio_CD4_CD8[3], delta_B
))
