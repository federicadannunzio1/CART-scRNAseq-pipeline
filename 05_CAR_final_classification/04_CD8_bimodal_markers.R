# ============================================================
#  9d – Analisi bimodalità CD8: CD8αβ vs CD8αα-like
#
#  RAZIONALE (da indicazione del prof):
#    Le cellule T CD8+ VERE esprimono il co-recettore come
#    ETERODIMERO CD8α/CD8β → entrambe le subunità CD8A e CD8B
#    devono essere alte.
#    Le cellule innate-like (MAIT, NKT, IEL-like) esprimono
#    solo CD8αα (OMODIMERO di CD8A) → CD8B è assente o molto bassa.
#
#    Conseguenza attesa: nella distribuzione di CD8B sulle
#    cellule annotate come CD8+ ci sarà BIMODALITÀ:
#      - picco a ~0  → CD8αα-like (MAIT, NKT, innate-like)
#      - picco a >0  → CD8αβ+ vere (citotossiche convenzionali)
#
#  FIRME GENICHE (indicate dal prof):
#    CD4 resting/naive: LTB, MAL, IL32, IL7R
#    CD8 vere:          CD8B, CTSW (Cathepsin W), GZMK
#
#  STRATEGIA:
#    1. Su tutti i campioni AB: estrai tutte le cellule T
#       (sia CAR+ che non-CAR) per avere una base ampia
#    2. Plotta le distribuzioni di densità di:
#       a. CD8B nelle cellule annotate CD8+ → bimodalità?
#       b. CD8A vs CD8B scatter (2D) → cloud CD8B-low?
#       c. CTSW e GZMK nelle CD8+ → confermano la bimodalità?
#       d. Firma CD4 (LTB+MAL+IL32+IL7R) in tutti i T cells
#          → alcune "CD8+" hanno firma CD4?
#    3. Classifica ogni cellula CD8+ annotata come:
#       - CD8_TRUE   : CD8A > 0 AND CD8B > 0
#       - CD8_ATYP   : CD8A > 0 AND CD8B == 0  (probabile CD8αα)
#       - CD8_LOW    : entrambi ≈ 0 (dropout)
#    4. Separa la visualizzazione CAR+ vs non-CAR per capire
#       se il bias è specifico delle CAR-T
#
#  OUTPUT: 9_CAR_final/res/
#    P10_CD8B_density_bimodal.png
#    P11_CD8A_vs_CD8B_scatter.png
#    P12_CD8_signature_CTSW_GZMK.png
#    P13_CD4_signature_in_Tcells.png
#    P14_classification_CD8true_vs_atyp.png
#    P15_UMAP_key_markers.png
#    CD8_bimodal_markers.xlsx
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

# ── Tipi cellulari T (per includere tutto il compartimento T) ─
CD8_TYPES <- c("Cytotoxic CD8+ T cells",
               "Naive CD8+ T cells",
               "Proliferating CD8+ T cells")

CD4_TYPES <- c("Naive CD4+ T cells", "Th1 cells", "Th2 cells",
               "Th17 cells", "Tfh cells", "Effector CD4+ T cells",
               "Tregs", "Proliferating CD4+ T cells")

ALL_T_TYPES <- c(CD8_TYPES, CD4_TYPES, "Memory T cells",
                 "NKT cells", "MAIT cells", "gamma-delta T cells")

# Geni da analizzare
GENES_LINEAGE <- c("CD4", "CD8A", "CD8B")
GENES_CD8_SIG <- c("CD8B", "CTSW", "GZMK")          # firma CD8 vera
GENES_CD4_SIG <- c("LTB", "MAL", "IL32", "IL7R")    # firma CD4 resting

ALL_GENES <- unique(c(GENES_LINEAGE, GENES_CD8_SIG, GENES_CD4_SIG))

AB_NAMES <- c("Ca_blood_AB", "Ca_bone_AB",
              "Bo_blood_AB", "Bo_bone_AB", "Me_bone_AB")

# ============================================================
# 1. CARICAMENTO
# ============================================================
section("Caricamento dati")

cat(sprintf("RDS esiste: %s\n", file.exists(rds_path)))
all_samples <- readRDS(rds_path)

# ============================================================
# 2. ESTRAZIONE ESPRESSIONE GENICA PER TUTTE LE CELLULE T (AB)
# ============================================================
section("Estrazione espressione genica (cellule T, campioni AB)")

expr_list <- list()

for (nm in AB_NAMES) {
  obj <- all_samples[[nm]]
  if (is.null(obj)) { cat(sprintf("  [SKIP] %s mancante\n", nm)); next }

  # JoinLayers Seurat v5
  if (length(grep("^counts\\.", Layers(obj), value = TRUE)) > 0)
    obj <- JoinLayers(obj)

  meta <- obj@meta.data

  # Seleziona cellule T
  mask_T <- !is.na(meta$cell_type) & meta$cell_type %in% ALL_T_TYPES
  bc_T   <- rownames(meta)[mask_T]
  if (length(bc_T) == 0) next

  # Geni disponibili
  genes_ok <- ALL_GENES[ALL_GENES %in% rownames(obj)]

  # FetchData (log-normalizzata, layer "data")
  expr_df <- FetchData(obj, vars = genes_ok, cells = bc_T,
                       layer = "data")

  # Aggiungi colonne mancanti come 0
  for (g in setdiff(ALL_GENES, colnames(expr_df)))
    expr_df[[g]] <- 0L

  # Aggiungi metadati utili
  expr_df$sample           <- nm
  expr_df$cell_type         <- meta[bc_T, "cell_type"]
  expr_df$IS_CAR            <- ifelse(
    !is.na(meta[bc_T, "IS_CAR_ALLIN_scREP"]) &
    meta[bc_T, "IS_CAR_ALLIN_scREP"] == "YES",
    "CAR+", "CAR-"
  )
  expr_df$barcode           <- rownames(expr_df)

  # Gruppo CD8 / CD4 / Altro
  expr_df$lineage_group <- dplyr::case_when(
    expr_df$cell_type %in% CD8_TYPES ~ "CD8+",
    expr_df$cell_type %in% CD4_TYPES ~ "CD4+",
    TRUE                              ~ "Altro T"
  )

  n_car <- sum(expr_df$IS_CAR == "CAR+")
  cat(sprintf("  %-15s | T cells: %4d | CAR+: %3d | geni ok: %s\n",
              nm, nrow(expr_df), n_car,
              paste(genes_ok, collapse=", ")))

  expr_list[[nm]] <- expr_df
}

expr_all <- bind_rows(expr_list)
cat(sprintf("\nTotale cellule T nei campioni AB: %d\n", nrow(expr_all)))
cat(sprintf("di cui CAR+: %d\n", sum(expr_all$IS_CAR == "CAR+")))

# Sottoinsiemi utili
expr_cd8     <- expr_all %>% filter(lineage_group == "CD8+")
expr_cd8_car <- expr_cd8 %>% filter(IS_CAR == "CAR+")
cat(sprintf("\nCD8+ totali: %d | CD8+ CAR+: %d\n",
            nrow(expr_cd8), nrow(expr_cd8_car)))

# ============================================================
# 3. CLASSIFICAZIONE CD8+ VERE vs ATIPICHE
#
#    La distinzione CD8αβ vs CD8αα si basa su CD8B:
#      CD8B > 0  → CD8αβ+ (eterodimero, cellule convenzionali)
#      CD8B == 0 → CD8αα-like (omodimero, innate-like: MAIT/NKT)
#
#    Nota: nella realtà dello scRNA-seq ci sono molti dropout.
#    Usiamo soglie conservative:
#      CD8B_HIGH: CD8B > 0  (qualsiasi espressione rilevata)
#      CD8B_ZERO: CD8B == 0 (nessuna molecola rilevata = dropout o vera assenza)
# ============================================================
section("Classificazione CD8αβ vs CD8αα-like")

expr_cd8 <- expr_cd8 %>%
  mutate(
    CD8_class = case_when(
      CD8A > 0 & CD8B > 0  ~ "CD8αβ+\n(CD8A+, CD8B+)",
      CD8A > 0 & CD8B == 0 ~ "CD8αα-like\n(CD8A+, CD8B=0)",
      CD8A == 0 & CD8B > 0 ~ "Solo CD8B\n(CD8A=0, CD8B+)",  # raro
      TRUE                  ~ "CD8A=0, CD8B=0\n(dropout)"
    )
  )

# Distribuzione
class_summary <- expr_cd8 %>%
  count(CD8_class) %>%
  mutate(pct = round(100 * n / sum(n), 1))
cat("\nClassificazione CD8+ (tutti):\n")
print(as.data.frame(class_summary))

class_summary_car <- expr_cd8_car %>%
  mutate(
    CD8_class = case_when(
      CD8A > 0 & CD8B > 0  ~ "CD8αβ+\n(CD8A+, CD8B+)",
      CD8A > 0 & CD8B == 0 ~ "CD8αα-like\n(CD8A+, CD8B=0)",
      CD8A == 0 & CD8B > 0 ~ "Solo CD8B\n(CD8A=0, CD8B+)",
      TRUE                  ~ "CD8A=0, CD8B=0\n(dropout)"
    )
  ) %>%
  count(CD8_class) %>%
  mutate(pct = round(100 * n / sum(n), 1))
cat("\nClassificazione CD8+ CAR+:\n")
print(as.data.frame(class_summary_car))

# ============================================================
# 4. GRAFICI
# ============================================================
section("Generazione grafici")

# Palette comune
COL_ALL  <- "#264653"   # CD8+ tutti  – teal scuro
COL_CAR  <- "#E63946"   # CD8+ CAR+  – rosso
COL_CD4  <- "#F4A261"   # CD4+        – arancione
COL_MEM  <- "#2A9D8F"   # Memory/Altro T

# ── helper densità sovrapposta ────────────────────────────────
dens_overlay <- function(df_all, df_car, gene,
                         title, xlab, binwidth = NULL) {

  vals_all <- df_all[[gene]]
  vals_car <- if (!is.null(df_car)) df_car[[gene]] else numeric(0)

  n_all <- sum(!is.na(vals_all))
  n_car <- sum(!is.na(vals_car))

  subtitle <- sprintf("CD8+ totali (blu): n=%d  |  CD8+ CAR+ (rosso): n=%d",
                      n_all, n_car)

  p <- ggplot(data.frame(val = vals_all[!is.na(vals_all)]),
              aes(x = val)) +
    geom_density(fill = COL_ALL, alpha = 0.35,
                 color = COL_ALL, linewidth = 0.9) +
    geom_rug(color = COL_ALL, alpha = 0.25,
             length = unit(0.03, "npc"))

  if (n_car > 1) {
    p <- p +
      geom_density(data = data.frame(val = vals_car[!is.na(vals_car)]),
                   aes(x = val),
                   fill = COL_CAR, alpha = 0.45,
                   color = COL_CAR, linewidth = 0.9) +
      geom_rug(data = data.frame(val = vals_car[!is.na(vals_car)]),
               aes(x = val),
               color = COL_CAR, alpha = 0.5,
               length = unit(0.03, "npc"))
  }

  p +
    labs(title = title, subtitle = subtitle,
         x = xlab, y = "Densità") +
    theme_classic(base_size = 12) +
    theme(plot.title    = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, color = "gray40"))
}

# ── P10: Densità CD8B nelle cellule CD8+ (BIMODALITÀ) ─────────
# Questo è il plot principale: se CD8B è bimodale → due popolazioni
p10a <- dens_overlay(expr_cd8, expr_cd8_car, "CD8B",
  title = "Distribuzione CD8B nelle cellule annotate CD8+",
  xlab  = "CD8B espressione (log-norm)")

p10b <- dens_overlay(expr_cd8, expr_cd8_car, "CD8A",
  title = "Distribuzione CD8A nelle cellule annotate CD8+",
  xlab  = "CD8A espressione (log-norm)")

# Annotazione percentuali CD8αα-like
pct_atyp_all <- class_summary %>%
  filter(grepl("CD8αα", CD8_class)) %>% pull(pct)
pct_atyp_car <- class_summary_car %>%
  filter(grepl("CD8αα", CD8_class)) %>% pull(pct)
if (length(pct_atyp_all) == 0) pct_atyp_all <- 0
if (length(pct_atyp_car) == 0) pct_atyp_car <- 0

p10c <- ggplot(class_summary,
               aes(x = reorder(CD8_class, -n), y = pct,
                   fill = CD8_class)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = paste0(pct, "%\n(n=", n, ")")),
            vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = c(
    "CD8αβ+\n(CD8A+, CD8B+)"       = COL_ALL,
    "CD8αα-like\n(CD8A+, CD8B=0)"  = COL_CAR,
    "Solo CD8B\n(CD8A=0, CD8B+)"   = "#AAAAAA",
    "CD8A=0, CD8B=0\n(dropout)"    = "#CCCCCC"
  )) +
  labs(title    = "Classificazione CD8+ (tutti)",
       subtitle = "CD8αα-like = CD8A rilevato, CD8B assente",
       x = NULL, y = "% cellule") +
  ylim(0, max(class_summary$pct) * 1.15) +
  theme_classic(base_size = 12) +
  theme(plot.title    = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
        axis.text.x   = element_text(size = 9))

p10d <- ggplot(class_summary_car,
               aes(x = reorder(CD8_class, -n), y = pct,
                   fill = CD8_class)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = paste0(pct, "%\n(n=", n, ")")),
            vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = c(
    "CD8αβ+\n(CD8A+, CD8B+)"       = COL_ALL,
    "CD8αα-like\n(CD8A+, CD8B=0)"  = COL_CAR,
    "Solo CD8B\n(CD8A=0, CD8B+)"   = "#AAAAAA",
    "CD8A=0, CD8B=0\n(dropout)"    = "#CCCCCC"
  )) +
  labs(title    = "Classificazione CD8+ CAR+",
       subtitle = sprintf("CD8αα-like nelle CAR: %.1f%%", pct_atyp_car),
       x = NULL, y = "% cellule") +
  ylim(0, max(class_summary_car$pct) * 1.15) +
  theme_classic(base_size = 12) +
  theme(plot.title    = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
        axis.text.x   = element_text(size = 9))

p10 <- (p10a | p10b) / (p10c | p10d) +
  plot_annotation(
    title = "P10 – CD8B bimodalità: CD8αβ+ (vere) vs CD8αα-like (innate-like)",
    subtitle = paste0(
      "Le cellule CD8αα-like (CD8A+, CD8B=0) sono MAIT/NKT/IEL-like, non CD8 convenzionali\n",
      sprintf("CD8αα-like: %.1f%% di tutte le CD8+  |  %.1f%% delle CD8+ CAR+",
              pct_atyp_all, pct_atyp_car)
    ),
    theme = theme(
      plot.title    = element_text(face = "bold", hjust = 0.5, size = 13),
      plot.subtitle = element_text(hjust = 0.5, color = "gray40", size = 10)
    )
  )
ggsave(file.path(out_dir, "P10_CD8B_density_bimodal.png"),
       p10, width = 14, height = 10, dpi = 300, bg = "white")
cat("  P10 salvato\n")

# ── P11: Scatter CD8A vs CD8B ─────────────────────────────────
# Aspettativa:
#   CD8αβ+:      nuvola in alto a destra (CD8A alto, CD8B alto)
#   CD8αα-like:  nuvola sull'asse X (CD8A alto, CD8B = 0)

make_scatter_cd8 <- function(df, title_str, col_pt = COL_ALL,
                              alpha_pt = 0.4) {
  # Calcola % per ogni quadrante
  n_tot <- nrow(df)
  pct_hetero <- round(100 * mean(df$CD8A > 0 & df$CD8B > 0), 1)
  pct_homo   <- round(100 * mean(df$CD8A > 0 & df$CD8B == 0), 1)
  pct_none   <- round(100 * mean(df$CD8A == 0 & df$CD8B == 0), 1)

  ggplot(df, aes(x = CD8A, y = CD8B)) +
    geom_point(color = col_pt, alpha = alpha_pt, size = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed",
               color = "firebrick", linewidth = 0.6) +
    annotate("text", x = max(df$CD8A) * 0.7,
             y = max(df$CD8B) * 0.85,
             label = sprintf("CD8αβ+: %.1f%%", pct_hetero),
             color = COL_ALL, size = 3.5, fontface = "bold") +
    annotate("text", x = max(df$CD8A) * 0.7,
             y = -0.12,
             label = sprintf("CD8αα-like: %.1f%%", pct_homo),
             color = COL_CAR, size = 3.5, fontface = "bold") +
    labs(title    = title_str,
         subtitle = sprintf("n = %d | La linea tratteggiata = CD8B = 0", n_tot),
         x = "CD8A espressione (log-norm)",
         y = "CD8B espressione (log-norm)") +
    theme_classic(base_size = 12) +
    theme(plot.title    = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, color = "gray40"))
}

p11a <- make_scatter_cd8(expr_cd8, "CD8A vs CD8B – tutte le CD8+ (AB)", COL_ALL)
p11b <- make_scatter_cd8(expr_cd8_car,
                          "CD8A vs CD8B – solo CD8+ CAR+", COL_CAR)

p11 <- (p11a | p11b) +
  plot_annotation(
    title = "P11 – Scatter CD8A vs CD8B: identificazione CD8αα-like",
    subtitle = paste0(
      "Le cellule sulla linea orizzontale (CD8B = 0) sono CD8αα-like\n",
      "(MAIT cells, NKT cells, IEL-like – non CD8αβ convenzionali)"
    ),
    theme = theme(
      plot.title    = element_text(face = "bold", hjust = 0.5, size = 13),
      plot.subtitle = element_text(hjust = 0.5, color = "gray40", size = 10)
    )
  )
ggsave(file.path(out_dir, "P11_CD8A_vs_CD8B_scatter.png"),
       p11, width = 14, height = 6, dpi = 300, bg = "white")
cat("  P11 salvato\n")

# ── P12: Firma CD8 vera (CTSW e GZMK) ────────────────────────
# CTSW: Cathepsin W → quasi esclusiva delle CD8+ citotossiche vere
# GZMK: Granzyme K  → CD8+ effector memory
#
# Aspettativa:
#   Se CD8αβ+ e CD8αα-like coesistono, CTSW sarà bimodale:
#   - alta nelle CD8αβ+, bassa/assente nelle CD8αα-like

p12a <- dens_overlay(expr_cd8, expr_cd8_car, "CTSW",
  title = "CTSW (Cathepsin W) nelle CD8+",
  xlab  = "CTSW espressione (log-norm)")

p12b <- dens_overlay(expr_cd8, expr_cd8_car, "GZMK",
  title = "GZMK (Granzyme K) nelle CD8+",
  xlab  = "GZMK espressione (log-norm)")

# Confronto CTSW tra CD8αβ+ e CD8αα-like
expr_cd8_class <- expr_cd8 %>%
  mutate(
    CD8_class2 = ifelse(CD8B > 0, "CD8αβ+ (CD8B > 0)",
                                   "CD8αα-like (CD8B = 0)")
  )

p12c <- ggplot(expr_cd8_class, aes(x = CD8_class2, y = CTSW,
                                    fill = CD8_class2)) +
  geom_violin(alpha = 0.7, scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.1, outlier.size = 0.3, fill = "white",
               alpha = 0.9) +
  scale_fill_manual(values = c(
    "CD8αβ+ (CD8B > 0)"   = COL_ALL,
    "CD8αα-like (CD8B = 0)" = COL_CAR
  ), guide = "none") +
  labs(title    = "CTSW: CD8αβ+ vs CD8αα-like",
       subtitle = "CTSW è specifica per le CD8 vere",
       x = NULL, y = "CTSW (log-norm)") +
  theme_classic(base_size = 12) +
  theme(plot.title    = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, color = "gray40"))

p12d <- ggplot(expr_cd8_class, aes(x = CD8_class2, y = GZMK,
                                    fill = CD8_class2)) +
  geom_violin(alpha = 0.7, scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.1, outlier.size = 0.3, fill = "white",
               alpha = 0.9) +
  scale_fill_manual(values = c(
    "CD8αβ+ (CD8B > 0)"   = COL_ALL,
    "CD8αα-like (CD8B = 0)" = COL_CAR
  ), guide = "none") +
  labs(title    = "GZMK: CD8αβ+ vs CD8αα-like",
       subtitle = "GZMK è espressa nelle CD8 effector memory",
       x = NULL, y = "GZMK (log-norm)") +
  theme_classic(base_size = 12) +
  theme(plot.title    = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, color = "gray40"))

p12 <- (p12a | p12b) / (p12c | p12d) +
  plot_annotation(
    title = "P12 – Firma CD8 vera: CTSW e GZMK",
    subtitle = "CTSW bassa nelle CD8αα-like conferma che non sono CD8 convenzionali",
    theme = theme(
      plot.title    = element_text(face = "bold", hjust = 0.5, size = 13),
      plot.subtitle = element_text(hjust = 0.5, color = "gray40", size = 10)
    )
  )
ggsave(file.path(out_dir, "P12_CD8_signature_CTSW_GZMK.png"),
       p12, width = 14, height = 10, dpi = 300, bg = "white")
cat("  P12 salvato\n")

# ── P13: Firma CD4 resting (LTB, MAL, IL32, IL7R) ─────────────
# Queste 4 molecole caratterizzano le CD4+ naive/resting.
# Se delle cellule annotate come CD8+ le esprimono ad alti livelli,
# probabilmente sono CD4+ con dropout di CD4 (abbastanza comune).

# Module score firma CD4 su tutti i T cells
# (non possiamo usare AddModuleScore perché abbiamo già dati estratti,
#  usiamo la somma normalizzata dei geni disponibili)
cd4_sig_genes <- c("LTB", "MAL", "IL32", "IL7R")
cd4_sig_avail <- cd4_sig_genes[cd4_sig_genes %in% colnames(expr_all)]

cat(sprintf("\n  Geni firma CD4 disponibili: %s\n",
            paste(cd4_sig_avail, collapse=", ")))

if (length(cd4_sig_avail) >= 2) {
  expr_all <- expr_all %>%
    mutate(
      CD4_sig_score = rowMeans(
        across(all_of(cd4_sig_avail)), na.rm = TRUE
      )
    )
  expr_cd8 <- expr_cd8 %>%
    mutate(
      CD4_sig_score = rowMeans(
        across(all_of(cd4_sig_avail)), na.rm = TRUE
      )
    )
  expr_cd8_car <- expr_cd8_car %>%
    mutate(
      CD4_sig_score = rowMeans(
        across(all_of(cd4_sig_avail)), na.rm = TRUE
      )
    )

  # Distribuzione del CD4 score per tipo cellulare
  # (ordine biologico: CD4 types first, then CD8, then other)
  type_order <- c(
    "Naive CD4+ T cells", "Th1 cells", "Th2 cells", "Th17 cells",
    "Tfh cells", "Effector CD4+ T cells", "Tregs",
    "Proliferating CD4+ T cells", "Memory T cells",
    "Cytotoxic CD8+ T cells", "Naive CD8+ T cells",
    "Proliferating CD8+ T cells", "MAIT cells", "NKT cells",
    "gamma-delta T cells"
  )
  type_order <- type_order[type_order %in% unique(expr_all$cell_type)]

  expr_all_ord <- expr_all %>%
    filter(cell_type %in% type_order) %>%
    mutate(cell_type = factor(cell_type, levels = rev(type_order)),
           group_col = case_when(
             lineage_group == "CD4+"   ~ COL_CD4,
             lineage_group == "CD8+"   ~ COL_ALL,
             TRUE                       ~ COL_MEM
           ))

  # Mediana per cell_type per ordinare le barre
  median_order <- expr_all_ord %>%
    group_by(cell_type) %>%
    summarise(med = median(CD4_sig_score, na.rm = TRUE), .groups = "drop") %>%
    arrange(med) %>%
    pull(cell_type)
  expr_all_ord <- expr_all_ord %>%
    mutate(cell_type = factor(cell_type, levels = median_order))

  p13a <- ggplot(expr_all_ord,
                 aes(x = cell_type, y = CD4_sig_score,
                     fill = lineage_group)) +
    geom_violin(alpha = 0.7, scale = "width", trim = TRUE) +
    geom_boxplot(width = 0.12, outlier.size = 0.2,
                 fill = "white", alpha = 0.85) +
    coord_flip() +
    scale_fill_manual(
      values = c("CD4+" = COL_CD4, "CD8+" = COL_ALL,
                 "Altro T" = COL_MEM),
      name = NULL
    ) +
    geom_hline(yintercept = 0, linetype = "dashed",
               color = "gray50", linewidth = 0.5) +
    labs(title    = "Distribuzione firma CD4 (LTB + MAL + IL32 + IL7R)",
         subtitle = "Se CD8+ ha score alto → probabile CD4+ con dropout",
         x = NULL,
         y = "CD4 signature score (media log-norm)") +
    theme_classic(base_size = 11) +
    theme(plot.title    = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
          legend.position = "bottom")

  # Singoli geni (violin) per i tipi CD8+
  cd4_sig_long <- expr_cd8 %>%
    select(cell_type, IS_CAR, all_of(cd4_sig_avail)) %>%
    pivot_longer(cols = all_of(cd4_sig_avail),
                 names_to = "gene", values_to = "expr") %>%
    mutate(
      label = paste0(cell_type, "\n(", IS_CAR, ")")
    )

  p13b <- ggplot(cd4_sig_long, aes(x = cell_type, y = expr,
                                    fill = IS_CAR)) +
    geom_violin(alpha = 0.65, scale = "width", trim = TRUE,
                position = position_dodge(0.8)) +
    geom_boxplot(width = 0.1, outlier.size = 0.2,
                 position = position_dodge(0.8), fill = "white") +
    facet_wrap(~ gene, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = c("CAR+" = COL_CAR, "CAR-" = COL_ALL),
                      name = NULL) +
    labs(title    = "Geni CD4 resting nelle CD8+ (CAR+ vs CAR-)",
         subtitle = "Valori alti in CD8+ → possibili CD4+ mislabeled",
         x = NULL, y = "Espressione (log-norm)") +
    theme_classic(base_size = 10) +
    theme(plot.title    = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
          axis.text.x   = element_text(angle = 30, hjust = 1, size = 8),
          strip.text    = element_text(face = "bold"),
          legend.position = "bottom")

  p13 <- p13a / p13b +
    plot_annotation(
      title = "P13 – Firma CD4 resting (LTB, MAL, IL32, IL7R) tra tutti i T cells",
      theme = theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 13)
      )
    )

  ggsave(file.path(out_dir, "P13_CD4_signature_in_Tcells.png"),
         p13, width = 14, height = 12, dpi = 300, bg = "white")
  cat("  P13 salvato\n")
} else {
  cat("  [SKIP P13] Geni firma CD4 non disponibili a sufficienza\n")
}

# ── P14: Summary classificazione CAR+ per campione ──────────
section("P14 – Classificazione CD8 per campione")

# Aggiungi la classe CD8 all'intero dataframe CD8+
expr_cd8_full <- expr_all %>%
  filter(lineage_group == "CD8+") %>%
  mutate(
    CD8_class2 = ifelse(CD8B > 0, "CD8αβ+", "CD8αα-like")
  )

# Conta per campione e status CAR
class_per_sample <- expr_cd8_full %>%
  group_by(sample, IS_CAR, CD8_class2) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(sample, IS_CAR) %>%
  mutate(pct = round(100 * n / sum(n), 1))

p14 <- ggplot(class_per_sample,
              aes(x = sample, y = pct, fill = CD8_class2)) +
  geom_col(position = "stack", width = 0.6) +
  geom_text(aes(label = ifelse(pct > 5, paste0(pct, "%"), "")),
            position = position_stack(vjust = 0.5),
            size = 3, color = "white", fontface = "bold") +
  facet_wrap(~ IS_CAR, nrow = 1) +
  scale_fill_manual(
    values = c("CD8αβ+" = COL_ALL, "CD8αα-like" = COL_CAR),
    name = NULL
  ) +
  labs(
    title    = "P14 – % CD8αβ+ vs CD8αα-like per campione",
    subtitle = "Sinistra: CD8+ CAR-  |  Destra: CD8+ CAR+",
    x = NULL, y = "% cellule CD8+") +
  theme_classic(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
        strip.text    = element_text(face = "bold", size = 12),
        axis.text.x   = element_text(angle = 35, hjust = 1),
        legend.position = "bottom")

ggsave(file.path(out_dir, "P14_classification_CD8true_vs_atyp.png"),
       p14, width = 12, height = 6, dpi = 300, bg = "white")
cat("  P14 salvato\n")

# ── P15: FeaturePlot UMAP – CD8B, CTSW, GZMK per campione ────
section("P15 – UMAP marcatori chiave")

UMAP_GENES_FINAL <- c("CD8B", "CTSW", "GZMK", "LTB", "IL7R")

for (nm in AB_NAMES) {
  obj <- all_samples[[nm]]
  if (is.null(obj)) next
  if (length(grep("^counts\\.", Layers(obj), value = TRUE)) > 0)
    obj <- JoinLayers(obj)

  genes_ok <- UMAP_GENES_FINAL[UMAP_GENES_FINAL %in% rownames(obj)]
  if (length(genes_ok) == 0) next

  cat(sprintf("  [%s] FeaturePlot geni: %s\n",
              nm, paste(genes_ok, collapse=", ")))

  ncol_fp <- min(length(genes_ok), 5)

  p_fp <- FeaturePlot(
    obj,
    features   = genes_ok,
    reduction  = "umap",
    ncol       = ncol_fp,
    min.cutoff = "q05",
    max.cutoff = "q95",
    order      = TRUE,
    pt.size    = 0.5
  ) &
    theme_classic(base_size = 9) &
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          plot.title = element_text(size = 10, face = "bold", hjust = 0.5))

  p_dim <- DimPlot(
    obj, reduction = "umap", group.by = "cell_type",
    label = TRUE, label.size = 2.5, repel = TRUE, pt.size = 0.5
  ) +
    ggtitle(paste0(nm, " – cell types")) +
    theme_classic(base_size = 9) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          legend.position = "none",
          axis.text = element_blank(), axis.ticks = element_blank())

  p_combined <- p_dim / p_fp +
    plot_layout(heights = c(1, 1)) +
    plot_annotation(
      title = paste0(nm, " – CD8B, CTSW, GZMK, LTB, IL7R"),
      theme = theme(plot.title = element_text(face = "bold",
                                              hjust = 0.5, size = 12))
    )

  out_fp <- file.path(out_dir, paste0("P15_UMAP_", nm, "_markers.png"))
  ggsave(out_fp, p_combined,
         width  = 4 * ncol_fp,
         height = 10,
         dpi = 300, bg = "white")
  cat(sprintf("  → %s\n", basename(out_fp)))
}

# ============================================================
# 5. TABELLE RIASSUNTIVE + EXCEL
# ============================================================
section("Export Excel")

# Tabella 1: classificazione globale CD8+
tab1 <- bind_rows(
  class_summary     %>% mutate(gruppo = "CD8+ tutti AB"),
  class_summary_car %>% mutate(gruppo = "CD8+ CAR+")
) %>% select(gruppo, CD8_class, n, pct)

# Tabella 2: per campione
tab2 <- class_per_sample %>%
  rename(classe = CD8_class2, pct_cellule = pct) %>%
  arrange(sample, IS_CAR, classe)

# Tabella 3: media espressione geni chiave per classe CD8
tab3 <- expr_cd8_class %>%
  group_by(CD8_class2) %>%
  summarise(
    n           = n(),
    CD8A_mean   = round(mean(CD8A,  na.rm=TRUE), 4),
    CD8B_mean   = round(mean(CD8B,  na.rm=TRUE), 4),
    CTSW_mean   = round(mean(CTSW,  na.rm=TRUE), 4),
    GZMK_mean   = round(mean(GZMK,  na.rm=TRUE), 4),
    .groups = "drop"
  )

# Tabella 4: firma CD4 nelle CD8+
if ("CD4_sig_score" %in% colnames(expr_cd8)) {
  tab4 <- expr_cd8 %>%
    mutate(CD8_class2 = ifelse(CD8B > 0, "CD8αβ+", "CD8αα-like")) %>%
    group_by(IS_CAR, CD8_class2) %>%
    summarise(
      n               = n(),
      CD4_sig_mean    = round(mean(CD4_sig_score, na.rm=TRUE), 4),
      LTB_mean        = round(mean(LTB,   na.rm=TRUE), 4),
      MAL_mean        = round(mean(MAL,   na.rm=TRUE), 4),
      IL32_mean       = round(mean(IL32,  na.rm=TRUE), 4),
      IL7R_mean       = round(mean(IL7R,  na.rm=TRUE), 4),
      .groups = "drop"
    )
} else {
  tab4 <- data.frame(note = "Geni firma CD4 non disponibili")
}

out_excel <- file.path(out_dir, "CD8_bimodal_markers.xlsx")
writexl::write_xlsx(
  list(
    "Classificazione_globale"  = as.data.frame(tab1),
    "Per_campione"             = as.data.frame(tab2),
    "Espressione_per_classe"   = as.data.frame(tab3),
    "Firma_CD4_nelle_CD8"      = as.data.frame(tab4)
  ),
  path = out_excel
)
cat(sprintf("\nExcel → %s\n", out_excel))

# ============================================================
# 6. RIEPILOGO
# ============================================================
section("RIEPILOGO FINALE")

cat(sprintf(
  "\nCD8+ totali (AB): %d
  di cui CD8αβ+ (CD8A+, CD8B+):    %.1f%%
  di cui CD8αα-like (CD8A+, CD8B=0): %.1f%%
\nCD8+ CAR+ (AB): %d
  di cui CD8αβ+ (CD8A+, CD8B+):    %.1f%%
  di cui CD8αα-like (CD8A+, CD8B=0): %.1f%%
",
  nrow(expr_cd8),
  class_summary$pct[grepl("CD8αβ", class_summary$CD8_class)],
  pct_atyp_all,
  nrow(expr_cd8_car),
  class_summary_car$pct[grepl("CD8αβ", class_summary_car$CD8_class)],
  pct_atyp_car
))

cat(paste0(
  "\nINTERPRETAZIONE ATTESA:\n",
  "  P10: densità CD8B bimodale → picco a 0 (CD8αα-like) + picco positivo (CD8αβ+)\n",
  "  P11: scatter CD8A vs CD8B → cloud sull'asse X = cellule CD8αα-like\n",
  "  P12: CTSW bassa nelle CD8αα-like → non sono CD8 citotossiche vere\n",
  "  P13: se CD8+ (specie CAR+) hanno firma CD4 alta → alcune sono CD4 mislabeled\n",
  "  P14: % CD8αα-like per campione → quale campione ha più cellule atipiche?\n",
  "  P15: UMAP → dove si localizzano le cellule CD8B-low?\n",
  "\nOutput in: ", out_dir, "\n"
))
