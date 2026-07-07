# ============================================================
#  9c – Confronto tre definizioni di cellule CAR-T
#
#  Tre definizioni a confronto per ogni campione:
#    DEF1 = CAR == "YES"          (riallineamento genomico diretto)
#    DEF2 = IS_CAR_ALLIN_scREP == "YES"  (clonotype VDJ, definizione Alessio)
#    DEF3 = DEF1 OR DEF2          (unione, la più completa)
#
#  Per ogni definizione: distribuzione CD4+ / CD8+ / Memory T / Other
#
#  Nota: per campioni AB l'annotazione cell_type è cluster-level
#  (PIPELINE_2 auto_annotate). Per campioni I è manuale (più affidabile).
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(writexl)
  library(scales)
})

# ── PERCORSI ──────────────────────────────────────────────────
base_dir <- path.expand("~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/CART/Code/Seurat_analysis/")
rds_path <- file.path(base_dir, "2_annotation", "all_samples_annotated_COMPLETE.rds")
out_dir  <- file.path(base_dir, "9_CAR_final", "res")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

section <- function(title)
  cat(paste0("\n", strrep("=", 65), "\n  ", title,
             "\n", strrep("=", 65), "\n"))

# ── Classificazione CD4/CD8 ────────────────────────────────────
cd4_types <- c("Naive CD4+ T cells","Th1 cells","Th2 cells","Th17 cells",
               "Tfh cells","Effector CD4+ T cells","Tregs",
               "Proliferating CD4+ T cells")
cd8_types <- c("Cytotoxic CD8+ T cells","Naive CD8+ T cells",
               "Proliferating CD8+ T cells")

group_label <- function(lbl) {
  if (is.na(lbl))          return("Unknown")
  if (lbl %in% cd4_types)  return("CD4+")
  if (lbl %in% cd8_types)  return("CD8+")
  if (lbl == "Memory T cells") return("Memory T")
  return("Other")
}

# ============================================================
# 1. CARICAMENTO
# ============================================================
section("Caricamento dati")
cat(sprintf("rds esiste: %s\n", file.exists(rds_path)))
all_samples <- readRDS(rds_path)

ab_names <- c("Ca_blood_AB","Ca_bone_AB","Bo_blood_AB","Bo_bone_AB","Me_bone_AB")
i_names  <- c("Bo_bone_I","Ca_bone_I","Me_bone_I")
all_names <- c(i_names, ab_names)

# ============================================================
# 2. ANALISI PER CAMPIONE E DEFINIZIONE
# ============================================================
section("Distribuzione CAR-T per definizione")

PALETTE_GROUP <- c(
  "CD4+"     = "#E63946",
  "CD8+"     = "#264653",
  "Memory T" = "#2A9D8F",
  "Other"    = "#AAAAAA",
  "Unknown"  = "#DDDDDD"
)

results <- list()

for (nm in all_names) {
  obj  <- all_samples[[nm]]
  if (is.null(obj)) { cat(sprintf("  [SKIP] %s\n", nm)); next }
  meta <- obj@meta.data

  # Vettori logici per le tre definizioni
  def1 <- !is.na(meta$CAR) & meta$CAR == "YES"
  def2 <- !is.na(meta$IS_CAR_ALLIN_scREP) & meta$IS_CAR_ALLIN_scREP == "YES"
  def3 <- def1 | def2

  cat(sprintf("\n%s | DEF1(CAR): %d | DEF2(ALLIN): %d | DEF3(union): %d\n",
              nm, sum(def1), sum(def2), sum(def3)))

  for (def_name in c("DEF1_CAR_direct","DEF2_IS_CAR_ALLIN","DEF3_union")) {
    mask <- switch(def_name,
      DEF1_CAR_direct  = def1,
      DEF2_IS_CAR_ALLIN = def2,
      DEF3_union        = def3
    )
    cells <- meta[mask, ]
    n     <- nrow(cells)
    if (n == 0) {
      results[[paste0(nm,"_",def_name)]] <- data.frame(
        sample=nm, definition=def_name, cell_type="(nessuna)", group="Other",
        n_CAR=0, stringsAsFactors=FALSE)
      next
    }

    tbl <- cells %>%
      group_by(cell_type) %>%
      summarise(n_CAR = n(), .groups = "drop") %>%
      mutate(
        sample     = nm,
        definition = def_name,
        group      = vapply(cell_type, group_label, character(1L))
      ) %>%
      arrange(desc(n_CAR))

    results[[paste0(nm,"_",def_name)]] <- tbl

    # Stampa sintetica
    g <- tbl %>% group_by(group) %>% summarise(n=sum(n_CAR),.groups="drop")
    n_cd4 <- sum(g$n[g$group=="CD4+"])
    n_cd8 <- sum(g$n[g$group=="CD8+"])
    n_mem <- sum(g$n[g$group=="Memory T"])
    n_oth <- sum(g$n[g$group=="Other"])
    tot_t <- n_cd4 + n_cd8 + n_mem
    cat(sprintf("  %-22s | n=%4d | CD4+: %3d (%4.1f%%) | CD8+: %3d (%4.1f%%) | Mem: %2d | Other: %2d\n",
                def_name, n, n_cd4, 100*n_cd4/max(tot_t,1),
                n_cd8, 100*n_cd8/max(tot_t,1), n_mem, n_oth))
  }
}

# ============================================================
# 3. TABELLE AGGREGATE
# ============================================================
section("Aggregazione")

all_res <- bind_rows(results)

# Sintesi per definizione (tutti campioni)
summary_def <- all_res %>%
  group_by(definition, group) %>%
  summarise(n = sum(n_CAR), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = n, values_fill = 0) %>%
  mutate(
    Total_T = rowSums(across(any_of(c("CD4+","CD8+","Memory T")))),
    pct_CD4 = round(100 * `CD4+` / pmax(Total_T,1), 1),
    pct_CD8 = round(100 * `CD8+` / pmax(Total_T,1), 1)
  )

cat("\n── Tutti i campioni aggregati ──\n")
print(as.data.frame(summary_def))

# Sintesi per definizione separata I vs AB
for (group_type in list(i=i_names, ab=ab_names)) {
  label <- if (identical(group_type, i_names)) "Campioni I" else "Campioni AB"
  cat(sprintf("\n── %s ──\n", label))
  s <- all_res %>%
    filter(sample %in% group_type) %>%
    group_by(definition, group) %>%
    summarise(n = sum(n_CAR), .groups = "drop") %>%
    pivot_wider(names_from = group, values_from = n, values_fill = 0) %>%
    mutate(
      Total_T = rowSums(across(any_of(c("CD4+","CD8+","Memory T")))),
      pct_CD4 = round(100 * `CD4+` / pmax(Total_T,1), 1),
      pct_CD8 = round(100 * `CD8+` / pmax(Total_T,1), 1)
    )
  print(as.data.frame(s))
}

# Dettaglio per campione (DEF3 union, la più completa)
cat("\n── Dettaglio per campione – DEF3 (unione) ──\n")
detail_def3 <- all_res %>%
  filter(definition == "DEF3_union") %>%
  group_by(sample, group) %>%
  summarise(n = sum(n_CAR), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = n, values_fill = 0) %>%
  mutate(
    Total_T = rowSums(across(any_of(c("CD4+","CD8+","Memory T")))),
    pct_CD4 = round(100 * `CD4+` / pmax(Total_T,1), 1),
    pct_CD8 = round(100 * `CD8+` / pmax(Total_T,1), 1)
  )
print(as.data.frame(detail_def3))

# Export Excel
writexl::write_xlsx(
  list(
    "Sintesi_per_definizione"  = as.data.frame(summary_def),
    "I_samples_per_def"        = as.data.frame(
      all_res %>% filter(sample %in% i_names) %>%
        group_by(definition,group) %>%
        summarise(n=sum(n_CAR),.groups="drop") %>%
        pivot_wider(names_from=group, values_from=n, values_fill=0)),
    "AB_samples_per_def"       = as.data.frame(
      all_res %>% filter(sample %in% ab_names) %>%
        group_by(definition,group) %>%
        summarise(n=sum(n_CAR),.groups="drop") %>%
        pivot_wider(names_from=group, values_from=n, values_fill=0)),
    "Per_campione_DEF3_union"  = as.data.frame(detail_def3),
    "Dettaglio_celltype"       = as.data.frame(
      all_res %>% filter(definition=="DEF3_union") %>%
        arrange(sample, desc(n_CAR)) %>%
        select(sample, cell_type, group, n_CAR))
  ),
  path = file.path(out_dir, "CAR_three_definitions_comparison.xlsx")
)
cat(sprintf("\nExcel → %s\n",
            file.path(out_dir, "CAR_three_definitions_comparison.xlsx")))

# ============================================================
# 4. GRAFICI
# ============================================================
section("Grafici")

# ── P10: Confronto 3 definizioni – tutti campioni ─────────────
p10_data <- all_res %>%
  group_by(definition, group) %>%
  summarise(n = sum(n_CAR), .groups = "drop") %>%
  filter(group %in% c("CD4+","CD8+","Memory T")) %>%
  mutate(definition = recode(definition,
    DEF1_CAR_direct   = "DEF1\nCAR diretto\n(genomico)",
    DEF2_IS_CAR_ALLIN = "DEF2\nIS_CAR_ALLIN\n(clonotype VDJ)",
    DEF3_union        = "DEF3\nUnione\n(DEF1 + DEF2)"))

p10 <- ggplot(p10_data, aes(x = definition, y = n, fill = group)) +
  geom_col(position = "stack", width = 0.6) +
  scale_fill_manual(values = PALETTE_GROUP) +
  labs(
    title    = "Distribuzione CAR-T per definizione – tutti i campioni",
    x = NULL, y = "Numero cellule", fill = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 10))

ggsave(file.path(out_dir, "P10_three_definitions_all.png"),
       p10, width = 10, height = 6, dpi = 300, bg = "white")
cat("  P10 salvato\n")

# ── P11: I campioni vs AB per definizione ──────────────────────
p11_data <- all_res %>%
  mutate(type = ifelse(sample %in% i_names, "I (infusione)", "AB (post-infusione)")) %>%
  group_by(type, definition, group) %>%
  summarise(n = sum(n_CAR), .groups = "drop") %>%
  filter(group %in% c("CD4+","CD8+","Memory T")) %>%
  mutate(definition = recode(definition,
    DEF1_CAR_direct   = "DEF1\n(genomico)",
    DEF2_IS_CAR_ALLIN = "DEF2\n(clonotype)",
    DEF3_union        = "DEF3\n(unione)"))

p11 <- ggplot(p11_data, aes(x = definition, y = n, fill = group)) +
  geom_col(position = "stack", width = 0.6) +
  facet_wrap(~ type, scales = "free_y") +
  scale_fill_manual(values = PALETTE_GROUP) +
  labs(
    title = "Distribuzione CAR-T: I vs AB per definizione",
    x = NULL, y = "Numero cellule", fill = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(size = 9))

ggsave(file.path(out_dir, "P11_I_vs_AB_three_definitions.png"),
       p11, width = 12, height = 6, dpi = 300, bg = "white")
cat("  P11 salvato\n")

# ── P12: % CD4 vs CD8 per definizione e tipo campione ─────────
p12_data <- all_res %>%
  mutate(type = ifelse(sample %in% i_names, "I", "AB")) %>%
  group_by(type, definition, group) %>%
  summarise(n = sum(n_CAR), .groups = "drop") %>%
  filter(group %in% c("CD4+","CD8+")) %>%
  group_by(type, definition) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  ungroup() %>%
  mutate(definition = recode(definition,
    DEF1_CAR_direct   = "DEF1 (genomico)",
    DEF2_IS_CAR_ALLIN = "DEF2 (clonotype)",
    DEF3_union        = "DEF3 (unione)"))

p12 <- ggplot(p12_data, aes(x = definition, y = pct, fill = group)) +
  geom_col(position = "stack", width = 0.6) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "gray60") +
  facet_wrap(~ type) +
  scale_fill_manual(values = PALETTE_GROUP) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(
    title    = "% CD4+ vs CD8+ per definizione CAR-T",
    subtitle = "Linea tratteggiata = 50%",
    x = NULL, y = "% di T cells CAR", fill = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 20, hjust = 1))

ggsave(file.path(out_dir, "P12_pct_CD4_CD8_by_definition.png"),
       p12, width = 11, height = 6, dpi = 300, bg = "white")
cat("  P12 salvato\n")

# ── P13: DEF1 (CAR dirette) – dettaglio cell_type per campione ─
p13_data <- all_res %>%
  filter(definition == "DEF1_CAR_direct", n_CAR > 0,
         group %in% c("CD4+","CD8+","Memory T","Other")) %>%
  mutate(
    sample_type = ifelse(sample %in% i_names, "I", "AB"),
    cell_type = factor(cell_type,
      levels = rev(unique(cell_type[order(group, n_CAR)])))
  )

p13 <- ggplot(p13_data,
              aes(x = reorder(cell_type, n_CAR), y = n_CAR, fill = group)) +
  geom_col(width = 0.7) +
  coord_flip() +
  facet_wrap(~ sample, scales = "free", ncol = 3) +
  scale_fill_manual(values = PALETTE_GROUP) +
  labs(
    title = "DEF1 – CAR dirette (CAR==YES): distribuzione per cell_type e campione",
    x = NULL, y = "Numero cellule", fill = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        strip.text = element_text(face = "bold", size = 9))

ggsave(file.path(out_dir, "P13_DEF1_direct_CAR_celltype.png"),
       p13, width = 16, height = 10, dpi = 300, bg = "white")
cat("  P13 salvato\n")

# ============================================================
# 5. RIEPILOGO FINALE
# ============================================================
section("RIEPILOGO FINALE")

cat("
╔══════════════════════════════════════════════════════════════╗
║  CONFRONTO TRE DEFINIZIONI – TUTTI I CAMPIONI               ║
╚══════════════════════════════════════════════════════════════╝
")
print(as.data.frame(summary_def))

cat("
─────────────────────────────────────────────────────────────
INTERPRETAZIONE:
  DEF1 (CAR==YES, genomico):
    Campioni I: alta copertura (1103 cellule totali), prevalenza CD4+
    Campioni AB: bassa copertura (74 cellule) – CAR expression
                 più bassa dopo infusione → meno reads mappanti

  DEF2 (IS_CAR_ALLIN, clonotype VDJ):
    Campioni I: SOTTOSTIMA rispetto a DEF1 perché perde le cellule
                CAR+ senza dati TCR (Bo: 779→142)
    Campioni AB: alta copertura via propagazione clonale, ma
                 biased CD8+ perché i seed clonotipi catturati
                 sono prevalentemente CD8+

  DEF3 (unione DEF1+DEF2):
    La più completa: recupera le CD4+ perse dalla DEF2 nei
    campioni I e mantiene la copertura clonale nei campioni AB
─────────────────────────────────────────────────────────────
")

cat(sprintf("\nOutput in: %s\n", out_dir))
for (f in c("CAR_three_definitions_comparison.xlsx",
            "P10_three_definitions_all.png",
            "P11_I_vs_AB_three_definitions.png",
            "P12_pct_CD4_CD8_by_definition.png",
            "P13_DEF1_direct_CAR_celltype.png"))
  cat(sprintf("  - %s\n", f))
