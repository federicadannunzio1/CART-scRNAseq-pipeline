# ============================================================
#  Q1: Dinamica della composizione CD4+ e CD8+
#      nel tempo (I → AB), stratificata per stato CAR
#
#  Domanda biologica:
#    Come cambia la distribuzione dei sottotipi CD4 e CD8
#    tra il prodotto di infusione (I) e dopo l'infusione (AB)?
#    Il pattern differisce tra cellule CAR+ e CAR-?
#
#  Nota metodologica:
#    N=3 pazienti → nessun test statistico formale,
#    analisi descrittiva + coerenza tra pazienti.
#    Le proporzioni sono calcolate DENTRO ogni compartimento
#    (CD4 o CD8), non rispetto a tutte le cellule del campione.
#
#  Prerequisiti:
#    all_samples_annotated_COMPLETE_IS_CAR_REVISED.rds
#
#  Output in out_dir/:
#    Q1_<paziente>_CD4_composition_CARpos_vs_CARneg.png
#    Q1_<paziente>_CD8_composition_CARpos_vs_CARneg.png
#    Q1_ALL_patients_composition_heatmap.png
#    Q1_composition_summary.xlsx
# ============================================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(scales)
library(openxlsx)

# ── UNICO PUNTO DA MODIFICARE ────────────────────────────────
rds_path <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/CART/Code/Seurat_analysis/2_annotation/all_samples_annotated_COMPLETE_IS_CAR_REVISED.rds"
out_dir  <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/CART/Code/Seurat_analysis/10_CART_functional_analysis/Q1_composition/"
# ─────────────────────────────────────────────────────────────

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

section <- function(title)
  cat(paste0("\n", strrep("=", 65), "\n  ", title, "\n",
             strrep("=", 65), "\n"))

# ============================================================
# DEFINIZIONE COMPARTIMENTI CD4 e CD8
# ============================================================

# Tipi cellulari che appartengono al compartimento CD4+
CD4_TYPES <- c(
  "Naive CD4+ T cells",
  "Effector CD4+ T cells",
  "Th1 cells",
  "Th2 cells",
  "Th17 cells",
  "Tfh cells",
  "Memory T cells",   # In questa pipeline = CD8 con basso naive/cytotox
                      # ATTENZIONE: se presenti in campioni CD4-enriched
                      # potrebbero essere ambigui. Verificare manualmente.
  "Tregs",
  "Proliferating CD4+ T cells"
)

# Tipi cellulari che appartengono al compartimento CD8+
CD8_TYPES <- c(
  "Cytotoxic CD8+ T cells",
  "Naive CD8+ T cells",
  "Proliferating CD8+ T cells"
)

# Ordine biologico per i plot
CD4_ORDER <- c(
  "Naive CD4+ T cells",
  "Th1 cells", "Th2 cells", "Th17 cells", "Tfh cells",
  "Effector CD4+ T cells",
  "Memory T cells",
  "Tregs",
  "Proliferating CD4+ T cells"
)

CD8_ORDER <- c(
  "Naive CD8+ T cells",
  "Cytotoxic CD8+ T cells",
  "Proliferating CD8+ T cells"
)

PALETTE_CD4 <- c(
  "Naive CD4+ T cells"         = "#E63946",
  "Th1 cells"                  = "#C1121F",
  "Th2 cells"                  = "#FF99C8",
  "Th17 cells"                 = "#FB5607",
  "Tfh cells"                  = "#800F2F",
  "Effector CD4+ T cells"      = "#F4A261",
  "Memory T cells"             = "#2A9D8F",
  "Tregs"                      = "#E9C46A",
  "Proliferating CD4+ T cells" = "#023E8A"
)

PALETTE_CD8 <- c(
  "Naive CD8+ T cells"         = "#577590",
  "Cytotoxic CD8+ T cells"     = "#264653",
  "Proliferating CD8+ T cells" = "#6A0572"
)

# Mappa paziente → campioni
# Bo: prodotto = Bo_bone_I | post-infusione = Bo_blood_AB (sangue) + Bo_bone_AB (midollo)
# Ca: prodotto = Ca_bone_I | post-infusione = Ca_blood_AB (sangue) + Ca_bone_AB (midollo)
# Me: prodotto = Me_bone_I | post-infusione = Me_bone_AB (solo midollo)
PATIENT_MAP <- list(
  Bo = list(
    I    = "Bo_bone_I",
    AB   = c("Bo_blood_AB", "Bo_bone_AB"),
    bone = "Bo_bone_AB",
    blood = "Bo_blood_AB"
  ),
  Ca = list(
    I    = "Ca_bone_I",
    AB   = c("Ca_blood_AB", "Ca_bone_AB"),
    bone = "Ca_bone_AB",
    blood = "Ca_blood_AB"
  ),
  Me = list(
    I    = "Me_bone_I",
    AB   = "Me_bone_AB",
    bone = "Me_bone_AB",
    blood = NULL
  )
)

# ============================================================
# CARICAMENTO
# ============================================================
section("Caricamento dati")

all_samples <- readRDS(rds_path)
cat("Campioni caricati:\n")
for (nm in names(all_samples))
  cat(sprintf("  %-20s | %5d celle\n", nm, ncol(all_samples[[nm]])))

# ── Funzione: rileva colonna CAR robustamente ───────────────
get_car_status <- function(obj, sample_name) {
  meta    <- obj@meta.data
  car_col <- NULL

  # Priorità: IS_CAR_ALLIN_scREP > IS_CAR > CAR
  for (col in c("IS_CAR_ALLIN_scREP", "IS_CAR", "CAR")) {
    if (col %in% colnames(meta)) { car_col <- col; break }
  }

  if (is.null(car_col)) {
    cat(sprintf("[WARN] %s: nessuna colonna CAR trovata → tutto CAR-\n",
                sample_name))
    return(rep("CAR-", ncol(obj)))
  }

  vals <- as.character(meta[[car_col]])
  car_pos <- grepl("^(YES|TRUE|yes|true|1)$", vals)
  cat(sprintf("  %s: colonna CAR = '%s' | CAR+ = %d (%.1f%%)\n",
              sample_name, car_col,
              sum(car_pos), 100 * mean(car_pos)))

  ifelse(car_pos, "CAR+", "CAR-")
}

# ============================================================
# STEP 1: ESTRAZIONE DATI COMPOSIZIONE
# ============================================================
section("STEP 1 | Calcolo proporzioni compartimenti CD4/CD8")

extract_compartment_props <- function(obj, sample_name, patient_id,
                                      timepoint, tissue) {
  meta        <- obj@meta.data
  car_status  <- get_car_status(obj, sample_name)
  cell_type   <- as.character(meta$cell_type)

  df <- data.frame(
    sample    = sample_name,
    patient   = patient_id,
    timepoint = timepoint,    # "I" oppure "AB"
    tissue    = tissue,       # "bone", "blood"
    cell_type = cell_type,
    car       = car_status,
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      compartment = case_when(
        cell_type %in% CD4_TYPES ~ "CD4",
        cell_type %in% CD8_TYPES ~ "CD8",
        TRUE                     ~ "Other"
      )
    ) %>%
    filter(compartment != "Other")

  # Calcola proporzioni DENTRO ogni compartimento × car_status
  df %>%
    group_by(sample, patient, timepoint, tissue, compartment, car, cell_type) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(sample, patient, timepoint, tissue, compartment, car) %>%
    mutate(
      n_compartment = sum(n),
      pct           = round(100 * n / n_compartment, 2)
    ) %>%
    ungroup()
}

# Costruisce dataset completo
all_props <- bind_rows(lapply(names(PATIENT_MAP), function(pid) {
  pm  <- PATIENT_MAP[[pid]]
  res <- list()

  # Campione I
  if (pm$I %in% names(all_samples)) {
    res[["I"]] <- extract_compartment_props(
      all_samples[[pm$I]], pm$I, pid, "I", "bone")
  } else {
    cat(sprintf("[WARN] %s: campione I '%s' non trovato\n", pid, pm$I))
  }

  # Campioni AB
  for (nm in pm$AB) {
    if (!nm %in% names(all_samples)) next
    tissue_lbl <- if (grepl("bone", nm)) "bone" else "blood"
    res[[nm]]  <- extract_compartment_props(
      all_samples[[nm]], nm, pid, "AB", tissue_lbl)
  }

  bind_rows(res)
}))

cat("\nProprietà estratte:\n")
print(table(all_props$patient, all_props$timepoint, all_props$compartment))

# Controlla celle per categoria (alcune potrebbero essere 0)
cat("\nCellule CAR+ per compartimento e timepoint:\n")
all_props %>%
  filter(car == "CAR+") %>%
  group_by(patient, timepoint, compartment) %>%
  summarise(n_total = sum(n), .groups = "drop") %>%
  print()

# ============================================================
# STEP 2: VISUALIZZAZIONE – BARPLOT PER PAZIENTE
# ============================================================
section("STEP 2 | Barplot composizione per paziente")

# Helper: crea barplot stacked per un paziente e un compartimento
make_composition_barplot <- function(df_patient, compartment_type,
                                     palette, order_vec, patient_id) {

  df_c <- df_patient %>%
    filter(compartment == compartment_type) %>%
    mutate(
      cell_type  = factor(cell_type, levels = rev(intersect(order_vec, unique(cell_type)))),
      # Etichetta asse x: sample_name → timepoint + tissue + CAR
      x_label    = paste0(timepoint, "\n", tissue, "\n", car),
      x_label    = factor(x_label, levels = unique(x_label[order(timepoint, tissue, car)]))
    )

  if (nrow(df_c) == 0) {
    cat(sprintf("  [SKIP] %s – compartimento %s vuoto\n",
                patient_id, compartment_type))
    return(NULL)
  }

  # Colori per i tipi presenti
  types_present <- as.character(unique(df_c$cell_type))
  cols <- palette[types_present]
  # Se un tipo non è in palette, colore grigio
  cols[is.na(cols)] <- "grey60"

  ggplot(df_c, aes(x = x_label, y = pct, fill = cell_type)) +
    geom_col(position = "stack", width = 0.80, color = "white",
             linewidth = 0.2) +
    geom_text(
      data = df_c %>%
        filter(pct >= 5),  # mostra % solo se >= 5%
      aes(label = paste0(round(pct, 0), "%")),
      position = position_stack(vjust = 0.5),
      size = 2.8, color = "white", fontface = "bold"
    ) +
    scale_fill_manual(values = cols, name = "Cell type") +
    scale_y_continuous(
      labels = function(x) paste0(x, "%"),
      limits = c(0, 105), expand = c(0, 0)
    ) +
    # Linea verticale tra I e AB
    geom_vline(
      xintercept = sum(unique(df_c$timepoint) == "I") + 0.5,
      linetype = "dashed", color = "gray40", linewidth = 0.6
    ) +
    annotate("text", x = 0.6, y = 103, label = "↑ Infusion product",
             size = 2.5, color = "gray40", hjust = 0) +
    annotate("text",
             x    = sum(unique(df_c$timepoint[order(df_c$timepoint)]) == "I") + 1,
             y    = 103,
             label = "Post-infusion ↑",
             size = 2.5, color = "gray40", hjust = 0) +
    labs(
      title    = paste0(patient_id, " – Compartimento ",
                        compartment_type, "+"),
      subtitle = paste0("Proporzione sottotipi DENTRO il compartimento ",
                        compartment_type, " (CAR+ vs CAR-)"),
      x = "Campione | Tessuto | Stato CAR",
      y = paste0("% dentro compartimento ", compartment_type, "+")
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title    = element_text(face = "bold", hjust = 0.5, size = 12),
      plot.subtitle = element_text(hjust = 0.5, color = "gray40", size = 9),
      axis.text.x   = element_text(size = 8),
      axis.title.x  = element_text(size = 9),
      legend.text   = element_text(size = 8),
      legend.key.size = unit(0.35, "cm")
    )
}

# Genera un pannello per ogni paziente: CD4 + CD8
for (pid in names(PATIENT_MAP)) {
  cat(paste0("\n── Paziente ", pid, " ──\n"))
  df_p <- all_props %>% filter(patient == pid)

  if (nrow(df_p) == 0) { cat("  Nessun dato.\n"); next }

  p_cd4 <- make_composition_barplot(df_p, "CD4", PALETTE_CD4,
                                     CD4_ORDER, pid)
  p_cd8 <- make_composition_barplot(df_p, "CD8", PALETTE_CD8,
                                     CD8_ORDER, pid)

  panels <- Filter(Negate(is.null), list(p_cd4, p_cd8))
  if (length(panels) == 0) next

  p_combined <- wrap_plots(panels, ncol = 2) +
    plot_annotation(
      title  = paste0("Paziente ", pid, " – Composizione CD4/CD8 nel tempo"),
      theme  = theme(plot.title = element_text(face = "bold",
                                               hjust = 0.5, size = 14))
    )

  out_path <- paste0(out_dir, "Q1_", pid, "_CD4CD8_composition.png")
  ggsave(out_path, p_combined, width = 14, height = 7,
         dpi = 300, bg = "white")
  cat(paste0("  → ", out_path, "\n"))
}

# ============================================================
# STEP 3: HEATMAP CONFRONTO INTER-PAZIENTE
#         (proporzione media per tipo, per timepoint e CAR status)
# ============================================================
section("STEP 3 | Heatmap confronto inter-paziente")

# Calcola proporzione media per tipo × timepoint × CAR × paziente
heatmap_df <- all_props %>%
  mutate(
    group_label = paste0(patient, " | ", timepoint, " | ", car)
  ) %>%
  group_by(group_label, compartment, cell_type) %>%
  # Media delle proporzioni tra campioni dello stesso paziente/timepoint
  # (es. bone_AB e blood_AB vengono mediati; se si vuole solo midollo,
  #  filtrare: filter(tissue == "bone") prima di questo blocco)
  summarise(pct_mean = mean(pct), .groups = "drop")

# Heatmap per compartimento CD4
make_heatmap <- function(df_heat, comp, palette, order_vec, title_str) {
  df_c <- df_heat %>%
    filter(compartment == comp) %>%
    mutate(
      cell_type   = factor(cell_type,
                           levels = intersect(order_vec, unique(cell_type))),
      group_label = factor(group_label,
                           levels = sort(unique(group_label)))
    )

  if (nrow(df_c) == 0) return(NULL)

  ggplot(df_c, aes(x = group_label, y = cell_type, fill = pct_mean)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = paste0(round(pct_mean, 0), "%")),
              size = 3, color = "black") +
    scale_fill_gradientn(
      colors = c("white", "#FFFDE7", "#FFD600", "#F57F17", "#B71C1C"),
      name   = "% media\nnell'compartimento",
      limits = c(0, 100)
    ) +
    scale_y_discrete(limits = rev(levels(df_c$cell_type))) +
    labs(
      title    = title_str,
      subtitle = paste0("Compartimento ", comp, "+\n",
                        "I=prodotto infusione | AB=post-infusione | ",
                        "valori = media tra campioni dello stesso gruppo"),
      x = NULL, y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title    = element_text(face = "bold", hjust = 0.5, size = 12),
      plot.subtitle = element_text(hjust = 0.5, color = "gray40", size = 9),
      axis.text.x   = element_text(angle = 40, hjust = 1, size = 9),
      axis.text.y   = element_text(size = 10),
      panel.grid    = element_blank()
    )
}

p_heat_cd4 <- make_heatmap(heatmap_df, "CD4", PALETTE_CD4, CD4_ORDER,
                            "Tutti i pazienti – Compartimento CD4+")
p_heat_cd8 <- make_heatmap(heatmap_df, "CD8", PALETTE_CD8, CD8_ORDER,
                            "Tutti i pazienti – Compartimento CD8+")

heat_panels <- Filter(Negate(is.null), list(p_heat_cd4, p_heat_cd8))
if (length(heat_panels) > 0) {
  p_heat_all <- wrap_plots(heat_panels, ncol = 1) +
    plot_annotation(
      title  = "Panoramica inter-paziente: proporzione sottotipi CD4/CD8",
      theme  = theme(plot.title = element_text(face = "bold",
                                               hjust = 0.5, size = 14))
    )
  out_heat <- paste0(out_dir, "Q1_ALL_patients_composition_heatmap.png")
  ggsave(out_heat, p_heat_all,
         width = max(14, length(unique(heatmap_df$group_label)) * 1.2 + 4),
         height = 14, dpi = 300, bg = "white")
  cat(paste0("  → ", out_heat, "\n"))
}

# ============================================================
# STEP 4: TABELLA QUANTITATIVA + SALVATAGGIO XLSX
# ============================================================
section("STEP 4 | Salvataggio tabella")

# Tabella per paziente/compartimento: ampia (pivoted)
summary_wide <- all_props %>%
  select(sample, patient, timepoint, tissue, compartment, car,
         cell_type, n, n_compartment, pct) %>%
  arrange(patient, compartment, timepoint, tissue, car,
          match(cell_type, c(CD4_ORDER, CD8_ORDER)))

wb <- createWorkbook()

addWorksheet(wb, "Proporzioni_complete")
writeData(wb, "Proporzioni_complete", summary_wide)

# Tabella pivot: cell_type vs sample (per CAR+)
pivot_car_pos <- all_props %>%
  filter(car == "CAR+") %>%
  select(sample, compartment, cell_type, pct) %>%
  pivot_wider(names_from = sample, values_from = pct, values_fill = 0)
addWorksheet(wb, "Pivot_CARpos")
writeData(wb, "Pivot_CARpos", pivot_car_pos)

# Tabella pivot: cell_type vs sample (per CAR-)
pivot_car_neg <- all_props %>%
  filter(car == "CAR-") %>%
  select(sample, compartment, cell_type, pct) %>%
  pivot_wider(names_from = sample, values_from = pct, values_fill = 0)
addWorksheet(wb, "Pivot_CARneg")
writeData(wb, "Pivot_CARneg", pivot_car_neg)

saveWorkbook(wb, paste0(out_dir, "Q1_composition_summary.xlsx"),
             overwrite = TRUE)
cat(paste0("  → Q1_composition_summary.xlsx\n"))

cat(paste0(
  "\n", strrep("=", 65), "\n",
  "  Q1 COMPLETATA\n\n",
  "  INTERPRETAZIONE ATTESA:\n",
  "  - Confronta CD4+ CAR+ in I: sono prevalentemente Naive?\n",
  "    Effector? Se c'è uno shift I→AB verso Effector/Cytotoxic,\n",
  "    suggerisce attivazione in vivo delle CART cells.\n",
  "  - Compartimento CD8+ CAR+: presenza di Cytotoxic CD8+ in I\n",
  "    vs espansione in AB → quantifica l'impatto dell'infusione.\n",
  "  - Confronto CAR+ vs CAR-: il prodotto I contiene cellule CD4\n",
  "    naive sia CAR+ che CAR-? O solo i CAR+ sono naive?\n",
  "  NOTA: N=3 pazienti → interpretare pattern, non test formali.\n",
  strrep("=", 65), "\n"
))
