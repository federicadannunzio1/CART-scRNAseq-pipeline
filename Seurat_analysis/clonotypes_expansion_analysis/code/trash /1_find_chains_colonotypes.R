library(purrr)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(ggplot2)
library(writexl)

# ==============================================================================
# 1. CONFIGURAZIONE PERCORSI
# ==============================================================================
seurat_path <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/CART/Data/seurat_obj_list/seurat_samples_sctype_azimuth_pbmc_bonemarrow_clonalvdj_CAR.rds"
base_path    <- "/Users/federicadannunzio/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/caruana-project/CART/Data/output_allineamento_original_no_CAR_si_VDJ"
output_dir   <- "/Users/federicadannunzio/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/caruana-project/CART/Code/Seurat_analysis/4_clonotypes_expansion_analysis/res"

if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ==============================================================================
# 2. STEP 1: CENSIMENTO CELLULE CAR+ (DA SEURAT)
# ==============================================================================
message("\n--- STEP 1: CENSIMENTO CELLULE CAR+ ---")

if(!file.exists(seurat_path)) stop("ERRORE: File Seurat non trovato!")
seurat_list <- readRDS(seurat_path)
flat_list   <- unlist(seurat_list, recursive = FALSE)

car_cells_map <- map_df(names(flat_list), function(obj_name) {
  meta  <- flat_list[[obj_name]]@meta.data
  is_car <- rep(FALSE, nrow(meta))
  
  if ("IS_CAR_ALLIN_scREP" %in% colnames(meta)) {
    is_car <- grepl("YES|TRUE", as.character(meta$IS_CAR_ALLIN_scREP), ignore.case = TRUE)
  } else if ("CAR" %in% colnames(meta)) {
    is_car <- grepl("YES|TRUE", as.character(meta$CAR), ignore.case = TRUE)
  }
  
  if (sum(is_car) == 0) return(NULL)
  
  meta %>%
    filter(is_car) %>%
    mutate(
      obj_name_r    = obj_name,
      full_barcode  = rownames(.),
      folder_name   = str_extract(full_barcode, "^.*(?=_[ACGT]+-[0-9]+)"),
      clean_barcode = str_extract(full_barcode, "[ACGT]+-[0-9]+")
    ) %>%
    select(obj_name_r, folder_name, clean_barcode)
})

cell_counts <- car_cells_map %>% group_by(folder_name) %>% summarise(Cells_Found = n())
message("Cellule CAR+ trovate per cartella:")
print(cell_counts)

# ==============================================================================
# 3. STEP 2: RECUPERO VDJ CON CDR3 COMPLETO
# ==============================================================================
message("\n--- STEP 2: RECUPERO FILE VDJ ---")

unique_folders <- unique(na.omit(car_cells_map$folder_name))

read_vdj_robust <- function(folder, base_dir) {
  f1 <- file.path(base_dir, "1", folder, "vdj_t", "filtered_contig_annotations.csv")
  f2 <- file.path(base_dir, "2", folder, "vdj_t", "filtered_contig_annotations.csv")
  target <- if(file.exists(f1)) f1 else if(file.exists(f2)) f2 else NULL
  
  if(is.null(target)) { message(paste("  MANCANTE:", folder)); return(NULL) }
  message(paste("  TROVATO:", folder))
  
  tryCatch({
    l1  <- readLines(target, n = 1)
    sep <- if(grepl(";", l1)) ";" else ","
    df  <- read_delim(target, delim = sep, show_col_types = FALSE)
    if(ncol(df) <= 1) df <- read_delim(target, delim = ",", show_col_types = FALSE)
    
    df %>%
      filter(tolower(productive) == "true") %>%
      group_by(barcode) %>%
      summarise(
        # V-genes
        TRA_V   = paste(unique(na.omit(v_gene[chain == "TRA"])), collapse = "/"),
        TRB_V   = paste(unique(na.omit(v_gene[chain == "TRB"])), collapse = "/"),
        # J-genes
        TRA_J   = paste(unique(na.omit(j_gene[chain == "TRA"])), collapse = "/"),
        TRB_J   = paste(unique(na.omit(j_gene[chain == "TRB"])), collapse = "/"),
        # CDR3 - LA VERA FIRMA DEL CLONE
        TRA_CDR3 = paste(unique(na.omit(cdr3[chain == "TRA"])), collapse = "/"),
        TRB_CDR3 = paste(unique(na.omit(cdr3[chain == "TRB"])), collapse = "/"),
        .groups = "drop"
      ) %>%
      mutate(
        folder_name    = folder,
        Gene_Label     = paste(ifelse(TRA_V == "", "?", TRA_V), "+", ifelse(TRB_V == "", "?", TRB_V)),
        Clone_ID_CDR3  = paste(TRA_CDR3, TRB_CDR3, sep = "_"),
        Clone_ID_Full  = paste0(TRA_V, ":", TRA_J, ":", TRA_CDR3, "_", TRB_V, ":", TRB_J, ":", TRB_CDR3),
        Has_Complete_TRA = (TRA_V != "" & TRA_CDR3 != ""),
        Has_Complete_TRB = (TRB_V != "" & TRB_CDR3 != ""),
        Clone_Quality  = case_when(
          Has_Complete_TRA & Has_Complete_TRB  ~ "Complete",
          Has_Complete_TRB & !Has_Complete_TRA ~ "TRB_only",
          Has_Complete_TRA & !Has_Complete_TRB ~ "TRA_only",
          TRUE                                 ~ "Incomplete"
        )
      )
  }, error = function(e) { message(paste("  ERRORE:", folder, e$message)); NULL })
}

vdj_database <- map_df(unique_folders, ~read_vdj_robust(.x, base_path))

# ==============================================================================
# 4. STEP 3: UNIONE SEURAT + VDJ
# ==============================================================================
message("\n--- STEP 3: UNIONE DATI ---")

full_data <- car_cells_map %>%
  inner_join(vdj_database, by = c("folder_name", "clean_barcode" = "barcode")) %>%
  mutate(
    patient = case_when(
      grepl("S151|S393", folder_name) ~ "Ca",
      grepl("S429|S435", folder_name) ~ "Bo",
      grepl("S345|S431", folder_name) ~ "Me",
      TRUE ~ "Unknown"
    ),
    stage = case_when(
      grepl("_I", folder_name) ~ "I",
      grepl("_A", folder_name) ~ "A",
      grepl("_B", folder_name) ~ "B",
      TRUE ~ "Unknown"
    )
  )

message("Distribuzione Pazienti x Stage:")
print(table(full_data$patient, full_data$stage))
message("Qualita cloni:")
print(table(full_data$Clone_Quality, full_data$patient))

# ==============================================================================
# 5. STEP 4: TOP 10 CLONI PER PAZIENTE
# ==============================================================================
message("\n--- STEP 4: TOP 10 CLONI ---")

# --- METODO CDR3 (CORRETTO) ---
top_clones_CDR3 <- full_data %>%
  filter(Clone_Quality == "Complete") %>%
  group_by(patient, Clone_ID_CDR3, Gene_Label, TRA_CDR3, TRB_CDR3) %>%
  summarise(total_n = n(), .groups = "drop") %>%
  arrange(patient, desc(total_n)) %>%
  group_by(patient) %>%
  slice_head(n = 10) %>%
  mutate(Rank = row_number())

# --- METODO V-GENE (CONFRONTO) ---
top_clones_Vgene <- full_data %>%
  group_by(patient, Gene_Label) %>%
  summarise(total_n = n(), .groups = "drop") %>%
  arrange(patient, desc(total_n)) %>%
  group_by(patient) %>%
  slice_head(n = 10) %>%
  mutate(Rank_ID = paste0("Clone ", row_number()))

# ==============================================================================
# 6. STEP 5: PREPARAZIONE DATI PER I PLOT
# ==============================================================================
message("\n--- STEP 5: PREPARAZIONE DATI PLOT ---")

# --- PLOT DATA V-GENE ---
plot_ready_Vgene <- full_data %>%
  inner_join(
    top_clones_Vgene %>% select(patient, Gene_Label, Rank_ID),
    by = c("patient", "Gene_Label")
  ) %>%
  group_by(patient, Rank_ID, Gene_Label, stage) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  complete(
    nesting(patient, Rank_ID, Gene_Label),
    stage = c("I", "A", "B"),
    fill = list(n_cells = 0)
  ) %>%
  filter(!(patient == "Me" & stage == "A")) %>%
  mutate(
    stage   = factor(stage,   levels = c("I", "A", "B")),
    Rank_ID = factor(Rank_ID, levels = paste0("Clone ", 10:1))
  )

# --- PLOT DATA CDR3 --
plot_ready_CDR3 <- full_data %>%
  filter(Clone_Quality == "Complete") %>%
  select(patient, Clone_ID_CDR3, stage, Gene_Label, TRA_CDR3, TRB_CDR3) %>%
  inner_join(
    top_clones_CDR3 %>% select(patient, Clone_ID_CDR3, Rank),
    by = c("patient", "Clone_ID_CDR3")
  ) %>%
  group_by(patient, Rank, Gene_Label, stage, TRB_CDR3) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  complete(
    nesting(patient, Rank, Gene_Label, TRB_CDR3),
    stage = c("I", "A", "B"),
    fill = list(n_cells = 0)
  ) %>%
  filter(!(patient == "Me" & stage == "A")) %>%
  mutate(
    stage      = factor(stage, levels = c("I", "A", "B")),
    Rank_Label = factor(paste0("Clone ", Rank), levels = paste0("Clone ", 10:1))
  )

message("plot_ready_Vgene: ", nrow(plot_ready_Vgene), " righe")
message("plot_ready_CDR3:  ", nrow(plot_ready_CDR3),  " righe")

# ==============================================================================
# 7. STEP 6: GENERAZIONE PLOT
# ==============================================================================
message("\n--- STEP 6: GENERAZIONE PLOT ---")

# --- PLOT 1: V-GENE ---
p1 <- ggplot(plot_ready_Vgene, aes(x = Rank_ID, y = n_cells, fill = stage)) +
  geom_bar(stat = "identity", position = position_dodge(preserve = "single"),
           width = 0.8, color = "black", linewidth = 0.2) +
  geom_text(
    data = plot_ready_Vgene %>%
      group_by(patient, Rank_ID) %>%
      slice(1) %>%
      group_by(patient) %>%
      mutate(y_pos = -max(n_cells) * 0.08),
    aes(x = Rank_ID, y = y_pos, label = Gene_Label),
    hjust = 1, size = 3.5, color = "grey20", inherit.aes = FALSE
  ) +
  facet_wrap(~patient, scales = "free", ncol = 1) +
  coord_flip(clip = "off", ylim = c(0, NA)) +
  scale_fill_manual(values = c("I" = "#619CFF", "A" = "#F8766D", "B" = "#00BA38")) +
  theme_minimal() +
  theme(
    plot.margin          = margin(10, 20, 10, 200),
    axis.text.y          = element_blank(),
    axis.ticks.y         = element_blank(),
    axis.text.x          = element_text(size = 10),
    legend.position      = "bottom",
    panel.grid.minor     = element_blank(),
    panel.grid.major.y   = element_blank(),
    strip.text           = element_text(face = "bold", size = 12)
  ) +
  labs(
    title    = "Top 10 Cloni CAR T (Metodo V-Gene)",
    subtitle = "cloni diversi con stesso V-gene",
    y = "Numero di Cellule", x = "", fill = "Stage"
  )

# --- PLOT 2: CDR3 ---
p2 <- ggplot(plot_ready_CDR3, aes(x = Rank_Label, y = n_cells, fill = stage)) +
  geom_bar(stat = "identity", position = position_dodge(preserve = "single"),
           width = 0.8, color = "black", linewidth = 0.2) +
  geom_text(
    data = plot_ready_CDR3 %>%
      group_by(patient, Rank_Label) %>%
      slice(1) %>%
      group_by(patient) %>%
      mutate(y_pos = -max(n_cells) * 0.08),
    aes(x = Rank_Label, y = y_pos, label = Gene_Label),
    hjust = 1, size = 3.5, color = "grey20", inherit.aes = FALSE
  ) +
  facet_wrap(~patient, scales = "free", ncol = 1) +
  coord_flip(clip = "off", ylim = c(0, NA)) +
  scale_fill_manual(values = c("I" = "#619CFF", "A" = "#F8766D", "B" = "#00BA38")) +
  theme_minimal() +
  theme(
    plot.margin          = margin(10, 20, 10, 200),
    axis.text.y          = element_blank(),
    axis.ticks.y         = element_blank(),
    axis.text.x          = element_text(size = 10),
    legend.position      = "bottom",
    panel.grid.minor     = element_blank(),
    panel.grid.major.y   = element_blank(),
    strip.text           = element_text(face = "bold", size = 12)
  ) +
  labs(
    title    = "Top 10 Cloni CAR T (Metodo CDR3)",
    subtitle = "Ogni clone e' definito da V+J+CDR3 univoco. Dal piu' espanso in alto al meno espanso in basso",
    y = "Numero di Cellule", x = "", fill = "Stage"
  )

print(p1)
print(p2)

# ==============================================================================
# 8. STEP 7: SALVATAGGIO
# ==============================================================================
message("\n--- STEP 7: SALVATAGGIO ---")

ggsave(file.path(output_dir, "Top10_Cloni_Vgene_vertical.png"),  p1, width = 12, height = 12, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "Top10_Cloni_CDR3_vertical.png"),   p2, width = 12, height = 12, dpi = 300, bg = "white")

write_xlsx(full_data,          file.path(output_dir, "RISULTATI_Cloni_Dati_Completi_con_CDR3.xlsx"))
write_xlsx(top_clones_CDR3,    file.path(output_dir, "RISULTATI_Top10_Cloni_CDR3.xlsx"))
write_xlsx(top_clones_Vgene,   file.path(output_dir, "RISULTATI_Top10_Cloni_Vgene.xlsx"))
write_xlsx(plot_ready_CDR3,    file.path(output_dir, "RISULTATI_Plot_CDR3.xlsx"))

message("Script completato! File salvati in: ", output_dir)
