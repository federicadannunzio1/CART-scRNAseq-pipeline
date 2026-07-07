# V Gene Family Analysis
# Converted from 06_vgene_families.Rmd

# ======================================================================
# LIBS
# ======================================================================
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(writexl)
library(stringr)
library(DT)
library(patchwork)

# ======================================================================
# PATHS
# ======================================================================
BASE    <- "/Users/federicadannunzio/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/caruana-project/CART/Code/07_clonotypes"
OUT_DIR <- file.path(BASE, "results")
dir.create(file.path(OUT_DIR, "figures"), showWarnings=FALSE)

# Legge i dati post-decontaminazione già calcolati in 05_conserved_families.Rmd
clones <- read_xlsx(file.path(OUT_DIR, "tables", "final_clone_sequences.xlsx"))
cat("Clonotipi unici (paziente × stage × TRA+TRB):", nrow(clones), "\n")
cat("Pazienti:", paste(unique(clones$patient), collapse=", "), "\n")
cat("Stage:   ", paste(unique(clones$stage),   collapse=", "), "\n")

# Aggiungi famiglia V (es. TRBV7-2 → TRBV7)
clones <- clones %>%
  mutate(
    TRB_family = str_extract(TRB_v_gene, "^TRBV[0-9]+"),
    TRA_family = str_extract(TRA_v_gene, "^TRAV[0-9]+"),
    pair_gene   = paste(TRA_v_gene,  TRB_v_gene,  sep=" + "),
    pair_family = paste(TRA_family,  TRB_family,  sep=" + ")
  )

# ======================================================================
# VBETA USAGE
# ======================================================================
vb_freq <- clones %>%
  group_by(patient, TRB_v_gene) %>%
  summarise(n_cells = sum(n_cells), .groups="drop") %>%
  group_by(patient) %>%
  mutate(freq = n_cells / sum(n_cells)) %>%
  ungroup() %>%
  # top 15 Vbeta per paziente per leggibilità
  group_by(patient) %>%
  slice_max(freq, n=15) %>%
  ungroup()

ggplot(vb_freq, aes(x=reorder(TRB_v_gene, freq), y=freq, fill=patient)) +
  geom_col(show.legend=FALSE) +
  facet_wrap(~patient, scales="free_y", ncol=3) +
  coord_flip() +
  scale_fill_manual(values=c(Bo="#E64B35", Ca="#4DBBD5", Me="#00A087")) +
  scale_y_continuous(labels=scales::percent_format(accuracy=1)) +
  labs(title="Top 15 famiglie Vβ — frequenza relativa per paziente",
       x=NULL, y="% cellule CAR+") +
  theme_classic(base_size=11) +
  theme(strip.text=element_text(face="bold", size=12))

# ======================================================================
# VALPHA USAGE
# ======================================================================
va_freq <- clones %>%
  group_by(patient, TRA_v_gene) %>%
  summarise(n_cells = sum(n_cells), .groups="drop") %>%
  group_by(patient) %>%
  mutate(freq = n_cells / sum(n_cells)) %>%
  ungroup() %>%
  group_by(patient) %>%
  slice_max(freq, n=15) %>%
  ungroup()

ggplot(va_freq, aes(x=reorder(TRA_v_gene, freq), y=freq, fill=patient)) +
  geom_col(show.legend=FALSE) +
  facet_wrap(~patient, scales="free_y", ncol=3) +
  coord_flip() +
  scale_fill_manual(values=c(Bo="#E64B35", Ca="#4DBBD5", Me="#00A087")) +
  scale_y_continuous(labels=scales::percent_format(accuracy=1)) +
  labs(title="Top 15 famiglie Vα — frequenza relativa per paziente",
       x=NULL, y="% cellule CAR+") +
  theme_classic(base_size=11) +
  theme(strip.text=element_text(face="bold", size=12))

# ======================================================================
# SHARED VGENE PAIRS
# ======================================================================
# Frequenza di ogni coppia per paziente
pair_per_patient <- clones %>%
  group_by(patient, TRA_v_gene, TRB_v_gene, pair_gene) %>%
  summarise(n_cells = sum(n_cells), .groups="drop") %>%
  group_by(patient) %>%
  mutate(freq = n_cells / sum(n_cells)) %>%
  ungroup()

# Coppie presenti in ≥2 pazienti
shared_pairs <- pair_per_patient %>%
  group_by(pair_gene, TRA_v_gene, TRB_v_gene) %>%
  summarise(
    n_pazienti   = n_distinct(patient),
    pazienti     = paste(sort(unique(patient)), collapse=" & "),
    n_cells_tot  = sum(n_cells),
    freq_Bo      = sum(freq[patient=="Bo"]),
    freq_Ca      = sum(freq[patient=="Ca"]),
    freq_Me      = sum(freq[patient=="Me"]),
    .groups = "drop"
  ) %>%
  filter(n_pazienti >= 2) %>%
  arrange(desc(n_pazienti), desc(n_cells_tot))

cat("Coppie Vα+Vβ presenti in ≥2 pazienti:", nrow(shared_pairs), "\n\n")
print(shared_pairs %>% select(pair_gene, pazienti, n_cells_tot,
                               freq_Bo, freq_Ca, freq_Me) %>%
        mutate(across(starts_with("freq"), ~round(.x*100,2))))

write_xlsx(shared_pairs,
           file.path(OUT_DIR, "tables", "shared_vgene_pairs.xlsx"))

# ======================================================================
# BUBBLE PAIRS
# ======================================================================
# Prendi le top coppie condivise + top coppie private per ciascun paziente
top_private <- pair_per_patient %>%
  group_by(patient) %>%
  slice_max(freq, n=8) %>%
  ungroup()

focus_pairs <- union(shared_pairs$pair_gene, top_private$pair_gene)

plot_bubble <- pair_per_patient %>%
  filter(pair_gene %in% focus_pairs) %>%
  mutate(
    is_shared = pair_gene %in% shared_pairs$pair_gene,
    pair_label = if_else(is_shared,
                         paste0("★ ", pair_gene),   # stella per condivisi
                         pair_gene)
  )

# Ordina per n_pazienti poi freq totale
pair_order <- plot_bubble %>%
  group_by(pair_label) %>%
  summarise(n_paz=n_distinct(patient), tot=sum(n_cells)) %>%
  arrange(desc(n_paz), desc(tot)) %>%
  pull(pair_label)

plot_bubble$pair_label <- factor(plot_bubble$pair_label, levels=rev(pair_order))

ggplot(plot_bubble, aes(x=patient, y=pair_label,
                         size=freq, color=patient)) +
  geom_point(alpha=0.8) +
  scale_size_continuous(range=c(2, 14),
                        labels=scales::percent_format(accuracy=0.1),
                        name="% cellule CAR+") +
  scale_color_manual(values=c(Bo="#E64B35", Ca="#4DBBD5", Me="#00A087"),
                     guide="none") +
  geom_text(aes(label=ifelse(freq>=0.02, scales::percent(freq,accuracy=0.1), "")),
            color="white", size=2.8, fontface="bold") +
  labs(title="Coppie Vα+Vβ nei clonotipi CAR+ per paziente",
       subtitle="★ = coppia presente in ≥2 pazienti  |  dimensione bolla = frequenza relativa",
       x=NULL, y=NULL) +
  theme_classic(base_size=11) +
  theme(axis.text.y = element_text(size=9,
                                    face=ifelse(grepl("★", rev(pair_order)),
                                                "bold","plain")),
        legend.position="right",
        panel.grid.major.y = element_line(color="grey92"))

ggsave(file.path(OUT_DIR, "figures", "bubble_vgene_pairs.png"),
       width=11, height=8, dpi=150, bg="white")

# ======================================================================
# FAMILY USAGE
# ======================================================================
fam_freq <- clones %>%
  group_by(patient, TRB_family) %>%
  summarise(n_cells=sum(n_cells), .groups="drop") %>%
  group_by(patient) %>%
  mutate(freq=n_cells/sum(n_cells)) %>%
  ungroup() %>%
  filter(!is.na(TRB_family))

# Famiglie presenti in ≥2 pazienti con freq >1%
shared_fam <- fam_freq %>%
  filter(freq >= 0.01) %>%
  group_by(TRB_family) %>%
  filter(n_distinct(patient) >= 2) %>%
  ungroup()

ggplot(fam_freq %>% semi_join(shared_fam, by="TRB_family"),
       aes(x=TRB_family, y=freq, fill=patient)) +
  geom_col(position="dodge", alpha=0.85) +
  scale_fill_manual(values=c(Bo="#E64B35", Ca="#4DBBD5", Me="#00A087")) +
  scale_y_continuous(labels=scales::percent_format(accuracy=1)) +
  labs(title="Famiglie Vβ presenti in ≥2 pazienti (freq >1%)",
       x=NULL, y="% cellule CAR+", fill="Paziente") +
  theme_classic(base_size=11) +
  theme(axis.text.x=element_text(angle=35, hjust=1))

ggsave(file.path(OUT_DIR, "figures", "vbeta_families_shared.png"),
       width=10, height=5, dpi=150, bg="white")

# ======================================================================
# FAMILY PAIRS SHARED
# ======================================================================
# Coppie di famiglie condivise tra ≥2 pazienti
fam_pair_freq <- clones %>%
  filter(!is.na(TRA_family), !is.na(TRB_family)) %>%
  group_by(patient, pair_family) %>%
  summarise(n_cells=sum(n_cells), .groups="drop") %>%
  group_by(patient) %>%
  mutate(freq=n_cells/sum(n_cells)) %>%
  ungroup()

shared_fam_pairs <- fam_pair_freq %>%
  filter(freq >= 0.005) %>%
  group_by(pair_family) %>%
  filter(n_distinct(patient) >= 2) %>%
  summarise(
    pazienti    = paste(sort(unique(patient)), collapse=" & "),
    n_pazienti  = n_distinct(patient),
    n_cells_tot = sum(n_cells),
    freq_Bo     = sum(freq[patient=="Bo"]),
    freq_Ca     = sum(freq[patient=="Ca"]),
    freq_Me     = sum(freq[patient=="Me"]),
    .groups="drop"
  ) %>%
  arrange(desc(n_pazienti), desc(n_cells_tot))

cat("Coppie di FAMIGLIE Vα+Vβ in ≥2 pazienti (freq ≥0.5%):\n")
print(shared_fam_pairs %>%
        mutate(across(starts_with("freq"), ~round(.x*100,1))))

write_xlsx(shared_fam_pairs,
           file.path(OUT_DIR, "tables", "shared_vgene_families.xlsx"))

# ======================================================================
# FAMILY HEATMAP
# ======================================================================
# Heatmap frequenza famiglie Vβ × paziente × stage
fam_stage <- clones %>%
  filter(!is.na(TRB_family)) %>%
  group_by(patient, stage, TRB_family) %>%
  summarise(n_cells=sum(n_cells), .groups="drop") %>%
  group_by(patient, stage) %>%
  mutate(freq=n_cells/sum(n_cells)) %>%
  ungroup() %>%
  mutate(pt_stage=factor(paste(patient, stage, sep="-"),
                         levels=c("Bo-I","Bo-A","Bo-B",
                                  "Ca-I","Ca-A","Ca-B",
                                  "Me-I","Me-B")))

# Mantieni solo famiglie presenti con >2% in almeno un campione
top_fam <- fam_stage %>%
  group_by(TRB_family) %>%
  filter(any(freq>0.02)) %>%
  ungroup() %>%
  pull(TRB_family) %>% unique()

ggplot(fam_stage %>% filter(TRB_family %in% top_fam),
       aes(x=pt_stage, y=reorder(TRB_family, freq), fill=freq)) +
  geom_tile(color="white", linewidth=0.4) +
  geom_text(aes(label=ifelse(n_cells>0, n_cells, "")),
            size=3, color="grey20") +
  scale_fill_gradientn(
    colors=c("white","#FEE8C8","#FC8D59","#D73027"),
    labels=scales::percent_format(accuracy=1),
    name="Freq. relativa"
  ) +
  labs(title="Frequenza famiglie Vβ per paziente e stage",
       subtitle="Numero = cellule assolute",
       x=NULL, y="Famiglia Vβ") +
  theme_classic(base_size=11) +
  theme(axis.text.x=element_text(angle=40, hjust=1, size=10),
        panel.grid=element_blank())

ggsave(file.path(OUT_DIR, "figures", "heatmap_vbeta_family_stage.png"),
       width=10, height=7, dpi=150, bg="white")

# ======================================================================
# SUMMARY DT
# ======================================================================
print(DT::datatable(
  shared_pairs %>%
    mutate(across(starts_with("freq"), ~scales::percent(.x, accuracy=0.1))) %>%
    rename(`Coppia Vα+Vβ`=pair_gene, Pazienti=pazienti,
           `N cellule tot`=n_cells_tot,
           `Freq Bo`=freq_Bo, `Freq Ca`=freq_Ca, `Freq Me`=freq_Me),
  caption="Coppie Vα+Vβ presenti in ≥2 pazienti",
  filter="top", options=list(pageLength=30, scrollX=TRUE), rownames=FALSE
))
