library(Seurat)
library(readxl)
library(ggplot2)
library(patchwork)


seurat_obj <- readRDS("/Users/federicadannunzio/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/caruana-project/Data/seurat_obj_list/seurat_samples_sctype_azimuth_pbmc_bonemarrow_clonalvdj_CAR.rds")

marker_file <- "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/Data/featureplots_cell_types/ignazio_markers.xlsx"   
markers_df <- read_excel(marker_file, sheet = 1, col_names = "markers")

markers_vector <- markers_df$markers

plot_markers_umap <- function(seurat_obj, 
                              markers, 
                              sample_name, 
                              group_name,
                              outdir = "~/federica.dannunzio@uniroma1.it - Google Drive/Drive condivisi/caruana-project/Code/Seurat_analysis/1_make_markers_featureplots/featureplots_ignazio_markers/NK-T_markers") {
  
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  markers_present <- markers[markers %in% rownames(seurat_obj)]
  
  if (length(markers_present) == 0) {
    message(paste("No markers found in:", sample_name, group_name))
    return(NULL)
  }
  
  p <- FeaturePlot(
    seurat_obj,
    features = markers_present,
    reduction = "umap",
    cols = c("lightgrey", "red"),
    ncol = 4
  ) + plot_annotation(
    title = paste("Sample:", sample_name, "| Group:", group_name)
  )
  

  file_path <- file.path(outdir, paste0(sample_name, "_", group_name, "_UMAP_markers.png"))
  ggsave(
    filename = file_path,
    plot = p,
    width = 16, 
    height = 4 * ceiling(length(markers_present) / 4),
    limitsize = FALSE
  )
  
  message(paste("Saved plot to:", file_path))
}

# struttura --> seurat_obj[[sample_name]][[group_name]]
for (sample_name in names(seurat_obj)) {

  current_sample_list <- seurat_obj[[sample_name]]
  
  for (group_name in names(current_sample_list)) {
    
    cat("Processing:", sample_name, "-", group_name, "\n")
    
    plot_markers_umap(
      seurat_obj = current_sample_list[[group_name]], 
      markers = markers_vector,
      sample_name = sample_name,
      group_name = group_name
    )
  }
}
