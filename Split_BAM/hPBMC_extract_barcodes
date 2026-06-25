```r 
######################################################################
#Extarct cell barcodes 
######################################################################
library(Signac)
library(Seurat)
H3K27ac <- readRDS("/home/wahid/project_scHMTF/GSE195725_processed_data/H3K27ac.rds")
H3K27me3 <- readRDS("/home/wahid/project_scHMTF/GSE195725_processed_data/H3K27me3.rds")
H3K4me1 <- readRDS("/home/wahid/project_scHMTF/GSE195725_processed_data/H3K4me1.rds")
H3K4me2 <- readRDS("/home/wahid/project_scHMTF/GSE195725_processed_data/H3K4me2.rds")
H3K4me3 <- readRDS("/home/wahid/project_scHMTF/GSE195725_processed_data/H3K4me3.rds")
H3K9me3 <- readRDS("/home/wahid/project_scHMTF/GSE195725_processed_data/H3K9me3.rds")
# =====================================================================
# Export cell type wise barcode files with cleaned cell type names
# =====================================================================
# List all loaded Seurat objects
seurat_objects <- list(
  H3K27ac = H3K27ac,
  H3K27me3 = H3K27me3,
  H3K4me1 = H3K4me1,
  H3K4me2 = H3K4me2,
  H3K4me3 = H3K4me3,
  H3K9me3 = H3K9me3
)

# Output base directory
output_base <- "/home/wahid/project_scHMTF/GSE195725_processed_data/bookchapter_figure/barcodes_by_celltype"
if (!dir.exists(output_base)) dir.create(output_base, recursive = TRUE)
# Loop over each histone mark
for (mark_name in names(seurat_objects)) {
  obj <- seurat_objects[[mark_name]]
  
  if (!"predicted.celltype.l1" %in% colnames(obj@meta.data)) {
    warning(paste(mark_name, "missing 'predicted.celltype.l1' – skipping"))
    next
  }
  # Create subfolder for this histone mark
  mark_folder <- file.path(output_base, mark_name)
  if (!dir.exists(mark_folder)) dir.create(mark_folder)
  # Get unique cell types (excluding NA)
  cell_types <- unique(obj$predicted.celltype.l1)
  cell_types <- cell_types[!is.na(cell_types)]
  for (ct in cell_types) {
    # Remove all spaces from cell type name (e.g., "CD8 T" -> "CD8T")
    ct_clean <- gsub(" ", "", ct)  
    barcodes <- rownames(obj@meta.data[obj$predicted.celltype.l1 == ct, , drop = FALSE])
    out_file <- file.path(mark_folder, paste0(ct_clean, "_barcodes.txt"))
    writeLines(barcodes, con = out_file)
    cat("Saved", length(barcodes), "barcodes for", mark_name, "->", ct_clean, "\n")
  }
}
cat("All barcode files exported to:", output_base, "\n")
```
