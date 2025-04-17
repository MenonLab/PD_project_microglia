# ---- Load Required Libraries ----
require(Matrix)
library(Seurat)
library(stringr)
library(scCustomize)  # For Read_CellBender_h5_Mat

# ---- STEP 1: Load Sample List & Exclude PSP Samples ----
# Load file names from CSV
allfiles <- read.csv(file = "converted_filenames.csv", sep = "\n", header = FALSE)
samples <- allfiles$V1

# Manually exclude PSP samples based on Serge's review
psp_samples <- c("EZ085", "EZ114")
samples <- samples[!substr(samples, 1, 5) %in% psp_samples]

# ---- STEP 2: Load & Preprocess Each Sample from CellBender Output ----
obj_list <- list()

for (sample in samples) {
  sample_name <- substr(sample, 1, 5)
  message("Processing sample: ", sample_name)

  location <- paste0("/mnt/vast/hpc/MenonLab/snRNAseq/datasets/PDF_przedborski/", sample)

  # Load matrix using CellBender-specific function
  cell_bender_mat <- Read_CellBender_h5_Mat(file_name = location)

  # Filter out Ensembl-style gene names with ":" and noisy genes
  sample_mRNA <- cell_bender_mat[grep(":", rownames(cell_bender_mat), invert = TRUE), ]
  sample_filtered <- sample_mRNA[grep("\\.|^RP[0-9]|-PS", rownames(sample_mRNA), invert = TRUE), ]

  # Filter cells with high mitochondrial content (>10%)
  mtpct <- colSums(sample_filtered[grep("^MT-", rownames(sample_filtered)), ]) / colSums(sample_filtered) * 100
  sample_filtered <- sample_filtered[grep("^MT-", rownames(sample_filtered), invert = TRUE), which(mtpct < 10)]

  # Filter cells with low total UMI count (<100)
  umicount <- colSums(sample_filtered)
  sample_filtered <- sample_filtered[, which(umicount > 100)]
  umicount <- colSums(sample_filtered)  # Recompute after filtering

  # Rename cells with sample ID prefix
  colnames(sample_filtered) <- paste(colnames(sample_filtered), sample_name, sep = "_")

  # Create Seurat object
  seurat_obj <- CreateSeuratObject(sample_filtered, project = sample_name)
  seurat_obj$sample_ID <- sample_name

  obj_list[[sample_name]] <- seurat_obj

  # Logging
  write(sample_name, file = "myfile_100.txt", append = TRUE)
  write(ncol(seurat_obj), file = "myfile_100.txt", append = TRUE)
}

# Save all individual sample objects
save(obj_list, file = "obj_list_100.rda")

# ---- STEP 3: Merge All Objects ----
EZ086 <- obj_list$EZ086
obj_list$EZ086 <- NULL

Merged_object <- merge(EZ086, y = obj_list, project = "PD_analysis_CellBender")
saveRDS(Merged_object, file = "Merged_new_obj.rds")

# ---- STEP 4: Top-Level Clustering ----
Merged_object <- Merged_object %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(verbose = FALSE)

Merged_object <- FindNeighbors(Merged_object, dims = 1:15)
Merged_object <- FindClusters(Merged_object, resolution = c(0.2, 0.3, 0.4, 0.5))
Merged_object <- RunUMAP(Merged_object, reduction = "pca", dims = 1:15)

# Final Save
saveRDS(Merged_object, file = "batch2_merged_toplevel_clustering.rds")
