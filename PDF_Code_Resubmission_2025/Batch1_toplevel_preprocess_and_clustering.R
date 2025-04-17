# ---- Load Libraries ----
library(Seurat)
library(stringr)
library(dplyr)

# ---- STEP 1: Load & Preprocess Each Sample ----
# Read file list (one filename per line)
allfiles <- read.csv(file = "file_list.csv", sep = "\n", header = FALSE)
samples <- allfiles$V1

obj_list <- list()

for (sample in samples) {
  sample_name <- substr(sample, 1, 5)
  message("Processing sample: ", sample_name)

  location <- paste0("path_to_data/", sample)  # <-- Replace with actual base path
  raw_data <- Read10X_h5(location)

  # Extract expression matrix (handle format differences)
  sample_mRNA <- if (str_detect(sample, "^PZ")) raw_data else raw_data$`Gene Expression`

  # Filter out ribosomal/pseudogenes and cells with high mito % or low UMI count
  sample_filtered <- sample_mRNA[grep("\\.|^RP[0-9]|-PS", rownames(sample_mRNA), invert = TRUE), ]
  mtpct <- colSums(sample_filtered[grep("^MT-", rownames(sample_filtered)), ]) / colSums(sample_filtered) * 100
  sample_filtered <- sample_filtered[grep("^MT-", rownames(sample_filtered), invert = TRUE), which(mtpct < 10)]

  umicount <- colSums(sample_filtered)
  sample_filtered <- sample_filtered[, which(umicount > 100)]
  umicount <- colSums(sample_filtered)  # update again

  # Rename cells with sample ID
  colnames(sample_filtered) <- paste(colnames(sample_filtered), sample_name, sep = "_")

  # Create Seurat object
  sample_obj <- CreateSeuratObject(sample_filtered, project = sample_name)
  sample_obj$sample_ID <- sample_name
  obj_list[[sample_name]] <- sample_obj

  # Optional logging
  write(sample_name, file = "sample_qc_log.txt", append = TRUE)
  write(ncol(sample_obj), file = "sample_qc_log.txt", append = TRUE)
}

# Save list of objects
save(obj_list, file = "obj_list.rda")

# ---- STEP 2: Merge and Preprocess All Samples ----
# Optionally remove problematic sample
obj_list$EZ020 <- NULL  # High variation; failed to cluster well even with UMI threshold

# Designate one sample as base, merge all others
PZ001 <- obj_list$PZ001
obj_list$PZ001 <- NULL

Merged_object <- merge(PZ001, y = obj_list, project = "PD_analysis")
rm(obj_list, PZ001)

# Save merged object before normalization
save(Merged_object, file = "Merged.raw.obj.rda")

# ---- STEP 3: Clustering Pipeline ----
Merged_object <- Merged_object %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(verbose = FALSE)

Merged_object <- FindNeighbors(Merged_object, dims = 1:15)
Merged_object <- FindClusters(Merged_object, resolution = c(0.2, 0.3, 0.4, 0.5))
Merged_object <- RunUMAP(Merged_object, reduction = "pca", dims = 1:15)

# Final save
save(Merged_object, file = "batch1_merged_toplevel_clustering_20220502.RData")
