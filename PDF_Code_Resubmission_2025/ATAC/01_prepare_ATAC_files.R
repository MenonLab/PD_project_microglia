
library(BSgenome.Hsapiens.UCSC.hg38)
library(ArchR)
library(tidyverse)

addArchRGenome("hg38")

# List all IDs 
ids <- c(paste0("EZ",c(paste0("00",1:9),paste0("0",10:84)),
                "PZ001","PZ002",
                paste0("0",85:99),paste0("",100:172) ))
files <- paste0("/path/to/cell_ranger_output/", ids,"-GEX/",ids,"out2/outs/atac_fragments.tsv.gz")

# Create Arrow Files
set.seed(54)
print("Creating Files")
ArrowFiles <- createArrowFiles(
  inputFiles = files,
  sampleNames = ids,
  minTSS = 4, 
  minFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  subThreading=F
)

ArrowFiles <- paste0(ids,".arrow")

## Remove files with missing information
for(i in seq_along(ArrowFiles)){
  if(length(list.files(pattern=ArrowFiles[i]))==0){
    print(i)
    ArrowFiles[i] <- NA
  }
}

ArrowFiles <- ArrowFiles[!is.na(ArrowFiles)]

print("Doublet Scores")
doubScores <- addDoubletScores(input = ArrowFiles, k = 10,knnMethod = "UMAP", LSIMethod = 1)

print("Create directory")
proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "PD",
  copyArrows = TRUE
)

proj <- saveArchRProject(ArchRProj = proj)