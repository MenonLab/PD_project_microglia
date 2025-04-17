setwd("path/to/home")

packages <- c("tidyverse","ArchR","BSgenome.Hsapiens.UCSC.hg38")
sapply(packages, function(package) library(package,character.only = T))

pathToMacs2 <- "/path/to/bin/macs2"

#################
### Add metadata into ArchR object
proj <-loadArchRProject(path="PD")
proj <- filterDoublets(proj)

## Adding metadata to ArchR object

new_clusters <- readRDS("path/to/micro_defined_taxonomy.rda")
cluster_ids <-  new_clusters@meta.data %>% select(sample_ID,final_clusters,final_names,Region,
                                                  Condition,Source,Source_L2) %>%
  mutate(cellid = str_remove(colnames(new_clusters),"_[A-Z]+[0-9]+$")) %>%
  mutate(cellid = paste0(sample_ID,"#", str_remove(cellid,"_[0-9]+"))) 

atac_ids <- data.frame(cellid = proj@cellColData %>% rownames()) %>%
  left_join(cluster_ids %>% select(-sample_ID), by = "cellid")
proj@cellColData$Region <- atac_ids$Region
proj <- addCellColData(proj, data = atac_ids$final_names,name =  "final_names",cells = atac_ids$cellid)

proj2 <- proj[!is.na(proj@cellColData$final_names),]
proj <- saveArchRProject(proj2, "EZ_Final_Micro",load=T,dropCells=T)


proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI", 
                        iterations = 4, 
                        clusterParams = list(
                          resolution = c(0.1, 0.2, 0.4), sampleCells = 10000, n.start = 10), 
                        varFeatures = 15000, dimsToUse = 1:30,force=T)

proj<- addHarmony(ArchRProj = proj, reducedDims = "IterativeLSI", name = "Harmony", groupBy = "Batch",force=T)
#proj <- addClusters(input = proj, reducedDims = "Harmony", method = "Seurat", name = "Clusters", resolution = 0.4,force=T)
proj <- addUMAP(ArchRProj = proj, reducedDims = "Harmony", name = "UMAP",
                nNeighbors = 30, minDist = 0.5, metric = "cosine",force=T)

proj <- addGroupCoverages(ArchRProj = proj, groupBy = "final_names",threads=1)
proj <- addReproduciblePeakSet(
  ArchRProj = proj, 
  groupBy ="final_names", 
  pathToMacs2 = pathToMacs2, force=T)

proj <- addPeakMatrix(proj,threads=1,force=T)
proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif",force=T)
proj <- addBgdPeaks(proj,force=T)
proj <- addDeviationsMatrix(proj, peakAnnotation = "Motif",force = TRUE, threads=1)

########

proj <- addGeneIntegrationMatrix(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = new_clusters,
  addToArrow = T,
  groupRNA ="final_names",
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore",
  threads=1,force=T
)

rm(microglia_se)

proj <- addImputeWeights(proj)

proj <- addCoAccessibility(
  ArchRProj = proj,
  reducedDims = "IterativeLSI"
)
proj <- addPeak2GeneLinks(
  ArchRProj = proj, corCutOff = .25,
  useMatrix = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI"
)

p2g <- getPeak2GeneLinks(
  ArchRProj = proj,
  corCutOff = 0.25,
  resolution = 1000,
  FDRCutOff = 5e-02,
  returnLoops = F
)

proj <- saveArchRProject(proj, "EZ_Final_Micro")
