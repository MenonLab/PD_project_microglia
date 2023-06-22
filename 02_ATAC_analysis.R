setwd("path/to/home")

packages <- c("tidyverse","stringr","RColorBrewer","limma","patchwork",
              "ArchR","BSgenome.Hsapiens.UCSC.hg38")
sapply(packages, function(package) library(package,character.only = T))

#################
### Add metadata into ArchR object
proj <-loadArchRProject(path="PD")
proj <- filterDoublets(proj)

annot <- readRDS(file = "/path/to/PD_RNA_top_level_metadata.RDS")
annot <- annot %>%
  rownames_to_column("id") %>%
  mutate(original = id) %>%
  separate(id, sep = "_", into = c("id", "batch")) %>%
  select(sample_ID, id, batch, original, TopLevel_CellType, seurat_clusters)

cell_ids <- data.frame(id = rownames(proj@cellColData),
                       original = rownames(proj@cellColData)) %>%
  separate(id, sep="#", into = c("sample_ID","id")) %>%
  left_join(annot, by = c("sample_ID","id"))

## Adding metadata to ArchR object
proj$sample_ID <- cell_ids$sample_ID
proj$id <- cell_ids$id
proj$ID  <- cell_ids$original.y
proj$TopLevel_CellType = cell_ids$TopLevel_CellType
proj <- proj[!is.na(proj$TopLevel_CellType),]

ez_id <- read_csv("/path_to/EZ_Meta.csv")
cell_ids <- data.frame(id = rownames(proj@cellColData), SeqID = proj$sample_ID) %>%
  left_join(ez_id, by = "SeqID")
proj$Condition  <- cell_ids$Condition
proj$Source <- cell_ids$Source
proj$Source_L2 <-cell_ids$Source_L2

rm(cell_ids,ez_id,annot)

####################
### Create subset ArchR object for microglia
### Add RNA clusters in Microglia ArchR object

proj_1 <- proj[proj$TopLevel_CellType == "Microglia",]
rm(proj)

microglia_meta <- read.csv("path/to/microgliametadata.csv")

load("/path/to/merged_microglia_clustering_metadata_v4_20221219.RData")
microglia_se <- microglia_merged_obj_v4
rm(microglia_merged_obj_v4)
microglia_se <- RenameCells(microglia_se, old.names = rownames(microglia_se@meta.data),
                            new.names = paste0(microglia_se$orig.ident,"#", str_remove(rownames(microglia_se@meta.data),"_[0-9]+")))

proj_1 <- addCellColData(
  ArchRProj = proj_1,
  data = as.character(microglia_se$RNA_snn_res.0.3)[colnames(microglia_se) %in% rownames(proj_1@cellColData)],
  name = "RNA_res.0.3",
  cells = colnames(microglia_se)[colnames(microglia_se) %in% rownames(proj_1@cellColData)]
)

proj_1$RNA_res.0.3[is.na(proj_1$RNA_res.0.3)] <- -1
proj_1 <- proj_1[proj_1$RNA_res.0.3!=-1,]

saveArchRProject(ArchRProj = proj_1, outputDirectory = paste0("Microglia_PD_NPD_0.3"), load = FALSE)   
rm(proj_1)

ct <- "Microglia"
suffix <- "PD_NPD_0.3"
group_var <- "Condition"  

proj <- loadArchRProject(path=paste0("/path/to/",ct,"_",suffix))
proj@cellColData$RNA_snn_res.0.3 <- factor(proj@cellColData$RNA_snn_res.0.3,
                                           levels=c(0:14),
                                           labels = c(paste0("M",1:7),"PVM",paste0("M",8:13),"Mono"))
proj@cellColData$RNA_snn_res.0.3 <- factor(proj@cellColData$RNA_snn_res.0.3,
                                           levels=c(paste0("M",1:13),"PVM","Mono"),
                                           labels = c(paste0("M",1:13),"PVM","MONO"))


proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI",
                        iterations = 4,
                        clusterParams = list(
                          resolution = c(0.1, 0.2, 0.4), sampleCells = 10000, n.start = 10),
                        varFeatures = 15000,
                        dimsToUse = 1:30)

proj<- addHarmony(ArchRProj = proj, reducedDims = "IterativeLSI", name = "Harmony", groupBy = "sample_ID")
proj <- addClusters(input = proj, reducedDims = "Harmony", method = "Seurat", name = "Clusters", resolution = 0.8)
proj <- addUMAP(ArchRProj = proj, reducedDims = "Harmony", name = "UMAP", nNeighbors = 30, minDist = 0.5, metric = "cosine")
saveArchRProject(ArchRProj = proj, outputDirectory = paste0(ct,"_",suffix), load = FALSE)

###########
####Add Peak Matrix and Deviations Matrix

pathToMacs2 <-  "/path/to/macs2"
    

proj <- loadArchRProject(path=paste0("/path/to/",ct,"_",suffix))
proj$Condition <- as.character(proj$Condition)
proj@genomeAnnotation$genome <- "BSgenome.Hsapiens.UCSC.hg38"
proj2 <- addGroupCoverages(ArchRProj = proj, groupBy = group_var,threads=1)
proj2 <- addReproduciblePeakSet(
  ArchRProj = proj2, 
  groupBy = group_var, 
  pathToMacs2 = pathToMacs2)
proj3 <- addPeakMatrix(proj2,threads=3)
proj3 <- addMotifAnnotations(ArchRProj = proj3, motifSet = "cisbp", name = "Motif")
proj3 <- addBgdPeaks(proj3)
proj3 <- addDeviationsMatrix(
  ArchRProj = proj3, 
  peakAnnotation = "Motif",
  force = TRUE, threads=1
)
saveArchRProject(ArchRProj = proj3, outputDirectory = paste0(ct,"_",suffix), load = FALSE)

proj3 <- loadArchRProject(path=paste0("/path/to/",ct,"_",suffix))

markersPeaks <- getMarkerFeatures(
  ArchRProj = proj3,
  useMatrix = "PeakMatrix", groupBy =group_var,k=75,
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",normBy="ReadsInTSS", threads=1, bgd="NPD", useGroups="PD")

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.05 & abs(Log2FC) >= .2", returnGR = TRUE)

list_df <- data.frame()
for(i in seq_along(markerList)){
  df <- markerList[[i]] %>% as.data.frame()
  if(nrow(df) >0){
    df$cluster <- paste0("cluster",i-1)
    list_df <- rbind(list_df,df)
  }
}
list_df %>% write_csv(paste0(ct,"_",suffix,"_differential_peaks.csv"))


motif_enrichment_in_markers <- function(direction){
  if(!(direction %in% c("M","P"))){
    print("Direction must be M or P")
    stop()
  }
  if(direction =="P"){
    cutoff_input <- "FDR <= 0.05 & Log2FC >= 0.2"
  } else{
    cutoff_input <- "FDR <= 0.05 & Log2FC <= -0.2"
  }
  motifs <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = proj3,
    peakAnnotation = "Motif",background ="all",
    cutOff = cutoff_input)
  
  motif_enrich_df <- data.frame()
  for(i in 1:ncol(motifs)){

    sub_df <- data.frame(assays(motifs[,i])[[1]])
    
    for(j in 2:length(assays(motifs[,1]))){
      col <- assays(motifs[,i])[[j]]
      sub_df <- cbind(sub_df, col)
    }
    colnames(sub_df) <- names(assays(motifs))
    sub_df$Cluster <- paste0("cluster",i-1) #index 0
    rownames(sub_df) <- NULL
    motif_enrich_df <- rbind(motif_enrich_df, sub_df)
  }
  motif_enrich_df %>%
    relocate(feature) %>%
    filter(mlog10Padj >= -log10(0.05)) %>%
    write_csv(paste0(ct,"_",suffix,"_motif_enrichment_",direction,"_0.05.csv"))
}
sapply(c("M","P"), function(direction) motif_enrichment_in_markers(direction))   
      
###########
#### Add Gene to Peak Links Matrix
  
load("/path/to/merged_microglia_clustering_metadata_v4_20221219.RData")
      
microglia_se <- microglia_merged_obj_v4
ids <- paste0(microglia_se$sample_ID,"#", str_remove(colnames(microglia_se),"_[0-9]+"))
rm(microglia_merged_obj_v4)

microglia_se <- RenameCells(microglia_se, new.names = ids)
ez_id <- read_csv("/path/to/EZ_Meta.csv") %>% select(SeqID,Condition) %>% distinct()

microglia_se@meta.data <- microglia_se@meta.data %>% left_join(ez_id, by = c("orig.ident"="SeqID"))
rownames(microglia_se@meta.data) <- colnames(microglia_se)

proj3 <- loadArchRProject(path=paste0("path/to/",ct,"_",suffix))
     
proj3@genomeAnnotation$genome <- "BSgenome.Hsapiens.UCSC.hg38"

group_var <- "Condition"

proj3 <- addGeneIntegrationMatrix(
  ArchRProj = proj3, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = microglia_se,
  addToArrow = T,
  groupRNA = group_var,
  nameCell = "predictedCell_Un",
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un",threads=1,force=T
)
rm(microglia_se)

proj3 <- addImputeWeights(proj3)

proj3 <- addPeak2GeneLinks(
  ArchRProj = proj3, corCutOff = .3,
  useMatrix = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI", threads=1
)
      
p2g <- getPeak2GeneLinks(
  ArchRProj = proj3,
  corCutOff = 0.30,
  resolution = 1000,
  returnLoops = TRUE
)
    
p <- plotPeak2GeneHeatmap(ArchRProj = proj3, groupBy = group_var, k=10, corCutOff = 0.30)

saveArchRProject(ArchRProj = proj3, outputDirectory = paste0(ct,"_",suffix), load = FALSE)

    
    
  