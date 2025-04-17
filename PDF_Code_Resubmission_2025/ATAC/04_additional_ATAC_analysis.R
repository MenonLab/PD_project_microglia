setwd("path/to/home")

packages <- c("tidyverse","stringr","RColorBrewer","limma","patchwork",
              "ArchR","plyranges","BSgenome.Hsapiens.UCSC.hg38")
sapply(packages, function(package) library(package,character.only = T))

pathToMacs2 <- "/path/to/bin/macs2"

#################
proj <-loadArchRProject(path="EZ_Final_Micro")
peak_set <- getPeakSet(proj)
markerList <- peak_set %>% as_tibble %>%
  mutate(name = paste0(seqnames,"_",start,"_",end))

all_res <- read.csv("pseudobulk_peak_results.csv") #pseudobulk peak results

new_clusters <- readRDS("path/to/micro_defined_taxonomy.rda") # snRNA-seq object
matches <- readRDS("EZ_Final_Micro/Annotations/Motif-Matches-In-Peaks.rds")
key_marks <- read.csv("path/to/micro_final_markers_wilcox_names.csv") #Sub-cluster marker genes

#######

#### Peak type enrichment, up-regulated peaks
p_peaks <- all_res %>% filter(cluster == y,p_adj <= .05,logFC > 0.3)
p_peaks <- markerList %>% filter(name %in% p_peaks$peak)
p_sig <- markerList %>% filter(name %in% p_peaks$name)
p_ns <- markerList %>% filter(!(name %in% p_peaks$name)) %>% filter(name %in% all_peaks$name)
p_table <- table(markerList$peakType[markerList$name %in% all_peaks$name], c(rep("sig",length(p_sig)), rep("ns",length(p_ns))))

res <- map_dfr(1:4, function(i){
  marginal_table <- as.matrix(p_table)
  marginal_table <- rbind(marginal_table[i,], apply(marginal_table[-i,], 2,sum))
  result <- fisher_test(marginal_table, alternative="greater") %>% mutate(peakType = rownames(p_table)[i]) %>%
    relocate(peakType)
  return(result)
})
res %>% write.csv(paste0(y,"_peakType_up-regulated-markers_fisher_test.csv"),quote=F,row.names=F)


#### Peak type enrichment, down-regulated peaks
n_peaks <- all_res %>% filter(cluster == y,p_adj <= .05,logFC < -0.3)
n_peaks <- markerList %>% filter(name %in% n_peaks$peak)
n_sig <- markerList %>% filter(name %in% n_peaks$name)
n_ns <- markerList %>% filter(!(name %in% n_peaks$name)) %>% filter(name %in% all_peaks$name)
n_table <- table(markerList$peakType[markerList$name %in% all_peaks$name], c(rep("sig",length(n_sig)), rep("ns",length(n_ns))))

res <- map_dfr(1:4, function(i){
  marginal_table <- as.matrix(n_table)
  marginal_table <- rbind(marginal_table[i,], apply(marginal_table[-i,], 2,sum))
  result <- fisher_test(marginal_table, alternative="greater") %>% mutate(peakType = rownames(n_table)[i]) %>%
    relocate(peakType)
  return(result)
})
res %>% write.csv(paste0(y,"_peakType_down-regulated-markers_fisher_test.csv"),quote=F,row.names=F)


##########

#Motif enrichment of peaks correlated to marker genes
p2g <- getPeak2GeneLinks(
  ArchRProj = proj,
  corCutOff = 0.25,
  resolution = 1000,
  FDRCutOff = 5e-02,
  returnLoops = F
)

idxRNA <- data.frame(idxRNA = 1:length(rownames(new_clusters)), gene = rownames(new_clusters))
idxATAC <- getPeakSet(proj) %>%
  mutate(peak = paste0(seqnames,"_",start, "_",start+500)) %>%
  as_tibble() %>% as.data.frame
idxATAC$idxATAC <-1:nrow(idxATAC)
idxATAC <- idxATAC %>% dplyr::select(idxATAC,peak)

p2g <- p2g %>%
  as.data.frame() %>%
  left_join(idxRNA, by = "idxRNA") %>%
  left_join(idxATAC, by = "idxATAC")


poi_input1 = markerList$name[p2g %>% filter(gene =="CD83") %>%pull(idxATAC)]
poi_input2 = markerList$name[p2g %>% filter(gene %in% c("CD83","BCAS2")) %>%pull(idxATAC)]
poi_input3 = markerList$name[p2g %>% filter(gene %in% (key_marks %>%filter(p_val<=.05, grepl("CD83",final_names))%>% 
                                                         top_n(9,avg_log2FC) %>% pull(gene))) %>%pull(idxATAC) %>% unique()]
poi_input4 = markerList$name[p2g %>% filter(gene %in% (key_marks %>%filter(p_val<=.05, grepl("CD83",final_names))%>% 
                                                         top_n(9,pct.1-pct.2) %>% pull(gene))) %>%pull(idxATAC) %>% unique()]
poi_input5 = markerList$name[p2g %>% filter(gene %in% (key_marks %>%filter(p_val<=.05, grepl("CD83",final_names))%>%
                                                         pull(gene))) %>%pull(idxATAC) %>% unique()]

pois_list <- list(poi_input1,poi_input2,poi_input3,poi_input4,poi_input5)
names(pois_list) <- c("CD83","CD83+BCAS2","Top 9 FC","Top 9 PCT","ALL")
names(pois_list) <- c("TMEM163","TMEM163+HAMP","Top 9 FC","Top 9 PCT","ALL")

motifs_marks <- map_dfr(1:5, function(x){
  test <- peakAnnoEnrichment2(seMarker = markerList,poi=pois_list[[x]],
                              peakAnnotation = "Motif",peakSet = peak_set,
                              matches = matches, background = "all") %>% 
    mutate(feature = str_remove(feature,"_[0-9]+$")) %>%
    arrange(p_adj) %>% as_tibble() %>% as.data.frame() %>%
    mutate(list = names(pois_list)[x])
  return(test)
})
colnames(motifs_marks)[ncol(motifs_marks)] <- "cluster"

motifs_marks %>% write.csv("CD83_BCAS2_P2G_LINKS_motifs_enrichment.csv", quote=F,row.names=F)
