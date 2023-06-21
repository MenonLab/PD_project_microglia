library(Seurat)
library(dplyr)

load("../obj_list.rda")

#removing EZ020 as it was showing high variation when looked independently. Even cutoff>1000nUMI couldnt cluster it. 
obj_list$EZ020 <- NULL

PZ001 <- obj_list$PZ001
obj_list$PZ001 <- NULL
Merged_object <- merge(PZ001, y = c(obj_list),  project = "PD analysis")
rm(obj_list,PZ001)

save(Merged_object, file ="Merged.obj.rda") 
#load("Merged.obj.rda")
Merged_object <- NormalizeData(Merged_object) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = FALSE)
head(Merged_object@meta.data)

Merged_object <- FindNeighbors(object =Merged_object, dims = 1:15)
Merged_object <- FindClusters(object =Merged_object, resolution = (0.2,0.3,0.4,0.5))
Merged_object <- RunUMAP(Merged_object, reduction = "pca", dims = 1:15)

save(Merged_object, file="merged_toplevel_clustering_20220502.RData")
head(Merged_object@meta.data)

head(Merged_object@meta.data)
pdf("toplevel_clustering_20220502.pdf",useDingbats=F)
p1= ElbowPlot(Merged_object,ndims=30);plot(p1)
p1=DimPlot(Merged_object,reduction = "umap",group.by="sample_ID")+NoLegend();plot(p1)
p1=DimPlot(Merged_object,reduction = "umap",group.by="sample_ID");plot(p1)
p1=DimPlot(Merged_object,reduction = "umap",label = TRUE);plot(p1)
p1=DotPlot(Merged_object,features=c("SNAP25","SLC17A7","SLC17A6","SYT1","TH","GAD1","GAD2","AQP4","AIF1","P2RYDC1","CSF1R","MOBP","MBP","OPALIN","OLIG1","PD$
p1=FeaturePlot(Merged_object,features=c("SNAP25","SLC17A6","GAD1","AQP4","OPALIN","PDGFRA"));plot(p1)
dev.off()

