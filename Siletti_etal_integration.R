library(Seurat)
library(dplyr)
library(harmony)
library(ggplot2)
##load refined inhouse neuron object
load("../../2_subclust_exEX20/micro/2_from_zena/merged_microglia_clustering_metadata_v4_20221219.RData")

## subset only SNCTRTL
Idents(microglia_merged_obj_v4) <- "region"
zena_micro_subsetted <- subset(microglia_merged_obj_v4, idents= "SNCTRL")
#head(zena_micro_subsetted@meta.data)

zena_micro_subsetted$zena_cluster_name <- zena_micro_subsetted$RNA_snn_res.0.3
zena_micro_subsetted$zena_cluster_name <- gsub("^", "Zena_", zena_micro_subsetted$zena_cluster_name)
zena_micro_subsetted$sex <- gsub("Woman", "female", zena_micro_subsetted$sex)
zena_micro_subsetted$sex <- gsub("Man", "male", zena_micro_subsetted$sex)


meta_Zena <- zena_micro_subsetted@meta.data[,c("sample_ID", "region","sex","patient_ID", "zena_cluster_name","age")]
write.csv(meta_Zena, "meta_Zena_microglia.csv")
Zena_counts <- zena_micro_subsetted@assays$RNA@counts
Zena_obj <- CreateSeuratObject(Zena_counts, meta.data=meta_Zena)
Zena_obj$donor_id <- Zena_obj$patient_ID
Zena_obj$patient_ID <- NULL

Zena_obj$Status <- Zena_obj$region

Zena_obj$Status <- gsub("SN", "Zena_PD", Zena_obj$Status)
Zena_obj$Status <- gsub("SNCTRL", "Zena_CTRL", Zena_obj$Status)

Zena_obj$dataset <- "Columbia"
Zena.list <- SplitObject(Zena_obj, split.by = "donor_id")  

head(Zena_obj@meta.data)
#head(Zena.list)
#Zena_obj$dataset <- "Columbia"

rm(microglia_merged_obj_v4,zena_micro_subsetted)

print("Zena_obj prepared")


##load Siletti etal data

Siletti_obj1 <- readRDS("../input/siletti_gene_symbols_microglia.rds")
Siletti_obj1$age <- Siletti_obj1$development_stage
Siletti_obj1$age <- as.integer(gsub("-year-old human stage", "", Siletti_obj1$age))
Siltti_meta <- Siletti_obj1@meta.data[,c("tissue","dissection","age","sex", "donor_id")]

Siletti_obj <- CreateSeuratObject(Siletti_obj1@assays$RNA@data, meta.data=Siltti_meta)
rm(Siletti_obj1,Siltti_meta)

Siletti_obj$region <- Siletti_obj$tissue
Siletti_obj$tissue <- NULL
Siletti_obj$dataset <- "Siletti_etal"

Siletti.list <- SplitObject(Siletti_obj, split.by = "donor_id")

Siletti.list$H18.30.001  <- NULL

head(Siletti_obj@assays$RNA@counts@Dimnames[[1]])




obj.list <- c(Zena_obj, Siletti_obj)
#obj.list <- c(Zena.list, Siletti.list)

obj.list
obj.list <- lapply(X = obj.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})  
  
rm(Zena.list,Siletti.list)

print("Normalisation done")

obj.features <- SelectIntegrationFeatures(object.list = obj.list)
print("obj.features done")

obj.anchors <- FindIntegrationAnchors(object.list = obj.list,anchor.features = obj.features)#, k.filter=150)
print("obj.anchors done")
rm(obj.list)
integrated_microglia <- IntegrateData(anchorset = obj.anchors, k.weight=25)
print("integration complete")
#save(integrated_microglia, file="cross_species_microglia_integrated_log_norm.RData")
integrated_microglia <- ScaleData(integrated_microglia,vars.to.regress=c("nCount_RNA"))

DefaultAssay(integrated_microglia) <- "integrated"
integrated_microglia <- RunPCA(object =integrated_microglia, npcs = 30, verbose = FALSE)
integrated_microglia <- FindNeighbors(object =integrated_microglia, dims = 1:22)
integrated_microglia <- FindClusters(object =integrated_microglia, resolution = c(0.2,0.3,0.4,0.5))
integrated_microglia <- RunUMAP(object =integrated_microglia, reduction = "pca", dims = 1:22)


saveRDS(integrated_microglia, file="clustered_bydataset_Siletti_zena_lognorm_microglia_integrated.rds")

pdf("bydataset_22pcs_lognorm_microglia_AY_Siletti.pdf",useDingbats=F)
p1=ElbowPlot(integrated_microglia, ndims= 30);plot(p1)
p1=DimPlot(integrated_microglia,reduction = "umap",group.by="donor_id")+NoLegend();plot(p1)
p1=DimPlot(integrated_microglia,reduction = "umap",group.by="donor_id");plot(p1)
p1=DimPlot(integrated_microglia,reduction = "umap",group.by="dataset");plot(p1)
p1=DimPlot(integrated_microglia,reduction = "umap",group.by="region");plot(p1)
#p1=DimPlot(integrated_microglia,reduction = "umap",group.by="Status");plot(p1)
p1=DimPlot(integrated_microglia,reduction = "umap",label = TRUE);plot(p1)
p1=DimPlot(integrated_microglia,reduction = "umap",group.by="integrated_snn_res.0.4");plot(p1)

dev.off()









