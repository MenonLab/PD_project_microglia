library(Seurat)
library(stringr)

allfiles <- read.csv(file ="file_list.csv", sep = "\n", header = F)
samples <- allfiles$V1

obj_list <- list()
for (sample in samples){
  #sample_name <- sapply(strsplit(sample,"-"), `[`, 1)
  sample_name <- substr(sample, 1, 5)

  print(sample_name)
  location = paste("/mnt/mfs/hgrcgrid/shared/MenonLab/snRNAseq/datasets/PDF_przedborski/",sample, sep="")
  raw_data <- Read10X_h5(location)
#str_detect(sample, "^EZ")
if ((str_detect(sample, "^PZ"))==TRUE) {
sample_mRNA <- raw_data}

else{
  sample_mRNA <- raw_data$`Gene Expression`
  }


  sample_filtered=sample_mRNA[grep("\\.|^RP[0-9]|-PS", rownames(sample_mRNA),invert=T),]  ###optional: remove RP and pseudogenes because of high noise
  mtpct=colSums(sample_filtered[grep("^MT-",rownames(sample_filtered)),])/colSums(sample_filtered)*100
  sample_filtered=sample_filtered[grep("^MT-",rownames(sample_filtered),invert=T),which(mtpct<10)]  ###remove most mitochondrial genes
  umicount=colSums(sample_filtered)
  keepcols=which(umicount>100)
  sample_filtered=sample_filtered[,keepcols]
  umicount=colSums(sample_filtered)
  print(sample_name)
  colnames(sample_filtered)=paste(colnames(sample_filtered),sample_name,sep="_")
  head(colnames(sample_filtered))
  sample_obj <- CreateSeuratObject(sample_filtered, project = sample_name)
  sample_obj
print("object_created")

head(sample_obj@meta.data)


  sample_obj$sample_ID <- sample_name
  
  obj_list[[sample_name]] <- sample_obj

  write(sample_name,file="myfile.txt",append=TRUE)
  write(ncol(sample_obj),file="myfile.txt",append=TRUE)
}


rm(raw_data,sample_filtered,sample_obj,sample_mRNA)
#head(obj_list$EZ001@meta.data)
save(obj_list, file="obj_list.rda")



