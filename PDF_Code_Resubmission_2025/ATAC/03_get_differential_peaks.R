
library(edgeR)
library(SingleCellExperiment)
library(Matrix)
library(Matrix.utils)
library(limma)

peak_set <- getPeakSet(proj)
markerList <- peak_set %>% as_tibble %>% mutate(name = paste0(seqnames,"_",start,"_",end))

counts <- getMatrixFromProject(proj,"PeakMatrix")
counts <- assays(counts)[[1]]
counts[is.na(counts)] <- 0
counts <- counts[,apply(counts,2, function(row) sum(row!=0)!=0)]

metadata <- proj@cellColData %>% as.data.frame
metadata$cluster_id <- factor(paste0(metadata$Sample,"_",metadata$Clusters))
colnames(counts) <- NULL

sce <- SingleCellExperiment(assays = list(counts = counts[,]),colData = metadata[,],
                            rowData = data.frame(markerList$name,row.names=markerList$name))
rownames(assays(sce)[[1]]) <- markerList$name
sce <- sce[rowSums(counts(sce)> 1) >= 10, ]
sample_variable <- "cluster_id"

kids <- purrr::set_names(levels(sce$cluster_id))
nk <- length(kids)
sids <- purrr::set_names(levels(as.factor(colData(sce)[[sample_variable]])))

sce$sampleInfo <- factor(colData(sce)[[sample_variable]],
                         labels = unique(colData(sce)[[sample_variable]]),
                         levels= unique(colData(sce)[[sample_variable]]))
ns <- length(sids)

n_cells <- table(sce$sampleInfo) %>%  as.vector()
names(n_cells) <- names(table(sce$sampleInfo))

m <- match(names(n_cells), sce$sampleInfo)
ei <- data.frame(colData(sce)[m, ], 
                 n_cells, row.names = NULL) %>% 
  dplyr::select( "cluster_id", "n_cells")  
groups <- colData(sce)[, c("cluster_id"),drop=F]

pb <- aggregate.Matrix(t(counts(sce)), 
                       groupings = groups, fun = "sum") 

ei <- ei[order(match(ei$cluster_id,rownames(pb))),]
ei$Intercept = 1
ei$sample <- str_extract(as.character(ei$cluster_id),"^EZ[0-9]+")
ei$final_names <- str_remove(as.character(ei$cluster_id),"^EZ[0-9]+_")
tab <- table(proj@cellColData$Sample, proj@cellColData$Batch) %>% as.data.frame() %>% filter(Freq!=0)

ei <- ei %>% left_join(tab %>% dplyr::select(-Freq), by = c("sample"="Var1")) %>%
  mutate(Var2 = as.numeric(factor(Var2, levels = c("batch1","batch2"), labels = 0:1))-1)
# EZ085-EZ172 is batch 2

deg_analysis <- function(final_name){
  print(final_name)
  ei$final_name <- 0
  ei$final_name[ei$final_names %in% final_name] <- 1 
  design <- ei %>%
    dplyr::select(any_of(c("Intercept","final_name","Var2"))) 
  rownames(design) <- ei$cluster_id
  design <- as.matrix(design)
  storage.mode(design) <- "numeric"
  pb1 <- pb
  y <- DGEList(t(pb1))
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  fit <- glmQLFTest(fit, coef = "final_name")
  fit <- topTags(fit, n = Inf, sort.by = "none")$table %>% 
    dplyr::mutate(gene = rownames(.), cluster_id = final_name) %>% 
    dplyr::rename(p_val = PValue, p_adj = FDR) %>%
    mutate(Cell_Type = final_name)
  res <- fit 
  out <- res[,c("Cell_Type","gene", "logFC", "logCPM", "F", "p_val", "p_adj")]
  rownames(out) <- NULL
  return(out)
  
}

all_res <- map_dfr(c("CD83_BCAS2","TMEM163_HAMP"), deg_analysis)
all_res %>% write.csv("pseudobulk_peak_results.csv",row.names=T)

