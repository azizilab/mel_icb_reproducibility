#!/usr/bin/env Rscript

### Seurat analysis for cellbender output provided as an argument
#### Author: Jana Biermann, PhD

print(paste('Start:',Sys.time()))

library(dplyr)
library(gplots)
library(ggplot2)
library(ggrastr)
library(DropletUtils)
library(SingleR)
library(celldex)
library(SingleCellExperiment)
library(scater)
library(pheatmap)
library(Matrix)
library(stringi)
library(Seurat)

pat <- commandArgs()[6]
doublet_rate <-0.0644
print(pat)

##### Loading, merging, QC, dimension reduction #####
### Load dataset
system(paste0("aws s3 sync s3://melanoma-ribas/cellbender/", pat,"/ data/",pat,"/ ","--exclude '*' --include '",pat,"_filtered.h5' --quiet"))
seu.data <- Read10X_h5(paste0("data/",pat,"/",pat,'_filtered.h5'), 
                       use.names = TRUE, unique.features = TRUE)

### Initialize the Seurat object with the raw (non-normalized data)
seu_raw <- CreateSeuratObject(counts = seu.data, project = pat, 
                              min.cells = 1, min.features = 1)

# Annotate MT genes
seu_raw <- PercentageFeatureSet(seu_raw, pattern = "^MT-", col.name = "percent.mt")
seu_raw <- PercentageFeatureSet(seu_raw, pattern = "^RPS", col.name = "percent.rps")
seu_raw <- PercentageFeatureSet(seu_raw, pattern = "^RPL", col.name = "percent.rpl")
seu_raw$percent.rp=seu_raw$percent.rps + seu_raw$percent.rpl

# Annotate
seu_raw$sample<-pat
seu_raw$patient<-stri_split_fixed(str = pat, pattern = "_", n = 2)[[1]][1]
seu_raw$condition<-stri_split_fixed(str = pat, pattern = "_", n = 2)[[1]][2]
seu_raw$barcode_orig<-rownames(seu_raw@meta.data)
seu_raw$barcode_pat<-paste0(seu_raw$barcode_orig,'_',pat)

clin<-read.csv('data/clin.csv',na.strings = '')
seu_raw@meta.data<-left_join(seu_raw@meta.data,clin,by='sample')
rownames(seu_raw@meta.data)<-seu_raw$barcode_orig

# Identify doublets using scrublet
#doublet_rate<-doublet_rate[doublet_rate$sample==pat,2]
writeMM(seu_raw@assays$RNA@counts, paste0('data/',pat,'/matrix_',pat,'_raw.mtx'))
system(paste('python3 brain_mets/scrublet/scrublet_code.py', pat, doublet_rate))
doublets <- read.table(paste0('data/',pat,'/doublets_',pat,'_raw.txt'),header = T)
seu_raw[['predicted_doublets']]<-doublets$predicted_doublets
seu_raw[['doublet_scores']]<-doublets$doublet_scores
system(paste0('rm data/',pat,'/matrix_',pat,'_raw.mtx'))
system(paste0('rm data/',pat,'/doublets_',pat,'_raw.txt'))


### subset 
minFeature<-200
maxFeature<- 12000
minCount<- 800
maxCount<- 70000
maxMT<-15

seu <- subset(seu_raw, subset = nFeature_RNA > minFeature & nFeature_RNA < maxFeature & 
                nCount_RNA > minCount & nCount_RNA < maxCount & 
                percent.mt < maxMT & predicted_doublets ==F)

### Workflow RNA
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu)
seu <- RunUMAP(seu, dims = 1:30)
seu <- FindNeighbors(seu, dims = 1:30)
seu <- FindClusters(seu)

### Add TCR info
# inspired by https://github.com/satijalab/seurat/issues/1546
system(paste0("aws s3 sync s3://melanoma-ribas/cellranger/v7.0.0/",pat,"_TCR/ data/",pat,"/ --exclude '*' --include 'filtered_contig_annotations.csv' --include 'clonotypes.csv' "))
tcr <- read.csv(paste0('data/', pat, '/filtered_contig_annotations.csv'))
clono <- read.csv(paste0('data/', pat, '/clonotypes.csv'))

# Subsets so only the first line of each barcode is kept,
# as each entry for given barcode will have same clonotype.
tcr <- tcr[!duplicated(tcr$barcode), ]
tcr$barcode_pat<-paste0(tcr$barcode,'_',pat)

# Only keep the barcode and clonotype columns. 
# We'll get additional clonotype info from the clonotype table.
tcr <- tcr[,c("barcode", 'barcode_pat',"raw_clonotype_id",'chain', 
              'v_gene', 'd_gene', 'j_gene', 'c_gene')]
names(tcr)[names(tcr) == "raw_clonotype_id"] <- "clonotype_id"

# Add info onto our original table by clonotype_id.
tcr <- merge(tcr, clono)

# Add to Seurat object
tcr$barcode <- NULL
seu@meta.data <- left_join(seu@meta.data,tcr,by='barcode_pat')
rownames(seu@meta.data)<-seu$barcode_orig
seu$tcr <- ifelse(is.na(seu$cdr3s_aa) == T, 'no_TCR', 'TCR')
table(seu$tcr)

# mait
table(seu$mait_evidence)
seu$mait <- ifelse(seu$mait_evidence %in% c('TRA:gene', 'TRA:gene;TRB:gene', 'TRA:gene+junction', 
                                            'TRA:gene+junction;TRB:gene', 'TRB:gene'), 'MAIT',
                   ifelse(seu$tcr=='TCR','no_MAIT',NA))
table(seu$mait)

# inkt
table(seu$inkt_evidence)
seu$inkt <- ifelse(seu$inkt_evidence %in% c('TRB:gene'), 'iNKT',
                   ifelse(seu$tcr=='TCR','no_iNKT',NA))
table(seu$inkt)

# both chains
seu$both_chains <- ifelse(grepl('TRB:', seu$cdr3s_aa) == T & 
                            grepl('TRA:', seu$cdr3s_aa) == T, 'both', 
                          ifelse(seu$tcr=='TCR','not_both',NA))
table(seu$both_chains)

# clone_size
seu$clone_size <- ifelse(seu$frequency > 1, 'Expanded', 
                         ifelse(seu$frequency == 1, 'Non-expanded', NA))
table(seu$clone_size)

# mait_inkt
seu$mait_inkt <- ifelse(seu$tcr=='TCR',paste0(seu$mait, ';', seu$inkt),NA)
table(seu$mait_inkt)



### cell type identification
seu_sce <- as.SingleCellExperiment(seu)

bped<-BlueprintEncodeData()
pred_bped_main <- SingleR(test = seu_sce, ref = bped, labels = bped$label.main)
pruneScores(pred_bped_main)
seu[['celltype_bped_main']]<-pred_bped_main$pruned.labels
pred_bped_fine <- SingleR(test = seu_sce, ref = bped, labels = bped$label.fine)
pruneScores(pred_bped_fine)
seu[['celltype_bped_fine']]<-pred_bped_fine$pruned.labels

hpca<-HumanPrimaryCellAtlasData()
pred_hpca_main <- SingleR(test = seu_sce, ref = hpca, labels = hpca$label.main)
pruneScores(pred_hpca_main)
seu[['celltype_hpca_main']]<-pred_hpca_main$pruned.labels
pred_hpca_fine <- SingleR(test = seu_sce, ref = hpca, labels = hpca$label.fine)
pruneScores(pred_hpca_fine)
seu[['celltype_hpca_fine']]<-pred_hpca_fine$pruned.labels


### stats
stats<-as.data.frame(matrix(data=NA,nrow = 1, ncol = 11))
colnames(stats)<-c('sample','n_raw_features','n_raw_cells','n_predicted_doublets','n_features','n_cells','median_features','median_counts','cutoff_features','cutoff_counts','cutoff_mt')
rownames(stats)<-pat
stats$sample<-pat
stats$n_raw_features<-dim(seu_raw@assays$RNA@counts)[1]
stats$n_raw_cells<-dim(seu_raw@assays$RNA@counts)[2]
stats$n_predicted_doublets <-length(which(seu_raw@meta.data$predicted_doublets ==T))
stats$n_features<-dim(seu@assays$RNA@counts)[1]
stats$n_cells<-dim(seu@assays$RNA@counts)[2]
stats$median_features<-round(median(seu@meta.data$nFeature_RNA))
stats$median_counts<-round(median(seu@meta.data$nCount_RNA))
stats$cutoff_features<-paste(minFeature,maxFeature)
stats$cutoff_counts<-paste(minCount,maxCount)
stats$cutoff_mt<-paste(maxMT)

### Save objects
ifelse(!dir.exists(file.path(paste0("data/",pat,'/'))), 
       dir.create(file.path(paste0("data/",pat,'/'))), FALSE)
saveRDS(seu, file = paste0("data/",pat,'/data_',pat,'_cb.rds'))

### write pdf reports
pdf(file = paste0("data/",pat,"/plots_", pat,"_cb.pdf"))

# stats
textplot(t(stats),cex=1.2,halign='left')

# plots raw data
ggplot(seu_raw@meta.data, aes(x=seu_raw$nCount_RNA,y = seu_raw$nFeature_RNA, col=seu_raw$percent.mt)) + 
  rasterise(geom_point(size=0.5,alpha=0.5),dpi=300)+ scale_colour_gradient(low="blue", high="green") + 
  labs(color = "Percent MT") + theme_classic() + ggtitle('Raw object')

ggplot(seu_raw@meta.data, aes(x=seu_raw$nCount_RNA,y = seu_raw$nFeature_RNA, col=seu_raw$doublet_scores)) + 
  rasterise(geom_point(size=0.5,alpha=0.5),dpi=300)+ scale_colour_gradient(low="lightgrey", high="darkviolet") + 
  labs(color = "doublet_scores") + theme_classic()+ ggtitle('Raw object')

print(VlnPlot(seu_raw, features = c("nFeature_RNA", "nCount_RNA", "percent.mt",'percent.rps','percent.rpl','doublet_scores'),
              ncol = 3,group.by = 'sample',pt.size = 0))

# QC plot filtered
ggplot(seu@meta.data, aes(x=seu$nCount_RNA,y = seu$nFeature_RNA, col=seu$percent.mt)) + 
  rasterise(geom_point(size=0.5,alpha=0.5),dpi=300)+ scale_colour_gradient(low="blue", high="green") + 
  labs(color = "Percent MT") + theme_classic()+ ggtitle('Filtered object')

ggplot(seu@meta.data, aes(x=seu$nCount_RNA,y = seu$nFeature_RNA, col=seu$doublet_scores)) + 
  rasterise(geom_point(size=0.5,alpha=0.5),dpi=300)+ scale_colour_gradient(low="lightgrey", high="darkviolet") + 
  labs(color = "doublet_scores") + theme_classic()+ ggtitle('Filtered object')

print(VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt",'percent.rps','percent.rpl','doublet_scores'),
              ncol = 3,group.by = 'patient',pt.size = 0))

FeaturePlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt",'percent.rps','percent.rpl','doublet_scores'), 
            min.cutoff = "q05", max.cutoff = 'q95',order=T, raster = T)

# PCA
print(ElbowPlot(seu))
DimPlot(seu, reduction = "pca",group.by = 'ident',raster = T,shuffle = T)

# UMAP
DimPlot(seu, reduction = "umap",label = T,group.by = 'ident',raster = T,shuffle = T)

# features
FeaturePlot(seu, features = c('MITF',"MLANA", 'PMEL', "DCT"), min.cutoff = "q05",
            max.cutoff = 'q95',order=T, raster = T)

FeaturePlot(seu, features = c('PTPRC',"CD8A", 'CD68', "MKI67"), min.cutoff = "q05",
            max.cutoff = 'q95',order=T, raster = T)

DimPlot(seu,group.by = c('clone_size','mait_inkt','both_chains','tcr'),
        shuffle = T,raster = T)
FeaturePlot(seu,'frequency',order = T,raster = T)


## bped
plotScoreHeatmap(pred_bped_fine, clusters=seu_sce@colData$ident,fontsize = 6,main='pred_bped_fine')
DimPlot(seu, reduction = "umap",label = T,group.by = 'celltype_bped_fine',repel = T,label.size = 2.5,raster = T,shuffle = T) + 
  guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
  theme(legend.text=element_text(size=6))

FeatureScatter(seu,feature1 ='nCount_RNA',feature2 = 'nFeature_RNA',shuffle = T,
               group.by = 'celltype_bped_main',raster = T)

VlnPlot(seu, features = c("nFeature_RNA"),group.by = 'celltype_bped_main',pt.size = 0)

## hpca
plotScoreHeatmap(pred_hpca_main, clusters=seu_sce@colData$ident,fontsize = 6,main='pred_hpca_main')
DimPlot(seu, reduction = "umap",label = T,group.by = 'celltype_hpca_main',repel = T,label.size = 2.5,raster = T,shuffle = T) + 
  guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
  theme(legend.text=element_text(size=6))

dev.off()

system(paste0("aws s3 sync data/",pat,"/ s3://melanoma-ribas/Seurat/",pat,"/ --exclude '*' --include '*_cb*' --exclude '.*' --quiet"))

print(paste('End:',Sys.time()))
