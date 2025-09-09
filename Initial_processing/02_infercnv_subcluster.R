#!/usr/bin/Rscript

#### inferCNV pipeline with sample name as argument
#### Author: Jana Biermann, PhD

print(paste("Start:",Sys.time()))

library(dplyr)
library(Seurat)
library(infercnv)
library(stringr)
library(gplots)
library(ggplot2)
library(scales)
library(viridis)

pat <- commandArgs()[6]

# Read Seurat object
system(paste0("aws s3 sync s3://melanoma-ribas/Seurat/", pat,"/ data/",pat,"/ ","--exclude '*' --include 'data_",pat,"_cb.rds'"))
seu <- readRDS(file = paste0('/home/ubuntu/data/',pat,'/data_',pat,'_cb.rds'))

# counts matrix
counts_matrix <- as.data.frame(GetAssayData(object = seu, slot = 'counts'))

# annotation
annot<-as.data.frame(seu$celltype_bped_fine)
annot[is.na(annot)] <- 'unknown'
immune<-c('CD4+ T-cells','CD4+ Tcm','CD4+ Tem','CD8+ T-cells','CD8+ Tcm','CD8+ Tem','Class-switched memory B-cells','DC','Eosinophils','Macrophages','Macrophages M1','Macrophages M2','Memory B-cells','Monocytes','naive B-cells','Neutrophils','NK cells','Plasma cells','Tregs')
annot$cell<-ifelse(annot$`seu$celltype_bped_fine` %in% immune, 'immune','non-immune')
annot$`seu$celltype_bped_fine`<-NULL

# gene order
gene_order<-read.table("/home/ubuntu/brain_mets/MBM/inferCNV/refdata-gex-GRCh38-2020-A_gen_pos.txt",header = F,row.names = 1)

# create the infercnv object
options(scipen = 100)
options(expressions=10000)
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts_matrix,
                                    annotations_file=annot,
                                    gene_order_file=gene_order,
                                    ref_group_names = 'immune')

# perform infercnv operations to reveal cnv signal
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir=paste0("vol/inferCNV_subcluster_",pat),  # dir is auto-created for storing outputs
                             cluster_by_groups=F,   # cluster
                             denoise=T,
                             HMM=T,
                             analysis_mode='subclusters',
                             tumor_subcluster_partition_method='leiden',
                             leiden_resolution=0.01,
                             #tumor_subcluster_partition_method='random_trees',
                             output_format='pdf',
                             num_threads = 4)

# identify malignant cells
immune<-c('CD4+ T-cells','CD4+ Tcm','CD4+ Tem','CD8+ T-cells','CD8+ Tcm','CD8+ Tem','Class-switched memory B-cells','DC','Eosinophils','Macrophages','Macrophages M1','Macrophages M2','Memory B-cells','Monocytes','naive B-cells','Neutrophils','NK cells','Plasma cells','Tregs')
seu[['immune']]<-ifelse(seu$celltype_bped_fine %in% immune, 'immune','non-immune')
seu<-add_to_seurat(seurat_obj=seu,infercnv_output_path=paste0("vol/inferCNV_subcluster_",pat))
cnv_cols<-grep('proportion_scaled_cnv_chr',names(seu@meta.data),value=T)
cnvs<-seu@meta.data[,cnv_cols]
seu$cnv_avg<-rowMeans(cnvs)
q10<- 0.1
seu$malignant<-ifelse(seu$cnv_avg > q10, 'malignant','non-malignant')

# add CNV metrics
cnv_cols<-grep('proportion_scaled_cnv_chr',names(seu@meta.data),value=T)
cnvs<-seu@meta.data[,cnv_cols]
seu$proportion_scaled_cnv_avg<-rowMeans(cnvs)

cnv_cols<-grep('proportion_cnv_chr',names(seu@meta.data),value=T)
cnvs<-seu@meta.data[,cnv_cols]
seu$proportion_cnv_avg<-rowMeans(cnvs)

cnv_cols<-grep('has_cnv_chr',names(seu@meta.data),value=T)
cnvs<-seu@meta.data[,cnv_cols]
seu$has_cnv_avg<-rowMeans(cnvs)


# save objects
saveRDS(seu,paste0("/home/ubuntu/vol/inferCNV_subcluster_", pat, '/', pat, "_cnv.rds"))


# results
pdf(paste0("/home/ubuntu/vol/inferCNV_subcluster_", pat, '/', pat, "_cnv_cut=0.1.pdf"))
textplot(addmargins(table(seu$immune,seu$malignant)),cex=1.2,halign='left')
textplot(table(seu$celltype_bped_fine,seu$malignant),cex=0.9,halign='left')
hist(seu$cnv_avg,breaks=100, main='Average absolute CNV level; all cells',
     xlab = 'Average absolute CNV proportion')
abline(v=q10,col="red")
DimPlot(seu, reduction='pca',group.by='malignant')
DimPlot(seu, reduction='umap',group.by='malignant')
DimPlot(seu, reduction='umap',group.by='immune')

FeaturePlot(seu, features =  c("proportion_scaled_cnv_avg"), 
            min.cutoff = "q05",max.cutoff = 'q95',order=T, raster = T)+ 
  scale_color_viridis(direction = -1)
FeaturePlot(seu, features =  c("proportion_cnv_avg"), 
            min.cutoff = "q05",max.cutoff = 'q95',order=T, raster = T)+ 
  scale_color_viridis(direction = -1)

FeaturePlot(seu, features = c('MLANA',"MITF", 'PMEL', "MKI67"), min.cutoff = "q05",
            max.cutoff = 'q95',order=T, raster = T)
FeaturePlot(seu, features = c('PTPRC',"CD8A", 'CD68','CD79A' ), min.cutoff = "q05",
            max.cutoff = 'q95',order=T, raster = T)
FeaturePlot(seu, features = c('ALB',"KRT15", 'COL1A1', "VWF"), min.cutoff = "q05",
            max.cutoff = 'q95',order=T, raster = T)

DimPlot(seu, reduction = "umap",label = T,group.by = 'celltype_bped_fine',repel = T,
        label.size = 2.5,raster = T,shuffle = T) + 
  guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
  theme(legend.position="none")

DimPlot(seu, reduction = "umap",label = T,group.by = 'celltype_hpca_main',repel = T,
        label.size = 2.5,raster = T,shuffle = T) + 
  guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
  theme(legend.position="none")

dev.off()

system(paste0("aws s3 sync /home/ubuntu/vol/inferCNV_subcluster_", pat, "/ s3://melanoma-ribas/inferCNV/inferCNV_subcluster_", pat,"/ --exclude '*._*' --exclude '*DS_S*' --quiet"))

print(paste("End:",Sys.time()))
