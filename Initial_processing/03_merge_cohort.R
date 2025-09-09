#!/usr/bin/env Rscript

#### Merge individual samples
#### Author: Jana Biermann, PhD

library(dplyr)
library(Seurat)
library(gplots)
library(ggplot2)
library(ggrastr)
library(viridis)
'%notin%' <- Negate('%in%')

label <- 'melanoma'
s3folder<-'melanoma-ribas'

pats<-c('F01_on','F01_pre','F02_on','F02_pre','F03_post1_on2','F03_post1_pre2','F04_pre',
        'F05_pre','F06_post1_pre2','F07_post1_pre2','F08_post','F09_post','F10_post',
        'F12_post','F12_pre','F15_pre','F16_post1_pre2','F16_pre','F17_post',
        'F18_post','F20_post1_pre2','F22_post','F23_post','F24_post','F25_post',
        'F26_post','F27_post1_pre2','F28_post1_pre2','F29_post1_pre2','F30_post',
        'F31_post','R204_pre','R294_on','R308_pre','R310_on1','R310_on2','R310_pre',
        'R319_on','R319_pre','R328_on','R329_on','R334_pre','R354_pre')

system(paste0("aws s3 sync s3://",s3folder,"/inferCNV/inferCNV_subcluster_", pats[1], 
              "/ data/",pats[1],"/ --exclude '*' --include '*_cnv.rds' "))
seu<- readRDS(paste0("data/",pats[1],"/",pats[1],"_cnv.rds"))

for(pat in pats[2:length(pats)]){
  print(paste('Processing:', pat))
  system(paste0("aws s3 sync s3://",s3folder,"/inferCNV/inferCNV_subcluster_", pat, 
                "/ data/",pat,"/ --exclude '*' --include '*_cnv.rds' "))
  tmp<- readRDS(paste0("data/",pat,"/",pat,"_cnv.rds"))
  seu<-merge(seu,tmp)
}

# Seurat workflow
seu <- NormalizeData(seu) %>% FindVariableFeatures()
seu <- ScaleData(seu)
seu <- RunPCA(seu, npcs = 100)
seu <- RunUMAP(seu, dims = 1:30)
seu <- FindNeighbors(seu, dims = 1:30)
seu <- FindClusters(seu)

# Cell cycle
seu <- CellCycleScoring(seu, s.features = cc.genes.updated.2019$s.genes, 
                        g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = F)
seu$cell_cycle <- ifelse(seu$G2M.Score > 0.1 | seu$S.Score > 0.1, 'cycling', 'non-cycling')

# Save merged object
if(!dir.exists(file.path('data/merged/'))){dir.create(file.path('data/merged/'))}
saveRDS(seu, file = paste0('data/merged/data_',label,'_merged.rds'))


### Plots
stats<-as.data.frame(matrix(data=NA,nrow = 1, ncol = 11))
colnames(stats)<-c('sample','n_features','n_cells','median_features','median_counts')
rownames(stats)<-label
stats$sample<-label
stats$n_features<-dim(seu@assays$RNA@data)[1]
stats$n_cells<-dim(seu@assays$RNA@data)[2]
stats$median_features<-round(median(seu@meta.data$nFeature_RNA))
stats$median_counts<-round(median(seu@meta.data$nCount_RNA))

pdf(file = paste0("data/merged/plots_",label,"_merged.pdf"))
textplot(t(stats),cex=1.2,halign='left')
ElbowPlot(seu,ndims = 100)
print(DimPlot(seu, reduction = "pca",group.by = 'ident',raster = T,shuffle = T,pt.size = 0.01))
print(DimPlot(seu, reduction = "pca",group.by = 'sample',raster = T,shuffle = T,pt.size = 0.01))
print(DimPlot(seu, reduction = "pca",group.by = 'patient',raster = T,shuffle = T,pt.size = 0.01))
print(DimPlot(seu, reduction = "pca",group.by = 'tissue',raster = T,shuffle = T,pt.size = 0.01))
print(DimPlot(seu, reduction = "pca",group.by = 'sample_group',raster = T,shuffle = T,pt.size = 0.01))
print(DimPlot(seu, reduction = "umap",label = T,group.by = 'ident',raster = T,shuffle = T,pt.size = 0.01))
print(DimPlot(seu, reduction = "umap",label = F,group.by = 'sample',raster = T,shuffle = T,pt.size = 0.01))
print(DimPlot(seu, reduction = "umap",label = F,group.by = 'patient',raster = T,shuffle = T,pt.size = 0.01))
print(DimPlot(seu, reduction = "umap",label = F,group.by = 'tissue',raster = T,shuffle = T,pt.size = 0.01))
print(DimPlot(seu, reduction = "umap",label = T,group.by = 'sample_group',raster = T,shuffle = T,pt.size = 0.01))
print(DimPlot(seu, reduction = "umap",label = F,group.by = 'orig.ident',raster = T,shuffle = T,pt.size = 0.01))
print(FeatureScatter(seu,'G2M.Score','S.Score',group.by = 'cell_cycle',pt.size = 0.01,shuffle = T,raster = T))
print(DimPlot(seu, reduction = "umap",label = T,group.by = 'cell_cycle',pt.size = 0.01,raster = T,shuffle = T))
print(DimPlot(seu, reduction = "umap",label = T,group.by = 'malignant',raster = T,shuffle = T))
DimPlot(seu, label = T,group.by = 'clone_size',raster = T,shuffle = T,
        label.size = 3,repel = T,pt.size=0.1)

DimPlot(seu, reduction = "umap",label = T,group.by = 'celltype_bped_fine',repel = T,
        label.size = 2.5,pt.size = 0.01) + 
  ggtitle('celltype_bped_fine') +
  guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
  theme(legend.text=element_text(size=6))
DimPlot(seu, reduction = "umap",label = T,group.by = 'celltype_hpca_main',repel = T,
        label.size = 2.5,raster = T,shuffle = T,pt.size = 0.01) + 
  ggtitle('celltype_hpca_main') +
  guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
  theme(legend.text=element_text(size=6))

FeaturePlot(seu, features = c('rna_MLANA',"rna_MITF", 'rna_PMEL', "rna_AXL"), min.cutoff = "q05",
            max.cutoff = 'q95',order=T, raster = T)

FeaturePlot(seu, features = c('rna_PTPRC',"rna_CD8A", 'rna_CD68', "rna_MKI67"), min.cutoff = "q05",
            max.cutoff = 'q95',order=T, raster = T)

FeaturePlot(seu, features = c('rna_ALB',"rna_CYP2E1", 'rna_COL1A1', "rna_VWF"), min.cutoff = "q05",
            max.cutoff = 'q95',order=T, raster = T)

print(FeaturePlot(seu, features = c("nCount_RNA", "nFeature_RNA", "percent.mt",'proportion_scaled_cnv_avg'), 
                  order=T,raster = T,pt.size = 0.01))
FeaturePlot(seu, features =  c("proportion_scaled_cnv_avg"), 
            min.cutoff = "q05",max.cutoff = 'q95',order=T, raster = T)+ 
  scale_color_viridis(direction = -1)
dev.off()

system(paste0("aws s3 sync data/merged/ s3://",s3folder,"/Seurat/merged/ --exclude '*' --include '*",label,"*' --quiet"))
