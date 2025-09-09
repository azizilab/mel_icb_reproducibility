#!/usr/bin/env Rscript

#### Seurat integration of CUIMC melanoma cohort
#### Author: Jana Biermann, PhD

print(Sys.time())

library(dplyr)
library(Seurat)
library(gplots)
library(ggplot2)
library(RColorBrewer)
library(viridis)
'%notin%' <- Negate('%in%')

label <- 'melanoma'
d<-30

seu<-readRDS('data/merged/data_melanoma_merged.rds')
DefaultAssay(seu)<-'RNA'
DietSeurat(seu)

obj.list<-SplitObject(seu,split.by = 'sample')

obj.list <- lapply(X = obj.list, FUN = function(x) {
  x <- x[grep("^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA|^MT-", rownames(x), value = T, invert = T), ]
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose=F)
})

features <- SelectIntegrationFeatures(object.list = obj.list)
obj.list <- lapply(X = obj.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

anchors <- FindIntegrationAnchors(object.list = obj.list,
                                  reduction = "rpca", 
                                  dims = 1:d)
seu <- IntegrateData(anchorset = anchors, dims = 1:d)
rm(anchors, obj.list)

# Seurat workflow
seu <- ScaleData(seu)
seu <- RunPCA(seu, npcs = 75)
seu <- RunUMAP(seu, dims = 1:d)
seu <- FindNeighbors(seu, dims = 1:d)
seu <- FindClusters(seu)


# Save object
ifelse(!dir.exists(file.path('data/integrated/')), 
       dir.create(file.path('data/integrated/')), FALSE)
saveRDS(seu, file = paste0('data/integrated/data_melanoma_integrated_RPCA_v4.rds'))

### stats
stats<-as.data.frame(matrix(data=NA,nrow = 1, ncol = 11))
colnames(stats)<-c('sample','n_features','n_cells','median_features','median_counts')
rownames(stats)<-label
stats$sample<-label
stats$n_features<-dim(seu@assays$integrated@data)[1]
stats$n_cells<-dim(seu@assays$integrated@data)[2]
stats$median_features<-round(median(seu@meta.data$nFeature_RNA))
stats$median_counts<-round(median(seu@meta.data$nCount_RNA))

pdf(paste0("data/integrated/plots_melanoma_integrated_RPCA_v4.pdf"))
textplot(t(stats),cex=1.2,halign='left')
ElbowPlot(seu,ndims = 75)
print(DimPlot(seu, reduction = "pca",group.by = 'ident',raster = T,shuffle = T)+coord_fixed())
print(DimPlot(seu, reduction = "pca",group.by = 'patient',raster = T,shuffle = T)+coord_fixed())
print(DimPlot(seu, reduction = "pca",group.by = 'tissue',raster = T,shuffle = T)+coord_fixed())
print(DimPlot(seu, reduction = "pca",group.by = 'treated',raster = T,shuffle = T)+coord_fixed())
print(DimPlot(seu, label = T,group.by = 'ident',raster = T,shuffle = T)+coord_fixed())
print(DimPlot(seu, label = F,group.by = 'patient',raster = T,shuffle = T)+coord_fixed())
print(DimPlot(seu, label = F,group.by = 'sample',raster = T,shuffle = T)+coord_fixed())
print(DimPlot(seu, label = F,group.by = 'tissue',raster = T,shuffle = T)+coord_fixed())
print(DimPlot(seu, label = T,group.by = 'cell_type_fine',shuffle = T,raster = T,
              repel = T,label.size = 2)+NoLegend())
print(DimPlot(seu, label = T,group.by = 'cell_type_main',shuffle = T,raster = T,
              repel = T)+coord_fixed())
print(DimPlot(seu, label = F,group.by = 'cell_type_main',shuffle = T,raster = T,repel = T,
              split.by = 'cell_type_main',ncol = 3))
print(DimPlot(seu, label = T,group.by = 'cell_type_int',shuffle = T,raster = T,
              repel = T)+coord_fixed())
print(DimPlot(seu, label = F,group.by = c('treated','responder','mucosal'),ncol = 2,
              raster = T,shuffle = T)&coord_fixed())
print(FeatureScatter(seu,'G2M.Score','S.Score',group.by = 'cell_cycle',shuffle = T,raster = T))
print(DimPlot(seu, label = F,group.by = c('cell_cycle','malignant','clone_size'),
              raster = T,shuffle = T,ncol = 2)&coord_fixed())
DimPlot(seu, label = T,group.by = 'celltype_bped_fine',repel = T,raster = T,shuffle = T,
        label.size = 2.5) + 
  guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
  theme(legend.text=element_text(size=6))
DimPlot(seu, label = T,group.by = 'celltype_hpca_main',repel = T,raster = T,shuffle = T,
        label.size = 2.5) + 
  guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
  theme(legend.text=element_text(size=6))

FeaturePlot(seu,features = c('rna_MLANA','rna_PMEL','rna_VWF','rna_COL1A1'),
            order = T,raster = T,max.cutoff = 'q95',min.cutoff = 'q05')& 
  scale_colour_gradientn(colours = rev(brewer.pal(11,'Spectral')))
FeaturePlot(seu, features = c('rna_PTPRC',"rna_CD8A", 'rna_CD68', "rna_MKI67"), 
            min.cutoff = "q05",max.cutoff = 'q95',order=T, raster = T)& 
  scale_colour_gradientn(colours = rev(brewer.pal(11,'Spectral')))
print(FeaturePlot(seu, features = c("rna_MLANA",'rna_PMEL','rna_MITF','rna_AXL'),
                  min.cutoff = "q05",max.cutoff = "q95", order = T,raster = T)& 
        scale_colour_gradientn(colours = rev(brewer.pal(11,'Spectral'))))

FeaturePlot(seu, features = c("nCount_RNA", "nFeature_RNA", 'percent.mt', "doublet_scores"), 
            min.cutoff = "q05",max.cutoff = "q95", order = T,raster = T)& 
  scale_colour_gradientn(colours = rev(brewer.pal(11,'Spectral')))
FeaturePlot(seu, features = c("percent.rps", "percent.rpl", 'G2M.Score', "S.Score"), 
            min.cutoff = "q05",max.cutoff = "q95", order = T,raster = T)& 
  scale_colour_gradientn(colours = rev(brewer.pal(11,'Spectral')))
FeaturePlot(seu, features = c('proportion_scaled_cnv_avg'), order=T,raster = T,
            pt.size=0.5, min.cutoff = "q05",max.cutoff = "q95")& 
  scale_colour_gradientn(colours = rev(brewer.pal(11,'Spectral')))
VlnPlot(seu,features = c("nCount_RNA", "nFeature_RNA", 'proportion_scaled_cnv_avg', "doublet_scores"),
        pt.size = 0,ncol = 2)&
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,size=5))
dev.off()


print(Sys.time())