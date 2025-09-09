#!/usr/bin/env Rscript

#### Main cell type annotation of integrated and merged Seurat object
#### Author: Jana Biermann, PhD

library(dplyr)
library(Seurat)
library(ggplot2)
library(gplots)
library(viridis)
library(RColorBrewer)
library(ggpubr)
'%notin%' <- Negate('%in%')

#setwd('~/Documents/melanoma')
label<-'melanoma'
directory<-paste0("data/annotation/main/")
directory_tum<-paste0("data/annotation/tumor/")

ifelse(!dir.exists(file.path(directory)), 
       dir.create(file.path(directory)), FALSE)
ifelse(!dir.exists(file.path(directory_tum)), 
       dir.create(file.path(directory_tum)), FALSE)

#### Merged object cluster overview ####
# Read in object
seu<-readRDS('data/merged/data_melanoma_merged.rds')

pdf('data/merged/plots_melanoma_merged_clusters.pdf')
print(DimPlot(seu, reduction = "umap",label = T,group.by = 'RNA_snn_res.0.3',raster = T,
              shuffle = T,pt.size = 0.01))
DimPlot(seu, reduction = "umap",label = F,group.by = 'sample',raster = T,
        shuffle = T,pt.size = 0.01)+ 
  guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
  theme(legend.text=element_text(size=6))

DimPlot(seu, reduction = "umap",label = T,group.by = 'sample',raster = T,
        shuffle = T,pt.size = 0.01)+ NoLegend()

ggplot(seu@meta.data, aes(x=RNA_snn_res.0.3, fill=sample))+
  geom_bar(position="fill",size=0.3,color='black') +
  theme_classic() + 
  ggtitle('Sample distribution among clusters')+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggplot(seu@meta.data, aes(x=RNA_snn_res.0.3, fill=patient))+
  geom_bar(position="fill",size=0.3,color='black') +
  theme_classic() + 
  ggtitle('Patient distribution among clusters')+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggplot(seu@meta.data, aes(x=RNA_snn_res.0.3, fill=malignant))+
  geom_bar(position="fill",size=0.3,color='black') +
  theme_classic() + 
  ggtitle('malignant distribution among clusters')+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

print(DimPlot(seu, reduction = "umap",label = F,group.by = 'orig.ident',raster = T,
              shuffle = T,pt.size = 0.01)+ 
        guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
        theme(legend.text=element_text(size=6)))
print(DimPlot(seu, reduction = "umap",label = F,group.by = 'tissue',raster = T,
              shuffle = T,pt.size = 0.01))
print(DimPlot(seu, reduction = "umap",label = T,group.by = 'sample_group',raster = T,
              shuffle = T,pt.size = 0.01)+ 
        guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
        theme(legend.text=element_text(size=6)))
#print(DimPlot(seu, reduction = "umap",label = T,group.by = 'cell_cycle',pt.size = 0.01,raster = T,shuffle = T))
DimPlot(seu, label = T,group.by = 'clone_size',raster = T,shuffle = T,
        label.size = 3,repel = T,pt.size=0.1)
DimPlot(seu, label = T,group.by = 'malignant',raster = T,shuffle = T,
        label.size = 3,repel = T,pt.size=0.1)
FeaturePlot(seu, features = 'proportion_scaled_cnv_avg',min.cutoff = "q05", 
            max.cutoff = "q95",order=T,raster=T,pt.size = 0.1)& 
  scale_colour_gradientn(colours = rev(brewer.pal(11,'Spectral')))

DimPlot(seu, reduction = "umap",label = T,group.by = 'celltype_bped_fine',repel = T,
        label.size = 2.5,pt.size = 0.01) + 
  guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
  theme(legend.text=element_text(size=6))
dev.off()


#### Annotation main (integrated object) ####
# Read in object
seu<-readRDS(paste0('data/integrated/data_',label,'_integrated_RPCA.rds'))

# Main cell type assignment after manual annotation based on DGE
seu$cell_type_main <- ifelse(seu$integrated_snn_res.0.8 %in% c(5,24,25,32,35), 'Myeloid cells', 'NA')
seu$cell_type_main <- ifelse(seu$integrated_snn_res.0.8 %in% c(0,1,2,3,4,6,7,8,11,14,15,16,18,19,21,26,27,29,37), 'Tumor cells', seu$cell_type_main)
seu$cell_type_main <- ifelse(seu$integrated_snn_res.0.8 %in% c(10,12,28,31,39), 'Stromal cells', seu$cell_type_main)
seu$cell_type_main <- ifelse(seu$integrated_snn_res.0.8 %in% c(9,13,17,30), 'T/NK cells', seu$cell_type_main)
seu$cell_type_main <- ifelse(seu$integrated_snn_res.0.8 %in% c(20,23,34), 'B/Plasma cells', seu$cell_type_main)
seu$cell_type_main <- ifelse(seu$integrated_snn_res.0.8 %in% c(22,36), 'Endothelial cells', seu$cell_type_main)
seu$cell_type_main <- ifelse(seu$integrated_snn_res.0.8 %in% c(33), 'Epithelial cells', seu$cell_type_main)
seu$cell_type_main <- ifelse(seu$integrated_snn_res.0.8 %in% c(38), 'Mast cells', seu$cell_type_main)


#### Update metadata in merged object and plot cell type distribution ####
# Read in object
seu<-readRDS('data/merged/data_melanoma_merged.rds')

clin<-read.csv('data/clin.csv',na.strings = '')
name_clin<-names(clin)[2:ncol(clin)]
seu@meta.data[,name_clin]<-NULL

seu$barcode_all<-row.names(seu@meta.data)
seu@meta.data<-left_join(seu@meta.data,clin,by='sample')
row.names(seu@meta.data)<-seu$barcode_all
saveRDS(seu,'data/merged/data_melanoma_merged_v2.rds')

# Read in object
seu<-readRDS('data/merged/data_melanoma_merged_v2.rds')

ct<-read.csv('data/annotation/main/data_melanoma_celltype_main_nontum.csv')
ct$X<-NULL
seu@meta.data<-left_join(seu@meta.data,ct,by='barcode_pat')
row.names(seu@meta.data)<-seu$barcode_all
seu$cell_type_main<-ifelse(is.na(seu$cell_type_main)==T,'Tumor cells',seu$cell_type_main)
seu$melanoma<-'melanoma'

df <- data.frame(sample = seu$sample, 
                 cell_type_main = seu$cell_type_main,
                 group=seu$group) 

# cell_type_main
df_main = df %>%
  group_by(sample, cell_type_main, group) %>%
  tally() %>%
  group_by(sample) %>%
  mutate(freq = n/sum(n))

pdf('data/merged/plots_melanoma_merged_v2.pdf')
print(DimPlot(seu, reduction = "umap",label = T,group.by = 'RNA_snn_res.0.3',raster = T,
              shuffle = T,pt.size = 0.01))
DimPlot(seu, reduction = "umap",label = F,group.by = 'sample',raster = T,
        shuffle = T,pt.size = 0.01)+ 
  guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
  theme(legend.text=element_text(size=6))

DimPlot(seu, reduction = "umap",label = T,group.by = 'sample',raster = T,
        shuffle = T,pt.size = 0.01)+ NoLegend()

DimPlot(seu, reduction = "umap",label = T,group.by = 'cell_type_main',raster = T,
        shuffle = T,pt.size = 0.01)

ggplot(seu@meta.data, aes(x=sample, fill=cell_type_main))+
  geom_bar(position="fill",size=0.3,color='black') +
  theme_classic() + 
  ggtitle('Cell type distribution')+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggplot(seu@meta.data, aes(x=group, fill=cell_type_main))+
  geom_bar(position="fill",size=0.3,color='black') +
  theme_classic() + 
  ggtitle('Cell type distribution')+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggplot(seu@meta.data, aes(x=tissue, fill=cell_type_main))+
  geom_bar(position="fill",size=0.3,color='black') +
  theme_classic() + 
  ggtitle('Cell type distribution')+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggplot(seu@meta.data, aes(x=sample, fill=cell_type_main))+
  geom_bar(position="fill",size=0.3,color='black') +
  theme_classic() + 
  ggtitle('Cell type distribution')+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  facet_grid(group~tissue)

ggplot(seu@meta.data, aes(x=melanoma, fill=cell_type_main))+
  geom_bar(position="fill",size=0.3,color='black') +
  theme_classic() + 
  ggtitle('Cell type distribution')+
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        strip.text.x=element_text(size=6,angle = 90))+
  facet_grid(tissue~group)

ggplot(seu@meta.data, aes(x=tissue, fill=cell_type_main))+
  geom_bar(position="fill",size=0.3,color='black') +
  theme_classic() + 
  ggtitle('Cell type distribution')+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        strip.text.x=element_text(size=5))+
  facet_grid(~group)

ggplot(seu@meta.data, aes(x=group, fill=cell_type_main))+
  geom_bar(position="fill",size=0.3,color='black') +
  theme_classic() + 
  ggtitle('Cell type distribution')+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  facet_grid(~tissue)

ggboxplot(df_main, x = 'cell_type_main', y = 'freq', color = 'group', add = 'jitter', 
          order = sort(unique(df_main$cell_type_main))) + ylim(0, 1) + 
  #stat_compare_means(aes(group = group),label = 'p.format',method = 'wilcox.test', size = 3) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.text = element_text(size=6)) +
  ggtitle('cell_type_main') #+ color_palette(palette = colSeq)

print(DimPlot(seu, reduction = "umap",label = F,group.by = 'tissue',raster = T,
              shuffle = T,pt.size = 0.01))
print(DimPlot(seu, reduction = "umap",label = T,group.by = 'group',raster = T,
              shuffle = T,pt.size = 0.01)+ 
        guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
        theme(legend.text=element_text(size=6)))
#print(DimPlot(seu, reduction = "umap",label = T,group.by = 'cell_cycle',pt.size = 0.01,raster = T,shuffle = T))
DimPlot(seu, label = T,group.by = 'clone_size',raster = T,shuffle = T,
        label.size = 3,repel = T,pt.size=0.1)
DimPlot(seu, label = T,group.by = 'malignant',raster = T,shuffle = T,
        label.size = 3,repel = T,pt.size=0.1)
FeaturePlot(seu, features = 'proportion_scaled_cnv_avg',min.cutoff = "q05", 
            max.cutoff = "q95",order=T,raster=T,pt.size = 0.1)& 
  scale_colour_gradientn(colours = rev(brewer.pal(11,'Spectral')))

DimPlot(seu, reduction = "umap",label = T,group.by = 'celltype_bped_fine',repel = T,
        label.size = 2.5,pt.size = 0.01) + 
  guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
  theme(legend.text=element_text(size=6))
dev.off()


#### Subset to non-tumor cells for further annotation ####
nontum <- subset(seu, RNA_snn_res.0.3 %in% c(0,1,9,25,32,33,36,39,41))
DimPlot(nontum, reduction = "umap",label = T,group.by = 'sample',raster = T,
        shuffle = T,pt.size = 0.01)+ NoLegend()

# Seurat workflow
nontum <- NormalizeData(nontum) %>% ScaleData()
nontum<-FindVariableFeatures(nontum)
nontum <- RunPCA(nontum, npcs = 50)
ElbowPlot(nontum,ndims = 50)
nontum <- RunUMAP(nontum, dims = 1:30)
nontum <- FindNeighbors(nontum, dims = 1:30)
nontum <- FindClusters(nontum,resolution = 0.1)

saveRDS(nontum,'data/integrated/data_nontum_merged.rds')

### Remove tumor clusters
nontum<-readRDS('data/integrated/data_nontum_merged.rds')
nontum <- subset(nontum, RNA_snn_res.0.1 %notin% c(3,9,11))

nontum <- NormalizeData(nontum) %>% ScaleData()
nontum<-FindVariableFeatures(nontum)
nontum <- RunPCA(nontum, npcs = 50)
ElbowPlot(nontum,ndims = 50)
nontum <- RunUMAP(nontum, dims = 1:30)
nontum <- FindNeighbors(nontum, dims = 1:30)
nontum <- FindClusters(nontum,resolution = 0.3)

nontum <- CellCycleScoring(nontum, s.features = cc.genes.updated.2019$s.genes,
                           g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = F)
nontum$cell_cycle <- ifelse(nontum$G2M.Score > 0.1 | nontum$S.Score > 0.1, 'cycling',
                            'non-cycling')

# Differential gene expression 
markers <- FindAllMarkers(nontum, only.pos = TRUE, assay = "RNA", min.pct = 0.25,
                          logfc.threshold = 0.25)
write.csv(markers, paste0(directory,'markers_nontum_res0.3.csv'), row.names = F)

# Main cell type assignment after manual annotation based on DGE
nontum$cell_type_main <- ifelse(nontum$RNA_snn_res.0.3 %in% c(0,5,10,11), 'Myeloid cells', 'NA')
nontum$cell_type_main <- ifelse(nontum$RNA_snn_res.0.3 %in% c(2,9,15,17), 'Stromal cells', nontum$cell_type_main)
nontum$cell_type_main <- ifelse(nontum$RNA_snn_res.0.3 %in% c(1,3,8), 'T/NK cells', nontum$cell_type_main)
nontum$cell_type_main <- ifelse(nontum$RNA_snn_res.0.3 %in% c(4,7), 'B/Plasma cells', nontum$cell_type_main)
nontum$cell_type_main <- ifelse(nontum$RNA_snn_res.0.3 %in% c(6,14), 'Endothelial cells', nontum$cell_type_main)
nontum$cell_type_main <- ifelse(nontum$RNA_snn_res.0.3 %in% c(12), 'Epithelial cells', nontum$cell_type_main)
nontum$cell_type_main <- ifelse(nontum$RNA_snn_res.0.3 %in% c(16), 'Mast cells', nontum$cell_type_main)
nontum$cell_type_main <- ifelse(nontum$RNA_snn_res.0.3 %in% c(13), ' pDCs', nontum$cell_type_main)


# Save object
saveRDS(nontum,'data/merged/data_nontum_merged_v2.rds')
celltype <- nontum@meta.data %>% dplyr::select('barcode_pat', 'cell_type_main')
write.csv(celltype, paste0(directory,'data_',label,'_celltype_main_nontum.csv'))
write.csv(nontum$barcode_pat, paste0(directory,'data_nontum_merged_v2_barcodes.csv'))

# Plots
pdf('data/merged/plots_nontum_merged_clusters_v2.pdf')
DimPlot(nontum, reduction = "umap",label = T,group.by = 'ident',raster = T,
        shuffle = T,pt.size = 0.01)
DimPlot(nontum, reduction = "umap",label = F,group.by = 'sample',raster = T,
        shuffle = T,pt.size = 0.01)+ 
  guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
  theme(legend.text=element_text(size=6))

DimPlot(nontum, reduction = "umap",label = T,group.by = 'sample',raster = T,
        shuffle = T,pt.size = 0.01)+ NoLegend()
DimPlot(nontum, reduction = "umap",label = T,group.by = 'cell_type_main',raster = T,
        shuffle = T,pt.size = 0.01)

ggplot(nontum@meta.data, aes(x=RNA_snn_res.0.3, fill=patient))+
  geom_bar(position="fill",size=0.3,color='black') +
  theme_classic() + 
  ggtitle('Patient distribution among clusters')+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggplot(nontum@meta.data, aes(x=RNA_snn_res.0.3, fill=malignant))+
  geom_bar(position="fill",size=0.3,color='black') +
  theme_classic() + 
  ggtitle('malignant distribution among clusters')+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

print(DimPlot(nontum, reduction = "umap",label = F,group.by = 'orig.ident',raster = T,
              shuffle = T,pt.size = 0.01)+ 
        guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
        theme(legend.text=element_text(size=6)))
print(DimPlot(nontum, reduction = "umap",label = F,group.by = 'tissue',raster = T,
              shuffle = T,pt.size = 0.01))
print(DimPlot(nontum, reduction = "umap",label = T,group.by = 'sample_group',raster = T,
              shuffle = T,pt.size = 0.01)+ 
        guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
        theme(legend.text=element_text(size=6)))
print(DimPlot(nontum, reduction = "umap",label = T,group.by = 'cell_cycle',pt.size = 0.01,raster = T,shuffle = T))
DimPlot(nontum, label = T,group.by = 'clone_size',raster = T,shuffle = T,
        label.size = 3,repel = T,pt.size=0.1)
DimPlot(nontum, label = T,group.by = 'malignant',raster = T,shuffle = T,
        label.size = 3,repel = T,pt.size=0.1)
FeaturePlot(nontum, features = 'proportion_scaled_cnv_avg',min.cutoff = "q05", 
            max.cutoff = "q95",order=T,raster=T,pt.size = 0.1)& 
  scale_colour_gradientn(colours = rev(brewer.pal(11,'Spectral')))

FeaturePlot(nontum, features = c("nCount_RNA", "nFeature_RNA", 'percent.mt', "doublet_scores"), 
            min.cutoff = "q05",max.cutoff = "q95", order = T,raster = T)& 
  scale_colour_gradientn(colours = rev(brewer.pal(11,'Spectral')))
FeaturePlot(nontum, features = c("percent.rps", "percent.rpl", 'G2M.Score', "S.Score"), 
            min.cutoff = "q05",max.cutoff = "q95", order = T,raster = T)& 
  scale_colour_gradientn(colours = rev(brewer.pal(11,'Spectral')))

DotPlot(nontum,assay = 'RNA',features = c('XBP1','MS4A1','CD19', #plasma/b
                                          'CD3E','TCF7', 'CD8A', #t cell
                                          'VCAN','LYZ','MERTK','F13A1','CLEC9A','CD1C', #myeloid
                                          'ZFAT','IL3RA','PLD4', #pDC
                                          'CPA3','TPSAB1', #mast
                                          'COL1A1','COL3A1','LAMA2', #FB
                                          'KRT15','KRT17','CCL19', #epithelial
                                          'VWF',
                                          'MKI67','TOP2A','PTPRC','MLANA','PMEL'
),
group.by = 'RNA_snn_res.0.3', dot.scale = 5)+ scale_color_viridis() + coord_flip()+
  theme(axis.text.y = element_text(size=7))

DotPlot(nontum,assay = 'RNA',features = c('XBP1','MS4A1','CD19', #plasma/b
                                          'CD3E','TCF7', 'CD8A', #t cell
                                          'VCAN','LYZ','MERTK','F13A1','CLEC9A','CD1C', #myeloid
                                          'ZFAT','IL3RA','PLD4', #pDC
                                          'CPA3','TPSAB1', #mast
                                          'COL1A1','COL3A1','LAMA2', #FB
                                          'KRT15','KRT17','CCL19', #epithelial
                                          'VWF',
                                          'MKI67','TOP2A','PTPRC','MLANA','PMEL'
),
group.by = 'cell_type_main', dot.scale = 5)+ scale_color_viridis() + coord_flip()+
  theme(axis.text.y = element_text(size=7))+RotatedAxis()

DimPlot(nontum, reduction = "umap",label = T,group.by = 'celltype_bped_fine',repel = T,
        label.size = 2.5,pt.size = 0.01,raster = T,shuffle = T) + 
  guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
  theme(legend.text=element_text(size=6))
dev.off()



#### Subset to tumor cells ####
tum <- subset(seu, barcode_pat %notin% nontum$barcode_pat)
DimPlot(tum, reduction = "umap",label = T,group.by = 'sample',raster = T,
        shuffle = T,pt.size = 0.01)+ NoLegend()

# Seurat workflow
tum <- NormalizeData(tum) %>% ScaleData()
tum<-FindVariableFeatures(tum)
tum <- RunPCA(tum, npcs = 100)
ElbowPlot(tum,ndims = 100)
tum <- RunUMAP(tum, dims = 1:50)
tum <- FindNeighbors(tum, dims = 1:50)
tum <- FindClusters(tum,resolution = 0.4)

tum <- CellCycleScoring(tum, s.features = cc.genes.updated.2019$s.genes,
                        g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = F)
tum$cell_cycle <- ifelse(tum$G2M.Score > 0.1 | tum$S.Score > 0.1, 'cycling',
                         'non-cycling')

tum$tissue<-ifelse(is.na(tum$tissue)==T,'skin',tum$tissue)

#saveRDS(tum,'data/merged/data_tum_merged.rds')
tum<-readRDS('data/merged/data_tum_merged.rds')

pdf('data/merged/plots_tum_merged_clusters.pdf')
DimPlot(tum, reduction = "umap",label = T,group.by = 'ident',raster = T,
        shuffle = T,pt.size = 0.01)
DimPlot(tum, reduction = "umap",label = F,group.by = 'sample',raster = T,
        shuffle = T,pt.size = 0.01)+ 
  guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
  theme(legend.text=element_text(size=6))

DimPlot(tum, reduction = "umap",label = T,group.by = 'patient',raster = T,
        shuffle = T,pt.size = 0.01)

ggplot(tum@meta.data, aes(x=RNA_snn_res.0.4, fill=sample))+
  geom_bar(position="fill",size=0.3,color='black') +
  theme_classic() + 
  ggtitle('Sample distribution among clusters')+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggplot(tum@meta.data, aes(x=RNA_snn_res.0.4, fill=patient))+
  geom_bar(position="fill",size=0.3,color='black') +
  theme_classic() + 
  ggtitle('Patient distribution among clusters')+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggplot(tum@meta.data, aes(x=RNA_snn_res.0.4, fill=malignant))+
  geom_bar(position="fill",size=0.3,color='black') +
  theme_classic() + 
  ggtitle('malignant distribution among clusters')+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

print(DimPlot(tum, reduction = "umap",label = F,group.by = 'orig.ident',raster = T,
              shuffle = T,pt.size = 0.01)+ 
        guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
        theme(legend.text=element_text(size=6)))
print(DimPlot(tum, reduction = "umap",label = F,group.by = 'tissue',raster = T,
              shuffle = T,pt.size = 0.01))
print(DimPlot(tum, reduction = "umap",label = T,group.by = 'sample_group',raster = T,
              shuffle = T,pt.size = 0.01)+ 
        guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
        theme(legend.text=element_text(size=6)))
print(DimPlot(tum, reduction = "umap",label = T,group.by = 'cell_cycle',
              pt.size = 0.01,raster = T,shuffle = T))
DimPlot(tum, label = T,group.by = 'clone_size',raster = T,shuffle = T,
        label.size = 3,repel = T,pt.size=0.1)
DimPlot(tum, label = T,group.by = 'malignant',raster = T,shuffle = T,
        label.size = 3,repel = T,pt.size=0.1)
FeaturePlot(tum, features = 'proportion_scaled_cnv_avg',min.cutoff = "q05", 
            max.cutoff = "q95",order=T,raster=T,pt.size = 0.1)& 
  scale_colour_gradientn(colours = rev(brewer.pal(11,'Spectral')))

FeaturePlot(tum, features = c("nCount_RNA", "nFeature_RNA", 'percent.mt', "doublet_scores"), 
            min.cutoff = "q05",max.cutoff = "q95", order = T,raster = T,pt.size = 0.5)& 
  scale_colour_gradientn(colours = rev(brewer.pal(11,'Spectral')))
FeaturePlot(tum, features = c("percent.rps", "percent.rpl", 'G2M.Score', "S.Score"), 
            min.cutoff = "q05",max.cutoff = "q95", order = T,raster = T,pt.size = 0.5)& 
  scale_colour_gradientn(colours = rev(brewer.pal(11,'Spectral')))

print(FeaturePlot(tum, features = c("MLANA",'PMEL','MITF','AXL'), pt.size = 0.5,
                  min.cutoff = "q05",max.cutoff = "q95", order = T,raster = T)& 
        scale_colour_gradientn(colours = rev(brewer.pal(11,'Spectral'))))

print(FeaturePlot(tum, features = c("MKI67",'PTPRC','VWF','COL1A1'), pt.size = 0.5,
                  min.cutoff = "q05",max.cutoff = "q95", order = T,raster = T)& 
        scale_colour_gradientn(colours = rev(brewer.pal(11,'Spectral'))))

DimPlot(tum, reduction = "umap",label = T,group.by = 'celltype_bped_fine',repel = T,
        label.size = 2.5,pt.size = 0.01,raster = T) + 
  guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
  theme(legend.text=element_text(size=6))
dev.off()


#### Clean tumor cells v2 (remove TCR doublets) ####

tum<-readRDS('data/merged/data_tum_merged.rds')

tum <- subset(tum, tcr=='no_TCR')

FeatureScatter(tum,feature1 = 'PTPRC',feature2 = 'doublet_scores',raster = T,
               group.by = 'malignant',shuffle = T)

# Seurat workflow
tum <- NormalizeData(tum) %>% ScaleData()
tum<-FindVariableFeatures(tum)
tum <- RunPCA(tum, npcs = 100)
ElbowPlot(tum,ndims = 100)
tum <- RunUMAP(tum, dims = 1:50)
tum <- FindNeighbors(tum, dims = 1:50)
tum <- FindClusters(tum,resolution = 0.5)

DimPlot(tum, reduction = "umap",label = T,group.by = 'sample',raster = T,
        shuffle = T,pt.size = 0.01)+ NoLegend()


pdf('data/merged/plots_tum_merged_clusters_v2.pdf')
DimPlot(tum, reduction = "umap",label = T,group.by = 'ident',raster = T,
        shuffle = T,pt.size = 0.01)
DimPlot(tum, reduction = "umap",label = F,group.by = 'sample',raster = T,
        shuffle = T,pt.size = 0.01)+ 
  guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
  theme(legend.text=element_text(size=6))

DimPlot(tum, reduction = "umap",label = T,group.by = 'patient',raster = T,
        shuffle = T,pt.size = 0.01)

ggplot(tum@meta.data, aes(x=RNA_snn_res.0.3, fill=sample))+
  geom_bar(position="fill",size=0.3,color='black') +
  theme_classic() + 
  ggtitle('Sample distribution among clusters')+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggplot(tum@meta.data, aes(x=RNA_snn_res.0.3, fill=patient))+
  geom_bar(position="fill",size=0.3,color='black') +
  theme_classic() + 
  ggtitle('Patient distribution among clusters')+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggplot(tum@meta.data, aes(x=RNA_snn_res.0.3, fill=malignant))+
  geom_bar(position="fill",size=0.3,color='black') +
  theme_classic() + 
  ggtitle('malignant distribution among clusters')+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggplot(tum@meta.data, aes(x=sample))+
  geom_bar(stat="count",size=0.3,color='black') +
  theme_classic() + 
  ggtitle('#cells per sample')+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

print(DimPlot(tum, reduction = "umap",label = F,group.by = 'orig.ident',raster = T,
              shuffle = T,pt.size = 0.01)+ 
        guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
        theme(legend.text=element_text(size=6)))
print(DimPlot(tum, reduction = "umap",label = F,group.by = 'tissue',raster = T,
              shuffle = T,pt.size = 0.01))
print(DimPlot(tum, reduction = "umap",label = T,group.by = 'sample_group',raster = T,
              shuffle = T,pt.size = 0.01)+ 
        guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
        theme(legend.text=element_text(size=6)))
print(DimPlot(tum, reduction = "umap",label = T,group.by = 'cell_cycle',
              pt.size = 0.01,raster = T,shuffle = T))
DimPlot(tum, label = T,group.by = 'malignant',raster = T,shuffle = T,
        label.size = 3,repel = T,pt.size=0.1)
FeaturePlot(tum, features = 'proportion_scaled_cnv_avg',min.cutoff = "q05", 
            max.cutoff = "q95",order=T,raster=T,pt.size = 0.1)& 
  scale_colour_gradientn(colours = rev(brewer.pal(11,'Spectral')))

FeaturePlot(tum, features = c("nCount_RNA", "nFeature_RNA", 'percent.mt', "doublet_scores"), 
            min.cutoff = "q05",max.cutoff = "q95", order = T,raster = T,pt.size = 0.5)& 
  scale_colour_gradientn(colours = rev(brewer.pal(11,'Spectral')))
FeaturePlot(tum, features = c("percent.rps", "percent.rpl", 'G2M.Score', "S.Score"), 
            min.cutoff = "q05",max.cutoff = "q95", order = T,raster = T,pt.size = 0.5)& 
  scale_colour_gradientn(colours = rev(brewer.pal(11,'Spectral')))

print(FeaturePlot(tum, features = c("MLANA",'PMEL','MITF','AXL'), pt.size = 0.5,
                  min.cutoff = "q05",max.cutoff = "q95", order = T,raster = T)& 
        scale_colour_gradientn(colours = rev(brewer.pal(11,'Spectral'))))

print(FeaturePlot(tum, features = c("MKI67",'PTPRC','VWF','COL1A1'), pt.size = 0.5,
                  min.cutoff = "q05",max.cutoff = "q95", order = T,raster = T)& 
        scale_colour_gradientn(colours = rev(brewer.pal(11,'Spectral'))))

FeatureScatter(tum,feature1 = 'PTPRC',feature2 = 'doublet_scores',raster = T,
               group.by = 'malignant',shuffle = T)

print(VlnPlot(tum,features = c("MLANA",'PTPRC','proportion_scaled_cnv_avg','doublet_scores'),
              pt.size = 0,ncol = 2) &
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,size=5)))

DimPlot(tum, reduction = "umap",label = T,group.by = 'celltype_bped_fine',repel = T,
        label.size = 2.5,pt.size = 0.01,raster = T) + 
  guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
  theme(legend.text=element_text(size=6))
dev.off()

saveRDS(tum,'data/merged/data_tum_merged_v2.rds')


#### Clean tumor cells v3 (remove samples with <50 cells) & v4 (remove 36 and non-malignant) ####

tum<-readRDS('data/merged/data_tum_merged_v2.rds')
tum <- subset(tum, orig.ident %notin% c('F25_post','F10_post'))

# Seurat workflow
tum <- NormalizeData(tum) %>% ScaleData()
tum<-FindVariableFeatures(tum)
tum <- RunPCA(tum, npcs = 100)
ElbowPlot(tum,ndims = 100)
tum <- RunUMAP(tum, dims = 1:50)
tum <- FindNeighbors(tum, dims = 1:50)
tum <- FindClusters(tum,resolution = 0.3)

saveRDS(tum,'data/merged/data_tum_merged_v3.rds')


tum<-readRDS('data/merged/data_tum_merged_v3.rds')
tum <- subset(tum, RNA_snn_res.0.3 != 36 & malignant=='malignant')

# Seurat workflow
tum <- NormalizeData(tum) %>% ScaleData()
tum<-FindVariableFeatures(tum)
tum <- RunPCA(tum, npcs = 100)
ElbowPlot(tum,ndims = 100)
tum <- RunUMAP(tum, dims = 1:50)
tum <- FindNeighbors(tum, dims = 1:50)
tum <- FindClusters(tum,resolution = 0.3)

tum <- CellCycleScoring(tum, s.features = cc.genes.updated.2019$s.genes,
                           g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = F)
tum$cell_cycle <- ifelse(tum$G2M.Score > 0.15 | tum$S.Score > 0.2, 'cycling',
                            'non-cycling')

# cell type assignment after manual annotation based on DGE
tum$cell_type_main <- 'Tumor cells'
tum$cell_type_fine <- 'Tumor cells'

celltype2 <- tum@meta.data %>% dplyr::select('barcode_pat', 'cell_type_main')
write.csv(celltype2, paste0(directory_tum,'data_',label,'_celltype_main_tum.csv'))
celltype3 <- tum@meta.data %>% dplyr::select('barcode_pat', 'cell_type_fine','cell_cycle')
write.csv(celltype3, paste0(directory_tum,'data_',label,'_celltype_fine_tum.csv'))

saveRDS(tum,'data/merged/data_tum_merged_v4.rds')

pdf('data/merged/plots_tum_merged_clusters_v4.pdf')
DimPlot(tum, reduction = "umap",label = T,group.by = 'ident',raster = T,
        shuffle = T,pt.size = 0.01)
DimPlot(tum, reduction = "umap",label = F,group.by = 'sample',raster = T,
        shuffle = T,pt.size = 0.01)+ 
  guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
  theme(legend.text=element_text(size=6))

DimPlot(tum, reduction = "umap",label = T,group.by = 'patient',raster = T,
        shuffle = T,pt.size = 0.01)

ggplot(tum@meta.data, aes(x=RNA_snn_res.0.3, fill=patient))+
  geom_bar(position="fill",size=0.3,color='black') +
  theme_classic() + 
  ggtitle('Patient distribution among clusters')+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggplot(tum@meta.data, aes(x=RNA_snn_res.0.3, fill=malignant))+
  geom_bar(position="fill",size=0.3,color='black') +
  theme_classic() + 
  ggtitle('malignant distribution among clusters')+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

print(DimPlot(tum, reduction = "umap",label = F,group.by = 'orig.ident',raster = T,
              shuffle = T,pt.size = 0.01)+ 
        guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
        theme(legend.text=element_text(size=6)))
print(DimPlot(tum, reduction = "umap",label = F,group.by = 'tissue',raster = T,
              shuffle = T,pt.size = 0.01))
print(DimPlot(tum, reduction = "umap",label = T,group.by = 'sample_group',raster = T,
              shuffle = T,pt.size = 0.01)+ 
        guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
        theme(legend.text=element_text(size=6)))
FeatureScatter(tum,'G2M.Score','S.Score',group.by = 'cell_cycle',raster = T)
print(DimPlot(tum, reduction = "umap",label = T,group.by = 'cell_cycle',
              raster = T,shuffle = T))
DimPlot(tum, label = T,group.by = 'malignant',raster = T,shuffle = T,
        label.size = 3,repel = T,pt.size=0.1)
FeaturePlot(tum, features = 'proportion_scaled_cnv_avg',min.cutoff = "q05", 
            max.cutoff = "q95",order=T,raster=T,pt.size = 0.1)& 
  scale_colour_gradientn(colours = rev(brewer.pal(11,'Spectral')))

FeaturePlot(tum, features = c("nCount_RNA", "nFeature_RNA", 'percent.mt', "doublet_scores"), 
            min.cutoff = "q05",max.cutoff = "q95", order = T,raster = T,pt.size = 0.5)& 
  scale_colour_gradientn(colours = rev(brewer.pal(11,'Spectral')))
FeaturePlot(tum, features = c("percent.rps", "percent.rpl", 'G2M.Score', "S.Score"), 
            min.cutoff = "q05",max.cutoff = "q95", order = T,raster = T,pt.size = 0.5)& 
  scale_colour_gradientn(colours = rev(brewer.pal(11,'Spectral')))

print(FeaturePlot(tum, features = c("MLANA",'PMEL','MITF','AXL'), pt.size = 0.5,
                  min.cutoff = "q05",max.cutoff = "q95", order = T,raster = T)& 
        scale_colour_gradientn(colours = rev(brewer.pal(11,'Spectral'))))

print(FeaturePlot(tum, features = c("MKI67",'PTPRC','VWF','COL1A1'), pt.size = 0.5,
                  min.cutoff = "q05",max.cutoff = "q95", order = T,raster = T)& 
        scale_colour_gradientn(colours = rev(brewer.pal(11,'Spectral'))))

FeatureScatter(tum,feature1 = 'PTPRC',feature2 = 'doublet_scores',raster = T,
               group.by = 'malignant',shuffle = T)

VlnPlot(tum,features = c("MLANA",'PTPRC','proportion_scaled_cnv_avg','doublet_scores'),
        pt.size = 0,ncol = 2) &
  theme(axis.text.x = element_text(size=5,angle = 90, vjust = 0.5, hjust = 1))

DimPlot(tum, reduction = "umap",label = T,group.by = 'celltype_bped_fine',repel = T,
        label.size = 2.5,pt.size = 0.01,raster = T) + 
  guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
  theme(legend.text=element_text(size=6))
dev.off()

