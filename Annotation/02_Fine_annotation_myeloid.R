#!/usr/bin/env Rscript

#### Fine cell type annotation of myeloid cell subset
#### Author: Jana Biermann, PhD

library(dplyr)
library(Seurat)
library(ggplot2)
library(gplots)
library(viridis)
library(hypeR)
library(RColorBrewer)
library(ggpubr)
library(scales)
'%notin%' <- Negate('%in%')

setwd('~/Documents/melanoma')
label<-'melanoma'
celltype<-'myeloid'
filename<-paste0(label,'_',celltype)
path.ct <- paste0('data/annotation/',celltype,'/')
ifelse(!dir.exists(file.path(path.ct)), 
       dir.create(file.path(path.ct)), FALSE)


seu<-readRDS('data/merged/data_melanoma_merged_v2.rds')

# Subset to myeloid cells
seu <- subset(seu, cell_type_main == 'Myeloid cells')
seu<-NormalizeData(seu)
seu <- ScaleData(object = seu)
seu<-FindVariableFeatures(seu)
seu <- RunPCA(object = seu,npcs = 75)
ElbowPlot(seu,ndims = 75)
seu <- FindNeighbors(seu, dims = 1:30)
seu <- FindClusters(seu, resolution = 0.4)
seu <- RunUMAP(object = seu, dims = 1:30)


## Integrate
sort(table(seu$sample))
DefaultAssay(seu)<-'RNA'
obj.list <- SplitObject(seu, split.by = "sample")
# Keep only patient with >50 cells
for (i in 1:length(obj.list)) {
  if (dim(obj.list[[i]]@assays$RNA@counts)[2] < 50) {
    obj.list[i] <- NA
  } else {
    obj.list[[i]] <- NormalizeData(obj.list[[i]])
    obj.list[[i]] <- FindVariableFeatures(obj.list[[i]])
  }
}
# Remove NAs
obj.list <- obj.list[!is.na(obj.list)]

anchors <- FindIntegrationAnchors(object.list = obj.list)
seu <- IntegrateData(anchorset = anchors,k.weight = 50)
#rm(anchors,obj.list)

#workflow
seu <- ScaleData(object = seu)
seu <- RunPCA(object = seu,npcs = 75)
ElbowPlot(seu,ndims = 75)
seu <- FindNeighbors(seu, dims = 1:20)
seu <- FindClusters(seu, resolution = 0.4)
seu <- RunUMAP(object = seu, dims = 1:20)

FeaturePlot(seu, features = c('CD3D','CD3E','TOX','MLANA','CD163','LYZ'),
            min.cutoff = "q05", max.cutoff = "q95",order=T,raster=T)& 
  scale_colour_gradientn(colours = rev(brewer.pal(11,'Spectral')))

# Cell cycle assignment
DefaultAssay(seu) <- 'RNA'
seu <- CellCycleScoring(seu, s.features = cc.genes.updated.2019$s.genes, 
                        g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = F)
DefaultAssay(seu) <- 'integrated'
seu$cell_cycle <- ifelse(seu$G2M.Score > 0.1 | seu$S.Score > 0.1, 'cycling', 'non-cycling')

seu$cell_type_main<-ifelse(seu$integrated_snn_res.0.3 %in% c(5,6,7,13),'Doublets','Myeloid cells')

# Save
ct <- seu@meta.data %>% dplyr::select('barcode_pat', 'cell_type_main')
write.csv(ct, paste0(path.ct,'data_',celltype,'_pt1.csv'))
ifelse(!dir.exists(file.path(path.ct)), 
       dir.create(file.path(path.ct),recursive = T), FALSE)
saveRDS(seu,paste0(path.ct,'data_',filename,'_pt1.rds'))
seu<-readRDS(paste0(path.ct,'data_',filename,'_pt1.rds'))


# get only myeloid
seu <- subset(seu, cell_type_main == 'Myeloid cells')

# reintegrate
sort(table(seu$sample))
DefaultAssay(seu)<-'RNA'
obj.list <- SplitObject(seu, split.by = "sample")
# Keep only patient with >50 cells
for (i in 1:length(obj.list)) {
  if (dim(obj.list[[i]]@assays$RNA@counts)[2] < 50) {
    obj.list[i] <- NA
  } else {
    obj.list[[i]] <- NormalizeData(obj.list[[i]])
    obj.list[[i]] <- FindVariableFeatures(obj.list[[i]])
  }
}
# Remove NAs
obj.list <- obj.list[!is.na(obj.list)]

anchors <- FindIntegrationAnchors(object.list = obj.list)
seu <- IntegrateData(anchorset = anchors,k.weight = 50)

seu <- ScaleData(object = seu)
seu <- RunPCA(object = seu,npcs = 50)
ElbowPlot(seu,ndims = 50)
seu <- FindNeighbors(seu, dims = 1:20)
seu <- FindClusters(seu, resolution = 0.4)
seu <- RunUMAP(object = seu, dims = 1:20)

# get module scores for sigs
sigs1<-read.csv('~/brain_mets/signatures/myeloid_signatures.csv',na.strings = '')
for(c in 1:ncol(sigs1)){
  seu<-AddModuleScore(object = seu,features = list(na.omit(sigs1[,c])),
                      name = colnames(sigs1)[c],assay = 'RNA',search=F)
}
sigs2<-read.csv('~/brain_mets/signatures/macrophages_ICI.csv',na.strings = '')
for(c in 1:ncol(sigs2)){
  seu<-AddModuleScore(object = seu,features = list(na.omit(sigs2[,c])),
                      name = colnames(sigs2)[c],assay = 'RNA',search=F)
}


# Save obj
ifelse(!dir.exists(file.path(path.ct)), 
       dir.create(file.path(path.ct),recursive = T), FALSE)
saveRDS(seu,paste0(path.ct,'data_',filename,'_pt2_reint.rds'))

seu<-readRDS(paste0(path.ct,'data_',filename,'_pt2_reint.rds'))


# Differential gene expression (DGE) based on clusters
markers <- FindAllMarkers(seu, only.pos = TRUE,assay='RNA', min.pct = 0.25, 
                          logfc.threshold = 0.25,test.use = 'MAST')
write.csv(markers, paste0(path.ct,'markers_',filename,'.csv'), row.names = F)
markers %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% top_n(15, wt=avg_log2FC) -> tops

DefaultAssay(seu)<-'RNA'
seu<-NormalizeData(seu) %>% ScaleData()
DefaultAssay(seu)<-'integrated'

pdf(file = paste0(path.ct,'markers_',filename,'_heatmap.pdf'),width = 10, height = 18)
DoHeatmap(seu, features = tops$gene, group.by = 'ident',raster = T,assay = 'RNA')
dev.off()


# Fine cell type assignment after manual annotation based on DGE reint
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(5), 'Macrophages (M1-like)', 'NA')
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(4), 'Macrophages (M2-like)', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(3), 'Monocytes',seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(10), 'cDC1', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(8), 'cDC2', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(9), 'DC3', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(7), 'Macrophages (cycling)',seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(0,2), 'Macrophages (phagocytic)', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(1,6,11), 'Macrophages (antigen-pres.)', seu$cell_type_fine)




### stats
stats<-as.data.frame(matrix(data=NA,nrow = 1, ncol = 11))
colnames(stats)<-c('sample','n_features','n_cells','median_features','median_counts')
rownames(stats)<-filename
stats$sample<-filename
stats$n_features<-dim(seu@assays$integrated@data)[1]
stats$n_cells<-dim(seu@assays$integrated@data)[2]
stats$median_features<-round(median(seu@meta.data$nFeature_RNA))
stats$median_counts<-round(median(seu@meta.data$nCount_RNA))

pdf(paste0(path.ct,'plots_',filename,'_reint2.pdf'))
textplot(t(stats),cex=1.2,halign='left')

print(DimPlot(seu, reduction = "umap",label = T,group.by = 'ident',shuffle = T,raster = T))
print(DimPlot(seu, reduction = "umap",label = F,group.by = 'patient',shuffle = T,raster = T))
print(DimPlot(seu, reduction = "umap",label = F,group.by = 'tissue',raster = T,
              shuffle = T))
print(DimPlot(seu, reduction = "umap",label = T,group.by = 'treated',raster = T,
              shuffle = T)+ 
        guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
        theme(legend.text=element_text(size=6)))
print(DimPlot(seu, reduction = "umap",label = T,group.by = 'cell_type_fine',
              shuffle = T,raster = T,repel = T))
print(DimPlot(seu, reduction = "umap",label = T,group.by = 'cell_type_main',
              shuffle = T,raster = T,repel = T))

print(DimPlot(seu, reduction = "umap",label = F,group.by = 'malignant',shuffle = T,raster = T))
print(DimPlot(seu, reduction = "umap",label = F,group.by = 'tcr',shuffle = T,raster = T))


FeatureScatter(seu,'G2M.Score','S.Score',group.by = 'cell_cycle',shuffle = T,raster = T)
print(DimPlot(seu, reduction = "umap",label = F,group.by = 'cell_cycle',shuffle = T,raster = T))

DimPlot(seu, reduction = "umap",label = T,group.by = 'celltype_bped_fine',
        repel = T,label.size = 3,shuffle = T,raster = T) + 
  guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
  theme(legend.text=element_text(size=7))

DimPlot(seu, reduction = "umap",label = T,group.by = 'celltype_bped_main',
        repel = T,label.size = 3,shuffle = T,raster = T) + 
  guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
  theme(legend.text=element_text(size=7))

DimPlot(seu, reduction = "umap",label = T,group.by = 'celltype_hpca_main',
        repel = T,label.size = 2.5,shuffle = T,raster = T) + 
  guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
  theme(legend.text=element_text(size=7))

print(FeaturePlot(seu, features = c("nCount_RNA", "nFeature_RNA", 'percent.mt'), 
                  min.cutoff = "q05",max.cutoff = "q95", order = T,raster = T)& 
        scale_colour_gradientn(colours = rev(brewer.pal(11,'Spectral'))))
print(FeaturePlot(seu, features = c("percent.rps", "percent.rpl", 'G2M.Score', "S.Score",'frequency','doublet_scores'), 
                   order = T,raster = F)& 
        scale_colour_gradientn(colours = rev(brewer.pal(11,'Spectral'))))

print(FeaturePlot(seu, features = c("FCGR3A", "CD14"), min.cutoff = "q05",max.cutoff = "q95",
                  order = T,raster = F)& 
        scale_colour_gradientn(colours = rev(brewer.pal(11,'Spectral'))))

VlnPlot(seu,c("M1_macrophage_polarization1","M1_bi1","M11",
              "M2_macrophage_polarization1","M2_bi1","M21"),
        pt.size = 0,group.by = 'integrated_snn_res.0.4')+NoLegend()

VlnPlot(seu,c("DC1_CD141.CLEC9A1","cDC11","DC2_CD1C_A1","cDC21",'DC3_CD1C_B1','DC41'),pt.size = 0,
        group.by = 'integrated_snn_res.0.4',flip = T,ncol = 3)+NoLegend()

VlnPlot(seu,c('DC51','DC6_pDC1',"pDC1","DC3_zilionis1","pDC_zilionis1","monoDC_zilionis1"),
        pt.size = 0,group.by = 'integrated_snn_res.0.4',flip = T,ncol=3)+NoLegend()


DotPlot(seu,assay = 'RNA',features = c(paste0(names(sigs1),1),paste0(names(sigs2),1)),
        group.by = 'integrated_snn_res.0.4', dot.scale = 3)+ 
  scale_color_viridis() + coord_flip()

DotPlot(seu,assay = 'RNA',features = c('FCN1','VCAN','LYZ', #monocyte
                                       'IL2RA','CD86','KLF6','HTRA1','PDK4','PLD4','IRF1','JAK2','ICAM1', #M1
                                       'PPARG','MERTK','FAM20C','F13A1','STAB1','SELENOP','CD163','FTL',
                                       'NR4A2','MRC1','CD163L1','CD5L','TMSB10','GPNMB','HS3ST2','IGF1','MAFB', 
                                       'TBXAS1',#M2
                                       'CLEC7A','CD81','IL10','CLEC10A','SPP1','VEGFA','CCL2' #TAM
),group.by = 'integrated_snn_res.0.4', dot.scale = 5)+ scale_color_viridis() + coord_flip() +RotatedAxis()

DotPlot(seu,assay = 'RNA',features = c('FCN1','VCAN','LYZ','CD14','FCGR3A','FCGR3B', #monocyte
                                       'IL2RA','CD86',#M1
                                       'PPARG','MERTK','FAM20C','F13A1','STAB1','SELENOP','FTL', #M2
                                       'CLEC9A','XCR1','BATF3','THBD',#cDC1
                                       'CD1C','FCGR2B', #cDC2
                                       'IL3RA', 'CLEC4C','TCF4','LILRA4','CD83','GZMB','JCHAIN','ITM2C', #pDC
                                       'IRF4','ITGAX','PF4',#mo-DC
                                       'CD5L','VCAM1','MARCO','VSIG4','CPVL','C1QA','CD163',
                                       'C3',#mg
                                       'rna_ASCL1'
),group.by = 'integrated_snn_res.0.4', dot.scale = 5)+ 
  scale_color_viridis() + coord_flip() +RotatedAxis()

DotPlot(seu,assay = 'RNA',features = c('CLEC9A','XCR1', #cDC1
                                       'CD1C','FCGR2B', #cDC2
                                       'BIRC3','LAMP3', #DC3
                                       'HLA-A','HLA-B','B2M', #antigen-pres
                                       'MKI67','TOP2A', #cycling
                                       'CD86','ICAM1','JAK2','IRF1', #M1
                                       'F13A1','MERTK','SELENOP', #M2
                                       'MSR1','CTSL','FCHO2', #phago
                                       'FCN1','VCAN' #monocyte
),group.by = 'cell_type_fine', dot.scale = 5)+ scale_color_viridis() + coord_flip() +RotatedAxis()

VlnPlot(seu,c('nCount_RNA', 'nFeature_RNA', 'percent.mt', 'doublet_scores','percent.rps', 'percent.rpl'),
        pt.size = 0,group.by = 'integrated_snn_res.0.4')+NoLegend()

for(i in seq(1, ncol(sigs1), 4)){
  four<-paste0(colnames(sigs1)[i:(i+3)],'1')
  print(FeaturePlot(seu, features = c(four),
                    min.cutoff = "q05", max.cutoff = "q95",order=T,raster=T)& 
          scale_colour_gradientn(colours = rev(brewer.pal(11,'Spectral'))))
}

for(i in seq(1, ncol(sigs2), 4)){
  four<-paste0(colnames(sigs2)[i:(i+3)],'1')
  print(FeaturePlot(seu, features = c(four),
                    min.cutoff = "q05", max.cutoff = "q95",order=T,raster=T)& 
          scale_colour_gradientn(colours = rev(brewer.pal(11,'Spectral'))))
}
VlnPlot(seu,c("jerby_mac_in_cold_tumors_up1","jerby_mac_in_cold_tumors_down1","jerby_mac_post_ICI_up1",
              "jerby_mac_post_ICI_down1"),
        pt.size = 0,group.by = 'treated',ncol = 2,split.by = 'responder')+
  theme(legend.position = 'bottom')
VlnPlot(seu,c("jerby_mac_in_cold_tumors_up1","jerby_mac_in_cold_tumors_down1","jerby_mac_post_ICI_up1",
              "jerby_mac_post_ICI_down1"),
        pt.size = 0,group.by = 'cell_type_fine',ncol = 2,split.by = 'responder')+
  theme(legend.position = 'right')
VlnPlot(seu,c("jerby_mac_in_cold_tumors_up1","jerby_mac_in_cold_tumors_down1","jerby_mac_post_ICI_up1",
              "jerby_mac_post_ICI_down1"),
        pt.size = 0,group.by = 'responder',ncol = 2)
VlnPlot(seu,c("jerby_mac_in_cold_tumors_up1","jerby_mac_in_cold_tumors_down1","jerby_mac_post_ICI_up1",
              "jerby_mac_post_ICI_down1"),
        pt.size = 0,group.by = 'integrated_snn_res.0.4',ncol = 2)
dev.off()


# Save
ct <- seu@meta.data %>% dplyr::select('barcode_pat', 'cell_type_fine','cell_cycle')
write.csv(ct, paste0(path.ct,'data_',filename,'.csv'))


# Differential gene expression (DGE) based on cell types
Idents(seu) <- seu$cell_type_fine
markers <- FindAllMarkers(seu, only.pos = TRUE, assay = 'RNA', min.pct = 0.25,
                          logfc.threshold = 0.25, test.use = 'MAST')
write.csv(markers, paste0(path.ct, 'markers_', filename, '_celltype.csv'), row.names = F)
markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  top_n(15, wt = avg_log2FC) -> tops

DefaultAssay(seu) <- 'RNA'
seu <- NormalizeData(seu) %>% ScaleData(features = rownames(seu@assays$RNA@data))
DefaultAssay(seu) <- 'integrated'

pdf(file = paste0(path.ct, 'markers_', filename, '_heatmap_celltype.pdf'), width = 10,
    height = 18)
DoHeatmap(seu, features = tops$gene, group.by = 'cell_type_fine', raster = T,
          assay = 'RNA')
dev.off()

#### quantification ####
df_treated<-seu@meta.data %>%
  group_by(sample, cell_type_fine, treated) %>%
  tally() %>%
  group_by(sample) %>%
  mutate(freq = n/sum(n))

df_mucosal<-seu@meta.data %>%
  group_by(sample, cell_type_fine, mucosal) %>%
  tally() %>%
  group_by(sample) %>%
  mutate(freq = n/sum(n))

df_responder<-seu@meta.data %>%
  group_by(sample, cell_type_fine, responder) %>%
  tally() %>%
  group_by(sample) %>%
  mutate(freq = n/sum(n))

df_time<-seu@meta.data %>%
  group_by(sample, cell_type_fine, time) %>%
  tally() %>%
  group_by(sample) %>%
  mutate(freq = n/sum(n))

df_tissue<-seu@meta.data %>%
  group_by(sample, cell_type_fine, tissue) %>%
  tally() %>%
  group_by(sample) %>%
  mutate(freq = n/sum(n))


pdf(paste0(path.ct,'plots_',filename,'_reint2_quantification.pdf'))
## boxplots
ggplot(df_treated,aes(x=cell_type_fine, freq,col=treated)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position=position_jitterdodge(jitter.width = 0.2),size=1,alpha=0.7)+
  stat_compare_means(method = "wilcox.test", label='p.format',size=3) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))

ggplot(df_mucosal,aes(x=cell_type_fine, freq,col=mucosal)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position=position_jitterdodge(jitter.width = 0.2),size=1,alpha=0.7)+
  stat_compare_means(method = "wilcox.test", label='p.format',size=3) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))

ggplot(df_responder,aes(x=cell_type_fine, freq,col=responder)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position=position_jitterdodge(jitter.width = 0.2),size=1,alpha=0.7)+
  stat_compare_means(method = "wilcox.test", label='p.format',size=3) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))

ggplot(df_time,aes(x=cell_type_fine, freq,col=time)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position=position_jitterdodge(jitter.width = 0.2),size=1,alpha=0.7)+
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(c("pre", "on"), c("pre", "post"), c("on", "post")))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))

ggplot(df_tissue,aes(x=cell_type_fine, freq,col=tissue)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position=position_jitterdodge(jitter.width = 0.2),size=1,alpha=0.7)+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))

## barplots
ggplot(seu@meta.data, aes(x = sample, fill = cell_type_fine)) + 
  geom_bar(stat = 'count', size = 0.3, color = 'black') + 
  theme_classic() + 
  theme(legend.position = 'right', 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggplot(seu@meta.data, aes(x = sample, fill = cell_type_fine)) + 
  geom_bar(position = 'fill', size = 0.3, color = 'black') + 
  theme_classic() + 
  theme(legend.position = 'right', 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggplot(seu@meta.data, aes(x = treated, fill = cell_type_fine)) + 
  geom_bar(stat = 'count', size = 0.3, color = 'black') + 
  theme_classic() + 
  theme(legend.position = 'right', 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.margin = margin(0.5, 8, 0.5, 0.5, "cm"))

ggplot(seu@meta.data, aes(x = treated, fill = cell_type_fine)) + 
  geom_bar(position = 'fill', size = 0.3, color = 'black') + 
  theme_classic() + 
  theme(legend.position = 'right', 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.margin = margin(0.5, 8, 0.5, 0.5, "cm"))

ggplot(seu@meta.data, aes(x = mucosal, fill = cell_type_fine)) + 
  geom_bar(position = 'fill', size = 0.3, color = 'black') + 
  theme_classic() + 
  theme(legend.position = 'right', 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.margin = margin(0.5, 8, 0.5, 0.5, "cm"))

ggplot(seu@meta.data, aes(x = responder, fill = cell_type_fine)) + 
  geom_bar(position = 'fill', size = 0.3, color = 'black') + 
  theme_classic() + 
  theme(legend.position = 'right', 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.margin = margin(0.5, 8, 0.5, 0.5, "cm"))

ggplot(seu@meta.data, aes(x = time, fill = cell_type_fine)) + 
  geom_bar(position = 'fill', size = 0.3, color = 'black') + 
  theme_classic() + 
  theme(legend.position = 'right', 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.margin = margin(0.5, 8, 0.5, 0.5, "cm"))

dev.off()


#### myeloid sigs vs t cell ratio ####
meta<-readRDS('data/annotation/main/data_melanoma_metadata.rds')

df_meta<-meta %>%
  group_by(sample, cell_type_main) %>%
  tally() %>%
  group_by(sample) %>%
  mutate(freq = n/sum(n))

df_meta_t<- df_meta %>%
  filter(cell_type_main=='T/NK cells') %>%
  select(sample, freq)

colnames(df_meta_t)[2]<-'t_cell_fraction'

seu$barcode_all<-rownames(seu@meta.data)
seu@meta.data<-left_join(seu@meta.data,df_meta_t, by='sample')
rownames(seu@meta.data)<-seu$barcode_all

seu<-subset(seu,cell_type_fine %in% c('Macrophages (M1-like)','Macrophages (M2-like)',
                                      'Macrophages (phagocytic)','Macrophages (antigen-pres.)'))

pdf(paste0(path.ct,'plots_',filename,'_correlation_Tcell_fraction.pdf'))
# jerby_mac_in_cold_tumors_up1
df<-seu@meta.data %>%
  select(sample,t_cell_fraction,jerby_mac_in_cold_tumors_up1,responder,treated) %>%
  group_by(sample) %>%
  summarise(t_cell_fraction=mean(t_cell_fraction),
            jerby_mac_in_cold_tumors_up1=mean(jerby_mac_in_cold_tumors_up1),
            responder=unique(responder),
            treated=unique(treated))

corr_res<-cor.test(df$t_cell_fraction[df$responder=='Responder'], 
                   df$jerby_mac_in_cold_tumors_up1[df$responder=='Responder'],
                   method = 'spearman')
corr_nonres<-cor.test(df$t_cell_fraction[df$responder=='Non-responder'], 
                   df$jerby_mac_in_cold_tumors_up1[df$responder=='Non-responder'],
                   method = 'spearman')

ggplot(df, aes(x=t_cell_fraction , y=jerby_mac_in_cold_tumors_up1,col=treated)) +
  geom_point()+
  geom_smooth(method=lm,  col='darkgray',lty=2)+
  theme_bw()+
  ggtitle(paste0('t_cell_fraction vs jerby_mac_in_cold_tumors_up1\nResonders: Rho = ',
                 round(corr_res$estimate[[1]],3),
                 ', P = ', scientific(corr_res$p.value,digits = 3),
                 '\nNon-resonders: Rho = ',
                 round(corr_nonres$estimate[[1]],3),
                 ', P = ', scientific(corr_nonres$p.value,digits = 3)))+
  theme(legend.position = 'right')+
  facet_wrap(~responder)

# jerby_mac_in_cold_tumors_down1
df<-seu@meta.data %>%
  select(sample,t_cell_fraction,jerby_mac_in_cold_tumors_down1,responder,treated) %>%
  group_by(sample) %>%
  summarise(t_cell_fraction=mean(t_cell_fraction),
            jerby_mac_in_cold_tumors_down1=mean(jerby_mac_in_cold_tumors_down1),
            responder=unique(responder),
            treated=unique(treated))

corr_res<-cor.test(df$t_cell_fraction[df$responder=='Responder'], 
                   df$jerby_mac_in_cold_tumors_down1[df$responder=='Responder'],
                   method = 'spearman')
corr_nonres<-cor.test(df$t_cell_fraction[df$responder=='Non-responder'], 
                      df$jerby_mac_in_cold_tumors_down1[df$responder=='Non-responder'],
                      method = 'spearman')

ggplot(df, aes(x=t_cell_fraction , y=jerby_mac_in_cold_tumors_down1,col=treated)) +
  geom_point()+
  geom_smooth(method=lm,  col='darkgray',lty=2)+
  theme_bw()+
  ggtitle(paste0('t_cell_fraction vs jerby_mac_in_cold_tumors_down1\nResonders: Rho = ',
                 round(corr_res$estimate[[1]],3),
                 ', P = ', scientific(corr_res$p.value,digits = 3),
                 '\nNon-resonders: Rho = ',
                 round(corr_nonres$estimate[[1]],3),
                 ', P = ', scientific(corr_nonres$p.value,digits = 3)))+
  theme(legend.position = 'right')+
  facet_wrap(~responder)

# jerby_mac_post_ICI_up1
df<-seu@meta.data %>%
  select(sample,t_cell_fraction,jerby_mac_post_ICI_up1,responder,treated) %>%
  group_by(sample) %>%
  summarise(t_cell_fraction=mean(t_cell_fraction),
            jerby_mac_post_ICI_up1=mean(jerby_mac_post_ICI_up1),
            responder=unique(responder),
            treated=unique(treated))

corr_res<-cor.test(df$t_cell_fraction[df$responder=='Responder'], 
                   df$jerby_mac_post_ICI_up1[df$responder=='Responder'],
                   method = 'spearman')
corr_nonres<-cor.test(df$t_cell_fraction[df$responder=='Non-responder'], 
                      df$jerby_mac_post_ICI_up1[df$responder=='Non-responder'],
                      method = 'spearman')

ggplot(df, aes(x=t_cell_fraction , y=jerby_mac_post_ICI_up1,col=treated)) +
  geom_point()+
  geom_smooth(method=lm,  col='darkgray',lty=2)+
  theme_bw()+
  ggtitle(paste0('t_cell_fraction vs jerby_mac_post_ICI_up1\nResonders: Rho = ',
                 round(corr_res$estimate[[1]],3),
                 ', P = ', scientific(corr_res$p.value,digits = 3),
                 '\nNon-resonders: Rho = ',
                 round(corr_nonres$estimate[[1]],3),
                 ', P = ', scientific(corr_nonres$p.value,digits = 3)))+
  theme(legend.position = 'right')+
  facet_wrap(~responder)

# jerby_mac_post_ICI_down1
df<-seu@meta.data %>%
  select(sample,t_cell_fraction,jerby_mac_post_ICI_down1,responder,treated) %>%
  group_by(sample) %>%
  summarise(t_cell_fraction=mean(t_cell_fraction),
            jerby_mac_post_ICI_down1=mean(jerby_mac_post_ICI_down1),
            responder=unique(responder),
            treated=unique(treated))

corr_res<-cor.test(df$t_cell_fraction[df$responder=='Responder'], 
                   df$jerby_mac_post_ICI_down1[df$responder=='Responder'],
                   method = 'spearman')
corr_nonres<-cor.test(df$t_cell_fraction[df$responder=='Non-responder'], 
                      df$jerby_mac_post_ICI_down1[df$responder=='Non-responder'],
                      method = 'spearman')

ggplot(df, aes(x=t_cell_fraction , y=jerby_mac_post_ICI_down1,col=treated)) +
  geom_point()+
  geom_smooth(method=lm,  col='darkgray',lty=2)+
  theme_bw()+
  ggtitle(paste0('t_cell_fraction vs jerby_mac_post_ICI_down1\nResonders: Rho = ',
                 round(corr_res$estimate[[1]],3),
                 ', P = ', scientific(corr_res$p.value,digits = 3),
                 '\nNon-resonders: Rho = ',
                 round(corr_nonres$estimate[[1]],3),
                 ', P = ', scientific(corr_nonres$p.value,digits = 3)))+
  theme(legend.position = 'right')+
  facet_wrap(~responder)
dev.off()
