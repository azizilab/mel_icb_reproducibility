#!/usr/bin/env Rscript

#### Fine cell type annotation of B cells
#### Author: Jana Biermann, PhD

library(dplyr)
library(Seurat)
library(ggplot2)
library(gplots)
library(viridis)
library(RColorBrewer)
'%notin%' <- Negate('%in%')

setwd('~/Documents/melanoma/')
label<-'melanoma'
celltype<-'bcells'
filename<-paste0(label,'_',celltype)
path.ct <- paste0('data/annotation/',celltype,'/')


seu<-readRDS('data/annotation/main/data_nontum_merged_v3.rds')
#seu$tissue<-ifelse(is.na(seu$tissue)==T,'skin',seu$tissue)

# Subset to B cells
seu <- subset(seu, cell_type_main == 'B/Plasma cells')

# get module scores for sigs
sigs<-read.csv('~/brain_mets/signatures/bcell_signatures.csv',na.strings = '')
for(c in 1:ncol(sigs)){
  seu<-AddModuleScore(object = seu,features = list(na.omit(sigs[,c])),
                      name = colnames(sigs)[c],assay = 'RNA',search=F)
}

# Cell cycle assignment
#DefaultAssay(seu) <- 'RNA'
seu <- CellCycleScoring(seu, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = F)
#DefaultAssay(seu) <- 'integrated'
seu$cell_cycle <- ifelse(seu$G2M.Score > 0.1 | seu$S.Score > 0.1, 'cycling', 'non-cycling')

sort(table(seu$sample))
obj.list <- SplitObject(seu, split.by = 'sample')

has_enough_cells <- function(seurat_obj) {
  return(dim(seurat_obj)[2] >= 20)
}

# Filter the list to keep only objects with at least 20 cells
obj.list <- Filter(has_enough_cells, obj.list)

# keep only patients with >n cells
n<-40
for (i in 1:length(obj.list)) {
  obj.list[[i]] <- NormalizeData(obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(obj.list[[i]])
  if (dim(obj.list[[i]]@assays$RNA@counts)[2] < n) {
    obj.list[i] <- NA
  }
}
# remove NAs
obj.list <- obj.list[!is.na(obj.list)]

# Find integration anchors
anchors <- FindIntegrationAnchors(object.list = obj.list, dims = 1:20, k.filter = n)

# Integrate data sets
seu <- IntegrateData(anchorset = anchors, dims = 1:20,k.weight=n)

# Rerun Seurat workflow
seu <- ScaleData(seu) %>% 
  RunPCA() %>% 
  RunUMAP(dims=1:20) %>%
  FindNeighbors(dims = 1:20) %>% 
  FindClusters(resolution = 0.8)

# Save obj
ifelse(!dir.exists(file.path(path.ct)), 
       dir.create(file.path(path.ct),recursive = T), FALSE)
saveRDS(seu,paste0(path.ct,'data_',filename,'_integrated_v1.rds'))

seu<-readRDS(paste0(path.ct,'data_',filename,'_integrated_v1.rds'))


# re-integrate
seu<-subset(seu, integrated_snn_res.0.8 %notin% c(1,9,10,12,13))
DefaultAssay(seu)<-'RNA'
obj.list <- SplitObject(seu, split.by = 'sample')
n<-40
for (i in 1:length(obj.list)) {
  obj.list[[i]] <- NormalizeData(obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(obj.list[[i]])
  if (dim(obj.list[[i]]@assays$RNA@counts)[2] < n) {
    obj.list[i] <- NA
  }
}
# remove NAs
obj.list <- obj.list[!is.na(obj.list)]

# Find integration anchors
anchors <- FindIntegrationAnchors(object.list = obj.list, dims = 1:20, k.filter = n)

# Integrate data sets
seu <- IntegrateData(anchorset = anchors, dims = 1:20,k.weight=n)

# Rerun Seurat workflow
seu <- ScaleData(seu) %>% 
  RunPCA() %>% 
  RunUMAP(dims=1:20) %>%
  FindNeighbors(dims = 1:20) %>% 
  FindClusters(resolution = 0.3)

seu <- FindClusters(seu, resolution = 0.3)


# Save obj
ifelse(!dir.exists(file.path(path.ct)), 
       dir.create(file.path(path.ct),recursive = T), FALSE)
saveRDS(seu,paste0(path.ct,'data_',filename,'_integrated_v2.rds'))

seu<-readRDS(paste0(path.ct,'data_',filename,'_integrated_v2.rds'))

# remove patient-specific artifacts
seu<-subset(seu, integrated_snn_res.0.3 %notin% c(7,8))
seu <- ScaleData(seu) %>% 
  RunPCA() %>% 
  RunUMAP(dims=1:20) %>%
  FindNeighbors(dims = 1:20) %>% 
  FindClusters(resolution = 0.3)
seu <- FindClusters(seu, resolution = 0.3)
saveRDS(seu,paste0(path.ct,'data_',filename,'_integrated_v3.rds'))

seu<-readRDS(paste0(path.ct,'data_',filename,'_integrated_v3.rds'))
seu$treated<-ifelse(seu$group %in% c('pre_naive','pre_naive_muc'),'Untreated','Treated')
seu$mucosal<-ifelse(seu$group %in% c('post_PD-1_muc','pre_naive_muc'),'Mucosal','Non-mucosal')
seu$responder<-ifelse(seu$best_response %in% c('CR','PR'),'Responder','Non-responder')


# Differential gene expression (DGE) based on clusters
markers <- FindAllMarkers(seu, only.pos = TRUE,assay='RNA', min.pct = 0.25, 
                          logfc.threshold = 0.25,test.use = 'MAST')
write.csv(markers, paste0(path.ct,'markers_',filename,'_integrated.csv'), row.names = F)
markers %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% top_n(15, wt=avg_log2FC) -> tops

DefaultAssay(seu)<-'RNA'
seu<-NormalizeData(seu) %>% ScaleData()
DefaultAssay(seu)<-'integrated'

pdf(file = paste0(path.ct,'markers_',filename,'_heatmap_integrated.pdf'),width = 10, height = 18)
DoHeatmap(seu, features = tops$gene, group.by = 'ident',raster = T,assay = 'RNA')
dev.off()


# Fine cell type assignment after manual annotation based on DGE
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.3 %in% c(4), 'Class-switched memory B cells', 'NA')
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.3 %in% c(1,8), 'Memory B cells', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.3 %in% c(0,2,3,7), 'Plasma cells', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.3 %in% c(5), 'Naive B cells', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.3 %in% c(6), 'Cycling cells', seu$cell_type_fine)



### stats
stats<-as.data.frame(matrix(data=NA,nrow = 1, ncol = 11))
colnames(stats)<-c('sample','n_features','n_cells','median_features','median_counts')
rownames(stats)<-celltype
stats$sample<-celltype
stats$n_features<-dim(seu@assays$RNA@data)[1]
stats$n_cells<-dim(seu@assays$RNA@data)[2]
stats$median_features<-round(median(seu@meta.data$nFeature_RNA))
stats$median_counts<-round(median(seu@meta.data$nCount_RNA))

pdf(paste0(path.ct,'plots_',filename,'_integrated_v03.pdf'))
textplot(t(stats),cex=1.2,halign='left')
ElbowPlot(seu,ndims = 50)
print(DimPlot(seu, reduction = "umap",label = T,group.by = 'ident',shuffle = T,raster = T))
print(DimPlot(seu, reduction = "umap",label = F,group.by = 'patient',shuffle = T,raster = T))
print(DimPlot(seu, label = T,group.by = 'cell_type_fine',shuffle = T,raster = T))
print(DimPlot(seu, reduction = "umap",label = F,group.by = 'tissue',raster = T,
              shuffle = T))
print(DimPlot(seu, reduction = "umap",label = T,group.by = 'sample_group',raster = T,
              shuffle = T)+ 
        guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
        theme(legend.text=element_text(size=6)))

FeatureScatter(seu,'G2M.Score','S.Score',group.by = 'cell_cycle', pt.size = 0.01)
print(DimPlot(seu, reduction = "umap",label = F,group.by = 'cell_cycle',shuffle = T,raster = T))

DimPlot(seu, reduction = "umap",label = T,group.by = 'celltype_bped_fine',repel = T,
        label.size = 3,shuffle = T,raster = F) + 
  guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
  theme(legend.text=element_text(size=6))

DimPlot(seu, reduction = "umap",label = T,group.by = 'celltype_hpca_main',repel = T,
        label.size = 2.5,shuffle = T,raster = T) + 
  guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
  theme(legend.text=element_text(size=6))

DimPlot(seu, reduction = "umap",label = T,group.by = 'celltype_hpca_fine',repel = T,
        label.size = 2.5,shuffle = T,raster = T) + 
  guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
  theme(legend.text=element_text(size=6))

print(FeaturePlot(seu, features = c("nCount_RNA", "nFeature_RNA", 'percent.mt', "doublet_scores"), 
                  min.cutoff = "q05",max.cutoff = "q95", order = T,raster = T)& 
        scale_colour_gradientn(colours = rev(brewer.pal(11,'Spectral'))))
print(FeaturePlot(seu, features = c("percent.rps", "percent.rpl", 'G2M.Score', "S.Score"), 
                  min.cutoff = "q05",max.cutoff = "q95", order = T,raster = T)& 
        scale_colour_gradientn(colours = rev(brewer.pal(11,'Spectral'))))

DimPlot(seu, label = T,group.by = 'clone_size',raster = T,shuffle = T,
        label.size = 3,repel = T)
DimPlot(seu, label = T,group.by = 'malignant',raster = T,shuffle = T,
        label.size = 3,repel = T)
FeaturePlot(seu, features = 'proportion_scaled_cnv_avg',min.cutoff = "q05", 
            max.cutoff = "q95",order=T,raster=T)& 
  scale_colour_gradientn(colours = rev(brewer.pal(11,'Spectral')))

DotPlot(seu,assay = 'RNA',features = c('CD27','SLAMF7','IRF4','CD38','PRDM1', #plasma cell
                                       'IGHM',# immature B cell (before entering lymph node)
                                       'IGHD', # mature naive B cell (IGM and IGD; in lymph node)
                                       'IGHG1','IGHG2','IGHG3','IGHG4', 'IGHE', # class-switched/ plasma cell
                                       'IGHA1','IGHA2',
                                       'SDC1','XBP1', 'CD207','TUBB','STMN1','TYMS', #plasmablast
                                       'CXCR5','TNFRSF13B', #Follicular B cells
                                       'POU2AF1', 'CD40', 'SUGCT', #Germinal center B cells
                                       'CD24','MS4A1','CD19', #b cell
                                       'CD22','CCR6','CR2', # class-switched/ memory B
                                       'TCL1A','LAIR1', #naive
                                       'IL2RA','BLK', #activated
                                       'MLANA','CD8A'
),
group.by = 'integrated_snn_res.0.3', dot.scale = 5)+ scale_color_viridis() + coord_flip()+
  theme(axis.text.y = element_text(size=7))

DotPlot(seu,assay = 'RNA',features = c('CD27','SLAMF7','IRF4','CD38','PRDM1', #plasma cell
                                       'IGHM',# immature B cell (before entering lymph node)
                                       'IGHD', # mature naive B cell (IGM and IGD; in lymph node)
                                       'IGHG1','IGHG2','IGHG3','IGHG4', 'IGHE', # class-switched/ plasma cell
                                       'IGHA1','IGHA2',
                                       'SDC1','XBP1', 'CD207','TUBB','STMN1','TYMS', #plasmablast
                                       'CXCR5','TNFRSF13B', #Follicular B cells
                                       'POU2AF1', 'CD40', 'SUGCT', #Germinal center B cells
                                       'CD24','MS4A1','CD19', #b cell
                                       'CD22','CCR6','CR2', # class-switched/ memory B
                                       'TCL1A','LAIR1', #naive
                                       'IL2RA','BLK' #activated
),
group.by = 'cell_type_fine', dot.scale = 5)+ scale_color_viridis() + coord_flip()+RotatedAxis()+
  theme(axis.text.y = element_text(size=7))

DotPlot(seu,assay = 'RNA',features = c('CD27','SLAMF7','IRF4','CD38','PRDM1', #plasma cell
                                       'IGHM',# immature B cell (before entering lymph node)
                                       'IGHD', # mature naive B cell (IGM and IGD; in lymph node)
                                       'IGHG1','IGHG2','IGHG3','IGHG4', 'IGHE', # class-switched/ plasma cell
                                       'IGHA1','IGHA2',
                                       'SDC1','XBP1', 'CD207','TUBB','STMN1','TYMS', #plasmablast
                                       'CXCR5','TNFRSF13B', #Follicular B cells
                                       'POU2AF1', 'CD40', 'SUGCT', #Germinal center B cells
                                       'CD24','MS4A1','CD19', #b cell
                                       'CD22','CCR6','CR2', # class-switched/ memory B
                                       'TCL1A','LAIR1', #naive
                                       'IL2RA','BLK' #activated
),
group.by = 'cell_type_fine', dot.scale = 5)+ scale_color_viridis() + coord_flip()+RotatedAxis()+
  theme(axis.text.y = element_text(size=7))

seu$cell_type_fine<-factor(seu$cell_type_fine,levels = c("Naive B cells","Memory B cells",
                                                         "Class-switched memory B cells",
                                                         "Plasma cells","Cycling cells"))
DotPlot(seu,assay = 'RNA',features = c('TCL1A','IGHD',
                                       'IGHM','CD22','CCR6','CR2','IL2RA',
                                       'BLK','MS4A1',
                                       'IGHG1','IGHG2','IGHG3','IGHG4', 'IGHE', # class-switched/ plasma cell
                                       'IGHA1','IGHA2','CD38','PRDM1',
                                       'TUBB','MKI67'
),
group.by = 'cell_type_fine', dot.scale = 5)+ scale_color_viridis() + coord_flip()+RotatedAxis()+
  theme(axis.text.y = element_text(size=7))


DotPlot(seu,assay = 'RNA',features = c(paste0(names(sigs),'1')),
        group.by = 'integrated_snn_res.0.3', dot.scale = 5)+ scale_color_viridis() + coord_flip()+
  theme(axis.text.y = element_text(size=7))
DotPlot(seu,assay = 'RNA',features = c(paste0(names(sigs),'1')),
        group.by = 'cell_type_fine', dot.scale = 5)+ scale_color_viridis() + coord_flip()+
  theme(axis.text.y = element_text(size=7))+RotatedAxis()

print(FeaturePlot(seu, features = c('IGHM','IGHD','IGHA1','IGHA2'), 
                  min.cutoff = "q05",max.cutoff = "q95", order = T,raster = T)& 
        scale_colour_gradientn(colours = rev(brewer.pal(11,'Spectral'))))
print(FeaturePlot(seu, features = c('IGHG1','IGHG2','IGHG3','IGHG4'), 
                  min.cutoff = "q05",max.cutoff = "q95", order = T,raster = T)& 
        scale_colour_gradientn(colours = rev(brewer.pal(11,'Spectral'))))
print(FeaturePlot(seu, features = c('IGHE'), 
                  min.cutoff = "q05",max.cutoff = "q95", order = T,raster = T)& 
        scale_colour_gradientn(colours = rev(brewer.pal(11,'Spectral'))))

print(FeaturePlot(seu, features = c('MLANA','VWF','CD3E','LYZ'), 
                  min.cutoff = "q05",max.cutoff = "q95", order = T,raster = T)& 
        scale_colour_gradientn(colours = rev(brewer.pal(11,'Spectral'))))

for(i in seq(1, ncol(sigs), 4)){
  four<-paste0(colnames(sigs)[i:(i+3)],'1')
  print(FeaturePlot(seu, features = c(four),min.cutoff = "q05", 
                    max.cutoff = "q95",order=T,raster=T)& 
          scale_colour_gradientn(colours = rev(brewer.pal(11,'Spectral'))))
}

VlnPlot(seu,c("nCount_RNA", "nFeature_RNA", 'percent.mt', "percent.rps", "percent.rpl", 'G2M.Score', "S.Score"),
        pt.size = 0)+NoLegend()

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


pdf(paste0(path.ct,'plots_',filename,'_v3_quantification.pdf'))
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

pdf(file = paste0(path.ct, 'markers_', filename, '_heatmap_celltype.pdf'), width = 10,
    height = 18)
DoHeatmap(seu, features = tops$gene, group.by = 'cell_type_fine', raster = T,
          assay = 'RNA')
dev.off()

