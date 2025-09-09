#!/usr/bin/env Rscript

#### Fine cell type annotation of endothelial cells
#### Author: Jana Biermann, PhD

library(dplyr)
library(Seurat)
library(ggplot2)
library(gplots)
library(viridis)
library(RColorBrewer)
library(ggpubr)
'%notin%' <- Negate('%in%')

my_palette<-c("lightgray" ,"#3288BD", "#66C2A5", "#ABDDA4", "#E6F598" ,"yellow", "#FEE08B" ,
              "#FDAE61", "#F46D43","#D53E4F" ,"#9E0142")

setwd('~/Documents/melanoma/')
label<-'melanoma'
celltype<-'endothelial'
filename<-paste0(label,'_',celltype)
path.ct <- paste0('data/annotation/',celltype,'/')
ifelse(!dir.exists(file.path(path.ct)), dir.create(file.path(path.ct),recursive = T), FALSE)

#system("aws s3 sync s3://melanoma-ribas/Seurat/annotation/main/ data/annotation/main/ --exclude '*' --include 'data_nontum_merged_v3.rds' ")
seu<-readRDS('data/annotation/main/data_nontum_merged_v3.rds')
seu$cell_type_main<-ifelse(seu$cell_type_main == ' pDCs','pDCs',seu$cell_type_main)

# Subset to neuronal cells
seu <- subset(seu, cell_type_main %in% c('Endothelial cells'))
seu<-NormalizeData(seu) %>% FindVariableFeatures()
seu <- ScaleData(object = seu)
seu <- RunPCA(object = seu,npcs = 75)
ElbowPlot(seu,ndims = 75)
seu <- FindNeighbors(seu, dims = 1:30)
seu <- FindClusters(seu, resolution = 0.4)
seu <- RunUMAP(object = seu, dims = 1:30)

# get module scores for sigs
sigs<-read.csv('~/brain_mets/signatures/endothelial.csv',na.strings = '')
for(c in 1:ncol(sigs)){
  seu<-AddModuleScore(object = seu,features = list(na.omit(sigs[,c])),
                      name = colnames(sigs)[c],assay = 'RNA',search=F)
}

# Cell cycle assignment
DefaultAssay(seu) <- 'RNA'
seu <- CellCycleScoring(seu, s.features = cc.genes.updated.2019$s.genes, 
                        g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = F)
DefaultAssay(seu) <- 'integrated'
seu$cell_cycle <- ifelse(seu$G2M.Score > 0.12 | seu$S.Score > 0.12, 'cycling', 'non-cycling')

# Integration
sort(table(seu$sample))
obj.list <- SplitObject(seu, split.by = 'sample')
n<-30
d<-20
k<-30
has_enough_cells <- function(seurat_obj) {
  return(dim(seurat_obj)[2] >= n)
}

# Filter the list to keep only objects with at least 20 cells
obj.list <- Filter(has_enough_cells, obj.list)

# Find integration anchors
anchors <- FindIntegrationAnchors(object.list = obj.list, dims = 1:d, k.filter = k)

# Integrate data sets
seu <- IntegrateData(anchorset = anchors, dims = 1:d,k.weight=k)

# Rerun Seurat workflow
seu <- ScaleData(seu) %>% 
  RunPCA() %>% 
  RunUMAP(dims=1:d) %>%
  FindNeighbors(dims = 1:d) %>% 
  FindClusters(resolution = 0.2)

seu<-FindClusters(seu,resolution = 0.2)

seu$treated<-ifelse(seu$group %in% c('pre_naive','pre_naive_muc'),'Untreated','Treated')
seu$mucosal<-ifelse(seu$group %in% c('post_PD-1_muc','pre_naive_muc'),'Mucosal','Non-mucosal')
seu$responder<-ifelse(seu$best_response %in% c('CR','PR'),'Responder','Non-responder')

#seu<-subset(seu,integrated_snn_res.0.4 !=7 & tcr != 'TCR') #rerun int
#seu<-subset(seu,integrated_snn_res.0.2 !=7) 
seu <- ScaleData(seu) %>% 
  RunPCA() %>% 
  RunUMAP(dims=1:d)

# Save obj
saveRDS(seu,paste0(path.ct,'data_',filename,'_v2.rds'))

seu<-readRDS(paste0(path.ct,'data_',filename,'_v2.rds'))


### stats
stats<-as.data.frame(matrix(data=NA,nrow = 1, ncol = 11))
colnames(stats)<-c('sample','n_features','n_cells','median_features','median_counts')
rownames(stats)<-filename
stats$sample<-filename
stats$n_features<-dim(seu@assays$integrated@data)[1]
stats$n_cells<-dim(seu@assays$integrated@data)[2]
stats$median_features<-round(median(seu@meta.data$nFeature_RNA))
stats$median_counts<-round(median(seu@meta.data$nCount_RNA))

pdf(paste0(path.ct,'plots_',filename,'.pdf'))
textplot(t(stats),cex=1.2,halign='left')

print(DimPlot(seu, label = T,group.by = 'ident',
              shuffle = T,raster = T,pt.size = 2)+coord_fixed())
print(DimPlot(seu, label = F,group.by = 'patient',
              shuffle = T,raster = T,pt.size = 2)+coord_fixed())
print(DimPlot(seu, label = F,group.by = 'tissue',
              shuffle = T,raster = T,pt.size = 2)+coord_fixed())
print(DimPlot(seu, label = F,group.by = c('responder','mucosal','treated'),
              shuffle = T,raster = T,pt.size = 4,ncol = 2))
FeatureScatter(seu,'G2M.Score','S.Score',group.by = 'cell_cycle', pt.size = 0.01)
print(DimPlot(seu, label = F,group.by = 'cell_cycle',
              shuffle = T,raster = T,pt.size = 2)+coord_fixed())
DimPlot(seu, label = F, group.by = 'cell_type_fine', repel = T, 
        raster=T,pt.size = 2)+coord_fixed()
DimPlot(seu, label = F, group.by = 'cell_type_main', repel = T, 
        raster=T,pt.size = 2)+coord_fixed()

DimPlot(seu, label = F, group.by = 'malignant', repel = T, 
        raster=T,pt.size = 2)+coord_fixed()

DimPlot(seu, label = T,group.by = 'celltype_bped_fine',repel = T,
        label.size = 3,shuffle = T,raster = T,pt.size = 2) + 
  guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
  theme(legend.text=element_text(size=7))+coord_fixed()

DimPlot(seu, label = T,group.by = 'celltype_hpca_main',repel = T,
        label.size = 2.5,shuffle = T,raster = T,pt.size = 2) + 
  guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
  theme(legend.text=element_text(size=7))+coord_fixed()

print(FeaturePlot(seu, features = c("nCount_RNA", "nFeature_RNA", 'percent.mt', "doublet_scores"), 
                  min.cutoff = "q05",max.cutoff = "q95", order = T,raster = T,pt.size = 4))
print(FeaturePlot(seu, features = c("percent.rps", "percent.rpl", 'G2M.Score', "S.Score"), 
                  min.cutoff = "q05",max.cutoff = "q95", order = T,raster = T,pt.size = 4))

print(FeaturePlot(seu, features = c("proportion_scaled_cnv_avg", "MLANA", 'PTPRC', "PMEL"), 
                  min.cutoff = "q05",max.cutoff = "q95", order = T,raster = T,pt.size = 2))

print(FeaturePlot(seu, features = c("rna_VEGFC", "rna_KIAA1217", 'rna_MMRN1', "rna_IL1R1"), 
                  min.cutoff = "q05",max.cutoff = "q95", order = T,raster = T,pt.size = 4))

# endothelial
print(FeaturePlot(seu, features = c("rna_PECAM1", 'rna_VWF', "rna_CD34", 'rna_CDH5'), 
                  min.cutoff = "q05",max.cutoff = "q95", order = T,raster = T,pt.size = 2))

DotPlot(seu,assay = 'RNA',features = c('VEGFC', 'DLL4', 'EFNB2','NOTCH4','HEY2', #arterial
                                       'RGCC','KIAA1217','ARHGAP18','PLVAP','ANGPT2', #capillary
                                       'MMRN1', 'FLT4','SEMA3D','RELN', #lymphatic
                                       'IL1R1','NR2F2','CDH11','EPHB4','VCAN', #venous
                                       'MKI67',
                                       'COL1A1','COL3A1','COL6A1','LAMA2','MLANA'), #FB
        group.by = 'integrated_snn_res.0.2', dot.scale = 5)+ scale_color_viridis() + coord_flip()+
  theme(axis.text.y = element_text(size=7))

DotPlot(seu,assay = 'RNA',features = c('VEGFC', 'DLL4', 'EFNB2','NOTCH4','HEY2', #arterial
                                       'RGCC','KIAA1217','ARHGAP18','PLVAP','ANGPT2', #capillary
                                       'MMRN1', 'FLT4','SEMA3D','RELN', #lymphatic
                                       'IL1R1','NR2F2','CDH11','EPHB4','VCAN', #venous
                                       'MKI67',
                                       'COL1A1','COL3A1','COL6A1','LAMA2','MLANA'), #FB
        group.by = 'cell_type_fine', dot.scale = 5)+ scale_color_viridis() + 
  coord_flip()+RotatedAxis()+
  theme(axis.text.y = element_text(size=7))

VlnPlot(seu,c("nCount_RNA", "nFeature_RNA", 'percent.mt', "doublet_scores","percent.rps", "percent.rpl", 'G2M.Score', "S.Score"),
        pt.size = 0)+NoLegend()
VlnPlot(seu,c("nCount_RNA", "nFeature_RNA", 'percent.mt', "doublet_scores","percent.rps", "percent.rpl"),
        pt.size = 0,group.by = 'cell_type_fine')+NoLegend()

DotPlot(seu,assay = 'RNA',features = paste0(names(sigs),1),
        group.by = 'integrated_snn_res.0.2', dot.scale = 3)+ scale_color_viridis() + coord_flip()
DotPlot(seu,assay = 'RNA',features = paste0(names(sigs),1),
        group.by = 'cell_type_fine', dot.scale = 3)+ scale_color_viridis() + coord_flip()+RotatedAxis()

for(i in seq(1, ncol(sigs), 4)){
  four<-paste0(colnames(sigs)[i:(i+3)],'1')
  print(FeaturePlot(seu, features = c(four),min.cutoff = "q05", 
                    max.cutoff = "q95",order=T,raster=T,pt.size = 3)&
          scale_colour_gradientn(colours = my_palette))
}
dev.off()


# Fine cell type assignment after manual annotation based on DGE
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.2 %in% c(0,6), 'Capillary EC', 'NA')
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.2 %in% c(5), 'Arterial EC', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.2 %in% c(2), 'Lymphatic EC', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.2 %in% c(1,3,4), 'Venous EC', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.2 %in% c(8), 'Cycling cells', seu$cell_type_fine)


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


# Differential gene expression (DGE) based on types
Idents(seu)<-seu$cell_type_fine
markers <- FindAllMarkers(seu, only.pos = TRUE,assay='RNA', min.pct = 0.25, 
                          logfc.threshold = 0.25,test.use = 'MAST')
write.csv(markers, paste0(path.ct,'markers_',filename,'_celltype.csv'), row.names = F)
markers %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% top_n(15, wt=avg_log2FC) -> tops
DefaultAssay(seu)<-'RNA'
seu<-NormalizeData(seu) %>% ScaleData()
DefaultAssay(seu)<-'integrated'
pdf(file = paste0(path.ct,'markers_',filename,'_heatmap_celltype.pdf'),width = 10, height = 18)
DoHeatmap(seu, features = tops$gene, group.by = 'cell_type_fine',raster = T,assay = 'RNA')
dev.off()


# Save object and labels
#saveRDS(seu,paste0(path.ct,'data_',filename,'.rds'))
ct <- seu@meta.data %>% select('barcode_pat', 'cell_cycle', 'cell_type_fine')
write.csv(ct, paste0(path.ct,'data_',filename,'_celltype.csv'))

