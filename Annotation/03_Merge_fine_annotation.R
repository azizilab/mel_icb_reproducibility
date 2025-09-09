#!/usr/bin/Rscript

#### Title: Merge cell type info
#### Author: Jana Biermann, PhD

library(dplyr)
library(ggplot2)
library(ggrepel)
library(scales)
library(Seurat)
library(reshape2)
library(viridis)
library(RColorBrewer)
library(ggpubr)
library(scater)
library(tidyr)
'%notin%' <- Negate('%in%')

#setwd('~/Documents/melanoma/')

directory<-'data/integrated/'
label<-'melanoma'

system("aws s3 cp s3://melanoma-ribas/Seurat/integrated/data_melanoma_integrated_RPCA.rds data/integrated/")
seu<-readRDS('data/integrated/data_melanoma_integrated_RPCA.rds')

system("aws s3 cp s3://melanoma-ribas/Seurat/annotation/endothelial/data_melanoma_endothelial_celltype.csv data/annotation/endothelial/")
system("aws s3 cp s3://melanoma-ribas/Seurat/annotation/stromal/data_melanoma_stromal_celltype.csv data/annotation/stromal/")
system("aws s3 cp s3://melanoma-ribas/Seurat/annotation/bcells/data_melanoma_bcells.csv data/annotation/bcells/")
system("aws s3 cp s3://melanoma-ribas/Seurat/annotation/main/data_melanoma_celltype_main_nontum.csv data/annotation/main/")
system("aws s3 cp s3://melanoma-ribas/Seurat/annotation/myeloid/data_melanoma_myeloid.csv data/annotation/myeloid/")
system("aws s3 cp s3://melanoma-ribas/Seurat/annotation/tcells/data_melanoma_tcells.csv data/annotation/tcells/")
system("aws s3 cp s3://melanoma-ribas/Seurat/annotation/tumor/data_melanoma_celltype_fine_tum.csv data/annotation/tumor/")
system("aws s3 cp s3://melanoma-ribas/Seurat/annotation/tumor/data_melanoma_celltype_main_tum.csv data/annotation/tumor/")


#### cell_type_main and fine ####

# Combine ct info
ct_ec<-read.csv('data/annotation/endothelial/data_melanoma_endothelial_celltype.csv',row.names = 1)
ct_fb<-read.csv('data/annotation/stromal/data_melanoma_stromal_celltype.csv',row.names = 1)
ct_b<-read.csv('data/annotation/bcells/data_melanoma_bcells.csv',row.names = 1)
ct_my<-read.csv('data/annotation/myeloid/data_melanoma_myeloid.csv',row.names = 1)
ct_t<-read.csv('data/annotation/tcells/data_melanoma_tcells.csv',row.names = 1)
ct_tum<-read.csv('data/annotation/tumor/data_melanoma_celltype_fine_tum.csv',row.names = 1)
ct_fine<-rbind.data.frame(ct_ec,ct_fb,ct_b,ct_my,ct_t,ct_tum)

ct_main<-read.csv('data/annotation/main/data_melanoma_celltype_main_nontum.csv',row.names = 1)
ct_main$cell_type_main<-ifelse(ct_main$cell_type_main == ' pDCs','pDCs',ct_main$cell_type_main)
ct_main_tum<-read.csv('data/annotation/tumor/data_melanoma_celltype_main_tum.csv',row.names = 1)[,1:2]
ct_main<-rbind.data.frame(ct_main,ct_main_tum)

# add main ct
seu$cell_type_main<-NULL
seu$barcode_all<-row.names(seu@meta.data)
seu@meta.data<-left_join(seu@meta.data,ct_main,by='barcode_pat')
row.names(seu@meta.data)<-seu$barcode_all

table(seu$cell_type_main,useNA='always')
print(DimPlot(seu, label = T,group.by = 'cell_type_main',shuffle = T,raster = T,repel = T))
print(DimPlot(seu, label = F,group.by = 'cell_type_main',shuffle = T,raster = T,repel = T,
              split.by = 'cell_type_main',ncol = 3))

# add fine ct
seu$cell_cycle<-NULL
seu$cell_type_fine<-NULL
seu$barcode_all<-row.names(seu@meta.data)
seu@meta.data<-left_join(seu@meta.data,ct_fine,by='barcode_pat')
row.names(seu@meta.data)<-seu$barcode_all
seu$cell_type_fine <- ifelse(seu$cell_type_main %in% c('Mast cells'),'Mast cells', seu$cell_type_fine)
seu$cell_type_main <- ifelse(seu$cell_type_main %in% c('Mast cells') & seu$cell_type_fine %in% c('Mast cells'),'Myeloid cells', seu$cell_type_main)
seu$cell_type_fine <- ifelse(seu$cell_type_main %in% c('Epithelial cells'),'Epithelial cells', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$cell_type_main %in% c('pDCs'),'pDCs', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$cell_type_main %in% c('Tumor cells') & seu$cell_cycle =='cycling','Cycling tumor cells', seu$cell_type_fine)
seu$cell_type_main<- ifelse(is.na(seu$cell_type_fine)==T,NA,seu$cell_type_main)
table(seu$cell_type_main,useNA='always')

table(seu$cell_type_fine,useNA='always')
print(DimPlot(seu, label = T,group.by = 'cell_type_fine',shuffle = T,raster = T,repel = T,label.size = 2)+
        NoLegend())


# add clin info
system('aws s3 cp s3://melanoma-ribas/Seurat/clin.csv data/')
clin<-read.csv('data/clin.csv',na.strings = '')
name_clin<-names(clin)[2:ncol(clin)]
seu@meta.data[,name_clin]<-NULL
seu$barcode_all<-row.names(seu@meta.data)
seu@meta.data<-left_join(seu@meta.data,clin,by='sample')
row.names(seu@meta.data)<-seu$barcode_all
seu$treated<-ifelse(seu$group %in% c('pre_naive','pre_naive_muc'),'Untreated','Treated')
seu$mucosal<-ifelse(seu$group %in% c('post_PD-1_muc','pre_naive_muc'),'Mucosal','Non-mucosal')
seu$responder<-ifelse(seu$best_response %in% c('CR','PR'),'Responder','Non-responder')


#### cell_type_int ####
# T cells
seu$cell_type_int <- ifelse(seu$cell_type_main %in% c('T/NK cells') & 
                              seu$cell_type_fine %in% c('Naive CD8+ T cells','CD8+ T cells'),
                            'CD8+ T cells', NA)
seu$cell_type_int <- ifelse(seu$cell_type_main %in% c('T/NK cells') & 
                              seu$cell_type_fine %in% c('CD4+ T cells','Naive CD4+ T cells'),
                            'CD4+ T cells', seu$cell_type_int)
seu$cell_type_int <- ifelse(seu$cell_type_main %in% c('T/NK cells') & 
                              seu$cell_type_fine %in% c('Tregs','Naive Tregs'),
                            'Tregs', seu$cell_type_int)
seu$cell_type_int <- ifelse(seu$cell_type_main %in% c('T/NK cells') & 
                              seu$cell_type_fine %in% c('NK cells 1','NK cells 2'),
                            'NK cells', seu$cell_type_int)
seu$cell_type_int <- ifelse(seu$cell_type_main %in% c('T/NK cells') & 
                              seu$cell_type_fine %in% c('Cycling cells'),
                            'Cycling cells', seu$cell_type_int)

# myeloid cells
seu$cell_type_int <- ifelse(seu$cell_type_main %in% c('Myeloid cells') & 
                              seu$cell_type_fine %in% c('Monocytes'),
                            'Monocytes', seu$cell_type_int)
seu$cell_type_int <- ifelse(seu$cell_type_main %in% c('Myeloid cells') & 
                              seu$cell_type_fine %in% c('cDC1', 'cDC2','DC3'),
                            'Dendritic cells', seu$cell_type_int)
seu$cell_type_int <- ifelse(seu$cell_type_main %in% c('Myeloid cells') & 
                              seu$cell_type_fine %in% c('Macrophages (antigen-pres.)','Macrophages (M1-like)','Macrophages (M2-like)','Macrophages (phagocytic)'),
                            'Macrophages', seu$cell_type_int)
seu$cell_type_int <- ifelse(seu$cell_type_main %in% c('Myeloid cells') & 
                              seu$cell_type_fine %in% c('Macrophages (cycling)'),'Cycling cells', seu$cell_type_int)
seu$cell_type_int <- ifelse(seu$cell_type_main %in% c('Myeloid cells') & 
                              seu$cell_type_fine %in% c('Mast cells'),
                            'Mast cells', seu$cell_type_int)

# b cells
seu$cell_type_int <- ifelse(seu$cell_type_main %in% c('B/Plasma cells') & 
                              seu$cell_type_fine %in% c('Plasma cells'),'Plasma cells', seu$cell_type_int)
seu$cell_type_int <- ifelse(seu$cell_type_main %in% c('B/Plasma cells') & 
                              seu$cell_type_fine %in% c('Class-switched memory B cells','Memory B cells','Naive B cells'),
                            'B cells', seu$cell_type_int)
seu$cell_type_int <- ifelse(seu$cell_type_main %in% c('B/Plasma cells') & 
                              seu$cell_type_fine %in% c('Cycling cells'),'Cycling cells', seu$cell_type_int)


# stromal, epithelial, endothelial, tumor
seu$cell_type_int <- ifelse(seu$cell_type_main %in% c('Stromal cells') & 
                              seu$cell_type_fine %in% c('Pericytes'),
                            'Pericytes', seu$cell_type_int)
seu$cell_type_int <- ifelse(seu$cell_type_main %in% c('Stromal cells') & 
                              seu$cell_type_fine %in% c('CAFs','CAFs (antigen-pres.)','CAFs (cytokines)','CAFs (ECM)'),
                            'CAFs', seu$cell_type_int)
seu$cell_type_int <- ifelse(seu$cell_type_main %in% c('Stromal cells') & 
                              seu$cell_type_fine %in% c('Cycling cells'),
                            'Cycling cells', seu$cell_type_int)

seu$cell_type_int <- ifelse(seu$cell_type_main %in% c('Epithelial cells') & 
                              seu$cell_type_fine %in% c('Epithelial cells'),
                            'Epithelial cells', seu$cell_type_int)

seu$cell_type_int <- ifelse(seu$cell_type_main %in% c('Endothelial cells'),
                            'Endothelial cells', seu$cell_type_int)
seu$cell_type_int <- ifelse(seu$cell_type_main %in% c('Endothelial cells')& 
                              seu$cell_type_fine %in% c('Cycling cells'),
                            'Cycling cells', seu$cell_type_int)

seu$cell_type_int <- ifelse(seu$cell_type_main %in% c('pDCs'),
                            'Dendritic cells', seu$cell_type_int)
seu$cell_type_int <- ifelse(seu$cell_type_main %in% c('Tumor cells'),
                            'Tumor cells', seu$cell_type_int)

print(DimPlot(seu, label = T,group.by = 'cell_type_int',shuffle = T,raster = T,repel = T))


# add immune info
seu$ct_immune<-ifelse(seu$cell_type_main %in% c('T/NK cells','B/Plasma cells','Myeloid cells',
                                                'pDCs'),'Immune','Non-immune')
print(DimPlot(seu, label = T,group.by = 'ct_immune',shuffle = T,raster = T,repel = T))

print(DimPlot(seu, label = T,group.by = 'cell_cycle',shuffle = T,raster = T,repel = T))
print(DimPlot(seu, label = F,group.by = 'cell_cycle',shuffle = T,raster = T,repel = T,split.by = 'cell_cycle'))

# Save object
saveRDS(seu, 'data/integrated/data_melanoma_integrated_RPCA_v2.rds')
ct <- seu@meta.data %>% dplyr::select('barcode_pat', 'cell_type_main', 'cell_type_int', 'cell_type_fine','cell_cycle')
write.csv(ct, 'data/annotation/main/data_melanoma_celltype_annotation.csv',row.names = T)


### stats
stats<-as.data.frame(matrix(data=NA,nrow = 1, ncol = 11))
colnames(stats)<-c('sample','n_features','n_cells','median_features','median_counts')
rownames(stats)<-label
stats$sample<-label
stats$n_features<-dim(seu@assays$integrated@data)[1]
stats$n_cells<-dim(seu@assays$integrated@data)[2]
stats$median_features<-round(median(seu@meta.data$nFeature_RNA))
stats$median_counts<-round(median(seu@meta.data$nCount_RNA))

pdf(paste0("data/integrated/plots_melanoma_integrated_RPCA_v2.pdf"))
textplot(t(stats),cex=1.2,halign='left')
ElbowPlot(seu,ndims = 100)
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

system("aws s3 sync data/integrated/ s3://melanoma-ribas/Seurat/integrated/ --dryrun")
system("aws s3 sync data/integrated/ s3://melanoma-ribas/Seurat/integrated/ ")

#### v3 (remove NAs) ####
# remove NAs
table(seu$cell_type_fine,useNA='always')
table(seu$cell_type_fine,seu$cell_type_main,useNA='always')
seu$cell_type_fine<-ifelse(is.na(seu$cell_type_fine) ==T, 'remove',seu$cell_type_fine)
seu<-subset(seu,cell_type_fine!='remove')

seu <- ScaleData(object = seu)
seu <- RunPCA(object = seu,npcs = 75)
ElbowPlot(seu,ndims = 75)
seu <- RunUMAP(object = seu, dims = 1:30)
seu <- FindNeighbors(seu, dims = 1:30)
seu <- FindClusters(seu)

# Save object
saveRDS(seu, 'data/integrated/data_melanoma_integrated_RPCA_v3.rds')


pdf(paste0("data/integrated/plots_melanoma_integrated_RPCA_v3.pdf"))
textplot(t(stats),cex=1.2,halign='left')
ElbowPlot(seu,ndims = 100)
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



