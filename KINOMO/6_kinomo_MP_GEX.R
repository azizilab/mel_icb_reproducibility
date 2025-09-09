#!/usr/bin/env Rscript

#### Title: Sanity check of gene expression of meta programs (MPs)
#### Author: Jana Biermann, PhD
#### Based on concept by Somnath Tagore, PhD


# Normalized Gene Contribution
# —————————————————
# 1. Sort the metagenes high to low based on stouffer integrated values (use top 100)
# 2. From the H matrices generated sample wise (using kinomo), identify the cell barcodes associated with each factor
# 3. Subset the barcodes from integrated single cell data based on factor barcodes (per metaprogram)
# 4. run a signature built using the ranked metagenes per metaprogram on the single cell data
#   a. on the individual factor
#   b. on all the factors associated with a metaprogram
# 5. Run EM-GMM model on the single cell (tpm) and identify the modality of the distribution per gene
# 6. Based on the modality (could be unimodal, bimodal or multimodal), identify the peak center, which we define as "normalized gene contribution"
# 7. Generate a matrix using the normalized gene contribution for each metagene per metaprogram
# 8. Generate the checkerboard


library(dplyr)
library(Seurat)
library(ggplot2)
library(gplots)
library(hypeR)
library(viridis)
library(pheatmap)
library(scales)
library(RColorBrewer)
library(hypeR)
library(msigdbr)
library(maptree)
library(reshape2)
'%notin%' <- Negate('%in%')

directory<-'data/KINOMO/'
label<-'melanoma'

setwd('~/Documents/melanoma')

ranks_tab<-read.csv('data/KINOMO/kinomo_melanoma_samples.csv')

pats<-c('F01_on','F01_pre','F02_on','F02_pre','F03_post1_on2','F03_post1_pre2','F04_pre',
        'F05_pre','F06_post1_pre2','F07_post1_pre2','F08_post','F09_post',
        'F12_post','F12_pre','F15_pre','F16_post1_pre2','F16_pre','F17_post',
        'F18_post','F20_post1_pre2','F22_post','F23_post',
        'F26_post','F27_post1_pre2','F28_post1_pre2','F29_post1_pre2','F30_post',
        'F31_post','R204_pre','R294_on','R308_pre','R310_on1','R310_on2','R310_pre',
        'R319_on','R319_pre','R328_on','R329_on','R334_pre','R354_pre')



#### Apply to data #####
# Download H files
for(pat in pats){
  ranks<-as.numeric(unlist(ranks_tab[ranks_tab$sample==pat,2:3]))
  for(rank in ranks){
    system(paste0("aws s3 sync s3://melanoma-ribas/KINOMO/",pat,"/ data/KINOMO/",pat,"/ --exclude '*' --include '*KINOMO_nmf_rank_",rank,"_H.csv' "))
  }
}

# Assign each barcode/cell to factor in individual sample
ranklabel<-'best'
top_number<-100
clustering_method<-'ward.D2'
meta_number<-11
MP_factors<-read.csv(paste0(directory,ranklabel, '/top',top_number,'/',clustering_method,
                            '/table_MP_factors_',ranklabel,'_top',top_number,"_",clustering_method,"_",meta_number,'MPs.csv'))

for(pat in pats){
  rank<-as.numeric(unlist(ranks_tab[ranks_tab$sample==pat,2]))
  H <- read.csv(paste0(directory,pat,"/",pat,"_KINOMO_nmf_rank_",rank,"_H.csv"),row.names = 1)
  H <- as.data.frame(t(H))
  top_H<-data.frame(barcode=gsub(pattern = '.',replacement = '-',x = rownames(H),fixed=T))
  top_H$top_factor <- max.col(H, "first")
  write.csv(top_H,paste0(directory,pat,"/",pat,"_KINOMO_nmf_rank_",rank,"_topH.csv"),row.names = F)
  print(pat)
  print(table(top_H$top_factor))
}


# Merge barcodes from individual samples into MPs
barcodes<-c()
for(m in 1:meta_number){
  tmp_factors<- subset(MP_factors, MP==m,sample_factor)
  tmp_barcodes<-c()
  for(factor in unlist(tmp_factors)){
    pat<-substr(factor,1,(nchar(factor)-8))
    rank<-as.numeric(unlist(ranks_tab[ranks_tab$sample==pat,2]))
    top_H<-read.csv(paste0(directory,pat,"/",pat,"_KINOMO_nmf_rank_",rank,"_topH.csv"))
    factor_number<-as.numeric(substr(factor,(nchar(factor)),nchar(factor)))
    tmp_barcodes<-c(tmp_barcodes,subset(top_H,top_factor==factor_number,barcode))
  }
  tmp_barcodes<-unlist(tmp_barcodes,recursive = T,use.names = F)
  barcodes[[m]]<-tmp_barcodes
  names(barcodes)[m]<-paste0('MP',m)
}
barcodes_df<-sapply(barcodes, function(x){
  c(x, rep(NA, max(sapply(barcodes,length)) - length(x)))})
write.csv(barcodes_df,paste0(directory,ranklabel, '/top',top_number,'/',clustering_method,
                       "/table_barcodes_",ranklabel,"_top",top_number,"_",clustering_method,"_",
                       meta_number,"MPs.csv"),row.names = F)


# Check gene expression of top 100 MP genes
seu<-readRDS('data/tumor/data_tum_merged_v5.rds')
MP_top100_df<-read.csv(paste0(directory,ranklabel, '/top',top_number,'/',clustering_method,
                       "/table_metaprograms_",ranklabel,"_top",top_number,"_",clustering_method,"_",
                       meta_number,"MPs_stouffer.csv"))

# Add MP groups to barcodes
seu$bc_group<-NA
for(m in 1:meta_number){
  seu$bc_group<-ifelse(seu$barcode_all %in% barcodes[[m]],paste0('MP',m),seu$bc_group)
}

# Apply MPs as signatures
for (m in 1:meta_number) {
  tmp_genes<-MP_top100_df[,m]
  seu <- AddModuleScore(object = seu, features = list(tmp_genes), 
                        name = paste0('MP',m), assay = 'RNA', search = F)
}


# Plot violins and dotplots
pdf(paste0(directory,ranklabel, '/top',top_number,'/',clustering_method,
           '/plots_DotPlots_',ranklabel,'_top',top_number,"_",clustering_method,"_",meta_number,'MPs.pdf'))

VlnPlot(seu, features=c(paste0('MP',seq(1:meta_number), '1')),pt.size=0,
        group.by = 'bc_group',stack=T,flip=T)+NoLegend()

for(m in 1:meta_number){
  tmp_genes<-MP_top100_df[,m]
  seu$bc_binary<-ifelse(seu$barcode_all %in% barcodes[[m]],paste0('MP',m),paste0('not_MP',m))
  print(DotPlot(seu,features = tmp_genes,group.by = 'bc_group',dot.scale = 2)+
    coord_flip()+
    scale_color_viridis()+
    ggtitle(paste0('MP',m))+
    theme(axis.text.y = element_text(colour = 'black',face = 'italic',size=5),
          plot.margin = margin(0.3, 4, 0.3, 0.3, 'cm')))
  
  print(DotPlot(seu,features = tmp_genes,group.by = 'bc_binary',dot.scale = 2,scale = F)+
    coord_flip()+
    scale_color_viridis()+
    ggtitle(paste0('MP',m,' binary'))+
    theme(axis.text.y = element_text(colour = 'black',face = 'italic',size=5),
          plot.margin = margin(0.3, 7, 0.3, 0.3, 'cm')))
  
  print(FeaturePlot(seu,paste0('MP',m,'1'))+
          scale_colour_gradientn(colours = rev(brewer.pal(11,'Spectral')))
  )
}
dev.off()


# Plot heatmaps
set.seed(1)
for(m in 1:meta_number){
  tmp_genes<-MP_top100_df[,m]
  sel_genes<-unlist(tmp_genes,recursive = T,use.names = F)
  seu_sel<-subset(seu,features = sel_genes,cells = sample(rownames(seu@meta.data),25000))
  seu_sel<-ScaleData(seu_sel,features=sel_genes)
  
  datamat<-as.matrix(seu_sel@assays$RNA@scale.data)
  datamat[datamat>3]<-3
  datamat[datamat< -3]<- (-3)
  datamat<- rescale(datamat, to=c(-1,1))
  datamat <- datamat[, names(sort(seu_sel$bc_group,decreasing = F))]
  
  # cluster anno
  annotation_col = data.frame(MPs = sort(seu_sel$bc_group,decreasing = F))
  rownames(annotation_col) = colnames(datamat)
  
  pdf(paste0(directory,ranklabel, '/top',top_number,'/',clustering_method,
             '/plots_heatmap_',ranklabel,'_top',top_number,"_",clustering_method,"_",
             meta_number,'MPs_MP',m,'.pdf'),height = 14,width = 10)
  print(pheatmap(datamat,
                 main = paste0('MP',m),
                 cluster_rows = F, cluster_cols = F, 
                 color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100),
                 show_colnames = F, show_rownames = T,
                 annotation_col = annotation_col))
  dev.off()
}  


sel_genes<-unique(unlist(MP_top100_df,recursive = T,use.names = F))
seu_sel<-subset(seu,features = sel_genes)

datamat<- rescale(seu_sel@assays$RNA@scale.data, to=c(-1,1))
datamat = datamat[match(sel_genes, rownames(datamat)), ]
datamat = datamat[, order(seu_sel$bc_group,decreasing = F)]

# marker anno
annotation_row = data.frame(MP = c(sort(rep(paste0('MP',seq(1:meta_number)),100))))
rownames(annotation_row) = rownames(datamat)

# cluster anno
cluster_anno = seu_sel$bc_group
MP1_sig = seu_sel$MP11
MP1_sig = MP1_sig[order(seu_sel$bc_group,decreasing = F)]
MP1_sig<- rescale(MP1_sig, to=c(0,1))

annotation_col = data.frame(cluster = factor(cluster_anno),MP1_sig=MP1_sig)
rownames(annotation_col) = colnames(datamat)

pdf(paste0(directory,ranklabel, '/top',top_number,'/',clustering_method,
           '/plots_heatmaps_',ranklabel,'_top',top_number,"_",clustering_method,"_",meta_number,'MPs.pdf'),height=17)
pheatmap(datamat,cluster_rows = F, cluster_cols = F, color = PurpleAndYellow(),
         show_colnames = F, show_rownames = F,annotation_col = annotation_col)
dev.off()

