#!/usr/bin/env Rscript

#### Title: Run EM-GMM model on the single cell and identify the modality of the distribution per gene
#### Author: Somnath Tagore, PhD; Jana Biermann, PhD

library(tidyverse)
library(reshape2)
library(matrixStats)
library(dplyr)
library(Seurat)
library(gplots)
library(ggplot2)
library(Matrix)
library(mixtools)
library(gridExtra)
library(pheatmap)
library(RColorBrewer)
library(hypeR)

directory<-'data/KINOMO/'
label<-'melanoma'

setwd('~/Documents/melanoma')

ranklabel<-'best'
top_number<-100
clustering_method<-'ward.D2'
meta_number<- 11


#### Run EM-GMM to get normalized GE ####

# load gene expression data
seu<-readRDS('data/tumor/data_tum_merged_v5.rds')
gene_exp<-seu@assays$RNA@data

#Mixture model implementation
#Select one gene at a time
MP_top100_df<-read.csv(paste0(directory,ranklabel, '/top',top_number,'/',clustering_method,
                              "/table_metaprograms_",ranklabel,"_top",top_number,"_",clustering_method,
                              "_",meta_number,"MPs_stouffer.csv"))
barcodes_df<-read.csv(paste0(directory,ranklabel, '/top',top_number,'/',clustering_method,
                             "/table_barcodes_",ranklabel,"_top",top_number,"_",clustering_method,
                             "_",meta_number,"MPs.csv"))

genes<-unique(unlist(MP_top100_df))
norm_exp<-matrix(NA,nrow = length(genes),ncol = meta_number)
rownames(norm_exp)<-genes
colnames(norm_exp)<-paste0('MP',1:meta_number)
for (m in 1:meta_number) {
  # set.seed(123456)
  # barcodes_down<-sample(na.omit(barcodes_df[,m]),size = 1000)
  # MP_gex<-gene_exp[,barcodes_down]
  MP_gex<-gene_exp[,na.omit(barcodes_df[,m])]
  
  pdf(paste0(directory,ranklabel, '/top',top_number,'/',clustering_method,
             '/plots_GMM_',ranklabel,'_top',top_number,"_",clustering_method,"_",
             meta_number,'MPs_MP',m,'.pdf'))
  for(g in genes){
    gex_g<-as.vector(MP_gex[g,])
    if(sum(gex_g==0)==length(gex_g)){
      norm_exp[g,m]<-0

    }else{
      set.seed(123456)
      gene_exp_gmm <- normalmixEM(gex_g, lambda=NULL, mu=NULL, sigma=sd(gex_g))
      
      x       <- with(gene_exp_gmm,seq(min(x),max(x),len=1000))
      pars    <- with(gene_exp_gmm,data.frame(comp=colnames(posterior), mu, sigma,lambda))
      em.df   <- data.frame(x=rep(x,each=nrow(pars)),pars)
      em.df$y <- with(em.df,lambda*dnorm(x,mean=mu,sd=sigma))
      print(ggplot(data.frame(x=gene_exp_gmm$x),aes(x,y=..density..)) + 
              geom_histogram(fill=NA,color="black")+
              geom_polygon(data=em.df,aes(x,y,fill=comp),color="grey50", alpha=0.5)+
              scale_fill_discrete("Component\nMeans",labels=format(em.df$mu,digits=3))+
              theme_bw()+
              ggtitle(g))
      
      #Select the component with highest Mean to be the average expression of that gene
      em.aggregate<-aggregate(em.df, list(em.df$comp), mean)
      norm_exp[g,m]<-max(em.aggregate$mu)
    }
  }
  dev.off()
}  
write.csv(norm_exp,paste0(directory,ranklabel, '/top',top_number,'/',clustering_method,
                 '/table_GMM_',ranklabel,'_top',top_number,"_",clustering_method,"_",meta_number,'MPs.csv'))


# Plot normalized heatmaps
MP_top100_df<-read.csv(paste0(directory,ranklabel, '/top',top_number,'/',clustering_method,
                              "/table_metaprograms_",ranklabel,"_top",top_number,"_",
                              clustering_method,"_",meta_number,"MPs_stouffer.csv"))
colnames(MP_top100_df)<-paste0('MP',1:meta_number)


pdf(paste0(directory,ranklabel, '/top',top_number,'/',clustering_method,
           '/plots_heatmapNormalized_',ranklabel,'_top',top_number,"_",clustering_method,"_",
           meta_number,'MPs.pdf'),height = 17,width = 4)
print(pheatmap(norm_exp,
               main = 'Normalized Gene Contribution',
               cluster_rows = F, cluster_cols = F, 
               fontsize_row=3,
               color = colorRampPalette(brewer.pal(n = 7, name ="Greys"))(100),
               show_colnames = T, show_rownames = T))
dev.off()


#### Generate list of final genes ####
head(norm_exp)
max_norm_MP<-data.frame(gene=rownames(norm_exp),max_norm_MP=paste0('MP',max.col(norm_exp)))
head(max_norm_MP)

# Keep final genes for each MP based on intersection with top normalized values
final_genes<-NULL
for (m in 1:meta_number) {
  final_genes[[paste0('MP',m)]]<-intersect(MP_top100_df[,m],max_norm_MP$gene[max_norm_MP$max_norm_MP==paste0('MP',m)])
}
final_genes_df<-sapply(final_genes, function(x){
  c(x, rep(NA, max(sapply(final_genes,length)) - length(x)))})
write.csv(final_genes_df,paste0(directory,ranklabel,'/top',top_number,'/',clustering_method,
                              '/table_finalgenes_',ranklabel,
                              '_top',top_number,"_",
                              clustering_method,"_",
                              meta_number,'MPs.csv'),row.names = F)

pdf(paste0(directory,ranklabel, '/top',top_number,'/',clustering_method,
           '/plots_heatmapNormalizedFinal_',ranklabel,'_top',top_number,"_",clustering_method,"_",
           meta_number,'MPs.pdf'),height = 8,width = 4)
print(pheatmap(norm_exp[unlist(final_genes),],
               main = 'Normalized Gene Contribution\n(final genes)',
               cluster_rows = F, cluster_cols = F, 
               border_color=NA,
               #cellwidth=8,cellheight = 8,
               scale = 'row',
               breaks = seq(0, 1, length.out = 10),
               fontsize_row=3,
               color = colorRampPalette(c("white", "black"))(10),
               show_colnames = T, show_rownames = T))
dev.off()

# Manually selected short list
final_genes_short<-list(MP1=c('MKI67','CIT','TOP2A'),
                              MP2=c('ST3GAL6','NELL1','DLC1','LRMDA'),
                              MP3=c('GAPDH','ALDOA','HLA-A','HLA-B','HLA-C'),
                              MP4=c('RACK1','GPNMB','GNAS','S100A6','S100B'),
                              MP5=c('NLGN1','MAML3','REV3L'),
                              MP6=c('GAS7','KLF6','ANXA2'),
                              MP7=c('SOX5','NRG3','NAV2'))

pdf(paste0(directory,ranklabel, '/top',top_number,'/',clustering_method,
           '/plots_heatmapNormalizedFinal_',ranklabel,'_top',top_number,"_",clustering_method,"_",
           meta_number,'MPs_short.pdf'))
print(pheatmap(norm_exp[unlist(final_genes_short),],
               main = 'Normalized Gene Contribution\n(final genes; short list)',
               cluster_rows = F, cluster_cols = F, 
               border_color=NA,
               #cellwidth=15,cellheight = 15,
               scale = 'row',
               breaks = seq(0, 1, length.out = 10),
               #fontsize_row=3,
               #color = colorRampPalette(brewer.pal(n = 7, name ="Greys"))(10),
               color = colorRampPalette(c("white", "black"))(10),
               show_colnames = T, show_rownames = T))
dev.off()


# Check pathways of final genes (long list)
HALLMARK <- msigdb_gsets(species = 'Homo sapiens', category = 'H')
KEGG <- msigdb_gsets(species = 'Homo sapiens', category = 'C2', subcategory = 'CP:KEGG')
REACTOME <- msigdb_gsets(species = 'Homo sapiens', category = 'C2', subcategory = 'CP:REACTOME')
CP <- msigdb_gsets(species = 'Homo sapiens', category = 'C2', subcategory = 'CP')
C5_BP <- msigdb_gsets(species = 'Homo sapiens', category = 'C5', subcategory = 'BP')
genesets_w <- as.list(as.data.frame(read.csv('~/brain_mets/signatures/wouters_mel_sigs.csv', na.strings = c('','NA', NA))))
genesets_w<-gsets$new(genesets_w,name = "Melanoma")

final_genes_df<-read.csv(paste0(directory,ranklabel,'/top',top_number,'/',clustering_method,
                                '/table_finalgenes_',ranklabel,
                                '_top',top_number,"_",
                                clustering_method,"_",
                                meta_number,'MPs.csv'))

pdf(paste0(directory,ranklabel, '/top',top_number,'/',clustering_method,
           '/plots_MP_pathwaysFinal_',ranklabel,"_top",top_number,"_",clustering_method,"_",meta_number,
           'MPs.pdf'),
    width = 15,height = 10)
for(m in 1:meta_number){
  h1<-hyp_dots(hypeR(as.character(na.omit(final_genes_df[,m])), genesets = genesets_w), 
               title = paste0('Wouters pathway selection\n',ranklabel,"; top",top_number,"; ",meta_number,' MPs (total)\nMP ',m), abrv = 80) + 
    theme_bw() + 
    theme(axis.title.y = element_blank(), 
          axis.text.y = element_text(color = 'black',size=5), 
          axis.text.x = element_text(color = 'black'), 
          panel.grid.minor = element_blank(),
          title = element_text(size=8)) + 
    scale_color_gradient(high = 'black', low = '#C51B8A')
  
  h2<-hyp_dots(hypeR(as.character(na.omit(final_genes_df[,m])), genesets = HALLMARK), 
               title = paste0('Hallmarks\n',ranklabel,"; top",top_number,"; ",meta_number,' MPs (total)\nMP ',m), abrv = 80) + 
    theme_bw() + 
    theme(axis.title.y = element_blank(), 
          axis.text.y = element_text(color = 'black',size=5), 
          axis.text.x = element_text(color = 'black'), 
          panel.grid.minor = element_blank(),
          title = element_text(size=8)) + 
    scale_color_gradient(high = 'black', low = '#C51B8A')
  
  h3<-hyp_dots(hypeR(as.character(na.omit(final_genes_df[,m])), genesets = KEGG), 
               title = paste0('KEGG\n',ranklabel,'\nMP ',m), abrv = 80) + 
    theme_bw() + 
    theme(axis.title.y = element_blank(), 
          axis.text.y = element_text(color = 'black',size=5), 
          axis.text.x = element_text(color = 'black'), 
          panel.grid.minor = element_blank(),
          title = element_text(size=8)) + 
    scale_color_gradient(high = 'black', low = '#C51B8A')
  
  h4<-hyp_dots(hypeR(as.character(na.omit(final_genes_df[,m])), genesets = REACTOME), 
               title = paste0('REACTOME\n',ranklabel,"; top",top_number,"; ",meta_number,' MPs (total)\nMP ',m), abrv = 80) + 
    theme_bw() + 
    theme(axis.title.y = element_blank(), 
          axis.text.y = element_text(color = 'black',size=5), 
          axis.text.x = element_text(color = 'black'), 
          panel.grid.minor = element_blank(),
          title = element_text(size=8)) + 
    scale_color_gradient(high = 'black', low = '#C51B8A')
  
  h5<-hyp_dots(hypeR(as.character(na.omit(final_genes_df[,m])), genesets = C5_BP), 
               title = paste0('C5_BP\n',ranklabel,'\nMP ',m), abrv = 80) + 
    theme_bw() + 
    theme(axis.title.y = element_blank(), 
          axis.text.y = element_text(color = 'black',size=5), 
          axis.text.x = element_text(color = 'black'), 
          panel.grid.minor = element_blank(),
          title = element_text(size=8)) + 
    scale_color_gradient(high = 'black', low = '#C51B8A')
  
  h6<-hyp_dots(hypeR(as.character(na.omit(final_genes_df[,m])), genesets = CP), 
               title = paste0('CP\n',ranklabel,"; top",top_number,"; ",meta_number,' MPs (total)\nMP ',m), abrv = 80) + 
    theme_bw() + 
    theme(axis.title.y = element_blank(), 
          axis.text.y = element_text(color = 'black',size=5), 
          axis.text.x = element_text(color = 'black'), 
          panel.grid.minor = element_blank(),
          title = element_text(size=8)) + 
    scale_color_gradient(high = 'black', low = '#C51B8A')
  
  print((h1+h2+h3)/(h4+h5+h6))
}
dev.off()


#### pathway dotplot to compare pathways between factors based on top 100 MP genes ####
ranklabel<-'best'
top_number<-100
clustering_method<-'ward.D2'
meta_number<- 7

MP_top100_df<-read.csv(paste0(directory,ranklabel,'/top',top_number,'/',clustering_method,
                              '/table_metaprograms_',ranklabel,
                              '_top',top_number,"_",
                              clustering_method,"_",
                              meta_number,'MPs_stouffer.csv'))

# genesets_merged<-gsets$new(
#   c(HALLMARK$list(), KEGG$list(),C5_BP$list(),
#     genesets_w$list()[setdiff(names(genesets_w$list()),c(names(HALLMARK$list()), names(KEGG$list()),names(C5_BP$list())))]),
#   name = "merged")

genesets_merged<-gsets$new(
  c(HALLMARK$list(), KEGG$list(),C5_BP$list(),
    genesets_w$list()[grep('Tirosh|jerby|rambow',x = names(genesets_w$list()),value = T)]),
  name = "merged")

## select top10 pathways per MP and collect output for all MP
collect_pws<-NULL
collect_res<-data.frame(label=NA,fdr=NA,overlap_fraction=NA,MP=NA)
for(m in 1:meta_number){
  h1<-hypeR(as.character(na.omit(MP_top100_df[,m])), genesets = genesets_merged)
  collect_pws_tmp<-h1$data %>% filter(fdr<0.05) %>% select(label) %>% top_n(n = 10)
  collect_pws<-c(collect_pws,collect_pws_tmp)
  collect_res_tmp<-h1$data %>% #filter(label %in% collect_pws) %>% 
    mutate(overlap_fraction=overlap/geneset, MP=paste0('MP',m)) %>% 
    select(label,fdr,overlap_fraction,MP)
  collect_res<-rbind.data.frame(collect_res,collect_res_tmp)
  
}
collect_pws<-unique(unlist(collect_pws))
collect_res<-collect_res[-1,]
collect_res$label<-substr(collect_res$label,1,70)

pdf(paste0(directory,ranklabel, '/top',top_number,'/',clustering_method,
           '/plots_MP_pathwaysFinalCombined_',ranklabel,"_top",top_number,"_",clustering_method,
           "_",meta_number,'MPs.pdf'))
print(ggplot(collect_res[collect_res$label %in% collect_pws,],aes(x=MP,y=label)) +
        geom_point(aes(size=-log10(fdr),color=overlap_fraction)) +
        #scale_color_gradientn(colors=viridis(n=100,option = 'B',direction = 1)) +
        scale_color_distiller(palette = 'YlOrRd',direction = 1)+
        labs(size="-log10(FDR)", color="Gene Ratio") +
        theme_bw() + 
        ggtitle(paste0('Top pathways\n',ranklabel,"; top",top_number,"; ",meta_number,' MPs (total)')) +
        theme(panel.grid.minor = element_blank(),
              #panel.grid.major = element_blank(),
              axis.text=element_text(size=10, colour = "black"),
              axis.text.y = element_text(size=6.5, colour = "black"),
              axis.text.x = element_text(angle = 90,size = 8),
              axis.title=element_blank(),
              panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black")))
dev.off()


