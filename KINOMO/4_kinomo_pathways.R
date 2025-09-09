#!/usr/bin/env Rscript

#### Title: Pathway analysis for NMF markers
#### Author: Jana Biermann, PhD

library(dplyr)
library(Seurat)
library(ggplot2)
library(gplots)
library(hypeR)
library(viridis)
'%notin%' <- Negate('%in%')

directory<-'data/KINOMO/'
label<-'melanoma'

setwd('~/Documents/melanoma')

ranks_tab<-read.csv('data/KINOMO/kinomo_melanoma_samples.csv')

pats<-c('F01_on','F01_pre','F02_on','F02_pre','F03_post1_on2','F03_post1_pre2','F04_pre',
        'F05_pre','F06_post1_pre2','F07_post1_pre2','F08_post','F09_post','F10_post',
        'F12_post','F12_pre','F15_pre','F16_post1_pre2','F16_pre','F17_post',
        'F18_post','F20_post1_pre2','F22_post','F23_post','F24_post','F25_post',
        'F26_post','F27_post1_pre2','F28_post1_pre2','F29_post1_pre2','F30_post',
        'F31_post','R204_pre','R294_on','R308_pre','R310_on1','R310_on2','R310_pre',
        'R319_on','R319_pre','R328_on','R329_on','R334_pre','R354_pre')

# Download gene lists
for(pat in pats){
  ranks<-as.numeric(unlist(ranks_tab[ranks_tab$sample==pat,2:3]))
  for(rank in ranks){
    system(paste0("aws s3 sync s3://melanoma-ribas/KINOMO/",pat,"/ data/KINOMO/",pat,"/ --exclude '*' --include '*KINOMO_nmf_rank_",rank,"_top100_W.csv' "))
  }
}


HALLMARK <- msigdb_gsets(species="Homo sapiens", category="H")
C5_BP <- msigdb_gsets(species="Homo sapiens", category="C5", subcategory="BP")
genesets_w <- as.list(as.data.frame(read.csv('~/brain_mets/signatures/wouters_mel_sigs.csv', na.strings = c('', NA))))

# Run pw analysis on individual samples, ranks, and factors
for(pat in pats){
  ranks<-as.numeric(unlist(ranks_tab[ranks_tab$sample==pat,2:3]))
  for(rank in ranks){
    genes<-read.csv(paste0("data/KINOMO/",pat,"/",pat,"_KINOMO_nmf_rank_",rank,"_top100_W.csv"),row.names = 1)
    pdf(paste0(directory, pat, '/plots_',pat,'_rank_',rank,'_top100_pathways.pdf'),width = 15,height = 7)
    for(i in 1:ncol(genes)){
      subs<-genes[,i]
      h1<-hyp_dots(hypeR(as.character(subs), genesets=HALLMARK),
                     title = paste0('Hallmarks:\n',colnames(genes)[i],';\nMax ranks: ',rank),
                     abrv=80)+ theme_bw() + 
              theme(axis.title.y = element_blank(), 
                    axis.text.y = element_text(color = 'black',size=5), 
                    axis.text.x = element_text(color = 'black'), 
                    panel.grid.minor = element_blank()) + 
              scale_color_gradient(high = 'black', low = '#C51B8A')
      
      h2<-hyp_dots(hypeR(as.character(subs), genesets=C5_BP),
                     title = paste0('GO BP:\n',colnames(genes)[i],';\nMax ranks: ',rank),
                     abrv=80)+ theme_bw() + 
              theme(axis.title.y = element_blank(), 
                    axis.text.y = element_text(color = 'black',size=5), 
                    axis.text.x = element_text(color = 'black'), 
                    panel.grid.minor = element_blank()) + 
              scale_color_gradient(high = 'black', low = '#C51B8A')
      
      h3<-hyp_dots(hypeR(as.character(subs), genesets=genesets_w),
                     title = paste0('Wouters pathways:\n',colnames(genes)[i],';\nMax ranks: ',rank),
                     abrv=80)+ theme_bw() + 
              theme(axis.title.y = element_blank(), 
                    axis.text.y = element_text(color = 'black',size=5), 
                    axis.text.x = element_text(color = 'black'), 
                    panel.grid.minor = element_blank()) + 
              scale_color_gradient(high = 'black', low = '#C51B8A')
      
      print(h1 + h2 + h3)
    }
    dev.off()
  }
}

system("aws s3 sync data/KINOMO/ s3://melanoma-ribas/KINOMO/ --exclude '*' --include '*_top100_pathways.pdf' --exclude '*._*' ")


# Run pw analysis on metaprograms
markers<-read.csv('data/NMF/gsea.metagenes.csv',row.names = 1)
markers<-read.csv('data/NMF/gsea.metagenes.0.3.csv',row.names = 1)

genesets_w <- as.list(as.data.frame(read.csv('~/brain_mets/signatures/wouters_mel_sigs.csv', na.strings = c('', NA))))
genesets_emt <- as.list(as.data.frame(read.csv('~/brain_mets/signatures/emt_sigs.csv', na.strings = c('', NA))))


pdf(paste0(directory,'plots_',label,'_pathway.pdf'),width = 10,height = 7)
for(i in 1:ncol(markers)){
  subs<-markers[,i]
  print(hyp_dots(hypeR(as.character(subs), genesets=genesets_w),
                 title = paste0('Wouters pathways; FB',i),abrv=80)+ theme_bw() + 
          theme(axis.title.y = element_blank(), 
                axis.text.y = element_text(color = 'black'), 
                axis.text.x = element_text(color = 'black'), 
                panel.grid.minor = element_blank()) + 
          scale_color_gradient(high = 'black', low = '#C51B8A'))
  
  print(hyp_dots(hypeR(as.character(subs), genesets=genesets_emt),
                 title = paste0('EMT pathways; FB',i),abrv=80)+ theme_bw() + 
          theme(axis.title.y = element_blank(), 
                axis.text.y = element_text(color = 'black'), 
                axis.text.x = element_text(color = 'black'), 
                panel.grid.minor = element_blank()) + 
          scale_color_gradient(high = 'black', low = '#C51B8A'))
}
dev.off()

