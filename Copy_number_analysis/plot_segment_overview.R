#!/usr/bin/Rscript

#### Title: Segment overview plot for numbat ichorCNA inferCNV
#### Author: Jana Biermann, PhD

library(dplyr)
library(ggplot2)
library(ggrepel)
library(scales)
library(Seurat)
library(reshape2)
library(grid)
library(rlist)
library(matrixStats)
library(copynumber)
'%notin%' <- Negate('%in%')

setwd('~/Documents/melanoma/')


pats_wgs<-c('F01_on','F01_pre','F02_on','F02_pre','F03_post1_on2','F03_post1_pre2','F04_pre',
            'F05_pre','F06_post1_pre2','F07_post1_pre2','F08_post','F09_post',
            'F12_post','F12_pre','F15_pre','F16_post1_pre2','F16_pre','F17_post',
            'F18_post','F20_post1_pre2','F22_post','F23_post',
            'F26_post','F27_post1_pre2','F28_post1_pre2','F29_post1_pre2','F30_post',
            'F31_post', 'R294_on', 'R310_pre', 'R310_on2', 'R319_pre', 'R319_on', 'R328_on',
            'R329_on', 'R334_pre', 'R354_pre')

pats<-c('F01_on','F01_pre','F02_on','F02_pre','F03_post1_on2','F03_post1_pre2','F04_pre',
        'F05_pre','F06_post1_pre2','F07_post1_pre2','F08_post','F09_post',
        'F12_post','F12_pre','F15_pre','F16_post1_pre2','F16_pre','F17_post',
        'F18_post','F20_post1_pre2','F22_post','F23_post',
        'F26_post','F27_post1_pre2','F28_post1_pre2','F29_post1_pre2','F30_post',
        'F31_post','R204_pre','R294_on','R308_pre','R310_on1','R310_on2','R310_pre',
        'R319_on','R319_pre','R328_on','R329_on','R334_pre','R354_pre')


#### Summary of inferCNV CNAs #####
chr_len<-read.csv('~/brain_mets/misc/GRCh38.p14_chr_length.csv')
chr_len<-chr_len[1:22,]
chr_len$Chromosome<-as.numeric(chr_len$Chromosome)
chr_len$Total_length_bp<-as.numeric(gsub(",", "", chr_len$Total_length_bp))
chr_len$Continuous_length_bp<-apply(chr_len,1,FUN = function(x){sum(chr_len[1:x,2])})

segs_comb<-NULL
for(pat in pats){
  print(paste0('Processing: ', pat))
  
  #system(paste0("aws s3 sync s3://melanoma-ribas/inferCNV/inferCNV_subcluster_", pat, "/ data/",pat,"/ --exclude '*' --include '*_inferCNV_geneloc.csv' "))
  df<- read.csv(paste0("data/",pat,"/table_",pat,"_inferCNV_geneloc.csv"))
  
  tmp<-data.frame(Sample=pat,
                  Chromosome=as.numeric(df$chrom),
                  Start=df$start.pos,
                  End=df$end.pos,
                  Num_Probes=NA,
                  Segment_Mean=log2(df$median))
  
  segs_comb<-rbind.data.frame(segs_comb,tmp)
}
write.csv(segs_comb,'data/tumor/table_cna_inferCNV_combined.csv',row.names = F)

segs_comb<-read.csv('data/tumor/table_cna_inferCNV_combined.csv')
df<-segs_comb

# Modify data to add xmin and xmax for each sample
sample_positions <- data.frame(Sample=unique(df$Sample),
                               xmin=seq(0.6, length(unique(df$Sample)) , by = 1),
                               xmax=seq(1.4, length(unique(df$Sample))+1 , by = 1))
df <- left_join(df,sample_positions,by='Sample')

# Get total bp position
for(chr in 1:22){
  tmp_segs<-df[df$Chromosome==chr,]
  if(chr==1){
    df$seg_start_continuous[df$Chromosome==1]<-tmp_segs$Start
    df$seg_end_continuous[df$Chromosome==1]<-tmp_segs$End
  }else{
    df$seg_start_continuous[df$Chromosome==chr]<-tmp_segs$Start+sum(chr_len[1:(chr-1),2])
    df$seg_end_continuous[df$Chromosome==chr]<-tmp_segs$End+sum(chr_len[1:(chr-1),2])
  }
}

hist(df$Segment_Mean,breaks = 500,ylim = c(0,1000))

df_clipped<-df
df_clipped$Segment_Mean[df_clipped$Segment_Mean>0.3]<-0.3
df_clipped$Segment_Mean[df_clipped$Segment_Mean< -0.3]<- -0.3

pdf(paste0("data/tumor/plots_melanoma_cna_overview_inferCNV.pdf"),width = 10)
print(ggplot(df_clipped,aes(x=Sample,y = seg_start_continuous))+
        geom_point(size=NA)+
        theme_bw()+
        xlab('Sample')+
        ylab('Genomic position (Chromosome)')+
        scale_y_continuous(breaks=chr_len$Continuous_length_bp, labels=chr_len$Chromosome,
                           expand = c(0,0))+
        geom_rect(aes(ymin = seg_start_continuous,
                      ymax = seg_end_continuous, 
                      xmin = xmin, 
                      xmax = xmax,
                      fill = Segment_Mean), 
                  lty=0)+
        scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
        theme(axis.text.x = element_text(size=5,hjust = 1),
              legend.position = 'bottom',
              panel.grid.major.x = element_blank(),
              panel.grid.major.y = element_blank(),
              panel.grid.minor.x = element_blank())+
        geom_hline(yintercept = chr_len$Continuous_length_bp, size=0.2,col='black',linetype=1)+
        coord_flip()+
        ggtitle(paste0('CNAs from inferCNV (all samples)')))
dev.off()



#### Summary of WGS CNAs #####
chr_len<-read.csv('~/brain_mets/misc/GRCh38.p14_chr_length.csv')
chr_len<-chr_len[1:22,]
chr_len$Chromosome<-as.numeric(chr_len$Chromosome)
chr_len$Total_length_bp<-as.numeric(gsub(",", "", chr_len$Total_length_bp))
chr_len$Continuous_length_bp<-apply(chr_len,1,FUN = function(x){sum(chr_len[1:x,2])})

# ichorCNA segs
segs_comb<-NULL
for(pat in pats_wgs){
  print(paste0('Processing: ', pat))
  
  #system(paste0("aws s3 sync s3://melanoma-ribas/ichorCNA/",pat,"/ data/",pat,"/ichorCNA/ --exclude '*' --include '*.seg' "))
  wgs<-read.table(paste0('data/',pat,'/ichorCNA/',
                         grep(pattern = 'tumor.cna.seg$',
                              x = list.files(paste0('data/',pat,'/ichorCNA/')),
                              value = T)),
                  header = T)
  
  wgs<-na.omit(wgs)
  wgs<-wgs[wgs$chr %notin% c('X','Y'),]

  tmp<-data.frame(Sample=pat,
                  Chromosome=as.numeric(wgs$chr),
                  Start=wgs$start,
                  End=wgs$end,
                  Segment_Mean=wgs[,6])
  
  segs_comb<-rbind.data.frame(segs_comb,tmp)
}

write.csv(segs_comb,'data/tumor/table_cna_WGS_combined.csv',row.names = F)


segs_comb<-read.csv('data/tumor/table_cna_WGS_combined.csv')
df<-segs_comb

# Modify data to add xmin and xmax for each sample
sample_positions <- data.frame(Sample=unique(df$Sample),
                               xmin=seq(0.6, length(unique(df$Sample)) , by = 1),
                               xmax=seq(1.4, length(unique(df$Sample))+1 , by = 1))
df <- left_join(df,sample_positions,by='Sample')

# Get total bp position
for(chr in 1:22){
  tmp_segs<-df[df$Chromosome==chr,]
  if(chr==1){
    df$seg_start_continuous[df$Chromosome==1]<-tmp_segs$Start
    df$seg_end_continuous[df$Chromosome==1]<-tmp_segs$End
  }else{
    df$seg_start_continuous[df$Chromosome==chr]<-tmp_segs$Start+sum(chr_len[1:(chr-1),2])
    df$seg_end_continuous[df$Chromosome==chr]<-tmp_segs$End+sum(chr_len[1:(chr-1),2])
  }
}

hist(df$Segment_Mean,breaks = 500,ylim = c(0,1000),xlim = c(-5,5))

df_clipped<-df
df_clipped$Segment_Mean[df_clipped$Segment_Mean>1.5]<-1.5
df_clipped$Segment_Mean[df_clipped$Segment_Mean< -1.5]<- -1.5

pdf(paste0("data/tumor/plots_melanoma_cna_overview_WGS.pdf"),width = 10)
print(ggplot(df_clipped,aes(x=Sample,y = seg_start_continuous))+
        geom_point(size=NA)+
        theme_bw()+
        xlab('Sample')+
        ylab('Genomic position (Chromosome)')+
        scale_y_continuous(breaks=chr_len$Continuous_length_bp, labels=chr_len$Chromosome,
                           expand = c(0,0))+
        geom_rect(aes(ymin = seg_start_continuous,
                      ymax = seg_end_continuous, 
                      xmin = xmin, 
                      xmax = xmax,
                      fill = Segment_Mean), 
                  lty=0)+
        scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
        theme(axis.text.x = element_text(size=5,hjust = 1),
              legend.position = 'bottom',
              panel.grid.major.x = element_blank(),
              panel.grid.major.y = element_blank(),
              panel.grid.minor.x = element_blank())+
        geom_hline(yintercept = chr_len$Continuous_length_bp, size=0.2,col='black',linetype=1)+
        coord_flip()+
        ggtitle(paste0('CNAs from WGS (all samples)')))
dev.off()



#### Summary of numbat CNAs #####
chr_len<-read.csv('~/brain_mets/misc/GRCh38.p14_chr_length.csv')
chr_len<-chr_len[1:22,]
chr_len$Chromosome<-as.numeric(chr_len$Chromosome)
chr_len$Total_length_bp<-as.numeric(gsub(",", "", chr_len$Total_length_bp))
chr_len$Continuous_length_bp<-apply(chr_len,1,FUN = function(x){sum(chr_len[1:x,2])})

# numbat segs
segs_comb<-NULL
for(pat in pats){
  print(paste('Processing:', pat))
  
  # Read-in numbat segs
  segs<-read.csv(paste0('data/',pat,'/table_',pat,'_segs.csv'))
  segs<-segs[segs$cnv_states!='neu',]
  
  clone<-read.csv(paste0('data/',pat,'/table_',pat,'_clone.csv'))
  clone$cell<-NULL
  clone<-clone[!duplicated(clone),]
  clone$sample<-paste0(pat,"_clone",clone$clone_opt)
  clone<-clone[clone$GT_opt != '',]
  
  # Add bulk
  tmp<-data.frame(Sample=paste0(pat,'_numbat_bulk'),
                  Chromosome=as.numeric(segs$CHROM),
                  Start=segs$seg_start,
                  End=segs$seg_end,
                  Segment_Mean=segs$log2FC_mean)
  segs_comb<-rbind.data.frame(segs_comb,tmp)
  
  # Add individual clones
  for(gt in clone$GT_opt){
    cl<-clone[clone$GT_opt==gt,3]
    for(CNA in strsplit(gt,',')[[1]]){
      if(grepl('_',CNA)==T){
        CNA_cut<-strsplit(CNA,'_')[[1]][1]
        tmp<-data.frame(Sample=cl,
                        Chromosome=as.numeric(segs[segs$seg_cons==CNA_cut,3]),
                        Start=segs[segs$seg_cons==CNA_cut,4],
                        End=segs[segs$seg_cons==CNA_cut,5],
                        Segment_Mean=segs[segs$seg_cons==CNA_cut,9])
        segs_comb<-rbind.data.frame(segs_comb,tmp)
      }else{
        tmp<-data.frame(Sample=cl,
                        Chromosome=as.numeric(segs[segs$seg_cons==CNA,3]),
                        Start=segs[segs$seg_cons==CNA,4],
                        End=segs[segs$seg_cons==CNA,5],
                        Segment_Mean=segs[segs$seg_cons==CNA,9])
        segs_comb<-rbind.data.frame(segs_comb,tmp)
      }
    }
  }
}

write.csv(segs_comb,'data/tumor/table_cna_numbat_combined.csv',row.names = F)


segs_comb<-read.csv('data/tumor/table_cna_numbat_combined.csv')
df<-segs_comb %>%
  filter(grepl("bulk", Sample))

# Modify data to add xmin and xmax for each sample
sample_positions <- data.frame(Sample=unique(df$Sample),
                               xmin=seq(0.6, length(unique(df$Sample)) , by = 1),
                               xmax=seq(1.4, length(unique(df$Sample))+1 , by = 1))
df <- left_join(df,sample_positions,by='Sample')

# Get total bp position
for(chr in 1:22){
  tmp_segs<-df[df$Chromosome==chr,]
  if(chr==1){
    df$seg_start_continuous[df$Chromosome==1]<-tmp_segs$Start
    df$seg_end_continuous[df$Chromosome==1]<-tmp_segs$End
  }else{
    df$seg_start_continuous[df$Chromosome==chr]<-tmp_segs$Start+sum(chr_len[1:(chr-1),2])
    df$seg_end_continuous[df$Chromosome==chr]<-tmp_segs$End+sum(chr_len[1:(chr-1),2])
  }
}

summary(df)
hist(df$Segment_Mean,breaks = 500,ylim = c(0,50))

df_clipped<-df
df_clipped$Segment_Mean[df_clipped$Segment_Mean>1.5]<-1.5
df_clipped$Segment_Mean[df_clipped$Segment_Mean< -1.5]<- -1.5

pdf(paste0("data/tumor/plots_melanoma_cna_overview_numbat.pdf"),width = 10)
print(ggplot(df_clipped,aes(x=Sample,y = seg_start_continuous))+
        geom_point(size=NA)+
        theme_bw()+
        xlab('Sample')+
        ylab('Genomic position (Chromosome)')+
        scale_y_continuous(breaks=chr_len$Continuous_length_bp, labels=chr_len$Chromosome,
                           expand = c(0,0))+
        geom_rect(aes(ymin = seg_start_continuous,
                      ymax = seg_end_continuous, 
                      xmin = xmin, 
                      xmax = xmax,
                      fill = Segment_Mean), 
                  lty=0)+
        scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
        theme(axis.text.x = element_text(size=5,hjust = 1),
              legend.position = 'bottom',
              panel.grid.major.x = element_blank(),
              panel.grid.major.y = element_blank(),
              panel.grid.minor.x = element_blank())+
        geom_hline(yintercept = chr_len$Continuous_length_bp, size=0.2,col='black',linetype=1)+
        coord_flip()+
        ggtitle(paste0('CNAs from numbat (all samples)')))
dev.off()

