#!/usr/bin/env Rscript

#### CNA frequency plots generated from inferCNV, ichorCNA and numbat outputs
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
library(patchwork)
'%notin%' <- Negate('%in%')

setwd('~/Documents/melanoma/')

pats<-c('F01_on','F01_pre','F02_on','F02_pre','F03_post1_on2','F03_post1_pre2','F04_pre',
        'F05_pre','F06_post1_pre2','F07_post1_pre2','F08_post','F09_post',
        'F12_post','F12_pre','F15_pre','F16_post1_pre2','F16_pre','F17_post',
        'F18_post','F20_post1_pre2','F22_post','F23_post',
        'F26_post','F27_post1_pre2','F28_post1_pre2','F29_post1_pre2','F30_post',
        'F31_post','R204_pre','R294_on','R308_pre','R310_on1','R310_on2','R310_pre',
        'R319_on','R319_pre','R328_on','R329_on','R334_pre','R354_pre')

pats_wgs<-c('F01_on','F01_pre','F02_on','F02_pre','F03_post1_on2','F03_post1_pre2','F04_pre',
            'F05_pre','F06_post1_pre2','F07_post1_pre2','F08_post','F09_post',
            'F12_post','F12_pre','F15_pre','F16_post1_pre2','F16_pre','F17_post',
            'F18_post','F20_post1_pre2','F22_post','F23_post',
            'F26_post','F27_post1_pre2','F28_post1_pre2','F29_post1_pre2','F30_post',
            'F31_post', 'R294_on', 'R310_pre', 'R310_on2', 'R319_pre', 'R319_on', 'R328_on',
            'R329_on', 'R334_pre', 'R354_pre')

pats_wgs_treated<-c('F01_on','F02_on','F03_post1_on2','F03_post1_pre2',
                    'F06_post1_pre2','F07_post1_pre2','F08_post','F09_post',
                    'F12_post','F17_post','F18_post','F20_post1_pre2','F22_post','F23_post',
                    'F26_post','F27_post1_pre2','F28_post1_pre2','F29_post1_pre2','F30_post',
                    'F31_post', 'R294_on', 'R310_on2', 'R319_pre', 'R319_on', 'R328_on',
                    'R329_on', 'R334_pre', 'R354_pre')

pats_wgs_untreated<-c('F01_pre','F02_pre','F04_pre','F05_pre',
                      'F12_pre','F15_pre','F16_post1_pre2','F16_pre',
                      'R310_pre')

pats_numbat_treated<-c('F01_on','F02_on','F03_post1_on2','F03_post1_pre2',
                    'F06_post1_pre2','F07_post1_pre2','F08_post','F09_post',
                    'F12_post','F17_post','F18_post','F20_post1_pre2','F22_post','F23_post',
                    'F26_post','F27_post1_pre2','F28_post1_pre2','F29_post1_pre2','F30_post',
                    'F31_post', 'R294_on', 'R310_on2', 'R319_pre', 'R319_on', 'R328_on',
                    'R329_on', 'R334_pre', 'R354_pre','R204_pre','R310_on1')

pats_numbat_untreated<-c('F01_pre','F02_pre','F04_pre','F05_pre',
                      'F12_pre','F15_pre','F16_post1_pre2','F16_pre',
                      'R310_pre','R308_pre')


#### inferCNV ######

# set up df for collecting all samples
all_samples<- read.csv(paste0("data/",pats[1],"/table_",pats[1],"_inferCNV_geneloc.csv"))
all_samples<-data.frame(chr_pos=paste0(all_samples$chrom,'_',all_samples$start.pos),
                    median=(all_samples$median - 1))
colnames(all_samples)[2]<-pats[1]

for(pat in pats[2:length(pats)]){
  print(paste('Processing:', pat))
  
  system(paste0("aws s3 sync s3://melanoma-ribas/inferCNV/inferCNV_subcluster_", pat, 
                "/ data/",pat,"/ --exclude '*' --include '*_inferCNV_geneloc.csv' "))
  
  if(!file.exists(paste0("data/",pat,"/table_",pat,"_inferCNV_geneloc.csv"))){
    # Get inferCNV gene table
    system("aws s3 sync s3://melanoma-ribas/Seurat/annotation/main/ data/annotation/main/ --exclude '*' --include 'data_nontum_merged_v2_barcodes.csv' ")
    nontum_cells<-read.csv('data/annotation/main/data_nontum_merged_v2_barcodes.csv')
    nontum_cells$barcode_orig<-substr(nontum_cells$x,1,18)
    nontum_cells$sample<-substr(nontum_cells$x,20,50)
    
    for(pat in pats){
      print(paste0('Starting ',pat))
      system(paste0("aws s3 sync s3://melanoma-ribas/inferCNV/inferCNV_subcluster_", pat, 
                    "/ data/",pat,"/ --exclude '*' --include 'run.final.infercnv_obj' "))
      infercnv_obj <- readRDS(paste0("data/",pat,"/run.final.infercnv_obj"))
      
      specific_tumor_cells = colnames(infercnv_obj@expr.data)[colnames(infercnv_obj@expr.data) %notin% 
                                                                nontum_cells[nontum_cells$sample==pat,3]]
      all(specific_tumor_cells %in% colnames(infercnv_obj@expr.data))
      
      if(length(specific_tumor_cells)>0){
        ischrarr = grep("chr",infercnv_obj@gene_order$chr)
        rowmediansarr = rowMedians(infercnv_obj@expr.data[ischrarr, specific_tumor_cells])
        atable = data.frame(sampleID = rep(pat, length(rowmediansarr)), 
                            chrom = substring(infercnv_obj@gene_order$chr[ischrarr], 4), 
                            start.pos = infercnv_obj@gene_order$start[ischrarr], 
                            end.pos = infercnv_obj@gene_order$stop[ischrarr], 
                            median = rowmediansarr)
        write.csv(atable, paste0("data/",pat,"/table_",pat,"_inferCNV_geneloc.csv"),row.names = F)
      }
      system(paste0("aws s3 sync data/",pat,"/ s3://melanoma-ribas/inferCNV/inferCNV_subcluster_", pat, 
                    "/  --exclude '*' --include '*_inferCNV_geneloc.csv' "))
      system(paste0("rm data/",pat,"/run.final.infercnv_obj"))
    }
  }
  
  loc_tmp<- read.csv(paste0("data/",pat,"/table_",pat,"_inferCNV_geneloc.csv"))
  loc_tmp<-data.frame(chr_pos=paste0(loc_tmp$chrom,'_',loc_tmp$start.pos),
                      median=(loc_tmp$median - 1))
  colnames(loc_tmp)[2]<-pat
  all_samples<-full_join(all_samples,loc_tmp,by = "chr_pos")
  
}

all_samples$chr<-as.numeric(sapply(strsplit(all_samples$chr_pos,"_"), `[`, 1))
all_samples$pos<-as.numeric(sapply(strsplit(all_samples$chr_pos,"_"), `[`, 2))
all_samples<-all_samples[,c('chr','pos',names(all_samples)[2:(ncol(all_samples)-2)])]

table(all_samples$chr)

only12<-all_samples[all_samples$chr == 12,]
reduced_samples<-all_samples[all_samples$chr != 12,]
reduced_samples<-rbind.data.frame(reduced_samples,
                                  only12[round(runif(1000,min=1,max=8389472)),])

infer_win <- winsorize(data=reduced_samples,pos.unit = "bp",verbose=T,method = 'pcf')
infer_seg<-pcf(infer_win,Y = reduced_samples,verbose = T)
plotGamma(infer_win,sample = 1,chrom = 1,gammaRange = c(1,1000))
#print(plotGenome(data = NULL, segments = infer_seg, pos.unit = "bp", sample = NULL))
infer_seg$mean[infer_seg$mean=='NaN']<-0
plotFreq(segments=infer_seg,thres.gain=0.02)
plotAberration(segments=infer_seg,thres.gain=0.02)
plotHeatmap(segments=infer_seg,upper.lim=0.02)

pdf(paste0("data/tumor/plots_melanoma_CNA_frequency_inferCNV.pdf"),width = 16)
hist(infer_seg$mean,breaks = 200)
hist(infer_seg$mean,breaks = 200,ylim = c(0,50))
plotFreq(segments=infer_seg,thres.gain=0.02)
plotFreq(segments=infer_seg,thres.gain=0.01)
plotFreq(segments=infer_seg,thres.gain=0.015)
plotFreq(segments=infer_seg,thres.gain=0.03)
plotFreq(segments=infer_seg,thres.gain=0.06)
plotAberration(segments=infer_seg,thres.gain=0.02)
plotAberration(segments=infer_seg,thres.gain=0.015)
plotHeatmap(segments=infer_seg,upper.lim=0.02)
plotHeatmap(segments=infer_seg,upper.lim=0.015)
dev.off()


#### Numbat all ######

# set up df for collecting all samples
num_seg<-NULL
for(pat in pats){
  print(paste('Processing:', pat))
  
  system(paste0("aws s3 sync s3://melanoma-ribas/numbat/",pat,"/ data/",pat,"/ --exclude '*' --include '*table_*_segs.csv' --include '*table_*_clone.csv' "))
  num_seg_tmp<-read.csv(paste0('data/',pat,'/table_',pat,'_segs.csv'))
  num_seg_tmp<-data.frame(sampleID=pat,
                          chrom=num_seg_tmp$CHROM,
                          arm=NA,
                          start.pos=num_seg_tmp$seg_start,
                          end.pos=num_seg_tmp$seg_end,
                          n.probes=NA,
                          mean=num_seg_tmp$log2FC_median)

  num_seg<-rbind.data.frame(num_seg,num_seg_tmp)
}

table(num_seg$chrom)
num_seg<-num_seg[which(!is.na(num_seg$mean)),]
hist(num_seg$mean,breaks = 200)

pdf(paste0("data/tumor/plots_melanoma_CNA_frequency_Numbat.pdf"),width = 16)
hist(num_seg$mean,breaks = 200)
plotFreq(segments=num_seg,thres.gain=0.4)
plotFreq(segments=num_seg,thres.gain=0.5,thres.loss = -0.4)
plotFreq(segments=num_seg,thres.gain=0.5,thres.loss = -0.3)
plotFreq(segments=num_seg,thres.gain=0.3)
plotFreq(segments=num_seg,thres.gain=0.5)
plotFreq(segments=num_seg,thres.gain=0.7)
plotAberration(segments=num_seg,thres.gain=0.5,thres.loss = -0.4)
plotHeatmap(segments=num_seg,upper.lim=0.5,lower.lim = -0.4)
plotHeatmap(segments=num_seg,upper.lim=0.4)
dev.off()




#### Numbat treated ######

# set up df for collecting all samples
num_seg<-NULL
for(pat in pats_numbat_treated){
  print(paste('Processing:', pat))
  
  system(paste0("aws s3 sync s3://melanoma-ribas/numbat/",pat,"/ data/",pat,"/ --exclude '*' --include '*table_*_segs.csv' --include '*table_*_clone.csv' "))
  num_seg_tmp<-read.csv(paste0('data/',pat,'/table_',pat,'_segs.csv'))
  num_seg_tmp<-data.frame(sampleID=pat,
                          chrom=num_seg_tmp$CHROM,
                          arm=NA,
                          start.pos=num_seg_tmp$seg_start,
                          end.pos=num_seg_tmp$seg_end,
                          n.probes=NA,
                          mean=num_seg_tmp$log2FC_median)
  
  num_seg<-rbind.data.frame(num_seg,num_seg_tmp)
}

table(num_seg$chrom)
num_seg<-num_seg[which(!is.na(num_seg$mean)),]
hist(num_seg$mean,breaks = 200)

pdf(paste0("data/tumor/plots_melanoma_CNA_frequency_Numbat_treated.pdf"),width = 16)
hist(num_seg$mean,breaks = 200)
plotFreq(segments=num_seg,thres.gain=0.4)
plotFreq(segments=num_seg,thres.gain=0.5,thres.loss = -0.4)
plotFreq(segments=num_seg,thres.gain=0.5,thres.loss = -0.3)
plotFreq(segments=num_seg,thres.gain=0.3)
plotFreq(segments=num_seg,thres.gain=0.5)
plotFreq(segments=num_seg,thres.gain=0.7)
dev.off()


dat<-read.csv('data/tumor/table_melanoma_chr_arm_corr.csv')
summary(dat)



#### Numbat untreated ######

# set up df for collecting all samples
num_seg<-NULL
for(pat in pats_numbat_untreated){
  print(paste('Processing:', pat))
  
  system(paste0("aws s3 sync s3://melanoma-ribas/numbat/",pat,"/ data/",pat,"/ --exclude '*' --include '*table_*_segs.csv' --include '*table_*_clone.csv' "))
  num_seg_tmp<-read.csv(paste0('data/',pat,'/table_',pat,'_segs.csv'))
  num_seg_tmp<-data.frame(sampleID=pat,
                          chrom=num_seg_tmp$CHROM,
                          arm=NA,
                          start.pos=num_seg_tmp$seg_start,
                          end.pos=num_seg_tmp$seg_end,
                          n.probes=NA,
                          mean=num_seg_tmp$log2FC_median)
  
  num_seg<-rbind.data.frame(num_seg,num_seg_tmp)
}

table(num_seg$chrom)
num_seg<-num_seg[which(!is.na(num_seg$mean)),]
hist(num_seg$mean,breaks = 200)

pdf(paste0("data/tumor/plots_melanoma_CNA_frequency_Numbat_untreated.pdf"),width = 16)
hist(num_seg$mean,breaks = 200)
plotFreq(segments=num_seg,thres.gain=0.4)
plotFreq(segments=num_seg,thres.gain=0.5,thres.loss = -0.4)
plotFreq(segments=num_seg,thres.gain=0.5,thres.loss = -0.3)
plotFreq(segments=num_seg,thres.gain=0.3)
plotFreq(segments=num_seg,thres.gain=0.5)
plotFreq(segments=num_seg,thres.gain=0.7)
dev.off()



#### ichorCNA all ######

# set up df for collecting all samples
wgs_seg<-NULL
for(pat in pats_wgs){
  print(paste('Processing:', pat))
  
  system(paste0("aws s3 sync s3://melanoma-ribas/ichorCNA/",pat,"/ data/",pat,"/ichorCNA/ --exclude '*' --include '*.seg' "))
  wgs_seg_tmp<-read.table(paste0('data/',pat,'/ichorCNA/',
                            grep(pattern = 'tumor.seg$',
                                 x = list.files(paste0('data/',pat,'/ichorCNA/')),
                                 value = T)),
                     header = T)
  wgs_seg_tmp<-data.frame(sampleID=pat,
                          chrom=wgs_seg_tmp$chr,
                          arm=NA,
                          start.pos=wgs_seg_tmp$start,
                          end.pos=wgs_seg_tmp$end,
                          n.probes=NA,
                          mean=wgs_seg_tmp$median)
  
  wgs_seg<-rbind.data.frame(wgs_seg,wgs_seg_tmp)
}

table(wgs_seg$chrom)
wgs_seg<-wgs_seg[which(!is.na(wgs_seg$mean)),]
hist(wgs_seg$mean,breaks = 200)

pdf(paste0("data/tumor/plots_melanoma_CNA_frequency_ichorCNA.pdf"),width = 16)
hist(wgs_seg$mean,breaks = 200)
plotFreq(segments=wgs_seg,thres.gain=0.1)
plotFreq(segments=wgs_seg,thres.gain=0.2)
plotFreq(segments=wgs_seg,thres.gain=0.3)
plotFreq(segments=wgs_seg,thres.gain=0.4)
plotFreq(segments=wgs_seg,thres.gain=0.5)
plotFreq(segments=wgs_seg,thres.gain=0.7)
dev.off()


#### ichorCNA treated ######

# set up df for collecting all samples
wgs_seg<-NULL
for(pat in pats_wgs_treated){
  print(paste('Processing:', pat))
  
  system(paste0("aws s3 sync s3://melanoma-ribas/ichorCNA/",pat,"/ data/",pat,"/ichorCNA/ --exclude '*' --include '*.seg' "))
  wgs_seg_tmp<-read.table(paste0('data/',pat,'/ichorCNA/',
                                 grep(pattern = 'tumor.seg$',
                                      x = list.files(paste0('data/',pat,'/ichorCNA/')),
                                      value = T)),
                          header = T)
  wgs_seg_tmp<-data.frame(sampleID=pat,
                          chrom=wgs_seg_tmp$chr,
                          arm=NA,
                          start.pos=wgs_seg_tmp$start,
                          end.pos=wgs_seg_tmp$end,
                          n.probes=NA,
                          mean=wgs_seg_tmp$median)
  
  wgs_seg<-rbind.data.frame(wgs_seg,wgs_seg_tmp)
}

table(wgs_seg$chrom)
wgs_seg<-wgs_seg[which(!is.na(wgs_seg$mean)),]
hist(wgs_seg$mean,breaks = 200)

pdf(paste0("data/tumor/plots_melanoma_CNA_frequency_ichorCNA_treated.pdf"),width = 16)
hist(wgs_seg$mean,breaks = 200)
plotFreq(segments=wgs_seg,thres.gain=0.1)
plotFreq(segments=wgs_seg,thres.gain=0.2)
plotFreq(segments=wgs_seg,thres.gain=0.3)
plotFreq(segments=wgs_seg,thres.gain=0.4)
plotFreq(segments=wgs_seg,thres.gain=0.5)
plotFreq(segments=wgs_seg,thres.gain=0.7)
dev.off()


#### ichorCNA untreated ######

# set up df for collecting all samples
wgs_seg<-NULL
for(pat in pats_wgs_untreated){
  print(paste('Processing:', pat))
  
  system(paste0("aws s3 sync s3://melanoma-ribas/ichorCNA/",pat,"/ data/",pat,"/ichorCNA/ --exclude '*' --include '*.seg' "))
  wgs_seg_tmp<-read.table(paste0('data/',pat,'/ichorCNA/',
                                 grep(pattern = 'tumor.seg$',
                                      x = list.files(paste0('data/',pat,'/ichorCNA/')),
                                      value = T)),
                          header = T)
  wgs_seg_tmp<-data.frame(sampleID=pat,
                          chrom=wgs_seg_tmp$chr,
                          arm=NA,
                          start.pos=wgs_seg_tmp$start,
                          end.pos=wgs_seg_tmp$end,
                          n.probes=NA,
                          mean=wgs_seg_tmp$median)
  
  wgs_seg<-rbind.data.frame(wgs_seg,wgs_seg_tmp)
}

table(wgs_seg$chrom)
wgs_seg<-wgs_seg[which(!is.na(wgs_seg$mean)),]
hist(wgs_seg$mean,breaks = 200)

pdf(paste0("data/tumor/plots_melanoma_CNA_frequency_ichorCNA_untreated.pdf"),width = 16)
hist(wgs_seg$mean,breaks = 200)
plotFreq(segments=wgs_seg,thres.gain=0.1)
plotFreq(segments=wgs_seg,thres.gain=0.2)
plotFreq(segments=wgs_seg,thres.gain=0.3)
plotFreq(segments=wgs_seg,thres.gain=0.4)
plotFreq(segments=wgs_seg,thres.gain=0.5)
plotFreq(segments=wgs_seg,thres.gain=0.7)
dev.off()

