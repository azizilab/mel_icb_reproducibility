#!/usr/bin/Rscript

#### Title: Chromosome arm comparison (numbat ichorCNA inferCNV)
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

#setwd('~/Documents/melanoma/')

pat<-'R294_on'

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


#### functions ####
# Annotate chr arm and split segments spanning both arms
annotate_chr_arms<-function(x){
  chr_arms_anno<-NULL
  for(chr in 1:22){
    df_tmp<-x[x$chr==chr,]
    df_tmp$arm<-NA
    
    df_tmp_cut<-NULL
    for(segm in 1:nrow(df_tmp)){
      if( df_tmp[segm,3] <= cyto[cyto$chrom==chr & cyto$arm == 'p',4]$arm_end){
        # entire segment on p
        df_tmp$arm[segm]<-'p'
        df_tmp_cut<-rbind.data.frame(df_tmp_cut,df_tmp[segm,])
        
      }else if(df_tmp[segm,2] > cyto[cyto$chrom==chr & cyto$arm == 'q',3]$arm_start){
        # entire segment on q
        df_tmp$arm[segm]<-'q'
        df_tmp_cut<-rbind.data.frame(df_tmp_cut,df_tmp[segm,])
        
      }else{
        #split p
        df_tmp_cut<-rbind.data.frame(df_tmp_cut,
                                         data.frame(chr=chr,
                                                    seg_start=df_tmp$seg_start[segm],
                                                    seg_end=cyto[cyto$chrom==chr & cyto$arm == 'p',4]$arm_end,
                                                    median=df_tmp$median[segm],
                                                    seg_length= cyto[cyto$chrom==chr & cyto$arm == 'p',4]$arm_end - df_tmp$seg_start[segm],
                                                    arm='p'))
        # split q
        df_tmp_cut<-rbind.data.frame(df_tmp_cut,
                                         data.frame(chr=chr,
                                                    seg_start=cyto[cyto$chrom==chr & cyto$arm == 'p',4]$arm_end +1,
                                                    seg_end=df_tmp$seg_end[segm],
                                                    median=df_tmp$median[segm],
                                                    seg_length= df_tmp$seg_end[segm] - cyto[cyto$chrom==chr & cyto$arm == 'p',4]$arm_end +1,
                                                    arm='q'))
      }
    }
    chr_arms_anno<-rbind.data.frame(chr_arms_anno,df_tmp_cut)
  }
  return(chr_arms_anno)
}


#### Code ####
# Read-in cytoband info
cyto<-read.delim('~/brain_mets/misc/GRCh38_cytoBand.txt',header = F)
cyto<- cyto %>% filter(nchar(V1)<=5)
cyto$arm<-substr(cyto$V4,1,1)
cyto <- cyto %>% group_by(V1,arm) %>% summarise(arm_start=min(V2),arm_end=max(V3) )
colnames(cyto)[1]<-'chrom_name'
cyto <- cyto[cyto$chrom_name %notin% c('chrX','chrY','chrM'),]
cyto$chrom<-as.numeric(substr(cyto$chrom_name,4,7))


pdf(paste0("data/tumor/plots_melanoma_chr_arm_corr_v2.pdf"),width = 16)
collect_corrs<-NULL
# Loop through pats (calculate correlations and plot)
for(pat in pats){
  print(paste('Processing:', pat))
  wgs<-FALSE
  inf<-FALSE
  p_num_inf<-NULL
  p_num_wgs<-NULL
  p_inf_wgs<-NULL
  
  #### numbat ####
  # Read-in numbat segs
  system(paste0("aws s3 sync s3://melanoma-ribas/numbat/",pat,"/ data/",pat,"/ --exclude '*' --include '*table_*_segs.csv' --include '*table_*_clone.csv' "))
  num<-read.csv(paste0('data/',pat,'/table_',pat,'_segs.csv'))
  df_num<-data.frame(chr=num$CHROM,
                     seg_start=num$seg_start,
                     seg_end=num$seg_end,
                     seg_length=num$seg_end - num$seg_start,
                     median=num$log2FC_median)
  df_num<-na.omit(df_num)
  
  # Annotate chr arm and split segments spanning both arms
  chr_arms_num<-annotate_chr_arms(df_num)
  chr_arms_num<-chr_arms_num %>% 
    mutate(chr_arm=paste0(chr,arm)) %>% 
    group_by(chr_arm) %>% 
    mutate(arm_length=sum(seg_length)) %>%
    mutate(arm_median_numbat=sum(median/arm_length*seg_length)) %>% 
    select(chr_arm,arm_median_numbat) %>%
    filter(!duplicated(chr_arm))
  write.csv(chr_arms_num,paste0('data/',pat,'/table_',pat,'_chr_arm_num.csv'))
  
  
  #### inferCNV ####
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
  
  atable<- read.csv(paste0("data/",pat,"/table_",pat,"_inferCNV_geneloc.csv"))
  atable<-data.frame(chr=atable$chrom,
                     pos=atable$start.pos,
                     median=(atable$median - 1))
  infer_cn <- winsorize(data=atable,pos.unit = "bp",verbose=F)
  infer_cn<-pcf(infer_cn,Y = atable,verbose = F)
  #print(plotGenome(data = NULL, segments = infer_cn, pos.unit = "bp", sample = NULL))
  infer_cn$sampleID<-paste0(pat,"_inferCNV")
  df_inf<-infer_cn[infer_cn$chrom %notin% c('X','Y'),]
  
  if(nrow(df_inf)>0){
    inf=TRUE
    df_inf<-data.frame(chr=df_inf$chrom,
                       seg_start=df_inf$start.pos,
                       seg_end=df_inf$end.pos,
                       seg_length=df_inf$end.pos - df_inf$start.pos,
                       median=df_inf$mean)
    df_inf<-na.omit(df_inf)
    
    # Annotate chr arm and split segments spanning both arms
    chr_arms_inf<-annotate_chr_arms(df_inf)
    chr_arms_inf<-chr_arms_inf %>% 
      mutate(chr_arm=paste0(chr,arm)) %>% 
      group_by(chr_arm) %>% 
      mutate(arm_length=sum(seg_length)) %>%
      mutate(arm_median_inferCNV=sum(median/arm_length*seg_length)) %>% 
      select(chr_arm,arm_median_inferCNV) %>%
      filter(!duplicated(chr_arm))
    write.csv(chr_arms_inf,paste0('data/',pat,'/table_',pat,'_chr_arm_inf.csv'))
    
    ## Plot
    chrom_corr<-inner_join(chr_arms_num, chr_arms_inf,by='chr_arm')
    corr_tmp<-cor.test(chrom_corr[,2][[1]], chrom_corr[,3][[1]],method = 'spearman')
    p_num_inf<-chrom_corr %>% ggplot(aes(x=.[[2]], y=.[[3]],col=chr_arm)) + 
            geom_point()+
            theme_bw()+
            geom_smooth(method=lm,  col='darkgray',lty=2)+
            ggtitle(paste0(pat,': ',names(chrom_corr)[2],' vs ',names(chrom_corr)[3],'\nRho: ',
                           round(corr_tmp$estimate[[1]],3),
                           ', cor P = ',
                           scientific(corr_tmp$p.value,digits = 3)))+
            ggrepel::geom_text_repel(aes(label = chr_arm),size=4,max.overlaps = 44)+
            xlab(names(chrom_corr)[2])+
            ylab(names(chrom_corr)[3])+
            theme(legend.position = 'none')
    
    collect_corrs<-rbind.data.frame(collect_corrs,
                                    data.frame(sample=pat,
                                               num_inf_Rho=corr_tmp$estimate[[1]],
                                               num_inf_P=corr_tmp$p.value,
                                               num_wgs_Rho=NA,
                                               num_wgs_P=NA,
                                               inf_wgs_Rho=NA,
                                               inf_wgs_P=NA))
  }
  
  #### ichorCNA ####
  if(pat %in% pats_wgs){
    wgs<-TRUE
    system(paste0("aws s3 sync s3://melanoma-ribas/ichorCNA/",pat,"/ data/",pat,"/ichorCNA/ --exclude '*' --include '*.seg' "))
    df_wgs<-read.table(paste0('data/',pat,'/ichorCNA/',
                           grep(pattern = 'tumor.seg$',
                                x = list.files(paste0('data/',pat,'/ichorCNA/')),
                                value = T)),
                    header = T)
    df_wgs<-data.frame(chr=df_wgs$chr,
                       seg_start=df_wgs$start,
                       seg_end=df_wgs$end,
                       seg_length=df_wgs$end - df_wgs$start,
                       median=df_wgs$median)
    
    # Annotate chr arm and split segments spanning both arms
    chr_arms_wgs<-annotate_chr_arms(df_wgs)
    chr_arms_wgs<-chr_arms_wgs %>% 
      mutate(chr_arm=paste0(chr,arm)) %>% 
      group_by(chr_arm) %>% 
      mutate(arm_length=sum(seg_length)) %>%
      mutate(arm_median_wgs=sum(median/arm_length*seg_length)) %>% 
      select(chr_arm,arm_median_wgs) %>%
      filter(!duplicated(chr_arm))
    write.csv(chr_arms_wgs,paste0('data/',pat,'/table_',pat,'_chr_arm_wgs.csv'))
    
    ## Plot
    chrom_corr<-inner_join(chr_arms_num, chr_arms_wgs,by='chr_arm')
    corr_tmp<-cor.test(chrom_corr[,2][[1]], chrom_corr[,3][[1]],method = 'spearman')
    p_num_wgs<-chrom_corr %>% ggplot(aes(x=.[[2]], y=.[[3]],col=chr_arm)) + 
      geom_point()+
      theme_bw()+
      geom_smooth(method=lm,  col='darkgray',lty=2)+
      ggtitle(paste0(pat,': ',names(chrom_corr)[2],' vs ',names(chrom_corr)[3],'\nRho: ',
                     round(corr_tmp$estimate[[1]],3),
                     ', cor P = ',
                     scientific(corr_tmp$p.value,digits = 3)))+
      ggrepel::geom_text_repel(aes(label = chr_arm),size=4,max.overlaps = 44)+
      xlab(names(chrom_corr)[2])+
      ylab(names(chrom_corr)[3])+
      theme(legend.position = 'none')
    
    collect_corrs<-rbind.data.frame(collect_corrs,
                                    data.frame(sample=pat,
                                               num_inf_Rho=NA,
                                               num_inf_P=NA,
                                               num_wgs_Rho=corr_tmp$estimate[[1]],
                                               num_wgs_P=corr_tmp$p.value,
                                               inf_wgs_Rho=NA,
                                               inf_wgs_P=NA))
  }
  
  if(inf==TRUE & wgs==TRUE){
    ## Plot
    chrom_corr<-inner_join(chr_arms_inf, chr_arms_wgs,by='chr_arm')
    corr_tmp<-cor.test(chrom_corr[,2][[1]], chrom_corr[,3][[1]],method = 'spearman')
    p_inf_wgs<-chrom_corr %>% ggplot(aes(x=.[[2]], y=.[[3]],col=chr_arm)) + 
      geom_point()+
      theme_bw()+
      geom_smooth(method=lm,  col='darkgray',lty=2)+
      ggtitle(paste0(pat,': ',names(chrom_corr)[2],' vs ',names(chrom_corr)[3],'\nRho: ',
                     round(corr_tmp$estimate[[1]],3),
                     ', cor P = ',
                     scientific(corr_tmp$p.value,digits = 3)))+
      ggrepel::geom_text_repel(aes(label = chr_arm),size=4,max.overlaps = 44)+
      xlab(names(chrom_corr)[2])+
      ylab(names(chrom_corr)[3])+
      theme(legend.position = 'none')
    
    collect_corrs<-rbind.data.frame(collect_corrs,
                                    data.frame(sample=pat,
                                               num_inf_Rho=NA,
                                               num_inf_P=NA,
                                               num_wgs_Rho=NA,
                                               num_wgs_P=NA,
                                               inf_wgs_Rho=corr_tmp$estimate[[1]],
                                               inf_wgs_P=corr_tmp$p.value))
  }
  
  if(inf==TRUE && wgs==TRUE){
    print(p_num_inf + p_num_wgs + p_inf_wgs)
  }else if(inf==TRUE){
    print(p_num_inf + plot_spacer())
  }else if(wgs==TRUE){
    print(p_num_wgs + plot_spacer())
  }
}

# Summarize Rhos and plot
collect_corrs<-collect_corrs %>%
  group_by(sample) %>% 
  summarise_all(list(~first(na.omit(.))))
write.csv(collect_corrs,paste0("data/tumor/table_melanoma_chr_arm_corr.csv"),row.names = F)

plot_corrs<-collect_corrs %>% 
  select(sample, grep('Rho',names(collect_corrs),T)) %>%
  melt()

plot_corrs2<-collect_corrs %>% 
  select(sample, grep('P',names(collect_corrs),T)) %>%
  melt()
plot_corrs$pvalue<-plot_corrs2$value

print(ggplot(plot_corrs,aes(x=sample,y=value,col=variable,size=-log10(pvalue)))+
  geom_point(alpha=0.9)+
  theme_bw() +
  ylab('Rho') +
  ylim(min(plot_corrs$value),1)+
  ggtitle('Rho per sample and comparison') +
  geom_hline(yintercept = 0, size=0.2,col='black',linetype=2)+
  labs(color='Comparison') +
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5),
        legend.position = 'bottom'))

print(ggplot(plot_corrs,aes(x=variable,y=value,col=variable))+
  geom_boxplot()+
  geom_jitter()+
  theme_bw() +
  ylab('Rho') +
  xlab('Comparison')+
  ggtitle('Rho by comparison') +
  geom_hline(yintercept = 0, size=0.2,col='black',linetype=2)+
  labs(color='Comparison') + plot_spacer())

dev.off()

