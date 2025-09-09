#!/usr/bin/Rscript

#### Title: Get filtered original barcodes per melanoma sample (input for numbat)
#### Author: Jana Biermann, PhD

library(dplyr)
library(Seurat)


seu_all<-readRDS('data/merged/data_melanoma_merged_v2.rds')
tmp<-seu_all@meta.data

for(pat in unique(tmp$sample)){
  sub<-tmp%>% filter(sample==pat) %>% select(barcode_orig) %>% unlist
  ifelse(!dir.exists(file.path(paste0('data/',pat))), 
         dir.create(file.path(paste0('data/',pat)),recursive = T), FALSE)
  write.table(sub, paste0('data/',pat,'/barcodes_',pat,'_filtered.csv'),
               sep=",",  row.names = F, col.names=FALSE, quote = FALSE)
}
