#!/usr/bin/env Rscript

#### Adjusted by: Jana Biermann, PhD
#### Run KINOMO for in first round using ranks 2-10 and plot for rank 4


# Author: Somnath Tagore, Ph.D. Title: Running KINOMO for performing Non-negative Matrix Factorization using gene expression data 
# Script Name: kinomo_run.R
# Last Updated: 01/24/2022

#Instructions
#The KINOMO repsository can be accessed via https://github.com/IzarLab/KINOMO.git
#Convert the gene expression data (raw counts) to Seurat object and perform necessary normalization, scaling etc.
#Run the following R script


## NMF

print(Sys.time())

library(dplyr)
library(Seurat)
library(ggplot2)
library(gplots)
library(NMF)
library(registry)
library(rngtools)
library(purrr)
library(cowplot)
library(stringr)
library(pkgmaker)
library(cluster)

#system('git clone https://github.com/IzarLab/KINOMO.git')


setwd('~/KINOMO')
pat <- commandArgs()[6]

directory<-paste0("~/KINOMO/output/",pat,"/")

seu<-readRDS('~/data/tumor/data_tum_merged_v4.rds')
seu<-subset(seu,orig.ident==pat)
DefaultAssay(seu) <- "RNA"
seu<-DietSeurat(seu)

# Remove MT & ribosomal genes
mito.genes <- grep(pattern = "^MT-", x = rownames(seu), value = TRUE)
rbl.genes <- grep(pattern = "^RB-", x = rownames(seu), value = TRUE)
rsl.genes <- grep(pattern = "^RS-", x = rownames(seu), value = TRUE)
rpl.genes <- grep(pattern = "^RPL-", x = rownames(seu), value = TRUE)
rbl.genes <- grep(pattern = "^RBL-", x = rownames(seu), value = TRUE)
rps.genes <- grep(pattern = "^RPS-", x = rownames(seu), value = TRUE)
rbs.genes <- grep(pattern = "^RBS-", x = rownames(seu), value = TRUE)
rbl1.genes <- grep(pattern = "^RB", x = rownames(seu), value = TRUE)
rsl1.genes <- grep(pattern = "^RS", x = rownames(seu), value = TRUE)
rpl1.genes <- grep(pattern = "^RPL", x = rownames(seu), value = TRUE)
rbl1.genes <- grep(pattern = "^RBL", x = rownames(seu), value = TRUE)
rps1.genes <- grep(pattern = "^RPS", x = rownames(seu), value = TRUE)
rbs1.genes <- grep(pattern = "^RBS", x = rownames(seu), value = TRUE)

counts <- GetAssayData(seu, assay = "RNA")
counts <- counts[-(which(rownames(counts) %in% c(mito.genes,rbl.genes,rsl.genes,rpl.genes,rbl.genes,rps.genes,rbs.genes,rbl1.genes,rsl1.genes,rpl1.genes,rbl1.genes,rps1.genes,rbs1.genes))),]
seu <- subset(seu, features = rownames(counts))

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu, features = rownames(seu))
seu <- RunPCA(seu, features = VariableFeatures(object = seu))
seu <- RunUMAP(seu, dims = 1:20)

mat<-as.matrix(seu@assays$RNA@data)
mat<-mat[rowSums(mat)>0,]

###Estimating the factorization rank
#A critical parameter in KINOMO-NMF is the factorization rank r. 
#It defines the number of metagenes used to approximate the target matrix. Given a 
#NMF method and the target matrix, a common way of deciding on r is to try different values, 
#compute some quality measure of the results, and choose the best value 
#according to this quality criteria.
#Several approaches have then been proposed to choose the optimal value of r. 
#For example, (Brunet2004) proposed to take the first value of r 
#for which the cophenetic coefficient starts decreasing, 
#(Hutchins2008) suggested to choose the first value where the RSS curve presents an 
#inflection point, and (Frigyesi2008) considered the smallest value at which 
#the decrease in the RSS is lower than the decrease of the RSS obtained 
#from random data.

ifelse(!dir.exists(file.path(directory)), 
       dir.create(file.path(directory),recursive = T), FALSE)
pdf(file = paste0(directory,pat,"_KINOMO_nmf.pdf"),height = 20,width=20)
#nmf.estimate.rank <- KINOMO(mat, 2:10, nrun=10)
nmf.estimate.rank <- nmf(mat, 2:10, nrun=10)
plot(nmf.estimate.rank)
dev.off()

saveRDS(nmf.estimate.rank,paste0(directory,pat,"_KINOMO_nmf_estimate.rds"))


#After the script successfully runs, the following files would be created in the current directory:
# 1. Rank plot
# 2. Meta-gene signature UMAPS
# 3. Meta-gene heatmaps
# 4. Associated .csv files
# 5. Associated .rds files
