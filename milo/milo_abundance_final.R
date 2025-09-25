#!/usr/bin/Rscript

#### Title: Milo with cell type as input argument
#### Author: Kevin Hoffer-Hawlik (adapted from Jana Biermann, PhD)

#BiocManager::install("miloR")

library(dplyr)
library(ggplot2)
library(ggrepel)
library(scales)
library(Seurat)
library(reshape2)
library(viridis)
library(RColorBrewer)
library(ggpubr)
library(miloR)
library(patchwork)
library(SingleCellExperiment)
library(tidyr)
library(ggrastr)
library(grid)
'%notin%' <- Negate('%in%')

setwd('~/Documents/khh/melanoma/')
directory<-'~/Documents/khh/melanoma/'

print(paste0('##### Start: ',Sys.time()))

# celltype <- commandArgs()[6]
celltype <- 'fine_T_mac_subsampled'

#### Create Milo object myeloid ####
# seu<-readRDS('data_melanoma_integrated_RPCA_v4.rds')
# sce <- as.SingleCellExperiment(seu)
# my_milo <- Milo(sce)
# celltype<-'myeloid'
# k<-40
# d<- 30
# prop <- 0.1
# reduced.dim <- "PCA"

#### Create Milo object main ####
# if(celltype=='main'){
#   seu<-readRDS('data_melanoma_integrated_RPCA_v4.rds')
#   seu$cell_type_fine<- seu$cell_type_main
#   sce <- as.SingleCellExperiment(seu, assay = "integrated")
#   my_milo <- Milo(sce)
#   celltype<-'main'
#   k<-65
#   d<- 30
#   prop <- 0.2
#   reduced.dim <- "PCA"
# }


#### Create Milo object fine ####
# if(celltype=='fine'){
#   seu<-readRDS('data_melanoma_integrated_RPCA_v4.rds')
#   seu$cell_type_fine[seu$cell_type_fine=='Naive Tregs'] <- "Tregs"
#   seu$cell_type_fine<- paste0(seu$cell_type_main, ' - ',seu$cell_type_fine)
#   seu$FvM<-ifelse(seu$Gender =='F','F','M')
#   sce <- as.SingleCellExperiment(seu, assay = "integrated")
#   my_milo <- Milo(sce)
#   celltype<-'fine'
#   k<-65
#   d<- 30
#   prop <- 0.2
#   reduced.dim <- "PCA"
#   # condis<-'FvM'
# }

#### Create Milo object with fine annotations just for T/NK and myeloid cells ####
# if(celltype=='fine_T_mac'){
#   seu<-readRDS('data_melanoma_integrated_RPCA_v4.rds')
#   
#   # fix labels and collapse less interesting fine cell type labels
#   seu$cell_type_fine[seu$cell_type_fine=='Naive Tregs'] <- "Tregs"
#   seu$cell_type_fine<- paste0(seu$cell_type_main, ' - ',seu$cell_type_fine)
#   seu$cell_type_fine[seu$cell_type_fine %in% c("B/Plasma cells - Class-switched memory B cells", "B/Plasma cells - Cycling cells", "B/Plasma cells - Memory B cells", "B/Plasma cells - Naive B cells", "B/Plasma cells - Plasma cells")] <- "B/Plasma cells"
#   seu$cell_type_fine[seu$cell_type_fine %in% c("Endothelial cells - Arterial EC", "Endothelial cells - Capillary EC", "Endothelial cells - Cycling cells", "Endothelial cells - Lymphatic EC", "Endothelial cells - Venous EC")] <- "Endothelial cells"
#   seu$cell_type_fine[seu$cell_type_fine %in% c("Stromal cells - CAFs", "Stromal cells - CAFs (antigen-pres.)", "Stromal cells - CAFs (cytokines)", "Stromal cells - CAFs (ECM)", "Stromal cells - Cycling cells", "Stromal cells - Pericytes")] <- "Stromal cells"
#   
#   sce <- as.SingleCellExperiment(seu, assay = "integrated")
#   my_milo <- Milo(sce)
#   celltype<-'fine'
#   k<-65
#   d<- 30
#   prop <- 0.2
#   reduced.dim <- "PCA"
#   # condis<-'FvM'
# }

#### Create Milo object tcells only ####
# if(celltype=='tcell'){
#   system("aws s3 sync s3://melanoma-ribas/Seurat/annotation/tcells/ data/annotation/tcells/ --exclude '*' --include 'data_melanoma_tcells_v4.rds' --quiet ")
#   seu<-readRDS('data/annotation/tcells/data_melanoma_tcells_v4.rds')
#   sce <- as.SingleCellExperiment(seu)
#   my_milo <- Milo(sce)
#   celltype<-'tcells'
#   k<-70
#   d<- 30
#   prop <- 0.1
#   reduced.dim <- "PCA"
# }

#### Create Milo object with fine annotations just for T/NK and myeloid cells, while subsampling to median patient counts ####
if(celltype=='fine_T_mac_subsampled'){
  seu<-readRDS('data_melanoma_integrated_RPCA_v4.rds')
  
  # fix labels
  seu$treated[seu$sample=='F16_post1_pre2'] <- 'Treated'
  seu$cell_type_fine[seu$cell_type_fine=='Naive Tregs'] <- "Tregs"
  myeloid_annotations <- read.csv('myeloid_scvi_annotations_032025.csv')
  seu$cell_type_fine[myeloid_annotations$X] <- myeloid_annotations$cell_type_fine
  
  # collapse less interesting fine cell type labels
  seu$cell_type_fine<- paste0(seu$cell_type_main, ' - ',seu$cell_type_fine)
  seu$cell_type_fine[seu$cell_type_fine %in% c("B/Plasma cells - Class-switched memory B cells", "B/Plasma cells - Cycling cells", "B/Plasma cells - Memory B cells", "B/Plasma cells - Naive B cells", "B/Plasma cells - Plasma cells")] <- "B/Plasma cells"
  seu$cell_type_fine[seu$cell_type_fine %in% c("Endothelial cells - Arterial EC", "Endothelial cells - Capillary EC", "Endothelial cells - Cycling cells", "Endothelial cells - Lymphatic EC", "Endothelial cells - Venous EC")] <- "Endothelial cells"
  seu$cell_type_fine[seu$cell_type_fine %in% c("Stromal cells - CAFs", "Stromal cells - CAFs (antigen-pres.)", "Stromal cells - CAFs (cytokines)", "Stromal cells - CAFs (ECM)", "Stromal cells - Cycling cells", "Stromal cells - Pericytes")] <- "Stromal cells"
  
  # make summary of cell types across patients
  pdf(paste0(directory,'cell_type_main_plot.pdf'))
  df_prop <- seu@meta.data %>%
    group_by(patient, cell_type_main) %>%
    summarise(count=n(), .groups="drop") %>%
    group_by(patient) %>%
    mutate(proportion = count / sum(count))
  df_totals <- seu@meta.data %>%
    group_by(patient) %>%
    summarise(total_cells=n())
  g <- ggplot(data=df_prop, aes(x=patient, y=proportion, fill=cell_type_main)) + 
    geom_bar(stat="identity") + labs(y="Proportion", x="Patient", fill="Cell Type") + 
    geom_text(data = df_totals,
              aes(x=patient, y=1.02, label=total_cells),
              inherit.aes = F, vjust=0.5, size=3.5, angle=90) + 
    theme(
      panel.background = element_blank(), panel.grid = element_blank(), 
      axis.line = element_line(color="black"), axis.ticks = element_line(color="black"), 
      axis.text.x = element_text(angle=90, vjust=0.5)
    )
  print(g)
  g <- ggplot(data=df_totals, aes(x=patient, y=total_cells)) + 
    geom_bar(stat="identity") + labs(y="Total Cell Counts", x="Patient") + 
    theme(
      panel.background = element_blank(), panel.grid = element_blank(), 
      axis.line = element_line(color="black"), axis.ticks = element_line(color="black"), 
      axis.text.x = element_text(angle=90, vjust = 0.5)
    )
  print(g)
  dev.off()
  
  # compare cell type proportions in R vs. NR, treated vs. untreated, mucosal vs. non-mucosal
  df_prop <- seu@meta.data %>%
    group_by(sample, cell_type_fine, treated, mucosal, responder) %>%
    summarise(count=n(), .groups="drop") %>%
    group_by(sample) %>%
    mutate(proportion = count / sum(count))
  write.csv(df_prop, 'cell_type_fine_proportions.csv')

  # remove tumors and lowly represented myeloid cell types
  remove <- c("Tumor cells - Cycling tumor cells", "Tumor cells - Tumor cells",
              "Myeloid cells - Langerhans Cell", "Myeloid cells - Mast cells",
              "Myeloid cells - Microglial Cell", "Myeloid cells - Osteoclast-like Cell"
              )
  mask <- seu@meta.data$cell_type_fine %notin% remove
  seu <- seu[,mask]
  
  # subsample so cells from one patient is not higher than median counts
  cell_counts_per_patient <- seu@meta.data %>% group_by(patient) %>% summarise(cell_count = n())
  median_cells_per_patient <- median(cell_counts_per_patient$cell_count)
  mask <- seu@meta.data %>% group_by(patient) %>% mutate(selected=row_number() <= min(median_cells_per_patient, n())) %>% ungroup %>% .$selected
  seu <- seu[,mask]
  
  DefaultAssay(seu) <- "RNA"
  seu <- seu %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% 
    RunPCA(npcs=30) %>% FindNeighbors(dims=1:30, k.param=30) %>% RunUMAP(dims=1:30)
  
  sce <- as.SingleCellExperiment(seu)
  my_milo <- Milo(sce)
  celltype<-'fine'
  k<- 90
  d<- 30
  prop <- 0.1
  reduced.dim <- "PCA"
}

#### Plots part 1 ####
pdf(paste0(directory,'plots_milo_',celltype,'_k=',k,'_d=',d,'_prop=',prop,'_pt1.pdf'))

DimPlot(seu, label = F, group.by = 'cell_type_main', repel = T,
        raster=T,pt.size = 2)+coord_fixed()

DimPlot(seu, label = F, group.by = 'cell_type_fine', repel = T,
        raster=T,pt.size = 2)+coord_fixed()

### Construct KNN graph
traj_milo <- buildGraph(my_milo, k = k, d = d, reduced.dim = reduced.dim)

### Defining representative neighbourhoods
# prop=0.1 for datasets of less than 30k cells
traj_milo <- makeNhoods(traj_milo, prop = prop, k = k, d=d, refined = T, reduced_dims = reduced.dim)

# want to have an average neighbourhood size over 5 x N_samples
# Otherwise consider rerunning makeNhoods increasing k and/or prop
text <- paste0('Average neighbourhood size should be over 5 x N_samples',
               '\n5 x N_samples: ', length(unique(seu$sample)) * 5,
              '\nNumber of samples: ', length(unique(seu$sample)),
              '\nNumber of cells: ', dim(seu)[2],
              '\nCurrent k: ',k,
              '\nCurrent prop: ',prop,
              '\nCurrent average nhood size: ', round(mean(colSums(nhoods(traj_milo))), digits = 2))
grid.newpage()
print(grid.text(label = text,just = "centre"))

print(plotNhoodSizeHist(traj_milo)+
        ggtitle(paste0('Histogram of #cells belonging to each neighbourhood')))
dev.off()

### Calculate within neighbourhood distances
# adds to the Milo object a n×m matrix, n = number of nhoods, m = number of samples
# Values indicate the number of cells from each sample counted in a neighbourhood
traj_milo <- calcNhoodDistance(traj_milo, d=d)

### Counting cells in neighbourhoods
traj_milo <- countCells(traj_milo, 
                        meta.data = data.frame(colData(traj_milo)), 
                        samples="sample")
head(nhoodCounts(traj_milo))
print(paste0('##### Finished part 1: ',Sys.time()))

rm(sce)
rm(seu)
save.image(paste0(directory,'image_milo_',celltype,'_k=',k,'_d=',d,'_prop=',prop,'.RData'))
print(paste0('##### Saved image: ',Sys.time()))

#load(paste0(directory,'image_milo_',celltype,'_k=',k,'_d=',d,'_prop=',prop,'.RData'))


#### Plots part 2 & 3 ####
# 'treated','responder','mucosal'
if(!exists('condis')){condis<-c('responder','mucosal','treated')}

for(condi in condis){
  print(paste0('##### Processing ',condi,': ',Sys.time()))
  
  pdf(paste0(directory,'plots_milo_',celltype,'_k=',k,'_d=',d,'_prop=',prop,'_pt2_',condi,'.pdf'))
  
  ### Differential abundance testing
  # hypothesis testing in a generalized linear model (GLM) framework
  # using the Negative Binomial GLM implementation
  # Design dataframe
  traj_design <- data.frame(colData(traj_milo))[,c("sample", condi)]
  traj_design <- distinct(traj_design)
  rownames(traj_design) <- traj_design$sample
  traj_design
  
  # Reorder rownames to match columns of nhoodCounts(milo)
  traj_design <- traj_design[colnames(nhoodCounts(traj_milo)), , drop=FALSE]
  traj_design
  
  # Perform differential neighbourhood abundance testing
  rownames(traj_design) <- traj_design$sample
  if(condi=='mucosal'){traj_design[,2]<-factor(traj_design[,2],levels = c('Non-mucosal','Mucosal'))}
  if(condi=='treated'){traj_design[,2]<-factor(traj_design[,2],levels = c('Untreated','Treated'))}
  if(condi=='responder'){traj_design[,2]<-factor(traj_design[,2],levels = c('Non-responder','Responder'))}
  if(condi=='CTvMT'){traj_design[,2]<-factor(traj_design[,2],levels = c('treated_individual','treated_combination'))}
  if(condi=='FvM'){traj_design[,2]<-factor(traj_design[,2],levels = c('M','F'))}
  
  da_results <- testNhoods(traj_milo, design = ~ traj_design[,2], design.df = traj_design)
  print(da_results %>%
    arrange(SpatialFDR) %>%
    head())
  
  print(ggplot(da_results, aes(PValue)) + 
          geom_histogram(bins=50)+
          theme_classic()+
          ggtitle(paste0('P-value distribution for condition: ',condi)))
  
  print(ggplot(da_results, aes(logFC, -log10(SpatialFDR))) + 
          geom_point() +
          theme_classic()+
          geom_hline(yintercept = 1)+
          ggtitle(paste0('Volcano plot for condition: ',condi)))
  
  ### Annotate nhoods
  da_results <- annotateNhoods(traj_milo, da_results, coldata_col = "cell_type_fine")
  head(da_results)
  
  write.csv(da_results,
            paste0(directory,'table_milo_',celltype,'_k=',k,'_d=',d,'_prop=',prop,'_',condi,'.csv'),
            row.names = F)
  
  print(ggplot(da_results, aes(logFC, -log10(SpatialFDR), col=da_results[,8])) +
            geom_point(alpha=0.8) +
            theme_bw() +
            geom_hline(yintercept = 1)+
            ggtitle(paste0('Volcano plot for condition: ',condi)))
  
  ### Visualize neighbourhoods displaying DA
  traj_milo <- buildNhoodGraph(traj_milo)
  
  print(plotNhoodGraphDA(traj_milo, da_results, alpha=0.05,show_groups='cell_type_fine')+
          ggtitle(paste0('Neighborhood graph\nCondition: ',condi)))
  dev.off()
  
  
  #### Plots part 3
  #system(paste0("aws s3 sync s3://melanoma-ribas/Seurat/global/ data/global/ --exclude '*' --include '*",celltype,"*.csv' "))
  #condi<-'treated'
  #condi<-'responder'
  #condi<-'CTvMT'
  #condi<-'mucosal'
  
  da_results<-read.csv(paste0(directory,'table_milo_',celltype,'_k=',k,'_d=',d,'_prop=',prop,'_',condi,'.csv'))
  
  pdf(paste0(directory,'plots_milo_',celltype,'_k=',k,'_d=',d,'_prop=',prop,'_pt3_',condi,'.pdf'))
  print(ggplot(da_results, aes(cell_type_fine_fraction)) + 
          geom_histogram(bins=50)+
          theme_classic()+
          ggtitle('Cell type assignment to neighborhoods'))
  
  signif_level<-ifelse(min(da_results$SpatialFDR)> 0.1,1,0.1)
  print(plotDAbeeswarm(da_results, group.by = "cell_type_fine",alpha = signif_level)+
          ggtitle(paste0('Condition: ',condi)))
  
  df <- tibble(cell_type = da_results$cell_type_fine,
               logFC = da_results$logFC,
               SpatialFDR = da_results$SpatialFDR,
               significant=ifelse(-log10(da_results$SpatialFDR)>1,0.6,0.1),
               direction=ifelse(-log10(da_results$SpatialFDR)>1 & da_results$logFC > 0,'up',
                                ifelse(-log10(da_results$SpatialFDR)>1 & da_results$logFC <= 0,'down',
                                       'NS')))
  if('down' %notin% df$direction){
    colorvector<-c('gray','#BC3D41')
  }else{
    colorvector<-c('#4F7EBB','gray','#BC3D41')
  }
  
  avg <- df %>%
    group_by(cell_type) %>%
    summarise(avg = mean(logFC),
              median=median(logFC))
  
  print(ggplot(df, aes(logFC, cell_type)) +
    geom_vline(xintercept = 0,lty='dashed',col='gray')+
    geom_boxplot(col='darkgray',fill=NA,outlier.shape = NA) +
    geom_jitter(aes(size = -log10(SpatialFDR), 
                    col = cell_type,
                    alpha = -log10(SpatialFDR)), 
                height = 0.2) +
    geom_point(data = avg, 
               aes(x = avg, y = cell_type)) +
    theme_bw() +
    theme(panel.grid=element_blank(),
          axis.text = element_text(color='black'))+
    ggtitle(paste0('Cell type: ',celltype,'\nCondition: ',condi)))
  
  print(ggplot(df %>% filter(-log10(SpatialFDR) >1), aes(logFC, cell_type)) +
          geom_vline(xintercept = 0,lty='dashed',col='gray')+
          geom_jitter(aes(size = -log10(SpatialFDR), 
                          col = cell_type,
                          alpha = -log10(SpatialFDR)), 
                      height = 0.2) +
          theme_bw() +
          theme(panel.grid=element_blank(),
                axis.text = element_text(color='black'))+
          ggtitle(paste0('Cell type: ',celltype,' - only significant nhoods','\nCondition: ',condi)))
  
  print(ggplot(df, aes(logFC, cell_type)) +
          geom_vline(xintercept = 0,lty='dashed',col='gray')+
          geom_jitter(aes(size = -log10(SpatialFDR), 
                          col = direction),
                      alpha = df$significant, 
                      height = 0.2) +
          scale_color_manual(values = colorvector)+
          geom_point(data = avg, 
                     aes(x = avg, y = cell_type)) +
          geom_point(data = avg, 
                     aes(x = median, y = cell_type),shape='|',size=4) +
          theme_bw() +
          theme(panel.grid=element_blank(),
                axis.text = element_text(color='black'))+
          ggtitle(paste0('Cell type: ',celltype,'\nCondition: ',condi)))
  
  print(ggplot(df, aes(logFC, cell_type)) +
          geom_vline(xintercept = 0,lty='dashed',col='gray')+
          geom_jitter(aes(size = -log10(SpatialFDR), 
                          col = cell_type),
                          alpha = df$significant, 
                      height = 0.2) +
          geom_point(data = avg, 
                     aes(x = avg, y = cell_type)) +
          geom_point(data = avg, 
                     aes(x = median, y = cell_type),shape='|',size=4) +
          theme_bw() +
          theme(panel.grid=element_blank(),
                axis.text = element_text(color='black'))+
          ggtitle(paste0('Cell type: ',celltype,'\nCondition: ',condi)))
  
  print(ggplot(df, aes(logFC, cell_type)) +
          geom_vline(xintercept = 0,lty='dashed',col='gray')+
          geom_jitter(aes(size = -log10(SpatialFDR), 
                          col = cell_type),
                      alpha = df$significant, 
                      height = 0.2) +
          theme_bw() +
          theme(panel.grid=element_blank(),
                axis.text = element_text(color='black'))+
          ggtitle(paste0('Cell type: ',celltype,'\nCondition: ',condi)))
  
  
  da_results$cell_type_fine <- ifelse(da_results$cell_type_fine_fraction < 0.5,
                                      "Mixed", da_results$cell_type_fine)
  print(plotDAbeeswarm(da_results, group.by = "cell_type_fine",alpha = signif_level)+
          ggtitle(paste0('Condition: ',condi)))
  
  df <- tibble(cell_type = da_results$cell_type_fine,
               logFC = da_results$logFC,
               SpatialFDR = da_results$SpatialFDR,
               significant=ifelse(-log10(da_results$SpatialFDR)>1,0.6,0.1),
               direction=ifelse(-log10(da_results$SpatialFDR)>1 & da_results$logFC > 0,'up',
                                ifelse(-log10(da_results$SpatialFDR)>1 & da_results$logFC <= 0,'down',
                                       'NS')))
  if('down' %notin% df$direction){
    colorvector<-c('gray','#BC3D41')
  }else{
    colorvector<-c('#4F7EBB','gray','#BC3D41')
  }
  
  avg <- df %>%
    group_by(cell_type) %>%
    summarise(avg = mean(logFC),
              median=median(logFC))
  
  print(ggplot(df, aes(logFC, cell_type)) +
    geom_boxplot(col='darkgray',fill=NA,outlier.shape = NA) +
    geom_jitter(aes(size = -log10(SpatialFDR), 
                    col = cell_type,
                    alpha = -log10(SpatialFDR)), 
                height = 0.15) +
    geom_point(data = avg, 
               aes(x = avg, y = cell_type)) +
    theme_bw() +
    theme(panel.grid=element_blank(),
          axis.text = element_text(color='black'))+
    geom_vline(xintercept = 0,lty='dashed',col='gray')+
    ggtitle(paste0('Cell type: ',celltype,'\nCondition: ',condi)))
  
  print(ggplot(df, aes(logFC, cell_type)) +
          geom_vline(xintercept = 0,lty='dashed',col='gray')+
          geom_jitter(aes(size = -log10(SpatialFDR), 
                          col = cell_type),
                      alpha = df$significant, 
                      height = 0.2) +
          geom_point(data = avg, 
                     aes(x = avg, y = cell_type)) +
          geom_point(data = avg, 
                     aes(x = median, y = cell_type),shape='|',size=4) +
          theme_bw() +
          theme(panel.grid=element_blank(),
                axis.text = element_text(color='black'))+
          ggtitle(paste0('Cell type: ',celltype,'\nCondition: ',condi)))
  
  print(ggplot(df, aes(logFC, cell_type)) +
          geom_vline(xintercept = 0,lty='dashed',col='gray')+
          geom_jitter(aes(size = -log10(SpatialFDR), 
                          col = direction),
                      alpha = df$significant, 
                      height = 0.2) +
          scale_color_manual(values = colorvector)+
          geom_point(data = avg, 
                     aes(x = avg, y = cell_type)) +
          geom_point(data = avg, 
                     aes(x = median, y = cell_type),shape='|',size=4) +
          theme_bw() +
          theme(panel.grid=element_blank(),
                axis.text = element_text(color='black'))+
          ggtitle(paste0('Cell type: ',celltype,'\nCondition: ',condi)))
  
  print(ggplot(df, aes(logFC, cell_type)) +
          geom_vline(xintercept = 0,lty='dashed',col='gray')+
          geom_jitter(aes(size = -log10(SpatialFDR), 
                          col = cell_type),
                      alpha = df$significant, 
                      height = 0.2) +
          theme_bw() +
          theme(panel.grid=element_blank(),
                axis.text = element_text(color='black'))+
          ggtitle(paste0('Cell type: ',celltype,'\nCondition: ',condi)))
  
  # df_tmp<-seu_meta %>%
  #   group_by(sample, cell_type_fine, !!sym(condi)) %>%
  #   tally() %>%
  #   group_by(sample) %>%
  #   mutate(freq = n/sum(n))
  # 
  # print(ggplot(df_tmp,aes(x=cell_type_fine, freq,col=!!sym(condi))) + 
  #         geom_boxplot(outlier.shape = NA) + 
  #         geom_point(position=position_jitterdodge(jitter.width = 0.2),size=1,alpha=0.7)+
  #         stat_compare_means(method = "wilcox.test", label='p.format',size=3) + 
  #         theme_classic() +
  #         theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5),
  #               legend.position = 'bottom'))
  dev.off()
  
  
  tiff(paste0(directory,'plots_milo_',celltype,'_k=',k,'_d=',d,'_prop=',prop,'_',condi,'_net.tiff'),
       width = 8, height = 8, units = "in",res = 300)
  print(plotNhoodGraphDA(traj_milo, da_results, alpha=0.05,show_groups='cell_type_fine')+
          ggtitle(paste0('Neighborhood graph\nCondition: ',condi)))
  dev.off()
  
}

# system(paste0("aws s3 sync data/global/ s3://melanoma-ribas/Seurat/global/ --exclude '*' --include '*",celltype,"*' --quiet "))

print(paste0('##### Done: ',Sys.time()))

#### Bespoke refined plots ####

directory<-'~/Documents/khh/melanoma/'
celltype<-'fine'
k<-90
d<-30
prop<-0.1
# condi<-'treated' # 'responder' 'mucosal' 'treated'
for(condi in c('responder', 'mucosal', 'treated')){

  da_results<-read.csv(paste0(directory,'table_milo_',celltype,'_k=',k,'_d=',d,'_prop=',prop,'_',condi,'.csv'))
  
  pdf(paste0(directory,'plots_milo_',celltype,'_k=',k,'_d=',d,'_prop=',prop,'_pt3_',condi,'_refined.pdf'))
  
  signif_level<-ifelse(min(da_results$SpatialFDR)> 0.1,1,0.1)
  
  df <- tibble(cell_type = da_results$cell_type_fine,
               logFC = da_results$logFC,
               SpatialFDR = da_results$SpatialFDR,
               significant=ifelse(-log10(da_results$SpatialFDR)>1,0.6,0.1),
               direction=ifelse(-log10(da_results$SpatialFDR)>1 & da_results$logFC > 0,'up',
                                ifelse(-log10(da_results$SpatialFDR)>1 & da_results$logFC <= 0,'down',
                                       'NS')))
  if('down' %notin% df$direction){
    colorvector<-c('gray','#BC3D41')
  }else{
    colorvector<-c('#4F7EBB','gray','#BC3D41')
  }
  
  avg <- df %>%
    group_by(cell_type) %>%
    summarise(avg = mean(logFC),
              median=median(logFC))
  
  print(ggplot(df %>% filter(-log10(SpatialFDR) >1), aes(logFC, cell_type)) +
          geom_vline(xintercept = 0,lty='dashed',col='gray')+
          geom_jitter(aes(size = -log10(SpatialFDR), 
                          col = cell_type,
                          alpha = -log10(SpatialFDR)), 
                      height = 0.2) +
          theme_bw() +
          theme(panel.grid=element_blank(),
                axis.text = element_text(color='black'))+
          ggtitle(paste0('Cell type: ',celltype,' - only significant nhoods','\nCondition: ',condi)))
  
  print(ggplot(df, aes(logFC, cell_type)) +
          geom_vline(xintercept = 0,lty='dashed',col='gray')+
          geom_jitter(aes(size = -log10(SpatialFDR), 
                          col = direction),
                      alpha = df$significant, 
                      height = 0.2) +
          scale_color_manual(values = colorvector)+
          geom_point(data = avg, 
                     aes(x = avg, y = cell_type)) +
          geom_point(data = avg, 
                     aes(x = median, y = cell_type),shape='|',size=4) +
          theme_bw() +
          theme(panel.grid=element_blank(),
                axis.text = element_text(color='black'))+
          ggtitle(paste0('Cell type: ',celltype,'\nCondition: ',condi)))
  
  
  da_results$cell_type_fine <- ifelse(da_results$cell_type_fine_fraction < 0.5,
                                      "Mixed", da_results$cell_type_fine)
  da_results <- da_results[da_results$cell_type_fine!='Mixed',] # remove Mixed
  
  df <- tibble(cell_type = da_results$cell_type_fine,
               logFC = da_results$logFC,
               SpatialFDR = da_results$SpatialFDR,
               significant=ifelse(-log10(da_results$SpatialFDR)>1,0.6,0.1),
               direction=ifelse(-log10(da_results$SpatialFDR)>1 & da_results$logFC > 0,'up',
                                ifelse(-log10(da_results$SpatialFDR)>1 & da_results$logFC <= 0,'down',
                                       'NS')))
  if('down' %notin% df$direction){
    colorvector<-c('gray','#2A286B')
  }else{
    colorvector<-c('#B0C4DE','gray','#2A286B')
  } # R: #2A286B; NR: #B0C4DE
  
  avg <- df %>%
    group_by(cell_type) %>%
    summarise(avg = mean(logFC),
              median=median(logFC))
  
  print(ggplot(df, aes(logFC, cell_type)) +
          geom_vline(xintercept = 0,lty='dashed',col='gray')+
          geom_jitter(aes(col = direction),
                      alpha = df$significant, 
                      height = 0.2) +
          scale_color_manual(values = colorvector)+
          theme_bw() +
          theme(panel.grid=element_blank(),
                axis.text = element_text(color='black'))+
          ggtitle(paste0('Cell type: ',celltype,' - no Mixed nhoods','\nCondition: ',condi)))
  
  print(ggplot(df %>% filter(-log10(SpatialFDR) >1), aes(logFC, cell_type)) +
          geom_vline(xintercept = 0,lty='dashed',col='gray')+
          geom_jitter(aes(size = -log10(SpatialFDR),
                          col = direction),
                      # alpha = df$significant,
                      height = 0.2) +
          scale_color_manual(values = colorvector)+
          theme_bw() +
          theme(panel.grid=element_blank(),
                axis.text = element_text(color='black'))+
          ggtitle(paste0('Cell type: ',celltype,' - no Mixed nhoods, only significant nhoods','\nCondition: ',condi)))
  
  dev.off()
}
