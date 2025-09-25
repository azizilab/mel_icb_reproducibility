#### Title: ContactTracing CUIMC downstream analysis using clone-labelled samples and T cells labelled with dpt
#### Checks for associated pathways using fgsea and hypergeometric GSEA
#### Author: Kevin Hoffer-Hawlik (adapted from Jana Biermann)

# renv::init()
# setwd("/Users/khofferhawlik/interactions_comparison/GSEA (R)")

library(dplyr)
library(Seurat)
library(ggplot2)
library(ggrastr)
library(gplots)
library(hypeR)
library(viridis)
library(scales)
library(RColorBrewer)
library(hypeR)
library(reshape2)
library(ggVennDiagram)
library(ggrepel)
library(UpSetR)
library(psych)
library(patchwork)
library(grid)
library(ggpubr)
library(ggbeeswarm)
# library(raster)
library(Cairo)
'%notin%' <- Negate('%in%')

#### GSEA pathway analysis #### 

library(fgsea)
source('gsea_functions.R')

directory<-'../data/cell_type_ct_dpt/' 

inters_paired<-read.delim2(paste0(directory,'ranked_interactions.txt'))
inters_paired<-inters_paired %>% filter(ligandDE=='True')
inters_paired<-inters_paired %>% filter(fdr_ligand<0.05)
inters_paired<-inters_paired %>% filter(numSigI1_fdr05_receptor>0)

inters_paired$ID<-paste0(inters_paired$ligand,';',inters_paired$cell_type_ligand,';',inters_paired$receptor,';',inters_paired$cell_type_receptor)
inters_paired$ligand_receptor<-paste0(inters_paired$ligand,'_',inters_paired$receptor)
inters_paired$DE<-ifelse(inters_paired$log2FC_ligand>0,'Up','Down')
inters_paired$log2FC_ligand<-as.numeric(inters_paired$log2FC_ligand)
inters_paired$fdr_ligand<-as.numeric(inters_paired$fdr_ligand)
inters_paired$ctLigand_ctReceptor<-paste0(inters_paired$cell_type_ligand,'::',inters_paired$cell_type_receptor)

inters_paired_DE<-inters_paired %>% filter(ligandDE=='True')
inters_paired_DE_up<-inters_paired %>% filter(ligandDE=='True' & DE=='Up')
inters_paired_DE_down<-inters_paired %>% filter(ligandDE=='True' & DE=='Down')

# Run GSEA (select comp for appropriate output file name)
# comp<-'all'

# comp<-'growing_sending'
# comp<-'shrinking_sending'
# comp<-'growing_receiving'
# comp<-'shrinking_receiving'

# comp<-'cd8_dpt1_sending'
# comp<-'cd8_dpt2_sending'
# comp<-'cd8_dpt3_sending'
# comp<-'cd4_dpt1_sending'
# comp<-'cd4_dpt2_sending'
# comp<-'cd4_dpt3_sending'

# comp<-'growing_sending_cd8_dpt1'
# comp<-'growing_sending_cd8_dpt2'
# comp<-'growing_sending_cd8_dpt3'
# comp<-'growing_sending_cd4_dpt1'
# comp<-'growing_sending_cd4_dpt2'
# comp<-'growing_sending_cd4_dpt3'
# comp<-'shrinking_sending_cd8_dpt1'
# comp<-'shrinking_sending_cd8_dpt2'
# comp<-'shrinking_sending_cd8_dpt3'
# comp<-'shrinking_sending_cd4_dpt1'
# comp<-'shrinking_sending_cd4_dpt2'
comp<-'shrinking_sending_cd4_dpt3'

pdf(paste0(directory,'R_analysis/plots_gsea_treated_paired_',comp,'.pdf'),width = 15,height = 16)
gene_list_full<-inters_paired %>% 
  filter(fdr_ligand < 0.05 & ligandDE=='True' & cell_type_ligand=='shrinking_Tumor_clones' & cell_type_receptor=='CD4_T_dpt_3') #change ct!
gene_list <- gene_list_full$log2FC_ligand
names(gene_list) = gene_list_full$ligand
gene_list = sort(gene_list, decreasing = TRUE)
gene_list = gene_list[!duplicated(names(gene_list))]
head(gene_list)

if(length(gene_list) > 5){
  res1 <- GSEA(gene_list, GO_files[1], pval = 0.05, scoreType = 'std')$Plot
  res2 <- GSEA(gene_list, GO_files[2], pval = 0.05, scoreType = 'std')$Plot
  res3 <- GSEA(gene_list, GO_files[3], pval = 0.05, scoreType = 'std')$Plot
  res4 <- GSEA(gene_list, GO_files[4], pval = 0.05, scoreType = 'std')$Plot
  res5 <- GSEA(gene_list, GO_files[5], pval = 0.05, scoreType = 'std')$Plot
  res6 <- GSEA(gene_list, GO_files[6], pval = 0.05, scoreType = 'std')$Plot
  res7 <- GSEA(gene_list, GO_files[7], pval = 0.05, scoreType = 'std')$Plot
  res8 <- GSEA(gene_list, GO_files[8], pval = 0.05, scoreType = 'std')$Plot
  
  print((res1 +res2 +res3 +res4 )/(res5 +res6 +res7 +res8 )&labs(subtitle = paste0('Cells included ',comp)))
}
dev.off()

#### hypeR pathway analysis ####

directory<-'../data/cell_type_ct_dpt_032425/'

# use genes from CT visualization output
inters_paired<-read.delim2(paste0(directory,'circos_myeloid_t_CUIMC_only/circle_plot_tabular.tsv'))
inters_paired<-inters_paired %>% filter(type=='ligand')

inters_paired$DE<-ifelse(inters_paired$log2FC>0,'Up','Down')
inters_paired$log2FC_ligand<-as.numeric(inters_paired$log2FC)
inters_paired$fdr<-as.numeric(inters_paired$fdr)

inters_paired_DE_up<-inters_paired %>% filter(DE=='Up')
inters_paired_DE_down<-inters_paired %>% filter(DE=='Down')

# Run hypergeometric GSEA (select comp for appropriate output file name)

# comp<-'ct_viz'
comp<-'ct_viz_myeloid_immune'
# comp<-'ct_viz_hem'

HALLMARK <- msigdb_gsets(species="Homo sapiens", category="H")
KEGG     <- msigdb_gsets(species="Homo sapiens", category="C2", subcategory="CP:KEGG")
REACTOME <- msigdb_gsets(species="Homo sapiens", category="C2", subcategory="CP:REACTOME")
PID <-msigdb_gsets(species="Homo sapiens", category="C2", subcategory="CP:PID") # signaling events mediated by stem cell factor receptor
CP <- msigdb_gsets(species="Homo sapiens", category="C2", subcategory="CP")
C5_BP <- msigdb_gsets(species="Homo sapiens", category="C5", subcategory="BP") # biological process
C5_MF <- msigdb_gsets(species="Homo sapiens", category="C5", subcategory="MF") # molecular function
genesets_w <- as.list(as.data.frame(read.csv('genesets/wouters_mel_sigs.csv', na.strings = c('', NA))))

# run based on CT visualization output
background<-inters_paired$target
markers<-inters_paired %>%
  group_by(DE) %>%
  mutate(gene=target)

pdf(paste0(directory,'R_analysis/plots_hyper_treated_paired_',comp,'.pdf'),width = 15,height = 10)
for(m in unique(markers$DE)){
  subs<-markers %>% 
    dplyr::filter(DE==m) %>% 
    mutate(abs=abs(log2FC_ligand)) %>%
    arrange(desc(abs)) %>%
    top_n(100) %>%
    magrittr::use_series(gene)
  
  h1<-hyp_dots(hypeR(as.character(subs), genesets = genesets_w, background=background), 
               title = paste0('Wouters pathway selection\n',m, ' ',comp), abrv = 80) + 
    theme_bw() + 
    theme(axis.title.y = element_blank(), 
          axis.text.y = element_text(color = 'black',size=5), 
          axis.text.x = element_text(color = 'black'), 
          panel.grid.minor = element_blank(),
          title = element_text(size=8)) + 
    scale_color_gradient(high = 'black', low = '#C51B8A')
  
  h2<-hyp_dots(hypeR(as.character(subs), genesets = HALLMARK, background=background), 
               title = paste0('Hallmarks\n',m, ' ',comp), abrv = 80) + 
    theme_bw() + 
    theme(axis.title.y = element_blank(), 
          axis.text.y = element_text(color = 'black',size=5), 
          axis.text.x = element_text(color = 'black'), 
          panel.grid.minor = element_blank(),
          title = element_text(size=8)) + 
    scale_color_gradient(high = 'black', low = '#C51B8A')
  
  h3<-hyp_dots(hypeR(as.character(subs), genesets = KEGG, background=background), 
               title = paste0('KEGG\n',m, ' ',comp), abrv = 80) + 
    theme_bw() + 
    theme(axis.title.y = element_blank(), 
          axis.text.y = element_text(color = 'black',size=5), 
          axis.text.x = element_text(color = 'black'), 
          panel.grid.minor = element_blank(),
          title = element_text(size=8)) + 
    scale_color_gradient(high = 'black', low = '#C51B8A')
  
  h4<-hyp_dots(hypeR(as.character(subs), genesets = REACTOME, background=background), 
               title = paste0('REACTOME\n',m, ' ',comp), abrv = 80) + 
    theme_bw() + 
    theme(axis.title.y = element_blank(), 
          axis.text.y = element_text(color = 'black',size=5), 
          axis.text.x = element_text(color = 'black'), 
          panel.grid.minor = element_blank(),
          title = element_text(size=8)) + 
    scale_color_gradient(high = 'black', low = '#C51B8A')
  
  h5<-hyp_dots(hypeR(as.character(subs), genesets = C5_BP, background=background), 
               title = paste0('C5_BP\n',m, ' ',comp), abrv = 80) + 
    theme_bw() + 
    theme(axis.title.y = element_blank(), 
          axis.text.y = element_text(color = 'black',size=5), 
          axis.text.x = element_text(color = 'black'), 
          panel.grid.minor = element_blank(),
          title = element_text(size=8)) + 
    scale_color_gradient(high = 'black', low = '#C51B8A')
  
  h6<-hyp_dots(hypeR(as.character(subs), genesets = CP, background=background), 
               title = paste0('CP\n',m, ' ',comp), abrv = 80) + 
    theme_bw() + 
    theme(axis.title.y = element_blank(), 
          axis.text.y = element_text(color = 'black',size=5), 
          axis.text.x = element_text(color = 'black'), 
          panel.grid.minor = element_blank(),
          title = element_text(size=8)) + 
    scale_color_gradient(high = 'black', low = '#C51B8A')
  
  print((h1+h2+h3)/(h4+h5+h6))
  
  # export results to a separate table for other visualizations
  export <- rbind(h1$data, h2$data, h3$data, h4$data, h5$data, h6$data)
  export <- export[!duplicated(export$label), ]
  export$fdr <- p.adjust(export$pval, method = "fdr")
  write.csv(export, paste0(directory,'R_analysis/table_hyper_treated_paired_',comp,'_',m,'.csv'))
}
dev.off()

#### Visualize stacked hematopoietic GSEA ####
# comp<-'ct_viz'
comp<-'ct_viz_myeloid_immune'
# comp<-'ct_viz_hem'

export_up <- read.csv(paste0(directory,'R_analysis/table_hyper_treated_paired_',comp,'_Up.csv'))
export_up$fill <- -log10(export_up$fdr) / max(-log10(export_up$fdr))
export_dn <- read.csv(paste0(directory,'R_analysis/table_hyper_treated_paired_',comp,'_Down.csv'))
export_dn$fill <- -log10(export_dn$fdr) / max(-log10(export_dn$fdr))
export_up$sign <- 1
export_dn$sign <- -1
export <- rbind(export_up,export_dn)
export$x <- export$sign * -log10(export$fdr)
export$fill <- export$sign * export$fill
export <- export[export$fdr < 0.05,] # remove insignificant

keep <- read.csv(paste0(directory,'R_analysis/keep.csv'), header=FALSE)[[1]] # remove broad/general pathways
export <- export[export$label %in% keep,]

export <- export[order(export$pval, decreasing=F),] # sort ascending by pval, to remove duplicates
export <- export[!duplicated(export$label), ] # remove duplicate hits, based on which has greater signif.
export <- export[order(export$x, decreasing=F),]
export$label <- factor(export$label, levels=export$label)
rownames(export) <- export$label
# df <- rbind(head(export,10),tail(export,10)) # only top 10 up and top 10 down
df <- export # all rows
pdf(paste0(directory,'R_analysis/plots_hyper_treated_paired_',comp,'_pretty.pdf'),width = 10,height = 15)
g<-ggplot(df, aes(x=label, y=x, fill=fill)) + geom_col() + coord_flip() + scale_fill_gradient2(
  low = "blue", 
  mid = "white", 
  high = "red", 
  midpoint = 0
) + theme_bw() + geom_vline(xintercept = 0, color='black') + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
print(g)
dev.off()