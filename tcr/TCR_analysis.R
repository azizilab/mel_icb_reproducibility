#!/usr/bin/env Rscript

### title: TCR Analysis
### author: Jana Biermann, PhD

print(Sys.time())

library(Seurat)
library(dplyr)
library(ggplot2)
library(gplots)
library(viridis)
library(scales)
library(stringr)
library(reshape2)
library(patchwork)
library(RColorBrewer)
'%notin%' <- Negate('%in%')


#setwd('~/Documents/melanoma/')
label<-'melanoma'
celltype<-'tcells'
filename<-paste0(label,'_',celltype)
path.ct <- paste0('data/annotation/',celltype,'/')
directory<-'data/tcells/'
ifelse(!dir.exists(file.path(directory)), dir.create(file.path(directory),recursive = T), FALSE)

colExp<-c('#C82F71','#7473A6')
colSha<-c('#a67473','#2fc886','#2f71c8')
colMAIT<-c('#009F3F','#E0C126','#F08080','#2fbec8')
colExpPub<- c('#d85a90','#a8285f','#807fae','#504f7c')


seu<-readRDS(paste0(path.ct, 'data_', filename, '_v2.rds'))
seu<-subset(seu, cell_type_fine != 'NK cells')


#### Plot TCR data ####

### Count TCR overlap in CD8+
seu_both<-subset(seu,both_chains =='both' & cell_type_fine %in% c('CD8+ T cells TCF7+','CD8+ T cells TOX+'))
res <- data.frame(sample = unique(seu_both$sample))
for (pat in unique(seu_both$sample)) {
  if(table(seu_both$sample,seu_both$cell_type_fine)[pat,1]>0 & 
     table(seu_both$sample,seu_both$cell_type_fine)[pat,2]>0){
    res$shared_TCR[res$sample == pat] <- 
      length(intersect(subset(seu_both, sample == pat & 
                                cell_type_fine == 'CD8+ T cells TCF7+')$cdr3s_aa, 
                       subset(seu_both, sample == pat & 
                                cell_type_fine == 'CD8+ T cells TOX+')$cdr3s_aa))
    res$unique_TOX[res$sample == pat] <- 
      length(unique(subset(seu_both, sample == pat & 
                             cell_type_fine == 'CD8+ T cells TOX+')$cdr3s_aa)) - res$shared_TCR[res$sample == pat]
    res$unique_TCF7[res$sample == pat] <- 
      length(unique(subset(seu_both, sample == pat & 
                             cell_type_fine == 'CD8+ T cells TCF7+')$cdr3s_aa)) - res$shared_TCR[res$sample == pat]
  }
}

# Plots
df <- seu@meta.data[, c('barcode_pat', 'patient','sample', 'frequency', 'clone_size', 'cdr3s_aa', 'cell_type_fine', 
                        'mait', 'mait_evidence', 'inkt', 'inkt_evidence', 'both_chains', 'mait_inkt','tissue','group','sample_group','time')]
df_both <- df %>% filter(both_chains == 'both')

pdf(file = paste0(directory,'plots_',filename,'_TCR.pdf'))
# Stacked bar plot both_chains
ggplot(df, aes(x = sample, fill = both_chains)) + 
  geom_bar(colour = 'black', position = 'fill', size = 0.25) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  ggtitle('TCR chains recovered across samples') + 
  ylab('Fraction')

# Stacked bar plot mait_inkt
ggplot(df_both, aes(x = sample, fill = mait_inkt)) + 
  geom_bar(colour = 'black', position = 'fill', size = 0.25) + 
  ggtitle('MAIT + iNKT in T cells') + 
  theme_classic() + 
  xlab('') + 
  ylab('Fraction of T cells') + 
  scale_fill_manual(values = colMAIT)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Stacked bar plot clone_size
df_clone_size <- df_both %>% group_by(sample,clone_size) %>% summarise(n = n())
ggplot(df_clone_size, aes(x = sample,y=n, fill = clone_size)) + 
  geom_bar(stat="identity",position=position_dodge(), colour = 'black',size = 0.25) + theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  xlab('') + ylab('# T cells') + scale_fill_manual(values = colExp)

# Stacked bar plot clone_size CD8
df_clone_size_cd8 <- df_both %>% filter(grepl('CD8',cell_type_fine)) %>% group_by(sample,clone_size) %>% summarise(n = n())
ggplot(df_clone_size_cd8, aes(x = sample,y=n, fill = clone_size)) + 
  geom_bar(stat="identity",position=position_dodge(), colour = 'black',size = 0.25) + theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  xlab('') + ylab('# CD8+ T cells') + scale_fill_manual(values = colExp)

# Bar plot clones
df_unique <- df_both %>% group_by(sample) %>% summarise(nClones = n_distinct(cdr3s_aa))
p1<-ggplot(df_unique, aes(x = sample, y = nClones)) + 
  geom_bar(stat = 'identity') + 
  theme_classic() + 
  theme(axis.text.x = element_blank(),axis.title.x = element_blank()) + 
  ylab('# Clones')

# Bar plot cells
df_cells <- df_both %>% group_by(sample) %>% summarise(nCells = n_distinct(barcode_pat))
p2<- ggplot(df_cells, aes(x = sample, y = nCells)) + 
  geom_bar(stat = 'identity') + theme_classic() + 
  theme(axis.text.x = element_blank(),axis.title.x = element_blank()) + 
  ylab('# Cells')

# Stacked bar plot clone_size
p3<-ggplot(df_both, aes(x = sample, fill = clone_size)) + 
  geom_bar(colour = 'black', position = 'fill', size = 0.25) + 
  theme_classic() + scale_fill_manual(values = colExp)+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  ylab('Fraction of T cells')

p1/p2/p3

# Circle stacked bar plot stratified clonotype + tissue
df_ct_strat <- df_both %>%
  group_by(tissue, clone_size, cell_type_fine) %>%
  summarise(n = n()) %>%
  mutate(freq = n/sum(n))

ggplot(df_ct_strat, aes(x = '', y = freq, group = cell_type_fine, 
                        fill = cell_type_fine, color = cell_type_fine)) + 
  geom_bar(width = 1, stat = 'identity', color = 'white', size = 0.1) + 
  coord_polar('y', start = 0) + 
  theme_void() + 
  facet_grid(clone_size ~ tissue)+
  theme(legend.position = 'bottom',
        strip.text.x = element_text(angle = 90))

# Circle stacked bar plot stratified clonotype + group
df_ct_strat <- df_both %>%
  group_by(group, clone_size, cell_type_fine) %>%
  summarise(n = n()) %>%
  mutate(freq = n/sum(n))

ggplot(df_ct_strat, aes(x = '', y = freq, group = cell_type_fine, 
                        fill = cell_type_fine, color = cell_type_fine)) + 
  geom_bar(width = 1, stat = 'identity', color = 'white', size = 0.1) + 
  coord_polar('y', start = 0) + 
  theme_void() + 
  facet_grid(clone_size ~ group)+
  theme(legend.position = 'bottom',
        strip.text.x = element_text(angle = 90))

# Circle stacked bar plot stratified clonotype + time
df_ct_strat <- df_both %>%
  group_by(time, clone_size, cell_type_fine) %>%
  summarise(n = n()) %>%
  mutate(freq = n/sum(n))

ggplot(df_ct_strat, aes(x = '', y = freq, group = cell_type_fine, 
                        fill = cell_type_fine, color = cell_type_fine)) + 
  geom_bar(width = 1, stat = 'identity', color = 'white', size = 0.1) + 
  coord_polar('y', start = 0) + 
  theme_void() + 
  facet_grid(clone_size ~ time)+
  theme(legend.position = 'bottom',
        strip.text.x = element_text(angle = 90))

# Circle stacked bar plot stratified cell type + tissue
df_clonesize_strat <- df_both %>%
  group_by(tissue, cell_type_fine, clone_size) %>%
  summarise(n = n()) %>%
  mutate(freq = n/sum(n))

ggplot(df_clonesize_strat, aes(x = '', y = freq, group = clone_size, fill = clone_size, color = clone_size)) + 
  geom_bar(width = 1, stat = 'identity', color = 'white', size = 0.1) + 
  coord_polar('y', start = 0) + 
  theme_void() + scale_fill_manual(values = colExp)+
  facet_grid(cell_type_fine ~ tissue)+
  theme(legend.position = 'bottom',
        strip.text.x = element_text(angle = 90))

# Circle stacked bar plot stratified cell type + group
df_clonesize_strat <- df_both %>%
  group_by(group, cell_type_fine, clone_size) %>%
  summarise(n = n()) %>%
  mutate(freq = n/sum(n))

ggplot(df_clonesize_strat, aes(x = '', y = freq, group = clone_size, fill = clone_size, color = clone_size)) + 
  geom_bar(width = 1, stat = 'identity', color = 'white', size = 0.1) + 
  coord_polar('y', start = 0) + 
  theme_void() + scale_fill_manual(values = colExp)+
  facet_grid(cell_type_fine ~ group)+
  theme(legend.position = 'bottom',
        strip.text.x = element_text(angle = 90))

# Circle stacked bar plot stratified cell type + time
df_clonesize_strat <- df_both %>%
  group_by(time, cell_type_fine, clone_size) %>%
  summarise(n = n()) %>%
  mutate(freq = n/sum(n))

ggplot(df_clonesize_strat, aes(x = '', y = freq, group = clone_size, fill = clone_size, color = clone_size)) + 
  geom_bar(width = 1, stat = 'identity', color = 'white', size = 0.1) + 
  coord_polar('y', start = 0) + 
  theme_void() + scale_fill_manual(values = colExp)+
  facet_grid(cell_type_fine ~ time)+
  theme(legend.position = 'bottom',
        strip.text.x = element_text(angle = 90))

# Circle stacked bar plot stratified tissue + group
df_clonesize_strat <- df_both %>%
  group_by(group, tissue, clone_size) %>%
  summarise(n = n()) %>%
  mutate(freq = n/sum(n))

ggplot(df_clonesize_strat, aes(x = '', y = freq, group = clone_size, fill = clone_size, color = clone_size)) + 
  geom_bar(width = 1, stat = 'identity', color = 'white', size = 0.1) + 
  coord_polar('y', start = 0) + 
  theme_void() + scale_fill_manual(values = colExp)+
  facet_grid(tissue ~ group)+
  theme(legend.position = 'bottom',
        strip.text.x = element_text(angle = 90))

# Stacked bar plot TCR overlap
df_overlap <- melt(res) %>% group_by(sample, variable)
ggplot(df_overlap, aes(x = sample, y = value, fill = variable)) + 
  geom_bar(colour = 'black', position = 'fill', stat = 'identity') + theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  ggtitle('Distribution of TCRs between CD8+ T cells') + 
  xlab('sample') + 
  ylab('Fraction of CD8+ T cells')+
  scale_fill_manual(values = colSha)
dev.off()



##### updated plots TCR #####
#system("aws s3 sync s3://melanoma-ribas/Seurat/annotation/tcells/ data/annotation/tcells/ --exclude '*' --include 'data_melanoma_tcells_v4.rds' ")
seu <- readRDS('data/annotation/tcells/data_melanoma_tcells_v4.rds')
table(seu$cell_type_fine,useNA='always')
table(seu$tcr,seu$sample)

#system("aws s3 sync s3://melanoma-ribas/Seurat/ data/ --exclude '*' --include 'clin.csv' --quiet")
clin<-read.csv('data/clin.csv',na.strings = '')
clin$PFS<-as.numeric(clin$PFS)
seu$barcode_all<-rownames(seu@meta.data)
seu@meta.data <- seu@meta.data[ , !names(seu@meta.data) %in% names(clin)[2:ncol(clin)]]
seu@meta.data<-left_join(seu@meta.data,clin,by='sample')
rownames(seu@meta.data)<-seu$barcode_all
seu$treat_resp<-paste0(seu$treated,'_',seu$responder)
seu$combination<-ifelse(seu$treated=='Untreated', 'Untreated',
                        ifelse(seu$group %in% c('on_combination','post_combination'),
                               'treated_combination','treated_individual'))

df <- seu@meta.data[, c('barcode_pat', 'patient','sample', 'frequency', 'clone_size', 'cdr3s_aa', 'cell_type_fine', 
                        'mait', 'mait_evidence', 'inkt', 'inkt_evidence', 'both_chains', 'mait_inkt',
                        'treated', 'responder','Tx_phase','tissue','Gender',
                        'genotype_simplified','best_response','mucosal','treat_resp','combination')]

df_both <- df %>% 
  filter(both_chains == 'both')

df_both2 <- df %>% 
  filter(both_chains == 'both' & cell_type_fine %notin% c('NK_1','NK_2'))


pdf(file = paste0(directory,'plots_',filename,'_TCR_update.pdf'))

# Circle stacked bar plot ct
df_ct_strat <- df_both %>%
  group_by(Gender, treat_resp, cell_type_fine) %>%
  summarise(n = n()) %>%
  mutate(freq = n/sum(n))

ggplot(df_ct_strat, aes(x = '', y = freq, group = cell_type_fine, 
                        fill = cell_type_fine, color = cell_type_fine)) + 
  geom_bar(width = 1, stat = 'identity', color = 'white', size = 0.1) + 
  coord_polar('y', start = 0) + 
  theme_void() + 
  facet_grid(treat_resp ~ Gender)+
  theme(legend.position = 'bottom',
        strip.text.x = element_text(angle = 90))

ggplot(df_both, aes(x = sample,  group = cell_type_fine, 
                     fill = cell_type_fine, color = cell_type_fine)) + 
  geom_bar(stat = 'count', position = 'fill',color = 'white', size = 0.1) + 
  theme_classic()+
  facet_grid(treat_resp ~ Gender,scales = 'free_x')+
  theme(legend.position = 'bottom',
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        strip.text.y = element_text(angle = 0))

df_ct_strat <- df_both2 %>%
  group_by(clone_size, treat_resp, cell_type_fine) %>%
  summarise(n = n()) %>%
  mutate(freq = n/sum(n))
ggplot(df_ct_strat, aes(x = '', y = freq, group = cell_type_fine, 
                        fill = cell_type_fine, color = cell_type_fine)) + 
  geom_bar(width = 1, stat = 'identity', color = 'white', size = 0.1) + 
  coord_polar('y', start = 0) + 
  theme_void() + 
  facet_grid(treat_resp ~ clone_size)+
  theme(legend.position = 'bottom',
        strip.text.x = element_text(angle = 90))

df_ct_strat <- df_both2 %>%
  group_by(clone_size, Gender, cell_type_fine) %>%
  summarise(n = n()) %>%
  mutate(freq = n/sum(n))
ggplot(df_ct_strat, aes(x = '', y = freq, group = cell_type_fine, 
                        fill = cell_type_fine, color = cell_type_fine)) + 
  geom_bar(width = 1, stat = 'identity', color = 'white', size = 0.1) + 
  coord_polar('y', start = 0) + 
  theme_void() + 
  facet_grid(Gender ~ clone_size)+
  theme(legend.position = 'bottom',
        strip.text.x = element_text(angle = 90))

# Circle stacked bar plot clone_size
df_ct_strat <- df_both2 %>%
  group_by(Gender, treat_resp, clone_size) %>%
  summarise(n = n()) %>%
  mutate(freq = n/sum(n))

ggplot(df_ct_strat, aes(x = '', y = freq, group = clone_size, 
                        fill = clone_size, color = clone_size)) + 
  geom_bar(width = 1, stat = 'identity', color = 'white', size = 0.1) + 
  coord_polar('y', start = 0) + 
  theme_void() + 
  scale_fill_manual(values = colExp)+
  facet_grid(treat_resp ~ Gender)+
  theme(legend.position = 'bottom',
        strip.text.x = element_text(angle = 90))


ggplot(df_both2, aes(x = sample,  group = clone_size, 
                        fill = clone_size, color = clone_size)) + 
  geom_bar(stat = 'count', position = 'fill',color = 'white', size = 0.1) + 
  scale_fill_manual(values = colExp)+
  theme_classic()+
  facet_grid(treat_resp ~ Gender,scales = 'free_x')+
  theme(legend.position = 'bottom',
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        strip.text.x = element_text(angle = 90))

ggplot(df_both2, aes(x = Gender,  group = clone_size, 
                     fill = clone_size, color = clone_size)) + 
  geom_bar(stat = 'count', position = 'fill',color = 'white', size = 0.1) + 
  scale_fill_manual(values = colExp)+
  theme_classic()+
  facet_grid(~treat_resp,scales = 'free_x')+
  theme(legend.position = 'bottom',
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        strip.text.x = element_text(angle = 90))
dev.off()
