#!/usr/bin/env Rscript

#### Cohort overview plots
#### Author: Jana Biermann, PhD

library(dplyr)
library(ggplot2)
library(gplots)
library(viridis)
library(RColorBrewer)

setwd('~/Documents/melanoma/')


drug_cols<-c( "#A65628",#adop
              "#49006A", #ipo
              "#7A0177","#AE017E",
              "#999999", #L
              "#E41A1C", #naive
              "#08519C", # nivo
              "#FF7F00", #combi
              "#238B45", #pembro
              "#41AB5D", "#74C476","#A1D99B","#C7E9C0",
              "#4292C6", #sp
              "#FFD92F", #tra
              "#F1B6DA") #vermu

#### clin overview v5 ####
df<-read.csv('data/clin_graph_v5.csv',na.strings = '')
df2 <- df %>% select(patient,Sample_J,Drug,survival_status,tissue,
                     best_response,OS,PFS,group,pair,
                     treatment_duration,treatment_start,treatment_drug)
df2<- df2 %>% 
  arrange(group) %>% 
  mutate(patient =ordered(patient, levels = rev(unique(patient))))

df_tmp<-df2 %>% select(patient) %>% unique
df_tmp$xmin<-rev(seq(0.6,nrow(df_tmp),1))
df_tmp$xmax<-rev(seq(1.4,nrow(df_tmp)+1,1))
df2<-left_join(df2,df_tmp,by='patient')

pdf('data/plot_patient_overview_v5.pdf',width = 10)
df2 %>%
  ggplot(aes(x=patient, y=Sample_J, color=Drug, shape=tissue)) + 
  geom_point(size=2.5) +
  scale_shape_manual(values=c(18,16,3,4,17))+
  scale_color_manual(values = drug_cols)+
  scale_fill_manual(values = drug_cols)+
  geom_point(data = df2, aes(x=patient, y=OS),col='black',shape='#',size=2) +
  geom_point(data = df2, aes(x=patient, y=PFS),col='black',shape='*',size=2.5) +
  geom_abline(intercept = 0, slope=0, linetype="dashed",color='black')+
  ylim(-800, 2500) +
  coord_flip() + 
  theme_bw() + 
  ylab('Days since ICI treatment')+
  geom_rect(aes(ymin = treatment_start, ymax = (treatment_start+treatment_duration), 
                xmin = xmin, xmax = xmax, fill = treatment_drug), 
            alpha = 0.4,lty=0)+
  geom_text(aes(label=best_response),size=2.5,color='black',y=2300)+
  geom_text(aes(label=survival_status),size=2.5,color='black',y=2500)+
  geom_text(aes(label=group),size=2,color='black',y=-700)+
  theme(legend.position = 'bottom',
        legend.key.size = unit(0.3, 'cm'),
        legend.text = element_text(size=7),
        legend.direction = 'vertical')
dev.off()

#### clinical overview ####
clin<-read.csv('data/clin.csv',na.strings = '')
pdf('data/plot_patient_clinical_overview.pdf')
ggplot(data=clin, aes(x=group, fill=sample)) +
  geom_bar(stat="count",color='black') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1))+
  ggtitle('Samples per group')

ggplot(data=clin, aes(x=group, fill=tissue)) +
  geom_bar(stat="count",color='black') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1))+
  ggtitle('Tissues per group')

ggplot(data=clin, aes(x=tissue, fill=sample)) +
  geom_bar(stat="count",color='black') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1))+
  ggtitle('Samples per tissue')

ggplot(data=clin, aes(x=group, fill=sample)) +
  geom_bar(stat="count",color='black') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1),
        legend.position = 'bottom')+
  ggtitle('Samples per group & tissue') +
  facet_grid(~tissue)

ggplot(data=clin, aes(x=Drug, fill=sample)) +
  geom_bar(stat="count",color='black') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1))+
  ggtitle('Samples per drug')

ggplot(data=clin, aes(x=time_drug, fill=sample)) +
  geom_bar(stat="count",color='black') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1))+
  ggtitle('Samples per time point & drug')
dev.off()


clin<-read.csv('data/clin.csv',na.strings = '')
clin$treated<-ifelse(clin$group %in% c('pre_naive','pre_naive_muc'),'Untreated','Treated')
clin$responder<-ifelse(clin$best_response %in% c('CR','PR'),'Responder','Non-responder')
clin$mucosal<-ifelse(clin$group %in% c('post_PD-1_muc','pre_naive_muc'),'Mucosal','Non-mucosal')
#write.csv(clin,'data/clin.csv',row.names = F)

pdf('data/plots_time_since_tx.pdf')
ggplot(data=clin, aes(x=reorder(sample,Tx_start_days), y=Tx_start_days,col=mucosal)) +
  geom_point() +
  theme_bw() +
  theme(legend.position = 'bottom')+
  geom_hline(yintercept = 0, linetype='dashed',col='darkgrey')+
  coord_flip()+
  ggtitle('Time from ICI start to sample')+
  facet_grid(treated~responder,scales = 'free_y')

ggplot(data=clin, aes(x=reorder(sample,Tx_start_days), y=Tx_start_days,col=mucosal)) +
  geom_point() +
  theme_bw() +
  theme(legend.position = 'bottom')+
  geom_hline(yintercept = 0, linetype='dashed',col='darkgrey')+
  coord_flip()+
  ggtitle('Time from ICI start to sample')+
  facet_grid(responder~.,scales = 'free_y')

ggplot(data=clin, aes(x=reorder(sample,Tx_start_days), y=Tx_start_days,col=mucosal)) +
  geom_bar(aes(y=treatment_duration,fill=time),stat = 'identity',alpha=0.5)+
  geom_point() +
  scale_color_manual(values = c('#E6AB02','black'))+
  theme_bw() +
  theme(legend.position = 'bottom')+
  coord_flip()+
  ggtitle('Time from ICI start to sample')+
  facet_grid(responder~.,scales = 'free_y')
dev.off()
