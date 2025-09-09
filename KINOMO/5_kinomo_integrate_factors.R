#!/usr/bin/env Rscript

#### Title: Combination of NMF factors to metafactors
#### Author: Jana Biermann, PhD
#### Based on code by Somnath Tagore, PhD

library(dplyr)
library(Seurat)
library(ggplot2)
library(gplots)
library(hypeR)
library(viridis)
library(pheatmap)
library(scales)
library(RColorBrewer)
library(hypeR)
library(msigdbr)
library(maptree)
library(reshape2)
'%notin%' <- Negate('%in%')

directory<-'data/KINOMO/'
label<-'melanoma'

setwd('~/Documents/melanoma')

ranks_tab<-read.csv('data/KINOMO/kinomo_melanoma_samples.csv')

pats<-c('F01_on','F01_pre','F02_on','F02_pre','F03_post1_on2','F03_post1_pre2','F04_pre',
        'F05_pre','F06_post1_pre2','F07_post1_pre2','F08_post','F09_post',
        'F12_post','F12_pre','F15_pre','F16_post1_pre2','F16_pre','F17_post',
        'F18_post','F20_post1_pre2','F22_post','F23_post',
        'F26_post','F27_post1_pre2','F28_post1_pre2','F29_post1_pre2','F30_post',
        'F31_post','R204_pre','R294_on','R308_pre','R310_on1','R310_on2','R310_pre',
        'R319_on','R319_pre','R328_on','R329_on','R334_pre','R354_pre')


#### Functions ####

# Stouffer integration without weight
stouffer.integ <- apply(my_matric,1,function(x){
  y <- sum(x)/sqrt(length(x))
  #y <- sum(w*x)/sqrt(sum(w^2))
  return(y)
})

# Identify z-score
zscore_data <- t(apply(my_matric,1,function(x){
  y <- (x - mean(x))/sd(x)
  return(y)
}))

# weighted stouffer using z-score
current.stouffer <- t(apply(zscore_data,1,function(x,clust,weight){
  temp.data <- data.frame(nes = x, cluster = clust, weight = weight)
  temp.data$w_z <- temp.data$nes*temp.data$weight
  temp.data$w_w <- temp.data$weight*temp.data$weight
  temp.agg <- aggregate(.~cluster, temp.data,sum)
  temp.agg$res <- temp.agg$w_z/sqrt(temp.agg$w_w)
  temp.res <- temp.agg$res
  return(temp.res)
}
clust = annot.col$Cluster, weight = annot.col$SilScore))
colnames(current.stouffer) <- levels(annot.col$Cluster)

# Rank sort using stouffer weights
rank.stouffer <- apply(current.stouffer,2,rank)
keep.num <- 100
keep.index <- apply(rank.stouffer,1,function(x){
  y <- ((max(x) > (nrow(rank.stouffer) + 1 - keep.num)) | min(x) <= keep.num)
  return(y)
})
my_matric <- my_matric[keep.index,]

# Cluster = the clusters identified by kinomo co-correlation
# weight = the average SilScore or Ward distance (whatever you used)

entropy = function(cat.vect){
  px  = table(cat.vect)/length(cat.vect)
  lpx = log(px, base=2)
  ent = -sum(px*lpx)
  return(ent)
}

entropy(MP_factors$time[MP_factors$MP==8])


#### Apply to data #####
# Download gene lists + rank files
for(pat in pats){
  ranks<-as.numeric(unlist(ranks_tab[ranks_tab$sample==pat,2:3]))
  for(rank in ranks){
    system(paste0("aws s3 sync s3://melanoma-ribas/KINOMO/",pat,"/ data/KINOMO/",pat,"/ --exclude '*' --include '*KINOMO_nmf_rank_",rank,"_top100_W.csv' "))
    system(paste0("aws s3 sync s3://melanoma-ribas/KINOMO/",pat,"/ data/KINOMO/",pat,"/ --exclude '*' --include '*KINOMO_nmf_rank_",rank,"_top200_W.csv' "))
    system(paste0("aws s3 sync s3://melanoma-ribas/KINOMO/",pat,"/ data/KINOMO/",pat,"/ --exclude '*' --include '*KINOMO_nmf_rank_",rank,"_top50_W.csv' "))
    system(paste0("aws s3 sync s3://melanoma-ribas/KINOMO/",pat,"/ data/KINOMO/",pat,"/ --exclude '*' --include '*KINOMO_nmf_rank_",rank,"_W.csv' "))
    #system(paste0("aws s3 sync s3://melanoma-ribas/KINOMO/",pat,"/ data/KINOMO/",pat,"/ --exclude '*' --include '*KINOMO_nmf_rank_",rank,".rds' "))
  }
}


### Part 1: Combining all NMF factors from each sample
for(top_number in c(50,100,200)){
  for(ranklabel in c('best','secondbest')){
    ifelse(!dir.exists(file.path(paste0(directory,ranklabel,"/top",top_number))), 
           dir.create(file.path(paste0(directory,ranklabel,"/top",top_number)),recursive = T), FALSE)
    
    # Get top genes from all samples
    genes_combined<-c()
    for(pat in pats){
      # Select rank
      if(ranklabel=='best'){
        rank<-as.numeric(unlist(ranks_tab[ranks_tab$sample==pat,2]))
      }else{
        rank<-as.numeric(unlist(ranks_tab[ranks_tab$sample==pat,3]))
      }
      
      # Read-in top n genes per factor
      W_topn<-read.csv(paste0(directory,pat,"/",pat,"_KINOMO_nmf_rank_",rank,"_top",top_number,"_W.csv"),row.names = 1)
      
      genes_combined<-unique(c(genes_combined,unlist(W_topn)))
    }
    
    # Get Ws for top genes from all samples
    W_combined<-NULL
    for(pat in pats){
      # Select rank
      if(ranklabel=='best'){
        rank<-as.numeric(unlist(ranks_tab[ranks_tab$sample==pat,2]))
      }else{
        rank<-as.numeric(unlist(ranks_tab[ranks_tab$sample==pat,3]))
      }
      
      # Read-in top n genes per factor
      W_topn<-read.csv(paste0(directory,pat,"/",pat,"_KINOMO_nmf_rank_",rank,"_top",top_number,"_W.csv"),row.names = 1)
      
      # Read-in W file for rank
      W_rank<-read.csv(paste0(directory,pat,"/",pat,"_KINOMO_nmf_rank_",rank,"_W.csv"))
      colnames(W_rank)<-c('gene',colnames(W_topn))
      
      # Subset W to genes in top_number
      W_rank_sub<- W_rank[W_rank$gene %in% genes_combined,]  
      
      if(is.null(W_combined)==T){
        W_combined<-W_rank_sub
      }else{
        W_combined<-full_join(W_combined,W_rank_sub,by='gene')
      }
    }
    
    write.csv(W_combined,paste0(directory,ranklabel, '/top',top_number,
                                '/table_W_combined_',ranklabel,'_top',top_number,'.csv'),row.names = F)
  }  
}

### Part 2: Correlating Ws and selecting MPs
clin<-read.csv('data/clin.csv')
clin<- clin %>% select('sample','group','best_response','time')
ranks_tab<-read.csv('data/KINOMO/kinomo_melanoma_samples.csv')

for(top_number in c(50,100,200)){
  for(ranklabel in c('best','secondbest')){
    W_combined<-read.csv(paste0(directory,ranklabel, '/top',top_number,
                                '/table_W_combined_',ranklabel,'_top',top_number,'.csv'),row.names = 1)
    
    for(clustering_method in c('ward.D2','complete')){
      ifelse(!dir.exists(file.path(paste0(directory,ranklabel,"/top",top_number,'/',clustering_method))), 
             dir.create(file.path(paste0(directory,ranklabel,"/top",top_number,'/',clustering_method)),recursive = T), FALSE)
      # Correlation and heatmap annotation
      W_corr <- cor(W_combined, method = "pearson",use = "pairwise.complete.obs")
      write.csv(W_corr,paste0(directory,ranklabel,'/top',top_number,"/table_W_correlation_",ranklabel,'_top',top_number,".csv"),row.names = F)
      hc<-hclust(dist(1-W_corr),method = clustering_method)
      heat_labels<-data.frame(sample_factor=colnames(W_corr)[hc$order])
      heat_labels$sample<-substr(heat_labels$sample_factor,1,(nchar(heat_labels$sample_factor)-8))
      heat_labels<-left_join(heat_labels,clin,by='sample')
      rownames(heat_labels)<-heat_labels$sample_factor
      df<-heat_labels
      heat_labels$sample_factor<-NULL
      
      ann_colors = list(
        sample = setNames(c(sample(hue_pal()(length(unique(heat_labels$sample))))),
                          unique(heat_labels$sample)),
        group = setNames(c(hue_pal()(length(unique(heat_labels$group)))),
                         unique(heat_labels$group)),
        best_response = setNames(c(hue_pal()(length(unique(heat_labels$best_response)))),
                                 unique(heat_labels$best_response)),
        time = setNames(c(hue_pal()(length(unique(heat_labels$time)))),
                          unique(heat_labels$time)))
      
      # Optimal number of clusters
      op_k <- kgs(hclust(dist(1-W_corr),method = clustering_method), dist(1-W_corr), maxclus = 20)
      #meta_number<- as.integer(names(op_k[which(op_k == min(op_k))]))
      #meta_number<-7
      
      for(meta_number in c(5:12)){
        MP<-data.frame(MP=cutree(hc,k=meta_number))
        MP$MP<-as.character(MP$MP)
        MP$sample_factor<-rownames(MP)
        MP_factors<-left_join(df,MP,by='sample_factor')
        write.csv(MP_factors,paste0(directory,ranklabel, '/top',top_number,'/',clustering_method,
                            '/table_MP_factors_',ranklabel,'_top',top_number,"_",clustering_method,"_",meta_number,'MPs.csv'),row.names = F)
        
        
        ## Save metaprogram genes to table
        MP_list<-c()
        for(m in 1:meta_number){
          df_sub<-MP_factors[MP_factors$MP==m,]
          MP_list_tmp<-c()
          
          for(pat in unique(df_sub$sample)){
            # Select rank
            if(ranklabel=='best'){
              rank<-as.numeric(unlist(ranks_tab[ranks_tab$sample==pat,2]))
            }else{
              rank<-as.numeric(unlist(ranks_tab[ranks_tab$sample==pat,3]))
            }
            
            # Read-in top n genes per factor
            W_topn<-read.csv(paste0(directory,pat,"/",pat,"_KINOMO_nmf_rank_",rank,"_top",top_number,"_W.csv"),row.names = 1)
            W_topn_sub<-W_topn[,intersect(df_sub$sample_factor,colnames(W_topn))]
            MP_list_tmp<-c(MP_list_tmp,unique(unlist(W_topn_sub)))
          }
          
          MP_list_tmp<-unique(MP_list_tmp)
          MP_list[[m]]<-MP_list_tmp
          names(MP_list)[m]<-paste0('MP',m)
        }
        MP_df<-sapply(MP_list, function(x){
          c(x, rep(NA, max(sapply(MP_list,length)) - length(x)))})
        write.csv(MP_df,paste0(directory,ranklabel, '/top',top_number,'/',clustering_method,
                               "/table_metaprograms_",ranklabel,"_top",top_number,"_",clustering_method,"_",meta_number,"MPs.csv"),row.names = F)
        
        # Get unique and shared genes across MPs
        blacklist<-rownames(W_combined)[lapply(rownames(W_combined),
                                               function(x){grep(x,MP_df) %>% length()})==ncol(MP_df)]
        whitelist<-rownames(W_combined)[lapply(rownames(W_combined),
                                               function(x){grep(x,MP_df) %>% length()})==1]
        MP_unique<-data.frame(MP=colnames(MP_df))
        MP_unique$unique<-apply(MP_df,2,
                                function(x){length(intersect(x,whitelist))})
        MP_unique$shared<-apply(MP_df,2,
                               function(x){length(na.omit(x))-length(intersect(x,whitelist))})
        MP_unique<-melt(MP_unique)
        
        ### Plots
        pdf(paste0(directory,ranklabel, '/top',top_number,'/',clustering_method,
                   '/plots_W_correlation_',ranklabel,'_top',top_number,"_",clustering_method,"_",meta_number,'MPs.pdf'),width = 10,height = 10)
        # Heatmap
        print(pheatmap(W_corr,
                       color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100),
                       breaks = seq(-1, 1, length.out = 100),
                       cutree_rows = meta_number,
                       cutree_cols = meta_number,
                       clustering_distance_rows = dist(1-W_corr),
                       clustering_distance_cols = dist(1-W_corr),
                       clustering_method = clustering_method,
                       annotation_col = heat_labels,
                       annotation_colors = ann_colors,
                       fontsize = 7,
                       fontsize_row = 4, 
                       fontsize_col = 4))
        
        # Dendrogram
        print(plot(hc,hang=-1,cex=0.3))
        print(rect.hclust(hc, border = "purple",k = meta_number))
        
        # KGS
        print(plot(names(op_k), op_k, xlab="# clusters", ylab="penalty",
                   main='KGS Measure for Pruning Hierarchical Clusters'))
        
        print(ggplot(MP_unique,aes(x=MP,y=value,fill=variable))+
                geom_bar(stat="identity")+
                theme_classic()+
                ggtitle('Unique genes per MP')+
                ylab('Number of unique genes per MP'))
        
        # Bar plots
        p1<-ggplot(MP_factors,aes(x=MP,fill=group))+
          geom_bar(stat="count",col='black')+
          theme_classic()+
          ylab('Number of factors')
      
        p2<-ggplot(MP_factors,aes(x=MP,fill=sample))+
          geom_bar(stat="count",col='black')+
          theme_classic()+
          ylab('Number of factors')+
          theme(legend.text = element_text(size=4),
                legend.key.size = unit(0.3,'cm'))
      
        p3<-ggplot(MP_factors,aes(x=MP,fill=best_response))+
          geom_bar(stat="count",col='black')+
          theme_classic()+
          ylab('Number of factors')
      
        p4<-ggplot(MP_factors,aes(x=MP,fill=time))+
          geom_bar(stat="count",col='black')+
          theme_classic()+
          ylab('Number of factors')
        print((p1 + p2) / (p3 +p4))
        
        p5<-ggplot(MP_factors,aes(x=MP,fill=group))+
          geom_bar(position="fill", stat="count",col='black')+
          theme_classic()+
          ylab('Percentage of factors')
      
        p6<-ggplot(MP_factors,aes(x=MP,fill=sample))+
          geom_bar(position="fill", stat="count",col='black')+
          theme_classic()+
          ylab('Percentage of factors')
        theme(legend.text = element_text(size=4),
                legend.key.size = unit(0.3,'cm'))
      
        p7<-ggplot(MP_factors,aes(x=MP,fill=best_response))+
          geom_bar(position="fill", stat="count",col='black')+
          theme_classic()+
          ylab('Percentage of factors')
      
        p8<-ggplot(MP_factors,aes(x=MP,fill=time))+
          geom_bar(position="fill", stat="count",col='black')+
          theme_classic()+
          ylab('Percentage of factors')
        print((p5 + p6) / (p7 +p8))
        
        dev.off()
        
      }
    }
  }
}

# Reduce genes in MPs
# for(top_number in c(50,100,200)){
#   for(ranklabel in c('best','secondbest')){
#     for(clustering_method in c('ward.D2','complete')){
#       for(meta_number in c(5:12)){
        
for(top_number in c(100)){
  for(ranklabel in c('best')){
    for(clustering_method in c('ward.D2')){
      for(meta_number in c(11)){

        MP_df<-read.csv(paste0(directory,ranklabel, '/top',top_number,'/',clustering_method,
                               "/table_metaprograms_",ranklabel,"_top",top_number,"_",clustering_method,"_",meta_number,"MPs.csv"))
        W_combined<-read.csv(paste0(directory,ranklabel, '/top',top_number,
                                    '/table_W_combined_',ranklabel,'_top',top_number,'.csv'),row.names = 1)
        MP_factors<-read.csv(paste0(directory,ranklabel, '/top',top_number,'/',clustering_method,
                                    '/table_MP_factors_',ranklabel,'_top',top_number,"_",clustering_method,"_",meta_number,'MPs.csv'))
        MP_list_top100<-c()
        # blacklist<-rownames(W_combined)[lapply(rownames(W_combined),
        #                                        function(x){grep(x,MP_df) %>% length()})==ncol(MP_df)]
        whitelist<-rownames(W_combined)[lapply(rownames(W_combined),
                                               function(x){grep(x,MP_df) %>% length()})==1]
        
        for(m in 1:ncol(MP_df)){
          W_sub<- as.data.frame(W_combined[na.omit(MP_df[,m]),MP_factors[MP_factors$MP==m,1]])
          #W_sub<-W_sub[which(rownames(W_sub) %notin% blacklist),]
          #W_sub<-W_sub[which(rownames(W_sub) %in% whitelist),]
          
          # Stouffer integration without weight
          stouffer.integ <- apply(W_sub,1,function(x){
            y <- sum(x)/sqrt(length(x))
            #y <- sum(w*x)/sqrt(sum(w^2))
            return(y)
          })
          
          if(is.null(names(sort(stouffer.integ,decreasing = T)[1:100]))){
            MP_list_top100[[m]]<- rep(NA,100)
          }else{
            MP_list_top100[[m]]<- names(sort(stouffer.integ,decreasing = T)[1:100])
          }
        } 
        MP_top100_df<-sapply(MP_list_top100, function(x){
          c(x, rep(NA, max(sapply(MP_list_top100,length)) - length(x)))})
        write.csv(MP_top100_df,paste0(directory,ranklabel,'/top',top_number,'/',clustering_method,
                                      '/table_metaprograms_',ranklabel,
                                      '_top',top_number,"_",
                                      clustering_method,"_",
                                      meta_number,'MPs_stouffer.csv'),row.names = F)
        
      }
    }
  }
}

### Part 3: Pathway analysis of MPs
# Read-in genesets from hypeR
HALLMARK <- msigdb_gsets(species = 'Homo sapiens', category = 'H')
KEGG <- msigdb_gsets(species = 'Homo sapiens', category = 'C2', subcategory = 'CP:KEGG')
REACTOME <- msigdb_gsets(species = 'Homo sapiens', category = 'C2', subcategory = 'CP:REACTOME')
CP <- msigdb_gsets(species = 'Homo sapiens', category = 'C2', subcategory = 'CP')
C5_BP <- msigdb_gsets(species = 'Homo sapiens', category = 'C5', subcategory = 'BP')
genesets_w <- as.list(as.data.frame(read.csv('~/brain_mets/signatures/wouters_mel_sigs.csv', na.strings = c('', NA))))

# for(top_number in c(50,100,200)){
#   for(ranklabel in c('best','secondbest')){
#     for(clustering_method in c('ward.D2','complete')){
#       for(meta_number in c(5:12)){

for(top_number in c(100)){
  for(ranklabel in c('best')){
    for(clustering_method in c('ward.D2')){
      for(meta_number in c(11)){
        MP_top100_df<-read.csv(paste0(directory,ranklabel,'/top',top_number,'/',clustering_method,
                            '/table_metaprograms_',ranklabel,
                            '_top',top_number,"_",
                            clustering_method,"_",
                            meta_number,'MPs_stouffer.csv'))
          
      pdf(paste0(directory,ranklabel, '/top',top_number,'/',clustering_method,
                 '/plots_MP_pathways_',ranklabel,"_top",top_number,"_",clustering_method,"_",meta_number,
                 'MPs.pdf'),
          width = 15,height = 10)
      for(m in 1:ncol(MP_top100_df)){
        h1<-hyp_dots(hypeR(as.character(MP_top100_df[,m]), genesets = genesets_w), 
                     title = paste0('Wouters pathway selection\n',ranklabel,"; top",top_number,"; ",meta_number,' MPs (total)\nMP ',m), abrv = 80) + 
          theme_bw() + 
          theme(axis.title.y = element_blank(), 
                axis.text.y = element_text(color = 'black',size=5), 
                axis.text.x = element_text(color = 'black'), 
                panel.grid.minor = element_blank(),
                title = element_text(size=8)) + 
          scale_color_gradient(high = 'black', low = '#C51B8A')
        
        h2<-hyp_dots(hypeR(as.character(MP_top100_df[,m]), genesets = HALLMARK), 
                     title = paste0('Hallmarks\n',ranklabel,"; top",top_number,"; ",meta_number,' MPs (total)\nMP ',m), abrv = 80) + 
          theme_bw() + 
          theme(axis.title.y = element_blank(), 
                axis.text.y = element_text(color = 'black',size=5), 
                axis.text.x = element_text(color = 'black'), 
                panel.grid.minor = element_blank(),
                title = element_text(size=8)) + 
          scale_color_gradient(high = 'black', low = '#C51B8A')
        
        h3<-hyp_dots(hypeR(as.character(MP_top100_df[,m]), genesets = KEGG), 
                     title = paste0('KEGG\n',ranklabel,'\nMP ',m), abrv = 80) + 
          theme_bw() + 
          theme(axis.title.y = element_blank(), 
                axis.text.y = element_text(color = 'black',size=5), 
                axis.text.x = element_text(color = 'black'), 
                panel.grid.minor = element_blank(),
                title = element_text(size=8)) + 
          scale_color_gradient(high = 'black', low = '#C51B8A')
        
        h4<-hyp_dots(hypeR(as.character(MP_top100_df[,m]), genesets = REACTOME), 
                     title = paste0('REACTOME\n',ranklabel,"; top",top_number,"; ",meta_number,' MPs (total)\nMP ',m), abrv = 80) + 
          theme_bw() + 
          theme(axis.title.y = element_blank(), 
                axis.text.y = element_text(color = 'black',size=5), 
                axis.text.x = element_text(color = 'black'), 
                panel.grid.minor = element_blank(),
                title = element_text(size=8)) + 
          scale_color_gradient(high = 'black', low = '#C51B8A')
        
        h5<-hyp_dots(hypeR(as.character(MP_top100_df[,m]), genesets = C5_BP), 
                     title = paste0('C5_BP\n',ranklabel,'\nMP ',m), abrv = 80) + 
          theme_bw() + 
          theme(axis.title.y = element_blank(), 
                axis.text.y = element_text(color = 'black',size=5), 
                axis.text.x = element_text(color = 'black'), 
                panel.grid.minor = element_blank(),
                title = element_text(size=8)) + 
          scale_color_gradient(high = 'black', low = '#C51B8A')
        
        h6<-hyp_dots(hypeR(as.character(MP_top100_df[,m]), genesets = CP), 
                     title = paste0('CP\n',ranklabel,"; top",top_number,"; ",meta_number,' MPs (total)\nMP ',m), abrv = 80) + 
          theme_bw() + 
          theme(axis.title.y = element_blank(), 
                axis.text.y = element_text(color = 'black',size=5), 
                axis.text.x = element_text(color = 'black'), 
                panel.grid.minor = element_blank(),
                title = element_text(size=8)) + 
          scale_color_gradient(high = 'black', low = '#C51B8A')
        
        print((h1+h2+h3)/(h4+h5+h6))
      }
      dev.off()
      }
    }
  }
}

