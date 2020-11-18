#library(ssh)
library(mnem)
library(Cairo)
setwd('/Volumes/Macintosh HD - Data/study/yr2sem2/rotation1/result/mnem')
source('/Volumes/Macintosh HD - Data/study/yr2sem2/rotation1/code/functions.R')
dir.create('combined mnem')

# load the dimension reduction results 
logodd.umap = readRDS('../dim reduction/logoddumap.rds')
logodd.tsne = readRDS('../dim reduction/logoddtsne.rds')
count.umap = readRDS('../dim reduction/countumap.rds')
count.tsne = readRDS('../dim reduction/counttsne.rds')

# analysis the mnem results for different parameters

# keep track of parameter sets that don't 
lack.list <- list()

# record the best log likelihood and k among all runs for each parameter setting
best.ll.k<- NULL

for (type in c('random','networks','cluster','cluster2','cluster3')){
  for (ksel1 in c('kmeans','hc')){
    for(ksel2 in c('silhouette','AIC','BIC')){
      # skip combinations when hc is used in ksel1, but silhouette is not used in ksel2
      if (!(ksel1 == 'hc' && ksel2 != 'silhouette')){
        for (ksel3 in c('cor','euclidean')){
          if (type %in% c('random','networks')){
            for (completeness in c(T,F)){
              # analyze individual run results with T or F
              parameters <- paste(type,completeness,ksel1,ksel2,ksel3,sep='_')
              res <- analyze.mnem (parameters,lack.list)
              best.ll.k <- rbind(best.ll.k,res[[1]])
              lack.list <- res[[2]]
            }
          }
          else{
            # analyze individual run results with only T 
            parameters <- paste(type,'TRUE',ksel1,ksel2,ksel3,sep='_')
            res <- analyze.mnem (parameters,lack.list)
            best.ll.k <- rbind(best.ll.k,res[[1]])
            lack.list <- res[[2]]
          }
        }
      }
    }
  }
}

colnames(best.ll.k) <- c('parameter','ll','k')

saveRDS(lack.list,file='lack.rds')
saveRDS(best.ll.k,file='bestllk.rds')


