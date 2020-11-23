library(mnem)
library(Cairo)
library(ggplot2)
library(ggfortify)
setwd('/Volumes/Macintosh HD - Data/study/yr2sem2/rotation1/local_result')
dir.create('data')

# combine the mnem results for different parameters
for (type in c('random','networks','cluster','cluster2','cluster3')){
  for (ksel1 in c('kmeans','hc')){
    for(ksel2 in c('silhouette','AIC','BIC')){
      # skip combinations when hc is used in ksel1, but silhouette is not used in ksel2
      if (!(ksel1 == 'hc' && ksel2 != 'silhouette')){
        for (ksel3 in c('cor','euclidean')){
          # analyze individual run results with T or F
          parameters <- paste(type,completeness,ksel1,ksel2,ksel3,sep='_')
          files <- list.files(pattern = parameters)
          mnem <- combine.mnem (parameters,files)
          saveRDS(mnem,file = paste0('./data/',parameter,'.rds'))
        }
      }
    }
  }
}

# record the best log likelihood and k among all runs for each parameter setting
result.summary <- make.summary('./mnem/data')
saveRDS(result.summary,file='summary.rds')

# draw graph to compare various type, k and complete
plot.summary(result.summary, path = './graphs')


