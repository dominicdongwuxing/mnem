library(umap)
library(Rtsne)
library(ggfortify)
library(ggplot2)
library(gridExtra)
library(mnem)

setwd('/Volumes/Macintosh HD - Data/study/yr2sem2/rotation1/result')
dir.create('./logodd umap')

R <- as.matrix(read.table('R_kcde.csv'), sep = '\t')

# explore how different min_dist affect the result 
for (d in c(0.0001,0.001,0.01,0.2,0.5)){
  custom.config <- umap.defaults
  custom.config$random_state <- 42
  custom.config$min_dist <- d
  result <- umap(t(R), config = custom.config)
  saveRDS(result, file=paste0('./logodd umap/min_dist_',d,'.rds'))
  graph.umap (result, readRDS('bestk.rds')$best, paste0('UMAP min_dist = ', umap$config$min_dist))
}



dir.create('./dim reduction/manhattan')

# try using manhattan distance
custom.config <- umap.defaults
custom.config$random_state <- 42
custom.config$metric <- 'manhattan'
result <- umap(t(R), config = custom.config)
saveRDS(result, file='./dim reduction/manhattan/logodd.umap.manhattan.rds')
graph.umap (result, readRDS('bestk.rds')$best, "UMAP Manhattan distance")


# use manhattan distance on tSNE
logodd.tsne.manhattan <- Rtsne(dist(t(R), method = 'manhattan'), is_distance = T)
saveRDS(logodd.tsne.manhattan, file='./dim reduction/manhattan/logodd.tsne.manhattan.rds')
graph.tSNE(logodd.tsne.manhattan,readRDS('bestk.rds')$best,"tSNE Manhattan distance")

# now compare PCA, UMAP and tSNE results
graph.dim.reduction.loggodd (
  result = readRDS('bestk.rds'),
  logodd.pca = readRDS('./dim reduction/logoddpca.rds'), 
  logodd.umap = readRDS('./dim reduction/manhattan/logodd.umap.manhattan.rds'),
  logodd.tsne = readRDS('./dim reduction/manhattan/logodd.tsne.manhattan.rds'))

# do umap with manhattan distance with normalized count data
result <- umap(t(as.matrix(read.table('normalized_transofrmed_data.txt',sep = '\t'))), config = custom.config)
saveRDS(result, file='./dim reduction/manhattan/count.umap.manhattan.rds')

# graph out the result 
S.genes.count <- as.vector(unlist(read.table('sgene_included.txt', sep = '\n')))
ggplot(result$layout, aes('umap1','umap2'))  + 
  geom_point(mapping = aes(x = result$layout[,1], y = result$layout[,2], color = S.genes.count))+
  xlab('umap1') + ylab('umap2')



