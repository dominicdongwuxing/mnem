setwd('/Volumes/Macintosh HD - Data/study/yr2sem2/rotation1/result')
library(mnem)
library(ggplot2)
library(ggfortify)
library(umap)
library(Rtsne)
library(gridExtra)
source('/Volumes/Macintosh HD - Data/study/yr2sem2/rotation1/code/functions.R')

# load best k result 
result <- readRDS('bestk.rds')
count.data <- as.matrix(read.table('normalized_transofrmed_data.txt',sep = '\t'))
log.data <- as.matrix(read.table('R_kcde.csv'), sep = '\t')
colnames(log.data) <- as.vector(unlist(read.table('sgene_R.txt', sep = '\n')))

# get histogram of responsibilities
graph.responsibilities (result)

# graph the raw and penalized log likelihood ratio along k values
graph.log.likelihood (result)

# inspect the responsibilities of cells and relationship to counts and log odds
check.sureness(result,count.data,log.data)

#  plots the convergence of the different EM iterations 
par(mfrow=c(2,2))
par(mar =c(5.1, 5, 4.1, 2.1))
plotConvergence(result$best, pch = 16, lwd = 3)

# plot the actual gene regulation network 
plot(result$best,oma = c(3,1,3,3), egene = F,cells = F,bestCell = F)

# do umap on logodds and normalized counts
logodd.umap <- umap(t(log.data))
count.umap <- umap(t(count.data))
saveRDS(logodd.umap,file = 'logoddumap.rds')
saveRDS(count.umap,file = 'countumap.rds')

# do pca on logodds and normalized counts
logodd.pca <- prcomp(t(log.data))
count.pca <- prcomp(t(count.data))
saveRDS(logodd.pca,file = 'logoddpca.rds')
saveRDS(count.pca,file = 'countpca.rds')

# do tSNE on logodds and normalized counts
logodd.tsne <- Rtsne(t(log.data))
count.tsne <- Rtsne(t(count.data))
saveRDS(logodd.tsne, file='logoddtsne.rds')
saveRDS(count.tsne, file='counttsne.rds')

# plot the dimension reduction results, but with only soft clustering
graph.dim.reduction.loggodd (
  result = result,
  logodd.pca = logodd.pca, 
  logodd.umap = logodd.umap,
  logodd.tsne = logodd.tsne)

# plot the umap and tSNE results, on count data
count.tsne <- readRDS('./dim reduction/counttsne.rds')
count.umap <- readRDS('./dim reduction/countumap.rds')



# discover how different parameter setting affects k 
ks <- explore.k (mnem:::modData(log.data))
saveRDS(ks, file = 'learnk.rds')


