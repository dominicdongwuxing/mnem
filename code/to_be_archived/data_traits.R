library(umap)
library(ggfortify)
library(ggplot2)

setwd('/Volumes/Macintosh HD - Data/study/yr2sem2/rotation1/result')
kde <- as.matrix(read.table('R_kde.csv'), sep = '\t')
kcde <- as.matrix(read.table('R_kcde.csv'), sep = '\t')

# check correlation between each cell over the two methods
cor.cells.pearson <- lapply(c(1:ncol(kde)), function(c) cor(kde[,c],kcde[,c]))
hist(xlim = (c(-1,1)), unlist(cor.cells.pearson), xlab = 'correlation', ylab = 'amount',
     main = 'log odds Pearson correlation between cells, kde vs kcde')

cor.cells.spearman <- lapply(c(1:ncol(kde)), function(c) cor(kde[,c],kcde[,c], method = 'spearman'))
hist(xlim = (c(-1,1)), unlist(cor.cells.spearman), xlab = 'correlation', ylab = 'amount',
     main = 'log odds Spearman correlation between cells, kde vs kcde')

cor.cells.kendall <- lapply(c(1:ncol(kde)), function(c) cor(kde[,c],kcde[,c], method = 'kendall'))
hist(xlim = (c(-1,1)), unlist(cor.cells.kendall), xlab = 'correlation', ylab = 'amount',
     main = 'log odds Kendall correlation between cells, kde vs kcde')


# check histogram of the two results
kde.val <- as.vector(kde)
kcde.val <- as.vector(kcde)

hist(xlim = c(-1,1),kcde.val, xlab = 'log odd value', ylab = 'amount', 
     main = 'log odd values of kcde vs kde result', breaks = 500, col = rgb(1,0,0,0.5))
hist(xlim = c(-1,1),kde.val, xlab = 'log odd value', ylab = 'amount', 
     main = 'log odd values of kde result', breaks = 1000, add = T, col = rgb(0,0,1,0.5))
legend('topleft', legend = c('kcde','kde'), fill = c('blue', 'red'), cex = 0.8)

# colnames are the S gene names for each cell
colnames(kcde) <- as.vector(unlist(read.table('sgene_R.txt', sep = '\n')))
S.genes <- colnames(kcde)

# do UMAP on kcde result
logodds.umap <- umap(t(kcde))

# do PCA on kcde result
logodds.pca <- prcomp(t(kcde))

# graph the dim reduction results
autoplot(logodds.pca) + geom_point(aes(color = S.genes))

ggplot(data = as.data.frame(cbind(c(1:100),cumsum(logodds.pca$sdev[1:100]))), aes(x = c(1:100), y = cumsum(logodds.pca$sdev[1:100])))  + 
        geom_point() + xlab('number of PCA component') + ylab('Cumulative variance explained in %')

ggplot(logodds.umap$layout, aes('umap1','umap2'))  + 
        geom_point(mapping = aes(x = logodds.umap$layout[,1], y = logodds.umap$layout[,2], color = S.genes))+
        xlab('umap1') + ylab('umap2') +xlim(-10,13) + ylim(-7.7,7.7)

# do dimension reduction in normalized transformed data (counts)
nor.trans.data <- as.matrix(read.table('normalized_transofrmed_data.txt',sep = '\t'))
nor.trans.pca <- prcomp(t(nor.trans.data))
nor.trans.umap <- umap(t(nor.trans.data))

S.genes.count <- as.vector(unlist(read.table('sgene_included.txt', sep = '\n')))
        
autoplot(nor.trans.pca) + geom_point(aes(color = S.genes.count))

ggplot(data = as.data.frame(cbind(c(1:100),cumsum(nor.trans.pca$sdev[1:100]))), aes(x = c(1:100), y = cumsum(nor.trans.pca$sdev[1:100])))  + 
        geom_point() + xlab('number of PCA component') + ylab('Cumulative variance explained in %')

ggplot(nor.trans.umap$layout, aes('umap1','umap2'))  + 
        geom_point(mapping = aes(x = nor.trans.umap$layout[,1], y = nor.trans.umap$layout[,2], color = S.genes.count))+
        xlab('umap1') + ylab('umap2') +xlim(-3,3) + ylim(-7.7,5)
