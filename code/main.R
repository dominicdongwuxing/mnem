library(Matrix,Linnorm,ggplot2,ggfortify,gridExtra,ks,snowfall,mnem,umap,ggpubr,Rtsne)

# set working directory to where the following raw data files are stored:
# GSM2396856_dc_3hr.mtx.txt, GSM2396856_dc_3hr_cellnames.csv
# path <- ...

# obtain normalized transformed data while plotting distribution of number of sgRNA in cells and S-genes 
normalized.transformed.data <- preprocessData (data.path = path, plot.distribution = T)

# do PCA on processed count data and visualize
nor.trans.pca <- prcomp(t(normalized.transformed.data))
autoplot(nor.trans.pca) + geom_point(aes(color=colnames(normalized.transformed.data))) + 
        theme_bw() + scale_color_discrete(name='S-gene')

# do UMAP on processed count data and visualize
custom.config <- umap.defaults
# first customize the UMAP with random seed and Manhattan distance metrics
custom.config$random_state <- 42
custom.config$metric <- 'manhattan'
nor.trans.umap <- umap(t(normalized.transformed.data), config = custom.config)
ggplot(nor.trans.umap$layout, aes('umap1','umap2')) +
        geom_point(mapping = aes(x = nor.trans.umap$layout[,1], y = nor.trans.umap$layout[,2], color = colnames(normalized.transformed.data)))+
        xlab('umap1') + ylab('umap2') + theme_bw() + scale_color_discrete(name='S-gene')

# compute log odd data from normalized transformed data
# use kernel cumulative density estimation (kcde) and kernel density estimation (kde) 
# to compute log odd matrix from normalized.transformed.data
kde <- computeLogodd (data = normalized.transformed.data, parallel = 4, kcde = F)
kcde <- computeLogodd (data = normalized.transformed.data, parallel = 4, kcde = T)

# check correlation between each cell over the two methods using 
# Kendall rank correlation coefficient
cor.cells.kendall <- unlist(lapply(c(1:ncol(kde)), function(c) cor(kde[,c],kcde[,c], method = 'kendall')))
p1 <- ggplot(as.data.frame(cor.cells.kendall),aes(x=cor.cells.kendall)) +
        geom_histogram(binwidth = 0.01,color = 'grey') +
        xlab('Kendall rank correlation coefficient') +
        ylab('amount of cells') + theme_bw() + theme(
                axis.text.x = element_text(size = 12),
                axis.text.x = element_text(size = 12))

# check histogram of the two results
kde.val <- as.vector(kde)
kcde.val <- as.vector(kcde)
p2 <- ggplot() + geom_histogram(aes(x=kcde.val,fill=rgb(0,0,1,0.5)),binwidth = 0.1,color=rgb(0,0,1,0.5)) + 
        geom_histogram(aes(x=kde.val,fill=rgb(1,0,0,0.5)),binwidth = 0.1,color=rgb(1,0,0,0.5)) +
        xlim(-1,1) + xlab('log odd value') + 
        scale_fill_manual(name='',values = c(rgb(0,0,1,0.5),rgb(1,0,0,0.5)), labels=c('kde','kcde'))+
        theme_bw() + theme(legend.position = c(0.85,0.85))+ theme(
                axis.text.x = element_text(size = 12),
                axis.text.x = element_text(size = 12))

grid.arrange(p1,p2,nrow=1)

# use kcde result as log odd matrix R
R <- kcde

# do UMAP on kcde result
logodds.umap <- umap(t(R), config = custom.config)

# do tSNE on kcde result
logodds.tsne <- Rtsne(dist(t(R), method = 'manhattan'), is_distance = T)

# do PCA on kcde result
logodds.pca <- prcomp(t(R))

# to do M&NEM on log odd data, first explore the best k 
# the first method is to set k = {1,2,3,4,5} and use the mnemk function
# which returns the best model among the 5 k values
mnemk.result <- mnemk(R,ks=1:5,starts=10)

# graph the raw and penalized log likelihood ratio along k values
graph.log.likelihood (mnemk.result)



# using mnemk and try all 5 ks takes a long time,
# alternatively, try to infer k by clustering the data first,
# with the mnem:::learnk function, which takes the set of arguments:
# ksel = ('kmeans'/'hc','silhouette'/'AIC'/'BIC','euclidean'/'cor'/'manhattan')
ks <- explore.k(R)

# visualize the result of k under different ksels 
ggballoonplot(ks, x = "ksel2", y = "k", size = 2, fill='black',
              facet.by = c("ksel1",'ksel3'),
              ggtheme = theme_bw()) + theme(axis.text=element_text(size=12)) 

# we know that ksel = c('kmeans','AIC'/'BIC'/'silhouette','euclidean') will give
# us K=2, and we want to compare the mnem result with such ksel setting with mnemk
# result when K=2 (the best K)
# run mnem by setting ksel = c('kmeans','AIC'/'BIC'/'silhouette','euclidean')
# with 10 replications each and save the files in the current directory
# compare mnemk (type = 'netowrks', complete = F) with the same setting using
# mnem 
plot.mnemk.mnem.comparison (mnemk.F,path = '.')

# now, use mnem function with R, and explore some other hyper-parameters:
# in addition to all ksel combinations, also try different initialization methods: 
# set type = 'networks'/'random'/'cluster'/'cluster2'/'cluster3'
# set complete = T/F, and do 100 replications for each combination by using HPC
# (ETH Euler). Note that for "cluster2" and "cluster3" only one replication is needed
# because they will produce the same result
# for type=networks/random, I set starts=1 in mnem and run for 100 times separately
# (each run/start is stochastic) and combine them back into one file so that it is faster 
# for type=cluster I set starts=100 and run for 1 time, this is because it is stochastic
# but the first start is always the same, so I would get same results if I set starts=1
# and run 100 times; for type=cluster2/3, it uses the same start (deterministic), 

# now download all the results into the working directory
# and create a new directory called "data" to store the combined files
# if it needs to be combined
dir.create('data')

# combine the mnem results for different parameters
for (type in c('random','networks','cluster','cluster2','cluster3')){
        for (ksel1 in c('kmeans','hc')){
                for(ksel2 in c('silhouette','AIC','BIC')){
                        # skip combinations when hc is used in ksel1, but silhouette is not used in ksel2
                        if (!(ksel1 == 'hc' && ksel2 != 'silhouette')){
                                for (ksel3 in c('cor','euclidean')){
                                        for (completeness in c('TRUE','FALSE')){
                                                # analyze individual run results with T or F
                                                parameters <- paste(type,completeness,ksel1,ksel2,ksel3,sep='_')
                                                files <- list.files(pattern = parameters)
                                                if (length(files) != 0){
                                                        mnem <- combine.mnem (parameters,files)
                                                        saveRDS(mnem,file = paste0('./data/',parameter,'.rds'))
                                                }
                                        }
                                }
                        }
                }
        }
}

# record the best log likelihood, average log likelihood with sd, 
# and number of k among all runs for each parameter setting
result.summary <- make.summary('./data')

# draw graph to compare various type, k and complete
plot.summary(result.summary)

# do wilcox test to compare different methods 
test.methods('./data')



# finally we conclude that the best set of hyper-parameters for mnem are:
# type=networks, k=2, complete=F, and we analyze the result
mnem.final.result <- readRDS('./data/networks_FALSE_kmeans_silhouette_euclidean.rds')

# we visualize the result of the best MNEM model
findAndVisualizeComponents(
        data = R, 
        mnem.result = mnem.final.result, 
        umap.result = logodds.umap, 
        save.path = '.')

# if there are cells with maximal confidence < 50%, inspect the responsibilities 
# of cells and relationship to counts and log odds
# and plot the relationship between max confidence of a cell to its mean (over E-genes)
# log odd value in R, and mean count value in normalized transformed data
# check.sureness(mnem.final.result,normalized.transformed.data,R)

# we can also visualize the result in both UMAP and t-SNE
graph.dim.reduction(mnem.final.result,logodds.umap,logodds.tsne)
