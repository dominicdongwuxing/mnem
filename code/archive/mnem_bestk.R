setwd('/Volumes/Macintosh HD - Data/study/yr2sem2/rotation1/result')
#setwd('/cluster/home/dongwu/rotation1')
library(mnem)

#####################################################################################

data <- as.matrix(read.table('R_kcde.csv'), sep = '\t')
# column names are S genes perturbed for each cell
colnames(data) <- as.vector(unlist(read.table('sgene_R.txt', sep = '\n')))

# row names are each E gene for each observation, but cutting the long prefix 'ENSMUSG000000...'
rownames(data) <- unlist(lapply(as.vector(unlist(read.table('egene_included.txt', sep = '\n'))),
                                function (x) strsplit(x,'_')[[1]][2]))

bestk <- mnemk(data, ks = 1:5)

saveRDS(bestk, file = 'bestk.rds')

