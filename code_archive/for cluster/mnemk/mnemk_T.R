#setwd('/Volumes/Macintosh HD - Data/study/yr2sem2/rotation1/result')
#setwd('/cluster/home/dongwu/rotation1')
library(mnem)
library(snowfall)
#####################################################################################


data <- as.matrix(read.table('/cluster/home/dongwu/rotation1/R_kcde.csv'), sep = '\t')
# column names are S genes perturbed for each cell
colnames(data) <- as.vector(unlist(read.table('/cluster/home/dongwu/rotation1/sgene_R.txt', sep = '\n')))

# row names are each E gene for each observation, but cutting the long prefix 'ENSMUSG000000...'
rownames(data) <- unlist(lapply(as.vector(unlist(read.table('/cluster/home/dongwu/rotation1/egene_included.txt', sep = '\n'))),
                                function (x) strsplit(x,'_')[[1]][2]))

result <- mnemk(data, starts = 10, parallel = 10, complete = T)
saveRDS(result, file = 'mnemk_T.rds')

