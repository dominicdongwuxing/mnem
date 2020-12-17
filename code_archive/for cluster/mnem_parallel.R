#setwd('/Volumes/Macintosh HD - Data/study/yr2sem2/rotation1/result')
#setwd('/cluster/home/dongwu/rotation1')
library(mnem)
library(snowfall)
#####################################################################################

type <- commandArgs(TRUE)[1]
paste('type:',type,typeof(type))
starts <- as.numeric(commandArgs(TRUE)[2])
paste('starts:',starts,typeof(starts))
completeness <- as.logical(commandArgs(TRUE)[3])
paste('completeness:',completeness,typeof(completeness))
run <- commandArgs(TRUE)[4]
paste('run:',run,typeof(run))
ksel1 <-commandArgs(TRUE)[5]
paste('ksel1:',ksel1,typeof(ksel1))
ksel2 <-commandArgs(TRUE)[6]
paste('ksel2:',ksel2,typeof(ksel2))
ksel3 <- commandArgs(TRUE)[7]
paste('ksel3:',ksel3,typeof(ksel3))
cores <- as.numeric(commandArgs(TRUE)[8])

data <- as.matrix(read.table('/cluster/home/dongwu/rotation1/R_kcde.csv'), sep = '\t')
# column names are S genes perturbed for each cell
colnames(data) <- as.vector(unlist(read.table('/cluster/home/dongwu/rotation1/sgene_R.txt', sep = '\n')))

# row names are each E gene for each observation, but cutting the long prefix 'ENSMUSG000000...'
rownames(data) <- unlist(lapply(as.vector(unlist(read.table('/cluster/home/dongwu/rotation1/egene_included.txt', sep = '\n'))),
                                function (x) strsplit(x,'_')[[1]][2]))

result <- mnem(data, starts = starts, type = type, complete = completeness, ksel = c(ksel1,ksel2,ksel3), parallel = cores)
saveRDS(result, file = paste0('/cluster/scratch/dongwu/mnem/',type,'_',completeness,'_',ksel1,'_',ksel2,'_',ksel3,'_',run,'.rds'))

