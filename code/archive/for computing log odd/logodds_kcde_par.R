setwd('/cluster/home/dongwu/rotation1')
source('/Volumes/Macintosh HD - Data/study/yr2sem2/rotation1/code/functions.R')
library('snowfall')
# number of cores request on the cluster
numCores <- 32
# use parallel computing
sfInit(parallel = TRUE, cpus = numCores)
sfLibrary(ks)

# load data and other previous results
data <- as.matrix(read.table('normalized_transofrmed_data.txt',sep = '\t'))
sgene.included <- as.vector(unlist(read.table("sgene_included.txt", sep='\n')))
sgene.unique <- as.vector(unlist(read.table('sgene_unique.txt', sep='\n')))

# get column index of control cells: cells with "MouseNTC" 
control.columns <- which(sgene.included == 'MouseNTC')

# get column indices of cells perturbed with each of the S gene
sgene.columns <- lapply(sgene.unique, function(x) which(sgene.included == x))


# export all global variable and function to slaves
sfExportAll()


Logoddslist <- sfLapply(seq_len(nrow(data)), log.odds.row)
sfStop()

Logodds <- do.call("rbind",Logoddslist)
write.table(Logodds,'R_kcde.csv',quote = F, sep = '\t', row.names = F, col.names = F)


