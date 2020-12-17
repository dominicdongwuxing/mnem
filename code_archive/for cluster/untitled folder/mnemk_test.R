#setwd('/Volumes/Macintosh HD - Data/study/yr2sem2/rotation1/result')
#setwd('/cluster/home/dongwu/rotation1')
library(mnem)
library(snowfall)
#####################################################################################

seed <- 2
Sgenes <- 3
Egenes <- 10
nCells <- 100
uninform <- 1
mw <- c(0.4, 0.3, 0.3)
Nems <- 3
noise <- 0.5
set.seed(seed)
simmini <- simData(Sgenes = Sgenes, Egenes = Egenes, Nems = Nems, mw = mw, nCells = nCells, 
                   uninform = uninform)
data <- simmini$data
ones <- which(data == 1)
zeros <- which(data == 0)
data[ones] <- rnorm(length(ones), 1, noise)
data[zeros] <- rnorm(length(zeros), -1, noise)

result <- mnemk(data, starts = 10,complete = T,parallel=10,ks = 1:5)
saveRDS(result, file = 'mnemk_sim.rds')

