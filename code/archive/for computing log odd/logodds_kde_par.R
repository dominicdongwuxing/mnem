setwd('/cluster/home/dongwu/rotation1')

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



# compute the log odds for each row, namely, each E gene i
log.odds.row <- function(i){
  # initialize vector
  R <- c(rep(0,length(sgene.included)))
  # F for null model on E gene i
  F.ic <- kde(data[i,control.columns],gridsize = 1000)
  # F for all S genes by the order of sgene.unique
  F.ij <- lapply(sgene.columns, function(x) kde(data[i,x],gridsize = 1000))
  for (k in 1:length(sgene.included)){
    # only compute non control cells
    if (sgene.included[k]!='MouseNTC'){
      # find which s gene this cell k is perturbed in terms of index in sgene.unique
      sgene.k <- which(sgene.included[k] == sgene.unique)
      estimates.D <- F.ij[[sgene.k]]$estimate
      eval.points.D <- F.ij[[sgene.k]]$eval.points
      # find out the estimate probability of cell k, which is the estimate of
      # the corresponding eval.point that is closest to cell k value
      D.idx <- which.min(abs(eval.points.D-data[i,k]))
      D <- estimates.D[D.idx]
        
      # do the same with N
      estimates.N <- F.ic$estimate  
      eval.points.N <- F.ic$eval.points
      N.idx <- which.min(abs(eval.points.N-data[i,k]))
      N <- estimates.N[N.idx]
        
      # deal with possible 0 
      if (D == 0 & N != 0) {
        D <- min(0.000001,N)
      } 
      else if (N == 0 & D != 0) {
        N <- min(0.000001,D)
      }
      else if (D == 0 & N == 0){
        D <- 1
        N <- 1
      }
      R[k] <- log(D/N)
    }
  }
  return (R[-control.columns])
}



# export all global variable and function to slaves
sfExportAll()


###################################################################################################


Logoddslist <- sfLapply(seq_len(nrow(data)), log.odds.row)
sfStop()

Logodds <- do.call("rbind",Logoddslist)
write.table(Logodds,'R_kde.csv',quote = F, sep = '\t', row.names = F, col.names = F)


