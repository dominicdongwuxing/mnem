library("ks")
#setwd('/cluster/home/dongwu/rotation1')


compute.log.odds <- function(data, sgene.included,sgene.unique){
  # initialize matrix 
  R <- matrix(rep(0,dim(data)[1]*length(sgene.included)), nrow = dim(data)[1], ncol = length(sgene.included))
  # get column index of control cells: cells with "MouseNTC" 
  control.columns <- which(sgene.included == 'MouseNTC')
  # get column indices of cells perturbed with each of the S gene
  sgene.columns <- lapply(sgene.unique, function(x) which(sgene.included == x))
  for (i in 1:dim(data)[1]){
    # F for null model on E gene i
    F.ic <- kcde(data[i,control.columns],gridsize = 1000)
    # F for all S genes by the order of sgene.unique
    F.ij <- lapply(sgene.columns, function(x) kcde(data[i,x],gridsize = 1000))
    for(k in 1:length(sgene.included)){
      # only compute non control cells
      if (sgene.included[k]!='MouseNTC'){
        # find which s gene this cell k is perturbed in terms of index in sgene.unique
        sgene.k <- which(sgene.included[k] == sgene.unique)
        estimates.D <- F.ij[[sgene.k]]$estimate
        eval.points.D <- F.ij[[sgene.k]]$eval.points
        D.idx <- which.min(lapply(eval.points.D,function(x) abs(x-data[i,k])))
        
        estimates.N <- F.ic$estimate  
        eval.points.N <- F.ic$eval.points
        N.idx <- which.min(lapply(eval.points.N, function(x) abs(x-data[i,k])))
        
        D <- min(estimates.D[D.idx],1-estimates.D[D.idx])
        if (D == 0){D <- 0.000001}
        N <- min(estimates.N[N.idx],1-estimates.N[N.idx])
        if (N == 0){N <- 0.000001}
        R[i,k] <- log(D/N)
      }
    }
  }
  
  return (R[,-control.columns])
}


###################################################################################################


data <- as.matrix(read.table('normalized_transofrmed_data.txt',sep = '\t'))
sgene.included <- as.vector(unlist(read.table("sgene_included.txt", sep='\n')))
sgene.unique <- as.vector(unlist(read.table('sgene_unique.txt', sep='\n')))
logodd.result <- compute.log.odds (data,sgene.included,sgene.unique)
write.table(logodd.result,'R.csv',quote = F, sep = '\t', row.names = F, col.names = F)
