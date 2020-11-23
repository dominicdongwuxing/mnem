# This script contains all functions for the project
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
# this is the first section of the script, which contains functions that performs 
# a particular section of the project,
# each section is written in a function that depends on lower order functions 
# from the second section.


# This function pre-processes the data in the following way:
# it reads 1, the raw data matrix, whose rows are E gene counts and columns individual cells
# 2, tag, whose first column contains the names of sgRNAs, second column contains the 
# corresponding cells that are found with this sgRNA
# 3, a cell name file indicating the names of the each cells in raw data column
# the 3 files need to be downloaded into the data path
# it then only includes E genes (row of the raw data) with median >= 1, 
# and cells (column of the raw data) with exactly one sgRNA detected (1 S gene perturbed)
# and performs normalization and transformation using the Linnorm package and returns
# the processed data
# it depends on packages : Matrix, Linnorm
preprocessData <- function(data.path = '.'){
  # get raw expression data, row: E gene counts; column: individual cell
  raw.data <- readMM(paste0(data.path,'/GSM2396856_dc_3hr.mtx.txt'))
  
  # get perturbation information 
  # column 1: sgRNA that knocks out a particular S gene name;
  # column 2: names of the cells that were perturbed by this sgRNA
  tag <- read.csv(file = paste0(data.path, "/GSM2396856_dc_3hr_cbc_gbc_dict_strict.csv"), header = F)
  
  # names of the 32777 cells as columns of raw.data
  cell.names <- read.csv(file = paste0(data.path,"/GSM2396856_dc_3hr_cellnames.csv"))
  
  # find which sgRNA was detected in raw data by matching each cell to its detected sgRNA
  sgRNA.indices <- lapply(cell.names[,2], function(x) grep(x,tag[,2]))
  
  # fill empty entries --- cells with no detected sgRNA into "$"
  sgRNA.indices[unlist(lapply(sgRNA.indices,function(x) length(unlist(x))==0))] <- '$'
  
  # boolean vector of whether each cell is a valid cell (cells with exactly 1 sgRNA match)
  sgRNA.inclusion.bool <- unlist(lapply(sgRNA.indices,function(x) length(x) == 1 && x != '$'))
  
  # get which sgRNA is detected in the valid cells 
  sgRNA.included <- sgRNA.indices[unlist(sgRNA.inclusion.bool)]
  
  # convert sgRNA indices (numbers) from sgRNA.included into s gene (names) 
  # by matching this number, which is row number in tag, to a unique S gene
  sgene.included <- unlist(lapply(sgRNA.included, function(x) unlist(strsplit(tag[x,1],'_'))[2]))
  
  # get a boolean vector of whether each E gene/row should be included (median >= 1)
  egene.inclusion.bool <- apply(raw.data,1,function(r) median(r) >= 1)
  
  # only include E genes/rows of raw.data with E gene median <= 1 and
  # only include cells/column of raw.data with exactly 1 sgRNA detected (1 S gene perturbed)
  included.data <- raw.data[egene.inclusion.bool,sgRNA.inclusion.bool]
  
  # normalize and transform data with Linnorm
  normalized.transformed.data <- Linnorm(included.data)
  # and the names of each cell is the names of the perturbed S gene for this cell
  colnames(normalized.transformed.data) <- sgene.included
  
  return (normalized.transformed.data)
}

# This function takes normalized.transformed.data, which is the output from preprocess.data,
# computes and outputs the log odd matrix using kcde function from library ks.
# by default, it uses parallel computing with 4 cores. 
# it depends on packages: ks, snowfall
computeLogodd <- function(data = NULL, parallel = 4) {
  # initialize parallel computing
  sfInit(parallel = TRUE, cpus = parallel)
  sfLibrary(ks)
  
  # get column names of data matrix, which is the perturbed S gene for each cell
  sgene.included <- colnames(data)
  
  # get a vector of all 24 S genes, exclude MouseNTC, which is control
  sgene.unique <- unique(sgene.included)
  sgene.unique <- sgene.unique[sgene.unique != 'MouseNTC']
  
  # get column index of control cells: cells perturbed with "MouseNTC" sgRNA
  control.columns <- which(sgene.included == 'MouseNTC')
  
  # get indices of columns/cells perturbed with each of the S gene
  sgene.columns <- lapply(sgene.unique, function(x) which(sgene.included == x))
  
  # compute the log odds for each row, namely, each E gene i
  log.odds.row <- function(i){
    # initialize vector
    R <- c(rep(0,length(sgene.included)))
    # F for null model on E gene i
    F.ic <- kcde(data[i,control.columns],gridsize = 1000)
    # F for all S genes by the order of sgene.unique
    F.ij <- lapply(sgene.columns, function(x) kcde(data[i,x],gridsize = 1000))
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
        D <- min(estimates.D[D.idx],1-estimates.D[D.idx])
        
        # do the same with N
        estimates.N <- F.ic$estimate  
        eval.points.N <- F.ic$eval.points
        N.idx <- which.min(abs(eval.points.N-data[i,k]))
        N <- min(estimates.N[N.idx],1-estimates.N[N.idx])
        
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
  
  # export all global variables and functions 
  sfExportAll()
  
  # compute log odd for each row in a parallel way
  Logoddslist <- sfLapply(seq_len(nrow(data)), log.odds.row)

  sfStop()
  
  # combine rows to reconstruct the log odd matrix
  Logodds <- do.call("rbind",Logoddslist)
  
  # column (cell) names are still the pertrubed S gene for the cell
  # the control cells are excluded from log odd matrix
  colnames(Logodds) <- sgene.included[sgene.included != 'MouseNTC']
  
  return(Logodds)
}

# This function takes the log odd matrix, which is the output from compute.logodd,
# performs dimension reduction with umap (using manhattan distance), and uses mnem 
# or mnemk to cluster the cells, the arguments for mnem or mnemk can be customized
# and then it colors the umap result by using the soft clustering result from mnem
# or mnemk with gradient coloring;
# uses mnemk with default setting
# to cluster the cells;
# the mnem result (as an mnem object), and umap can be given as arguments 
# (therefore omitting those steps)
# it returns the mnem object, plots the umap (colored by clusters) result,
# plots the network, the convergence, the histogram of responsibilities, 
# the distribution of max confidence among cells, the plots are saved as svg 
# it depends on packages: mnem, snowfall, umap, ggfortify, ggplot2, gridExtra
findAndVisualizeComponents <- function(
  data = NULL, 
  mnem.result = NULL, 
  umap.result = NULL, 
  strategy = 'mnemk', 
  ksel = NULL,
  starts = NULL,
  type = NULL,
  complete = NULL,
  parallel = NULL,
  save.path = '.',
  ks = NULL,
  ...){
  # do umap on data if not given
  if (is.null(umap.result)){
    custom.config <- umap.defaults
    custom.config$random_state <- 42
    custom.config$metric <- 'manhattan'
    umap.result <- umap(t(data), config = custom.config)
  }
  
  # do mnem or mnemk if not given
  if (is.null(mnem.result)){
    # set starts is default 3 if starts not given
    if (is.null(starts)){starts = 3}
    
    # set type is default networks if not given
    if (is.null(type)){type = 'networks'}
    
    # set completeness is default false if not given
    if (is.null(complete)){complete = F}
    
    # set parallel as number of starts if not given
    if (is.null(parallel)){parallel = starts}
    
    if (strategy == 'mnemk'){
      # set the default ksel for mnemk if ksel if not given
      if (is.null(ksel)){ksel = c("kmeans", "silhouette", "cor")}
      
      # try k from 1 to 5 is ks not given
      if (is.null(ks)){ks = seq_len(5)}
      
      # do mnemk
      mnem.result <- mnemk(data, ksel=ksel,starts=starts,type=type,complete=complete,
                           parallel=parallel, ks=ks, ...)

    }
    else if (strategy == 'mnem'){
      # set a good ksel combination for mnem if ksel is not given
      if (is.null(ksel)){ksel = c("kmeans", "silhouette", "euclidean")}
      
      # do mnem
      mnem.result <- mnem(data, ksel = ksel,starts=starts,type=type,complete=complete,
                          parallel=parallel, ...)
    }
  }
  
  # if it is an mnemk result then also graph the log likelihood vs k
  # and then only extract the best result (an mnem object) to work with 
  if (!is.null(mnem.result$best)){
    graph.log.likelihood (mnem.result, save.path)
    mnem.object <- mnem.result$best
  }
  # otherwise keep working with the mnem result, which is already an mnem object
  else {mnem.object <- mnem.result}
  
  # pick the best run and plot umap with the soft clustering, if k is 2 or 3
  graph.umap (umap.result, mnem.object, 'UMAP', save.path)
  
  # pick the best run and plot convergence, network, histogram of responsibilities and
  # distribution of max confidence 
  
  # get responsibilities, k by l
  gamma <- getAffinity(x = mnem.object$probs, mw = mnem.object$mw)

  # draw and save histogram of responsibilities
  svg(filename=paste0(save.path, '/Histogram of responsibilities.svg'))
  hist(gamma, xlab = "responsibilities", ylab = 'count',
       main = "Histogram of responsibilities")
  dev.off()

  # draw and save distribution of max confidence among cells
  svg(filename=paste0(save.path, '/Distribution of max confidence among cells.svg'))
  confidence <- apply(gamma, 2, max) 
  h <- hist(confidence)
  h$density <- h$counts/sum(h$counts)
  par(mfrow=c(1,1))
  plot(h, freq = F, xlab = 'max confidence among components', ylab = 'fraction of cells',
       main = 'Distribution of max confidence among cells', col = 'grey')
  dev.off()
  
  # draw and save plot convergence
  svg(filename=paste0(save.path, '/convergence.svg'))
  par(mfrow=c(2,2))
  par(mar =c(5.1, 5.1, 4.1, 2.1))
  plotConvergence(mnem.object, pch = 20, cex.main = 0.8, cex.lab = 0.8)
  dev.off()
  
  # draw and save network if k < 10
  if (dim(gamma)[1] < 10){
    svg(filename=paste0(save.path, '/network.svg'), width = 9, height = 9)
    par(mar=c(1,1,1,1))
    plot(mnem.object,egene = F,cells = F,bestCell = F)
    dev.off()
  }

  return (mnem.result)
}

# this function combines the different stages of the project
# data from the Perturb Seq paper needs to be downloaded in the data path
# it takes an argument to stop at 3 different stages of the study, data preprocessing,
# computing log odd, or performing the mnem and then analyze the result
# default is to perform all stages 
# it accepts additional parameters for the mnem function
# it depends on packages: 
# Matrix, Linnorm if stops at preprocess;
# + ks, snowfall if stops at log odd;
# + mnem, snowfall, umap, ggfortify, ggplot2, gridExtra if stops at analyze 
createStudy <- function(stages = c('preprocess','log odd','analyze'), 
                        data.path = '.',
                        save.path = '.', ...){
  if ('preprocess' %in% stages){
    normalized.transformed.data <- preprocessData(data.path=data.path)
  }
  else {
    return ('No Result')
  }
  
  if ('log odd' %in% stages){
    R <- computeLogodd (data = normalized.transformed.data, parallel = 4)
  }
  else {
    return (normalized.transformed.data)
  }
  
  if ('analyze' %in% stages){
    mnem.result <- findAndVisualizeComponents(
      data = NULL, 
      mnem.result = NULL, 
      umap.result = NULL, 
      strategy = 'mnemk', 
      ksel = NULL,
      starts = NULL,
      type = NULL,
      completeness = NULL,
      parallel = NULL,
      save.path = save.path,
      ks = NULL,
      ...) 
    return (mnem.result)
  }
  else {
    return (R)
  }
  
}
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
# This is the second section of the script, which contains all lower order functions
# that the functions from the first section depends on, and extra functions that 
# are used in my project, but are not used in the first section functions.

# For data preprocessing:

# this function looks at the distribution of number of perturbed S-genes over all cells
show.sgRNA.distribution <- function (indices, title = "Distribution of detected sgRNAs") {
  # firstly exlcude cells with no sgRNA matched
  indices <- indices[unlist(sapply(indices,function (x) x!='$'))]
  # find number of all cells
  sum.cells <- length(indices)
  # find frequency of cells with control sgRNA (m_MouseNTC_100_A_67005 as index 50) detected
  control <- sum(unlist(sapply(indices,function (x) x == 50)))/sum.cells
  # find maximum number of S genes perturbed in a single cell, over all cells
  max.no <- max(sapply(indices,function (x) length(x)))
  # record frequencies of all different amount of non-trivial pertubation 
  freqs <- lapply(c(1:max.no), function (n) sum(sapply(indices,function (x) length(x)==n))/sum.cells)
  # for the cells with 1 sgRNA, subtract those with control sgRNA
  freqs[[1]] <- freqs[[1]] - control
  freqs <- c(control, freqs)
  print(freqs)
  barplot(unlist(freqs), names.arg = c(0:max.no),main = title, ylim = c(0,0.8), ylab = 'Frequency',xlab = 'no. of sgRNA detected per cell')
}

show.sgRNA.percentage <- function(indices, title = "Percentage of sgRNA matching"){
  # find number of all cells
  sum.cells <- length(indices)
  # find frequency of cells with no sgRNA matched
  empty <- sum(unlist(sapply(indices,function (x) x=='$')))/sum.cells
  # find frequency of cells with no and 1 sgRNA matched
  one <- sum(unlist(sapply(indices,function (x) length(x) == 1)))/sum.cells
  # subtract to get frequency of cells with no and 1 sgRNA matched
  one <- one - empty
  freqs <- c(empty, one, 1 - empty - one)
  barplot(unlist(freqs), names.arg = c('0 sgRNA', '1 sgRNA', '>1 sgRNA'), main = title, ylab = 'Frequency')
}

##########################################################################################

# For mnemk result analysis 

# extract legend from a graph
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

# graphs the histogram of responsibilities of mnemk, soft: cluster = 0 hard: cluster = 1
graph.responsibilities <- function (result, cluster = 0){
  svg(filename = 'Histogram of responsibilities.svg')
  # get k (no. of components) by l (number of cells) matrix of soft clustering pr.
  prob <- getAffinity(x = result$best$probs, mw = result$best$mw, affinity = cluster)
  k <- dim(prob)[1]
  hist(prob, xlab = "responsibilities", ylab = 'count',
       main = paste("Histogram of responsibilities under best k =",k ))
  dev.off()
}

# input: the mnemk result
# graph the raw and penalized log likelihood ratio along k values
graph.log.likelihood <- function (result, save.path = '.'){
  svg(filename = paste0(save.path, '/log likelihood ratio under various k.svg'))
  k <- seq(length(result$ics))
  raw <- result$lls
  penalized <- result$ics
  par = (mar=c(5,4,4,6) )
  par(oma=c(0,0,0,1))
  plot(k, penalized, ylim = c(min(penalized),max(penalized)), pch = 16, type = 'b', col = 'red',
       xlab = '', ylab = '', main = 'log likelihood ratio under various k', axes = F, xaxt = "n")
  axis(2, pretty(range(penalized)),col = 'red', col.axis = 'red')
  mtext("penalized",side = 2, line = 2.5, col = 'red')
  box()
  
  par(new = TRUE)
  plot(k, raw, ylim = c(min(raw),max(raw)), pch = 16, type = 'b', col = 'blue', axes = F,
       xlab = '', ylab = '', xaxt = "n")
  axis(4, pretty(range(raw)),col = 'blue', col.axis = 'blue')
  mtext('raw', col = 'blue', line = 0, outer = T, side = 4)
  
  axis(1, at = 1:k)
  mtext('k', side = 1, line = 2.5)
  legend('bottomright',legend=c('penalized','raw'), text.col = c('red','blue'), 
         pch = c(16,16), col = c('red','blue'), cex = 0.7)
  dev.off()
}

# inspect the responsibilities of cells and relationship to counts and log odds
check.sureness <- function (result, count.data, log.data, cluster = 0){
  # get responsibilities, k by l
  gamma <- getAffinity(x = result$best$probs, mw = result$best$mw, affinity = cluster)
  # get the maximum responsibility out of the ks for each cell
  # if cluster = 1 (hard clustering) it would be all 1s
  confidence <- apply(gamma, 2, max) 
  # plot histogram of Distribution of max confidence among cells
  h <- hist(confidence)
  h$density <- h$counts/sum(h$counts)
  par(mfrow=c(1,1))
  plot(h, freq = F, xlab = 'max confidence among components', ylab = 'fraction of cells',
       main = 'Distribution of max confidence among cells', col = 'grey')
  # inspect the cells that are not confident (pr. < 0.5) of being in any component
  unsure.cells.idx <- which(confidence < 0.5)
  sgene.included <- as.vector(unlist(read.table('sgene_included.txt', sep = '\n')))
  # get the mean of normalized count over E genes for each cell
  count.mean <- unlist(apply(count.data[,sgene.included != 'MouseNTC'],2,mean))
  # get the mean of log odds over E genes for each cell
  logodd.mean <- unlist(apply(log.data,2,mean))
  # get count means of confident and doubtful cells
  confident.count <- count.mean[-unsure.cells.idx]
  doubtful.count <- count.mean[unsure.cells.idx]
  # do wilcox test and extract p value
  p.count <- format(wilcox.test(confident.count,doubtful.count)$p.value,scientific = TRUE, digits = 2)
  # get log odd means of confident and doubtful cells
  confident.logodd <- logodd.mean[-unsure.cells.idx]
  doubtful.logodd <- logodd.mean[unsure.cells.idx]
  p.logodd <- format(wilcox.test(confident.logodd,doubtful.logodd)$p.value,scientific = TRUE, digits = 2)
  par(mfrow=c(1,2))
  plot(count.mean,confidence, pch = '.',
       xlab = 'mean of normalized counts over E genes', ylab = 'assignment confidence of individual cell')
  abline(v=c(mean(confident.count),mean(doubtful.count)), col = c('blue','red'), lty = 2, lwd = 3)
  legend('bottomleft', title = paste('p value:',p.count), title.col = 'black', text.font = 2,
         legend=c('confident mean','doubtful mean'), y.intersp = 0.6, x.intersp = 1,
         text.col = c('blue','red'), 
         lty = c(2,2), col = c('blue','red'),cex = 0.5,text.width = 0.1)
  plot(logodd.mean,confidence, pch = '.',
       xlab = 'mean of log odds over E genes', ylab = 'assignment confidence of individual cell')
  abline(v=c(mean(confident.logodd),mean(doubtful.logodd)), col = c('blue','red'), lty = 2, lwd = 3)
  legend('bottomright', title = paste('p value:',p.logodd), title.col = 'black', text.font = 2,
         legend=c('confident mean','doubtful mean'),y.intersp = 0.6, x.intersp = 1,
         text.col = c('blue','red'), 
         lty = c(2,2), col = c('blue','red'),cex = 0.5,text.width = 0.1)
  par(mfrow=c(1,1))
}

# graph a umap result with a mnem result clustering
graph.umap <- function (umap, mnem, title, save.path = '.'){
  soft <- getAffinity(x = mnem$probs, mw = mnem$mw, affinity = 0) 
  k <- dim(soft)[1]
  if (k == 2){
    svg(filename = paste0(save.path, '/UMAP with clustering.svg'), width = 9.4, height = 9.4)
    # make a fake graph just to get the legend
    a <- data.frame(x=seq(2),y=seq(2))
    p0 <- ggplot(a)+geom_point(mapping=aes(x=x,y=y,col=c('red','green')))+
      scale_color_manual(name='',values = c('red','green'),labels=c('100% k=1','100% k=2'))
    
    mylegend <- g_legend(p0)
    
    p <- ggplot(umap$layout, aes('umap1','umap2'))  +
      geom_point(mapping = aes(x = umap$layout[,1], y = umap$layout[,2]),color = rgb(soft[1,],soft[2,],0,0.35), size = 0.1) +
      ggtitle(title) + xlab('umap1') + ylab('umap2') + theme(plot.title = element_text(hjust = 0.5,size = 13))
    
    grid.arrange(p, mylegend, nrow = 1,widths = c(10,1))
    dev.off()
  }
  
  else if (k == 3){
    svg(filename = paste0(save.path,'/UMAP with clustering.svg'), width = 9.4, height = 9.4)
    # make a fake graph just to get the legend
    a <- data.frame(x=seq(3),y=seq(3))
    p0 <- ggplot(a)+geom_point(mapping=aes(x=x,y=y,col=c('red','green','blue')))+
      scale_color_manual(name='',values = c('red','green','blue'),labels=c('100% k=1','100% k=2','100% k=3'))
    
    mylegend <- g_legend(p0)
    
    p <- ggplot(umap$layout, aes('umap1','umap2'))  +
      geom_point(mapping = aes(x = umap$layout[,1], y = umap$layout[,2]),color = rgb(soft[1,],soft[2,],soft[3,],0.35), size = 0.1) +
      ggtitle(title) + xlab('umap1') + ylab('umap2') + theme(plot.title = element_text(hjust = 0.5,size = 13))
    
    grid.arrange(p, mylegend, nrow = 1,widths = c(10,1))
    dev.off()
  }
}

graph.tSNE <- function (tsne, mnem, title){
  soft <- getAffinity(x = mnem$probs, mw = mnem$mw, affinity = 0) 
  k <- dim(soft)[1]
  # make a fake graph just to get the legend
  a <- data.frame(x=seq(3),y=seq(3))
  p0 <- ggplot(a)+geom_point(mapping=aes(x=x,y=y,col=c('red','green','blue')))+
    scale_color_manual(name='',values = c('red','green','blue'),labels=c('100% k=1','100% k=2','100% k=3'))
  
  mylegend <- g_legend(p0)
  
  p <- ggplot(tsne$Y, aes('tsne1','tsne2'))  +
    geom_point(mapping = aes(x = tsne$Y[,1], y = tsne$Y[,2]),color = rgb(soft[1,],soft[2,],soft[3,],0.35), size = 0.1) +
    ggtitle(title) + xlab('tsne1') + ylab('tsne2') 
  
  grid.arrange(p, mylegend, nrow = 1,widths = c(10,1))
}

# discover how different parameter setting affects k 
explore.k <- function(data){
  ksel1 <- c('kmeans','hc')
  ksel2 <- c('silhouette','BIC','AIC')
  ksel3 <- c('cor','euclidean','manhattan')
  ks <- list()
  
  for (i in ksel1){
    for (j in ksel2){
      for (k in ksel3){
        if (!(i == 'hc' && (j == 'BIC' || j == 'AIC'))){
          ksel <- paste(i,j,k)
          result <- mnem:::learnk(data, ksel = c(i,j,k))
          ks[[ksel]] <- result$k
          saveRDS(result, file = paste0(i,j,k,'.rds'))
        }
      }
    }
  }
  
  return (ks)
}


# The next 3 functions are not used

# graph the pca, umap and tSNE on log odd data, using the mnemk result
graph.dim.reduction.loggodd <- function(
  result = result,
  logodd.pca = logodd.pca, 
  logodd.umap = logodd.umap,
  logodd.tsne = logodd.tsne){
  
  # get soft assignments of the cells
  hard <- getAffinity(x = result$best$probs, mw = result$best$mw, affinity = 1)
  soft <- getAffinity(x = result$best$probs, mw = result$best$mw, affinity = 0)
  
  # make a fake graph just to get the legend
  a <- data.frame(x=seq(3),y=seq(3))
  p0 <- ggplot(a)+geom_point(mapping=aes(x=x,y=y,col=c('red','green','blue')))+
    scale_color_manual(name='',values = c('red','green','blue'),labels=c('100% k=1','100% k=2','100% k=3'))
  
  mylegend <- g_legend(p0)
  
  # pca plot
  embedding <- logodd.pca$x[,1:2]
  # get variance explained by the first 2 PCs
  var <- logodd.pca$sdev^2
  pc1.xlab <- paste('PC1 (',format(var[1]/sum(var)*100, digits = 3),'%)')
  pc2.xlab <- paste('PC2 (',format(var[2]/sum(var)*100, digits = 3),'%)')
  
  p1 <- ggplot(embedding, aes(pc1.xlab,pc1.xlab))  +
    geom_point(mapping = aes(x = embedding[,1], y = embedding[,2]),color = rgb(soft[1,],soft[2,],soft[3,],0.35), size = 0.1)+
    xlab(pc1.xlab) + ylab(pc2.xlab) + coord_fixed(ratio=0.8) +  
    ggtitle("PCA") + theme(plot.title = element_text(hjust = 0.5,size = 13))
  
  # umap plot
  p2 <- ggplot(logodd.umap$layout, aes('umap1','umap2'))  +
    geom_point(mapping = aes(x = logodd.umap$layout[,1], y = logodd.umap$layout[,2]),color = rgb(soft[1,],soft[2,],soft[3,],0.35), size = 0.1)+
    xlab('umap1') + ylab('umap2') + coord_fixed(ratio=1) + 
    ggtitle("UMAP") + theme(plot.title = element_text(hjust = 0.5,size = 13))
  # if doing euclidean distance result, then add  + xlim(c(-5,5)) + ylim(c(-7,3.5)) +
  
  # visualize tSNE
  p3 <- ggplot(logodd.tsne$Y, aes('tSNE1','tSNE2'))  +
    geom_point(mapping = aes(x = logodd.tsne$Y[,1], y = logodd.tsne$Y[,2]),color = rgb(soft[1,],soft[2,],soft[3,],0.35), size = 0.1)+
    xlab('tSNE1') + ylab('tSNE2') + coord_fixed(ratio=0.9) + 
    ggtitle("tSNE") + theme(plot.title = element_text(hjust = 0.5,size = 13))
  
  grid.arrange(arrangeGrob(p1,p2,p3,nrow = 1), mylegend, nrow = 1,widths = c(10,1))
  
}

# takes an mnem result (not mnemk result) and graphs the umap and tSNE clustering
# of log odd and count data
# only graph when k=2 or 3
graph.dim.reduction <- function (
  result = result,
  logodd.umap = logodd.umap,
  logodd.tsne = logodd.tsne,
  count.umap = count.umap,
  count.tsne = count.tsne,
  path = '.'
){
  # this is the s gene names of all columns (cells) of the count data
  sgene.included <- as.vector(unlist(read.table('../sgene_included.txt', sep = '\n')))
  # boolean of all the control and non control cells
  control <- sgene.included == 'MouseNTC'
  soft <- getAffinity(x = result$probs, mw = result$mw, affinity = 0) 
  
  k <- dim(soft)[1]
  if (k == 2){
    # make fake graphs just to get the legends
    a.count <- data.frame(x=seq(3),y=seq(3))
    p.count <- ggplot(a.count) + geom_point(mapping=aes(x=x,y=y,color = c('red','blue','black'))) +
      scale_color_manual(values = c('red','blue','black'),
                         name = '',
                         labels = c('100% k=1','100% k=2','control'))  + 
      theme(legend.position = "bottom", legend.box = "horizontal")
    
    mylegend.count <- g_legend(p.count)
    
    a.logodd <- data.frame(x=seq(2),y=seq(2))
    p.logodd <- ggplot(a.logodd)+geom_point(mapping=aes(x=x,y=y,col=c('red','blue')))+
      scale_color_manual(name='',values = c('red','blue'),labels=c('100% k=1','100% k=2')) +
      theme(legend.position = "bottom", legend.box = "horizontal")
    
    mylegend.logodd <- g_legend(p.logodd)
    
    p1<- ggplot(count.umap$layout,aes('umap1','umap2')) +
      geom_point(mapping=aes(x=count.umap$layout[,1],y=count.umap$layout[,2]),
                 color=ifelse(control,'black',rgb(soft[1,],0,soft[2,],0.35)), size = 0.1)+
      
      xlab('umap1') + ylab('umap2') +  xlim(c(-3.5,3)) + coord_fixed(ratio=0.45) +
      ggtitle("UMAP") + theme(plot.title = element_text(hjust = 0.5,size = 10),plot.margin=unit(c(1,1,-0.5,1), "cm"))
    
    
    p2 <- ggplot(count.tsne$Y,aes('tSNE1','tSNE2')) +
      geom_point(mapping=aes(x=count.tsne$Y[,1],y=count.tsne$Y[,2]),
                 color = ifelse(control,'black',rgb(soft[1,],0,soft[2,],0.35)), size = 0.1)+
      xlab('tSNE1') + ylab('tSNE2') + coord_fixed() +
      ggtitle("tSNE") + theme(plot.title = element_text(hjust = 0.5, size = 10),plot.margin=unit(c(1,1,-0.5,1), "cm"))
    
    p3<- ggplot(logodd.umap$layout,aes('umap1','umap2')) +
      geom_point(mapping=aes(x=logodd.umap$layout[,1],y=logodd.umap$layout[,2]),
                 color=rgb(soft[1,],0,soft[2,],0.35), size = 0.1)+
      
      xlab('umap1') + ylab('umap2') +  xlim(c(-5,5)) + ylim(c(-7,3.5)) + coord_fixed(ratio=1) +
      ggtitle("UMAP") + theme(plot.title = element_text(hjust = 0.5,size = 10),plot.margin=unit(c(1,1,-0.5,1), "cm"))
    
    
    p4 <- ggplot(logodd.tsne$Y,aes('tSNE1','tSNE2')) +
      geom_point(mapping=aes(x=logodd.tsne$Y[,1],y=logodd.tsne$Y[,2]),
                 color = rgb(soft[1,],0,soft[2,],0.35), size = 0.1)+
      xlab('tSNE1') + ylab('tSNE2') + coord_fixed() +
      ggtitle("tSNE") + theme(plot.title = element_text(hjust = 0.5, size = 10),plot.margin=unit(c(1,1,-0.5,1), "cm"))
    
    svg(filename = paste0(path,'/clustering on count.svg'))
    grid.arrange(arrangeGrob(p1,p2,nrow = 1), mylegend.count,nrow = 2, heights = c(15,1))
    dev.off()
    
    svg(filename = paste0(path,'/clustering on logodd.svg'))
    grid.arrange(arrangeGrob(p3,p4,nrow = 1), mylegend.logodd, nrow = 2,heights=c(15,1))
    dev.off()
    
  }
  else{
    # make fake graphs just to get the legends
    a.count <- data.frame(x=seq(4),y=seq(4))
    p.count <- ggplot(a.count) + geom_point(mapping=aes(x=x,y=y,color = c('red','green','blue','black'))) +
      scale_color_manual(values = c('red','green','blue','black'),
                         name = '',
                         labels = c('100% k=1','100% k=2','100% k=3','control')) +
      theme(legend.position = "bottom", legend.box = "horizontal")
    mylegend.count <- g_legend(p.count)
    
    a.logodd <- data.frame(x=seq(3),y=seq(3))
    p.logodd <- ggplot(a.logodd)+geom_point(mapping=aes(x=x,y=y,col=c('red','green','blue')))+
      scale_color_manual(name='',values = c('red','green','blue'),labels=c('100% k=1','100% k=2','100% k=3')) +
      theme(legend.position = "bottom", legend.box = "horizontal")
    
    mylegend.logodd <- g_legend(p.logodd)
    
    p1<- ggplot(count.umap$layout,aes('umap1','umap2')) +
      geom_point(mapping=aes(x=count.umap$layout[,1],y=count.umap$layout[,2]),
                 color=ifelse(control,'black',rgb(soft[1,],soft[2,],soft[3,],0.35)), size = 0.1)+
      
      xlab('umap1') + ylab('umap2') +  xlim(c(-3.5,3)) + coord_fixed(ratio=0.42) +
      ggtitle("UMAP") + theme(plot.title = element_text(hjust = 0.5,size = 10),plot.margin=unit(c(1,1,-0.5,1), "cm"))
    
    
    p2 <- ggplot(count.tsne$Y,aes('tSNE1','tSNE2')) +
      geom_point(mapping=aes(x=count.tsne$Y[,1],y=count.tsne$Y[,2]),
                 color = ifelse(control,'black',rgb(soft[1,],soft[2,],soft[3,],0.35)), size = 0.1)+
      xlab('tSNE1') + ylab('tSNE2') + coord_fixed() +
      ggtitle("tSNE") + theme(plot.title = element_text(hjust = 0.5, size = 10),plot.margin=unit(c(1,1,-0.5,1), "cm"))
    
    svg(filename = paste0(path,'/clustering on count.svg'))
    grid.arrange(arrangeGrob(p1,p2,nrow = 1), mylegend.count, nrow = 2,heights = c(15,1))
    dev.off()
    
    p3<- ggplot(logodd.umap$layout,aes('umap1','umap2')) +
      geom_point(mapping=aes(x=logodd.umap$layout[,1],y=logodd.umap$layout[,2]),
                 color=rgb(soft[1,],soft[2,],soft[3,],0.35), size = 0.1)+
      
      xlab('umap1') + ylab('umap2') +  xlim(c(-5,5)) + ylim(c(-7,3.5)) + coord_fixed(ratio=1) +
      ggtitle("UMAP") + theme(plot.title = element_text(hjust = 0.5,size = 10),plot.margin=unit(c(1,1,-0.5,1), "cm"))
    
    
    p4 <- ggplot(logodd.tsne$Y,aes('tSNE1','tSNE2')) +
      geom_point(mapping=aes(x=logodd.tsne$Y[,1],y=logodd.tsne$Y[,2]),
                 color = rgb(soft[1,],soft[2,],soft[3,],0.35), size = 0.1)+
      xlab('tSNE1') + ylab('tSNE2') + coord_fixed() +
      ggtitle("tSNE") + theme(plot.title = element_text(hjust = 0.5, size = 10),plot.margin=unit(c(1,1,-0.5,1), "cm"))
    
    svg(filename = paste0(path,'/clustering on logodd.svg'))
    grid.arrange(arrangeGrob(p3,p4,nrow = 1), mylegend.logodd, nrow = 2,heights = c(15,1))
    dev.off()
  }
}


#graph the umap and tSNE on count data, using the mnemk result
graph.dim.reduction.count <- function(
  result = result,
  count.umap = count.umap,
  count.tsne = count.tsne
){
  # this is the s gene names of all columns (cells) of the count data
  sgene.included <- as.vector(unlist(read.table('sgene_included.txt', sep = '\n')))
  # boolean of all the control and non control cells
  control <- sgene.included == 'MouseNTC'
  non.control <- sgene.included != 'MouseNTC'
  soft <- getAffinity(x = result$best$probs, mw = result$best$mw, affinity = 0) 
  
  # make a fake graph just to get the legend
  a <- data.frame(x=seq(4),y=seq(4))
  p0 <- ggplot(a) + geom_point(mapping=aes(x=x,y=y,color = c('red','green','blue','black'))) +
    scale_color_manual(values = c('red','green','blue','black'),
                       name = '',
                       labels = c('100% k=1','100% k=2','100% k=3','control'))  
  mylegend <- g_legend(p0)

  p1<- ggplot(count.umap$layout,aes('umap1','umap2')) +
    geom_point(mapping=aes(x=count.umap$layout[,1],y=count.umap$layout[,2]),
               color=ifelse(control,'black',rgb(soft[1,],soft[2,],soft[3,],0.35)), size = 0.1)+

    xlab('umap1') + ylab('umap2') +  xlim(c(-3.5,3)) + coord_fixed(ratio=0.42) +
    ggtitle("UMAP") + theme(plot.title = element_text(hjust = 0.5,size = 10))


  p2 <- ggplot(count.tsne$Y,aes('tSNE1','tSNE2')) +
    geom_point(mapping=aes(x=count.tsne$Y[,1],y=count.tsne$Y[,2]),
               color = ifelse(control,'black',rgb(soft[1,],soft[2,],soft[3,],0.35)), size = 0.1)+
    xlab('tSNE1') + ylab('tSNE2') + coord_fixed() +
    ggtitle("tSNE") + theme(plot.title = element_text(hjust = 0.5, size = 10))
    
  grid.arrange(arrangeGrob(p1,p2,nrow = 1), mylegend, nrow = 1,widths = c(10,1),
               top = 'Dimension reduction on count data with soft clustering')
  
}



##########################################################################################

# For mnem result analysis

# for each set of parameters 
# extract mnem files and combine it, save the combined version
combine.mnem <- function(parameters, files){
  mnems <- readRDS(files[1])
  # if there is only 1 file then it doesn't need to be combined
  if (length(files) > 1){
    for (i in 2:length(files)){
      mnems$limits[[i]] <- readRDS(files[i])$limits[[1]]
    }
    # get the latest log likelihood for each run and find out the index of the best run
    idx <- which.max(sapply(seq(length(mnems$limits)), 
                            function (x) mnems$limits[[x]]$ll[length(mnems$limits[[x]]$ll)]))
    # convert the objects (other than limits) of mnem into that of the best run
    # don't need to change the "complete" since it is the same 
    mnems$mw <- mnems$limits[[idx]]$mw
    mnems$probs <- mnems$limits[[idx]]$probs
    mnems$lls <- mnems$limits[[idx]]$ll
    mnems$phievo <- mnems$limits[[idx]]$phievo
    mnems$thetaevo <- mnems$limits[[idx]]$thetaevo
    mnems$mwevo <- mnems$limits[[idx]]$mwevo
    mnems$ll <- mnems$lls[length(mnems$lls)] 
    mnems$comp <- readRDS(files[idx])$comp
  }

  # save the combined file
  saveRDS(mnems, file = paste0('./combined mnem','/',parameters,'.rds'))
  return(mnems)
}

# this function reads the mnem results for different parameter settings
# and summarises results into a matrix containing 
# the parameter setting itself, the log likelihood and its sd (if there are replicates),
# number of k and the molecular weights
make.summary <- function(path){
  result.summary<- NULL
  # get all the rds files 
  files <- list.files(path = path, pattern = '.rds')
  for (file in files){
    mnem.result <- readRDS(paste0(path,'/',file))
    parameters <- strsplit(strsplit(file,'.rds')[[1]],'_')[[1]]
    type <- parameters[1]
    complete <- parameters[2]
    ksel1 <- parameters[3]
    ksel2 <- parameters[4]
    ksel3 <- parameters[5]
    ll <- mnem.result$ll
    k <- length(mnem.result$mw)
    mw <- paste(unlist(lapply(mnem.result$mw, function(x) round(x,2))),collapse = ',')

    # calculate the standard deviation of the log likelihood at the each of each convergence 
    #sd <- sd(unlist(lapply(seq_len(length(mnem.result$limits)), function (x) mnem.result$limits[[x]]$ll[length(mnem.result$limits[[x]]$ll)])))

    res.vec <- c(type,complete,ksel1,ksel2,ksel3,ll,k,mw)
    result.summary <- rbind(result.summary,res.vec)
  }
  colnames(result.summary) <- c('type','complete','clustering method',
                                'model selection criteria','distance metric',
                                'best loglikelihood','number of k','molecular weights')
  rownames(result.summary) <- NULL
  return (result.summary)
}

# plot the summary of mnem result and save to path 
# the graph compares various type, k and complete
plot.summary <- function(result.summary, path = '.'){
  svg(filename = paste0(path,'/parameter comparison.svg'))
  result.summary <- as.data.frame(result.summary)
  types <- c('random','cluster','cluster2','cluster3','networks')
  # for k = 2 or 3, and complete = true or false, we build four vectors of length 
  # 5 = number of types
  # each vector at each position is the best log likelihood under the particular type,
  # complete, and number of k
  k2T <- unlist(lapply(types, 
                        function (x) 
                          as.numeric(max(result.summary[which(result.summary['type'] == x 
                                         & result.summary['number of k'] == 2 
                                         & result.summary['complete'] == 'TRUE'),
                                   'best loglikelihood']))))
  k2F <- unlist(lapply(types, 
                        function (x) 
                          as.numeric(max(result.summary[which(result.summary['type'] == x 
                                         & result.summary['number of k'] == 2 
                                         & result.summary['complete'] == 'FALSE'),
                                   'best loglikelihood']))))
  k3T <- unlist(lapply(types, 
                        function (x) 
                          as.numeric(max(result.summary[which(result.summary['type'] == x 
                                         & result.summary['number of k'] == 3 
                                         & result.summary['complete'] == 'TRUE'),
                                   'best loglikelihood']))))
  k3F <- unlist(lapply(types, 
                       function (x) 
                         as.numeric(max(result.summary[which(result.summary['type'] == x 
                                        & result.summary['number of k'] == 3 
                                        & result.summary['complete'] == 'FALSE'),
                                  'best loglikelihood']))))
  
  p <- ggplot() + 
    geom_line(aes(x=types,y=k2T,col='k = 2', linetype='observed'), group='k=2, observed',  size = 0.5) +
    geom_line(aes(x=types,y=k2F,col='k = 2', linetype='expected'),group='k=2, expected',  size = 0.5) +
    geom_line(aes(x=types,y=k3T,col='k = 3', linetype='observed'), group='k=3, observed', size = 0.5) +
    geom_line(aes(x=types,y=k3F,col='k = 3', linetype='expected'), group='k=3, expected',  size = 0.5) + xlab('') + ylab('log likelihood') +
    theme_bw() + scale_color_discrete('number of k') + scale_linetype_discrete('log likelihood')
    
  dev.off()
}
