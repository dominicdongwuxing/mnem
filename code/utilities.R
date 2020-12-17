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
# it depends on packages : Matrix, Linnorm and ggplot2 if plot.distribution is set TRUE
preprocessData <- function(data.path = '.', plot.distribution = F){
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
  
  # plot detected sgRNA distribution and bar plot of S-gene distribution 
  if (plot.distribution){
    dist <- as.data.frame(table(unlist(lapply(sgRNA.indices, function(x) length(unlist(x))))))
    dist[,2] <- dist[,2]/sum(dist[,2])
    dist[6,2]<-sum(dist[6:9,2])
    dist <- dist[1:6,]
    # p1 is the distribution of number of sgRNA detected in cells
    ggplot(dist, aes(x=Var1,y=Freq)) + geom_bar(stat='identity') + 
      scale_x_discrete(name = 'number of sgRNAs detected in cells', labels = c('0','1','2','3','4','>=5')) + 
      scale_y_continuous(name = 'Fraction of cells', breaks = seq(0,0.6,by=0.1)) +theme_bw()
  }
  
  # fill empty entries --- cells with no detected sgRNA into "$"
  sgRNA.indices[unlist(lapply(sgRNA.indices,function(x) length(unlist(x))==0))] <- '$'
  
  # boolean vector of whether each cell is a valid cell (cells with exactly 1 sgRNA match)
  sgRNA.inclusion.bool <- unlist(lapply(sgRNA.indices,function(x) length(x) == 1 && x != '$'))
  
  # get which sgRNA is detected in the valid cells 
  sgRNA.included <- sgRNA.indices[unlist(sgRNA.inclusion.bool)]
  
  # convert sgRNA indices (numbers) from sgRNA.included into s gene (names) 
  # by matching this number, which is row number in tag, to a unique S gene
  sgene.included <- unlist(lapply(sgRNA.included, function(x) unlist(strsplit(tag[x,1],'_'))[2]))
  
  # plot distribution of S-genes over included cells
  if (plot.distribution){
    dist.s <- as.data.frame(table(sgene.included))
    ggplot(dist.s,aes(x=sgene.included,y=Freq)) + geom_bar(stat='identity') + 
      xlab('S-gene') + ylab('amount of cells') + theme_bw() + theme(axis.text.x = element_text(angle = 45,size=11))
    
  }
  
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
# can choose between using kcde (cumulative density estimation) or kde (density estimation)
# it depends on packages: ks, snowfall
computeLogodd <- function(data = NULL, parallel = 4, kcde = T) {
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
  
  if (kcde){
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
  }
  else {
    log.odds.row <- function(i){
      # initialize vector
      R <- c(rep(0,length(sgene.included)))
      # F for null model on E gene i
      F.ic <- kde(data[i,control.columns],gridsize = 1000)
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
# it depends on packages: mnem, snowfall, umap, ggfortify, ggplot2, gridExtra,Cairo
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
      if (is.null(ksel)){ksel = c("hc", "silhouette", "euclidean")}
      
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
       main = 'Distribution of max confidence among cells', col = 'grey',xlim = c(0,1))
  dev.off()
  
  # draw and save plot convergence
  svg(filename=paste0(save.path, '/convergence.svg'))
  par(mfrow=c(2,2))
  par(mar =c(5.1, 5.1, 4.1, 2.1))
  plotConvergence(mnem.object, pch = 20, cex.main = 0.8, cex.lab = 0.8, cex = 0.3)
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
# + mnem, snowfall, umap, ggfortify, ggplot2, gridExtra, Cairo if stops at analyze 
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



##########################################################################################

# For exploring the effect of ksel parameters on number of cluster k

# discover how different parameter setting affects k 
explore.k <- function(data){
  ksel1 <- c('kmeans','hc')
  ksel2 <- c('silhouette','BIC','AIC')
  ksel3 <- c('cor','euclidean','manhattan')
  ks <- data.frame()
  
  for (i in ksel1){
    for (j in ksel2){
      for (k in ksel3){
        if (!(i == 'hc' && (j == 'BIC' || j == 'AIC'))){
          ksel <- paste(i,j,k)
          result <- mnem:::learnk(data, ksel = c(i,j,k))
          if (i=='hc'){i<-'hierarchical clustering'}
          if (k=='cor'){k<-'correlation'}
          ks <- rbind(ks,c(i,j,k,result$k))
          saveRDS(result, file = paste0(i,j,k,'.rds'))
        }
      }
    }
  }
  
  colnames(ks) <- c('ksel1','ksel2','ksel3','k')
  ks[,4] <- as.numeric(ks[,4])
  
  return (ks)
}

##########################################################################################

# For mnemk result analysis 

# extract legend from a graph
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}


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
  
  axis(1, at = 1:length(result$ics))
  mtext('k', side = 1, line = 2.5)
  legend('bottomright',legend=c('penalized','raw'), text.col = c('red','blue'), 
         pch = c(16,16), col = c('red','blue'), cex = 0.7)
  dev.off()
}

# inspect the responsibilities of cells and relationship to counts and log odds
check.sureness <- function (result, count.data, log.data, cluster = 0){
  # get responsibilities, k by l
  gamma <- getAffinity(x = result$probs, mw = result$mw, affinity = cluster)
  # get the maximum responsibility out of the ks for each cell
  # if cluster = 1 (hard clustering) it would be all 1s
  confidence <- apply(gamma, 2, max) 
  # inspect the cells that are not confident (pr. < 0.5) of being in any component
  unsure.cells.idx <- which(confidence < 0.5)
  sgene.included <- colnames(log.data)
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
      xlab('umap1') + ylab('umap2') + theme(plot.title = element_text(hjust = 0.5,size = 13))+theme_bw()
    
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
      xlab('umap1') + ylab('umap2') + theme(plot.title = element_text(hjust = 0.5,size = 13)) +theme_bw()+
      theme(axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20),
            axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20))
    
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
    ggtitle(title) + xlab('tsne1') + ylab('tsne2') + theme_bw() +
    theme(axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20),
          axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20))
  
  grid.arrange(p, mylegend, nrow = 1,widths = c(10,1))
}

graph.dim.reduction <- function (
  result = result,
  logodd.umap = logodd.umap,
  logodd.tsne = logodd.tsne,
  path = '.'
){

  soft <- getAffinity(x = result$probs, mw = result$mw, affinity = 0) 
  
  k <- dim(soft)[1]
  if (k == 2){
    # make fake graphs just to get the legends
    a.logodd <- data.frame(x=seq(2),y=seq(2))
    p.logodd <- ggplot(a.logodd)+geom_point(mapping=aes(x=x,y=y,col=c('red','blue')))+
      scale_color_manual(name='',values = c('red','blue'),labels=c('100% k=1','100% k=2')) +
      theme(legend.position = "bottom", legend.box = "horizontal",
            legend.key = element_rect(colour = "transparent", fill = "white"))
    
    mylegend.logodd <- g_legend(p.logodd)
    
    p1<- ggplot(logodd.umap$layout,aes('umap1','umap2')) +
      geom_point(mapping=aes(x=logodd.umap$layout[,1],y=logodd.umap$layout[,2]),
                 color=rgb(soft[1,],0,soft[2,],0.35), size = 0.1)+
      
      xlab('umap1') + ylab('umap2') +  theme_bw()+coord_fixed(ratio=1.1)+
      ggtitle("UMAP") + theme(plot.title = element_text(hjust = 0.5,size = 12),plot.margin=unit(c(1,1,-0.5,1), "cm"))
    
    
    p2 <- ggplot(logodd.tsne$Y,aes('tSNE1','tSNE2')) +
      geom_point(mapping=aes(x=logodd.tsne$Y[,1],y=logodd.tsne$Y[,2]),
                 color = rgb(soft[1,],0,soft[2,],0.35), size = 0.1)+
      xlab('tSNE1') + ylab('tSNE2') + theme_bw()+coord_fixed()+
      ggtitle("tSNE") + theme(plot.title = element_text(hjust = 0.5, size = 12),plot.margin=unit(c(1,1,-0.5,1), "cm"))
    
    #svg(filename = paste0(path,'/clustering on logodd.svg'))
    grid.arrange(arrangeGrob(p1,p2,nrow = 1), mylegend.logodd, nrow = 2,heights=c(15,1))
    #dev.off()
    
  }
  else{
    # make fake graphs just to get the legends
    a.logodd <- data.frame(x=seq(3),y=seq(3))
    p.logodd <- ggplot(a.logodd)+geom_point(mapping=aes(x=x,y=y,col=c('red','green','blue')))+
      scale_color_manual(name='',values = c('red','green','blue'),labels=c('100% k=1','100% k=2','100% k=3')) +
      theme(legend.position = "bottom", legend.box = "horizontal",
            legend.key = element_rect(colour = "transparent", fill = "white"))
    
    mylegend.logodd <- g_legend(p.logodd)
    
    p1<- ggplot(logodd.umap$layout,aes('umap1','umap2')) +
      geom_point(mapping=aes(x=logodd.umap$layout[,1],y=logodd.umap$layout[,2]),
                 color=rgb(soft[1,],soft[2,],soft[3,],0.35), size = 0.1)+
      
      xlab('umap1') + ylab('umap2') + theme_bw()+coord_fixed(ratio=1.1)+
      ggtitle("UMAP") + theme(plot.title = element_text(hjust = 0.5,size = 12),plot.margin=unit(c(1,1,-0.5,1), "cm"))
    
    
    p2 <- ggplot(logodd.tsne$Y,aes('tSNE1','tSNE2')) +
      geom_point(mapping=aes(x=logodd.tsne$Y[,1],y=logodd.tsne$Y[,2]),
                 color = rgb(soft[1,],soft[2,],soft[3,],0.35), size = 0.1)+
      xlab('tSNE1') + ylab('tSNE2') + theme_bw()+coord_fixed()+
      ggtitle("tSNE") + theme(plot.title = element_text(hjust = 0.5, size = 12),plot.margin=unit(c(1,1,-0.5,1), "cm"))
    
    #svg(filename = paste0(path,'/clustering on logodd.svg'))
    grid.arrange(arrangeGrob(p1,p2,nrow = 1), mylegend.logodd, nrow = 2,heights = c(15,1))
    #dev.off()
  }
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

  return(mnems)
}

# this function reads the mnem results for different parameter settings
# and summarises results into a matrix containing 
# the parameter setting itself, the log likelihood and its sd (if there are replicates),
# number of k and the molecular weights
make.summary <- function(path){
  result.summary<- NULL
  # get all the rds files, which should be the mnem result for each parameter setting
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

    # calculate the mean and standard deviation of the log likelihood at the each of each convergence 
    mean <- mean(unlist(lapply(seq_len(length(mnem.result$limits)), function (x) mnem.result$limits[[x]]$ll[length(mnem.result$limits[[x]]$ll)])))
    sd <- sd(unlist(lapply(seq_len(length(mnem.result$limits)), function (x) mnem.result$limits[[x]]$ll[length(mnem.result$limits[[x]]$ll)])))
    res.vec <- c(type,complete,ksel1,ksel2,ksel3,ll,mean,sd,k,mw)
    result.summary <- rbind(result.summary,res.vec)
  }
  colnames(result.summary) <- c('type','complete','clustering method',
                                'model selection criteria','distance metric',
                                'best loglikelihood','mean loglikelihood','sd loglikelihood','number of k','molecular weights')
  rownames(result.summary) <- NULL
  return (result.summary)
}

# plot the summary of mnem result 
# the graph compares various type, k and complete
plot.summary <- function(result.summary){
  result.summary <- as.data.frame(result.summary)
  types <- c('random','cluster','cluster2','cluster3','networks')
  # for k = 2 or 3, and complete = true or false, we build four vectors of length 
  # 5 = number of types
  # each vector at each position is the mean log likelihood under the particular type,
  # complete, and number of k
  k2T.mean <- unlist(lapply(types, 
                        function (x) 
                          as.numeric(result.summary[which(result.summary['type'] == x 
                                         & result.summary['number of k'] == 2 
                                         & result.summary['complete'] == 'TRUE'),
                                   'mean loglikelihood'])))
  
  k2T.max <- unlist(lapply(types, 
                            function (x) 
                              as.numeric(result.summary[which(result.summary['type'] == x 
                                                              & result.summary['number of k'] == 2 
                                                              & result.summary['complete'] == 'TRUE'),
                                                        'best loglikelihood'])))
  k2T.sd <- unlist(lapply(types, 
                          function (x) 
                            as.numeric(result.summary[which(result.summary['type'] == x 
                                          & result.summary['number of k'] == 2 
                                          & result.summary['complete'] == 'TRUE'),
                                   'sd loglikelihood'])))
  
  k2F.mean <- unlist(lapply(types, 
                        function (x) 
                          as.numeric(result.summary[which(result.summary['type'] == x 
                                         & result.summary['number of k'] == 2 
                                         & result.summary['complete'] == 'FALSE'),
                                   'mean loglikelihood'])))
  
  k2F.max <- unlist(lapply(types, 
                            function (x) 
                              as.numeric(result.summary[which(result.summary['type'] == x 
                                                              & result.summary['number of k'] == 2 
                                                              & result.summary['complete'] == 'FALSE'),
                                                        'best loglikelihood'])))
  k2F.sd <- unlist(lapply(types, 
                       function (x) 
                         as.numeric(result.summary[which(result.summary['type'] == x 
                                                         & result.summary['number of k'] == 2 
                                                         & result.summary['complete'] == 'FALSE'),
                                                   'sd loglikelihood'])))
  
  k3T.mean <- unlist(lapply(types, 
                        function (x) 
                          as.numeric(result.summary[which(result.summary['type'] == x 
                                         & result.summary['number of k'] == 3 
                                         & result.summary['complete'] == 'TRUE'),
                                   'mean loglikelihood'])))
  
  k3T.max <- unlist(lapply(types, 
                            function (x) 
                              as.numeric(result.summary[which(result.summary['type'] == x 
                                                              & result.summary['number of k'] == 3 
                                                              & result.summary['complete'] == 'TRUE'),
                                                        'best loglikelihood'])))
  k3T.sd <- unlist(lapply(types, 
                       function (x) 
                         as.numeric(result.summary[which(result.summary['type'] == x 
                                                         & result.summary['number of k'] == 3 
                                                         & result.summary['complete'] == 'TRUE'),
                                                   'sd loglikelihood'])))
  k3F.mean <- unlist(lapply(types, 
                       function (x) 
                         as.numeric(result.summary[which(result.summary['type'] == x 
                                        & result.summary['number of k'] == 3 
                                        & result.summary['complete'] == 'FALSE'),
                                  'mean loglikelihood'])))
  
  k3F.max <- unlist(lapply(types, 
                            function (x) 
                              as.numeric(result.summary[which(result.summary['type'] == x 
                                                              & result.summary['number of k'] == 3 
                                                              & result.summary['complete'] == 'FALSE'),
                                                        'best loglikelihood'])))
  k3F.sd <- unlist(lapply(types, 
                       function (x) 
                         as.numeric(result.summary[which(result.summary['type'] == x 
                                                         & result.summary['number of k'] == 3 
                                                         & result.summary['complete'] == 'FALSE'),
                                                   'sd loglikelihood'])))
  
   ggplot() + 
    geom_line(aes(x=types,y=k3T.mean,col='k = 3', linetype='expectation'), group='k3Tl', size = 1) +
    geom_point(aes(x=types,y=k3T.mean,shape='expectation mean',col='k = 3'), group='k3Ta',  size = 3) +
    geom_point(aes(x=types,y=k3T.max,shape='expectation max',col='k = 3'), group='k3Tm',  size = 3) +
    geom_errorbar(aes(x=types,ymin=k3T.mean-k3T.sd, ymax=k3T.mean+k3T.sd,col='k = 3'),width=0.2,size=0.7)+
    
    geom_line(aes(x=types,y=k3F.mean,col='k = 3', linetype='observed'), group='k3Fl',  size = 1) + 
    geom_point(aes(x=types,y=k3F.mean,shape='observed mean',col='k = 3'), group='k3Fa', size = 3) +
    geom_point(aes(x=types,y=k3F.max,shape='observed max',col='k = 3'), group='k3Fm',  size = 3) +
    geom_errorbar(aes(x=types,ymin=k3F.mean-k3F.sd, ymax=k3F.mean+k3F.sd,col='k = 3'),width=0.2,size=0.7)+
  
    geom_line(aes(x=types,y=k2T.mean,col='k = 2',linetype='expectation'), group='k2Tl',  size = 1) +
    geom_point(aes(x=types,y=k2T.mean,shape='expectation mean',col='k = 2'), group='k2Ta',  size = 3) +
    geom_point(aes(x=types,y=k2T.max,shape='expectation max',col='k = 2'), group='k2Tm',  size = 3) +
    geom_errorbar(aes(x=types,ymin=k2T.mean-k2T.sd, ymax=k2T.mean+k2T.sd,col='k = 2'),width=0.2,size=0.7)+
     
    geom_line(aes(x=types,y=k2F.mean,col='k = 2', linetype='observed'),group='k2Fl',  size = 1) +
    geom_point(aes(x=types,y=k2F.mean,shape='observed mean',col='k = 2'), group='k2Fa',  size = 3) +
    geom_point(aes(x=types,y=k2F.max,shape='observed max',col='k = 2'), group='k2Fm',  size = 3) +
    geom_errorbar(aes(x=types,ymin=k2F.mean-k2F.sd, ymax=k2F.mean+k2F.sd,col='k = 2'),width=0.2,size=0.7)+
     
    xlab('') + ylab('log likelihood') +
    theme_bw() + 
    scale_shape_manual(name='',values=c(2,1,17,16),labels=c('expectation max','expectation mean','observed max','observed mean'))+
    scale_color_manual(name='',values=c('red','blue'),labels=c('k=2','k=3')) + 
    scale_linetype_manual(name='',values=c('dotted','solid'),labels=c('expectation','observed')) +
    theme(legend.position = 'right') + theme(axis.text=element_text(size=14,face='bold'),
                                              axis.title.y=element_text(size=14,face='bold'))
    

}

# compare mnemk (type = 'netowrks', complete = F) with the same setting using
# mnem. Both mnemk and mnem should have 10 replications 
# It takes the mnemk result and the path that stores the mnem result 
# and output the graph
plot.mnemk.mnem.comparison <- function(mnemk.F,path){
  # read mnem files with k=2 (when using ksel = c('kmeans',*,'euclidean'))
  # where * is AIC/BIC/silhouette, and extract the log likelihoods for the
  # 10 replications 
  aic.file <- readRDS(paste0(path,'/networks_FALSE_kmeans_AIC_euclidean.rds'))
  bic.file <- readRDS(paste0(path,'/networks_FALSE_kmeans_BIC_euclidean.rds'))
  sil.file <- readRDS(paste0(path,'/networks_FALSE_kmeans_silhouette_euclidean.rds'))
  
  aic.values <- unlist(lapply(c(1:10), 
                            function(x) aic.file$limits[[x]]$ll[length(aic.file$limits[[x]]$ll)]))
  bic.values <- unlist(lapply(c(1:10), 
                              function(x) bic.file$limits[[x]]$ll[length(bic.file$limits[[x]]$ll)]))
  sil.values <- unlist(lapply(c(1:10), 
                              function(x) sil.file$limits[[x]]$ll[length(sil.file$limits[[x]]$ll)]))
  mnemk.values <- unlist(lapply(c(1:10), 
                              function(x) mnemk.F$best$limits[[x]]$ll[length(mnemk.F$best$limits[[x]]$ll)]))
  
  types <- c('k-AIC-e','k-BIC-e','k-silhouette-e','k=2')
  ggplot() + geom_boxplot(aes(x=types[1],y=aic.values)) +
    geom_boxplot(aes(x=types[2],y=bic.values)) +
    geom_boxplot(aes(x=types[3],y=sil.values)) +
    geom_boxplot(aes(x=types[4],y=mnemk.values)) + xlab('') + ylab('log likelihood')+
    ylim(470000,474500)+
    theme_bw()+ theme(axis.text=element_text(size=12,face='bold'),
                      axis.title.y=element_text(size=14,face='bold')) 
  
  
  # organize data into a data frame
  ll <- c(aic.values,bic.values,sil.values,mnemk.values)
  group <- unlist(lapply(types, function (x) rep(x,10)))
  my.data <- data.frame(ll,group)
  # check outliners 
  outliner.test <- my.data %>% group_by(group) %>% identify_outliers(ll)
  outliners <- TRUE %in% outliner.test[,4]
  # check normality
  model  <- lm(ll ~ group, data = my.data)
  normality.test <- shapiro_test(residuals(model))
  # check variance homogeneity
  homogeneity.test <- my.data %>% levene_test(ll ~ group)
  # if the data has no extreme outliners, is normally distributed and has
  # homogeneous variance, then we do one-wau ANOVA to see if any group
  # has a statistically significantly differnet log likelihood 
  if (!outliners && (normality.test['p.value'] > 0.05) && homogeneity.test['p'] > 0.05){
    res.aov <- aov(ll ~ group, data = my.data)
    return (list(outliner=outliner.test,normality=normality.test,homogeneity=homogeneity.test,aov=summary(res.aov)))
  }
}

# do wilcox test to compare different methods 
test.methods<-function(path){
  data <- NULL
  files <- list.files(path = path, pattern = '.rds')
  for (file in files){
    mnem.result <- readRDS(paste0(path,'/',file))
    parameters <- strsplit(strsplit(file,'.rds')[[1]],'_')[[1]]
    type <- parameters[1]
    complete <- parameters[2]
    k <- length(mnem.result$mw)
    replication <- length(mnem.result$limits)
    lls <- unlist(lapply(seq_len(length(mnem.result$limits)), 
                         function (x) mnem.result$limits[[x]]$ll[length(mnem.result$limits[[x]]$ll)]))
    res.matrix <- cbind(lls,rep(type,replication),rep(complete,replication),rep(k,replication))
    data <- rbind(data,res.matrix)
  }
  colnames(data) <- c('ll','type','complete','k')
  # compare complete = T/F (observed/expectation)
  complete <- wilcox.test(as.numeric(data[data[,'complete'] == 'TRUE','ll']),
                          as.numeric(data[data[,'complete'] == 'FALSE','ll']))$p.value
  # compare k=2/3 under complete = F ('observed')
  #k <- wilcox.test(as.numeric(data[data[,'complete'] == 'FALSE' & data[,'k'] == '2','ll']),
                   #as.numeric(data[data[,'complete'] == 'FALSE' & data[,'k'] == '3','ll']))$p.value
  # compare k=2/3 under cluster, networks, random and 'observed'
  cluster <- wilcox.test(as.numeric(data[data[,'complete'] == 'FALSE' & data[,'k'] == '2' & data[,'type'] == 'cluster','ll']),
                                  as.numeric(data[data[,'complete'] == 'FALSE' & data[,'k'] == '3' & data[,'type'] == 'cluster','ll']))$p.value
  networks <- wilcox.test(as.numeric(data[data[,'complete'] == 'FALSE' & data[,'k'] == '2' & data[,'type'] == 'networks','ll']),
                                 as.numeric(data[data[,'complete'] == 'FALSE' & data[,'k'] == '3' & data[,'type'] == 'networks','ll']))$p.value
  random <- wilcox.test(as.numeric(data[data[,'complete'] == 'FALSE' & data[,'k'] == '2' & data[,'type'] == 'random','ll']),
                                   as.numeric(data[data[,'complete'] == 'FALSE' & data[,'k'] == '3' & data[,'type'] == 'random','ll']))$p.value
  networks_random <- wilcox.test(as.numeric(data[data[,'complete'] == 'FALSE' & data[,'k'] == '2' & data[,'type'] == 'networks','ll']),
                                 as.numeric(data[data[,'complete'] == 'FALSE' & data[,'k'] == '2' & data[,'type'] == 'random','ll']))$p.value
  networks_cluster <- wilcox.test(as.numeric(data[data[,'complete'] == 'FALSE' & data[,'k'] == '2' & data[,'type'] == 'networks','ll']),
                                  as.numeric(data[data[,'complete'] == 'FALSE' & data[,'k'] == '2' & data[,'type'] == 'cluster','ll']))$p.value
  return(list(complete=complete,k=k,cluster.networks=cluster.networks,
              networks.random=networks.random,random.cluster=random.cluster))
}




