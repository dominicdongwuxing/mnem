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

