
superTargetDTI <- function(Y, K, cutoffSuper) {
  ### INPUT:
  ## Y: matrix, 0-1 interaction adjacency matrix
  ## K: matrix, kenel matrix
  ## cutoffSuper: scalar
  
  ### OUTPUT:
  
  ## agglomerative hierarchical clustering
  
  ## convert to 'dist' matrix
  d <- 1 - K
  ## match to matlab in SuperTarget code  
  d <- d - diag(diag(d))
  
  ## convert to 'dist' object
  d <- as.dist(d)
  
  ## cluster package
  agn <- cluster::agnes(d, method = "ward")

  agn <- as.hclust(agn)

  C <- cutree(agn, h = 1.1)

  uC <- unique(C)
  
  ## number of cluster
  nC <- length(uC)
  
  Yshrink <- Y
  ## YY: used for do.call()
  YY <- as.data.frame(Y)
  
  for (i in 1:nC) {
    curCluster <- which(C %in% i)
    ## do.call: pairwise max for a data.frame
    ## do.call: second parameter is data.frame
    Yshrink[, curCluster] <- do.call(pmax.int, YY[, curCluster])
  }
  
  return(Yshrink)
}


