
inferZeros <- function(inMat, simRow, K = NULL) {
  # INPUT
  # inMat: input matrix
  # simRow: similarity matrix between rows
  # K: number of neighbor used to infer zeros

  # OUTPUT
  # a matrix of complete input matrix
  
  if (!is.matrix(inMat)) {
    stop("inMat MUST be matrix class! \n")
  }
  
  if (!is.matrix(simRow)) {
    stop("simRow MUST be matrix class! \n")
  }
  
  
  numRow <- nrow(inMat)
  numCol <- ncol(inMat)
  
  if (is.null(K)) {
    K <- numCol
  }  
  
  indexZeros <- which(apply(inMat, 1, sum) == 0)
  numIndexZeros <- length(indexZeros)
  
  if (numIndexZeros > 0) {
    simRow[, indexZeros] <- 0
  }
  
  simRow <- simRow - diag(diag(simRow))
  
  if (numIndexZeros > 0) {
    # cat("Inferring zeros in the inMat! \n")
    # flush.console()
    for (i in indexZeros) {
      currSimForZeros <- simRow[i, ]
      indexRank <- rank(currSimForZeros)
      indexNeig <- which(indexRank > (numRow - K))
      simCurr <- currSimForZeros[indexNeig]
      inMat_Known <- inMat[indexNeig, , drop = FALSE]
      inMat[i, ] <- simCurr %*% inMat_Known / sum(simCurr)
    }
  }
  return(inMat)
}