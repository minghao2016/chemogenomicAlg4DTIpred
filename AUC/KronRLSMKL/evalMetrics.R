
evalMetrics <- function(Ypred,
                        testSet
                        ) {

  ##     MPR
  percentileRank <- function(x) {
    ### INPUT:
    ## x: vector, unsorted predicted scores for user i
    
    ### OUTPUT:
    ## pr: vector, with percentile rank
    
    rx <- rle(sort(x))
    smaller <- cumsum(c(0, rx$lengths))[seq(length(rx$lengths))]
    larger <- rev(cumsum(c(0, rev(rx$lengths))))[-1]
    rxpr <- smaller / (smaller + larger)
    res <- rxpr[match(x, rx$values)]
    pr <- 1 - res
    names(pr) <- names(x)
    return(pr)
  }
  
  testedCol <- sort(unique(testSet[, 2]))
  MPR <- 0
  
  for (i in testedCol) {
    ypredi <- Ypred[, i]
    ## percentile ranking result: prr
    prr <- percentileRank(ypredi)
    idxTest <- testSet[testSet[, 2] == i, 1]
    MPR <- MPR + mean(prr[idxTest])
  }
  MPR <- MPR / length(testedCol)
  
  ## result
  metrics <- matrix(NA, nrow = 1, ncol = 1)
  colnames(metrics) <- c("MPR")
  metrics <- as.data.frame(metrics)
  metrics[1, ] <- MPR
  
  return(metrics)
}





