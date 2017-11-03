
MLKNNtest <- function(trainY, testSim, testY, K = 3, Prior, PriorN, Cond, CondN) {
  ### INPUT:
  ## trainY: 
  
  ### OUTPUT:
  ##
  
  
  numClass <- nrow(trainY)
  numTrain <- ncol(trainY)
  numTest <- ncol(testY)
  
  ## neighbor positions of 'test set' data
  ## matrix to List by row
  L <- split(testSim, row(testSim)) ## good! row(X)
  ## neighbor positions, top K
  ## order(), decreasing = TRUE
  neig <- Map(order, L, decreasing = TRUE)
  neig <- lapply(X = neig, "[", 1:K)
  
  ## compute outputs
  Outputs <- matrix(0, numClass, numTest)
  
  for (i in 1:numTest) {
    ## neighbor labels for current instance
    neigLabel <- trainY[, neig[[i]]]
    
    ## number of neighbors with label == 1
    temp <- rowSums(neigLabel == 1)
    
    ## probability in and out 
    ProbIn <- Prior * Cond[cbind(1:numClass, temp + 1)]
    ProbOut <- PriorN * CondN[cbind(1:numClass, temp + 1)]
    
    P2 <- ProbIn + ProbOut
    
    ## are all zeros for 'ProbIn + ProbOut'?
    isZero <- (P2 == 0)
    
    isAllZero <- (sum(isZero) == length(isZero))
    isNotAllZero <- (sum(isZero) == 0)
    
    if (isAllZero) {
      Outputs[, i] <- Prior
    } else if (isNotAllZero) {
      Outputs[, i] <- ProbIn / P2
    } else {
      ## zero position
      zeroPos <- which(isZero)
      Outputs[zeroPos, i] <- Prior[zeroPos]
      
      ## not zero position
      notZeroPos <- which(!isZero)
      Outputs[notZeroPos, i] <- ProbIn[notZeroPos] / P2[notZeroPos]
    }
  }
  
  return(Outputs)
}








