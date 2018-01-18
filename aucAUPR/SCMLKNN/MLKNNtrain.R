
MLKNNtrain <- function(trainSim, trainY, K = 3, smooth = 1) {

  
  numClass <- nrow(trainY)
  numTrain <- ncol(trainY)
  
  ## compute Prior and PriorN
  tempCi <- rowSums(trainY == 1)
  Prior <- (smooth + tempCi) / (smooth * 2 + numTrain)
  PriorN <- 1 - Prior
  

  ## matrix to List by row
  L <- split(trainSim, row(trainSim)) 
  ## neighbor positions, top K
  ## order(), decreasing = TRUE
  neig <- Map(order, L, decreasing = TRUE)
  neig <- lapply(X = neig, "[", 1:K)
  
  tempCi <- matrix(0, nrow = numClass, ncol = K + 1)
  tempNCi <- tempCi
  
  
  for (i in 1:numTrain) {
    ## neighbor labels for current instance
    neigLabel <- trainY[, neig[[i]]]
    ## number of neighbors with label == 1
    temp <- rowSums(neigLabel == 1)
    
    isOne <- (trainY[, i] == 1)
    ## use: mat[mat] to access the elements of a matrix, and modify the matrix elements
    curEleOne <- cbind(which(isOne), temp[isOne] + 1)
    tempCi[curEleOne] <- tempCi[curEleOne] + 1
    
    isZero <- !isOne
    curEleZero <- cbind(which(isZero), temp[isZero] + 1)
    tempNCi[curEleZero] <- tempNCi[curEleZero] + 1
  }
  
  temp1 <- rowSums(tempCi)
  temp2 <- rowSums(tempNCi)
  
  Cond <- (smooth + tempCi) / (smooth * (K + 1) + temp1)
  CondN <- (smooth + tempNCi) / (smooth * (K + 1) + temp2)
  
  res <- list(Prior = Prior, PriorN = PriorN,
              Cond = Cond, CondN = CondN)
  
  return(res)
}


