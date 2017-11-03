




calcPredScoreAUCAUPR <- function(U,
                                 V,
                                 simDrug,
                                 simTarget,
                                 knownDrugIndex,
                                 knownTargetIndex,
                                 testIndexRow,
                                 testIndexCol,
                                 K = 5,
                                 thisAlpha,
                                 thisBeta,
                                 thisGamma,
                                 testSet,
                                 Y) {
  
  
  
  if (K < 0) {
    stop("K MUST be '>=' 0! \n")
  }
  
  if (K > 0) {
    ## cat("with K smoothing! \n")
    ## for drug
    indexTestD <- unique(testIndexRow)
    testD <- U[indexTestD, ]
    testD <- cbind(indexTestD, testD)
    numTest <- length(indexTestD)
    numColTestD <- ncol(testD)
    simDrugKnown <- simDrug[, knownDrugIndex]
    numDrugKnown <- length(knownDrugIndex)
    
    for (i in 1:numTest) {
      indexCurr <- indexTestD[i]
      isNewDrug <- !(indexCurr %in% knownDrugIndex)
      if (isNewDrug) {
        cat("Have New Drug", "\n")
        flush.console()
        simDrugNew <- simDrugKnown[indexCurr, ] # vector
        indexRank <- rank(simDrugNew) # vector
        indexNeig <- which(indexRank > (numDrugKnown - K))
        simCurr <- simDrugNew[indexNeig] # vector
        # index for U
        index4U <- knownDrugIndex[indexNeig]
        U_Known <- U[index4U, , drop = FALSE] # force to matrix
        # vec %*% matrix => matrix
        testD[i, 2:numColTestD] <- (simCurr %*% U_Known) / sum(simCurr)
      }
    }
    
    Unew <- U
    Unew[indexTestD, ] <- testD[, -1]
    
    ## for target
    # unique index for test target
    indexTestT <- unique(testIndexCol)
    testT <- V[indexTestT, ]
    # add first column as labels
    testT <- cbind(indexTestT, testT) # 1st column is unique test label
    # number of unique test set
    numTest <- length(indexTestT)
    # number of column for testT
    numColTestT <- ncol(testT)
    # known similarity matrix for targets
    simTargetKnown <- simTarget[, knownTargetIndex]
    # number of known targets
    numTargetKnown <- length(knownTargetIndex)
    
    for (i in 1:numTest) {
      indexCurr <- indexTestT[i]
      isNewTarget <- !(indexCurr %in% knownTargetIndex)
      if (isNewTarget) {
        cat("Have New Target", "\n")
        flush.console()
        simTargetNew <- simTargetKnown[indexCurr, ]
        indexRank <- rank(simTargetNew) 
        indexNeig <- which(indexRank > (numTargetKnown - K))
        simCurr <- simTargetNew[indexNeig] 
        index4V <- knownTargetIndex[indexNeig]
        V_Known <- V[index4V, , drop = FALSE]
        testT[i, 2:numColTestT] <- (simCurr %*% V_Known) / sum(simCurr)
      }
    }

    Vnew <- V
    Vnew[indexTestT, ] <- testT[, -1]
    
    Vnewt <- t(Vnew)
    UnewVnewt <- Unew %*% Vnewt
    val <- thisAlpha * UnewVnewt +
      thisBeta * (simDrug %*% UnewVnewt) + thisGamma * (UnewVnewt %*% simTarget)
    Ypred <- sigmoid(val)
    
    ## for mpr
    ## result <- evalMetrics(Ypred = Ypred, testSet = testSet)
    
    ## calculate auc and aupr
    
    ## columns in the test set
    # colInTest <- unique(testSet[, "colIndex"])
    # yLabel <- Y[, colInTest]
    # yScore <- Ypred[, colInTest]
    # 
    # yLabel <- as.vector(yLabel)
    # yScore <- as.vector(yScore)
    # result <- calAUPR(yLabel, yScore)
	
  } else {  
    
    Vt <- t(V)
    UVt <- U %*% Vt
    val <- thisAlpha * UVt + thisBeta * (simDrug %*% UVt) + 
      thisGamma * (UVt %*% simTarget) 
    Ypred <- sigmoid(val)
    
    ## result <- evalMetrics(Ypred = Ypred, testSet = testSet)
    
    ## calculate auc and aupr
    
    ## columns in the test set
    # colInTest <- unique(testSet[, "colIndex"])
    # yLabel <- Y[, colInTest]
    # yScore <- Ypred[, colInTest]
    # 
    # yLabel <- as.vector(yLabel)
    # yScore <- as.vector(yScore)
    # result <- calAUPR(yLabel, yScore)
  }
  return(Ypred)
}




