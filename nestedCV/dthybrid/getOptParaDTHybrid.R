
getOptParaDTHybrid <- function(Yfold, innerFold = 5, numSplitInner = 1,
           paraList = list(lambda = seq(0.2, 0.9, by = 0.1),
                           alpha = seq(0.2, 0.3, by = 0.1)),
           sd,
           st) {
    
  ##
  Folds <- doCVPositiveOnly3(Yfold, kfold = innerFold, numSplit = numSplitInner)
  
  lenL <- length(paraList)
  para <- expand.grid(paraList$lambda, paraList$alpha)
  colnames(para) <- c("lambda", "alpha")
  nCombin <- nrow(para)
  
  res <- matrix(NA, nrow = nCombin, ncol = lenL + 1)
  colnames(res) <- c("lambda", "alpha", paste0("mpr", innerFold, "CV"))
  
  mprTotal <- matrix(NA, nrow = innerFold, ncol = 1)
  colnames(mprTotal) <- "mprValue"
  
  for (jj in 1:nCombin) {
    lambda <- para[jj, "lambda"]
    alpha <- para[jj, "alpha"]
    for (i in 1:innerFold) {
      cat("...inner loop...", "lambda =", lambda, "; alpha =", alpha, "; fold:", i,
          "/", innerFold, "...doing nestedCV...", "\n")
      flush.console()
      
      Yfold <- Folds[[1]][[i]][["foldMat"]]
      testSet <- Folds[[1]][[i]][["testSet"]]
      
      Ypred <- computeRecommendation(A = Yfold, lambda = lambda, alpha = alpha,
                                     S = sd, S1 = st)
      result <- evalMetrics(Ypred = Ypred, testSet = testSet)
      mprTotal[i, 1] <- result[, "MPR"]
    }
    res[jj, "lambda"] <- lambda
    res[jj, "alpha"] <- alpha
    res[jj, "mpr5CV"] <- mean(mprTotal[, "mprValue"])
  }
  
  cutValue <- min(res[, "mpr5CV"])
  isBest <- res[, "mpr5CV"] == cutValue
  bestRes <- res[isBest, ]
  ## if multiple results, then just need 1st one
  if (!is.vector(bestRes)) {
    bestRes <- bestRes[1, ]
  }
  
  return(bestRes)
}









