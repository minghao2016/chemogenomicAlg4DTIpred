



calcPredScore <- function(Y,
                          U,
                          V,
                          simDrug,
                          simTarget,
                          thisAlpha,
                          thisBeta,
                          thisGamma,
                          testSet) {
  
  
  
  
  Vt <- t(V)
  UVt <- U %*% Vt
  val <- thisAlpha * UVt + thisBeta * (simDrug %*% UVt) + thisGamma * (UVt %*% simTarget)
  Ypred <- sigmoid(val)
  
  testLabel <- Y[testSet]
  score <- Ypred[testSet]
  
  result <- calAUPR(testLabel, score)
  
  return(result)
}




