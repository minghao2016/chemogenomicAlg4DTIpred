

blm <- function(Yfold = Yfold,
                Kd = Kd,
                Kt = Kt,
                cc = 1) {
  ### INPUT:
  ##
  
  
  ### OUTPUT:
  ##
  
  
  
  ## prediction by target
  ## number of target
  nDrug <- nrow(Yfold)
  nTarget <- ncol(Yfold)
  YpredByTarget <- matrix(NA, nrow = nDrug, ncol = nTarget)

  ## pre-computing
  trainK <- as.kernelMatrix(Kd)  
  
  for (i in 1:nTarget) {
    ###################### change to factor  #################
    curY <- factor(Yfold[, i])

    model <- ksvm(trainK, curY, kernel = "matrix", C = cc, cross = 0)

    testK <- as.kernelMatrix(Kd[, SVindex(model), drop = F])
    YpredByTarget[, i] <- predict(model, testK, type = "decision")
  }
  
  ## prediction by drug
  YfoldT <- t(Yfold)
  nTarget <- nrow(YfoldT)
  nDrug <- ncol(YfoldT)
  YpredByDrug <- matrix(NA, nrow = nTarget, ncol = nDrug)
  
  ## pre-computing
  trainK <- as.kernelMatrix(Kt)
  
  for (i in 1:nDrug) {
    ##############################################
    curY <- factor(YfoldT[, i]) ## must factor()
    ##############################################
    model <- ksvm(trainK, curY, kernel = "matrix", C = cc, cross = 0)
    testK <- as.kernelMatrix(Kt[, SVindex(model), drop = F])
    YpredByDrug[, i] <- predict(model, testK, type = "decision")
  }
  
  ## get maximum
  Ypred <- pmax(YpredByTarget, t(YpredByDrug))
  
  return(Ypred)
}








