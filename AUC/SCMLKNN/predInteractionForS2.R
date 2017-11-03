
predInteractionForS2 <- function(Yfold, simD, simT) {
  ### INPUT:
  ## Yfold: matrix, drug-target interaction matrix
  ## simD: list, for drug similarity
  ## simT: list, for target similarity
  
  ### OUTPUT:
  ##
  
  ## Reduce() handles components of one List
  nDrug <- length(simD)
  simD <- Map("*", 1/nDrug, simD)
  Kd <- Reduce("+", simD)
  
  ## Map() handles two List
  nTarget <- length(simT)
  simT <- Map("*", 1/nTarget, simT)
  Kt <- Reduce("+", simT)
  

  cutoffSuper <- 1.1
  
  ## super target Y
  # source("superTargetDTI.R")
  Yst <- superTargetDTI(Y = Yfold, K = Kt, cutoffSuper = cutoffSuper)
  
  ## transpose Yfold
  YfoldT <- t(Yfold)
  
  # source("MLKNNtrain.R")
  res <- MLKNNtrain(trainSim = Kd, trainY = YfoldT, K = 3, smooth = 1)
  Prior <- res$Prior
  PriorN <- res$PriorN
  Cond <- res$Cond
  CondN <- res$CondN
  # source("MLKNNtest.R")
  YpredOrg <- MLKNNtest(trainY = YfoldT, testSim = Kd, testY = YfoldT, K = 3,
                      Prior = Prior, PriorN = PriorN, Cond = Cond, CondN = CondN)
  ## output 1
  YpredOrg <- t(YpredOrg)
  
  ## super target
  YstT <- t(Yst)
  ## prediction based on Yst
  res <- MLKNNtrain(trainSim = Kd, trainY = YstT, K = 3, smooth = 1)
  Prior <- res$Prior
  PriorN <- res$PriorN
  Cond <- res$Cond
  CondN <- res$CondN
  YpredSt <- MLKNNtest(trainY = YstT, testSim = Kd, testY = YstT, K = 3,
                        Prior = Prior, PriorN = PriorN, Cond = Cond, CondN = CondN)
  ## output 2
  YpredSt <- t(YpredSt)
  
  ## final results
  Ypred <- YpredOrg * YpredSt
  
  return(Ypred)
}



