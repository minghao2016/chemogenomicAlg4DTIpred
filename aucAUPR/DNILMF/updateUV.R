updateUV <- function(cc = 5,
                     inMat,
                     thisAlpha = 0.8,
                     thisBeta = 0.1,
                     Sd,
                     thisGamma = NULL,
                     St,
                     lamU,
                     lamV,
                     numLat = 50,
                     initMethod = "useNorm",
                     thisSeed = 123,
                     maxIter = 100) {
  
  # INPUT
  # cc:
  # inMat:
  # thisAlpha:
  # thisBeta:
  # Sd:
  # thisGamma:
  # St:
  # lamU:
  # lamV:
  # numLat:
  # initMethod:
  # thisSeed:
  # maxIter:
  
  # OUTPUT:
  # a list with two elements: U and V
  
  
  if ((thisAlpha > 1) | (thisAlpha < 0)) {
    stop("thisAlpha should be [0, 1]! \n")
  }
  
  if ((thisBeta > 1) | (thisBeta < 0)) {
    stop("thisBeta should be [0, 1]! \n")
  }
  
  if (is.null(thisGamma)) {
    thisGamma <- 1 - thisAlpha - thisBeta
  }
  
  numRow <- nrow(inMat)
  numCol <- ncol(inMat)
  
  if (initMethod == "useNorm") {
    U <- matrix(NA, nrow = numRow, ncol = numLat)
    U <- apply(U, 2, function(x) {
      rnorm(x, mean = 0, sd = 1)
    })
    U <- 1 / numLat * U
    
    V <- matrix(NA, nrow = numCol, ncol = numLat)
    V <- apply(V, 2, function(x) {
      rnorm(x, mean = 0, sd = 1)
    })
    V <- 1 / numLat * V
    
  } else if (initMethod == "useSeed") {
    set.seed(thisSeed)
    U <- matrix(NA, nrow = numRow, ncol = numLat)
    U <- apply(U, 2, function(x) {
      rnorm(x, mean = 0, sd = 1)
    })
    U <- 1 / numLat * U

    
    V <- matrix(NA, nrow = numCol, ncol = numLat)
    V <- apply(V, 2, function(x) {
      rnorm(x, mean = 0, sd = 1)
    })
    V <- 1 / numLat * V

    
  } else {
    stop("initMethod should be one of {useNorm, useSeed}\n")  
  }
  
  sumGradU <- matrix(0, nrow = numRow, ncol = numLat)
  sumGradV <- matrix(0, nrow = numCol, ncol = numLat)

  # last log-likelihood
  lastLog <- calcLogLik(
      cc = cc,
      inMat = inMat,
      thisAlpha = thisAlpha,
      U = U,
      V = V,
      thisBeta = thisBeta,
      Sd = Sd,
      thisGamma = thisGamma,
      St = St,
      lamU = lamU,
      lamV = lamV)
  
  currDeltaLL <- 1000
  # main loop
  for (i in 1:maxIter) {
    # gradU
    gradU <- calcDeriv(
      cc = cc,
      inMat = inMat,
      thisAlpha = thisAlpha,
      U = U,
      V = V,
      thisBeta = thisBeta,
      Sd = Sd,
      thisGamma = thisGamma,
      St = St,
      lamU = lamU,
      lamV = lamV,
      isGradU = TRUE)
    sumGradU <- sumGradU + (gradU ^ 2)
    stepSize <- 1 / sqrt(sumGradU)
    U <- U + stepSize * gradU
    
    # gradV
    gradV <- calcDeriv(
      cc = cc,
      inMat = inMat,
      thisAlpha = thisAlpha,
      U = U,
      V = V,
      thisBeta = thisBeta,
      Sd = Sd,
      thisGamma = thisGamma,
      St = St,
      lamU = lamU,
      lamV = lamV,
      isGradU = FALSE
    )
    sumGradV <- sumGradV + (gradV ^ 2)
    stepSize <- 1 / sqrt(sumGradV)
    V <- V + stepSize * gradV
    
    currLog <- calcLogLik(
      cc = cc,
      inMat = inMat,
      thisAlpha = thisAlpha,
      U = U,
      V = V,
      thisBeta = thisBeta,
      Sd = Sd,
      thisGamma = thisGamma,
      St = St,
      lamU = lamU,
      lamV = lamV)
    
    # delta log-likelihood
    deltaLog <- (currLog - lastLog) / abs(lastLog)
    
    # stop earlier
    if (abs(deltaLog) < 1e-5) {
      break
    }
    
    if ((i > 50) & (deltaLog > currDeltaLL)) {
      break
    }
    
    currDeltaLL <- deltaLog
    lastLog <- currLog
  }
  
  UV <- list(U = U, V = V)
  return(UV)
}
