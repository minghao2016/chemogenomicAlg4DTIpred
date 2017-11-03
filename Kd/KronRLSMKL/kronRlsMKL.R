
kronRlsMKL <- function(K1, K2, Yfold, lmd, sgm = NULL, maxiter = NULL) {
  ## INPUT:
  # K1: list, multiple kernels
  # K1: list, multiple kernels
  # Yfold: matrix, training set
  # lmd: scalar, for RLS regularized coefficient
  # sgm: scalar, for 'weight' regularized coefficient
  # maxiter: scalar, maximum iteraction number
  
  ## OUTPUT:
  
  
  if (is.null(maxiter)) {
    maxiter <- 20
  }
  
  if (is.null(sgm)) {
    sgm <- 0.2
  }
  
  ## number of kernels
  nA <- length(K1)
  nB <- length(K2)
  
  ## initialize kernel weights (uniform)
  alph <- rep(1 / nA, times = nA)
  bta <- rep(1 / nB, times = nB)
  
  iter <- 0
  incr <- 1000
  limitIncr <- 0.01
  combFval <- rep(0, times = maxiter)
  
  ## iterative steps: optimize weights
  while ((iter < maxiter) && (abs(incr) > limitIncr)) {
    
    iter <- iter + 1
    
 
    K1comb <- combineKernels(alph, K1)
    K2comb <- combineKernels(bta, K2)
    

    A <- kronRlsC(K1comb, K2comb, Yfold, lmd)
    

    u <- (Yfold - lmd * A / 2)
    

    Ma <- vector("list", length = nA)
    for (i in 1:nA) {
      ## fix K2, and compute K1 without K1 weights
      ## $Ma \in R^{n \times m}$, while $Yfold \in R^{m \times n}
      Ma[[i]] <- t(K2comb) %*% t(A) %*% K1[[i]]
    }
    
    ## optimal alpha
    # source("optWeights.R")
    ## nonlinear programming
    NP <- optWeights(w0 = alph, u = u, Ma = Ma, lmd = lmd, sgm = sgm)
    xAlpha <- NP$optW
    fvalueAlpha <- NP$optV
    

    Mb <- vector("list", length = nB)
    for (i in 1:nB) {
      Mb[[i]] <- t(K2[[i]]) %*% t(A) %*% K1comb
    }
    
    ## optimal beta
    NP <- optWeights(w0 = bta, u = u, Ma = Mb, lmd =lmd, sgm = sgm)
    xBeta <- NP$optW
    fvalueBeta <- NP$optV
    
    ## update values for next iteration
    alph <- xAlpha
    bta <- xBeta
    
    combFval[iter] <- fvalueAlpha + fvalueBeta
    
    if (iter > 1) {
      incr <- combFval[iter - 1] - combFval[iter]
    }
  }
  ## perform prediction with optimal weights
  K1comb <- combineKernels(alph, K1)
  K2comb <- combineKernels(bta, K2)
  

  Yhat <- kronRls(K1comb, K2comb, Yfold, lmd)
  
  res <- list(Yhat = Yhat, alph = alph, bta = bta)
  
  return(res)
}





