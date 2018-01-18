

optWeights <- function(w0, u, Ma, lmd, sgm) {
  ## INPUT:
  # w0: vector, initial weight to be optimized
  # u: matrix
  # Ma: matrix, dimension equal to t(Yfold)
  # lmd: scalar, regularized coefficient of RLS
  # sgm: scalar, regularized coefficient of weight
  
  ## OUTPUT:
  
  # library(nloptr)
  
  ## creat objective function
  fn <- function(w, u, Ma, lmd, sgm) {
    ## INPUT:
    # w: vector, weight to be optimized
    # u: matrix 
    # Ma: matrix
    # lmd: scalar
    # sgm: scalar
    
    ## OUTPUT:
    
    # source("combineKernels.R")
    V <- combineKernels(w, Ma)
    num <- nrow(V) * ncol(V)
    J <- norm(u - t(V), "F") / (2 * lmd * num) + sgm * (norm(w, "2") ^ 2)
    return(J)
  }
  
  n <- length(w0)
  LB <- rep(0, times = n)
  UB <- rep(1, times = n)
  
  ## sum(w) = 1
  heq <- function(w) {
    return(sum(w) - 1)
  }
  
  NP <- nloptr::auglag(
    w0,
    fn,
    gr = NULL,
    lower = LB,
    upper = UB,
    hin = NULL,
    heq = heq,
    localsolver = "LBFGS",
    u = u,
    Ma = Ma,
    lmd = lmd,
    sgm = sgm
  )
  ## optimal weights
  optW <- NP$par
  ## optimal values
  optV <- NP$value
  ## convergence?
  conv <- NP$convergence
  
  if (conv < 0) {
    ## stop("NOT convergence! \n")
    warning("NOT convergence! \n")
  }  
  
  L <- list(optW = optW, optV = optV)
  
  return(L)
}





