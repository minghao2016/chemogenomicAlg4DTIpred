
kronRlsC <- function(K1, K2, Yfold, lmd = NULL) {

  
  ## INPUT:
  # K1: kernel of drugs
  # K2: kernel of targets
  # Yfold: Y matrix of training set
  # lmd: lambda, for (K + lambda * I)^(-1) %*% y
  
  ## OUTPUT:
  # C: matrix, for coefficient C, that is matrix C
  
  if (is.null(lmd)) {
    lmd <- 1
  }
  

  r1 <- eigDecomp(K1)
  v1 <- r1$v
  vec1 <- r1$vec
  
  r2 <- eigDecomp(K2)
  v2 <- r2$v
  vec2 <- r2$vec
  
  L <- kronecker(t(v2), v1)
  

  L <- 1 / (L + lmd) 

  ## m1
  m1 <- t(vec1) %*% Yfold %*% vec2 
  

  ## m2
  m2 <- m1 * L
  #------------------------
  
  ## Y predicted values, for kronRLS
  C <- vec1 %*% m2 %*% t(vec2)
  
  
  return(C)
}



