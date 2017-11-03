

eigDecomp <- function(K) {
  ## INPUT:
  # K: matrix, PSD kernel
  
  res <- base::eigen(K)
  eigValue <- res$values
  eigVec <- res$vectors
  
  res <- list(v = eigValue, vec = eigVec)
  
  return(res)
}




