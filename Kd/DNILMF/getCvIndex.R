getCvIndex <- function(totNum, nfold) {
  # OBJECTIVE
  # get indices for nfold cross-validation
  
  # INPUT
  # totNum: total number of samples
  # nfold: number of folds you need
  
  # OUTPUT
  # list: with nfold element
  
  # Cross-validation folds
  lenSeg <- ceiling(totNum / nfold)
  incomplete <- nfold * lenSeg - totNum     
  complete <- nfold - incomplete                  
  inds <- matrix(c(sample(1:totNum), rep(NA, incomplete)), nrow = lenSeg, byrow = TRUE)
  folds <- lapply(as.data.frame(inds), function(x) c(na.omit(x)))
  
  return(folds)
}

