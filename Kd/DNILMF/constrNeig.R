constrNeig <- function(simMatDrug, simMatTarget, K = 0) {
  # INPUT
  # simMatDrug: similarity matrix for drug
  # simMatTarget: similarity matrix for target
  # K: number of neighbors
  
  # OUTPUT
  # a list with laplacian matrix for drugs (lapD) and laplacian matrix for targets (lapT)
  # and simD and simT, which have already put diagonal elements to zeros
  
  
  numDrug <- nrow(simMatDrug)
  numTarget <- nrow(simMatTarget)
  
  simD <- simMatDrug - diag(diag(simMatDrug))
  simT <- simMatTarget - diag(diag(simMatTarget))
 
  if (K < 0) {
    stop("K MUST be greater or equal to zero! \n")
  }
  
  if (K > 0) {
    # for drug
    rankIndex <- t(apply(simD, 1, rank))
    neigSimD <- simD * (rankIndex > (numDrug - K))
    # calc laplacian matrix
    D1 <- rowSums(neigSimD)
    D2 <- colSums(neigSimD)
    lapD <- 0.5 * (diag(D1 + D2) - (neigSimD + t(neigSimD)))
    # for target
    rankIndex <- t(apply(simT, 1, rank))
    neigSimT <- simT * (rankIndex > (numTarget - K))
    # calc laplacian matrix
    D1 <- rowSums(neigSimT)
    D2 <- colSums(neigSimT)
    lapT <- 0.5 * (diag(D1 + D2) - (neigSimT + t(neigSimT)))
    # return
    lapList <- list()
    lapList$lapD <- lapD
    lapList$lapT <- lapT
    lapList$simD <- simD
    lapList$simT <- simT
    lapList$simDOrg <- simMatDrug
    lapList$simTOrg <- simMatTarget
    return(lapList)
  } else { # without neighbors
    # for drug
    D1 <- rowSums(simD)
    D2 <- colSums(simD)
    lapD <- 0.5 * (diag(D1 + D2) - (simD + t(simD)))    
    # for target
    D1 <- rowSums(simT)
    D2 <- colSums(simT)
    lapT <- 0.5 * (diag(D1 + D2) - (simT + t(simT)))    
    # return 
    lapList <- list()
    lapList$lapD <- lapD
    lapList$lapT <- lapT
    lapList$simD <- simD
    lapList$simT <- simT
    lapList$simDOrg <- simMatDrug
    lapList$simTOrg <- simMatTarget
    return(lapList)
  }
}
