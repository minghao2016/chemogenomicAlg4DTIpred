


calcLogLik <- function(
  cc = 5,
  inMat,
  thisAlpha = 0.8,
  U,
  V,
  thisBeta = 0.1,
  Sd,
  thisGamma = 0.1,
  St,
  lamU,
  lamV) {
  
  # INPUT
  # cc: default 5
  # inMat: input interaction matrix
  # thisAlpha: weight coefficient for U %*% t(V)
  # U: latent matrix for rows
  # V: latent matrix for cols
  # thisBeta: weight coefficient for Sd %*% U %*% t(V)
  # Sd: similarity for drugs
  # thisGamma: weighte coefficient for U %*% t(V) %*% St
  # St: similarity for target
  # lamU: lambda for U
  # lamV: lambda for V
  ## lamLD: lambda for laplacian of drug
  ## LD: laplacin matrix of drug
  ## lamLT: lambda for laplacian of target
  ## LT: laplacin matrix of target
  
  # OUTPUT
  # a scalar of log-likelihood
    
  
  Y <- inMat
  cY <- cc * Y
  
  Vt <- t(V)
  
  UVt <- U %*% Vt
  
  M <- thisAlpha * UVt + thisBeta * (Sd %*% UVt) + thisGamma * (UVt %*% St)


  log1pexpRes <- log1pexp(M)
  
  
  LL <- sum(cY * M) - sum((1 + cY - Y) * log1pexpRes) -
    0.5 * lamU * (base::norm(U, "F") ^ 2) - 0.5 * lamV * (base::norm(V, "F") ^ 2)
  
  return(LL)
}
