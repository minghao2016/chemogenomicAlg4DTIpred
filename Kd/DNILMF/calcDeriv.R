


calcDeriv <- function(
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
  lamV,
  isGradU = TRUE) {
  
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
  # isGradU: TRUE or FALSE
  
  # OUTPUT
  # gradient of U or V
  
  
  Y <- inMat
  cY <- cc * Y
  
  Vt <- t(V)
  
  UVt <- U %*% Vt
  
  A <- thisAlpha * UVt + thisBeta * (Sd %*% UVt) + thisGamma * (UVt %*% St)

  M <- 1 + cY - Y
  # sourceCpp("sigmoid.cpp")
  sigmoidRes <- sigmoid(A)
  P <- M * sigmoidRes
  
  
  
  if (isGradU) {
    YV <- Y %*% V
    PV <- P %*% V
    Dt <- t(Sd)
    TtV <- t(St) %*% V
    
    gradU <-
      thisAlpha * cc * YV + thisBeta * cc * (Dt %*% YV) + thisGamma * cc * (Y %*% TtV) -
      thisAlpha * PV - thisBeta * (Dt %*% PV) - thisGamma * (P %*% TtV) -
      lamU * U
    # return
    return(gradU)
  } else {
    Yt <- t(Y)
    YtU <- Yt %*% U
    DU <- Sd %*% U
    Pt <- t(P)
    PtU <- Pt %*% U
    gradV <-
      thisAlpha * cc * YtU + thisBeta * cc * (Yt %*% DU) + thisGamma * cc * (St %*% YtU) -
      thisAlpha * PtU - thisBeta * (Pt %*% DU) - thisGamma * (St %*% PtU) -
      lamV * V
    # return
    return(gradV)
  }
}




