
setwd("your\\dir")

rm(list = ls())

## current data set name
db <- "nr"

switch (db,
        en = {
          cat("en data\n")
          flush.console()
          sd <- read.table("e_simmat_dc.txt")
          sd <- as.matrix(sd)
          st <- read.table("e_simmat_dg.txt")
          st <- as.matrix(st)
          Y <- read.table("e_admat_dgc.txt")
          Y <- as.matrix(Y) 
          Y <- t(Y)         
        },
        ic = {
          cat("ic data\n")
          flush.console()
          sd <- read.table("ic_simmat_dc.txt")
          sd <- as.matrix(sd)
          st <- read.table("ic_simmat_dg.txt")
          st <- as.matrix(st)
          Y <- read.table("ic_admat_dgc.txt")
          Y <- as.matrix(Y)
          Y <- t(Y)
        },
        gpcr = {
          cat("gpcr data\n")
          flush.console()
          sd <- read.table("gpcr_simmat_dc.txt")
          sd <- as.matrix(sd)
          st <- read.table("gpcr_simmat_dg.txt")
          st <- as.matrix(st)
          Y <- read.table("gpcr_admat_dgc.txt")
          Y <- as.matrix(Y)
          Y <- t(Y)
        },
        nr = {
          cat("nr data\n")
          flush.console()
          sd <- read.table("nr_simmat_dc.txt")
          sd <- as.matrix(sd)
          st <- read.table("nr_simmat_dg.txt")
          st <- as.matrix(st)
          Y <- read.table("nr_admat_dgc.txt")
          Y <- as.matrix(Y)
          Y <- t(Y)
        },
        stop("db should be one of the follows: 
             {en, ic, gpcr, nr}\n")
        )

## load required packages
pkgs <- c("matrixcalc", "data.table", "Rcpp", "ROCR", "Bolstad2", "MESS")
rPkgs <- lapply(pkgs, require, character.only = TRUE)

## source required R files
rSourceNames <- c(
  "doCVPositiveOnly.R",
  "doCVPositiveOnly3.R",
  "constrNeig.R",
  "inferZeros.R",
  "calcLogLik.R",
  "calcDeriv.R",
  "updateUV.R",
  "evalMetrics.R",
  "calcPredScore.R"
)
rSN <- lapply(rSourceNames, source, verbose = FALSE)

## sourceCPP required C++ files
cppSourceNames <- c("fastKF.cpp", "fastKgipMat.cpp", "log1pexp.cpp", "sigmoid.cpp")
cppSN <- lapply(cppSourceNames, sourceCpp, verbose = FALSE)


## convert to kernel
isKernel <- FALSE    ## default=TRUE
if (isKernel) {
  if (!isSymmetric(sd)) {
    sd <- (sd + t(sd)) / 2
  }
  epsilon <- 0.1
  while (!is.positive.semi.definite(sd)) {
    sd <- sd + epsilon * diag(nrow(sd))
  }
  if (!isSymmetric(st)) {
    st <- (st + t(st)) / 2
  }
  epsilon <- 0.1
  while (!is.positive.semi.definite(st)) {
    st <- st + epsilon * diag(nrow(st))
  }
}

Y <- t(Y)
tmp <- sd
sd <- st
st <- tmp

Y[1:3, 1:3]


## hyper-parameters

isDefaultPara <- TRUE
if (isDefaultPara) {
  if (db == "nr") {
    numLat = 10
  } else {
    numLat <- 50    
  }
  cc <- 7           
  thisAlpha <- 0.1  
  lamU <- 2         
  lamV <- 2        
  K1 <- 5           
} else {
  ## best for gpcr data
  numLat <- 90
  cc <- 6
  thisAlpha <- 0.4
  lamU <- 2
  lamV <- 2
  K1 <- 2
}

## values according to hyper-parameters
thisBeta <- (1 - thisAlpha)/2
thisGamma <- 1 - thisAlpha - thisBeta


nNeig <- 3
nIter <- 3

## sim diffusion or linear combine?
isLinearCom <- FALSE

## do cross-validation
kfold <- 10
numSplit <- 5

## DT-Hybrid method
savedFolds <- doCVPositiveOnly3(Y, kfold = kfold, numSplit = numSplit)

## saving results
resMetrics <- matrix(NA, nrow = kfold, ncol = 1)
colnames(resMetrics) <- c("MPR")
resMetrics <- as.data.frame(resMetrics)
finalResult <- vector("list", length = numSplit)

# main loop
for (i in 1:numSplit) {
  for (j in 1:kfold) {
    cat("numSplit:", i, "/", numSplit, ";", "kfold:", j, "/", kfold, "\n")
    flush.console()
    
    Yfold <- savedFolds[[i]][[j]][[6]]
    
    Yr <- inferZeros(Yfold, sd, K = K1)
    Yc <- inferZeros(t(Yfold), st, K = K1)
    
    KgipD <- fastKgipMat(Yr, 1)
    KgipT <- fastKgipMat(Yc, 1)
    
    if (isLinearCom) {
      ## linear combination
      cat("linear combine \n")
      flush.console()
      sd_temp <- 0.5 * (KgipD + sd)
      st_temp <- 0.5 * (KgipT + st)
    } else {
      cat("sim diffusion \n")
      flush.console()
      ## default: nNeig = 3, nIter = 2
      sd_temp <- fastKF(KgipD, sd, nNeig, nIter) 
      st_temp <- fastKF(KgipT, st, nNeig, nIter)
    }
    
    lap <- constrNeig(sd_temp, st_temp, K = K1)
    ########################################################
    
    lapD <- lap$lapD
    lapT <- lap$lapT
    simD <- lap$simD       
    simT <- lap$simT
    simDOrg <- lap$simDOrg 
    simTOrg <- lap$simTOrg
    
    ## use AdaGrid to update U and V
    UV <- updateUV(
      cc = cc,
      inMat = Yfold,
      thisAlpha = thisAlpha,
      thisBeta = thisBeta,
      Sd = simD,
      thisGamma = thisGamma,
      St = simT,
      lamU = lamU,
      lamV = lamV,
      numLat = numLat,
      initMethod = "useNorm",
      thisSeed = 123,
      maxIter = 100)
    
    U <- UV$U
    V <- UV$V
    
    testSet <- savedFolds[[i]][[j]][[1]]
    knownDrugIndex <- savedFolds[[i]][[j]][[4]]
    knownTargetIndex <- savedFolds[[i]][[j]][[5]] 
    testIndexRow <- savedFolds[[i]][[j]][[2]]      
    testIndexCol <- savedFolds[[i]][[j]][[3]]      
    ## result
    result <- calcPredScore(
      U = U,
      V = V,
      #simDrug = simDOrg,
      simDrug = simD,
      #simTarget = simTOrg,
      simTarget = simT,
      knownDrugIndex = knownDrugIndex,
      knownTargetIndex = knownTargetIndex,
      testIndexRow = testIndexRow,
      testIndexCol = testIndexCol,
      K = K1,
      thisAlpha = thisAlpha, 
      thisBeta = thisBeta,   
      thisGamma = thisGamma,
      testSet = testSet)
    resMetrics[j, ] <- result
  }
  finalResult[[i]] <- resMetrics
}

# combine result
resCom <- as.data.frame(data.table::rbindlist(finalResult))
resMean <- colMeans(resCom)
se <- sqrt(var(resCom[, 1]) / length(resCom[, 1]))
cat("DNILMF:", "MPR =", round(resMean, 3), "+\\-", round(se, 3),  "\n")
flush.console()
