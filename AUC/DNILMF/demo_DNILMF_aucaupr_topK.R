
setwd("your dir")

rm(list = ls())

## current data set name
db <- "en"

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
pkgs <- c("matrixcalc", "data.table", "Rcpp", "ROCR", 
  "Bolstad2", "MESS", "reshape2")
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
  "calcPredScore.R",
  "calAUPR.R",
  "calcPredScoreAUCAUPR.R"
)
rSN <- lapply(rSourceNames, source, verbose = FALSE)

## sourceCPP required C++ files
cppSourceNames <- c("fastKF.cpp", "fastKgipMat.cpp", "log1pexp.cpp", "sigmoid.cpp")
cppSN <- lapply(cppSourceNames, sourceCpp, verbose = FALSE)

## C++ functions
sourceCpp("sortIdxByCol2.cpp")



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
if (db == "nr") {
  topK = c(10, 20)
} else {
  topK <- c(20, 30, 40, 50, 60, 70)
}

LenK <- length(topK)
resMetrics <- matrix(NA, nrow = LenK, ncol = 3)
colnames(resMetrics) <- c("topK", "auc", "aupr")
resMetrics[, "topK"] <- topK
resMetrics <- as.data.frame(resMetrics)

foldResult <- vector("list", length = kfold)

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
      sd <- 0.5 * (KgipD + sd)
      st <- 0.5 * (KgipT + st)
    } else {
      cat("sim diffusion \n")
      flush.console()
      ## default: nNeig = 3, nIter = 2
      sd <- fastKF(KgipD, sd, nNeig, nIter) 
      st <- fastKF(KgipT, st, nNeig, nIter)
    }
    
    lap <- constrNeig(sd, st, K = K1)
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
    Ypred <- calcPredScoreAUCAUPR(
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
      testSet = testSet,
      Y = Y)
    
    ## columns in the test set
    colInTest <- unique(testSet[, "colIndex"])
    yLabel <- Y[, colInTest, drop = FALSE]
    yScore <- Ypred[, colInTest, drop = FALSE]
    
    for (ii in 1:LenK) {
      curK <- topK[ii]
      res <- sortIdxByCol2(yScore, yLabel, curK)
      yLabelSorted <- res$labelSorted
      yScoreSorted <- res$scoreSorted
      ## result
      yLabelSorted <- as.vector(yLabelSorted)
      yScoreSorted <- as.vector(yScoreSorted)
      ## auc and aupr
      result <- calAUPR(yLabelSorted, yScoreSorted)
      resMetrics[ii, c("auc", "aupr")] <- result[1, ]
    }
    foldResult[[j]] <- resMetrics
  }
  finalResult[[i]] <- foldResult
}

resDT <- NULL
for (i in 1:numSplit) {
  resDT <- rbind(resDT, data.table::rbindlist(finalResult[[i]]))
}

print(resDT[, mean(auc), by = topK])
print(resDT[, sqrt(var(auc) / length(auc)), by = topK])

print(resDT[, mean(aupr), by = topK])
print(resDT[, sqrt(var(aupr) / length(aupr)), by = topK])


# save to file
curDate <- format(Sys.time(), format = "%Y%m%d")
curTime <- format(Sys.time(), format =  "%Hh%Mm%Ss")
savedFileName <- paste0(db, "_topK_", curDate, ".", curTime, ".RData")
cat("\n\n")
print(savedFileName)
save.image(file = savedFileName)

