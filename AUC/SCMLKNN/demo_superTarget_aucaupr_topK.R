

setwd("D:\\000\\chemogenomicAlg4DTIpred-master\\aucaupr\\SCMLKNN")

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
          "Bolstad2", "MESS", "nloptr", "cluster")
rPkgs <- lapply(pkgs, require, character.only = TRUE)

## source required R files
rSourceNames <- c(
  "doCVPositiveOnly.R",
  "doCVPositiveOnly3.R",
  "evalMetrics.R",
  "getCvIndex.R",
  "superTargetDTI.R",
  "MLKNNtrain.R",
  "MLKNNtest.R",
  "predInteractionForS2.R",
  "calAUPR.R"
)
rSN <- lapply(rSourceNames, source, verbose = FALSE)


## sourceCPP required C++ files
cppSourceNames <- c("fastKF.cpp", "fastKgipMat.cpp", 
                    "log1pexp.cpp", "sigmoid.cpp")
cppSN <- lapply(cppSourceNames, sourceCpp, verbose = FALSE)



## convert to kernel
isKernel <- TRUE     
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

## do cross-validation
kfold <- 10
numSplit <- 5

## DT-Hybrid method
savedFolds <- doCVPositiveOnly3(Y, kfold = kfold, numSplit = numSplit)

## C++ functions
sourceCpp("sortIdxByCol2.cpp")

## saving results
if (db == "nr") {
  topK = c(10, 20)
} else {
  topK <- c(20, 30, 40, 50, 60, 70)
}
## saving results
LenK <- length(topK)
resMetrics <- matrix(NA, nrow = LenK, ncol = 3)
colnames(resMetrics) <- c("topK", "auc", "aupr")
resMetrics[, "topK"] <- topK
resMetrics <- as.data.frame(resMetrics)

foldResult <- vector("list", length = kfold)

finalResult <- vector("list", length = numSplit)

## main loop
for (i in 1:numSplit) {
  for (j in 1:kfold) {
    cat("numSplit:", i, "/", numSplit, ";", "kfold:", j, 
        "/", kfold, "\n")
    flush.console()
    
    ## training set with the test set links removed
    Yfold <- savedFolds[[i]][[j]][[6]]
    
    ## extract test set
    testSet <- savedFolds[[i]][[j]][[1]]
    knownDrugIndex <- savedFolds[[i]][[j]][[4]]
    knownTargetIndex <- savedFolds[[i]][[j]][[5]] 
    testIndexRow <- savedFolds[[i]][[j]][[2]]      
    testIndexCol <- savedFolds[[i]][[j]][[3]]      

    Ypred <- predInteractionForS2(Yfold = Yfold, simD = list(sd), simT = list(st)) 

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

