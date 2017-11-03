

setwd("D:\\000\\chemogenomicAlg4DTIpred-master\\aucaupr\\dthybrid")

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
pkgs <- c("matrixcalc", "data.table", "Rcpp", "ROCR", "Bolstad2", "MESS")
rPkgs <- lapply(pkgs, require, character.only = TRUE)

## source required R files
rSourceNames <- c(
  "doCVPositiveOnly.R",
  "doCVPositiveOnly3.R",
  "evalMetrics.R",
  "Recommendation.R",
  "calAUPR.R"
)
rSN <- lapply(rSourceNames, source, verbose = FALSE)

## C++ functions
sourceCpp("sortIdxByCol2.cpp")

## final Y, sd, st
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

## saving results
resMetrics <- matrix(NA, nrow = kfold, ncol = 2)
colnames(resMetrics) <- c("auc", "aupr")
resMetrics <- as.data.frame(resMetrics)
finalResult <- vector("list", length = numSplit)

topK <- c(20, 30, 40, 50, 60, 70)
LenK <- length(topK)
resMetrics <- matrix(NA, nrow = LenK, ncol = 3)
colnames(resMetrics) <- c("topK", "auc", "aupr")
resMetrics[, "topK"] <- topK

# main loop
for (i in 1:numSplit) {
  for (j in 1:kfold) {
    cat("numSplit:", i, "/", numSplit, ";", "kfold:", j, "/", kfold, "\n")
    flush.console()
    
    Yfold <- savedFolds[[i]][[j]][[6]]

    ## row is drug, col is target
    if (db == "en") {
      theAl <- 0.4
    } else if (db == "ic") {
      theAl <- 0.3
    } else if (db == "gpcr") {
      theAl <- 0.2
    } else {
      theAl <- 0.4
    }
    Ypred <- computeRecommendation(A = Yfold, lambda = 0.5, 
                                   alpha = theAl, S = sd, S1 = st)
    
    testSet <- savedFolds[[i]][[j]][[1]]
    
    ## result
    ## result2 <- evalMetrics(Ypred = Ypred, testSet = testSet)
    
    
    ## columns in the test set
    colInTest <- unique(testSet[, "colIndex"])
    yLabel <- Y[, colInTest]
    yScore <- Ypred[, colInTest]
    
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
    
    
    
    
   
    
    resMetrics[j, ] <- result[1, ]
  }
  finalResult[[i]] <- resMetrics
}

# combine result
resCom <- as.data.frame(data.table::rbindlist(finalResult))

auc <- resCom[, "auc"]
aupr <- resCom[, "aupr"]

aucMean <- mean(auc)
aucSe <- sqrt(var(auc) / length(auc))

auprMean <- mean(aupr)
auprSe <- sqrt(var(aupr) / length(aupr))

aucMeanSig <- round(aucMean, 3)
aucSeSig <- round(aucSe, 3)

auprMeanSig <- round(auprMean, 3)
auprSeSig <- round(auprSe, 3)

cat("DTHybrid:", "AUC =", aucMeanSig, "+\\-", aucSeSig,  "\n",
    "    AUPR =", auprMeanSig, "+\\-", auprSeSig, "\n")
flush.console()

# save to file
curDate <- format(Sys.time(), format = "%Y%m%d")
curTime <- format(Sys.time(), format =  "%Hh%Mm%Ss")
savedFileName <- paste0(db, "_", "auc", aucMeanSig, "+-", aucSeSig, "_", 
                        "aupr", auprMeanSig, "+-", auprSeSig, "_", curDate, ".", curTime, ".RData")
cat("\n\n")
print(savedFileName)
save.image(file = savedFileName)

