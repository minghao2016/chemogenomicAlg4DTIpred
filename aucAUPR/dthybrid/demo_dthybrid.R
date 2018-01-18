
setwd("your dir")

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
        kd = {
          cat("kd data\n")
          Y <- read.table("drug-target_interaction_affinities_Kd__Davis_et_al.2011.txt")
          Y[Y <= 30] <- 1
          Y[Y > 30] <- 0
          Y <- as.matrix(Y)
          sd <- read.table("drug-drug_similarities_2D.txt")
          sd <- as.matrix(sd)
          st <- read.table("target-target_similarities_WS_normalized.txt")
          st <- as.matrix(st)
        },
        stop("db should be one of the follows: 
             {en, ic, gpcr, nr}\n")
        )


if (db == "kd") {
  dim(Y) ## 68 * 442
  dim(sd) ## 68 * 68
  dim(st) ##  442 * 442
  
  
  idxZeroCol <- which(colSums(Y) == 0)
  Y <- Y[, -idxZeroCol]
  st <- st[-idxZeroCol, -idxZeroCol]
  
  ## which(colSums(Y) == 0)
  
  idxZeroRow <- which(rowSums(Y) == 0)
  Y <- Y[-idxZeroRow, ]
  sd <- sd[-idxZeroRow, -idxZeroRow]
  
  which(rowSums(Y) == 0)
  which(colSums(Y) == 0)
  
  dim(Y)  ## 65 373
  dim(sd) ## 65 65
  dim(st) ## 373 373
  
  sd[1:3, 1:3]
  st[1:3, 1:3]
  Y[1:3, 1:3]
}

## load required packages
pkgs <- c("matrixcalc", "data.table", "Rcpp", "ROCR", "Bolstad2", "MESS")
rPkgs <- lapply(pkgs, require, character.only = TRUE)

## source required R files
rSourceNames <- c(
  "doCVPositiveOnly.R",
  "doCVPositiveOnly3.R",
  "doCrossVal.R",
  "calAUPR.R",
  "evalMetrics.R",
  "Recommendation.R"
)
rSN <- lapply(rSourceNames, source, verbose = FALSE)


Y <- t(Y)
tmp <- sd
sd <- st
st <- tmp

Y[1:3, 1:3]

## do cross-validation
kfold <- 10
numSplit <- 5

## DT-Hybrid method
savedFolds <- doCrossVal(Y, nfold = kfold, nsplit = numSplit)

## for saving results
AUPRVec <- vector(length = kfold)
AUCVec <- vector(length = kfold)
finalResult <- matrix(NA, nrow = numSplit, ncol = 2)
colnames(finalResult) <- c("AUPR", "AUC")


# main loop
for (i in 1:numSplit) {
  for (j in 1:kfold) {
    cat("numSplit:", i, "/", numSplit, ";", "kfold:", j, "/", kfold, "\n")
    flush.console()
    
    Yfold <- savedFolds[[i]][[j]][[1]]

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
    
    testSet <- savedFolds[[i]][[j]][[2]]
    
    testLabel <- Y[testSet]
    score <- Ypred[testSet]
    
    result <- calAUPR(testLabel, score)
    
    AUPRVec[j] <- result[1, "aupr"]
    AUCVec[j] <- result[1, "auc"]
  }
  AUPR <- mean(AUPRVec)
  AUC <- mean(AUCVec)
  finalResult[i, "AUPR"] <- AUPR
  finalResult[i, "AUC"] <- AUC
}


auc <- round(mean(finalResult[, "AUC"]), 3)
aucSD <- round(sd(finalResult[, "AUC"]), 3)

aupr <- round(mean(finalResult[, "AUPR"]), 3)
auprSD <- round(sd(finalResult[, "AUPR"]), 3)



# save to file
curDate <- format(Sys.time(), format = "%Y-%m-%d")
curTime <- format(Sys.time(), format =  "%H.%M.%S")
savedFileName <- paste0(db, "_", curDate, "_", curTime, "_auc", auc, "+-", aucSD, "_aupr", aupr, "+-", auprSD, ".RData")
cat("\n\n")
print(savedFileName)
save.image(file = savedFileName)