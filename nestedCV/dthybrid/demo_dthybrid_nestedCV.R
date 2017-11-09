

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
pkgs <- c("matrixcalc", "data.table", "Rcpp", "ROCR", "Bolstad2", "MESS")
rPkgs <- lapply(pkgs, require, character.only = TRUE)

## source required R files
rSourceNames <- c(
  "doCVPositiveOnly.R",
  "doCVPositiveOnly3.R",
  "evalMetrics.R",
  "Recommendation.R",
  "getOptParaDTHybrid.R"
)
rSN <- lapply(rSourceNames, source, verbose = FALSE)


Y <- t(Y)
tmp <- sd
sd <- st
st <- tmp

Y[1:3, 1:3]

## do cross-validation
kfold <- 10     ## default = 10
numSplit <- 1   ## default = 5

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
    cat("outer loop...", "numSplit:", i, "/", numSplit, ";", "kfold:", j, "/", 
        kfold, "\n")
    flush.console()
    
    Yfold <- savedFolds[[i]][[j]][[6]]
    
    ## perform nested CV for best parameter selection
    bestPara <- getOptParaDTHybrid(
      Yfold, innerFold = 5, numSplitInner = 1,
      paraList = list(lambda = seq(0.2, 0.9, by = 0.1),
                      alpha = seq(0.2, 0.9, by = 0.1)),
      sd,
      st)
    ## best parameters
    bestLambda <- bestPara["lambda"]
    bestAlpha <- bestPara["alpha"]
    
    testSet <- savedFolds[[i]][[j]][[1]]
    
    Ypred <- computeRecommendation(A = Yfold, lambda = bestLambda, 
                                   alpha = bestAlpha, S = sd, S1 = st)
    ## result
    result <- evalMetrics(Ypred = Ypred, testSet = testSet)
    resMetrics[j, ] <- result[, "MPR"]
  }
  finalResult[[i]] <- resMetrics
}

# combine result
resCom <- as.data.frame(data.table::rbindlist(finalResult))

mpr <- resCom[, "MPR"]

mprMean <- mean(mpr) 

mprSe <- sqrt(var(mpr) / length(mpr))

mprMeanSig <- round(mprMean, 3)
mprSeSig <- round(mprSe, 3)

cat("DTHybrid:", "MPR =", mprMeanSig, "+\\-", mprSeSig,  "\n")
flush.console()

# save to file
curDate <- format(Sys.time(), format = "%Y%m%d")
curTime <- format(Sys.time(), format =  "%Hh%Mm%Ss")
savedFileName <- paste0(db, "_", "mpr", mprMeanSig, "+-", mprSeSig, "_nestedCV_", curDate, ".", curTime, ".RData")
cat("\n\n")
print(savedFileName)
save.image(file = savedFileName)



