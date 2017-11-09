

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
  "splitData.R"
)
rSN <- lapply(rSourceNames, source, verbose = FALSE)


Y <- t(Y)
tmp <- sd
sd <- st
st <- tmp

Y[1:3, 1:3]

####################
## sub-sample

Yorg <- Y
sdOrg <- sd
stOrg <- st


if (db = "en") {
  ## en
  subNum <- 150
  numRmSmp <- sample(subNum)
  Y <- Yorg[-numRmSmp,]
  sd <- sd[-numRmSmp,-numRmSmp]
  zeroCol <- which(colSums(Y) == 0)
  if (length(zeroCol) == 0) {
    st <- st
  } else {
    Y <- Y[, -zeroCol]
    st <- st[-zeroCol,-zeroCol]
  }
}

if (db == "ic") {
  subNum <- 130
  ## ic
  numRmSmp <- sample(subNum)
  Y <- Yorg[-numRmSmp, ]
  sd <- sdOrg[-numRmSmp, -numRmSmp]
  zeroCol <- which(colSums(Y) == 0)
  if (length(zeroCol) == 0) {
    st <- stOrg
  } else {
    Y <- Y[, -zeroCol]
    st <- stOrg[-zeroCol, -zeroCol]
  }
}

if (db == "gpcr") {
  subNum <- 70
  ## ic
  numRmSmp <- sample(subNum)
  Y <- Yorg[-numRmSmp, ]
  sd <- sdOrg[-numRmSmp, -numRmSmp]
  zeroCol <- which(colSums(Y) == 0)
  if (length(zeroCol) == 0) {
    st <- stOrg
  } else {
    Y <- Y[, -zeroCol]
    st <- stOrg[-zeroCol, -zeroCol]
  }
}



## do cross-validation
kfold <- 10
numSplit <- 5

## DT-Hybrid method
##savedFolds <- doCVPositiveOnly3(Y, kfold = kfold, numSplit = numSplit)
savedFolds <- splitData(Y, kfold, numSplit)


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
    
    ## look at splitData()
    Yfold <- savedFolds[[i]][[j]][["trainFold"]]

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
    
    ## look at splitData()
    testSet <- savedFolds[[i]][[j]][["testIndex"]]
    
    ## result
    result2 <- evalMetrics(Ypred = Ypred,testSet = testSet)
    resMetrics[j, ] <- result2
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

nr <- nrow(Y)
nc <- ncol(Y)
# save to file
curDate <- format(Sys.time(), format = "%Y%m%d")
curTime <- format(Sys.time(), format =  "%Hh%Mm%Ss")
savedFileName <- paste0(db, "_", "mpr", mprMeanSig, "+-", mprSeSig, "_", "numRmSmp", subNum, "_", "nr", nr, "nc", nc, "_", curDate, ".", curTime, ".RData")
cat("\n\n")
print(savedFileName)
save.image(file = savedFileName)
print(dim(Y))
