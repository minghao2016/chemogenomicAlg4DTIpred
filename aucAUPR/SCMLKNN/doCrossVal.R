## training and test set split
doCrossVal <- function(adjMat = en, nfold = 10, nsplit = 5) {
  ##### GOAL
  ## ensure each row and each column have at least link after cv
  
  
  ##### INPUT
  ## adjMat: matrix, with binary values {0, 1}
  ## nfold: scalar, number of cross-validation
  ## nsplit: scalar, number of splits
  
  ##### OUTPUT
  
  
  ## all one positions
  onePos <- which(adjMat == 1, arr.ind = TRUE)
  onePosOrg <- onePos
  
  ## row sums
  rs <- rowSums(adjMat) ## row is target
  ## col sums
  cs <- colSums(adjMat)
  
  ##########################################################
  ##   import to select rational row and column
  ##########################################################
  ## From which row/column to select links to do cross-validation? 
  ## you can set to:
  ## 1: means row/column at least have two links 
  ## 2: means row/column at least have three links
  ## 3: means row/column at least has four links 
  ## more: means row/column at least has (more + 1) links
  nLink <- 1   ## default 2
  cat("nLink =", nLink, "\n")
  okRow <- which(rs > nLink)
  okCol <- which(cs > nLink)
  ##########################################################
  ##########################################################
  
  onePos <- onePos[onePos[, "row"] %in% okRow, ]
  onePos <- onePos[onePos[, "col"] %in% okCol, ]
  # nrow(onePos)
  
  
  # library(data.table)
  onePos <- as.data.table(onePos)
  
  ## good idea
  ## randomly remove [one row], when 'col' > 1
  dupCol <- onePos[, .N, by = col][N > 1][, col]
  if (length(dupCol) > 0) {
    tmp <- onePos[col %in% dupCol]
    ## randomly select one row: sample()
    tmp <- tmp[, sample(row, 1), by = col] ## [by = col]
    tmp <- cbind(row = tmp[, V1], tmp)
    tmp[, V1 := NULL]
    ## way to match unwanted rows: paste0()
    isKept <- !(onePos[, paste0(row, ".", col)] %in% tmp[, paste0(row, ".", col)])
    onePos <- onePos[isKept, ]
  }
  
  ## randomly remove [one row], when 'row' > 1
  dupRow <- onePos[, .N, by = row][N > 1][, row]  
  if (length(dupRow) > 0) {
    tmp <- onePos[row %in% dupRow]
    tmp <- tmp[, sample(col, 1), by = row]
    setnames(tmp, "V1", "col")
    isKept <- !(onePos[, paste0(row, ".", col)] %in% tmp[, paste0(row, ".", col)])
    onePos <- onePos[isKept, ]
  }
  ###########################################################
  ## here, [onePos] can 'safely' do cross-validation
  
  
  ## deal with ZEROs position
  zeroPos <- which(adjMat == 0, arr.ind = TRUE)
  zeroPosOrg <- zeroPos
  zeroPos <- as.data.table(zeroPos)
  
  
  ########################################################
  #####################   nested list   ##################
  ## save final list
  ## nested list, should define two lists
  savedFolds <- vector(mode = "list", length = nsplit)
  names(savedFolds) <- paste0("split_", 1:nsplit)
  
  # save nfold list for each split
  cvFolds <- vector("list", length = nfold)
  names(cvFolds) <- paste0("fold_", 1:nfold)
  ########################################################
  ########################################################
  
  ## num is based on 'onePos'
  ## kevin's way to generate CV folds
  numOne <- nrow(onePos)
  if (numOne < nfold) {
    stop("available number of ONEs is less than [nfold]!\n")
  }
  
  numZero <- nrow(zeroPos)
  if (numZero < nfold) {
    stop("available number of ZEROs is less than [nfold]!")
  }
  
  ######################################################
  
  for (ii in 1:nsplit) {
    ## must use: ceiling() to get enough entries in matrix
    foldSize <- ceiling(numOne / nfold)
    ## good idea: each row is one-fold
    foldsOne <- matrix(NA, nrow = nfold, ncol = foldSize)
    #############################
    ##   random sample  #########
    foldsOne[1:numOne] <- sample(numOne)
    ##   random sample  #########
    #############################
    ## List
    ## 'matrix to List', fast
    ## row(folds), where 'folds' should be 'matrix'
    foldsOne <- split(foldsOne, row(foldsOne))
    ## Map(na.omit, folds) OR lapply(folds, na.omit), but with attribute
    foldsOne <- lapply(foldsOne, function(x) x[!is.na(x)])
    
    ## names(folds) <- paste0("fold", 1:nfold)
    
    ## for ZEROs
    foldSize <- ceiling(numZero / nfold)
    foldsZero <- matrix(NA, nrow = nfold, ncol = foldSize)
    foldsZero[1:numZero] <- sample(numZero)
    foldsZero <- split(foldsZero, row(foldsZero))
    foldsZero <- lapply(foldsZero, function(x) x[!is.na(x)])
    
    for (jj in 1:nfold) {
      ## for ONEs
      curFoldOne <- foldsOne[[jj]]
      testIndexOne <- onePos[curFoldOne, ]
      
      ## trainIndexOne <- onePos[setdiff(1:numOne, curFoldOne), ]
      
      ## for ZEROs
      curFoldZero <- foldsZero[[jj]]
      testIndexZero <- zeroPos[curFoldZero, ]
      
      ## trainIndexZero <- zeroPos[setdiff(1:numZero, curFoldZero), ]
      
      ## pick up these rows as 'test set'
      ## matrix
      testIndex <- rbind(testIndexOne, testIndexZero)
      testIndex <- as.matrix(testIndex)
      ###################################
      ## same: doCVPosiiveOnly3.R
      ## testIndexRow <- testIndex[, row]
      ## testIndexCol <- testIndex[, col]
      ###################################
      
      
      ## split by row
      ## testIndexByRow <- split(testIndex[, col], f = testIndex[, row])
      # length(testIndex) 
      # testIndexByRow[1:3]
      
      ## put test test as ZEROs, then extract training index by row
      ## matrix-matrix operation
      adjMatTemp <- adjMat
      adjMatTemp[testIndex] <- 0
      
      ## which(x, arr.ind = TRUE)
      ## matrix
      ## trainIndex <- which(adjMatTemp == 1, arr.ind = TRUE) ## used in cascadeSim
      ## trainIndex <- rbind(trainIndexOne, trainIndexZero)
      ## List
      ## trainIndexByRow <- split(trainIndex[, col], f = trainIndex[, row])
      ## trainIndexByRow <- lapply(trainIndexByRow, function(x) unname(x))
      # length(trainIndexByRow)
      # trainIndexByRow[1:3]
      # which(adjMat[1, ] == 1); which(adjMat[2, ] == 1)
      
      ## prepare for saving
      ## rownames(trainIndex) <- NULL ## used in cascadeSim
      ## trainIndex <- as.matrix(trainIndex)
      
      ################################
      ## same: doCVPositiveOnly3.R
      # vector
      ## knownDrugIndex <- sort(unique(trainIndex[, "row"]))
      # vector
      ## knownTargetIndex <- sort(unique(trainIndex[, "col"]))
      ################################
      
      ## results into List
      res <- list(
        ## training fold
        trainFold = adjMatTemp,
        ## 2-column matrix, i = row, j = column
        ## trainIndex = trainIndex, ## 2-column matrix
        ## List
        ## trainIndexByRow = trainIndexByRow,
        ## 2-column matrix
        testIndex = testIndex ## 2-column matrix
        ## List
        ## testIndexByRow = testIndexByRow,
        ## below same: doCVPositiveOnly3.R
        ## testIndexRow = testIndexRow,
        ## testIndexCol = testIndexCol,
        ## knownDrugIndex = knownDrugIndex,
        ## knownTargetIndex = knownTargetIndex
      )
      # length(res)
      # str(res)
      cvFolds[[jj]] <- res
    }
    savedFolds[[ii]] <- cvFolds
  }
  ## save to file  
  # save(savedFolds, file = "savedFolds.RData")
  return(savedFolds)
}





