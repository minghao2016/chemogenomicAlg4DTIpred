doCrossValidation <- function(inMat, kfold = 10, numSplit = 5) {
  # INPUT
  # inMat: input matrix
  # kfold: k-fold cross-validation
  # numSplit: number of splits
  
  # OUTPUT
  # a list of multiple components
  
  
  # library(reshape2)
  # library(data.table)
  
  

  rownames(inMat) <- NULL
  colnames(inMat) <- NULL
  
  
  # width -> length, stack by column using "melt()"
  # i, j, v format
  triplet <- as.data.table(reshape2::melt(inMat))
  setnames(triplet, c("rowIndex", "colIndex", "value"))
  
  numTri <- nrow(triplet)
  
  #####################   nested list   ############
  # save final list
  # nested list, should define two lists
  savedFolds <- vector(mode = "list", length = numSplit)
  names(savedFolds) <- paste0("split_", 1:numSplit)
  
  # save kfold list for each split
  cvFolds <- vector("list", length = kfold)
  names(cvFolds) <- paste0("fold_", 1:kfold)
  ##################################################
  
  
  source("getCvIndex.R")
  
  for (i in 1:numSplit) {
    #######
    folds <- getCvIndex(totNum = numTri, nfold = kfold)
    #######
    for (j in 1:kfold) {
      
      currIndex <- folds[[j]]
      
      testData <- triplet[currIndex]
      
      testIndex <- testData[, 1:2, with = FALSE]
      testIndex <- as.matrix(testIndex)
      
      ###############################
      testLabel <- inMat[testIndex]
      # test labels: used for calculating AUPR and AUC
      testLabel2 <- testData[, value]
      ###############################
      
      
      testIndexRow <- testIndex[, "rowIndex"]
      testIndexCol <- testIndex[, "colIndex"]
      
      # known information for drug-target matrix
      tmpTriplet <- triplet
      tmpTriplet[currIndex, "value"] <- 0
      knownInteraction <- tmpTriplet[value > 0]
      knownDrugIndex <- unique(knownInteraction[, rowIndex])
      knownTargetIndex <- unique(knownInteraction[, colIndex])
       
      # fold matrix in the test set
      tmp <- inMat
      tmp[testIndex] <- 0
      
      # source("getInteractType.R")
      interactType <- getInteractType(tmp)
      
      # 1
      cvFolds[[j]]$testLabel <- testLabel
      # 2
      cvFolds[[j]]$testIndex <- testIndex
      # 3
      cvFolds[[j]]$testIndexRow <- testIndexRow
      # 4
      cvFolds[[j]]$testIndexCol <- testIndexCol
      # 5
      cvFolds[[j]]$knownDrugIndex <- knownDrugIndex
      # 6
      cvFolds[[j]]$knownTargetIndex <- knownTargetIndex
      # 7
      cvFolds[[j]]$foldMat <- tmp
      # 8
      cvFolds[[j]]$interactType <- interactType
      # 9
      cvFolds[[j]]$currIndexTest <- currIndex
      # 10
      cvFolds[[j]]$testLabel2 <- testLabel2
    }
    savedFolds[[i]] <- cvFolds
  }
  
  cat("save 'savedFolds.RData' to disk! \n")
  flush.console()
  
  save(savedFolds, file = "savedFolds.RData")
  
  return(savedFolds)
}

