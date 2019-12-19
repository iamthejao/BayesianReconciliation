source("reconciliationMethods.R")
library(foreach)
library(doMC)
library(bigmemory)
options(bigmemory.allow.dimnames=TRUE)
registerDoMC(detectCores()-1)
library(tictoc)

# dset="infantgts"
# h=1
# fmethod="ets"
# iTest=49
# seed=0
# testSize=50

#hierRecBayesianExperiment(dset, h, fmethod, iTest)

hierRecBayesianExperiment <- function(dset, h, fmethod="ets", iTest=1, testSize=50,
                        seed=0, synth_n=100, synthCorrel=0.5, # Parameters only for synthetic TS
                        printResults=FALSE, testProbability=0.05)
  {
  
  # Setting seed and loading dataset
  set.seed(seed)
  hierTs = loadDataset(dset)
  allTs = allts(hierTs)
  columnNames = colnames(allTs)
  numTs <- ncol(allTs)
  alpha <- 0.2
  
  #if h=1, the possible preds are the whole test size lenght; 
  #if h=2, the possible preds are the (test size lenght -1); etc.
  possiblePreds <- testSize - h + 1
  if (iTest>possiblePreds){
    stop("iTest>possiblePreds")
  }
  
  # Slicing train and test
  timeIdx             <- time(hierTs$bts[,1])
  startTrain          <- timeIdx[1]
  endTrain            <- length(timeIdx) - h - (iTest - 1)
  train               <- window(hierTs, start = startTrain, end = timeIdx[endTrain] )
  test                <- window(hierTs, start = timeIdx[endTrain +1], end=timeIdx[endTrain + h])
  allTsTrain <- allts(train)
  
  # Creating matrices to save results
  # Making vector and matrixes for storaging values in shared memory for parallel calculation
  # Tourism or ARIMA, otherwise it is quicker to do single core due to overhead
  parallel=FALSE
  if ((dset == 'tourism') || (fmethod == 'arima')){
    parallel=TRUE
  }
  
  # Predictions from model
  preds <- big.matrix(ncol=numTs, nrow=1, type='double', shared=parallel)
  
  # Fit from model
  fitted <- big.matrix(nrow=dim(allTsTrain)[1], ncol = numTs, type='double', shared=parallel)
  colnames(fitted) <- colnames(allTsTrain)
  
  # Iteration on parallel or single core
  if (parallel){
    '%iter%' <- foreach::'%dopar%'
  } else {
    '%iter%' <- foreach::'%do%'
  }
    
  print("Starting series fits")
  tic("fit models")
  foreach (i = 1:numTs) %iter% {
    if (parallel==FALSE){
      i #does nothing
      #print(paste("Series num: ", i))
    }
    if (fmethod=="ets"){
      model <- ets(allTsTrain[,i], additive.only = TRUE)
    }
    else if (fmethod=="arima"){
      model <- auto.arima(allTsTrain[,i])
    }
    tmp <- forecast(model, h=h, PI=FALSE)
    fitted[,i] <- model$fitted 
    preds[i] <- tmp$mean[h]
  }
  toc()
  
  # Sum matrix and indices
  mSumMatrix = smatrix(train)
  mSumMatrixT = t(mSumMatrix)
  bottomIdx <- seq( nrow(mSumMatrix) - ncol(mSumMatrix) +1, nrow(mSumMatrix))
  upperIdx <- setdiff(1:nrow(mSumMatrix),bottomIdx)
  
  # Recovering matrices from shared mem
  fitted = as.matrix(fitted)
  preds = as.matrix(preds)
  y = aggts(train)
  residuals = (y - fitted)
  residuals[, upperIdx] = residuals[, upperIdx] * -1
  
  # Base and BottomUp MSE
  rmseBase =  sqrt(mean  ( (allts(test)[h,] - preds)^2 ))
  predsBottomUp = mSumMatrix %*% preds[bottomIdx]
  rmseBottomUp = sqrt(mean  ( (allts(test)[h,] - predsBottomUp)^2 ))
  
  # Calculating covariance matrix using shrink
  mCov_shrmint = estimateCovariance(residuals, method='shrmint')
  
  # Calculating outputs
  # Full covariance bayes with and wo positivity constraint
  # Diagonal covariance bayes with and wo positivity constraint
  outBayesFull = bayesReconFull(preds, mSumMatrix, mCov_shrmint, positivity=FALSE)
  outBayesFullPositive = bayesReconFull(preds, mSumMatrix, mCov_shrmint, positivity=TRUE)
  outBayes = bayesReconNotFull(preds, mSumMatrix, mCov_shrmint, positivity=FALSE)
  outBayesPositive = bayesReconNotFull(preds, mSumMatrix, mCov_shrmint, positivity=TRUE)
  
  # Taking full samples
  target = allts(test)[h,]
  sampleFull = outBayesFull$sample %*% t(mSumMatrix)
  sampleFullPos = outBayesFullPositive$truncSample %*% t(mSumMatrix)
  sample = outBayes$sample %*% t(mSumMatrix)
  samplePos = outBayesPositive$truncSample %*% t(mSumMatrix)
  
  # Calculating Statistics
  statsBayesFull = calculateStatistics(target,
                                       outBayesFull$coherentPreds, outBayesFull$coherentPredsMedianSample,
                                       sampleFull, upperIdx, bottomIdx,
                                       suffix="_Full")
  
  statsBayesFullPos = calculateStatistics(target,
                                          outBayesFullPositive$coherentPredsTrunc, outBayesFullPositive$coherentPredsMedianTrunc,
                                          sampleFullPos, upperIdx, bottomIdx,
                                          suffix="_FullPos")
  
  statsBayes = calculateStatistics(target,
                                   outBayes$coherentPreds, outBayes$coherentPredsMedianSample,
                                   sample, upperIdx, bottomIdx,
                                   suffix="")
  
  statsBayesPos = calculateStatistics(target,
                                   outBayesPositive$coherentPredsTrunc, outBayesPositive$coherentPredsMedianTrunc,
                                   samplePos, upperIdx, bottomIdx,
                                   suffix="_Pos")
  
  # Checking MinT library implementation with a certain probability
  rmseMint = NaN
  if (runif(1) < testProbability){
    print("Performing MinT on library for double check")
    tic("Running MinT")
    if (fmethod == "ets"){
      fcastMintShr <-forecast.gts(train, h = h, method = "comb", weights="mint",
                                  fmethod=fmethod, covariance="shr", additive.only = TRUE,
                                  parallel=parallel, num.cores=10)
    } else {
      fcastMintShr <-forecast.gts(train, h = h, method = "comb", weights="mint",
                                  fmethod=fmethod, covariance="shr",
                                  parallel=parallel, num.cores=10)
    }
    fcastMintShr = allts(fcastMintShr)[h,]
    rmseMint  <- sqrt(mean((allts(test)[h,] - fcastMintShr)^2))
    diff = fcastMintShr - outBayesFull$coherentPreds
    if (norm(diff, "2") > 1.0){
      print(statsBayesFull$rmseMean_Full)
      print(rmseMint)
      stop("MinT and Bayes result diverged.")
    }
    toc()
  }
  
  if (dset=="synthetic"){
    dset <- paste0(dset,"_correl",synthCorrel,"_n",synth_n,"_seed",seed)
  }
  else if (dset=="syntheticLarge"){
    dset <- paste0("largeSynthetic_n",synth_n,"_seed",seed,"_correl", synthCorrel)
  }
  
  statistics = c(statsBayesFull, statsBayesFullPos, statsBayes, statsBayesPos)
  output = list(dset=dset, h=h, iTest=iTest, fmethod=fmethod, rmseBase=rmseBase, rmseBottomUp=rmseBottomUp, rmseMint=rmseMint)
  dataFrame = data.frame(c(output, statistics))
  return(dataFrame)
}