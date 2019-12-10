source("reconciliationMethods.R")
library(foreach)
library(doMC)
library(bigmemory)
options(bigmemory.allow.dimnames=TRUE)
registerDoMC(detectCores()-1)
library(tictoc)

dset="infantgts"
h=1
fmethod="ets"
iTest=1
seed=0
testSize=50

#hierRecBayesianExperiment(dset, h, fmethod, iTest)

hierRecBayesianExperiment <- function(dset, h, fmethod="ets", iTest=1, testSize=50,
                        seed=0, synth_n=100, synthCorrel=0.5, # Parameters only for synthetic TS
                        printResults=FALSE, testProbability=0.1)
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
  mseBase =  mean  ( (allts(test)[h,] - preds)^2 )
  predsBottomUp = mSumMatrix %*% preds[bottomIdx]
  mseBottomUp = mean  ( (allts(test)[h,] - predsBottomUp)^2 )
  
  # # Calculating initial MSEs
  mCov_shrmint = estimateCovariance(residuals, method='shrmint')
  #mCov_glasso = estimateCovariance(residuals, method='glasso')
  
  outBayesShr = bayesReconFull(preds, mSumMatrix, mCov_shrmint)
  fcastBayesShr = outBayesShr$coherentPreds
  
  outBayesShrNotFull = bayesReconNotFull(preds, mSumMatrix, mCov_shrmint)
  fcastBayesShrNotFull = outBayesShrNotFull$coherentPreds
  
  mseBayesShrFull =  mean((allts(test)[h,] - fcastBayesShr)^2)
  mseBayesShrNotFull =  mean((allts(test)[h,] - fcastBayesShrNotFull)^2)
  
  # Energy Score
  target = allts(test)[h,]
  covarianceFull = outBayesShr$coherentCovariance
  meanFull = as.numeric(fcastBayesShr)
  
  covarianceNotFull = outBayesShrNotFull$coherentCovariance
  meanNotFull = as.numeric(fcastBayesShrNotFull)
  
  esBayesFull = energyScore(target, meanFull, covarianceFull)$score
  esBayesNotFull = energyScore(target, meanNotFull, covarianceNotFull)$score
  
  #fcastBayesGlasso = bayesReconFull(preds, mSumMatrix, mCov_glasso)
  #fcastBayesGlassoNotFull = bayesReconNotFull(preds, mSumMatrix, mCov_glasso)
  mseBayesGlassoFull =  NaN#mean((allts(test)[h,] - fcastBayesGlasso)^2)
  mseBayesGlassoNotFull =  NaN#mean((allts(test)[h,] - fcastBayesGlassoNotFull)^2)
  
  # Checking MinT library implementation with a certain probability
  mseMintShr = NaN
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
    mseMintShr  <- mean((allts(test)[h,] - fcastMintShr)^2)
    diff = fcastMintShr - fcastBayesShr
    if (norm(diff, "2") > 1.0){
      print(mseBayesShrFull)
      print(mseMintShr)
      stop("MinT and Bayes result diverged.")
    }
    toc()
  }
  
  if (printResults){
    print(c("MseBase: ",mseBase))
    print(c("MseBottomUp: ",mseBottomUp))
    print(c("MseBayesShr: ",mseBayesShrFull))
    print(c("MseMinTShr: ",mseMintShr))
    print(c("MseBayesShrNotFull: ",mseBayesShrNotFull))
    print(c("MseBayesGlasso: ",mseBayesGlassoFull))
    print(c("MseBayesGlassoNotFull: ",mseBayesGlassoNotFull))
  }
  
  if (dset=="synthetic"){
    dset <- paste0(dset,"_correl",synthCorrel,"_n",synth_n,"_seed",seed)
  }
  else if (dset=="syntheticLarge"){
    dset <- paste0("largeSynthetic_n",synth_n,"_seed",seed,"_correl", synthCorrel)
  }
  
  cols = c("dset", "h", "iTest", "fmethod",
           "mseBase", "mseBottomUp",
           "mseBayesShr", "mseBayesShrFull",
           "mseMinT",
           "mseBayesGlasso", "mseBayesGlassoFull",
           "esBayesShrNotFull", "esBayesShrFull")
  
  dataFrame <- data.frame(dset, h, iTest, fmethod, mseBase, mseBottomUp, mseBayesShrNotFull, mseBayesShrFull,
                          mseMintShr, mseBayesGlassoNotFull, mseBayesGlassoFull,
                          esBayesNotFull, esBayesFull)
  colnames(dataFrame) <- cols
  return(dataFrame)
}