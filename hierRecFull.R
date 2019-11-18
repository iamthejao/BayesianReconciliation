source("reconciliationMethods.R")
library(foreach)
library(doMC)
library(bigmemory)
options(bigmemory.allow.dimnames=TRUE)
registerDoMC(detectCores()-1)
library(tictoc)

# dset="tourism"
# h=1
# fmethod="ets"
# iTest=1
# seed=0

hierRecFull <- function (dset, h=1, fmethod="ets", iTest=1, 
                     seed=0, synth_n=100, synthCorrel=0.5, printResults=FALSE, saveFiles=TRUE){
  
  testSize <- 50
  set.seed(seed)
  hierTs = loadDataset(dset)
  
  #if h=1, the possible preds are the whole test size lenght; 
  #if h=2, the possible preds are the (test size lenght -1); etc.
  possiblePreds <- testSize - h + 1
  if (iTest>possiblePreds){
    stop("iTest>possiblePreds")
  }
  #here the experiment starts
  timeIdx <- time(hierTs$bts[,1])
  startTrain          <- timeIdx[1]
  endTrain            <- length(timeIdx) - h - (iTest - 1)
  train               <- window(hierTs, start = startTrain, end = timeIdx[endTrain] ) #timeIdx[startTrain]
  test                <- window(hierTs, start = timeIdx[endTrain +1], end=timeIdx[endTrain + h])
  
  # Calculating using Mint for comparison
  #  sometimes the sample matrix is not positive definite and minT crashes
  # the matrix is computed internally by.libPaths() minT and cannot be controlled from here.
  # How can I make this faster?
  
  print("Starting MinT")
  tic("MinT")
  startMint = Sys.time()
  mseCombMintSample <- NA
  fcastCombMintShr <-
    forecast(train, h = h, method = "comb", weights="mint", fmethod=fmethod, 
             covariance="shr")
  mseCombMintShr  <- hierMse(fcastCombMintShr, test,  h)
  toc()
  
  # Organizing variables for BayesianReconciliation
  # Recompute predictions to be  accessed by the Bayesian method
  allTsTrain <- allts(train)
  numTs <- ncol(allTsTrain)
  alpha <- 0.2
  
  
  # Making vector and matrixes for storaging values in shared memory for parallel calculation
  sigma <- big.matrix(ncol=numTs, nrow=1, type='double') #vector(length = numTs)
  preds <- big.matrix(ncol=numTs, nrow=1, type='double') #vector(length = numTs)
  
  # Residuals from model
  residuals <- big.matrix(nrow=dim(allTsTrain)[1], ncol = numTs, type='double')#matrix(nrow=dim(allTsTrain)[1], ncol = numTs)
  colnames(residuals) <- colnames(allTsTrain)
  # Fit from model
  fitted <- big.matrix(nrow=dim(allTsTrain)[1], ncol = numTs, type='double')#matrix(nrow=dim(allTsTrain)[1], ncol = numTs)
  colnames(fitted) <- colnames(allTsTrain)
  # Same as allTsTrain
  actual <- big.matrix(nrow=dim(allTsTrain)[1], ncol = numTs, type='double')#matrix(nrow=dim(allTsTrain)[1], ncol = numTs)
  colnames(actual) <- colnames(allTsTrain)
  
  # Fit models on training data
  print("Starting Bayes")
  tic("Bayes")
  startBayes = Sys.time()
  foreach (i = 1:numTs) %dopar% {
    print(paste("Current TS:", i))
    if (fmethod=="ets"){
      model <- ets(allTsTrain[,i], additive.only = TRUE)
    }
    else if (fmethod=="arima"){
      model <- auto.arima(allTsTrain[,i])
    }
    tmp <- forecast(model, h=h, level=1-alpha)
    residuals[,i] <- model$residuals #we could store model$residuals if we allowed  multiplicative errors
    fitted[,i] <- model$fitted 
    actual[,i] <- model$x
    preds[i] <- tmp$mean[h]
    sigma[i] <- abs ( (tmp$mean[1] - tmp$upper[1])  / (qnorm(alpha / 2)) )
  }
  
  # Recovering
  sigma = as.matrix(sigma)
  preds = as.matrix(preds)
  residuals = as.matrix(residuals)
  fitted = as.matrix(fitted)
  actual = as.matrix(actual)
  
  # Calculate MSE for all covariance types
  mSumMatrix = smatrix(train)
  mseBase =  mean  ( (allts(test)[h,] - preds)^2 )
  
  calibration50 <- NA#checkCalibration(preds, sigma, test, coverage = 0.5)
  calibration80 <- NA#checkCalibration(preds, sigma, test, coverage = 0.8)
  
  bottomIdx <- seq( nrow(mSumMatrix) - ncol(mSumMatrix) +1, nrow(mSumMatrix))
  predsIncoherent = mSumMatrix %*% preds[bottomIdx]
  mseIncoherent =  mean  ( (allts(test)[h,] - predsIncoherent)^2 )
  
  mseBayesDiag =  mean  ( (allts(test)[h,] - bayesReconcFull(preds, residuals, mSumMatrix, "diagonal", sigmas=sigma))^2 )
  mseBayesSample =  mean  ( (allts(test)[h,] - bayesReconcFull(preds, residuals, mSumMatrix, "sam"))^2 )
  mseBayesGlasso =  mean  ( (allts(test)[h,] - bayesReconcFull(preds, residuals, mSumMatrix, "glasso"))^2 )
  mseBayesShr =  mean  ( (allts(test)[h,] - bayesReconcFull(preds, residuals, mSumMatrix, "shr"))^2 )
  
  toc()
  
  if (printResults){
    print(c("MseBase: ",mseBase))
    print(c("MseMint: ", mseCombMintShr))
    print(c("MseIncoherent: ",mseIncoherent))
    print(c("MseCoherentDiagonal: ",mseBayesDiag))
    print(c("MseCoherentSam: ",mseBayesDiag))
    print(c("MseCoherentGlasso: ",mseBayesGlasso))
    print(c("MseCoherentShr: ",mseBayesShr))
  }
  
  if (saveFiles){
    # Saving Files with results
    # Save to file the results, at every iteration
    if (dset=="synthetic"){
      dataFrame <- data.frame(h, iTest, fmethod, synth_n, synthCorrel, corrB2_U, mseBase,mseCombMintSample,
                              mseCombMintShr, mseBayesDiag, mseBayesSample, mseBayesGlasso, mseBayesShr)
      colnames(dataFrame) <- c("h","iTest","fmethod","sampleSize","correlB1_U","correlB2_U",
                               "mseBase","mseMintSample","mseCombMintShr","mseBayesDiag","mseBayesSample",
                               "mseBayesGlasso", "mseBayesShr")
      dset <- paste0(dset,"_correl",synthCorrel,"_n",synth_n)
    }
    else if (dset=="syntheticLarge"){
      dataFrame <- data.frame(h, iTest, fmethod, synth_n, mseBase,mseCombMintSample,
                              mseCombMintShr, mseBayesDiag, mseBayesSample, mseBayesGlasso, mseBayesShr)
      colnames(dataFrame) <- c("h","iTest","fmethod","sampleSize",
                               "mseBase","mseMintSample","mseCombMintShr","mseBayesDiag","mseBayesSample",
                               "mseBayesGlasso","mseBayesShr")
      dset <- paste0("largeSynthetic_n",synth_n)
    }
    else
    {
      dataFrame <- data.frame(h, iTest, fmethod, dset, calibration50, calibration80, 
                              mseBase,mseCombMintSample,mseCombMintShr,mseBayesDiag,mseBayesSample,mseBayesGlasso,mseBayesShr)
    }
    return(dataFrame)
  }
}