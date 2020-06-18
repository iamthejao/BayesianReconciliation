list.of.packages <- c("foreach", "doMC", "tictoc", "forecast")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

source("reconciliationMethods.R")
library(foreach)
library(doMC)
library(forecast)
registerDoMC(floor(detectCores()/2))
library(tictoc)

print(paste0("CPUS Available: ", detectCores()))
print(paste0("CPUS Used: ", floor(detectCores()/2)))

# Set of variables usefull for debugging
# dset = 'infantgts'
# synth_n = 10
# synthCorrel = 0.5
# h = 1
# fmethod = 'ets'
# iTest = 10
# seed = 0
# savePredictions=FALSE
# saveSamples = FALSE
# runPositive = FALSE
# enforceKhOne = FALSE
# sampleSize = 100000

hierRecBayesianExperiment <- function(dset, h, fmethod="ets", iTest=1, testSize=50,
                        seed=0, synth_n=100, savePredictions=TRUE, saveSamples=FALSE,
                        enforceKhOne=FALSE, sampleSize=100000) # Reasonable sampling point according to convergence test
  {
  
  # Setting seed and loading dataset
  set.seed(seed)
  hierTs = loadDataset(dset, synth_n=synth_n, synthCorrel=synthCorrel)
  allTs = allts(hierTs)
  columnNames = colnames(allTs)
  numTs <- ncol(allTs)
  alpha <- 0.2
  
  # Update dset name to full specification in case of synthetic data
  if (dset=="syntheticLarge") {
    dset <- paste0(dset, "_n", synth_n)
  }
  
  # Printing experiment
  print("=========================")
  print(paste0("DSET: ",dset))
  print(paste0("H: ",h))
  print(paste0("Kh is 1: ",enforceKhOne))
  print(paste0("METHOD: ",fmethod))
  print(paste0("iTest: ",iTest))
  print(paste0("SAMPLE SIZE: ",sampleSize))
  print(paste0("SEED: ",seed))
  print(paste0("SAVE PREDS: ",savePredictions))
  print(paste0("SAVE SAMPLES: ",saveSamples))
  print("=========================")
  print("")

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
  
  # Setting target
  target = allts(test)[h,]
  
  # Creating matrices to save results
  # Making vector and matrixes for storaging values in shared memory for parallel calculation
  # Tourism or ARIMA, otherwise it is quicker to do single core due to overhead
  parallel=FALSE
  if ((dset == 'tourism') || (fmethod == 'arima')){
   parallel=TRUE
  }
  
  # Predictions and sigma from model
  preds <- matrix(ncol=numTs, nrow=1)
  sigma <- matrix(ncol=numTs, nrow=1)
  
  # Fit from model
  fitted <- matrix(nrow=dim(allTsTrain)[1], ncol = numTs)
  colnames(fitted) <- colnames(allTsTrain)
  
  # Iteration on parallel or single core
  if (parallel){
    '%iter%' <- foreach::'%dopar%'
  } else {
    '%iter%' <- foreach::'%do%'
  }
    
  print("Starting series fits")
  tic("fit models")
  
  # Fitting methods with ARIMA or ETS
  # fitOrReadTS function fits and tries to save model in a file for future uses
  # or reads if it already exists
  # Here we fit a model for each timeseries in the tree
  outputs <- foreach (i = 1:numTs) %iter% {
    
    ans = fitOrReadTS(allTsTrain[,i], fmethod)
    model = ans$model
    
    tmp <- forecast(model, h=h, PI=TRUE, level=1-alpha)
    
    f = model$fitted
    p = tmp$mean[h]
    
    # Saving forecast std, it is not currently used but could be used in future work
    s = abs((tmp$mean[h] - tmp$upper[h])/(qnorm(alpha/2)))
    
    out = list(fitted=f, preds=p, sigma=s, flag=ans$flag)
    return(out)
  }
  
  # Amount loaded from disk
  loadCount = sum(as.numeric(lapply(outputs, function(l){l$flag})))
  print(paste0("Timeseries loaded from disk: ", (loadCount/numTs)*100, "%"))
  toc()
  
  # Populating matrix of predictions and fits
  for (i in 1:numTs) {
    fitted[,i] <- outputs[[i]]$fitted
    preds[i] <- outputs[[i]]$preds
    sigma[i] <- outputs[[i]]$sigma
  }
  
  tic("Calculating statistics") 
  
  # Kh assumption
  if (enforceKhOne){
    kh = 1
  } else {
    kh = h 
  }
  
  # Sum matrix, indices and kh
  mSumMatrix = smatrix(train)
  mSumMatrixT = t(mSumMatrix)
  bottomIdx <- seq( nrow(mSumMatrix) - ncol(mSumMatrix) +1, nrow(mSumMatrix))
  upperIdx <- setdiff(1:nrow(mSumMatrix),bottomIdx)

  # Preds BottomUp
  predsBottomUp = mSumMatrix %*% preds[bottomIdx]
  
  # Recovering residuals
  y = aggts(train)
  residuals = (y - fitted)
  
  # Calculating covariance matrix using shrink
  # Using the method 'shrmint' should yield the **same** (up to a certain numerical precision) 
  # result for the coherent predictions when comparing pMinT (our implementation)
  # and MinT from the forecast.gts library.
  mCov_shrmint = estimateCovariance(residuals, method='shrmint')
  
  # Covariance for the base case
  mCov_base = matrix(0, nrow=numTs, ncol=numTs)
  diag(mCov_base) = diag(mCov_shrmint) * kh
  
  # Saving predictions
  tic("Generating Data")
  
  # Calculating outputs
  # Corr stands for the correlated case considering cross-covariance. In the paper, it is called pMinT
  # Indep stands for the independent case where cross-covariance = 0. In the paper, it is called Linear Gaussian

  # Correlated upper and bottom residuals
  outBayesCorr = bayesRecon(preds, mSumMatrix, mCov_shrmint, reconType="pmint",
                            sampleSize=sampleSize,
                            kh=kh)
  
  # Independent (no cross-correlation) upper and bottom residuals
  outBayesIndep = bayesRecon(preds, mSumMatrix, mCov_shrmint, reconType="lg",
                             sampleSize=sampleSize,
                             kh=kh)
  
  # In order to calculate the energyScore of the base
  # Sample from a MVN with diagonal covariance matrix
  baseSample = sampleMVN(preds, cbind(mCov_base), sampleSize=sampleSize)
  
  # BottomUp sampling, sample from the bottom and apply BottomUp method.
  mSigmaB = cbind(mCov_shrmint[bottomIdx, bottomIdx])
  bottomUpSample = sampleMVN(preds[bottomIdx], kh*mSigmaB, sampleSize=sampleSize)
  bottomUpSample = cppMatMult(bottomUpSample, mSumMatrixT)
  outBottomUp = list("preds"=predsBottomUp)
  toc()
  
  tic("Calculating statistics")
  # Calculating Statistics
  
  # Groupings for statistics
  # A group is a list of indexes with a specific label
  # For example, "Upper" and the indexes for the upper timeseries
  groupings = getGroupingsIndex(hierTs)
  
  # RMSE, MAE, ES for every group
  statsBase = calculateStatistics(target,
                                  preds, preds,
                                  baseSample, groupings,
                                  suffix="_Base")
  
  statsBottomUp = calculateStatistics(target,
                                      predsBottomUp, predsBottomUp,
                                      bottomUpSample, groupings,
                                      suffix="_BottomUp")
  
  statsBayesCorr = calculateStatistics(target,
                                       outBayesCorr$coherentPreds, outBayesCorr$coherentPredsMedianSample,
                                       outBayesCorr$sample, groupings,
                                       suffix="_Corr")
  
  statsBayesIndep = calculateStatistics(target,
                                        outBayesIndep$coherentPreds, outBayesIndep$coherentPredsMedianSample,
                                        outBayesIndep$sample, groupings,
                                        suffix="_Indep")
  toc()
  
  # Takes a lot of space depending on sampleSize
  if (saveSamples) {
    saveSamples = paste0("samples/", paste(dset, fmethod, h, iTest, sep="_"),".RData")
    corrSample = outBayesCorr$sample
    corrPosSample = outBayesCorrPos$sample
    indepSample = outBayesIndep$sample
    indepPosSample = outBayesIndepPos$sample
    save(residuals,bottomUpSample, bottomUpPosSample, baseSample, corrSample, corrPosSample, indepSample, indepPosSample,
         file=saveSamples)
  }
  
  # Saving predictions, posterior mean, posterior covariance, etc from bayesRecon method.
  # All but samples for storage reasons
  if (savePredictions){
    # Removing samples in order to take less space
    outBayesCorr$sample = NULL
    outBayesIndep$sample = NULL
    savePredictions = paste0("predictions/", paste(dset, fmethod, h, kh, iTest, sep="_"),".RData")
    save(outBayesCorr, outBayesIndep, preds, outBottomUp, target, file=savePredictions)
  }
  
  # Saving outputs in a dataframe and sending out
  statistics = c(statsBase, statsBottomUp, statsBayesCorr, statsBayesIndep)
  output = list(dset=dset, h=h, iTest=iTest, fmethod=fmethod, khOne=enforceKhOne, seed=seed)
  dataFrame = data.frame(c(output, statistics))
  
  return(dataFrame)
}