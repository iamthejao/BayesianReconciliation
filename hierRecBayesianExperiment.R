list.of.packages <- c("foreach", "doMC", "tictoc")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

source("reconciliationMethods.R")
library(foreach)
library(doMC)
registerDoMC(floor(detectCores()/2))
library(tictoc)

print(paste0("CPUS Available: ", detectCores()))
print(paste0("CPUS Used: ", floor(detectCores()/2)))

# Set of variables usefull for debugging
dset = 'tourism'
synth_n = 10
synthCorrel = 0.5
h = 2
fmethod = 'ets'
iTest = 18
seed = 0
savePredictions=FALSE
saveSamples = FALSE
runPositive = FALSE
enforceKhOne = FALSE
sampleSize = 100000

hierRecBayesianExperiment <- function(dset, h, fmethod="ets", iTest=1, testSize=50,
                        seed=0, synth_n=100, synthCorrel=0.5, # Parameters only for synthetic TS
                        testProbability=0.10,
                        savePredictions=TRUE, saveSamples=FALSE,
                        runPositive=FALSE, enforceKhOne=FALSE, sampleSize=100000) # Reasonable sampling point according to convergence test
  {
  
  # Setting seed and loading dataset
  set.seed(seed)
  hierTs = loadDataset(dset, synth_n=synth_n, synthCorrel=synthCorrel)
  allTs = allts(hierTs)
  columnNames = colnames(allTs)
  numTs <- ncol(allTs)
  alpha <- 0.2
  
  # Update dset name to full specification in case of synthetic data
  if (dset=="synthetic"){
    dset <- paste0(dset,"_correl",synthCorrel,"_seed",seed)
  } else if (dset=="syntheticLarge") {
    dset <- paste0("largeSynthetic_n",synth_n,"_seed",seed)
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
  print(paste0("RUN POSITIVITY CONSTRAINT: ",runPositive))
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
  
  # Flipping upper hierarchy due to residual definition
  residuals[, upperIdx] = residuals[, upperIdx] * -1
  
  # Calculating covariance matrix using shrink
  mCov_shrmint = estimateCovariance(residuals, method='shrmint')
  
  # Covariance for the base case
  mCov_base = matrix(0, nrow=numTs, ncol=numTs)
  diag(mCov_base) = diag(mCov_shrmint) * kh
  #diag(mCov_base) = sigma^2
  
  # Saving predictions
  tic("Generating Data")
  
  # Calculating outputs
  # Corr stands for the correlated case considering cross-covariance. In the paper, it is called pMinT
  # Indep stands for the independent case where cross-covariance = 0. In the paper, it is called Linear Gaussian

  # Correlated upper and bottom residuals
  outBayesCorr = bayesRecon(preds, mSumMatrix, mCov_shrmint, positivity=FALSE, noiseType="correlated",
                            sampleSize=sampleSize,
                            kh=kh)
  
  # Independent (no cross-correlation) upper and bottom residuals
  outBayesIndep = bayesRecon(preds, mSumMatrix, mCov_shrmint, positivity=FALSE, noiseType="independent",
                             sampleSize=sampleSize,
                             kh=kh)
  
  # Future work on positivity constraint
  if (runPositive){
    # Correlated upper and bottom residuals and truncated distribution
    outBayesCorrPos = bayesRecon(preds, mSumMatrix, mCov_shrmint, positivity=TRUE, noiseType="correlated",
                                 sampleSize=sampleSize,
                                 kh=kh)
    
    # Independent upper and bottom residuals and truncated distribution
    outBayesIndepPos = bayesRecon(preds, mSumMatrix, mCov_shrmint, positivity=TRUE, noiseType="independent",
                                  sampleSize=sampleSize,
                                  kh=kh)
  } else {
    outBayesCorrPos = NULL
    outBayesIndepPos = NULL
  }
  
  # In order to calculate the energyScore of the base
  # Sample from a MVN with diagonal covariance matrix
  baseSample = sampleMVN(preds, cbind(mCov_base), sampleSize=sampleSize)
  
  # BottomUp sampling, sample from the bottom and apply BottomUp method.
  mSigmaB = cbind(mCov_shrmint[bottomIdx, bottomIdx])
  bottomUpSample = sampleMVN(preds[bottomIdx], kh*mSigmaB, sampleSize=sampleSize)
  bottomUpSample = cppMatMult(bottomUpSample, mSumMatrixT)
  outBottomUp = list("preds"=predsBottomUp)
  
  bottomUpPosSample = NULL
  if (runPositive){
    # BottomUp Positivity sampling
    bottomUpPosSample = sampleMVN(preds[bottomIdx], kh*mSigmaB, sampleSize=sampleSize, positivity=TRUE)
    bottomUpPosSample = cppMatMult(bottomUpPosSample, mSumMatrixT)
    bottomUpPosMeanSample = colMeans(bottomUpPosSample)
    bottomUpPosMedianSample = colMedians(bottomUpPosSample)
    outBottomUp["predsMeanPosSample"] = bottomUpPosMeanSample
    outBottomUp["predsMedianPosSample"] = bottomUpPosMedianSample
  }
  
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
  
  statsBottomUpPos = NULL
  statsBayesCorrPos = NULL
  statsBayesIndepPos = NULL
  if (runPositive){
    statsBottomUpPos = calculateStatistics(target,
                                           bottomUpPosMeanSample, bottomUpPosMedianSample,
                                           bottomUpPosSample, groupings,
                                           suffix="_BottomUpPos")
    statsBayesCorrPos = calculateStatistics(target,
                                            outBayesCorrPos$coherentPredsSample, outBayesCorrPos$coherentPredsMedianSample,
                                            outBayesCorrPos$sample, groupings,
                                            suffix="_CorrPos")
    statsBayesIndepPos = calculateStatistics(target,
                                             outBayesIndepPos$coherentPredsSample, outBayesIndepPos$coherentPredsMedianSample,
                                             outBayesIndepPos$sample, groupings,
                                             suffix="_IndepPos")
    
  }
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
    outBayesCorrPos$sample = NULL
    outBayesIndep$sample = NULL
    outBayesIndepPos$sample = NULL
    savePredictions = paste0("predictions/", paste(dset, fmethod, h, kh, iTest, sep="_"),".RData")
    save(outBayesCorr, outBayesCorrPos, outBayesIndep, outBayesIndepPos, preds, outBottomUp, target, file=savePredictions)
  }
  
  # At this point, we double check if the bayesian reconciliation method pMinT matches
  # the results of MinT from Hyndman's library
  # With a certain probability (this step can be very costly for big datasets since we dont have the models cached)
  # normDiff is saved to be evaluated later
  rmseMint = NaN
  normDiff = NaN
  set.seed(Sys.time())
  if (runif(1) <= testProbability){
    fcastMintShr = tryCatch(
      {
        set.seed(seed)
        print("Performing MinT on library for double check")
        tic("Running MinT")
        if (fmethod == "ets"){
          fcastMintShr <-forecast.gts(train, h = h, method = "comb", weights="mint",
                                      fmethod=fmethod, covariance="shr", additive.only = TRUE,
                                      parallel=parallel, num.cores=5)
        } else {
          fcastMintShr <-forecast.gts(train, h = h, method = "comb", weights="mint",
                                      fmethod=fmethod, covariance="shr",
                                      parallel=parallel, num.cores=5)
        }
        toc()
        fcastMintShr = allts(fcastMintShr)[h,]
        fcastMintShr},
      error=function(cond){
      print("MinT ERROR")
      print(cond)
      return(NULL)}
      , finally={}
      )
    
    if (is.null(fcastMintShr)){
      rmseMint = -1
      normDiff = -1
    } else {
      rmseMint  <- sqrt(mean((allts(test)[h,] - fcastMintShr)^2))
      diff = (fcastMintShr - outBayesCorr$coherentPreds)
      normDiff = norm(diff, "2")
      normFcast = norm(fcastMintShr, "2")
      if ((normDiff/normFcast) > 0.01){
        print("Bayes and MinT slightly different.")
        print(statsBayesCorr$rmseMeanAll_Corr)
        print(rmseMint)
        #stop("MinT and Bayes result diverged.")
      }
    }
  }
  
  # Saving outputs in a dataframe and sending out
  statistics = c(statsBase, statsBottomUp, statsBottomUpPos, statsBayesCorr, statsBayesCorrPos, statsBayesIndep, statsBayesIndepPos)
  output = list(dset=dset, h=h, iTest=iTest, fmethod=fmethod, rmseMint=rmseMint, normDiffMint=normDiff)
  dataFrame = data.frame(c(output, statistics))
  
  return(dataFrame)
}