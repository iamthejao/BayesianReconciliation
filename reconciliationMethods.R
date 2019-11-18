library(hts)
library(huge) #covariance matrix via glasso
library(SHIP) #shrinkage of covarianca matrix
library(matlib)
source("loadTourism.R")

estimateCovariance <- function(residuals, method="diagonal", diagonal=NULL, labels=NULL){
  #prior covariance for the bottom time series
  #bottomVar <- sigma[bottomIdx]^2
  #bottomResiduals <- residuals[,bottomIdx]
  if (method=="diagonal"){
    if (is.null(diagonal) == FALSE){
      covar = diag(diagonal)
    }
    else{
      covar <- diag(apply(residuals, 2, var))  
    }
  }
  else if (method=="sam"){
    #the covariances are the covariances of the time series
    #the variances are the variances of the forecasts, hence the variances of the residuals
    covar <- cov(residuals)
  }
  else if (method=="glasso"){
    #the covariances are the covariances of the time series
    #the variances are the variances of the forecasts, hence the variances of the residuals
    out.glasso <- huge(residuals, method = "glasso", cov.output = TRUE)
    out.select <- huge.select(out.glasso, criterion = "ebic")
    covar <- out.select$opt.cov
  }
  else if (method=="shr"){
    if (is.null(diagonal) == FALSE){
      sigmaDiag = diag(diagonal)
    }
    else{
      sigmaDiag <- diag(apply(residuals, 2, var))  
    }
    covar <-  shrink.estim(residuals, tar=build.target(residuals,type="D"))[[1]]
  }
  # Adding labels
  if (is.null(labels) == FALSE){
    colnames(covar) <- labels
    rownames(covar) <- labels
  }
  return(covar)
}

bayesReconcFull <- function(preds, residuals, mSumMatrix, covarMethod="shr", sigmas=NULL){
  
  # Defining sumMatrix and idxs
  bottomIdx <- seq( nrow(mSumMatrix) - ncol(mSumMatrix) +1, nrow(mSumMatrix))
  upperIdx <- setdiff(1:nrow(mSumMatrix),bottomIdx)
  
  # Defining covariance
  mCovar = estimateCovariance(residuals, method=covarMethod, labels=colnames(residuals), diagonal=sigmas)
  mA = mSumMatrix[upperIdx, ]
  mAtr = t(mA)
  
  # Taking predictions
  vBottomPreds = preds[bottomIdx]
  vTopPreds = preds[upperIdx]
  
  # Priors
  vMeans = mSumMatrix %*% vBottomPreds
  vPriorMean = vMeans[bottomIdx]
  vIncoherence = (vTopPreds - vMeans[upperIdx])
  
  # Base covariance
  mSigmaB = mCovar[bottomIdx, bottomIdx]
  
  # Upper covariance
  mSigmaU = mCovar[upperIdx, upperIdx]
  
  # Cross covariance
  mM = mCovar[bottomIdx, upperIdx]
  mMtr = mCovar[upperIdx, bottomIdx]
  
  # Update rule
  
  # Gain is n_base by + n_base
  mP1Gain = ((mSigmaB %*% mAtr) + mM)
  mP2Gain = (mSigmaU + (mA %*% mSigmaB %*% mAtr) + (mA%*%mM) + (mMtr%*%mAtr))
  mGain = mP1Gain %*% solve(mP2Gain + 1e-6*diag(mP2Gain))
  
  # Posterior mean
  vPosteriorMean = vPriorMean + mGain %*% vIncoherence
  
  # Posterior Cov
  mSigmaBp = mSigmaB - mGain %*% (mA%*%mSigmaB+mMtr)
  
  # Coherent predictions
  vCoherentPreds = mSumMatrix %*% vPosteriorMean
  
  return(vCoherentPreds)
}

#covariance can be "diagonal", "sam" or "glasso"
bayesRecon <- function (preds, residuals, S, covariance="shr", sigma=NULL){
  bottomIdx <- seq( nrow(S) - ncol(S) +1, nrow(S))
  upperIdx <- setdiff(1:nrow(S),bottomIdx)
  
  #prior mean and covariance of the bottom time series
  priorMean <- preds[bottomIdx]
  Y_vec <- preds[upperIdx]
  
  
  #prior covariance for the bottom time series
  bottomVar <- sigma[bottomIdx]^2
  bottomResiduals <- residuals[,bottomIdx]
  if (covariance=="diagonal"){
    priorCov <- diag(bottomVar)
  }
  else if (covariance=="sam"){
    #the covariances are the covariances of the time series
    #the variances are the variances of the forecasts, hence the variances of the residuals
    priorCov <- cov(bottomResiduals)
  }
  else if (covariance=="glasso"){
    #the covariances are the covariances of the time series
    #the variances are the variances of the forecasts, hence the variances of the residuals
    out.glasso <- huge(bottomResiduals, method = "glasso", cov.output = TRUE)
    out.select <- huge.select(out.glasso, criterion = "ebic")
    priorCov <- out.select$opt.cov
  }
  else if (covariance=="shr"){
    sigmaDiag <- diag(bottomVar)
    priorCov <-  shrink.estim(bottomResiduals, tar=build.target(bottomResiduals,type="D"))[[1]]
  }
  
  upperVar <- sigma[upperIdx]^2
  #covariance for the upper time series; we need managing separately the case where only a single time series is present
  #as diag will try to create a matrix of size upperVar instead.
  upperResiduals <- residuals[,upperIdx]
  if (length(upperIdx)==1) {
    Sigma_y <- upperVar
  }
  
  else if (covariance=="diagonal"){
    Sigma_y <- diag(upperVar)
  }
  #if we only one upper time series, there is no covariance matrix to be estimated. 
  else if (covariance=="glasso") {
    #get variance and covariance of the residuals
    out.glasso <- huge(upperResiduals, method = "glasso", cov.output = TRUE)
    out.select <- huge.select(out.glasso, criterion = "ebic")
    Sigma_y <- out.select$opt.cov
  }
  else if (covariance=="sam") {
    #get variance and covariance of the residuals
    Sigma_y <- cov(upperResiduals)
  }
  else if (covariance=="shr") {
    sigma_y_diag <- diag(upperVar)
    Sigma_y <-  shrink.estim(upperResiduals, tar=build.target(upperResiduals,type="D"))[[1]]
  }
  #==updating
  #A explains how to combin the bottom series in order to obtain the
  # upper series
  
  #if upperIdx contains a single row, R behaves oddily; hence we need to manually manage that case.
  if (length(upperIdx)==1){
    A <- cbind(S[upperIdx,])
  }
  else {
    A <- t(S[upperIdx,])
  }
  
  M <- ncol ( t(A) %*% priorCov %*% A + Sigma_y )
  correl <- priorCov %*% A %*%
    solve (t(A) %*% priorCov %*% A + Sigma_y + 1e-6*diag(M))
  
  postMean <- priorMean + correl  %*%
    (Y_vec - t(A) %*% priorMean)
  bayesPreds <- buReconcile(postMean, S, predsAllTs = FALSE)
  return(bayesPreds)
}

loadDataset <- function(dset, synth_n=100, synthCorrel=0.5){
  
  feasibleDset <- c("infantgts", "tourism", "synthetic", "syntheticLarge")
  if (! (dset %in% feasibleDset)){
    print("feasible dset are:")
    print(feasibleDset)
    stop ("wrong dset supplied" )
  }
  
  if (dset=="tourism"){
    hierTs <- loadTourism()
  }
  else if (dset=="infantgts"){
    hierTs <- infantgts
  }
  else if (dset=="synthetic"){
    source("drawSmallHierarchy.R")
    #we generate the hierarchy with two bottom time series
    #training and the test
    
    #synth_n=100
    #synthCorrel=0.5
    #h=1
    
    listSynth <- artificialTs(n=synth_n + h, correl = synthCorrel)
    synthTs <- listSynth$bottomTs
    corrB2_U <- listSynth$corrB2_U
    y=ts(synthTs, frequency = 1)
    hierTs <- hts(y, bnames = colnames(y))
  }
  else if (dset=="syntheticLarge"){
    source("drawLargeHierarchy.R")
    #we generate the hierarchy with *four* bottom time series
    synthTs <- simulFourBottom(n=synth_n)
    y=ts(synthTs, frequency = 1)
    hierTs <- hts(y, bnames = colnames(y), characters=c(1,1))
  }
  return(hierTs)
}

hierMse <- function (htsPred, htsActual, h) {
  #receives two hts objects, containing  forecast and actual value.
  #computes the mse for the whole hierarchy.
  mse <- mean  ( (allts(htsPred)[h,] - allts(htsActual)[h,])^2 )
  return (mse)
}

#The buReconcile function computes the bu prediction given the predictions (1 x tot time series) and the S matrix
#(tot time series X bottom time series)
#predsAllTs is a flag: is set to true, the input preds contains predictions for all the hierarchy
#and the function retrieves the bottom series; if set to false, this is not needed
#as preds only contains only bottom time series
buReconcile <- function (preds,S, predsAllTs = FALSE) {
  bottomPreds <- preds
  if (predsAllTs) {
    #retrieves the bottom prediction from all predictions
    upperIdx <- 1 : (nrow(S) - ncol(S))
    bottomIdx <- setdiff (1:nrow(S), upperIdx)
    bottomPreds <- preds [,bottomIdx]
  }
  
  buPreds <- preds
  
  #nrow(S) is the total number of time series
  for (i in 1:nrow(S)){
    buPreds[i] <- S[i,] %*% bottomPreds
  }
  return (buPreds)
}

#check the calibration of the prediction interval with coverage (1-currentAlpha)
checkCalibration <- function(h, preds,sigmas,htsActual,coverage){
  stdQuant <- abs(qnorm((1-coverage)/2))
  included <- vector(length = length(preds))
  actual <- allts(htsActual)[h,]
  for (ii in seq_along(preds)){
    upper <- preds[ii] + stdQuant * sigmas[ii]
    lower <- preds[ii] - stdQuant * sigmas[ii]
    included[ii] <- (actual[ii] > lower) &  (actual[ii] < upper)
  }
  return (mean(included))
}