library(hts)
library(huge) #covariance matrix via glasso
library(SHIP) #shrinkage of covarianca matrix
library(matlib)
library(MASS)
library(scoringRules)
library(tmvtnorm)
source("loadTourism.R")

rowVars <- function(x, ...) {
  rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
}

colVars <- function(x, ...) {
  x = t(x)
  tmp = rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
  return(t(tmp))
}

# THIS IS WRONG FOR TRUNCATED, I HAVE TO SAMPLE FROM THE TRUNCATED BASE DIST AND THEN CONTINUE
energyScore <- function(target, mean, covariance, sample=50000, truncated=FALSE){
  
  # This gives an error saying that the covariance matrix is not SPD
  #mCholCov = chol(covariance)
  #dims = dim(covariance)[1]
  #rnormal = matrix(rnorm(dims * twoSample), twoSample, dims)
  #mMean = t(do.call(cbind, foreach(i=1:(2*sample))%do%{mean}))
  #samples = (rnormal %*% mCholCov) + mMean
  #rtmvnorm
  
  twoSample = 2 * sample
  #mTarget = t(do.call(cbind, foreach(i=1:sample)%do%{target})) slow
  if (truncated){
    samples = rtmvnorm(twoSample, mean=as.numeric(mean), sigma=covariance,
             lower=rep(0, length(mean)), algorithm="rejection")
  }else{
    samples = mvrnorm(n=twoSample, mu=mean, Sigma=covariance)
  }
  
  samples1 = samples[(1:sample),]
  samples2 = samples[((sample+1):twoSample), ]
  
  term1 = t(apply(samples1, 1, function(x){x-target}))#(samples1 - mTarget) slow
  normst1 = apply(term1, 1, function(x){norm(x,"2")})
  
  term2 = (samples1 - samples2)
  normst2 = apply(term2, 1, function(x){norm(x,"2")})
  
  score = (sum(normst1) - 0.5 * sum(normst2))/(sample-1)
  scoreLib = NaN#scoringRules::es_sample(target, t(samples))
  out = list(score=score, scoreLib=scoreLib)
  return(out)
}

estimateCovariance <- function(residuals, method="diagonal", diagonal=NULL, labels=NULL){
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
  else if (method == "sammint") {
    n <- nrow(residuals)
    covar <- crossprod(residuals)/n
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
      mTar = diag(diagonal)
    }
    else{
      # Both expressions lead to the same result
      mTar = build.target(residuals,type="D")#diag(apply(residuals, 2, var))  
    }
    covar = shrink.estim(residuals, tar=mTar)[[1]]
  }
  else if (method=="shrmint"){
    
    lowerD <-function (x) 
    {
      n <- nrow(x)
      return(diag(apply(x, 2, crossprod)/n))
    }
    
    tar <- lowerD(residuals)
    covar <- hts:::shrink.estim(residuals, tar)[[1]]
  }
  else {
    print("Method not known.")
  }
  # Adding labels
  if (is.null(labels) == FALSE){
    colnames(covar) <- labels
    rownames(covar) <- labels
  }
  return(covar)
}

bayesReconFull <- function(preds, mSumMatrix, mCovar, positivity=FALSE){
  
  # Defining sumMatrix and idxs
  bottomIdx <- seq( nrow(mSumMatrix) - ncol(mSumMatrix) +1, nrow(mSumMatrix))
  upperIdx <- setdiff(1:nrow(mSumMatrix),bottomIdx)
  
  # Defining covariance
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
  mSigmaBp = round(mSigmaBp, 2)
  
  # Coherent predictions
  vCoherentPreds = mSumMatrix %*% vPosteriorMean
  mCoherentVariance = mSumMatrix %*% mSigmaBp %*% mSumMatrixT
  colnames(mCoherentVariance) = colnames(mCovar)
  rownames(mCoherentVariance) = rownames(mCovar)
  
  out = list(posteriorMean=vPosteriorMean, posteriorVariance=mSigmaBp,
             coherentPreds=vCoherentPreds, coherentCovariance=mCoherentVariance)
  
  if (positivity){
    vPosteriorMeanTrunc = mtmvnorm(mean=as.numeric(vPosteriorMean),
                                   sigma=mSigmaBp,
                                   lower=rep(0, length(vPosteriorMean)),
                                   doComputeVariance=FALSE)$tmean
    out$posteriorMeanTrunc = vPosteriorMeanTrunc
    out$coherentPredsTrunc = mSumMatrix %*% vPosteriorMeanTrunc
  }
  return(out)
}

bayesReconNotFull <- function(preds, mSumMatrix, mCovar, positivity=FALSE){
  
  # Defining sumMatrix and idxs
  bottomIdx <- seq( nrow(mSumMatrix) - ncol(mSumMatrix) +1, nrow(mSumMatrix))
  upperIdx <- setdiff(1:nrow(mSumMatrix),bottomIdx)
  
  # Defining covariance
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
  
  # Update rule
  # Gain is n_base by + n_base
  mP1Gain = (mSigmaB %*% mAtr)
  mP2Gain = (mSigmaU + (mA %*% mSigmaB %*% mAtr))
  mGain = mP1Gain %*% solve(mP2Gain + 1e-6*diag(mP2Gain))
  
  # Posterior mean
  vPosteriorMean = vPriorMean + mGain %*% vIncoherence
  
  # Posterior Cov
  mSigmaBp = mSigmaB - mGain %*% (mSigmaU + mA %*% mSigmaB %*% mAtr) %*% t(mGain)
  # BOTH LEAD TO SAME RESULT
  #test <- mSigmaB - mGain %*% mA %*% mSigmaB
  
  # Coherent predictions
  vCoherentPreds = mSumMatrix %*% vPosteriorMean
  # Coherent predictions
  vCoherentPreds = mSumMatrix %*% vPosteriorMean
  mCoherentVariance = mSumMatrix %*% mSigmaBp %*% mSumMatrixT
  colnames(mCoherentVariance) = colnames(mCovar)
  rownames(mCoherentVariance) = rownames(mCovar)
  
  out = list(posteriorMean=vPosteriorMean, posteriorVariance=mSigmaBp,
             coherentPreds=vCoherentPreds,
             coherentCovariance = mCoherentVariance)
  if (positivity){
    vPosteriorMeanTrunc = mtmvnorm(mean=as.numeric(vPosteriorMean),
                                   sigma=mSigmaBp,
                                   lower=rep(0, length(vPosteriorMean)),
                                   doComputeVariance=FALSE)$tmean
    out$posteriorMeanTrunc = vPosteriorMeanTrunc
    out$coherentPredsTrunc = mSumMatrix %*% vPosteriorMeanTrunc
  }
  return(out)
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

# #covariance can be "diagonal", "sam" or "glasso"
# # old function for reconciliation, leaving here to check if I forget something
# old_bayesReconNotFull2 <- function (preds, S, mCovar){
#   
#   bottomIdx <- seq( nrow(S) - ncol(S) +1, nrow(S))
#   upperIdx <- setdiff(1:nrow(S),bottomIdx)
#   
#   #prior mean and covariance of the bottom time series
#   priorMean <- preds[bottomIdx]
#   Y_vec <- preds[upperIdx]
#   
#   #prior covariance for the bottom time series
#   priorCov = mCovar[bottomIdx, bottomIdx]
#   
#   #covariance for the upper time series
#   Sigma_y = mCovar[upperIdx, upperIdx]
#   
#   #if upperIdx contains a single row, R behaves oddily; hence we need to manually manage that case.
#   if (length(upperIdx)==1){
#     A <- cbind(S[upperIdx,])
#   }
#   else {
#     A <- t(S[upperIdx,])
#   }
#   
#   M <- ncol ( t(A) %*% priorCov %*% A + Sigma_y )
#   G <- priorCov %*% A %*%
#     solve (t(A) %*% priorCov %*% A + Sigma_y + 1e-6*diag(M))
#   
#   postMean <- priorMean + G  %*%
#     (Y_vec - t(A) %*% priorMean)
#   bayesPreds <- buReconcile(postMean, S, predsAllTs = FALSE)
#   
#   posteriorCov = NULL#priorCov - G %*% (Sigma_y + A %*% )
#   
#   out = list(posteriorMean=postMean, posteriorVariance=posteriorCov, coherentPreds=bayesPreds)
#   return(bayesPreds)
# }
# 
# #covariance can be "diagonal", "sam" or "glasso"
# # Old function for reconciliation
# old_bayesRecon <- function (preds, residuals, S, covariance="shr", sigma=NULL){
#   bottomIdx <- seq( nrow(S) - ncol(S) +1, nrow(S))
#   upperIdx <- setdiff(1:nrow(S),bottomIdx)
#   
#   #prior mean and covariance of the bottom time series
#   priorMean <- preds[bottomIdx]
#   Y_vec <- preds[upperIdx]
#   
#   
#   #prior covariance for the bottom time series
#   bottomVar <- sigma[bottomIdx]^2
#   bottomResiduals <- residuals[,bottomIdx]
#   if (covariance=="diagonal"){
#     priorCov <- diag(bottomVar)
#   }
#   else if (covariance=="sam"){
#     #the covariances are the covariances of the time series
#     #the variances are the variances of the forecasts, hence the variances of the residuals
#     priorCov <- cov(bottomResiduals)
#   }
#   else if (covariance=="glasso"){
#     #the covariances are the covariances of the time series
#     #the variances are the variances of the forecasts, hence the variances of the residuals
#     out.glasso <- huge(bottomResiduals, method = "glasso", cov.output = TRUE)
#     out.select <- huge.select(out.glasso, criterion = "ebic")
#     priorCov <- out.select$opt.cov
#   }
#   else if (covariance=="shr"){
#     sigmaDiag <- diag(bottomVar)
#     priorCov <-  shrink.estim(bottomResiduals, tar=build.target(bottomResiduals,type="D"))[[1]]
#   }
#   
#   upperVar <- sigma[upperIdx]^2
#   #covariance for the upper time series; we need managing separately the case where only a single time series is present
#   #as diag will try to create a matrix of size upperVar instead.
#   upperResiduals <- residuals[,upperIdx]
#   if (length(upperIdx)==1) {
#     Sigma_y <- upperVar
#   }
#   
#   else if (covariance=="diagonal"){
#     Sigma_y <- diag(upperVar)
#   }
#   #if we only one upper time series, there is no covariance matrix to be estimated. 
#   else if (covariance=="glasso") {
#     #get variance and covariance of the residuals
#     out.glasso <- huge(upperResiduals, method = "glasso", cov.output = TRUE)
#     out.select <- huge.select(out.glasso, criterion = "ebic")
#     Sigma_y <- out.select$opt.cov
#   }
#   else if (covariance=="sam") {
#     #get variance and covariance of the residuals
#     Sigma_y <- cov(upperResiduals)
#   }
#   else if (covariance=="shr") {
#     sigma_y_diag <- diag(upperVar)
#     Sigma_y <-  shrink.estim(upperResiduals, tar=build.target(upperResiduals,type="D"))[[1]]
#   }
#   #==updating
#   #A explains how to combin the bottom series in order to obtain the
#   # upper series
#   
#   #if upperIdx contains a single row, R behaves oddily; hence we need to manually manage that case.
#   if (length(upperIdx)==1){
#     A <- cbind(S[upperIdx,])
#   }
#   else {
#     A <- t(S[upperIdx,])
#   }
#   
#   M <- ncol ( t(A) %*% priorCov %*% A + Sigma_y )
#   correl <- priorCov %*% A %*%
#     solve (t(A) %*% priorCov %*% A + Sigma_y + 1e-6*diag(M))
#   
#   postMean <- priorMean + correl  %*%
#     (Y_vec - t(A) %*% priorMean)
#   bayesPreds <- buReconcile(postMean, S, predsAllTs = FALSE)
#   return(bayesPreds)
# }