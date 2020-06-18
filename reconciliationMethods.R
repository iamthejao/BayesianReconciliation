list.of.packages <- c("hts", "huge", "SHIP", "matlib", "MASS", "tmvtnorm", "truncnorm", "matrixStats", "wordspace",
                      "digest", "Rcpp", "corpcor")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(hts)
library(huge) #covariance matrix via glasso
library(SHIP) #shrinkage of covariance matrix
library(matlib)
library(MASS)
library(tmvtnorm)
library(truncnorm)
library(matrixStats)
library(wordspace)
library(digest)
library(Rcpp)
library(corpcor)
sourceCpp("cppAlgebra.cpp")
source("loadTourism.R")

# Make groups from hierarchical timeseries
getGroupingsIndex <- function(hierTs){
  out = tryCatch(
    {gm=get_groups(hierTs)
    levels = apply(gm, 1, function(x)length(unique(x)))
    levelLabels = names(levels)
    list(levels=levels, labels=levelLabels)},
    error=function(cond){
      levels = lapply(hierTs$labels, function(i){length(i)})
      levelLabels = names(levels)
      levelLabels[length(levels)] = "Bottom"
      return(list(levels=levels, labels=levelLabels))})
  
  levels = out$levels
  levelLabels = out$labels
  curr = 1
  depth = 1
  levelIndexes = foreach (i = cumsum(levels)) %do% {
    serie = curr:i
    curr = i+1
    size = length(serie)
    
    name = levelLabels[depth]
    if (name != "Bottom"){
      name = paste0("Depth", depth, "(", name,")")
    }
    levelLabels[depth] = name
    
    depth = depth + 1
    serie
  }
  
  mSumMatrix = smatrix(hierTs)
  allIdx = 1:ncol(allts(hierTs))
  bottomIdx <- seq( nrow(mSumMatrix) - ncol(mSumMatrix) +1, nrow(mSumMatrix))
  upperIdx <- setdiff(1:nrow(mSumMatrix),bottomIdx)
  grp = list()
  grp["indexes"]=list(c(list(allIdx, upperIdx), levelIndexes))
  grp["labels"]= list(c(c("All", "Upper"), levelLabels))
  return(grp)
}

# Receives a TS object (series), fit, hash and store or read
# In order to read, function must receive the EXACT same object, including
# Precision, col and row names, type
# Otherwise wont recognize
fitOrReadTS <- function(series, fmethod, location="storage/") {
  
  key = digest(series, algo="md5")
  file = paste0(location, key, "_", fmethod, ".model")
  
  FOUND=0
  if(file.exists(file)){
    model = readRDS(file)
    FOUND=1
  } else {
    if (fmethod=="ets"){
      model <- ets(series, additive.only = TRUE)
    }
    else if (fmethod=="arima"){
      model <- auto.arima(series)
    }
    saveRDS(model, file=file)
  }
  return(list(model=model, flag=FOUND))
}

# Calculate rmse, mae and es for each group returns a list with the calculations
calculateStatistics <- function(target, mean, median, sample, groupings, suffix=""){
  
  idxs = groupings[['indexes']] #list(1:length(target), upperIdx, bottomIdx)
  labels = groupings[['labels']] #c("All", "Upper", "Bottom")
  out = list()
  
  for(i in 1:length(labels)){
    
    label = labels[i]
    idx = idxs[[i]]
    
    rmseMean = sqrt(mean((target[idx] - mean[idx])^2))
    rmseMedian = sqrt(mean((target[idx] - median[idx])^2))
    maeMean = mean(abs(target[idx] - mean[idx]))
    maeMedian = mean(abs(target[idx] - median[idx]))
    es = energyScore(target[idx], sample[,idx, drop=FALSE])
    
    out[paste0("rmseMean",label,suffix)] = rmseMean
    out[paste0("rmseMedian",label,suffix)] = rmseMedian
    out[paste0("maeMean",label,suffix)] = maeMean
    out[paste0("maeMedian",label,suffix)] = maeMedian
    out[paste0("es",label,suffix)] = es
  }
  return(out)
}

# Forces symmetry in a covariance matrix
makeSymmetric <- function(covariance){
  
  upperIdx = upper.tri(covariance)
  lowerIdx = lower.tri(covariance)
  
  # This should not make any difference since the matrix is (or very close to) symmetric
  if (sum(covariance[upperIdx]) > sum(covariance[lowerIdx])){
    sym = forceSymmetric(covariance, "U")
  } else {
    sym = forceSymmetric(covariance, "L")
  }
  return(sym)
}

# Fast Vectorized implementation following
# "Assessing probabilistic forecasts of multivariate quantities, with an application to ensemble predictions of surface winds"
# https://www.stat.washington.edu/sites/default/files/files/reports/2008/tr537.pdf
# Double checked against the library scoringRules es_sample function and other implementations
# Instead of dividing the sample in two I take a second one as a random permutation of the first
energyScore <- function(target, sample){
  size = dim(sample)[1]
  randPermutation = sample(size)
  
  samples1 = sample
  samples2 = sample[randPermutation, ]
  
  mTarget = as.matrix(rep(1, size)) %*% as.vector(target)
  term1 = (samples1 - mTarget)
  normst1 = rowNorms(term1)
  
  term2 = (samples1 - samples2)
  normst2 = rowNorms(term2)
  
  score = mean(normst1) - 0.5 * mean(normst2)
  return(score)
}

# Pay attention
# The sammint and shrmint only make sense if data is centered around 0
# It is the same method used in MinT and that is why it is implemented here
# So the results would match exactly
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
      mTar = build.target(residuals,type="D")
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

# Function to sample from multivariate normals using functions from many packages
sampleMVN <- function(mean, sigma, sampleSize, positivity=FALSE, seed=0, fromMarginals=FALSE,
                      truncAlgorithm='gibbs') {
  
  set.seed(seed)
  
  if (isSymmetric.matrix(sigma) == FALSE){
    sigma = as.matrix(makeSymmetric(sigma))
  }
  
  if (is.positive.definite(sigma) == FALSE){
    print("Approximating positive definite")
    for (i in 12:1){
      sigma_prime = make.positive.definite(sigma, tol=10^-i)
      print(paste0("L1 norm diff approximation: ", sum(abs(sigma-sigma_prime))))
      if (is.positive.definite(sigma_prime)){
        sigma = sigma_prime
        break
      }
    }
  }
  
  if (positivity){
    if (fromMarginals == FALSE){
      sample = tmvtnorm::rtmvnorm(sampleSize, mean=as.numeric(mean), sigma=sigma,
                                  lower=rep(0, length(mean)), algorithm=truncAlgorithm)
      # sample = TruncatedNormal::rtmvnorm(n=sampleSize, mu=as.numeric(mean), sigma=sigma, lb=rep(0, length(mean)))
    } else {
      sample = do.call(cbind, foreach(i = 1:length(mean)) %dopar% {truncnorm::rtruncnorm(n=sampleSize, a=0, mean=mean[i], sd=sigma[i,i])})
    }
  } else {
    if (fromMarginals == FALSE){
      sample = MASS::mvrnorm(n=sampleSize, mu=mean, Sigma=sigma)  
    } else {
      sample = do.call(cbind, foreach(i = 1:length(mean)) %dopar% {rnorm(n=sampleSize, mean=mean[i], sd=sigma[i,i])})
    }
  }
  return(sample)
}

# Bayesian Reconciliation function
bayesRecon <- function(preds, mSumMatrix, mCovar, sampleSize=100000, reconType="pmint",
                       seed=0, kh=1){
  
  mCovar = cbind(mCovar)
  recons <- c("pmint", "lg")
  if (! (reconType %in% recons)){
    print("method types are:")
    print(recons)
    stop ("Wrong method supplied." )
  }
  
  # Defining sumMatrix and idxs
  bottomIdx <- seq( nrow(mSumMatrix) - ncol(mSumMatrix) +1, nrow(mSumMatrix))
  upperIdx <- setdiff(1:nrow(mSumMatrix),bottomIdx)
  
  # Defining covariance
  mA = mSumMatrix[upperIdx, ]
  mAtr = t(mA)
  mSumMatrixT = t(mSumMatrix)
  
  # Taking predictions
  vBottomPreds = preds[bottomIdx]
  vTopPreds = preds[upperIdx]
  
  # BottomUp Prior
  vPriorMeans = mSumMatrix %*% vBottomPreds
  vBottomPriorMean = vPriorMeans[bottomIdx]
  vIncoherence = (vTopPreds - vPriorMeans[upperIdx])
  
  # Bottom covariance
  mSigmaB = mCovar[bottomIdx, bottomIdx]
  
  # Upper covariance
  mSigmaU = mCovar[upperIdx, upperIdx]
  
  # Update rule
  # Gain is n_base by n_base
  # cppMatMult is faster to bigger matrices
  # Smaller matrices I use %*% since the overhead would not make it worth
  
  if (reconType=="lg"){
    
    # mSigma11 is mSigmaB
    mP1Gain = cppMatMult(mSigmaB, mAtr) # mSigma12 and mSigma21
    mP2Gain = (mSigmaU + cppMatMult(cppMatMult(mA, mSigmaB), mAtr)) #mSigma22
    
    # # mSigma11 is mSigmaB
    # mP1Gain = (mSigmaB %*% mAtr) # mSigma12 and mSigma21
    # mP2Gain = (mSigmaU + (mA %*% mSigmaB %*% mAtr)) #mSigma22
    
  } else if (reconType=="pmint") {
    
    # Cross covariance
    mM = -1 * mCovar[bottomIdx, upperIdx]
    mMtr = -1 * mCovar[upperIdx, bottomIdx]
    
    # Gain is n_base by + n_base
    # mSigma11 is mSigmaB
    mP1Gain = (cppMatMult(mSigmaB,mAtr) + mM) # mSigma12 and mSigma21
    mP2Gain = (mSigmaU + cppMatMult(cppMatMult(mA,mSigmaB),mAtr) + cppMatMult(mA,mM) + cppMatMult(mMtr,mAtr)) #mSigma22
    
    # mP1Gain = ((mSigmaB %*% mAtr) + mM) # mSigma12 and mSigma21
    # mP2Gain = (mSigmaU + (mA %*% mSigmaB %*% mAtr) + (mA%*%mM) + (mMtr%*%mAtr)) #mSigma22
  } else {
    stop("Wrong method in bayesian reconciliation")
  }
  
  # Adjusting for kh
  # mP1Gain is Cov(Bt+h, Ut+h | It,b)
  mP1Gain = mP1Gain * kh
  # mP2Gain is Cov(Ut+h | It,b)
  mP2Gain = mP2Gain * kh
  # mSigmaB * kh
  mSigmaB_kh = mSigmaB * kh
  
  # Gain, kh disappear in the Gain formula due to multiplying kh * (1/kh) from the inverse part
  mGain = cppMatMult(mP1Gain, cppInverse(mP2Gain + 1e-6*diag(mP2Gain)))
  
  # Prior Covariance
  mPriorCovJoint = rbind(cbind(mP2Gain, t(mP1Gain)), cbind(mP1Gain, mSigmaB_kh))
  
  # Posterior mean
  vBottomPosteriorMean = vBottomPriorMean + mGain %*% vIncoherence
  
  # Posterior Cov
  # Same as kh x (mSigmaB - mGain x t(mP1Gain))
  # in practice is (kh mSigmaB) - mGain (t(kh Cov(Bt+h, Ut+h | It,b)))
  # kh can be put out
  mSigmaBPosterior = mSigmaB_kh - cppMatMult(mGain, t(mP1Gain))
  
  # it is not symmetric due to matrix multiplication error accumulating, so we force symmetry
  # and mirror the part with highest sum of variance (the difference should be close to zero anyway)
  mSigmaBPosterior = as.matrix(makeSymmetric(mSigmaBPosterior))
  
  # Coherent predictions and full covariance matrix
  vCoherentPreds = mSumMatrix %*% vBottomPosteriorMean
  mCoherentCovariance = cppMatMult(cppMatMult(mSumMatrix, mSigmaBPosterior), mSumMatrixT)
  mCoherentCovariance = as.matrix(makeSymmetric(mCoherentCovariance))
  
  colnames(mCoherentCovariance) = colnames(mCovar)
  rownames(mCoherentCovariance) = rownames(mCovar)
  
  out = list(posteriorMean=vBottomPosteriorMean, posteriorVariance=mSigmaBPosterior,
             coherentPreds=vCoherentPreds, coherentCovariance=mCoherentCovariance,
             incoherentPreds=preds, priorCovJoint=mPriorCovJoint, priorMeanJoint=vPriorMeans,
             seed=seed, reconType=reconType)
  
  sample = sampleMVN(vBottomPosteriorMean, mSigmaBPosterior, sampleSize, seed=seed)
  
  # Mean & Median
  meanSample = colMeans(sample)
  medianSample = colMedians(sample)
  
  # Returning full bottom up sample
  out$sample = cppMatMult(sample, mSumMatrixT) #sample %*% t(mSumMatrix)
  out$posteriorMeanSample = meanSample
  out$posteriorMedianSample = medianSample
  out$coherentPredsSample = mSumMatrix %*% meanSample
  out$coherentPredsMedianSample = mSumMatrix %*% medianSample
  
  return(out)
}

loadDataset <- function(dset, synth_n=100, synthCorrel=0.5){
  
  feasibleDset <- c("infantgts", "tourism", "syntheticLarge")
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
  else if (dset=="syntheticLarge"){
    source("drawLargeHierarchy.R")
    #we generate the hierarchy with *four* bottom time series
    synthTs <- simulFourBottom(n=synth_n, noisy=TRUE, factor=2)
    y=ts(synthTs, frequency = 1)
    hierTs <- hts(y, bnames = colnames(y), characters=c(1,1))
  }
  return(hierTs)
}