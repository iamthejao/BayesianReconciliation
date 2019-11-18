library(hts)
library(huge) #covariance matrix via glasso
library(SHIP) #shrinkage of covarianca matrix
library(car)
library(matlib)
source("drawSmallHierarchy.R")
source("drawLargeHierarchy.R")
source("reconciliationMethods.R")

# Settings
h = 1
iTest = 1
testSize <- 50

#if h=1, the possible preds are the whole test size lenght; 
#if h=2, the possible preds are the (test size lenght -1); etc.
possiblePreds <- testSize - h + 1

# Drawing ArtificialTs

# Those are the bottom time series
# for the big synth we have 4 on bottom that builds
# the format 1 2 4 (Total, (A,B), (Aa, Ab, Ba, Bb))

synth_n=100
#synthCorrel=0.5
#listSynth <- artificialTs(n=synth_n + h, correl = synthCorrel)
#synthTs <- listSynth$bottomTs
#corrB2_U <- listSynth$corrB2_U
synthTs <- simulFourBottom(n=synth_n)

# Making timeseries object
y=ts(synthTs, frequency = 1)
# Making hierarchical timeseries Object with 1 level (sum of basis)
hierTs <- hts(y, bnames = colnames(y), characters = c(1,1))
plot.gts(hierTs)

# Train and test
timeIdx             <- time(hierTs$bts[,1])
startTrain          <- timeIdx[1]
endTrain            <- length(timeIdx) - h - (iTest - 1)
train               <- window(hierTs, start = timeIdx[startTrain], end = timeIdx[endTrain] )
test                <- window(hierTs, start =timeIdx[endTrain +1], end=timeIdx[endTrain + h])

# Recompute predictions to be  accessed by the Bayesian method
allTsTrain <- allts(train)
numTs <- ncol(allTsTrain)
alpha <- 0.2

# Making vector and matrixes for storaging values
sigma <- vector(length = numTs)
preds <- vector(length = numTs)

# Residuals from model
residuals <- matrix(nrow=dim(allTsTrain)[1], ncol = numTs)
colnames(residuals) <- colnames(allTsTrain)

# Fit from model
fitted <- matrix(nrow=dim(allTsTrain)[1], ncol = numTs)
colnames(fitted) <- colnames(allTsTrain)

# Same as allTsTrain
actual <- matrix(nrow=dim(allTsTrain)[1], ncol = numTs)
colnames(actual) <- colnames(allTsTrain)

# Fit models on training data
fmethod = "arima" #ets
for (i in 1:numTs){
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

# Plotting basis and fitted basis
par(mar=c(1,1,1,1))
par(mfrow=c(4,1))
for (i in colnames(y)){
  plot.default(fitted[,i], main=i, col="red", ylim=c(min(y[,i]*0.95), max(y[,i]*1.05)))
  lines(fitted[,i])
  points(y[,i], main=i, pch=12)
}

mSumMatrix = smatrix(train)
mseBase =  mean  ( (allts(test)[h,] - preds)^2 )

bottomIdx <- seq( nrow(mSumMatrix) - ncol(mSumMatrix) +1, nrow(mSumMatrix))
predsIncoherent = mSumMatrix %*% preds[bottomIdx]
mseIncoherent =  mean  ( (allts(test)[h,] - predsIncoherent)^2 )

mseBayesDiag =  mean  ( (allts(test)[h,] - bayesReconcFull(preds, residuals, mSumMatrix, "diagonal", sigmas=sigma))^2 )
mseBayesSam =  mean  ( (allts(test)[h,] - bayesReconcFull(preds, residuals, mSumMatrix, "sam"))^2 )
mseBayesGlasso =  mean  ( (allts(test)[h,] - bayesReconcFull(preds, residuals, mSumMatrix, "glasso"))^2 )
mseBayesShr =  mean  ( (allts(test)[h,] - bayesReconcFull(preds, residuals, mSumMatrix, "shr"))^2 )

print(c("MseBase: ",mseBase))
print(c("MseIncoherent: ",mseIncoherent))
print(c("MseCoherentDiagonal: ",mseBayesDiag))
print(c("MseCoherentSam: ",mseBayesDiag))
print(c("MseCoherentGlasso: ",mseBayesGlasso))
print(c("MseCoherentShr: ",mseBayesShr))

if (dset=="synthetic"){
  dataFrame <- data.frame(h, fmethod, synth_n, synthCorrel, corrB2_U, mseBase,mseCombMintSample,
                          mseCombMintShr, mseBayesDiag, mseBayesSample, mseBayesGlasso, mseBayesShr)
  colnames(dataFrame) <- c("h","fmethod","sampleSize","correlB1_U","correlB2_U",
                           "mseBase","mseMintSample","mseCombMintShr","mseBayesDiag","mseBayesSample",
                           "mseBayesGlasso", "mseBayesShr")
  dset <- paste0(dset,"_correl",synthCorrel,"_n",synth_n)
}

else if (dset=="syntheticLarge"){
  dataFrame <- data.frame(h, fmethod, synth_n, mseBase,mseCombMintSample,
                          mseCombMintShr, mseBayesDiag, mseBayesSample, mseBayesGlasso, mseBayesShr)
  colnames(dataFrame) <- c("h","fmethod","sampleSize",
                           "mseBase","mseMintSample","mseCombMintShr","mseBayesDiag","mseBayesSample",
                           "mseBayesGlasso","mseBayesShr")
  dset <- paste0("largeSynthetic_n",synth_n)
}

else
{
  dataFrame <- data.frame(h, fmethod, dset, calibration50, calibration80, 
                          mseBase,mseCombMintSample,mseCombMintShr,mseBayesDiag,mseBayesSample,mseBayesGlasso,mseBayesShr)
}


filename <- paste("results/mse_",dset,".csv",sep="")

writeNames <- TRUE
if(file.exists(filename)){
  writeNames <- FALSE
}
write.table(dataFrame, file=filename, append = TRUE, sep=",", row.names = FALSE, col.names = writeNames)






  