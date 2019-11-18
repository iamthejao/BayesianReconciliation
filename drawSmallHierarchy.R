#draws the autoreg and the noise correlation
drawPhiSigma <- function(correl=0.01){  
  bottomTs <- 2
  sigmaFound = FALSE
  while (sigmaFound == FALSE){
    phi <- runif(bottomTs, min=-0.98, max=0.98)
    #covariance of the noise of the two bottom time series(it is also a correlation as we assume the noise
    # to have variance 1)
    sigma_e12 <- get_rho12_exact(phi, correl)
    
    #we need to be not na and to have at least one abs root smaller than 1
    absSmaller_1 <- abs(sigma_e12) < 1
    sigmaFound  <- !is.na(sigma_e12[1]) & sum(absSmaller_1)>0
    if (sigmaFound) {
      sigma_e12 <- sigma_e12[absSmaller_1]
      #select the first root in case they are both valid
      if (length(sigma_e12)>1){
        sigma_e12 <- sigma_e12[1]
      }
    }
  }
  return(list("phi"=phi, "sigma_e12"= sigma_e12))
}


#determines the value of sigma12 in order to obtain the desired correlation between B1 and U
get_rho12_numerical <- function (phi, corB1_U) {
  
  #we assume everywhere noise with variance 1
  inv_sigmaB1 <- sqrt(1 - phi[1]^2)
  a <- 1/(1-phi[1]^2)
  b <- 1 - prod(phi)
  c <- 1/(1-phi[1]^2)
  d <- 1/(1-phi[2]^2)
  
  #range is +-1 becasue we are talking about a correlation
  inverse = function (f, lower = -1, upper = 1) {
    function (y) uniroot((function (x) f(x) - y), lower = lower, upper = upper)[1]
  }
  
  
  correl_inverse = inverse(function (x) inv_sigmaB1 * (a + x/(1-prod(phi))) * 
                             1/(sqrt(c + d + 2 * x/b)) )
  
  rho12 = correl_inverse(corB1_U)
  return(rho12)
}

get_rho12_exact <- function(phi, corB1_U){
  a <- -1 * (1 - corB1_U^2) / (1-prod(phi)) 
  b <- (1 - corB1_U^2)^2 / (1-prod(phi))^2
  c <- (1-phi[1]^2)/(1-prod(phi))^2
  d <- 1/(1-phi[1]^2) * (1-corB1_U^2) - corB1_U^2 * 1/(1-phi[2]^2)
  e <- (1 - phi[1]^2) / (1-prod(phi))^2
  
  
  rho12 <- vector(length = 2)
  
  rho12[1] <- a + sqrt (b - c * d)
  rho12[2] <- a - sqrt (b - c * d)
  
  
  rho12 <- rho12/e
  return(rho12)
} 

#draw correlated noise for n time steps, return a vector (sampleSize x 2)
#assumes the variance to be 1, the bottomCorrel is set by the user
#in this very specific case, covariance = correlation.
#by default 2 bottom time series are assumed
drawNoise <- function(sampleSize, bottomCovar, bottomTs=2){
  library(MASS)
  #generate an excess of noise observations
  #non-generic code, but might be enough for some limited simulation
  # if (bottomTs==2)
  noise <- mvrnorm(n=2*sampleSize, mu = c(0,0), Sigma= rbind(c(1,bottomCovar),c(bottomCovar,1)))
  # if (bottomTs==4)
  #   noise <- mvrnorm(n=2*sampleSize, mu = c(0,0,0,0), Sigma= rbind(c(1,bottomCovar,bottomCovar, bottomCovar),
  #                                                                  c(bottomCovar,1, bottomCovar, bottomCovar),
  #                                                                  c(bottomCovar, bottomCovar, 1, bottomCovar),
  #                                                                  c(bottomCovar, bottomCovar, bottomCovar, 1)
  #   ))
  
  return (noise)
}


simulVAR1 <- function(n, phi, noise){
  howManyTs <- length(phi)
  simul <- matrix( nrow = (n+1), ncol = howManyTs)
  simul[1,] = rep(0,howManyTs)
  c <- 1
  for (t in 2:(n+1)){
    for (ts in 1:howManyTs ){
      simul[t,ts] <- c + phi[ts] * simul[t-1,ts] + noise[t,ts]
    }
  }
  return (simul[-1,])
}


#draws phi; compute the noise correlation in order to obtain the desired correlation B1-U;
#simulate the two bottom time series with correlated noise and the upper time series as a sum
#return the hierarchy as a data frame
#n is the sample size, while correl is the required correlation between B1 and U
artificialTs <- function(n, correl){
  pars <- drawPhiSigma(correl)
  noise <- drawNoise(n, pars$sigma_e12)
  
  bottomTs <- simulVAR1(n, pars$phi, noise)
  # upperTs <- rowSums(bottomTs)
  # hierarchy <- data.frame(cbind(bottomTs, upperTs))
  #we only return the bottom time series, the upper is computed by hts
  bottomTs <- data.frame(bottomTs)
  colnames(bottomTs) <- c("b1","b2")
  
  # debug
  # print ("phi:")
  # print(phi)
  # expectedVar_B1 <- 1/(1-phi[1]^2)
  # print("expected and empirical var of B1:")
  # print (c(expectedVar_B1, var(hierarchy[,1])))
  phi <- pars$phi
  sigma_e12 <- pars$sigma_e12
  expectedCovarB1B2 <- sigma_e12 / (1 - prod(phi))
  # empiricalCovarB1B2 <- cov(hierarchy)[1,2]
  # 
  # print ("expected and empirical covar of B1 and B2:")
  # print (c(expectedCovarB1B2, empiricalCovarB1B2))
  # 
  # print ("expected and empirical correl of B1 and B2:")
  # #we assume sigma^2=1, which disappears from below
  
  expectedCorrelB1B2 = expectedCovarB1B2 * sqrt( (1-phi[1]^2) * (1-phi[2]^2) )  
  # empiricalCorrelB1B2 = cor(hierarchy)[1,2]
  # print (c(expectedCorrelB1B2, empiricalCorrelB1B2))
  # 
  # print ("expected and empirical variance of U:")
  # #we assume sigma^2=1, which disappears from below
  varBottom = 1 / (1 - phi^2)
  expectedVar_U = varBottom[1] + varBottom[2] + 2 * expectedCorrelB1B2 * 
    sqrt(prod (varBottom))
  # expectedVar_U_second = varBottom[1] + varBottom[2] + 2 * correl / 
  #   (1-prod (phi))
  # 
  # empiricalVar_U = cov(hierarchy)[3,3]
  # print (c(expectedVar_U, expectedVar_U_second, empiricalVar_U))
  # 
  # print ("expected and empirical covariance (B1,U):")
  
  
  expectedCovB1_U = 1/ (1-phi[1]^2) + sigma_e12/(1-prod(phi))
  expectedCovB2_U = 1/ (1-phi[2]^2) + sigma_e12/(1-prod(phi))
  # empiricalCovB1_U = cov(hierarchy)[1,3]
  # print (c(expectedCovB1_U, empiricalCovB1_U))
  
  # print ("empirical correl of B1 and U:")
  expectedCorB1_U =  expectedCovB1_U / sqrt  ( 1/(1-phi[1]^2) * expectedVar_U )
  expectedCorB2_U =  expectedCovB2_U / sqrt  ( 1/(1-phi[2]^2) * expectedVar_U )
  
  # hierarchy = cbind(bottomTs, rowSums(bottomTs))
  # empiricalCorB1_U = cor(hierarchy)[1,3]
  # empiricalCorB2_U = cor(hierarchy)[2,3]
  # print("actual correlation B1 U")
  # print (expectedCorB1_U)
  # print("actual correlation B2 U")
  # print (expectedCorB2_U)
  
  
  
  return (list("bottomTs"=bottomTs, "corrB2_U"=expectedCorB2_U))
}




