#a simple function to draw synthetic data using for the noise the covariance matrix
#of the minT paper, sec 3.2 



simulFourBottom <- function(n){
  library(MASS)
  bottomTs <- 4
  #draw pars
  phi <- runif(bottomTs)
  
  #draw noise
  Sigma <- matrix(data = rbind(c(5,3,2,1), 
                               c(3,4,2,1), 
                               c(2,2,5,3), 
                               c(1,1,3,4)), nrow = 4,ncol = 4)
  #generate an excess of noise
  noise <- mvrnorm(n=2*n, mu = c(0,0,0,0), Sigma= Sigma)
  simul <- matrix( nrow = (n+1), ncol = bottomTs)
  simul[1,] = rep(0,bottomTs)
  c <- 1
  for (t in 2:(n+1)){
    for (ts in 1:bottomTs ){
      simul[t,ts] <- c + phi[ts] * simul[t-1,ts] + noise[t,ts]
    }
  }
  colnames(simul) <- c("aa","ab","ba","bb")
  return (simul[-1,])
}

