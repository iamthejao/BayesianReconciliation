#a simple function to draw synthetic data using for the noise the covariance matrix
#of the minT paper, sec 3.2 

simulFourBottom <- function(n, diag=FALSE, noisy=FALSE, factor=1){
  library(MASS)
  bottomTs <- 4
  #draw pars
  phi <- runif(bottomTs)
  
  #draw noise
  if (diag==TRUE){
    Sigma <- matrix(data = rbind(c(5,0,0,0), 
                                 c(0,4,0,0), 
                                 c(0,0,5,0), 
                                 c(0,0,0,4)), nrow = 4,ncol = 4)
  } else {
    Sigma <- matrix(data = rbind(c(5,3,2,1), 
                                 c(3,4,2,1), 
                                 c(2,2,5,3), 
                                 c(1,1,3,4)), nrow = 4,ncol = 4)
  }
  #generate an excess of noise
  noise <- mvrnorm(n=2*n, mu = c(0,0,0,0), Sigma= Sigma)
  epsNoise <- rnorm(n=2*n, mean = 0, sd=sqrt(5*factor))
  
  simul <- matrix( nrow = (n+1), ncol = bottomTs)
  simul[1,] = rep(0,bottomTs)
  c <- 1
  for (t in 2:(n+1)){
    for (ts in 1:bottomTs ){
      
      if (noisy){
        if ((ts %% 2) == 0){
          sign = 1
        } else {
          sign = -1
        }
      } else {
        sign = 0
      }
      simul[t,ts] <- c + phi[ts] * simul[t-1,ts] + noise[t,ts] + (sign * epsNoise[t])
    }
  }
  colnames(simul) <- c("aa","ab","ba","bb")
  return (simul[-1,])
}