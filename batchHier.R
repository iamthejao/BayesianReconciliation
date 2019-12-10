library(foreach)
library(doMC)
#registerDoMC(10)
source("hierRecBayesianExperiment.R")

batchHierAndParse <- function(dset, folder, seed=0, hs=c(1,2,3,4), models=c("ets", "arima")){
  
  for (h in hs){
    for (model in models){
      print(paste(as.character(h)," ",model))
      tic("iTest from to 50-h")
      limit = 50-h
      dsets  = foreach (iTest = (1:limit), .export=c(dset, model, h)) %do% 
        {
          ans = hierRecBayesianExperiment(get("dset"), 1, fmethod=get("model"), iTest=iTest, seed=seed)
          return(ans)
        }
      toc()
      #return(dsets)
      # Appending par results together and writing
      #dataFrame = do.call("rbind", dsets)
      #filename <- paste(folder, "/bayesian_mse_",dset,".csv",sep="")
      #writeNames <- TRUE
      #if(file.exists(filename)){
      # writeNames <- FALSE
      #}
      #write.table(dataFrame, file=filename, append = TRUE, sep=",", row.names = FALSE, col.names = writeNames)
    }
  }
}

t = batchHierAndParse("synthetic", "results/bayesian_results")