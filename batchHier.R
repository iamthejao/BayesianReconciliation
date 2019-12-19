library(foreach)
library(doMC)
registerDoMC(10)
source("hierRecBayesianExperiment.R")

batchHierAndParse <- function(dset, folder, file, seed=0, hs=c(1,2,3,4), models=c("ets", "arima")){
  
  for (model in models){
    for (h in hs){
      print(paste(as.character(h)," ",model))
      tic("iTest from to 50-h")
      limit = 50-h
      dsets  = foreach (iTest = (1:limit)) %dopar% 
        {
          ans = hierRecBayesianExperiment(dset, h, fmethod=model, iTest=iTest, seed=seed)
          return(ans)
        }
      toc()
      gc()
      # Appending par results together and writing
      dataFrame = do.call("rbind", dsets)
      filename <- paste(folder, "/",file,dset,".csv",sep="")
      writeNames <- TRUE
      if(file.exists(filename)){
      writeNames <- FALSE
      }
      write.table(dataFrame, file=filename, append = TRUE, sep=",", row.names = FALSE, col.names = writeNames)
    }
  }
}

batchHierAndParse("infantgts", "results/bayesian_results", "arimaOnly", hs=c(1,2,3,4), models=c("ets"))



