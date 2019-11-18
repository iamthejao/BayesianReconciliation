library(foreach)
library(doMC)
registerDoMC(5)

batchHierAndParse <- function(dset, version, folder, seed=0, hs=c(1,2,3,4), models=c("ets", "arima")){
  source("parseHierResults_aggregatedH.R")
  if (version == "new"){
    source("hierRecFull.R")
  } else{
    source("hierRecOld.R")
  }
  
  for (h in hs){
    for (model in models){
      print(paste(as.character(h)," ",model))
      tic("iTest from to 50-h")
      dsets <- foreach (iTest = (1:(50-h))) %dopar% {
        print(iTest)
        if (version == "new"){
          ans = hierRecFull(dset,h=h, fmethod = model, iTest = iTest, seed=seed)
        } else  {
          ans = hierRecOld(dset,h=h, fmethod = model, iTest = iTest, seed=seed)
        }
        return(ans)
      }
      toc()
      
      # Appending par results together and writing
      dataFrame = do.call("rbind", dsets)
      filename <- paste(folder, "/mse_",dset,".csv",sep="")
      writeNames <- TRUE
      if(file.exists(filename)){
       writeNames <- FALSE
      }
      write.table(dataFrame, file=filename, append = TRUE, sep=",", row.names = FALSE, col.names = writeNames)
    }
  }
  parseHierResults(dset, folder)
}

#batchHierAndParse("infantgts", "old", "results/current_method")
#batchHierAndParse("infantgts", "new", "results/full_covariance")
