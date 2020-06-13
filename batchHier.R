list.of.packages <- c("foreach", "optparse")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(optparse)

# Here it starts the main
parser = OptionParser()
parser <- add_option(parser, c("-d", "--dataset"), type="character", default=NULL,
                     help="Insert dataset name (infantgts, tourism, synthetic, syntheticLarge)",
                     metavar="string")

parser <- add_option(parser, c("-s", "--hstep"), type="integer", default=NULL,
                     help="Insert h steps ahead (1,2,3,4)",
                     metavar="numeric")

parser <- add_option(parser, c("-m", "--method"), type="character", default=NULL,
                     help="Insert method ARIMA or ETS",
                     metavar="string")

parser <- add_option(parser, c("-i", "--init"), type="integer", default=1,
                     help="Test range initial value",
                     metavar="numeric")

parser <- add_option(parser, c("-e", "--end"), type="integer", default=50,
                     help="Test range end value",
                     metavar="numeric")

parser <- add_option(parser, c("-l", "--size"), type="integer", default=100,
                     help="Synthetic large number of timesteps",
                     metavar="numeric")

parser <- add_option(parser, c("-c", "--corr"), type="double", default=0.5,
                     help="Synthetic small timeseries correlation",
                     metavar="numeric")

parser <- add_option(parser, c("-p", "--path"), type="character", default="results/bayesian_results",
                     help="Output path",
                     metavar="string")

parser <- add_option(parser, c("-k", "--khone"), type="logical", default=TRUE,
                     help="Enforce kh = 1 or kh = h",
                     metavar="numeric")

args = commandArgs(trailingOnly=TRUE)
parsed = parse_args(parser, args=args)

print(parsed)

source("hierRecBayesianExperiment.R")
library(foreach)

batchHierAndParse <- function(dset, h, model, iTestFrom=1, iTestTo=50, seed=0, folder="results/bayesian_results_kh_basekh",
                              synth_n=100, synthCorrel=0.5, enforceKhOne=TRUE){
  
  # Guaranteeing arguments are parsed correctly
  limit = 50 - h + 1
  iTestTo = min(limit, iTestTo)
  if (iTestFrom >= iTestTo){
    return(NULL)
  }
  
  if (dset=="synthetic"){
    dset_name <- paste0(dset,"_correl",synthCorrel,"_seed",seed)
  } else if (dset=="syntheticLarge") {
    dset_name <- paste0("largeSynthetic_n",synth_n,"_seed",seed)
  } else {
    dset_name = dset
  }
  
  filename <- paste(folder, "/",h,"_", dset_name,"_",model,"_", iTestFrom,"to", iTestTo,".csv",sep="")
  print(paste("Filename: ", filename))
  
  # Save table at every run
  for (iTest in (iTestFrom:iTestTo)) {
    
    print(paste0("iTest: ", iTest))
    
    # Probabiliy of calling MinT function to actually checking if the results match (and save)
    if (dset != 'tourism'){
      p = 0.05
    } else {
      p = 0.01
    }
    
    print(paste0("Prob of checking MinT: ",p))
    
    # Calculate
    dataFrame = hierRecBayesianExperiment(dset, h, fmethod=model, iTest=iTest, seed=seed, testProbability=p,
                                          savePredictions=TRUE, synth_n=synth_n, synthCorrel=synthCorrel,
                                          saveSamples=FALSE, enforceKhOne=enforceKhOne)
    # Collect garbage
    gc()
    
    # Write
    writeNames <- TRUE
    if(file.exists(filename)){
      writeNames <- FALSE
    }
    write.table(dataFrame, file=filename, append = TRUE, sep=",", row.names = FALSE, col.names = writeNames)
  }
  return(NULL)
}

tic("Batch Run")
if (is.null(parsed$hstep)){
  for(h in 1:4){
    batchHierAndParse(parsed$dataset, h, parsed$method,
                      parsed$init, parsed$end,
                      synth_n=parsed$size, synthCorrel=parsed$corr,
                      folder=parsed$path, enforceKhOne=parsed$khone, seed=123321)
  }
}else{
  batchHierAndParse(parsed$dataset, parsed$hstep, parsed$method,
                    parsed$init, parsed$end,
                    synth_n=parsed$size, synthCorrel=parsed$corr,
                    folder=parsed$path, enforceKhOne=parsed$khone)
}