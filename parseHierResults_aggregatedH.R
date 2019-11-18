parseHierResults <- function (dset, folder){
  prevWd = getwd()
  setwd(paste(prevWd, folder, sep="/"))
  #parse the results of hierarchical non-temporal reconciliation
  #readt the mse, extract the proportion of favorable signs and the produces the boxplot
  library(readr)
  #source('bayesianSignedRank.R')
  results <- read_csv(paste("mse_",dset,".csv",sep=""))
  results <- unique(results) #because some experiements on the cluster are duplicated
  fmethods <- unique(results$fmethod)
  horizons <- unique(results$h)
  #we realize that we are only interested in h up to 4
  horizons <- horizons[horizons<5]
  
  configs <- length(fmethods) * length(horizons) 
  
  #we need first to instantiate the data frame with placeholder values, and then we fill the correct values
  comparison <- data.frame(
    cases = rep(fmethods[1],configs),
    h = rep(1,configs),
    fmethod=rep(fmethods[1],configs),
    # medianBaseMint=rep(-1,configs),
    medianBaseBayesShr=rep(-1,configs),
    # medianBaseBayesGlasso=rep(-1,configs),
    medianMintBayesShr =rep(-1,configs),
    # medianMintBayesGlasso =rep(-1,configs),
    # pValMedianMintBayesShr=rep(-1,configs),
    # pValMedianMintBayesGlasso=rep(-1,configs),
    stringsAsFactors = FALSE
  )
  
  aggrComparison <- comparison[1:length(fmethods),]
  
  
  #analysis for each h
  counter <- 1
  for (fmethod in fmethods){
    for (h in horizons){
      print(paste(fmethod,dset))
      comparison$fmethod[counter] <- fmethod
      idx = results$fmethod==fmethod & results$h==h 
      if (sum(idx)>0){
        subresults <- results[idx,]
        comparison$cases[counter] <- sum(idx)
        comparison$fmethod[counter] <- fmethod
        comparison$h[counter] <- h
        # comparison$medianBaseMint[counter] <- median(subresults$mseBase / subresults$mseCombMintShr)
        comparison$medianBaseBayesShr[counter] <- round ( median(subresults$mseBase / subresults$mseBayesShr), digits = 2)
        # comparison$medianBaseBayesGlasso[counter] <- median(subresults$mseBase / subresults$mseBayesGlasso)
        comparison$medianMintBayesShr[counter] <- round ( median(subresults$mseCombMintShr / subresults$mseBayesShr), digits = 2)
        # comparison$medianMintBayesGlasso[counter] <- median(subresults$mseCombMintShr / subresults$mseBayesGlasso)
        # comparison$pValMedianMintBayesShr[counter] <- wilcox.test(log(subresults$mseCombMintShr/ subresults$mseBayesShr),
                                                                  # alternative="less")$p.value
        # comparison$pValMedianMintBayesGlasso[counter] <- wilcox.test(log(subresults$mseCombMintShr / subresults$mseBayesGlasso), 
                                                                     # alternative="less")$p.value
      }
      counter <- counter + 1
    }
    
    #generate the bplot with ggplot2  - eventually commented out
    library(ggplot2)
    pdfname <- paste("plot","_",dset,"_",fmethod,".pdf",sep = "")
    denom <- subresults$mseBase 
    resLenght <- length(subresults$mseBase)
    
    
    
    relMse <- rbind(matrix(subresults$mseCombMintShr/denom), matrix(subresults$mseBayesShr/denom), matrix(subresults$mseBayesGlasso/denom))
    label <-  factor(rbind(matrix(rep("MinT",resLenght)),matrix(rep("Bayes-shr",resLenght)),
                     matrix(rep("Bayes-glasso",resLenght))),
                     levels = c("MinT","Bayes-shr","Bayes-glasso"))
    
    dataPlot <- as.data.frame(relMse)
    dataPlot$label <- label
    currentPlot <- ggplot(dataPlot, aes(x = label, y = log10(relMse))) + geom_boxplot()  +
      stat_boxplot(geom = "errorbar", width = 0.5) +  #draw the whiskers
      scale_x_discrete(name = "") +
      scale_y_continuous(name = "Log (relative mse)") + 
      ggtitle(paste (dset, fmethod))
    
    scaling <- 1.8 #to avoid large outliers that make the boxplot unreadable
    if (dset=="tourism"){
      scaling<- 1.1  
    }
    else if (fmethod=="ets"){
      scaling<- 3 
    }
    
    
    ylim1 = boxplot.stats(log(dataPlot$V1))$stats[c(1, 5)]
    currentPlot = currentPlot + coord_cartesian(ylim = ylim1*scaling)  + geom_hline(yintercept = 0, color='darkblue', linetype="dashed")
    print(currentPlot)
    ggsave(pdfname, width = 4, height = 3)
}
filename=paste("summaryEachH_",dset,".csv",sep="")
write.table(comparison,file=filename,sep=",",row.names = FALSE)

#analysis aggregated over h
counter <- 1
for (fmethod in fmethods){
    aggrComparison$fmethod[counter] <- fmethod
    idx = results$fmethod==fmethod 
    if (sum(idx)>0){
      subresults <- results[idx,]
      aggrComparison$cases[counter] <- sum(idx)
      aggrComparison$fmethod[counter] <- fmethod
      aggrComparison$medianBaseMint[counter] <- median(subresults$mseBase / subresults$mseCombMintShr)
      aggrComparison$medianBaseBayesShr[counter] <- median(subresults$mseBase / subresults$mseBayesShr)
      # aggrComparison$medianBaseBayesGlasso[counter] <- median(subresults$mseBase / subresults$mseBayesGlasso)
      aggrComparison$medianMintBayesShr[counter] <- median(subresults$mseCombMintShr / subresults$mseBayesShr)
      # aggrComparison$medianMintBayesGlasso[counter] <- median(subresults$mseCombMintShr / subresults$mseBayesGlasso)
      # aggrComparison$pValMedianMintBayesShr[counter] <- wilcox.test(log(subresults$mseCombMintShr/ subresults$mseBayesShr),alternative="less")$p.value
      # aggrComparison$pValMedianMintBayesGlasso[counter] <- wilcox.test(log(subresults$mseCombMintShr / subresults$mseBayesGlasso),alternative="less")$p.value
    }
    counter <- counter + 1
  }

filename=paste("summary_",dset,".csv",sep="")
write.table(aggrComparison,file=filename,sep=",",row.names = FALSE)
setwd(prevWd)
}
