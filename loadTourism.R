loadTourism  <- function() {
  tourism <- read.csv("TourismData_v3.csv", 
                      na = "empty")
  
  #the first two columns contain the time
  tourism <- tourism[,-1:-2]
  
  #add the zone code, constituted by the first two characters
  # of each label
  tsnames <- colnames(tourism)
  zones <- vector(length = length(tsnames))
  for (i in 1:length(tsnames)) {
    tsnames[i] <- paste(substr(tsnames[i],1,2),tsnames[i],sep="")
    zones[i] <- substr(tsnames[i],1,2)
  }


  #add the state code, constituted by the first character
  # of each label
  for (i in 1:length(tsnames)) {
    tsnames[i] <- paste(substr(tsnames[i],1,1),tsnames[i],sep="")
  }
  #sets the labels containing the full geographic information
  colnames(tourism) <- tsnames
  
  #we known that it start in Jan 1998
  y=ts(tourism, frequency = 12, start = c(1998,1))
  
  
  groupedTourism <- gts(y=y, characters = list(c(1,2, 3), 3))
  #at this point,   
  #Number of groups at each level: 1 7 27 76 4 28 108 304 
  #Total number of series: 555 
  #Number of observations per series: 228 
  
  return (groupedTourism)
}