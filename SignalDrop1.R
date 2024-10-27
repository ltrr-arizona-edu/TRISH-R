SignalDrop1 <- function(N,r){
  # Drop in  maximum climate signal in recent years as chronologies drop out 
  # D. Meko; last revised 2022-09-20
  #
  # Have set of single-site reconstructions, each with a statistic (e.g., R-squared) and
  # and ending year that may differ. Want to know how that signal strength drops year by
  # year over the interval the chronologies end.
  #
  #--- IN
  # N: vector of ending years
  # r: vector of the statistics
  #
  #---OUT
  #
  # Output: named list:
  #   yrx1, x1:  vectors of years and number of available series (from N)
  #   yrx2, x1:  vectors of years and maximum statistic for those years

  yr1 <- min(N); yr2 <- max(N); yrx <- yr1:yr2 # yrx is vector of years covering ending years in N
  yrx1 <- yrx; yrx2 <- yrx # year vectors for the output
  mx <- length(yrx) # length of output vectors yrx1, yrx2, x1, x2
  
  A <- t(replicate(mx,N)) # replicate N into matrix with identical rows like N
  B <- replicate(length(N),yrx) # replicate year vector covering x1, x2 to identical columns
  
  L <- A >= B # logical matrix true of ending year in A at least as high as corresp. year in B
  x1 <- rowSums(L) # vector of the total number of series in each year

  R  <- t(replicate(mx,r)) # row-replicate the vector of statistics to mx rows
  R[!L] <- NA   # convert to NA the elements representing series not present
  x2 <- apply(R,1,max,na.rm=TRUE) # maximum statistic in each year (represents maximum for the 
  # available series)
  Output <- list(yrx1=yrx1,x1=x1,yrx2=yrx2,x2=x2)
  return(Output)
}