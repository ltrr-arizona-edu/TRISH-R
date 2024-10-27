PrewhitenChrons <- function(X,p,outputDir){
  # Convert tsm of chronologies to prewhitened matrix using AR model order p
  # D. Meko; last revised 2023-01-21
  #
  #--- INPUT
  #  
  # X: [data frame] tsm of chronologies, year as column 1. All series need not
  #   cover all years of X
  # p: scalar:  prewhiten using AR model of order p (allowable are p=1, 2, or 3)
  # outputDir: any error message generated will go to this system folder
  #
  #
  #--- OUTPUT
  #
  # Output: named list with revised data frame. Will start p years later than X
  # and data columns with AR(p)-whitened versions of original columns.
  #   Xwhitened -- data frame of AR-whitened versions of input X
  #
  #
  #---NOTES
  #
  # For each series, mean is subtracted first, AR modeling is done, and residuals are
  # shifted to have same mean as original data over the years covered by residuals
  #
  # Output data frame inherits column names of input data frame, and has row names 
  # as a character vector of years. This even though first data column is also the year.
  
  source(paste(code_dir,"emssgUNH.R",sep="")) # write error file to system, specified output folder
  
  L <- !is.data.frame(X) | !any(p==c(1,2,3))
  if (L){
    emssg <- 'PrewhitenChrons: X must be data frame, and AR order p must be one of {1,2,3}'
    ResTemp<-emssgUNH(emssg,outputDir)
    stop(emssg)
  }
  
  # ALLOCATING
  # Allocate Y same row-size as X, but without year column
  nY = dim(X)[2]-1 # number of time series in X
  mY = dim(X)[1] # number of years in matrix X
  Y <- matrix(NA,nrow=mY,ncol=nY)
  
  k <- 1:mY   # nominal year of X and Y, equivalent to row number
  
  # Loop over time series in X, fitting AR model, getting residuals, 
  # restoring mean; place in column of Y
  for (n in 1:(nY)){
    j1 <- n+1
    x <- X[,j1]
    L <- complete.cases(x)
    ix <- k[L] # in X, these would be rows with data
    ilast<- ix[length(ix)] # last value of whitened series
    # should be in this row of X and Y
    x <- x[L] #  the time series, w/o flanking NA

    # SET UP REGRESSION
    #
    # Put predictand in v, predictor(s) on W 
    mx <-length(x)
    j=(p+1):mx
    v <- x[j] # predictand for AR modeling
    mv <- length(v)
    J = matrix(rep(j,p), ncol = p)
    mJ <- dim(J)[1]
    
    # Build matrix to subtract
    j1 = 1:p
    J1 <- matrix(j1, nrow=mJ, ncol=length(j1), byrow=TRUE) # row-dupe to mJ rows
    
    # Subtract indices
    J2 <- J-J1
    # each col of J2 is a vector index into x
    
    W = matrix(NA,nrow =mv,ncol=p)
    for (kk in 1:p){
      jthis <- J2[,kk]
      w = x[jthis]
      W[,kk]=w
    }
    
    # Regress
    M <- lm(v ~ W)
    y <- M$residuals+mean(v) # the prewhitened chronology
    
    # Re-insert
    i2 <- ilast; i1 <- ilast - length(y)+1 # start and end row
    Y[i1:i2,n]=y
  }
  yrY <- X[,1]
  
  # Trim off all-NA rows
  N1 <-rowSums(!is.na(Y)) # vector with number of non-NA in each row of Y
  L <- N1>0
  Y <- Y[L,]; yrY <- yrY[L]
 
   # Add year col, give col names, make as data frame
  Y <- cbind(yrY,Y)
  rownames(Y) <- as.character(yrY)
  colnames(Y) <- names(X)
  Y <- as.data.frame(Y)

  #   STORE OUTPUT
  Output <- list("Xwhitened"=Y)
  return(Output)
}

