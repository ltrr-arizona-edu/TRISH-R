firstLast <- function(X,yrX){
  # First and last year of data in columns of a time series matrix
  # Dave Meko, February 2025
  #
  #---INPUT
  #
  # X [matrix] one or more time series in matrix; no internal NA allowed, but
  #   can have different start and end years
  # yrX [vector] years of X
  #
  #
  #---OUTPUT
  # 
  # Output: a two column matrix of first and last years
  #
  #---REFERENCES
  #
  #---REVISIONS
  #
  # revised 2025-02-01:  first working version
  #
  #---NOTES
  #
  # Built on Meko's firstlst.m Matlab function
  # Stops with error if X not a matrix, yrX not a vector, column dimension of 
  # X not equal to length of yrX, or yrX does not increment by 1
  
  
  #--- SOURCE FUNCTIONS
  
  source(paste(code_dir,"emssgUNH.R",sep="")) # write error file to system, specified output folder
 
  
  #---CHECK INPUTS
  
  # X and yrX must be compatible in row-size; yrX must increment by 1
  if (!is.matrix(X)) stop('X must be a matrix')
  if (!is.vector(yrX)) stop('yrX must be a vector')
  L1 <- dim(X)[1] == length(yrX)
  if (!L1)  stop('Row size of X not equal to length of yrX')
  d = diff(yrX)
  if (!all(d == 1)) stop('yrX does not increment by 1')
  
  
  #--- DUPICATE YEAR VECTOR INTO MATRIX SAME COL-SIZE AS X
  
  ncols <- dim(X)[2]
  V <- as.matrix(replicate(ncols,yrX)) # year matrix
  
  
  #--- COMPUTE FIRST AND LAST YEARS
  
  L1 <- is.na(X) # flag missing values in time series matrix
  V[L1] <- NA   # replace year-matrix values with NA if no data
  
  # Tmin and Tmax will be numeric variables
  Tmin <- apply(V,2,min, na.rm = TRUE) 
  Tmax <- apply(V,2,max, na.rm = TRUE) 
  T <- cbind(Tmin,Tmax) # matrix, start year in col1, end year in col 2
  rownames(T) <- colnames(X)
  colnames(T) <- c('YearGo','YearStop')
  
  #-- OUTPUT

  Output <- T
  return (Output)  
}


