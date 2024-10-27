TrimTsm1<- function(X,yrgo,yrsp,nNeg,nPos){
  # Trim a time series matrix with constraints for single-site reconstruction (SSR)
  # D. Meko; last revised 2021-12-28
  #
  # X [matrix]r time series matrix, year as col 1
  # yrgo [1x1]r nNeg + desired start year of output matrix. For example if yrgo=1750 and
  #   nNeg=2, desired start year is 1748. Output tsm includes only those columns with data in
  #   year yrgo-mlead. If yrgo=NA, output tsm is trimmed to begin with first year with data for
  #   all series in input X (first year of common period)
  # yrsp [1x1]r desired end year - nPos of output tsm. For example, if yrsp=1990 and nPos=2,
  #   desired end year is 1992. Output tsm is truncated to end in yrsp+nPos or in the most recent
  #   year with data for any of the series in input X. If yrsp=NA, output X ends in the year of
  #   most recent data in X after column-screening for yrgo.
  #
  # Returns list Output, with fields:
  #   X: [matrix] trimmed time series matrix of input X; numbers of rows and columns generally
  #     reduced from input, X, but otherwise the same form.
  #   ix: [matrix] one-column matrix indicating which columns of the input X are columns of
  #     the output X (disregarding the year column)
  #
  # Input matrix X assumed to have year in column 1, values for time series (e.g., tree-ring indices) in
  # remaining columns.
  
  nX <- dim(U)[2]-1 # number of time series in U
  mX <- dim(U)[1] # number of years in U
  
  yrX = as.matrix(U[,1])
  X<-as.matrix(U[,2:(nX+1)])
 
  L <- is.na(X) 
  n1<-rowSums(L) # numberic vector, number of NA in each row of X
  ifull<-which(n1==0) # rows of X with no NA at any site
  iFirstFull=min(ifull) # first row in X with no data missing
  
  iany <- which(n1<nX)
  iLastAny<-max(iany) # last row with data at any sight
    

  #--- ROW AND COLUMN TRIMMING INDICES FOR START YEAR
  
  if (is.na(yrgo)){
    i1<-iFirstFull # start with first year for which data at all sites
    jcols<-1:nX # include all original columns
  } else {
    i1 <- max((yrgo-yrX[1]+1-2),1) # need tree-ring data two years before yrgo, if possible
    xthis <- X[i1,]
    L<-!is.nan(xthis)
    jcols<-which(L) # want columns with data in that year
  }
  
  #---COL-TRIM THE MATRIX FOR FIRST-YEAR SPECS
  X<-X[,jcols]
  
  # revision 2023Feb09a
  if (length(jcols)==1)
    X <- as.matrix(X)
  end
  
  nX2 <- dim(X)[2] # number of series in col-trimmed version of X
  
  
  #---FIND LAST YEAR WITH ANY DATA IN COL-TRIMMED MATRIX
  
  L <- is.na(X) 
  n1<-rowSums(L) # numberic vector, number of NA in each row of X
  iany <- which(n1<nX2)
  iLastAny<-max(iany) # last row with data at any sight

  #--- ROW  TRIMMING INDEX FOR LAST YEAR

  if (is.na(yrsp)){
    i2<-iLastAny
  } else {
   itemp <- (yrsp+nPos)-yrX[1]+1
   i2<- min(mX,itemp)
  }

  #--- TRIM THE MATRIX
  
  X <- X[i1:i2,]
  yrX <- yrX[i1:i2]
  
  
  #--- REBUILD DATA FRAME
  
  c1 <-colnames(U)[1]
  c2<-colnames(U)[-1]
  c3 <- c(c1,c2[jcols])
  Y <-cbind(yrX,X)
  Y<-data.frame(Y)
  colnames(Y)<-c3
  
  Output<-list(X=Y,ix=jcols)
  return(Output)
}
