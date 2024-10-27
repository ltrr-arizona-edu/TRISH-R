SeasClim<-function(X,begmo,endmo,kopt) {
  # Seasonalize monthly climate data in a 3-column (year-month-value) input format
  # D Meko
  # Last revised 2022 MaY 02
  #
  # Seasonalize monthly climate data in a 3-column (year-month-value) input format
  #
  #--- INPUT
  #
  # X [matrix]: monthly data in 3 columns: year-month-value
  # begmo [numeric]: start month of season (1=Jan, 12=Dec)
  # endmo [numeric]: end month of season (1=Jan, 12=Dec).
  # kopt [numeric]: indicator of how data should be seasonalized
  #   =1: sum over months
  #   =2: average over months
  #
  #---OUTPUT
  #
  # Output: list with fields
  #   F [matrix, 2-col] year and seasonalized climate variable
  #   eflag: error flag
  #     = 0 no problem
  #     = 1 input data not consistent in having 12 months per year, with year incrementing by 1
  #
  #---NOTES
  #
  # Seasons may cross year boundary, but may not exceed 12 months in length
  
  #--- NUMBER OF MONTHS IN SEASON
  
  if (endmo>=begmo){
    nmos <- endmo-begmo+1
  } else {
    nmos <- endmo + (12-begmo)+1
  }
  
  
  #--- GET INDEX TO ROWS OF ALL END MONTH, START MONTH
  L <- X[,2]==endmo
  i2 <-which(L)
  L <- X[,2]==begmo
  i1 <-which(L)
  
  #---ADJUST LAST INDICES SO THAT START MONTH PRECEDES END MONTH
  if (i1[length(i1)]>i2[length(i2)]){
    i1 <- i1[-length(i1)]
  }
  
  #---ADJUST FIRST INDICES SO THAT START MONTH PRECEDES END MONTH
  if (i1[1]>i2[1]){
    i2 <- i2[-1]
  }
  
  #--- CHECK INDICES
  #
  # Nummber of start months must equal number of end months; intervals must all be equal to expected number
  # of months in season; year of end month must increment by 1
  L1 <- length(i1)==length(i2)
  n1 <- i2-i1+1
  L2 <- all(n1==nmos)
  yr <- X[i2,1] # years of end month
  d <-diff(yr)
  L3<-FALSE
  if (all(d==1)){
    L3<-TRUE
  }
  L <- L1 && L2 && L3
  if (!L){
    eFlag<-1
    Output <- list('tsm'=NA,'eFlag'=eFlag,'begmo'=begmo,'endmo'=endmo)
    return(Output)
  }
  rm(L1,L2,L3,L,n1,d)
  
  
  #---BUILD INDEX MATRIX TO PULL SUBSET OF MONTHLY DATA
  #
  # Index to rows of input matrix X. Goal is to use index matrix to make subset matrix of values in X and
  # then sum or average over the submatrix to get the seasonalized data
  
  I2 <- t(t(rep(1, nmos))) %*% t(i2) # row-dupe end-month indices into matrix with nmos rows
  
  # Build increment matrix
  a = (nmos-1):0
  a=-1.0*a
  A <- t(t(rep(1, length(i2)))) %*% t(a)
  A <- t(A)
  
  I3 <-I2+A   # Combine into final extraction index matris
  nrow <- dim(I3)[1]
  ncol <-dim(I3)[2]
  
  
  #---PULL SUBMATRX MONTHLY DATA, SUM OR AVERAGE OVER ROWS, BIND WITH YEAR
  
  x <- X[,3] # vector of data values
  Y <- x[I3] # vector of subsets of monthly data 
  dim(Y)=c(nrow,ncol) # mxn m is # months in season and n years
  
  # Treat as P or T
  if (kopt==1){
    y <- colSums(Y)
  }else{
    y <-colMeans(Y)
  }
  y <- colMeans(Y) # vector of seasonalized data
  
  F <-cbind(yr,y) # 2-col matrix, year and value
  eFlag<-0
  
  #--- ORGANIZE OUTPUT
  
  Output <- list('tsm'=F,'eFlag'=eFlag,'begmo'=begmo,'endmo'=endmo)
  return(Output)
}