PeriodCommon <- function(X, Y) {
  # Common period of a time series with time series matrix
  # D Meko, Last revised 2021-05-17
  #
  # X [matrix] multiple time series; time vector as column 1
  # Y [matrix] single time series; time vector as column 1
  #
  # Returns list with parts
  #   X, Y [matrix] like input X, Y, but for years in which no data are
  #     missing in Y or in any of the series in X
  #   Beware that time (e.g., year) is first column of these matrices
  #   tgo, tsp [numeric]: start and ending times (e.g., years) of X and Y
  #
  # Why? Written to facilitate organization of predictors and predicand for
  # calls to regression functions.
  X<-as.matrix(X)
  Y<-as.matrix(Y)
  
  # Y ; Separate time column; 
  y1 <- Y[,-1]
  yry1 <-Y[,1]
  yrX1<-X[,1]
  X1<-X[,-1]
  
  # In case of X having just one series
  X1<-as.matrix(X1)
  yrX1<-as.matrix(yrX1)
  
  #--- BUILD NA MTX TO HOLD ALL OF AND X AND Y
  nA<-ncol(X1)+1
  ton<-min(yry1[1],yrX1[1])
  toff<-max(yry1[length(yry1)],yrX1[dim(X1)[1]])
  mA<-toff-ton+1
  
  A<-matrix(NA,mA,nA)
  yrA<-(ton:toff)
  
  
  #--- FILL COLS OF A, PREDICTAND IN COL 1 if applicable
  
  irow<-yry1-ton+1  # row indices of tarter slots in A for y
  A[irow,1]=y1;
  irow<-yrX1-ton+1  # row indices of tarter slots in A for X1
  A[irow,2:nA]=X1
  
  
  #-- FIND ROWS WITH NO NA'S AND PULL SEGMENT
       
  L=complete.cases(A)
  A<-A[L,]
  yrA<-yrA[L]
  mA<-nrow(A)
  nA<-ncol(A)
  
  # error message if time vector does not increment by 1
  L<-all(diff(yrA)==1)
  if(!L) stop('year column of A does not increment by 1')
  
  Ynew <- cbind(yrA,A[,1])
  Xnew <- cbind(yrA,A[,2:nA])
  tgo<-yrA[1]
  tsp<-yrA[mA]
  Output<-list(X=Xnew,Y=Ynew,tgo=tgo,tsp=tsp)
}


