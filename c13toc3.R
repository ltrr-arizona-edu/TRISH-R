c13toc3 <- function(V) {
  # 13-column monthly climate matrix to three-column matrix (year, month value)
  # D. Meko
  # Last revised 2021-05-30
  #
  #--- IN
  #
  # V [matrix]  13 columns, year and Jan-Dec data
  #
  #--- OUT
  #
  # X [matrix]  3 columns, year, month day
  
  mV <- nrow(V) # number of rows (years) in the 13-col matrix
  
  yrV<-V[,1]
  V<-V[,-1]
  
  # Reshaped data
  v <- as.vector(t(V))
  x<-as.matrix(v) # data as 1-col matrix
  
  # Build year 1-col matrix
  yr<-yrV[1]:yrV[mV]
  yr <- t(replicate(12,yr))
  yr <- as.vector(yr)
  yr <- as.matrix(yr) # year as 1col matrix
  
  # Build month matrix
  u <- 1:12
  U <- replicate(mV,u)
  u <- as.vector(U)
  month <- as.matrix(u)
  #u <- as.matrix(u) # 1 col matrix
  #u <-  replicate(u,mV)
  
  X <-cbind(yr,month,x)
  
  return(X)
}
