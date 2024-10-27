LeaveOut<-function(nNeg, nPos, mA) {
  # Build pointer matrix for use in leave-m-out cross-validation
  # D Meko
  # Last revised 2021 May 25
  #
  # Returns logical pointer matrix whose jth column points to the predictor cols
  # for calibration of a model to supply the jth observation
  # 
  # nNeg [numberic]: maximum negative lag to be considered on tree rings
  # nPos [numeric]: maximum positive lag ...
  # mA [numeric] # number of of observations in calibration period
  #
  # Note that nNeg should be a positive number (e.g., 2 indicates 
  # lage upt to -2)
  # 
  # Leave-m-out cross-validation ensures that if lags are in the regression model
  # no tree-ring observations used to predict the "left-out" observation are also
  # used in calibration of the model that supplied the prediction (Meko 1997)
  
  m <- 1+ 4*max(nNeg,nPos) # leave this many out
  mhalf <- (m-1)/2
  i1 <- mhalf+1
  i2 <- mA-mhalf
  
  A<-matrix(1,mA,mA)
  
  for (j in 1:mA){
    a <- rep(1,mA)
    if (j<i1){
      igo<-1
      isp<-j+mhalf
      a[igo:isp]<-0
      A[,j]<-a
    } else if (j>i2) {
      igo<-j-mhalf
      isp<-mA
      a[igo:isp]<-0
      A[,j]<-a
    } else {
      igo<-(j-mhalf)
      isp<-(j+mhalf)
      a[igo:isp]<-0
      A[,j]<-a
    }
  }
  Lin<-A
  Leave <- list("Lin"=Lin,"NumberLeftOut"=m)
  return(Leave)
}