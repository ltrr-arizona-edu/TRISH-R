ssValid<- function(y,X,ical,ival,i1) {
  # Split-sample calibration-validation of regression model 
  # D Meko
  # last revised 2024-03-08
  #
  # Does one half of the split sample validation/calibration. Generally would be called
  # twice, first time with ical and ival pointing to first and second halves of data, and then
  # to second and first halves of data. Written because needed by function resonsw4.
  #
  # y [matrix] single-col matrix of predictand
  # X [matrix]  predictors, not all of which may be in model model
  # ical [vector]i  vector of rows for calibration part of y,X
  # ival [vector]i vector of rows of validation...
  # i1 [vector]i  columns of X that are to be used as predictors in regression
  #
  # Returns named list Output with fields:
  #   RE reduction of error statistic
  #   PearsonRcalib: correlation of predicted values with observed for cal period
  #   PearsonRvalid: correlation of predicted values with observed for val period
  #   MeanObsCalibPd:  mean obs predictand for calib period
  #   MeanObsValidPd:  mean obs predictand for validation period
  #   MeanRecCalibPd: mean recon predctand for calibration period
  #   MeanRecValidPd:  mean recon predctand for validation period
  #   RsquaredCalib: Calibration R squared of regression
  #   RsquaredValid: Prediction R squared; as computed here this is equivalent to 
  #     the reduction of error statistic
  #   MeanAbsErrorCalibPd:  mean absolute error of reconstruction for calib period
  #   MeanAbsErrorValidPd:  mean absolute error of reconstruction for validation period  
  #   SampleSizeCalib:number of observations in calib period
  #   SampleSizeValid:number of observations in validation period
  #   SummaryStatisticsMatrix (includes all of the above, with value for calib
  #       period in first col and for validation period in second col
  #   ssPred: vector of predictions for validation period
  #
  # Output.A adds nothing new to the other fields of Outlook, but may be useful for a quick look 
  # at statistics of calibration and validation (cols 1 and 2 of A). 
  #
  # revised 2024-03-08: cosmetic. Correction of typos in comments
  
  #--- ALLOCATE

  A<-matrix(NA,nrow=6,ncol=2) # to hold for calib (col1) and validation (col2):
  # Pearson r observed with predicted
  # Observed mean
  # Reconstructed mean
  # R squared
  # Mean absolute error
  # Nunber of observaions
  
  #--- CALIBRATION
  
  y <- as.matrix(y) # in case y happened to be passes as vector
  yc<-as.matrix((y[ical,1]))
  X<-as.matrix(X)
  Xc<-as.matrix(X[ical,i1]);

  G<-lm(yc~Xc)
  
  A[1,1]<-cor(G$fitted.values,as.vector(y[ical])) # Pearson r obs with recon
  A[2,1]<-mean(y[ical]) # mean observed y
  A[3,1]<-mean(G$fitted.values) # mean recon y
  A[4,1]<-summary(G)$r.squared # R squared of calibration
  A[5,1]<-mean(abs(G$residuals)) # mean absolute error
  A[6,1]<-length(ical)
  
  PredNull<-rep(A[2,1],length(ival))# The "null" prediction, defined as the observed
  # calib-period mean, expanded to a vector the length of the calibration period
  
  #--- VALIDATION
  
  mValid<-length(ival)# number of obs in validation period
  eNull <-  y[ival]-PredNull # null-prediction  errors, computed as observed y minus 
  # prediction (vector)
  
  # Use calib-pd model to generate prediction for validation period. 
  X1 <- cbind((matrix(1,nrow=mValid,ncol=1)),X[ival,i1]) # predictor matrix, with leading col of ones
  b<-as.matrix(G$coefficients) # matrix, 1 col
  yhat<-X1 %*% b # predictions for validation period (1-col matrix)
  ev <-  y[ival]-yhat # reconstruction errors, computed as observed y minus predicted
  
  # Store some results
  A[1,2]<- cor(as.matrix(y[ival]),yhat) # Pearson r calibration-period prediction with observed
  A[2,2]<- mean(y[ival]) # mean observed
  A[3,2]<- mean(yhat) # mean reconstructed
  
  # Compute RE and the prediction R squared, which as I am computing them are identical
  SSEv<-sum(ev*ev) # sum of squares of validation errors
  SSEnull  <- sum(eNull*eNull) # sum of squares of null-prediction residuals
  RE<-1-(SSEv/SSEnull) # reduction of error statistic
  A[4,2]<- RE # prediction R squared
  A[5,2]<-mean(abs(ev))
  A[6,2]<-length(ival)
  

  #--- RETURN OUTPUT AS LIST
  
  Output <- list("RE"=RE,"PearsonRcalib"=A[1,1],"PearsonRvalid"=A[1,2],"MeanObsCalibPd"=A[2,1],
  "MeanObsValidPd"=A[2,2],"MeanRecCalibPd"=A[3,1],"MeanRecValidPd"=A[3,2],
  "RsquaredCalib"=A[4,1],"RsquaredValid"=A[4,2],"MeanAbsErrorCalibPd"=A[5,1],
  "MeanAbsErrorValidPd"=A[5,2],"SampleSizeCalib"=A[6,1],"SampleSizeValid"=A[6,2],
  "SummaryStatisticsMatrix"=A,"ssPreds"=yhat)

  return(Output)
}