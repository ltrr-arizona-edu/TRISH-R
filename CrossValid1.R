CrossValid1<- function(X, y, nNeg, nPos, i2) {
  # Leave-m-out cross-validation of a previously estimated stepwise regression model. 
  # D Meko
  # last revised 2024-03-08
  #
  # IN
  # X is matrix of potential predictors
  # y is vector of predictand
  # nNeg and nPos are the maximum negative and positive lags that were allowed in
  #   the stepwise regression model of predictand on chronologies. These settings
  #   are used by function LeaveOut to compute m=1+4*max(nNeg,nPos), the nunber
  #   of observations to be omitted in each iteration of leave-m-out cross-
  #   validation
  # i2 is order that columns of matrix of potential predictors entered model
  #   in calling script or function, which is assumed to have applied stepwise
  #   regression; that order is repeated here in the stepwise cross-validation
  # 
  # OUT
  # CV [named list] statistics for maximim-RE model from forward stepwise
  #   Names self-describing. Includes step of maximum RE; columns of the original
  #   matrix of potential predictors in the maximum-RE model; reduction of error
  #   statistic, validation root-mean-square errot, cross-validation predictions,
  #   cross-validation residuals of the final model; and number of observations
  #   left out in leave-m-out cross-validation. Also in CV is the numerical
  #   vector REcvAll of RE at each step of the modeling, and corresponding vector
  #   RMSEVall or RMSE at each step
  #
  # NOTES
  #
  # Differs from CrossValid2 in that computes vector of statistics (e.g. REall) for 
  # each step of a previously run stepwise regression. In contrast, CrossValid2
  # deals just with the simpler case of one model (e.g. scalar RE). CrossValid1
  # also is more complicated than CrossValid1 because CrossValid2 has to identify
  # the "best" model (by max RE)
  #
  # Absolutely important that input column pointer i1 indicates columns of X 
  # in order (left to right) as they entered stepwise in the regression assumed
  # to have been done before calling this function.
  #
  # revised 2024-03-08: cosmetic; expansion and clarification of comments
  
  library("pracma") # needed for emulation of Matlab "backslash" operator through
  # QR decomposition
  
  source(paste(code_dir,"LeaveOut.R",sep="")) # form pointer matrix for leave-m-out cross-validation
  
  
  #--- Build pointer matrix for predictor sets
  
  X<-as.matrix(X)
  mX <-dim(X)[1]
  y <- as.matrix(y)
  H<-LeaveOut(nNeg,nPos,mX)
  Lin <- H$Lin
  mOut <- H$NumberLeftOut
  
  
  #### Cross-validation modeling
  
  #... Storage for models for all steps
  E1 <- matrix(NA,mX,length(i2)) # to store cross-valid residuals for all models
  P1 <- E1 # to hold cross-validation predictions ...
  RMSEvAll <- rep(NA,length(i2)) # to store RMSEv for all models
  REall <- RMSEvAll  # ... RE ...
  
  for (k in 1:length(i2)){
    
    #--- Storage for various 1-col matrices specific to k-step model
    w1 <- matrix(NA,mX,1) # cv predictions
    w2 <-  w1 # null predictions (equal to calibration means)
    ithis <- i2[1:k]
    
    #--- long-term predictor  matrix
    Xthis <- as.matrix(X[,ithis])
    a1this <- matrix(1,mX,1)
    Xthis <- cbind(a1this,Xthis)
    
    
    for (n in 1:mX){
      
      #--- Build predictor matrix
      Lthis <- Lin[,n]
      Lthis <- as.logical(as.matrix(Lthis))
      nthis <- sum(Lthis)
      u <- as.matrix(X[Lthis,ithis])
      a1 <- matrix(1,nthis,1)
      U <- cbind(a1,u) # predictor matrix
      
      # Build predictand as 1-col matrix
      v <- as.matrix(y[Lthis])
      
      #--- Matrix left division to estimate regression parameters
      b <- mldivide(U,v) # [matrix, 1 col, with coefficients, constant term first]
      
      #--- Estimated predictand for central "left-out" observation
      vhat1 <- Xthis[n,] %*% b
      w1[n] <- vhat1;
      w2[n] <- mean(v)
      
    }
    
    #--- Residuals time series
    e1<-y-w1 # cross-validation 
    e2 <- y-w2 # null-model  (using calib means as predictions)
    
    #--- Validaton statistics
    SSE1  <- sum(e1*e1) # sum of squares of cross-validation errors
    MSE1 <- SSE1/mX  # mean square error of cross-validation
    RMSEv <- sqrt(MSE1) # root-mean-square error of cross-validation
    SSE2 <- sum(e2*e2) # sum of square of null-model residuals
    RE <- 1 - SSE1/SSE2 # reduction of error statistic
    
    #--- Store statistics for model step
    E1[,k] <- e1
    P1[,k] <- w1
    RMSEvAll[k]<- RMSEv
    REall[k] <-RE
  }
  

  #--- Find maximum-RE model and its statistics
  kmax <- which.max(REall)  # at this step
  i2cv <- i2[1:kmax] 
  REwinner <- REall[kmax]
  Pcv <- as.matrix(P1[,kmax]) # cv predictions
  Ecv <- as.matrix(E1[,kmax]) # cv errors
  RMSEcv <- RMSEvAll[kmax]
  
  
  CV <- list("REmaxStep"=kmax,"ColumnsIn"=i2cv,"REcvAll"=REall,"REcv"=REwinner,
             "CVpredictions"=Pcv,"CVresiduals"=Ecv,"RMSEvall"=RMSEvAll,"RMSEcv"=RMSEcv,
             "LeftOut"=mOut)
  return(CV)
}