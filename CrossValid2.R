CrossValid2<- function(X, y, nNeg,nPos) {
  # Leave-m-out cross-validation of a regression model. 
  # D Meko
  # last revised 20220104
  #
  # IN
  # X is matrix of predictors
  # y is vector of predictand
  # nNeg and nPos are the maximum negative and positive lags that were considered when
  #   the model was fit. 
  #
  # OUT
  # Output [named list] cross-validation statistics: 
  #   REcv (1x1)r cross-validation reduction of error
  #   CVpredictions [m 1 col]r cross-validation predictions
  #   CVresidual [m, 1 col]r cross-validation residuals
  #   RMSEcv (1x1)r root-mean-square error of cross-validation
  #   LeftOut (1x1)i how many obs left out in each cross-validation model
  #
  # This is simplified from CrossValid1(), which handles various steps of a model
  # previously fit by forward stepwise regression.
  
  library("pracma") # needed for emulation of Matlab "backslash" operator through
  # QR decomposition
  
  source(paste(code_dir,"LeaveOut.R",sep="")) # form pointer matrix for leave-m-out cross-validation

  #--- Build pointer matrix for cross-validation predictor sets
  
  X<-as.matrix(X)
  mX <-dim(X)[1]
  y <- as.matrix(y)
  H<-LeaveOut(nNeg,nPos,mX)
  Lin <- H$Lin # logical pointer matrix; each col marks obs to use as 1 
  mOut <- H$NumberLeftOut
  
  
  #### Cross-validation modeling
  
  #--- Storage 
  w1 <- matrix(NA,mX,1) # to hold cv predictions
  w2 <-  w1 # to hold null predictions (equal to calibration means)
  
  #--- long-term predictor matrix, with 1's in first col
  a1this <- matrix(1,mX,1)
  Xthis <- cbind(a1this,X)
  
  
  for (n in 1:mX){ # Loop over observations
    
    #--- Build predictor matrix
    Lthis <- Lin[,n]
    Lthis <- as.logical(as.matrix(Lthis))
    nthis <- sum(Lthis)
    u <- as.matrix(X[Lthis,])
    a1 <- matrix(1,nthis,1)
    U <- cbind(a1,u) # predictor matrix, this cv model
    
    # Build predictand as 1-col matrix
    v <- as.matrix(y[Lthis])

    #--- Matrix left division to estimate regression parameters
    b <- mldivide(U,v) # [matrix, 1 col, with coefficients, constant term first]

    #--- Estimated predictand for central "left-out" observation
    vhat1 <- Xthis[n,] %*% b
    w1[n] <- vhat1;
    w2[n] <- mean(v)
  }

  #--- Time series of residuals
  e1<-y-w1 # cross-validation 
  e2 <- y-w2 # null-model  (using calib means as predictions)
  
  #--- Validaton statistics
  SSE1  <- sum(e1*e1) # sum of squares of cross-validation errors
  MSE1 <- SSE1/mX  # mean square error of cross-validation
  RMSEcv <- sqrt(MSE1) # root-mean-square error of cross-validation
  SSE2 <- sum(e2*e2) # sum of square of null-model residuals
  REcv <- 1 - SSE1/SSE2 # reduction of error statistic

  Output <- list("REcv"=REcv,"CVpredictions"=w1,"CVresiduals"=e1,"RMSEcv"=RMSEcv,
             "LeftOut"=mOut)
  return(Output)
}