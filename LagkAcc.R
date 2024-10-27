LagkAcc <- function(X,k) {
  # D. Meko 
  # Last revised 2022-09-5
  # Lag-k autocorrelation coefficient(s) of a time series matrix or vector.
  #
  # Converts matrix or vector to z-scores, shifts one version k slots relative to the other,
  # and computes the average product of the first and last mX-k observations of each series,
  # where mX is number of rows, or observations, in X
  #
  #---IN
  #
  # X: matrix or vector without any NA; number of observations mX
  # k: the lag of the desired autocorrelation; require k< mX/4
  #
  #---OUT
  #
  # Output: list with fields
  #   rk: lag-k autocorrelation coefficient(s) [scalar or vector] 
  #   k: lag
  #   Lflag: flag (logical, length 2)
  #     (1) X has at least one NA 
  #     (2) k >= 1/4 the number of observations in X
  #   ErrorMessage [vector]c : error message associated with Lflag
  
  Lflag <-c(FALSE,FALSE) # initialize as no error flags
  ErrorMessage <- "No problems"
  
  # Error message; No NA allowed
  if (!all(complete.cases(X))){
    Lflag[1]<-TRUE
    ErrorMessage <- 'X contains a NA'
    Output <- list(Lflag=Lflag,ErrorMessage=ErrorMessage)
    return(Output)
  }
  
  if (is.vector(X)){
    X <- matrix(X)
  } else {
    # no action needed; X is matrix
  }
  Z <- scale(X,center=TRUE,scale<-TRUE) # time series matrix to zscores
  mZ<- dim(Z)[1]
  
  # Error message; lag k must be less than 1/4 mZ
  if (k>=mZ/4){
    Lflag[2]<-TRUE
    ErrorMessage <- 'Lag k must be less than 1/4 the number of obs in X'
    Output <- list(Lflag=Lflag,ErrorMessage=ErrorMessage)
    return(Output)
  }
  nZ<- dim(Z)[2]
  Zhead <- head(Z,mZ-k)
  Ztail <- tail(Z,mZ-k)
  Z1 <- Zhead * Ztail # products of standardized departures;
  # These standardized departures were computed using the means and standard
  # deviations for the full mZ observation in X
  rk <- colMeans(Z1)
  txt1 <- paste('Lag-',as.character(k),'autocorrelations',sep='')
  Output <- list(rk=rk,lag=k,Lflag=Lflag,ErrorMessage=ErrorMessage)
}