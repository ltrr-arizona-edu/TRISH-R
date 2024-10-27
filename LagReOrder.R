LagReOrder<- function(X){
  # Reorder columnms of lagged tree-ring index
  # D. Meko; last revised 2021-12-31
  #
  # X, [m]r lagged matrix of time series, assumed L to R to have unlagged, followed by 
  #   negative lags 1,2,3, ... followed by positive lags 1,2, 3, For example;
  #   [t t-1 t-2 t+1 t+2]
   #
  # Returns a matrix X:
  #   X [m]r reordered version of input X. From L to R have highest negative to highest
  #     positive lag. For example, [t-2  t-1  t   t+1   t+2]. 
  #
  # Why? Utility function to organize matrix of lagged tree-ring chronology to form needed by reconsw4().
  # This utility function intended for input X representing lagged values of a single time series, not
  # of, say, multiple chronologies. The context is single-site reconstruction (SSR).
  
  #--- CHECK INPUT
  
  if (!is.matrix(X)) {stop('X should be matrix')}
  ncols  <- dim(X)[2]
  if ((ncols %% 2)==0 ) {stop('X must have odd number of columns')}
  
  #--- REORDER
  
  jc <- 1 # col number in original of of central col in target
  jLeft  <-  rev(1 + (1:nNeg)) # negative lags ordered L to R highest to lowest
  jRight <- (2+nNeg):ncols  # positive lags ordered L to R lowest to highest
  j <- c(jLeft,jc,jRight)
  

  Output <- X[,j]

  return(Output)
}


 