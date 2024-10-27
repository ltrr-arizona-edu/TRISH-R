LagModel2Sign<- function(i1,k,b) {
  # Build string representation of signs of regression coefficient in lagged stepwise regression model. 
  # D Meko
  # last revised 2024-03-04
  # IN
  # i1 [v]j<k: columns of lagged tree-ring matrix in model; some subset of 12345,
  #   in the order from low (1) to high (5) column number (e.g., 234 for cols 2-4). Note this is NOT
  #   arranged to reflect order of entry into the stepwise regression. Note also the regression
  #   coefficients in b (below) go along with columns in i1. Thus b(1) is coefficient on chronology
  #   in i1(1), etc. 
  # k [v]1:  the number of columns in the lagged tree-ring matrix used by calling
  #   function to do the regression. Assumed odd length and symmetric around lag zero. 
  #   For example, if k=5, the matrix would have been at lags {t-2 t-2 t t+1 t+2}
  #   from the predictand. 
  # b[v] regression coefficients, constant term first, then other terms corresponding to cols
  #   i1 of k-column lagged matrix
  #
  # OUT
  # s1 is a character string of length k, with 0's if corresponding lag is not
  # in the model, "'P' if lag is in model with positive coefficient, and "N" if in model
  # with negative coefficient. For example,
  # s1="0NP0P" means lags 0 with positive coefficient, lag t-1 with negative coefficient,
  # and lag t+2 with positive coefficient in model.
  #
  # Why. Lagged single-site regression of flow on a lagged chronology, forward stepwise.
  # Wanted the output tables of SSR models from ReconAnalog to indicate the sign of the
  # coefficients in models (already had the information on which lags are in model and in
  # what order they entered, from function LagModel2Char.R)
  #
  # revised 2024-03-04: opening comments revised for clarity
  
  #---Check input
  L1 <- !(k==floor(k)) ||  !(all(i1==floor(i1))) || k==0
  L2 <- any(i1>k) || any(i1<0)
  L3 <- (k %% 2)==0
  if (L1 || L2 || L3) {
    stop("i1 and k must be integers>0; i1 must be odd length with all elements <=k")
  }
  
  #--- Initialize length-k vector with zeros, meaning lag not in model
  # Fill string with "0"
  si <- rep('0',k); 
  
  # Change some elements of si to "P" or "N" depending on sign of coefficient
  # Coefficients in b are in col-order of s1 according to pointer i1
  b <- b[-1] # strip constant term off model, leaving coefficients
  
  # Positive
  L <- b>0
  iP <- i1[L]
  si[iP] <- 'P'

  # Negative
  L <- b<0
  iN <- i1[L]
  si[iN] <- 'N'
  
  si<-paste(si, sep = '', collapse = '') # concatenate 
  return(si)
}