trimnan<- function(X){
  # Get row indices of 1-col matrix after trimming of leading and trailing NA
  # D. Meko; last revised 2021-12-24
  #
  # X, [matrix]r the time series; assumed single column
  #
  # Returns vector:
  #   Output: a vector of row indices of X marking the row after any trimming
  #   off of leading and trailing NAs. Internal NAs still allowed.
  #
  # Why? Utility function for eliminating leading and trailing missing values 
  # of a time series
  #
  # Checks that X is a 1-col matrix and not all NA. Internal NAs may still exist in 
  # the values of X pointed to by Output. This function does not deal with the 
  # actual time variable (e.g., year). That is handled at an upper level.
  
  L1 <- is.matrix(X)
  if (!L1){
    stop('X not a matrix')
  }
  nx<- dim(X)[2]
  if (!nx==1){
    stop('X is a matrix, but not 1-column')
  }
  L3<- complete.cases(X)
  if (!any(L3)){
    stop('X cannot be all-NA')
  }
  i3<-which(L3)
  igo<-min(i3)
  isp<-max(i3)
  Output<-igo:isp
}

