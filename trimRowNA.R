trimRowNA<- function(X){
  # Row index for trimming trailing and leading all-NA rows from a matrix
  # D. Meko; last revised 2022-01-07
  #
  # X, [matrix]r the time series; assumed nX columns, no "time" column
  # Returns vector:indices of X marking the row rows remaining 
  #   after any trimming off leading and trailing NAs. But if any of the
  #   columns after trimming has and internal NA, returns a string with
  #   message instead of the vector of row indices.
  #
  # Why? Utility function for manipulating time series matrices of tree-ring
  # prdictors for a reconstruction model

  # Return an error string if X is not matrix
  if (!is.matrix(X)){return("X is not a matrix in trimRowNA()")}
  nX<- dim(X)[2] # number of colunns in X

  # Identify first and last rows that are not all-NA
  i1<- which(apply(X, 1, function(x) !all(is.na(x)))) # identify rows all-NA
  igo1 <- i1[1]; isp1 <- i1[length(i1)]
  
  # Check that no internal NA in any columns after trimming X
  X1<-X[igo1:isp1,,drop=FALSE]
  for (n in 1:nX){
    x<- X[,n,drop=FALSE]
    i2 <- which(complete.cases(x))
    L<-all(diff(i2)==1)
    if (!L){return("Internal NA in one or more cols of X n call to trimRowNA()")}
  }
  
  Output<-igo1:isp1
  return(Output)
}

