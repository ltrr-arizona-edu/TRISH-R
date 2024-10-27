TranFlow<- function(x,ktran){
  # Transformation of flows or some other single time series.
  # D. Meko; last revised 2021-12-26
  #
  # x, [matrix]r the time series (e.g., flows); assumed single column
  # ktran [scalar]i type of transform requested
  #   =1: no tranform requested
  #   =2: square root (valid iff x non-negative)
  #   =3: log10 (valid iff x all-positive)
  #
  # Returns Output, a named list, with parts:
  #   x [matrix] the (possibly) transformed verstion of input x. If not possible
  #     to apply the requested transform, returns original x
  #   flag [scalar]i  flag indicating if transform applied, and why not if not
  #     =0 No problem; transform applied
  #     =1 Square root transform rejected because some x non-negative
  #     =2 Log10 transform rejected becase x not all-positive
  #   Transformed [logical]1x1  whether flow transformed or not (T=yes,F=No)
  #
  # Why? For flow reconstruction, a transformation sometimes gives a stronger reconstruction
  # as measured by calibration statistics, or avoids the problem of highly skewed residuals, or
  # residuals whose variance strongly depends on the size of the predicted value (non-constant variance).
  #
  # Bombs if input x is all-NA or is not a 1-col matrix. Checks whether desired
  # transform mathematically sensible (e.g., log10(0) is minus infinity) and if not,
  # returns original x, and flag indicating why could not transform. 
  #
  # Missing values: input x may conatin some NAs. If so, corresponding elements of 
  # transformed x are NA.
  #
  # Only two types of transformatons are supported: square root and log10 transform

  # Check input series
  L1 <- is.matrix(x)
  if (!L1){
    stop('x not a matrix')
  }
  nx<- dim(x)[2]
  if (!nx==1){
    stop('x is a matrix, but not 1-column')
  }
  
  if (ktran==1){
    sTran<-''
  } else if (ktran==2) {
    sTran<- "(sqrt-transformed)"
  } else if (ktran==3) {
    sTran <- '(log10-transformed)'
  } else {
    sTran<- '(Invalid ktran choice)'
  }


  #--- TRANSFORM
  v<-x
  Transformed<-FALSE
  flagTr<-0
  if (ktran==1){
    v<-v
  } else if (ktran==2) {
    if (!any(v<0)){
      v<-sqrt(v)
      Transformed<-TRUE
    } else {
      flagTr<-1
      sTran<-''
    }
  } else  {
    if(!any(v<=0)){
      v=log10(v)
      Transformed<-TRUE
    } else {
      # Write error file
      flagTr<-2
      sTran<-''
    }
  }
  
  Output<-list(x=v,flag=flagTr,Transformed=Transformed,sTran=sTran)
  return(Output)
}

