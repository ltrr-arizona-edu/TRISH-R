xyCI <- function(X){
  # Upper and lower CI into x,y for CI polygon
  # D. Meko
  # Revised 2022-09-23
  #
  # Upper and lower CI into x,y for CI polygon. Utility function to get x,y polygon
  # coordinates for, say, plotting a shaded confidence interval
  #
  #---IN
  #
  # X: [matrix] 3-column time series matrix with year, lower CI and upper CI
  #
  #---OUT
  #
  # Output: a list with elements:
  #   x [vector] x coordinates for CI polygon
  #   y [vector] y coordinates for CI polygon
  
  #---CHECK THAT INPPUT X IS 3-COL MATRIX WITH VALUES IN COL 1 INCREMENTING BY 1
  L1 <- is.matrix(X)
  L2 <- dim(X)[2]==3
  
  if (!L1 && L2){
    stop('X must be 3-col matrix as input to xyCR')
  }
  
  
  #---UNLOAD INPUT
  
  x1 <- X[,1] # time vector
  # check that time increments by 1
  d <- diff(x1)
  if (!all(d==1)){
    stop('Time column in input X to xyCI must increment by 1')
  }
  u1 <- X[,2] # lower CI (e.g., lower 50% CI)
  u2 <- X[,3]  # upper CI (e.g., upper 50% CI)

  #--- MAKE POLYGON COORDINATES
  
  # reorganize the data values for the lower and upper Ci into a vector whose end value equals the first value
  yP  <- cbind(u1,rev(u2)) # matrix with lower CI as col 1, reversed upper CI as col 2
  yP <- c(yP) # convert matrix to vector:  col 2 appended to end of col 1
  yP <- append(yP,yP[1]) # append first year's value of lower CI as the last value in the vector
  
  # Similarly reorganize the time (e.g., year) vector
  xP  <- cbind(x1,rev(x1))
  xP <- c(xP)
  xP <- append(xP,xP[1])
  Output <- list('x'=xP,'y'=yP)
  return(Output)
}