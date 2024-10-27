EffectSS <- function(x,y) {
  # D. Meko 
  # Last revised 2022-09-09
  # Effective sample size -- effective number of "independent" observations.
  #
  # Effective sample size, Nprime, is computed from the lag-1 autocorrelation of a time series
  # or pair of series. In univariate mode, Nprime can be applied for adjustment of significance
  # of univariate statistics (e.g., uncertainty of the sample mean or variance). In bivariate
  # mode, the effecitive sample size can be applied to adjust signficance of the correlation
  # coefficient for the two series.
  #
  #---IN
  #
  # x: time series matrix or vector without any NA; number of observations mx, number of series nx
  # y: ditto; but if passing vector and matrix, make sure y is the vector and x the matrix
  #
  #---OUT
  #
  # Output: list with fields
  #   Nprime: scalar or vector of effective sample size 
  #   Lflag: flag (logical, length 2)
  #     (1) x or y have at least one NA 
  #     (2) y is not vector and not same col-size as x
  #   ErrorMessage [vector]c : error message associated with Lflag
  #
  #---NOTES
  #
  # If input argument y is NA, Nprime is computed in univariate mode:
  #   If x is vector, Nprime is the scalar effective sample size
  #   If x is a matrix, Nprime is the vector of the effective sample sizes of the individual series
  # If input argument y is not NA, Nprime is computed in bivariate mode:
  #   If y is a vector, Nprime is the effective sample size for correlation of y with all
  #     series in x (Nprime can be scalar or vector, depending on whether x is scalar or vector)
  #   If y is a matrix with ny>1 columns, ny must equal nx, and Nprime is the effective sample size
  #     for correlation of each column of x with same column of y (Nprime is a vector)
  # Method. In univariate mode, if original sample size is N, effective sample size is
  #   Nprime = N(1-r1)/(1+r1), where r1 is the lag-1 autocorrelation
  # Method. In bivariate mode, for pair of series, x and y, effective sample size is
  #   Nprime = N(1-r1r2)/(1+r1r2), where r1 is lag-1 autocorrelation of x and r2 is lag-1
  #   autocorrelation of y

  source(paste(code_dir,"LagkAcc.R",sep="")) # optional transformation of flows
  Lflag <-c(FALSE,FALSE) # initialize as no error flags
  ErrorMessage <- "No problems"
  Nprime <- NA

  klag <-1 # will only need lag-1 autocorrelation
  
  if (!all(complete.cases(x)) | (!all(is.na(y)) & !all(complete.cases(y)))){
    # ERROR MESSAGE
    Lflag[1]<-TRUE
    ErrorMessage <- 'x or y  contain a NA'
    Output <- list(Lflag=Lflag,ErrorMessage=ErrorMessage,Nprime=Nprime)
    return(Output)
  }
  if (is.vector(x)){
    N <- length(x)
    nx <-1
  } else {
    N <- dim(x)[1]
    nx <- dim(x)[2]
  }
  if (all(is.na(y))){
    # Univariate mode for effective sample size
    ResTemp <- LagkAcc(x,klag)
    a <- 1-ResTemp$rk
    b <- 1+ResTemp$rk
    f <- a/b
    Nprime <- floor(f*N)
    L <- ResTemp$rk <= 0
    Nprime[L] <- N # if lag-1 r of either series lag-1 autocorr non-positive, Nprime equals N
  } else {
    # bivariate mode (e.g., for significane adjustment for correlation)
    if (is.vector(y)){
      if (nx ==1){
        # both x and y are vectors
        ResTemp<- rkGet(x,y,klag)
        r1x <- ResTemp$r1x; r1y <- ResTemp$r1y
      } else {
        # y vector, x matrix
        y = matrix(replicate(nx,y),nrow=N) # replicate y to same col-size as x
        ResTemp<- rkGet(x,y,klag)
        r1x <- ResTemp$r1x; r1y <- ResTemp$r1y
      }
    } else {
      # y and x both matrix
      ny <- dim(y)[2]
      if (ny != nx){
        # ERROR MESSAGE
        Lflag[2]<-TRUE
        ErrorMessage <- 'y is matrix and not same col-size as x'
        Output <- list(Lflag=Lflag,ErrorMessage=ErrorMessage,Nprime=Nprime)
        return(Output)
      }
      ResTemp<- rkGet(x,y,klag)
      r1x <- ResTemp$r1x; r1y <- ResTemp$r1y
    }
    rr <- r1x * r1y # products of lag-1 autocorrelation
    f <- (1-rr)/(1+rr)
    Nprime <- floor(f*N)
    L <- r1x <=0 | r1y <=0
    Nprime[L] <- N
  }
  
  Output <- list(Nprime=Nprime,Lflag=Lflag,ErrorMessage=ErrorMessage)
}
rkGet <- function(x,y,klag){
  ResTempx <- LagkAcc(x,klag)
  ResTempy <- LagkAcc(y,klag)
  r1x <- ResTempx$rk
  r1y <- ResTempy$rk
  Out1 <- list(r1x=r1x,r1y=r1y)
}