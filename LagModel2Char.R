LagModel2Char<- function(i1,k) {
  # Build string represntation of a lagged stepwise regression model. 
  # D Meko
  # last revised 20220105
  #
  # IN
  # i1 [v]j<k: the colunns ordered left-to-right of the lagged tree-ring matrix
  #   in the order that they entered forward stepwise
  # k [v]1:  the number of columns in the lagged tree-ring matrix used by calling
  #   function to do the regression. Assumed odd length and symmetric around lag zero. 
  #   For example, if k=5, the matrix would have been at lags {t-2 t-2 t t+1 t+2}
  #   from the predictand. 
  #
  # OUT
  # s1 is a character string of length k, with 0's if corresponding lag is not
  # in the model. The number at each position gives the order that a particular lag
  # entered the forward stepwise regression.  For example,
  # s1="02103" means lag 0 entered first, followed by lag t-1, followed by lag t+2
  #
  # Why. Lagged single-site regression of flow on a lagged chronology, forward stepwise.
  # Needed a simple way to show lags in the model and the order of their entry. The
  # string output s1 can be used in a table reporting results. 

  
  #---Check input
  L1 <- !(k==floor(k)) ||  !(all(i1==floor(i1))) || k==0
  L2 <- any(i1>k) || any(i1<0)
  L3 <- (k %% 2)==0
  if (L1 || L2 || L3) {
    stop("i1 and k must be integers>0; i1 must be odd length with all elements <=k")
  }
  
  #---Fill a numeric vector version of desired string
  ni1 <- length(i1)
  ModelCoded<-vector(mode<-"numeric",length=k) # initialize as vector of zeros
  for (n in (1:ni1)){
    islot<-i1[n]
    ModelCoded[islot]<-n # ModelCoded, still numeric vector
  }
  s1<-as.character(ModelCoded) # convert numweric vector to character
  s1<-paste(s1, sep = '', collapse = '') # concatenate 
  return(s1)
}