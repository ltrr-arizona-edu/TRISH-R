Tsm2Scores1 <- function(D) {
  # Time series matrix to PC scores
  # D. Meko 
  # Last revised 2022-08-28
  #
  # Time series matrix to PC scores
  #
  #---IN
  #
  # D is list with members:
  #   X: time series matrix (data frame); no NA
  #   yrX: years for X (vector)
  #   nmsX: names of variables in X (vector)
  #   khow: how PCA done (1=correl mtx, 2=covariance mtx)
  #  
  #--- OUT -- named list "Output" with elements:
  #
  # Scores: PC scores (data frame)
  # yrScores: year for Scores (vector)
  # Loadings: loading matrix from PCA (data frame)
  # EigValues: eigenvalues (vector)
  # PctVar: percentage of variance of X accounted for by each PC (vector)
  # CumPctVar: cumulative of PctVar (vector)
  # Eig1Cutoff: last most important PC, by equivalent eigenvalue-of-1 threshold (scalar)
  #
  #--- NOTES
  # 
  # Written for TRISH as function to be called by ReconAnalog-RecPCR1. Assumes X 
  # has at least two columns.
  #
  # Uses singular value decomposition of the time series matrix, and R function prcomp.
  # Eig1Cutoff is strictly applicable to PCA on the correlation matrix, or singular
  # value decomposition (SVD) of a data matrix that is centered and scaled such that
  # each column has unit variance. If the number of variables in X is p, and khow=1, 
  # the eigenvalues sum to p, and the "average" eigenvalue is 1. Any eigenvalue greater
  # than 1 represents "greater than average" variance. 
  #
  # If PCA is done on the covariance matrix, the matrix X is centered, but not scaled,
  # prior to PCA. The centered variables (columns of X) retain their original variance.
  # The sum of the variances of the columns of X is identically equal to the sum of 
  # the variance of the PC scores derived by PCA on the covariance matrix of X. With
  # khow=2, Output$Eig1Cutoff is set to however many PCs have greater variance than
  # the average variance of the PC scores. 
  
  
  library(resample) # function colVars
  
  #=== UNLOAD
  
  X <- D$X; yrX <- D$yrX
  nmsX <- D$nmsX; khow <-D$khow
  rm(D)

  #=== PRINCIPAL COMPONENTS ANALYSIS
  
  # User function for centering a matrix
  center_colmeans <- function(x) {
    xcenter = colMeans(x)
    x - rep(xcenter, rep.int(nrow(x), ncol(x)))
  }
  
  # PCA and compute scores, which are not an element returned by prcomp
  if (khow==2){
    # PCA and compute scores
    # PCA on covariance matrix
    P <- prcomp("x"=X,"retx"=TRUE,center=TRUE,scale=FALSE)
    Xc <- center_colmeans(X)
    # To compute scores, need to multiply centered data times loadings
    F <- Xc %*% P$rotation; yrF <-yrX   # scores
  } else {
    # PCA on correlation matrix
    P <- prcomp("x"=X,"retx"=TRUE,center=TRUE,scale=TRUE)
    # Scores are z-scored data time loading
    df <- data.frame(X)
    Xz <-sapply(df, function(df) (df-mean(df))/sd(df))
    F <- Xz %*% P$rotation; yrF <-yrX   # scores
  }
  

  # P$rotation: loadings are returned by prcomp
  
  # Compute eigenvalues. These should be the columnn variances of the scores, F.
  # Checked with my Matlab eigen: results match to 5 digits; some differences 
  # after that. Conclude OK. 
  EV <- colVars(F) # eigenvalues
  Pct1 <- 100*EV/sum(EV) # pctg variance in each PC
  CumPct1 <- cumsum(Pct1) # cumulative pctage variance
  
  # "Eigenvalue of 1" cutoff
  # For PCA on covariance matrix, I use mean of EVs as the threshold
  EigThresh <- mean(EV)
  L <- EV>EigThresh
  EigCut <- max(which(L))
  
  Output <- list('Scores'=F,'yrScores'=yrF,'Loadings'=P$rotation,'EigValues'=EV,
                 'PctVar'=Pct1,'CumPctVar'=CumPct1,'Eig1Cutoff'=EigCut)
  return(Output)
}