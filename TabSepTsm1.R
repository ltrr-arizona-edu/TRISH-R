TabSepTsm1 <- function(D) {
  # Write tab-separated file of reconstruction with CI, plus observed predictand
  # D. Meko 
  # Last revised 2022-06-27
  #
  # Called from a script or function. For example, was called by function RecLR1
  #
  # D is list with members:
  #   outDir: string tells where to write file (e.g.,  "/home/dave/AAAtrish2/test_out/")
  #   filename: name of output file (without txt; for example, ObservedAndReconstructedTimeSeries)
  #   header: vector of 5 strings to be headers for columns. For example,
  #     "Year","Obs RO (mm)","Reconstruction","Lower 50% CI","Upper 50% CI"
  #   observed: matrix with two columns (year and observed predictand)
  #   recon: 5-col matrix with year, recon, lower 50% CI and upper 50% CI
  #   fmtsH: vector of string formats for the column headers (e.g., c[1]="Year")
  #   fmtsD: vector of string formats for a row of the data matrix
  #
  #
  # Notes
  # 
  # In general, observed could have data extending closer to present then recon. 
  # Not a problem!
  
  
  #=== UNLOAD
  textH <- D$header
  Y <- D$observed
  V <- D$recon
  fnm <- D$filename
  fmtsD <- D$fmtsD; fmtsH<-D$fmtsH
  outDir <- D$outDir
  
  #---Trim to complete cases
  L <- complete.cases(Y)
  Y <- Y[L,]
  L <- complete.cases(V)
  V <- V[L,]
  rm(D)
  
  #--- Allocate storage
  
  yrgo <- min(V[1,1],Y[1,1])
  yrmax1 <- max(V[,1]) ;  yrmax2 <- max(Y[,1])
  yrsp <- max(yrmax1,yrmax2)
  
  yrX <- yrgo:yrsp
  mX <- length(yrX)
  X = matrix(NA,nrow=mX,ncol=5)
  
  
  #--- Fill X
  
  X[,1]=yrX
  
  # obs
  yr <- Y[,1]; y <- Y[,2]
  irow = yr - yrgo + 1
  X[irow,2] <- y
  
  # recon
  yr <- V[,1];  V <- V[,-1]
  irow = yr - yrgo + 1
  X[irow,3:5] <- V
  
  
  #--- OUTPUT FILE
  
  pf = paste(outDir,fnm,'.txt',sep="")
  
  #--- Title
  fprintf('%s\n\n',fnm,file=pf,append="FALSE")
  
  # Header line
  for (n in 1:length(textH)){
    xthis <- textH[n]
    fmt <- fmtsH[n]
    fprintf(fmt,xthis,file=pf,append="TRUE")
  }
 
  # Data lines
  X1 = t(X)
  fprintf(fmtsD,X1,file=pf,append=TRUE)
  
  Output <- NA
  return(Output)
}