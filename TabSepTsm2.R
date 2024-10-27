TabSepTsm2 <- function(D) {
  # Write tab-separated file of a time series matrix, with heading line
  # D. Meko 
  # Last revised 2022-08-25
  #
  # Called from a script or function. For example, was called by function RecLR1
  #
  # D is list with members:
  #   outDir: string tells where to write file (e.g.,  "/home/dave/AAAtrish2/test_out/")
  #   filename: name of output file (without txt; for example, ObservedAndReconstructedTimeSeries)
  #   header: vector of strings to be headers for columns. For example,
  #     "Year","RO (mm)","meanSSD"
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
  textH <- D$textH
  Y <- D$dataB # year, predictand, predictor
  fnm <- D$filename
  fmtsB <- D$fmtsB; fmtsH<-D$fmtsH
  outDir <- D$outDir
  textT <-D$textT # Tail
  
  rm(D)

  #--- Output file
  pf = paste(outDir,fnm,'.txt',sep="")
  
  #--- Title and header line
  fprintf('%s\n\n',fnm,file=pf,append="FALSE")
  for (n in 1:length(textH)){
    xthis <- textH[n]
    fmt <- fmtsH[n]
    fprintf(fmt,xthis,file=pf,append="TRUE")
  }
  fprintf('%s\n','',file=pf,append="TRUE")

  #--- BODY
  X1 = t(Y)
  fprintf(fmtsB,X1,file=pf,append=TRUE)
  
  #--- TAIL
  fprintf('%s\n','',file=pf,append=TRUE)
  fprintf('%s\n',textT,file=pf,append=TRUE)
  
  Output <- NA
  return(Output)
}