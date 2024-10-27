TabSepTsm3 <- function(D) {
  # Write general tab-separated file of a time series matrix, with heading line
  # D. Meko 
  # Last revised 2022-09-09
  #
  # Called from a script or function. For example, was called by function RecMLR1, first
  # write PC scores to tab-sep file
  #
  # D is list with members:
  #   outDir: string tells where to write file (e.g.,  "/home/dave/AAAtrish2/test_out/")
  #   filename: name of output file (without txt; for example, PCscoresTimeSeries)
  #   textH: vector of strings to be headers for columns. For example,
  #     "Year","PC1","PC2", etc
  #   dataB: matrix of data (all numeric) for the body. Time (year) in column 1 and other variables in remaining columns
  #   textT: text to be written below the time series listing. The "T" stands for "tail."  This text
  #     can be prepared using "paste," with \n to change lines
  #   fmtsH: vector of string formats for the column headers (e.g., c[1]='%5d')
  #   fmtsB: vector of string formats for a row of the data matrix
  #
  #
  #=== NOTES 
  
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
  fprintf(fmtsB,t(Y),file=pf,append=TRUE)
  
  #--- TAIL
  fprintf('%s\n','',file=pf,append=TRUE)
  fprintf('%s\n',textT,file=pf,append=TRUE)
  Output <- NA
  return(Output)
}