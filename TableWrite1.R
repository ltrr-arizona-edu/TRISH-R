TableWrite1 <- function(D) {
  # Write table with any number of columns and rows
  # D. Meko 
  # Last revised 2022-09-08
  #
  # Called from a script or function. For example, was called by function RecMLR1 
  # to write table of correlations of PC scores with some variable y.
  # Table has (H)ead, (B)ody, and (T)ail.
  #
  # D is list with members:
  #   outputDir: string tells where to write file (e.g.,  "/home/dave/AAAtrish2/test_out/")
  #   textH: header text; vector of strings whose elements are lines of text to go above the
  #     the table. 
  #       textH[1] is special, will be the first line, and is also (without ".txt"), the
  #           file name of the output table. For example, "Table4-PCA2"
  #       textH[2] is a short descriptive caption, for example ("Correlation of PCs with RO")
  #       textH[3-?] are column headings
  #   dataB: named list with the table data. Columns must be compatible with textH[3]. Number
  #       of rows in table is computed from contents of DataB.
  #   fmtHB: named list with formats for header row (Head) and table row (Body)
  #   textT: text to go below table (T="tail"). Multi-row easily written by using 
  #         "paste" in combination with \n.
  #  
  
  #=== UNLOAD
  
  outputDir <- D$outputDir
  textT <- D$textT
  textH <- D$textH
  dataB<- D$dataB
  fmtH <- D$fmtHB[1]; fmtB <- D$fmtHB[2]
  BunnyTrack <- D$BunnyTrack
  
  ncol <- length(textH)-2 # number of columns in table (first two elements title a)
  H <- textH[3:length(textH)]
  

  #=== BUILD FILE NAME AND WRITE TITLE AND HEADER
 
  fnm <- textH[1]
  pf = paste(outputDir,fnm,'.txt',sep="")
  xthis <- paste(fnm,': ', textH[2],sep='')
  fprintf('%s\n\n',xthis,file=pf,append="FALSE")
  fprintf('%s\n',BunnyTrack,file=pf,append="TRUE")
  fprintf(fmtH$Head,H,file=pf,append="TRUE")
  fprintf('%s\n',' ',file=pf,append="TRUE")
  
  #=== BUILD BODY
  
  fmtBody <- fmtB$Body
  var1<-dataB[[1]]
  for (n in 1:length(dataB[[1]])){
    # PC number 
    fmt <- fmtBody[1]
    xthis <- var1[n]
    fprintf(fmt,xthis,file=pf,append="TRUE")
    # data vector
    fmt <- fmtBody[2:length(fmtBody)]
    xthis  <- c(dataB$r[n],dataB$Thresh1[n],dataB$Thresh2[n],dataB$r1PC[n])
    fprintf(fmt,xthis,file=pf,append="TRUE")
  }
  fprintf('%s\n\n',BunnyTrack,file=pf,append="TRUE")
  fprintf('%s\n',textT,file=pf,append="TRUE")

  Output <- NA
  return(Output)
}