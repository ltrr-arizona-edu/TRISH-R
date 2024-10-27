Table1Column <- function(D) {
  # Write table with just 1 column of data after column of variable names
  # D. Meko 
  # Last revised 2022-06-26
  #
  # Called from a script or function. For example, was called by function RecLR1 to gene
  # tables of statistics for calibration, validation and analysis of residuals of
  # reconstruction model. Table has (H)ead, (B)ody, and (T)ail.
  #
  # D is list with members:
  #   outputDir: string tells where to write file (e.g.,  "/home/dave/AAAtrish2/test_out/")
  #   textH: header text; vector of strings whose elements are lines of text to go above the
  #     the table. textH[1] is special, will be the first line, and is also (without ".txt"), the
  #     file name of the output table. For example, "Table3-MSR-LR1-Calibration". textH[2] is 
  #     text that, after " - ", will be printed after the filename at the top of the table.
  #     Remaining textH[?] are the column headings. 
  #   textB: variable names; vector of strings to be printed in first column
  #   TfmtB: string with the format for column 1 (e.g., '%10s\t). Be sure that 
  #     1) the number before "s" is at least as long as any string in textB
  #     2) you use "\t" to make tab separation of cols 1 and 2
  #   dataB: corresponding variables as named in calling program; a vector of data variables
  #   DfmtsB: data formats for body; vector of strings that are fprint formats. For example,
  #     "c('%-8g','%-6.2f). Keep in mind:
  #       1) formats in DfmtsB go 1-1 with data in dataB; consider formats accordingly
  #       2) the "-" in the format strings left justify column 1
  #   textT: tail text; vector of strings, each of which is a line in the tail under the table
  #       If do not want a tail, make textT and empty vector (textT <-c()) in calling function
  #
  # Notes
  # Title line is followed by a blank line. Vector textH might have only one element. If more, those
  # lines of head text are followed by another blank line.
  # Table body has a line of "===================" above and below
  # Body is followed by optional "tail," which, if it exists, is separated by a blank line from
  # the "=================" below the body.
  
  
  #=== UNLOAD
  
  textH <- D$textH
  textB <- D$textB
  TfmtB <- D$TfmtB
  dataB <- D$dataB
  DfmtsB <- D$DfmtsB;
  textT <- D$textT
  outDir <-D$outDir
  rm (D)
  
  #--- BUILD OUTPUT FILENAME 

  # Build path/filename for output file
  TableTitle <- textH[1] # will be (with suffix txt) the outfile name, such as
  fnm <- paste(TableTitle,'.txt',sep="") 
  pf1 <- paste(outDir,fnm,sep="")
 
  if (file.exists(pf1)){file.remove(pf1)} # must remove old version of file  
  
  # Line to go above and below table
  baseLine <- "========================================="

  #=== HEADER
  
  for (n in 1:length(textH)){
    if (n==1 | n==length(textH)){
      fmtthis <- '%s\n\n'
    } else {
      fmtthis <- '%s\n'
    }
    vThis = textH[n]
    if (n==1){
      vThis <- paste(vThis,' - ',textH[2])
      fprintf(fmtthis,vThis,file=pf1,append=FALSE)
    } else if (n==2){
      # nothing
    } else {
      fprintf(fmtthis,vThis,file=pf1,append=TRUE)
    }
  }
  fprintf('%s\n',baseLine,file=pf1,append=TRUE)
  rm(vThis, fmtthis, textH, n)
  
  
  #=== BODY
  
  nT <- length(dataB)
  for (n in 1:length(dataB)){
    DfmtB <- DfmtsB[n]
    Tdata  <- dataB[n]
    vName <- textB[n]
    fprintf(TfmtB,vName,file=pf1,append=TRUE)
    fprintf(DfmtB,Tdata,file=pf1,append=TRUE)
  }
  
  #=== TAIL
  
  fprintf('%s\n\n',baseLine,file=pf1,append=TRUE)
  nTail <- length(textT)
  fmtT <- {'%s\n'} # for lines of tail
  for (n in 1:nTail){
    vthis <- textT[n]
    fprintf(fmtT,vthis,file=pf1,append=TRUE)
  }
  
  Output <- NA
  return(Output)
}