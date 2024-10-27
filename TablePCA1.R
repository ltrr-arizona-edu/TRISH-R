TablePCA1 <- function(D) {
  # Write PCA summary table, including list of loadings
  # D. Meko 
  # Last revised 2022-09-07
  #
  # Called from a script or function. For example, was called by function RecLR1
  #
  # D is list with members:
  #   outDir: string tells where to write file (e.g.,  "/home/dave/AAAtrish2/test_out/")
  #   filename: name of output file (without txt; for example, Table3-PCA1.txt
  #   header: vector of strings to be headers for columns. For example,
  #     "N","Site#","SiteID","PC1", "PC2", etc
  #   P: list with loadings and other data needed
  #   fmtsH: vector of string formats for the column headers (e.g., c[1]="Year")
  #   fmtsB: vector of string formats for a row of the body of table
  #     1 applies to cols 1 and 2) "N" and "Site#"
  #     2 applies to col 3 (site ID)
  #     3 applies (after replicating) to remaing columns (e.g., '%-12.8g\t)
  #
  #
  # Notes
  # 
  # Column header and numeric format widths should be compatible; moreover,
  # Col 1 should have width no smaller than 4 (to accomdate "Cum%')
  # Col 2 should have width no smaller than 5 (to accomodate "Site#")
  # Col 3 should be sized as large as length of longest site ID
  
  #=== UNLOAD
  textH <- D$textH$Heading
  DfmtsB <- D$DfmtB
  DfmtsH <- D$DfmtH
  fnm <- D$textH$Title
  tit1 <- D$textH$SubTitle
  fmtsB <- D$fmtB; fmtsH<-D$fmtH
  TfmtB <-D$TfmtB
  outDir <- D$outDir
  textT <-D$textT # Tail
  Y <- D$dataB$ResPCA$Loadings
  EV <- D$dataB$ResPCA$EigValues
  PctVar <- D$dataB$ResPCA$PctVar
  CumPctVar <- D$dataB$ResPCA$CumPctVar
  jScreened <- D$dataB$jScreened
  SiteID <- D$dataB$SiteID
  textT <- D$textT
  BunnyTrack=D$BunnyTrack
  rm(D)
  
  
  nPC <- dim(Y)[2] # how many PCs
  nSites <- nPC;  # same number of variables
  jPC <- 1:nSites
  
  #--- Output file
  pf = paste(outDir,fnm,'.txt',sep="")

  #--- Title and header line
  fprintf('%s\n\n',paste(fnm,' - ',tit1,sep=""),file=pf,append="FALSE")
  fprintf('%s\n',BunnyTrack,file=pf,append="TRUE")
  
    # first three col headers
  for (n in 1:3){
    xthis <- textH[n]
    fmt <- TfmtB$Left[n]
    fprintf(fmt,xthis,file=pf,append="TRUE")
  }
  # headers for PCs
  for (n in 4:length(textH)){
    xthis <- textH[n]
    fmt <- TfmtB$Right[[n-3]]
    fprintf(fmt,xthis,file=pf,append="TRUE")
  }

#--- BODY, PART 1: LOADNGS

for (n in 1:nSites){
  # First 3 columns
  xthis  <-c(jPC[n],jScreened[n])
  # The sequential number and the database site numer
  for (m in 1:2) {
    fmt <- DfmtsB$Left[m]
    fprintf(fmt,xthis[m],file=pf,append="TRUE")
  }
  # Site id
  fmt <-DfmtsB$Left[3]
  fprintf(fmt,SiteID[n],file=pf,append="TRUE")
  # Loadings on this site
  xthis <- Y[n,]
  fmt  <- DfmtsB$Right
  fprintf(fmt,xthis,file=pf,append="TRUE")
}

fprintf('%s\n','',file=pf,append="TRUE") # blank line

# Eigenvalue line
fprintf('%s\t\t','Eigenvalue         ',file=pf,append="TRUE") 
fmt  <- DfmtsB$Right
fprintf(fmt,EV,file=pf,append="TRUE")

# Lines for %  Var. and Cum % Var.Eigenvalue line
fmt <- DfmtsB$Pctg
fprintf('%s\t\t','Pctg Variance     ',file=pf,append="TRUE") 
fprintf(fmt,PctVar,file=pf,append="TRUE")

fprintf('%s\t\t','Cum. Pctg Variance',file=pf,append="TRUE") 
fmt <- DfmtsB$Pctg
fprintf(fmt,CumPctVar,file=pf,append="TRUE")

fprintf('%s\n',BunnyTrack,file=pf,append="TRUE")

  #--- TAIL

  fprintf('%s\n','',file=pf,append=TRUE)
  fprintf('%s\n',textT,file=pf,append=TRUE)
  Output <- NA
  return(Output)
}