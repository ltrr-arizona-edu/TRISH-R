crn2tsv <- function(Din){
  # Multiple crn files to tsv files of chronologies and sample size
  # Dave Meko, January 2025
  #
  #---INPUT
  #
  # Input is a list (Din) with following fields:
  #
  # pathin: path to input crn files
  # suff.pattern: [string] suffix pattern regular expression <"\\.crn$">
  # pathout: path to output file of tsm of chronologies
  # fileoutFull: name of output file of chronologies, all chrons, full length
  # fileoutCut: name of truncated file of chronologies; starts in year
  #   yrGo and includes only those chronologies with data in yrGo
  # read.crn.arg [list] following input argument to dplR function read.crn (see it)
  #   header <NULL>
  #   encoding <getOption("encoding")>
  #   long <TRUE>
  # yrGo [scalar] desired start year of output time series matrix of
  #   chronologies truncated on early end. That matrix will include only
  #   the chronologies with data extending back to yrGo
  # LNaN  [logical] missing values as NaN (TRUE) or as NA (FALSE)
  # LcolNumbers [logical] <- leading row of column number (TRUE) or not (FALSE)
  # maskChron [vector] chronologies to exclude from analysis
  #
  #---OUTPUT
  # 
  # pfout: tab-separated tsm with site chronologies
  # Output: a list with elements:
  #   Xfull, Xcut:  data frames of time series matrices of chronologies
  #   Sfull, Scut:  data frames with sample size (number of cores) for above
  #   codesFull, codesCut:  character codes of the columns of the matrices
  #     after year column
  #
  #---REFERENCES
  #
  #---REVISIONS
  #
  # revised 2025-01-17: to allow NaN or NA missing value code and to allow
  #   a leading row of column numbers in output tsv files
  # revised 2025-01-31: add Din$maskChron input arg to allow excluding us of specific
  #   crn files in the input data folder (e.g., known bad species)
  # revised 2025-02-01: added storing first and last years non-missing data in each
  #   series, and writing tsv files of those for full and truncated matrix
  #
  #---NOTES
  #
  # Assumes "code_dir" in global space. That is where user-written functions
  # are stored.
  #
  # Libraries needed (calling script should load these)
  #   dplR, treeclim, data.table
  #
  # Example driver script: chron2_UF2025b.R 
  # That script has details of preparation of data and reason for the various
  # options in Din. 
  
  #--- SOURCE FUNCTIONS
  
  source(paste(code_dir,"emssgUNH.R",sep="")) # write error file to system, specified output folder
  source(paste(code_dir,"trimRowNA.R",sep="")) # trim mtx so that no all-NA rows
  source(paste(code_dir,"firstLast.R",sep="")) # years of first and last non-missing values in tsm
  
  
  #--- UNLOAD INPUT
  
  pathin <- Din$pathin
  suff.pattern <- Din$suff.pattern
  pathout <- Din$pathout
  fileoutFull <- Din$fileoutFull
  fileoutCut <- Din$fileoutCut
  read.crn.arg <- Din$read.crn.arg # some args for call to read.crn
  yrGo  <- Din$yrGo
  LNaN  <- Din$LNaN
  LcolNumbers <- Din$LcolNumbers
  maskChron <- Din$maskChron
  
  # Missing value code for output tsv files 
  if (LNaN){
    mvCode <- "NaN"
  } else {
    mvCode <- "NA"
  }

  # logical for writing with col names. 
  # In tsv output file, if write column numbers as row 1, will
  # need to have a separate write.table call to write column names, followed
  # by another write.table call to write the table with col.names==FALSE
  # This is a quirk in write.table that needed workaround
  if (LcolNumbers){
    LcnCode <- FALSE # if writing line of column numbers, want col.names==FALSE
    # in call to write matrix
  } else {
    LcnCode <- TRUE  # otherwise
  }
  
  #--- MAKE LIST OF CRN FILES; STORE CHRONS AND SAMPLE SIZES
  
  # call a function to read names of all files in pathin with specified suffix 
  # into a character, "files", that has length equal to the number of read-in files. 
  files <- list.files(path = pathin, pattern = "\\.crn$", full.names = TRUE)
  
  # Mask (otionally) some of the crn files so that they will not be read
  if (any(is.na(maskChron))){
    # do nothing if maskChron set to NA
  } else {
    L1 <- anyDuplicated(files) || anyDuplicated(maskChron)
    L2 <- any(maskChron>length(files))
    if (L1) stop('Duplicate crn filename in character files, or duplicate maskChron')
    if (L2) stop('Some member of maskChron greater than length(files)')
    files <- files[-maskChron] # omit the chronologies in files pointed to by maskChron
    rm(L1,L2)
  }
  
  data_list <- list() # initialize empty list
  codes <- files # initialize col names 
  codesFull <- codes
 
  # Loop over files, storing data and sample sizes in list, and keeping
  # track of earliest start and most recent end year
  j <- 0
  for (file in files) {
    j <- j+1
    
    # reduce file name to col heading
    g <- files[[j]]
    g1 <- strsplit(g,'/')
    g2 <- g1[[1]]
    ng <- length(g2)
    fname1 <- g2[ng]
    L1 <-  grepl( '_VarStab',fname1, fixed = TRUE)
    if (L1){
      codesFull[j]=strsplit(fname1,'_')[[1]][1]
    } else {
      codesFull[j] <- strsplit(fname1, ".")[[1]][1]
    }
    rm(L1)
    
    # Read crn file and get dimension
    data <- read.crn(file)
    d <- dim(data)
    
    if (j==1){
      FirstYear <- as.numeric(row.names(data)[1])
      LastYear <- as.numeric(row.names(data)[d[1]])
    } else {
      FirstYear <- min(FirstYear,as.numeric(row.names(data)[1]))
      LastYear <- max(LastYear,as.numeric(row.names(data)[d[1]]))
    }
    data_list[[file]] <- data
  }
  
  # Allocate for a data frame with all series
  
  nser1 <- j # total number of read-in chrons
  FirstYear <- as.numeric(FirstYear)
  LastYear <- as.numeric(LastYear)
  nyr1 <- LastYear-FirstYear+1
  A <- matrix(NaN,nrow=nyr1,ncol=nser1) # data frame for chronologies
  mA <- dim(A)[1]
  nA <- dim(A)[2]

  yrA <- FirstYear:LastYear; yrB = yrA
  B <- A # to store sample size
  mB <- mA; nB <- nA
  
  
  #--- BUILD FULL TIME SERIES MATRICES OF VALUES AND SAMPLE SIZE
  
  for (n  in 1:nser1){
    d <- data_list[[n]]
    k = as.numeric(row.names(d))
    irow <- k-FirstYear+1;
    A[irow,n] <- d[,1];
    B[irow,n] <- d[,2];
  }

  # Trim rows so that no leading or trailing all-NA
  # also check that no internal NA in any series
  ikeep <- trimRowNA(A)
  
  # vector row pointer to keepers
  A <- A[ikeep,]; B<- B[ikeep,]
  yrA <- yrA[ikeep];  yrB <- yrB[ikeep]
  nser <- dim(A)[2]

  # Compute first and last year with real data for each column of A (and B)
  colnames(A) <- codesFull
  Res_firstLast.Full <- firstLast(A,yrA)

  # Bind year column to matrices of chronologies and sample size
  A1 <- cbind(yrA,A); B1 <- cbind(yrB,B) # add year as first column
  colnames(A1) <- c('Year',codesFull); colnames(B1) <- c('Year',codesFull)

  #--- Optionally write row of column numbers to start file of full-length chronologies
  if (LcolNumbers){
    v <- 0:nser
    v <- as.matrix(v); v <- t(v)
    pf <- paste(pathout,fileoutFull,sep='')
    write.table(v, pf, append = FALSE, sep = "\t", dec = ".",
                row.names = FALSE, col.names = FALSE,quote=FALSE,
                na=mvCode) 
    h <- colnames(A1)
    h <- as.matrix(h)
    h <- t(h)
    
    # append column names 
    write.table(h, pf, append = TRUE, sep = "\t", dec = ".",
                row.names = FALSE, col.names = FALSE,quote=FALSE,
                na=mvCode) 
  } else {
  }

  #--- Write  remainder of file of full-length chronologies
  pf <- paste(pathout,fileoutFull,sep='')
  write.table(A1, pf, append = !LcnCode, sep = "\t", dec = ".",
              row.names = FALSE, col.names = LcnCode,quote=FALSE,
              na=mvCode)
  

  #--- SAMPLE SIZE OUTPUT TSV FILE, FULL CHRONOLOGY MATRIX
  
  # Optional row with column numbers, col 0 as year
  if (LcolNumbers){
    SSfile <- paste('SS_',fileoutFull,sep='')
    pf <- paste(pathout,SSfile,sep='')
    write.table(v, pf, append = FALSE, sep = "\t", dec = ".",
                row.names = FALSE, col.names = FALSE,quote=FALSE,
                na=mvCode) 
    # append column names 
    write.table(h, pf, append = TRUE, sep = "\t", dec = ".",
                row.names = FALSE, col.names = FALSE,quote=FALSE,
                na=mvCode) 
  } else {
  }
  #--- Write  remainder of file of full-length chronologies
  SSfile <- paste('SS_',fileoutFull,sep='')
  pf <- paste(pathout,SSfile,sep='')
  write.table(B1, pf, append = !LcnCode, sep = "\t", dec = ".",
              row.names = FALSE, col.names = LcnCode,quote=FALSE,
              na=mvCode)

  #--- ROW-TRUNCATED AND POSSIBLY COL-TRUNCATED VERSION OF CHRON MATRIX
  
  # trim rows at start
  i1 <- yrA >= yrGo
  A <- A[i1,]; yrA <- yrA[i1]
  B <- B[i1,]; yrB <- yrB[i1]
  
  # index to trim columns to chrons with data in yrGo
  i2 <- is.na(A[1,])
  i2 <- !i2
  
  # stop with error message if no remaining chrons have data in yrGo
  if (!any(i2)) stop('No remaining series have data in yrGo')
  
  # trim columns to series with data in yrGo, and ensure that if only 1 series
  # A and B are matrices
  A <- A[,i2]; B <- B[,i2]
  codesCut <- codesFull[i2]
  A <- as.matrix(A); B <- as.matrix(B)

  # trim off leading or trailing all-NA rows
  ikeep <- trimRowNA(A) # vector row pointer to keepers
  A <- as.matrix(A[ikeep,]); B<- as.matrix(B[ikeep,])
  yrA <- yrA[ikeep];  yrB <- yrB[ikeep]
  nser <- dim(A)[2] # revised nser
  
  # Compute first and last year with real data for each column of A (and B)
  colnames(A) <- codesCut
  Res_firstLast.Cut <- firstLast(A,yrA)
  
  # Bind year column to matrices of chronologies and sample size
  A2 <- cbind(yrA,A); B2 <- cbind(yrB,B) # add year as first column
  colnames(A2) <- c('Year',codesCut); colnames(B2) <- c('Year',codesCut)
  
  #--- Optionally write row of column numbers to start file of trimmed chronologies
  if (LcolNumbers){
    v <- 0:nser
    v <- as.matrix(v); v <- t(v)
    pf <- paste(pathout,fileoutCut,sep='')
    write.table(v, pf, append = FALSE, sep = "\t", dec = ".",
                row.names = FALSE, col.names = FALSE,quote=FALSE,
                na=mvCode) 
    h <- colnames(A2) # col names of trimmed matrix
    h <- as.matrix(h)
    h <- t(h)

    # append column names 
    write.table(h, pf, append = TRUE, sep = "\t", dec = ".",
                row.names = FALSE, col.names = FALSE,quote=FALSE,
                na=mvCode) 
  } else {
  }
  
  #--- Write file of truncated and time-screened chronologies
  pf <- paste(pathout,fileoutCut,sep='')
  write.table(A2, pf, append = !LcnCode, sep = "\t", dec = ".",
              row.names = FALSE, col.names = LcnCode,quote=FALSE,
              na=mvCode)

  
  #--- SAMPLE SIZE OUTPUT TSV FILE, CUT CHRONOLOGY MATRIX
  
  # Optional row with column numbers, col 0 as year
  if (LcolNumbers){
    SSfile <- paste('SS_',fileoutCut,sep='')
    pf <- paste(pathout,SSfile,sep='')
    write.table(v, pf, append = FALSE, sep = "\t", dec = ".",
                row.names = FALSE, col.names = FALSE,quote=FALSE,
                na=mvCode) 
    # append column names 
    write.table(h, pf, append = TRUE, sep = "\t", dec = ".",
                row.names = FALSE, col.names = FALSE,quote=FALSE,
                na=mvCode) 
  } else {
  }
  #--- Write  remainder of file of full-length chronologies
  SSfile <- paste('SS_',fileoutCut,sep='')
  pf <- paste(pathout,SSfile,sep='')
  write.table(B2, pf, append = !LcnCode, sep = "\t", dec = ".",
              row.names = FALSE, col.names = LcnCode,quote=FALSE,
              na=mvCode)
  
  
  #--- WRITE TAB-SEP FILES OF THE FIRST AND LAST YEAR TABLES
  # 
  # full matrix
  Yearsfile <- paste('Years_',fileoutFull,sep='')
  pf <- paste(pathout,Yearsfile,sep='')
  write.table(Res_firstLast.Full, file = pf, append = FALSE, quote = FALSE, sep = "\t",
              eol = "\n", na = "NA", dec = ".", row.names = TRUE,
              col.names = FALSE, qmethod = c("escape", "double"),
              fileEncoding = "")
  
  # cut matrix
  Yearsfile <- paste('Years_',fileoutCut,sep='')
  pf <- paste(pathout,Yearsfile,sep='')
  write.table(Res_firstLast.Cut, file = pf, append = FALSE, quote = FALSE, sep = "\t",
              eol = "\n", na = "NA", dec = ".", row.names = TRUE,
              col.names = FALSE, qmethod = c("escape", "double"),
              fileEncoding = "")
  
  
  
  
  
  
  
  
  #--- OUTPUT LIST

  Output <- list('Xfull'=A1,'Xcut'=A2,'Sfull'=B1,'Scut'=B2,
                 'namesFull'=codesFull,'namesCut'=codesCut,
                 'yearsFull' = Res_firstLast.Full,
                 'yearsCut' = Res_firstLast.Cut)
  return (Output)  
}


