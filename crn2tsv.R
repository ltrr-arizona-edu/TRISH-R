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
  
  files <- list.files(path = pathin, pattern = "\\.crn$", full.names = TRUE)
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
  ikeep <- trimRowNA(A) # vector row pointer to keepers
  A <- A[ikeep,]; B<- B[ikeep,]
  yrA <- yrA[ikeep];  yrB <- yrB[ikeep]
  nser <- dim(A)[2]
 
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
  
  # trim columns to chrons with data in yrGo
  i2 <- is.na(A[1,])
  i2 <- !i2
  A <- A[,i2]; B <- B[,i2]
  codesCut <- codesFull[i2]

  # trim off leading or trailing all-NA rows
  ikeep <- trimRowNA(A) # vector row pointer to keepers
  A <- A[ikeep,]; B<- B[ikeep,]
  yrA <- yrA[ikeep];  yrB <- yrB[ikeep]
  nser <- dim(A)[2] # revised nser
  

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
  
  #--- OUTPUT LIST
  
  Output <- list('Xfull'=A1,'Xcut'=A2,'Sfull'=B1,'Scut'=B2,
                 'namesFull'=codesFull,'namesCut'=codesCut)
  return (Output)  
}


