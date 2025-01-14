crn2tsv <- function(Din){
  # Convert one or more tree-ring "crn" files to tab-sep time series matrix
  # D. Meko; last revised 2025-01-11 -- still being written
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
  #
  #---OUTPUT
  # 
  # pfout: tab-separated tsm with site chronologies
  # Output: a list with elements:
  #   Xfull, Xcut:  data frames of time series matrices of chronologies
  #   Sfull, Scut:  data frames with sample size (number of cores) for above
  #   CodesFull, CodesCut:  character codes of the columns of the matrices
  #     after year column
  #
  #---REFERENCES
  #
  #---REVISIONS
  #
  # revised 2025-01-12.... 
  #
  #---NOTES
  #
  # Assumes "code_dir" in global space. That is where user-written functions stored.
  
  # Libraries needed (calling script should load these)
  #   dplR, treeclim, data.table
  
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
  A <- matrix(NA,nrow=nyr1,ncol=nser1) # data frame for chronologies
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
  
  
  # Bind year column to matrices of chronologies and sample size
  A1 <- cbind(yrA,A); B1 <- cbind(yrB,B) # add year as first column
  colnames(A1) <- c('Year',codesFull); colnames(B1) <- c('Year',codesFull)
  
  #--- Write file of full-length chronologies
  pf <- paste(pathout,fileoutFull,sep='')
  write.table(A1, pf, append = FALSE, sep = "\t", dec = ".",
              row.names = FALSE, col.names = TRUE,quote=FALSE)
  
  
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

  
    # Bind year column to matrices of chronologies and sample size
  A2 <- cbind(yrA,A); B2 <- cbind(yrB,B) # add year as first column
  colnames(A2) <- c('Year',codesCut); colnames(B2) <- c('Year',codesCut)
  
  #--- Write file of truncated and time-screened chronologies
  pf <- paste(pathout,fileoutCut,sep='')
  write.table(A2, pf, append = FALSE, sep = "\t", dec = ".",
              row.names = FALSE, col.names = TRUE,quote=FALSE)
  
  #--- OUTPUT LIST
  
  Output <- list('Xfull'=A1,'Xcut'=A2,'Sfull'=B1,'Scut'=B2,
                 'namesFull'=codesFull,'namesCut'=codesCut)
  return (Output)  
}


