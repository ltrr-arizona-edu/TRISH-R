KnnAnalog <- function(U,yrU,v,yrv,kNN) {
  # Analog reconstruction from Euclidean distance
  # D. Meko
  # Last revised 2022-09-17
  #
  #--- IN
  #
  # U [matrix] time series matrix for which analog values wanted (typically PC scores)
  # yrU: vector or 1-col time series matrix of years
  # v, yrv: vectors of data and years to supply analogs for each row of U
  #`kNN: [scalar]: find the kNN nearest analogs 
  #
  #--- OUT
  #
  # Output: list with three fields:
  #   YearsAnalog [matrix]: covering same years as input U, with year as column 1 and 
  #     analog years (from yrv) for kNN nearest neighbors in nex kNN columns. Ordered left to
  #     right from nearest neighbor.
  #   DataAnalog [matrix] same size as YearsAnalog, with the corresponding data from v
  #   Distance [matrix] same size as YearsAnalog, with the Euclidean distances
  #   Recon [matrix] three-column matrix, same years as yrU, with year and the analog reconstruction
  #     Col 1: year
  #     Col 2: analog year, from nearest neighbor (col 2 YearsAnalog) for that part of U not
  #       overlapping v, and from the second-nearest neighbor (col 3 of YearsAnalog) for the years
  #       of U overlapping v. 
  #     Col 3: corresponding analog data values from columns 2 or 3 of DataAnalog
  #     Col 4: either 1 or 2, indicating that analog year and data comes from nearest neighbor (1)
  #       or second nearest neighbor (2)
  # NearestObserved: nearest observed v as predictor for v
  #     Col 1-2: input yrv and v repeated
  #     Col 3: nearest-observed v 
  #     Col 4: year of nearest-observed v
  #
  #=== NOTES
  #
  # Written to get "analog" reconstruction for TRISH by calls from ReconAnalog-RecMLR1. In that
  # application, U are PC scores of a time series matrix of screened single-site reconstructions
  # of some hydrologic variable v.
  #
  # It is possible for U to have more recent data than v. If so, the most recent analog reconstructed
  # values will come from the nearest neighbor rather than the second nearest neighbor
  # 
  #---- Method
  #
  # Find Euclidean distance from each row or U to all other rows
  # Use the subset of Euclidean distances to rows in the overlap of U and v to assign kNN
  #   nearest neighbors (similar rows in U from that overlap)
  # For years of U not overlapping v, assign the analog reconstruction as the nearest neighbor
  #   value of v
  # For years of U overlapping, assign the analog reconstruction as the sencond nearest neighbor
  #   value of v. This step is needed because the nearest neighbor for years in that overlap would
  #   by definition be the years itself -- yielding a perfect reconstruction for years yrv
  # Store matrices of results, including the kNN years, values and distances
  U <- as.matrix(U) # in case only one time series in U
  mU <- dim(U)[1]; nU<- dim(U)[2];  # number of years and variables (e.g., PCs) in U
  mv <- length(v) # number of years in v
  
  # Compute Euclidean distance of each row of U to all other rows
  H <- dist(U,method="euclidean") 
  H <- as.matrix(H) # want matrix; this will be mU x mU
  
  #--- GET DISTANCE SUB-MATRIX FOR OVERLAP
  # overlap of U and v
  L1 <- yrU >= yrv[1] & yrU <= yrv[length(v)]
  X <- H[L1,] # Distance matrix for overlap of U and v
  yrX <- yrU[L1]
 
  #--- SORT THE DISTANCE SUB-MATRIX FROM NEAREST TO FARTHEST 
  Sort1 <- function(x){
    y <- sort(x,method='quick',index.return=TRUE)
  }
  A <- apply(X,2,Sort1)
  # A is a list, with A[[1]]$x a vector of the sorted neighbor values of v for yrU[1]
  # and A[[1]]$ix a vector of the corresponding year 
  # yrv[A[[1]]$ix] is a vector of the corresponding analog year s for yrU[1]
  
  #-- ALLOCATE SOME MATRICES FOR RESULTS
  Y1  <- matrix(data=NA,nrow=mU,ncol=kNN+1) # for analog data
  Y2 <- Y1 #  for analog years
  Y3 <- Y1 # for analog distances
  Y4 <- matrix(data=NA,nrow=mU,ncol=4) # for year, blended analog recon, analog year, and whether
  #   from nearest or second-nearest neighbor
  nY1 <- kNN+1; nY2 <- nY1; nY3 <- nY1 # number of columns in Y1, Y2, Y3
  Y1[,1]=yrU;  Y2[,1] <- yrU; Y3[,1] <- yrU; Y4[,1] <- yrU;  # first col of these tsm's is year
  
  
  #---  FILL MATRICES Y1, Y2 AND Y3
  
  for (n in 1:mU){
    y1  <- v[A[[n]]$ix[1:kNN]] #  vector of analog values from v
    y2  <- yrv[A[[n]]$ix[1:kNN]] # vector of analog years
    y3  <- A[[n]]$x[1:kNN] # vector of kNN distances
    Y1[n,2:nY1] <- y1
    Y2[n,2:nY2] <- y2
    Y3[n,2:nY3] <- y3
  }
  
  #--- FILL MATRIX Y4
  
  # Initialize with all years using nearest neighbor
  Y4[,2] <- Y1[,2];
  Y4[,3] <- Y2[,2]
  Y4[,4] <- 1
  
  # But second nearest for overlap of U and v
  Y4[L1,2] <- Y1[L1,3]
  Y4[L1,3] <- Y2[L1,3]
  Y4[L1,4] <- 2
  
  
  #=== NEAREST-OBSERVED ANALYSIS
  #
  # Fill a vector w with the values of v closest to those in each yrv.
  V <- matrix(replicate(length(v),v),nrow=length(v)) #replicate cols of v
  W <- t(V)
  D <- abs(V-W)
  A <- apply(D,1,Sort1) # sort rows
  # A is a list, with A[[1]]$x a vector of the sorted neighbor values of v for yrv[1]
  # and A[[1]]$ix a vector of pointer to the corresponding year in yrv
  # yrv[A[[1]]$ix] is a vector of the corresponding nearest-observed years for yrv[1]
  
  W <- matrix(data=NA,nrow=mv,ncol=4) # for year, v, nearest-observed v,
  # and year of nearest observd v
  W[,1]=yrv; W[,2]=v
  for (n in 1:mv){
    X <- A[[n]] 
    irow <- X$ix[2]
    W[n,3] <- v[irow]; W[n,4] <- yrv[irow]
  }
 
  #=== OUTPUT LIST
  
  Output <- list(YearsAnalog=Y2,DataAnalog=Y1,Distance=Y3,Recon=Y4,
                 NearestObserved=W)
  return(Output)
}
