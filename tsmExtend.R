tsmExtend<- function(X,yrX,yrsp,N1,N2){
  # Extend time series matrix on recent end by quantile method
  # D. Meko; last revised 2022-01-26
  #
  # X, [matrix]r the time series; assumed nX columns, no "time" column
  # yrX [matrix]i 1-col year matrix for X
  # yrsp: desired stop year of extended time series matrix (tsm)
  #   If yrsp=NA, extend X so that all series have data up to year of most recent series in X
  #   If yrsp=nnnn, where nnnn is any year, extend all series to be full through year nnnn, as long as
  #     nnn is not later than yrX(length(yrX)). If n>yrX(length(yrX)), treat as if yrsp=NA
  #   If yrsp=-1, truncate X to end with the last year with data at all series in X.
  # N1: common period to all series in X should be at least this many years (e.g., 50)
  # N2: each series in X must overlap somewhere with each other series by at least N2 years (e.g., 100)
  #
  # Returns list Output, with fields:
  #   Y: the extended version of input X
  #   yrY: the extended version of input yrX
  #   khow: how the extension was made
  #     =1 no NA in any column of X; no action needed; return Y as X
  #     =2 successful extension, but note that input yrsp later than last year with
  #        data for any series
  #     =3 successful extension as requested
  #     =4 aborted because common period of all X <N1 years
  #     =5 aborted because for some year i of series A needing extension, no
  #         other series that has data in year i overlaps A by N>=N2 years
  #     =6 aborted because yrX not continuous
  #     =7 aborted because first year of X has a NA in some column
  #
  # N1: the common period of N>=N1 years is used for the spearman correlation (r) matrix to set
  #     order of series examined for a quantile prediction
  # N2: a period N>=N2 of overlap of series A and B is examined to find the quantile from that
  #     period of N years from which to pull an extension value
  # Method. Overall common period is used to get spearman r matrix. That matrix allows 
  # for each column of X and ordering of the other series from most-correlated to least-correlated.
  # For given year i of  series A needing "filling" all other series are examined to find which of 
  # the other series has data in year i. Then, in order of decreasing r, the two-series overlap
  # of N>=N2 years is used to pull the overlap. The most highly correlated series that also has data
  # in year i provides the estimate. The value of B in year i is compared with the values of 
  # B in the N>=N2 overlap to get a quantile, or non-exceedance probability. The "prediction" is linearly
  # interpolated from the sorted values of A in the N>=N2 overlap such that the estimate has the same
  # non-exceedance probability in the overlap.
  
  # Debugging stressors
  #X <- X[216:316,,drop=FALSE]; yrX <- yrX[216:316,,drop=FALSE]
  # X <- X[-250,,drop=FALSE]; yrX <- yrX[-250,,drop=FALSE]
  # X [1,4]<- NA
  
  mX <- dim(X)[1] # number of rows in X
  nX <- dim(X)[2] # number of cols in X
  Ifill<-matrix(0,mX,nX) # to hold either 0 or the column of X used as predictor
  
  # Input X should be a tsm that increments by 1
  L1 <- !all(diff(yrX)==1)
  if (L1) {
    khow=6; 
    Output <- list("Y"=NA,"yrY"=NA,"khow"=khow)
    return(Output) # abort year vector, yrX, not continuous
  }
  
  #--- All columns of input tsm X should have valid data in first years
  xtemp<-X[1,]
  if (any(is.na(xtemp))){
    khow<-7
    warning('First row of input tsm X must not have any NA')
    Output <- list("Y"=NA,"yrY"=NA,"khow"=khow)
    return(Output) # abort; a NA in first year of X, some column
  }
  
  
  
  # Check that matrix already not without NA and that yrsp not later than last
  # tree-ring data in any series. 
  khow<-3 # initialize as successful extension, no unusual circumstances
  if (!any(is.na(X))){
    khow<-1
    Y<-X; yrY<-yrX
  } else if (yrsp> yrX[length(yrX),1]) {
    yrsp<-NA
    khow=2
  }
  yrab <- yrX
  if (khow==1){
    # no action needed; Y, yrY, khow already set
  } else {
    
    #--- FIND COMMON PERIOD OF X AND COMPUTE SPEARMAN r MATRIX
    L <- complete.cases(X)
    U <- X[L,,drop="FALSE"]
    yrU <- yrX[L,,drop="FALSE"]
    mU <- dim(U)[1]; nU <- dim(U)[2]
    L1 <- dim(U)[1]<N1
    if (L1) {
      khow=4; 
      Output <- list("Y"=NA,"yrY"=NA,"khow"=khow)
      return(Output) # abort because common period too short
    }
    
    
    R <- cor(U,method="spearman")   # Spearman r matrix, common period
    
    
    #--- MAKE POINTER MATRIX THAT FOR EACH COL OF INPUT X GIVES THE HIGHEST TO
    # LOWEST SISTER SERIES FROM LEFT TO RIGHT 
    
    H <- matrix(NA,nrow=nX,ncol=(nX-1)) # allocate for index
    
    # Loop over rows of r matrix, sorting high to low and filling rows of X
    j1 <- 1:nX
    for (n in 1:nX){
      
      jthis <- j1[-n]
      r <-R[n,]
      r <- r[-n]
      tmp<-sort(r,decreasing=TRUE,index.return=TRUE)
      H[n,]<-jthis[tmp$ix]
    }
    
    
    #--- IDENTIFY ANALOG AND FILL IN
    #
    # Loop over all series (cols) of U. For a series, start with most recent year of missing
    # data and work backwards, filling in values using as predictor series that most highly
    # correlated chronology that has data for tha year. 
    Y<-X; yrY<-yrX   # to hold filled in version of X
    
    for (n in 1:nX){
      ncheck<-n
      ix <- H[n,] # cols of predictor series, L-R from most to least preferred
      xthis <- X[,n,drop=FALSE] # 1-col matrix of key series (needing filling in)
      if (!is.na(xthis[length(xthis)])){
        # no action needed
      } else {
        
        # work back until a year in which data not missing
        kwh<-TRUE
        j<-mX
        while (kwh){
          u <- X[j,ix] # vector of row of data for the other series for year needing filling in series n
          jgood <- !is.na(u);
          jfind <- which(jgood);
          jfind <- jfind[1]; 
          ifind <-ix[jfind]; # predictor col of X   
          
          # Get common period of columns n and ifind of X
          a <- X[,n]
          b <- X[,ifind]
          L <- !is.na(a) & !is.na(b) # mark rows with data for both series
          a <- a[L]; b<-b[L]; yrab <- yrX[L,1,drop=FALSE]
          nab <- length(a)
          
          if (nab<N2){          
            warning('Insufficient overlap (N<N2) of good data at predictor and predictand')
            khow<-5; 
            Output <- list("Y"=NA,"yrY"=NA,"khow"=khow)
            return(Output) # abort because too short overlap of predictor and predictand 
          }
          
          # Get non-exceedance probabiiity of sorted series of length nab
          p = (1:nab)/(nab+1);
          asort<-sort(a) # ascending sort
          bsort<-sort(b)
          
          # Get current value of predictor series; if outside the range of predictor
          # series in common period, set estimate as the max or min of predictand series
          # for the common period
          xhave <- X[j,ifind]
          if (xhave >= max(bsort)){
            apred <- max(asort)
            
          } else if (xhave<=min(bsort)){
            apred<-min(asort)
            
          } else {
            # interpolate by non-exceedance probability
            #phave<- interp1(bsort,p,xhave,method="linear")
            #apred <- interp1(bsort,asort,xhave,method="linear")
            
            apredList <- approx(bsort,asort,xhave,method="linear",ties="ordered")
            apred <- apredList$y
          } # end if (xhave >= max(bsort))
          
          Y[j,n] <- apred # filled in value
          
          Ifill[j,n] <- ifind # corresponding column of predictor
            j<-j-1 # will check next year back
            # if (j==305){
            #   print(j)
            #   browser()
            # }

            if (is.na(X[j,n])){
              # proceed with the estimation
            } else {
              kwh <- FALSE # will exit the while loop
            }
        } # end while
      } # if (!is.nan(uthis[length(uthis)]))

    } #   for (n in 1:nX)
  } # if khow==1
  Output <- list("Y"=Y,"yrY"=yrY,"khow"=khow)
                 
                 
}