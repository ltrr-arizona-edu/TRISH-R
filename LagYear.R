LagYear<- function(X,tGo,tSp,nNeg=0,nPos=0,ktrim){
  # Create a lagged verson of a time series matrix
  # D. Meko; last revised 2021-12-31
  #
  # X, [matrix] the time series; no time column; could be a single series, or many
  # tGo,tSp:[numeric] first and last times (e.g., years)
  # nNeg [numeric] : maximum negative lag, >=0
  # nPos[numeric]  : maximum positive lag,>=0
  # ktrim (1x1)i:  option for trimming of leading and trailing rows of lagged tsm
  #   ==1 lagged tsm trimmed to exclude any leading or trailing years with a NA
  #         (all-NA are trimmed off regardless of ktrim)
  #   ==2 lagged tsm not so trimmed 
  #
  # Returns a named list, with parts:
  #   X [matrix] lagged version of input X; first the unlagged variables,
  #     then negative lags, then positive lags
  #   tGo, tSp [numeric] start and end times of lagged matrix
  #   ids [character] variable names, coded by column in input X and by lag
  #      Exampple: "X1"   "X1N1" "X1N2" "X1P1" "X1P2"
  #
  # Why? Utility function to prepare matrix of chronologies, or a single 
  # chronology as predictors for reconstruction modeling. Climate can influence
  # growth in multiple years, and so the information on the current 
  # year's climate from the current ring may be conditional on past rings.
  # Likewise, climate can directly influence growth in multiple years, and so 
  # we expect that positive lags of the tree-ring index might improve the 
  # signal on current year's climate beyond what is available from just the 
  # current ring.
  #
  # Accepts X with one or more time series
  # Allows posive and/or negative lags, and accepts no lags
  #
  # revised 2021-12-30: include new input arg ktrim, which will affect year coverage of output tsm
  # revised 2021-12-31: tgo and tsp to "tGo" & "tSp": tsp is a built in R function
  
  #--- Code specific to whether one or more time series in input matrix
  if (dim(X)[2]==1) {
    mX<-nrow(X)
    nX<-ncol(X)
    X0<-X
    mX0<-mX
    nX0<-nX
    
    a<-NA
    # Negative lags
    if (nNeg>0){
      jset<-1:nNeg
      for (j in jset){
        X1<-matrix(a,j,byrow=TRUE)
        X2<-X0[1:(mX-j),]
        X2<-as.matrix(X2)
        X3<-rbind(X1,X2)
        X<-cbind(X,X3)
      }
    } 
    # Positive lags
    if (nPos>0){
      jset<-1:nPos
      for (j in jset){
        X2<-matrix(a,j,byrow=TRUE)
        X1<-X0[(j+1):(mX),]
        X1<-as.matrix(X1)
        X3<-rbind(X1,X2)
        X<-cbind(X,X3)
      }
    }
  } else {
    # Mulitple serie in input matrix
    mX=dim(X)[1]
    nX=dim(X)[2]
    X0<-X
    mX0=mX
    nX0=nX
    
    # Negative lags
    if (nNeg>0){
      jset<-1:nNeg
      for (j in jset){
        X1<-matrix(NA,j,nX)
        X2<-X0[(1):(mX-j),]
        X3<-rbind(X1,X2)
        X<-cbind(X,X3)
      }
    }
    # Positive lags
    if (nPos>0){
      jset<-1:nPos
      for (j in jset){
        X2<-matrix(NA,j,nX)
        X1<-X0[(j+1):(mX),]
        X3<-rbind(X1,X2)
        X<-cbind(X,X3)
      }
    }
  }
  
  #--- Trim leading and trailing rows, and refresh start and end year

  if (ktrim==1){
  if (nNeg>0){
    X<-X[-(1:nNeg),];
  }
  if (nPos>0){
    n1<-nrow(X)
    X<-X[1:(n1-nPos),]
  }
  tGo<-tGo+nNeg
  tSp <- tSp-nPos
  
  } else {
    # no action needed
  }
  mX=nrow(X)
  nX=ncol(X)

  
  #--- Build series ids (column headings
  
  # lag 0
  c1<-c('a')
  jset <- 1:nX0    
  for (j in jset){
    c1[j]<- sprintf("X%s", j)
  }
  
  # negative lags
  if (nNeg>0){
    for (k in 1:nNeg){
      for (j in jset){
        kslot<-nX0+(k-1)*nX0+j
        c1[kslot]<- sprintf("X%sN%s", j,k)
      }
    } 
  }
  
  # positive lags
  if (nPos>0){
    for (k in 1:nPos){
      for (j in jset){
        kslot<-nX0+nNeg*nX0+(k-1)*nX0+j
        c1[kslot]<- sprintf("X%sP%s", j,k)
      }
    } 
  }

  #--- TRIM OFF ANY ALL-NA ROWS
  

  L <- (rowSums(is.na(X))) == nX
  X <- X[!L,]
  yrX <- tGo:tSp
  yrX <- yrX[!L]
  tGo=yrX[1]; tSp <- yrX[length(yrX)]


  
  Output<-list(X=X,tGo=tGo,tSp=tSp,ids=c1)
  return(Output)
}


 