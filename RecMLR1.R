RecMLR1 <- function(D) {
  # Multi-site reconstruction (MSR) by multiple linear regression on SSRs or their PCs
  # D. Meko 
  # Last revised 2024-03-06
  #
  # Called from a script or function (e.g., ReconAnalog1) that has generated the SSRs.
  # Consider the predictand, y, and predictor matrix, X. Here, y is regressed on 
  # either the matrix of screened SSRs or PCs of that matrix. The method uses
  # stepwise forward multiple linear regression. 
  #
  # D is list with members:
  #   Text (1x4)s: symbol for y label; units of y, in parens; longer name of y; name of y season; Example:
  #    "RO","(mm)","Runoff", "12-month season ending in month 9"
  #   U, yrU: matrix of screened SSRs (matrix); matrix of years (1-col matrix)
  #   nmsU: names of screened SSRs (vector)
  #   jScreened: pointer from columns of screened SSRs to site number in original user network
  #   v, yrv: numeric and integer; predictand and years
  #   yrsC (1x2)d   first and last year of desired calibration period; if NaN, default to first or last
  #     available year of overlap of u and v
  #   yrEnd (1x1)d desired last year of provided time series output of reconstruction and of plot of 
  #     reconstruction with 50% error bar
  #   nNeg, nPos: integer (both positive) of max negative and possitive lag allowed by calling function in
  #     SSR modeling. This affects m in leave-m-out cross-validation
  #   incR2a: the stwp of "approximate" maximum adjusted R-squared is the before which the increase
  #     in adjusted R-squared is less than incR2a. I have used 0.01, figuring that, say, an
  #     increase in adjusted R-squared from 0.50 to 0.5099 is not worth a more complicated model
  #   kstop: stopping rule for forward stepwise. 
  #     =1 at approximate maximum adjusted R-squared
  #     =2 at maximum cross-validation RE, but at no higher step than that of approximate maximum adjusted R-squared
  #   NcMin: minimum acceptable number of years for calibration of MSR model
  #   PCoption: option for building pool of potential predictors from scores of PCs of SSRs
  #     1: specify k most important PCs, for example from viewing scree plot
  #     2: the m<fN PCs most highly correlated with predictand; N=number of years in calibration period
  #         of MSR model; f LE 0.20, defaul f=0.10
  #   f: factor that multiplied times number of calibration years puts cap on number of PCs allowed in 
  #       pool of potential predictors of MSR model (see PCoption)
  #   alphaR: If using analog method (methMSR=3), only those PCs whose correlation with y are significant
  #     at this alpha level (default 0.05) are retained for use in selection of analog years
  #   PCApredictors (TorF) whether to use PCs of SSRs as predictors
  #   methMSR: method of reconstruction (1=SLR1, 2=MLR1-PCA or MLR1-noPCA), 3=Analog
  #   PdfDescribe: string referring to pdf file with description of method; used in table notes
  #   nPCsKeep: user-specified number of "most important" PCs of SSRs to include in pool of potential predictors
  #   kHowPCA: option for PCA on correlation (1) or covariance (2) matrix
  #   ScreenAnalogPCs: Logical to screen the PCs used in analog reconstruction (methMSr=3)
  #   by correlation with predictand.If TRUE, only those PCs whose correlations with y
  #   are significant at alpha-level alphaR (see earlier) are used to identify analogs
  #
  #   NextFigNumber [i] # start naming figures as Figure0?.png, where ? is NextFigNumber
  #   outputDir: folder that output files are written to <"/home/dave/AAAtrish2/test_out/">
  #
  # Notation below uses "<>" to indicate "default" values of input specifications. For example, 
  # <f=010> means the default for variable f is 0.10.  Or, with options 1 and 2, <2> indicates that
  # the default option is "2."
  #
  # Revision
  # Rev 2023-02-09.  Minor, to avoid fatal error in in a listing annotated on a figure.
  # Rev 2023-03-29.  format '-3d' to '3.0f' in response to error with analog method call for calb table
  # Rev 2023-04-03.  To handle special case in which specified start and end years of calibration for
  #   SSR models are incompatible with data available for MSR model
  # Rev 2023-04-15. Typo "nPCskeep" corrected
  # Rev 2023-05-08. Cosmetic, to figure summarizing calibration statistics. Switched to using layout()
  #   for screen splitting to allow sub-windows to be of different width. Formerly used split.screen.
  # Rev 2023-05-16.  yrEnd newly provided input argument, for truncation of final reconstruction
  # Rev 2023-06-02. remove error message and automatic bail when user specifies a calibration period impossible
  #   for the time coverage of tree-ring and hydro data. Now the calibration period is simply truncated
  #   to be as long as possible, and program allowed to proceed. 
  # Rev 2023-11-26. (1) To used input arg yrEnd as last year of desired reconstruction output rather
  #   than compute yrEnd internally. (2) To fix a fatal error when running in "analog" mode
  # revised 2024-03-06. Cosmetic change for labeling of figure showing ACFs of reconstruction for calibration
  #   and earlier years
  
  #   (methMSR=3)
  library(car)
  library(nortest)
  library(gplots)
  library(pracma)
  source(paste(code_dir,"CrossValid2.R",sep="")) # leave-m-out cross-validation  
  source(paste(code_dir,"ssValid.R",sep="")) # split-sample validation 
  source(paste(code_dir,"mannken1.R",sep="")) # time plot and trend test of reg. resids 
  source(paste(code_dir,"stemACF.R",sep="")) # stem plot of acf, with CI & annotaton 
  source(paste(code_dir,"xyCI.R",sep="")) # compute polygon (for shaded CI) from lower and upper CI 
  source(paste(code_dir,"Table1Column.R",sep="")) # write a table file with just 1 data column
  source(paste(code_dir,"TabSepTsm1.R",sep="")) # write a tab-sep file with obs, recon, 50% CI
  source(paste(code_dir,"TabSepTsm2.R",sep="")) # write a tab-sep file with model input data, 50% CI
  source(paste(code_dir,"Tsm2Scores1.R",sep="")) # time series matrix to scores of PCsI
  source(paste(code_dir,"ForwStep3.R",sep="")) # stepwise regression to get stopping step
  source(paste(code_dir,"LagkAcc.R",sep="")) # lag-k autocorrelation(s) of vector or matrix
  source(paste(code_dir,"TablePCA1.R",sep="")) # write tailored tab-sep table of PCA loadings
  source(paste(code_dir,"TableWrite1.R",sep="")) # write multi-column table
  source(paste(code_dir,"EffectSS.R",sep="")) # effective sample size (adjusted for autocorrelation)
  source(paste(code_dir,"TabSepTsm3.R",sep="")) # for time series output of PC scores
  source(paste(code_dir,"KnnAnalog.R",sep="")) # nearest neighbor analogs
  source(paste(code_dir,"emssgUNH.R",sep="")) # write error file to system, specified output folder
  
  #====== HARD CODE
  
  cBlue1 <- "#3399FF"; cMagenta1 <- '#FF00FF'; cGreen <- '#00FF00'
 
  flagBail<-0 # flag for bailing out of function
  #   0 = no problems
  #   1 = fatal error; message returned and program aborts
  #   2 = no abortion, but calibration period had to be modified to suit coverage of climate 
  #       series and SSR matrix; message is returned to calling program
  flagMsg<-'No problems'
  minLength1 <-130 # if length of series > minLength1 in reconstruction time series plot, line
  # without plot characters is plotted; otherwise line with symbols
  
  # 1 specified calibration period too short or inconsistent with data coverage 
  # 
  #======================= UNLOAD LIST, AND RENAME SOME VARIABLES
  
  U <-D$U ; yrU<-D$yrU  # full tsm of screened pf SSRs

  nmsU <- D$nmsU # ids of screened SSRs
  jScreened <- D$jScreened # pointer from screened SSRs to user-database site number
  v <-D$v ; yrv<-D$yrv  # full length predictand
  yrgo1<-D$yrsC[1]; yrsp1<-D$yrsC[2] # desired start and end year of calib period of MSR model
  yrEnd = D$yrEnd; # final reconstruction to be truncated at this year
  nNeg<-D$nNeg; nPos<-D$nPos # max neg and pos lags allowed in SSR modeling
  incR2<-D$incR2a # stepwise in MSR will not choose a model whose increase in adjusted 
  #   R-squared from the previous step is less than incR2
  kstop <-D$kstop # stopping rule for stepwise
  N1 <- D$NcMin # mimimum acceptable number of years for calibration of MSR model
  PCoption <- D$PCoption # how pool of potential predictors to be filled:
  #  1 Use will specify to use most important K PCs after viewing scree plot
  # <2> The m<fN PCs most highly correlated with predictand, where N is number of 
  #   observations in calibration period of MSR model, and f is factor (see next)
  f <- D$f # If PCoption=2, f sets the maximum possible number of predictors that will
  #   be allowed in the pool of potential predictors for the MSR model. <f=0.10>
  alphaR <- D$alphaR # for analog method, PCs retained only if correlated with y at this alpha level
  ScreenAnalogPCs <- D$ScreenAnalogPCs # whether or not to screen the analog PCs with correlation
   PCApredictors <- D$PCApredictors # TorF for using PCs of SSRs as the predictors
  methMSR <- D$methMSR # method of MSR: 1=SLR1, 2=MLR1-PCA or MLR1-noPCA, 3= Analog
  PdfDescribe <- D$PdfDescribe # string referring to pdf file that describes recon. method
  nPCsKeep <- D$nPCsKeep # user-specified number of PCs to include in pool of potential predictors
  kHowPCA <- D$kHowPCA # 1= on correlation mtx, 2=on covariance mtx
  NextFigNumber <-D$NextFigNumber # start naming figures as Figure0?.png, where ? is NextFigNumber
  outputDir <-D$outputDir # outut to be written to this system folder
  # Note that outputDir is also defined in the global environment. So, I think that I would
  # not actually need to pass outputDir as argument to functions
  HydroName <- D$Text[3] # name of hydrologic variable (e.g., "Runoff")
  HydroLabel <- D$Text[1] # label of hydro variable, for plots (e.g., "RO)
  HydroUnits <- D$Text[2] # units of hydro varialbe with parens, (e.g., "(mm
  HydroUnits2 <- substr(HydroUnits,2,(nchar(HydroUnits)-1)) # units w/o pare s
  HydroSeason <- D$Text[4] # season of hydro variable (e.g., ""12-month season ending in month 9" )
  strTailPart1 <- PdfDescribe
  if (methMSR==3){
    RecMethod<- paste('Method: Analog years of observed',HydroLabel,'from PC scores of SSRs')
  } else {
    if (PCApredictors) {
      RecMethod<- paste('Method: stepwise regression of observed',HydroLabel,'on PC scores of SSRs')
    } else{
      RecMethod<- paste('Method: stepwise regression of observed',HydroLabel,'on SSRs')
    }
  }
  rm(D)
  # all of the above are numeric or integer, except that yru is 1-col matrix
  
  #============= GET CALIBRATION DATA AND CHECK THAT CALIBRATION PERIOD LONG ENOUGH
  
  # Compute longest possible calibration period given the overlap of U and v
  yrgo2 <- max(yrU[1],yrv[1]) # earliest overlap year of U and v
  yrsp2 <- min(yrU[dim(U)[1]],yrv[length(v)]) # earliest overlap year of U and v
  
  # Return error message and fail if overlap of screened SSRs (U) and available hydro series v <N1 nears
  if ((yrsp2 -yrgo2 +1)< N1){
    emmsg <- paste('Available overlap of matrix of screened SSRs with target hydro series less than',
                   ' the required minimum of ',as.character(N1),' years.',sep='')
    ResTemp<-emssgUNH(emssg,outputDir)
    stop(emssg)
  }
  # If specified calibration start or end is NA, replace with earliest possible start,
  # latest possible end.
  if (is.na(yrgo1)){
    yrgo1<-yrgo2  # start calibration as early as possible if yrgo1 is NA
  }
  if (is.na(yrsp1)){
    yrsp1<-yrsp2 # end calibration as late as possible if yrsp1 is NA
  }
 
  #---- If calibration period specified on TRISH screen 2 not possible for MSR model  not 
  # possible for MSR model because of insufficient time time coverage by matrix U of
  # screened SSR and vector of climate series v, force the calibration period to the full
  # available overlap of the matrix and vector.
  L1 <- yrgo1 < yrgo2
  if (L1){
    yrgo1 <- yrgo2
  }
  L2 <- yrsp1 > yrsp2
  if (L2){
    yrsp1 <- yrsp2
  }
  L <- L1 || L2
  
  # Following commented out because decided to just truncate the calibration period to years supported
  # by the data. User will have to live with it. 
  # if (L){
  #   flagBail<-2
  #   flagMsg<-paste('RecMLR1 message: MSR calibration period forced to ',as.character(yrgo1),'-',
  #                  as.character(yrsp1), ' in response to available time covrage of vector of ',
  #                  ' hydro data and matrix of screened SSRs',sep='')
  # }

  U0 <- U # save the original SSRs for debugging
  
  #========================= OPTIONAL CONVERSION OF SSRs to PCs
  
  if (methMSR==3 | isTRUE(PCApredictors)){
    D <- rep(NA,4) # make empty list; "D" otherwise ia some kind R function
    DinPCA <- list("X"=U,"yrX"=yrU,"nmsX"=nmsU,"khow"=kHowPCA)
    ResPCA <- Tsm2Scores1(DinPCA)
    
    # Loadings
    LoadPC <- ResPCA$Loadings
    row.names(LoadPC)<-nmsU # assign chronology ids as row names
    namesX <- colnames(LoadPC)
    
    # Replace SSRs of U with the PC scores
    U <- ResPCA$Scores
    
    #====================== FIGURE 1x1 SCREE PLOT
    #
    # Figure files are numbered within this function with first figure as figure01.png. where
    # In general, figure files are named Figure?? is a number built from NextFigNumber+jFigAdd. jFigAdd 
    # will start at 0 but increment for later figures
    jFigAdd <- 0
    FigNumber <- NextFigNumber+jFigAdd   # for naming this png
    if (FigNumber<10){
      fileOut <- paste(outputDir,'Figure0', as.character(FigNumber),'-PCA1','.png',sep="")
    } else {
      fileOut <- paste(outputDir,'Figure', as.character(FigNumber),'-PCA1','.png',sep="")
    }

    # Blocks for plot
    nPC <- length(namesX)
    j <- 1:nPC # for x axis, PC#
    
    # Pct variance annotation for up to first 7 PCs
    if (nPC<=7){
      nListPct <- nPC
      str1 <- paste('Pctg variance of the ',as.character(nPC),' SSRs explained by PCs')
    } else {
      nListPct <-7 # will list variance explained for only first 7 Pcs
      str1 <- paste('Pctg variance of the ',as.character(nPC),' SSRs explained by first 7 PCs')
    }
    str1 <-paste(str1,'\nPC   %Var  Cum%\n')
    B <- matrix(nrow=nListPct,ncol=2)
    B[,1] <- ResPCA$PctVar[1:nListPct]; B[,2] <- ResPCA$CumPctVar[1:nListPct] # rev 2023-02-09
    
    for (n in 1:nListPct){
      str1a  <- paste(as.character(n),'  ',sprintf('%5.0f %5.0f\n',B[n,1],B[n,2]))
      str1 <- paste(str1,str1a)
    }
    rm(str1a,B,nListPct)
    
    # Plot figure
    png(filename=fileOut, width = 960, height = 480)
    par(mar=c(5,5,5,1),cex.main=1.4,cex.axis=1.5, cex.lab=1.5, cex.main=1.5)
    plot(j,ResPCA$EigValues,type="b",pch=1,col="blue",
         xlab="PC Number",ylab="Eigenvalue",
         main=paste('Scree Plot of Eigenvalues',
                    '\n(dashed red line at mean)'))
    abline(h=mean(ResPCA$EigValues),col='red',lty=2)
    text(nPC/2,max(ResPCA$EigValues),str1,adj=c(0,1),cex=1.3)
    dev.off()
    
    
    
    #====================== FIGURE 1x1 HEAT MAP OF LOADING OF PCS ON SSRS
    #
    jFigAdd=jFigAdd+1
    FigNumber <- NextFigNumber+jFigAdd   # for naming this png
    if (FigNumber<10){
      fileOut <- paste(outputDir,'Figure0', as.character(FigNumber),'-PCA2','.png',sep="")
    } else {
      fileOut <- paste(outputDir,'Figure', as.character(FigNumber),'-PCA2','.png',sep="")
    }
    png(filename=fileOut, width = 960, height = 480)
    par(mar=c(5,6,5,2),cex.main=1.4,cex.axis=1.5, cex.lab=1.3, cex.main=1.5)
    
    heatmap.2(LoadPC,Rowv = NA, Colv = NA,trace='none',key=TRUE,
              dendrogram='none',
              margins=c(5,8),
              lhei=c(1.5,5),lwid=c(1,5),
              keysize=1.0,key.title=NA,key.ylab=NA,key.xlab='Loading',
              key.ytickfun = NULL,density.info="none",
              main='PC Loadings on SSRs')
    
    dev.off()
    rm(str1)
  } else {
    # mode np-PCA
    namesX <- nmsU
  }
  
  #========================= PULL CALIBRATION DATA 
  L <- yrU>=yrgo1 & yrU<=yrsp1
  U1 <- U[L,]; yrU1<- yrU[L]
  L <- yrv>=yrgo1 & yrv<=yrsp1
  v1 <- v[L]; yrv1<- yrv[L]
  ncalib <- length(v1) # length of MSR calibration period
  # Status. U1 is matrix; yrU1, v1 and yrv are vector.
  
  
  #========================= CORRELATION OF HYDRO SERIES WITH SSRs OR THEIR PCS; LAG-1 
  # AUTOCORRELATIONS OF SAME
  
  #--- Check on desired threshold correlation
  alphaRR <- c(0.10, 0.05, 0.01) # acceptable alpha levels (2 tailed)
  Siggys <- c('90%','95%','99%')
  alphaTemp <- 1-alphaRR/2
  ThreshsR <- qnorm(alphaTemp) 
  L <- alphaR == alphaRR
  if (!any(L)){
    emssg <- 'Invalid setting for critical alpha for correlation screening; must be one of {0.10,0.05 ,0.01}'
    ResTemp<-emssgUNH(emssg,outputDir)
    stop(emssg)
  }
  ThreshSiggy <- ThreshsR[L] # this over sqrt(N-2) gives 2=tailed confidence interval
  Siggy <- Siggys[L]
  rm(L,alphaRR,Siggys,alphaTemp,ThreshsR)
  rThresh <- ThreshSiggy/sqrt(length(v1)-2) # threshold level of correlation
  
  #--- Correlations
  rv1U1 <- cor(v1,U1)
  rCI = rThresh # rename for convenience; threshold correlation for confidence interval
  
  #--- Lag-1 autocorrelations
  ResLag1 <- LagkAcc(v1,1)
  r1v1 <- ResLag1$rk
  ResLag1 <- LagkAcc(U1,1)
  r1U1 <- ResLag1$rk
  rm(ResLag1)
  
  
  
  #====================== FIGURE 1x1 BAR CHART OF R OF HYDRO SERIES WITH SSRs or THEIR PCS
  
  tit1 <- c(paste('Correlation with ',Siggy,' CI (blue) and autocorrelation (magenta),',
                  as.character(yrgo1),'-',as.character(yrsp1)),
            paste('\n(magenta horizontal line is r(1) of ',HydroLabel,')',sep=""))
  if (methMSR==3 | PCApredictors){
    xlab1 <- 'PC #'
    legText <- c(paste('r of',HydroLabel,'with PCs'),
                 'r(1) of PCs')
    jFigAdd <- jFigAdd+1
  } else {
    xlab1 <- 'SSR # (in order as in matrix of screened SSRs)'
    legText <- c(paste('r of',HydroLabel,'with SSRs'),
                 'r(1) of SSR')
    jFigAdd <- 0
  } 
  
  FigNumber <- NextFigNumber+jFigAdd   # for naming this png
  if (FigNumber<10){
    fileOut <- paste(outputDir,'Figure0', as.character(FigNumber),'-MSRcalibration1','.png',sep="")
  } else {
    fileOut <- paste(outputDir,'Figure', as.character(FigNumber),'-MSRcalibration1','.png',sep="")
  }
  png(filename=fileOut, width = 960, height = 480)
  par(mar=c(5,6,5,2),cex.main=1.0,cex.axis=1.5, cex.lab=1.3, cex.main=1.5)
  h <- t(cbind(t(rv1U1),r1U1))
  
  barplot(h,ylim=c(-1,1),beside=TRUE,main=tit1,xlab=xlab1,ylab='r',
          width=c(1.0,0.5),names.arg=as.character(1:dim(U1)[2]),
          col=c(cBlue1,cMagenta1))
  legend(x=dim(U1)[2],y=1,xjust=1,yjust=1,legend=legText,fill=c(cBlue1,cMagenta1),
         cex=1.3)
  abline(h=0,col='black')
  abline(h=c(rCI, -rCI), col=cBlue1, lty=2, lwd=4)
  abline(h=r1v1,col=cMagenta1,lty=1,lwd=2)
  
  dev.off()
 # axis(1, at =(1:length(rv1U1)),labels = nmsU)

  #========== POSSIBLE REDUCTION OF POOL OF POTENTIAL PREDICTORS OR ANALOG SET OF PC'S
  # 
  # Input arg f gives a fraction, which if multiplied times the number of observations in
  # the calibration period of the MSR model, puts an upper limit on the allowable number
  # of variables in the pool of potential predictors. 
  #
  # The pool must be less than fN, where N is the length of the calibration period. How this
  # constraint is applied depends on if the MSR is done using the SSRs or their PCs as 
  # predictors. If using the SSRs, the pool is reduced so that it includes only the
  # m<fN SSRs most highly correlated with y. If using PCs of the SSRs, and if PCoption=2, the
  # same correlation screening is done, but on the PCs of the SSRs rather than on the SSRs.
  # If the using PCs and PCopton=1 (user specified to use a certain number of "nost important"
  # PCs), execution stops with an error message if the the specified number of PCs is too
  # large.
  #
  # This section also handles correlation screening for the analog reconstruction method
  # (methMSR=3). For that, screening is done according to the specified alphaR alpha-level of 
  # required signficance of correlation of PCs with the predictand.
   
  nmax <- floor(f*length(v1)) # maximum allowable size of pool of potential predictors
  if ((methMSR==3 | PCApredictors) && PCoption==1){
    if (nPCsKeep > nmax) {
      emssg <- paste('You instructed to retain more than',sprintf('%g',nmax),'PCs; too many!')
      ResTemp<-emssgUNH(emssg,outputDir)
      stop(emssg)
    } else {
      # Assign as pool the first nPCsKeep PCs of the SSRs
      jU2toU1 <- 1:nPCsKeep # pointer to cols of U1
      U2 <- U1[,(1:nPCsKeep)] # pool of potential predictors
      yrU2 <- yrU1;
      namesU2 <- namesX[jU2toU1]
    }
  } else {
    # Predictors could be SSRs or their PCs, and correlation with y to be used, if needed,
    # to keep pool small.
    # Correlate columns of U1 with y and compute absolute correlation for retaining predictors
    # Threshold is two-tailed test, no adjustment for autcorrelation; input alphaR must be one of
    #  [0.10 0.05 0.01]
    rThese <- cor(v1,U1) # vector of correlations of SSRs or PCs with predictand
    if (methMSR==3){
      if (ScreenAnalogPCs){
        # Analog method; retain only the PCs with correlation significance greater than threshold
        Lkeep <- abs(rThese)>rThresh 
        if (!any(Lkeep)){
          # ERROR MESSAGE
          emssg <- paste('No PCs passed screening for correlation with',HydroLabel,
                         'at alpha=',sprintf('%g',alphaR))
          ResTemp<-emssgUNH(emssg,outputDir)
          stop(emssg)
        } else {
          U2 <- U1[,Lkeep]; yrU2 <- yrU1
          jU2toU1 <- which(Lkeep) # pointer of U2 back to U1, cols
          namesU2 <- namesX[Lkeep]
        }
      } else {
        # no screening of analog PCs
        U2 <- U1; yrU2 <- yrU1; namesU2 <- namesX; jU2toU1 <- 1:dim(U1)[2]
      }
    } else {
      # Not analog method; will screen to avoid too large a pool
      if (dim(U1)[2]>nmax) {   #if pool too large
        # --- Correlation screening
        rSort <- sort(abs(rThese),decreasing=TRUE,index.return=TRUE)
        rCut <- rSort$x[nmax] # PCs or SSRs with smaller absolute r will be dropped
        Lkeep <- abs(rThese)>=rCut
        U2 <- U1[,Lkeep]; yrU2 <- yrU1
        jU2toU1 <- which(Lkeep)
        namesU2 <- namesX[Lkeep]
      } else {
        # No problem with size of pool; use them all
        U2 <- U1;  yrU2 <- yrU1  # pool of potential predictors
        jU2toU1 <- 1:dim(U1)[2] # pointer of U2 back to U1, cols
        namesU2 <- namesX
      }
    }
  }
  
  # In case U2 has just 1 column, make sure U2 is a matrix
  U2 <- as.matrix(U2)


  #================  EUCLIDEAN DISTANCE 
  #
  # If methMSR==3, find instrumental-period first-k nearest neighbors in yrv1 for
  # every year of to every year in U. Compute Euclidean distances for vectors of
  # PC scores. Build 2-dim time series matrices of the years and of the "analog y"
  #
  # Also build time series of "nearest-observed" y and compute statistics of 
  # accuracy if such a series were used as the reconstruction
  
  
  if (methMSR == 3){
    kNN <- 2 # Want kNN nearest neighbors from yrv1, v1 for every year of U
    ResKNN <- KnnAnalog(U[,jU2toU1],yrU,v1,yrv1,kNN)
    # Among outputs is ResKNN$Recon, with "reconstruction" in col 2 and the corresponding
    # analog years in col 3. Col 4 tells if recon is nearest or 2nd nearest neighbor
    # Also has "NearestObserved," which is a null reconstruction using as the 
    # reconstruction in any year of yrv1 the nearest v to the value of v in yrv1 
    # 
    
    #--- ASSESSMENT OF ACCURACY OF ANALOG AND "NEAREST-OBSERVED" RECONS
    L <- ResKNN$Recon[,4] == 2
    # Analog predictions
    Fits <- ResKNN$Recon[L,2] 
    yrFits <- yrv1
    # Nearest-observed predictions
    FitsNO <- ResKNN$NearestObserved[,3] 
    yrFitsNO <- yrv1
    
    #--- SOS terms-- Analog
    Efits <- v1-Fits; # error, or obs minus recon
    SSE  <- sum(Efits*Efits) # sum of squares of errors
    SSv1 <-  sum((v1-mean(v1)) * (v1-mean(v1))) # sos of departures of v1 from its mean
    REtemp <- 1-SSE/SSv1
    rtemp <- cor(v1,Fits)
    
    #--- SOS terms-- Nearest-observed
    EfitsNO <- v1-FitsNO; # error, or obs minus nearest-observed
    SSENO  <- sum(EfitsNO*EfitsNO) # sum of squares of errors
    REtempNO <- 1-SSENO/SSv1
    rtempNO <- cor(v1,FitsNO)
    
    #--- STORE
    Analog <- list(Fits=Fits,yrFits=yrFits,Resids=Efits,Correl=rtemp,RE=REtemp,
                   DataKnn=ResKNN$DataAnalog,YearsKnn=ResKNN$YearsAnalog,
                   DistanceKnn=ResKNN$Distance,ReconNN12=ResKNN$Recon,
                   RMSE=sqrt(mean(SSE)),
                   FitsNO=FitsNO,yrFitsNO=yrFitsNO,ResidsNO=EfitsNO,CorrelNO=rtempNO,RENO=REtempNO,
                   RMSENO = sqrt(mean(SSENO)))
    rm(rtemp,REtemp,SSv1,SSE,Efits,Fits,yrFits,L,ResKNN,
       rtempNO,REtempNO,SSENO,EfitsNO,FitsNO)
  }


  if (methMSR != 3){
    #========================= STEPWISE REGRESSION
    # 
    # To select forward stepwise from pool of potential predctors and do an initial
    # fit of model
    #
    # Status. 
    # U2, yrU2: pool of potential predictors
    # iU1: pointer (cols) of U2 back to U1 and to long matrix U
    # v1, yrv1: predictand vector
    ResMLR1a <- ForwStep3(U2,namesU2, v1,kstop,nNeg,nPos, incR2a)
    inmodelU3 <- ResMLR1a$ColsInModel
    npred <- length(inmodelU3) # number of predictors in MSR model
    U3 <- U2[,inmodelU3]; yrU3 <-yrU2
    namesU3 <- namesU2[inmodelU3] # names of predictors in model
    jU3toU1 <- jU2toU1[inmodelU3] # pointer (col) back to full U and to original calibration U1
    Fits <- ResMLR1a$Model$fitted.values # calib period predictions
    ModelSummary <- summary(ResMLR1a$Model)
    RMSEc <- ModelSummary$sigma

    #========================= RE-FIT REGRESSION WITH LM FUNCTION (for debugging only)
    # 
    # Function lm gives comprehensive statistics. Can uncomment the three lines, and check that the
    # statistics returned exactly match those in ResMLR1a
    # M <- lm(v1~U3) # model object
    # M1<-summary(M)
    
    
    #========================= STORE CALIBRATION STATISTICS
    
    # Significance of overall F
    # M1$fstatistic has F, dfnum, dfdenom in 1-3
    pF <-ResMLR1a$Fpvalue
    
    
    OutputCal<-list('flag'=flagBail,'Msg'=flagMsg,'lmModel'=ResMLR1a$Model,'yearGoCal'=yrgo1,'yearSpCal'=yrsp1,
                    'coefficients'=ResMLR1a$Coefficients,'Rsquared'=ResMLR1a$Rsquared,'F'=ResMLR1a$Foverall,'pF'=pF,
                    'RsquaredAdj'=ResMLR1a$RsquaredAdj,'RMSEc'=RMSEc)
    
  } else {
    # Analog method
  }
  
  #====================== FIGURE CALIBRATION2
  #

  #--- Build some strings for use in plots
  strCalPd <- paste(' Calibration period: ',sprintf('%d',yrgo1),'-',sprintf('%d',yrsp1),
                    ' (N=', sprintf('%d',length(yrv1)),' yr)',sep='')
  
  #--- Build Tit1 figure png filename

  if (methMSR==3){
    Fits <- Analog$Fits
    strPart1 <- '\nNumber of: screened SSRs / PCs for analogs'
    strPart2 <- paste('\nPCA on years ',as.character(yrU[1]),'-',as.character(yrU[length(yrU)]),
                      paste('\nAnalogs from years ',as.character(yrv1[1]),'-',as.character(yrv1[length(v1)]),sep=''),sep='')
    txtCommon <- c('Common','common')
  } else {
    txtCommon <- c('Calibration','calbration')
    if (PCApredictors){
      strPart1 <- '\n   # screened PCs in pool / # PCs in final model'
      strPart2 <- ' (see Table 8 for model coefficients)'
    } else {
      strPart1 <- '\n   # screened SSRs in pool / # in final model'
      strPart2 <- ' (see Table 6 for model coefficients)'
    }
  }
  jFigAdd<-jFigAdd+1
  FigNumber <- NextFigNumber+jFigAdd   # for naming this png
  if (FigNumber<10){
    fileOut <- paste(outputDir,'Figure0', as.character(FigNumber),'-MSRcalibration2.png',sep="")
  } else {
    fileOut <- paste(outputDir,'Figure', as.character(FigNumber),'-MSRcalibration2.png',sep="")
  }
  # Build figure
  png(filename=fileOut, width = 960, height = 480)
  layout.matrix <- matrix(c(1,2,3,3),nrow=2,byrow=TRUE)
  layout(layout.matrix,heights=c(2,2),widths=c(1,2))

  # Scatter
  par(mar=c(4,5,2,3))
  r<-cor(Fits,v1)
  strTit <- paste('Recon vs Obs ',HydroName,', r=',as.character(round(r,digits=2)))
  plot(v1,Fits,ylab=paste('Recon',HydroLabel,HydroUnits),
       xlab=paste('Obs',HydroLabel,HydroUnits),main=strTit)
  abline(lm(Fits~v1),col="red")

    
  # Stats table
  par(mar=c(0,0,0,0))
  plot(1,1,xaxt="n",yaxt="n",bty="n",pch="",ylab="",xlab="", main="", sub="",
       xlim=c(0,1),ylim=c(0,1))
  if (methMSR==3){
    strText <-paste(RecMethod,strPart1,
                    '\n      ',sprintf('%g',dim(U1)[2]),
                    ' / ',sprintf('%g',dim(U2)[2]),strPart2,sep='')
    text(x=0.05,y=0.95,'Model Statistics',adj=c(0,1),cex=1.5,font=2)
    text(x=0.05,y=0.8,strText,adj=c(0,1),cex=1.2)
  } else {
    if (npred==1){
      U3 <- as.matrix(U3) # if happen to be only 1 predictor in final model, still want U3 as matrix
    }
    strText <-paste(RecMethod,strPart1,
                    '\n      ',sprintf('%g',dim(U2)[2]),
                    ' / ',sprintf('%g',dim(U3)[2]),strPart2,
                    '\n',strCalPd,
                    '\n   R-squared =',sprintf('%5.2f',ResMLR1a$Rsquared),
                    '\n   F=',sprintf('%.5g',ResMLR1a$Foverall),'(p=',sprintf('%.4g',pF),')',
                    '\n   RMSE=',sprintf('%g',RMSEc),HydroUnits2,' (equivalent to std error of the estimate)')
    if (methMSR==3){
      text(x=0.05,y=0.95,'Analog Model Statistics',adj=c(0,1),cex=1.5,font=2)
    } else {
      text(x=0.05,y=0.95,'Calibration Statistics',adj=c(0,1),cex=1.5,font=2)
    }
    text(x=0.05,y=0.8,strText,adj=c(0,1),cex=1.4)
  }

  # Time plots, obs and rec
  zx <- c(yrgo1,yrsp1)
  zy <- c(mean(v1),mean(v1))
  
  ylims <- c(min(v1,Fits),max(v1,Fits))
  ylimsInc <- 0.05 * diff(ylims)
  ylims <- c(ylims[1]-ylimsInc,ylims[2]+ylimsInc)
  rm (ylimsInc)
  
  par(mar=c(4,5,3,2),cex.axis=1.5, cex.lab=1.5, cex.main=1.5)
  plot(yrv1,v1,type="b",pch=1,col="blue",xlim=c(yrgo1,yrsp1),ylim=ylims,
       ylab=paste(HydroLabel,HydroUnits),xlab="Year",
       main=paste(txtCommon[1],'Period Observed (blue) and Reconstructed (red)',HydroName))
  lines(yrv1,Fits,type="b",pch=2,col="red")
  lines(zx,zy)
  h<-par("usr"); yoffset<- (h[4]-h[3])/100; ytop <- h[4]-yoffset
  # legend(yrgo1+1,ytop,legend=c("Obs", "Recon"),
  #        col=c("blue", "red"), lty=1, cex=1.2)
  
  dev.off()
  rm (strPart1,strPart2,strText,strTit,r,h)
  
  
  #====================== FIGURE:  CALIBRATION3
  #
  # CW from LL: histogram of yhat;  hist of obs y; acf of obs y, acf of yhat
  
  #--- same xlims for all histograms
  xlo = min(c(min(v1),min(Fits)))
  xhi = max(c(max(v1),max(Fits)))
  xinc <- 0.05*(xhi-xlo)
  xlims <- c(xlo-xinc,xhi+xinc)
  rm(xlo,xhi,xinc)
  #xlims <-c((c(min(v1),min(Fits))), (c(max(v1),max(Fits))))
  
  
  #--- Build figure png filename
  jFigAdd <- jFigAdd+1
  FigNumber <- NextFigNumber+jFigAdd   # for naming this png
  if (FigNumber<10){
    fileOut <- paste(outputDir,'Figure0', as.character(FigNumber),'-MSRcalibration3.png',sep="")
  } else {
    fileOut <- paste(outputDir,'Figure', as.character(FigNumber),'-MSRcalibration3.png',sep="")
  }
  # Build figure
  png(filename=fileOut, width = 960, height = 480)
  # Create the layout
  nf <- layout( matrix(c(1,2,3,4), ncol=2) )
  
  nbin <- floor(5*log10(length(v1))) # Panofsky rule of thumb for number of bins
  xtweak1 <- 0.00001* (max(v1)-min(v1));
  xtweak2 <- 0.00001* (max(Fits)-min(Fits));
  brks1 = seq(min(v1)-xtweak1,max(v1)+xtweak1,length.out=(nbin+1))
  brks2 = seq(min(Fits)-xtweak2,max(Fits)+xtweak2,length.out=(nbin+1))
  xlab1 <- paste(HydroLabel,HydroUnits)
  MaxLag <- floor(length(v1)/4)
  
  par(mar=c(5,6,2,2),cex.axis=1.5, cex.lab=1.5, cex.main=1.5)
  Tit1 <- paste('Histogram, Observed',HydroLabel,'(N=',sprintf('%d',length(v1)),'yr)')
  hist(v1,main=Tit1,breaks=brks1,xlim=xlims,xlab=xlab1)
  
  par(mar=c(5,6,2,2))
  Tit1 <- paste('Histogram, Reconstructed',HydroLabel,'(N=',as.character(length(v1)),'yr)')
  hist(Fits,main=Tit1,breaks=brks2,xlim=xlims,xlab=xlab1)
  
  par(mar=c(5,5,4,2))
  Tit1 <- paste('ACF, Observed',HydroLabel,',with 95% CI')
  acf(v1,lag.max=MaxLag,type='correlation',main=Tit1)
  
  par(mar=c(5,5,4,2))
  Tit1 <- paste('ACF, Reconstructed',HydroLabel,', with 95% CI')
  acf(Fits,lag.max=MaxLag,type='correlation',main=Tit1)
  dev.off()

  
  #========================= ANALYSIS OF RESIDUALS (local 3,4,5)
  #
  # Will need the regression residuasl. First figure (1x2) will be histogram and
  # scatter of residuals on predicted values. Second figure (1x1) will be time plot of
  # the residuals with a fitted (non-parametric fit) trend line and annotated result
  # of Mann-Kendall trend test. The significance will be adjusted as needed for autocorrelation
  # of the residuals over and above that in a linear trend. The third figure (1x1) will be
  # the acf of residuals, with annotated DW test results
  
  
  #============= FIGURE: ANALYSI OF RESIDUALS: HISTORGRAM & CONSTANCY OF VARIANCE 
  #
  # 1st of 3 analysis of residuals plots
  
  if (methMSR==3){
    Resids <- Analog$Resids
  } else {
    Resids <- ResMLR1a$Model$residuals
  }
  
  # Lilliefors test of normality of residuals
  hLillie <- lillie.test(Resids);
  if (methMSR==3){
    Tit1 <- paste('Residuals, E, of analog reconstruction of',HydroName,
                  '\n  (p=',sprintf('%.2g',hLillie$p.value),'for H0 that E normal, from Lilliefors Test)')
  } else {
    Tit1 <- paste('Residuals, E, of regression of',HydroName,'on tree rings',
                  '\n  (p=',sprintf('%.2g',hLillie$p.value),'for H0 that E normal, from Lilliefors Test)')
  }
  
  
  # Breusch-Pagan test for heterogeneity of regression residuals
  if (methMSR==3){
    BP <-NA
    Tit2 <- paste('Scatter of Residuals on Fitted Values','\n')
  } else {
    BP <- ncvTest(ResMLR1a$Model)
    Tit2 <- paste('Scatter of Residuals on Fitted Values',
                  '\n  (p=',sprintf('%.2g',BP$p),' for H0 that E homoscedastic)',
                  '\n(from Breusch-Pagan Test)')
  }
 
  # Buld filename for plot
  jFigAdd <- jFigAdd+1
  FigNumber <- NextFigNumber+jFigAdd   # for naming this png
  if (FigNumber<10){
    fileOut <- paste(outputDir,'Figure0', as.character(FigNumber),'-AnalysisResiduals1.png',sep="")
  } else {
    fileOut <- paste(outputDir,'Figure', as.character(FigNumber),'-AnalysisResiduals1.png',sep="")
  }
  
  # Figure size and margins
  png(filename=fileOut, width = 960, height = 480)
  layout.matrix <- matrix(c(1,2), nrow = 1, ncol = 2)
  layout(layout.matrix,heights=2,widths=c(1,1)) 
  
  # Left plot

  par(mar = c(5.1, 4.5, 5.1, 2.1),cex.axis=1.1, cex.lab=1.5, cex.main=1.3)
  hist(Resids,xlab=paste('Residual',HydroUnits),ylab='Frequency',main=Tit1)
   
  # # right plot
  par(mar = c(5.1, 4.1, 6.1, 2.1),cex.axis=1.1, cex.lab=1.5, cex.main=1.3)
  plot(Fits,Resids,xlab=paste('Predicted',HydroUnits),
       ylab=paste('Residual',HydroUnits),
       main=Tit2)
  abline(h=0,lty=2,col='#808080') # dash gray
  dev.off()

  
  #--- TIME PLOT OF REGRESSION RESIDUALS, WITH MANN-KENDAL TREND TEST (1X1) 
  
  if (methMSR==3){
    RegWord <- 'Model'
  } else {
    RegWord <- 'Regression'
  }
  
  jFigAdd <- jFigAdd+1
  FigNumber <- NextFigNumber +jFigAdd 
  
  # Prepare input for mannken1
  X  <- cbind(yrv1,Resids) # matrix with year and regression or analog residuals

  kopt<- c(2,1) # want plot; want adjustment of significance of Mann-Kendall statistic
  # for autocorrelation if warranted
  kplot <-2 # for TRISH, the time plot of residuals, with annotated MK test results
  # and non-parametric-fit straight line fit to trend
  ylabTemp1 <- paste('Residual',HydroUnits)
  ylabTemp2 <- paste('Detrended Residual',HydroUnits)
  textPlot <- c('Regression Residuals with Nonparametric-Fit Trend Line,','Year',ylabTemp1,
                ylabTemp2)
  Din <- list(X=X,kopt=kopt,kplot=kplot,NextFigNumber=FigNumber,textPlot=textPlot,outputDir=outputDir)
  
  # mannken1 to get statistics and plot
  ResMK <- mannken1(Din)
  rm(Din)
  
  
  #--- ACF OF REGRESSION RESIDUALS, INCLUDING 95% Ci 
  
  jFigAdd <- jFigAdd+1
  FigNumber <- NextFigNumber +jFigAdd 

    # acf and its 95% CI
  lagsPlot <- floor(min(ncalib/4,20))
  acfMy <- acf(Resids, lag.max=lagsPlot, type = "correlation", 
               plot = FALSE)
  k <- acfMy$lag # lags
  w <- acfMy$acf # acf
  
  # DW statistic of the regression residuals
  if (methMSR==3){
    strDW <- 'Durbin-Watson not computed (inappropriate for analog method)'
  } else {
    DW <-  durbinWatsonTest(ResMLR1a$Model)
    strDW <- paste('p=',sprintf('%g',DW$p),
                   ': Durbin-Watson test (2-sided) of H0 that population lag-1 autocorrelation is zero',
                   sep='')
  }
  # Text for plot
  Tit1 <- paste('ACF of Residuals with 95% CI (N=',sprintf('%g',ncalib),' yr)')
  textPlot <- c(Tit1,'Lag(yr)','r',strDW,'-AnalysisResiduals3')
  
  # Store inputs required by stemACF in a list (see opening comments there)
  Din <- list(x=k,y=w,nsize=ncalib,kAlpha=1,FigNumber=FigNumber,
              outputDir=outputDir,linecol1='#0022CC',linecol2='#696969',
              linecol3='#E60000',textPlot=textPlot)
  ResNull <- stemACF(Din)
  
  #========================= VALIDATE AND STORE STATISTICS
  
  #ResCV <- CrossValid2(u1, v1, nNeg,nPos) # cross-validation
  
  # Split-sample validation -- not if analog method
  if (methMSR==3){
    # Analog method: no split-sample validation.
    # For analog method, the model residuals can be considered validation residuals, because
    # no model has been fit and no tuning to improve estimates for the overlap period
    # wit climate. The stored Analog$RMSE can therefore be used in place of RMSEcv.
    # For analog method there is also no "RE" statistic. But a pseudo-RE is computed as
    # 1-SSE1/SSE2, where SSE is the sum-of-squares of residuals for the overlap period
    # with y, and SSE2 is the sum-of-squares of departures of observed y from its mean.
    # This RE is therefore a skill statistic comparing errors for the reconstruction with errors
    # of a null reconstruction consisting of the observed mean for each year. 
    OutputVal<-list('RE'=Analog$RE,'RMSE'=Analog$RMSE,'Correlation'=Analog$Correl)
  } else {
    
    iAstop <- ceiling(length(v1)/2) # end row index in v1 of first half of data, assumed longer than
    # longer of the two halves if length of v1 odd
    iBgo <- iAstop+1 # start row of second half
    iA <- 1:iAstop # row indices of first half of full calib period
    iB <- iBgo:length(v1) # ... of second half
    
    #--- Calibrate on early, validate on late, then reverse
    ical<-iA; ival<-iB
    i1 <- 1:dim(U3)[2]; # all columns of U3 are the predictors in the final model
    ResSS1=ssValid(v1,U3,ical,ival,i1);
    REa1<-ResSS1$RE # RE for calib on early, valid on late
    ical<-iB; ival<-iA
    ResSS2=ssValid(v1,U3,ical,ival,i1);
    REb1<-ResSS1$RE # RE for calib on late, valid on early
    OutputVal<-list('mLeaveOutCV'=ResMLR1a$CrossValidStorage$LeftOut,
                    'REcv'=ResMLR1a$CrossValidStorage$REcv,'RMSEcv'=ResMLR1a$CrossValidStorage$RMSEcv,
                    'REcalEarlyValLate'=REa1,'REcalLateValEarly'=REb1)
  }
  
  
  #====================== FIGURE Validation: time plots of obs y, recon, y, cv predictions of y; with observed mean line
  #
  # Completely different figures is made for analog vs other methods (methMSR)
  #
  #--- Build figure png filename
  jFigAdd <- jFigAdd+1
  FigNumber <- NextFigNumber+jFigAdd   # for naming this png
  if (FigNumber<10){
    fileOut <- paste(outputDir,'Figure0', as.character(FigNumber),'-Validation1.png',sep="")
  } else {
    fileOut <- paste(outputDir,'Figure', as.character(FigNumber),'-Validation1.png',sep="")
  }

  ylims <- c(min(v,Fits),max(v1,Fits))
  ylims <- c(ylims[1]-0.05*diff(ylims),ylims[2]+0.05*diff(ylims))
  if (methMSR==3){
    # Recall:  v1, yrv1 is observed predictand; Fits, yrFits are analog-reconstructed
    # predictand, which is 2nd nearest neighbor;
    # Build figure, time plot at left, window at right

        # Figure size and margins
    png(filename=fileOut, width = 960, height = 480)
    layout.matrix <- matrix(c(1,2), nrow = 1, ncol = 2)
    layout(layout.matrix,heights=2,widths=c(4,1)) 
    
    par(mar=c(4,4,5,1),cex.main=1.4)
    plot(yrv1,v1,ylim=ylims,type="b",pch=1,col="blue",xlim=c(yrgo1,yrsp1),
         ylab=paste(HydroLabel,HydroUnits),xlab="Year",
         main=paste('Observed (blue) and Analog-Predicted (red)',HydroName,
        '\n(Years=',sprintf('%d',yrgo1),'-',sprintf('%d',yrsp1),
                    '; black line at observed mean; green * at nearest observed)'))
    lines(yrv1,Fits,type="b",pch=2,col="red")
    lines(yrv1,Analog$FitsNO,type="p",pch=8,col=cGreen)
    abline(h=mean(v1))
    
    #--- Stats window
    par(mar = c(0,0,0,0))
    xlims <- c(0,1); ylims <- c(0,1)
    plot(0,type='n',axes=FALSE,ann=FALSE,xlim=xlims,ylim=ylims)
    
    # Analog
    strText <-paste('\n   Analog',
                    '\n \n  RMSE=',sprintf('%g',Analog$RMSE),HydroUnits2,
                    '\n  RE=',sprintf('%.2g',Analog$RE),
                    '\n  r=',sprintf('%.2g',Analog$Correl))
    text(x=0.01,y=0.95,'Validation statistics',adj=c(0,1),cex=1.3,font=2)
    text(x=0.01,y=0.9,strText,adj=c(0,1),cex=1.2)
    
    # Nearest Observed
    strText <-paste('\n  Best Possible','\n (Nearest Observed)',
                    '\n\n  RMSE=',sprintf('%g',Analog$RMSENO),HydroUnits2,
                    '\n  RE=',sprintf('%.4g',Analog$RENO),
                    '\n  r=',sprintf('%.4g',Analog$CorrelNO))
    text(x=0.01,y=0.65,strText,adj=c(0,1),cex=1.2)
    dev.off()
  } else {
    # Build figure
    png(filename=fileOut, width = 960, height = 480)
    par(mar=c(4,4,5,1),cex.main=1.4)
    plot(yrv1,v1,ylim=ylims,type="b",pch=1,col="blue",xlim=c(yrgo1,yrsp1),
         ylab=paste(HydroLabel,HydroUnits),xlab="Year",
         main=paste('Time Plots of Observed, Reconstructed, and Cross-Validation-Predicted',HydroName,'(',sprintf('%d',yrgo1),'-',sprintf('%d',yrsp1),')',
                    '\n(black line at observed mean)'))
    lines(yrv1,Fits,type="b",pch=2,col="red")
    lines(yrv1,ResMLR1a$CrossValidStorage$CVpredictions,type="b",pch=17,col="#990099")
    lines(zx,zy)
    h<-par("usr"); yoffset<- (h[4]-h[3])/100; ytop <- h[4]-yoffset
    legend(yrgo1+1,ytop,legend=c("Obs", "Recon","CVpred"),
           col=c("blue", "red","#990099"),pch=c(1,2,17),lty=1, cex=1.2)
    dev.off()
  }
  
  #====================== FIGURE VALIDATION2
  #
  # Does not apply to analog method
  # Text strings for cross-validation
  strAnPd <- strCalPd # built for earlier plot: gives calib period and length
  if (methMSR!=3){
    jFigAdd <- jFigAdd+1
    FigNumber <- NextFigNumber+jFigAdd   # for naming this png
    if (FigNumber<10){
      fileOut <- paste(outputDir,'Figure0', as.character(FigNumber),'-Validation2.png',sep="")
    } else {
      fileOut <- paste(outputDir,'Figure', as.character(FigNumber),'-Validation2.png',sep="")
    }
    # Lilliefors test of normality of CV residuals
    hLillie <- lillie.test(ResMLR1a$CrossValidStorage$CVresiduals);
    Tit1 <- paste('Cross-Validation Residuals',
                  '\n(p=',sprintf('%.2g',hLillie$p.value),', H0: normally distributed',
                  ' [Lilliefors Test])',sep='')
    
    # Text strings for split-sample validation
    yrtemp <- yrv1[iA] # year vector, early split
    yrgoA <- yrtemp[1];  yrspA <- yrtemp[length(yrtemp)]
    yrtemp <- yrv1[iB] # year vector, late split
    yrgoB <- yrtemp[1];  yrspB <- yrtemp[length(yrtemp)]
    rm(yrtemp)
    strSplitA <- paste('   A: ',sprintf('%d',yrgoA),'-',sprintf('%d',yrspA),
                       ' (N=', sprintf('%d',length(iA)),' yr)',sep='')
    strSplitB <- paste('   B: ',sprintf('%d',yrgoB),'-',sprintf('%d',yrspB),
                       ' (N=', sprintf('%d',length(iB)),' yr)',sep='')
    

    # Figure size and margins
    png(filename=fileOut, width = 960, height = 480)
    layout.matrix <- matrix(c(1,2), nrow = 1, ncol = 2)
    layout(layout.matrix,heights=2,widths=c(1,1))
    strText <-paste('\nCross-validation (cv) method: leave-',as.character(ResMLR1a$CrossValidStorage$LeftOut),'out',
                    '\n',strAnPd,
                    '\n   RMSEcv=',sprintf('%g',ResMLR1a$CrossValidStorage$RMSEcv),HydroUnits2,
                    '\n   REcv=',sprintf('%.2g',ResMLR1a$CrossValidStorage$REcv),
                    '\n\nSplit-sample validation',
                    '\n',strSplitA,
                    '\n',strSplitB,
                    '\n    RE{A}=',sprintf('%.2g',ResSS1$RE),' (calibrated on A, validated on B)',
                    '\n    RE{B}=',sprintf('%.2g',ResSS2$RE),' (calibrated on B, validated on A)')

    # Left plot histogram
    par(mar = c(5.1, 4.5, 5.1, 2.1),cex.axis=1.1, cex.lab=1.5, cex.main=1.3)
    hist(ResMLR1a$CrossValidStorage$CVresiduals,xlab=paste('Residual',HydroUnits),ylab='Frequency',
         main=Tit1)

        # right plot, stats
    par(mar = c(0,0,0,0))
    xlims <- c(0,1); ylims <- c(0,1)
    plot(0,type='n',axes=FALSE,ann=FALSE,xlim=xlims,ylim=ylims)
    text(x=0.05,y=0.95,'Validation statistics',adj=c(0,1),cex=1.5,font=2)
    text(x=0.05,y=0.8,strText,adj=c(0,1),cex=1.2)
    #plot(1,1,pch=1,xlim=xlims,ylim=ylims)
    dev.off()
  }

  #========================= RECONSTRUCTION WITH 50% CI
  if (methMSR ==3){
    Xr <- as.matrix(U[,jU2toU1]) # matrix,long-term scores of PCs used to find analogs
    L<-complete.cases(Xr)
    Xr <- as.matrix(Xr[L,]) ; yrXr <- as.matrix(yrU[L,])
    mXr <- dim(Xr)[1]
    yrgo3 <- yrXr[1,1]; yrsp3 <- yrXr[mXr,1] # start and end year of recon
    yh <- Analog$ReconNN12[,1:2] # year and reconstructed y, full recon
    
    # Delta y for 50% CI
    xStdNorm75 <-qnorm(0.75, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
    deltaRec50 <- xStdNorm75 * Analog$RMSE
    yhLo <- yh[,2]-deltaRec50; yhHi <-yh[,2]+deltaRec50
    yh <- cbind(yh,yhLo,yhHi) # matrix with year, recon, lower 50 upper 50
  } else {
    #Xr <- as.matrix(U[,inmodelU3]) # matrix, long-term predictors; INCORRECT
    Xr <- as.matrix(U[,jU3toU1]) # matrix, long-term predictors
    L<-complete.cases(Xr)
    yrXr <- as.matrix(yrU[L,])
    # Add ones column and reconstruct
    mXr <- dim(Xr)[1]
     Xones<-matrix(1,nrow=mXr,ncol=1)
    Xr <- cbind(Xones,Xr)
    yrgo3 <- yrXr[1,1]; yrsp3 <- yrXr[mXr,1] # start and end year of recon
    yh <- Xr %*% ResMLR1a$Coefficients # reconstruction as 1-col matrix
    yh <- cbind(yrXr,yh) # year and recon y, as generated from the reg coefficients

    # Delta y for 50% CI
    xStdNorm75 <-qnorm(0.75, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
    deltaRec50 <- xStdNorm75 * ResMLR1a$CrossValidStorage$RMSEcv
    yhLo <- yh[,2]-deltaRec50; yhHi <-yh[,2]+deltaRec50
    yh <- cbind(yh,yhLo,yhHi) # matrix with year, recon, lower 50 upper 50
  }

  # Truncate reconstruction matrix yh to end with yrEnd, which is the last year of the
  #   quantile-extended matrix of SSRs.  Beware that some of the most recent years of
  #   the reconstruction may be based on extended SSRs
  Ltemp <- yh[,1] <= yrEnd # marks rows of yh to be retained
  yh <- yh [Ltemp,] # truncate tsm yh
  rm(Ltemp)
  
  
  # Quality control check that "fitted values" from calibration period match
  # the reconstruction for that period arrived at by applying regression coefficients to
  # selected columns of predictor matrix. This check not applicable for analog method.
  if (methMSR!=3){
    yfit1<-ResMLR1a$Model$fitted.values # from the regression model
  L<- yh[,1]>=yrgo1 & yh[,1]<=yrsp1
  yfit2 <- yh[L,2]
  dtemp = abs(yfit2 - yfit1)
  bigTemp <- max(c(max(abs(yfit1)),max(abs(yfit1)))); # largest absolute value of recon
  # during calib period, by either way of generating (regression output of manual application
  # of coefficients to columns (correct, presumabley) of predictor matrix)
  L= dtemp>1E-9 *bigTemp # maximum difference in the two versions of calib-period recon
  #     different by more than 1 billionth the largest value?
  
  if (any(L)){
    emmsg <- paste('ResMLR1a$Model$fitted.values for at least one year of calibration',
                   'period differ from reconstruction stored in yh by more than 1E-8 times largest absolute',
                   'reconstructed value in either of those series -- error from RecMLR1.R')
    ResTemp<-emssgUNH(emssg,outputDir)
    stop(emssg)
  }
} else {
  # no check for analog method
}
  
  #====================== FIGURE (local 8): 1x1 TSP OF FULL RECON WITH 50% CI
  
  # time series for plot and CI
  y <- yh[,2]; yry <- yh[,1] # for time series

  # Compute shaded polygon x and y
  Xtemp <- yh[,-2] # matrix with year as col 1, lower CI as col 2, upper CI as col 2
  ResTemp <- xyCI(Xtemp)
  xP <- ResTemp$x;   yP <- ResTemp$y
 
  # Limits for plot
  yLo <- min(yh[,2:4])
  yHi <- max(yh[,2:4])
  ynudge <- 0.02 * (yHi-yLo)
  ylims = c(yLo-ynudge, yHi+ynudge)
  xlims = c(yrXr[1]-1,yrEnd+1)
  
  # Strings for plot
  strRecYrs <- paste(sprintf('%d',yrXr[1]),'-',sprintf('%d',yrEnd),sep='')
  Tit1 <- paste('Reconstructed ',HydroName,', ',strRecYrs,
                '\n(50% CI shaded; dashed line at reconstructed mean)',sep='')
  ylab1 <- paste(HydroLabel,HydroUnits)

  #--- Build figure filename
  jFigAdd <- jFigAdd+1
  FigNumber <- NextFigNumber+jFigAdd   # for naming this png
  if (FigNumber<10){
    fileOut <- paste(outputDir,'Figure0', as.character(FigNumber),'-Reconstruction1.png',sep="")
  } else {
    fileOut <- paste(outputDir,'Figure', as.character(FigNumber),'-Reconstruction1.png',sep="")
  }
  
  # Build figure
  png(filename=fileOut, width = 960, height = 480)
  par(mar=c(5,5,5,1),cex.main=1.4,cex.axis=1.2,cex.lab=1.2)
  if (length(y)>minLength1){
    plot(yry,y,type="l",col="blue",xlim=xlims,ylim=ylims,
         ylab=ylab1,xlab='Year',main=Tit1)
  } else {
    plot(yry,y,type="b",pch=1,col="blue",xlim=xlims,ylim=ylims,
         ylab=ylab1,xlab='Year',main=Tit1)
  }
  abline(h=mean(y),lty=2,col='#808080') # dash gray
  adjustcolor("red",alpha.f=0.5) 
  #polygon(yryP,yP,col='#FFEE99') # flavescent
  polygon(xP,yP, col=rgb(1.00,0,0,0.1),border=NA) # mustard
  dev.off()

  
  
  #====================== FIGURE: Reconstruction analysis : 2x2.
  #
  # At left, top and bottom are ACFs of recon for calib years and earlier
  # At right is single frame with box plots or recon for same
  
  
  #--- Pull recon for calib perod and for earlier
  # w1, w2 will the reconsturction for those period
  # Already have y, yry as full length recon
  L <- yry >= yrgo1; # calib pd
  w1 <-y[L]; yrw1<- yry[L]
  L <- yry < yrgo1
  w2 <- y[L]; yrw2 <-yry[L]
  
  #---Make some strings for use in plots
   strPd1<- paste(sprintf('%d',yrgo1),'-',sprintf('%d',yrsp1),
                  ' (N=', sprintf('%d',length(w1)),' yr)',sep='')
   strPd2<- paste(sprintf('%d',yrgo3),'-',sprintf('%d',yrgo1-1),
                  ' (N=', sprintf('%d',length(w2)),' yr)',sep='')
   strAnnote1 <- paste('A: ',strPd1)
   strAnnote2 <- paste('B: ',strPd2)
  
  Tit1 <-'ACF of Reconstruction with 95% CI, Calibration Period'
  Tit2 <- 'ACF of Reconstruction with 95% CI, Earlier Years'
  Tit3 <- paste('Distribution of Reconstructed',HydroLabel)
  
  #--- Build figure filename
  jFigAdd <- jFigAdd+1
  FigNumber <- NextFigNumber+jFigAdd   # for naming this png
  if (FigNumber<10){
    fileOut <- paste(outputDir,'Figure0', as.character(FigNumber),'-Reconstruction2.png',sep="")
  } else {
    fileOut <- paste(outputDir,'Figure', as.character(FigNumber),'-Reconstruction2.png',sep="")
  }


    # Figure size and margins
  png(filename=fileOut, width = 960, height = 480)
  layout.matrix <- matrix(c(1,2,3,3), nrow = 2,byrow=FALSE)
  #layout(layout.matrix,heights=2,widths=c(1,1))
  layout(layout.matrix)#,heights=2,widths=c(1,1))
  par(mar=c(5,5,5,1),cex.main=1.4,cex.axis=1.2,cex.lab=1.2)
  
  #--- Upper left, acf of recon for calib pd
  
  MaxLag <- floor(length(w1)/4)
  MaxLag  <- min(c(MaxLag),20)
  
  par(mar=c(5,5,4,2))
  acf(w1,lag.max=MaxLag,ylim=c(-1,1),type='correlation',main=Tit1)
  text(MaxLag,1,strAnnote1,adj=c(1,1),cex=1.5)
  
  #---Lower left, acf of recon for years before start of calib pd
  
   MaxLag <- floor(length(w2)/4)
   MaxLag  <- min(c(MaxLag),20)
  
  par(mar=c(5,5,4,2))
   acf(w2,lag.max=MaxLag,ylim=c(-1,1),type='correlation',main=Tit2)
   text(MaxLag,1,strAnnote2,adj=c(1,1),cex=1.5)
   
  
  #---Right -- side by side box plots of recon for calib period and prior
  
  par(mar=c(5,8,4,1),cex.lab=1.3)
  namesBP<-c('Period A','Period B')
  boxplot(w1,w2,notch=FALSE,ylab=xlab1,
          main=Tit3,names=namesBP)
   dev.off()

   
  #=== TABLE: CALIBRATION1

  #--- Header
  if (PCApredictors | methMSR==3){
    TableTitle <- "Table5-Calibration1"
  } else {
    TableTitle <- "Table3-Calibration1"
  }
  if (methMSR==3){
    TitleAdd <- 'Statistics of analog MSR model'
  } else {
    TitleAdd <- 'Calibration statistics of MSR model'
  }
  SSRdef <- '   (SSR: "single-site reconstruction")'
  textH<- c(TableTitle,TitleAdd,
            paste("Predictand:",HydroName,"for",HydroSeason),
            RecMethod,SSRdef)
  # --- Body (this includes the variable labels that go in first column)
  if (methMSR==3){
    PCsUsed <- namesU2[1]
    if (length(PCsUsed)>1){
      for (n in 2:length(namesU2)){
        PCsUsed <- paste(PCsUsed,namesU2[n])
      }
    } else {
    }
    textB <-c("YearGo","YearStop","Npool","alphaR","Npredictors","RSME","RE","r")
    TfmtB <- '%-11s\t' # format for name of variable; size for longest
    # debug DfmtsB <- c('%-4d\n','%-4d\n','%-3d\n','%-6.2f\n','%-3d\n','%-8.3g\n','%-6.2f\n','%-6.2f\n')
    DfmtsB <- c('%-4.0f\n','%-4.0f\n','%-3.0f\n','%-6.2f\n','%-3.0f\n','%-8.3g\n','%-6.2f\n','%-6.2f\n')
    dataB <- c(yrgo1,yrsp1,dim(U1)[2],alphaR,dim(U2)[2],Analog$RMSE,
               Analog$RE,Analog$Correl)
    textT <- c(paste('Units of RMSE: ',HydroUnits2), 
               paste('PCs used for analogs: ',PCsUsed),
               strTailPart1)
  } else {
    textB <-c("YearGo","YearStop","Npool","Npredictors","R2","F","pF","R2adj","RMSEc")
    TfmtB <- '%-10s\t' # format for name of variable; size for longest
    DfmtsB <- c('%-4.0f\n','%-4.0f\n','%-3.0f\n','%-3.0f\n','%-6.2f\n','%-8.3g\n','%-8.3g\n','%-6.2f\n','%-10g\n')
    dataB <- c(yrgo1,yrsp1,dim(U2)[2],npred,OutputCal$Rsquared,OutputCal$F,OutputCal$pF,
               OutputCal$RsquaredAdj,OutputCal$RMSEc)
    #---Tail
    textT <- c(paste('Units of RMSEc: ',HydroUnits2),strTailPart1)
  }
  D1 <- list(textH=textH,textB=textB,TfmtB=TfmtB,dataB=dataB,DfmtsB=DfmtsB,
             textT=textT,outDir =outputDir)
  
  #---Function call for table
  
  ResTemp <- Table1Column(D1)
  rm(D1,textT,dataB,DfmtsB,TfmtB,textB,textH,TableTitle)

  

  #=== TABLE AnalysisOfResiduals

  #--- Header
  if (PCApredictors | methMSR==3){
    TableTitle <- "Table6-AnalysisResiduals1"
  } else {
    TableTitle <-  "Table4-AnalysisResiduals1"
  }

  if (methMSR==3){
    TitleAdd <- 'Normality and trend'
  } else {
    TitleAdd <- 'Normality, autocorrelation, trend, heteroskedasticity'
  }
  textH<- c(TableTitle,TitleAdd)

  if (methMSR==3){
    # --- Body
    textB <-c("YearGo","YearStop","pNormal","TrendSlope","pTrend")
    TfmtB <- '%-10s\t' # format for name of variable; size for longest
    DfmtsB <- c('%-4.0f\n','%-4.0f\n','%-8.3g\n','%-8.5g\n','%-8.3g\n')
    dataB <- c(yrgo1,yrsp1,hLillie$p.value,ResMK$b,ResMK$pvalue)
    #---Tail
    textT <- c('Tests applied include: Lilliefors for trend and Mann-Kendall, for trend',
               paste('Units of TrendSlope: ',HydroUnits2,' per year',sep=""),
               strTailPart1)
  } else {
    #I SCREWED UP HERE AND DELETED CODE
    # --- Body
    textB <-c("YearGo","YearStop","pNormal","DW","pDW","TrendSlope","pTrend","BP","dfBP","pBP")
    TfmtB <- '%-10s\t' # format for name of variable; size for longest
    DfmtsB <- c('%-4.0f\n','%-4.0f\n','%-8.3g\n','%-8.4g\n','%-8.5g\n','%-8.5g\n','%-8.3g\n',
                '%-8.3g\n','%-8.5g\n','%-8.3g\n')
    dataB <- c(yrgo1,yrsp1,hLillie$p.value,DW$dw,DW$p,ResMK$b,ResMK$pvalue,BP$ChiSquare,BP$Df,BP$p)
    #---Tail
    textT <- c('Tests applied include: Lilliefors for trend; Durbin-Watson for autocorrelation, Mann-Kendall for trend',
               '\n   and Breusch-Pagan for constancy of variance\n',
               paste('Units of TrendSlope: ',HydroUnits2,' per year',sep=""),
               strTailPart1)

  }
  D1 <- list(textH=textH,textB=textB,TfmtB=TfmtB,dataB=dataB,DfmtsB=DfmtsB,
             textT=textT,outDir =outputDir)

  #---Function call for table
  ResTemp <- Table1Column(D1)
  rm(D1,textT,dataB,DfmtsB,TfmtB,textB,TableTitle,textH)



  #=== TABLE Validation1
  
  if (methMSR==3){
    # Table not applicable
  } else {
    
    #--- Header
    if (PCApredictors){
      TableTitle <- "Table7-Validation1"
    } else {
      TableTitle <- "Table5-Validation1"
    }
    
    TitleAdd <- 'Cross-validation and split-sample validation'
    textH<- c(TableTitle,TitleAdd)
    
    # --- Body
    textB <-c("NleaveOut","RMSEcv","REcv","YearGoA","YearStopA","YearGoB","YearStopB",
              "REsplitA","REsplitB")
    TfmtB <- '%-9s\t' # format for name of variable; size for longest
    DfmtsB <- c('%-3.0f\n','%-8.5g\n','%-5.2f\n','%-4.0f\n','%-4.0f\n',
                '%-4.0f\n','%-4.0f\n','%-5.2f\n','%-5.2f\n')
    dataB <- c(ResMLR1a$CrossValidStorage$LeftOut,ResMLR1a$CrossValidStorage$RMSEcv,
               ResMLR1a$CrossValidStorage$REcv,yrgoA,yrspA,yrgoB,yrspB,ResSS1$RE,ResSS2$RE)
    
    #---Tail
    
    textT <- c('"NleaveOut" is number of observatations left out in cross-validation (cv).',
               'RMSE and RE refer to root-mean-square error and reduction-of-error.',
               'Start and end years are listed for split-sample early (A) and late (B) parts.',
               paste('Units of RMSEcv: ',HydroUnits2,sep=""),
               strTailPart1)
    
    D1 <- list(textH=textH,textB=textB,TfmtB=TfmtB,dataB=dataB,DfmtsB=DfmtsB,
               textT=textT,outDir =outputDir)
    
    #---Function call for table
    ResTemp <- Table1Column(D1)
    rm(D1,textT,dataB,DfmtsB,TfmtB,textB,TableTitle,textH,ResTemp)
  }


  
    #=== TABLE MSR coefficients
  
  if (methMSR==3){
    # table not applicable for analog method
  } else {
    #--- Title
    if (PCApredictors){
      TableTitle <- "Table8-CoefficientsMSR"
    } else {
      TableTitle <- "Table6-CoefficientsMSR"
    }
    
    TitleAdd <- 'Coefficients of MSR regression model'
    textH<- c(TableTitle,TitleAdd)
    
    # --- Body
    
    textB <- names(ResMLR1a$Coefficients) # vector of strings with names of predictors
    # if only one predictors names(ResMLR1a$Coefficients) from lm does not list the name
    # of the predictor; first entry is "(Intercept)" and second is ""
    if (npred==1){
      textB[2] <- namesU3[1]
    }
    TfmtB <- '%-12s\t' # format for name of variable; size for longest
    DfmtsB <- rep('%-12.8g\n',npred+1)
    dataB <- ResMLR1a$Coefficients
    
    #dataB <- c(M$coefficients[1],M$coefficients[2])
    
    #---Tail
    textT <- c(strTailPart1)
    
    D1 <- list(textH=textH,textB=textB,TfmtB=TfmtB,dataB=dataB,DfmtsB=DfmtsB,
               textT=textT,outDir=outputDir)
    
    #---Function call for table
    ResTemp <- Table1Column(D1)
    rm(D1,textT,dataB,DfmtsB,TfmtB,textB,TableTitle,textH,ResTemp)
  }

  

  #=== TABLE PCA loadings

  if (PCApredictors | methMSR==3){

    #--- Title
    TableTitle <- "Table3-PCA1"

    TableSubTitle <- "Loadings of principal components of single-site reconstructions"

    #--- Heading
    c1=vector()
    for (n in 1:nPC){
      c1[n] <- paste('PC',as.character(n),sep='')
    }
    c1 <- c('N','Site#','SiteID',c1)
    textH<- list(Title=TableTitle,SubTitle=TableSubTitle,Heading=c1)

    # longest id & fmt to accommodate
    nbig1 <- 8
    nbig <- max(c(nchar(nmsU),nbig1))
    fmtID <- paste('%',as.character(nbig),'s\t',sep='')
    rm(c1,nbig,nbig1)
    
    # --- Body
    nwide1=12;
    nwide2 =5;
    nwide3=3;
    
    fmtPCa <- paste('%',as.character(nwide1),'s\t',sep='') # for headers of PC cols, except last
    fmtPCb <- paste('%',as.character(nwide1),'s\n\n',sep='') # ... last
    fmtPC1 <- paste('%',as.character(nwide1),'.5g\t',sep='') # for loadings, except last
    fmtPC2 <- paste('%',as.character(nwide1),'.5g\n',sep='') # for loadings,  last
    fmtPct1 <- paste('%',as.character(nwide1),'.3g\t',sep='') # for loadings, except last
    fmtPct2 <- paste('%',as.character(nwide1),'.3g\n',sep='') # for loadings,  last
    
    textB <- list(SiteID=rownames(LoadPC),Lower=c('EV','%','Cum%')) # text for body
    TfmtB1 <- c('%4s\t','%5s\t',fmtID) # for first 3 cols
    #TfmtB2 <- c(rep('%8s\t',nPC-1),'%8s\n\n') # for the rest -- vectors of loadngs
    TfmtB2 <- c(rep(fmtPCa,nPC-1),fmtPCb) # for the rest -- vectors of loadngs
    TfmtB <- list(Left=TfmtB1,Right=TfmtB2)
    DfmtB1 <- c('%4d\t','%5d\t',fmtID)
    #DfmtB2 <-  c(rep('%8.5g\t',nPC-1),'%8.5g\n') # for the rest -- vectors of loadngs
    DfmtB2 <-  c(rep(fmtPC1,nPC-1),fmtPC2) # for the rest -- vectors of loadngs
    #DfmtB3 <-  c(rep('%8.3g\t',nPC-1),'%8.3g\n') # for pctg amd cum pctg
    DfmtB3 <-  c(rep(fmtPct1,nPC-1),fmtPct2) # for pctg amd cum pctg


    DfmtB <- list(Left=DfmtB1,Right=DfmtB2,Pctg=DfmtB3)
    dataB <- list(ResPCA=ResPCA,jScreened=jScreened,SiteID=rownames(LoadPC))

    #---Tail
    textT <- paste('Loadings are graphically shown by heat map in Figure 06\n',
                   strTailPart1,sep='')

    # Line to go above and below table (cosmetic only). No need to change this, but you
    # can if you want the width of line to perfectly match width of table as viewed in
    # some text editor. Would first need to view tabel in the text editor, draw the
    # track to desired width, and replace the line in the following statement.
    BunnyTrack <- strrep('=',(nPC+3)*(nwide1+4)) # 4 for the tabs; 3 for the cols before loading


    D1 <- list(textH=textH,textB=textB,TfmtB=TfmtB,dataB=dataB,DfmtB=DfmtB,
               textT=textT,BunnyTrack=BunnyTrack,outDir=outputDir)

    #---Function call for table

    ResTemp <- TablePCA1(D1)
    rm(D1,textT,dataB,DfmtB,TfmtB,textB,TableTitle,textH,ResTemp,BunnyTrack)
    


    #=== TABLE PCA CORRELATION WITH HYDRO
    #
    # Recall: rv1U1, r1U1, r1V1 are: vector of correlations of hydro with PCs,
    # vector of lag-1 autocorrelations of U1, and scalar of lag-1 autocorrelation
    # of v1. Refer bar plots in Figure 7
    #
    # This table will, for each PC, list the following
    # 1) PC#, 2) r with y, 3) r(1) of PC, 4) r(1) of y, threshold1, threshold2
    # Threshold1 is alpha=alphaR (e.g., 0.05) two-tailed, significance threshold, neglecting autocorrelation
    # Threshold2 is similar threshold, but using an adjustment to effective sample size:
    #    Nprime = N(1-r1r2)/(1+r1r2), where N is samples size, Nprime is effective
    #     sample size, r1, r2 are lag-1 autocorrelations of the pair of variables. If
    #     either r1 or r2 are non-positve, no adjustment is made (Nprime=N)
    #--- Title to Headings

    # Compute autocorrelation-adjusted critical levels for correlations of v1 with U1
    ResEff <- EffectSS(U1,v1) # if matrix and vector, matrix must be the first arg
    Nprime <- ResEff$Nprime # vector of effective sample sizes for correlaitons
    rCI1 = ThreshSiggy/sqrt(Nprime-2) # adjusted Siggy (e.g. 95%) CI
    rm(ResEff,Nprime)
    textH <- c('Table4-PCA2',
               paste('Correlation of PC scores with ',HydroName,' (',as.character(yrgo1),
                     '-',as.character(yrsp1),')',sep=''),
               'PC#','Corr','Thresh1','Thresh2','r(1)')
    # Formats
    fmt1<-c('%4s\t','%6s\t','%7s\t','%7s\t','%7s\n')
    fmt2<-c('%4s\t','%6.2f\t','%7.2f\t','%7.2f\t','%7.2f\n')
    fmtHB <- list(Head=fmt1,Body=fmt2)

    # Body and tail
    dataB <- list(PC=colnames(U),r=rv1U1,Thresh1=rep(rCI,nPC),Thresh2=rCI1,
                  r1PC=r1U1)
    textT <- paste('Listed are the correlations (Corr) of ',HydroName,' with PC scores,',
                   paste('\nthresholds for',Siggy,'significance of correlation neglecting (Thresh1)'),
                   '\nand accounting for (Thresh2) lag-1 autcorrelation in the two series,',
                   '\nand the lag-1 autocorrelations (r(1)) of the PCs. The two thresholds',
                   '\ndiffer only if both series have positive autocorrelation. For ',HydroName,',',
                   '\nlag-1 autocorrelation is r(1)=',sprintf('%5.2f',r1v1),'. A correlation of ',HydroName,
                   '\nwith a PC is judged "significant" only if the absolute value of',
                   '\ncorrelation is greater than the threshold.',
                   '\n\nSee Wilks (2019) for autocorrelation adjustment.',sep='',
                   '\n',PdfDescribe)

    #--- Function call to write table
    BunnyTrack <- '============================================================='
    D = list(outputDir=outputDir,textH=textH,dataB=dataB,fmtHB=fmtHB,textT=textT,
             BunnyTrack=BunnyTrack)
    ResTemp <- TableWrite1(D)
  }


  
  #=== TIME SERIES DATA: FULL-LENGTH SCREENED SSRS OR THEIR PCS
  #
  # Recall that ResPCA$Scores and ResPCA$yrScores are the data needed if PCApredictors
  # Recall than for noPCA mode, long-terms screened SSRs are in U, yrU, and ids in nmsU2

  # Title and data
  if (PCApredictors | methMSR==3){
    TableTitle <- "PCscoresTimeSeries"
    dataB <- cbind(ResPCA$yrScores,ResPCA$Scores) # Data
    fmtsB <- c('%6g\t',rep('%8.5g\t',(nPC-1)),'%8.5g\n')
    # Headings
    textH <- c('Year',colnames(LoadPC)) # text
    germFmt  <- '%8s\t'
    c1 <- rep(germFmt,(nPC-1)); c1 <-c (c1,'%8s\n')
    fmtsH <- c('%6s\t',c1) # format  for text of headings
    # Tail
    if (methMSR==3){
      textT <- paste('Scores of principal components of single-site reconstructions (SSRs).',
                     '\nNot all of these these PCs may have been used in selecting analog',
                     '\nyears (see table note for Table 5).',sep='')
    } else {
    textT <- paste('Scores of principal components of single-site reconstructions (SSRs).',
                   '\nNot all of these these PCs may have been selected as predictors for the ',
                   '\nMSR model (see Table 8)',sep='')
    }
  } else {
    TableTitle <- "ScreenedSSRtimeSeries"
    dataB <- cbind(yrU,U) # Data
    fmtsB <- c('%6g\t',rep('%8.5g\t',(dim(U)[2]-1)),'%8.5g\n')
    textH <- c('Year',namesX) # text
    germFmt  <- '%8s\t'
    c1 <- rep(germFmt,(dim(U)[2]-1)); c1 <-c (c1,'%8s\n')
    fmtsH <- c('%6s\t',c1) # format  for text of headings
    # Tail
    textT <- paste('Screened single-site reconstructions (SSRs).',
                   '\nNot all of these these SSRs may have been selected as predictors for the ',
                   '\nMSR model (see Table 6)',sep='')
  }
  D1 <- list(outDir=outputDir,filename=TableTitle,textH=textH,dataB=dataB,
             textT=textT,fmtsH=fmtsH,fmtsB=fmtsB)
  ResTemp <- TabSepTsm3(D1)
  
  

  #=== TIME SERIES DATA"  REGRESSION MODEL INPUT (DOES NOT APPLY TO ANALOG MODEL)
  if (methMSR != 3){
    
    # ResMLR1a$Model$model -- model inputs, y first, then predictors in order of how coefficients listed
    # yrv1 -- years (a vector)
    
    #--- Title
    TableTitle <- "RegressionInputTimeSeries"
    textTitle<- c(TableTitle)
    
    #--- Head
    textH <- c('Year',
               paste(HydroLabel,HydroUnits),
               namesU3)
    nmaxH <- max(nchar(textH)) # length of longest header string
    if (nmaxH<12){
      nmaxH <-12
    }
    nH <- length(textH) # number of headers
    fmtH1  <- paste('%-',as.character(nmaxH),'s\t',sep='')
    fmtH2  <- paste('%-',as.character(nmaxH),'s\n',sep='')
    fmtsH  <- rep(fmtH1,(nH-1))
    fmtsH <- c(fmtsH,fmtH2)
    
    # --- Body
    fmtsB <- c('%-12.0f\t','%-12.8g\t',rep('%-12.8g\t',(npred-1)),'%-12.8g\n')
    dataB <- cbind(as.matrix(yrv1),ResMLR1a$Model$model)
    
    #---Tail
    textT <- c("In the MSR model, second column is regressed on remaining columns.",
               strTailPart1)
    D1 <- list(filename=TableTitle,textH=textH,fmtsH=fmtsH,dataB=dataB,
               fmtsB=fmtsB,textT=textT,outDir=outputDir)
    
    #---Function call for write
    ResTemp <- TabSepTsm2(D1)
    rm(D1,textTitle,TableTitle,textH,fmtH1,fmtH2,fmtsH,dataB,fmtsB,textT,ResTemp)
  }
  

  
  #=== TIME SERIES DATA: AnalogYearsTimeSeries (applies only to analog method)
  
  if (methMSR==3){
    W <- Analog$ReconNN12 # year, reon y, analog year, indicator if nearest neighbor or second-nearest
    
    #--- Title
    TableTitle <- "AnalogYearsTimeSeries"
    textTitle<- c(TableTitle)
    
    #--- Head
    textH <- c('Year',
               paste(HydroLabel,HydroUnits),
               'Analog Year','Neighbor')
    nmaxH <- max(nchar(textH)) # length of longest header string
    if (nmaxH<12){
      nmaxH <-12
    }
    nH <- length(textH) # number of headers
    fmtH1  <- paste('%-',as.character(nmaxH),'s\t',sep='')
    fmtH2  <- paste('%-',as.character(nmaxH),'s\n',sep='')
    fmtsH  <- rep(fmtH1,(nH-1))
    fmtsH <- c(fmtsH,fmtH2)
    
    # --- Body
    fmtsB <- c('%-12.0f\t','%-12.8g\t','%-12.5g\t','%-12.0g\n')
    dataB <- W
    
    #---Tail
    textT <- c(paste(' Analog year (col 3) observed ',HydroName,' supplies reconstruction for year in col 1.',sep=''),
    'Col 4 indicates if analog year is nearest (1) or second-nearest (2) neighbor to year of reconstruction',
         strTailPart1)
    D1 <- list(filename=TableTitle,textH=textH,fmtsH=fmtsH,dataB=dataB,
               fmtsB=fmtsB,textT=textT,outDir=outputDir)
    
    #---Function call for write
    ResTemp <- TabSepTsm2(D1)
    rm(D1,textTitle,TableTitle,textH,fmtH1,fmtH2,fmtsH,dataB,fmtsB,textT,ResTemp)
  }
  

  
  #=== TIME SERIES DATA: Reconstruction with 50% confidence interval
  #
  # Status.
  # yh: 4 column tsm with year, recon, lower50%, upper 50% (matrix)
  # v, yrv: observed predictand and years (vectors)
  # v might have more recent data than yh because SSRs by lagged regression cannot extend
  #   beyond m years before the end of the tree-ring data, where m is maximum lag allowed
  #   in reconstruction model.
  # In formats, make sure that field lengths for header and data match and that the length of
  # non-year columns is at least as long as the longest header element
  fmtsH <- c('%6s\t','%10s\t','%10s\t','%10s\t','%10s\n') # for header lin
  fmtsD <- c('%6g\t','%10g\t','%10g\t','%10g\t','%10g\n') # for data matrix
  
  D1 <- list(header=c("Year",paste("Obs",HydroLabel,HydroUnits),
                      "Reconstructed","Lower 50% CI","Upper 50% CI"),
             observed=cbind(yrv,v), recon=yh,outDir=outputDir,
             fmtsH=fmtsH, fmtsD=fmtsD,
             filename="ReconstructionWithConfidenceIntervalTimeSeries")
  ResTemp <-  TabSepTsm1(D1)
  
  
  #=== ORGANIZE DATA FOR RETURN TO CALLING FUNCTION
  
  if (methMSR==3){
    ResCV <- NA
    CalData <- list('Year'=yrv1,'y'=v1,'yhat'=Fits,'PCscores' <- U2, 'Residuals'=Analog$Resids)
    CalMtx <- cbind(Analog$ReconNN12) # year, recon, analog year, whether analog year 
    OutputCal <- list('flag'=flagBail,'Msg'=flagMsg)
  } else {
    
    # calib period:  year, obs, rec, CVpredictions, e_cal, e_cv
    
    ResCV <- ResMLR1a$CrossValidStorage # Renaming list for convenience so that can
    # use code copied from RecSLR1
    
    CalData <- list('Year'=yrv1,'y'=v1,'yhat'=Fits,'yhatCV'=ResCV$CVpredictions,
                    'Residuals'=ResMLR1a$residuals,'ResidualsCV'=ResCV$CVresiduals,
                    'PredictorMtx'=Xr)
    CalMtx <- cbind(yrv1,v1,Fits,ResCV$CVpredictions,ResMLR1a$residuals,ResCV$CVresiduals)
  }
  OutputRec <- list('yearGoRec'=yrgo3,'yearSpRec'=yrsp3,'xStdNorm75'=xStdNorm75,
                    'deltaRec50'=deltaRec50,'yhat'=yh,'CalibrationData'=CalData,'CalibrationMtx'=CalMtx)
 
  
   #=== OUTPUT BACK TO CALLING FUNCTIONS
  
  Output <- c(OutputCal,OutputVal,OutputRec)
  
  return(Output)
}