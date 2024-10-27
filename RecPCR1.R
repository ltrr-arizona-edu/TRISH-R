RecPCR1 <- function(D) {
  # Multi-site reconstruction (MSR) by regression on PCs of single-site reconstructions (SSRs)
  # D. Meko 
  # Last revised 2022-08-28
  #
  # Called from a script or function (e.g., ReconAnalog1) that has generated the SSRs.
  # Consider the predictand, y, and predictor matrix, X. Here, y is regressed on the mean of the
  # SSRs in X. The method therefore is simple linear regression. 
  #
  # D is list with members:
  #   Text (1x4)s: symbol for y label; units of y, in parens; longer name of y; name of y season; Example:
  #    "RO","(mm)","Runoff", "12-month season ending in month 9"
  #   U, yrU: matrix of screened SSRs (matrix); matrix of years (1-col matrix)
  #   nmsU: names of screened SSRs (vector)
  #   v, yrv: numeric and integer; predictand and years
  #   yrsC (1x2)d   first and last year of desired calibraiton period; if NaN, default to first or last
  #     available year of overlap of u and v
  #   nNeg, nPos: integer (both positive) of max negative and possitive lag allowed by calling function in
  #     SSR modeling. This affects m in leave-m-out cross-validation
  #   NcMin: minimum acceptable number of years for calibration of MSR model
  #   PCoption: option for building pool of potential predictors from scores of PCs of SSRs
  #     1: specify k most important PCs, for example from viewing scree plot
  #     2: the m<fN PCs most highly correlated with predictand; N=number of years in calibration period
  #         of MSR model; f LE 0.20, default f=0.10
  #   f: factor that multiplied times number of calibration years puts cap on number of PCs allowed in 
  #       pool of potential predictors of MSR model (see PCoption)
  #   NextFigNumber [i] # start naming figures as Figure0?.png, where ? is NextFigNumber
  #   outputDir: folder that output files are written to <"/home/dave/AAAtrish2/test_out/">
  #
  # Notation below uses "<>" to indicate "default" values of input specifications. For example, 
  # <f=010> means the default for variable f is 0.10.  Or, with options 1 and 2, <2> indicates that
  # the default option is "2."
  
  library(car)
  library(nortest)
  source(paste(code_dir,"CrossValid2.R",sep="")) # leave-m-out cross-validation  
  source(paste(code_dir,"ssValid.R",sep="")) # split-sample validation 
  source(paste(code_dir,"mannken1.R",sep="")) # time plot and trend test of reg. resids 
  source(paste(code_dir,"stemACF.R",sep="")) # stem plot of acf, with CI & annotaton 
  source(paste(code_dir,"xyCI.R",sep="")) # compute polygon (for shaded CI) from lower and upper CI 
  source(paste(code_dir,"Table1Column.R",sep="")) # write a table file with just 1 data column
  source(paste(code_dir,"TabSepTsm1.R",sep="")) # write a tab-sep file with obs, recon, 50% CI
  source(paste(code_dir,"TabSepTsm2.R",sep="")) # write a tab-sep file with model input data, 50% CI
  
  flagBail<-0 # flag for bailing out of function
  flagMsg<-'No problems'
  minLength1 <-130 # if length of series > minLength1 in reconstruction time series plot, line
  # without plot characters is plotted; otherwise line with symbols
  
  # 1 specified calibration period too short or inconsistent with data coverage 
  # 
  #======================= UNLOAD LIST, AND RENAME SOME VARIABLES
  
  U <-D$U ; yrU<-D$yrU  # full tsm of screened pf SSRs
  v <-D$v ; yrv<-D$yrv  # full length predictand
  yrgo1<-D$yrsC[1]; yrsp1<-D$yrsC[2] # desired start and end year of calib period of MSR model
  nNeg<-D$nNeg; nPos<-D$nPos # max neg and pos lags allowed in SSR modeling
  N1 <- D$NcMin # mimimum acceptable number of years for calibration of MSR model
  PCoption <- D$PCoption # how pool of potential predictors to be filled:
  #  1 Use will specify to use most important K PCs after viewing scree plot
  # <2> The m<fN PCs most highly correlated with predictand, where N is number of 
  #   observations in calibration period of MSR model, and f is factor (see next)
  f <- D$f # If PCoption=2, f sets the maximum possible number of predictors that will
  #   be allowed in the pool of potential predictors for the MSR model. <f=0.10>
  NextFigNumber <-D$NextFigNumber # start naming figures as Figure0?.png, where ? is NextFigNumber
  outputDir <-D$outputDir # outut to be written to this system folder
  # Note that outputDir is also defined in the global environment. So, I think that I would
  # not actually need to pass outputDir as argument to functions
  HydroName <- D$Text[3] # name of hydrologic variable (e.g., "Runoff")
  HydroLabel <- D$Text[1] # label of hydro variable, for plots (e.g., "RO)
  HydroUnits <- D$Text[2] # units of hydro varialbe with parens, (e.g., "(mm
  HydroUnits2 <- substr(HydroUnits,2,(nchar(HydroUnits)-1)) # units w/o pare s
  HydroSeason <- D$Text[4] # season of hydro variable (e.g., ""12-month season ending in month 9" )
  RecMethod<- paste('Method: stepwise regression of observed',HydroLabel,'on PC scores of SSRs')
  rm(D)
  # all of the above are numeric or integer, except that yru is 1-col matrix
  
  
  #============= GET CALIBRATION DATA AND CHECK THAT CALIBRATION PERIOD LONG ENOUGH
  
  # Compute longest possible calibration period given the overlap of U and v
  yrgo2 <- max(yrU[1],yrv[1]) # earliest overlap year of U and v
  yrsp2 <- min(yrU[dim(U)[1]],yrv[length(v)]) # earliest overlap year of U and v
  
  # If specified calibration start or end is NA, replace with earliest possible start,
  # latest possible end.
  if (is.na(yrgo1)){
    yrgo1<-yrgo2  # start calibration as early as possible if yrgo1 is NA
  }
  if (is.na(yrsp1)){
    yrsp1<-yrsp2 # end calibration as late as possible if yrgo1 is NA
  }
  
  #---- Bail with error if calibration period too short or if specified calibration period not
  # possible with time coverage of U and v
  L <- (yrgo1<yrgo2) || (yrsp1>yrsp2) ||  ((yrsp1-yrgo1+1)<N1)
  if (L){
    flagBail<-1
    flagMsg<-paste('RecLR1 error: specified calibration years,',as.character(yrgo1),'-',
                   as.character(yrsp1), 'inconsistent with data coverage or less than allowable',
                   'minimum length (',as.character(N1),'yr) of calibration period')
    Output<-list('flag'=flagBail,'Msg'=flagMsg)
    return(Output)
  }
  
  
  
  
  browser()
  
  #========================= PULL CALIBRATION DATA AND REGRESS
  
  L <- yru>=yrgo1 & yru<=yrsp1
  u1 <- u[L]; yru1<- yru[L]
  
  L <- yrv>=yrgo1 & yrv<=yrsp1
  v1 <- v[L]; yrv1<- yrv[L]
  
  M <- lm(v1~u1) # model object
  M1<-summary(M)
  ncalib <- length(v1)
  
  #========================= STORE CALIBRATION STATISTICS
  
  # Significance of overall F
  # M1$fstatistic has F, dfnum, dfdenom in 1-3
  pF <-pf(M1$fstatistic[1],M1$fstatistic[2],M1$fstatistic[3],lower.tail=FALSE) 
  
  OutputCal<-list('flag'=flagBail,'Msg'=flagMsg,'lmModel'=M,'yearGoCal'=yrgo1,'yearSpCal'=yrsp1,
                  'coefficients'=M$coefficients,'Rsquared'=M1$r.squared,'F'=M1$fstatistic[1],'pF'=pF,
                  'RsquaredAdj'=M1$adj.r.squared,'RMSEc'=M1$sigma)
  
  #====================== FIGURE (local 1): 1x2 CALIBRATION PERIOD TIME PLOTS AND STATISTICS
  #
  # Figure files are numbered within this function with first figure as figure01.png. where
  # In general, figure files are named Figure?? is a number built from NextFigNumber+jFigAdd. jFigAdd 
  # will start at 0 but increment for later figures
  
  #--- Build some strings for use in plots
  strCalPd <- paste('  Calibration period: ',sprintf('%d',yrgo1),'-',sprintf('%d',yrsp1),
                   ' (N=', sprintf('%d',length(yrv1)),' yr)',sep='')
  
  #--- BuildTit1 figure png filename
  jFigAdd <- 0
  FigNumber <- NextFigNumber+jFigAdd   # for naming this png
  if (FigNumber<10){
    fileOut <- paste(outputDir,'Figure0', as.character(FigNumber),'.png',sep="")
  } else {
    fileOut <- paste(outputDir,'Figure', as.character(FigNumber),'.png',sep="")
  }
  
  # Build figure
  png(filename=fileOut, width = 960, height = 480)
  split.screen(c(2,1)) # Makes Screen 1 and 2
  split.screen(c(1,2), screen=1) # Makes Screen 3 and 4
  
  # time plots, obs and rec
  screen(2)
  zx <- c(yrgo1,yrsp1)
  zy <- c(mean(v1),mean(v1))
  par(mar=c(4,4,3,1))
  plot(yrv1,v1,type="b",pch=1,col="blue",xlim=c(yrgo1,yrsp1),
       ylab=paste(HydroLabel,HydroUnits),xlab="Year",
       main=paste('Calibration Period Observed and Reconstructed',HydroName))
  lines(yrv1,M$fitted.values,type="b",pch=2,col="red")
  lines(zx,zy)
  h<-par("usr"); yoffset<- (h[4]-h[3])/100; ytop <- h[4]-yoffset
  legend(yrgo1+1,ytop,legend=c("Obs", "Recon"),
         col=c("blue", "red"), lty=1, cex=1.2)
  
  screen(3)
  par(mar=c(4,4,2,8))
  r<-cor(M$fitted.values,v1)
  strTit <- paste('Recon vs Obs ',HydroName,', r=',as.character(round(r,digits=2)))
  plot(v1,M$fitted.values,ylab=paste('Recon',HydroLabel,HydroUnits),
       xlab=paste('Obs',HydroLabel,HydroUnits),main=strTit)
  abline(lm(M$fitted.values~v1),col="red")
  
  
  screen(4)
  par(mar=c(0,0,0,0))
  strText <-paste(RecMethod,'\n',strCalPd,
                  '\n   R-squared =',sprintf('%.2g',M1$r.squared),
                  '\n   F=',sprintf('%.5g',M1$fstatistic[1]),'(p=',sprintf('%.4g',pF),')',
                  '\n   RMSE=',sprintf('%g',M1$sigma),HydroUnits2,' (equivalent to std error of the estimate)')
  text(x=0.05,y=0.95,'Calibration Statistics',adj=c(0,1),cex=1.5,font=2)
  text(x=0.05,y=0.8,strText,adj=c(0,1),cex=1.2)
  
  dev.off()

  #====================== FIGURE (local 2): 2x2 DISTRIBUTIONS AND ACFS
  #
  # CW from LL: histogram of yhat;  hist of obs y; acf of obs y, acf of yhat
  #--- Uniform xlims for histograms
  xlo = min(c(min(v1),min(M$fitted.values)))
  xhi = max(c(max(v1),max(M$fitted.values)))
  xinc <- 0.05*(xhi-xlo)
  xlims <- c(xlo-xinc,xhi+xinc)
  rm(xlo,xhi,xinc)
  xlims <-c((c(min(v1),min(M$fitted.values))), (c(max(v1),max(M$fitted.values))))
  
  #--- Build figure png filename
  jFigAdd <- jFigAdd+1
  FigNumber <- NextFigNumber+jFigAdd   # for naming this png
  if (FigNumber<10){
    fileOut <- paste(outputDir,'Figure0', as.character(FigNumber),'.png',sep="")
  } else {
    fileOut <- paste(outputDir,'Figure', as.character(FigNumber),'.png',sep="")
  }
  # Build figure
  png(filename=fileOut, width = 960, height = 480)
  # Create the layout
  nf <- layout( matrix(c(1,2,3,4), ncol=2) )
  
  nbin <- floor(5*log10(length(v1))) # Panofsky rule of thumb for number of bins
  brks1 = seq(min(v1)*.99999,max(v1)*1.00001,length.out=(nbin+1))
  brks2 = seq(min(M$fitted.values)*.99999,max(M$fitted.values)*1.00001,length.out=(nbin+1))
  xlab1 <- paste(HydroLabel,HydroUnits)
  MaxLag <- floor(length(v1)/4)
  
  par(mar=c(5,6,2,2),cex.axis=1.5, cex.lab=1.5, cex.main=1.5)
  Tit1 <- paste('Histogram, Observed',HydroLabel,'(N=',sprintf('%d',length(v1)),'yr)')
  hist(v1,main=Tit1,breaks=brks1,xlim=c(150,450),xlab=xlab1)
  
  par(mar=c(5,6,2,2))
  Tit1 <- paste('Histogram, Reconstructed',HydroLabel,'(N=',as.character(length(v1)),'yr)')
  hist(M$fitted.values,main=Tit1,breaks=brks2,xlim=c(150,450),xlab=xlab1)
  
  par(mar=c(5,5,4,2))
  Tit1 <- paste('ACF, Observed',HydroLabel,',with 95% CI')
  acf(v1,lag.max=MaxLag,type='correlation',main=Tit1)
  
  par(mar=c(5,5,4,2))
  Tit1 <- paste('ACF, Reconstructed',HydroLabel,', with 95% CI')
  acf(M$fitted.values,lag.max=MaxLag,type='correlation',main=Tit1)
  
  dev.off()

  #========================= ANALYSIS OF RESIDUALS (local 3,4,5)
  #
  # Will need the regression residuasl. First figure (1x2) will be histogram and
  # scatter of residuals on predicted values. Second figure (1x1) will be time plot of
  # the residuals with a fitted (non-parametric fit) trend line and annotated result
  # of Mann-Kendall trend test. The significance will be adjusted as needed for autocorrelation
  # of the residuals over and above that in a linear trend. The third figure (1x1) will be
  # the acf of residuals, with annotated DW test results
  
  
  #--- HISTOGRAM AND SCATTER OF RESIDS ON PREDICTED (1X)
  #
  # 1st of 3 analysis of residuals plots
  
  # Lilliefors test of normality of residuals
  hLillie <- lillie.test(M$residuals);
  Tit1 <- paste('Residuals, E, of regression of',HydroName,'on tree rings',
  '\n  (p=',sprintf('%.2g',hLillie$p.value),' for H0 that E normal, from Lilliefors Test)')
  
  # Breusch-Pagan test for heterogeneity of regresson residuals
  BP <- ncvTest(M)
  Tit2 <- paste('Scatter of Residuals on Fitted Values',
                '\n  (p=',sprintf('%.2g',BP$p),' for H0 that E homoscedastic)',
                '\n(from Breusch-Pagan Test)')
 
  # Buld filename for plot
  jFigAdd <- jFigAdd+1
  FigNumber <- NextFigNumber+jFigAdd   # for naming this png
  if (FigNumber<10){
    fileOut <- paste(outputDir,'Figure0', as.character(FigNumber),'.png',sep="")
  } else {
    fileOut <- paste(outputDir,'Figure', as.character(FigNumber),'.png',sep="")
  }
  
  # Figure size and margins
  png(filename=fileOut, width = 960, height = 480)
  layout.matrix <- matrix(c(1,2), nrow = 1, ncol = 2)
  layout(layout.matrix,heights=2,widths=c(1,1)) 
  
  # Left plot

  par(mar = c(5.1, 4.5, 5.1, 2.1),cex.axis=1.1, cex.lab=1.5, cex.main=1.3)
  hist(M$residuals,xlab=paste('Residual',HydroUnits),ylab='Frequency',main=Tit1)
   
  # # right plot
  par(mar = c(5.1, 4.1, 6.1, 2.1),cex.axis=1.1, cex.lab=1.5, cex.main=1.3)
  plot(M$fitted.values,M$residuals,xlab=paste('Fitted Values',HydroUnits),
       ylab=paste('Residual',HydroUnits),
       main=Tit2)
  abline(h=0,lty=2,col='#808080') # dash gray
  dev.off()
  

  #--- TIME PLOT OF REGRESSION RESIDUALS, WITH MANN-KENDAL TREND TEST (1X1) 
  
  jFigAdd <- jFigAdd+1
  FigNumber <- NextFigNumber +jFigAdd 
  
  # Prepare input for mannken1
  X  <- cbind(yrv1,M$residuals) # matrix with year and regression residuals

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
  acfMy <- acf(M$residuals, lag.max=lagsPlot, type = "correlation", 
               plot = FALSE)
  k <- acfMy$lag # lags
  w <- acfMy$acf # acf
  
  # DW statistic of the regression residuals
  DW <-  durbinWatsonTest(M)
  strDW <- paste('p=',sprintf('%g',DW$p),
                 ': Durbin-Watson test (2-sided) of H0 that population lag-1 autocorrelation is zero',
                 sep='')
  
    # Text for plot
  Tit1 <- paste('ACF of Residuals with 95% CI (N=',sprintf('%g',ncalib),' yr)')
  textPlot <- c(Tit1,'Lag(yr)','r',strDW)
  
  # Store inputs required by stemACF in a list (see opening comments there)
  Din <- list(x=k,y=w,nsize=ncalib,kAlpha=1,FigNumber=FigNumber,
              outputDir=outputDir,linecol1='#0022CC',linecol2='#696969',
              linecol3='#E60000',textPlot=textPlot)
  ResNull <- stemACF(Din)
  

  #========================= VALIDATE AND STORE STATISTICS
  
  ResCV <- CrossValid2(u1, v1, nNeg,nPos) # cross-validation
  
  # Split-sample validation
  
  iAstop <- ceiling(length(v1)/2) # end row index in v1 of first half of data, assumed longer than
  # longer of the two halves if length of v1 odd
  iBgo <- iAstop+1 # start row of second half
  iA <- 1:iAstop # row indices of first half of full calib period
  iB <- iBgo:length(v1) # ... of second half
  
  #--- Calibrate on early, validate on late, then reverse
  ical<-iA; ival<-iB
  i1 <- 1; # col 1 of u1 is predictor; moot because here u1 is vector
  ResSS1=ssValid(v1,u1,ical,ival,i1);
  REa1<-ResSS1$RE # RE for calib on early, valid on late
  ical<-iB; ival<-iA
  ResSS2=ssValid(v1,u1,ical,ival,i1);
  REb1<-ResSS1$RE # RE for calib on late, valid on early
  
  OutputVal<-list('mLeaveOutCV'=ResCV$LeftOut,'REcv'=ResCV$REcv,'RMSEcv'=ResCV$RMSEcv,
                  'REcalEarlyValLate'=REa1,'REcalLateValEarly'=REb1)
  
  #====================== FIGURE (local 6): time plots of obs y, recon, y, cv predictions of y; with observed mean line
  #
  #--- Build figure png filename
  jFigAdd <- jFigAdd+1
  FigNumber <- NextFigNumber+jFigAdd   # for naming this png
  if (FigNumber<10){
    fileOut <- paste(outputDir,'Figure0', as.character(FigNumber),'.png',sep="")
  } else {
    fileOut <- paste(outputDir,'Figure', as.character(FigNumber),'.png',sep="")
  }
  
  # Recall
  # v1,yrv1 is observed;  M$fittedvalues is recon;  ResCV$CVpredictions is cv-prdicted
  # zx,zy are 1x2s that allow for line at observd mean; yrgo1,yrsp1 are 1st & last of calib years
  
  # Build figure
  png(filename=fileOut, width = 960, height = 480)
  par(mar=c(4,4,5,1),cex.main=1.4)
  plot(yrv1,v1,type="b",pch=1,col="blue",xlim=c(yrgo1,yrsp1),
       ylab=paste(HydroLabel,HydroUnits),xlab="Year",
       main=paste('Time Plots of Observed, Reconstructed, and Cross-Validation-Predicted',HydroName,'(',sprintf('%d',yrgo1),'-',sprintf('%d',yrsp1),')',
                                   '\n(black line at observed mean)'))
  lines(yrv1,M$fitted.values,type="b",pch=2,col="red")
  lines(yrv1,ResCV$CVpredictions,type="b",pch=17,col="#990099")
  lines(zx,zy)
  h<-par("usr"); yoffset<- (h[4]-h[3])/100; ytop <- h[4]-yoffset
  legend(yrgo1+1,ytop,legend=c("Obs", "Recon","CVpred"),
         col=c("blue", "red","#990099"),pch=c(1,2,17),lty=1, cex=1.2)
  dev.off()
  

  #====================== FIGURE (local 7): 1x2 VALIDATION PLOTS AND STATISTICS
  
  # Lilliefors test of normality of CV residuals
  hLillie <- lillie.test(ResCV$CVresiduals);
  Tit1 <- paste('Cross-Validation Residuals',
                '\n(p=',sprintf('%.2g',hLillie$p.value),', H0: normally distributed',
                ' [Lilliefors Test])',sep='')
  
  # Text strings for cross-validation
  strAnPd <- strCalPd # built for earlier plot: gives calib period and length

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

  
  #--- Build figure filename
  jFigAdd <- jFigAdd+1
  FigNumber <- NextFigNumber+jFigAdd   # for naming this png
  if (FigNumber<10){
    fileOut <- paste(outputDir,'Figure0', as.character(FigNumber),'.png',sep="")
  } else {
    fileOut <- paste(outputDir,'Figure', as.character(FigNumber),'.png',sep="")
  }
  
  # Figure size and margins
  png(filename=fileOut, width = 960, height = 480)
  layout.matrix <- matrix(c(1,2), nrow = 1, ncol = 2)
  layout(layout.matrix,heights=2,widths=c(1,1))
  
  # Left plot histogram
  
  par(mar = c(5.1, 4.5, 5.1, 2.1),cex.axis=1.1, cex.lab=1.5, cex.main=1.3)
  hist(ResCV$CVresiduals,xlab=paste('Residual',HydroUnits),ylab='Frequency',main=Tit1)
  
  # right plot, stats
  
  par(mar = c(0,0,0,0))
  xlims <- c(0,1); ylims <- c(0,1)
  plot(0,type='n',axes=FALSE,ann=FALSE,xlim=xlims,ylim=ylims)
  #plot(1,1,pch=1,xlim=xlims,ylim=ylims)
  strText <-paste('\nCross-validation (cv) method: leave-',as.character(ResCV$LeftOut),'out',
                  '\n',strAnPd,
                  '\n   RMSEcv=',sprintf('%g',ResCV$RMSEcv),HydroUnits2,
                  '\n   REcv=',sprintf('%.2g',ResCV$REcv),
                  '\n\nSplit-sample validation',
                  '\n',strSplitA,
                  '\n',strSplitB,
                  '\n    RE{A}=',sprintf('%.2g',ResSS1$RE),' (calibrated on A, validated on B)',
                  '\n    RE{B}=',sprintf('%.2g',ResSS2$RE),' (calibrated on B, validated on A)')
  text(x=0.05,y=0.95,'Validation statistics',adj=c(0,1),cex=1.5,font=2)
  text(x=0.05,y=0.8,strText,adj=c(0,1),cex=1.2)
  
  
  dev.off()
  

  #========================= RECONSTRUCTION AND 50% CI
  
  Xr <- as.matrix(u) # matrix, 1 col, of long-term mean of SSRs
  L<-complete.cases(Xr)
  Xr <- as.matrix(Xr[L,]) ; yrXr <- as.matrix(yru[L,])
  mXr <- dim(Xr)[1]
  yrgo3 <- yrXr[1,1]; yrsp3 <- yrXr[mXr,1] # start and end year if recon
  
  # Add ones column and reconstruct
  Xones<-matrix(1,nrow=mXr,ncol=1)
  Xr <- cbind(Xones,Xr)
  yh <- Xr %*% M$coefficients # reconstruction as 1-col matrix
  yh <- cbind(yrXr,yh)
  
  # Delta y for 50% CI
  xStdNorm75 <-qnorm(0.75, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
  deltaRec50 <- xStdNorm75 * ResCV$RMSEcv
  yhLo <- yh[,2]-deltaRec50; yhHi <-yh[,2]+deltaRec50
  yh <- cbind(yh,yhLo,yhHi) # matrix with year, recon, lower 50 upper 50
  
  
  #====================== FIGURE (local 8): 1x1 TSP OF FULL RECON WITN 50% CI
  
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
  xlims = c(yrXr[1]-1,yrXr[length(yrXr)]+1)
  
  # Strings for plot
  strRecYrs <- paste(sprintf('%d',yrXr[1]),'-',sprintf('%d',yrXr[length(yrXr)]),sep='')
  Tit1 <- paste('Reconstructed ',HydroName,', ',strRecYrs,
                '\n(50% CI shaded; dashed line at reconstructed mean)',sep='')
  ylab1 <- paste(HydroLabel,HydroUnits)
  
  #--- Build figure filename
  jFigAdd <- jFigAdd+1
  FigNumber <- NextFigNumber+jFigAdd   # for naming this png
  if (FigNumber<10){
    fileOut <- paste(outputDir,'Figure0', as.character(FigNumber),'.png',sep="")
  } else {
    fileOut <- paste(outputDir,'Figure', as.character(FigNumber),'.png',sep="")
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

  
  #====================== FIGURE (local 9): 2x2. 
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
    fileOut <- paste(outputDir,'Figure0', as.character(FigNumber),'.png',sep="")
  } else {
    fileOut <- paste(outputDir,'Figure', as.character(FigNumber),'.png',sep="")
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
  
  
  #=== TABLE Table3-MSR-LR1, summarizing calibration of model
  
  #--- Header
  TableTitle <- "Table3-MSR-LR1-Calibration"
  SSRdef <- '   (SSR: "single-site reconstruction")'
  textH<- c(TableTitle,
            paste("Predictand:",HydroName,"for",HydroSeason),
            RecMethod,SSRdef)
  
  # --- Body
  textB <-c("YearGo","YearStop","R2","F","pF","R2adj","RMSEc")
  TfmtB <- '%-10s\t' # format for name of variable; size for longest
  DfmtsB <- c('%-4d\n','%-4d\n','%-6.2f\n','%-8.3g\n','%-8.3g\n','%-6.2f\n','%-10g\n')
  dataB <- c(yrgo1,yrsp1,OutputCal$Rsquared,OutputCal$F,OutputCal$pF,OutputCal$RsquaredAdj,
              OutputCal$RMSEc)

  #---Tail
  textT <- c(paste('Units of RMSEc: ',HydroUnits2),
             "See TrishOutputDescribeLR1.pdf for definitions of variables in column 1")
  D1 <- list(textH=textH,textB=textB,TfmtB=TfmtB,dataB=dataB,DfmtsB=DfmtsB,
          textT=textT,outDir =outputDir)
  
  #---Function call for table
  ResTemp <- Table1Column(D1)
  
  rm(D1,textT,dataB,DfmtsB,TfmtB,textB,textH,TableTitle)

  
  
  #=== TABLE Table4-MSR-LR1-AnalysisOfResiduals
  
  #--- Header
  TableTitle <- "Table4-MSR-LR1-AnalysisOfResiduals"
  textH<- c(TableTitle)

  # --- Body
  textB <-c("YearGo","YearStop","pNormal","DW","   pDW","TrendSlope","   pTrend",
            "BP Test ChiSq","   dfBP","   pBP")
  TfmtB <- '%-13s\t' # format for name of variable; size for longest
  DfmtsB <- c('%-4d\n','%-4d\n','%-8.3g\n','%-5.2g\n','%-8.3g\n','%-8.5g\n','%-8.3g\n',
              '%-8.5g\n','%-4d\n','%-8.3g\n')
  dataB <- c(yrgo1,yrsp1,hLillie$p.value,DW$dw,DW$p,ResMK$b,ResMK$pvalue,
             BP$ChiSquare,BP$Df,BP$p)

  #---Tail
  
  textT <- c(paste('Units of TrendSlope: ',HydroUnits2,' per year',sep=""),
             "See TrishOutputDescribeLR1.pdf for definitions of variables in column 1")
             
  D1 <- list(textH=textH,textB=textB,TfmtB=TfmtB,dataB=dataB,DfmtsB=DfmtsB,
             textT=textT,outDir =outputDir)
  
  #---Function call for table
  ResTemp <- Table1Column(D1)
  rm(D1,textT,dataB,DfmtsB,TfmtB,textB,TableTitle,textH)
  
  
  
  #=== TABLE Table5-MSR-LR1-Validation
  
  #--- Header
  TableTitle <- "Table5-MSR-LR1-Validation"
  textH<- c(TableTitle)
  
  # --- Body
  textB <-c("NleaveOut","RMSEcv","REcv","YearGoA","YearStopA","YearGoB","YearStopB",
  "REsplitA","REsplitB")
  TfmtB <- '%-9s\t' # format for name of variable; size for longest
  DfmtsB <- c('%-3d\n','%-8.5g\n','%-5.2f\n','%-4d\n','%-4d\n',
              '%-4d\n','%-4d\n','%-5.2f\n','%-5.2f\n')
  dataB <- c(ResCV$LeftOut,ResCV$RMSEcv,ResCV$REcv,
             yrgoA,yrspA,yrgoB,yrspB,ResSS1$RE,ResSS2$RE)
  
  #---Tail
  
  textT <- c(paste('Units of RMSEcv: ',HydroUnits2,sep=""),
             "See TrishOutputDescribeLR1.pdf for definitions of variables in column 1")
  
  D1 <- list(textH=textH,textB=textB,TfmtB=TfmtB,dataB=dataB,DfmtsB=DfmtsB,
             textT=textT,outDir =outputDir)
  
  #---Function call for table
  ResTemp <- Table1Column(D1)
  rm(D1,textT,dataB,DfmtsB,TfmtB,textB,TableTitle,textH,ResTemp)

  
  #=== TABLE Table6-MSR-LR1-Coefficients
  
  #--- Title
  TableTitle <- "Table6-MSR-LR1-Coefficients"
  textH<- c(TableTitle)

  # --- Body
  textB <- c("Intercept","meanSSR") # headings of cols
  TfmtB <- '%-12s\t' # format for name of variable; size for longest
  DfmtsB <- c('%-12.8g\n','%-12.8g\n')
  dataB <- c(M$coefficients[1],M$coefficients[2])
  
  #---Tail
    textT <- c("See TrishOutputDescribeLR1.pdf for definitions of variables in column 1")
  
  D1 <- list(textH=textH,textB=textB,TfmtB=TfmtB,dataB=dataB,DfmtsB=DfmtsB,
             textT=textT,outDir=outputDir)
  
  #---Function call for table
  ResTemp <- Table1Column(D1)
  rm(D1,textT,dataB,DfmtsB,TfmtB,textB,TableTitle,textH,ResTemp)
  
  
  #=== TIME SERIES DATA"  RegressionInput
  
   #--- Title
  TableTitle <- "RegressionInputTimeSeries"
  textTitle<- c(TableTitle)
  
  #--- Head
  textH <- c('Year',
             paste(HydroLabel,HydroUnits),'meanSSD')
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
  fmtsB <- c('%-12.0f\t','%-12.8g\t','%-12.8g\n')
  dataB <- cbind(as.matrix(yrv1),M$model)
  
  #---Tail
  textT <- c("For the LR1 regression model, second column is regressed on third",
               "See TrishOutputDescribeLR1.pdf for definitions of columns")
  D1 <- list(filename=TableTitle,textH=textH,fmtsH=fmtsH,dataB=dataB,
             fmtsB=fmtsB,textT=textT,outDir=outputDir)
  
  #---Function call for write
  ResTemp <- TabSepTsm2(D1)
  rm(D1,textTitle,TableTitle,textH,fmtH1,fmtH2,fmtsH,dataB,fmtsB,textT,ResTemp)
  
  
  #=== TIME SERIES DATA: Reconstruction with 50% confidence interval
  
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
             filename="ObservedAndReconstructedTimeSeries")
  
  ResTemp <-  TabSepTsm1(D1)
  
  
  
  
  
  #=== ORGANIZE DATA FOR RETURN TO CALLING FUNCTION
  #
  # calib period:  year, obs, rec, CVpredictions, e_cal, e_cv
  
  CalData <- list('Year'=yrv1,'y'=v1,'yhat'=M$fitted.values,'yhatCV'=ResCV$CVpredictions,
                  'Residuals'=M$residuals,'ResidualsCV'=ResCV$CVresiduals,
                  'PredictorMtx'=Xr)
  CalMtx <- cbind(yrv1,v1,M$fitted.values,ResCV$CVpredictions,M$residuals,ResCV$CVresiduals)
  
  OutputRec <- list('yearGoRec'=yrgo3,'yearSpRec'=yrsp3,'xStdNorm75'=xStdNorm75,
                    'deltaRec50'=deltaRec50,'yhat'=yh,'CalibrationData'=CalData,'CalibrationMtx'=CalMtx)
  

  
  
  
  
  
  #=== OUTPUT BACK TO CALLING FUNCTIONS
  
  Output <- c(OutputCal,OutputVal,OutputRec)

  return(Output)
}