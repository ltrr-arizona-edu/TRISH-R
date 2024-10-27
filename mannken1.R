mannken1 <- function(Din){
  # X,kopt,outputDir,kplot,NextFigNumber
  # Mann-Kendall trend test for a time series
  # D. Meko
  # Last revised 2023-11-24
  #
  #---IN
  #
  # Input is a list with following fields:
  #
  # X [matrix]: year as col 1, time series as col 2 (see Notes)
  # kopt [vector]1x2;  options
  #   [1] Plot of time series, trend line, and detrended series
  #     1 skip plotting
  #     2 plot
  #   [2] Adjust signficance of Mann-Kendall statistic for lag-1 autocorrelations
  #     1 Yes
  #     2 No
  # kplot [scalar] which version of plots to produce
  #   1 Plot time series with fitted trend line & horiz line at mean; figure file
  #     named "mannken1_F1.png"
  #   2 Likewise plot time series with fitted trend line;n but this is TRISH-specific,
  #     intended for looking at trend in regression residuals in context of other
  #     calls; horiz line at 0, and file named like "Figure??-AnalysisResduals2.png", with name 
  #     part ?? controlled by input arg NextFigNumber
  #   3 The GEOS485A version:  2x1 plot of trend-fit to time series (top) and detrended series (bottom)
  # outputDir [char] folder to which any plots go (e.g., '/home/dave/test_out/')
  #      If "Null", this plotting of pngs to an output folder is ignored
  # textPlot [vector, char]3  title, xlabel, ylabel for time series plot
  # NextFigNumber [integer]: if called for TRISH plot (kplot=2), the figure file
  #   will be named Figure0?.png, where ? is NextFigNumber. Ignored otherwise.
  #
  #---OUT
  #
  # Output is list with elements:
  #   What: list with three elements telling (1) which function created Output, 
  #     (1) user to refer to comment section of that function to get details on the list items, and 
  #     (2) the date on which the Output list was created.
  #   statistic: Mann-Kendal statistic (see Notes and References)
  #   AnalysisPeriod: string first and last times of analysis period (e.g., '1950-2020')
  #   Lacf(L) request possible adjustment of signficance for positive significant lag-1 autocorrelation
  #     (only enacted if the residuals to a least-squares straight line fit to the series have
  #     significant positive lag-1 autocorrelation (0.05, one-tailed test)
  #   Ladjusted (L): whether autcorrelation adjustment ended up being applied
  #   vif [scalar] variance inflation factor (set to 1 if have not requested autcorrelation adjustment
  #     or if that adjustment not warranted by the data); see references on Mann Kendall test
  #   pvalue:  p-value for significance of MK test; two-tailed null hypothesis of no trend
  #   ngrp: number of groups of ties
  #   nties: total number of ties
  #   Lflag (1x2)L   flag
  #     (1) inadequate sample size (need 10 or more obs)
  #       T: sample size too small (fewer than 10 obs)
  #       F: sample size not too small
  #     (2) identically 0 slope (summation needed for test statistic dentically 0; this would 
  #       result in an infinite test statistic; See Haan,2002)
  #         T: identically 0 slope, returns original as detrended without going through non-parametric fit
  #         F:   slope either positive or negative (not necessarily significant); non-parametric fit proceeds
  #   b (1x1)r  slope (nonparametric estimate)
  #   a (1x1)r  intercept ...
  #   equation: equation for trend line (nonparametric fit)
  #   X (mx x 2)r:  time series matrix (time as col 1) of original series for specified analysis period
  #   xhat (mx x 1)r   trend line 
  #   xdetrended (mx x 1)r detrended x (by non-parametric trend line)
  #   ErrorMessage [vector]c : error message associated with Lflag
  #
  #--=NOTES
  #
  # X: assumed 2-col matrix, year as col 1, data as col2. Assumed that x has not missing values
  #   and that yrx increments by 1
  # Lflag[1]: If time series has fewer than 10 observations, Lflag(1)=T, and Lflag and ErrorMessage are the only
  #   output list element returned; here Lflag[2] is set to FALSE
  # Lflag[2]: set to T if essentially no slope in trend line. If so, the test statistic cannot 
  #   be computed because a denominator in the equation for the test statistic is zero. In this case
  #   Lflag and ErrorMessage are the only list elements returned.
  # pvalue: this is for a two-tailed test. H0 is no trend. A small p-value indicates reject H0. 
  #   For example, if p-value==0.09, we reject H0 at alpha=0.10 level. If p-value==0.0499, we 
  #   reject H0 at alpha=0.05.  
  # nties, ngrp: handling of ties follows Salas(1993)
  # b, a, equation: Nonparametric trend line fit following equations in Haan (2002). The detrended series 
  #   Result.xdetrended is shifted such that has same median as input series x. 
  # Autocorrelation adjustmentment. If requested, applied only if residuals rom an least squared fit 
  #   straight line fit to the original time series have significant (0.05 alpha) lag-1 positive 
  #   autocorrelation by a 1-tailed test. 
  # vif, or variance inflation factor: annotation of variance inflation factor (vif) at upper left
  #   of time series with fitted line. If kopt[2]==2 (you do not enable autocorrelation djustment), vif 
  #   is not annotated. Otherwise, vif is annotated, but will be "vif=1.0" if the autocorrelation  is
  #   not justified (no significant lag-1 autocorrelation in residuals from least-squares-fit line)
  # Horizontal gray dashed gray line on plots: This is at the median of x if kplot= 1 or 3, and at 0 if kplot=2.
  #   The rationale for 0 is that with kplot=2 we are chacking for trend in regression residuals, which should
  #   should have a mean of zero for the calibraton period of the regression model. The median residual can differ
  #   from zero.
  #
  # Test data from Haan (2002) used to originally check results
  #   x<- c(25.56, 33.28, 34.03, 35.72, 39.33, 32.21, 30.76, 44.45, , 42.69)
  #   t <- 1978:1987;
  #   Cbind t and x into X and call mannken1.R; compare results to those in Haan(2002)
  #
  #---REFERENCES
  # Haan, C. T. (2002). Statistical methods in hydrology (Second ed.). Iowa State University Press. (496 pp)
  # Helsel, D. R., & Hirsch, R. M. (1992). Statistical methods in water research. Amsterdan, The Netherlands: Elsevier.
  # Salas, J. D. (1992). Analysis and modeling of hydrologic time series. In D. R. Maidment (Ed.), Handbook of hydrology (p. 1-72). New York: McGraw-Hill, Inc.
  # Wilks, D. S. (2019). Statistical methods in the atmospheric sciences (Fourth ed.). Cambridge, MA: Elsevier. (818 pp)
  # 
  # Algorithms for Mann-Kendall statistic and adjustment of its signficance for autocorrelation from 
  # Wilks (2019). Handling of ties in time series as recommended by Salas (1992). Test data from Haan (2002). Detreding
  # is done following Haan (2002, p. 345), who gives equations for non-parametric estimation of slope and intercept of a linear
  # trend line. Haan (2002) got the equations for the estimation from Helse and Hirsch (1992).
  #
  # revised 2023-11-24. minor correction to labeling of plots
  
  source(paste(code_dir,"ties1.R",sep="")) # optional transformation of flows
  

  #--- UNLOAD
  
  X <- Din$X; kopt <- Din$kopt; kplot <- Din$kplot
  outputDir <- Din$outputDir; NextFigNumber <- Din$NextFigNumber 
  textPlot <- Din$textPlot
  
    # Hard code, needed if kplot==2 (special case of TRISH plot)
  jFigAdd <- 0 # increment this for any plot after first
  ErrorMessage <- "No problems"
  Lflag <- c(FALSE,FALSE) # initialize error flags
  
  # Check input
  L1  <- is.matrix(X)
  L2 <- dim(X)[2]==2
  x <- X[,2];   yrx<- X[,1]  # vectgrs
  L3 <- !any(is.na(x)) && all(diff(yrx)==1)
  L = L1 & L2 & L3
  if (!L){
    stop('Something amiss with in put X; should be 2-col matrux with year as col 1; no missing values')
  }

  # Need at least 10 observations
  mx <-length(x)
  if (mx<10){
    Lflag[1]<-TRUE
    Lflag[2]<-FALSE
    ErrorMessage <- 'Few than 10 observations in series; cannot run function mannken1'
    Output <- list(Lflag=Lflag,ErrorMessage=ErrorMessage)
    return(Output)
  }


  # Build string label of time coverage (e.g., '1890-1989')
  strFL <- paste(as.character(yrx[1]),'-',as.character(yrx[length(yrx)]),sep='')

  #---Optional check for lag-1 autocorrelation of residuals from least-squares straight-line  trend fit to x
  Ladjusted <-FALSE
  Lacf <- FALSE
  vif <-1 # initialize to effectively make no autocorrelaton adjustment to variance of the statistic 
  

  if (kopt[2]==1){
 
    # adjustment to be explored
    Lacf<-TRUE
    M <- lm(x~yrx) # regress x on yrx
    e <- M$residuals # residuals from straight line fit to x vs yrx (trend line)
    r <- acf(e,plot=FALSE) # autocorrelation function of those residuals; acf object
    r1 <- r$acf[2] # lag-1 autocorrelation
    #rm(M,e,r,acf)

    # threhsold for statistically significant positive lag-1 autocorrelation (alpha=0.05)
    r95 <- (-1+1.645*sqrt(mx-2))/(mx-1)
    if (r1>r95){
      Ladjusted=TRUE
      f <- (1+r1)/(1-r1)
      Nprime <- floor(mx/f) # effective sample size
    }else{
      f <-1
      Nprime <- mx
    }
    vif <-f; rm(f)
  } else {
    # do not consider adjustment  
  }

  #---BUILD SUMS FOR COMPUTATION OF STATISTIC
  
  # Col-dupe, then col-dupe x
  A <- matrix(x,nrow=mx,ncol=mx,byrow=FALSE)
  B <- matrix(x,ncol=mx,nrow=mx,byrow=TRUE)
  
  # Difference matrix, C
  C <- A-B
  # Consider elements below the diagonal of C:
  #   col 1 is difference of all succeeding values and x(1) 
  #   col 2 is difference of all succeeding values and x(2) 
  #... col (mx-1) is difference of x(mx ) and x(mx-1)

  # Because interested only in elements below the diagonal, convert 
  # elements above diagonal to 0; then lop off first row and last col
  L<-upper.tri(C,diag=FALSE)
  C[L]=0; 
  D <- C[-1,]; D <- D[,-mx]

  # Logicals for positive and negative differences
  Lp <- D>0;
  Ln <- D<0;
  Lz <- D==0;
  
  # Quantities for the test statistic
  E <- D
  E[Lp] <- (-1)
  E[Ln] <- 1
  E[Lz] <- 0
  
  # Test statistic is based on difference of number of positive and negative differences
  s1 <- sum(E)
  s <- sum(s1)
  

  # There is and adjustent for ties in values of x; deal with that
  T <- ties1(x)
  if (length(T$ngroups)==0){
    vties <- 0
    nties <-0 
    ngrp <-0
  } else {
    ngrp <- T$ngroups
    nties <- sum(T$nties) # total number of x involved in ties
    w<-0
    for (k in 1:ngrp){
      e <- T$nties[k]
      h = e*(e-1)*(2*e+5)
      w <- w+h
    }
    vties <-w
  }

  #======== COMPUTE TEST STATISTIC, u
  
  N <- mx
  if (s>0){
    m <- (-1)
  } else if (s<0){
    m <- 1
  } else {
    # special case of exactly as many positive as negative differences. This indicates no trend, and
    # also causes problem n computation of u because s occurs in a denominator
    m <-0
    Lflag <- c(FALSE,TRUE)
    ErrorMessage <- 'In mannken1, sums of neg and pos s are equal. For sure no trend'
    Output <- list(Lflag=Lflag,ErrorMessage=ErrorMessage)
    return(Output)
  }  

  if (m==0){
    # For case of no slope, do not try to compute Mann-Kendall statistic
  } else {
    v <- ((N*(N-1)*(2*N+5))-vties)/18 # "variance of the sampling distribution of S" Wilks, 2019, p 173
    if (kopt[2]==1){
      v <- vif*v # adjust variance with variance inflation factor, if enabled and warranted
    }  else {
    }
    u <- (s+m)/sqrt(v); # test statistic for Mann Kendall trend test
    Tstatistic=u;  # eq 14.10 in Haan (2002)
    
    w <- pnorm(abs(u), mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
    pvalue <- 2*(1-w)
    pstring <- as.character(round(pvalue,digits=6))
    strP <- paste('p = ',pstring,', Mann-Kendall test (H0: no trend)',sep='')
    
    #================== NON-PARAMETRIC TREND LINE
    #Haan 2002, p. 345

    # matrix with (mx-1) rows, each of which is duped vector (mx-1):1
    N1 <- N-1
    k <- N1:1
    K <- matrix(k,ncol=N1,nrow=N1,byrow=TRUE)
    rm(k)

    # matrix with (mx-1) cols, each of which is duped vector (mx-2):0
    N2 <- N-2
    j <- N2:0
    J <- matrix(j,ncol=N1,nrow=N1,byrow=FALSE)
    rm(j)
    
    # More matrices
    H <- K-J
    L<-upper.tri(H,diag=TRUE)
    Q1 <- D/(lower.tri(H,diag=TRUE)*H)
    
    b <- median(Q1,na.rm=TRUE) # slope
    tmed <-median(yrx)

    a <- median(x)-b*tmed
    
    # Build a string for trend line
    if (b<0){
      bb <- ' - '
    } else {
      bb <- ' + '
    }
    eqn1 = paste('y = ', sprintf('%g',a), bb, sprintf('%g',b), 't, ','trend line\n',sep='')
    
    # Generate the prediction of x from the non-parametric trend model
    xhat <- a + b*yrx  # prediction by trend model

    # Generate a detrended version of x; shift that to have the same median as x
    xdetrended <- x-xhat # before shift
    xmed1 <- median(xdetrended)
    xmed2 <- median(x)
    d1 <- xmed1-xmed2
    xdetrended <- xdetrended-d1
    
    #============= OPTIONAL PLOT, TO BE RETURNED AS A PNG
    # see these inputs:
    # kopt ....whether want plot
    # kplot... which plots
    
    if (kopt[1]==2){
      # You want plots

      # Next setting apply to plot with trend line in it, regardless of kplot setting
      Tit1 <- textPlot[1]; xlab <- textPlot[2]; ylab <- textPlot[3]
      Tit1 <- paste(Tit1,strFL)
      Tit2 <- textPlot[4]
      yhi <- max(c(max(x),max(xhat),max(xdetrended))); ylo = min(c(min(x),min(xhat),min(xdetrended)))
      yexpand = 0.10*(yhi-ylo)
      ylims <- c(ylo,yhi+yexpand)
      
      yrgo1 <- yrx[1]; yrsp1 <- yrx[length(yrx)]
      if (kplot==3){
        # 1x2 with seris and trend line at top, detrended and original series at bottom
        
        #--- Build figure png filename
        fileOut <- paste(outputDir,'mannken1_F1a','.png',sep="")
        zm <- c(xmed2,xmed2) # horizontal line will be at median x, which should also be the
        # median of xdetrended

        
        # Figure size and margins
        png(filename=fileOut, width = 960, height = 480)
        layout.matrix <- matrix(c(1,2), nrow = 2, ncol = 1)
        layout(layout.matrix,heights=1,widths=1) 
        
        # Top plot
        par(mar=c(4.0,8,2,8),cex.axis=1.1, cex.lab=1.5, cex.main=1.3)
        plot(yrx,x,type="b",pch=1,col="blue",xlim=c(yrgo1,yrsp1),xlab=xlab,
             ylab=ylab, main=Tit1,ylim=ylims)
        lines(yrx,xhat,type="l",col="red") # non-parametic-fit trend line
        abline(h=zm,lty=2,col='#808080') # dash gray
        
        # annote test info
        f1 <- 1.1
        ySep <- yexpand # y-separate from above line of text
        text(yrx[1],ylims[2],eqn1,adj=c(0,1),cex=1.0) # line eqn
        text(yrx[1],ylims[2]-f1*ySep,strP,adj=c(0,1),cex=1.0) # pvalue for MK test
        # conditional annotated text on variance inflation factor
        if (Lacf){
          # You asked to look at the acf of residuals of a least-squares-fit straight line fit to 
          # the time series. If the residuals from this line have positive lag-1 autocorrelation signficant 
          # at p<0.05 by a one-tailed test, vif will be computed and vif>1.0. If no significant
          # lag-1 autocorrelation, vif is set to 1.0.
           strVIF <- paste('VIF=',sprintf('%g',vif))
           text(yrx[1],ylims[2]-2*f1*ySep,strVIF,adj=c(0,1),cex=1.0) # variance inflation factor
        }
        
        # Bottom plot
        par(mar=c(4.0,8,2,8),cex.axis=1.1, cex.lab=1.5, cex.main=1.3)
        plot(yrx,xdetrended,type="b",pch=1,col="blue",xlim=c(yrgo1,yrsp1),xlab=xlab,
             ylab=ylab, main=Tit2,ylim=ylims)
        abline(h=zm,lty=2,col='#808080') # dash gray
        dev.off()
      } else if ((kplot==2) |(kplot==1)){
        # 1x1  with the series and fitted trend line, with special naming of figure file for TRISH
        
        #--- Build figure png filename
        if (kplot==2){
          # TRISH-special
          zm <- c(0,0)
          FigNumber <- NextFigNumber+jFigAdd   # for naming this png
          if (FigNumber<10){
            fileOut <- paste(outputDir,'Figure0', as.character(FigNumber),'-AnalysisResiduals2.png',sep="")
          } else {
            fileOut <- paste(outputDir,'Figure', as.character(FigNumber),'-AnalysisResiduals2.png',sep="")
          }
        } else {
          zm <- c(xmed2,xmed2)
          # Most general
          if (outputDir=="Null"){
          } else {
            fileOut <- paste(outputDir,'mannken1_F1a','.png',sep="")
          }
        }
        #--- Build time plot of time series with trend line and annotation
        if (outputDir=="Null"){
        } else {
          #fileOut <- paste(outputDir,'mannken1_F1a','.png',sep="")
          png(filename=fileOut, width = 960, height = 480)
          par(mar=c(5,6,2,2),cex.axis=1.5, cex.lab=1.5, cex.main=1.5)
        }
        plot(yrx,x,type="b",pch=1,col="blue",xlim=c(yrgo1,yrsp1),xlab=xlab,
             ylab=ylab, main=Tit1,ylim=ylims)
        lines(yrx,xhat,type="l",col="red") # non-parametic-fit trend line
        abline(h=zm,lty=2,col='#808080') # dash gray
        
        # annote test info
        ySep <- yexpand/2 # y-separate from above line of text
        
        # Simplify annotation if outputDir="Null"
        if (outputDir=="Null"){
        } else {
          text(yrx[1],ylims[2],eqn1,adj=c(0,1),cex=1.2) # line eqn
          text(yrx[1],ylims[2]-ySep,strP,adj=c(0,1),cex=1.2) # pvalue for MK test
        }
        # conditional annotated text on variance inflation factor
        if (Lacf){
          # You asked to look at the acf of residuals of a least-squares-fit straight line fit to 
          # the time series. If the residuals from this line have positive lag-1 autocorrelation signficant 
          # at p<0.05 by a one-tailed test, vif will be computed and vif>1.0. If no significant
          # lag-1 autocorrelation, vif is set to 1.0.
          strVIF <- paste('VIF=',sprintf('%g',vif))
          
          if (outputDir=="Null"){
            nullTxt  <- paste(strP,'VIF=',sprintf('%g',vif))
            text(yrx[1],ylims[2],nullTxt,adj=c(0,1),cex=1.0)
          } else {
            text(yrx[1],ylims[2]-2*ySep,strVIF,adj=c(0,1),cex=1.2) # variance inflation factor
          }
         }
        if (outputDir=="Null"){
        } else {
          dev.off()
        }
      } else {
        stop ('kplot must be 1 or 2')
      }
    }
  }
  
  #--- BUILD OUTPUT LEST
  
  creation <- 'List created by function mannkenn'
  definitions <- 'See opening comment section of the creation function'
  dateCreated <- Sys.Date()
  What <- list("creation"=creation,"definitions"=definitions,"dateCreated"=dateCreated)
 
  Output <- list("What"=What,"statistic"=Tstatistic,"AnalysisPeriod"=strFL,"Lacf"=Lacf,"Ladjusted"=Ladjusted,
                 "vif"=vif,"pvalue"=pvalue,"ngrp"=ngrp,"nties"=nties,"Lflag"=Lflag,"b"=b,"a"=a,
                 "equation"=eqn1,"X"=X,"xhat"=xhat,"xdetrended"=xdetrended)
  return(Output)
}

