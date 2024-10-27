NashSutt<- function(y, yh, yry, ycMean, kplot, outputDir,gFileType) {
  # Nash-Suttcliffe efficiency 
  # D Meko
  # last revised 20230726
  #
  # IN
  # y: vector of observed
  # yry: year vector for y and yh
  # yh: vector or reconstructed (or predicted) for same observations
  # ycMean: (1x1): calibration-period mean of y
  # kplot (1x1): plot option
  #   1= plot on current device
  #   2= plot a graphics file of with suffix "gFileType" to folder "outputDir"
  #   3= skip plot
  # outputDir: folder to write plot file to (ignored unless kplot=2)
  # gFileType: type of graphics file (e.g., "png") -- ignored unless kplot=2
  #
  # OUT
  # Output [named list] 
  #   NSE (1x1)r Nash-Suttcliffe efficiency
  #   RE (1x1)r reduction of error statistic
  #   rTest: list with result of correlation (#Pearson) test for IV period:
  #     $estimate = correlation
  #     $parameter = df
  #     $p,vaue = pvalue for test
  #   e: vector of errors (observed minus predicted)
  #   yMean: mean y (for Nash-Suttcliffe period)
  #   ycMean: mean y for calibration period on which prediction based (echo of input)
  #   SSE: sum of squares of e
  #   SSEnull1: sum of squares of departures (y-yMean)
  #   SSEnull2: sum of squares of departures (y-ycMean)
  #
  # Notes
  #
  # The NSE and RE are similarly computed from ratios of sum-of-squares terms. NSE uses the validation period
  # mean as the null prediction, giving SSEnull1. RE uses the calibration period means as the null prediction,
  # giving SSEnull2. Equations are
  #
  # SSE = 1 - (SSE/SSEnull1)
  # RE = 1 - (SSE/SSEnull2)
  library("pracma") # needed for emulation of Matlab "backslash" operator through
  # QR decomposition
  
  source(paste(code_dir,"LeaveOut.R",sep="")) # form pointer matrix for leave-m-out cross-validation

    #--- CHECK INPUT
  # y and yh should be same-length numeric vectors
  
  if (!is.vector(y) || (!is.vector(yh))){
    stop('y and yh must be vectors')
  } else {
  }

  if (!is.numeric(y) || (!is.numeric(yh))){
    stop('y and yh must be numeric')
  } else {
  }

  L <- length(y) == length(yh)
  if (!L){
    stop('y and yh must be same length')
  }

  #--- COMPUTE STATISTICS
  
  yMean <- mean(y) # verification-period observed mean (null prediction for NSE)
  e <- y-yh # validation period errors
  e1 <- y-yMean # errors for null prediction of validation mean
  e2 <- y-ycMean # errors for null prediction of calibration mean
  
  SSE <- sum(e*e) # sum of squares, recon errors
  SSEnull1 <- sum(e1*e1) # sum of squares, reconstruction consisting of validation mean each year
  SSEnull2 <- sum(e2*e2) # sum of squares, reconstruction consisting of calibration mean each year
  
  NSE <- 1 - (SSE/SSEnull1) # Nash-Suttcliffe efficiency
  RE <- 1 - (SSE/SSEnull2) # Reduction of error statistic
  
  rTest <- cor.test(y,yh)
  r <- rTest$estimate # Pearson r
  df <- rTest$parameter
  p <- rTest$p.value

  
  #--- PLOT
  
  # Expand right side of clipping rect to make room for the legend
  par(xpd=T, mar=par()$mar+c(0,0,0,6))

    # ylim
  yhi <- max(c(max(yh),max(y)))
  ylo <- min(c(min(yh),min(y)))
  ylims <- c(ylo,yhi) 
  yhi <- yhi+0.1*diff(ylims)
  ylims <- c(ylo,yhi) 
  xlims <- c(yry[1]-1,yry[length(yry)]+1)
  
  str_r <- paste('NSE=',sprintf('%5.2f',NSE),';  RE=',sprintf('%5.2f',RE),'; r=',sprintf('%5.2f',r), ' (N=',as.character(df),
                 ', p=',sprintf('%6g',p),')',sep='')
  str_tit <- paste('Independent verification, ',as.character(yry[1]),'-',
                   as.character(yry[length(yry)]),'\n',str_r)
  plot(yry,y,type='b',ylim=ylims,xlim=xlims,
       xlab='Year',ylab='Data value',main=str_tit) # main plot of obs
  lines(yry,yh,col='red') # Line for reconstruction
  points(yry,yh,col='red',pch=2)
  # Null prediction -- the observed mean for indep period
  lines(xlims,c(yMean,yMean),col='black')
  # Null prediction by calib mean
  lines(xlims,c(ycMean,ycMean),col='black',lty='dashed',lwd=2)
  #abline(h=ycMean,col='black',lty='dashed')
  # The reconstructed mean for independent period
  lines(xlims,c(mean(yh),mean(yh)),col='red')
  #abline(h=mean(yh),col='red')
  
  # In legend below, "26" ignored as a plot character. This is not an error.
  legend(xlims[2]+1,ylims[2], legend=c("Observed","Recon","ObsMeanV","ObsMeanC","RecMeanV"),
         col=c("black","red", "black","black","red"),
         lty=c(1,1,1,2,1),pch=c(1,2,26,26,26),lwd=c(1,1,1,2,1), cex=0.9)

  # Restore default clipping rect
  par(mar=c(5, 4, 4, 2) + 0.1)
  
  
  #--- OUTPUT
  
  Output <- list("NSE"=NSE,"RE"=RE,"rTest"=rTest,"Errors"=e)
  return(Output)
}