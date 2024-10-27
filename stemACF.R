stemACF <- function(Din){
  # Stem plot for ACF, ala Matlab
  # D. Meko
  # Revised 2022-05-17
  #
  # Stem plot for ACF, ala Matlab.
  #
  #---IN
  #
  # Din: list with elemets as follows:
  #   x [vector] lags for plot (assume start at lag 0, where r(0)==1)
  #   nsize [integer] numbe observations in time series on which acf was computed
  #   y [vector] acf at those lags, x
  #   kAlpha [integer] alpha level for CI (two-sided test of H0 that acf(k)=0)
  #   =1 0.05 level (95% CI)
  #   =2 0.01 level (99% CI)
  #   FigNumber [integer]; png with name Figure??.png will be written to outputDir,
  #     where ?? (e.g., 1) is the number of the figure. 
  #   outputDir [char] output directory (e.g., )
  #   linecol1, 2, 3: line colors for stems, zero line, and CI
  #   textPlot [list]  text for xlabel, ylabel, title, and optional upper left annotation
  #       If textPlot[4]=='null', nothing is annotated; also, textPlot[5] is either 'null' or some
  #       string, such as '-AnalysisResiduals3', which is to be built into the filename. Thus, might
  #       have output "Figure12-AnalysisResiduals3.png"
  #
  #
  #---OUT
  #
  #---NOTES
  #
  # Little checking of inputs for proper class, etc. 
  # Motivation:I got the idea for this plot function from this page:
  #   https://www.r-bloggers.com/2009/11/matlab-style-stem-plot-with-r/
  # textPlot: typically something like c('Lag k (yr)', 'r(k))','ACF of ...',xxx)
  #   where xxx might be some string giving additional information, such as DW statistic
  # txtPlot[5]: set this to 'Null' for general use (outside of TRISH)
  
  #---HARD CODE
  
  ylims <- c(-1.05,1.05) # acf can range from -1 to +1

  #---UNLOAD INPUT
  
  x <- Din$x; y <- Din$y; kAlpha <- Din$kAlpha; FigNumber <- Din$FigNumber
  outputDir <- Din$outputDir; nsize <- Din$nsize
  linecol1 <- Din$linecol1; linecol2 <- Din$linecol2; linecol3 <- Din$linecol3
  
  xlab1 <- Din$textPlot[2]
  ylab1 <- Din$textPlot[3]
  Tit1 <- Din$textPlot[1]
  txtAnn <- Din$textPlot[4]
  txtFnm <- Din$textPlot[5]
  
  
  #---COMPUTE CRITICAL R FOR CONFIDENCE INTERVAL
  
  if (kAlpha==1){
    rcrit <- 1.96/sqrt(nsize)
  } else if (kAlpha==2){
    rcrit <- 2.5758/sqrt(nsize)
  } else {
    stop('stem1 accepts kAlpha of either 1 or 2')
  }
  
  
  #---SET UP FIGURE
  
  #--- Build figure png filename
  if (FigNumber<10){
    if (txtFnm=='Null'){
      fileOut <- paste(outputDir,'Figure0', as.character(FigNumber),'.png',sep="")
    } else {
      fileOut <- paste(outputDir,'Figure0', as.character(FigNumber),txtFnm,'.png',sep="")
    }
  } else {
    if (txtFnm=='Null'){
      fileOut <- paste(outputDir,'Figure', as.character(FigNumber),'.png',sep="")
    } else {
      fileOut <- paste(outputDir,'Figure', as.character(FigNumber),txtFnm,'.png',sep="")
    }
  }
  
  # Figure size and margins
  png(filename=fileOut, width = 960, height = 480)
  par(mar = c(5.1, 4.5, 5.1, 2.1),cex.axis=1.1, cex.lab=1.5, cex.main=1.3,cex=1.2)
  
  
  #---MAKE PLOT

  plot(x,y,pch=16,type='p',col=linecol1,xlab=xlab1,ylab=ylab1,
       main=Tit1,ylim=ylims)
  for (i in 1:length(x)){
    lines(c(x[i],x[i]), c(0,y[i]),col=linecol1)
  }
  # horiz lines at 0 and upper and lower CI
  abline(h=0,lty=1,col=linecol2) 
  abline(h=rcrit,lty=5,col=linecol3) 
  abline(h=-rcrit,lty=5,col=linecol3) 
  
  # optional annotation 
  if (txtAnn=='null'){
    # no action needed
  } else {
    text(1.1,ylims[2],txtAnn,adj=c(0,1),cex=1.2)
  }
  dev.off()
  
  Output <-NA
  return(Output)
}