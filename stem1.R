stem1 <- function(Din){
# Stem plot for ACF, ala Matlab
# D. Meko
# Revised 2022-05-17
#
# Stem plot, useful for plotting ACF with CI
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
#   nextFigNumber [integer]; png with name graph-??.png will be written to outputDir,
#     where ?? (e.g., 1) is the number of the figure
#   outputDir [char] output directory (e.g., )
#   linecol1, 2, 3: line colors for stems, zero line, and CI
#   textPlot [list]  text for xlabel, ylabel, title, and upper left annotation
#
#---NOTES
#
# I got the idea for this plot function from this page:
# https://www.r-bloggers.com/2009/11/matlab-style-stem-plot-with-r/
  


#---UNLOAD INPUT

x <- Din$x; y <- Din$y; Kalpha <- Din$kAlpha; nextFigNumber <- Din$NextFigNumber
outputDir <- Din$outputDr; nsize <- Din$nsize
linecol1 <- Din$linecol1; linecol2 <- Din$linecol2; linecol3 <- Din$linecol3

xlab1 <- Din$xtextPlot[1]
ylab1 <- Din$xtextPlot[2]
Tit1 <- Din$xtextPlot[3]
txtAnn <- Din$textPLot[4]

#---COMPUTE CRITICAL R FOR CONFIDENCE INTERVAL

if (kAlpha==1){
  rcrit <- 1.96/sqrt(nsize)
} else if (kAlpha==2){
  rcrit <- 2.5758/sqrt(nsize)
} else {
  stop('stem1 accepts kAlpha of either 1 or 2')
}


#---SET UP FIGURE

#--- Build graphics png filename
FigNumber <- NextFigNumber   # for naming this png
if (FigNumber<10){
  fileOut <- paste(outputDir,'graph-0', as.character(FigNumber),'.png',sep="")
} else {
  fileOut <- paste(outputDir,'graph-', as.character(FigNumber),'.png',sep="")
}

# Figure size and margins
png(filename=fileOut, width = 960, height = 480)
par(mar = c(5.1, 4.5, 5.1, 2.1),cex.axis=1.1, cex.lab=1.5, cex.main=1.3)


#---MAKE PLOT

ylims <- c(-1.05,1.05)

plot(x,y,pch=16,type='p',col=linecol1,xlab=xlab1,ylab=ylab1,
     main=Tit1)
for (i in 1:length(x)){
  lines(c(x[i],x[i]), c(0,y[i]),col=linecol)
}
abline(h=0,lty=1,col=linecol2) 
abline(h=rcrit,lty=1,col=linecol3) 
abline(h=-rcrit,lty=1,col=linecol3) 

dev.off()
}