# chron1_UF2025b.R
# Chronology development sample script
# D Meko
# Last revised 2025-01-10
#
# Loops over rwl files (see nms1) and develops site chronologies with and without 
# variance stabilization. In the specified output folder (see outputDir), you will 
# find many files, all beginning with a code that identifies the chronology. 
# Plots show the differences in the chronologies with and without variance
# stabilization. A summary plot shows the range shrinkage due to variance
# stabilation. Range shrinkage defines as the ratio of range of time series after
# variance stabilization to range before variance stabilization. 
#
# chron1_UF2025b.R was modified from script chron1_ba1.R, written to standardize a
# set of tree-ring collections that have 4 data type (not just total-width)
#
# .crn --- the standard site chronology
# VarStab.crn -- the variance-stabilized chronology
# .txt -- a tab-separated table with, among other things, the windowed EPS
# .png -- chronology plot (dplR style), of standard site chronology with shaded sample depth
# VarStab.png -- time plots of standard chronology, before and after variance stabilization


######## LOAD LIBRARIES

library(dplR)
library(treeclim)
library(data.table)

######## CLEANING AND PREP

rm(list = setdiff(ls(), lsf.str())) # Remove existing variables except functions

# Hard code some system folders (need to tailor to your laptop)
# For ultimate simplicity you can make all of these folders the same!
code_dir <- "/home/dave/GitWork/TRISH-R/" # any user-written functions here
outputDir <-"/home/dave/AAA_sahar/test_out/" # any output files go here
dataInDir <- "/home/dave/AAA_sahar/data/" # rwl input files here
workingDir <- "/home/dave/AAA_sahar/" # folder for this R project 

# set current working directory
setwd(workingDir)

# Hard code lists of filenames.
nms1 <- c('Earth001','Earth002','Earth003','Earth004') # codes for sites

# Detrending spline and spline smoother for smoothed time series plot
NyrSpline <-100 # use NyrSpline-year spline for detrending
NyrSmooth=10; # this spline for smoothing annual serie

# File name for bar plot showing change in variance due to variance adjustment
flnm4 <- "RangeShrink.png"


######## SOURCE NEEDED FUNCTIONS

source(paste(code_dir,"emssgUNH.R",sep="")) # writes error message to specified output folder


######## LOOP OVER SITES 

nsites  <- length(nms1)
F1 <- matrix(NA,nrow=nsites,ncol=1) # to hold ratio of range (annual data) after to before variance-stabilization
F1s <- matrix(NA,nrow=nsites,ncol=1) # to hold ratio of range (smoothed) after to before variance-stabilization
rownames(F1) <- nms1
colnames(F1) <- "Ratio"
rownames(F1s) <- nms1
colnames(F1s) <- "Ratio"

for (n in 1:nsites){
  nm1 <- nms1[n]
  pf  <- paste(dataInDir,nm1,'.rwl',sep='')
  a1 <- read.rwl(pf,format="tucson") # read in the rwl object as  data frame
  a2 <-detrend(a1, y.name = names(a1), make.plot = FALSE,
               method = "Spline",
               nyrs = NyrSpline, f = 0.5, pos.slope = FALSE,
               constrain.nls = c("never", "when.fail", "always"),
               verbose = FALSE, return.info = FALSE,
               difference = FALSE)
  
  
  # Build path/filenames for ouput
  fnm1 <- paste(nm1,'',sep='') 
  fnm1s <- paste(fnm1,'_VarStab',sep='');
  fnm1 <- paste(outputDir,fnm1,sep='')
  fnm1s <- paste(outputDir,fnm1s,sep='')

  # Call function to compute standardized chronology, without and with variance-
  # stabilization
  a3 <- chron(a2, biweight = TRUE, prewhiten = FALSE) 
  a3s <- chron.stabilized(a2, winLength=50, biweight = TRUE, running.rbar = FALSE)
  
  
  #---WRITE CRN TO OUTPUT FOLDER
  
  pf <- paste(fnm1,'.crn',sep='')
  write.crn(a3, pf, header = NULL, append = FALSE) # standard site chron
  pf <- paste(fnm1s,'.crn',sep='')
  write.crn(a3s, pf, header = NULL, append = FALSE) # var-stab. standard site chron
  
  #--- CHRONOLOGY RUNNING STATISTICS
  
  s1 <- rwi.stats.running(a2, ids = NULL, period = "max",
                          method = "spearman",
                          prewhiten=FALSE,n=NULL,
                          running.window = TRUE,
                          window.length = 50,
                          window.overlap = 25,
                          first.start = NULL,
                          min.corr.overlap = 30,
                          round.decimals = 3,
                          zero.is.missing = TRUE)
  s2 <- s1[,c(1,3,4,5,14,15)] # reduced table of running statistics
  pf <- paste(fnm1,'.txt',sep='')
  
  fwrite(s2, file = pf, append = FALSE, quote = "auto",
         sep = "\t", sep2 = c("","|",""),
         eol = if (.Platform$OS.type=="windows") "\r\n" else "\n",
         na = "", dec = ".", row.names = FALSE, col.names = TRUE,
         qmethod = "double",
         logical01 = getOption("datatable.logical01", FALSE),  # due to change to TRUE; see NEWS
         scipen = getOption('scipen', 0L),
         dateTimeAs = "ISO",
         buffMB = 8L, 
         showProgress = getOption("datatable.showProgress", interactive()),
         compress = "none",
         yaml = FALSE,
         bom = FALSE)
  
  #write.csv(s2, pf, row.names=TRUE)
  #---  PLOT CHRONOLOGY   
  
  # Will be adding line for variance stabilized chron. Need year and value.
  yrx <- as.numeric (rownames(a3s))
  x <- a3s$vsc
  # No variance stabilization
  fnmPlot <- paste(fnm1,'.png',sep='')
  png(filename=fnmPlot, width = 960, height = 480)
  par(mar=c(5,5,4,5))
  # Create the layout
  crn.plot(a3, add.spline = TRUE, nyrs = 10, f = 0.5,
           crn.line.col='grey50',spline.line.col='darkred',
           samp.depth.col='grey90',
           samp.depth.border.col='grey80',
           crn.lwd=1,spline.lwd=1.75,
           abline.pos=1,abline.col='black',
           abline.lty=1,abline.lwd=1,
           xlab="Time",ylab="Index")
  dev.off()
  
  
  # SPLINE SMOOTHED VERSIONS OF CHRONOLOGY
  
  y1 <- ffcsaps(a3$std,nyrs=NyrSmooth)  # original
  y2 <- ffcsaps(x,nyrs = NyrSmooth) # variance stabilized
  
  # TIME PLOTS COMPARING STANDARD CHRONOLOGY WITHOUT AND WITH VARIANCE STABILIZATION
  #
  # 2x1, with annual in 1, 10-year spline in 2
  
  fnmPlotVS <- paste(fnm1s,'.png',sep='')
  
  # Figure size and margins
  png(filename=fnmPlotVS, width = 960, height = 480)
  layout.matrix <- matrix(c(1,2), nrow = 2, ncol = 1)
  layout(layout.matrix,heights=c(1,1),widths=2) 
  
  # Limits, keep same for A, B
  xLo <- min(yrx); xHi <- max(yrx); yLo <- min(a3$std); yHi <- max(a3$std)
  
  # Create empty plot, top - A
  Title1 <- paste(nm1,': effect of variance-stabilization (red) on site chronology (black)',
                  sep='')
  par(mar = c(2.1,4.1,2.1,2.1),cex.axis=1.1, cex.lab=1.5, cex.main=1.3)
  plot(1,type='n',xlim=c(xLo,xHi),ylim=c(yLo,yHi),
       xlab='Year', ylab='Index',main=Title1)
  
  # Add lines
  lines(yrx,a3$xxxstd,col='black')
  lines(yrx,x,col='red')
  
  # Create empty plot, bottom, B
  Title1 <- paste(nm1,': smoothed by 10-year spline',
                  sep='')
  par(mar = c(2.1,4.1,2.1,2.1),cex.axis=1.1, cex.lab=1.5, cex.main=1.3)
  plot(1,type='n',xlim=c(xLo,xHi),ylim=c(yLo,yHi),
       xlab='Year', ylab='Index',main=Title1)
  
  # Add lines
  lines(yrx,y1,col='black')
  lines(yrx,y2,col='red')
  
  dev.off()
  
  f1 <- c(diff(range(x))/diff(range(a3$std)) ,  diff(range(y2))/diff(range(y1)))
  F1[n,1]<- f1[1]; F1s[n,1]<- f1[2]
  F1 <- round(F1,2); F1s  <- round(F1s,2)
}

# BAR CHARTS OF % COMPRESSION OF RANGE DUE TO VARIANCE STABILIZATION
#
# 1X2, annual in 1, smoothed in 2

pf <- paste(outputDir,flnm4,sep='')
# Figure size and margins
png(filename=pf, width = 960, height = 480)
layout.matrix <- matrix(c(1,2), nrow = 1, ncol = 2)
layout(layout.matrix,widths=c(1,1),heights=2) 
ylims <-c(0,1.1)

# annual
par(mar = c(4.1,5.1,4.1,4.1),cex.axis=1.1, cex.lab=1.5, cex.main=1.3)
barplot(t(F1), main="Annual",
        xlab="Chronology", ylab="Ratio of Range",
        col=c("green"),
        ylim=ylims, beside=TRUE)
legend(x="topright",legend=colnames(F1),fill=c("green"),
       inset=c(0.30,0))
       
# smoothed
par(mar = c(4.1,5.1,4.1,4.1),cex.axis=1.1, cex.lab=1.5, cex.main=1.3)
barplot(t(F1s), main="Smoothed",
        xlab="Chronology", ylab="Ratio of Range",
        col=c("green"),
        ylim=ylims, beside=TRUE)
legend(x="topright",legend=colnames(F1),fill=c("green"),
       inset=c(0.30,0))
dev.off()
