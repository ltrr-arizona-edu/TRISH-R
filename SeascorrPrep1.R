# Prepare climate and tree-ring data for call to treeclim/seascorr
# 
# D Meko
# Last revised 20250526
#
# Written because needed to convert format of monthly P and T data from what
# came from KNMI explorer to what treeclim's seascorr function wants
#
# This script prepares climate and tree-ring data for call to
#   treeClim_seascorr.Rmd.  Seascorr accepts the monthly climate data
#   as a data frame with years
#   as the row names. Two time series of monthly climate data are required, a
#   primary and secondary. Seascorr can accept that data as a 4-col matrix or 
#   frame, with cols: year, month, primary climate variable, secondary ...
# This script assumes you have in your working folder:
#   1) A tab-sep or space-sep file of each of the climate variables, with year
#     as col 1 and jan-Dec data as cols 2-13
#   2) A tab-sep or space-sep file of tree-ring chronologies with year as col-1
#     and one or more chrons as the remaining columns. You edit the script to say
#     which of the chronologies to use
#
# The script finally saves the two data frames for later loading by the script 
# that calls seascorr. Last statement saves the data as 'SeascorrTahoeInput.RData'. You 
# would revise output filename as well as the input filenames, as needed, to be consistent
# with your data.
#
# Search for "TailorMe" to find those lines that need to be tailored to your data

source("c13toc3.R") # converts 13-col to 3 col matrix

#--- Read input climate data files 
V1 <- read.table('Pmonthly.dat',sep="",header=FALSE) # read 13-col space-sep file TailorMe
V2 <- read.table('Tmonthly.dat',sep="",header=FALSE) # read 13-col space-sep file TailorMe

# Read table of tree-ring chronologies; pull desired, store as data frame
T1 <- read.table('ChronsStdTahoe5.txt',sep='\t',header=TRUE) # TailorMe
T2 <- T1[,-1] # strip off year col
rownames(T2) <- T1[,1]
rnms  <- rownames(T2)
c2 <- colnames(T2)
jwant <- 4 # which of the data columns will supply desired chron TailorMe
T2 <- as.data.frame(T2[,jwant])
c2 <- c2[jwant]
colnames(T2) <- c2
A <- T2 #  This is the data frame with the desired chronology
rownames(A) <- rnms
L <- complete.cases(A)
A <- as.data.frame(A[L,])
rownames(A) <- rnms[L]
colnames(A) <- c2
rm(c2,L,T2,T1)


#--- PUT 13-COL CLIMATE DATA INTO MATRICES, AND TRIM TO COVER SAME YEARS

V1 <- as.matrix(V1)
yrV1 <- V1[,1]
V2 <- as.matrix(V2)
yrV2 <- V2[,1]
mV1 <- nrow(V1)
mV2 <- nrow(V2)

yrgo <- max(yrV1[1],yrV2[1])
yrsp <- min(yrV1[mV1],yrV2[mV2])

L <- (yrV1 >= yrgo) & (yrV1 <= yrsp)
yrV1 <- yrV1[L]
V1 <- as.matrix(V1[L,])

L <- (yrV2 >= yrgo) & (yrV2 <= yrsp)
yrV2 <- yrV2[L]
V2 <- as.matrix(V2[L,])           


#--- CONVERT PRIMARY AND SECONDARY CLIMATE SERIES TO 3-COL AND COMBINE THEM

X1 <- c13toc3(V1) # 3-col primary clim variable
X2 <- c13toc3(V2) # 3-col secondary clim variable
X <- cbind(X1,X2[,3])


#--- CONVERT 4-COL CLIMATE MATRIX INTO A DATAFRAME 

X <- as.data.frame(X)
colnames(X) <- c('Year','Month','P','T') # TailorMe


#--- STATUS
#
# X [data frame] climate data: year, month, primary secondary
# A [data frame] tree ring chronology: years as rownames

save(X,A,file='SeascorrTahoeInput.RData') # TailorMe