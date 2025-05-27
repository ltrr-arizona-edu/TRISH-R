# Prepare climate and tree-ring data for call to treeclim/seascorr
# 
# D Meko
# Last revised 20250527
#
# Written because I needed to convert format of monthly P and T downloaded from 
# KNMI Explorer to format readable as input to package treeclim's seascorr
# function.
#
# This script prepares climate and tree-ring data for call to 
# treeClim_seastcorr.Rmd. Monthly monthly climate data must be a data frame
# with years as the row names. Two time series of monthly climate data are required: a
# primary and secondary. Seascorr reads that data as a 4-col matrix or with 
# cols: year, month, primary climate variable, secondary variable.
#
# This script assumes you have in your working folder:
#   1) A tab-sep or space-sep file of each of the climate variables, with year
#     as col 1 and Jan-Dec data as cols 2-13
#   2) A tab-sep or space-sep file of tree-ring chronologies with year as col-1
#     and one or more chronilogies as the remaining columns. You edit the script 
#     to tell which of the chronologies to use
#
# The script finally saves the two data frames for later loading by the script 
# TreeClim_seascorr.Rmd. Last statement saves the data as 'SeascorrTahoeInput.RData'. 
# You should edit to match your desired filename for ".RData" file, and should also
# revise the input filenames to be consistent with your data. 
#
# To assist in making those edits, I have added a comment containing the text "TailorMe"
# before lines needing revision. Search for "TailorMe" to find those lines.

# Edit following two lines to match where you have stored user-written functions (UWFs)
# and input files
code_dir <- "/home/dave/GitWork/TRISH-R/" # path to UWFs, TailorMe
data_dir <- "/home/dave/GitWork/TRISH-R/" # path to input data, TailorMe

source(paste(code_dir,'c13toc3.R',sep=''))  # converts 13-col to 3 col matrix



#--- Read input climate data files as data.frame 
# Edit next two lines as needed for path & filname
V1 <- read.table(paste(data_dir,'Pmonthly_seascorr.dat',sep=''),sep="",header=FALSE) # read 13-col space-sep file TailorMe
V2 <- read.table(paste(data_dir,'Tmonthly_seascorr.dat',sep=''),sep="",header=FALSE) # read 13-col space-sep file TailorMe

# Read table of tree-ring chronologies; pull desired, store as data frame
T1 <- read.table(paste(data_dir,'ChronsTahoe5.txt',sep=''),sep='\t',header=TRUE) # TailorMe
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

save(X,A,file=paste(data_dir,'SeascorrTahoeInput.RData',sep='')) # TailorMe