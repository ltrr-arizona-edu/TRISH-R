# chron2_UF2025b.R
# Chronology time series matrix development script
# D Meko
# Last revised 2025-01-17
#
# Reads a file of crn filenames; loops over the filenames, reading in site
# chronologies and sample sizes and building time series matrices of those, 
# with year as column 1.
#
# Input specs allow control over the time coverage of the resulting matrix and
# which chronologies it contains. Writes two tab-separated-value (tsv) files of
# time series chronologies, and corresponding two files of sample size.
#
# One version of tsv file has all chronologies or sample sizes covering interval
# from the first year with data in any series to the most recent year with
# data in any series. The other version of tsv file has the matrix truncated 
# to start at specified year.  Output filenames are coded like this, for example:
#
#   EarthChronsFull.txt, SS_EarthChronsFul.txt
#   EarthChrons1466.txt, SS_EarthChrons1466.txt
#
# where 1466 in this case is the specified start year for the "cut" version of
# chronologies. 
#
# You have control through edited lines (search this file for occurrences of
# "Edit next) of following:
#  1) start year of truncated or "cut" version of output time series matrix
#  2) whether missing values are coded as "NaN" or "NA"
#  3) whether the output files of time series matrix should have a leading row with
#     column number (0 above the Year column) indicating sequential number of
#     chronology in the matrix.
#
# Besides the output tsv files, you will find a list "ChronData" in your 
# environment after running this script. ChronData is a list with elements:
#  Xfull, Xcut, Sfull, Scut: array matrices of the chronologies and sample sizes,
#     with year as column 1
#  namesFull, namesCut:  the corresponding series codes, or names
#
#
#--- NOTE
#
# The option for a leading row with the column number (0  1 2 ...) is included
# because web-based tool requires that row for an uploaded data matrix. Likewise
# TRISH requires "NaN" for the missing value code. 
#
# The option for "cut" or "truncated" time series matrix output is for 
# compatibility of running ReconAnalog.R reconstruction script standalong outside
# of TRISH. In that case, the all series in the matrix are required to have
# data in the first year of the matrix. For TRISH, the full matrix can be uploaded
# and TRISH interactively does any required row and column trimming.
#
#------ PREPARATION STEPS FOR TESTING ON SAMPLE DATA
#
# 1) store the four provided sample crn files
#     (Earth001_VarStat.rwl,  ...002, 003, 004)
#   in a previously empty folder, your "dataInDir"
# 2) copy this script to your R project folder, and function crn2tsv.R to
#   your "codeDir" (see below)
# 3) RStudio: start it, and open your R project or create new project
# 4) Install R libraries listed under "LOAD LIBRARIES" in your "codeDir"
# 5) search (using ctl-F) for all occurrences of the string "Edit next" in
#   this R script and edit any paths accordingly
# 8) try running this script; if successful, two new files will appear in
#   your specified (see below) "outputDir)
#
# To run on your own crn files, follow directions above, except will have to
# edit lines dealing with paths/filenames and with the truncation start year "yrGo",
# and have to comment out two lines to get the desired mode of output -- NaN vs NA
# and leading column-number row or not.


######## LOAD LIBRARIES

library(dplR)
library(treeclim)
library(data.table)

######## CLEANING AND PREP

rm(list = setdiff(ls(), lsf.str())) # Remove existing variables except functions

# Hard code some system folders (need to tailor to your laptop)
# For ultimate simplicity you can make all of these folders the same!
# Edit next 4 lines
code_dir <- "/home/dave/GitWork/TRISH-R/" # any user-written functions here
outputDir <-"/home/dave/AAA_sahar/data/crns/" # any output files go here
dataInDir <- "/home/dave/AAA_sahar/data/crns" # rwl input files here
workingDir <- "/home/dave/AAA_sahar/" # folder for this R project 

# set current working directory
setwd(workingDir)

#--- SOURCE

source(paste(code_dir,"crn2tsv.R",sep=""))

#--- HARD CODE

pathin <- dataInDir   # input crns here
suff.pattern <- "\\.crn$"   # search for this pattern for crns
pathout <- outputDir # tab-sep file of chronologies goes here

# Edit next 2 lines
fileoutFull <- "EarthChronsFull.txt"  # EarthChrons is a sample dataset
fileoutPrefix <- 'EarthChrons'

# args for read.crn
read.crn.arg <- list('NULL','getOption("encoding")','TRUE')
names(read.crn.arg) <- c('heading','encoder','long')


# Edit next  line
yrGo <- 1466

# Edit next two lines by commenting out 1 line
LNaN <- 'TRUE'  # missing data as NaN
#LNaN < 'FALSE'  # missing data as NA

# Edit next two lines by commenting out 1 line
#LcolNumbers <- 'TRUE'  # add a leading row of column numbers to tsv output file
LcolNumbers <- 'FALSE' # do not add such a leading row

# Build filename for matrix truncated in yrGo
fileoutCut <- paste(fileoutPrefix,as.character(yrGo),'.txt',sep='')

#--- END HARD CODE


#--- CALL FUNCTION TO READ CRNS, ORGANIZE AS MATRIX, AND WRITE OUTPUT FILE

Din <- list(pathin,suff.pattern,pathout,fileoutFull,fileoutCut,
            read.crn.arg,yrGo,LNaN,LcolNumbers)
names(Din) <- c('pathin','suff.pattern','pathout','fileoutFull','fileoutCut',
  'read.crn.arg','yrGo','LNaN','LcolNumbers')

ChronData <- crn2tsv(Din)

