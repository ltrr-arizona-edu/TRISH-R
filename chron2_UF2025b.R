# chron2_UF2025b.R
# Chronology time series matrix development script
# D Meko
# Last revised 2025-01-12

#
# Reads a file of crn filenames; loops over the filenames, reading in site
# chronologies and building a time series matrix of those, with year as column 1.
# Input specs allow control over the time coverage of the resulting matrix and
# which chronologies it contains. Writes two tab-separated-value (tsv) files of
# time series chronologies. One tsv file has all chronologies covering interval
# from the first year with data in any series to the most recent year with
# data in any series. The other tsv file has the matrix truncated at a specified
# start year. 
#
# Preparation
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
# edit lines dealing with filenames and with truncation start year "yrGo"


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
fileoutFull <- "EarthChronsFull.txt" 
fileoutPrefix <- 'EarthChrons'

# args for read.crn
read.crn.arg <- list('NULL','getOption("encoding")','TRUE')
names(read.crn.arg) <- c('heading','encoder','long')


# Edit next  line
yrGo <- 1466

# Build filename for matrix truncated in yrGo
fileoutCut <- paste(fileoutPrefix,as.character(yrGo),'.txt',sep='')

#--- END HARD CODE


#--- CALL FUNCTION TO READ CRNS, ORGANIZE AS MATRIX, AND WRITE OUTPUT FILE

Din <- list(pathin,suff.pattern,pathout,fileoutFull,fileoutCut,
            read.crn.arg,yrGo)
names(Din) <- c('pathin','suff.pattern','pathout','fileoutFull','fileoutCut',
  'read.crn.arg','yrGo')

Jacky <- crn2tsv(Din)

