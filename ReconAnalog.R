################################################################################
#
# ReconAnalog.R
# Runoff reconstruction by analog and non-analog methods, all two-stage
# D Meko
# Last revised 2024-04-19
#
# This script and functions it calls are the reconstruction framework for TRISH 
# (Tree-Ring Integrated System for Hydrology), a web-based tool developed by a 
# team of researchers from the University of New Hampshire (UNH) and University 
# of Arizona (UA). TRISH interactively builds the reconstruction predictand from
# a global water balance model using climatic inputs from various alternative climate
# products. ReconAnalog can also be run standalone from the R prompt or in Rstudio
# with a predictand time series (e.g., a gaged flow record) supplied by the user.
#
# ReconAnalog run standalone requires a javascript object notation (json) initialization file, 
# recon.init, that specifies inputs, directories, and program settings. TRISH builds this json
# file from user input at web screens. The typical user will not need to change any code in
# ReconAnalog. 
#
# ReconAnalog can optionally do reconstructions by four different statistical methods, each of which
# by a pdf file written by ReconAnalog to a system output folder. 
#
#================ TREE RING INPUT
#
# ReconAnalog requires input time series of tree-ring chronologies and a predictand
# seasonal hyrdrologic or climatological variable. ReconAnalog also requires metadata
# for the chronologies and specifications for the reconstruction. How ReconAnalog gets
# all this data depends on whether ReconAnalog is run within TRISH or standalone from
# your laptop in Rstudio. TRISH culls the tree-ring data from its internally stored time 
# series and metadata, and prepares the predictand seasonal hydroclimatic predictand from
# monthly data via the UNM water balance model. In standalone mode, you supply the 
# data files of tree-ring chronologies and metadata, and the time series of predictand 
# hydroclimate data in tab-separated files with format described below. For example, 
# the tree-ring time series matrix and metadata read by statements:
#
#   U <- read.table(tr_file,sep="\t",header=TRUE) # input chronologies
#   Tmeta <- read.table(trM_file,sep="\t",header=TRUE) # input chronologies
#
# The tree-ring time series data provided by TRISH will already have been screened for
# time coverage and geographic domain before ReconAnalog is called. Outside of TRISH, all 
# chronologies in the provided file of chronologies are available for use. The time series 
# and metadata used outside of TRISH should have the following form.
#
# 1) Time series. Tab-separared matrix of chronologies. The year is column 1 and the time
#   series are in remaining columns. All chronologies should have data in the first year,#
#   but chronologies are allowed to end in different years (e.g., NaN-filled on recent end).
#   Row 1 must contain tab-separated headers, beginning with "Year" for column 1. Remaining
#   headers are site ids with maximum allowed length of 12 characters, Avoid spaces or 
#   hyphens ("-") in site codes.
# 2) Metadata. Tab-separated metadata for the sites in the data matrix. Rows, after header row, 
#   should correspond to the columns of the data matrix, in the same order. 
#   1 sequential number (1 to however many are in the time series matrix)
#   2 site number corresponding to the column of the chronology in the network the
#     the user uploaded to UNH
#   3 site code (1-12 characters, no spaces, and maybe some other rules we should specify)
#   4 longitude east (decimal degrees; negative if west)
#   5 latitude north (decimal degrees; negatve if southern hemispher)
#   6 elevation (m above msl)
#   7 species code (4-letter code, following ITRDB convention)
#   8 data-type (1 letter code): R=total ring width, E=earlywood width, L=latewood width,
#     X=maximum density. 
#   9 First and last year that chronology had valid data.
#   Here is and example line:
#       7	 45	BUT    	  93.367000	  64.283000	  113	LAGM	R 1723   1999
# 
#
#
############# JSON INITIALIZATION FILE
# 
# ReconAnalog was revised in Nov 2022 so that user-settable inputs are no longer modified
# by revising lines of code, but instead by changing settings of an input JSON file (e.g., 
# Recon.init). The JSON file has 28 inputs, described in some detail here. Default 
# settings from Meko's trial run outside of TRISH on a Linux laptop are included. Below, for
# brevity "flow" refers to the predictand hydroclimate time series. Actually, this could be
# some other type of variable (e.g., precipitation). "SSR" and "MSR" refer to "single-site"
# reconstruction (reconstruction of flow separately from each chronology) and "multi-site"
# reconstruction (combination of the SSRs into a single final reconstruction). TRISH users do
# not need to be concerned about preparing a json init file because the json file is generated
# by TRISH by user inputs at the TRISH screens.
#
# "code_dir" :  "/home/dave/Data/RlibraryMeko/",
#   Directory with user-written R functions that ReconAnalog must be able to access
# "pdf_dir" :  /home/dave/Projects/ba2/TRISHvisual/",
#   Directory with four pdfs describing in detail the alternative reconstruction methods.
#   Pdf files are copied from this directory to the output directory on the system running TRISH
# "tr_file" : "treeData.txt",
#   Name of the tab-separated file with time series of tree-ring chronologies. 
# "trM_file" : "treeMeta.txt",
#   Name of file with tree-ring metadata. Likewise, this file is generated by TRISH from user-uploaded
#   data.
# "cl_file" : "hydroData.txt",
#   Name of tab-separated file with time series (year and value) of predictand climatic or 
#   hydroclimatic predictand.
# "outputDir" :  "/home/dave/AAAtrish2/test_out/",
#   The system folder to which the output will be written
# "NameNetwork":  "Kyzyl",
#   The name of the tree-ring network. Standalone users can set this to whatever name
#   desireable; it has now effect on the computations or labeling.
# "PrewhitenOrder" : 0,
#   Whether or not to prewhiten chronologies with AR model before use in regression modeling.
#   If "0", do not prewhiten. If "1", "2", or "3," prewhiten with that order of AR model. Note
#   that this is prewhitening at the site-chronology level, which is different from "residual"
#   chronologies, which apply AR modeling to index time series from individual cores and then 
#   combine those.
# "LallowLags" : true
#   Whether to allow lags in the models. If "false," only lag-0 models are allowed.
#   If "true," lags t-2 to t+2 from the predictand are allowed in pool of potential predictors
# "NsitesUserNetwork" : 274,
#   The number of sites in the network provided to TRISH by the user. Standalone, set this equal
#   to the number of chronologies in the input matrix "trm_file" (see above). The column "N2" in 
#   file "trm_file" is a convenience allowing the user to cross-reference the time series in 
#   "tr_file" to an outside database of user chronologies.
# "YearScreen" : [1700,1997],
#   Start and end year of mandatory common period of coverage by all chronologies
#   to be used in the reconstruction. Others will be ignored by TRISH in culling the
#   user's supplied tree-ring data. All culled chronologies are required to have 
#   complete (no missing) data for the inclusive period bracketed by "YearScreen".
#   TRISH lets user input these two years and does the screening. Standalone, "YearScreen"
#   has no effect, and can be set to the first and last year of the matrix in "tr_file". 
# "NafterYearScreen" : 36,
#   Number of chronologies remaining after screening for "YearScreen." Standalone, 
#   this should be set to the number of time series in the tab-separated time series
#   matrix "tr_file"
# "NafterPolygon" : 36,
#   Number of chronologies remaining after screening for sites being in the specified
#   geographic domain (e.g., polygon drawn interactively to mark tree-ring site domain 
#   in TRISH ). Outside of TRISH, "NafterPolygon" is ignored, and can be set identical
#   to "NafterYearScreen" and "NafterPolygon"
# "HydroVariable" :  "RO",
#   Code for the hydrclimatic variable represented by the predictand. Must be a member of
#   a recognized set of codes. TRISH uses this code for determining how to compute seasonalized
#   data (e.g., precip as a sum, temperature as a mean). The code is also used for building
#   labels for figures and assigning units. From time to time, new hydrological variables
#   are supported. Search for line starting "Dtypes <-" to see the currently accepted codes.
# "ClimDatSet" :  "CRU",
#   The source of the climate data used by UNH water balance model and TRISH to compute the 
#   seasonalized predictand. Must be a member of a recognized set of codes. TRISH lets the
#   user select this from a dropdown menu. Outside of TRISH, "ClimDatSet" has no effect,
#   and is ignored; you can set it to "CRU", for example.
# "HydroSeason" : [9,12],
#   Ending month (1=jan, 12=dec) and number of months in season of predictand. "HydroSeason"
#   tells TRISH how to seasonalize the target predictand from monthly data. Outside of
#   of TRISH, "HydroSeason" defines the season of the input hydrologic time series predictand,
#   and, should match whatever the season is of your input data in "cl_file"
# "yrgoc" : -99999, 
#   Desired start of y data to be used for calibrating SSR models and final (MSR)
#   model. If -99999, use earliest possible year. Here, -99999 is used instead of
#   NA because JSON does not handle NA. "yrgoc" should be in the range of years 
#   covered by the input data in "cl_file". If outside the range, "yrgoc" is forced 
#   to the first year of that data.
# "yrspc": -99999
#   Desired stop of y data to be used for calibrating SSR models and MSR model; If -99999,
#   use latest possible year. Analogous conditions regarding NA, etc., apply as for yrgoc 
#   Setting both "yrgoc" and "yrspc" to -99999 lets the SSR models adapt to the
#   variable overlap of chronologies with flow. Programs gives an error message if
#   your setting for yrspc is inconsistent with the ending year of available predictand, the
#   ending year of the most-recent-ending tree-ring chronology SSR-modeled, and the setting for
#   LallowLags. Of course, yrspc cannot be set later than the ending year of the predictand.
#   But it also cannot be set later than the ending year of the most-recent ending
#   chronology if lags are not allowed, and no later than two years before that if lags
#   are allowed. Note that the chronology constraint on yrspc is base only on the most-recent-
#   ending chronology. The SSR modeling automatically truncates the calibration period to
#   end earlier if some chronology has ends before the most-recent-ending chronology.
# "ktran" : 1,
#   Optional transformation of the predictand, y, before reconstruction.
#   None (1), square root(2) or log10 (3). The resulting reconstruction will  be in units 
#   of transformed flow. TRISH issues error messages if the selected transformation
#   is inconsistent with the data (e.g., log10 transformaton when there are flows 
#   of zero)
# "methMSR" :  2,
#   Method for MSR reconstruction. This and PCApredictors effectively specifies the 
#   method to be used. (1) Simple Linear regression (SLR), (2) MLR (multiple linear regression)
#   on SSRs or their PCs, (3)) Analog nearest-neighbor PCA 
# "PCApredictors" : true,
#   Whether PCs of the SSRs (true) or the SSRs themselves (false) comprise the pool of
#   potential predictors for the MSR.
#  kHowPCA <-2  # if PCA, is it done on the correlation (1) or covariance (2) matrix
#   of the SSRs. Covariance matrix makes more sense if differences in the variances 
#   of individual time series on which the PCA is run are important. Because ReconAnalog
#   runs PCA on SSRs, and because the variance of an SSR reflects its calibration accuracy,
#   the best selection here is kHowPCA=2. However, the option is available to do the PCA
#   on the correlation matrix.
# "PCoption" : 2,
#   If reconstruction method is MLR of y on PCs of the SSRs (methMSR=2,
#   PCApredictors=TRUE), PCoption specifies how many and which PCs should be in the 
#   pool of potential predictors. If a different reconstruction method, "PCoption"
#   has no effect. Options are 1) directly specify the number of PCs, and 2) use the
#   m<fN PCs whose scores are highest-correlated with y. Here, N is the number of years 
#   for calibration of the MSR model, and f is a factor less than or equal to 0.2. 
#   By default, f=0.10, meaning that the pool of potential predictors must be less than
#   1/10 N, where N is the length of the calibration period.
# "nPCsKeep" :  1,
#   If PCoption=1, this is the number of PCs to keep in the pool of potential predictors
#   of the MSR model (inclusion by priority of size of eigenvalue). If PCoption!=1, 
#   nPCsKeep is ignored 
# "f" :   0.10,
#   A decimal fraction limiting the sized of the pool of potential predictors of the
#   MSR model to less than this specified decimal fraction of number of years in the 
#   calibration period of the MSR model. For example if 100 years in the calibration period,
#   f=0.10 dictates fewer than 10 potential predictors in that pool.
# "alphaR"  :  0.05,
#   Critical alpha-level for correlation (r) of y with SSRs or their PCs when using the
#   analog (methMSR=3 method for MSR reconstruction. Otherwise, alphaR is ignored.
#   The only acceptable values for alphaR are 0.10, 0.05, and .01, corresponding to
#   90%, 95% and 99% significance by a two-tailed test of # the null hypothesis of 
#   zero population correlation. In analog MSR reconstruction, only those PCs 
#   significantly correlated with flow at level alphaR are used to identify analog
#   years.
# "Lcausal" : true
#   Whether SSR models should be rejected if the only predictors of flow in year t are
#   tree-ring variables in earlier years (e.g., t-1 or t-2). It is sensible to set
#   Lcausal=true if the chronologies are standard chronologies. With residual chronologies,
#   this is still probably the more reasonable choice. 
# "RequireStable" : true
#   Whether the non-rejection of a SSR model should depend on the temporal stability of the 
#   model. "TRUE" instructs to reject model if the RE for either half of the split-sample
#   validation is negative. "FALSE" instructs the computed split-sample RE to be ignored in
#   deciding whether to reject. "TRUE" is the more stringent requirement, and it seems reasonable that
#   a model fit to the second half, say, of the overlap of flow and tree rings should do a
#   better job of reconstructing the first half of flow than a "null" reconstruction equalt
#   to the second-half mean of observed flow in every year. But, plausibly, one half of the
#   flow record itself might be flawed (e.g., less reliable gage, or very little year to year
#   climate fluctuation. The option is therefore available for "RequireStable"=false.
#
#--- REVISIONS
#
# revised 2023-03-24: to check specified json inputs yrgoc and yrspc against start and end year of 
#   available seasonalized hydro series. Error message if inconsistent. Similar check against end year
#   of most recent tree-ring series, which must (for lagged predictors) extend at least two years beyond
#   yrspc. 
# revised 2023-04-09: to allow screening of SSRs to optionally be relaxed so that temporal 
#   stability (positive split RE on both halves) not required
# revised 2023-04-17: incorporate multiple edits from AlexP, including optional reading of json filename from
#   command line; otherwise, if no command line args, assumes "init01.json"
# revised 2023-04-19: write progress percentage (see ) for use by TRISH in progress bar
# revised 2023-05-08: 1) add "seasonal" to the labeling of hydrovariables in plots,
#   2) add "Soil Moisture" to the allowable date types, 3) add annotation of the season end month and
#     number of months to the SSR Figure that shows decrease in signal strengh as chronologies drop out
# revised 2023-05-16: to ensure that the provided reconstruction does not extend beyond (on recent end) 
#   the more recent of 1) specified end year of calibration,  and 2) last year in which all screened
#   SSR have data. This is to avoid flakey reconstruction results that can result in last few years
#   when screened SSR matrx us extended statistically to cover throught the last year of the most
#   recent SSR
# revised 2023-05-25: for better appearance of the multi-line text annotation on Figure 01 - SSR. 
#   Revision ensures that lines of that annotation are not crunched together, which could formerly
#   happen depending on user's screen resolution.
# revised 2023-05-26: to add a column "Sign" to the two output tables summarizing SSRs. Already 
#   had a column "Model" that indicated (e.g., "00120") which lags (L to R, t-2 to t+2) are in
#   model and in which order they entered (here, lag 0 at step 1, lag t+1 at step 2).  Now have 
#   column "Sign" telling the sign of those regression coefficient. For this example, 
#   "00NP0" indicate negative lag-0 coefficient and positive lag-1 coefficient. This 
#   revision involved modifying lists and indexing. Function reconsw4 also had to be
#   revised, and a new function, LagModel2SignR, written.
# revised 2023-05-26: "check.names=FALSE" to read.table to overcome problem with conversion 
#   of hyphens and underscores in column headings.
# revised 2023-11-28.  (1) To use as the last year of the generated final reconstruction (MSR)
#   the last year of the quantile-extended matrix of screened SSRs. See code
#   "yrEnd <- yrY3[length(yrY3)]". (2) To add annotation of the specified end year of calibration
#   period to SSR plot summarizing drop in accuracy as chronologies drop out on recent end.
#   (3) Revision of comments for clarity and correction of typos. (4) Added check that seasonalized
#   climate time series has no internal NA and a year vector that increments by 1 year
# revised 2024-03-04: (1) added code (commented out) that can be uncommented for testing reconsw4
#   (search "#debug on"), (2) revision of comments for typos and clarity, (2) Figure01-SSR1.png 
#   cosmetic change to fix annotation sometimes overlaying bars of bar chart
# revised 2024-04-18: (1) corrected code "(!is.na(yrgoc) & !is.na(yrspc))" that would
#   have made program bombing when input yrspc is set to NA and yrgoc is set to a specific year.
#   (2) extended SSR Figure 3 (search "strEndF") to include annotation of the latest feasible
#   SSR calibration period ending year given the time coverage of tree-ring data, time coverage
#   of predictand, and logical setting for whether lags to be allowed in models. (3) cosmetic 
#   correction of typos a syntax in comments.



######## LOAD LIBRARIES

library(car) # for durbinWatsonTest()
library(rjson) # for reading json input files
library(pracma) # for isempty)
library(ggplot2)

################################################################################

# Remove all variables except functions
rm(list = setdiff(ls(), lsf.str()))

################################################################################
#
# JSON initializatioon input. Input data as list from json; in for loop, store the list elements
inputArgs = commandArgs(trailingOnly = TRUE)	### AlexP addition to use input from command line argument list 
if (length(inputArgs) > 0) {			### AlexP addition to use input from command line argument list
  json_file	= as.character(inputArgs[1])	### AlexP addition to use input from command line argument list
  myData <- fromJSON(file=json_file)		### AlexP addition to use input from command line argument list
}	else {
  myData <- fromJSON(file="Recon.init")  ### DaveM so that assumes this json name if no args in command line
}				

naJSON <- -99999 # regard -99999 as NA
X<- myData
for (j in 1:length(X)){
  a <- names(X)[j]
  b <- paste(a,'<-X$',a,sep='')
  s <- paste(a,'<-',b)
  eval(parse(text=s))
}
rm(a,b,s)
if (yrgoc == naJSON){
  yrgoc <- NA
}
if (yrspc == naJSON){
  yrspc <- NA
}

### PROGRESS BAR -- FUNCTION AND FIRST CALL
#
# TRISH online tool writes a progress bar, which requires a one-line output message at various
# stages of the analysis reporting the percentage of run completed. For now, these messages are
# all writting in this main calling script ReconAnalog. We estimate that 15% of the time is taken up
# by TRISH before calling ReconAnalog -- so, the starting 15% complete.
#{ "message": "Starting the analysis...", "percent": 15 }

# For UNH server, must use a different target path/filename for the progress update file than in Rstudio
# on laptop
if (length(inputArgs) > 0){
  pfProg <- paste(outputDir,'../status.js',sep='') # this if on UNH server
} else {
  pfProg <- paste(outputDir,'ProgressTRISH.txt',sep='') # this if not on UNH server
}
pctDone <- 15; pctInc <- 0 # starting pctg done and increment of pctg to be apply next
mssgProg <- "Starting the analysis..."
ProgTrack <- function(pfProg,mssgProg,pctDone,pctInc) {
  pctDone <- pctDone + pctInc
  mssgThis <- paste('{ \"message\": \"',mssgProg,'\", \"percent\": ',as.character(pctDone),' }',sep='')
  fprintf('%s',mssgThis,file=pfProg,append="FALSE")
  return (pctDone)
}
pctDone <- ProgTrack(pfProg,mssgProg,pctDone,pctInc)


######## SOURCE THE FUNCTIONS ReconAnalog DEPENDS ON

source(paste(code_dir,"TranFlow.R",sep="")) # optional transformation of flows
source(paste(code_dir,"trimnan.R",sep="")) # get indices to allow trimming of leading and trailing NA
source(paste(code_dir,"emssgUNH.R",sep="")) # write error file to system, specified output folder
source(paste(code_dir,"reconsw4.R",sep=""))  # single-site reconstruction (SSR) by distributed-lag stepwise regres
source(paste(code_dir,"TrimTsm1.R",sep=""))  # trim tree-ring time series matrix in preparation for SSR
source(paste(code_dir,"trimRowNA.R",sep=""))  # row indices of a matrix after trimming off leading and trailing all-NA rows
source(paste(code_dir,"tsmExtend.R",sep="")) # extend time series matrix on recent end
source(paste(code_dir,"RecSLR1.R",sep="")) # Reconstruction by simple linear regression of y on mean of SSRs
source(paste(code_dir,"RecMLR1.R",sep="")) # Reconstruction by multiple linear regression or analog method
source(paste(code_dir,"SignalDrop1.R",sep="")) # drop in maximum accuracy as tree-ring network thins in recent years
source(paste(code_dir,"PrewhitenChrons.R",sep="")) # convert tsm of chronologies to prewhitened chronologies

MinCalibLength <- 30 # hard-coded minimum allowable length of calibration period for SSR models

########## READ FILES OF PREDICTAND, TREE-RING TIME SERIES, TREE-RING METADATA

V <- read.table(cl_file,sep="\t",header=TRUE) # input flow; UNH will not read this in because
# user picks the hydroData variable interactively in TRISH
U <- read.table(tr_file,sep="\t",header=TRUE,check.names=FALSE) # input chronologies
Tmeta <- read.table(trM_file,sep="\t",header=TRUE) # input chronologies
# Status. U is tsm of chronologies and Tmeta is table of corresponding metadata. So far
# the chronologies in the tsm and metadata table are those satisfying (1) in polygon, and (2) complete

# Update progress bar
mssgProg <- "Data read in..."


########## OPTIONAL PREWHITENING OF CHRONOLOGIES
#
# Option to prewhiten (remove autocorrelation from) chronologies using autoregressive model
# of order 1, 2, or 3 (AR(1),AR(2), AR(3)). Order 0 instructs to not prewhiten. 
# Prewhitening essentially allows user to check whether "residual" chronologies might have
# stronger signal than standard chronologies for y. If so, user may want to develop
# residual chronologies as separate network and submit those to TRISH for testing. 
# Prewhitened chronologies are not exactly the same as residual chronologies because with
# prehitening the autoregressive modeling is done on the site chronology, but with 
# residual chronologies (as defined in dendro literature) the modeling is done on the 
# standard core indices before averaging those into a site chronology. But experience
# has shown the two versions (prewhitened and residual) are very similar.

if (PrewhitenOrder==0){
  # No action needed
  PWtext <- 'Not prewhitened'
  } else {
    if (any(PrewhitenOrder==c(1,2,3))){
      # Call function to covert U to prewhitened chronologies
      # If order p, will end up converting first p values of chronology to NA (lose p leading years)
      PWtext <- paste('Prewhitened with AR(',as.character(PrewhitenOrder),') model',sep='')
      ResTemp <- PrewhitenChrons(U,PrewhitenOrder,outputDir)
      U <- ResTemp$Xwhitened
      mssgProg <- "Data read in and chronologies whitened..."
    } else {
      # Invalid order (must be 1,2,or 4); stop with error message
      # Write error file
      emssg <- 'Allowable AR(p) prewhitening models are p=0,1,2 or 3'
      ResTemp<-emssgUNH(emssg,outputDir)
      stop(emssg)
    }
  }

# Update progress bar
pctInc <- 5
pctDone <- ProgTrack(pfProg,mssgProg,pctDone,pctInc)
# pctDone is 20% to her


# Get first and last year of the tree-ring TSM. It is assumed that this matrix has 
# no all-NA rows
yrTree <- U[,1]
yrgoTree <- yrTree[1]; yrspTree <- yrTree[length(yrTree)]  


### Check that ending year of most-recent-ending tree-ring chronology consistent with input specified
# end year of calibration period. If lags allowed, tree-ring series must extend at least 2 years
# beyond end of specified calibration period. If most-recent-ending chronology has too early an end year, 
# all of the chronologies will similarly have too early an end year. Best in that case to abort here. Later,
# checks are done in function reconsw4 on individual tree ring chronologies.

nExtra <- 0 # initialize how many years tree ring series must extend beyond yrspc
# If lags enabled, it is assumed that maximum negative and positive lag is 2 years.
# Will therefore need 2 years of tree-ring chronology after specified yrspc ending 
# year of the calibration period to satisfy the +2 lag on predictors
if (LallowLags){
  nExtra <- 2
}

if (!is.na(yrspc)){
  L1 = yrspc > (yrspTree-nExtra)
  if (L1){
  emssg <- paste('For tree matrix ending in ', as.character(yrspTree), 
  ', the specified last year of calibration cannot be later than ', as.character(yrspTree-nExtra),
  '\n  (two additional years needed for lagged model)',sep='')
  ResTemp<-emssgUNH(emssg,outputDir)
  stop(emssg)
  }
  rm(L1)
}


# Check that input yrgoc and yrspc -- when both are specified rather than NA -- 
# are consistent with hard-coded minimum allowable length of calibration period of SSR models.
# If one or both of yrgoc and yrspc are input as NA (-99999 in Recon.init), the
# check is done series by series within function reconsw4.
L1 <- (!is.na(yrgoc) & !is.na(yrspc))
if (L1){
  ntemp <- yrspc-yrgoc+1
  if (ntemp < MinCalibLength) {
    emssg <- paste('Specified calibration years ',as.character(yrgoc), ' and ', as.character(yrspc),
       ' mark a calibration period fewer than the minimum allowable ',as.character(MinCalibLength),
       ' years.', sep='')
    ResTemp<-emssgUNH(emssg,outputDir)
    stop(emssg)
  }
}


##########  CHECK NUMBER OF SERIES IN TREE-RING MATRIX

nms1<-colnames(U[-1]) # chronology ids; if using TRISH, for all chronologies from polygon screening
nchron<-length(nms1) # number of chrons from polygon screening
if (NafterYearScreen != nchron){
  # Write error file
  emssg <- 'Number of data columns of input tree-ring matrix inconsistent with nAfterPolygon'
  ResTemp<-emssgUNH(emssg,outputDir)
  stop(emssg)
}

######### BUILD SOME STRINGS TO BE USED LATER IN LABELS

ClimDataSet <- c('CRU','Delaware','Reanalyis') # input climate data used by UNH WBM in TRISH
# This variable is included here for information only. When ReconAnalog is used
# within TRISH, the user has a dropdown window that allows selection of the input
# data. But ReconAnalog does not care about the input data set and works only with
# whatever processed WBM output or seasonalized climate series is provided. Also,
# seasonalizing monthly data is done ahead of calling ReconAnalog (either in TRISH or
# outside of TRISH), there is no need to specify that a particular variable is summed (e.g., P)
# vs averaged (T) to generate seasonal data from monthly data.
#
# Note that TRISH so far deals with only the first four Dtypes below. Other types are
# included for use of ReconAnalog outside of TRISH. For example, F1 is volume flow of a
# river in thousand acre feet (kaf). These extra types are needed so that graphics 
# have properly labeled axes. Axes labels use DtypesLabel rather than Dtypes 
# (e.g., "Flow (kaf)" rather than "Flow1 (kaf)".

# Discharge, Runoff, Soil Moisture, Temperature, Precip
Dtypes <- c('Q','RO','SM','T','P','Flow1')
DtypesLabel <- c('Seasonal Q','Seasonal RO','Seasonal SM','Seasonal T','Seasonal P','Seasonal Flow')
LabsUnits <- c('(cms)','(mm)','(mm)','(Deg C)','(mm)','(kaf)')
HydNames  <- c('Discharge','Runoff','Soil Moisture','Temperature','Precipitation','Flow')
ithis = which(Dtypes==HydroVariable)
if (isempty(ithis)){
  emssg<-paste(HydroVariable, ' is not one of: ', paste(Dtypes,collapse=' '),sep='')
  ResTemp<-emssgUNH(emssg,outputDir)
  stop(emssg)
}
Dtype<-DtypesLabel[ithis]; LabUnits<-LabsUnits[ithis]; HydName<-HydNames[ithis]

# Text to be used in reminding user of the selected season; annotated on one of the SSR figures
MonthsOfYear <- c('Jan','Feb','Mar','Apr','May','June','July','Aug','Sept','Oct','Nov','Dec')
LabelSeason1 <- paste('Ending Month of Season = ',MonthsOfYear[HydroSeason[1]],sep='')
LabelSeason2 <- paste('Season length = ',as.character(HydroSeason[2]),' Months',sep='')



########## RENAME A VARIABLE; SET UP MENU OF METHODS

n1<- NsitesUserNetwork # number of sites in the full user-supplied network; 
# system must provide this; renamed here for convenience
msrMethods <- c('SLR1','MLR1-PCA','MLR1-noPCA','Analog')


########## BUILD FILENAME FOR PDF DESCRIBING MSR METHOD TO BE USED
#
# This pdf will we written to system output folder. User and refer to the pdf for
# details on the method used, including definitions of downloaded output files

if (methMSR==1){
  jMethod <- 1
} else if (methMSR==3) {
  jMethod <-4
} else {
  if (isTRUE(PCApredictors)) {
    jMethod <-2
  } else {
    jMethod <- 3
  }
}
PdfFile <- paste('TrishOutputDescribe',msrMethods[jMethod],'.pdf',sep='')
PdfDescribe <- paste('See',PdfFile,'for more',sep=' ')
rm(jMethod)


########## CLEAR OUTPUT DIRECTORY AND COPY PDF DESCRIBING MSR METHOD THERE

PathClean <-  paste(outputDir,'*',sep='')
unlink(PathClean,recursive=F,force=F)
file.copy(from=paste(pdf_dir,PdfFile,sep=''),to=paste(outputDir))
rm(PathClean)

############# HARD CODED SETTINGS
# 
# These setting not to be changed by casual user. Can be used by developer to explore 
# possibilities for extending or modifying ReconAnalog. TRISH users cannot change
# these setting from TRISH menus. Standalone users could change setttings, but run the
# risk of unintended consequences. 

nNeg<-2 #$ maximum negative lag allowed on chronologies for SSR models
nPos<-2 # maximum positive lag allowed on chronologies for SSR models
# ReconAnalog() was written to specifically apply lagged regression with maximum of
# 2 negative and positive lags on the chronology. Accordingly, leaving 9 out in 
# cross-validation guarantees that cross-validation estimated do not depend on
# any of the tree-ring data that are used also in calibrating the cross-validation
# model. 
# yrgo1 and yrsp1 best both set to NA. This lets the time coverage of the tree-ring data itself
# determine the time comverage of the reconstruction. In some early trialsl I played with varying
# yrgo1 and yrsp1, but decided makes more sense to go with NA
yrgo1<-yrgoTree+nNeg # desired start year of reconstruction; actual reconstruction coverage will depend on 
#   coverage of tree-ring chronologies in final model
yrsp1<-yrspTree-nPos  # desired end year of reconstruction (including through calibration period)
# yrgo1 and yrsp1 are the desired start and end years of reconstructed flow.
# The tree-ring matrix supplied by TRISH will be trimmed in ReconAnalog() to include only
# those chronologies with data in year yrgo1-2 (allows for two negative lags in 
# single-site regression).
# This will eliminate from consideration any chronologies that start at a later year. 
# The tree-ring matrix will be trimmed to end in the earlier of (yrsp1+2) and (the last year
# of data for any site passing the screening for yrgo1). If yrgo1=NA, the matrix is trimmed to
# start with the first year for which all sites in the basin domain have data. If yrsp1=NA,
# the matrix is trimmed to end with the last year of data at any one of the sites.
N1 <- 50 # in forward extension of SSR matrix, common period of all SSRs must be
# at least N1 year (e.g., 50
N2 <- 100 # in forward extension of SSR matrix, series needing an extension in year
# i must overlap by N>=N2 years with some other series that does have a value in 
# year i.
N3 <- 30 # minimum acceptable number of years for calibration on MSR. Allows 15/15, which is
# ridiculously low, for split-sample validation. 
incR2a<-0.01 # Critical increment in adjusted R-squared of MSR model. Stepwise model is assumed to reach
# "approximate maximum" when next step would yield increase in adjusted R-squared less than inR2a.
# Stepwise models in SSRs and MSR are not allowed to proceed beyond the approximate maximum of
# adjusted R-squared. Further, depending on "kstop" (see nex), the model may stop an an even
# earlier step to satisfy constraints on cross-validation. 
kstop <-2 # Stepwise forced to stop at either the maximum adjusted R-squared (kstop=1) or at some earlier
# step if cross-validation RE reaches a maximum at an earlier step (kstop=2)
ScreenAnalogPCs <- TRUE # Screen the PCs used in analog reconstruction (methMSr=3)
# by correlation with predictand. If TRUE, only those PCs whose correlations with y
# are significant at alpha-level alphaR (see earlier) are used to identify analogs

######### MAKE FLOWS MATRIX; TRIM OFF ANY LEADING OR TRAILING NA; STOP IF INTERNAL NA

V<-as.matrix(V)
v<-as.matrix(V[,2])
v <- as.numeric(v) # in case any NA in v
v <- as.matrix(v)
yrv <- V[,1,drop=FALSE]
yrv <- as.matrix(as.numeric(yrv))

i1<-trimnan(v) # row indices of v with leading and trailing NAs stripped
vTemp <- v[i1]
yrvTemp <- yrv[i1]

# Check that no internal NA in the nan-trimmed time series of seasonalized climate 
L<- any(is.na(vTemp))
if (L) {
  emssg<-'Internal NA (missing value) in the vector of seasonalized climate time series'
  ResTemp<-emssgUNH(emssg,outputDir)
  stop(emssg)
}

# Check that the year increments by 1 in the seasonalized climate series 
#  (i.e., in case there is a row missing, but no internal NA)
d = diff(yrvTemp)
L<-!all(d==1)
if (L) {
  emssg<-'Year for seasonal climate input does not increment by 1'
  ResTemp<-emssgUNH(emssg,outputDir)
  rm(ResTemp)
  stop(emssg)
}

# Row index i1 OK; re-grab vectors of the climate series and its year
V<-V[i1,,drop=FALSE] # trim the matrix of leading and trailing NAs
v <- v[i1,,drop=FALSE]; yrv<-yrv[i1,,drop=FALSE]
rm(L,d,i1,vTemp,yrvTemp)

# Check that specified desired start and end year of calibration period consistent with time coverage of
# predictand, v. 
if (!is.na(yrgoc)){
  L <- yrgoc < yrv[1]
  if (L){
    emssg<-paste('Specified yrgoc=',as.character(yrgoc),' is earlier than start (',as.character(yrv[1]),
    ') of seasonalized hydro series v',sep='')
    ResTemp<-emssgUNH(emssg,outputDir)
    rm(ResTemp)
    stop(emssg)
  }
  rm(L)
}
if (!is.na(yrspc)){
  L <- yrspc > yrv[length(yrv)]
  if (L){
    emssg<-paste('Specified yrspc=',as.character(yrspc),' is later than end (',as.character(yrv[length(yrv)]),
                 ') of seasonalized hydro series v',sep='')
    ResTemp<-emssgUNH(emssg,outputDir)
    rm(ResTemp)
    stop(emssg)
  }
  rm(L)
}


############### TRANSFORM FLOWS (OPTIONAL)
# Allowed are square root or log10. Call a function to do the transform. If transform
# not reasonable physically, function called returns a flag that prompts this script
# to abort and also prints a message to outputDir

kBogus<-FALSE
Transformed<-FALSE
if (ktran==1){
  # If not call TranFlow, set ResTf to empty list
  ResTf<-vector(mode = "list", length = 0) 
} else if (ktran==2 || ktran==3) {
  ResTf <- TranFlow(v,ktran)
}  else {
  kbogus<-TRUE
}

#--- UNH HANDLING OF TRANSFORM MESSAGE
emssg<-'None' # initialize error message
if (kBogus){
  emssg<-'Invalid specified ktran: option must be 1, 2 or 3'
}
if (length(ResTf)==0){
  sTran<-''
} else {
  if (ResTf$flag==0){
    sTran<-ResTf$sTran
    V[,2]<-ResTf$x
    
  } else if (ResTf$flag==1) {
    sTran<-''
    emssg<-'Sqrt transform invalid;<br />series has negative values'
  } else if (ResTf$flag==2) {
    emssg<- 'Log10 transform invalid;<br />series not all-positive'
  }
}

# Conditional error message to OutputDir
if (emssg=='None'){
} else {
  # Write error file
  ResTemp<-emssgUNH(emssg,outputDir)
  stop(emssg)
}

# Depending on transform, modify y units, and store with other information applicable to 
# all four MSR methods. First might have to modify label for units of y
if (ktran==1){
  # no transform, units OK
} else if (ktran==2) {
  LabUnits<-paste('[sqrt',LabUnits,']')
} else if (ktran==3){
  LabUnits<-paste('[log10',LabUnits,']')
}
txtSeas <- paste0(as.character(HydroSeason[2]),'-month season ending in month ',as.character(HydroSeason[1]))
RecListx<-c(Dtype,LabUnits,HydName,txtSeas) # general list for calls to MSR recon


###### CHECK FOR INTERNAL NA IN TREE-RING SERIES, ABORT WITH MESSAGE IF FOUND

emssg<-'None'
for (n in 1:nchron){
  j<-n+1
  u <- as.matrix(U[,j])
  ResTemp<-trimnan(u) # index to rows of nan-trimmed u
  u1<-u[ResTemp]
  L<-any(is.na(u1))
  if (L){
    emssg<-paste('Internal NA in tree-rng series ',as.character(n),': ',nms1[n])
    ResTemp<-emssgUNH(emssg,outputDir)
    stop(emssg)
  }
}

################################################################################
#
# TRIM TREE-RING MATRIX TO COVER JUST THE PERIOD NEEDED FOR SSR'S; COL-TRIM AS
# NEEDED SO THAT TSM INCLUDES ONLY THE CHRONS THAT COULD PROVIDE THE RECONSTRUCTION
# TIME COVERAGE; ALSO ROW-TRIM THE TREE-RING METADATA TABLE ACCORDINGLY
#
# See earlier comments for inputs yrgo1 and yrsp1
#
# The reconstruction method allows lags t-2 to t+2 on tree ring. Therefore, if you
# specify you want the reconstruction to start in year t0, all series must have 
# data in year t0-2. This next section may therefore reduce the number of colunns in 
# the tree ring matrix by removing those without data in year to-2. The matrix will
# be row-truncated to begin in year yrgo1-2 and to stop in the earlier of
# [year yrsp1+2 and the last year with data for any chronology]

mlead<-2; # this many leading tree-ring values are needed to produce first reconstructed value
ResTrim1 <- TrimTsm1(U,yrgo1,yrsp1,nNeg,nPos)
X<-ResTrim1$X # tsm of tree-ring chronologies to be converted by SSR
ix2 <- ResTrim1$ix # index of series in X to the columns of tsm of map-screened chronologies
# provided by system
nms2 <- nms1[ix2] # column names of the (possibly) column-reduced tree-ring matrix
# These ids are identically equal to the colnames (after Year) of X, as returned
# by TrimTsm1()

# Metadata trimming,and check for exact match of column headers of chronology tsm with
# Ids in the metadata
Tmeta <- Tmeta[ix2,]  # row-trim
IdMeta=Tmeta$Id
IdMeta <- gsub(" ", "", IdMeta, fixed = TRUE)
L<-identical(IdMeta,nms2)
#--- Bomb message if not exact match
emsg1 <- c('Ids in header row of time series matrix of chronologies do not exactly match Ids in metadata')
if (!L){
  eBomb<-emssgUNH(emsg1,outputDir)
  stop(eBomb)
}

#---  Compute last allowable end year for SSR calibration. Could be limited by last
# yeqr of available tree rings or last year of available predictand. And must consider if
# lags allowed on tree rings. 
yrT <- X[,1] # year vector for tree-ring matrix of those chronologies to be SSR modeled
yrTend <- yrT[length(yrT)]# last year of the SSR-capable tree-ring matrix
if (LallowLags){
  yrspcLimitT  <- yrTend-2 # limit on end year of SSR calibration imposed by tree-ring coverage
} else {
  yrspcLimitT  <- yrTend
}
yrspcLimitv <- yrv[length(yrv)] # limit on end year of SSR calibration imposed by predictand coverage
end
yrspcLimit <- min(c(yrspcLimitT,yrspcLimitv)) # SSR calib period cannot end later than this year

################################################################################
#
#  CONVERT TREE-RING SERIES TO SINGLE-SITE RECONSTRUCTIONS (SSR's) OF FLOW
#
# All series in X will be converted. The regression statistics and the SSRs will be
# stored. Only a subset of those SSRs, depending on the calibration and validation
# statistics, will ultimately be used in the PCA analysis and final reconstruction.
# The SSR statistics will include a "reject" flag, which will be used to screen out
# bad tree-ring sites. A site will be rejected if any of the following are true:
# 1) p-value of overall F of regression model >=0.05
# 2) REcv (Reduction of error statistic from leave-9-out cross-validation) <=0
# 3) RE from either half of split-sample cross-validation <=0
# 4) Model uses just past tree-ring index to reconstruct current y  (optional, see Lcausal input)
#
# SSR modeling is done by reconsw4(), and comments there have more details on the method.
# Essentially, y(t) (flow in year t) is regressed stepwise forward on tree-ring
# {x(t-2), x(t-1), x(t), x(t+1),x(t+2)}. Entry of variables is stopped if no increase
# in validation skill, or very small in increase in adjusted R-squared (see input ktop)

# SSR time series will be in matrix Y1; these will be the full set of estimated SSRs
nchron2 <- dim(X)[2]-1 # number of chronologies after screening for time coverage
nY1<-nchron2;
mY1 <- dim(X)[1] # Y1 will store the SSRs, and will initially be same size as X
Y1 <- matrix(data=NA,nrow=mY1,ncol=nchron2) # to hold the SSRs
yrY1 <- matrix(X[,1])
yrgoY1<-yrY1[1,1]
yrspY1<-yrY1[mY1,1]


#--- Table of statistics of all models passing the spatial-temporal screening.
# The data frame SSRdf1 is the summary table, described below. A second table,
# SSRdf2, has a similar structure, but includes only those site passing statistical
# screening for a strong, validated, temporally stable signal for flow.

SSRdf1 <- data.frame(matrix(nrow = nchron2, ncol = 16))
names(SSRdf1) <- c("N1", "N2","Site","StartC","EndC","Model","Sign","R2adj","Fp",
"REcv","REa","REb","Refit","StartR","EndR","Reject")
# N1, N2, are sequential number in this table, and in the user's full network.
# Site is a site id; StartC & EndC are first and last years of calibration period;
#   Model is a coded string indicated which of lags t-2 to t+2 are in final model
#   and the order of entry; Sign is a similar string indicating sign of regression
#   coefficients
# R2adj is regression adjusted R-square; Fp is the p-value of overall F of the
#   regression (pF<0.05 means significant model); 
# REcv is th reduction-of-error statistic from leave-9-out cross-validation; REa 
#   and REb are reduction-of-error for split-sample calibration/validation (e.g., 
#   REa is for validation on first half when model is fit to second half; 
# Refit is a logical (0 0r 1) indicating whether the final selectio of lags
#   allowed refitting of the model to a slightly longer calibration period than 
#   possible if all lags t-2 to t+2 are in model; 
# StartR and EndR are the start and end years of the reconstruction (SSR for this site)
# Reject is a logical (0 or 1) indicating if site was rejected for further use in 
#   in reconstruction. Note that only sites not rejected are used in the later
#   PCA step. Numbering n3 therefore goes with sites having a 0 in this column.


#--- Generate SSRs
#
# Function reconsw4() is called in a loop over all chronologies in the drawn polygon
# that covers the minimum acceptable reconstruction interval input at the TRISH window.
# Each chronology is converted to an individual estimate of flow by lagged stepwise
# forward regression. Five lagged values, lagged t-2 to t+2 from the year of flow, are
# considered as potential predictors. The forward stepwise process is stopped when an
# an additional step would result in less than c1 increase in adjusted R squared.
# Leave-9-out cross-validation and split-sample calibration/validation are then applied
# to test the skill of prediction and to possibly simplify the model. Input option
# "LallowLags" lets you forgo the lags, and in that case also uses leave-1-out 
# cross-validation,
# 
# Status: V and X are data frames with prepared flows and tree-ring chronologies

timeGo<-proc.time() # start timer to check cpu time and clock time for running 
yrsCalWindow <- c(yrgoc,yrspc) # calibration of model will consider flows only within this
# window

SSRprelim<-vector(mode="list",nchron2) # to hold lists from preliminary SSR modeling for each chronology
SSRlags1 <- vector(); SSRlags2<-vector() # initialize empty vectors to hold concatenated lags
# in models (all SSRs for1 ; those passing final screening for 2)
# Fill list SSRdf1 members with network site numbers and site IDs as supplied in list treeMeta
SSRdf1[,2]<-Tmeta$N2
SSRdf1[,3]<-IdMeta

# Progress bar strategy. Know have 20 % done and will allocate additional 65% to SSR modeling, which 
# can be time consuming, especially when there are many chronologies Will update
# progress every 20 chronologies. If 20 or few chronologies only one updated needed
nslabs <- ceil(nchron2/20) # will divide the available 45% into nslabs increments such than 
# percentage work required for SSR modeling will not exceed 45%
pctInc <- floor(65/nslabs)
rm(nslabs)

iprogress <- 0 # for updating progress bar every 20 SSRs


for (n in 1:nchron2){ # loop over tree-ring chronologies
  x <- as.matrix(X[,c(1,(n+1))]) # 2-col data frame with year and chronology
  
  # #debug on  test reconsw4
  # # Uncomment block after debugging
  # if (n==26){
  #   save(x,V,nNeg,nPos,yrsCalWindow,incR2a,Lcausal,LallowLags,MinCalibLength,RequireStable,
  #        file = "a1.RData")
  # }
  # #debug off
  
  ResSSR <- reconsw4(x,V,nNeg,nPos,yrsCalWindow,incR2a,Lcausal,LallowLags,
                     MinCalibLength,RequireStable)
  SSRprelim[[n]]<-ResSSR
  SSRdf1[n,1]<-n
  SSRdf1[n,4]<-ResSSR$yearsCal[1]
  SSRdf1[n,5]<-ResSSR$yearsCal[2]
  SSRdf1[n,6]<-ResSSR$ModelCoded
  SSRdf1[n,7]<-ResSSR$ModelSign
  SSRdf1[n,8]<-ResSSR$RsquaredAdj
  SSRdf1[n,9]<-ResSSR$Fp
  SSRdf1[n,10]<-ResSSR$REcv
  SSRdf1[n,11]<-ResSSR$REa
  SSRdf1[n,12]<-ResSSR$REb
  SSRdf1[n,13]<-ResSSR$Lrefit
  SSRdf1[n,14]<-ResSSR$yearsRec[1]
  SSRdf1[n,15]<-ResSSR$yearsRec[2]
  SSRdf1[n,16]<-ResSSR$Lreject
  
  # Build vector of lags in models (1=t-2... 5=t+2)
  # Make two versions: one for all SSRs, the other for those passing screeing.
  # Result (SSRlags1, SSRlags2) are vectors with the concatenated lags in models.
  # Values 1,2,3,4,5 correspond to lags t-2,t-2,t,t+1,t+2
  
  SSRlags1 <- c(SSRlags1,ResSSR$Model)  
  if (!ResSSR$Lreject){ 
    SSRlags2 <- c(SSRlags2,ResSSR$Model)
  }

  #---Store this SSR in matrix of SSRs
  irow <- ResSSR$yhat[,1]-yrgoY1[1]+1 # target rows in Y1
  Y1[irow,n]=ResSSR$yhat[,2,drop=FALSE]
  
  # Update progress if have run through 20 
  iprogress <- iprogress + 1
  if (iprogress >= 20){
    # Update progress bar
    mssgProg <- paste(as.character(n),'/',as.character(nchron2),' SSR models completed...',sep='')
    pctDone <- ProgTrack(pfProg,mssgProg,pctDone,pctInc)
    iprogress <- 0 # re-initialize counter
    } else {
  }
}

# Progress bar message about SSR modeling complete
if (pctDone>85){
  error('Programming error: pctDone should be less than 85% after SSR modeling')
}
mssgProg <- "SSR modeling completed..."
pctDone <- 85; pctInc <-0 # Move progress to 85% when all SSR modeling done
pctDone <- ProgTrack(pfProg,mssgProg,pctDone,pctInc)

#---Trim off any all-NA rows of SSR matrix Y1
i1<- trimRowNA(Y1)
if (is.vector(i1)){
 Y1<-Y1[i1,,drop=FALSE]
 yrY1<- yrY1[i1,,drop=FALSE]
}
rm(i1)

mSSRdf1<-nrow(SSRdf1) # number of rows in the table



#--- Write SSR results for all chronologies models to a file using fprintf
pf1<-file.path(outputDir,paste("Table1-SSR1",".txt",sep=""))
if (file.exists(pf1)){file.remove(pf1)} # must remove old versio of file
fmt1<-"%4s%4s%8s%7s%5s%8s%8s%6s%9s%8s%6s%6s%7s%6s%6s%7s\n"
fmt2<-"%4d%4d%8s%7d%5d%8s%8s%6.2f%9.2G%8.2f%6.2f%6.2f%7s%6d%6d%7s\n"
TitleLine <- 'Table1-SSR1 - Statistics of single site reconstruction (SSR) models'
fprintf('%s\n\n',TitleLine,file=pf1,append=FALSE)
fprintf(fmt1,"N1","N2","Site","Goc","Endc", "Model","Sign", "R2a",  "pF", "REcv", "REa",
        "REb", "Refit", "Gor", "Endr", "Reject",file=pf1,append=TRUE)
for (n in 1:mSSRdf1){
  fprintf(fmt2,SSRdf1[n,1],SSRdf1[n,2],SSRdf1[n,3],SSRdf1[n,4],SSRdf1[n,5],
          SSRdf1[n,6],SSRdf1[n,7],SSRdf1[n,8],SSRdf1[n,9],SSRdf1[n,10],
          SSRdf1[n,11],SSRdf1[n,12],SSRdf1[n,13],SSRdf1[n,14],SSRdf1[n,15],
          SSRdf1[n,16],file=pf1,append=TRUE)
}
fprintf('%s\n\n','',file=pf1,append=TRUE)
fprintf('%s\n','This table applies to chronologies before screening for hydrologic signal',
        file=pf1,append=TRUE)
fprintf('%s\n',paste('Chronologies:',PWtext),file=pf1,append=TRUE)
fprintf('%s\n',PdfDescribe,file=pf1,append=TRUE)


#--- Make second table, SSRdf2, a data frame that is a subset of SSRdf1 with just those
# chronologies not rejected according to the Lreject criterion.
Lreject<-SSRdf1[,16]

# If no chronologies have stable signal for flow, bail, with suggestions to user
if (all(Lreject)){
  emssgThis<- paste('ReconAnalog aborted:',
                    '\nNo chronology has a stable signal for ', HydName,'.',
                    '\nSome things you can try: ',
                    '\n1) turn off \"Reject non-stable\" in "Reconstruction Model Specifications" section;',	### AlexP change
                    '\n2) use a different climate variable or season;',
                    '\n3) try a different climate polygon or screening of the tree-ring network.',
                    sep='')
  eBomb<-emssgUNH(emssgThis,outputDir)
  stop(eBomb)
}

# If only one chronology has as stable signal and you have called for PCA, punt
nGood <- sum(!Lreject,na.rm=TRUE) 
if ((nGood == 1) && (methMSR==2 | methMSR==3) ){
  emssgThis<- paste('ReconAnalog aborted:',
                    '\nOnly one chronology has a stable signal for ', HydName,'.',
                    '\nThe selected multivariabe reconstruction method does not apply. You can try this: ',
                    '\n1) to get a reconstruction from this one site, select \"SLR\" as reconstruction \"Method\"',	
                    '\n2) try a different climate variable or season;',
                    '\n3) try a different climate polygon or screening of the tree-ring network.',
                    sep='')
  eBomb<-emssgUNH(emssgThis,outputDir)
  rm(nGood)
  stop(eBomb)
}


ix3 <- ix2[!Lreject] # col-pointer of "accepted" series to original tree-ring matrix
nms3 <- nms2[!Lreject] # ids of series passing screening
SSRdf2<-SSRdf1[!Lreject,]
mSSRdf2<-nrow(SSRdf2) # number of rows in the table
j1<- 1:mSSRdf2
SSRdf2[,1]=j1

pf2<-file.path(outputDir,paste("Table2-SSR2",".txt",sep=""))
if (file.exists(pf2)){file.remove(pf2)} # must remove old version of file
fmt1<-"%4s%4s%8s%7s%5s%8s%8s%6s%9s%8s%6s%6s%7s%6s%6s%7s\n"
fmt2<-"%4d%4d%8s%7d%5d%8s%8s%6.2f%9.2G%8.2f%6.2f%6.2f%7s%6d%6d%7s\n"

TitleLine <- 'Table2-SSR2 - Statistics of screened single site reconstruction (SSR) models'
fprintf('%s\n\n',TitleLine,file=pf2,append=FALSE)

fprintf(fmt1,"N1","N2","Site","Goc","Endc", "Model","Sign", "R2a",  "pF", "REcv", "REa",
        "REb", "Refit", "Gor", "Endr", "Reject",file=pf2,append=TRUE)
for (n in 1:mSSRdf2){
  fprintf(fmt2,SSRdf2[n,1],SSRdf2[n,2],SSRdf2[n,3],SSRdf2[n,4],SSRdf2[n,5],
          SSRdf2[n,6],SSRdf2[n,7],SSRdf2[n,8],SSRdf2[n,9],SSRdf2[n,10],
          SSRdf2[n,11],SSRdf2[n,12],SSRdf2[n,13],SSRdf2[n,14],SSRdf2[n,15],
          SSRdf2[n,16],file=pf2,append=TRUE)
}
fprintf('%s\n\n','',file=pf2,append=TRUE)
fprintf('%s\n','This table applies to chronologies passing screening for hydrologic signal',
        file=pf2,append=TRUE)
fprintf('%s\n',paste('Chronologies:',PWtext),file=pf2,append=TRUE)
fprintf('%s\n',PdfDescribe,file=pf2,append=TRUE)

#--- Make time series matrix of SSRs passing signal for screening (n0n-rejects)
Y2 <- Y1[,!Lreject,drop=FALSE]
yrY2 <-yrY1

#---Trim off any all-NA rows of SSR matrix Y2
i1<- trimRowNA(Y2)
if (is.vector(i1)){
  Y2<-Y2[i1,,drop=FALSE]
  yrY2<- yrY2[i1,,drop=FALSE]
}
mY2<-dim(Y2)[1] # number of rows in Y2
jScreened=SSRdf2$N2 # pointer from cols of screened SSRs
# to columns in original users tree-ring network

Tcpu<- (proc.time()-timeGo)[1] #...... processing time
Tclock<- (proc.time()-timeGo)[3] #...... clock time


################################################################################
#
# ANALOG EXTENSION OF SCREENED SSR'S ON RECENT END
#
# This extension is included to allow use of as long a calibration period as possible
# for the MSR model when tree-ring chronologies have variable ending year. For an 
# example of why this extension might be useful, consider a matrix of 10 tree-ring
# chronologies, 9 of which end in year 2022 and the other in 1998. Say you are using
# a PCA-based reconstruction method, which dictates that the ending year of the 
# reconstruction is no later than the earliest-ending chronology. The method described
# here works with the SSRs, or the single-site reconstructions generated from individual
# chronologies, and extends the SSR of the earliest-ending SSR so that it ends in 2022 instead 
# of 1995. This extension gives 28 additional years to the calibration period of the
# MSR model. 
#
# Have, say, N SSRs. Have a target end year that you want all series to cover.
# The first step is to find the common period for all N SSRs and compute the 
# correlation matrix (Spearman) for that common period. 
# Loop over all N SSRs, each time defining the current series as "key" and the
# rest as "others." Loop over key series: (1) Spearman r of key with others, and
# sorting of others from most correlated to least. (2) Proceed for next steps going
# from most to least correlated. (3) Pull full overlap of the two series -- this could
# generally be longer than the common period used for the Spearman correlation matrix.
# (4) Loop over the years of key series needing filled in. (5) Fill in all that are possible 
# from this member of others. (6) If still values to fill, proceed over more
# of the others, in sequence from the member of others most correlated with key to the
# member least correlated. The analog values is computed as folows.
#
# Analog method used, assuming have a key series and a predictor and the data 
# for the full overlap of the two series. Sort the two series from smallest to largest for
# that overlap. Have the value, x, of the predictor series for the year with data missing
# at key series. Compute the non-exceedance probability of that value in the predictor 
# series for the overlap. Assign the estimate as the interpolated value of the 
# kwy series with the same non-exceedance probability in the overlap. 

# Truncate matrix Y2, yrY2 on early end so that first row has no NA. After that truncation,
# Y2 should have no NAs on early end, but generally will have NAs in some cols on recent end
i1 <- which(complete.cases(Y2))
Y2 <- Y2[i1[1]:mY2,,drop=FALSE]
yrY2 <- yrY2[i1[1]:mY2,,drop=FALSE]

#--- Call function to extend tree-ring matrix on recent end
ResME <- tsmExtend(Y2,yrY2,yrsp1,N1,N2) # returns named list with Y, yrY

#--- Bomb out messages from tsmExtend()
emsgs2 <- c('No need to extend','OK, but yrsp later than last year of data in X',
            'No problem','tsmExtend aborted: common period of all series too short',
            'tsmExtend aborted: insufficient common period of predictor and predictand',
            'tsmExtend aborted: yrX not continuous','tsmExtend aborted: first year of X has some NA')
khow<-ResME$khow
if (khow>3){
  emsg2<-emsgs2[[khow]]
  eBomb<-emssgUNH(emsg2,outputDir)
  stop(eBomb)
}

Y3 <- ResME$Y; yrY3 <- ResME$yrY # forward-extend tree-ring matrix
if (any(is.na(Y3))){
  eBomb<-emssgUNH('ReconAnalog() aborted: matrix Y3 has a NA',outputDir)
  stop(eBomb)
} 

# Revision 2023-05-16, th en re-revised 2023-11-22, to set the end year for the reconstructed y as the last year
# of the quantile-extended matrix of SSRs. The user should be aware that the recent end of the 
# reconstruction could be be based on statistically extended tree-ring data for some chronologies.
# This could be very bad if, for example, the most recent chronology has a weak SSR signal. One of the
# SSR plots shows the drop in strongest SSR signal as chronologies drop out. This plot can be used
# as a guide for the last you should trust the reconstruction.
yrEnd <- yrY3[length(yrY3)]


################################################################################
#
# FIGURES FOR SINGLE-SITE RECONSTRUCTIONS
#
# Idea is that in TRISH user can use radio button choose from 1) bar chart of site-screening,
# 2) boxplot summaries of model statistics, 3) line plot of time coverage of SSRs, and
# 4) z-score time-series plot (annual and smoothed on one set of axes; mean of the SSRs 
# after converting to z-scores)
n2 <-nchron # number of sites returned by polygon screening
n3 <- length(ix2) # number of chronologies after screening for reconstruction window
n4 <- dim(Y3)[2] # number of chronologies passing final screening for hydrologic signal


#----FIGURE 01.  1x2 OF 1)BAR PLOT SUMMARIZING NUMBER OF CHRONOLOGIES, ORIGINAL AND SCREENED,
# AND 2) ADJUSTED R-SQUARED OF ALL FITTED SSR MODELS AND OF THOSE PASSING FINAL
# SCREENING FOR HYDROLOGIC SIGNAL


#--- Numbers of chronologies (bar chart)

# Preliminaries
xtemp <- c(n1,n2,n3,n4)
DeltaTemp<- 0.1*(max(xtemp)-min(xtemp))

#yhi <- 0.1*(max(xtemp)-min(xtemp)) + max(xtemp)
yhi <- DeltaTemp + max(xtemp)
ylo <- 0
ylims <- c(ylo,yhi)
ylims2 <- ylims
ylims2[2]<-ylims2[2]+0.4*ylims2[2]
DeltaTemp2 <- ylims2[2]/30

barnames <- c('N1','N2','N3','N4')

strAnnote1<-paste('\nN1: source network','\nN2: polygon-selected',
               '\nN3: SSR-modeled','\nN4: Passed screening')

png(filename=paste(outputDir,"Figure01-SSR1.png",sep=""), width = 960, height = 480)
layout.matrix <- matrix(c(1,2), nrow = 1, ncol = 2)
layout(layout.matrix,heights=2,widths=c(2,1))
par(mar=c(5,4,4,4),cex.main=1.3)

bp <- barplot(xtemp,ylim=ylims2,xlab='Screening Stage',col='Pink',border=TRUE,
              names.arg=barnames,main='Number of Chronologies',cex.lab=1.3)
text(bp,xtemp+DeltaTemp2,labels=xtemp)
abline(h=0)

# Annotate meaning of N1, N2, N3, N4
xtemp <- 3.65; ytemp <- ylims2[2]
text(xtemp,ytemp,strAnnote1,adj=c(0,1),cex=1.2)

rm(xtemp,barnames,ylims,ylims2,bp,DeltaTemp2)

#--- Adjusted R squared

#screen(2)
par(mar=c(5,4,4,1),cex.main=1.3)

namesBP<-c('N3','N4')
boxplot(SSRdf1[,8],SSRdf2[,8],notch=FALSE,
        ylab = "Adj R-squared of SSR Models",
        main="Adjusted R-squared",names=namesBP,
        cex.lab=1.2)
dev.off()

#----FIGURE 02.  1x2 OF HISTOGRAMS OF WHICH LAGS ARE IN THE SSR MODELS. 1) ALL SSRS, AND
# 2) THOSE SSRS PASSING SCREENING FOR HYDROLOGIC SIGNAL. EACH SSR MODEL MAY HAVE
# FROM 1-5 LAGS, RANGING FROM t-2 to t+2 YEARS FROM THE YEAR OF FLOW. THE HISTOGRAMS
# SUM OVER MODELS, SO THA THE GRAND TOTAL NUMBER OF LAGS IN THE HISTOGRAM IS
# GREATER THAN THE NUMBER OF MODELS.


#--- Left: Histogram of lags, all N3 SSRs 

png(filename=paste(outputDir,"Figure02-SSR2.png",sep=""), width = 960, height = 480)
layout.matrix <- matrix(c(1,2), nrow = 1, ncol = 2)
layout(layout.matrix,heights=2,widths=c(1,1))

# Left histpgram
par(mar=c(5,4,4,1),cex.main=1.4)
hBreaks<-c(0.5,1.5,2.5,3.5,4.5,5.5)
# Changing x axis
xTicks<-seq(1, 5, by=1)
xTickLabs <- c('-2','-1','0','+1','+2')
n3a<-length(SSRlags1)
title1<-paste('Histogram of Lags (',n3,'N3 Models,',n3a,'Total Lags)')
hist(SSRlags1,breaks=hBreaks,xlim=c(0.5,5.5),xaxt='n',
     main=title1,xlab='',cex.lab=1.2)
vtemp<-par('usr')
mtext('Lag (yr)',side=1,line=1.5,cex=1.2)
for (k in 1:5){
  mtext(xTickLabs[k], side = 1, line = 0, outer = FALSE, 
      at = k, adj = NA, padj = NA, cex = NA, col = NA, 
      font = NA)
}


#--- Right: Histogram of lags, SSRs passing screening for hydro signal

ylims<-c(vtemp[3],vtemp[4]) # same y limits as previuous
par(mar=c(5,4,4,1))
hBreaks<-c(0.5,1.5,2.5,3.5,4.5,5.5)
# Changing x axis
xTicks<-seq(1, 5, by=1)
xTickLabs <- c('-2','-1','0','+1','+2')
n4a<-length(SSRlags2)
title2<-paste('Histogram of Lags (',n4,'N4 Models,',n4a,'Total Lags)')
hist(SSRlags2,breaks=hBreaks,xlim=c(0.5,5.5),ylim=ylims,xaxt='n',
     main=title2,xlab='',yaxs='i',cex.lab=1.2)
mtext('Lag (yr)',side=1,line=1.5,cex=1.2)
for (k in 1:5){
  mtext(xTickLabs[k], side = 1, line = 0, outer = FALSE, 
        at = k, adj = NA, padj = NA, cex = NA, col = NA, 
        font = NA)
}
dev.off()

#----FIGURE 03. 1x1. TWO-Y-AXIS TIME PLOT TO HELP GUIDE USER IN CHOICE OF 
# CALIRATION PERIOD FOR MULTI-SITE-RECONSTRUCTION (MSR) MODEL. 
# 
# Before extension by tsmEndtend, the signal-screened SSRs generally end in
# different years. These time plots cover the tail years of the SSRs: from 
# the last year with data for all SSRs through the last year of the SSR with
# most recent data. One plot is the maximum adjusted R-squared of available
# SSRs. The other plot is the number of available SSRs. Annotated on the plot
# is also the ending year of the observed hydro variable. The user will eventually
# need to select the end year of calibration of the MSR model. This year
# cannot be later than the last year of the observed hydro series. The year could
# be as late as the last year SSR at any site, but this may not be a good
# idea if the adjusted R-squared of the most recent SSR is small. 
#
# It is also possible that the last available hydro data for calibration ends 
# before the ending year of the earliest-ending SSR. In that case this plot of drop
# in R-squared as chronologies drop out toward the present be just a curiosity because
# the end year of the calibration period cannot be later than the end year of the 
# hydro series

#--- Prepare data needed for the plots
#
# Status. Have full-length observed hydro in 1-col matrices, yrv, v. 
# Have the adjusted R-squared of signal-screened models in SSRdf2[,8].
# Have the corresponding end years of the SSRs in  SSRdf2$EndR. 
# Have time series of SSRs in matrices Y3,yrY3
ResSD <- SignalDrop1(SSRdf2$EndR,SSRdf2[,8])
x1 <- ResSD$x1;   x2 <- ResSD$x2;  yrx1 <- ResSD$yrx1; yrx2 <- ResSD$yrx1; 

# Want x axis to start year before and end year after the relevant period
yrsEnd <- unique(SSRdf2$EndR) # unique ending years of screened SSRs
xHead <- max(yrv) # head of arrow here; also is last year of available hydro series
#yrLo <- min(yrsEnd)
yrLo <- min(min(yrsEnd),xHead-5) # xHead-5 in case hydro series ends before end of any SSR
yrHi <- (xHead - yrLo) + xHead
yrHi <- max(yrHi,max(yrsEnd)+2)
xlims <- c(yrLo-0.05,yrHi+0.05) # limits for x axis
ylims <- c(min(x1)-0.2,max(x1)+0.2) # limits for y axis

# Arrow
yHead <- ylims[1]+ diff(ylims)/2
yTail <- yHead + diff(ylims)/10
xTail <- xHead+ (xlims[2]-xHead)/2 # x position of tail

#--- yyplot
png(filename=paste(outputDir,"Figure03-SSR3.png",sep=""), width = 960, height = 480)
par(mar=c(5,5,4,5),cex.main=1.4)
plot(yrx1,x1,xaxt='n',yaxt='n',type="b",pch=1,col='red3',xlim=xlims,ylim=ylims,cex=1,
     main='Drop in Signal Strength with Loss of SSRs (Chronologies) on Recent End',xlab='Year',
     ylab='Number of SSRs',cex.lab=1.2)
# vertical line at last year of hydro series
axis(1,at=seq(yrLo,yrHi,1))
axis(2,at=seq(min(x1),max(x1),1))
abline(v=max(yrv),col='Magenta',lty=2)

# arrow to the vertical line
arrows(xTail,yTail,xHead,yHead)
txtTemp<- paste('Last year of',Dtype)

# annotate seasons and end year of specified calibration of recon models
# Calibration period of SSR could end sooner if chronology ends sooner, but 
# period not allow to extend more recent than specified year. 

# String for calibration period end year. NA indicates let the data
# determindx the end year, which will depend on then end year of chronology, end year
# of y, and whether lags allowed. 
strEnd <- paste(as.integer(yrspc),'=','specified calib. end year',sep=' ') 

# String for latest feasible end year of calibration period. Depends on on end years of
# y and of chronologies, and also on whether lags allowed. The key chronology for this determination
# is the one ending most recently, because ReconAnalog automatically extends SSRs by a quantile 
# extension algorithm to the year of the most-recent ending SSR. The actual calibration periods for
# individual SSRs can end earlier than this last feasible end year, because the function reconsw4
# will trim back the calibration period end year if is not possible given the end year of the
# chronology and the optional lagging. 
if (LallowLags){
  yrspF <- yrspcLimit # last feasible end year for SSR calibration period; 
  strF1 <- '   (assuming lags allowed)'
}else{
  yrspF <- yrspcLimitv # last feasible end year for SSR calibration period; 
  strF1 <- '   (assuming lags allowed)'
}
strEndF <- paste(as.integer(yrspF),'=','latest feasible calib. end year',sep=' ') 


txtTemp2 <- paste('\n ',HydroVariable,'=',HydName,
   '\n ',LabelSeason1,
   '\n ',LabelSeason2,
   '\n ',strEnd,
   '\n ',strEndF,
   '\n ',strF1)
text(xTail,yTail,txtTemp,adj=c(0,0.5),cex=1.4)
text(xHead,yHead,txtTemp2,adj=c(0,1),cex=1.3)
  
par(new = T)
plot(yrx2,x2, pch=16, type='b',lty=0,axes=F, xlab=NA, ylab=NA, xlim=xlims,
     cex.lab=1.3,cex.axis=1.2)
axis(side = 4)
mtext(side = 4, line = 3, 'Maximum adjusted R-squared')
legend("topright",
       legend=c("N of SSRs","Max Adj R-squared"),
       lty=c(1,0), pch=c(1, 16), col=c("red3", "black"),cex=1.3)
dev.off()
rm(strEnd,strEndF,strF1)
#rm(txtTemp,txtTemp2,xTail,xHead,yTail,yHead,yrsEnd,yrLo,yrHi,xlims,ylims,x1,x2,yrx1,yrx2)



#----FIGURE 04. 1x1. SCATTERPLOT OF OBSERVED HYDRO (SEASONALIZED HYDROLOGIC SERIES) ON THE
# MEAN OF THE SIGNAL-SCREENED SSRS. THIS PLOT WILL HELP USER DECIDE WHAT TO CHOOSE AS METHOD
# FOR THE MULTI-SITE RECONSTRUCTION (MSR)
#
# Status. 
# Y3, yrY3 are matrices with the screened SSRs
# V is matrix of hydro series, with year as col 1 and data as col 2

# Prepare the two series for the scatterplot
w <-rowMeans(Y3, na.rm=TRUE) # average of the SSRs (vector)
yrw = yrY3 # 
W<- as.matrix(cbind(yrw,w)) # bind yrw and w into a time series matrix
ResPC<- PeriodCommon(W,V)# get common period of W and V
yrgo1 <- ResPC$tgo; yrsp1 <- ResPC$tsp # start and end years of common period 
w1 <- ResPC$X[,2]; yrw1 <- ResPC$X[,1]; # mean-SSR and its year, as vectors, for common period 
# with observed hydro
v1 <- ResPC$Y[,2]; yrv1 <- ResPC$Y[,1]; # observed hydro and its year, as vectors, for common period 

# Pearson correlation of observed hydro with mean of SSRs
r = cor(v1,w1)
rStr<- paste('r=',toString(round(r,digits=2)))
strMain=paste('Scatter of Mean of',as.character(n4),
              'Single-Site Reconstructions (SSRs) of',Dtype,'on Observed',Dtype,
              '\n(Fits are straight-line (red) and loess (blue))') # plot title
rm(r)

# Strings for plots
ylabTemp <- paste('Mean SSR',LabUnits)


#--- scatterplot
png(filename=paste(outputDir,"Figure04-SSR4.png",sep=""), width = 760, height = 480)
par(mar=c(5,5,6,4))
# Call function from package car for the scatterplot. This function as called will regress w1 on
# v1 and plot the lease-squared-fit straight line. Also plotted is a loess (local regression) curve
# and loess curves to represent variance of the loess estimate. The loess curves use a a span of 2/3 and
# are estimated by R function loess.R. For the error bars, the negative and positive departures from
# the loess estimate of the mean are squared and themselves fit wit a loess curve. The plotted lines are
# at the square root of those squared-departure fits. Because the two bordering smoothed line represent the
# typical positive and negative departure, the confidence interval can be heuristically considered 
# a 50% confidence interval.
scatterplot(v1,w1,boxplots=FALSE,
            regLine=list(method=lm, lty=1, lwd=2, col="red"),
            ylab=ylabTemp,
            xlab=paste('Observed',Dtype,LabUnits),
            main=strMain,cex.lab=1.2)
text(min(v1),max(w1),rStr,adj=c(0,1),cex=2)
dev.off()
rm (ylabTemp)


################################################################################
#
# COMBINE SSR'S INTO A FINAL MULTI-SITE RECONSTRUCTION (MSR
# Method depends on settings for methMSR and PCApredictors. If PCA is involved in
# reconstruction, method might also depend on settings of nkHowPCA, PCoption 
# and nPCsKeep
#
# methMSR has three possible values: (1) simple linear regression, 
# (2) stepwise multiple linear regression on SSRs (PCApredictors=false) or their 
# PCs (PCApredictors=true), and (3) PCA analog
#
# Simple linear regression is done by calling function RecSLR1
# The other methods ared done by calling function RecMLR1

ReconMethods <- c("Simple linear regression of y on mean of SSRs", 
                  "Multiple linear regression of y on SSRs or their PCs",
                  "Principal component analog nearest neighbor")
ReconMethod <- ReconMethods[methMSR]
NextFigNumber<-5 # because SSR has already produced figure files Figure01.png to Figure04.png

# Set calibraton period of MSR to start with latest of [yrgoc;d first available year of hydro series; first available year of mean SSR[
# Set period to end with earliest of [yrgoc; last available year of hydro series' last available year of mean SSR]
if (methMSR==1){
  # Recon by simple linear regression, using RecSLR1(). 
  SpecListSLR1 <- list("PdfDescribe"=PdfDescribe,"Text"=RecListx,"u"=w,"yru"=yrw,
                       "v"=V[,2],"yrv"=V[,1],"yrsC"=yrsCalWindow,"yrEnd"=yrEnd,"nNeg"=nNeg,"nPos"=nPos,
                       "NcMin"=N3,"NextFigNumber"=NextFigNumber,"outputDir"=outputDir)
  # Until RecSLR1() finished, will not code the call with the above inputs; will just have that
  # function read MyData from wd.
  #save(SpecListLR,file="MyData")
  #---- stop
  
  # Progress bar update
  mssgProg <- "SSRs prepared for multi-site reconstruction modeling..."
  pctInc <- 10;
  pctDone <- ProgTrack(pfProg,mssgProg,pctDone,pctInc)
  
  
  # Call to function for method RecSLR1
  save(SpecListSLR1,file=paste(outputDir,"MyDataSLR1.dat",sep="")) 
  Z <- RecSLR1(SpecListSLR1)
  if (Z$flag>0){
    emssg<-Z$Msg
    ResTemp<-emssgUNH(emssg,outputDir)
    stop(emssg)
  }
} else if (methMSR==2 | methMSR==3){
  # Recon by regression on sreened SSRs or their PCs, with call to RecMLR1 
  #save(SpecListMLR1,file="MyData"). Note that lags have been dealt with at the SSR step. No
  # lags are included in the MSR model. But, nPos and nNeg are used in the MSR model to set
  # m in leave-m-out cross-validation, on grounds that the SSRs did use lagging.
  SpecListMLR1 <- list("Text"=RecListx,"U"=Y3,"yrU"=yrY3,"nmsU"=nms3,"jScreened"=jScreened,"v"=V[,2],"yrv"=V[,1],
                     "yrsC"=yrsCalWindow,"yrEnd"=yrEnd,"nNeg"=nNeg,"nPos"=nPos,"incR2a"=incR2a,"kstop"=kstop,
                     "NcMin"=N3, "PCoption"=PCoption,"f"=f,"PCApredictors"=PCApredictors,
                     "methMSR"=methMSR,"PdfDescribe"=PdfDescribe, "nPCsKeep"=nPCsKeep,"alphaR"=alphaR,
                     "ScreenAnalogPCs"=ScreenAnalogPCs, "kHowPCA"=kHowPCA,"NextFigNumber"=NextFigNumber,
                     "outputDir"=outputDir)
  save(SpecListMLR1,file=paste(outputDir,"MyDataMLR1.dat",sep=""))
  Z <- RecMLR1(SpecListMLR1)
  if (Z$flag == 1){
    emssg<-Z$Msg
    ResTemp<-emssgUNH(emssg,outputDir)
    stop(emssg)
  }
}
# Progress bar update
mssgProg <- "Reconstruction modeling completed!"
pctInc <- 0;
pctDone <- 100
pctDone <- ProgTrack(pfProg,mssgProg,pctDone,pctInc)

print("Done!")			### AlexP addition to print status to the STDOUT
