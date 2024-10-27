reconsw4<- function(X,Y,nNeg,nPos,yrsCalWindow,c1,Lcausal,LallowLags,MinCalibLength,RequireStable){
  # reconsw4(): single-site reconstruction (SSR) by stepwise forward distributed-lag regression
  # D. Meko; last revised 2024-04-17
  #
  # When LallowLags=True, for predictand y and chronology x, when LallowLags=True, code 
  # assumes a potential-predictor pool of x(t-2) to x(t+2) for predicting y(t). If all 
  # predictors were to enter, the final model would have 5 predictors: x(t-2), x(t-1), x(t),
  # x(t+1) and x(t+2). Predictors enter forward stepwise until the adjusted R-squared
  # fails to increase by c1 (see below) from the previous step. The predictor most highly 
  # correlated with y(t) enters first. Additional predictors enter in order of decreasing
  # correlation with regression residuals from the previous step. At each step the model 
  # is cross-validated by leave-9-out cross-validation, and tested for stability by 
  # split-sample validation. The reduction-of-error statistic is computed as the skill
  # statistic in these validation exercises. 
  #
  # The maximum adjusted R-squared model is then simplified, or "cut back" to the step
  # of maximum validation skill, as judged by the maximum computed reduction-of-error
  # statistic. The model is then accepted or rejected based on several criteria. 
  # Rejection results from ANY of the following being true:
  #   1 REcv<=0:   no skill of validation by leave-9-out cross-validation
  #   2 REa <=0 or REb <=0:  no skill in either of the split-sample validation tests
  #   3 pF >=0.05: no significant calibration skill, as judged by p-value of overall F
  #   4 Illogical lagged model: final model implies that y(t) can be predicted from
  #     just the past values of x(t). Logic, by example, is that this years runoff logically
  #     should not be predictable from a model with only past years' tree rings as predictors.
  #
  # Model input parameters Lcauasal and LallowStable allow overriding of the rejection by
  # criteria 4 and 2 above (see comments on those input parameters)
  #
  # When LallowLags=false, the predictor pool for y(t) includes just x(t), and the above
  # stepwise description can be ignored. In this case, the model is simple linear regression,
  # and of course "Lcausal" has no effect. But the screening is still guided by skill 
  # requirements for cross-validation and split-sample validation.
  #
  # The broadest possible lagging of t-2 to t+2 is critical in this function because code builds
  # variables such as a 5-character pointer output coding which lags are in the final model. So,
  # the function will bomb if inputs nNeg and nPos inputs are not both set to 2. 
  #
  #=== INPUT
  #
  # X [matrix]r tree-ring index as col 2, year as col 1; may have leading and trailing NA
  # Y [matrix]r flow (or some other predicand) time series in col 2, year in col 1
  #   Time series in Y and X will need to overlap by at least 30 years
  # nNeg [1x1]i  maximum number of negative lags allowed. The only acceptable setting
  #   is nNeg=2, which means consider lags t-2 and t-1 relative to flow as possible predictors.
  # nPos [1x1]i  maximum number of positive lags allowed.  The only acceptable setting
  #   is nPos=2, which means consider lags t+2 and t+1 relative to flow as possible predictors.
  #   Note that you can over-ride the input settings for nNeg and nPos with the input
  #   argument LallowLags, such that effectively no lags are allowed in the model. 
  # yrsCalWindow (1x2)i  calibration will make use of flow data only within this period
  #   designated by first and last year
  # c1 (1x1)r critical necessary increment of adjusted R-square and cross-validation reduction of error (REcv)
  #   required to include an anddition step in stepwise regression. Models are first fit up to the step
  #   at which R2adj has increased by at least c1 from prveious step. After leave-9-out cross-validation, 
  #   model is possibly further simplified such that last step must yield an increase of at least
  #   c1 (e.g., 0.01) in REcv. This threshold is in interest of parsimony, to avoid a more complicated model
  #   if practical gain in accuracy and skill is negligible.
  # Lcausal (1x1)L  TRUE if reject any model that has negative lags only on the tree-ring series. Makes most
  #   sense to do this if using standard chronologies, but not for residual chronologies, as a negative lag might
  #   be compensating for over-whitening or under-whitening. 
  # LallowLags (1x1)L TRUE if allow lagged model; FALSE if force model to be lag-0 only
  # MinCalibLength (1x1)i minimum allowable length of calibration period; if overlap of X and Y is too short
  #   for this, the function bombs with an error message. The error message is also written to the output
  #   directory "/outputDir/", which must be in the global space (e.g., from calling script)
  # RequireStable (1x1)L TRUE if model tagged for "rejection" when the the RE statistic is negative for 
  #   either half of the split-sample validation. FALSE if the the split-sample RE is to be ignored in
  #   the decision to reject
  # MinCalibLength (1x1)i minimum allowable length of calibration period; if overlap of X and Y is too short
  #   for this, the function bombs with an error message. The error message is also written to the output
  #   directory "/outputDir/", which must be in the global space (e.g., from calling script)
  #
  #
  #=== OUTPUT
  #
  # Returns named list Output with following elements:
  #   Model [vector]i columns of [t-2 t-1 t t+1 t+2] lagged tree-ring matrix in final model, L-to-R
  #     in order as the variables entered stepwise. For example [5 2] means model has t+2 and t-1
  #     as predictors, and that t+2 (element 5) entered first.
  #   ModedCoded [string] 1x5 string showing which of the five potential predictors are in the
  #     final model, and order that they entered. This could be used in a supplemental table. For
  #     example, "02001" indicates lags t-1 and t+2 are in the model, and that t+2 entered before
  #     t-1
  #   ModelSign [string] 1x5 string with "0" if lag not in model, and "P" or "N" indicating lag
  #     in model with coefficient of positive or negative sign
  #   yearsCal(1x2)i  first and last year of calibration
  #   yearsRec(1x2)i first and last year of reconstruction
  #   MaxLagNegPos (1x2)i maximum negative and positive lags considered in modeling
  #   LeftOutCV (1x1)i  number of observations left out in cross-validation
  #   IncrementR2adj (1x1)r critical threshold for meaningful "increase" in adjusted R square
  #   RegCoefs (vector)r regression coefficients, constant term first
  #   Rsquared (1x1)r R-squared of calibration
  #   RsquaredAdj (1x1)r adjusted ....
  #   F (1x1)r overall F of regression
  #   pF (1x1)r p-value of overall F
  #   Lsig (1x1)L final model overall-F significant AT 0.05
  #   REcv (1x1)r reduction of error statistic for leave-9-out cross-validati0on
  #   REa (1x1)r RE for model calibrated on first half, validated on second
  #   REb (1x1)r RE for model calibrated on second half, validated on first
  #   RMSE (1x2)r  root mean square error of calibration (1) and cross-validation (2)
  #   LREcut (1x1)L final model was truncated, or cut off, at a lower lag than indicated
  #     by maximum adjusted R-squared because cross-validation RE was higher at a lower lag.
  #   Lrefit (1x1)L final model was re-fit to a longer calibration period
  #     because the calibration/validation procedure resulted in a model
  #     with fewer than two positive lags.
  #   LNegOnly(1x1)L TRUE is negative lags only in final model; FALSE otherwise
  #   Lreject (1x1)L   if true, tree-ring chronology rejected because at least one 
  #     screening criterion not satisfies. If false, chronology not rejected. Screening
  #     criteria always include a significant calibration signal (p-value of overall F) 
  #     a positive skill of cross-validation. Additional criteria apply if Lcausal=true
  #     or RequireStable=true (see above)
  #   Yh [matrix] reconstruction; 2-col matrix, with year in first column 
  #
  # revised 20221213: adding option LallowLags
  # revised 20230325: added input arg MinCalibLength and code that gives error message if calibration
  #   interval for model (in terms of y years) is fewer than MinCalibLength observations
  # revised 2023-04-09: to allow the requirement for temporal stability of models to be optional.
  #   See RequireStable
  # revised 2023-05-26: call to new function, LagModel2Sign, to get new list item, 
  # revised 2024-03-04: (1) fixed error in call to LagModel2Sign for coding the lags in model to
  #   sign of regression coefficients, (2) added commented-out debug code that can be uncommented
  #   for future debugging lagged regression [search "#edebug"], (3) fixed error in getting 
  #   lagged tree-ring matrix ready for regression refit (see Z<- Xorig[,ColsInModel,drop=FALSE] )
  # revised 2024-04-17: comments clarified and corrected for typos
 
  source(paste(code_dir,"LagYear.R",sep="")) # build lagged matrix of chrons
  source(paste(code_dir,"LagReOrder.R",sep="")) # build lagged matrix of chrons
  source(paste(code_dir,"CrossValid1.R",sep=""))  # leave-m-out cross-validation of stepwise models
  source(paste(code_dir,"CrossValid2.R",sep=""))  # leave-m-out cross-validation of stepwise models
  source(paste(code_dir,"PeriodCommon.R",sep="")) # common period of overlap of chrons and predictand
  source(paste(code_dir,"ForwStep2.R",sep="")) # forward stepwise regression
  source(paste(code_dir,"ssValid.R",sep="")) # split-sample cross-validation
  source(paste(code_dir,"LagModel2Char.R",sep="")) # Build string representation of lagged model estimate forward stepwise
  source(paste(code_dir,"LagModel2Sign.R",sep="")) # Build string representation of signs of coefs in model
  source(paste(code_dir,"emssgUNH.R",sep="")) # write error file to system, specified output folder

  
  #--- CHECK INPUT

  # X and Y must be matrix and 2-column
  L<-is.matrix(X) && is.matrix(Y)
  if (!L) {stop('X or Y not matrix')}
  L<- dim(X)[2]==2 && dim(Y)[2]==2
  if (!L) {stop('X or Y must be 2-column')}
  
  # Maximum lag of 2: both nNeg and nPos must be 2
  # Be aware that LallowLags can override these settings and allow no-lag
  # models, but nNeg and nPos must still be input as 2.
  L <- (nNeg==2 && nPos==2) 
  if (!L){
    emssg <- 'reconsw4 requires that input nNeg=2 and nPos=2}'
    ResTemp<-emssgUNH(emssg,outputDir)
    stop(emssg)
  }
  
  
  # yrsCalWindow should be length 2 vector
  if (!length(yrsCalWindow)==2) {stop('yrsCalWindow should be length 2 vector')}
  
  # c1 should be greater than or equal to zero but no greater  than 0.05
  L <- (c1>=0) && (c1<=0.05) 
  if (!L) {stop('c1 should be non-negative and not greater than 0.05')}
  

  ##############################################################################
  #
  # BUILD MATRIX OF LAGGED TREE-RING INDEX, AFTER TRIMMING TO COMPLETED CASES5 
  
  a5<-c(1,2,3,4,5) # vector needed later; refers to cols 1-5 of lagged tree-ring matrix

    # LagYear() requires X as single-column matrix
  mX <- dim(X)[1]
  tGo <- X[1,1]
  tSp <- X[mX,1]
  yrX <- tGo:tSp

  ktrim <-2 # do not trim lagged matrix to exclude any leading or trailing years containing one or more NA;
  # But all-NA rows are trimmed off
  x <- matrix(X[,2])
  L <- complete.cases(x)
  x <- as.matrix(x[L]) # trim x to non-NA
  mlength <- dim(x)[1]
  yrX <- yrX[L]
  tGo <- yrX[1]
  tSp <- yrX[mlength]
  yrX <- tGo:tSp
  
  ResLag<-LagYear(x,tGo,tSp,nNeg,nPos,ktrim) # to build lagged matrix

  X <- ResLag$X # lagged tree-ring matrix: col order is first unlagged, then
  # increasing negative lags, then increasing positive lags. So, if nNeg=2 and nPos=2,
  # cols would be: t t-1 t-2 t+1 t+2
  tGo <- ResLag$tGo
  tSp <- ResLag$tSp; # start and end year of  returned lagged matrix X
  yrX <- tGo:tSp
 
  
  #--- Call function to re-order columns of lagged matrix. LagYear returned matrix columns
  # ordered no-lag, followed by negative lags and then positive lags with increasing lag 
  # left to right. For consistency with other functions, want cols of X arranged L to R 
  # from highest negative to highest positive lag. This to be consistent with what is expected 
  # in functions yet to be called.
  # After call to LagReOrder, if nNeg and NPos are 2, for example, X has cols ordered
  # t-2 t-1 0 t+1 t+2
  X <- LagReOrder(X)


  ##############################################################################
  #
  # STORE VERSIONS OF LAGGED TREE-RING MATRIX: ORIGINAL, AND INCLUDING ONLY
  # THOSE ROWS WITHOUT ANY NAN
  #
  # Because will use the latter in first pass of stepwise regression, where need
  # data for all five potential predictors, and depending on the final selected lags
  # may be able to use additional rows of X for final calibration.
  nX <- dim(X)[2]  # cols in lagged matrix
  Xorig <- X # full lagged matrix
  yrXorig <- as.matrix(tGo:tSp) # year col-vector for Xorig

  L <- complete.cases(Xorig)
  X <- Xorig[L,]  # Lagged tree-ring matrix, any row with a NA deleted
  yrX <- as.matrix(yrXorig[L])


################################################################################
#
# PREPARE MATRICES OF FLOWS AND LAGGED TREE-RINGS FOR STEPWISE

#--- PeriodCommon() to trim lagged tree-ring matrix and flow time series to same 
# time coverage. Flow data are already in matrix of desired form in V.
U<- cbind(yrX,X) # matrix of lagged chronology, with year as col 1
namesU<-c('XN2','XN1','X0','XP1','XP2')
ResPC <- PeriodCommon(U,V) 
U1<-ResPC$X
V1<-ResPC$Y

#--- Optional truncation of calibration years; user not allowed to specify a start
# year of calibration earlier than first year of overlap of lagged predictor matrix with y
# 
if (is.na(yrsCalWindow[1])){
  tGo <- ResPC$tgo
  } else if (yrsCalWindow[1]<ResPC$tgo){
  tGo <- ResPC$tgo
  } else {
    tGo <- yrsCalWindow[1]
}
if (is.na(yrsCalWindow[2])){
  tSp <- ResPC$tsp
} else if (yrsCalWindow[2]>ResPC$tsp){
  tSp <- ResPC$tsp
} else {
  tSp <- yrsCalWindow[2]
}

#--- Truncation of flows and tree-ring matrix for calibration
yrU1<-as.matrix(U1[,1])
U1<-as.matrix(U1[,-1])
yrV1<-as.matrix(V1[,1])
V1<-as.matrix(V1[,-1])

# V1 and U1 Ok as-is if the desired calibration period matches the overlap of y and lagged X
# Otherwise, must truncate.
L<- tGo==ResPC$tgo && tSp==ResPC$tsp
if (!L){
  tGo<-max(tGo,ResPC$tgo)
  tSp<-min(tSp,ResPC$tsp)
  L1<-yrV1>=tGo & yrV1<=tSp
  U1 <- as.matrix(U1[L1,]) ; yrU1<-as.matrix(yrU1[L1])
  V1 <- as.matrix(V1[L1,]) ; yrV1<-as.matrix(yrV1[L1])
}else{
  # no need for action
}

# Minimum allowable length of calibration period is MinCalibLength. Return error message
# and bomb out if this minimum length of calibration period violated
nCalib <- tSp-tGo+1 
Lbad <- nCalib < MinCalibLength
if (Lbad){
  emssg <- paste('Function reconsw4 requires calibration period of at least ', as.character(MinCalibLength),' years.',
     '\nOverlap for at least one tree-ring chronology allows only ', as.character(nCalib),' years.', sep='')
      ResTemp<-emssgUNH(emssg,outputDir)
  stop(emssg)
}
rm(Lbad)


################################################################################
#
# PRELIMINARY REGRESSION TO SET ORDER OF ENTRY OF VARIABLES OR MODEL FOR WHICH 
# ADJUSTED R-SQUARE IS APPROXIMATELY MAXIMUM; DEPENDING ON INPUT LallowLags, MAY SIMPLIFY
# BECAUSE ONLY LAG-0 IS POSSIBLE.
#
# Status. Flow and lagged tree-ring matrix of calibration length are in matrices U1,V1, with
# 1-col year matrices yrU1 and yrV1. tGo and tSp are the corresponting start and end years 
# for calibration. To go back to longer data, have available longer calib data in ResPC fields. 
# Have full--length lagged tree-ring matrix in Xorig, yrXorig for later use in reconstruction.

#--- Forward stepwise of flow on lagged tree rings

namesU<-c('XN2','XN1','X0','XP1','XP2') # matrices have been organized L to R from lag-2 to lag+2

# #debug on  
# write.csv(U1,file='U1.csv')
# write.csv(V1,file='V1.csv')
# # In MATLAB, can troubleshoot with simple script (see text_reconsw4.m), such as
# #   y = csvread('V1.csv',1,1)
# #   X = csvread('U1.csv',1,1)
# #   stepwise(X,y)
# #debug off


if (LallowLags){
  ResFS1<- ForwStep2(U1,namesU,V1,c1)
} else {
  ResFS1<- ForwMoot(U1,namesU,V1)
}

# $StepMaxR2adj -- the indicated final step
# $ColsInOrderEntry -- where {1 2 3 4 5} are lags t-2 t-2 t t+1 t+2
i1<-ResFS1$ColsInOrderEntry # Cols entered into max-R2a model, in order as entered
ni1<-length(i1)

# Store a few statistics, which could change later in response to re-fit with one or two more years on
# recent end if data allow (fewer than 2 positive lags in SSE model and y not limiting)
ModelPicked<-i1
ColsInModel <- ResFS1$ColsInModel # cols in model in order of RegCoefs (after constant term)
RegCoefs<-ResFS1$Coefficients
R2<-ResFS1$Rsquared
R2a<-ResFS1$RsquaredAdj
Foverall<-ResFS1$Foverall
Fp<-ResFS1$Fpvalue
yearsCal1<-c(yrV1[1],yrV1[length(yrV1)])
if (Fp<=0.05){
  Lsig<-TRUE
} else {
  Lsig<-FALSE
}
H<-summary(ResFS1$Model)
rmseCal<-H$sigma

################################################################################
#
# LEAVE-9-OUT CROSS-VALIDATION OF MODELS AT EACH STEP UP TO THE LAST STEP OF THE
# APPROXIMATE MAX-ADJR2 MODEL. 
# At each step compute REcv (reduction of error from cross-validation).
# Select as last step in model the step with maximum REcv. Argument i1 
# indicates cols of the lagged tree-ring index (1=lag-2, 3=no lag, 5=lag+2), in order
# of entry in the preliminary stepwise regression.
#
ResCV1<-CrossValid1(U1, V1, nNeg, nPos, i1)
if (ResCV1$REmaxStep<ni1){ 
  # stepwise cross-validation results in maximum RE at an earlier step than the
  # step of maximum adjusted R-squared of calibration. Will "cut" back, or simplify,
  # the model
  LREcut<-TRUE # 
  ni1<-ResCV1$REmaxStep # number of lags in model, after simplifying
  i1 <- i1[1:ni1] # which cols are in simplified model [12345 maps as lag-2 to lag+2]
  ModelPicked<-i1  # note that from L to R these are in order of entry in stepwise
  
  # Revise 02 Mar 2024
  # Want to reorder the column numbers i1  L to R from lowest to highest. This for
  # proper collating of regression coefficients with cols of lagged predictor matrix
  L <- is.element(a5,ModelPicked)
  ColsInModel <- a5[L] 
  # Revise end
  
} else {
  LREcut<-FALSE
}

# Build a 5-digit string indicating order of entry of predictors. Positions in the
# 5 indicate lag (1=lag-2, 3=no lag, 5=lag+2, etc). The number in the position is
# the step of entry into the regression. A "0" means does not enter. For example,
# [00100] means only lag-0 in model, and lag-0 of course first to enter. Another 
# example [20100] means lags -2 and 0 in model, with no-lag entering first and 
# lag -2 second
ModelCoded<-LagModel2Char(ModelPicked,5) # code the lags and order of entry into a string
REcv<-ResCV1$REcv
rmseCV<-ResCV1$RMSEcv
LeftOutCV<-ResCV1$LeftOut


################################################################################
#
# SPLIT-SAMPLE VALIDATION TO CHECK TEMPORAL STABILITY OF THE CURRENT MODEL AS
# IDENTIFIED BY MAX R2AD AND POSSIBLY SIMPLIFIED BY CROSS-VALIDATION

iAstop <- ceiling(dim(V1)[1]/2) # end row index in V1 of first half of data, assumed
# longer of the two halves if row-size of V1 odd
iBgo <- iAstop+1 # start row of second half
iA <- 1:iAstop # row indices of first half of full calib period
iB <- iBgo:(dim(V1)[1]) # ... of second half

#--- Calibrate on early, validate on late, then reverse
ical<-iA; ival<-iB
ResSS1=ssValid(V1,U1,ical,ival,i1);
REa1<-ResSS1$RE # RE for calib on early, valid on late
ical<-iB; ival<-iA
ResSS1=ssValid(V1,U1,ical,ival,i1);
REb1<-ResSS1$RE # RE for calib on late, valid on early


################################################################################
#
# DECIDE WHETHER YOU NEED TO REFIT AND RE-VALIDATE FINAL MODEL
# 
# Depending on how many lags are in the identified "best" model, and on the time coverage
# of the observed flows, it may be possible to extend the calibration period forward one or
# two years later than the last year for calibration in the fit requiring lag +2 on tree rings.
# Also, if cross-validation forced a simplification from the near-maximum-R2adj model, the
# predictors will be different from those in ResFS1$. Either situation indicates you should 
# refit. 
#
# Status: Y is a matrix with year in first column and the full available flows in the
# second column. Xorig, yrXorig are the full-length tree-ring chronologies; Xorig is 5-column and 
# has already been organized so that L-to-R are lags t-2 to t+2.
# Input yrsCalWindow is specified first and last year for period outside of which
# flows are not to be used for calibration; if NA for either first or last of this length-2 vector, 
# it is assumed to be OK to use  any of the flows in Y. The vector i1 indicates the columns of the
# lagged tree-ring matrix in the final model. The start and end years of the model fit in 
# the preliminary stepwise (allowing lags -2 to +2 to enter) is yearsCal1.

# Get flow data eligible for use in calibration

Y2<-Y
L<-complete.cases(Y2)
yrY2 <- Y2[L,1,drop=FALSE]
Y2<-Y[L,2,drop=FALSE]
yrgo2<-yrY2[1]; yrsp2<-yrY2[length(yrY2)] # initially consider all years usable
if (!is.na(yrsCalWindow[1])){
  yrgo2 <- max(yrsCalWindow[1],yrY2[1])
}
if (!is.na(yrsCalWindow[2])){
  yrsp2 <- min(yrsCalWindow[2],yrY2[length(yrY2)])
}
L<-yrY2>=yrgo2 & yrY2<=yrsp2
Y2<-Y2[L,,drop=FALSE]
yrY2<-yrY2[L,,drop=FALSE]

# Get tree-ring matrix with complete cases for columns i1

# revised 2024-03-04
Z<- Xorig[,ColsInModel,drop=FALSE] # changed from Z<- Xorig[,i1,drop=FALSE]
L <- complete.cases(Z)

Z <- Z[L,,drop=FALSE]
yrZ <- yrXorig[L,1,drop=FALSE]

#--- Find common period of Y1 and Z, which will define the "revised" calibration period
Zm <- as.matrix(cbind(yrZ,Z))
Y2m <- as.matrix(cbind(yrY2,Y2))
ResPC2 <- PeriodCommon(Zm,Y2m)
yearsCal2<-c(ResPC2$tgo, ResPC2$tsp) # revised calibration period

#--- Model will be refit if LREcut or if yearsCal2[2]>yrsCal[2]. One condition is that 
# cross-validation resulted in a simpler model that was initially fit (results stored in
# ResFS1$). The other condition is that fewer than two positive lags in the model may allow
# a later year for the calibration period than used in generating ResFS1$

if (LREcut || (yearsCal2[2]>yearsCal1[2])){
  Lrefit<-TRUE
} else {
  Lrefit<-FALSE
}

  
################################################################################
#
# IF INDICATED, RE-FIT AND RE-VALIDATE REGRESSION MODEL

if (!Lrefit){
  yearsCal<-yearsCal1 # if no refitting of model, final calib period is same as 
  # for preliminary stepwise regression
} else {
  # Store the predictor(s) and predictand for revised model
  U1<-ResPC2$X
  yrU1<-U1[,1,drop=FALSE]
  U1 <- U1[,-1,drop=FALSE] # drop year column
  colnames(U1)<-namesU[ColsInModel]; # U1, like X, organized L to R as lag-2 to Lag+2, as 1-5
  V1<-ResPC2$Y
  yrV1<-V1[,1,drop=FALSE]
  V1 <- V1[,-1,drop=FALSE]
  
  # Regression and storage of revised calibration statistics
  G <- lm(V1 ~ U1)
  H<-summary(G)
  R2<-H$r.squared
  R2a<-H$adj.r.squared
  Foverall<-H$fstatistic[1]
  Fp<-Fpvalue(G)
  if (Fp<=0.05){
    Lsig<-TRUE
  } else {
    Lsig<-FALSE
  }
  yearsCal<-yearsCal2
  RegCoefs<-G$coefficients
  rmseCal<-H$sigma
  
  # #debug on: These csv files for Matlab stepwise check on entry and coefficients after refit
  # # Commment out this block after debug check
  # write.csv(U1,file='U1r.csv')
  # write.csv(V1,file='V1r.csv')
  # #debug off
  
  #--- Cross-validation of the re-fit model
  ResCVrefit<- CrossValid2(U1,V1,nNeg,nPos)
  REcv<-ResCVrefit$REcv
  rmseCV<-ResCVrefit$RMSEcv
  LeftOutCV<-ResCVrefit$LeftOut
  
  #--- Split-sample validation of the re-fit model
  InModel<- 1:length(ModelPicked) # because for this refit the predictor matrix has
  # already been culled; all columns are in the model
  
  # Need pointer to rows of V1, first and second "halves"
  iAstop <- ceiling(dim(V1)[1]/2) # end row index in V1 of first half of data, assumed
  # longer of the two halves if row-size of V1 odd
  iBgo <- iAstop+1 # start row of second half
  iA <- 1:iAstop # row indices of first half of full calib period
  iB <- iBgo:(dim(V1)[1]) # ... of second half

  #--- Calibrate on early, validate on late, then reverse
  ical<-iA; ival<-iB
  ResSS2=ssValid(V1,U1,ical,ival,InModel);
  REa2<-ResSS2$RE # RE for calib on early, valid on late
  ical<-iB; ival<-iA
  ResSS2=ssValid(V1,U1,ical,ival,InModel);
  REb2<-ResSS2$RE # RE for calib on late, valid on early
}

# Set the correct REa and REb
if (Lrefit){
  REa<-REa2
  REb<-REb2
} else {
  REa<-REa1
  REb<-REb1
}


################################################################################
#
# DECIDE IF CHRONOLOGY SHOULD BE REJECTED FOR RECONSTRUCTION
#
# Reject if any one of these is true: 1) overall F of regression not significant
# at p<0.05, 2) Cross-validation reduction of error (REcv) not greater than zero,
# 3) reduction of error computed from split-sample validation on either half
# (REa,REb) not greater than zero. Depending on input Lcausal, could also reject 
# if only lags in model are negative on tree rings. Depending on RequireStable,
# could dispense with requirement for temporal stability, such that OK if REa<0 
# or REb<0 or if both are <0.
Lreject<-FALSE
if (RequireStable){
  L<- !Lsig || REcv <=0 || REa <= 0 || REb<=0
} else {
  L<- !Lsig || REcv <=0 
}
if (L) {Lreject=TRUE}

# Optional rejection of model with negative lags only on tree rings
LNegOnly=FALSE
if (!any(i1>=3)){
  LNegOnly <-TRUE
}
if (Lcausal && LNegOnly) Lreject<-TRUE


################################################################################
#
# APPLY FITTED REGRESSION MODEL TO RECONSTRUCT. 

# Get sub-matrix of full-length lagged chronologies with complete cases for the
# lags in model
# Xr <- as.matrix(Xorig[,ModelPicked]) # code error corrected by next line on 2 Mar 2024
Xr <- as.matrix(Xorig[,ColsInModel]) 
L<-complete.cases(Xr)
Xr <- as.matrix(Xr[L,])
yrXr<-matrix(yrXorig[L,1])
mXr <- dim(Xr)[1]

# Add ones column and reconstruct
Xones<-matrix(1,nrow=mXr,ncol=1)
Xr <- cbind(Xones,Xr)
yh <- Xr %*% RegCoefs # reconstruction as 1-col matrix
yh <- cbind(yrXr,yh)
yearsRec<-c(yh[1,1],yh[mXr,1])


###### BUILD 5-CHARACTER CODE OF SIGN OF REGRESSION COEFFICIENT ON LAGS -2 TO +2
#
# Code like [00P0P], meaning positive coefficients on lags 0 and t+1, and other lags not in model
#ModelSign <- LagModel2Sign(i1,5,RegCoefs) ---- 
ModelSign <- LagModel2Sign(ColsInModel,5,RegCoefs) # revised 2024-03-04


################################################################################
#  
# MAKE NAMED LIST FOR RETURN 

Output<-list("ColsInModel"=ColsInModel,"Model"=ModelPicked,"ModelCoded"=ModelCoded,
  "ModelSign"=ModelSign,"yearsCal"=yearsCal,"yearsRec"=yearsRec,"MaxLagNegPos"=c(nNeg,nPos),
  "LeftOutCV"=LeftOutCV,"IncrementR2adj"=c1,"RegCoefs"=RegCoefs,"Rsquared"=R2,
  "RsquaredAdj"=R2a,"F"=Foverall,"Fp"=Fp,  "Lsig"=Lsig,"REcv"=REcv,"REa"=REa,"REb"=REb,
  "RMSE"=c(rmseCal,rmseCV),"LREcut"=LREcut,"Lrefit"=Lrefit,"LNegOnly"=LNegOnly,
  "Lreject"=Lreject,"yhat"=yh)
return(Output)
}

ForwMoot <- function(X,namesX,y) {
  # Moot-point forward stepwise regression for special case of one potential predictor.
  # D. Meko
  # Last revised 2022-12-13
  #
  # Tailored for reconsw4, such that works with middle column (lag 0) of a  matrix
  # of predictors assumed to be lagged t-1:t+2 from predictand
  # 
  # INPUT ARGUMENTS
  # y [matrix] 1-col  of predictand
  # X [matrix] one-col matrix of the single potential predictor
  # namesX [character] vector of id of potential predictor
  #
  # OUTPUT
  # H: named list, 
  # names(H)<-c('Model','StepMaxR2adj','ColsInModel','Coefficients',
  # 'ColsInOrderEntry','Rsquared','RsquaredAdj','Step',
  # 'RsquaredAllSteps','RsquaredAdjAllSteps','Foverall',
  # 'Fpvalue')
  #    Model [lm object] is a special type of R object, and holds many of the
  #     regression statistics and data, such as the residuals and predicted data.
  #     R has functions that operate on lm objects. For example,
  #     summary(H#Model) shows a summary table of the regression
  #   The Coefficients, combined with ColsInModel would allow a reconstruction to
  #     be generated from the long time series of potenial predictors X. 
  #     ColsInModel gives the columns of that matrix that the coefficients apply 
  #      to. By plotting RsquaredAllSteps or RsquaredAdjAllSteps agains step, you
  #     can disply how R-square and adjusted R-square changes with step in the 
  #     stepwise modelin. 
  #   Most of the other list items are obvious from their names. StepMaxR2adj is
  #     the step at which adjusted R-square reaches a maximum. 
  
  source(paste(code_dir,"Fpvalue.R",sep="")) # p-value of overall-F from lm() [not written by Meko]
  
  # Regression
  G <- lm(y ~ X[,3]) # regression on col 3 of X, which is no-lag from y
  kstep <- 1 # only one possible step in simple linear regression
  
  Foverall<-summary(G)$fstatistic[1]
  p<-Fpvalue(G)
    
  H  <- list('Model'=G,'StepMaxR2adj'=kstep,'ColsInModel'=3,'Coefficients'=G$coefficients,
  'ColsInOrderEntry'=3,'Rsquared'=summary(G)$r.squared,'RsquaredAdj'=summary(G)$adj.r.squared,
  'Step'=kstep,'RsquaredAllSteps'=summary(G)$r.squared,
  'RsquaredAdjAllSteps'=summary(G)$adj.r.squared,'Foverall'=Foverall,'Fpvalue'=p)
  
  return(H)

}

