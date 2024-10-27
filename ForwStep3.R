ForwStep3 <- function(X,namesX, y,kstop,nNeg,nPos,ic) {
  # Forward stepwise regression with validation/calibration stopping rules
  # D. Meko
  # Last revised 2024-03-09
  #
  # Initial forward stepwise entry is done until all variables are in. The final model
  # is selected optionally as the model for step at which adjusted R-square is "approximately"
  # maximum, or at that or a lower step at which cross-validation RE is maximum. 
  #
  # INPUT ARGUMENTS
  # y [matrix] 1-col  of predictand
  # X [matrix] one-or-more cols of pool of potential predictors
  # namesX [character] vector of ids of potential predictors (cols of X)
  # kstop [numeric] stopping criterion: 1= approximate maximimum adjusted R-squared,
  #   2=maximum cross-validation RE
  # nNeg [numberic] maximum negative lag on chronologies allowed in modeling
  # nPos [numeric] maximum positive...CrossValidStorage
  # ic [numeric] critical increment in adjusted R-squared. If kstop=1, entry stops when next step would 
  #   yield increase of adjusted R-squared less than ic. If kstop=2, model could stop at an even earlier
  #   step if cross-validation RE reaches a maximum before the step of approximate maximum 
  #   adjusted R-squared.
  #
  # OUTPUT
  # H: named list, with elements: 
  c('Model','StopCriterion','StoppingStep','ColsInModel','Coefficients',
    'ColsInOrderEntry','Rsquared','RsquaredAdj','Step',
    'RsquaredAllSteps','RsquaredAdjAllSteps','Foverall',
    'Fpvalue',
    'CrossValidStorage','NleftOutinCV','REcv','RMSEcv',
    'ModelTable')
  #    Model [lm object] special type of R object that holds many of the
  #     regression statistics and data, such as the residuals and predicted data.
  #     R has functions that operate on lm objects. For example,
  #     summary(H$Model) shows a summary table of the regression
  #   StopCriterion: criterion for picking best model (1=max adjusted R-squares; 
  #      2=maximum cross-validation RE)
  #   StoppingStep: the step at which stepwise regression stopped according to StopCriterion
  #   ColsInModel: which columns of input X are in the final model; the order is important,
  #      and matches the order of coefficients in H$coefficients (after the intercept)
  #   Coefficients: the regression coefficients, starting with intercept term and then
  #      in order as in H$ColsInModel.The Coefficients, combined with ColsInModel allow 
  #      a reconstruction when applied to matrix X. 
  #   ColsInOrderEntry: order in which the columns of X entered the stepwise model. In general,
  #      this differs from order in H$ColsInModel
  #   Rsquared,RsquaredAdj: R-squared and adjusted R-squared of the final model
  #   Step, RsquaredAllSteps,RsquaredAdjAllSteps: step is a numeric vector (1,2, ...) of 
  #      the step in the stepwise process; the other two variables are the R-squared and
  #      adjusted R-squared at each step
  #   Foverall,Fpvalue: overall-F of the final regression model and p-value of that F
  #   CrossValidStorage: a list of detailed cross-validation statistics returned by
  #      CrossValid1(). No need to dig into this unless curious. See the comments for
  #      Crossvalid1() for definition of list.
  #   NleftOutinCV,REcv,RMSEcv: cross-validation statistics. NleftOutCV is the number of 
  #      observations left out at each iteration of cross-validation. With lags, this is
  #      greater than 1, and set to assure no tree-ring data providing a cross-validation
  #      prediction are also used to calibrate the prediction model. REcv is the 
  #      reduction-of-error statistic. RMSEcv is the root-mean-square error of cross-validation.
  #   ModelTable: a table object listing the same coefficients as in H$coefficients, but also
  #      including the corresponding id of the predictor, cross-referenced to columns of
  #      input matrix X and indicating if lagged positively or negatively and by how much.
  #      For example, X7 is the seventh column of X, no lags; and X2P2 is the second column of
  #      of X lagged +2 years forward from the predictand.
  #
  # revised 2024-03-09: 1) fixed code error that under some circumstances could
  #   result too few stepwise steps being allow based on maximum REcv. Fix may
  #   give a stronger (higher R-squared) model, with more predictors; 
  #   2)cosmetic. Comments cleaned up for typos and clarity.

  source(paste(code_dir,"CrossValid1.R",sep="")) # leave-m-out cross-validation  
  
  #--- ALLOCATE AND INITIALIZE 
  Np<-dim(X)[2]  # size predictor pool
  i1<-1:Np # index to predictors in pool
  i2<-rep(NA,Np) # index in order of entry
  Lin=rep(FALSE,Np)  # initial for in model
  R2<-rep(NA,Np)
  R2a<-rep(NA,Np)
  
  
  #--- FOREWARD STEPWISE REGRESSION
  
  for (j in 1:Np){
    if (j==1){
      # First pass pick X correlated highest (absolute R) with y
      colnames(X)<-namesX # re-initialize these 
      r<-cor(y,X)
      H<-sort(abs(r),decreasing=TRUE,index.return=TRUE)
      iwinner<-H$ix[1]
      i2[j]<-iwinner
      Lin[iwinner]<-TRUE # Logical to cols of X in model to be estimated
      U<-X[,Lin,drop="FALSE"]  # cull matrix of those predictors
      G<-lm(y~U)  # estimate model (G is a "lm" object )
      R2this<-summary(G)$r.squared
      R2adjthis<-summary(G)$adj.r.squared
      R2[j]<-R2this  # store R-squared for this model
      R2a[j]<-R2adjthis # store adjusted R-squared ...
      e<-as.matrix(G$residuals) # model residuals; 
    } else {
      # iterations after the first. e are the residuals of the previous step.
      # At each step, e is correlated with the cols of X not yet in model.
      # "Chosen one" is the the col of X with highest correlation
      colnames(X)<-namesX # re-initialize these 
      ibullpen <-i1[!Lin] # cols of X NOT yet in models
      r<-cor(e,X[,ibullpen])
      H<-sort(abs(r),decreasing=TRUE,index.return=TRUE)
      iwinner<-H$ix[1] # pointer to highest correlated member of sub-matrix
      i2[j]<-ibullpen[iwinner] # corresponding column of that member in X
      ithis<-ibullpen[iwinner] #... renamed for simplicity
      Lin[ithis]<-TRUE # logical of cols of X in current model to be estimated
      # Next statements exactly same as those described for first iteration
      U<-X[,Lin,drop="FALSE"]
      G<-lm(y~U) # estimation
      R2this<-summary(G)$r.squared
      R2adjthis<-summary(G)$adj.r.squared
      R2[j]<-R2this
      R2a[j]<-R2adjthis
      e<-as.matrix(G$residuals)
    }
  }
  
  #--- FIND "BEST" MODEL AND RE-FIT IT
  rm(H)
  Lin<-rep(FALSE,Np) # re-initialize as no variables in model
  
  if (kstop==1){
    # MAXIMUM ADJUSTED R-SQUARED VERSION
    # Commented out code for using maximum adjusted R2 as criterion for best
    # model. For large number of variables in pool of potential predictos, and
    # with those possibly chosen by correlation screening from larger number
    # of possible predictors, this criterion tends to over-fit model. Decided
    # to stop entry instead FIRST maximimum of adjusted R2
    #
    s<-sort(R2a,decreasing=TRUE,index.return=TRUE) # sorted adjusted R-squared,
    #  from highest to lowest.
    iwinner<-s$ix[1] # step at which adjusted R-squared highest
  }else{
    # MAXIMUM RE VERSION
    CV <-CrossValid1(X, y, nNeg, nPos, i2) 
    iwinner <- CV$REmaxStep # revised 2024-03-09
  }
  
  i2a<-i2[1:iwinner] # cols of X in model to be fit
  Lin[i2a]<-TRUE # turn on logical for cols of X in model to be fit

  #---- PULL SUBMATRIX OF PREDICTORS AND ESTIMATE MODEL
  U<-X[,Lin,drop="FALSE"]
  colnames(X)<-namesX # re-initialize these 
  colnames(U)<-colnames(X)[Lin]
  G<-lm(y~U) # estimation
  
  # Fix names of coefficients so not prefixed with "U"
  # Cosmetic, so that when user types "G$coefficients" at terminal he
  # is not bothered by that "U"
  bb <- G$coefficients # Regression coefficients
  bbnms  <- names(bb)
  bbnms1 <- bbnms[-1] # names of terms, without "(Intercept)"
  bbnms2 <- substring(bbnms1,2)
  names(G$coefficients) <- c(bbnms[1],bbnms2)
  rm(bb,bbnms,bbnms1,bbnms2)
  

  #--- STATISTICS OF FINAL FITTED MODEL
  
  kstep<-sum(Lin) # final step (max adj R-square)
  
  if (kstop==1){
    stopHow<-'Approximate maximum of adsjusted R-squared'
  }else{
    stopHow<-'Maximum cross-validation RE'
  }
  
  H = list() 
  H[[1]]<-G
  H[[2]]<-stopHow
  
  H[[3]]<-kstep
  H[[4]]<-i1[Lin]
  H[[5]]<-G$coefficients
  H[[6]]<-i2[1:kstep]
  H[[7]]<-summary(G)$r.squared
  H[[8]]<-summary(G)$adj.r.squared
  jstep<-(1:Np)
  H[[9]]<-jstep
  H[[10]]<-R2
  H[[11]]<-R2a
  H[[12]]<-summary(G)$fstatistic[1]
  
  #  function for p-value of F; from
  # https://stackoverflow.com/questions/5587676/pull-out-p-values-and-r-squared-from-a-linear-regression
  lmp <- function (modelobject) {
    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(modelobject)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
  }
  p<-lmp(G)
  H[[13]]<-p # pvalue of overall F
  
  if (kstop==2){
    H[[14]] <- CV
    H[[15]] <- CV$LeftOut
    H[[16]] <- CV$REcv
    H[[17]] <- CV$RMSEcv
  }else{
    H[[14]] <- NA
    H[[15]] <- NA
    H[[16]] <- NA
    H[[17]] <- NA
  }
  # Build table with regression coeffients opposite variable names
  tab<-matrix(t(G$coefficients))
  colnames(tab)<-c('Model')
  rownames(tab)<-c('Intercept',namesX[H[[4]]])
  tab<- as.table(tab)
  H[[18]]<-tab
  
  names(H)<-c('Model','StopCriterion','StoppingStep','ColsInModel','Coefficients',
              'ColsInOrderEntry','Rsquared','RsquaredAdj','Step',
              'RsquaredAllSteps','RsquaredAdjAllSteps','Foverall',
              'Fpvalue',
              'CrossValidStorage','NleftOutinCV','REcv','RMSEcv',
              'ModelTable')
  return(H)
}




