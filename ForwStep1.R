ForwStep1 <- function(X,namesX, y,kstop,nNeg,nPos) {
  # Forward stepwise regression
  # D. Meko
  # Last revised 2023-08-16
  #
  # The forward stepwise entry is done until all variables are in. The final model
  # is selected optionally as the model for step at which adjusted R-square is maximum
  # or the cross-validation RE is a maximum
  #
  # INPUT ARGUMENTS
  # y [matrix] 1-col  of predictand
  # X [matrix] one-or-more cols of pool of potential predictors
  # namesX [character] vector of ids of potential predictors (what's in cols of X)
  # kstop [numeric] stopping criterion: 1=maximimum adjusted R-squared,
  #   2=maximum cross-validation RE
  # nNeg [numeric] maximum negative lag on chronologies allowed in modeling
  # nPos [numeric] maximum positive...
  #
  # OUTPUT
  # H: named list, with elements: 
  # c('Model','StopCriterion','StoppingStep','ColsInModel','Coefficients',
  #   'ColsInOrderEntry','Rsquared','RsquaredAdj','Step',
  #   'RsquaredAllSteps','RsquaredAdjAllSteps','Foverall',
  #   'Fpvalue',
  #   'CrossValidStorage','NleftOutinCV','REcv','RMSEcv',
  #   'ModelTable')
  #    Model [lm object] special type of R object that holds many of the
  #     regression statistics and data, such as the residuals and predicted data.
  #     R has functions that operate on lm objects. For example,
  #     summary(H$Model) shows a summary table of the regression
  #   StopCriterion: criterion for picking best model (1=max adjusted R-squared; 
  #      2=maximum cross-validation RE)
  #   StoppingStep: the step at which stepwise regression stopped according to StopCriterion
  #   ColsInModel: which columns of input X are in the final model; the order is important,
  #      and matches the order of coefficients in H$coefficients (after the intercept)
  #   Coefficients: the regression coefficients, starting with intercept term and then
  #      in order as in H$ColsInModel.The coefficients, combined with ColsInModel allow 
  #      a reconstruction to be generated from matrix X. 
  #   ColsInOrderEntry: order in which the columns of X entered the stepwise model. In general,
  #      this differs from order in H$ColsInModel
  #   Rsquared,RsquaredAdj: R-squared and adjusted R-squared of the final model
  #   Step, RsquaredAllSteps,RsquaredAdjAllSteps: step is a numeric vector (1,2, ...) of 
  #      the step in the stepwise process; the other two variables are the R-squared and
  #      adjusted R-squared at each step
  #   Foverall,Fpvalue: overall-F of the final regression model and p-value of that F
  #   CrossValidStorage: a list of detailed cross-validation statistics returned by
  #      CrossValid1(). See opening comments of function Crossvalid1 details
  #   NleftOutinCV,REcv,RMSEcv: cross-validation statistics. NleftOutCV is the number of 
  #      observations left out at each iteration of cross-validation. With lags, this is
  #      greater than 1, and set to ensure that no tree-ring data are used both to provide
  #      a cross-validation prediction in a specific and to calibrate the model giving 
  #      that prediction. REcv is the reduction-of-error statistic. RMSEcv is the 
  #      root-mean-square error.
  #   ModelTable: a table object listing the same coefficients as in H$coefficients, but also
  #      including the corresponding id of the predictor, cross-referenced to columns of
  #      input matrix X and indicating if lagged positively or negatively and by how much.
  #      For example, X7 is the seventh column of X, no lags; and X2P2 is the second column of
  #      of X lagged +2 years forward from the predictand.
  #
  # Rev 2023-08-16: cosmetic, to make sure that summary(G) retains  correct col 
  #    name when just on predictor, on call to lm()
  # Rev 2024-12-21: so that global "code_dir" used to access CrossValid1.R
  
  source(paste(code_dir,"CrossValid1.R",sep="")) 

  #--- ALLOCATE AND INITIALIZE 
  Np<-dim(X)[2]  # size predictor pool
  i1<-1:Np # index to predictors in pool
  i2<-rep(NA,Np) # index in order of entry
  Lin=rep(FALSE,Np)  # logical pointer to variables in model
  R2<-rep(NA,Np) # Regression R-squared at each step
  R2a<-rep(NA,Np) # Adjusted...
  
  
  #--- FORWARD STEPWISE REGRESSION
  
  for (j in 1:Np){
    if (j==1){
      
      # First, pass pick X correlated highest (absolute correlation) with y
      colnames(X)<-namesX # re-initialize these 
      r<-cor(y,X)
      H<-sort(abs(r),decreasing=TRUE,index.return=TRUE)
      iwinner<-H$ix[1]
      i2[j]<-iwinner
      Lin[iwinner]<-TRUE # turn on logical for variable just entered

      U<-X[,Lin,drop="FALSE"]  # cull matrix of those predictors
      # If only 1 predictor, need next to make sure lm receive col name
      if (dim(U)[2]==1){
        colnames(U) <- colnames(X)[Lin]
      }
      G<-lm(y~U)  # estimate model (G is a "lm" object )
      R2this<-summary(G)$r.squared
      R2adjthis<-summary(G)$adj.r.squared
      R2[j]<-R2this  # store R-squared for this model
      R2a[j]<-R2adjthis # store adjusted R-squared ...
      e<-as.matrix(G$residuals) # model residuals; 
    } else {
      # iteration after the first. e are the residuals of the previous step.
      # At each step, e is correlated with the cols of X not yet in model.
      # "Chosen one" is the the col of X with highest correlations
      colnames(X)<-namesX # re-initialize these 
      ibullpen <-i1[!Lin] # cols of X NOT yet in models
      r<-cor(e,X[,ibullpen])
      H<-sort(abs(r),decreasing=TRUE,index.return=TRUE)
      iwinner<-H$ix[1] # pointer to highes correlated member of sub-matrix
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
  
  #--- LEAVE-M-OUT CROSS-VALIDATION, STEPWISE
  #
  # Do regardless of stopping rule. CV will be a list with cross-validation results
  CV <-CrossValid1(X, y, nNeg, nPos, i2) 
  
  #--- FIND "BEST" MODEL AND RE-FIT IT
  rm(H)
  Lin<-rep(FALSE,Np) # re-initialize as no variables in model
  
  if (kstop==1){
    # MAXIMUM ADJUSTED R-SQUARED VERSION
    # Commented out code for using maximum adjusted R2 as criterion for best
    # model. For large number of variables in pool of potential predictors, and
    # with those possibly chosen by correlation screening from larger number
    # of possible predictors, this criterion tends to over-fit model. Decided
    # to stop entry instead FIRST maximum of adjusted R2
    #
    s<-sort(R2a,decreasing=TRUE,index.return=TRUE) # sorted adjusted R-squared,
    #  from highest to lowest.
    iwinner<-s$ix[1] # step at which adjusted R-squared highest
  }else{
    # MAXIMUM RE VERSION
    iwinner <- CV$REmaxStep # step at which RE highest
  }
  i2a<-i2[1:iwinner] # cols of X in model to be fit
  Lin[i2a]<-TRUE # turn on logical for cols of X in model to be fit

  
  #---- PULL SUBMATRIX OF PREDICTORS AND RE-ESTIMATE MODEL
  U<-X[,Lin,drop="FALSE"]
  colnames(X)<-namesX # re-initialize these 
  colnames(U)<-colnames(X)[Lin]
  G<-lm(y~U) # estimation
  
  # Do not want variables to be listed U1, U2, etc. Want real IDs
  bbnms <- names(G$coefficients)
  bbnms1 <- bbnms[-1] # names of terms, without "(Intercept)"
  bbnms2 <- colnames(X)[Lin]
  names(G$coefficients) <- c(bbnms[1],bbnms2)
  rm(bb,bbnms,bbnms1,bbnms2)
  
  #--- STATISTICS OF FINAL FITTED MODEL
  
  kstep<-sum(Lin) # final step (max adj R-square)
  
  if (kstop==1){
    stopHow<-'Max R-squared adjusted'
  }else{
    stopHow<-'Max RE'
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
  H[[14]] <- CV
  H[[15]] <- CV$LeftOut 
  if (kstop==2){
    H[[16]] <- CV$REcv
    H[[17]] <- CV$RMSEcv
  }else{
    # CrossValid1 returns CVREv and RMSEcv as the "best" according to maximum RE, but
    # if using maximum adjusted R-squared, want instead to use index iwinner to grab those
    # applicable statistics.
    H[[16]] <- CV$REcvAll[iwinner]
    H[[17]] <- CV$RMSEvall[iwinner]
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




