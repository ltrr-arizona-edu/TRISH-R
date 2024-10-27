ForwStep2 <- function(X,namesX, y, c1) {
  # Forward stepwise regression
  # D. Meko
  # Last revised 2024-03-04
  #
  # The forward stepwise entry is done until all variables are in. The final model
  # is selected as the model for step at which adjusted R-square is maximum
  #
  # INPUT ARGUMENTS
  # y [matrix] 1-col  of predictand
  # X [matrix] one-or-more cols of pool of potential predictors
  # namesX [character] vector of ids of potential predictors (what's in cols of X)
  # c1 [1x1]r required incremental increase in adjusted R-square to warrant additional 
  #   step in forward stepwise
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
  #
  # revised 2024-03-04: corrected error in call to ForwStep2 for stepwise variables entering
  #   after the first to to enter. At "r<-cor(e,X[,ibullpen])". Error could have resulted in
  #   MLR regression models picking next variable (not in model) most highlly correlated
  #   with predictand, y, rather than with residuals, e, from previous step. Change could
  #   result in model with stronger statistics (e.g., R-squared)
  
  
  
  source(paste(code_dir,"Fpvalue.R",sep="")) # p-value of overall-F from lm() [not written by Meko]
  
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
      # First, pass pick X correlated highest (absolute R) with y
      colnames(X)<-namesX # re-initialize these 
      r<-cor(y,X)
      H<-sort(abs(r),decreasing=TRUE,index.return=TRUE)
      iwinner<-H$ix[1]
      i2[j]<-iwinner
      Lin[iwinner]<-TRUE # Logical to cols of X in model to be estimated
      U<-X[,Lin]  # culll matrix of those predictors
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
      r<-cor(e,X[,ibullpen]) # rev 2024-03-04
      H<-sort(abs(r),decreasing=TRUE,index.return=TRUE)
      iwinner<-H$ix[1] # pointer to highes correlated member of sub-matrix
      i2[j]<-ibullpen[iwinner] # corresponding column of that member in X
      ithis<-ibullpen[iwinner] #... renamed for simplicity
      Lin[ithis]<-TRUE # logical of cols of X in current model to be estimated
      # Next statements exactly same as those described for first iteration
       U<-X[,Lin]
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
  
  # MAXIMUM ADJUSTED R-SQUARED VERSION
  # Commented out code for using maximum adjusted R2 as criterion for best
  # model. For large number of variables in pool of potential predictos, and
  # with those possibly chosen by correlation screening from larger number
  # of possible predictors, this criterion tends to over-fit model. Decided
  # to stop entry instead FIRST maximimum of adjusted R2
  #
  # s<-sort(R2a,decreasing=TRUE,index.return=TRUE) # sorted adjusted R-squared, dec.
  # iwinner<-s$ix[1] # step at which adjusted R-squared highest 
  # i2a<-i2[1:iwinner] # cols of X in model to be fit
  # Lin[i2a]<-TRUE # turn on logical for cols of X in model to be fit

  # FIRST LOCAL MAXIMUM OF ADJUSTED R-SQUARED
  
  if (length(Lin)==1){
    # only one predictor in pool; only one step
    iwinner<-1
  } else {
    d<-diff(R2a) # first difference of adjusted R-squared with step in model
    if (all(d>=c1)){
      # Adjusted R-squared always increasing by at least c1 with step
      iwinner<-length(Lin)
    } else {
      # adjusted R-squared increases by less than c1 at some step (step ithis+1)
      ithis<-which(d < c1) # find step of first "insufficient" increase  in adjusted R-squared
      iwinner<-ithis[1]
    }
  }
  i2a<-i2[1:iwinner] # cols of X in final model 
  Lin[i2a]<-TRUE # turn on logical for cols of X in model to be fit
  
  #---- PULL SUBMATRIX OF PREDICTORS AND ESTIMATE MODEL
  # if adjR2 begins dropping from 1st step, only 1 series enters matrix...
  if(length(i2a)==1){
    U <- matrix(X[,Lin])
  }
  else{
    U <- X[,Lin]
  }
  colnames(X)<-namesX # re-initialize these 
  colnames(U)<-colnames(X)[Lin]
  G<-lm(y~U) # estimation
  
  
  #--- STATISTICS OF FINAL FITTED MODEL
  
  kstep<-sum(Lin) # final step (max adj R-square)
  
  H = list() 
  H[[1]]<-G
  H[[2]]<-kstep
  H[[3]]<-i1[Lin]
  H[[4]]<-G$coefficients
  H[[5]]<-i2[1:kstep]
  H[[6]]<-summary(G)$r.squared
  H[[7]]<-summary(G)$adj.r.squared
  jstep<-(1:Np)
  H[[8]]<-jstep
  H[[9]]<-R2
  H[[10]]<-R2a
  H[[11]]<-summary(G)$fstatistic[1]
  p<-Fpvalue(G)
  H[[12]]<-p # pvalue of overall F
  names(H)<-c('Model','StepMaxR2adj','ColsInModel','Coefficients',
              'ColsInOrderEntry','Rsquared','RsquaredAdj','Step',
              'RsquaredAllSteps','RsquaredAdjAllSteps','Foverall',
              'Fpvalue')
  return(H)
}




