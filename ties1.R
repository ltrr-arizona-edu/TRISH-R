ties1 <- function(x){
# Identify ties in a time series
# D Meko
# Last revised 2022-05-10
#
#--- Input
#   x [numeric] time series
#--- Output is a list N with fields
#   ngroups = number of groups of ties
#   nties number of members each group
#--- Notes
#   Modeled after my matlab ties1.m function
#   Needed in context of Mann-Kendall trend analysis (mannken1.R)
#   Make sure that input x is a vector, and has no NAs
#   $nties: if a signle tie in a series, would make the nties(j)=2, since two
#     data values are involved in a single tie

#--- Check input
L1 <- is.vector(x)
L2 <- !any(is.na(x))
L <- L1 & L2
if (!L){
  stop('Input x must be a vector with no NAs')
}


mx <- length(x)

j1<-1:mx   # numbering vector same length as x

#--- Check that not all values in x are unique; if all are unique, no need for further
# work, and can return empty elements in list Output
c <- unique(x) # unique value in x; a 1-col matrix, even if x a vector
if (length(c)==mx) {
  Output <- list(ngroups=vector(),nties=vector())
  return(Output)
}

#--- find the dupicates
L <- duplicated(x)
d  <- x[L] # the duplicates, a vector
ndupe=length(d)
  
#--- find the unique duplicates, which will represent "groups" of ties
g <- unique(d)  
ngrp <- length(g) # number of groups

#--- Number of members for each group. Will dupe vectors to matrics. 

# Convert vector of unique non-uniques in x into a row vector and row-dupe it to row
# size equal to the number of duplicates, ndupe rows
V <- matrix(g, nrow=ndupe, ncol=length(g), byrow=TRUE)

#  col-dupe the vector of  values, d,  involved in any ties  
F <- matrix(d, nrow=length(d), ncol=ngrp, byrow=FALSE)

# associate each tie-value with a group; number of trues in each col are number
# of tie-values in each group; 
L <- V==F 

# sum over rows and add 1 to get number of obs in x associated with ties in the ngrp groups
nties <- colSums(L)+1  # a vector

#---- BUILD LIST
Output <-list(ngroups=ngrp,nties=nties)
}




