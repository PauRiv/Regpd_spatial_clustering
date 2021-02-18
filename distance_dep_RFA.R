### CONVEX COMBINATION OF DISTANCES (BOTH RFA AND DEPENDENCE)
library(Kendall)
source("RFA_function.R")#RFA clustering, see Le Gall et al., 2020
######
D_lambda <- function(lambda, x,y){
  ## COMPUTE A CONVEX COMBINATION OF KENDALL CORRELATION AND PW RATIO DIFFERENCE
  # ARGUMENTS
  # x, y, time series
  # lambda, in [0;1] weight on the RFA part
  # VALUE
  # lambda*|omega(x)-omega(y)| + (1-lambda)*tau(x,y)
  
  #pre-processing: x and y with same length
  if(lambda>1|lambda<0)stop("lambda must be between 0 and 1")
  x = x[which(!is.na(x*y))]; y = y[which(!is.na(x*y))]
  x = na.omit(x);y = na.omit(y)
  
  d_omega = abs(xi.Ratio(x, independence = FALSE)
                -xi.Ratio(y,independence = FALSE))
  KT = as.double(Kendall(x,y)$tau[1])
  dist = lambda*d_omega + (1-lambda)*(1-KT)#or/2 ?
  return(dist)
}


























