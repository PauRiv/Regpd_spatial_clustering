##################### FUNCTIONS FOR REGIONAL INFERENCE ########################
library(plyr)#rbind.fill.matrix
library(mev)#egp2.fit, qegd,...
library(evd)#qgpd
###############################################
######### FIT EGPD(ki,s,xi)

fitEGPDkSemiRegCensored <- function(M,method="pwm",cens_thres=c(1,Inf),ParReg = "mean", thres=0){
  #   DESCRIPTION
  # Fit a regional version of EGPD in Naveau et al. (2016).
  #   ARGUMENTS
  # M, matrix of positive precipitation 
  # Each row corresponds to a day, each col to a site
  # censor_thres, bounds to estimate parameters
  # ParReg, method for regionalizing shape parameter xi : 
  #         either "mean" (i.e. mean of the at-site estimates)
  #         or "pool" (pool normalized data to estimate regional xi)
  # thres = threshold of wet days (e.g. 2mm) if not null, provides conditional parameters
  
  #   VALUE
  # List of parameters: 
  #   - the fist element of the list is the list of semi-regional scale parameters
  #   - the second and last element of the list contains regional shape parameter kappa and the regional shape parameter xi estimated on pooled normalized data
  ###############################################
  #if only one temporal serie put in a matrix
  
  if(is.null(dim(M))){
    M = as.matrix(M,ncol=1)}
  nday = nrow(M)
  nstat = ncol(M)
  MN = matrix(NA,nrow = nday,ncol = nstat)
  as_mean = as_Lcensor_thres = rep(NA,nstat)
  as_Rcensor_thres = rep(Inf,nstat)
  #consider only days with data
  obs_rr = matrix(nrow=0,ncol=0)
  for (station in 1:nstat){
    y = M[,station]
    y = y[!is.na(y)]
    y = y[y>thres]-thres#consider only positive precipitation when thres =0
    as_mean[station] = mean(y,na.rm=TRUE)
    obs_rr = rbind.fill.matrix(obs_rr,t(as.matrix(y)))#and fill the matrix with NA if needed
  }
  M = t(obs_rr)
  fitreg = as.data.frame(matrix(NA,nrow = 4,ncol=nstat),
                         row.names = c("kappa","sigma","xi.site","xi.reg"))
  for (station in 1:nstat) {
    moy = as_mean[station]
    y = na.omit(M[,station])
    as_Lcensor_thres[station] = cens_thres[1]/moy
    as_Rcensor_thres[station] = cens_thres[2]/moy
    #y = y/mean(y)#test to fit on normalized data
    InitPar = c(.9,.5,.1)
    fitreg[1:3,station] = fit.extgp(y[y>0],method = method,init=InitPar,
                                    censoring = cens_thres,
                                    plots = FALSE)$fit$pwm
    MN[,station] = y/moy
  }
  yn = as.vector(MN)
  Xi.reg.pool = fit.extgp(yn[yn>0],method=method,init = InitPar,
                          censoring = c(max(as_Lcensor_thres),min(as_Rcensor_thres)),plots = FALSE)$fit$pwm[3]
  Xi.reg.mean = mean(data.matrix(fitreg[3,]))
  if(ParReg=="pool"){xi.reg = Xi.reg.pool}
  if(ParReg=="mean"){xi.reg = Xi.reg.mean}
  Kappa.site = as.vector(data.matrix(fitreg[1,]))
  SigmaSemiReg = as_mean*xi.reg/(Kappa.site*beta(Kappa.site,1-xi.reg)-1)
  fitreg[2,] = SigmaSemiReg
  fitreg[4,] = rep(xi.reg,nstat)
  return(fitreg)
}




































