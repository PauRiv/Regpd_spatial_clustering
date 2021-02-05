######### FIT EGPD(k,s,xi)
library(mev)
library(evd)#pgpd
library(plyr)#rbind.fill.matrix
library(zipfR)#incomplete beta function
############################
## PRELIMINARY FUNCTIONS
############################
IB <- function(x,y,a,b){
  #incomplete beta as defined in appendix of Naveau et al., 2016
  z = Ibeta(y,a,b)-Ibeta(x,a,b)
  return(z)
}
H <- function(x,sigma,xi){
  z = pgpd(x,scale=sigma,shape=xi)
  return(z)
}
pG = function(kappa,u){
  fct=u^kappa# = G(u)
  return(fct)
}
Fbar <- function(x,kappa,sigma,xi){
  z = 1-pG(kappa,H(x,sigma,xi))
  return(z)
}

fitEGPDkSemiRegCensoredIter <- function(M,method="pwm",cens_thres=c(0,Inf),round=0.1,
                                        ParReg = "mean", thres=0, precision=.0001,
                                        ParInit=c(0.5,0.5,0.2), loop.max = 10){
  #   DESCRIPTION
  # Fit a SEMI regional version of EGPD in Naveau et al. (2016).
  #   ARGUMENTS
  # M, matrix of positive precipitation 
  # Each row corresponds to a day, each col to a site
  # censor_thres, bounds to estimate parameters
  # ParReg, method for regionalizing shape parameter xi : 
  #         either "mean" (i.e. mean of the at-site estimates)
  #         or "pool" (pool normalized data to estimate regional xi)
  # thres = threshold of wet days (e.g. 2mm) if not null, provides conditional parameters
  
  #   VALUE
  # List of parameters (Theta,Theta_0)
  #   - the fist element of the list is the list of semi-regional flexibility parameter kappa,
  #     semi-regional scale parameters sigma and regional shape parameter xi.
  #   - the second element of the list contains initialization of parameters (i.e. at-site estimates)
  ###############################################
  #if only one temporal serie put in a matrix
  
  if(is.null(dim(M))){
    M = as.matrix(M,ncol=1)}
  nday = nrow(M)
  nstat = ncol(M)
  
  Sigma_names=Kappa_names=rep(NA,nstat)
  CensoredMean = moyN  = rep(NA,nstat)
  #consider only days with data
  obs_rr = matrix(nrow=0,ncol=0)
  obs_norm = matrix(0,nrow = 0,ncol =0)
  Theta = list()
  Theta_init = list()
  kappa_init = sigma_init = xisite_init= rep(NA,nstat)
  for (station in 1:nstat){
    Sigma_names[station] = paste("sigma",station,sep="_")
    Kappa_names[station] = paste("kappa",station,sep="_")
    
    y = M[,station]
    y = y[!is.na(y)]
    y = y[y>thres]-thres#consider only positive precipitation when thres =0
    moy = mean(y,na.rm=TRUE)
    obs_rr = rbind.fill.matrix(obs_rr,t(as.matrix(y)))#and fill the matrix with NA if needed
    obs_norm = rbind.fill.matrix(obs_norm,t(as.matrix(y/moy)))#and fill the matrix with NA if needed
  }#end loop on station for data formatting
  M = t(obs_rr)
  MN=t(obs_norm)
  #definir Theta et sa forme (matrice ? liste ? dataframe ?)
  
  for (station in 1:nstat) {
    y =na.omit(M[,station])
    fit_init = fit.extgp(y,model=1,method=method,init =ParInit, censoring = cens_thres, rounded = round,plots = FALSE)
    kappa_init[station] = fit_init$fit$pwm[1]
    sigma_init[station] = fit_init$fit$pwm[2]
    xisite_init[station] = fit_init$fit$pwm[3]
  }#end loop on station for initialization
  xi_init = mean(xisite_init)
  Theta_init = list("kappa"=kappa_init,"sigma"=sigma_init,"xi.reg"=rep(xi_init,nstat))
  
  Theta = Theta_init
  Theta_0 = list("kappa"=kappa_init,"sigma"=sigma_init,"xi.site"=xisite_init,"xi.reg"=rep(xi_init,nstat))
  loop <-0; increment <-precision +.1 
  log.init<-0
  Xi.old <- unique(Theta$xi.reg)
  u = cens_thres[1]
  while((loop<loop.max)&(increment>precision)){
    loop <- loop+1
    Sigma.old <- Theta$sigma
    Kappa.old <- Theta$kappa
    
    
    for (i in 1:nstat){ y = na.omit(M[,i])
    CensoredMean[i] <- mean(y[y>u],na.rm=TRUE)
    }#at-site censored mean 
    
    Sigma.new <- (Xi.old*CensoredMean)/(Kappa.old*IB(H(u,Sigma.old,Xi.old),
                                                     1,Kappa.old,1-Xi.old)/Fbar(u,Kappa.old,
                                                                                Sigma.old,Xi.old)-1)
    Theta$sigma <- Sigma.new
    for (i in 1:nstat){
      y = na.omit(M[,i])
      moyN[i] = mean(y[y/Sigma.new[i]>u/Sigma.new[i]]/Sigma.new[i],na.rm=TRUE)}# at site censored mean (data normalized by sigma.new)
    
    Kappa.new = (Xi.old*moyN)/(1*IB(H(u,Sigma.new,Xi.old),
                                    1,Kappa.old,1-Xi.old)/Fbar(u,Kappa.old,
                                                               Sigma.new,Xi.old)-1/Kappa.old)
    Theta$kappa = Kappa.new
    #update kappa and sigma in Theta
    
    increment<-max(abs(Sigma.old-Theta$sigma),abs(Kappa.old-Theta$kappa))
    
  }#end while
  
  
  return(list("Theta"=Theta,"Theta_0"=Theta_0))
}


fitEGPDkREGCensorIter <- function(M,method="pwm",cens_thres=c(0,Inf),round=0.1,
                                  ParReg = "mean", thres=0, precision=.0001,
                                  ParInit=c(0.5,0.5,0.2), loop.max = 20){
  #   DESCRIPTION
  # Fit a REGIONAL  version of EGPD in Naveau et al. (2016).
  #   ARGUMENTS
  # M, matrix of positive precipitation 
  # Each row corresponds to a day, each col to a site
  # censor_thres, bounds to estimate parameters
  # ParReg, method for regionalizing shape parameter xi : 
  #         either "mean" (i.e. mean of the at-site estimates)
  #         or "pool" (pool normalized data to estimate regional xi)
  # thres = threshold of wet days (e.g. 2mm) if not null, provides conditional parameters
  
  #   VALUE
  # List of parameters: (Theta,Theta_0)
  #   - the fist element of the list is the list of regional kappa, semi-regional scale parameters sigma_i, and regional shape parameter xi
  #   - the second element of the list contains initialization parameters (at-site estimates)
  ###############################################
  #if only one temporal serie put in a matrix
  
  if(is.null(dim(M))){
    M = as.matrix(M,ncol=1)}
  nday = nrow(M)
  nstat = ncol(M)
  
  Sigma_names=Kappa_names=rep(NA,nstat)
  CensoredMean = moyN  = rep(NA,nstat)
  #consider only days with data
  obs_rr = matrix(nrow=0,ncol=0)
  obs_norm = matrix(0,nrow = 0,ncol =0)
  Theta = list()
  Theta_init = list()
  kappasite_init = sigma_init = xisite_init= rep(NA,nstat)
  for (station in 1:nstat){
    Sigma_names[station] = paste("sigma",station,sep="_")
    
    y = M[,station]
    y = y[!is.na(y)]
    y = y[y>thres]-thres#consider only positive precipitation when thres =0
    moy = mean(y,na.rm=TRUE)
    obs_rr = rbind.fill.matrix(obs_rr,t(as.matrix(y)))#and fill the matrix with NA if needed
    obs_norm = rbind.fill.matrix(obs_norm,t(as.matrix(y/moy)))#and fill the matrix with NA if needed
  }#end loop on station for data formatting
  M = t(obs_rr)
  MN=t(obs_norm)
  #definir Theta et sa forme (matrice ? liste ? dataframe ?)
  
  for (station in 1:nstat) {
    y =na.omit(M[,station])
    fit_init = fit.extgp(y,model=1,method=method,init =ParInit, censoring = cens_thres, rounded = round,plots = FALSE)
    kappasite_init[station] = fit_init$fit$pwm[1]
    sigma_init[station] = fit_init$fit$pwm[2]
    xisite_init[station] = fit_init$fit$pwm[3]
  }#end loop on station for initialization
  Theta_loc=list("kappai" = kappasite_init,"sigmai"=sigma_init,"xii"=xisite_init)
  
  kappa_init = mean(kappasite_init)
  xi_init = mean(xisite_init)
  Theta_init = list("kappa.reg"=rep(kappa_init,nstat),"sigma"=sigma_init,"xi.reg"=rep(xi_init,nstat))
  
  Theta = Theta_init
  Theta_0 = list("kappa.reg"=rep(kappa_init,nstat),"sigma"=sigma_init,"xi.reg"=rep(xi_init,nstat),
                 "kappa.site"=kappasite_init,"xi.site"=xisite_init)
  loop <-0; increment <-precision +.1 
  log.init<-0
  Xi.old <- unique(Theta$xi.reg)
  Kappa.old = unique(Theta$kappa.reg)
  u = cens_thres[1]
  while((loop<loop.max)&(increment>precision)){
    loop <- loop+1
    Sigma.old <- Theta$sigma
    
    
    for (i in 1:nstat){ y = na.omit(M[,i])
    CensoredMean[i] <- mean(y[y>u],na.rm=TRUE)
    }#at-site censored mean 
    
    Sigma.new <- (Xi.old*CensoredMean)/(Kappa.old*IB(H(u,Sigma.old,Xi.old),
                                                     1,Kappa.old,1-Xi.old)/Fbar(u,Kappa.old,
                                                                                Sigma.old,Xi.old)-1)
    Theta$sigma <- Sigma.new
    
    increment<-max(abs(Sigma.old-Theta$sigma))
    
  }#end while
  
  
  return(list("Theta"=Theta,"Theta_loc"=Theta_loc))
}


































