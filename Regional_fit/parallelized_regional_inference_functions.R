######### FIT EGPD(k,s,xi)
library(mev)
library(evd)#pgpd
library(plyr)#rbind.fill.matrix
library(zipfR)#incomplete beta function
library(foreach);library(iterators);library(parallel);library(doParallel) # parallelization
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

fitEGPDkSemiRegCensoredIter_paral <- function(M, ncores, method="pwm", cens_thres=c(0,Inf),round=0.1,
                                        ParReg = "mean", thres=0, precision=.0001,
                                        ParInit=c(0.5,0.5,0.2), loop.max = 20){
  #   DESCRIPTION
  # Fit a (semi)regional version of EGPD in Naveau et al. (2016).
  #   ARGUMENTS
  # M, matrix of positive precipitation (1 col = 1 station)
  # ncores, number of cores to use in the parallelization
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
  
  registerDoParallel(cores=ncores)
  
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
    # obs_norm = rbind.fill.matrix(obs_norm,t(as.matrix(y/moy)))#and fill the matrix with NA if needed
  }#end loop on station for data formatting
  M = t(obs_rr)
  #definir Theta et sa forme (matrice ? liste ? dataframe ?)
  
  Theta_list <- foreach(station=1:nstat) %dopar% {
    library(mev)
    y =na.omit(M[,station])
    fit_init = fit.extgp(y,model=1,method=method,init =ParInit, censoring = cens_thres, rounded = round,plots = FALSE)
    kappa_init = fit_init$fit$pwm[1]
    sigma_init = fit_init$fit$pwm[2]
    xisite_init = fit_init$fit$pwm[3]
    
    return(c(kappa_init, sigma_init, xisite_init))
  } #end foreach station 
  
  Theta_mat0 <- matrix(unlist(Theta_list), ncol = length(Theta_list))
  
  kappa_init <- Theta_mat0[1,]
  sigma_init <- Theta_mat0[2,]
  xisite_init <- Theta_mat0[3,]
  
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
    
    CensoredMean_list <- foreach(station=1:nstat) %dopar% {
      y = na.omit(M[,station])
      return(mean(y[y>u],na.rm=TRUE))
    } #end foreach station 
    
    CensoredMean <- unlist(CensoredMean_list)
    
    Sigma.new <- (Xi.old*CensoredMean)/(Kappa.old*IB(H(u,Sigma.old,Xi.old),
                                                     1,Kappa.old,1-Xi.old)/Fbar(u,Kappa.old,
                                                                                Sigma.old,Xi.old)-1)
    Theta$sigma <- Sigma.new
    
    moyN_list <- foreach(station=1:nstat) %dopar% {
      y = na.omit(M[,station])
      return(mean(y[y/Sigma.new[station]>u/Sigma.new[station]]/Sigma.new[station],na.rm=TRUE))
    } #end foreach station
    
    moyN <- unlist(moyN_list) # at site censored mean (data normalized by sigma.new)
    
    
    Kappa.new = (Xi.old*moyN)/(1*IB(H(u,Sigma.new,Xi.old),
                                    1,Kappa.old,1-Xi.old)/Fbar(u,Kappa.old,
                                                               Sigma.new,Xi.old)-1/Kappa.old)
    Theta$kappa = Kappa.new
    #update kappa and sigma in Theta
    
    increment<-max(abs(Sigma.old-Theta$sigma),abs(Kappa.old-Theta$kappa))
    
  }#end while
  
  
  return(list("Theta"=Theta,"Theta_0"=Theta_0))
}# end function fitEGPDkSemiRegCensoredIter_paral




fitEGPDkSemiRegCensoredIter.boot_paral <- function(M, ncores, method="pwm", cens_thres=c(0,Inf),round=0.1,
                                                   sites = "default", thres=0, precision=.0001,
                                              ParInit=c(0.5,0.5,0.2), loop.max = 20){
  #   DESCRIPTION
  # Fit a (semi)regional version of EGPD in Naveau et al. (2016).
  #   ARGUMENTS
  # M, matrix of positive precipitation (1 col = 1 station)
  # ncores, number of cores to use in the parallelization
  # Each row corresponds to a day, each col to a site
  # censor_thres, bounds to estimate parameters
  # sites, indices of GP/stations where local parameters are computed (e.g. medoid, min/max silhouette)
  # thres = threshold of wet days (e.g. 2mm) if not null, provides conditional parameters
  
  #   VALUE
  # List of parameters (Theta,Theta_0)
  #   - the fist element of the list is the list of semi-regional flexibility parameter kappa,
  #     semi-regional scale parameters sigma and regional shape parameter xi.
  #   - the second element of the list contains initialization of parameters (i.e. at-site estimates)
  ###############################################
  #if only one temporal serie put in a matrix
  
  registerDoParallel(cores=ncores)
  
  if(is.null(dim(M))){
    M = as.matrix(M,ncol=1)}
  nday = nrow(M)
  nstat = ncol(M)
  if(is.character(sites)){if(sites=="default"){sites = 1:nstat}else{stop("undefined option for choice of sites")} }
  nb_stat = length(sites)
  
  Sigma_names=Kappa_names=rep(NA,nstat)
  CensoredMean = moyN  = rep(NA,nb_stat)
  #consider only days with data
  obs_rr = matrix(nrow=0,ncol=0)
  # obs_norm = matrix(0,nrow = 0,ncol =0)
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
    # obs_norm = rbind.fill.matrix(obs_norm,t(as.matrix(y/moy)))#and fill the matrix with NA if needed
  }#end loop on station for data formatting
  M = t(obs_rr)

  Theta_list <- foreach(station=1:nstat) %dopar% {
    library(mev)
    y =na.omit(M[,station])
    fit_init = fit.extgp(y,model=1,method=method,init =ParInit, censoring = cens_thres, rounded = round,plots = FALSE)
    kappa_init = fit_init$fit$pwm[1]
    sigma_init = fit_init$fit$pwm[2]
    xisite_init = fit_init$fit$pwm[3]
    
    return(c(kappa_init, sigma_init, xisite_init))
  } #end foreach station 
  
  Theta_mat0 <- matrix(unlist(Theta_list), ncol = length(Theta_list))
  
  kappa_init <- Theta_mat0[1,]
  sigma_init <- Theta_mat0[2,]
  xisite_init <- Theta_mat0[3,]
  
  xi_init = mean(xisite_init)
  Theta_init = list("kappa"=kappa_init[sites],"sigma"=sigma_init[sites],"xi.reg"=rep(xi_init,nb_stat))
  
  Theta = Theta_init
  Theta_0 = list("kappa"=kappa_init[sites],"sigma"=sigma_init[sites],"xi.site"=xisite_init[sites],"xi.reg"=rep(xi_init,nb_stat))
  loop <-0; increment <-precision +.1 
  log.init<-0
  Xi.old <- unique(Theta$xi.reg)
  u = cens_thres[1]
  while((loop<loop.max)&(increment>precision)){
    loop <- loop+1
    Sigma.old <- Theta$sigma
    Kappa.old <- Theta$kappa
    
    CensoredMean_list <- foreach(station=1:nb_stat) %dopar% {
      y = na.omit(M[,station])
      return(mean(y[y>u],na.rm=TRUE))
    } #end foreach station 
    
    CensoredMean <- unlist(CensoredMean_list)
    
    Sigma.new <- (Xi.old*CensoredMean)/(Kappa.old*IB(H(u,Sigma.old,Xi.old),
                                                     1,Kappa.old,1-Xi.old)/Fbar(u,Kappa.old,
                                                                                Sigma.old,Xi.old)-1)
    Theta$sigma <- Sigma.new
    
    moyN_list <- foreach(station=1:nb_stat) %dopar% {
      site = sites[i]
      y = na.omit(M[,site])
      return(mean(y[y/Sigma.new[station]>u/Sigma.new[station]]/Sigma.new[station],na.rm=TRUE))
    } #end foreach station
    
    moyN <- unlist(moyN_list) # at site censored mean (data normalized by sigma.new)
    
    
    Kappa.new = (Xi.old*moyN)/(1*IB(H(u,Sigma.new,Xi.old),
                                    1,Kappa.old,1-Xi.old)/Fbar(u,Kappa.old,
                                                               Sigma.new,Xi.old)-1/Kappa.old)
    Theta$kappa = Kappa.new
    #update kappa and sigma in Theta
    
    increment<-max(abs(Sigma.old-Theta$sigma),abs(Kappa.old-Theta$kappa))
    
  }#end while
  
  
  return(list("Theta"=Theta,"Theta_0"=Theta_0))
}# end function fitEGPDkSemiRegCensoredIter.boot_paral



fitEGPDkREGCensoredIter.boot_paral <- function(M, ncores, method="pwm", cens_thres=c(0,Inf),round=0.1,
                                                   sites = "default", thres=0, precision=.0001,
                                                   ParInit=c(0.5,0.5,0.2), loop.max = 20){
  #   DESCRIPTION
  # Fit a regional version of EGPD in Naveau et al. (2016).
  #   ARGUMENTS
  # M, matrix of positive precipitation (1 col = 1 station)
  # ncores, number of cores to use in the parallelization
  # Each row corresponds to a day, each col to a site
  # censor_thres, bounds to estimate parameters
  # sites, indices of GP/stations where local parameters are computed (e.g. medoid, min/max silhouette)
  # thres = threshold of wet days (e.g. 2mm) if not null, provides conditional parameters
  
  #   VALUE
  # List of parameters (Theta,Theta_0)
  #   - the fist element of the list is the regional flexibility parameter kappa,
  #     semi-regional scale parameters sigma and regional shape parameter xi.
  #   - the second element of the list contains initialization of parameters (i.e. at-site estimates)
  ###############################################
  #if only one temporal serie put in a matrix
  
  registerDoParallel(cores=ncores)
  
  if(is.null(dim(M))){
    M = as.matrix(M,ncol=1)}
  nday = nrow(M)
  nstat = ncol(M)
  if(is.character(sites)){if(sites=="default"){sites = 1:nstat}else{stop("undefined option for choice of sites")} }
  nb_stat = length(sites)
  
  Sigma_names=Kappa_names=rep(NA,nstat)
  CensoredMean = moyN  = rep(NA,nb_stat)
  #consider only days with data
  obs_rr = matrix(nrow=0,ncol=0)
  # obs_norm = matrix(0,nrow = 0,ncol =0)
  Theta = list()
  Theta_init = list()
  kappa_init = sigma_init = xisite_init= rep(NA,nstat)
  for (station in 1:nstat){
    Sigma_names[station] = paste("sigma",station,sep="_")
    Kappa_names[station] = paste("kappa",station,sep="_")
    
    y = M[,station]
    y = y[!is.na(y)]
    y = y[y>thres]-thres#consider only positive precipitation when thres =0
    # moy = mean(y,na.rm=TRUE)
    obs_rr = rbind.fill.matrix(obs_rr,t(as.matrix(y)))#and fill the matrix with NA if needed
    # obs_norm = rbind.fill.matrix(obs_norm,t(as.matrix(y/moy)))#and fill the matrix with NA if needed
  }#end loop on station for data formatting
  M = t(obs_rr)
  
  Theta_list <- foreach(station=1:nstat) %dopar% {
    library(mev)
    y =na.omit(M[,station])
    fit_init = fit.extgp(y,model=1,method=method,init =ParInit, censoring = cens_thres, rounded = round,plots = FALSE)
    kappa_init = fit_init$fit$pwm[1]
    sigma_init = fit_init$fit$pwm[2]
    xisite_init = fit_init$fit$pwm[3]
    
    return(c(kappa_init, sigma_init, xisite_init))
  } #end foreach station 
  
  Theta_mat0 <- matrix(unlist(Theta_list), ncol = length(Theta_list))
  
  kappasite_init <- Theta_mat0[1,]
  sigma_init <- Theta_mat0[2,]
  xisite_init <- Theta_mat0[3,]
  
  kappa_init = mean(kappasite_init)#regional version of kappa
  xi_init = mean(xisite_init)
  Theta_init = list("kappa.reg"=rep(kappa_init,nb_stat),"sigma"=sigma_init[sites],"xi.reg"=rep(xi_init,nb_stat))
  
  Theta = Theta_init
  Theta_0 = list("kappa"=kappa_init[sites],"sigma"=sigma_init[sites],"xi.site"=xisite_init[sites],
                 "kappa.reg" = rep(kappa_init,nb_stat),"xi.reg"=rep(xi_init,nb_stat))
  loop <-0; increment <-precision +.1 
  log.init<-0
  Xi.old <- unique(Theta$xi.reg)
  Kappa.old = unique(Theta$kappa.reg)
  u = cens_thres[1]
  while((loop<loop.max)&(increment>precision)){
    loop <- loop+1
    Sigma.old <- Theta$sigma
    
    CensoredMean_list <- foreach(station=1:nb_stat) %dopar% {
      y = na.omit(M[,station])
      return(mean(y[y>u],na.rm=TRUE))
    } #end foreach station 
    
    CensoredMean <- unlist(CensoredMean_list)
    
    Sigma.new <- (Xi.old*CensoredMean)/(Kappa.old*IB(H(u,Sigma.old,Xi.old),
                                                     1,Kappa.old,1-Xi.old)/Fbar(u,Kappa.old,
                                                                                Sigma.old,Xi.old)-1)
    Theta$sigma <- Sigma.new

    #update kappa and sigma in Theta
    
    increment<-max(abs(Sigma.old-Theta$sigma),abs(Kappa.old-Theta$kappa))
    
  }#end while
  
  
  return(list("Theta"=Theta,"Theta_0"=Theta_0))
}# end function fitEGPDkSemiRegCensoredIter.boot_paral










































