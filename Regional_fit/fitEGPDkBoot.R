library(plyr)#rbind.fill.matrix
library(mev)
library(evd)
library(evd)#pgpd
library(zipfR)#incomplete beta function
library(foreach);library(iterators);library(parallel);library(doParallel) # parallelization

#simulation of data below
fitEGPDk.boot <- function(M, sites = "default", ncores, cens_thres=c(1,Inf), round=0.1,
                          thres=0, precision=.0001,
                          ParInit=c(0.5,0.5,0.2), loop.max = 20){
  #   DESCRIPTION
  # Fit regional, semi regional and local versions of EGPD in Naveau et al. (2016).
  #   ARGUMENTS
  # M, matrix of positive precipitation (1 col = 1 station)
  # Each row corresponds to a day, each col to a site
  # ncores, number of cores to use in the parallelization
  # censor_thres, bounds to estimate parameters
  # thres = threshold of wet days (e.g. 2mm) if not null, provides conditional parameters
  
  #   VALUE
  # List of parameters (ThetaReg, ThetaSReg ,Theta_0)
  #   - the fist element of the list is the list of regional flexibility parameter kappa,
  #     semi-regional scale parameters sigma and regional shape parameter xi.
  #   - the second element of the list is the list of semi regional flexibility parameter kappa,
  #     and scale parameters sigma,  regional shape parameter xi.
  #   - the third element of the list contains initialization of parameters (i.e. at-site estimates)
  ###############################################
  #if only one temporal serie put in a matrix
  
  registerDoParallel(cores=ncores)
  
  if(is.null(dim(M))){
    M = as.matrix(M,ncol=1)}
  nday = nrow(M)
  nstat = ncol(M)
  if(is.character(sites)){if(sites=="default"){sites = 1:nstat}else{stop("undefined option for choice of sites")} }
  nb_stat = length(sites)
  
  SigmaS_names =  SigmaR_names = KappaS_names = rep(NA,nstat)
  CensoredMean = moyN  = rep(NA,nb_stat)
  #consider only days with data
  obs_rr = matrix(nrow=0,ncol=0)
  Theta = list()
  Theta_init = list()
  kappa_init = sigma_init = xisite_init= rep(NA,nstat)
  ## BUILD NAMES, REMOVE NA IN PRECIP MATRIX
  for (station in 1:nstat){
    SigmaR_names[station] = paste("sigmaR",station,sep="_")
    SigmaS_names[station] = paste("sigmaS",station,sep="_")
    KappaS_names[station] = paste("kappaS",station,sep="_")
    
    y = M[,station]
    y = y[!is.na(y)]
    y = y[y>thres]-thres#consider only positive precipitation when thres =0
    moy = mean(y,na.rm=TRUE)
    obs_rr = rbind.fill.matrix(obs_rr,t(as.matrix(y)))#and fill the matrix with NA if needed
  }#end loop on station for data formatting
  M = t(obs_rr)
  
  # INITIALIZATION PARAMETERS (provide at-site parameters)
  Theta_list <- foreach(station=1:nstat) %dopar% {
    library(mev)
    y =na.omit(M[,station])
    fit_init = fit.extgp(y,model=1,method="pwm",init = ParInit, censoring = cens_thres, rounded = round,plots = FALSE)
    kappa_init = fit_init$fit$pwm[1]
    sigma_init = fit_init$fit$pwm[2]
    xisite_init = fit_init$fit$pwm[3]
    return(c(kappa_init, sigma_init, xisite_init))
  } #end foreach station 
  Theta_mat0 <- matrix(unlist(Theta_list), ncol = length(Theta_list))
  
  kappa_init <- Theta_mat0[1,]
  sigma_init <- Theta_mat0[2,]
  xisite_init <- Theta_mat0[3,]
  #REGIONAL PARAMETERS
  xi_reg = mean(xisite_init)
  shape = list("xi.reg" = xi_reg, "xi.site" = xisite_init)
  kappa_reg = mean(kappa_init)
  
  # INIT OF THE LOOP ON CHOSEN SITES
  Theta_init = list("kappa"=kappa_init[sites],"sigma"=sigma_init[sites],"xi.reg"=rep(xi_reg,nb_stat))
  ThetaS = ThetaR = Theta_init
  
  loop <-0; increment <-precision +.1 
  log.init<-0
  Xi.old <- unique(ThetaS$xi.reg)
  u = cens_thres[1]
  
  
  while((loop<loop.max)&(increment>precision)){
    loop <- loop+1
    #########################
    ## SEMI REGIONAL FLEX AND SCALE PARAMETERS
    Sigma.old <- ThetaS$sigma
    Kappa.old <- ThetaS$kappa
    
    CensoredMean_list <- foreach(station=1:nb_stat) %dopar% {
      y = na.omit(M[,station])
      return(mean(y[y>u],na.rm=TRUE))
    } #end foreach station 
    
    CensoredMean <- unlist(CensoredMean_list)
    Sigma.new <- (Xi.old*CensoredMean)/(Kappa.old*IB(H(u,Sigma.old,Xi.old),
                                                     1,Kappa.old,1-Xi.old)/Fbar(u,Kappa.old,
                                                                                Sigma.old,Xi.old)-1)
    ThetaS$sigma <- Sigma.new
    
    moyN_list <- foreach(station=1:nb_stat) %dopar% {
      site = sites[station]
      y = na.omit(M[,site])
      return(mean(y[y/Sigma.new[station]>u/Sigma.new[station]]/Sigma.new[station],na.rm=TRUE))
    } #end foreach station
    
    moyN <- unlist(moyN_list) # at site censored mean (data normalized by sigma.new)
    
    
    Kappa.new = (Xi.old*moyN)/(1*IB(H(u,Sigma.new,Xi.old),
                                    1,Kappa.old,1-Xi.old)/Fbar(u,Kappa.old,
                                                               Sigma.new,Xi.old)-1/Kappa.old)
    ThetaS$kappa = Kappa.new
    #update kappa and sigma in ThetaS
    
    #########################
    ##  REGIONAL FLEX AND SCALE PARAMETERS
    
    SigmaR.old <- ThetaR$sigma
    
    SigmaR.new <- (xi_reg*CensoredMean)/(kappa_reg*IB(H(u,SigmaR.old,xi_reg),
                                                      1,kappa_reg,1-xi_reg)/Fbar(u,kappa_reg,
                                                                                 SigmaR.old,xi_reg)-1)
    
    
    ThetaR$sigma <- SigmaR.new
    
    
    
    
    
    
    increment<-max(abs(SigmaR.old-ThetaR$sigma),
                   abs(Kappa.old-ThetaS$kappa),
                   abs(Sigma.old-ThetaS$sigma))
    
  }#end while
  
  kappaS = ThetaS$kappa; sigmaS = ThetaS$sigma; sigmaR = ThetaR$sigma
  return(list("shape"=shape, 
              "scale"= list("sigma.local"=sigma_init, "sigma.semireg"=sigmaS,"sigma.reg"=sigmaR),
              "flexibility" = list("kappa.local"= kappa_init, "kappa.semireg" = kappaS, "kappa.reg" = kappa_reg)))
}#end function fitEGPDk.boot




simulEGPD<-function( kappa=1/10, sigma, xi,n=100){
  # Arguments
  #   xi = shape parameter with length= nb regions
  #   n = sample size
  #   
  # Valeur
  #  rEGPD matrix with n rows (a row = a wet day), and number of sites (a col= an EGPD(kappa,sigma,kappa) sample)
  
  m<-length(xi)  # nombre regions
  if(length(sigma)!=m) stop("xi and sigma with different dimensions")
  out<-matrix(NA,n,m)
  for(i in 1:m){
    if(xi[i]>0){out[,i]<-(sigma[i]/xi[i])*((1-(runif(n))^(1/kappa[i]))^(-xi[i])-1)}#temp serie with length n
    if(xi[i]==0){out[,i] <- -sigma[i]*log(1-(runif(n))^(1/kappa[i]))}
  }
  return(out)
}
xi = c(.1,.2)
sigma = c(1,16)
kappa = c(1,4)
M = simulEGPD(kappa = kappa, sigma=sigma,xi=xi,n=200)




























