##################### FUNCTIONS FOR REGIONAL INFERENCE ########################
library(plyr)#rbind.fill.matrix
library(mev)#egp2.fit, qegd,...
library(evd)#qgpd
###############################################
# function:  weights in Algo 1   from a sample U_1,...,U_n
###############################################
pBB = function(w,u){ 
  # ARGUMENTS 
  # w are Bernstein weights
  # m is the degree of Bernstein Polynomial
  m_w <- length(w)
  beta.weights <-matrix(NA,m_w, length(u))
  w.beta <- beta.weights
  
  for (i in 1:m_w){
    beta.weights[i,] <- pbeta(u,shape1 = i, shape2 = m_w-i+1)
    w.beta[i,] <- as.numeric(w[i]) * beta.weights[i,]
  }
  fct<-apply(w.beta, MARGIN = 2, FUN = "sum")
  return(fct)
}
weights <- function(m,x,xi, sigma){ #, kap=1){
  u <- evd:::pgpd(x, scale=sigma, shape=xi)
  w<- numeric(m);
  FnU <- ecdf(u)
  ii<-c(1:m)/m;iii<-c(0:(m-1))/m
  w<-FnU(ii) -FnU(iii)
  if(w[m]==0){
    # cat("\n  1-G_{m,m}(1-1/m)=",1-pBB(w,1-1/m))
    w[m] = 1-pBB(w,1-1/m);
    w = w/sum(w) 
  }
  # if(w[m]==0){
  #   w[m] = kap/m;
  #   w = w/sum(w) 
  #   cat("\t  kap/m=", round(kap/m,3), "\t w[m]=", round(w[m],3))
  # }
  return(w)
}
pwmGP<-function(x,xi.positive=T){
  # pwmGP
  # Description
  # 
  #   Inference of xi and sigma for a GP sample with PWM inference
  #   with the assumption that 0 \geq xi \gep 1 
  #   (if estimate xi is negative then it is forced to be equal to zero)
  #
  # Arguments
  # x vector of non-negative values 
  #
  # Value
  # Estimates of c("sigma","xi")  
  ###############################################
  x<-x[!is.na(x)]
  Fn <-ecdf(x)
  m0<-mean(x)
  m1<-mean(x*(1-Fn(x)))
  xi.pwm <- (m0-4*m1)/(m0-2*m1)
  if(xi.pwm<0){
    if(xi.positive) xi.pwm <-0
    # cat("\n xi set to zero because its  PWM estimate is negative")
  }
  sigma.pwm <- m0*(1-xi.pwm)
  out<-c(sigma.pwm,xi.pwm)
  names(out)<-c("sigma","xi")
  return(out)
}
#####################################
## Regional and at-site estimation of GPD parameters by PWM
pwmGPreg<-function(M,xi.positive=T){
  ## ARGUMENTS
  # M = matrix with temporal series in col, nrow=ndays and ncol= number of stations
  # series come from sites in a same homogeneous region, see Le Gall et al., 2020
  ## VALUE
  # par_reg = vector of length nstat+1 containing semi regional estimates of
  # sigma for the nstat first coordinates and regional xi for last coordinate
  if(is.null(dim(M))){
    #if M is a temporal serie, put in a matrix
    M = as.matrix(M,ncol=1)}
  nday = nrow(M)
  nstat = ncol(M)
  obs_rr = matrix(nrow=0,ncol=0)
  pwm0 = pwm1 = nb_rr =scale_names= rep(NA,nstat)
  
  for (i in 1:nstat){
    x = M[,i]
    x = x[!is.na(x)]
    #x = x[x > 1]#consider only positive precipitation
    obs_rr = rbind.fill.matrix(obs_rr,t(as.matrix(x)))#and fill the matrix with NA if needed
    nb_rr[i] = length(x)
    
    Fni <- ecdf(x)#at site cdf
    pwm0[i] = mean(x)#at-site pwm 0
    pwm1[i] = mean(x*(1-Fni(x)))#at-site pwm1
    scale_names[i] = paste("sigma",i,sep="_")
  }
  #regional estimation of xi
  x = as.vector(M)
  x<-x[!is.na(x)]
  Fn <-ecdf(x)
  m0reg<-mean(x)
  m1reg<-mean(x*(1-Fn(x)))
  xi.reg <- (m0reg-4*m1reg)/(m0reg-2*m1reg)
  if(xi.reg<0){
    if(xi.positive) xi.reg <- 0
    # cat("\n xi set to zero because its  PWM estimate is negative")
  }
  sigma.pwm <- pwm0*(1-xi.reg)
  out<-c(sigma.pwm,xi.reg)
  names(out)<-c(scale_names,"xi.reg")
  return(out)
}
############################
#### FIT EGPD(G,sigma,xi) REGIONALLY
#M is a positive matrix (transformation x[x>2]-2 done)
fitGPreg<-function(M,m=10, precision=.0001,loop.max=50, verba=F){
  #   DESCRIPTION
  # Fit the regional version of Bernstein mixture GPD in Tencaliec et al. (2019).
  #   ARGUMENTS
  # M, matrix of positive precipitation (matrix.obs[matrix.obs>detection_threshold]-detection_threshold)
  # Each row corresponds to a day, each col to a site
  # m, degree of Bernstein Polynomial
  # precision   the optimisation scheme stops when the parameters change is lower than precision
  #   or when the number of iteration is larger than 
  # loop.max 
  #   VALUE
  # Theta list of parameters: 
  #   - the fist element of the list is the list of Bernstein weights for each site
  #   - the second and last element of the list contains semi-regional scale (sigma_k) parameters and the regional shape (xi) parameter (mean of at-site shape parameters)
  ###############################################
  
  #if only one temporal serie put in a matrix
  if(is.null(dim(M))){
    M = as.matrix(M,ncol=1)}
  nday = nrow(M)
  nstat = ncol(M)
  
  #consider only days with data
  obs_rr = matrix(nrow=0,ncol=0)
  for (i in 1:nstat){
    y = M[,i]
    y = y[!is.na(y)]
    #y = y[y > 1]#consider only positive precipitation
    obs_rr = rbind.fill.matrix(obs_rr,t(as.matrix(y)))#and fill the matrix with NA if needed
  }
  M = t(obs_rr)
  M.new <- M
  Theta = rep(NA,nstat+1)
  Thetai = matrix(data = NA,nrow = 2,ncol = nstat)#matrix of the couple of at-site estimates (sigma_i, xi_i)
  Xi.site = rep(NA,nstat)
  Scale_names = rep(NA,nstat)
  Wi = list()#list of weight estimates for each site
  #initialization of EGPD parameters
  for (i in 1:nstat) {
    Para = egp2.fit(na.omit(M.new[,i]), model=1,method="pwm", init=c(1,1,.2),plots=F)$fit$pwm
    Theta[i]= Para[2]
    Xi.site[i] = Para[3]
    Scale_names[i] = paste("sigma", i, sep = "_")
  }
  Theta[nstat+1] = mean(Xi.site)
  names(Theta) = c(Scale_names,"xi.reg")
  loop <-0; increment <-precision +.1 
  log.init<-0
  #Patricia's loop with regional or semi-regional estimation of xi and sigma
  while((loop<loop.max)&(increment>precision)){
    loop <- loop+1
    Xi.old <- Theta["xi.reg"]
    for (i in 1:nstat) {
      Sigma <- paste("sigma",i,sep="_")
      Sigma.old <- Theta[Sigma]
      x <- na.omit(M.new[,i])
      z <- evd:::pgpd(x, scale=Theta[Sigma], shape=Theta["xi.reg"])
      wi <- weights(m,x,xi=Theta["xi.reg"], sigma=Theta[Sigma]) #, kap=2)
      Wi[[i]] <- wi
      v <- pBB(wi,z)
      x.new <- evd:::qgpd(v, scale=Theta[Sigma], shape=Theta["xi.reg"])
      atsite_param <- pwmGP(x.new)
      Thetai[,i] <- atsite_param
      #faire un vecteur des m_0
    }#end for
    #ESSAYER D'ESTIMER xi.reg EN POOLANT LES DONNEES PLUT?T QU'EN MOYENNANT
    #PoolPrecip <- na.omit(as.vector(M.new))
    # xi.reg <- pwmGP(PoolPrecip)$xi
    #Theta <- c(Thetai[1,],xi.reg)
    
    #Thetai[1,] <- m0i*(1-xi.reg)
    Theta <-c(Thetai[1,],mean(Thetai[2,]))#regional shape parameter is the mean of at-site estimate
    names(Theta) = c(Scale_names,"xi.reg")
    increment<-abs(Xi.old-Theta["xi.reg"])
    
    if (verba) {
      #NOT ADAPTED TO REGIONALISATION
      cat("\n theta=",round(Theta,3), "\t loop=", loop,"\t increment=",increment)
      log.new<-sum(log(dEGP.BB(x,c(wi,Theta))))
      log.old<-sum(log(dEGP.BB(x,c(wi,sigma.old,xi.old))))
      cat("\t diff log-likelihood=",round(log.new-log.old,5))
    }
  }#end while
  
  
  Param<-list(Wi,Theta)
  return(Param)
}





















