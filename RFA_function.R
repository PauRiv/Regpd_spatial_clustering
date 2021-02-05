#################### FUNCTIONS FOR RFA : CLUSTERING ON PWM RATIO ######################
##############################     LAST MODIFICATION 05/02/2021            ##############################

######################
## ESTIMATION
sampling<-function(data,thres=1){ 
  #   ARGUMENTS
  # data = matrix of observation. Each row corresponds to a day and each column to a site/station
  #   VALUE
  # X, matrix of STRICTLY POSITIVE observation. Each row corresponds to a day and each column to a site/station
  nday <- nrow(data) #number of observation per station
  nb_stat<-dim(thres)[2] #number of stations
  
  data[data<thres]=NA
  X=data-thres
  return(X)}

xi.Ratio <-function(x,independence=TRUE){
  # Arguments
  #   x = temporal serie for a station
  #   n = length of temporal serie
  #
  # Values
  #   xi.R = ratio estimator
  #   
  #   a0, a1, a2, 3 = Estimator of PWM 0,1,2 and 3
  #   See : Diebolt08, Improving PWM method for the GEV distribution; and Guillou, 2009
  #
  #   Apply to precip data matrix to obtain a vector
  #   R.vect of length=nb_stat, vector of the ratio values for each site
  
  
  x<-sort(x,na.last=NA) #remove Na and sort
  
  if(independence==TRUE){# weights for PWM estimators
    weight <- (c(0:(length(x)-1))-.5)/(length(x)-1)#for PWM with cdf
    surv_weight <- 1-weight#for PWM with survival function
  }
  if(independence==FALSE){F=ecdf(x)
  weight <- F(x)
  surv_weight <- 1-weight#for PWM with survival function
  }
  
  
  a0 <- mean(x)
  a1 <- mean(weight*x)
  a2 <- mean(weight^2*x)
  a3 <- mean(weight^3*x)
  
  xi.R <-(3*a2-2*a1)/(2*a1-a0)
  #if(F) xi.R <- (1/xi.R)
  return(xi.R)
}
R.vec<-function(X){
  #Argument
  #   X=matrix of daily precip, a row = a day, a col = a site
  # Value 
  #   R.vect = vector of R estimates, with length = nb of sites
  
  R.vect <- apply(X, 2,  xi.Ratio)
  return(R.vect)
}

#CLUSTERING FUNCTION
clustering_algo<-function(data,clustering_method, nb_clusters=NULL){
  # ARGUMENTS
  #   data, vector of points to cluster
  #   clustering_method, method used for clustering : "pam", "cah" or "kmeans"
  #   nb_clusters, final number of clusters
  
  # if(clustering_method=="cah"){
  #   ###   distance by pairs of X
  #   #function(u,v){abs(u-v)}) on R.vect vector of ratio values site to site
  #   d.R.mat<-dist(data)
  #   
  #   ###   Hierarchical Clustering
  #   hclust.out<-hclust(d.R.mat)
  #   classes = cutree(hclust.out, k=nb_clusters) 
  #   # OR : number of clusters depending on the criterion
  #   
  # }  #CAH
  # 
  
  # if(clustering_method=="kmeans"){
  #   
  #   ### initialisation and clustering
  #   init.centers <- quantile(data, probs=seq(0,1,length=nb_clusters))
  #   clustz <- kmeans(data,centers=init.centers,iter.max=20)
  #   classes <- clustz$cluster}#kmeans
  
  if(clustering_method=="pam"){
    clusters_info <- pamk(as.vector(data))
   
    #classes <- clusters_info$pam}
    if(length(nb_clusters)!=0){classes <- pam(as.vector(data),nb_clusters,cluster.only = TRUE)}}#pam
  if(clustering_method=="PAMfmado"){
    #extract max
    #precip_max <-
    classif<-PAMfmado.R(x=precip_max,K=nb_clusters)#,max.min = thres_precip)
    clusters<-classif$clustering
  }
  
  
  return(clusters_info)
  #return(classes)
  }#classif


################################
## END FUNCTIONS FOR CLUSTERING ON PWM RATIO
################################


## FUNCTIONS FOR SIMULATION 
# SIMULATION OF DATA

simulEGPD<-function(xi, sigma, kappa=1/10, n=100){
  # Arguments
  #   xi = shape parameter with length= nb regions
  #   n = sample size
  #   
  # Valeur
  #  rEGPD matrix with n rows (a row = a wet day), and number of "stations" (=n.x*n.y) col (a col= an EGPD(xi,sigma,kappa) sample)
  
  m<-length(xi)  # nombre regions
  if(length(sigma)!=m) stop("xi and sigma with different dimensions")
  out<-matrix(NA,n,m)
  for(i in 1:m){
    if(xi[i]>0){out[,i]<-(sigma[i]/xi[i])*((1-(runif(n))^(1/kappa[i]))^(-xi[i])-1)}#temp serie with length n
    if(xi[i]==0){out[,i] <- -sigma[i]*log(1-(runif(n))^(1/kappa[i]))}
  }
  return(out)
}

simulKAPPA<-function(xi, sigma,mu=0,flex=1/10,n=100){
  # Arguments
  #   xi (=-shape1 in hos05)= shape parameter of length nb_reg, sigma = scale, mu = location, flex=h (=shape parameter in hosking, 2005)
  #   n = taille echantillon
  #   h=-1 = GLogistic, h=0 = GEV (not defined here), h=1=GPD
  # Values
  #   rKAPPA matrix with n rows and n.x*n.y (number of "stations") col.
  #   A row = a day; a col=a col= a KAPPA(xi,sigma,mu,flex) sample
  #    de parametres xi et sigma =1, mu et flexibilit?
  
  m<-length(xi)
  if(length(sigma)!=m) stop("xi and sigma with different dimensions")
  out<-matrix(NA,n,m)
  for(i in 1:m){
    out[,i]<-rkappa4(n=n,shape1=-xi[i],shape2=flex[i],scale=sigma[i],location=0)
  }
  return(out)
  
}

# THEORETICAL VALUES OF RATIOS
R.EGPD<-function(kappa, xi){
  #ARGUMENTS 
  # kappa and xi, shape parameters of EGPD
  #VALUE
  # Theoretical value of ratio of PWM defined in Le Gall et al. 2020 for EGPD(kappa, sigma, xi)
  RR = (3*beta(3*kappa, 1-xi)-beta(kappa, 1-xi))/(2*beta(2*kappa, 1-xi)-beta(kappa, 1-xi))
  #if(F) RR<-1/RR
  return(RR)
}
R.EGPD0<-function(kappa){
  #ARGUMENTS 
  # kappa shape parameter of EGPD
  #VALUE
  # Theoretical value of ratio of PWM defined in Le Gall et al. 2020 for EGPD(kappa, sigma, xi=0)
  a0 = kappa*gamma(kappa)
  a1 = a0*(1-2^(-kappa -1))
  a2 = a0*(1-2^(-kappa)+3^(-kappa-1))
  RR = (3*a2-a0)/(2*a1-a0)
  return(RR)
}

R.GPD<-function(xi){
  #ARGUMENTS 
  # xi, shape parameter of GPD
  #VALUE
  # Theoretical value of ratio of PWM defined in Le Gall et al. 2020 for GPD(sigma, xi)
  RGPD<-(5-xi)/(3-xi)
  return(RGPD)
}
R.approx<-function(xi){
  #ARGUMENTS 
  #  xi, shape parameter of GEV
  #VALUE
  # Theoretical value of ratio of PWM defined in Le Gall et al. 2020 for GEV(xi) (and EGPD with kappa -> infty)
  Rlim<-(3^xi-1)/(2^xi-1)
  return(Rlim)
}
RL.th<-function(kap.grid, xi.grid,sigma.grid,p,normalised=TRUE){
  #ARGUMENTS
  #   kap.grid, kappa (shape) parameter of EGPD with G(u)=u^kappa
  #   xi.grid, strictly positive shape parameter
  #   sigma.grid, scale parameter
  #   proba, probability
  
  #VALUE
  #  return level x_p for EGPD(G,xi.grid,sigma.grid) with G(u)=u^kappa
  
  #if(xi.grid==0){y=-sigma.grid*log(1-p^(1/kap.grid))}#sigmaN=sigma if xi=0
  y=sigma.grid/xi.grid*((1-p^(1/kap.grid))^(-xi.grid)-1)
  if(normalised==TRUE){sigmaN<-1/(kap.grid*beta(kap.grid,1-xi.grid)-1)#sigmaNormalised=xi.grid*sigmaN
  y=sigmaN*((1-p^(1/kap.grid))^(-xi.grid)-1)#}
  }
  return(y)
}
RL.th0<-function(kap.grid,sigma.grid,p,normalised=TRUE){
  #ARGUMENTS
  #   kap.grid, kappa (shape) parameter of EGPD with G(u)=u^kappa
  #   sigma.grid, scale parameter
  #   proba, probability
  
  #VALUE
  #  (normalised) return level x_p for EGPD(G,xi.grid=0,sigma.grid) with G(u)=u^kappa
  
  y=-sigma.grid*log(1-p^(1/kap.grid))
  if(normalised==TRUE){y=-log(1-p^(1/kap.grid))}
  return(y)
}
#################
## NORMALIZED VERSIONS OF RATIO (PARAMETRIC APPROACH)
R.norm<-function(kappa, xi){
  
  RN<-(R.EGPD(kappa,xi)-R.approx(xi))/(R.GPD(xi)-R.approx(xi))
  return(RN)
}
R.norm2<-function(kappa, xi){
  RN2<-(R.EGPD(kappa,xi)-R.EGPD(kappa,0.8))/(R.EGPD(kappa,0.00001)-R.EGPD(kappa,0.8))
  return(RN2)
}
R.norm3<-function(kappa,xi){
  # approx of R.norm wth allowed value of xi
  RN3<-(R.EGPD(kappa,xi)-R.approx(0.05))/(R.GPD(0.05)-R.approx(0.05))
  return(RN3)
}
R.norm4<-function(kappa, xi){
  #approx of R.norm2 with allowed value of kappa
  RN4<-(R.EGPD(kappa,xi)-R.EGPD(2,0.8))/(R.EGPD(2,0.00001)-R.EGPD(2,0.8))
  #RN4<-(R.EGPD(kappa,xi)-R.GPD(0.8))/(R.GPD(0.00001)-R.GPD(0.8))
  return(RN4)
}

R.vec2<-function(x){
  # Arguments
  #   x = temporal serie for a station
  #   n = length of temporal serie
  #
  # Values
  #   xi.RN = normalised ratio estimator
  #   
  #   b0, b1, b2, b3 = Estimator of PWM 0,1,2 and 3
  #   See : Diebolt08, Improving PWM method for the GEV distribution; and Guillou, 2009
  #
  #   Apply to precip data matrix to obtain a vector
  #   R.vect of length=nb_stat, vector of the ratio values for each site
  #Argument
  #   X=matrix of daily precip, a row = a day, a col = a site
  # Value 
  #   R.vect2 = vector of R.norm2 estimates, with length = nb of sites
  
  R.vect2 <- apply(X, 2,  xi.Ratio.norm)
  return(R.vect)
}
























