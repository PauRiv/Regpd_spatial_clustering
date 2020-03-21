#Always clean before doing anything 
#rm(list=ls(all=TRUE))
#set.seed(2018)

######################
###  packages
######################
library(cluster) #PAM
library(rgdal)#map (fun : readOGR)
library(ggplot2)# graphes
library(fpc)#pamk : silhouette criterion
library(factoextra)#optimal number of cluster for kmeans algo
library(JLutils)#best.cutree
library(mev)#fit EGPD
######################
# Sampling
######################
#
# Data= Matrix of daily precip (in mm)
#A row= a day, a col= a site
#rownames= aaaammjj
#
###########################
echantillonnage<-function(donnees,seuil=1){ 
  nday <- nrow(donnees) #number of observation per station
  nb_stat<-dim(donnees)[2] #number of stations
  
  donnees[donnees<seuil]=NA
  X=donnees
  return(X)}

###################################################################
## Fonction du ratio
##################################################################
xi.Ratio <-function(x){
  # Arguments
  #   x = temporal serie for a station
  #   n = length of temporal serie
  #
  # Values
  #   xi.R = ratio estimator
  #   
  #   b0, b1, b2, b3 = Estimator of PWM 0,1,2 and 3
  #   See : Diebolt08, Improving PWM method for the GEV distribution; and Guillou, 2009
  #
  #   Apply to precip data matrix to obtain a vector
  #   R.vect of length=nb_stat, vector of the ratio values for each site
  
  
  x<-sort(x,na.last=NA) #remove Na and sort
  
  # weights for PWM estimators
  poids1 <- (c(0:(length(x)-1))-.5)/(length(x)-1)
  weight1 <- 1-poids1

  
  b0 <- mean(x)
  b1 <- mean(weight1*x)
  b2 <- mean(weight1^2*x)
  b3 <- mean(weight1^3*x)
  
  xi.R <-(3*b2-b0)/(2*b1-b0)
  # option xi.R<-(3*B3-B2)/(2*B2-B1)
  xi.R <- (1/xi.R)
 
  return(xi.R)
}
xi.ratio.alt <-function(x){
  # Arguments
  #   x = temporal serie for a station
  #   n = length of temporal serie
  # Values
  #   xi.R = ratio estimator
  #   
  #   b0, b1, b2, b3 = Estimator of PWM 0,1,2 and 3
  #   See : Diebolt08, Improving PWM method for the GEV distribution; and Guillou, 2009
  #
  #   Apply to precip data matrix to obtain a vector
  #   R.vect of length=nb_stat, vector of the ratio values for each site
  
  x<-sort(x,na.last=NA) #remove NA and sort
  
  #weights for PWM estimators
  poids1 <- (c(0:(length(x)-1))-.5)/(length(x)-1)
  weight1 <- 1-poids1
  
  
  b0 <- mean(x)
  b1 <- mean(weight1*x)
  b2 <- mean(weight1^2*x)
  b3 <- mean(weight1^3*x)
  
  xi.R <-(2*b2-b1)/b2
  # deuxieme ratio xi.R<-(3*B3-B2)/(2*B2-B1)
  #xi.R <- (1/xi.R)
  
  return(xi.R)
}


###################################################
#Estimation of ratio for each site
###################################################
R.vec<-function(X){
  R.vect <- apply(X, 2,  xi.Ratio)
  return(R.vect)
}
R.vec.alt<-function(X){
  R.vect.alt <- apply(X, 2,  xi.ratio.alt)
  return(R.vect.alt)
}

###################################################
## Clustering : return cluster for each station 
##################################################
# clustering of R to obtain cluster of sites with same distribution (except scale-factor)
#

classific<-function(X,mode_classif, k.clust){
  
  if(mode_classif=="cah"){
    ###   distance by pairs of X
    #function(u,v){abs(u-v)}) on R.vect vector of ratio values site to site
    d.R.mat<-dist(X)
    
    ###   Hierarchical Clustering
    hclust.out<-hclust(d.R.mat)
    classes = cutree(hclust.out, k=k.clust) 
    # OR : number of clusters depending on the criterion
    
  }  #CAH
  
  
  if(mode_classif=="kmeans"){
    
    ### initialisation and clustering
    init.centers <- quantile(X, probs=seq(0,1,length=k.clust))
    clustz <- kmeans(X,centers=init.centers,iter.max=20)
    classes <- clustz$cluster}#kmeans
  
  if(mode_classif=="pam"){
    classes <- pam(X,k.clust,cluster.only = TRUE)}#pam
  
  
  
  return(classes)}#classif
