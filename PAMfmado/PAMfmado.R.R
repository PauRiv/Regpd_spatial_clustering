library(cluster)
PAMfmado.R <-function(x,K,J=0,threshold=.99,max.min=0, distance="manhattan"){
  #
  # Description 
  #   This function performs the PAM algorithm based on the F-madogram distance
  #
  #
  # Aguments 
  #   x a matrix
  #   K number of clusters
  #   J number of resampling for which the stations are randomly moved to break the dependence 
  #                    J=0 means no resampling
  #   threshold corresponding to the quantile level for the resampling.
  #   the resulting quantile is stored via sil.final<<-quantile(sil.vec[sil.vec>0],probs=threshold)
  #                    J=0 means no resampling
  # 
  #
  Nnb = ncol(x) 
  Tnb =nrow(x)
  cat("\n Number of stations (Nnb)=",Nnb, "\t time series length (Tnb)=",Tnb,"\t number of cluster=",K)

  #--- DISTANCE MATRIX
  #--- F-MADOGRAM  
  cat("\n Computing the F-madogram distance")
  V = array(NaN, dim = c(Tnb,Nnb)) # ; VV = array(NaN, dim = c(Tnb,Nnb)) ; Z = array(NaN, dim = c(Tnb,Nnb))
  for(p in 1:Nnb) {
    x.vec = as.vector(x[,p])
    x.vec[x.vec < max.min]=NA # thresholding 
    Femp = ecdf(x.vec)(x.vec)
    V[,p] = Femp
  }
  # with Madogram
  DD = dist(t(V),method = distance,diag = TRUE, upper = TRUE)/(2*Tnb)
  # !!! NA are excluded for when computing the distance between two rows but the dist is then multiplied by (Tnb/ nb of non-NA) 
  # in order to have comparable distances (ie wrt the number of data points) among all rows 
  
  #---------- CLUSTERING WITH PAM -------#
  output = pam(DD,K,diss = TRUE,medoids = NULL)
  
  #################################################################################
  ####  resampling step to compute the silhouette coeff for independent case
  #################################################################################
  if(J>0){
    sample.R = function(x) sample(x)
    sil.vec=rep(0,J)
    for(j in 1:J){
      #x.sample=x[,sample(1:Nnb)]
      x.sample=apply(x,2,sample.R) # resampling within each column 
      for(p in 1:Nnb) {
        x.bis = as.vector(x.sample[,p])
        x.bis[x.bis < max.min]=NA # thresholding 
        Femp = ecdf(x.bis)(x.bis)
        V[,p] = Femp
#         Femp = ecdf(x.sample[,p])
#         V[,p] = Femp(x.sample[,p])
      }
      # with Madogram
      DD = dist(t(V),method = "manhattan",diag = TRUE, upper = TRUE)/(2*Tnb)
      sil.vec[j]=pam(DD,K,diss = TRUE,medoids = NULL,keep.diss=FALSE,keep.data=FALSE)$silinfo$avg.width
    }
    sil.final<<-quantile(sil.vec[sil.vec>0],probs=threshold)
  }
  
  return(output)
}
