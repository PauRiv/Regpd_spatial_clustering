library(cluster)
library(foreach);library(iterators);library(parallel);library(doParallel)
distance.fmado.R <-function(x, K, max.min=0, distance="manhattan", ncores = 24){
  #
  # Description 
  #   This function performs the PAM algorithm based on the F-madogram distance, prior to the PAM algorithm
  #
  #
  # Aguments 
  #   x a matrix with block maxima (each col=one station, each row=one day)
  #   distance a charactere specifying the name of the distance
  #   ncores an positive integer, the number of cores to use in the paralelization
  #
  registerDoParallel(cores=ncores)
  
  Nnb = ncol(x) 
  Tnb =nrow(x)
  cat("\n Number of stations (Nnb)=",Nnb, "\t time series length (Tnb)=",Tnb,"\t number of cluster=",K)
  
  #--- DISTANCE MATRIX
  #--- F-MADOGRAM  
  V = array(NaN, dim = c(Tnb,Nnb))
  
  list_V <- foreach(p = 1:Nnb) %dopar% {
    x.vec = as.vector(x[,p])
    x.vec[x.vec < max.min]=NA # thresholding 
    Femp = ecdf(x.vec)(x.vec)
    return(Femp)
  } #end foreach list_V
  
  for(p in 1:Nnb) {
    V[,p] = list_V[[p]]
  }
  # with Madogram
  DD = dist(t(V),method = distance,diag = TRUE, upper = TRUE)/(2*Tnb)
  # !!! NA are excluded for when computing the distance between two rows but the dist is then multiplied by (Tnb/ nb of non-NA) 
  # in order to have comparable distances (ie wrt the number of data points) among all rows 
}
