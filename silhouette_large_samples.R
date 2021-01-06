################# SILHOUETTE FOR LARGE SAMPLES 
###############   ONLY ON NEAREST POINTS OR ONLY BY CONSIDERING DISTANCE TO THE MEDOIDS #################
#################                17/12/2020           ##################
## PACKAGES
library(cluster)#pam
library(fpc)#pamk, silhouette criterion
## SILHOUETTE ON LESS POINTS
reduced_silhouette <- function(R_vector,clusters_vec, medoids, nb_points = 100){
  # COMPUTE THE SILHOUETTE COEFFICIENT ON FEWER POINTS THAN THE ORIGINAL SUBSET
  #   ARGUMENTS
  # R_vector, vector of points to be classified
  # clusters_vec, vector containing the cluster indice for each point of the dataset
  # medoids, matrix or vector of indices of medoids
  # nb_points, number of points in the subset used to compute silhouette coefficient
  #   VALUE
  # sil_reduced, silhouette coefficient computed on the subset of the nb_points points closest to their medoids

  ind_subset = matrix(NA,nrow = nb_points,ncol = length(medoids))
  dist_mat = outer(R_vector,R_vector, FUN=function(u,v){abs(u-v)})
  dist_mat = as.dist(dist_mat)
  d.R.mat = as.matrix(dist_mat)
  for (NoCluster in 1:length(c(medoids))) {
    medoid = medoids[NoCluster]
    minimal_dist = sort(c(d.R.mat[,medoid]))[1:nb_points]
    for (i in 1:nb_points) {
      ind_subset[i,NoCluster] = which(c(d.R.mat[,medoid])==minimal_dist[i])
    }#end loop nb_point
  
  }#end NoCluster
  sub_clust_vec = clusters_vec[sort(c(ind_subset))]
  R_vect_sub = R_vector[sort(c(ind_subset))]
  sub_dist_mat = outer(R_vect_sub,R_vect_sub, FUN=function(u,v){abs(u-v)})
  sub_dist_mat = as.dist(sub_dist_mat)
 
  sil_reduced = silhouette(sub_clust_vec,sub_dist_mat)#summary(sil_REGPD) for partition silhouette criterion
  return(list("asw" = mean(sil_reduced[,3]),"silhouette_values"= sil_reduced))
}
#silhouette_medoid <- function(R_vector,)






























