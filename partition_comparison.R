################### COMPARE PARTITIONS WITH SEVERAL CRITERIA
library(clusterCrit) # package with numerous clustering validation crits

compute_criteria <- function(R.vect, krange, partition){
  ## COMPUTE SEVERAL CRITERIA TO EVALUATE PARTITION
  #   ARGUMENTS
  # R.vect, vector of points to be classified
  # krange, the number of clusters in each partition
  # partition, a list of partition. The i-th vector of the list 
  #            is a partition with i different clusters
  #   VALUE
  # dataframe containing the values of each criterion (silhouette, dunn,
  #                                 xie_beni, s_dbw, davies_bouldin)
  #   for each partition

  cluster_crits = c("silhouette","dunn","xie_beni","s_dbw","davies_bouldin") #criteria
  crit_summary <- data.frame(array(NA, c(max(krange),5))) # to store the crits
  names(crit_summary) <- cluster_crits
  
  for (k in krange) {
    result=intCriteria(traj = as.matrix(R.vect), part = partition[[k]] , cluster_crits) # function in Cluster crit package
    crit_summary$silhouette[k] = result$silhouette
    crit_summary$dunn[k] = result$dunn
    crit_summary$xie_beni[k] = result$xie_beni
    crit_summary$s_dbw[k] = result$s_dbw
    crit_summary$davies_bouldin[k] = result$davies_bouldin
    
  }
  
  return(crit_summary)

  }# end function

#optimal number of clusters according to each criterion
# Maximize: silhouette, dunn  
which.max(crit_summary$silhouette)
which.max(crit_summary$dunn)
# Minimize: xie_beni, s_dbw and davies_bouldin
which.min(crit_summary$xie_beni)
which.min(crit_summary$s_dbw)
which.min(crit_summary$davies_bouldin)
#plot criteria
palette = c('#1f78b4','#33a02c','#e31a1c','#ff7f00','#6a3d9a')
plot(x=krange, y = summary_criterion$silhouette,ylim=c(0,max(summary_criterion)),
     ylab="criterion value",
     col=palette[1], type="b", pch = 16)
for (i in 2:5) {
  color = palette[i]
  lines(x=krange, y = summary_criterion[,i],col=color,type = "b", pch=16)
  
}
legend("topright",legend = cluster_crits,col = palette,pch=16)



































