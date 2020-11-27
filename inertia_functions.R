##################
## INERTIA
# PRELIMINARY FUNCTIONS : WITHIN SUM OF SQUARES, BETWEEN SUM OF SQUARES AND TOTAL SUM OF SQUARES

withinss<-function(data,NbCluster){
  # ARGUMENTS 
  #   data, dataframe containing  vector of points and vector of corresponding clusters 
  #   NbCluster, name/number of the cluster
  
  # VALUE
  #   wss, within cluster NbClust sum of squares
  #   G, isobarycenter of cluster NbCluster
  #   points, points of the cluster 
  #   cluster, vector of cluster
  x<-data[,1]
  cluster<-data[,2]
  G=mean(x[cluster==NbCluster])#isobarycenter of cluster NbCluster
  wss<-sum((x[cluster==NbCluster]-G)^2)
  return(within.ss=list("wss"=wss,"isobary"=G,"points"=x,"cluster"=cluster))
}
tot.withinss<-function(data){
  # ARGUMENTS 
  #   data, dataframe containing vector of points and vector of corresponding clusters
  # VALUE
  #   tot.wss, total within sum square
  #   partial.bary, vector of isobarycenters of clusters
  number_clusters=length(unique(data[,2]))
  partial.wss<-rep(0,number_clusters)
  partial.bary<-partial.wss
  for (NbCluster in 1:number_clusters) {
    partialwss<-withinss(data,NbCluster)
    partial.wss[NbCluster]<-partialwss$wss
    partial.bary[NbCluster]<-partialwss$isobary
  }
  tot.wss=sum(partial.wss)
  return(list("tot.wss"=tot.wss,"partial.bary"=partial.bary))
}
between.ss<-function(data){
  x=data[,1]#points
  G=mean(x)#isobary
  number_clusters=length(unique(data[,2]))
  n<-rep(0,number_clusters)
  partial.bary<-tot.withinss(data)$partial.bary
  for (NbCluster in 1:number_clusters) {
    n[NbCluster]=sum(data[,2]==NbCluster)
  }
  bss=sum(n*(G-partial.bary)^2)
  return(bss)
}
tot.ss<-function(data){
  # ARGUMENTS 
  #   data, dataframe with two columns containing vector of points and vector of corresponding clusters
  # VALUE
  #   ratio =  total sum of squares 
  tot.wss=tot.withinss(data)$tot.wss
  bet.ss=between.ss(data)
  totss=tot.wss+bet.ss
  return(totss)
}
# INERTIA FUNCTION : RATIO OF SUMS OF SQUARES
ratio_ss<-function(data){
  # ARGUMENTS 
  #   data, A dataframe with two columns:
  #         - column 1 contains vector of points (e.g. estimates of REGPD)
  #         - column 2 contains the vector of corresponding clusters
  # VALUE
  #   ratio =  total sum of squares 
  ratio<-tot.withinss(data)$tot.wss /tot.ss(data)
  return(ratio)
}
#EXAMPLE
# R.vect is the vector of REGPD estimates at sites n1, n2, ...
# classif is the vector of cluster (number of the cluster of site n1, number of the cluster of n2,..)

# REGPD_data<-t(data.frame(rbind(R.vect,classif),row.names = c("R","cluster")))

# inertia _REGPD is the inertia for the number of clusters provided in "classif" 
# inertia_REGPD<-ratio_ss(REGPD_data)

