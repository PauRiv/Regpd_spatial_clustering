#################### MAPPING CLUSTERS OF STATIONS OBTAINED BY PWM RATIO PAM CLUSTERING (margins) AND PAMfmado (dependence) #############
############################                     LAST MODIFICATION: 04/09/2020              #################
#Always clean before doing anything 
rm(list=ls(all=TRUE))
source("RFA_function.R")#functions for RFA clustering, see Le Gall et al., 2020
source("PAMfmado.R.R")# function for F-madogram clustering, see Bador et al., 2015
library(fpc)#pamk, silhouette criterion
# PLOT PACKAGES
#library(rgdal)#map (fun : readOGR)
#library(RColorBrewer)#map colors
#library(rworldmap)
#library(raster)#relief map
#library(dplyr)#colors for clusters
##################################################################
### LOAD DATA (PRECIPITATION AND COORDINATES OF SITES)
##################################################################
#switzerland
#load("~/Thèse/Codes/Suisse/RFA_REGPD/matricePluieTriee.RData")# precip data for 191 stations, 85years
#load("~/Thèse/Codes/Suisse/RFA_REGPD/TableauCoordTrie.RData")#coordinates

###################################################
# PARAMETERS FOR CLUSTERING
###################################################
thres_precip <- 1#only wet days
clustering_method<-"pam"  #"pam" # "cah" # ou "kmeans" # or "PAMfmado"
season <- "all"
ratio <-"classic"# "classic" for moments of order 0 to 2 or "alternative" for higher moments
independence=FALSE#if TRUE, F(X) is estimated by (i/n)
###############################################################
### EXTRACTING DATA
###############################################################
#Matrix of observation, each row corresponds to a day and each column to a site
Fint <- matricePluieTriee 
nb_stat<-ncol(Fint)
nday <- nrow(Fint)
# Coordinates table, dataframe with 2 column (x and y coordinates), row names = codes/names of stations
coord_table <- tableau.coord.trie

#EXTRACT SEASONAL PRECIPITATION
## METS PEUT ETRE TA FONCTION DE SEUILLAGE ET D'EXTRACTION SELON LES SAISONS ICI
vec_date = as.Date(row.names(Fint))
vec_month = format(vec_date,"%m")
if(season=="summer"){
  Fint<-Fint[vec_month%in%c("06","07","08"),]
}
if(season=="fall"){
  Fint<-Fint[vec_month%in%c("09","10","11"),]
}
if(season=="winter"){
  Fint<-Fint[vec_month%in%c("12","01","02"),]
}
if(season=="spring"){
  Fint<-Fint[vec_month%in%c("03","04","05"),]
}
vec_date_season = as.Date(row.names(Fint))
vec_year=format(vec_date_season,"%y")


if(clustering_method=="PAMfmado"){
  #extract annual maxima
  A<-matrix(data=NA,nrow=length(unique(vec_year)),ncol = nb_stat)
  tm_NA<- list()
  for (j in 1:nb_stat) {
    for (i in 1:length(unique(vec_year))) {
      NoYear<-unique(vec_year)[i]
      A[i,j]<-max(Fint[vec_year==NoYear,j],na.rm=TRUE)
    }
    #indices of station with year(s) without any max
    if(sum(is.infinite(A[,j]))>0){tm_NA=c(tm_NA,j)}
  }
  
  #tm_NA <- unlist(tm_NA,use.names = FALSE)
  #Fint <- A[,-tm_NA]
  #coord_table <- coord_table[-tm_NA,]
}

#nb_wet_days_a_year<-mean(apply(Fint>thres_precip, FUN =sum, MARGIN = 2 ),na.rm=TRUE)/(nday/365.25)
#nb_wet_days_summer<-mean(apply(Fint>thres_precip, FUN =sum, MARGIN = 2 ),na.rm=TRUE)/(nday/92)
#nb_wet_days_fall<-mean(apply(Fint>thres_precip, FUN =sum, MARGIN = 2 ),na.rm=TRUE)/(nday/91)

#############################################################
#SAMPLING: x = mat of precip. a row = a day, a col=a station 
x<- sampling(data=Fint,thres=thres_precip) #sampling of data
# ESTIMATING PWM RATIO 
if(ratio=="classic"){
  R.vect<-R.vec(x) #vector of ratio for each station
}else{R.vect<-R.vec.alt(x)}

#############################################################
## OPTIMAL NUMBER OF CLUSTERS
#############################################################
if(clustering_method=="PAMfmado"){
  nb_reg<-pamk(as.vector(x))$nc
}
if(clustering_method=="pam"){
  nb_reg<-pamk(as.vector(R.vect))$nc
}#end pam
if(clustering_method=="cah"){
  d.xi.mat<-outer(as.vector(R.vect),as.vector(R.vect),FUN=function(u,v){abs(u-v)})
  d.xi.mat<-as.dist(d.xi.mat)
  hclustxi.out<-hclust(d.xi.mat)
  #plot inertia 
  inertie <- sort(hclustxi.out$height, decreasing = TRUE)
  plot(inertie[1:6], type = "s", xlab = "Nombre de classes", ylab = "Inertie")
  nb_reg<-best.cutree(hclustxi.out)#nb classe optimal en terme d'inertie
}#end cah
if(clustering_method=="kmeans"){
  ratio_ss <- data.frame(cluster = seq(from = 1, to = 9, by = 1))
  for (k in 1:9) {
    km_model <- kmeans(na.omit(R.vect), k, nstart = 20)#code original = pas de na.omit
    ratio_ss$ratio[k] <- km_model$tot.withinss / km_model$totss #inertie intraclasse normalisee
  }
  
  ggplot(ratio_ss, aes(cluster, ratio)) +
    geom_line() +
    geom_point()
  #If inertia not plotted, re-do ggplot
  print("nb_reg ?")}# END Kmeans

#####################################################
## CLUSTERING
#####################################################
if(clustering_method=="PAMfmado"){
  classif<-PAMfmado.R(x=Fint,K=nb_reg)#,max.min = thres_precip)
  clusters<-classif$clustering
} else{
  clusters<-classific(R.vect,clustering_method, nb_reg)#vector of clusters for meth_classif method
  clusters_obj <- pam(R.vect,k=nb_reg)
  }
#nb_reg=2
row.names(coord_table) <- names(R.vect)

########################################################
## MAP OF CLUSTERED STATIONS
########################################################
# WORLD MAP (WITH RELIEF)
visu_clusters = function(clusters, loc, x, y, title, xlim = "auto", ylim = "auto", resuming = median,
                         interest = 'LAT', silhouette = FALSE, kls = NULL, palette = "Set1",
                         save = FALSE, path = NULL){
  #' Visualisation of clustering on an European map
  #'
  #' @param clusters clustering object, as returned by PAM or kls_clusters_extremes or clusters_extremes
  #' @param loc data frame containing the localisation of each object (at least two columns, one of latitudes and one of longitudes)
  #' @param x name of the column containing the longitudes (for the x-axis)
  #' @param y name of the column containing the latitudes (for the x-axis)
  #' @param title title to be given to the plot
  #' @param xlim x limits for the plot, if xlim = "auto" then there are automatically computed, otherwise should be of the form xlim = c(xmin, xmax)
  #' @param ylim y limits for the plot, if ylim = "auto" then there are automatically computed, otherwise should be of the form ylim = c(ymin, ymax)
  #' @param resuming function to be used to determining the ordering, with respect to the column "interest". If median, then the median of the "interest" of each cluster will be computed, and the colors will be assigned with respect to this.
  #' @param interest variable to consider for the ordering of colors. See more details on "resuming"
  #' @param silhouette wether to display the points with a size proportional to their silhouette width.
  #' @param kls KLs between each points, should only be giiven if silhouette = TRUE.
  #' @param palette which color palette to use (from RBrewer)
  #' @param save wether to save the resulting plot or not.
  #' @param path path of the file when you want to save it.
  #'
  #' @return plot of the clustering on an European map
  if(silhouette){
    if (sum(kls<0, na.rm = TRUE) >  0){
      kls = kls-min(kls, na.rm=TRUE)
    }
    dist = stats::as.dist(kls)
  }else{sizes = rep(1,length(clusters$clustering))}
  
  if(length(class(clusters)) == 2){
    unique_clustering = clusters
    clusters = list()
    clusters[[1]] = unique_clustering
  }
  
  for(i in 1:length(clusters)){
    
    clusters_item = clusters[[i]]
    
    if(silhouette){
      sil = cluster::silhouette(clusters_item, dist)
      widths = sil[order(as.integer(names(sil[,3]))),3]
      sizes = (0.85)*(widths+1)/2+0.15
    }
    
    title_item = title[i]
    if(!is.null(path)){
      path_item = path[i]
    }
    
    nb_clusters = length(clusters_item$medoids)
    list_colors = RColorBrewer::brewer.pal(n = nb_clusters, name = palette)
    
    clustered_loc = loc
    clustered_loc$clustering = clusters_item$clustering
    clustered_loc$var = loc[,interest]
    
    summary_by_cluster <- clustered_loc %>%
      group_by(clustering) %>%
      dplyr::summarize(position = resuming(var))
    spatial_order_clusters = summary_by_cluster[,"clustering"]
    temp = seq(1,nb_clusters)
    order = order(summary_by_cluster$position)
    ordering = rep(1,nb_clusters)
    ordering[order] = temp
    spatial_order_clusters$order = ordering
    clustered_loc = left_join(clustered_loc, spatial_order_clusters)
    list_colors = c("red", "black")#for visible sites in 2 clusters
    colors =list_colors[clustered_loc$order]
    
    reliefData <- stack("~/Th?se/Codes/Suisse/Mapping/HYP_HR_SR_OB_DR.tif")
    newext <- c(-10, 25, 40, 53)
    reliefData.c <- crop(reliefData, newext)
    
    newmap = rworldmap::getMap(resolution = "low")
    sp::plot(newmap, xlim = c(6, 11), ylim = c(45, 48), asp = 1, main=title_item)
    plotRGB(reliefData.c, add=TRUE)
    sp::plot(newmap, xlim = c(6, 11), ylim = c(45, 48), add = TRUE)
    graphics::points(loc[,x], loc[,y], pch=18, col=colors, cex=sizes)
    graphics::points(loc[as.integer(clusters_item$medoids),x], loc[as.integer(clusters_item$medoids),y], pch=5)
    if(silhouette){
      avg_sil = summary(sil)$avg.width
      #graphics::text(-5,51, paste("Average silhouette width:", round(avg_sil,4)),
      #                 font = 2)
      legend(6.1,48.2, legend=c("Average silhouette width:", round(avg_sil,4)), cex=0.8)
    }
    if (save){
      grDevices::dev.copy(png,path_item,width=1200,height=800,res=150)
      grDevices::dev.off()
    }
  }
}
visu_clusters(clusters_obj, MetaDataTriee[,], 'longitude', 'latitude', NULL,
              interest = "latitude", silhouette = FALSE, kls = kls)#, save = TRUE, path = path)

# LOCAL MAP (NO RELIEF)
coord_table$clustering = clusters
## ALWAYS THE SAME ORDER FOR COLORS
summary_by_cluster <- coord_table %>%
  group_by(clustering) %>%
  dplyr::summarize(position = median(y))
spatial_order_clusters = summary_by_cluster[,"clustering"]#contient la liste des No des clusters
temp = seq(1,nb_reg)
order = order(summary_by_cluster$position)
ordering = rep(1,nb_reg)
ordering[order] = temp#ordre des clusters selon les coordonn?es des pts des medoides
spatial_order_clusters$order = ordering #on ajoute la colone des num?ros de clusters
coord_table = left_join(coord_table, spatial_order_clusters) #concatener avec le tableau des metadata
Bassin<-readOGR("~/Th?se/Codes/Suisse/BASIN1_CH1903_LV03.shp") #st_read() in sf package
Swizerland<-readOGR("~/Th?se/Codes/Suisse/CHE_adm0_CH1903_LV03.shp")

palette="Set1"
list_colors = brewer.pal(n = nb_reg, name = palette)
#list_colors = c("#377EB8","#E41A1C")

#x11()
plot(Swizerland)#,main=paste("Season = ", season, ", clustering =", clustering_method,", ratio =", ratio ))#))#Switzerland
plot(Bassin, add=TRUE)
#plotRGB(reliefData,add=TRUE)#Basin 
#avec longitude et latitude
points(coord_table, col = list_colors[clusters],pch=18,cex=1.2)#clustered stations 
legend(x="topright", legend=unique(clusters), col=unique(list_colors[clusters]), pch=18,bg="white",bty="n")
#########################################
##     END CLUSTERING
#########################################

#########################################
##     INFERENCE AND HOMOGENEITY TESTS
#########################################

#########################################################
## REGIONAL ESTIMATION OF EGPD(kappa,sigma,xi) PARAMETERS
#########################################################
# Precipitation in each two clusters
MatPluieClust1<-x[,clusters==1]
MatPluieClust2<-x[,clusters==2]

#save(MatPluieClust2,file="~/Th?se/Codes/Suisse/precip_cluster2.RData")

#regional estimation of shape parameters kappa and xi
#parameter_clust1=fit.extgp(data=na.omit(mat_pluie_norm[,clusters=1]),method="pwm",init=c(0,0,0), model=1)#model=1 ie Gu)=u^\kappa
#parameter_clust1=fit.extgp(data=na.omit(Fint[,classif=1]),method="pwm",init=c(0,0,0), model=1)#model=1 ie Gu)=u^\kappa

#kappa       sigma          xi (normalised)
#0.612240515 2.791833664 0.002161607 

#parameter_clust2=fit.extgp(data=na.omit(mat_pluie_norm[,clusters=2]),method="pwm",init=c(0,0,0), model=1)#model=1 ie Gu)=u^\kappa
#parameter_clust2=fit.extgp(data=na.omit(Fint[,classif=2]),method="pwm",init=c(0,0,0), model=1)#model=1 ie Gu)=u^\kappa

#kappa      sigma         xi (normalised) 
#0.7382127 2.4648298  0.2557899 

#######################################################
#at-site parameters on a map
#######################################################

#parameter_point_winter=apply(X=na.omit(Fint), MARGIN =2, FUN = fit.extgp,method="pwm", init=c(0,0,0), model=1)#margin=2 ie applied to columns
#save(parameter_point_summer,file='~/Th?se/Codes/Suisse/ParametresStationEGPD_ete.RData')

season="all"
Bassin<-readOGR("~/Th?se/Codes/Suisse/BASIN1_CH1903_LV03.shp") #st_read() dans le package sf
Swizerland<-readOGR("~/Th?se/Codes/Suisse/CHE_adm0_CH1903_LV03.shp")

#plot(Swizerland, main="At-site xi parameter")
#plot(Bassin,add=TRUE)#Basin
para=3 #3 for kappa, #4 for sigma #5 for xi
K<-2
#old.par = par(no.readonly = TRUE) 
colorK=rainbow(K)

colorK <- terrain.colors(K) # rainbow(K)
colorK <- heat.colors(K)[K:1] # rainbow(K)
colorK <- rainbow(K+3)[K:1] # rainbow(K)
#z<-parametre_stations_spring[,para]
z<-val_station[,para]
col.z<-seq(from=min(z)-.0001,to=max(z), length.out = K+1)
seuil<-0.75

col.z<-c(0, seuil,2)

plot(Swizerland, main=paste0("At-site " ,colnames(val_station)[para], " parameter, "," season = ", season))#, seuil))

plot(Bassin,add=TRUE)#Basin

for(i in 1:K){
  #kappa sigma or xi
  selection<-subset(val_station, (col.z[i+1] >=kappa )& (kappa> col.z[i]))
  points(selection, pch= 21, col=0, bg=colorK[i],cex=1.5) #,cex=exp(strength))
}
legend(x="topleft",legend=round(col.z[-1],2), text.col=colorK)




#at-site return levels

return_level_EGPD<-function(p,vec_parametre){
  y=((1-p^(1/vec_parametre[1]))^(-vec_parametre[3])-1)*vec_parametre[2]/vec_parametre[3]
  return(y)
}
T=10
colorT <- rainbow(T+3)[T:1] # rainbow(K)

vec_kap=as.vector(val_station[,"kappa"])
vec_sig=as.vector(val_station[,"sigma"])
vec_ksi=as.vector(val_station[,"xi"])
plot(Swizerland, main="At-site return levels")
p=1-1/(30*nb_wet_days_a_year)#1-1/T

return_levels<-apply(X=parametre_stations[,-c(1,2)],return_level_EGPD,p=p, MARGIN=1 )

return_stat<-data.frame(cbind(parametre_stations,return_levels))

x<-return_stat[,6]
col.x<-seq(from=min(x)-.0001,to=max(x), length.out = T+1)


Bassin<-readOGR("~/Th?se/Codes/Suisse/BASIN1_CH1903_LV03.shp") #st_read() in sf package
Swizerland<-readOGR("~/Th?se/Codes/Suisse/CHE_adm0_CH1903_LV03.shp")

plot(Swizerland, main="At-site 30-years return levels")
plot(Bassin,add=TRUE)#Bassin d'etude

for(i in 1:T){
  #kappa sigma or xi
  selection<-subset(return_stat, (col.x[i+1] >=return_levels )& (return_levels> col.x[i]))
  points(selection, pch= 21, col=0, bg=colorT[i],cex=1.5) #,cex=exp(strength))
}
legend(x="left",legend=round(col.x[-1],2), text.col=colorT)


#unclustered ratio on a map

Bassin<-readOGR("~/Th?se/Codes/Suisse/BASIN1_CH1903_LV03.shp") #st_read() in sf package
Swizerland<-readOGR("~/Th?se/Codes/Suisse/CHE_adm0_CH1903_LV03.shp")

plot(Swizerland, main="At-site xi parameter")
plot(Bassin,add=TRUE)#Basin
#old.par = par(no.readonly = TRUE) 

#Choose your favorite color palette
colorK=rainbow(K)

colorK <- terrain.colors(K) 
colorK <- heat.colors(K)[K:1] 
colorK <- rainbow(K+3)[K:1] 
z<-val_station[,6]
col.z<-seq(from=min(z)-.0001,to=max(z), length.out = K+1)
col.z<-c(0,0.776,2)

plot(Swizerland, main=paste0("At-site unclustered ratio"))

plot(Bassin,add=TRUE)#Basin

for(i in 1:K){
  selection<-subset(val_station, (col.z[i+1] >=R.vect )& (R.vect> col.z[i]))
  points(selection, pch= 21, col=0, bg=colorK[i],cex=1.5) #,cex=exp(strength))
}
legend(x="topright",legend=round(col.z[-1],2), text.col=colorK)


###########################################################################
### To be deleted ; homogeneity tests on clusters
###########################################################################
library(lmomRFA)
atsite.lmoment<-data.frame(regsamlmu(x),classif)
atsite.lmom<-regsamlmu(x)

H_values<-matrix(NA,nrow=nb_reg,ncol = 3)
V_obs<-H_values#for three statistical tests (cf Hosking&Wallis 2005)
mu_V<-V_obs#matrix with nb_reg rows and 3 col, mean of H for each cluster and each statistic of test
sd_V<-mu_V#idem but standard deviation
Dcrit<-matrix(NA,nrow = nb_reg, ncol = 2)
D<-list()
info_test=list()
for (NbCluster in 1:nb_reg) {
  test_info_reg<-regtst(atsite.lmoment[classif==NbCluster,-dim(atsite.lmoment)[2]])
  #info_test[NbCluster]<-test_info_reg
  statistics<-test_info_reg$H
  mean_sim<-test_info_reg$vbar
  sd_sim<-test_info_reg$vsd
  Vobs<-test_info_reg$vobs
  Dcritical<-test_info_reg$Dcrit
  Disc_measure<-test_info_reg$D
  H_values[NbCluster,]<-statistics
  V_obs[NbCluster,]<-Vobs
  mu_V[NbCluster,]<-mean_sim
  sd_V[NbCluster,]<-sd_sim
  Dcrit[NbCluster,]<-Dcritical
  D[[NbCluster]]<-Disc_measure
  
}
ind_disc_stat<- c(which(D[[1]]>3),which(D[[2]]>3))#stations with discordancy measures above the 5% confidence level

points(coord_table[ind_disc_stat,], col = 3, pch=16,cex=1.2) #discordant stations

legend(x="topright", legend=c(unique(clusters),"disc"), col=c(unique(clusters),3), pch=16,bg="white",bty="n")
REGPD_Dcrit<-Dcrit
REGPD_D<-D
REGPD_H_values<-H_values
REGPD_mu_V<-mu_V
REGPD_sd_V<-sd_V
REGPD_Vobs<-V_obs
REGPD_test_info_reg<-test_info_reg

#save(REGPD_test_info_reg,file = '~/Th?se/Codes/Suisse/RFA_Hos05/REGPD_hos05/REGPD_homo_test.RData')

####################################################################################
## END
####################################################################################
