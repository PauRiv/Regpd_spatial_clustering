#Always clean before doing anything 
rm(list=ls(all=TRUE))
source("~/Thèse/Codes/Suisse/reg_suisse_EGPD_ratio.R")
source("~/Thèse/Codes/ClusterMax/R/PAMfmado.R.R")
##################################################################
### Chargement des donnees (precipitations et coordonnees)
##################################################################
#bassin 
load("~/Thèse/Codes/Suisse/MetadataPrec_SEL.RData")#coord,nb_an_obs,nom_stations
#switzerland
load("~/Thèse/Codes/Suisse/matricePluieTriee.RData")# precip data for 191 stations, 85years
load("~/Thèse/Codes/Suisse/TableauCoordTrie.RData")#tableau coord 
load("~/Thèse/Codes/Suisse/StationsSummerParEGPD.RData")#at site parameters
load("~/Thèse/Codes/Suisse/StationsWinterParEGPD.RData")#at site parameters
load("~/Thèse/Codes/Suisse/StationsFallParEGPD.RData")#at site parameters
load("~/Thèse/Codes/Suisse/StationsSpringParEGPD.RData")#at site parameters

###################################################
# Parametres : methode classif
###################################################
seuil_precip<-1#only wet days
mode_classif<-"pam" # "cah" # or "kmeans" # or "PAMfmado"
saison<-"toutes" #"ete" "printemps" "automne" "hiver" "toutes"
ratio<-"classic"#"alternative" or "classic"
###############################################################
### Extracting data
###############################################################
#Fint <- matrixPrec_SEL#basin
#Fint <- matrixPrec#switzerland tot, whole year
Fint<-matricePluieTriee
nb_stat<-dim(Fint)[2]

vec_date = as.Date(row.names(Fint))
vec_month = format(vec_date,"%m")

if(saison=="ete"){
  season="summer"
  Fint<-matricePluieTriee[vec_month%in%c("06","07","08"),]
}
if(saison=="automne"){
  season="fall"
  Fint<-matricePluieTriee[vec_month%in%c("09","10","11"),]
}
if(saison=="hiver"){
  season="winter"
  Fint<-matricePluieTriee[vec_month%in%c("12","01","02"),]
}
if(saison=="printemps"){
  season="spring"
  Fint<-matricePluieTriee[vec_month%in%c("03","04","05"),]
}
if(saison=="toutes"){
  season="all"
  Fint<-matricePluieTriee
}

#tableau.coord=data.frame(x=MetadataPrec[,5],y=MetadataPrec[,6])#coord of stations
nday <- nrow(Fint)
nb_stat<-dim(Fint)[2]#=nombre stations
nb_wet_days_a_year<-mean(apply(Fint>1, FUN =sum, MARGIN = 2 ),na.rm=TRUE)/(nday/365.25)
nb_wet_days_summer<-mean(apply(Fint>1, FUN =sum, MARGIN = 2 ),na.rm=TRUE)/(nday/92)
nb_wet_days_fall<-mean(apply(Fint>1, FUN =sum, MARGIN = 2 ),na.rm=TRUE)/(nday/91)
#############################################################
#Sampling : x = mat of precip. a row = a day, a col=a station [useless if EGPD]
x<- echantillonnage(donnees=Fint,seuil=seuil_precip)-1 #sampling of data

if(ratio=="classic"){
  R.vect<-R.vec(x) #vector of ratio for each station
}else{R.vect<-R.vec.alt(x)}

#############################################################
## Optimal number of clusters 
#############################################################
if(mode_classif=="PAMfmado"){
  nb_reg<-pamk(na.omit(as.vector(x)))$nc
}
if(mode_classif=="pam"){
  nb_reg<-pamk(na.omit(as.vector(R.vect)))$nc
}#fin pam

if(mode_classif=="cah"){
  d.xi.mat<-outer(na.omit(as.vector(R.vect)),na.omit(as.vector(R.vect)),FUN=function(u,v){abs(u-v)})
  d.xi.mat<-as.dist(d.xi.mat)
  hclustxi.out<-hclust(d.xi.mat)
  #tracer inertie 
  inertie <- sort(hclustxi.out$height, decreasing = TRUE)
  plot(inertie[1:6], type = "s", xlab = "Nombre de classes", ylab = "Inertie")
  nb_reg<-best.cutree(hclustxi.out)#nb classe optimal en terme d'inertie
}#fin cah

if(mode_classif=="kmeans"){
  ratio_ss <- data.frame(cluster = seq(from = 1, to = 9, by = 1))
  for (k in 1:9) {
    km_model <- kmeans(na.omit(R.vect), k, nstart = 20)#code original = pas de na.omit
    ratio_ss$ratio[k] <- km_model$tot.withinss / km_model$totss #inertie intraclasse normalisee
  }
  
  ggplot(ratio_ss, aes(cluster, ratio)) +
    geom_line() +
    geom_point()
  #If inertia not plotted, re-do ggplot
  print("nb_reg ?")}# FIN Kmeans

#####################################################
## Classif 
#####################################################
if(mode_classif=="PAMfmado"){
  classif<-PAMfmado.R(x=Fint,K=3)#,max.min = seuil_precip)
  classif<-classif$clustering
} else{
  classif<-classific(na.omit(R.vect),mode_classif, nb_reg)#vector of clusters ofr meth_classif method
  
}

row.names(tableau.coord.trie)<-names(R.vect)



########################################################
## Mapping of clustered stations
########################################################
#tableau.coord=data.frame(x=MetadataPrec[,5],y=MetadataPrec[,6])#coord of stations
Bassin<-readOGR("~/Thèse/Codes/Suisse/BASIN1_CH1903_LV03.shp") #st_read() in sf package

Swizerland<-readOGR("~/Thèse/Codes/Suisse/CHE_adm0_CH1903_LV03.shp")


#x11()
plot(Swizerland,main=paste("Season= ", season, ", clustering=", mode_classif,", ratio=", ratio ))#,saison))#Switzerland
plot(Bassin,add=TRUE)#Basin 
points(tableau.coord.trie, col = classif,pch=16,cex=1.2)#clustered stations 
legend(x="topright", legend=unique(classif), col=unique(classif), pch=16,bg="white",bty="n")


#normalisation of precip by their sigma (index-flood)
vect_moy=apply(Fint, MARGIN = 2,FUN = mean,na.rm=TRUE)
if(length(vect_moy)==nb_stat){
  mat_moy=matrix(data=rep(vect_moy,nday),ncol=nb_stat,byrow=TRUE)#mat nb_days*nb_stat in (i,j)=mean of precip at station j
  mat_pluie_norm=Fint/mat_moy #normalised (by at-site mean) rainfall intensities
}


#regional estimation of shape parameters kappa and xi
parameter_clust1=fit.extgp(data=na.omit(mat_pluie_norm[,classif=1]),method="pwm",init=c(0,0,0), model=1)#model=1 ie Gu)=u^\kappa
#parameter_clust1=fit.extgp(data=na.omit(Fint[,classif=1]),method="pwm",init=c(0,0,0), model=1)#model=1 ie Gu)=u^\kappa

#kappa       sigma          xi (normalised)
#0.612240515 2.791833664 0.002161607 

parameter_clust2=fit.extgp(data=na.omit(mat_pluie_norm[,classif=2]),method="pwm",init=c(0,0,0), model=1)#model=1 ie Gu)=u^\kappa
#parameter_clust2=fit.extgp(data=na.omit(Fint[,classif=2]),method="pwm",init=c(0,0,0), model=1)#model=1 ie Gu)=u^\kappa

#kappa      sigma         xi (normalised) 
#0.7382127 2.4648298  0.2557899 

#estim for sub-region1 of cluster 2
parameter_clust21=fit.extgp(data=na.omit(mat_pluie_norm[,(classif=2)&(sub_classif2=3)]),method="pwm",init=c(0,0,0), model=1)#model=1 ie Gu)=u^\kappa
#estimation station by station of kappa and xi

#estim for sub-region2 of cluster 2
parameter_clust22=fit.extgp(data=na.omit(mat_pluie_norm[,(classif=2)&(sub_classif2=4)]),method="pwm",init=c(0,0,0), model=1)#model=1 ie Gu)=u^\kappa

#estim for sub-region2 of cluster 3
parameter_clust23=fit.extgp(data=na.omit(mat_pluie_norm[,(classif=2)&(sub_classif2=5)]),method="pwm",init=c(0,0,0), model=1)#model=1 ie Gu)=u^\kappa

#estim for sub-region2 of cluster 4
parameter_clust24=fit.extgp(data=na.omit(mat_pluie_norm[,(classif=2)&(sub_classif2=6)]),method="pwm",init=c(0,0,0), model=1)#model=1 ie Gu)=u^\kappa

#estim for sub-region2 of cluster 5
parameter_clust25=fit.extgp(data=na.omit(mat_pluie_norm[,(classif=2)&(sub_classif2=7)]),method="pwm",init=c(0,0,0), model=1)#model=1 ie Gu)=u^\kappa

#estim for sub-region2 of cluster 6
parameter_clust26=fit.extgp(data=na.omit(mat_pluie_norm[,(classif=2)&(sub_classif2=8)]),method="pwm",init=c(0,0,0), model=1)#model=1 ie Gu)=u^\kappa



#######################################################
#at-site parameters on a map
#######################################################

#parameter_point_winter=apply(X=na.omit(Fint), MARGIN =2, FUN = fit.extgp,method="pwm", init=c(0,0,0), model=1)#margin=2 ie applied to columns
#save(parameter_point_summer,file='~/Thèse/Codes/Suisse/ParametresStationEGPD_ete.RData')

season="all"
Bassin<-readOGR("~/Thèse/Codes/Suisse/BASIN1_CH1903_LV03.shp") #st_read() dans le package sf
Swizerland<-readOGR("~/Thèse/Codes/Suisse/CHE_adm0_CH1903_LV03.shp")

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


Bassin<-readOGR("~/Thèse/Codes/Suisse/BASIN1_CH1903_LV03.shp") #st_read() in sf package
Swizerland<-readOGR("~/Thèse/Codes/Suisse/CHE_adm0_CH1903_LV03.shp")

plot(Swizerland, main="At-site 30-years return levels")
plot(Bassin,add=TRUE)#Bassin d'etude

for(i in 1:T){
  #kappa sigma or xi
  selection<-subset(return_stat, (col.x[i+1] >=return_levels )& (return_levels> col.x[i]))
  points(selection, pch= 21, col=0, bg=colorT[i],cex=1.5) #,cex=exp(strength))
}
legend(x="left",legend=round(col.x[-1],2), text.col=colorT)


#unclustered ratio on a map

Bassin<-readOGR("~/Thèse/Codes/Suisse/BASIN1_CH1903_LV03.shp") #st_read() in sf package
Swizerland<-readOGR("~/Thèse/Codes/Suisse/CHE_adm0_CH1903_LV03.shp")

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
atsite.lmoment<-data.frame(regsamlmu(Fint),classif)
atsite.lmom<-regsamlmu(Fint)

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
  mu_V<-mean_sim
  sd_V<-sd_sim
  Dcrit[NbCluster,]<-Dcritical
  D[[NbCluster]]<-Disc_measure
  
}
ind_disc_stat<- c(which(D[[1]]>3),which(D[[2]]>3))#stations with discordancy measures above the 5% confidence level

points(tableau.coord.trie[ind_disc_stat,], col = 3, pch=16,cex=1.2) #discordant stations

legend(x="topright", legend=c(unique(classif),"disc"), col=c(unique(classif),3), pch=16,bg="white",bty="n")
REGPD_Dcrit<-Dcrit
REGPD_D<-D
REGPD_H_values<-H_values
REGPD_mu_V<-mu_V
REGPD_sd_V<-sd_V
REGPD_Vobs<-V_obs



#save(REGPD_Dcrit,file = '~/Thèse/Codes/Suisse/RFA_Hos05/REGPD_hos05/REGPD_critical_measure.RData')
#save(REGPD_H_values,file = '~/Thèse/Codes/Suisse/RFA_Hos05/REGPD_hos05/REGPD_statistic_H.RData')
#save(REGPD_mu_V,file = '~/Thèse/Codes/Suisse/RFA_Hos05/REGPD_hos05/REGPD_theo_mean_V.RData')
#save(REGPD_sd_V,file = '~/Thèse/Codes/Suisse/RFA_Hos05/REGPD_hos05/REGPD_theo_sd_V.RData')
#save(REGPD_Vobs,file = '~/Thèse/Codes/Suisse/RFA_Hos05/REGPD_hos05/REGPD_empirical_V.RData')

####################################################################################
## END
####################################################################################
######
# Sub clustering (USELESS)
###############################"
# first : pam

# optimal number of clusters
clusterNo<-2
if(mode_classif=="PAMfmado"){
  nb_sub_reg<-pamk(na.omit(R.vect[classif==clusterNo]))$nc
}
if(mode_classif=="pam"){
  nb_sub_reg<-pamk(na.omit(R.vect[classif==clusterNo]))$nc
}#fin pam

if(mode_classif=="cah"){
  d.xi.mat<-outer(na.omit(R.vect[classif==clusterNo]),na.omit(R.vect[classif==clusterNo]),FUN=function(u,v){abs(u-v)})
  d.xi.mat<-as.dist(d.xi.mat)
  hclustxi.out<-hclust(d.xi.mat)
  #tracer inertie 
  inertie <- sort(hclustxi.out$height, decreasing = TRUE)
  plot(inertie[1:6], type = "s", xlab = "Nombre de classes", ylab = "Inertie")
  nb_sub_reg<-best.cutree(hclustxi.out)#nb classe optimal en terme d'inertie
}#fin cah

if(mode_classif=="kmeans"){
  ratio_ss <- data.frame(cluster = seq(from = 1, to = 9, by = 1))
  for (k in 1:9) {
    km_model <- kmeans(na.omit(R.vect[classif==clusterNo]), k, nstart = 20)#code original = pas de na.omit
    ratio_ss$ratio[k] <- km_model$tot.withinss / km_model$totss #inertie intraclasse normalisee
  }
  
  ggplot(ratio_ss, aes(cluster, ratio)) +
    geom_line() +
    geom_point()
  #Si pas de plot de l'inertie, refaire tourner le ggplot
  print("nb_sub_reg ?")}# FIN Kmeans

#####################################################
## sub-clustering (USELESS)
#####################################################
if(mode_classif=="PAMfmado"){
  sub_classif<-PAMfmado.R(x=Fint,K=3)#,max.min = seuil_precip)
  sub_classif<-classif$clustering
} else{
  sub_classif<-classific(na.omit(R.vect[classif==clusterNo]),mode_classif, nb_sub_reg)#creation des classes de R.vect pour nb_reg regions et pour la methode mode_classif
  
}
if(clusterNo==1){
  sub_classif1<-sub_classif
}
if(clusterNo==2){
  sub_classif2<-sub_classif+2
}
tableau.coord.plus<-data.frame(cbind(tableau.coord.trie,R.vect,classif))
###End sub-clustering
#plotting sub-regions
Bassin<-readOGR("~/Thèse/Codes/Suisse/BASIN1_CH1903_LV03.shp") #st_read() dans le package sf

Swizerland<-readOGR("~/Thèse/Codes/Suisse/CHE_adm0_CH1903_LV03.shp")


#sub-reg one
plot(Swizerland,main=paste("Season= ", season, ", clustering=", mode_classif,", ratio=", ratio ))#,saison))#Suisse
plot(Bassin,add=TRUE)#Bassin d'etude
points(tableau.coord.trie[classif==1,], col = 1,pch=16,cex=1.2)#stations classees
points(tableau.coord.trie[classif==2,], col = sub_classif2+3,pch=16,cex=1.2)#stations classees
legend(x="topright", legend=unique(sub_classif1+2), col=unique(sub_classif1+2), pch=16,bg="white",bty="n")

#sub_reg 2

plot(Swizerland,main=paste("Season= ", season, ", clustering=", mode_classif,", ratio=", ratio ))#,saison))#Suisse
plot(Bassin,add=TRUE)#Bassin d'etude
points(tableau.coord.trie[classif==2,], col = sub_classif2+3,pch=16,cex=1.2)#stations classees
legend(x="topright", legend=unique(sub_classif2+2), col=unique(sub_classif2+3), pch=16,bg="white",bty="n")


    