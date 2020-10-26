##############################################################
##### Spatial clustering on ERA5 precipitation data (CH) #####
##############################################################

message("Are you sure you want to run the next line?")
message("It will clean the environment... I mean, virtually :/")
rm(list=ls(all=TRUE))

# Set seed for reproductible work
set.seed(2020)


########## packages required ##########

library(ncdf4) # To handle NetCDF files
library(lubridate) # For easier date manipulations
library(fpc)#pamk, silhouette criterion
library(cluster)# PAM
##
#PACKAGES FOR PLOT
# library(rworldmap)
# library(raster)#relief map
# library(dplyr)#colors for clusters
######### Functions for clustering ######
source("RFA_function.R")#RFA clustering, see Le Gall et al., 2020
source("PAMfmado.R.R")# F-madogram clustering, see Bador et al., 2015

######### Open the netCDF file ##########
# Location on Pauline's Laptop:
folder_data <- "C:/Users/admin/Documents/Regpd_spatial_clustering_LOCAL/"

#Location on Philomène's Laptop
folder_data <- "~/Thèse/Collaboration_Pauline/data/"

# use the function nc_open() to open a netcdf file
nc_precip_ERA5 = nc_open(filename = paste0(folder_data, "era5_daily_precip_1979-01-02_2018-12-31_CH.nc"))
# use the function ncvar_get to extract a variable
precip_ERA5_CH <- ncvar_get(nc = nc_precip_ERA5, varid = "TP")
lon_ERA5 <- ncvar_get(nc = nc_precip_ERA5, varid = "lon")
lat_ERA5 <- ncvar_get(nc = nc_precip_ERA5, varid = "lat")
times_ERA5 <- ncvar_get(nc = nc_precip_ERA5, varid = "time")
# Change time to a more convinient format
date_ERA5 <- format(as.POSIXlt(times_ERA5*60*60, origin="1979-01-02 12:30:00"), format = "%Y-%m-%d")
# Clsoe the ncdf once you've extract all you wanted
nc_close(nc_precip_ERA5)
############################
## CLUSTERING ALGO
# ESTIMATING PWM RATIO
serie_temp = precip_ERA5_CH[1,1,]
threshold = 1

x <- sampling(data=serie_temp,thres=threshold) #sampling of data
R.vect<-xi.Ratio(x) #vector of ratio for each station

#extract positive precipitation ERA-5 switz
positive_precip <- apply(X = precip_ERA5_CH, FUN = sampling, MARGIN = c(1,2))
#compute PWM ratio on positive precipitation
R_CH <- apply(X=positive_precip, FUN = xi.Ratio, MARGIN = c(2,3))#c(2,3) because exchanges dimension

image(R_CH)

# CLUSTER ESTIMATED VALUES OF RATIO

spatial_clusters = clustering_algo(data=R_CH,clustering_method="pam")#later: clustering_method="PAMfmado" for clustering on spatial dependence only
image(matrix(spatial_clusters$pamobject$clustering, nrow=23))


# PLOT CLUSTERED SITES ON A MAP
#PYTHON ?
# reliefData <- stack("~/Thèse/Codes/Suisse/Mapping/HYP_HR_SR_OB_DR.tif")
# newext <- c(-10, 25, 40, 53)
# reliefData.c <- crop(reliefData, newext)
# 
# newmap = rworldmap::getMap(resolution = "low")
# sp::plot(newmap, xlim = c(6, 11), ylim = c(45, 48), asp = 1)
# plotRGB(reliefData.c, add=TRUE)
# sp::plot(newmap, xlim = c(6, 11), ylim = c(45, 48), add = TRUE)
# graphics::points(loc[,x], loc[,y], pch=18, col=colors, cex=sizes)


# Matrix with NAs:
R_with_NA <- matrix(c(NA, 1.669390, 1.665978 ,1.659580,
                      NA ,1.672934 ,1.700153 ,NA,
                      NA , NA , 1.731162, 1.711460), nrow = 4)
image(R_with_NA)
cluster_with_NA <- clustering_algo(data=R_with_NA,
                                   clustering_method="pam")
# PAM LOOP FOR SEVERAL NUMBER OF CLUSTERS
# list of pam objects (clustering vector, silhouette criteria, medoid centers)
list_partition <- list()
# number of clusters in the partition
nb_clusters = c(2:5)
sil_crit = rep(0, length(nb_clusters))
for (i in 1:length(nb_clusters)) {
  NbClusters = nb_clusters[i]
  list_partition[[i]] = pam(R_CH,NbClusters)
  sil_crit[i] = list_partition[[i]]$silinfo$avg.width
}
# plot silhouette criterion
plot(1:max(nb_clusters),c(0,sil_crit), xlab = "Number of clusters", ylab= "Average silhouette width")#optimal silhouette criterion is argument of highest silhouette criterion

## F-MADOGRAM CLUSTERING (see Bador et al., 2015, Bernard et al., 2013, Saunders, 2020)
# M is the matrix of annual/seasonal/weekly maxima
K=2 #3,4,5 # number of clusters
Spatial_clustering_Fmad = PAMfmado.R(M,K) #vector of clusters, medoids...
dep_clusters = Spatial_clustering_Fmad$clustering #vector of clusters

## REGIONAL INFERENCE

m = 30 #degree of BB 
M1 <- #Matrix of precip in Cluster 1 #idem for clusters 2..5

EGPDm_C1 = fitGPreg(M1,m=m) #idem for the others clusters








