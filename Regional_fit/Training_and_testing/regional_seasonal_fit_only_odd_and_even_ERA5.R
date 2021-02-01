######################################################################
##### Regional inference of EGPD on ERA5 precipitation data (EU) #####
######################################################################
rm(list=ls(all=TRUE))

########## packages required ##########
library(ncdf4) # To handle NetCDF files


seas <- "JJA"
print(seas)
# if (seas=="JJA") {nclust <- 3; var_name="cluster_num_3"} else {nclust <- 2; var_name="cluster_num_2"}

nclust <- 2; var_name="cluster_num_2"

########## load functions, cluster data & precip ##########

#Partition
nc_clusters <- nc_open(paste0("/scratch3/pauline/Regpd_spatial_clustering/Rvalue_clusterpam_16_ERA5_EU", seas, "V1.nc"))
optimal_partition <- ncvar_get(nc_clusters, varid = var_name)
lon_clus <- ncvar_get(nc_clusters, varid = "longitude")
lat_clus <- ncvar_get(nc_clusters, varid = "latitude")
nc_close(nc_clusters)

#regional fitting function
source("/scratch3/pauline/Regpd_spatial_clustering/Regio_Fit_halfdata/parallelized_regional_inference_functions.R")

#precip data
## Dimensions=lon*lat*time
nc_precip_ERA5 = nc_open("/scratch3/pauline/EGPD_on_ERA5/Data/era5_daily_precip_1979-01-02_2018-12-31_Europe_land.nc")

#Gridpoint corresponding to the station in ERA5
lon_ERA5 <- ncvar_get(nc = nc_precip_ERA5, varid = "lon")
num_lon <- length(lon_ERA5)
lat_ERA5 <- ncvar_get(nc = nc_precip_ERA5, varid = "lat")
num_lat <- length(lat_ERA5)
times_ERA5 <- ncvar_get(nc = nc_precip_ERA5, varid = "time")
date_ERA5 <- format(as.POSIXlt(times_ERA5*60*60, origin="1979-01-02 12:30:00"), format = "%Y-%m-%d")
precip_ERA5 <- ncvar_get(nc = nc_precip_ERA5, varid = "TP")
nc_close(nc_precip_ERA5)

###############################################
# extr_seas_precip
# Description
# 
# Extract from a precip time serie the non NA and >1 precip of a given season 
# 
# Arguments
# season               string in in c("SON", ""DJF", "MAM", "JJA")
# precip_timeserie  	 precipitation time serie to process
# precip_date          dates corresponding to the precip_timeserie (same length)
#
# Value
# vector with non NA and >0 precip of the given season 
###############################################

extr_seas_pos_precip <- function(season, precip_timeserie, precip_date, thshld = 0){
  stopifnot(length(precip_timeserie)==length(precip_date))
  SEASONS <- c("MAM", "JJA", "SON", "DJF")
  stopifnot(season %in% SEASONS)
  ref_seas <- which(SEASONS==season)
  months_in_season <- (((c(1,2,3)+(3*ref_seas-1))-1)%%12)+1 #gives for ex: c(12,1,2) if season == "DJF"
  
  return(precip_timeserie[(lubridate::month(precip_date)==months_in_season[1]
                           | lubridate::month(precip_date)==months_in_season[2]
                           | lubridate::month(precip_date)==months_in_season[3]) & 
                            (precip_timeserie>thshld)])
}#end for keep_seasonal_positive_precip function


#Perform the regional fitting
# dim(precip_ERA5)
# M=precip_ERA5[214:276,104,1:3500]
# Mpositive <- matrix(nrow=0,ncol=0)
# for (sta in 1:nrow(M)) {
#   y <- M[sta,]
#   y <- y[y>0]
#   Mpositive = rbind.fill.matrix(Mpositive,t(as.matrix(y)))
# }
# Mpos <- t(Mpositive)
# Q <- fitEGPDkSemiRegCensoredIter(Mpos)
# dim(Q)

List_reg_fit <- list()
List_keep_coord <- list()
for(clust_nb in 1:nclust){
  M_pos_precip <- matrix(nrow=0,ncol=0)
  keep_lon <- numeric()
  keep_lat <- numeric()
  
  counting <- 0
  for (LON in 1:num_lon) {
    for (LAT in 1:num_lat) {
      
      if(LON%%2==LAT%%2) {
        if(!is.na(optimal_partition[LON,LAT]) & optimal_partition[LON,LAT]==clust_nb){
            counting <- counting + 1
            print(counting)
            keep_lon[counting] <- LON ; keep_lat[counting] <- LAT
            y <- extr_seas_pos_precip(precip_timeserie = precip_ERA5[LON,LAT,], season = seas, precip_date = date_ERA5, thshld = 0)
            ind_1_third <- (1:length(y))[which(((1:length(y))%%3)==0)]
            M_pos_precip <- rbind.fill.matrix(M_pos_precip, t(as.matrix(y[ind_1_third])))
          }#end if in cluster
        
        }#end if ODD or EVEN
      }#end for LAT
    }#end for LON : took 1h
    rm(counting)
    
    List_keep_coord[[clust_nb]] <- list(keep_lon = keep_lon, keep_lat = keep_lat)
    List_reg_fit[[clust_nb]] <- fitEGPDkSemiRegCensoredIter_paral(M = t(M_pos_precip),
                                                                  ncores = floor(detectCores()/3 - 1), cens_thres = c(1,Inf))
        
      
}#end for clust_nb

save(List_reg_fit, file = paste0("/scratch3/pauline/Regpd_spatial_clustering/Regio_Fit_halfdata/res_reg_fit_only_even_OR_odd_",seas,"_1third.Rdata"))
save(List_keep_coord, file = paste0("/scratch3/pauline/Regpd_spatial_clustering/Regio_Fit_halfdata/keep_coord_reg_fit_only_even_OR_odd_",seas,"_1third.Rdata"))
