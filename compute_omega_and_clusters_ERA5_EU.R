##############################################################
##### Spatial clustering on ERA5 precipitation data (EU) #####
##############################################################

message("Are you sure you want to run the next line?")
message("It will clean the environment... I mean, virtually :/")
rm(list=ls(all=TRUE))

# Set seed for reproductible work
set.seed(2020)

########## packages required ##########

library(ncdf4) # To handle NetCDF files
library(lubridate) # For easier date manipulations
library(fpc) #pamk, silhouette criterion
library(cluster) # PAM
library(tictoc) #timer
library(foreach);library(iterators);library(parallel);library(doParallel) # parallelization

######### Functions for clustering ######
source("/scratch3/pauline/Regpd_spatial_clustering/RFA_function.R")#RFA clustering, see Le Gall et al., 2020
source("/scratch3/pauline/Regpd_spatial_clustering/PAMfmado.R.R")# F-madogram clustering, see Bador et al., 2015

###############################################
# extr_seas_pos_precip
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
# vector with non NA and >1 precip of the given season 
###############################################

extr_seas_pos_precip <- function(season, precip_timeserie, precip_date){
  stopifnot(length(precip_timeserie)==length(precip_date))
  SEASONS <- c("MAM", "JJA", "SON", "DJF")
  stopifnot(season %in% SEASONS)
  ref_seas <- which(SEASONS==season)
  months_in_season <- (((c(1,2,3)+(3*ref_seas-1))-1)%%12)+1 #gives for ex: c(12,1,2) if season == "DJF"
  
  return(precip_timeserie[(lubridate::month(precip_date)==months_in_season[1]
                           | lubridate::month(precip_date)==months_in_season[2]
                           | lubridate::month(precip_date)==months_in_season[3])
                          &(!is.na(precip_timeserie)) &(precip_timeserie>1)])
}#end for keep_seasonal_positive_precip function

########## Settings for parallelization ##########
nb_cores <- floor(detectCores()/3 - 1)

registerDoParallel(cores=nb_cores) #n_lon=464=16*29



########## Open the netCDF file ##########
## Dimensions=lon*lat*time
nc_precip_ERA5 = nc_open("/scratch3/pauline/EGPD_on_ERA5/Data/era5_daily_precip_1979-01-02_2018-12-31_Europe_land.nc")

#Gridpoint corresponding to the station in ERA5
lon_ERA5 <- ncvar_get(nc = nc_precip_ERA5, varid = "lon")
num_lon <- length(lon_ERA5)
lat_ERA5 <- ncvar_get(nc = nc_precip_ERA5, varid = "lat")
num_lat <- length(lat_ERA5)
times_ERA5 <- ncvar_get(nc = nc_precip_ERA5, varid = "time")
date_ERA5 <- format(as.POSIXlt(times_ERA5*60*60, origin="1979-01-02 12:30:00"), format = "%Y-%m-%d")

nc_close(nc_precip_ERA5)

R_value <- matrix(data = NA, nrow = num_lat, ncol = num_lon)
nc_precip_ERA5 = nc_open("/scratch3/pauline/EGPD_on_ERA5/Data/era5_daily_precip_1979-01-02_2018-12-31_Europe_land.nc")
tic()
precip_ERA5 <- ncvar_get(nc = nc_precip_ERA5, varid = "TP")
toc() #1.5min
nc_close(nc_precip_ERA5)

GP_kept <- !(apply(X=precip_ERA5[,,dim(precip_ERA5)[3]], MARGIN = c(1,2), FUN = is.na))
ref_in_list <- matrix(data = NA, nrow = nrow(GP_kept), ncol = ncol(GP_kept))
dim(precip_ERA5);dim(GP_kept)


nmin_wetdays <- 100
##### Select the season #####
seas = "JJA"

precip_positive_list <- list()

counting <- 0
tic()
for (LON in 1:num_lon) {
  for (LAT in 1:num_lat) {
    
    if(GP_kept[LON,LAT]){
      positive_precip <- extr_seas_pos_precip(season = seas, precip_timeserie = precip_ERA5[LON,LAT,], precip_date = date_ERA5)
      if (length(positive_precip)>nmin_wetdays){
        counting <- counting + 1
        precip_positive_list[[counting]] <- positive_precip
        ref_in_list[LON,LAT] <- counting
      } else {
        GP_kept[LON,LAT] <- FALSE
      }#end if long enough
    }#end if kept
    
  }#end for LAT
  print(paste("processed lon", LON, "out of", num_lon))
}#end for LON
toc()#50min

nb_GP_kept <- sum(GP_kept)

Rval_list <- foreach(GP=1:nb_GP_kept) %dopar% {
  
  RESU <- xi.Ratio(precip_positive_list[[GP]])
  
  attr(RESU, 'id') <- GP #keep track of where I am
  
  return(RESU)
} #end foreach GP 
# 40sec

Rval_vec <- unlist(Rval_list)


# # Store just Rvalue ####
# 
# R_matrix <- matrix(data = NA, ncol = ncol(GP_kept), nrow = nrow(GP_kept))
# counting <- 0
# tic()
# for (LON in 1:num_lon) {
#   for (LAT in 1:num_lat) {
# 
#     if(GP_kept[LON,LAT]){
#       counting <- counting + 1
#       if(counting==attr(Rval_list[[counting]], "id")){
#         R_matrix[LON,LAT] <- Rval_vec[counting]
#       }#end if good ID
#     }#end if kept
# 
#   }#end for LAT
#   print(paste("processed lon", LON, "out of", num_lon))
# }#end for LON
# toc()




# # create and write the netCDF file -- ncdf4 version
# ncfname <- "/scratch3/pauline/Regpd_spatial_clustering/Rvalue_ERA5.nc"
# # define dimensions
# londim <- ncdim_def(name = "longitude", units = "degrees_east", vals = lon_ERA5)
# latdim <- ncdim_def(name = "latitude", units = "degrees_north", vals = lat_ERA5)
# 
# # define variables
# fillvalue <- 1e32
# 
# R_val <- ncvar_def(name = "R_val", units = "no_unit", dim = list(londim, latdim),
#                    missval = fillvalue)
# 
# # create netCDF file and put arrays
# ncout <- nc_create(ncfname,R_val,force_v4=TRUE)
# 
# # put variables
# ncvar_put(ncout, R_val, R_matrix)
# 
# nc_close(ncout)


# compute clusters #####

nb_clusters = c(2:16)
list_partition <- list()
sil_crit = rep(0, length(nb_clusters))
tic()
list_partition <- foreach(ind=1:length(nb_clusters)) %dopar% {
  NbClusters = nb_clusters[ind]
  partition_am = pam(Rval_vec,NbClusters)
  return(partition_am)
}#end foreach ind
toc() #for 16 clusters: #17.7h SON # h DJF # 30h MAM # h JJA
save(list_partition, file = paste0("/scratch3/pauline/Regpd_spatial_clustering/list_partition_16_ERA5_EU",
                                   seas, ".Rdata"))

load(file = paste0("/scratch3/pauline/Regpd_spatial_clustering/list_partition_16_ERA5_EU",
                   seas, ".Rdata"))

for (ind in 1:length(nb_clusters)) {
  sil_crit[ind] <- list_partition[[ind]]$silinfo$avg.width
}#end for ind

# plot silhouette criterion
# optimal silhouette criterion is argument of highest silhouette criterion
plot(1:max(nb_clusters),c(0,sil_crit), xlab = "Number of clusters", ylab= "Average silhouette width",
     main=paste("ERA-5", seas))

# Store Rvalue and cluster #####

R_matrix <- matrix(data = NA, ncol = ncol(GP_kept), nrow = nrow(GP_kept))
cluster_matrix_2 <- matrix(data = NA, ncol = ncol(GP_kept), nrow = nrow(GP_kept))
cluster_matrix_3 <- matrix(data = NA, ncol = ncol(GP_kept), nrow = nrow(GP_kept))
cluster_matrix_4 <- matrix(data = NA, ncol = ncol(GP_kept), nrow = nrow(GP_kept))
cluster_matrix_5 <- matrix(data = NA, ncol = ncol(GP_kept), nrow = nrow(GP_kept))
cluster_matrix_6 <- matrix(data = NA, ncol = ncol(GP_kept), nrow = nrow(GP_kept))
cluster_matrix_7 <- matrix(data = NA, ncol = ncol(GP_kept), nrow = nrow(GP_kept))
cluster_matrix_8 <- matrix(data = NA, ncol = ncol(GP_kept), nrow = nrow(GP_kept))
cluster_matrix_9 <- matrix(data = NA, ncol = ncol(GP_kept), nrow = nrow(GP_kept))
cluster_matrix_10 <- matrix(data = NA, ncol = ncol(GP_kept), nrow = nrow(GP_kept))
cluster_matrix_11 <- matrix(data = NA, ncol = ncol(GP_kept), nrow = nrow(GP_kept))
cluster_matrix_12 <- matrix(data = NA, ncol = ncol(GP_kept), nrow = nrow(GP_kept))
cluster_matrix_13 <- matrix(data = NA, ncol = ncol(GP_kept), nrow = nrow(GP_kept))
cluster_matrix_14 <- matrix(data = NA, ncol = ncol(GP_kept), nrow = nrow(GP_kept))
cluster_matrix_15 <- matrix(data = NA, ncol = ncol(GP_kept), nrow = nrow(GP_kept))
cluster_matrix_16 <- matrix(data = NA, ncol = ncol(GP_kept), nrow = nrow(GP_kept))
counting <- 0
tic()
for (LON in 1:num_lon) {
  for (LAT in 1:num_lat) {
    
    if(GP_kept[LON,LAT]){
      counting <- counting + 1
      if(counting==attr(Rval_list[[counting]], "id")){
        R_matrix[LON,LAT] <- Rval_vec[counting]
        cluster_matrix_2[LON,LAT] <- list_partition[[1]]$clustering[counting]
        cluster_matrix_3[LON,LAT] <- list_partition[[2]]$clustering[counting]
        cluster_matrix_4[LON,LAT] <- list_partition[[3]]$clustering[counting]
        cluster_matrix_5[LON,LAT] <- list_partition[[4]]$clustering[counting]
        cluster_matrix_6[LON,LAT] <- list_partition[[5]]$clustering[counting]
        cluster_matrix_7[LON,LAT] <- list_partition[[6]]$clustering[counting]
        cluster_matrix_8[LON,LAT] <- list_partition[[7]]$clustering[counting]
        cluster_matrix_9[LON,LAT] <- list_partition[[8]]$clustering[counting]
        cluster_matrix_10[LON,LAT] <- list_partition[[9]]$clustering[counting]
        cluster_matrix_11[LON,LAT] <- list_partition[[10]]$clustering[counting]
        cluster_matrix_12[LON,LAT] <- list_partition[[11]]$clustering[counting]
        cluster_matrix_13[LON,LAT] <- list_partition[[12]]$clustering[counting]
        cluster_matrix_14[LON,LAT] <- list_partition[[13]]$clustering[counting]
        cluster_matrix_15[LON,LAT] <- list_partition[[14]]$clustering[counting]
        cluster_matrix_16[LON,LAT] <- list_partition[[15]]$clustering[counting]
      }#end if good ID
    }#end if kept
    
  }#end for LAT
  print(paste("processed lon", LON, "out of", num_lon))
}#end for LON
toc()#19 sec





# create and write the netCDF file -- ncdf4 version
ncfname <- paste0("/scratch3/pauline/Regpd_spatial_clustering/Rvalue_clusterpam_16_ERA5_EU",seas,".nc")
# define dimensions
londim <- ncdim_def(name = "longitude", units = "degrees_east", vals = lon_ERA5)
latdim <- ncdim_def(name = "latitude", units = "degrees_north", vals = lat_ERA5)

# define variables
fillvalue <- 1e32

R_val <- ncvar_def(name = "R_val", units = "no_unit", dim = list(londim, latdim),
                   missval = fillvalue)

cluster_PAM_2 <- ncvar_def(name = "cluster_num_2", units = "no_unit", dim = list(londim, latdim),
                           missval = fillvalue)
cluster_PAM_3 <- ncvar_def(name = "cluster_num_3", units = "no_unit", dim = list(londim, latdim),
                           missval = fillvalue)
cluster_PAM_4 <- ncvar_def(name = "cluster_num_4", units = "no_unit", dim = list(londim, latdim),
                           missval = fillvalue)
cluster_PAM_5 <- ncvar_def(name = "cluster_num_5", units = "no_unit", dim = list(londim, latdim),
                           missval = fillvalue)
cluster_PAM_6 <- ncvar_def(name = "cluster_num_6", units = "no_unit", dim = list(londim, latdim),
                           missval = fillvalue)
cluster_PAM_7 <- ncvar_def(name = "cluster_num_7", units = "no_unit", dim = list(londim, latdim),
                           missval = fillvalue)
cluster_PAM_8 <- ncvar_def(name = "cluster_num_8", units = "no_unit", dim = list(londim, latdim),
                           missval = fillvalue)
cluster_PAM_9 <- ncvar_def(name = "cluster_num_9", units = "no_unit", dim = list(londim, latdim),
                           missval = fillvalue)
cluster_PAM_10 <- ncvar_def(name = "cluster_num_10", units = "no_unit", dim = list(londim, latdim),
                            missval = fillvalue)
cluster_PAM_11 <- ncvar_def(name = "cluster_num_11", units = "no_unit", dim = list(londim, latdim),
                            missval = fillvalue)
cluster_PAM_12 <- ncvar_def(name = "cluster_num_12", units = "no_unit", dim = list(londim, latdim),
                            missval = fillvalue)
cluster_PAM_13 <- ncvar_def(name = "cluster_num_13", units = "no_unit", dim = list(londim, latdim),
                            missval = fillvalue)
cluster_PAM_14 <- ncvar_def(name = "cluster_num_14", units = "no_unit", dim = list(londim, latdim),
                            missval = fillvalue)
cluster_PAM_15 <- ncvar_def(name = "cluster_num_15", units = "no_unit", dim = list(londim, latdim),
                            missval = fillvalue)
cluster_PAM_16 <- ncvar_def(name = "cluster_num_16", units = "no_unit", dim = list(londim, latdim),
                            missval = fillvalue)

# create netCDF file and put arrays
ncout <- nc_create(ncfname,
                   list(R_val, cluster_PAM_2, cluster_PAM_3, cluster_PAM_4, cluster_PAM_5,
                        cluster_PAM_6, cluster_PAM_7, cluster_PAM_8, cluster_PAM_9,
                        cluster_PAM_10, cluster_PAM_11, cluster_PAM_12
                        , cluster_PAM_13, cluster_PAM_14, cluster_PAM_15, cluster_PAM_16
                   ),
                   force_v4=TRUE)

# put variables
ncvar_put(ncout, R_val, R_matrix)
ncvar_put(ncout, cluster_PAM_2, cluster_matrix_2)
ncvar_put(ncout, cluster_PAM_3, cluster_matrix_3)
ncvar_put(ncout, cluster_PAM_4, cluster_matrix_4)
ncvar_put(ncout, cluster_PAM_5, cluster_matrix_5)
ncvar_put(ncout, cluster_PAM_6, cluster_matrix_6)
ncvar_put(ncout, cluster_PAM_7, cluster_matrix_7)
ncvar_put(ncout, cluster_PAM_8, cluster_matrix_8)
ncvar_put(ncout, cluster_PAM_9, cluster_matrix_9)
ncvar_put(ncout, cluster_PAM_10, cluster_matrix_10)
ncvar_put(ncout, cluster_PAM_11, cluster_matrix_11)
ncvar_put(ncout, cluster_PAM_12, cluster_matrix_12)
ncvar_put(ncout, cluster_PAM_13, cluster_matrix_13)
ncvar_put(ncout, cluster_PAM_14, cluster_matrix_14)
ncvar_put(ncout, cluster_PAM_15, cluster_matrix_15)
ncvar_put(ncout, cluster_PAM_16, cluster_matrix_16)


nc_close(ncout)

