
rm(list = ls())
seas <- "JJA"

load(paste0("/scratch3/pauline/Regpd_spatial_clustering/Regio_Fit_halfdata/Test_fitting/test_AD_leftGP_",seas, ".Rdata"))
load(paste0("/scratch3/pauline/Regpd_spatial_clustering/Regio_Fit_halfdata/Test_fitting/test_AD_leftGP_", seas, "_LocalFit.Rdata"))
load(paste0("/scratch3/pauline/Regpd_spatial_clustering/Regio_Fit_halfdata/keep_coord_reg_fit_only_odd_", seas, "_1third.Rdata"))


library("ncdf4")


nc_precip_ERA5 = nc_open("/scratch3/pauline/EGPD_on_ERA5/Data/era5_daily_precip_1979-01-02_2018-12-31_Europe_land.nc")
lon_ERA5 <- ncvar_get(nc = nc_precip_ERA5, varid = "lon")
num_lon <- length(lon_ERA5)
lat_ERA5 <- ncvar_get(nc = nc_precip_ERA5, varid = "lat")
num_lat <- length(lat_ERA5)
nc_close(nc_precip_ERA5)



matrix_regio_pvalAD <- matrix(data = NA, nrow = num_lon, ncol = num_lat)
matrix_regio_pvalAD005 <- matrix(data = NA, nrow = num_lon, ncol = num_lat)
matrix_local_pvalAD <- matrix(data = NA, nrow = num_lon, ncol = num_lat)
matrix_local_pvalAD005 <- matrix(data = NA, nrow = num_lon, ncol = num_lat)
matrix_cluster <- matrix(data = NA, nrow = num_lon, ncol = num_lat)

n_clust <- length(List_keep_coord)

for (clust in 1:n_clust) {
  n_sta <- length(List_keep_coord[[clust]]$keep_lon)
  for (sta in 1:n_sta) {
    ref_lon <- List_keep_coord[[clust]]$keep_lon[sta]
    ref_lat <- List_keep_coord[[clust]]$keep_lat[sta]
    
    matrix_regio_pvalAD[ref_lon,ref_lat] <- (pval_test_ad_with_left[[clust]])[sta]
    matrix_local_pvalAD[ref_lon,ref_lat] <- (pval_test_ad_with_left_local[[clust]])[sta]
    
    matrix_cluster[ref_lon,ref_lat] <- clust
    
  }#end for sta
}#end for clust

matrix_regio_pvalAD005 <- matrix_regio_pvalAD>0.05
matrix_local_pvalAD005 <- matrix_local_pvalAD>0.05



ncfname <- paste0("/scratch3/pauline/Regpd_spatial_clustering/Regio_Fit_halfdata/Test_fitting/pvalue_AD_test_reg_fit_",seas,".nc")

# define dimensions
londim <- ncdim_def(name = "longitude", units = "degrees_east", vals = lon_ERA5)
latdim <- ncdim_def(name = "latitude", units = "degrees_north", vals = lat_ERA5)

# define variables
fillvalue <- 1e32

pval_Fit_AD <- ncvar_def(name = "pval", units = "no_unit", dim = list(londim, latdim), missval = fillvalue)
pval_005_AD <- ncvar_def(name = "pval_gt005", units = "no_unit", dim = list(londim, latdim), missval = fillvalue)
pval_Fit_AD_loc <- ncvar_def(name = "pval_loc", units = "no_unit", dim = list(londim, latdim), missval = fillvalue)
pval_005_AD_loc <- ncvar_def(name = "pval_loc_gt005", units = "no_unit", dim = list(londim, latdim), missval = fillvalue)

Partition <- ncvar_def(name = "partition", units = "no_unit", dim = list(londim, latdim), missval = fillvalue)

ncout <- nc_create(ncfname, list(pval_Fit_AD, pval_005_AD, pval_Fit_AD_loc, pval_005_AD_loc, Partition), force_v4=TRUE)


ncvar_put(ncout, pval_Fit_AD, matrix_regio_pvalAD)
ncvar_put(ncout, pval_005_AD, matrix_regio_pvalAD005)
ncvar_put(ncout, pval_Fit_AD_loc, matrix_local_pvalAD)
ncvar_put(ncout, pval_005_AD_loc, matrix_local_pvalAD005)
ncvar_put(ncout, Partition, matrix_cluster)

nc_close(ncout)
