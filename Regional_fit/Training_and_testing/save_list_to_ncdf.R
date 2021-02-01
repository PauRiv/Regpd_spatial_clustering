

seas <- "MAM"

load(paste0("/scratch3/pauline/Regpd_spatial_clustering/Regio_Fit_halfdata/res_reg_fit_only_even_OR_odd_", seas, "_1third.Rdata"))
load(paste0("/scratch3/pauline/Regpd_spatial_clustering/Regio_Fit_halfdata/keep_coord_reg_fit_only_even_OR_odd_", seas, "_1third.Rdata"))


library("ncdf4")


nc_precip_ERA5 = nc_open("/scratch3/pauline/EGPD_on_ERA5/Data/era5_daily_precip_1979-01-02_2018-12-31_Europe_land.nc")
lon_ERA5 <- ncvar_get(nc = nc_precip_ERA5, varid = "lon")
num_lon <- length(lon_ERA5)
lat_ERA5 <- ncvar_get(nc = nc_precip_ERA5, varid = "lat")
num_lat <- length(lat_ERA5)
nc_close(nc_precip_ERA5)



matrix_regio_kappa <- matrix(data = NA, nrow = num_lon, ncol = num_lat)
matrix_regio_sigma <- matrix(data = NA, nrow = num_lon, ncol = num_lat)
matrix_regio_xi <- matrix(data = NA, nrow = num_lon, ncol = num_lat)

matrix_local_kappa <- matrix(data = NA, nrow = num_lon, ncol = num_lat)
matrix_local_sigma <- matrix(data = NA, nrow = num_lon, ncol = num_lat)
matrix_local_xi <- matrix(data = NA, nrow = num_lon, ncol = num_lat)
matrix_local_xi_site <- matrix(data = NA, nrow = num_lon, ncol = num_lat)

matrix_cluster <- matrix(data = NA, nrow = num_lon, ncol = num_lat)

n_clust <- length(List_reg_fit)

for (clust in 1:n_clust) {
  n_sta <- length(List_keep_coord[[clust]]$keep_lon)
  for (sta in 1:n_sta) {
    ref_lon <- List_keep_coord[[clust]]$keep_lon[sta]
    ref_lat <- List_keep_coord[[clust]]$keep_lat[sta]
    
    matrix_regio_kappa[ref_lon,ref_lat] <- List_reg_fit[[clust]]$Theta$kappa[sta]
    matrix_regio_sigma[ref_lon,ref_lat] <- List_reg_fit[[clust]]$Theta$sigma[sta]
    matrix_regio_xi[ref_lon,ref_lat] <- List_reg_fit[[clust]]$Theta$xi.reg[sta]
    
    matrix_local_kappa[ref_lon,ref_lat] <- List_reg_fit[[clust]]$Theta_0$kappa[sta]
    matrix_local_sigma[ref_lon,ref_lat] <- List_reg_fit[[clust]]$Theta_0$sigma[sta]
    matrix_local_xi_site[ref_lon,ref_lat] <- List_reg_fit[[clust]]$Theta_0$xi.site[sta]
    matrix_local_xi[ref_lon,ref_lat] <- List_reg_fit[[clust]]$Theta_0$xi.reg[sta]
    
    matrix_cluster[ref_lon,ref_lat] <- clust
    
  }#end for sta
}#end for clust




ncfname <- paste0("/scratch3/pauline/Regpd_spatial_clustering/Regio_Fit_halfdata/regional_fitting_even_OR_odd_ERA5_EU_",seas,".nc")

# define dimensions
londim <- ncdim_def(name = "longitude", units = "degrees_east", vals = lon_ERA5)
latdim <- ncdim_def(name = "latitude", units = "degrees_north", vals = lat_ERA5)

# define variables
fillvalue <- 1e32

Regio_kappa <- ncvar_def(name = "kappa_regio", units = "no_unit", dim = list(londim, latdim), missval = fillvalue)
Regio_sigma <- ncvar_def(name = "sigma_regio", units = "no_unit", dim = list(londim, latdim), missval = fillvalue)
Regio_xi <- ncvar_def(name = "xi_regio", units = "no_unit", dim = list(londim, latdim), missval = fillvalue)

Local_kappa <- ncvar_def(name = "kappa_local", units = "no_unit", dim = list(londim, latdim), missval = fillvalue)
Local_sigma <- ncvar_def(name = "sigma_local", units = "no_unit", dim = list(londim, latdim), missval = fillvalue)
Local_xi <- ncvar_def(name = "xi_local_site", units = "no_unit", dim = list(londim, latdim), missval = fillvalue)
Local_xi_reg <- ncvar_def(name = "xi_local_reg", units = "no_unit", dim = list(londim, latdim), missval = fillvalue)

Partition <- ncvar_def(name = "partition", units = "no_unit", dim = list(londim, latdim), missval = fillvalue)

ncout <- nc_create(ncfname,
                   list(Regio_kappa, Regio_sigma, Regio_xi, Local_kappa, Local_sigma, Local_xi, Local_xi_reg, Partition),
                   force_v4=TRUE)


ncvar_put(ncout, Regio_kappa, matrix_regio_kappa)
ncvar_put(ncout, Regio_sigma, matrix_regio_sigma)
ncvar_put(ncout, Regio_xi, matrix_regio_xi)
ncvar_put(ncout, Local_kappa, matrix_local_kappa)
ncvar_put(ncout, Local_sigma, matrix_local_sigma)
ncvar_put(ncout, Local_xi, matrix_local_xi_site)
ncvar_put(ncout, Local_xi_reg, matrix_local_xi)
ncvar_put(ncout, Partition, matrix_cluster)

nc_close(ncout)
