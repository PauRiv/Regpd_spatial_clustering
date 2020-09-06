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
