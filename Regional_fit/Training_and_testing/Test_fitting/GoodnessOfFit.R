################################################################
#######                                                  #######
#######            Goodness of the fit tests:            #######
#######       Compare with the other half of the data    #######
#######                                                  #######
#######            Pauline Rivoire 20.01.2021            #######
#######                                                  #######
################################################################


rm(list = ls())
source("/scratch3/pauline/ExtendedGeneralizedPareto/mixtureEGPfit_April2019.R")

seas <- "SON"
# if(seas == "JJA"){print("/!\ 3 clusters in JJA")}



# Load Data ---------------------------------------------------------------

# Fitted EGPD
load(paste0("/scratch3/pauline/Regpd_spatial_clustering/Regio_Fit_halfdata/res_reg_fit_only_even_OR_odd_", seas, "_1third.Rdata"))
load(paste0("/scratch3/pauline/Regpd_spatial_clustering/Regio_Fit_halfdata/keep_coord_reg_fit_only_even_OR_odd_", seas, "_1third.Rdata"))

# Raw precip

nc_precip_ERA5 = nc_open("/scratch3/pauline/EGPD_on_ERA5/Data/era5_daily_precip_1979-01-02_2018-12-31_Europe_land.nc")

#Gridpoint corresponding to the station in ERA5
lon_ERA5 <- ncvar_get(nc = nc_precip_ERA5, varid = "lon")
num_lon <- length(lon_ERA5)
lat_ERA5 <- ncvar_get(nc = nc_precip_ERA5, varid = "lat")
num_lat <- length(lat_ERA5)
times_ERA5 <- ncvar_get(nc = nc_precip_ERA5, varid = "time")
date_ERA5 <- format(as.POSIXlt(times_ERA5*60*60, origin="1979-01-02 12:30:00"), format = "%Y-%m-%d")
precip_ERA5 <- ncvar_get(nc = nc_precip_ERA5, varid = "TP") ## Dimensions=lon*lat*time
nc_close(nc_precip_ERA5)


# List_reg_fit[[1]]$Theta$kappa[1]
# List_reg_fit[[1]]$Theta$sigma[1]
# List_reg_fit[[1]]$Theta$xi.reg[1]




library(mev)
# pextgp(q=c(3,6,9,15,17), kappa = List_reg_fit[[1]]$Theta$kappa[1],
# sigma = List_reg_fit[[1]]$Theta$sigma[1], xi = List_reg_fit[[1]]$Theta$xi.reg[1])

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



library(foreach);library(iterators);library(parallel);library(doParallel) # parallelization
registerDoParallel(cores=floor(detectCores()/3))




# AD test regional fit ------------------------------------------------

pval_test_ad_with_left <- list()

for (clust in 1:length(List_keep_coord)) {
  npoints <- length(List_reg_fit[[clust]]$Theta$xi.reg)
  pval_test_ad_with_left[[clust]] <- numeric(length = npoints)
  print(paste("cluster",clust))
  
  list_pval <- foreach(GP=1:npoints) %dopar% {
    LON <- List_keep_coord[[clust]]$keep_lon[GP]
    LAT <- List_keep_coord[[clust]]$keep_lat[GP]
    
    if(!is.na(sum(precip_ERA5[LON-1,LAT,]))){
      seas_precip_local <- extr_seas_pos_precip(season = seas, precip_timeserie = precip_ERA5[LON-1,LAT,], precip_date = date_ERA5)
    } else {
      if(!is.na(sum(precip_ERA5[LON+1,LAT,]))){
        seas_precip_local <- extr_seas_pos_precip(season = seas, precip_timeserie = precip_ERA5[LON+1,LAT,], precip_date = date_ERA5)
      } else {
        if(!is.na(sum(precip_ERA5[LON,LAT-1,]))){
          seas_precip_local <- extr_seas_pos_precip(season = seas, precip_timeserie = precip_ERA5[LON,LAT-1,], precip_date = date_ERA5)
        }else {
          if(!is.na(sum(precip_ERA5[LON,LAT+1,]))){
            seas_precip_local <- extr_seas_pos_precip(season = seas, precip_timeserie = precip_ERA5[LON,LAT+1,], precip_date = date_ERA5)
          } else {
            seas_precip_local <- NA
          }#end if else top
        }#end if else bottom
      }#end if else right
    }#and if else left
    
    
    if(!is.na(seas_precip_local)){
      ADtest <- goftest::ad.test(x=seas_precip_local[(1:length(seas_precip_local))[which((1:length(seas_precip_local)%%3)==1)]],
                                 null = pextgp, kappa = List_reg_fit[[clust]]$Theta$kappa[GP], sigma = List_reg_fit[[clust]]$Theta$sigma[GP],
                                 xi = List_reg_fit[[clust]]$Theta$xi.reg[GP])
      return(ADtest$p.value)
    } else {
      return(NA)
    }
  }#end foreach GP
  pval_test_ad_with_left[[clust]] <- unlist(list_pval)
}#end for clust

save(pval_test_ad_with_left, file = paste0("/scratch3/pauline/Regpd_spatial_clustering/Regio_Fit_halfdata/Test_fitting/test_AD_leftGP_",
                                           seas,".Rdata"))





# ADtest local fit ---------------------------------------------------
pval_test_ad_with_left_local <- list()

for (clust in 1:length(List_keep_coord)) {
  npoints <- length(List_reg_fit[[clust]]$Theta_0$xi.site)
  pval_test_ad_with_left_local[[clust]] <- numeric(length = npoints)
  print(paste("cluster",clust))
  
  list_pval <- foreach(GP=1:npoints) %dopar% {
    LON <- List_keep_coord[[clust]]$keep_lon[GP]
    LAT <- List_keep_coord[[clust]]$keep_lat[GP]
    
    if(!is.na(sum(precip_ERA5[LON-1,LAT,]))){
      seas_precip_local <- extr_seas_pos_precip(season = seas, precip_timeserie = precip_ERA5[LON-1,LAT,], precip_date = date_ERA5)
    } else {
      if(!is.na(sum(precip_ERA5[LON+1,LAT,]))){
        seas_precip_local <- extr_seas_pos_precip(season = seas, precip_timeserie = precip_ERA5[LON+1,LAT,], precip_date = date_ERA5)
      } else {
        if(!is.na(sum(precip_ERA5[LON,LAT-1,]))){
          seas_precip_local <- extr_seas_pos_precip(season = seas, precip_timeserie = precip_ERA5[LON,LAT-1,], precip_date = date_ERA5)
        }else {
          if(!is.na(sum(precip_ERA5[LON,LAT+1,]))){
            seas_precip_local <- extr_seas_pos_precip(season = seas, precip_timeserie = precip_ERA5[LON,LAT+1,], precip_date = date_ERA5)
          } else {
            seas_precip_local <- NA
          }#end if else top
        }#end if else bottom
      }#end if else right
    }#and if else left
    
    
    if(!is.na(seas_precip_local)){
      ADtest <- goftest::ad.test(x=seas_precip_local[(1:length(seas_precip_local))[which((1:length(seas_precip_local)%%3)==1)]],
                                 null = pextgp, kappa = List_reg_fit[[clust]]$Theta_0$kappa[GP], sigma = List_reg_fit[[clust]]$Theta_0$sigma[GP],
                                 xi = List_reg_fit[[clust]]$Theta_0$xi.site[GP])
      return(ADtest$p.value)
    } else {
      return(NA)
    }
  }#end foreach GP
  pval_test_ad_with_left_local[[clust]] <- unlist(list_pval)
}#end for clust

save(pval_test_ad_with_left_local, file = paste0("/scratch3/pauline/Regpd_spatial_clustering/Regio_Fit_halfdata/Test_fitting/test_AD_leftGP_",
                                                 seas,"_LocalFit.Rdata"))


# QQ-plots regio ------------------------------------------------------------------
load(paste0("/scratch3/pauline/Regpd_spatial_clustering/Compute_Clusters/list_partition_16_ERA5_EU_", seas,"_V3.Rdata"))
load(paste0("/scratch3/pauline/Regpd_spatial_clustering/Compute_Clusters/coord_in_list_",seas,".Rdata"))
status_med <- character(length = 2)
for (med in 1:length(list_partition[[1]]$medoids)) {
  
  LON <- keep_coord_list[list_partition[[1]]$medoids[med],"LON ref"]
  LAT <- keep_coord_list[list_partition[[1]]$medoids[med],"LAT ref"]
  
  if(LON%%2==LAT%%2){
    status_med <- "medoid=fitted distrib"
    coord_ref <- which(List_keep_coord[[med]]$keep_lon==LON & List_keep_coord[[med]]$keep_lat==LAT)
    
    theoritical_param <- as.matrix(c(List_reg_fit[[med]]$Theta$kappa[coord_ref], List_reg_fit[[med]]$Theta$sigma[coord_ref],
                                     List_reg_fit[[med]]$Theta$xi.reg[coord_ref]))
    theoritical_param_loc <- as.matrix(c(List_reg_fit[[cl]]$Theta_0$kappa[coord_ref], List_reg_fit[[cl]]$Theta_0$sigma[coord_ref],
                                         List_reg_fit[[cl]]$Theta_0$xi.site[coord_ref]))
    rownames(theoritical_param) <- c("kappa", "sigma", "xi");rownames(theoritical_param_loc) <- c("kappa", "sigma", "xi")
    
    
    if(!is.na(sum(precip_ERA5[LON-1,LAT,]))){
      seas_precip_local <- extr_seas_pos_precip(season = seas, precip_timeserie = precip_ERA5[LON-1,LAT,], precip_date = date_ERA5)
    } else {
      if(!is.na(sum(precip_ERA5[LON+1,LAT,]))){
        seas_precip_local <- extr_seas_pos_precip(season = seas, precip_timeserie = precip_ERA5[LON+1,LAT,], precip_date = date_ERA5)
      } else {
        if(!is.na(sum(precip_ERA5[LON,LAT-1,]))){
          seas_precip_local <- extr_seas_pos_precip(season = seas, precip_timeserie = precip_ERA5[LON,LAT-1,], precip_date = date_ERA5)
        }else {
          if(!is.na(sum(precip_ERA5[LON,LAT+1,]))){
            seas_precip_local <- extr_seas_pos_precip(season = seas, precip_timeserie = precip_ERA5[LON,LAT+1,], precip_date = date_ERA5)
          } else {
            seas_precip_local <- NA
          }#end if else top
        }#end if else bottom
      }#end if else right
    }#and if else left
    
    seas_precip1tird <- seas_precip_local[(1:length(seas_precip_local))[which((1:length(seas_precip_local)%%3)==1)]]
    
  } else {
    status_med <- "medoid=empirical distrib"
    
    seas_precip_local <- extr_seas_pos_precip(season = seas, precip_timeserie = precip_ERA5[LON,LAT,], precip_date = date_ERA5)
    seas_precip1tird <- seas_precip_local[(1:length(seas_precip_local))[which((1:length(seas_precip_local)%%3)==1)]]
    
    if(!is.na(sum(precip_ERA5[LON-1,LAT,]))){
      coord_ref <- which(List_keep_coord[[med]]$keep_lon==(LON-1) &
                           List_keep_coord[[med]]$keep_lat==(LAT))
      
      
      theoritical_param <- as.matrix(c(List_reg_fit[[med]]$Theta$kappa[coord_ref], List_reg_fit[[med]]$Theta$sigma[coord_ref],
                                     List_reg_fit[[med]]$Theta$xi.reg[coord_ref]))
      theoritical_param_loc <- as.matrix(c(List_reg_fit[[cl]]$Theta_0$kappa[coord_ref], List_reg_fit[[cl]]$Theta_0$sigma[coord_ref],
                                           List_reg_fit[[cl]]$Theta_0$xi.site[coord_ref]))
      
      rownames(theoritical_param) <- c("kappa", "sigma", "xi");rownames(theoritical_param_loc) <- c("kappa", "sigma", "xi")
      
    } else {
      if(!is.na(sum(precip_ERA5[LON+1,LAT,]))){
        coord_ref <- which(List_keep_coord[[med]]$keep_lon==(LON+1) &
                             List_keep_coord[[med]]$keep_lat==(LAT))
        
        
        theoritical_param <- as.matrix(c(List_reg_fit[[med]]$Theta$kappa[coord_ref], List_reg_fit[[med]]$Theta$sigma[coord_ref],
                                         List_reg_fit[[med]]$Theta$xi.reg[coord_ref]))
        theoritical_param_loc <- as.matrix(c(List_reg_fit[[cl]]$Theta_0$kappa[coord_ref], List_reg_fit[[cl]]$Theta_0$sigma[coord_ref],
                                             List_reg_fit[[cl]]$Theta_0$xi.site[coord_ref]))
        
        rownames(theoritical_param) <- c("kappa", "sigma", "xi");rownames(theoritical_param_loc) <- c("kappa", "sigma", "xi")
      } else {
        if(!is.na(sum(precip_ERA5[LON,LAT-1,]))){
          coord_ref <- which(List_keep_coord[[med]]$keep_lon==(LON) &
                               List_keep_coord[[med]]$keep_lat==(LAT-1))
          
          
          theoritical_param <- as.matrix(c(List_reg_fit[[med]]$Theta$kappa[coord_ref], List_reg_fit[[med]]$Theta$sigma[coord_ref],
                                           List_reg_fit[[med]]$Theta$xi.reg[coord_ref]))
          theoritical_param_loc <- as.matrix(c(List_reg_fit[[cl]]$Theta_0$kappa[coord_ref], List_reg_fit[[cl]]$Theta_0$sigma[coord_ref],
                                               List_reg_fit[[cl]]$Theta_0$xi.site[coord_ref]))
          
          rownames(theoritical_param) <- c("kappa", "sigma", "xi");rownames(theoritical_param_loc) <- c("kappa", "sigma", "xi")
        }else {
          if(!is.na(sum(precip_ERA5[LON,LAT+1,]))){
            coord_ref <- which(List_keep_coord[[med]]$keep_lon==(LON) &
                                 List_keep_coord[[med]]$keep_lat==(LAT+1))
            
            
            theoritical_param <- as.matrix(c(List_reg_fit[[med]]$Theta$kappa[coord_ref], List_reg_fit[[med]]$Theta$sigma[coord_ref],
                                             List_reg_fit[[med]]$Theta$xi.reg[coord_ref]))
            theoritical_param_loc <- as.matrix(c(List_reg_fit[[cl]]$Theta_0$kappa[coord_ref], List_reg_fit[[cl]]$Theta_0$sigma[coord_ref],
                                                 List_reg_fit[[cl]]$Theta_0$xi.site[coord_ref]))
            
            rownames(theoritical_param) <- c("kappa", "sigma", "xi");rownames(theoritical_param_loc) <- c("kappa", "sigma", "xi")
          } else {
            theoritical_param <- NA
          }#end if else top
        }#end if else bottom
      }#end if else right
    }#and if else left
    
  }#end if else LON LAT both even OR odd
  
  proba_quantiles=c(seq(0.1,0.8,by = 0.1), seq(0.81,0.95, by = 0.01), seq(0.951,0.99, by = 0.005))
  
  Q_emp <- quantile(ecdf(seas_precip1tird), probs=proba_quantiles)
  Q_th <- qEGP(x = proba_quantiles, param = c(theoritical_param["xi",], theoritical_param["sigma",], theoritical_param["kappa",]))
  Q_th_loc <- qEGP(x = proba_quantiles, param = c(theoritical_param_loc["xi",], theoritical_param_loc["sigma",], theoritical_param_loc["kappa",]))
  
  jpeg(filename = paste0("/scratch3/pauline/Regpd_spatial_clustering/Regio_Fit_halfdata/Test_fitting/qqplot_medoid_clust_",med,"_",seas,".jpg"),
       width = 500,height = 500)
  plot(Q_emp, Q_th, main=paste0(seas,": QQplot medoid cluster ", med, "\n",status_med),
       xlab="Empirical quantile [mm/wetday]", ylab="Theoritical quantile [mm/wetday]", pch=18, col="gray", cex=2)
  points(Q_emp, Q_th_loc, pch="+",cex=1.2,col="black")
  abline(b=1, a=0, lty=2)
  legend("topleft", pch=c(18, 3), col = c("gray", "black"), legend = c("Regional fitting", "Local fitting"))
  dev.off()
  
}#end for med


# hist(list_partition[[1]]$silinfo$widths[,"sil_width"])
# min(list_partition[[1]]$silinfo$widths[,"sil_width"])
# max(list_partition[[1]]$silinfo$widths[,"sil_width"])
# 
# rownames(list_partition[[1]]$silinfo$widths)[1]
# rownames(list_partition[[1]]$silinfo$widths)[length(list_partition[[1]]$silinfo$widths[,"sil_width"])]


SIL_CLUS_1 <- as.matrix(list_partition[[1]]$silinfo$widths[as.numeric(which(list_partition[[1]]$silinfo$widths[,"cluster"]==1)),"sil_width"])
SIL_CLUS_2 <- as.matrix(list_partition[[1]]$silinfo$widths[as.numeric(which(list_partition[[1]]$silinfo$widths[,"cluster"]==2)),"sil_width"])

max_sil_1 <- row.names(SIL_CLUS_1)[1]
#list_partition[[1]]$silinfo$widths[rownames(list_partition[[1]]$silinfo$widths)==max_sil_1,"sil_width"]
min_sil_1 <- row.names(SIL_CLUS_1)[nrow(SIL_CLUS_1)]
if(!length(which(List_keep_coord[[1]]$keep_lon==keep_coord_list[as.numeric(min_sil_1),"LON ref"] &
                 List_keep_coord[[1]]$keep_lat==keep_coord_list[as.numeric(min_sil_1),"LAT ref"]))){
  min_sil_1 <- row.names(SIL_CLUS_1)[nrow(SIL_CLUS_1)-1]
}
#list_partition[[1]]$silinfo$widths[rownames(list_partition[[1]]$silinfo$widths)==max_sil_2,"sil_width"]
max_sil_2 <- row.names(SIL_CLUS_2)[1]
#list_partition[[1]]$silinfo$widths[rownames(list_partition[[1]]$silinfo$widths)==min_sil_1,"sil_width"]
min_sil_2 <- row.names(SIL_CLUS_2)[nrow(SIL_CLUS_2)]
#list_partition[[1]]$silinfo$widths[rownames(list_partition[[1]]$silinfo$widths)==min_sil_2,"sil_width"]

point_names <- c("max silh cl.1","min silh cl.1","max silh cl.2","min silh cl.2")
for (GP in point_names) {
  
  if(GP=="max silh cl.1"){GP_ref <- as.numeric(max_sil_1); id=1; cl=1; filname="maxsilh_clus1"}
  if(GP=="min silh cl.1"){GP_ref <- as.numeric(min_sil_1); id=2; cl=1; filname="minsilh_clus1"}
  if(GP=="max silh cl.2"){GP_ref <- as.numeric(max_sil_2); id=3; cl=2; filname="maxsilh_clus2"}
  if(GP=="min silh cl.2"){GP_ref <- as.numeric(min_sil_2); id=4; cl=2; filname="minsilh_clus2"}
  
  LON <- keep_coord_list[GP_ref,"LON ref"]
  LAT <- keep_coord_list[GP_ref,"LAT ref"]
  
  status_GP <- character()
  
  if(LON%%2==LAT%%2){
    status_GP <- "GP=fitted distrib"
    coord_ref <- which(List_keep_coord[[cl]]$keep_lon==LON & List_keep_coord[[cl]]$keep_lat==LAT)
    
    theoritical_param <- as.matrix(c(List_reg_fit[[cl]]$Theta$kappa[coord_ref], List_reg_fit[[cl]]$Theta$sigma[coord_ref],
                                     List_reg_fit[[cl]]$Theta$xi.reg[coord_ref]))
    theoritical_param_loc <- as.matrix(c(List_reg_fit[[cl]]$Theta_0$kappa[coord_ref], List_reg_fit[[cl]]$Theta_0$sigma[coord_ref],
                                         List_reg_fit[[cl]]$Theta_0$xi.site[coord_ref]))
    
    if(length(theoritical_param>0 & length(theoritical_param_loc)>0)){
      rownames(theoritical_param) <- c("kappa", "sigma", "xi");rownames(theoritical_param_loc) <- c("kappa", "sigma", "xi")}
    
    if(!is.na(sum(precip_ERA5[LON-1,LAT,]))){
      seas_precip_local <- extr_seas_pos_precip(season = seas, precip_timeserie = precip_ERA5[LON-1,LAT,], precip_date = date_ERA5)
    } else {
      if(!is.na(sum(precip_ERA5[LON+1,LAT,]))){
        seas_precip_local <- extr_seas_pos_precip(season = seas, precip_timeserie = precip_ERA5[LON+1,LAT,], precip_date = date_ERA5)
      } else {
        if(!is.na(sum(precip_ERA5[LON,LAT-1,]))){
          seas_precip_local <- extr_seas_pos_precip(season = seas, precip_timeserie = precip_ERA5[LON,LAT-1,], precip_date = date_ERA5)
        }else {
          if(!is.na(sum(precip_ERA5[LON,LAT+1,]))){
            seas_precip_local <- extr_seas_pos_precip(season = seas, precip_timeserie = precip_ERA5[LON,LAT+1,], precip_date = date_ERA5)
          } else {
            seas_precip_local <- NA
          }#end if else top
        }#end if else bottom
      }#end if else right
    }#and if else left
    
    seas_precip1tird <- seas_precip_local[(1:length(seas_precip_local))[which((1:length(seas_precip_local)%%3)==1)]]
    
  } else {
    status_GP <- "GP=empirical distrib"
    
    seas_precip_local <- extr_seas_pos_precip(season = seas, precip_timeserie = precip_ERA5[LON,LAT,], precip_date = date_ERA5)
    seas_precip1tird <- seas_precip_local[(1:length(seas_precip_local))[which((1:length(seas_precip_local)%%3)==1)]]
    
    if(length(which(List_keep_coord[[cl]]$keep_lon==(LON-1) & List_keep_coord[[cl]]$keep_lat==(LAT)))){
      coord_ref <- which(List_keep_coord[[cl]]$keep_lon==(LON-1) &
                           List_keep_coord[[cl]]$keep_lat==(LAT))
      
      
      theoritical_param <- as.matrix(c(List_reg_fit[[cl]]$Theta$kappa[coord_ref], List_reg_fit[[cl]]$Theta$sigma[coord_ref],
                                       List_reg_fit[[cl]]$Theta$xi.reg[coord_ref]))
      theoritical_param_loc <- as.matrix(c(List_reg_fit[[cl]]$Theta_0$kappa[coord_ref], List_reg_fit[[cl]]$Theta_0$sigma[coord_ref],
                                           List_reg_fit[[cl]]$Theta_0$xi.site[coord_ref]))
      
      rownames(theoritical_param) <- c("kappa", "sigma", "xi")
      rownames(theoritical_param_loc) <- c("kappa", "sigma", "xi")
      
    } else {
      if(length(which(List_keep_coord[[cl]]$keep_lon==(LON+1) & List_keep_coord[[cl]]$keep_lat==(LAT)))){
        coord_ref <- which(List_keep_coord[[cl]]$keep_lon==(LON+1) &
                             List_keep_coord[[cl]]$keep_lat==(LAT))
        
        
        theoritical_param <- as.matrix(c(List_reg_fit[[cl]]$Theta$kappa[coord_ref], List_reg_fit[[cl]]$Theta$sigma[coord_ref],
                                         List_reg_fit[[cl]]$Theta$xi.reg[coord_ref]))
        theoritical_param_loc <- as.matrix(c(List_reg_fit[[cl]]$Theta_0$kappa[coord_ref], List_reg_fit[[cl]]$Theta_0$sigma[coord_ref],
                                             List_reg_fit[[cl]]$Theta_0$xi.site[coord_ref]))
        
        rownames(theoritical_param) <- c("kappa", "sigma", "xi")
        rownames(theoritical_param_loc) <- c("kappa", "sigma", "xi")
      } else {
        if(length(which(List_keep_coord[[cl]]$keep_lon==(LON) & List_keep_coord[[cl]]$keep_lat==(LAT-1)))){
          coord_ref <- which(List_keep_coord[[cl]]$keep_lon==(LON) &
                               List_keep_coord[[cl]]$keep_lat==(LAT-1))
          
          
          theoritical_param <- as.matrix(c(List_reg_fit[[cl]]$Theta$kappa[coord_ref], List_reg_fit[[cl]]$Theta$sigma[coord_ref],
                                           List_reg_fit[[cl]]$Theta$xi.reg[coord_ref]))
          theoritical_param_loc <- as.matrix(c(List_reg_fit[[cl]]$Theta_0$kappa[coord_ref], List_reg_fit[[cl]]$Theta_0$sigma[coord_ref],
                                               List_reg_fit[[cl]]$Theta_0$xi.site[coord_ref]))
          
          rownames(theoritical_param) <- c("kappa", "sigma", "xi")
          rownames(theoritical_param_loc) <- c("kappa", "sigma", "xi")
        }else {
          if(length(which(List_keep_coord[[cl]]$keep_lon==(LON) & List_keep_coord[[cl]]$keep_lat==(LAT+1)))){
            coord_ref <- which(List_keep_coord[[cl]]$keep_lon==(LON) &
                                 List_keep_coord[[cl]]$keep_lat==(LAT+1))
            
            
            theoritical_param <- as.matrix(c(List_reg_fit[[cl]]$Theta$kappa[coord_ref], List_reg_fit[[cl]]$Theta$sigma[coord_ref],
                                             List_reg_fit[[cl]]$Theta$xi.reg[coord_ref]))
            theoritical_param_loc <- as.matrix(c(List_reg_fit[[cl]]$Theta_0$kappa[coord_ref], List_reg_fit[[cl]]$Theta_0$sigma[coord_ref],
                                                 List_reg_fit[[cl]]$Theta_0$xi.site[coord_ref]))
            
            rownames(theoritical_param) <- c("kappa", "sigma", "xi")
            rownames(theoritical_param_loc) <- c("kappa", "sigma", "xi")
          } else {
            theoritical_param <- NA
          }#end if else top
        }#end if else bottom
      }#end if else right
    }#and if else left
    
  }#end if else LON LAT both even OR odd
  if(length(theoritical_param)>0 & length(theoritical_param_loc)>0){
    proba_quantiles=c(seq(0.1,0.8,by = 0.1), seq(0.81,0.95, by = 0.01), seq(0.951,0.99, by = 0.005))
  
    Q_emp <- quantile(ecdf(seas_precip1tird), probs=proba_quantiles)
    Q_th <- qEGP(x = proba_quantiles, param = c(theoritical_param["xi",], theoritical_param["sigma",], theoritical_param["kappa",]))
    Q_th_loc <- qEGP(x = proba_quantiles, param = c(theoritical_param_loc["xi",], theoritical_param_loc["sigma",], theoritical_param_loc["kappa",]))
    
    jpeg(filename = paste0("/scratch3/pauline/Regpd_spatial_clustering/Regio_Fit_halfdata/Test_fitting/qqplot",filname,"_",seas,".jpg"),
        width = 500,height = 500)
    plot(Q_emp, Q_th, main=paste0(seas,": QQplot ", GP, "\n",status_GP),
        xlab="Empirical quantile [mm/wetday]", ylab="Theoritical quantile [mm/wetday]", pch=18, col="gray", cex=2)
    points(Q_emp, Q_th_loc, pch="+",cex=1.2,col="black")
    abline(b=1, a=0, lty=2)
    legend("topleft", pch=c(18, 3), col = c("gray", "black"), legend = c("Regional fitting", "Local fitting"))
    dev.off()
  }#end if pb at GP
  
}#end for GP

