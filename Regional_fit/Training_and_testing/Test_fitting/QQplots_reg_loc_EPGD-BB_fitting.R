

library(mev);library(ncdf4)
source("/scratch3/pauline/ExtendedGeneralizedPareto/mixtureEGPfit_April2019.R")


nc_precip_ERA5 = nc_open("/scratch3/pauline/EGPD_on_ERA5/Data/era5_daily_precip_1979-01-02_2018-12-31_Europe_land.nc")
lon_ERA5 <- ncvar_get(nc = nc_precip_ERA5, varid = "lon") ; num_lon <- length(lon_ERA5)
lat_ERA5 <- ncvar_get(nc = nc_precip_ERA5, varid = "lat") ; num_lat <- length(lat_ERA5)
times_ERA5 <- ncvar_get(nc = nc_precip_ERA5, varid = "time")
date_ERA5 <- format(as.POSIXlt(times_ERA5*60*60, origin="1979-01-02 12:30:00"), format = "%Y-%m-%d")
precip_ERA5 <- ncvar_get(nc = nc_precip_ERA5, varid = "TP") ## Dimensions=lon*lat*time
nc_close(nc_precip_ERA5)

seas <- "DJF"

# Fitted EGPD.kappa
load(paste0("/scratch3/pauline/Regpd_spatial_clustering/Regio_Fit_halfdata/res_reg_fit_only_even_OR_odd_", seas, "_1third.Rdata"))
load(paste0("/scratch3/pauline/Regpd_spatial_clustering/Regio_Fit_halfdata/keep_coord_reg_fit_only_even_OR_odd_", seas, "_1third.Rdata"))

# Fitted EGPD.bb
nc_ERAfit <- nc_open(filename = paste0("/scratch3/pauline/Bootstrap_nonparametric_EGPD/ERA5_Europe/fit_1thirddata_EGPD_parameters_ERA5_m30_B100_",
                                       seas,"_1979_2018_wetdays1mm_minnb500.nc"))
w1 <- ncvar_get(nc=nc_ERAfit,varid = "w1");w2 <- ncvar_get(nc=nc_ERAfit,varid = "w2");
w3 <- ncvar_get(nc=nc_ERAfit,varid = "w3");w4 <- ncvar_get(nc=nc_ERAfit,varid = "w4");
w5 <- ncvar_get(nc=nc_ERAfit,varid = "w5");w6 <- ncvar_get(nc=nc_ERAfit,varid = "w6");
w7 <- ncvar_get(nc=nc_ERAfit,varid = "w7");w8 <- ncvar_get(nc=nc_ERAfit,varid = "w8");
w9 <- ncvar_get(nc=nc_ERAfit,varid = "w9");w10 <- ncvar_get(nc=nc_ERAfit,varid = "w10");
w11 <- ncvar_get(nc=nc_ERAfit,varid = "w11");w12 <- ncvar_get(nc=nc_ERAfit,varid = "w12");
w13 <- ncvar_get(nc=nc_ERAfit,varid = "w13");w14 <- ncvar_get(nc=nc_ERAfit,varid = "w14");
w15 <- ncvar_get(nc=nc_ERAfit,varid = "w15");w16 <- ncvar_get(nc=nc_ERAfit,varid = "w16");
w17 <- ncvar_get(nc=nc_ERAfit,varid = "w17");w18 <- ncvar_get(nc=nc_ERAfit,varid = "w18");
w19 <- ncvar_get(nc=nc_ERAfit,varid = "w19");w20 <- ncvar_get(nc=nc_ERAfit,varid = "w20");
w21 <- ncvar_get(nc=nc_ERAfit,varid = "w21");w22 <- ncvar_get(nc=nc_ERAfit,varid = "w22");
w23 <- ncvar_get(nc=nc_ERAfit,varid = "w23");w24 <- ncvar_get(nc=nc_ERAfit,varid = "w24");
w25 <- ncvar_get(nc=nc_ERAfit,varid = "w25");w26 <- ncvar_get(nc=nc_ERAfit,varid = "w26");
w27 <- ncvar_get(nc=nc_ERAfit,varid = "w27");w28 <- ncvar_get(nc=nc_ERAfit,varid = "w28");
w29 <- ncvar_get(nc=nc_ERAfit,varid = "w29");w30 <- ncvar_get(nc=nc_ERAfit,varid = "w30");
xi <- ncvar_get(nc=nc_ERAfit,varid = "xi");sigma <- ncvar_get(nc=nc_ERAfit,varid = "sigma")
nc_close(nc_ERAfit)


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
  SEASONS <- c("MAM", "JJA", "SON", "DJF")
  stopifnot((length(precip_timeserie)==length(precip_date)) & (season %in% SEASONS))
  ref_seas <- which(SEASONS==season)
  months_in_season <- (((c(1,2,3)+(3*ref_seas-1))-1)%%12)+1 #gives for ex: c(12,1,2) if season == "DJF"
  
  return(precip_timeserie[(lubridate::month(precip_date)==months_in_season[1]
                           | lubridate::month(precip_date)==months_in_season[2]
                           | lubridate::month(precip_date)==months_in_season[3]) & 
                            (precip_timeserie>thshld)])
}#end for keep_seasonal_positive_precip function


# QQ-plots medoid ------------------------------------------------------------------
load(paste0("/scratch3/pauline/Regpd_spatial_clustering/Compute_Clusters/list_partition_16_ERA5_EU_", seas,"_V3.Rdata"))
load(paste0("/scratch3/pauline/Regpd_spatial_clustering/Compute_Clusters/coord_in_list_",seas,".Rdata"))

for (med in 1:length(list_partition[[1]]$medoids)) {
  status_med <- character()
  LON <- keep_coord_list[list_partition[[1]]$medoids[med],"LON ref"]
  LAT <- keep_coord_list[list_partition[[1]]$medoids[med],"LAT ref"]
  
  if(LON%%2==LAT%%2){
    LON_bb <- LON;LAT_bb <- LAT
    status_med <- "medoid=fitted distrib"
    coord_ref <- which(List_keep_coord[[med]]$keep_lon==LON & List_keep_coord[[med]]$keep_lat==LAT)
    
    nwet_days_fit <- floor(length(extr_seas_pos_precip(season = seas, precip_timeserie = precip_ERA5[LON,LAT,], precip_date = date_ERA5))/3)
    
    theoritical_param <- as.matrix(c(List_reg_fit[[med]]$Theta$kappa[coord_ref], List_reg_fit[[med]]$Theta$sigma[coord_ref],
                                     List_reg_fit[[med]]$Theta$xi.reg[coord_ref]))
    theoritical_param_loc <- as.matrix(c(List_reg_fit[[med]]$Theta_0$kappa[coord_ref], List_reg_fit[[med]]$Theta_0$sigma[coord_ref],
                                         List_reg_fit[[med]]$Theta_0$xi.site[coord_ref]))
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
    nwet_days_fit <- length(seas_precip1tird)
    if(!is.na(sum(precip_ERA5[LON-1,LAT,]))){
      LON_bb <- LON - 1
      LAT_bb <- LAT
      coord_ref <- which(List_keep_coord[[med]]$keep_lon==(LON-1) & List_keep_coord[[med]]$keep_lat==(LAT))
      
      
      theoritical_param <- as.matrix(c(List_reg_fit[[med]]$Theta$kappa[coord_ref], List_reg_fit[[med]]$Theta$sigma[coord_ref],
                                       List_reg_fit[[med]]$Theta$xi.reg[coord_ref]))
      theoritical_param_loc <- as.matrix(c(List_reg_fit[[med]]$Theta_0$kappa[coord_ref], List_reg_fit[[med]]$Theta_0$sigma[coord_ref],
                                           List_reg_fit[[med]]$Theta_0$xi.site[coord_ref]))
      
      rownames(theoritical_param) <- c("kappa", "sigma", "xi");rownames(theoritical_param_loc) <- c("kappa", "sigma", "xi")
      
    } else {
      if(!is.na(sum(precip_ERA5[LON+1,LAT,]))){
        LON_bb <- LON + 1
        LAT_bb <- LAT
        coord_ref <- which(List_keep_coord[[med]]$keep_lon==(LON+1) & List_keep_coord[[med]]$keep_lat==(LAT))
        
        
        theoritical_param <- as.matrix(c(List_reg_fit[[med]]$Theta$kappa[coord_ref], List_reg_fit[[med]]$Theta$sigma[coord_ref],
                                         List_reg_fit[[med]]$Theta$xi.reg[coord_ref]))
        theoritical_param_loc <- as.matrix(c(List_reg_fit[[med]]$Theta_0$kappa[coord_ref], List_reg_fit[[med]]$Theta_0$sigma[coord_ref],
                                             List_reg_fit[[med]]$Theta_0$xi.site[coord_ref]))
        
        rownames(theoritical_param) <- c("kappa", "sigma", "xi");rownames(theoritical_param_loc) <- c("kappa", "sigma", "xi")
      } else {
        if(!is.na(sum(precip_ERA5[LON,LAT-1,]))){
          LON_bb <- LON
          LAT_bb <- LAT - 1
          coord_ref <- which(List_keep_coord[[med]]$keep_lon==(LON) & List_keep_coord[[med]]$keep_lat==(LAT-1))
          
          
          theoritical_param <- as.matrix(c(List_reg_fit[[med]]$Theta$kappa[coord_ref], List_reg_fit[[med]]$Theta$sigma[coord_ref],
                                           List_reg_fit[[med]]$Theta$xi.reg[coord_ref]))
          theoritical_param_loc <- as.matrix(c(List_reg_fit[[med]]$Theta_0$kappa[coord_ref], List_reg_fit[[med]]$Theta_0$sigma[coord_ref],
                                               List_reg_fit[[med]]$Theta_0$xi.site[coord_ref]))
          
          rownames(theoritical_param) <- c("kappa", "sigma", "xi");rownames(theoritical_param_loc) <- c("kappa", "sigma", "xi")
        }else {
          if(!is.na(sum(precip_ERA5[LON,LAT+1,]))){
            LON_bb <- LON
            LAT_bb <- LAT + 1
            coord_ref <- which(List_keep_coord[[med]]$keep_lon==(LON) & List_keep_coord[[med]]$keep_lat==(LAT+1))
            
            
            theoritical_param <- as.matrix(c(List_reg_fit[[med]]$Theta$kappa[coord_ref], List_reg_fit[[med]]$Theta$sigma[coord_ref],
                                             List_reg_fit[[med]]$Theta$xi.reg[coord_ref]))
            theoritical_param_loc <- as.matrix(c(List_reg_fit[[med]]$Theta_0$kappa[coord_ref], List_reg_fit[[med]]$Theta_0$sigma[coord_ref],
                                                 List_reg_fit[[med]]$Theta_0$xi.site[coord_ref]))
            
            rownames(theoritical_param) <- c("kappa", "sigma", "xi");rownames(theoritical_param_loc) <- c("kappa", "sigma", "xi")
          } else {
            theoritical_param <- NA
          }#end if else top
        }#end if else bottom
      }#end if else right
    }#and if else left
    
  }#end if else LON LAT both even OR odd
  
  qemp2 <- sort(seas_precip1tird)
  p=seq(from=1,to=length(qemp2)-0.1,length.out=length(qemp2))/length(qemp2)
  Q_th <- qextgp(p = p, kappa = theoritical_param["kappa",], sigma = theoritical_param["sigma",], xi = theoritical_param["xi",])
  Q_th_loc <- qextgp(p = p, kappa = theoritical_param_loc["kappa",], sigma = theoritical_param_loc["sigma",], xi = theoritical_param_loc["xi",])
  
  qempbb <- sort(seas_precip1tird[seas_precip1tird>1])
  pB=seq(from=1,to=length(qempbb)-0.1,length.out=length(qempbb))/length(qempbb)
  if(!is.na(w1[LON_bb,LAT_bb])){
    Q_BB <- qEGP.BB(x = pB, param=c(w1[LON_bb,LAT_bb],w2[LON_bb,LAT_bb],w3[LON_bb,LAT_bb],
                                    w4[LON_bb,LAT_bb],w5[LON_bb,LAT_bb],w6[LON_bb,LAT_bb],
                                    w7[LON_bb,LAT_bb],w8[LON_bb,LAT_bb],w9[LON_bb,LAT_bb],
                                    w10[LON_bb,LAT_bb],w11[LON_bb,LAT_bb],w12[LON_bb,LAT_bb],
                                    w13[LON_bb,LAT_bb],w14[LON_bb,LAT_bb],w15[LON_bb,LAT_bb],
                                    w16[LON_bb,LAT_bb],w17[LON_bb,LAT_bb],w18[LON_bb,LAT_bb],
                                    w19[LON_bb,LAT_bb],w20[LON_bb,LAT_bb],w21[LON_bb,LAT_bb],
                                    w22[LON_bb,LAT_bb],w23[LON_bb,LAT_bb],w24[LON_bb,LAT_bb],
                                    w25[LON_bb,LAT_bb],w26[LON_bb,LAT_bb],w27[LON_bb,LAT_bb],
                                    w28[LON_bb,LAT_bb],w29[LON_bb,LAT_bb],w30[LON_bb,LAT_bb],
                                    sigma[LON_bb,LAT_bb],xi[LON_bb,LAT_bb])) + 1
    max_x_y <- max(c(Q_th, Q_th_loc,Q_BB))} else {max_x_y <- max(c(Q_th, Q_th_loc))}
  jpeg(filename = paste0("/scratch3/pauline/Regpd_spatial_clustering/Regio_Fit_halfdata/Test_fitting/QQ_withBB/qqplot_medoid_clust_",
                         med,"_reg-loc-BB_",seas,".jpg"), width = 500,height = 500, quality = 85)
  ref_pos_precip <- which(qemp2>1)
  plot(qemp2[ref_pos_precip], Q_th[ref_pos_precip], main=paste0(seas,": QQplot medoid cluster ", med, "\n",status_med, "; nb. fitted wet days=", nwet_days_fit),
       xlim=c(0,max_x_y), ylim=c(0,max_x_y), xlab="Empirical quantile [mm/wetday]", ylab="Theoritical quantile [mm/wetday]", pch=18, col="gray", cex=2)
  if(!is.na(w1[LON_bb,LAT_bb])){points(qempbb, Q_BB, pch=4,cex=1.2,col="deepskyblue2")
  } else {text(x=2,y=max_x_y-1, labels = "BB not fitted", col = "deepskyblue2")}
  points(qemp2[ref_pos_precip], Q_th_loc[ref_pos_precip], pch="+",cex=1.2,col="black")
  abline(b=1, a=0, lty=2)
  legend("bottomright", pch=c(18, 3,4), col = c("gray", "black","deepskyblue2"), legend = c("Regional fitting", "Local fitting","Local BB fitting"))
  dev.off()
  
}#end for med



# QQ-plots maxmin silh regio ----------------------------------------------
load(paste0("/scratch3/pauline/Regpd_spatial_clustering/Compute_Clusters/list_partition_16_ERA5_EU_", seas,"_V3.Rdata"))
load(paste0("/scratch3/pauline/Regpd_spatial_clustering/Compute_Clusters/coord_in_list_",seas,".Rdata"))

SIL_CLUS_1 <- as.matrix(list_partition[[1]]$silinfo$widths[as.numeric(which(list_partition[[1]]$silinfo$widths[,"cluster"]==1)),"sil_width"])
SIL_CLUS_2 <- as.matrix(list_partition[[1]]$silinfo$widths[as.numeric(which(list_partition[[1]]$silinfo$widths[,"cluster"]==2)),"sil_width"])

max_sil_1 <- row.names(SIL_CLUS_1)[1]
min_sil_1 <- row.names(SIL_CLUS_1)[nrow(SIL_CLUS_1)]
if(!length(which(List_keep_coord[[1]]$keep_lon==keep_coord_list[as.numeric(min_sil_1),"LON ref"] &
                 List_keep_coord[[1]]$keep_lat==keep_coord_list[as.numeric(min_sil_1),"LAT ref"]))){
  min_sil_1 <- row.names(SIL_CLUS_1)[nrow(SIL_CLUS_1)-1]
}
max_sil_2 <- row.names(SIL_CLUS_2)[1]
min_sil_2 <- row.names(SIL_CLUS_2)[nrow(SIL_CLUS_2)]

point_names <- c("max silh cl.1","min silh cl.1","max silh cl.2","min silh cl.2")
for (GP in point_names) {
  
  if(GP=="max silh cl.1"){GP_ref <- as.numeric(max_sil_1); id=1; cl=1; filname="maxsilh_clus1"}
  if(GP=="min silh cl.1"){GP_ref <- as.numeric(min_sil_1); id=2; cl=1; filname="minsilh_clus1"}
  if(GP=="max silh cl.2"){GP_ref <- as.numeric(max_sil_2); id=3; cl=2; filname="maxsilh_clus2"}
  if(GP=="min silh cl.2"){GP_ref <- as.numeric(min_sil_2); id=4; cl=2; filname="minsilh_clus2"}
  
  LON <- keep_coord_list[GP_ref,"LON ref"] ; LAT <- keep_coord_list[GP_ref,"LAT ref"]
  nwet_days_fit <- floor(length(extr_seas_pos_precip(season = seas, precip_timeserie = precip_ERA5[LON,LAT,], precip_date = date_ERA5))/3)
  status_GP <- character()
  
  if(LON%%2==LAT%%2){
    LON_bb <- LON ; LAT_bb <- LAT
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
    nwet_days_fit <- length(seas_precip1tird)
    if(length(which(List_keep_coord[[cl]]$keep_lon==(LON-1) & List_keep_coord[[cl]]$keep_lat==(LAT)))){
      LON_bb <- LON - 1 ; LAT_bb <- LAT
      coord_ref <- which(List_keep_coord[[cl]]$keep_lon==(LON-1) & List_keep_coord[[cl]]$keep_lat==(LAT))
      
      theoritical_param <- as.matrix(c(List_reg_fit[[cl]]$Theta$kappa[coord_ref], List_reg_fit[[cl]]$Theta$sigma[coord_ref],
                                       List_reg_fit[[cl]]$Theta$xi.reg[coord_ref]))
      theoritical_param_loc <- as.matrix(c(List_reg_fit[[cl]]$Theta_0$kappa[coord_ref], List_reg_fit[[cl]]$Theta_0$sigma[coord_ref],
                                           List_reg_fit[[cl]]$Theta_0$xi.site[coord_ref]))
      
      rownames(theoritical_param) <- c("kappa", "sigma", "xi")
      rownames(theoritical_param_loc) <- c("kappa", "sigma", "xi")
      
    } else {
      if(length(which(List_keep_coord[[cl]]$keep_lon==(LON+1) & List_keep_coord[[cl]]$keep_lat==(LAT)))){
        LON_bb <- LON + 1 ; LAT_bb <- LAT
        coord_ref <- which(List_keep_coord[[cl]]$keep_lon==(LON+1) & List_keep_coord[[cl]]$keep_lat==(LAT))
        
        theoritical_param <- as.matrix(c(List_reg_fit[[cl]]$Theta$kappa[coord_ref], List_reg_fit[[cl]]$Theta$sigma[coord_ref],
                                         List_reg_fit[[cl]]$Theta$xi.reg[coord_ref]))
        theoritical_param_loc <- as.matrix(c(List_reg_fit[[cl]]$Theta_0$kappa[coord_ref], List_reg_fit[[cl]]$Theta_0$sigma[coord_ref],
                                             List_reg_fit[[cl]]$Theta_0$xi.site[coord_ref]))
        
        rownames(theoritical_param) <- c("kappa", "sigma", "xi")
        rownames(theoritical_param_loc) <- c("kappa", "sigma", "xi")
      } else {
        if(length(which(List_keep_coord[[cl]]$keep_lon==(LON) & List_keep_coord[[cl]]$keep_lat==(LAT-1)))){
          LON_bb <- LON ; LAT_bb <- LAT - 1
          coord_ref <- which(List_keep_coord[[cl]]$keep_lon==(LON) & List_keep_coord[[cl]]$keep_lat==(LAT-1))
          
          theoritical_param <- as.matrix(c(List_reg_fit[[cl]]$Theta$kappa[coord_ref], List_reg_fit[[cl]]$Theta$sigma[coord_ref],
                                           List_reg_fit[[cl]]$Theta$xi.reg[coord_ref]))
          theoritical_param_loc <- as.matrix(c(List_reg_fit[[cl]]$Theta_0$kappa[coord_ref], List_reg_fit[[cl]]$Theta_0$sigma[coord_ref],
                                               List_reg_fit[[cl]]$Theta_0$xi.site[coord_ref]))
          
          rownames(theoritical_param) <- c("kappa", "sigma", "xi")
          rownames(theoritical_param_loc) <- c("kappa", "sigma", "xi")
        }else {
          if(length(which(List_keep_coord[[cl]]$keep_lon==(LON) & List_keep_coord[[cl]]$keep_lat==(LAT+1)))){
            LON_bb <- LON ; LAT_bb <- LAT + 1
            coord_ref <- which(List_keep_coord[[cl]]$keep_lon==(LON) & List_keep_coord[[cl]]$keep_lat==(LAT+1))
            
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
    
    qemp2 <- sort(seas_precip1tird)
    p=seq(from=1,to=length(qemp2)-0.1,length.out=length(qemp2))/length(qemp2)
    Q_th <- qextgp(p = p, kappa = theoritical_param["kappa",], sigma = theoritical_param["sigma",], xi = theoritical_param["xi",])
    Q_th_loc <- qextgp(p = p, kappa = theoritical_param_loc["kappa",], sigma = theoritical_param_loc["sigma",], xi = theoritical_param_loc["xi",])
    
    qempbb <- sort(seas_precip1tird[seas_precip1tird>1])
    pB=seq(from=1,to=length(qempbb)-0.1,length.out=length(qempbb))/length(qempbb)
    Q_BB <- qEGP.BB(x = pB, param=c(w1[LON_bb,LAT_bb],w2[LON_bb,LAT_bb],w3[LON_bb,LAT_bb],
                                    w4[LON_bb,LAT_bb],w5[LON_bb,LAT_bb],w6[LON_bb,LAT_bb],
                                    w7[LON_bb,LAT_bb],w8[LON_bb,LAT_bb],w9[LON_bb,LAT_bb],
                                    w10[LON_bb,LAT_bb],w11[LON_bb,LAT_bb],w12[LON_bb,LAT_bb],
                                    w13[LON_bb,LAT_bb],w14[LON_bb,LAT_bb],w15[LON_bb,LAT_bb],
                                    w16[LON_bb,LAT_bb],w17[LON_bb,LAT_bb],w18[LON_bb,LAT_bb],
                                    w19[LON_bb,LAT_bb],w20[LON_bb,LAT_bb],w21[LON_bb,LAT_bb],
                                    w22[LON_bb,LAT_bb],w23[LON_bb,LAT_bb],w24[LON_bb,LAT_bb],
                                    w25[LON_bb,LAT_bb],w26[LON_bb,LAT_bb],w27[LON_bb,LAT_bb],
                                    w28[LON_bb,LAT_bb],w29[LON_bb,LAT_bb],w30[LON_bb,LAT_bb],
                                    sigma[LON_bb,LAT_bb],xi[LON_bb,LAT_bb])) + 1
    max_x_y <- max(c(Q_th, Q_th_loc,Q_BB))
    
    jpeg(filename = paste0("/scratch3/pauline/Regpd_spatial_clustering/Regio_Fit_halfdata/Test_fitting/QQ_withBB/qqplot",filname,"_reg-loc-BB_",seas,".jpg"),
         width = 500,height = 500, quality = 85)
    ref_pos_precip <- which(qemp2>1)
    plot(qemp2[ref_pos_precip], Q_th[ref_pos_precip], main=paste0(seas,": QQplot ", GP, "\n",status_GP, "; nb. fitted wet days=", nwet_days_fit),
         xlim = c(0,max_x_y), ylim=c(0,max_x_y), xlab="Empirical quantile [mm/wetday]", ylab="Theoritical quantile [mm/wetday]", pch=18, col="gray", cex=2)
    points(qempbb, Q_BB, pch="+",cex=1.2,col="deepskyblue2")
    points(qemp2[ref_pos_precip], Q_th_loc[ref_pos_precip], pch="+",cex=1.2,col="black")
    abline(b=1, a=0, lty=2)
    legend("bottomright", pch=c(18, 3,3), col = c("gray", "black","deepskyblue2"), legend = c("Regional fitting", "Local fitting", "Local BB fitting"))
    dev.off()
  }#end if pb at GP
  
}#end for GP
