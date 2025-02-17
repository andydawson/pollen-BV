library(raster)
library(VoCC)
library(ggplot2)
library(terra)
library(tidyterra)

locs_grid <- readRDS('data/grid_5.0.RDS')
# preds <- readRDS('/scratch/pollen_model_outputs/polya-gamma-predictions_5.0_overdispersed.RDS')
taxa = readRDS("data/taxa_5.0.RDS")

pbs_ll = readRDS('data/map-data/geographic/pbs_ll.RDS')

pbs = readRDS('data/map-data/geographic/pbs.RDS')

######################################################################################################################
##
######################################################################################################################

# read raster masks (provides masks for spatial domain/resolution of genetic + ENM data)
# stack <- stack('data/map-data/study_region_daltonIceMask_lakesMasked_linearIceSheetInterpolation.tif')
# stack[stack == 1] <- NA
# saveRDS(stack, 'data/map-data/study_region_mask_glacierNA.RDS')
stack <- readRDS('data/map-data/study_region_mask_glacierNA.RDS')
mask_bins = (700:0)*30
names(stack) <- mask_bins

# # convert rasterstack masks into list of rasters, reverse order of items to correspond with pollen time
# mask_list <- unstack(stack)
# n_times <- length(mask_list)
# mask_names <- NA
# for(i in 1:n_times){
#   mask_names[i] <- names(mask_list[[i]])
#   mask_names[i] <- as.integer(substr(mask_names[i], start = 2, stop = nchar(mask_names[i]))) * 30
# }
# names(mask_list) <- mask_names
# mask_list <- rev(mask_list)
# mask_names <- rev(mask_names)
# interp_times <- as.numeric(mask_names)
# 
# # create vector of pollen times
# time <- c(-70, seq(705, by = 990, length.out = 22))
# time <- data.frame(timeFrom = time[-1], timeTo = time[-23])
# time$time_mid <- (time$timeFrom + time$timeTo) / 2  # to get midpoint of time bins
# times <- time$time_mid
# 
# # interpolate pollen rasters to temporal resolution of 30yrs
# # convert list (iterations) of lists (times) to list of rasterstacks (stacked by time)
# rasterstack_list <- list()
# for(i in 1:n_iter){
#   rasterstack_list[[i]] <- stack(raster_list[[i]])
# }
# 
# # interpolate for each iteration
# interp_list <- list()
# for(i in 1:n_iter){
#   interp_list[[i]] <- interpolateRasters(rasterstack_list[[i]], interpFrom = times, interpTo = interp_times)
#   interp_list[[i]] <- unstack(interp_list[[i]])
# }
# saveRDS(interp_list, 'output/interp_preds_n50_v4.0.RDS')
# 
# 
# # mask out ice/water
# masked_interp_list <- interp_list
# for(i in 1:n_iter){
#   for(j in 1:n_times){
#     masked_interp_list[[i]][[j]] <- raster::resample(masked_interp_list[[i]][[j]], y = mask_list[[j]])
#     masked_interp_list[[i]][[j]] <- mask(masked_interp_list[[i]][[j]], mask = mask_list[[j]])
#   }
# }
# saveRDS(masked_interp_list, 'output/interp_masked_preds_n50_v4.0.RDS')

######################################################################################################################
##
######################################################################################################################

# In both paleoclimatic simulations, the timescale is expressed as time before present, where ‘present’
# follows radiocarbon dating conventions and is deﬁned as 0 ka BP (1950). Decade 0 is deﬁned as from 
# 1st January 1951 to 31st December 1960. CCSM3 data extends from decade −2200 to +3, which means it
# ends on December 31st, 1990,

# raster brick name suggests first layer is -2200
ccsm3_bins = rev(c(-3, -2, -1, seq(0, 2200, by=1)))*10
n_ccsm3_bins = length(ccsm3_bins)

# tmax_month = raster('data/GCM/ccsm3_22-0k_temp.nc', level=i, varname='tmax')
# clim_rast = brick('data/GCM/ccsm3_22-0k_temp.nc', level=i, varname='tmax')

clim_rast = brick('data/GCM/ccsm3_22-0k_temp.nc', level=6, varname='tmax')


# foo = brick('data/GCM/ccsm3_22-0k_temp.nc', level=i, varname='tmax')

# clim_rast_annual = brick('data/GCM/ccsm3_22-0k_temp.nc', level=1, varname='tmax')


# clim_mask_stack <- list()
# # for(i in 1:n_iter){
# # rasterstack_list[[i]] <- stack(raster_list[[i]])
# # }
# for (i in 1:n_ccsm3_bins){
#   print(i) 
#   
#   bin = ccsm3_bins[i]
#   
#   mask_bin_id = which.min(abs(mask_bins - bin))
#   
#   mask_proj = projectRaster(stack[[mask_bin_id]], clim_rast[[i]])
#   clim_mask = raster::mask(clim_rast[[i]], mask=mask_proj)
#   clim_mask_stack[[i]] <- stack(clim_mask)
# }
# 
# clim_mask_stack = brick(clim_mask_stack)
# 
# saveRDS(clim_mask_stack, 'data/clim_mask_stack.RDS')
# 
# clim_mask_stack = readRDS('data/clim_mask_stack.RDS')

######################################################################################################################
## HERE: climate velocity for range for each taxon
######################################################################################################################

# In both paleoclimatic simulations, the timescale is expressed as time before present, where ‘present’
# follows radiocarbon dating conventions and is deﬁned as 0 ka BP (1950). Decade 0 is deﬁned as from 
# 1st January 1951 to 31st December 1960. CCSM3 data extends from decade −2200 to +3, which means it
# ends on December 31st, 1990,

# raster brick name suggests first layer is -2200
ccsm3_bins = rev(c(-3, -2, -1, seq(0, 2200, by=1)))*10
n_ccsm3_bins = length(ccsm3_bins)

# tmax_month = raster('data/GCM/ccsm3_22-0k_temp.nc', level=i, varname='tmax')
# clim_rast = brick('data/GCM/ccsm3_22-0k_temp.nc', level=i, varname='tmax')

clim_rast = brick('data/GCM/ccsm3_22-0k_temp.nc', level=6, varname='tmax')
clim_rast = brick('data/GCM/ccsm3_22-0k_temp.nc', level=6, varname='tmin')
clim_rast = brick('data/GCM/ccsm3_22-0k_prcp.nc', level=6, varname='prcp')
clim_rast = brick('data/GCM/ccsm3_22-0k_et.nc', level=6, varname='pet')
clim_rast = brick('data/GCM/ccsm3_22-0k_et.nc', level=6, varname='aet')


# clim_rast = terra::rast(raster::stack(clim_rast))

# clim_rast = rast('data/GCM/ccsm3_22-0k_temp.nc', subds='tmax', lyrs=6)

# # clim_rast_annual = brick('data/GCM/ccsm3_22-0k_temp.nc', level=1, varname='tmax')
# 
# 
# 
# 
# # mask out ice/water
# 

name_pair = data.frame(suffix = c('temp', 'temp', 'prcp', 'et', 'et', 'gdd', 'gdd', 'vap'), 
                       varname = c('tmax', 'tmin', 'prcp', 'pet', 'aet', 'gdd0', 'gdd5', 'vap'))

months = c('jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec')

N_climvars = nrow(name_pair)

for (n in 4:N_climvars){
  
  suffix = name_pair$suffix[n]
  var = name_pair$varname[n]
  
  for (k in 1:12){
    
    print(paste0('Var ', var, '; Month ', k))
    
    clim_rast = brick(paste0('data/GCM/ccsm3_22-0k_', suffix, '.nc'), 
                      level=k, 
                      varname=var)
    
    
    clim_rast_stack <- list()
    
    for (i in 1:n_ccsm3_bins){
      # print(i)
      clim_rast_stack[[i]] <- stack(clim_rast[[i]])
    }
    
    clim_rast_brick = brick(clim_rast_stack)
    clim_rast_spat = rast(clim_rast_brick)
    
    writeRaster(clim_rast_spat, paste0('data/clim_rast/clim_rast_spat_', var, '_', months[k], '.tif'), overwrite=TRUE)
  }
}


# ######################################################################################################################
# ##
# ######################################################################################################################
# 
# # time periods: 16-14, 14-12, 12-10,10-7,7-4,4-1
# 
# # periods = data.frame(period_mid = c(15, 13, 11, 8.5, 5.5, 2.5),
# #                      period_old = c(16, 14, 12, 10, 7, 4),
# #                      period_young = c(14, 12, 10, 7, 4, 1))
# # periods = data.frame(period_mid = seq(19, 1, by=-2),
# #                      period_old = seq(20, 2, by=-2),
# #                      period_young =seq(18, 0, by=-2))
# # periods = periods * 1000
# 
# time = c(-70, seq(705, by = 990, length.out = 22))
# 
# periods = data.frame(period_old = time[-1], 
#                      period_young = time[-23])
# periods$period_mid = (periods$period_old + periods$period_young)/2
# 
# periods = periods[nrow(periods):1, ]
# n_periods = nrow(periods) - 1
# 
# gv_all = data.frame(x = numeric(0),
#                     y = numeric(0),
#                     voccMag = numeric(0),
#                     period = numeric(0))
# 
# for (i in 1:2){#n_periods){
#   
#   print(paste0("Period ", i, " of ", n_periods))
#   
#   period_ids = which((ccsm3_bins < periods[i,'period_old']) & (ccsm3_bins > periods[i+1, 'period_young']))
#   
#   clim_mask_stack_sub = clim_mask_stack[[period_ids]]
#   clim_period = brick(clim_mask_stack_sub)
#   
#   ccsm3_bins[period_ids]
#   period_ends = ccsm3_bins[period_ids]
#   
#   # temporal trend
#   # th minimum number of obs needed to calculate trend
#   vt <- tempTrend(clim_period, th = 10)
#   
#   # plot(vt)
#   
#   # spatial gradient
#   vg <- spatGrad(clim_period, th = 0.0001, projected = FALSE)
#   
#   # plot(vg)
#   
#   # climate velocity
#   gv <- gVoCC(vt, vg)
#   
#   # plot(gv)
#   
#   gv_df = as.data.frame(gv$voccMag, xy=TRUE)
#   # gv_df$voccMag = gv_df$voccMag/abs(diff(period_ends)/10)
#   # gv_df$voccMag = gv_df$voccMag/abs(diff(period_ends)/10)
#   
#   # hist(gv_df$voccMag, breaks=100)
#   
#   gv_all = rbind(gv_all, 
#                  data.frame(gv_df, period = (periods[i, 'period_mid'] + periods[i+1, 'period_mid'])/2))
# }
# 
# # gv_all[which(gv_all$voccMag>100), 'voccMag'] = 100
# # gv_all[which(gv_all$voccMag<(-100)), 'voccMag'] = -100
# 
# saveRDS(gv_all, 'output/climate_velocity.RDS')
# 
# gv_all$voccMag_binned = cut(gv_all$voccMag, breaks=c(-100, -50, -10, -1, 0, 1, 10, 50, 100), labels=FALSE)
# 
# ggplot() +
#   geom_histogram(data=gv_all, aes(x=sign(voccMag)*sqrt(abs(voccMag)), y=after_stat(density)), bins=30) +
#   # geom_density(data=log(gv_all), aes(x=voccMag)) +
#   geom_vline(xintercept=0, colour='red', lty=2) +
#   facet_wrap(~period) +
#   theme_bw()
# 
# # hist(gv$voccMag)
# 
# ggplot() +
#   geom_boxplot(data=gv_all, aes(x=factor(period), y=sign(voccMag)*sqrt(abs(voccMag)))) +
#   # geom_density(data=log(gv_all), aes(x=voccMag)) +
#   geom_hline(yintercept=0, colour='red', lty=2) +
#   # facet_wrap(~period) +
#   theme_bw()

######################################################################################################################
##
######################################################################################################################

# time periods: 16-14, 14-12, 12-10,10-7,7-4,4-1

# periods = data.frame(period_mid = c(15, 13, 11, 8.5, 5.5, 2.5),
#                      period_old = c(16, 14, 12, 10, 7, 4),
#                      period_young = c(14, 12, 10, 7, 4, 1))
# periods = data.frame(period_mid = seq(19, 1, by=-2),
#                      period_old = seq(20, 2, by=-2),
#                      period_young =seq(18, 0, by=-2))
# periods = periods * 1000

# clim_mask_stack = readRDS('data/clim_mask_stack.RDS')
# foo = rast(clim_mask_stack)

clim_rast = rast('data/clim_rast_spat.tif')

# taxon_range_mask_iter = readRDS(paste0('output/',taxa[1],'_range_mask_holocene_iter_200_v5.0.RDS'))
taxon_range_mask = readRDS(paste0('output/',taxa[17],'_range_mask_holocene_filled_iter_200_v5.0.RDS'))

presence_mask_poly = vect('output/taxa_presence_mask_pollen_hull_poly')

proj_mask = crs(presence_mask_poly)


#clim_rast = project(clim_rast, taxon_range_mask)

clim_rast = project(clim_rast, proj_mask)

N_taxa = length(taxa)

time_ends <- c(-70, seq(705, by = 990, length.out = 22))
# time <- data.frame(timeFrom = time[-1], timeTo = time[-23])
# time$time_mid <- (time$timeFrom + time$timeTo) / 2
# # specify time bins; use median of time bins for subsequent work
# time <- c(-70, seq(705, by = 990, length.out = 22))
# time <- data.frame(timeFrom = time[-1], timeTo = time[-23])
# time$time_mid <- (time$timeFrom + time$timeTo) / 2

periods = data.frame(period_old = time_ends[-1], 
                     period_young = time_ends[-23])
periods$period_mid = (periods$period_old + periods$period_young)/2

periods = periods[nrow(periods):1, ]

n_periods = nrow(periods) - 1


for (n in 1:N_climvars){
  
  # suffix = name_pair$suffix[n]
  var = name_pair$varname[n]
  
  gv_all = data.frame(var = character(0),
                      month = character(0),
                      x = numeric(0),
                      y = numeric(0),
                      voccMag = numeric(0),
                      period = numeric(0))
  
  for (k in 1:12){
    
    gv_month = data.frame(var = character(0),
                          month = character(0),
                          x = numeric(0),
                          y = numeric(0),
                          voccMag = numeric(0),
                          period = numeric(0))
    
    month = months[k]
    
    print(paste0('Var ', var, '; Month ', k))
    
    clim_rast = rast(paste0('data/clim_rast/clim_rast_spat_', var, '_', months[k], '.tif'))
    
    clim_rast = project(clim_rast, proj_mask)
    
    for (i in 1:n_periods){
      
      print(paste0("Period ", i, " of ", n_periods))
      
      period_ids = which((ccsm3_bins < periods[i,'period_old']) & (ccsm3_bins > periods[i+1, 'period_young']))
      
      # clim_mask_stack_sub = subset(clim_mask_stack, period_ids)
      clim_rast_sub = subset(clim_rast, period_ids)
      clim_rast_sub_mask = clim_rast_sub
      
      clim_period = brick(clim_rast_sub_mask)
      
      ccsm3_bins[period_ids]
      period_ends = ccsm3_bins[period_ids]
      
      # temporal trend
      # th minimum number of obs needed to calculate trend
      vt = tempTrend(clim_period, th = 10)
      # plot(vt)
      
      # spatial gradient
      vg = spatGrad(clim_period, th = 0.0001, projected = FALSE)
      # plot(vg)
      
      # climate velocity
      gv = gVoCC(vt, vg)
      # plot(gv)
      
      gv_rast = rast(gv)
      
      # gv_mask = terra::mask(gv_rast, taxon_mask)
      
      # gv_df = as.data.frame(gv$voccMag, xy=TRUE)
      gv_df = as.data.frame(gv_rast$voccMag, xy=TRUE)
      
      gv_month = rbind(gv_month, 
                       data.frame(var = var,
                                  month = month,
                                  gv_df, 
                                  period = (periods[i, 'period_mid'] + periods[i+1, 'period_mid'])/2))
      
      gv_all = rbind(gv_all, 
                     data.frame(var = var,
                                month = month,
                                gv_df, 
                                period = (periods[i, 'period_mid'] + periods[i+1, 'period_mid'])/2))
      
    }
    
    
    saveRDS(gv_month, paste0('output/climate_velocity/', var, '_', month, '_CV.RDS'))
    
    saveRDS(gv_month, paste0('output/climate_velocity/', var, '_', month, '_CV.RDS'))
    
    
  }
  
  saveRDS(gv_all, paste0('output/climate_velocity/', var, '_CV.RDS'))
}


n=2

# suffix = name_pair$suffix[n]
var = name_pair$varname[n]

gv_all = data.frame(var = character(0),
                    month = character(0),
                    x = numeric(0),
                    y = numeric(0),
                    voccMag = numeric(0),
                    period = numeric(0))

for (k in 1:12){
  
  
  month = months[k]
  
  print(paste0('Var ', var, '; Month ', k))
  
  gv_month = readRDS(paste0('output/climate_velocity/', var, '_', month, '_CV.RDS'))
  
  
  gv_all = rbind(gv_all, 
                 gv_month)
}

saveRDS(gv_all, paste0('output/climate_velocity/', var, '_CV.RDS'))

# for (i in 1:n_periods){
#   
#   print(paste0("Period ", i, " of ", n_periods))
#   
#   period_ids = which((ccsm3_bins < periods[i,'period_old']) & (ccsm3_bins > periods[i+1, 'period_young']))
#   
#   # clim_mask_stack_sub = subset(clim_mask_stack, period_ids)
#   clim_rast_sub = subset(clim_rast, period_ids)
#   clim_rast_sub_mask = clim_rast_sub
#   
#   clim_period = brick(clim_rast_sub_mask)
#   
#   ccsm3_bins[period_ids]
#   period_ends = ccsm3_bins[period_ids]
#   
#   # temporal trend
#   # th minimum number of obs needed to calculate trend
#   vt <- tempTrend(clim_period, th = 10)
#   # plot(vt)
#   
#   # spatial gradient
#   vg <- spatGrad(clim_period, th = 0.0001, projected = FALSE)
#   # plot(vg)
#   
#   # climate velocity
#   gv <- gVoCC(vt, vg)
#   # plot(gv)
#   
#   gv_rast = rast(gv)
#   
#   # gv_mask = terra::mask(gv_rast, taxon_mask)
#   
#   # gv_df = as.data.frame(gv$voccMag, xy=TRUE)
#   gv_df = as.data.frame(gv_rast$voccMag, xy=TRUE)
#   
#   gv_all = rbind(gv_all, 
#                  data.frame(gv_df, 
#                             period = (periods[i, 'period_mid'] + periods[i+1, 'period_mid'])/2))
# }
# 
# saveRDS(gv_all, 'output/TMAX_JUNE_CV.RDS')


# gv_taxon = data.frame(x = numeric(0),
#                       y = numeric(0),
#                       voccMag = numeric(0),
#                       period = numeric(0),
#                       taxon = character(0))
# for (tx in 1:N_taxa){
#   
#   print(paste0("Taxon ", tx, ": ", taxa[tx]))
#   
#   # taxon_range_mask_iter = readRDS(paste0('output/',taxa[tx],'_range_mask_holocene_iter_200_v5.0.RDS'))
#   # taxon_range_mask_iter_filled = readRDS(paste0('output/',taxa[tx],'_range_mask_holocene_filled_iter_200_v5.0.RDS'))
#   # 
#   # taxon_range_mask = lapply(taxon_range_mask_iter, rast)
#   # taxon_range_mask_iter = raster(taxon_range_mask_iter)
#   # taxon_range_mask_iter = subset(taxon_range_mask_iter, 1)
#   
#   # ggplot() +
#   #   geom_spatraster(data=taxon_range_mask_iter, aes(fill=max)) +
#   #   scale_fill_viridis_c(na.value = "pink")
#   
#   
#   # taxon_range_mask_iter_filled = subset(taxon_range_mask_iter_filled, 1)
#   
#   # ggplot() +
#   #   geom_spatraster(data=taxon_range_mask_iter_filled, aes(fill=focal_mean)) +
#   #   scale_fill_viridis_c(na.value = "pink")
#   
#   taxon_mask = presence_mask_poly[tx,]
#   
#   for (i in 1:n_periods){
#     
#     period_ids = which((ccsm3_bins < periods[i,'period_old']) & (ccsm3_bins > periods[i+1, 'period_young']))
#     # period_id_old = which.min(abs(ccsm3_bins - periods[i,'period_old'] ))
#     # period_id_young = which.min(abs(ccsm3_bins - periods[i,'period_young'] ))
#     # period_ids = c(period_id_old, period_id_young)
#     
#     # clim_mask_stack_sub = subset(clim_mask_stack, period_ids)
#     clim_rast_sub = subset(clim_rast, period_ids)
#     
#     
#     # ggplot() +
#     #   geom_spatraster(data=clim_rast_sub[[1]]) +
#     #   scale_fill_viridis_c(na.value = "pink")
#     
#     # clim_rast_sub_mask = terra::mask(clim_rast_sub, taxon_range_mask_iter)
#     clim_rast_sub_mask = clim_rast_sub
#     
#     # ggplot() +
#     #   geom_spatraster(data=clim_rast_sub_mask[[1]]) +
#     #   scale_fill_viridis_c(na.value = "pink")
#     # 
#     # clim_period = brick(clim_mask_stack_sub)
#     clim_period = brick(clim_rast_sub_mask)
#     
#     ccsm3_bins[period_ids]
#     period_ends = ccsm3_bins[period_ids]
#     
#     # temporal trend
#     # th minimum number of obs needed to calculate trend
#     vt <- tempTrend(clim_period, th = 10)
#     
#     # plot(vt)
#     
#     # spatial gradient
#     vg <- spatGrad(clim_period, th = 0.0001, projected = FALSE)
#     
#     # plot(vg)
#     
#     
#     # climate velocity
#     gv <- gVoCC(vt, vg)
#     
#     # plot(gv)
#     
#     gv_rast = rast(gv)
#     
#     gv_mask = terra::mask(gv_rast, taxon_mask)
#     
#     # gv_df = as.data.frame(gv$voccMag, xy=TRUE)
#     gv_df = as.data.frame(gv_mask$voccMag, xy=TRUE)
#     
#     # gv_df$voccMag = gv_df$voccMag/abs(diff(period_ends)/10)
#     # gv_df$voccMag = gv_df$voccMag/abs(diff(period_ends)/10)
#     
#     # hist(gv_df$voccMag, breaks=100)
#     
#     gv_taxon = rbind(gv_taxon, 
#                      data.frame(gv_df, 
#                                 period = (periods[i, 'period_mid'] + periods[i+1, 'period_mid'])/2, 
#                                 taxon = rep(taxa[tx]))
#     )
#   }
# }
# 
# 
# # gv_taxon[which(gv_taxon$voccMag>500), 'voccMag'] = 500
# # gv_taxon[which(gv_taxon$voccMag<(-100)), 'voccMag'] = -100
# 
# saveRDS(gv_taxon, 'output/climate_velocity_taxon.RDS')

gv_all = data.frame(x = numeric(0),
                    y = numeric(0),
                    voccMag = numeric(0),
                    period = numeric(0))


for (i in 1:n_periods){
  
  print(paste0("Period ", i, " of ", n_periods))
  
  period_ids = which((ccsm3_bins < periods[i,'period_old']) & (ccsm3_bins > periods[i+1, 'period_young']))
  
  # clim_mask_stack_sub = subset(clim_mask_stack, period_ids)
  clim_rast_sub = subset(clim_rast, period_ids)
  clim_rast_sub_mask = clim_rast_sub
  
  clim_period = brick(clim_rast_sub_mask)
  
  ccsm3_bins[period_ids]
  period_ends = ccsm3_bins[period_ids]
  
  # temporal trend
  # th minimum number of obs needed to calculate trend
  vt <- tempTrend(clim_period, th = 10)
  # plot(vt)
  
  # spatial gradient
  vg <- spatGrad(clim_period, th = 0.0001, projected = FALSE)
  # plot(vg)
  
  # climate velocity
  gv <- gVoCC(vt, vg)
  # plot(gv)
  
  gv_rast = rast(gv)
  
  # gv_mask = terra::mask(gv_rast, taxon_mask)
  
  # gv_df = as.data.frame(gv$voccMag, xy=TRUE)
  gv_df = as.data.frame(gv_rast$voccMag, xy=TRUE)
  
  gv_all = rbind(gv_all, 
                 data.frame(gv_df, 
                            period = (periods[i, 'period_mid'] + periods[i+1, 'period_mid'])/2))
}


# gv_taxon[which(gv_taxon$voccMag>500), 'voccMag'] = 500
# gv_taxon[which(gv_taxon$voccMag<(-100)), 'voccMag'] = -100

saveRDS(gv_all, 'output/TMAX_JUNE_CV.RDS')

# gv_taxon$voccMag_binned = cut(gv_taxon$voccMag, 
#                               breaks=c(-100, -50, -10, -1, 0, 1, 10, 50, 100, 200, 500), labels=FALSE)
# 
# ggplot() +
#   geom_histogram(data=gv_taxon, aes(x=sign(voccMag)*sqrt(abs(voccMag)), y=after_stat(density)), bins=30) +
#   # geom_density(data=log(gv_all), aes(x=voccMag)) +
#   geom_vline(xintercept=0, colour='red', lty=2) +
#   facet_wrap(~period) +
#   theme_bw()
# 
# # hist(gv$voccMag)
# 
# ggplot() +
#   geom_boxplot(data=gv_all, aes(x=factor(period), y=sign(voccMag)*sqrt(abs(voccMag)))) +
#   # geom_density(data=log(gv_all), aes(x=voccMag)) +
#   geom_hline(yintercept=0, colour='red', lty=2) +
#   # facet_wrap(~period) +
#   theme_bw()

# saveRDS(gv_all, paste0('output/climate_velocity/', var, '_CV.RDS'))



# gv_all = readRDS('output/TMAX_JUNE_CV.RDS')

presence_mask_poly = vect('output/taxa_presence_mask_pollen_hull_poly')

proj_mask = crs(presence_mask_poly)

clim_rast = project(clim_rast, proj_mask)

for (i in 4:N_climvars){
  
  varname = name_pair$varname[i]
  
  gv_all = readRDS(paste0('output/climate_velocity/', varname, '_CV.RDS'))
  
  gv_taxa = data.frame(var = character(0),
                       month = character(0),
                       x = numeric(0),
                       y = numeric(0),
                       voccMag = numeric(0),
                       period = numeric(0),
                       taxon = character(0))
  
  for (tx in 1:N_taxa){
    
    print(paste0("Taxon ", tx, ": ", taxa[tx]))
    
    taxon_mask = presence_mask_poly[tx,]
    
    gv_in = terra::extract(taxon_mask, gv_all[,c('x', 'y')])
    
    gv_taxon = gv_all[which(!is.na(gv_in$FID)),]
    
    
    gv_taxa = rbind(gv_taxa, 
                    data.frame(gv_taxon, 
                               # period = (periods[i, 'period_mid'] + periods[i+1, 'period_mid'])/2, 
                               taxon = rep(taxa[tx])))
    
  }
  
  saveRDS(gv_taxa, paste0('output/climate_velocity/', varname, '_CV_taxa.RDS'))
  
}

# gv_all = readRDS('output/TMAX_JUNE_CV.RDS')
# 
# presence_mask_poly = vect('output/taxa_presence_mask_pollen_hull_poly')
# 
# proj_mask = crs(presence_mask_poly)
# 
# clim_rast = project(clim_rast, proj_mask)
# 
# gv_taxa = data.frame(x = numeric(0),
#                      y = numeric(0),
#                      voccMag = numeric(0),
#                      period = numeric(0),
#                      taxon = character(0))
# 
# for (tx in 1:N_taxa){
#   
#   print(paste0("Taxon ", tx, ": ", taxa[tx]))
#   
#   taxon_mask = presence_mask_poly[tx,]
#   
#   gv_in = terra::extract(taxon_mask, gv_all[,c('x', 'y')])
#   
#   gv_taxon = gv_all[which(!is.na(gv_in$FID)),]
#   
#   
#   gv_taxa = rbind(gv_taxa, 
#                   data.frame(gv_taxon, 
#                              # period = (periods[i, 'period_mid'] + periods[i+1, 'period_mid'])/2, 
#                              taxon = rep(taxa[tx])))
#   
# }
# 
# saveRDS(gv_taxa, 'output/TMAX_JUNE_CV_taxa.RDS')

# ######################################################################################################################
# ## climate velocity for range for each taxon
# ######################################################################################################################
# 
# # In both paleoclimatic simulations, the timescale is expressed as time before present, where ‘present’
# # follows radiocarbon dating conventions and is deﬁned as 0 ka BP (1950). Decade 0 is deﬁned as from 
# # 1st January 1951 to 31st December 1960. CCSM3 data extends from decade −2200 to +3, which means it
# # ends on December 31st, 1990,
# 
# # raster brick name suggests first layer is -2200
# ccsm3_bins = rev(c(-3, -2, -1, seq(0, 2200, by=1)))*10
# n_ccsm3_bins = length(ccsm3_bins)
# 
# # tmax_month = raster('data/GCM/ccsm3_22-0k_temp.nc', level=i, varname='tmax')
# # clim_rast = brick('data/GCM/ccsm3_22-0k_temp.nc', level=i, varname='tmax')
# 
# clim_rast = brick('data/GCM/ccsm3_22-0k_temp.nc', level=6, varname='tmax')
# 
# # clim_rast_annual = brick('data/GCM/ccsm3_22-0k_temp.nc', level=1, varname='tmax')
# 
# # mask out ice/water
# 
# clim_mask_stack <- list()
# # for(i in 1:n_iter){
# # rasterstack_list[[i]] <- stack(raster_list[[i]])
# # }
# for (i in 1:n_ccsm3_bins){
#   print(i) 
#   
#   bin = ccsm3_bins[i]
#   
#   mask_bin_id = which.min(abs(mask_bins - bin))
#   
#   mask_proj = projectRaster(stack[[mask_bin_id]], clim_rast[[i]])
#   clim_mask = raster::mask(clim_rast[[i]], mask=mask_proj)
#   clim_mask_stack[[i]] <- stack(clim_mask)
# }
# 
# clim_mask_stack = brick(clim_mask_stack)
# clim_mask_stack = rast(clim_mask_stack)
# 
# writeRaster(clim_mask_stack, 'data/clim_mask_stack.tif', overwrite=TRUE)
# 
# # saveRDS(clim_mask_stack, 'data/clim_mask_stack.RDS')

######################################################################################################################
##
######################################################################################################################

# time periods: 16-14, 14-12, 12-10,10-7,7-4,4-1

# periods = data.frame(period_mid = c(15, 13, 11, 8.5, 5.5, 2.5),
#                      period_old = c(16, 14, 12, 10, 7, 4),
#                      period_young = c(14, 12, 10, 7, 4, 1))
# periods = data.frame(period_mid = seq(19, 1, by=-2),
#                      period_old = seq(20, 2, by=-2),
#                      period_young =seq(18, 0, by=-2))
# periods = periods * 1000

# clim_mask_stack = readRDS('data/clim_mask_stack.RDS')
# foo = rast(clim_mask_stack)

clim_rast = rast('data/clim_rast_spat.tif')

# taxon_range_mask_iter = readRDS(paste0('output/',taxa[1],'_range_mask_holocene_iter_200_v5.0.RDS'))
taxon_range_mask = readRDS(paste0('output/',taxa[17],'_range_mask_holocene_filled_iter_200_v5.0.RDS'))

presence_mask_poly = vect('output/taxa_presence_mask_pollen_hull_poly')

proj_mask = crs(presence_mask_poly)


#clim_rast = project(clim_rast, taxon_range_mask)

clim_rast = project(clim_rast, proj_mask)

N_taxa = length(taxa)

time <- c(-70, seq(705, by = 990, length.out = 22))
# time <- data.frame(timeFrom = time[-1], timeTo = time[-23])
# time$time_mid <- (time$timeFrom + time$timeTo) / 2
# # specify time bins; use median of time bins for subsequent work
# time <- c(-70, seq(705, by = 990, length.out = 22))
# time <- data.frame(timeFrom = time[-1], timeTo = time[-23])
# time$time_mid <- (time$timeFrom + time$timeTo) / 2

periods = data.frame(period_old = time[-1], 
                     period_young = time[-23])
periods$period_mid = (periods$period_old + periods$period_young)/2

periods = periods[nrow(periods):1, ]

n_periods = nrow(periods) - 1

gv_taxon = data.frame(x = numeric(0),
                      y = numeric(0),
                      voccMag = numeric(0),
                      period = numeric(0),
                      taxon = character(0))
for (tx in 1:N_taxa){
  
  print(paste0("Taxon ", tx, ": ", taxa[tx]))
  
  # taxon_range_mask_iter = readRDS(paste0('output/',taxa[tx],'_range_mask_holocene_iter_200_v5.0.RDS'))
  # taxon_range_mask_iter_filled = readRDS(paste0('output/',taxa[tx],'_range_mask_holocene_filled_iter_200_v5.0.RDS'))
  # 
  # taxon_range_mask = lapply(taxon_range_mask_iter, rast)
  # taxon_range_mask_iter = raster(taxon_range_mask_iter)
  # taxon_range_mask_iter = subset(taxon_range_mask_iter, 1)
  
  # ggplot() +
  #   geom_spatraster(data=taxon_range_mask_iter, aes(fill=max)) +
  #   scale_fill_viridis_c(na.value = "pink")
  
  
  # taxon_range_mask_iter_filled = subset(taxon_range_mask_iter_filled, 1)
  
  # ggplot() +
  #   geom_spatraster(data=taxon_range_mask_iter_filled, aes(fill=focal_mean)) +
  #   scale_fill_viridis_c(na.value = "pink")
  
  taxon_mask = presence_mask_poly[tx,]
  
  for (i in 1:n_periods){
    
    period_ids = which((ccsm3_bins < periods[i,'period_old']) & (ccsm3_bins > periods[i+1, 'period_young']))
    # period_id_old = which.min(abs(ccsm3_bins - periods[i,'period_old'] ))
    # period_id_young = which.min(abs(ccsm3_bins - periods[i,'period_young'] ))
    # period_ids = c(period_id_old, period_id_young)
    
    # clim_mask_stack_sub = subset(clim_mask_stack, period_ids)
    clim_rast_sub = subset(clim_rast, period_ids)
    
    
    # ggplot() +
    #   geom_spatraster(data=clim_rast_sub[[1]]) +
    #   scale_fill_viridis_c(na.value = "pink")
    
    # clim_rast_sub_mask = terra::mask(clim_rast_sub, taxon_range_mask_iter)
    clim_rast_sub_mask = clim_rast_sub
    
    # ggplot() +
    #   geom_spatraster(data=clim_rast_sub_mask[[1]]) +
    #   scale_fill_viridis_c(na.value = "pink")
    # 
    # clim_period = brick(clim_mask_stack_sub)
    clim_period = brick(clim_rast_sub_mask)
    
    ccsm3_bins[period_ids]
    period_ends = ccsm3_bins[period_ids]
    
    # temporal trend
    # th minimum number of obs needed to calculate trend
    vt <- tempTrend(clim_period, th = 10)
    
    # plot(vt)
    
    # spatial gradient
    vg <- spatGrad(clim_period, th = 0.0001, projected = FALSE)
    
    # plot(vg)
    
    
    # climate velocity
    gv <- gVoCC(vt, vg)
    
    # plot(gv)
    
    gv_rast = rast(gv)
    
    gv_mask = terra::mask(gv_rast, taxon_mask)
    
    # gv_df = as.data.frame(gv$voccMag, xy=TRUE)
    gv_df = as.data.frame(gv_mask$voccMag, xy=TRUE)
    
    # gv_df$voccMag = gv_df$voccMag/abs(diff(period_ends)/10)
    # gv_df$voccMag = gv_df$voccMag/abs(diff(period_ends)/10)
    
    # hist(gv_df$voccMag, breaks=100)
    
    gv_taxon = rbind(gv_taxon, 
                     data.frame(gv_df, 
                                period = (periods[i, 'period_mid'] + periods[i+1, 'period_mid'])/2, 
                                taxon = rep(taxa[tx]))
    )
  }
}


# gv_taxon[which(gv_taxon$voccMag>500), 'voccMag'] = 500
# gv_taxon[which(gv_taxon$voccMag<(-100)), 'voccMag'] = -100

saveRDS(gv_taxon, 'output/climate_velocity_taxon.RDS')

gv_taxon$voccMag_binned = cut(gv_taxon$voccMag, 
                              breaks=c(-100, -50, -10, -1, 0, 1, 10, 50, 100, 200, 500), labels=FALSE)

ggplot() +
  geom_histogram(data=gv_taxon, aes(x=sign(voccMag)*sqrt(abs(voccMag)), y=after_stat(density)), bins=30) +
  # geom_density(data=log(gv_all), aes(x=voccMag)) +
  geom_vline(xintercept=0, colour='red', lty=2) +
  facet_wrap(~period) +
  theme_bw()

# hist(gv$voccMag)

ggplot() +
  geom_boxplot(data=gv_all, aes(x=factor(period), y=sign(voccMag)*sqrt(abs(voccMag)))) +
  # geom_density(data=log(gv_all), aes(x=voccMag)) +
  geom_hline(yintercept=0, colour='red', lty=2) +
  # facet_wrap(~period) +
  theme_bw()

######################################################################################################################
## taxon cv: more vars and more efficient
######################################################################################################################

# time periods: 16-14, 14-12, 12-10,10-7,7-4,4-1

# periods = data.frame(period_mid = c(15, 13, 11, 8.5, 5.5, 2.5),
#                      period_old = c(16, 14, 12, 10, 7, 4),
#                      period_young = c(14, 12, 10, 7, 4, 1))
# periods = data.frame(period_mid = seq(19, 1, by=-2),
#                      period_old = seq(20, 2, by=-2),
#                      period_young =seq(18, 0, by=-2))
# periods = periods * 1000

# clim_mask_stack = readRDS('data/clim_mask_stack.RDS')
# foo = rast(clim_mask_stack)

clim_rast = rast('data/clim_rast_spat.tif')

# taxon_range_mask_iter = readRDS(paste0('output/',taxa[1],'_range_mask_holocene_iter_200_v5.0.RDS'))
taxon_range_mask = readRDS(paste0('output/',taxa[17],'_range_mask_holocene_filled_iter_200_v5.0.RDS'))

presence_mask_poly = vect('output/taxa_presence_mask_pollen_hull_poly')

proj_mask = crs(presence_mask_poly)


#clim_rast = project(clim_rast, taxon_range_mask)

clim_rast = project(clim_rast, proj_mask)

N_taxa = length(taxa)

time <- c(-70, seq(705, by = 990, length.out = 22))
# time <- data.frame(timeFrom = time[-1], timeTo = time[-23])
# time$time_mid <- (time$timeFrom + time$timeTo) / 2
# # specify time bins; use median of time bins for subsequent work
# time <- c(-70, seq(705, by = 990, length.out = 22))
# time <- data.frame(timeFrom = time[-1], timeTo = time[-23])
# time$time_mid <- (time$timeFrom + time$timeTo) / 2

periods = data.frame(period_old = time[-1], 
                     period_young = time[-23])
periods$period_mid = (periods$period_old + periods$period_young)/2

periods = periods[nrow(periods):1, ]

n_periods = nrow(periods) - 1

gv_taxon = data.frame(x = numeric(0),
                      y = numeric(0),
                      voccMag = numeric(0),
                      period = numeric(0),
                      taxon = character(0))
# for (tx in 1:N_taxa){
#   
#   print(paste0("Taxon ", tx, ": ", taxa[tx]))




for (i in 1:n_periods){
  
  print(paste0("Time period ", i, " of ", n_periods))
  
  period_ids = which((ccsm3_bins < periods[i,'period_old']) & (ccsm3_bins > periods[i+1, 'period_young']))
  
  # clim_mask_stack_sub = subset(clim_mask_stack, period_ids)
  clim_rast_sub = subset(clim_rast, period_ids)
  
  
  # clim_rast_sub_mask = terra::mask(clim_rast_sub, taxon_range_mask_iter)
  clim_rast_sub_mask = clim_rast_sub
  
  # 
  # clim_period = brick(clim_mask_stack_sub)
  clim_period = brick(clim_rast_sub_mask)
  
  ccsm3_bins[period_ids]
  period_ends = ccsm3_bins[period_ids]
  
  # temporal trend
  # th minimum number of obs needed to calculate trend
  vt <- tempTrend(clim_period, th = 10)
  
  # plot(vt)
  
  # spatial gradient
  vg <- spatGrad(clim_period, th = 0.0001, projected = FALSE)
  
  # plot(vg)
  
  
  # climate velocity
  gv <- gVoCC(vt, vg)
  
  # plot(gv)
  
  gv_rast = rast(gv)
  
  for (tx in 1:N_taxa){
    
    print(paste0("Taxon ", tx, ": ", taxa[tx]))
    
    taxon_mask = presence_mask_poly[tx,]
    
    
    gv_mask = terra::mask(gv_rast, taxon_mask)
    
    # gv_df = as.data.frame(gv$voccMag, xy=TRUE)
    gv_df = as.data.frame(gv_mask$voccMag, xy=TRUE)
    
    # gv_df$voccMag = gv_df$voccMag/abs(diff(period_ends)/10)
    # gv_df$voccMag = gv_df$voccMag/abs(diff(period_ends)/10)
    
    # hist(gv_df$voccMag, breaks=100)
    
    gv_taxon = rbind(gv_taxon, 
                     data.frame(gv_df, 
                                period = (periods[i, 'period_mid'] + periods[i+1, 'period_mid'])/2, 
                                taxon = rep(taxa[tx]))
    )
  }
}


# gv_taxon[which(gv_taxon$voccMag>500), 'voccMag'] = 500
# gv_taxon[which(gv_taxon$voccMag<(-100)), 'voccMag'] = -100

saveRDS(gv_taxon, 'output/TMAX_JUNE_CV_taxon.RDS')

gv_taxon$voccMag_binned = cut(gv_taxon$voccMag, 
                              breaks=c(-100, -50, -10, -1, 0, 1, 10, 50, 100, 200, 500), labels=FALSE)

ggplot() +
  geom_histogram(data=gv_taxon, aes(x=sign(voccMag)*sqrt(abs(voccMag)), y=after_stat(density)), bins=30) +
  # geom_density(data=log(gv_all), aes(x=voccMag)) +
  geom_vline(xintercept=0, colour='red', lty=2) +
  facet_wrap(~period) +
  theme_bw()

# hist(gv$voccMag)

ggplot() +
  geom_boxplot(data=gv_all, aes(x=factor(period), y=sign(voccMag)*sqrt(abs(voccMag)))) +
  # geom_density(data=log(gv_all), aes(x=voccMag)) +
  geom_hline(yintercept=0, colour='red', lty=2) +
  # facet_wrap(~period) +
  theme_bw()


# ######################################################################################################################
# ## annual averages at 500-year resolution
# ######################################################################################################################
# 
# # In both paleoclimatic simulations, the timescale is expressed as time before present, where ‘present’
# # follows radiocarbon dating conventions and is deﬁned as 0 ka BP (1950). Decade 0 is deﬁned as from 
# # 1st January 1951 to 31st December 1960. CCSM3 data extends from decade −2200 to +3, which means it
# # ends on December 31st, 1990,
# 
# # raster brick name suggests first layer is -2200
# ccsm3_bins = seq(0, 22000, by=500)
# n_ccsm3_bins = length(ccsm3_bins)
# 
# clim_rast = brick('data/GCM/ccsm3_22-0k_temp.nc', level=6, varname='tmax')
# 
# 
# clim_mask_stack <- list()
# # for(i in 1:n_iter){
# # rasterstack_list[[i]] <- stack(raster_list[[i]])
# # }
# for (i in 1:n_ccsm3_bins){
#   print(i) 
#   
#   bin = ccsm3_bins[i]
#   
#   clim_bin = raster(paste0('data/GCM/CCSM/', bin, 'BP/an_avg_TMAX.tif'))
#   crs(clim_bin) = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
#   
#   mask_bin_id = which.min(abs(mask_bins - bin))
#   
#   mask_proj = projectRaster(stack[[mask_bin_id]], clim_bin)
#   clim_mask = mask(clim_bin, mask=mask_proj)
#   clim_mask_stack[[i]] <- stack(clim_mask)
# }
# 
# clim_mask_stack = brick(clim_mask_stack)
# 
# saveRDS(clim_mask_stack, 'data/clim_mask_stack.RDS')
# 
# ######################################################################################################################
# ##
# ######################################################################################################################
# 
# # time periods: 16-14, 14-12, 12-10,10-7,7-4,4-1
# 
# # periods = data.frame(period_mid = c(15, 13, 11, 8.5, 5.5, 2.5),
# #                      period_old = c(16, 14, 12, 10, 7, 4),
# #                      period_young = c(14, 12, 10, 7, 4, 1))
# periods = data.frame(period_mid = seq(19, 1, by=-2),
#                      period_old = seq(20, 2, by=-2),
#                      period_young =seq(18, 0, by=-2))
# periods = periods * 1000
# 
# n_periods = nrow(periods)
# 
# gv_all = data.frame(x = numeric(0),
#                     y = numeric(0),
#                     voccMag = numeric(0),
#                     period = numeric(0))
# 
# for (i in 1:n_periods){
#   
#   period_ids = which((ccsm3_bins < periods[i,'period_old']) & (ccsm3_bins > periods[i, 'period_young']))
#   # period_id_old = which.min(abs(ccsm3_bins - periods[i,'period_old'] ))
#   # period_id_young = which.min(abs(ccsm3_bins - periods[i,'period_young'] ))
#   # period_ids = c(period_id_old, period_id_young)
#   
#   clim_mask_stack_sub = clim_mask_stack[[period_ids]]
#   # clim_period = brick(clim_mask_stack_sub)
#   
#   ccsm3_bins[period_ids]
#   period_ends = ccsm3_bins[period_ids]
#   
#   # temporal trend
#   # th minimum number of obs needed to calculate trend
#   vt <- tempTrend(clim_mask_stack, th = 10)
#   
#   plot(vt)
#   
#   # spatial gradient
#   vg <- spatGrad(clim_mask_stack, th = 0.0001, projected = FALSE)
#   
#   plot(vg)
#   
#   # climate velocity
#   gv <- gVoCC(vt, vg)
#   
#   plot(gv)
#   
#   gv_df = as.data.frame(gv$voccMag, xy=TRUE)
#   # gv_df$voccMag = gv_df$voccMag/abs(diff(period_ends)/10)
#   # gv_df$voccMag = gv_df$voccMag/abs(diff(period_ends)/10)
#   
#   hist(gv_df$voccMag, breaks=100)
#   
#   gv_all = rbind(gv_all, 
#                  data.frame(gv_df, period = periods[i, 'period_mid']))
# }
# 
# gv_all[which(gv_all$voccMag>100), 'voccMag'] = 100
# 
# saveRDS(gv_all, 'output/climate_velocity.RDS')
# 
# gv_all$voccMag_binned = cut(gv_all$voccMag, breaks=c(-100, -50, -10, -1, 0, 1, 10, 50, 100, 200), labels=FALSE)
# 
# ggplot() +
#   geom_histogram(data=gv_all, aes(x=sign(voccMag)*sqrt(abs(voccMag)), y=..density..), bins=30) +
#   # geom_density(data=log(gv_all), aes(x=voccMag)) +
#   geom_vline(xintercept=0, colour='red', lty=2) +
#   facet_wrap(~period) +
#   theme_bw()
# 
# # hist(gv$voccMag)
# 
# ggplot() +
#   geom_boxplot(data=gv_all, aes(x=factor(period), y=sign(voccMag)*sqrt(abs(voccMag)))) +
#   # geom_density(data=log(gv_all), aes(x=voccMag)) +
#   geom_hline(yintercept=0, colour='red', lty=2) +
#   # facet_wrap(~period) +
#   theme_bw()
# 
# 
# # ggplot() +
# #   geom_boxplot(data=gv_all, aes(x=factor(period), y=voccMag)) +
# #   # geom_density(data=log(gv_all), aes(x=voccMag)) +
# #   geom_hline(yintercept=0, colour='red', lty=2) +
# #   # facet_wrap(~period) +
# #   theme_bw()



# # Now the distance-based velocities
# # Take 1960-1970 as base period against 2000-2009
# r2 <- stack(mean(clim_rast [[1:10]], na.rm = TRUE), mean(clim_rast[[41:50]], na.rm = TRUE))
# # prepare the data frame with the necessary variables
# clim <- na.omit(data.frame(getValues(r2), cid = 1:ncell(clim_rast )))
# clim[,c("x","y")] <- xyFromCell(clim_rast, clim$cid)
# # 1965-2004 (40 yr), 500 km search radius
# v <- dVoCC(clim, n = 1, tdiff = 40, method = "Single", climTol = 0.1, geoTol = 500, distfun = "GreatCircle", trans = NA, lonlat = TRUE)
# 
# # Change sign as needed and create the distance-based velocity raster
# ind <-  which(r2[[1]][v$focal] > r2[[2]][v$target])
# v$velBis <- v$vel
# v$velBis[ind] <- v$vel[ind] * -1
# 
# # put output in raster format
# dv <- raster(gv)
# dv[v$focal] <- v$velBis
# marshift$GV <- with(marshift, raster::extract(abs(gv[[1]]), cbind(long, lat), buffer = (Shift* (timespan/10)*1000), fun = mean, na.rm = TRUE, small = TRUE))
# marshift$DV <- with(marshift, raster::extract(abs(dv), cbind(long, lat), buffer = (Shift* (timespan/10)*1000), fun = mean, na.rm = TRUE, small = TRUE))
# # retrieve the closest marine cell for those centroids falling on land
# marshift$GV <- with(marshift, ifelse(is.na(GV), gv[[1]][which.min(replace(distanceFromPoints(gv[[1]], cbind(long,lat)), is.na(gv[[1]]), NA))], GV))
# marshift$DV <- with(marshift, ifelse(is.na(DV), dv[which.min(replace(distanceFromPoints(dv, cbind(long, lat)), is.na(dv), NA))],DV))
# # fit the regression models
# Mgv <- lm(Shift^(1/4) ~ I((GV*10)^(1/4)), data = marshift, weights=years_data)
# Mdv <- lm(Shift^(1/4) ~ I((DV*10)^(1/4)), data = marshift, weights=years_data)
# summary(Mgv) 
# summary(Mdv)
# 
# 
# Produce the observed vs predicted scatterplots with regression lines (Fig. 2 in Garcia Molinos et al. 2019).
# 
# ```{r}
# # first compare both velocities
# my.at <- seq(-50, 50, by = 5)
# p1 <- rasterVis::levelplot(gv[[1]], par.settings = BuRdTheme, at=my.at, main = 'Gradient-based vocc', margin = FALSE)
# my.at <- seq(-20, 20, by = 2)
# p2 <- rasterVis::levelplot(dv, par.settings = BuRdTheme, at=my.at, main = 'Distance-based vocc', margin = FALSE)
# gridExtra::grid.arrange(p1, p2, ncol=1)
# # scatter plots with the resulting regression line
# p1 <- ggplot(na.omit(marshift), aes(x=I((GV*10)^(1/4)), y=Shift^(1/4))) + geom_point(color="grey") + geom_smooth(method=lm, se=FALSE) +
#   theme_classic() + scale_color_brewer(palette="Accent") + labs(x="Predicted shift (x^1/4; km/yr)", y = "Observed shift (y^1/4; km/yr)")
# p2 <- ggplot(na.omit(marshift), aes(x=I((DV*10)^(1/4)), y=Shift^(1/4))) + geom_point(color="grey") + geom_smooth(method=lm, se=FALSE) +
#   theme_classic() + scale_color_brewer(palette="Accent") + labs(x="Predicted shift (x^1/4; km/yr)", y = "Observed shift (y^1/4; km/yr)")
# grid.arrange(p1, p2, nrow = 1)
# ```



# 
# tmax_month_alb = projectRaster(tmax_month[[idx_keep[k]]], crs = crs(alb_proj))
# # tmax_month_alb = projectRaster(tmax_month, crs = crs(alb_proj))
# tmax_month_pm = raster::extract(tmax_month_alb, lct_paleo_spatial)

