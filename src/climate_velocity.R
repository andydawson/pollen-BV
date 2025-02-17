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

######################################################################################################################
## reformat climate rasters
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

# clim_rast = brick('data/GCM/ccsm3_22-0k_temp.nc', level=6, varname='tmax')
# clim_rast = brick('data/GCM/ccsm3_22-0k_temp.nc', level=6, varname='tmin')
# clim_rast = brick('data/GCM/ccsm3_22-0k_prcp.nc', level=6, varname='prcp')
# clim_rast = brick('data/GCM/ccsm3_22-0k_et.nc', level=6, varname='pet')
# clim_rast = brick('data/GCM/ccsm3_22-0k_et.nc', level=6, varname='aet')

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


######################################################################################################################
## CV: by var, for each month
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
# 
# clim_rast = rast('data/clim_rast_spat.tif')
# 
# # taxon_range_mask_iter = readRDS(paste0('output/',taxa[1],'_range_mask_holocene_iter_200_v5.0.RDS'))
# taxon_range_mask = readRDS(paste0('output/',taxa[17],'_range_mask_holocene_filled_iter_200_v5.0.RDS'))
# 
# presence_mask_poly = vect('output/taxa_presence_mask_pollen_hull_poly')
# 
# proj_mask = crs(presence_mask_poly)


#clim_rast = project(clim_rast, taxon_range_mask)
# 
# clim_rast = project(clim_rast, proj_mask)

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
    
  }
  
  saveRDS(gv_all, paste0('output/climate_velocity/', var, '_CV.RDS'))
}

######################################################################################################################
## CV: mask maps using taxon-specific distributions
######################################################################################################################

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
