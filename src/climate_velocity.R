library(raster)
library(VoCC)
library(ggplot2)


######################################################################################################################
##
######################################################################################################################

# read raster masks (provides masks for spatial domain/resolution of genetic + ENM data)
# stack <- stack('data/map-data/study_region_daltonIceMask_lakesMasked_linearIceSheetInterpolation.tif')
# stack[stack == 1] <- NA
saveRDS(stack, 'data/map-data/study_region_mask_glacierNA.RDS')
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

# foo = mask(clim_rast[[1]], stack[[1]])
# 
# foo = projectRaster(stack[[1]], clim_rast[[1]])
# 
# bar = mask(clim_rast[[1]], foo)

# mask out ice/water

clim_mask_stack <- list()
# for(i in 1:n_iter){
# rasterstack_list[[i]] <- stack(raster_list[[i]])
# }
for (i in 1:n_ccsm3_bins){
  print(i) 
  
  bin = ccsm3_bins[i]
  
  mask_bin_id = which.min(abs(mask_bins - bin))
  
  mask_proj = projectRaster(stack[[mask_bin_id]], clim_rast[[i]])
  clim_mask = mask(clim_rast[[i]], mask=mask_proj)
  clim_mask_stack[[i]] <- stack(clim_mask)
}

clim_mask_stack = brick(clim_mask_stack)

saveRDS(clim_mask_stack, 'data/clim_mask_stack.RDS')

######################################################################################################################
##
######################################################################################################################

# time periods: 16-14, 14-12, 12-10,10-7,7-4,4-1

# periods = data.frame(period_mid = c(15, 13, 11, 8.5, 5.5, 2.5),
#                      period_old = c(16, 14, 12, 10, 7, 4),
#                      period_young = c(14, 12, 10, 7, 4, 1))
periods = data.frame(period_mid = seq(19, 1, by=-2),
                     period_old = seq(20, 2, by=-2),
                     period_young =seq(18, 0, by=-2))
periods = periods * 1000

n_periods = nrow(periods)

gv_all = data.frame(x = numeric(0),
                    y = numeric(0),
                    voccMag = numeric(0),
                    period = numeric(0))

for (i in 1:n_periods){
  
  period_ids = which((ccsm3_bins < periods[i,'period_old']) & (ccsm3_bins > periods[i, 'period_young']))
  # period_id_old = which.min(abs(ccsm3_bins - periods[i,'period_old'] ))
  # period_id_young = which.min(abs(ccsm3_bins - periods[i,'period_young'] ))
  # period_ids = c(period_id_old, period_id_young)
  
  clim_mask_stack_sub = clim_mask_stack[[period_ids]]
  clim_period = brick(clim_mask_stack_sub)
  
  ccsm3_bins[period_ids]
  period_ends = ccsm3_bins[period_ids]
  
  # temporal trend
  # th minimum number of obs needed to calculate trend
  vt <- tempTrend(clim_period, th = 10)
  
  plot(vt)
  
  # spatial gradient
  vg <- spatGrad(clim_period, th = 0.0001, projected = FALSE)
  
  plot(vg)
  
  # climate velocity
  gv <- gVoCC(vt, vg)
  
  plot(gv)
  
  gv_df = as.data.frame(gv$voccMag, xy=TRUE)
  # gv_df$voccMag = gv_df$voccMag/abs(diff(period_ends)/10)
  # gv_df$voccMag = gv_df$voccMag/abs(diff(period_ends)/10)
  
  hist(gv_df$voccMag, breaks=100)
  
  gv_all = rbind(gv_all, 
                 data.frame(gv_df, period = periods[i, 'period_mid']))
}

gv_all[which(gv_all$voccMag>100), 'voccMag'] = 100

saveRDS(gv_all, 'output/climate_velocity.RDS')

gv_all$voccMag_binned = cut(gv_all$voccMag, breaks=c(-100, -50, -10, -1, 0, 1, 10, 50, 100, 200), labels=FALSE)

ggplot() +
  geom_histogram(data=gv_all, aes(x=sign(voccMag)*sqrt(abs(voccMag)), y=..density..), bins=30) +
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

# ggplot() +
#   geom_boxplot(data=gv_all, aes(x=factor(period), y=voccMag)) +
#   # geom_density(data=log(gv_all), aes(x=voccMag)) +
#   geom_hline(yintercept=0, colour='red', lty=2) +
#   # facet_wrap(~period) +
#   theme_bw()



# Now the distance-based velocities
# Take 1960-1970 as base period against 2000-2009
r2 <- stack(mean(clim_rast [[1:10]], na.rm = TRUE), mean(clim_rast[[41:50]], na.rm = TRUE))
# prepare the data frame with the necessary variables
clim <- na.omit(data.frame(getValues(r2), cid = 1:ncell(clim_rast )))
clim[,c("x","y")] <- xyFromCell(clim_rast, clim$cid)
# 1965-2004 (40 yr), 500 km search radius
v <- dVoCC(clim, n = 1, tdiff = 40, method = "Single", climTol = 0.1, geoTol = 500, distfun = "GreatCircle", trans = NA, lonlat = TRUE)

# Change sign as needed and create the distance-based velocity raster
ind <-  which(r2[[1]][v$focal] > r2[[2]][v$target])
v$velBis <- v$vel
v$velBis[ind] <- v$vel[ind] * -1

# put output in raster format
dv <- raster(gv)
dv[v$focal] <- v$velBis
marshift$GV <- with(marshift, raster::extract(abs(gv[[1]]), cbind(long, lat), buffer = (Shift* (timespan/10)*1000), fun = mean, na.rm = TRUE, small = TRUE))
marshift$DV <- with(marshift, raster::extract(abs(dv), cbind(long, lat), buffer = (Shift* (timespan/10)*1000), fun = mean, na.rm = TRUE, small = TRUE))
# retrieve the closest marine cell for those centroids falling on land
marshift$GV <- with(marshift, ifelse(is.na(GV), gv[[1]][which.min(replace(distanceFromPoints(gv[[1]], cbind(long,lat)), is.na(gv[[1]]), NA))], GV))
marshift$DV <- with(marshift, ifelse(is.na(DV), dv[which.min(replace(distanceFromPoints(dv, cbind(long, lat)), is.na(dv), NA))],DV))
# fit the regression models
Mgv <- lm(Shift^(1/4) ~ I((GV*10)^(1/4)), data = marshift, weights=years_data)
Mdv <- lm(Shift^(1/4) ~ I((DV*10)^(1/4)), data = marshift, weights=years_data)
summary(Mgv) 
summary(Mdv)


Produce the observed vs predicted scatterplots with regression lines (Fig. 2 in Garcia Molinos et al. 2019).

```{r}
# first compare both velocities
my.at <- seq(-50, 50, by = 5)
p1 <- rasterVis::levelplot(gv[[1]], par.settings = BuRdTheme, at=my.at, main = 'Gradient-based vocc', margin = FALSE)
my.at <- seq(-20, 20, by = 2)
p2 <- rasterVis::levelplot(dv, par.settings = BuRdTheme, at=my.at, main = 'Distance-based vocc', margin = FALSE)
gridExtra::grid.arrange(p1, p2, ncol=1)
# scatter plots with the resulting regression line
p1 <- ggplot(na.omit(marshift), aes(x=I((GV*10)^(1/4)), y=Shift^(1/4))) + geom_point(color="grey") + geom_smooth(method=lm, se=FALSE) +
  theme_classic() + scale_color_brewer(palette="Accent") + labs(x="Predicted shift (x^1/4; km/yr)", y = "Observed shift (y^1/4; km/yr)")
p2 <- ggplot(na.omit(marshift), aes(x=I((DV*10)^(1/4)), y=Shift^(1/4))) + geom_point(color="grey") + geom_smooth(method=lm, se=FALSE) +
  theme_classic() + scale_color_brewer(palette="Accent") + labs(x="Predicted shift (x^1/4; km/yr)", y = "Observed shift (y^1/4; km/yr)")
grid.arrange(p1, p2, nrow = 1)
```



# 
# tmax_month_alb = projectRaster(tmax_month[[idx_keep[k]]], crs = crs(alb_proj))
# # tmax_month_alb = projectRaster(tmax_month, crs = crs(alb_proj))
# tmax_month_pm = raster::extract(tmax_month_alb, lct_paleo_spatial)

