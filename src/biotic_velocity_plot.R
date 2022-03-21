#setwd('C:/Users/abrow/Documents/pg-pollen')
# require(tidyr)
require(ggplot2)
# require(rasterVis)
# require(fields)
# require(rgdal)
require(raster)
# require(enmSdm)
# require(rgeos)
# require(sp)
# require(dplyr)
# require(holoSimCell)
# require(gridExtra)
# require(ggrepel)
# require(parallel)

######################################################################################################################
## CALCULATE AND PLOT BIOTIC VELOCITIES
######################################################################################################################

# read prediction output (50 iterations randomly sampled from all 1000 model iterations)
locs_grid <- readRDS('data/grid_4.1.RDS')
# preds <- readRDS('output/polya-gamma-predictions_4.1_overdispersed.RDS')
taxa = readRDS("data/taxa_4.1.RDS")

# read raster masks (provides masks for spatial domain/resolution of genetic + ENM data)
stack <- stack('data/map-data/study_region_daltonIceMask_lakesMasked_linearIceSheetInterpolation.tif')
proj <- proj4string(stack)#do some plotting
#first make a mega dataset containing all taxa
abv  = NULL #going to increment rows in a df, slow but not many to do
for (tx in 1:length(taxa))
{
  bv = readRDS(paste0('output/',taxa[tx],'_bvs_n200_v4.1.RDS'))
  bv$taxon=taxa[tx]
  abv = rbind(abv,bv)
}


pdf ("figures/biotic_velocities.pdf",width=12,height=12)

ggplot(abv, aes(x=as.factor(timeFrom), y=centroidVelocity)) + geom_boxplot() + facet_wrap(~taxon) + theme(axis.text.x = element_text(angle=90,hjust=1))


ggplot(abv, aes(x=as.factor(timeFrom), y=quantile_quant0p05)) + geom_boxplot() + facet_wrap(~taxon) + theme(axis.text.x = element_text(angle=90,hjust=1))

ggplot(abv, aes(x=as.factor(timeFrom), y=quantile_quant0p95)) + geom_boxplot() + facet_wrap(~taxon) + theme(axis.text.x = element_text(angle=90,hjust=1))

dev.off()

#  nsQuantVelocity_quant0p05
#  nsQuantVelocity_quant0p95

######################################################################################################################
## northern range quantile
######################################################################################################################

library(dplyr)
abv_nq = abv %>% 
  group_by(timeFrom, taxon) %>% 
  summarize(nbv_mean = mean(nsQuantVelocity_quant0p95),
            nbv_median = median(nsQuantVelocity_quant0p95),
            nbv_sd = sd(nsQuantVelocity_quant0p95),
            nbv_lo = quantile(nsQuantVelocity_quant0p95, c(0.05)),
            nbv_hi = quantile(nsQuantVelocity_quant0p95, c(0.95)), .groups="keep")


abv_nq$sig = rep(NA)
abv_nq$sig = (sign(abv_nq$nbv_lo) + sign(abv_nq$nbv_hi))/2

nq_sig = abv_nq %>% 
  group_by(timeFrom) %>%
  summarize(n_sig = sum(abs(sig)))

ggplot(data = nq_sig) +
  geom_point(aes(x = timeFrom, y = n_sig)) +
  theme_bw() +
  xlab('year BP') +
  ylab('Number of taxa')


ggplot(data=abv_nq) +
  geom_boxplot(aes(x=factor(timeFrom), y=nbv_sd))

######################################################################################################################
## southern range quantile
######################################################################################################################

abv_sq = abv %>% 
  group_by(timeFrom, taxon) %>% 
  summarize(bv_mean = mean(nsQuantVelocity_quant0p05),
            bv_median = median(nsQuantVelocity_quant0p05),
            bv_sd = sd(nsQuantVelocity_quant0p05),
            bv_lo = quantile(nsQuantVelocity_quant0p05, c(0.05)),
            bv_hi = quantile(nsQuantVelocity_quant0p05, c(0.95)), .groups="keep")


abv_sq$sig = rep(NA)
abv_sq$sig = (sign(abv_sq$bv_lo) + sign(abv_sq$bv_hi))/2

sq_sig = abv_sq %>% 
  group_by(timeFrom) %>%
  summarize(n_sig = sum(abs(sig)))

ggplot(data = sq_sig) +
  geom_point(aes(x = timeFrom, y = n_sig)) +
  theme_bw() +
  xlab('year BP') +
  ylab('Number of taxa')

ggplot(data=abv_sq) +
  geom_boxplot(aes(x=factor(timeFrom), y=bv_sd))

######################################################################################################################
## centroid
######################################################################################################################

foo3 = abv %>% 
  group_by(timeFrom, taxon) %>% 
  summarize(cbv_mean = mean(centroidVelocity),
            cbv_median = median(centroidVelocity),
            cbv_sd = sd(centroidVelocity),
            cbv_lo = quantile(centroidVelocity, c(0.05)),
            cbv_hi = quantile(centroidVelocity, c(0.95)), .groups="keep")

ggplot(data=foo3) +
  geom_line(aes(x=timeFrom, y=cbv_median, group=taxon, colour=taxon))


cbv_all = abv %>% 
  group_by(timeFrom) %>% 
  summarize(cbv_mean_all = mean(centroidVelocity), .groups="keep")

cbv_merge = merge(abv[,c('timeFrom', 'centroidVelocity', 'taxon')], cbv_all, by=c('timeFrom'))

# foo_merge = merge(foo3, foo_all, by=c('timeFrom'))

cbv_merge$diff = cbv_merge$centroidVelocity - cbv_merge$cbv_mean_all 

# foo3$total = rep(NA)
# foo3$total = foo3$nbv_lo + foo3$nbv_hi

# ggplot(data=foo_merge) +
#   geom_point(aes(x = timeFrom, y = diff, colour = taxon))
# 
cbv_diff = cbv_merge %>% 
  group_by(timeFrom, taxon) %>% 
  summarize(diff_lo = quantile(diff, c(0.05)),
            diff_hi = quantile(diff, c(0.95)),
            diff_sd = sd(diff), .groups="keep")

cbv_diff$sig = rep(NA)
cbv_diff$sig = (sign(cbv_diff$diff_lo) + sign(cbv_diff$diff_hi))/2

# bar3 = foo3_merge %>%
#   group_by(timeFrom) %>%
#   summarize(total = mean(nbv_median))
# 
# ggplot(data = bar3) +
#   geom_point(aes(x = timeFrom, y = total)) +
#   theme_bw() +
#   xlab('year BP') +
#   ylab('Total movement')


cbv_sig = cbv_diff %>% 
  group_by(timeFrom) %>%
  summarize(n_sig = sum(abs(sig)))

ggplot(data = cbv_sig) +
  geom_point(aes(x = timeFrom, y = n_sig)) +
  theme_bw() +
  xlab('year BP') +
  ylab('Number of taxa')

ggplot(data=cbv_diff) +
  geom_boxplot(aes(x=factor(timeFrom), y=diff_sd))


######################################################################################################################
## all bv sig
######################################################################################################################


bv_all = rbind(data.frame(nq_sig, type='nq'),
data.frame(sq_sig, type='sq'))
bv_all = rbind(bv_all,
               data.frame(cbv_sig, type='c'))

ggplot(data = bv_all) +
  geom_point(aes(x = timeFrom, y = n_sig, colour=type)) +
  theme_bw() +
  xlab('year BP') +
  ylab('Number of taxa')

