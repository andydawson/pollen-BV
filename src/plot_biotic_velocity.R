#
#setwd('C:/Users/abrow/Documents/pg-pollen')
require(tidyr)
require(ggplot2)
# require(rasterVis)
# require(fields)
require(rgdal)
require(raster)
require(enmSdm)
require(rgeos)
# require(sp)
# require(dplyr)
# require(holoSimCell)
# require(gridExtra)
# require(ggrepel)
# require(parallel)
library(reshape2)

######################################################################################################################
## read in metadata
######################################################################################################################

# read prediction output (50 iterations randomly sampled from all 1000 model iterations)
locs_grid <- readRDS('data/grid_4.1.RDS')
# preds <- readRDS('output/polya-gamma-predictions_4.1_overdispersed.RDS')
taxa = readRDS("data/taxa_4.1.RDS")

# read raster masks (provides masks for spatial domain/resolution of genetic + ENM data)
# stack <- stack('data/map-data/study_region_daltonIceMask_lakesMasked_linearIceSheetInterpolation.tif')
# proj <- proj4string(stack)#do some plotting
stack = readRDS('data/map-data/study_region_mask_glacierNA.RDS')

######################################################################################################################
## read in biotic velocities
######################################################################################################################

#first make a mega dataset containing all taxa
abv  = NULL #going to increment rows in a df, slow but not many to do
for (tx in 1:length(taxa))
{
  bv = readRDS(paste0('output/',taxa[tx],'_bvs_n200_v4.1.RDS'))
  bv$taxon=taxa[tx]
  abv = rbind(abv,bv)
}
saveRDS(abv, 'output/pollen_BVs_all.RDS')


######################################################################################################################
## plot biotic velocities versus time as boxplots
######################################################################################################################

pdf ("figures/biotic_velocities.pdf",width=12,height=12)

ggplot(abv, aes(x=as.factor(timeFrom), y=centroidVelocity)) + geom_boxplot() + facet_wrap(~taxon) + theme(axis.text.x = element_text(angle=90,hjust=1))

ggplot(abv, aes(x=as.factor(timeFrom), y=quantile_quant0p05)) + geom_boxplot() + facet_wrap(~taxon) + theme(axis.text.x = element_text(angle=90,hjust=1))

ggplot(abv, aes(x=as.factor(timeFrom), y=quantile_quant0p95)) + geom_boxplot() + facet_wrap(~taxon) + theme(axis.text.x = element_text(angle=90,hjust=1))

dev.off()


# just centroidVelocity, subset to exclude very high values (for visualization only)
ggplot(subset(abv, centroidVelocity<2500), aes(x=as.factor(timeFrom), y=centroidVelocity/100)) + 
  geom_boxplot() + 
  facet_wrap(~taxon) + 
  xlab('Time (kYBP)') +
  ylab('Biotic Velocity (km/decade)') +
  theme_bw(14) +
  theme(axis.text.x = element_text(angle=90,hjust=1)) +
  scale_x_discrete(labels=c(round(-abv$timeFrom/1000)))
ggsave("figures/biotic_velocities_centroid.png")

# overall summary of the centroid velocities
# across all taxa, and all times, and all iterations
summary(abv$centroidVelocity/100)

# plot histograms of the distributions of centroid velocities by taxon
# each histogram includes all times, and all iterations
pdf ("figures/biotic_velocities_centroid_histogram.pdf",width=12,height=12)
ggplot(subset(abv, centroidVelocity<2500)) + 
  geom_histogram(aes(x=centroidVelocity/100, y=..density..)) + 
  facet_wrap(~taxon) + 
  # xlab('Time (kYBP)') +
  xlab('Biotic Velocity (km/decade)') +
  theme_bw(14) #+
dev.off()

# boxplot of centroid velocities by taxon
# over all times, and all iterations, by taxon
ggplot(subset(abv, centroidVelocity<2500), aes(x=taxon, y=centroidVelocity/100)) + 
  geom_boxplot(outlier.shape=NA) + 
  xlab('Taxon') +
  ylab('Biotic Velocity (km/decade)') +
  theme_bw(14) +
  ylim(c(0,10)) +
  theme(axis.text.x = element_text(angle=45,hjust=1)) 
ggsave("figures/biotic_velocities_centroid_boxplot.png")

#  nsQuantVelocity_quant0p05
#  nsQuantVelocity_quant0p95

######################################################################################################################
## summarize centroid velocities
######################################################################################################################

# calculate the mean, median, standard deviation, 5% quantile, and 95% quantiles
# summary stats for each taxon for each time
cbv_summary = abv %>% 
  group_by(timeFrom, taxon) %>% 
  summarize(cbv_mean = mean(centroidVelocity),
            cbv_median = median(centroidVelocity),
            cbv_sd = sd(centroidVelocity),
            cbv_lo = quantile(centroidVelocity, c(0.05)),
            cbv_hi = quantile(centroidVelocity, c(0.95)), .groups="keep")

# plot median centroid biotic velocity through time
# perhaps not useful
ggplot(data=cbv_summary) +
  geom_line(aes(x=timeFrom, y=cbv_median, group=taxon, colour=taxon))


ggplot(data=subset(cbv_summary, cbv_median<2000)) +
  geom_point(aes(x=cbv_median, y=timeFrom, colour=taxon)) +
  theme_bw()

## quantify the difference between BV for each taxon for each and overall mean BV
# calculate the overall mean centroid biotic velocity
cbv_summary_all = abv %>% 
  group_by(timeFrom) %>% 
  summarize(cbv_mean = mean(centroidVelocity),
            cbv_median = median(centroidVelocity),
            cbv_sd = sd(centroidVelocity),
            cbv_lo = quantile(centroidVelocity, c(0.05)),
            cbv_hi = quantile(centroidVelocity, c(0.95)), .groups="keep")

# merge with abv biotic velocities
# add column that is difference between cBV and overall mean cBV
cbv_merge = merge(abv[,c('timeFrom', 'centroidVelocity', 'taxon')], cbv_summary_all, by=c('timeFrom'))
cbv_merge$diff = cbv_merge$centroidVelocity - cbv_merge$cbv_mean 


# difference between mean cBV for each time and cBV for that taxon for that time
ggplot(data=subset(cbv_merge, diff<2000)) +
  geom_boxplot(aes(x=factor(timeFrom), y=diff)) +
  geom_hline(yintercept = 0, colour='red', lty=2) + 
  facet_wrap(~taxon) + 
  theme_bw()

# difference between mean cBV for each time and cBV for that taxon for that time
# same as previous figure but now panels for each time
ggplot(data=subset(cbv_merge, diff<2000)) +
  geom_boxplot(aes(x=factor(taxon), y=diff)) +
  geom_hline(yintercept = 0, colour='red', lty=2) + 
  facet_wrap(~timeFrom) + 
  theme_bw()

# not useful
ggplot(data=cbv_summary_all) +
  # geom_point(aes(y=cbv_median, x=factor(timeFrom))) +
  geom_linerange(aes(x=cbv_median, y=factor(timeFrom), xmin=cbv_lo, xmax=cbv_hi))+
  # geom_vline(xintercept = 0, colour='red', lty=2) + 
  theme_bw()

# calculate the overall mean centroid biotic velocity
# plot median centroid biotic velocity through time
# perhaps not useful
ggplot(data=cbv_summary_all) +
  geom_point(aes(x=cbv_median, y=factor(timeFrom))) +
  geom_linerange(aes(y=timeFrom, xmin=cbv_lo, xmax=cbv_hi))

# summarize differences in cBV mean and taxon cBV for each time
cbv_diff = cbv_merge %>% 
  group_by(timeFrom, taxon) %>% 
  summarize(diff_lo = quantile(diff, c(0.05)),
            diff_hi = quantile(diff, c(0.95)),
            diff_median = quantile(diff, c(0.5)),
            diff_sd = sd(diff), .groups="keep")

# determine if 0 is contained in the 90% interval of differences
# if interval includes the zero, then we can not conclude thar cBV for that taxon for that time 
# is different from the cbvMean
cbv_diff$sig = rep(NA)
cbv_diff$sig = (sign(cbv_diff$diff_lo) + sign(cbv_diff$diff_hi))/2

# standard deviations of differences between cBV mean and cBV by taxon by time
# maybe not useful, just shows more uncertainty father back in time (I think)
ggplot(data=cbv_diff) +
  geom_boxplot(aes(x=factor(timeFrom), y=diff_sd)) +
  facet_wrap(~taxon)

# median differences between cBV mean and cBV by taxon by time
# maybe not useful
ggplot(data=cbv_diff) +
  geom_hline(yintercept = 0, colour='red') + 
  geom_boxplot(aes(x=factor(timeFrom), y=diff_median)) +
  facet_wrap(~taxon)

ggplot(data=cbv_diff) +
  geom_boxplot(aes(x=factor(timeFrom), y=diff_median, colour=taxon))


cbv_sig = cbv_diff %>% 
  group_by(timeFrom) %>%
  summarize(n_sig = sum(abs(sig)))

ggplot(data = cbv_sig) +
  geom_point(aes(x = timeFrom, y = n_sig)) +
  theme_bw() +
  xlab('year BP') +
  ylab('Number of taxa')



## correlation among BVs for taxa
taxa = unique(abv$taxon)
N_taxa = length(taxa)

for (i in 1:N_taxa){
  taxon = taxa[i]
  
  for (iter in 1:200){
    abv[which(abv$taxon == taxon),'iter'][((iter-1)*21 + 1):((iter-1)*21+21)] =  iter
  }
}

# foo = abv[,c('timeFrom', 'centroidVelocity', 'taxon', 'iter')] %>% 
#   group_by(timeFrom, iter) %>%
#   pivot_wider(names_from=taxon, 
#                   values_from = centroidVelocity)

foo = cbv_summary[,c('timeFrom', 'taxon', 'cbv_median')] %>% 
  group_by(timeFrom) %>%
  pivot_wider(names_from = taxon, 
              values_from = cbv_median)

# cor_all = data.frame(taxon1 = character(0), 
# taxon2 = character(0), 
# matrix(NA, nrow=0, ncol=200))

dat_sub = foo[,2:ncol(foo)]
cor_mat = cor(dat_sub)
# cor_mat[lower.tri(cor_mat)] = NA
cor_melt = melt(cor_mat)
cor_all = cor_melt[!is.na(cor_melt$value),]
colnames(cor_all) = c('taxon1', 'taxon2', 'cor')

# for (iter in 2:200){
#   dat_sub = foo[which(foo$iter == iter),3:ncol(foo)]
#   cor_mat = cor(dat_sub)
#   cor_mat[lower.tri(cor_mat)] = NA
#   cor_melt = melt(cor_mat)
#   cor_melt = cor_melt[!is.na(cor_melt$value),]
#   cor_all[,paste0('X', as.character(iter))] = cor_melt$value
# }
# 
# cor_means = data.frame(cor_all[,1:2], cor = rowMeans(cor_all[,3:ncol(cor_all)]))
ggplot() +
  geom_tile(data=cor_all, aes(x=taxon1, y=taxon2, fill=cor)) +
  theme_bw()


breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)
cor_all$bin = cut(cor_all$cor, breaks, labels=FALSE)

summary(cor_all)

cor_high = cor_all[which((cor_all$cor >= 0.4)&(cor_all$cor<1)),]


ggplot() +
  geom_point(data=cor_high, aes(x=taxon1, y=taxon2, colour=factor(bin), size=cor))

ggplot() +
  geom_point(data=cor_high, aes(x=taxon1, y=taxon2, colour=factor(bin), size=factor(bin)))

ggplot() +
  geom_point(data=cor_high, aes(x=taxon1, y=taxon2, colour=cor, size=cor)) 

ggplot() +
  geom_tile(data=cor_high, aes(x=taxon1, y=taxon2, fill=factor(bin)))

ggplot() +
  geom_point(data=cor_high, aes(x=taxon1, y=taxon2)) 




# abv %>% 
#   group_by(taxon, iter) %>%
#   mutate()
# 
# 
# ######################################################################################################################
# ## northern range quantile
# ######################################################################################################################
# 
# library(dplyr)
# abv_nq = abv %>% 
#   group_by(timeFrom, taxon) %>% 
#   summarize(nbv_mean = mean(nsQuantVelocity_quant0p95),
#             nbv_median = median(nsQuantVelocity_quant0p95),
#             nbv_sd = sd(nsQuantVelocity_quant0p95),
#             nbv_lo = quantile(nsQuantVelocity_quant0p95, c(0.05)),
#             nbv_hi = quantile(nsQuantVelocity_quant0p95, c(0.95)), .groups="keep")
# 
# 
# ggplot(data=abv_nq) +
#   geom_boxplot(aes(x=factor(timeFrom), y=nbv_mean))
# 
# ggplot(data=abv_nq) +
#   geom_boxplot(aes(x=factor(timeFrom), y=nbv_median))
# 
# ggplot(data=abv_nq) +
#   geom_boxplot(aes(x=factor(timeFrom), y=nbv_sd))
# 
# 
# 
# abv_nq$sig = rep(NA)
# abv_nq$sig = (sign(abv_nq$nbv_lo) + sign(abv_nq$nbv_hi))/2
# 
# nq_sig = abv_nq %>% 
#   group_by(timeFrom) %>%
#   summarize(n_sig = sum(abs(sig)))
# 
# ggplot(data = nq_sig) +
#   geom_point(aes(x = timeFrom, y = n_sig)) +
#   theme_bw() +
#   xlab('year BP') +
#   ylab('Number of taxa')
# 
# 
# 
# ######################################################################################################################
# ## southern range quantile
# ######################################################################################################################
# 
# abv_sq = abv %>% 
#   group_by(timeFrom, taxon) %>% 
#   summarize(bv_mean = mean(nsQuantVelocity_quant0p05),
#             bv_median = median(nsQuantVelocity_quant0p05),
#             bv_sd = sd(nsQuantVelocity_quant0p05),
#             bv_lo = quantile(nsQuantVelocity_quant0p05, c(0.05)),
#             bv_hi = quantile(nsQuantVelocity_quant0p05, c(0.95)), .groups="keep")
# 
# 
# abv_sq$sig = rep(NA)
# abv_sq$sig = (sign(abv_sq$bv_lo) + sign(abv_sq$bv_hi))/2
# 
# sq_sig = abv_sq %>% 
#   group_by(timeFrom) %>%
#   summarize(n_sig = sum(abs(sig)))
# 
# ggplot(data = sq_sig) +
#   geom_point(aes(x = timeFrom, y = n_sig)) +
#   theme_bw() +
#   xlab('year BP') +
#   ylab('Number of taxa')
# 
# ggplot(data=abv_sq) +
#   geom_boxplot(aes(x=factor(timeFrom), y=bv_sd))
# 
# 
# ######################################################################################################################
# ## all bv sig
# ######################################################################################################################
# 
# 
# bv_all = rbind(data.frame(nq_sig, type='nq'),
#                data.frame(sq_sig, type='sq'))
# bv_all = rbind(bv_all,
#                data.frame(cbv_sig, type='c'))
# 
# ggplot(data = bv_all) +
#   geom_point(aes(x = timeFrom, y = n_sig, colour=type)) +
#   theme_bw() +
#   xlab('year BP') +
#   ylab('Number of taxa')
# 
# ######################################################################################################################
# ## centroid summary
# ######################################################################################################################
# 
# ggplot(data=abv) +
#   geom_histogram(aes(x=centroidVelocity)) +
#   facet_wrap(~taxon)
# 
# 
# 
# foo3 = abv %>% 
#   group_by(timeFrom, taxon) %>% 
#   summarize(cbv_mean = mean(centroidVelocity),
#             cbv_median = median(centroidVelocity),
#             cbv_sd = sd(centroidVelocity),
#             cbv_lo = quantile(centroidVelocity, c(0.05)),
#             cbv_hi = quantile(centroidVelocity, c(0.95)), .groups="keep")
# 
# ggplot(data=foo3) +
#   geom_point(aes(x=timeFrom, y=cbv_median, group=taxon, colour=taxon))
# 
# 
# cbv_all = abv %>% 
#   group_by(timeFrom) %>% 
#   summarize(cbv_mean_all = mean(centroidVelocity), .groups="keep")