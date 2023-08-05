library(dplyr)
library(ggplot2)
library(reshape2)
library(broom)
library(tidyr)
library(corrplot)
library(ggcorrplot)
library(FactoMineR)
library(factoextra)


######################################################################################################################
## read and summarize biotic velocities
######################################################################################################################

pollen_BVs = readRDS('output/pollen_BVs_all_v5.0.RDS')

pollen_BVs_summary = pollen_BVs %>% 
  group_by(timeFrom, timeTo, timeMid, taxon) %>% 
  summarize(bvc_mean = mean(centroidVelocity),
            bvc_median = median(centroidVelocity),
            bvc_sd = sd(centroidVelocity),
            bvc_lo = quantile(centroidVelocity, c(0.05)),
            bvc_hi = quantile(centroidVelocity, c(0.95)), 
            bvn_mean = mean(nsQuantVelocity_quant0p95),
            bvn_median = median(nsQuantVelocity_quant0p95),
            bvn_sd = sd(nsQuantVelocity_quant0p95),
            bvn_lo = quantile(nsQuantVelocity_quant0p95, c(0.05)),
            bvn_hi = quantile(nsQuantVelocity_quant0p95, c(0.95)), 
            bvs_mean = mean(nsQuantVelocity_quant0p05),
            bvs_median = median(nsQuantVelocity_quant0p05),
            bvs_sd = sd(nsQuantVelocity_quant0p05),
            bvs_lo = quantile(nsQuantVelocity_quant0p05, c(0.05)),
            bvs_hi = quantile(nsQuantVelocity_quant0p05, c(0.95)), 
            .groups="keep") %>%
  mutate(bvn_sig = ifelse(sign(bvn_lo) == sign(bvn_hi), 1, 0),
         bvs_sig = ifelse(sign(bvs_lo) == sign(bvs_hi), 1, 0))

# pollen_BVs_summary$period = abs(pollen_BVs_summary$timeFrom + 990/2)

pollen_BVs_summary$bvn_sig = factor(pollen_BVs_summary$bvn_sig,
                                    levels = c(0, 1),
                                    labels = c('Not sig', 'Sig')) 
pollen_BVs_summary$bvs_sig = factor(pollen_BVs_summary$bvs_sig,
                                    levels = c(0, 1),
                                    labels = c('Not sig', 'Sig')) 

head(pollen_BVs_summary)

dim(pollen_BVs_summary)

pollen_BVs_summary$period = abs(pollen_BVs_summary$timeMid)

# pollen_BVs_summary$period = floor(pollen_BVs_summary$period/1000)*1000

time_labels = unique(pollen_BVs_summary$timeFrom)
time_period = unique(pollen_BVs_summary$timeMid)

ggplot(data=pollen_BVs_summary) +
  geom_point(aes(x=timeFrom, y=bvc_median)) +
  geom_linerange(aes(x=timeFrom, ymin=bvc_lo, ymax=bvc_hi)) +
  facet_wrap(~taxon) +
  ylab('Biotic velocity (m/year)') +
  xlab('Time (YBP)') +
  scale_x_continuous(name = 'Time (YBP)', breaks = time_labels, labels = time_labels) +
  theme_bw(14) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=8)) 
ggsave('figures/BV_centroid_quantiles_fixedy.pdf', width = 14, height=8)
ggsave('figures/BV_centroid_quantiles_fixedy.png', width = 14, height=8)

ggplot(data=pollen_BVs_summary) +
  geom_point(aes(x=timeFrom, y=bvc_median)) +
  geom_linerange(aes(x=timeFrom, ymin=bvc_lo, ymax=bvc_hi)) +
  facet_wrap(~taxon, scales = 'free_y') +
  ylab('Biotic velocity (m/year)') +
  xlab('Time (YBP)') +
  theme_bw(14) +
  scale_x_continuous(name = 'Time (YBP)', breaks = time_labels, labels = time_labels) +
  theme_bw(14) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=7.5)) 
ggsave('figures/BV_centroid_quantiles_freey.pdf', width = 14, height=8)
ggsave('figures/BV_centroid_quantiles_freey.png', width = 14, height=8)

ggplot(data=pollen_BVs_summary) +
  geom_hline(yintercept = 0, colour='darkgrey') +
  geom_point(aes(x=timeFrom, y=bvn_median, colour=bvn_sig)) +
  geom_linerange(aes(x=timeFrom, ymin=bvn_lo, ymax=bvn_hi, colour=bvn_sig)) +
  facet_wrap(~taxon) +
  ylab('Biotic velocity (m/year)') +
  xlab('Time (YBP)') +
  scale_x_continuous(name = 'Time (YBP)', breaks = time_labels, labels = time_labels) +
  theme_bw(14) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=7.5)) 
ggsave('figures/BV_north_quantiles_fixedy.pdf', width = 14, height=8)
ggsave('figures/BV_north_quantiles_fixedy.png', width = 14, height=8)

ggplot(data=pollen_BVs_summary) +
  geom_hline(yintercept = 0, colour='darkgrey') +
  geom_point(aes(x=timeFrom, y=bvn_median, colour=bvn_sig)) +
  geom_linerange(aes(x=timeFrom, ymin=bvn_lo, ymax=bvn_hi, colour=bvn_sig)) +
  facet_wrap(~taxon, scales = 'free_y') +
  ylab('Biotic velocity (m/year)') +
  xlab('Time (YBP)') +
  scale_x_continuous(name = 'Time (YBP)', breaks = time_labels, labels = time_labels) +
  theme_bw(14) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=7.5)) 
ggsave('figures/BV_north_quantiles_freey.pdf', width = 14, height=8)
ggsave('figures/BV_north_quantiles_freey.png', width = 14, height=8)

ggplot(data=pollen_BVs_summary) +
  geom_hline(yintercept = 0, colour='darkgrey') +
  geom_point(aes(x=timeFrom, y=bvs_median, colour=bvs_sig)) +
  geom_linerange(aes(x=timeFrom, ymin=bvs_lo, ymax=bvs_hi, colour=bvs_sig)) +
  facet_wrap(~taxon) +
  ylab('Biotic velocity (m/year)') +
  xlab('Time (YBP)') +
  scale_x_continuous(name = 'Time (YBP)', breaks = time_labels, labels = time_labels) +
  theme_bw(14) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=7.5)) 
ggsave('figures/BV_south_quantiles_fixedy.pdf', width = 14, height=8)
ggsave('figures/BV_south_quantiles_fixedy.png', width = 14, height=8)

ggplot(data=pollen_BVs_summary) +
  geom_hline(yintercept = 0, colour='darkgrey') +
  geom_point(aes(x=timeFrom, y=bvs_median, colour=bvs_sig)) +
  geom_linerange(aes(x=timeFrom, ymin=bvs_lo, ymax=bvs_hi, colour=bvs_sig)) +
  facet_wrap(~taxon, scales = 'free_y') +
  ylab('Biotic velocity (m/year)') +
  xlab('Time (YBP)') +
  scale_x_continuous(name = 'Time (YBP)', breaks = time_labels, labels = time_labels) +
  theme_bw(14) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=7.5)) 
ggsave('figures/BV_south_quantiles_freey.pdf', width = 14, height=8)
ggsave('figures/BV_south_quantiles_freey.png', width = 14, height=8)

######################################################################################################################
## correlation of taxon-specific biotic velocity time series
######################################################################################################################

BVC_taxa_wide = pollen_BVs_summary[,colnames(pollen_BVs_summary) %in% 
                                     c('timeFrom', 'taxon', 'bvc_mean')] %>% 
  pivot_wider(id_cols = 'timeFrom', values_from = 'bvc_mean', names_from = 'taxon') 

BVC_cor = cor(BVC_taxa_wide[,2:ncol(BVC_taxa_wide)])

corrplot(BVC_cor)

BVC_cor_p = rcorr(as.matrix(BVC_taxa_wide[,2:ncol(BVC_taxa_wide)]))
BVC_r = BVC_cor_p$r
BVC_p = BVC_cor_p$P
BVC_p[which(is.na(BVC_p))] = 0

png('figures/BVC_corr_plot.png')#, width = 14, height=8, units='in')
corrplot(BVC_r, 
         order = 'FPC', 
         p.mat = BVC_p, 
         insig = 'blank',
         type = 'lower')
dev.off()

BVN_taxa_wide = pollen_BVs_summary[,colnames(pollen_BVs_summary) %in% c('timeFrom', 'taxon', 'bvn_mean')] %>% 
  pivot_wider(id_cols = 'timeFrom', values_from = 'bvn_mean', names_from = 'taxon') 

BVN_cor = cor(BVN_taxa_wide[,2:ncol(BVN_taxa_wide)])

corrplot(BVN_cor)


BVN_cor_p = rcorr(as.matrix(BVN_taxa_wide[,2:ncol(BVN_taxa_wide)]))
BVN_r = BVN_cor_p$r
BVN_p = BVN_cor_p$P
BVN_p[which(is.na(BVN_p))] = 0

png('figures/BVN_corr_plot.png')#, width = 14, height=8)
corrplot(BVN_r, 
         order = 'FPC', 
         p.mat = BVN_p, 
         insig = 'blank',
         type = 'lower')
dev.off()

BVS_taxa_wide = pollen_BVs_summary[,colnames(pollen_BVs_summary) %in% c('timeFrom', 'taxon', 'bvs_mean')] %>% 
  pivot_wider(id_cols = 'timeFrom', values_from = 'bvs_mean', names_from = 'taxon') 

BVS_cor = cor(BVS_taxa_wide[,2:ncol(BVS_taxa_wide)])

corrplot(BVS_cor)


BVS_cor_p = rcorr(as.matrix(BVS_taxa_wide[,2:ncol(BVS_taxa_wide)]))
BVS_r = BVS_cor_p$r
BVS_p = BVS_cor_p$P
BVS_p[which(is.na(BVS_p))] = 0

png('figures/BVS_corr_plot.png')#, width = 14, height=8)
corrplot(BVS_r, 
         order = 'FPC', 
         p.mat = BVS_p, 
         insig = 'blank',
         type = 'lower')
dev.off()

# foo = melt(BVC_r)
# colnames(foo) = c('taxon1', 'taxon2', 'cor')
# bar = melt(BVC_p)
# colnames(bar) = c('taxon1', 'taxon2', 'p')
# 
# foobar = merge(foo, bar)
# 
# ggcorrplot(BVC_r, 
#            hc.order = TRUE,
#            type = "lower",
#            p.mat = BVC_p,
#            insig = "blank")
#   
# ggcorrplot(BVC_r, 
#            hc.order = TRUE,
#            type = "lower",
#            p.mat = BVC_p)

# 
# pollen_BVs_summary[,colnames(pollen_BVs_summary) %in% c('timeFrom', 'taxon', 'bvc_mean')] %>%
#   pairwise_cor(taxon, timeFrom, bvc_mean)

# BVC_BVS_cor = pollen_BVs_summary %>%
#   group_by(taxon) %>%
#   summarize(cor = cor(bvc_mean, bvs_mean))
# 
# BVN_BVS_cor = pollen_BVs_summary %>%
#   group_by(taxon) %>%
#   summarize(cor = cor(bvn_mean, bvs_mean))
# 
# 
# cor_mat = cor(dat_sub)
# # cor_mat[lower.tri(cor_mat)] = NA
# cor_melt = melt(cor_mat)
# cor_all = cor_melt[!is.na(cor_melt$value),]
# colnames(cor_all) = c('taxon1', 'taxon2', 'X1')

######################################################################################################################
## principal component analysis
######################################################################################################################

# CENTROID
BVC_wide = BVC_taxa_wide[,2:ncol(BVC_taxa_wide)]

BVC_normalized = scale(BVC_wide)

BVC_pca = prcomp(BVC_normalized)
BVC_pca$x[,'PC1']

foo = princomp(BVC_normalized)

# 
# BVN_pca$loadings[,1:3]
# BVN_pca$scores

eigs <- BVC_pca$sdev^2
rbind(SD = sqrt(eigs),
  Proportion = eigs/sum(eigs),
  Cumulative = cumsum(eigs)/sum(eigs))

png('figures/BVC_scree_plot.png')#, width = 14, height=8)
fviz_eig(BVC_pca, addlabels = TRUE)
dev.off()

fviz_pca_var(BVC_pca, col.var = "black")
fviz_cos2(BVC_pca, choice = "var", axes = 1:2)
BVC_pca_time = data.frame(period=time_period, BVC_normalized %*% BVC_pca$rotation)

pca_time_melt = melt(BVC_pca_time, id.vars='period')
pca_time_melt$type = rep('BVC')

# NORTH
BVN_pca = prcomp(BVN_normalized)
BVN_pca$x[,'PC1']
# 
# BVN_pca$loadings[,1:3]
# BVN_pca$scores

png('figures/BVN_scree_plot.png')#, width = 14, height=8)
fviz_eig(BVN_pca, addlabels = TRUE)
dev.off()

fviz_pca_var(BVN_pca, col.var = "black")
fviz_cos2(BVN_pca, choice = "var", axes = 1:2)
BVN_pca_time = data.frame(period=time_period, BVN_normalized %*% BVN_pca$rotation)

pca_time_melt = rbind(pca_time_melt, data.frame(melt(BVN_pca_time, id.vars='period'), type = rep('BVN')))

# SOUTH
BVS_wide = BVS_taxa_wide[,2:ncol(BVS_taxa_wide)]

BVS_normalized = scale(BVS_wide)
BVS_pca = prcomp(BVS_normalized)
BVS_pca$x[,'PC1']
# 
# BVN_pca$loadings[,1:3]
# BVN_pca$scores

png('figures/BVS_scree_plot.png')#, width = 14, height=8)
fviz_eig(BVS_pca, addlabels = TRUE)
dev.off()

fviz_pca_var(BVS_pca, col.var = "black")
fviz_cos2(BVS_pca, choice = "var", axes = 1:2)
BVS_pca_time = data.frame(period=time_period, BVS_normalized %*% BVS_pca$rotation)

pca_time_melt = rbind(pca_time_melt, data.frame(melt(BVS_pca_time, id.vars='period'), type = rep('BVS')))


# 
# ggplot(data=subset(pca_time_melt, variable %in% c('PC1', 'PC2'))) +
#   geom_point(aes(x=period, y=value, colour=variable)) +
#   geom_line(aes(x=period,  y=value, colour=variable)) +
#   ylab('Principal component') +
#   xlab('Time (YBP)') +
#   scale_x_continuous(name = 'Time (YBP)', breaks = time_labels, labels = time_labels) +
#   theme_bw(14) +
#   theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=8)) 
# ggsave('figures/BVC_PCA_time_series.pdf', width = 14, height=8)
# ggsave('figures/BVC_PCA_time_series.pdf', width = 14, height=8)

climate_periods = data.frame(time_from = c(21, 18.5, 14.7, 12.9)*(-1000), 
                             time_to = c(18.5, 14.7, 12.9, 11.7)*(-1000),
                             period_name = c('LGM', 'Oldest Dryas', 'Bolling-Allerod', 'Younger Dryas'))
climate_periods$period_name = factor(climate_periods$period_name,
                                     levels = c('LGM', 'Oldest Dryas', 'Bolling-Allerod', 'Younger Dryas'))

PC_first_two = subset(pca_time_melt, variable %in% c('PC1', 'PC2'))

ggplot(data=PC_first_two) +
  geom_rect(data=climate_periods, inherit.aes=FALSE, aes(xmin=time_from, 
                                                         xmax=time_to, 
                                                         ymin=min(PC_first_two$value),
                                                         ymax=max(PC_first_two$value),
                                                         group=period_name,
                                                         fill=period_name),
            colour='transparent',
            alpha=0.3) +
  geom_point(aes(x=period, y=value, colour=type, shape=type), size=3) +
  geom_line(aes(x=period,  y=value, colour=type), linewidth=1) +
  ylab('Principal component') +
  xlab('Time (YBP)') +
  scale_x_continuous(name = 'Time (YBP)', breaks = time_labels, labels = time_labels) +
  theme_bw(20) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=16)) +
  facet_grid(variable~.) 
 #+
  # geom_text(data=climate_periods, 
  #   aes(x = time_from, y = 5.5, label = period_name),
  #   size = 3, vjust = 0, hjust = 0, nudge_x = 50, check_overlap = TRUE
  # )
ggsave('figures/BV_PCA_time_series.pdf', width = 14, height=10)
ggsave('figures/BV_PCA_time_series.png', width = 14, height=10)


######################################################################################################################
## climate velocity
######################################################################################################################

CVs = readRDS('output/climate_velocity.RDS')

CVs_summary = CVs %>% 
  group_by(period) %>% 
  summarize(cv_mean = mean(voccMag, na.rm=TRUE),
            cv_median = median(voccMag, na.rm=TRUE),
            cv_sd = sd(voccMag, na.rm=TRUE),
            cv_lo = quantile(voccMag, c(0.05), na.rm=TRUE),
            cv_hi = quantile(voccMag, c(0.95), na.rm=TRUE), .groups="keep")

ggplot(data=CVs_summary) +
  geom_point(aes(x=period, y=cv_median)) +
  geom_linerange(aes(x=period, ymin=cv_lo, ymax=cv_hi)) +
  ylab('Climate velocity (m/year)') +
  xlab('Time (YBP)') +
  theme_bw(14) +
  theme(axis.text.x = element_text(angle=90,hjust=1)) 


ggplot(data=CVs) +
  geom_boxplot(aes(x=factor(period), y=voccMag), outlier.shape = NA) +
  # geom_linerange(aes(x=period, ymin=cv_lo, ymax=cv_hi)) +
  ylab('Climate velocity (km/decade)') +
  xlab('Time (YBP)') +
  ylim(c(-3,3)) +
  theme_bw(16) +
  theme(axis.text.x = element_text(angle=45,hjust=1, size=14)) +
  geom_hline(yintercept=0, colour='grey', lty=2, lwd=2) +
  scale_x_discrete(limits=rev)
ggsave('figures/CV_vs_time_boxplot.pdf', width=10, height=6)
ggsave('figures/CV_vs_time_boxplot.png', width=10, height=6)


PC_first_two$period = -(PC_first_two$period)

PC_CV = merge(PC_first_two, CVs_summary, by.x = 'period', by.y = 'period')

PC_CV_wide = PC_CV[,colnames(PC_CV) %in% c('period', 'variable', 'type', 'value', 'cv_mean')] %>% 
  pivot_wider(id_cols = c('period', 'cv_mean'), values_from = 'value', names_from = c('variable', 'type')) 

PC_CV_cor = cor(PC_CV_wide[,2:ncol(PC_CV_wide)])

corrplot(BVS_cor)


BVS_cor_p = rcorr(as.matrix(BVS_taxa_wide[,2:ncol(BVS_taxa_wide)]))
BVS_r = BVS_cor_p$r
BVS_p = BVS_cor_p$P
BVS_p[which(is.na(BVS_p))] = 0

corrplot(BVS_r, 
         order = 'FPC', 
         p.mat = BVS_p, 
         insig = 'blank',
         type = 'lower')

# BVs between two consecutive time bins 990 years apart
# 
Vs = merge(pollen_BVs_summary, CVs_summary, by.x = 'period', by.y = 'period')

head(Vs)


# ggplot() +
#   geom_point(data=Vs, aes(x=cv_mean, y= bv_mean/100, colour=taxon)) +
#   geom_smooth(data=Vs, aes(x=cv_mean, y= bv_mean/100, colour=taxon), method=lm) +
#   # facet_wrap(~taxon) +
#   theme_bw(14) +
#   xlab('Climate velocity (km/decade)') +
#   ylab('Biotic velocity (km/decade)')
# ggsave('figures/BV_vs_CV_mean_scatter.pdf', width=10, height=8)
# 
# ggplot() +
#   geom_point(data=Vs, aes(x=cv_median, y= bv_median/100)) +
#   geom_smooth(data=Vs, aes(x=cv_median, y= bv_median/100), method=lm) +
#   facet_wrap(~taxon) +
#   theme_bw(14) +
#   xlab('Climate velocity (km/decade)') +
#   ylab('Biotic velocity (km/decade)')
# ggsave('figures/BV_vs_CV_median_scatter_taxon.pdf', width=10, height=8)

# ggplot(data=Vs) +
#   geom_point(aes(x=cv_mean, y= bvc_mean/100)) +
#   geom_smooth(aes(x=cv_mean, y= bvc_mean/100), method=lm) +
#   facet_wrap(~taxon) +
#   theme_bw(14) +
#   xlab('Climate velocity (km/decade)') +
#   ylab('Biotic velocity (km/decade)')
# ggsave('figures/BVC_vs_CV_mean_scatter_taxon.pdf', width=10, height=8)


ggplot() +
  geom_point(data=Vs, aes(x=cv_mean, y= bvc_mean/100)) +
  geom_smooth(data=Vs, aes(x=cv_mean, y= bvc_mean/100), method=lm) +
  facet_wrap(~taxon, scales = 'free_y') +
  theme_bw(14) +
  xlab('Climate velocity (km/decade)') +
  ylab('Biotic velocity (km/decade)')
ggsave('figures/BVC_vs_CV_mean_scatter_taxon.pdf', width=10, height=8)
ggsave('figures/BVC_vs_CV_mean_scatter_taxon.png', width=10, height=8)


ggplot(data=Vs) +
  geom_point(aes(x=cv_mean, y= bvc_mean/100)) +
  geom_linerange(aes(x=cv_mean, ymin=bvc_lo/100, ymax=bvc_hi/100)) +
  geom_smooth(aes(x=cv_mean, y= bvc_mean/100), method=lm) +
  facet_wrap(~taxon, scales = 'free_y') +
  theme_bw(14) +
  xlab('Climate velocity (km/decade)') +
  ylab('Biotic velocity (km/decade)')
ggsave('figures/BVC_vs_CV_mean_scatter_credible_taxon.pdf', width=10, height=8)
ggsave('figures/BVC_vs_CV_mean_scatter_credible_taxon.png', width=10, height=8)


ggplot(data=Vs) +
  geom_point(aes(x=cv_mean, y= bvn_mean/100)) +
  geom_linerange(aes(x=cv_mean, ymin=bvn_lo/100, ymax=bvn_hi/100)) +
  geom_smooth(aes(x=cv_mean, y= bvn_mean/100), method=lm) +
  facet_wrap(~taxon, scales = 'free_y') +
  theme_bw(14) +
  xlab('Climate velocity (km/decade)') +
  ylab('Biotic velocity (km/decade)')
ggsave('figures/BVN_vs_CV_mean_scatter_credible_taxon.pdf', width=10, height=8)
ggsave('figures/BVN_vs_CV_mean_scatter_credible_taxon.png', width=10, height=8)


ggplot(data=Vs) +
  geom_point(aes(x=cv_mean, y= bvs_mean/100)) +
  geom_linerange(aes(x=cv_mean, ymin=bvs_lo/100, ymax=bvs_hi/100)) +
  geom_smooth(aes(x=cv_mean, y= bvs_mean/100), method=lm) +
  facet_wrap(~taxon, scales = 'free_y') +
  theme_bw(14) +
  xlab('Climate velocity (km/decade)') +
  ylab('Biotic velocity (km/decade)')
ggsave('figures/BVS_vs_CV_mean_scatter_credible_taxon.pdf', width=10, height=8)
ggsave('figures/BVS_vs_CV_mean_scatter_credible_taxon.png', width=10, height=8)


ggplot(data=Vs) +
  geom_point(aes(x=bvc_mean, y= bvn_mean/100)) +
  geom_linerange(aes(x=bvc_mean, ymin=bvn_lo/100, ymax=bvn_hi/100)) +
  geom_smooth(aes(x=bvc_mean, y= bvn_mean/100), method=lm) +
  facet_wrap(~taxon, scales = 'free') +
  theme_bw(14) +
  xlab('Centroid biotic velocity (km/decade)') +
  ylab('Northern Biotic velocity (km/decade)')
ggsave('figures/BVN_vs_BVC_mean_scatter_credible_taxon.pdf', width=10, height=8)
ggsave('figures/BVN_vs_BVC_mean_scatter_credible_taxon.png', width=10, height=8)

ggplot(data=Vs) +
  geom_point(aes(x=bvc_mean/100, y= bvs_mean/100)) +
  geom_linerange(aes(x=bvc_mean/100, ymin=bvs_lo/100, ymax=bvs_hi/100)) +
  geom_smooth(aes(x=bvc_mean/100, y= bvs_mean/100), method=lm) +
  facet_wrap(~taxon, scales = 'free') +
  theme_bw(14) +
  xlab('Centroid biotic velocity (km/decade)') +
  ylab('Southern Biotic velocity (km/decade)')
ggsave('figures/BVS_vs_BVC_mean_scatter_credible_taxon.pdf', width=10, height=8)
ggsave('figures/BVS_vs_BVC_mean_scatter_credible_taxon.png', width=10, height=8)


ggplot(data=Vs) +
  geom_point(aes(x=bvs_mean/100, y= bvn_mean/100)) +
  geom_linerange(aes(x=bvs_mean/100, ymin=bvn_lo/100, ymax=bvn_hi/100)) +
  geom_smooth(aes(x=bvs_mean/100, y= bvn_mean/100), method=lm) +
  facet_wrap(~taxon, scales = 'free') +
  theme_bw(14) +
  xlab('Southern biotic velocity (km/decade)') +
  ylab('Northern Biotic velocity (km/decade)') 
ggsave('figures/BVN_vs_BVS_mean_scatter_credible_taxon.pdf', width=10, height=8)
ggsave('figures/BVN_vs_BVS_mean_scatter_credible_taxon.png', width=10, height=8)


Vs_melt_mean = melt(Vs[,c('period', 'taxon', 'bvc_mean', 'bvn_mean', 'bvs_mean', 'cv_mean')], 
     id.vars = c('period', 'taxon', 'cv_mean'))
Vs_melt_lo = melt(Vs[,c('period', 'taxon', 'bvc_lo', 'bvn_lo', 'bvs_lo', 'cv_mean')], 
                    id.vars = c('period', 'taxon', 'cv_mean'))
Vs_melt_hi = melt(Vs[,c('period', 'taxon', 'bvc_hi', 'bvn_hi', 'bvs_hi', 'cv_mean')], 
                  id.vars = c('period', 'taxon', 'cv_mean'))

Vs_melt = data.frame(Vs_melt_mean[,1:4], 
                     type = substr(Vs_melt_mean$variable, 1, 3),
                     bv_mean = Vs_melt_mean$value,
                     bv_lo = Vs_melt_lo$value,
                     bv_hi = Vs_melt_hi$value)


ggplot(data=Vs_melt) +
  geom_point(aes(x=cv_mean, y= bv_mean/100, colour=type)) +
  geom_linerange(aes(x=cv_mean, ymin=bv_lo/100, ymax=bv_hi/100, colour=type)) +
  geom_smooth(aes(x=cv_mean, y= bv_mean/100, colour=type, group=type), method=lm) +
  facet_wrap(~taxon, scales = 'free_y') +
  theme_bw(14) +
  xlab('Climate velocity (km/decade)') +
  ylab('Biotic velocity (km/decade)')
ggsave('figures/BVs_vs_CV_mean_scatter_credible_taxon.pdf', width=10, height=8)
ggsave('figures/BVs_vs_CV_mean_scatter_credible_taxon.png', width=10, height=8)


# ggplot(data=Vs_melt) +
  



ggplot(data=subset(Vs, period<14000)) +
  geom_point(aes(x=cv_mean, y= bvc_mean/100)) +
  geom_linerange(aes(x=cv_mean, ymin=bvc_lo/100, ymax=bvc_hi/100)) +
  geom_smooth(aes(x=cv_mean, y= bvc_mean/100), method=lm) +
  facet_wrap(~taxon, scales = 'free_y') +
  theme_bw(14) +
  xlab('Climate velocity (km/decade)') +
  ylab('Biotic velocity (km/decade)')
ggsave('figures/BVC_vs_CV_mean_scatter_credible_taxon_14000_mod.pdf', width=10, height=8)

ggplot(data=subset(Vs, period<12000)) +
  geom_point(aes(x=cv_mean, y= bvc_mean/100)) +
  geom_linerange(aes(x=cv_mean, ymin=bvc_lo/100, ymax=bvc_hi/100)) +
  geom_smooth(aes(x=cv_mean, y= bvc_mean/100), method=lm) +
  facet_wrap(~taxon, scales = 'free_y') +
  theme_bw(14) +
  xlab('Climate velocity (km/decade)') +
  ylab('Biotic velocity (km/decade)')
ggsave('figures/BVC_vs_CV_mean_scatter_credible_taxon_holocene.pdf', width=10, height=8)
ggsave('figures/BVC_vs_CV_mean_scatter_credible_taxon_holocene.png', width=10, height=8)

# ggplot(data=subset(Vs_melt, period<12000)) +
#   geom_point(aes(x=cv_mean, y= bv_mean/100, colour = type)) +
#   geom_linerange(aes(x=cv_mean, ymin=bv_lo/100, ymax=bv_hi/100, colour = type)) +
#   geom_smooth(aes(x=cv_mean, y= bv_mean/100, colour=type, group=type), method=lm) +
#   facet_wrap(~taxon, scales = 'free_y') +
#   theme_bw(14) +
#   xlab('Climate velocity (km/decade)') +
#   ylab('Biotic velocity (km/decade)')
# ggsave('figures/BVC_vs_CV_mean_scatter_credible_taxon_holocene.pdf', width=10, height=8)
# ggsave('figures/BVC_vs_CV_mean_scatter_credible_taxon_holocene.png', width=10, height=8)


ggplot() +
  geom_point(data=Vs, aes(x=cv_mean, y= bvc_mean/100)) +
  geom_smooth(data=Vs, aes(x=cv_mean, y= bvc_mean/100), method=lm) +
  theme_bw(14) +
  xlab('Climate velocity (km/decade)') +
  ylab('Biotic velocity (km/decade)')
ggsave('figures/BVC_vs_CV_mean_scatter.pdf', width=10, height=8)

# ggplot() +
#   geom_point(data=Vs, aes(x=cv_median, y= bv_median/100)) +
#   geom_smooth(data=Vs, aes(x=cv_median, y= bv_median/100), method=lm) #+
#   # facet_wrap(~taxon)

ggplot() +
  geom_point(data=Vs, aes(x=cv_median, y= bvc_median)) +
  facet_wrap(~period)

BV_CV_lm = Vs_melt %>%
  nest_by(taxon, type) %>%
  mutate(mod = list(lm(bv_mean ~ cv_mean, data = data))) %>%
  reframe(tidy(mod)) 
colnames(BV_CV_lm)[which(colnames(BV_CV_lm) == 'estimate')] = 'LGM'

BV_CV_lm_hol = subset(Vs_melt, period<12000) %>%
  nest_by(taxon, type) %>%
  mutate(mod = list(lm(bv_mean ~ cv_mean, data = data))) %>%
  reframe(tidy(mod))
colnames(BV_CV_lm_hol)[which(colnames(BV_CV_lm_hol) == 'estimate')] = 'Holocene'

BV_CV_slopes = data.frame(subset(BV_CV_lm, term == 'cv_mean'), 
                          Holocene = subset(BV_CV_lm_hol, term == 'cv_mean')$Holocene)
BV_CV_slopes$diff = BV_CV_slopes$LGM - BV_CV_slopes$Holocene

ggplot(data=BV_CV_slopes) +
  geom_histogram(aes(x = diff), breaks = seq(-60, 200, by=15)) + 
  facet_grid(type~.)
                          
# BV_CV_slopes = rbind(data.frame(subset(BV_CV_lm, term == 'cv_mean'), period = '21k'),
#                      data.frame(subset(BV_CV_lm_hol, term == 'cv_mean'), period = 'holocene'))

BV_CV_slopes$sig = 0
BV_CV_slopes$sig[which(BV_CV_slopes$p.value<0.05)] = 1




ggplot(data=BV_CV_slopes) +
  geom_point(aes(x = estimate, 
                 y = taxon, 
                 shape = factor(sig), 
                 colour=type), 
             size=2)

ggplot(data=slopes) +
  geom_histogram(aes(x = estimate, 
                 colour = type,
                 fill = type), 
               alpha=0.5)

ggplot(data=BV_CV_slopes) +
  geom_histogram(aes(x = LGM), bins = 15, breaks = seq(-60, 200, by = 15),
                 alpha=0.5, colour='blue', fill='lightblue', alpha=0.5) +
  geom_vline(xintercept = 0, colour = 'black', lty = 2) +
  facet_grid(type~.) +
  theme_bw() +
  theme(text = element_text(size=14))
ggsave('figures/BV_vs_CV_slopes_histogram.pdf', width=10, height=8)
ggsave('figures/BV_vs_CV_slopes_histogram.png', width=10, height=8)


# ggplot(data=slopes) +
#   geom_density(aes(x = estimate), 
#                  alpha=0.5) +
#   facet_grid(type~.)


######################################################################################################################
## NS quant 95 BV
######################################################################################################################

