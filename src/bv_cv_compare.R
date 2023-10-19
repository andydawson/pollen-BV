library(dplyr)

pollen_BVs = readRDS('output/pollen_BVs_all.RDS')

pollen_BVs_summary = pollen_BVs %>% 
  group_by(timeFrom, timeTo, taxon) %>% 
  summarize(bv_mean = mean(centroidVelocity),
            bv_median = median(centroidVelocity),
            bv_sd = sd(centroidVelocity),
            bv_lo = quantile(centroidVelocity, c(0.05)),
            bv_hi = quantile(centroidVelocity, c(0.95)), .groups="keep")

pollen_BVs_summary$period = abs(pollen_BVs_summary$timeFrom + 990/2)

dim(pollen_BVs_summary)

pollen_BVs_summary$period = floor(pollen_BVs_summary$period/1000)*1000


ggplot(data=pollen_BVs_summary) +
  geom_point(aes(x=timeFrom, y=bv_median)) +
  geom_linerange(aes(x=timeFrom, ymin=bv_lo, ymax=bv_hi)) +
  facet_wrap(~taxon) +
  ylab('Biotic velocity (m/year)') +
  xlab('Time (YBP)') +
  theme_bw(14) +
  theme(axis.text.x = element_text(angle=90,hjust=1)) 
          

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
  ylab('Biotic velocity (m/year)') +
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
ggsave('figures/CV_vs_time_boxplot.png', width=10, height=6)


Vs = merge(pollen_BVs_summary, CVs_summary)

head(Vs)

ggplot() +
  geom_point(data=Vs, aes(x=cv_median, y= bv_median/100)) +
  geom_smooth(data=Vs, aes(x=cv_median, y= bv_median/100), method=lm) +
  facet_wrap(~taxon) +
  theme_bw(14) +
  xlab('Climate velocity (km/decade)') +
  ylab('Biotic velocity (km/decade)')
ggsave('figures/BV_vs_CV_vs_scatter.pdf', width=10, height=8)

# ggplot() +
#   geom_point(data=Vs, aes(x=cv_median, y= bv_median/100)) +
#   geom_smooth(data=Vs, aes(x=cv_median, y= bv_median/100), method=lm) #+
#   # facet_wrap(~taxon)

ggplot() +
  geom_point(data=Vs, aes(x=cv_median, y= bv_median)) +
  facet_wrap(~period)
