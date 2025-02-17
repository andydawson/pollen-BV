library(dplyr)
library(ggplot2)
library(reshape2)
library(broom)
library(tidyr)
library(corrplot)
library(ggcorrplot)
library(FactoMineR)
library(factoextra)
library(Hmisc)

name_pair = data.frame(suffix = c('temp', 'temp', 'prcp', 'et', 'et'), 
                       varname = c('tmax', 'tmin', 'prcp', 'pet', 'aet'),
                       varlabel = c('maximum temperature', 
                                    'minimum temperature', 
                                    'precipitation', 
                                    'potential evapotranspiration',
                                    'actual evapotr'))

months = c('jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec')

N_climvars = nrow(name_pair)

Vs_recalc = FALSE

######################################################################################################################
## read and summarize biotic velocities
######################################################################################################################

# pollen_BVs = readRDS('output/pollen_BVs_all_v5.0.RDS')
pollen_BVs = readRDS('output/pollen_BVs_all_masked_n200_v5.0.RDS')

pollen_BVs_summary = pollen_BVs %>% 
  group_by(timeFrom, timeTo, timeMid, taxon) %>% 
  dplyr::summarize(bvc_mean = mean(centroidVelocity),
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
taxa = unique(pollen_BVs_summary$taxon)


BV_colnames = as.vector(sapply(c('bvc', 'bvn', 'bvs'), function(x) paste0(x, c('_mean', '_median', '_lo', '_hi'))))
pollen_BVs_summary[, BV_colnames] = pollen_BVs_summary[, BV_colnames]/100

######################################################################################################################
## BV temporal patterns
######################################################################################################################
# pollen_BVs_scale = pollen_BVs %>% 
#   group_by(taxon) %>% 
#   dplyr::summarize(bvc_mean = mean(centroidVelocity),
#                    bvc_median = median(centroidVelocity),
#                    bvc_sd = sd(centroidVelocity))
# 
# pollen_BVs_scale = pollen_BVs_summary %>% 
#   group_by(taxon) %>% 
#   dplyr::summarize(bvc_all_mean = mean(bvc_mean),
#                    # bvc_median = median(bv),
#                    bvc_all_sd = sd(bvc_mean, na.rm=TRUE))#,
#                    # bvc_scaled = (bvc_mean - bvc_all_mean)/bvc_all_sd)
# 
# pollen_BVs_summary= merge(pollen_BVs_summary, pollen_BVs_scale, by = 'taxon', all.x = TRUE)
# pollen_BVs_summary$bvc_scaled = (pollen_BVs_summary$bvc_mean - pollen_BVs_summary$bvc_all_mean) / pollen_BVs_summary$bvc_all_sd
# 
# 
# ggplot(data=pollen_BVs_summary) + 
#   geom_point(aes(x=timeMid, y=bvc_scaled)) +
#   facet_wrap(~taxon) +
#   geom_hline(yintercept = 0, linetype=2)+
#   geom_hline(yintercept = 2, linetype=3) +
#   geom_hline(yintercept = -2, linetype=3)




######################################################################################################################
## BV temporal patterns
######################################################################################################################
BVC_taxa_wide = pollen_BVs_summary[,colnames(pollen_BVs_summary) %in% 
                                     c('timeFrom', 'taxon', 'bvc_mean')] %>% 
  pivot_wider(id_cols = 'timeFrom', values_from = 'bvc_mean', names_from = 'taxon') 

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
BVN_taxa_wide = pollen_BVs_summary[,colnames(pollen_BVs_summary) %in% 
                                     c('timeFrom', 'taxon', 'bvn_mean')] %>% 
  pivot_wider(id_cols = 'timeFrom', values_from = 'bvn_mean', names_from = 'taxon') 

BVN_wide = BVN_taxa_wide[,2:ncol(BVN_taxa_wide)]

BVN_normalized = scale(BVN_wide)

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
BVS_taxa_wide = pollen_BVs_summary[,colnames(pollen_BVs_summary) %in% 
                                     c('timeFrom', 'taxon', 'bvs_mean')] %>% 
  pivot_wider(id_cols = 'timeFrom', values_from = 'bvs_mean', names_from = 'taxon') 

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

saveRDS(pca_time_melt, 'output/PCA_BVs.RDS')
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
## CVs
######################################################################################################################



######################################################################################################################
## climate velocity masked by taxon
######################################################################################################################

# name_pair

# add all CVs as columns to this data frame
# matching by period and taxon for n = 1
# then matching by period, taxon, and month for other vars
# Vs_all = pollen_BVs_summary

if (Vs_recalc){
  
  for (n in 1:N_climvars){
    
    varname = name_pair$varname[n]
    CVs_var = readRDS(paste0('output/climate_velocity/', varname, '_CV_taxa.RDS'))
    CVs_var$voccMag = CVs_var$voccMag / 100
    CVs_var$voccMag_pos = abs(CVs_var$voccMag)
    
    CVs_summary = CVs_var %>%
      group_by(period, taxon, month) %>%
      dplyr::summarize(cv_mean = mean(voccMag_pos, na.rm=TRUE),
                       cv_median = median(voccMag_pos, na.rm=TRUE),
                       cv_sd = sd(voccMag_pos, na.rm=TRUE),
                       cv_lo = quantile(voccMag_pos, c(0.05), na.rm=TRUE),
                       cv_hi = quantile(voccMag_pos, c(0.95), na.rm=TRUE), .groups="keep")
    CVs_summary$clim_var = rep(varname, nrow(CVs_summary))
    
    Vs = merge(pollen_BVs_summary, 
               CVs_summary, by = c('period', 'taxon'))
    Vs$era_name = NA
    Vs$era_name[which(Vs$period>18000)] = 'near LGM'
    Vs$era_name[which((Vs$period<18000)&(Vs$period>12000))] = 'post LGM'
    Vs$era_name[which((Vs$period<12000)&(Vs$period>2000))] = 'pre-industrial Holocene'
    Vs$era_name[which(Vs$period<2000)] = 'Anthropocene'
    
    Vs$split_holocene = NA
    Vs$split_holocene[which(Vs$period<12000)] = 'Holocene'
    Vs$split_holocene[which(Vs$period>=12000)] = 'pre-Holocene'
    
    head(Vs)
    
    
    pdf(paste0('figures/BVC_vs_CV_mean_scatter_credible_taxon_', varname, '.pdf'))
    
    for (i in 1:12){
      
      month= months[i]
      print(paste0(varname, '; ', month))
      
      p = ggplot(data=Vs[which(Vs$month == month),]) +
        geom_point(aes(x=cv_mean, y= bvc_mean)) +
        geom_linerange(aes(x=cv_mean, ymin=bvc_lo, ymax=bvc_hi)) +
        geom_smooth(aes(x=cv_mean, y= bvc_mean), method=lm) +
        # facet_grid(~taxon) +
        facet_wrap(~taxon, scales = 'free') +
        theme_bw(14) +
        xlab('Climate velocity (km/decade)') +
        ylab('Biotic velocity (km/decade)') +
        labs(title = paste0(varname, '; ', month))
      # ggsave(pasteo('figures/BVC_vs_CV_mean_scatter_credible_taxon_', varname, '.pdf'), width=10, height=8)
      # ggsave(pasteo('figures/BVC_vs_CV_mean_scatter_credible_taxon_', varname, '.pdf'), width=10, height=8)
      print(p) 
    }
    dev.off()
    
    # if (n==1){
    #   Vs_all = merge(Vs_all, 
    #                  CVs_summary, by = c('period', 'taxon'))
    # } else {
    #   Vs_all = merge(Vs_all, CVs_summary, by = c('period', 'taxon', 'month'))
    # }
    if (n==1){
      Vs_all = Vs
    } else {
      Vs_all = rbind(Vs_all, Vs)
    }
  }
  
  saveRDS(Vs_all, 'output/BV_CV_taxa_vars_all.RDS')
} else {
  Vs_all = readRDS('output/BV_CV_taxa_vars_all.RDS')
}

Vs = Vs_all

Vs_melt_mean = melt(Vs[,c('period', 'era_name', 'split_holocene', 'clim_var', 'month', 'taxon', 'bvc_mean', 'bvn_mean', 'bvs_mean', 'cv_mean')], 
                    id.vars = c('period', 'era_name', 'split_holocene', 'clim_var', 'month', 'taxon', 'cv_mean'))
Vs_melt_median = melt(Vs[,c('period', 'era_name', 'split_holocene', 'clim_var','month',  'taxon', 'bvc_median', 'bvn_median', 'bvs_median', 'cv_mean')], 
                      id.vars = c('period', 'era_name', 'split_holocene', 'clim_var', 'month', 'taxon', 'cv_mean'))
Vs_melt_lo = melt(Vs[,c('period', 'era_name', 'split_holocene', 'clim_var','month', 'taxon', 'bvc_lo', 'bvn_lo', 'bvs_lo', 'cv_mean')], 
                  id.vars = c('period', 'era_name', 'split_holocene', 'clim_var', 'month', 'taxon', 'cv_mean'))
Vs_melt_hi = melt(Vs[,c('period', 'era_name', 'split_holocene', 'clim_var','month',  'taxon', 'bvc_hi', 'bvn_hi', 'bvs_hi', 'cv_mean')], 
                  id.vars = c('period', 'era_name', 'split_holocene', 'clim_var', 'month', 'taxon', 'cv_mean'))

Vs_melt = data.frame(Vs_melt_mean[,1:7], 
                     type = substr(Vs_melt_mean$variable, 1, 3),
                     bv_mean = Vs_melt_mean$value,
                     bv_median = Vs_melt_median$value,
                     bv_lo = Vs_melt_lo$value,
                     bv_hi = Vs_melt_hi$value)

# Vs_melt$split_name = NA
# Vs_melt$split_name[which(Vs_melt$period <= 12000)] = 'Holocene'
# Vs_melt$split_name[which(Vs_melt$period > 12000)] = 'pre_Holocene'

Vs_melt$month = factor(Vs_melt$month, levels = months)

N_taxa = length(taxa)


# BVs between two consecutive time bins 990 years apart
Vs_taxon_cor = Vs %>%
  # pivot_wider(id_cols = c(period, taxon), 
  #             names_from = month, 
  #             values_from = cv_median) %>%
  group_by(taxon) %>%
  dplyr::summarize(BVN_BVC =  cor(bvn_median, bvc_median),
                   BVS_BVC =  cor(bvs_median, bvc_median),
                   BVN_BVS =  cor(bvn_median, bvs_median)) %>%
  pivot_longer(cols = c('BVN_BVC', 'BVS_BVC', 'BVN_BVS'))

colnames(Vs_taxon_cor) = c('taxon', 'type', 'All21k')



Vs_taxon_cor_split = Vs %>%
  # pivot_wider(id_cols = c(period, taxon), 
  #             names_from = month, 
  #             values_from = cv_median) %>%
  group_by(taxon, split_holocene) %>%
  dplyr::summarize(BVN_BVC =  cor(bvn_median, bvc_median),
                   BVS_BVC =  cor(bvs_median, bvc_median),
                   BVN_BVS =  cor(bvn_median, bvs_median)) %>%
  pivot_longer(cols = c('BVN_BVC', 'BVS_BVC', 'BVN_BVS'))

Vs_taxon_cor_split_long = Vs_taxon_cor_split %>%
  pivot_wider(id_cols = c(taxon, name),
              values_from = value,
              names_from = split_holocene)
colnames(Vs_taxon_cor_split_long) = c('taxon', 'type', 'Hol', 'preHol')

Vs_taxon_cor_all = merge(Vs_taxon_cor, Vs_taxon_cor_split_long, by = c('taxon', 'type'))
# Vs_taxon_cor_all$type = factor(Vs_taxon_cor_all$type, 
#                                levels = c('BVN_BC')
# Vs_taxon_cor_all$type = factor(Vs_taxon_cor_all$type, 
#                                levels = c('BVN_BVC_cor',
#                                           'BNS_BVC_cor',
#                                           'BVN_BVS_cor'))#,
# labels = c('BVN-BVC',
#            'BVS-BVC', 
#            'BVN-BVS'))
ggplot(data=Vs_taxon_cor_all) +
  geom_point(aes(x=preHol, y=Hol, colour=type)) +
  geom_smooth(method=lm, formula=y~x, aes(x=preHol, y=Hol, colour=type))

Vs_taxon_cor_all$diff_preHol_Hol = Vs_taxon_cor_all$preHol - Vs_taxon_cor_all$Hol

# not interesting
ggplot(data=Vs_taxon_cor_all) +
  geom_density(aes(x=diff_preHol_Hol), fill='lightgrey') + #, colour=variable, fill=variable)) + 
  # geom_vline(data=Vs_taxon_cor_all_mean, aes(xintercept=cor_mean), colour='dodgerblue') + 
  facet_grid(type~.) +
  theme_bw(16) +
  xlab('Correlation')
# ggsave('figures/BV_cor_by_split.pdf', width=12, height=8)

Vs_taxon_cor_all_melt = melt(Vs_taxon_cor_all, 
                             id.vars = c('taxon', 'type'))

Vs_taxon_cor_all_mean = Vs_taxon_cor_all_melt %>%
  group_by(type, variable) %>%
  dplyr::summarize(cor_mean = mean(value))

ggplot(data=Vs_taxon_cor_all_melt) +
  geom_density(aes(x=value), fill='lightgrey') + #, colour=variable, fill=variable)) + 
  geom_vline(data=Vs_taxon_cor_all_mean, aes(xintercept=cor_mean), colour='dodgerblue') + 
  facet_grid(variable~type) +
  theme_bw(16) +
  xlab('Correlation')
ggsave('figures/BV_cor_by_split.pdf', width=12, height=8)


ggplot(data=subset(Vs_taxon_cor_all_melt, variable != 'All21k')) +
  geom_point(aes(x=value, y=taxon, colour=variable)) + 
  theme_bw(16) +
  xlab('Correlation') +
  facet_grid(type~.)
# ggsave('figures/BV_cor_by_split.pdf', width=12, height=8)

ggplot(data=subset(Vs_taxon_cor_all_melt, variable != 'All21k')) +
  geom_col(aes(x=value, y=taxon, colour=variable, fill=variable), position="dodge") + 
  theme_bw(16)+
  xlab('Correlation') +
  facet_grid(type~.)

ggplot(data=subset(Vs_taxon_cor_all_melt, variable != 'All21k')) +
  geom_col(aes(y=value, x=taxon, colour=variable, fill=variable), position="dodge") + 
  theme_bw(16)+
  xlab('Correlation') +
  facet_grid(type~.)


foo = subset(Vs_taxon_cor_all_melt, variable != 'diff_preHol_Hol')
ggplot(data=subset(foo, variable != 'All21k')) +
  geom_col(aes(y=value, x=taxon, colour=variable, fill=variable), position="dodge") + 
  theme_bw(16)+
  xlab('Correlation') +
  facet_grid(type~.)
# BV_CV_lm = Vs_melt[, c('period', 'era_name', 'split_holocene', 'taxon', 'name', 
#                        'type', 'bv_median')] %>%
#   pivot_wider(id_cols = c(period, era_name, split_holocene, taxon, name),
#               values_from = value,
#               names_from = split_holocene) %>%

BVN_BVC_lm = Vs %>%
  nest_by(taxon) %>%
  mutate(mod = list(lm(bvn_median ~ bvc_median, data = data))) %>%
  reframe(tidy(mod)) %>%
  subset(., term=='bvc_median')
BVN_BVC_text <- data.frame(
  taxon=BVN_BVC_lm$taxon,
  xpos = rep(Inf, length(taxa)),
  ypos =  rep(Inf, length(taxa)),
  annotateText = round(BVN_BVC_lm$estimate,2),
  hjustvar = rep(1, length(taxa)),
  vjustvar = rep(1, length(taxa))) 

BVS_BVC_lm = Vs %>%
  nest_by(taxon) %>%
  mutate(mod = list(lm(bvs_median ~ bvc_median, data = data))) %>%
  reframe(tidy(mod))  %>%
  subset(., term=='bvc_median')
BVS_BVC_text <- data.frame(
  taxon=BVS_BVC_lm$taxon,
  xpos = rep(Inf, length(taxa)),
  ypos =  rep(Inf, length(taxa)),
  annotateText = round(BVS_BVC_lm$estimate,2),
  hjustvar = rep(1, length(taxa)),
  vjustvar = rep(1, length(taxa))) 

BVN_BVS_lm = Vs %>%
  nest_by(taxon) %>%
  mutate(mod = list(lm(bvn_median ~ bvs_median, data = data))) %>%
  reframe(tidy(mod)) %>%
  subset(., term=='bvs_median')
BVN_BVS_text <- data.frame(
  taxon=BVN_BVS_lm$taxon,
  xpos = rep(Inf, length(taxa)),
  ypos =  rep(Inf, length(taxa)),
  annotateText = round(BVN_BVS_lm$estimate,2),
  hjustvar = rep(1, length(taxa)),
  vjustvar = rep(1, length(taxa))) 

# geom_text(data=BVN_BVC_text, x=Inf,y=-Inf,hjust=1.5,vjust=-1.5, aes(label=annotateText)) +


# range_act <- range(range(results$act), range(results$pred))
# dummy <- data.frame(pred = range_act, value = range_act,
#                     variable = "act", stringsAsFactors=FALSE)

pdf('figures/BV_vs_BV_mean_scatter_credible_taxon.pdf', width=14, height=12)
p = ggplot() +
  geom_point(data=Vs, aes(x=bvc_mean, y= bvn_mean), alpha=0.5) +
  geom_linerange(data=Vs, aes(x=bvc_mean, ymin=bvn_lo, ymax=bvn_hi), alpha=0.5) +
  geom_linerange(data=Vs, aes(y=bvn_mean, xmin=bvc_lo, xmax=bvc_hi), alpha=0.5) +
  geom_smooth(data=Vs, aes(x=bvc_mean, y= bvn_mean), method=lm, fullrange=TRUE) +
  theme_bw(14) +
  theme(aspect.ratio=1) +
  xlab('Centroid biotic velocity (km/decade)') +
  ylab('Northern Biotic velocity (km/decade)') +
  geom_text(data=BVN_BVC_text, x=Inf,y=-Inf,hjust=1.5,vjust=-1.5, aes(label=annotateText)) +
  facet_wrap(~taxon, scales = 'free') +
  # facet_wrap(~taxon) +
  geom_abline(intercept=0, slope=1, linetype=2) +
  geom_hline(yintercept=0, linetype=2) 
# ggsave('figures/BVN_vs_BVC_mean_scatter_credible_taxon.pdf', width=10, height=8)
# ggsave('figures/BVN_vs_BVC_mean_scatter_credible_taxon.png', width=10, height=8)
print(p)
# dev.off()
# 
# pdf(paste0('figures/BVS_vs_BVC_mean_scatter_credible_taxon_', varname, '.pdf'))
# 
# for (i in 1:12){
#   
#   month= months[i]
#   
#   print(paste0(varname, '; ', month))

p = ggplot(data=Vs) +
  geom_point(aes(x=bvc_mean, y= bvs_mean), alpha=0.5) +
  geom_linerange(aes(x=bvc_mean, ymin=bvs_lo, ymax=bvs_hi), alpha=0.5) +
  geom_linerange(aes(y=bvs_mean, xmin=bvc_lo, xmax=bvc_hi), alpha=0.5) +
  geom_smooth(aes(x=bvc_mean, y= bvs_mean), method=lm, fullrange=TRUE) +
  facet_wrap(~taxon, scales = 'free') +
  theme_bw(14) +
  theme(aspect.ratio=1) +
  xlab('Centroid biotic velocity (km/decade)') +
  ylab('Southern Biotic velocity (km/decade)') +
  geom_text(data=BVS_BVC_text, x=Inf,y=-Inf,hjust=1.5,vjust=-1.5, aes(label=annotateText)) +
  facet_wrap(~taxon, scales = 'free') +
  geom_abline(intercept=0, slope=1, linetype=2) +
  geom_hline(yintercept=0, linetype=2)
print(p)
# ggsave('figures/BVS_vs_BVC_mean_scatter_credible_taxon.pdf', width=10, height=8)
# ggsave('figures/BVS_vs_BVC_mean_scatter_credible_taxon.png', width=10, height=8)
# }
# dev.off()

# pdf(paste0('figures/BVN_vs_BVS_mean_scatter_credible_taxon_', varname, '.pdf'))
# for (i in 1:12){
#   
#   month= months[i]
#   
#   print(paste0(varname, '; ', month))

p = ggplot(data=Vs) +
  geom_point(aes(x=bvs_mean, y= bvn_mean), alpha=0.5) +
  geom_linerange(aes(x=bvs_mean, ymin=bvn_lo, ymax=bvn_hi), alpha=0.5) +
  geom_linerange(aes(y=bvn_mean, xmin=bvs_lo, xmax=bvs_hi), alpha=0.5) +
  geom_smooth(aes(x=bvs_mean, y= bvn_mean), method=lm, fullrange=TRUE) +
  facet_wrap(~taxon, scales = 'free') +
  theme_bw(14) +
  theme(aspect.ratio=1) +
  xlab('Southern biotic velocity (km/decade)') +
  ylab('Northern Biotic velocity (km/decade)') +
  geom_text(data=BVN_BVS_text, x=Inf,y=-Inf,hjust=1.5,vjust=-1.5, aes(label=annotateText)) +
  facet_wrap(~taxon, scales = 'free') +
  geom_abline(intercept=0, slope=1, linetype=2) +
  geom_hline(yintercept=0, linetype=2)
print(p)
# ggsave('figures/BVN_vs_BVS_mean_scatter_credible_taxon.pdf', width=10, height=8)
# ggsave('figures/BVN_vs_BVS_mean_scatter_credible_taxon.png', width=10, height=8)
# } 
dev.off()


CVs_sub = Vs[, c('period', 'taxon', 'clim_var', 'month', 'cv_median')] 


pdf(paste0('figures/CV_monthly_cor_by_taxon.pdf'))
for (i in 1:N_climvars){
  
  varname = name_pair$varname[i]
  
  print(paste0('Variable: ', varname))
  
  for (k in 1:N_taxa){
    taxon = taxa[k]
    
    CVs_var_taxon = CVs_sub[which((CVs_sub$clim_var == varname) & 
                                    (CVs_sub$taxon == taxon)),]
    
    CVs_var_taxon_cor = CVs_var_taxon %>%
      pivot_wider(id_cols = c(period, taxon), 
                  names_from = month, 
                  values_from = cv_median) %>%
      group_by(taxon) %>%
      dplyr::summarize(clim_cor =  cor(.[, c(months)])) 
    # dplyr::reframe(clim_cor =  cor(.[, c(months)])) 
    
    
    CVs_var_taxon_cor_mat = as.matrix(CVs_var_taxon_cor[, 2])
    colnames(CVs_var_taxon_cor_mat) = months
    row.names(CVs_var_taxon_cor_mat) = months
    
    # p = ggcorrplot(CVs_var_taxon_cor_mat, 
    #                method = "circle", 
    #                type = "lower", 
    #                pch=8, 
    #                show.diag=TRUE) +
    #   labs(title = paste0(varname, '; ', taxon))
    
    # p = ggcorrplot(CVs_var_taxon_cor_mat, method = "circle", 
    p = ggcorrplot(CVs_var_taxon_cor_mat, 
                   method = "square", 
                   type = "lower", 
                   show.diag=TRUE) +
      labs(title = paste0(varname, '; ', taxon))
    print(p)
  }
}
dev.off()


pdf(paste0('figures/CV_monthly_cor_across_taxa.pdf'))
for (i in 1:N_climvars){
  
  varname = name_pair$varname[i]
  
  print(paste0('Variable: ', varname))
  
  for (j in 1:12){
    month = months[j]
    
    CVs_var_month = CVs_sub[which((CVs_sub$clim_var == varname) & 
                                    (CVs_sub$month == month)),]
    
    CVs_var_month_cor = CVs_var_month %>%
      pivot_wider(id_cols = c(period, month), 
                  names_from = taxon, 
                  values_from = cv_median) %>%
      # group_by(month) %>%
      dplyr::summarize(clim_cor =  cor(.[, c(taxa)])) 
    # dplyr::reframe(clim_cor =  cor(.[, c(months)])) 
    
    
    CVs_var_month_cor_mat = as.matrix(CVs_var_month_cor)
    colnames(CVs_var_month_cor_mat) = taxa
    row.names(CVs_var_month_cor_mat) = taxa
    
    p = ggcorrplot(CVs_var_month_cor_mat, 
                   method = "circle", 
                   type = "lower", 
                   pch=4, 
                   show.diag=TRUE) +
      labs(title = paste0(varname, '; ', month))
    
    print(p)
  }
}
dev.off()


# foo = subset(CVs_sub, clim_var == 'aet') %>%
#   # pivot_wider(id_cols = c(period, taxon), 
#   #             names_from = month, 
#   #             values_from = cv_median) %>%
#   group_by(taxon) %>%
#   dplyr::summarize(clim_cor =  cor(cv_median))

# subset(foo, clim_var == 'aet') %>% 
#   group_by(taxon, month) %>%
#   summarize(cor())

## Linear model: BV vs CV; All 21 k
BV_CV_lm = Vs_melt %>%
  nest_by(taxon, type, clim_var, month) %>%
  mutate(mod = list(lm(bv_mean ~ cv_mean, data = data))) %>%
  reframe(tidy(mod)) 
colnames(BV_CV_lm)[which(colnames(BV_CV_lm) == 'estimate')] = 'All_21k'
colnames(BV_CV_lm)[which(colnames(BV_CV_lm) == 'std.error')] = 'All_21k.std.error'
colnames(BV_CV_lm)[which(colnames(BV_CV_lm) == 'statistic')] = 'All_21k.statistic'
colnames(BV_CV_lm)[which(colnames(BV_CV_lm) == 'p.value')] = 'All_21k.p.value'

## Linear model: BV vs CV; pre-Holocene
BV_CV_lm_post = subset(Vs_melt, split_holocene == 'pre-Holocene') %>%
  nest_by(taxon, type, clim_var, month) %>%
  mutate(mod = list(lm(bv_mean ~ cv_mean, data = data))) %>%
  reframe(tidy(mod)) 
colnames(BV_CV_lm_post)[which(colnames(BV_CV_lm_post) == 'estimate')] = 'preHol'
colnames(BV_CV_lm_post)[which(colnames(BV_CV_lm_post) == 'std.error')] = 'preHol.std.error'
colnames(BV_CV_lm_post)[which(colnames(BV_CV_lm_post) == 'statistic')] = 'preHol.statistic'
colnames(BV_CV_lm_post)[which(colnames(BV_CV_lm_post) == 'p.value')] = 'preHol.p.value'

## Linear model: BV vs CV; Holocene
BV_CV_lm_hol = subset(Vs_melt, split_holocene == 'Holocene') %>%
  nest_by(taxon, type, clim_var, month) %>%
  mutate(mod = list(lm(bv_mean ~ cv_mean, data = data))) %>%
  reframe(tidy(mod))
colnames(BV_CV_lm_hol)[which(colnames(BV_CV_lm_hol) == 'estimate')] = 'Hol'
colnames(BV_CV_lm_hol)[which(colnames(BV_CV_lm_hol) == 'std.error')] = 'Hol.std.error'
colnames(BV_CV_lm_hol)[which(colnames(BV_CV_lm_hol) == 'statistic')] = 'Hol.statistic'
colnames(BV_CV_lm_hol)[which(colnames(BV_CV_lm_hol) == 'p.value')] = 'Hol.p.value'

## Linear model: BV vs CV; Merge period estimates
BV_CV_slopes = merge(BV_CV_lm, BV_CV_lm_hol, 
                     by = c('taxon', 'type', 'clim_var', 'month', 'term'))
BV_CV_slopes = merge(BV_CV_slopes, BV_CV_lm_post, 
                     by = c('taxon', 'type', 'clim_var', 'month', 'term'))

# only keep slopes (not intercepts)
BV_CV_slopes = subset(BV_CV_slopes, term == 'cv_mean')

# difference in slopes
BV_CV_slopes$diff_all21_Hol = BV_CV_slopes$All_21k - BV_CV_slopes$Hol
BV_CV_slopes$diff_preHol_Hol = BV_CV_slopes$preHol - BV_CV_slopes$Hol

BV_CV_slopes$month = factor(BV_CV_slopes$month, levels = months)

BV_CV_slopes$same_sign_all21_Hol = sign(BV_CV_slopes$All_21k) * sign(BV_CV_slopes$Hol)
BV_CV_slopes$same_sign_preHol_Hol = sign(BV_CV_slopes$preHol) * sign(BV_CV_slopes$Hol)

BV_CV_slopes$sign_case_all21_Hol = NA
BV_CV_slopes$sign_case_all21_Hol[which((sign(BV_CV_slopes$All_21k)>0) & (sign(BV_CV_slopes$Hol)>0))] = 'p_p'
BV_CV_slopes$sign_case_all21_Hol[which((sign(BV_CV_slopes$All_21k)>0) & (sign(BV_CV_slopes$Hol)<0))] = 'p_n'
BV_CV_slopes$sign_case_all21_Hol[which((sign(BV_CV_slopes$All_21k)<0) & (sign(BV_CV_slopes$Hol)<0))] = 'n_n'
BV_CV_slopes$sign_case_all21_Hol[which((sign(BV_CV_slopes$All_21k)<0) & (sign(BV_CV_slopes$Hol)>0))] = 'n_p'

BV_CV_slopes$sign_case_preHol_Hol = NA
BV_CV_slopes$sign_case_preHol_Hol[which((sign(BV_CV_slopes$preHol)>0) & (sign(BV_CV_slopes$Hol)>0))] = 'p_p'
BV_CV_slopes$sign_case_preHol_Hol[which((sign(BV_CV_slopes$preHol)>0) & (sign(BV_CV_slopes$Hol)<0))] = 'p_n'
BV_CV_slopes$sign_case_preHol_Hol[which((sign(BV_CV_slopes$preHol)<0) & (sign(BV_CV_slopes$Hol)<0))] = 'n_n'
BV_CV_slopes$sign_case_preHol_Hol[which((sign(BV_CV_slopes$preHol)<0) & (sign(BV_CV_slopes$Hol)>0))] = 'n_p'


BV_CV_slopes$sign_case_preHol_Hol = factor(BV_CV_slopes$sign_case_preHol_Hol,
                                           levels = c('p_p', 'p_n', 'n_p', 'n_n'),
                                           labels = c('+ / +', '+ / -', '- / +', '- / -'))


ggplot(data=BV_CV_slopes) +
  geom_histogram(aes(x = All_21k)) +#, breaks = seq(-60, 200, by=15)) + 
  facet_grid(type~.)

ggplot(data=BV_CV_slopes) +
  geom_point(aes(x = All_21k, y=Hol)) +#, breaks = seq(-60, 200, by=15)) + 
  facet_grid(type~.)

ggplot(data=BV_CV_slopes) +
  geom_point(aes(x = preHol, y=Hol)) +#, breaks = seq(-60, 200, by=15)) + 
  facet_grid(type~.)


# BVN_BVS_lm = Vs %>%
#   nest_by(taxon) %>%
#   mutate(mod = list(lm(bvn_median ~ bvs_median, data = data))) %>%
#   reframe(tidy(mod)) %>%
#   subset(., term=='bvs_median')

pdf(paste0('figures/BVC_vs_CV_mean_scatter_credible_taxon_split.pdf'), width=12, height=10)
for (i in 1:N_climvars){
  
  varname = name_pair$varname[i]
  
  print(paste0('Variable: ', varname))
  
  Vs_climvar = subset(Vs, clim_var == varname)
  
  for (k in 1:12){
    
    month= months[k]
    
    print(paste0('BVC vs ', varname, '; ', month, '; split'))
    
    BV_CV_slopes_climvar = BV_CV_slopes[which((BV_CV_slopes$clim_var == varname) & 
                                                (BV_CV_slopes$month==month) & 
                                                (BV_CV_slopes$type=='bvc')),]
    
    BV_CV_slopes_text = data.frame(
      taxon=BV_CV_slopes_climvar$taxon,
      xpos = rep(Inf, length(taxa)),
      ypos =  rep(Inf, length(taxa)),
      annotateText = paste0('pH: ', round(BV_CV_slopes_climvar$preHol,2), 
                            '; H: ', round(BV_CV_slopes_climvar$Hol,2)),
      hjustvar = rep(1, length(taxa)),
      vjustvar = rep(1, length(taxa))) 
    
    p = ggplot(data=Vs_climvar[which(Vs_climvar$month == month),]) +
      geom_point(aes(x=cv_mean, y= bvc_mean, colour=split_holocene)) +
      geom_linerange(aes(x=cv_mean, ymin=bvc_lo, ymax=bvc_hi, colour=split_holocene)) +
      geom_smooth(aes(x=cv_mean, y= bvc_mean, colour=split_holocene), method=lm) +
      # facet_grid(~taxon) +
      facet_wrap(~taxon, scales = 'free') +
      theme_bw(18) +
      xlab('Climate velocity (km/decade)') +
      ylab('Biotic velocity (km/decade)') +
      labs(title = paste0(varname, '; ', month)) +
      geom_text(data=BV_CV_slopes_text, x=Inf,y=Inf,hjust=1.2,vjust=1.5, 
                aes(label=annotateText), size=3)
    # ggsave(pasteo('figures/BVC_vs_CV_mean_scatter_credible_taxon_', varname, '.pdf'), width=10, height=8)
    # ggsave(pasteo('figures/BVC_vs_CV_mean_scatter_credible_taxon_', varname, '.pdf'), width=10, height=8)
    print(p) 
  }
  
  for (k in 1:12){
    
    month= months[k]
    
    print(paste0('BVN vs ', varname, '; ', month, '; split'))
    
    BV_CV_slopes_climvar = BV_CV_slopes[which((BV_CV_slopes$clim_var == varname) & 
                                                (BV_CV_slopes$month==month) & 
                                                (BV_CV_slopes$type=='bvn')),]
    
    BV_CV_slopes_text = data.frame(
      taxon=BV_CV_slopes_climvar$taxon,
      xpos = rep(Inf, length(taxa)),
      ypos =  rep(Inf, length(taxa)),
      annotateText = paste0('pH: ', round(BV_CV_slopes_climvar$preHol,2), 
                            '; H: ', round(BV_CV_slopes_climvar$Hol,2)),
      hjustvar = rep(1, length(taxa)),
      vjustvar = rep(1, length(taxa))) 
    
    p = ggplot(data=Vs_climvar[which(Vs_climvar$month == month),]) +
      geom_point(aes(x=cv_mean, y= bvn_mean, colour=split_holocene)) +
      geom_linerange(aes(x=cv_mean, ymin=bvn_lo, ymax=bvn_hi, colour=split_holocene)) +
      geom_smooth(aes(x=cv_mean, y= bvn_mean, colour=split_holocene), method=lm) +
      # facet_grid(~taxon) +
      facet_wrap(~taxon, scales = 'free') +
      theme_bw(18) +
      xlab('Climate velocity (km/decade)') +
      ylab('Biotic velocity (km/decade)') +
      labs(title = paste0(varname, '; ', month)) +
      geom_text(data=BV_CV_slopes_text, x=Inf,y=Inf,hjust=1.2,vjust=1.5, 
                aes(label=annotateText), size=3)
    print(p) 
  }
  
  for (k in 1:12){
    
    month= months[k]
    
    print(paste0('BVS vs ', varname, '; ', month, '; split'))
    
    BV_CV_slopes_climvar = BV_CV_slopes[which((BV_CV_slopes$clim_var == varname) & 
                                                (BV_CV_slopes$month==month) & 
                                                (BV_CV_slopes$type=='bvs')),]
    
    BV_CV_slopes_text = data.frame(
      taxon=BV_CV_slopes_climvar$taxon,
      xpos = rep(Inf, length(taxa)),
      ypos =  rep(Inf, length(taxa)),
      annotateText = paste0('pH: ', round(BV_CV_slopes_climvar$preHol,2), 
                            '; H: ', round(BV_CV_slopes_climvar$Hol,2)),
      hjustvar = rep(1, length(taxa)),
      vjustvar = rep(1, length(taxa))) 
    
    p = ggplot(data=Vs_climvar[which(Vs_climvar$month == month),]) +
      geom_point(aes(x=cv_mean, y= bvs_mean, colour=split_holocene)) +
      geom_linerange(aes(x=cv_mean, ymin=bvs_lo, ymax=bvs_hi, colour=split_holocene)) +
      geom_smooth(aes(x=cv_mean, y= bvs_mean, colour=split_holocene), method=lm) +
      # facet_grid(~taxon) +
      facet_wrap(~taxon, scales = 'free') +
      theme_bw(18) +
      xlab('Climate velocity (km/decade)') +
      ylab('Biotic velocity (km/decade)') +
      labs(title = paste0(varname, '; ', month)) +
      geom_text(data=BV_CV_slopes_text, x=Inf,y=Inf,hjust=1.2,vjust=1.5, 
                aes(label=annotateText), size=3)
    print(p) 
  }
  
}
dev.off()

for (i in 1:N_climvars){
  
  varname = name_pair$varname[i]
  
  print(paste0('Variable: ', varname))
  
  Vs_climvar = subset(Vs, clim_var == varname)
  
  pdf(paste0('figures/BVC_vs_CV_mean_scatter_credible_taxon_', varname, '.pdf'))
  
  for (k in 1:12){
    
    month= months[k]
    
    print(paste0('BVC vs ', varname, '; ', month))
    
    p = ggplot(data=Vs_climvar[which(Vs_climvar$month == month),]) +
      geom_point(aes(x=cv_mean, y= bvc_mean)) +
      geom_linerange(aes(x=cv_mean, ymin=bvc_lo, ymax=bvc_hi)) +
      geom_smooth(aes(x=cv_mean, y= bvc_mean), method=lm) +
      # facet_grid(~taxon) +
      facet_wrap(~taxon, scales = 'free') +
      theme_bw(14) +
      xlab('Climate velocity (km/decade)') +
      ylab('Biotic velocity (km/decade)') +
      labs(title = paste0(varname, '; ', month))
    # ggsave(pasteo('figures/BVC_vs_CV_mean_scatter_credible_taxon_', varname, '.pdf'), width=10, height=8)
    # ggsave(pasteo('figures/BVC_vs_CV_mean_scatter_credible_taxon_', varname, '.pdf'), width=10, height=8)
    print(p) 
  }
  dev.off()
  
  pdf(paste0('figures/BVC_vs_CV_mean_scatter_credible_taxon_split_', varname, '.pdf'), width=12, height=10)
  
  for (k in 1:12){
    
    month= months[k]
    
    print(paste0('BVC vs ', varname, '; ', month, '; split'))
    
    p = ggplot(data=Vs_climvar[which(Vs_climvar$month == month),]) +
      geom_point(aes(x=cv_mean, y= bvc_mean, colour=split_holocene)) +
      geom_linerange(aes(x=cv_mean, ymin=bvc_lo, ymax=bvc_hi, colour=split_holocene)) +
      geom_smooth(aes(x=cv_mean, y= bvc_mean, colour=split_holocene), method=lm) +
      # facet_grid(~taxon) +
      facet_wrap(~taxon, scales = 'free') +
      theme_bw(18) +
      xlab('Climate velocity (km/decade)') +
      ylab('Biotic velocity (km/decade)') +
      labs(title = paste0(varname, '; ', month))
    # ggsave(pasteo('figures/BVC_vs_CV_mean_scatter_credible_taxon_', varname, '.pdf'), width=10, height=8)
    # ggsave(pasteo('figures/BVC_vs_CV_mean_scatter_credible_taxon_', varname, '.pdf'), width=10, height=8)
    print(p) 
  }
  dev.off()
  
  pdf(paste0('figures/BVN_vs_CV_mean_scatter_credible_taxon_', varname, '.pdf'))
  
  for (k in 1:12){
    
    month= months[k]
    
    print(paste0('BVN vs ', varname, '; ', month))
    
    p = ggplot(data=Vs_climvar[which(Vs_climvar$month == month),]) +
      geom_point(aes(x=cv_mean, y= bvn_mean)) +
      geom_linerange(aes(x=cv_mean, ymin=bvn_lo, ymax=bvn_hi)) +
      geom_smooth(aes(x=cv_mean, y= bvn_mean), method=lm) +
      # facet_grid(~taxon) +
      facet_wrap(~taxon, scales = 'free') +
      theme_bw(14) +
      xlab('Climate velocity (km/decade)') +
      ylab('Biotic velocity (km/decade)') +
      labs(title = paste0(varname, '; ', month))
    # ggsave(pasteo('figures/BVC_vs_CV_mean_scatter_credible_taxon_', varname, '.pdf'), width=10, height=8)
    # ggsave(pasteo('figures/BVC_vs_CV_mean_scatter_credible_taxon_', varname, '.pdf'), width=10, height=8)
    print(p) 
  }
  dev.off()
  
  pdf(paste0('figures/BVN_vs_CV_mean_scatter_credible_taxon_split_', varname, '.pdf'), width=12, height=10)
  
  for (k in 1:12){
    
    month= months[k]
    
    print(paste0('BVN vs ', varname, '; ', month, '; split'))
    
    p = ggplot(data=Vs_climvar[which(Vs_climvar$month == month),]) +
      geom_point(aes(x=cv_mean, y= bvn_mean, colour=split_holocene)) +
      geom_linerange(aes(x=cv_mean, ymin=bvn_lo, ymax=bvn_hi, colour=split_holocene)) +
      geom_smooth(aes(x=cv_mean, y= bvn_mean, colour=split_holocene), method=lm) +
      # facet_grid(~taxon) +
      facet_wrap(~taxon, scales = 'free') +
      theme_bw(18) +
      xlab('Climate velocity (km/decade)') +
      ylab('Biotic velocity (km/decade)') +
      labs(title = paste0(varname, '; ', month))
    # ggsave(pasteo('figures/BVC_vs_CV_mean_scatter_credible_taxon_', varname, '.pdf'), width=10, height=8)
    # ggsave(pasteo('figures/BVC_vs_CV_mean_scatter_credible_taxon_', varname, '.pdf'), width=10, height=8)
    print(p) 
  }
  dev.off()
  
  
  pdf(paste0('figures/BVS_vs_CV_mean_scatter_credible_taxon_', varname, '.pdf'))
  
  for (k in 1:12){
    
    month= months[k]
    
    print(paste0('BVS vs ', varname, '; ', month))
    
    p = ggplot(data=Vs_climvar[which(Vs_climvar$month == month),]) +
      geom_point(aes(x=cv_mean, y= bvs_mean)) +
      geom_linerange(aes(x=cv_mean, ymin=bvs_lo, ymax=bvs_hi)) +
      geom_smooth(aes(x=cv_mean, y= bvs_mean), method=lm) +
      # facet_grid(~taxon) +
      facet_wrap(~taxon, scales = 'free') +
      theme_bw(14) +
      xlab('Climate velocity (km/decade)') +
      ylab('Biotic velocity (km/decade)') +
      labs(title = paste0(varname, '; ', month))
    # ggsave(pasteo('figures/BVC_vs_CV_mean_scatter_credible_taxon_', varname, '.pdf'), width=10, height=8)
    # ggsave(pasteo('figures/BVC_vs_CV_mean_scatter_credible_taxon_', varname, '.pdf'), width=10, height=8)
    print(p) 
  }
  dev.off()
  
  pdf(paste0('figures/BVS_vs_CV_mean_scatter_credible_taxon_split_', varname, '.pdf'), width=12, height=10)
  
  for (k in 1:12){
    
    month= months[k]
    
    print(paste0('BVS vs ', varname, '; ', month, '; split'))
    
    p = ggplot(data=Vs_climvar[which(Vs_climvar$month == month),]) +
      geom_point(aes(x=cv_mean, y= bvs_mean, colour=split_holocene)) +
      geom_linerange(aes(x=cv_mean, ymin=bvs_lo, ymax=bvs_hi, colour=split_holocene)) +
      geom_smooth(aes(x=cv_mean, y= bvs_mean, colour=split_holocene), method=lm) +
      # facet_grid(~taxon) +
      facet_wrap(~taxon, scales = 'free') +
      theme_bw(14) +
      xlab('Climate velocity (km/decade)') +
      ylab('Biotic velocity (km/decade)') +
      labs(title = paste0(varname, '; ', month))
    # ggsave(pasteo('figures/BVC_vs_CV_mean_scatter_credible_taxon_', varname, '.pdf'), width=10, height=8)
    # ggsave(pasteo('figures/BVC_vs_CV_mean_scatter_credible_taxon_', varname, '.pdf'), width=10, height=8)
    print(p) 
  }
  dev.off()
  
}

types = c('bvc', 'bvn', 'bvs')



# pdf('figures/BV_vs_CV_linear_slope_scatter_preHol_Hol_by_type1.pdf')#, width=10, height=8)
# 
# for (i in 1:N_climvars){
#   
#   varname = name_pair$varname[i]
#   
#   print(paste0('Variable: ', varname))
#   
#   BV_CV_slopes_climvar = subset(BV_CV_slopes, clim_var == varname)
#   
#   
#   for (j in 1:3){
#     
#     type = types[j]
#     
#     BV_CV_slopes_climvar_type = BV_CV_slopes_climvar[which(BV_CV_slopes_climvar$type == type),]
#     
#     p = ggplot(data=BV_CV_slopes_climvar)+
#       geom_point(aes(x = preHol, y=Hol, colour=sign_case_preHol_Hol)) +#, breaks = seq(-60, 200, by=15)) + 
#       # facet_grid(month~type, scales = 'free') +
#       # facet_grid(month~type) +
#       # facet_wrap(~month) + #, scales='free') +
#       geom_vline(xintercept = 0, linetype = 1, colour='dodgerblue', alpha=0.5) + 
#       geom_hline(yintercept = 0, linetype = 1, colour='dodgerblue', alpha=0.5) + 
#       theme_bw(16) +
#       xlab('Pre-Holocene linear slope') +
#       ylab('Holocene linear slope') +
#       labs(title = paste0(toupper(type), ' vs ', toupper(varname)), colour='PH / H') +
#       coord_equal()
#     print(p)
#   }
# }
# dev.off()

pdf('figures/BV_vs_CV_linear_slope_scatter_preHol_Hol_by_type.pdf')#, width=10, height=8)

for (j in 1:3){
  
  type = types[j]
  
  BV_CV_slopes_type = BV_CV_slopes[which(BV_CV_slopes$type == type),]
  
  for (i in 1:N_climvars){
    
    varname = name_pair$varname[i]
    
    print(paste0('Variable: ', varname))
    
    BV_CV_slopes_climvar_type = subset(BV_CV_slopes_type, clim_var == varname)
    
    
    
    # p = ggplot(data=BV_CV_slopes_climvar_type)+
    #   geom_point(aes(x = preHol, y=Hol)) +#, breaks = seq(-60, 200, by=15)) + 
    #   # facet_grid(month~type, scales = 'free') +
    #   # facet_grid(month~type) +
    #   facet_wrap(~month) +
    #   geom_vline(xintercept = 0, linetype = 1, colour='dodgerblue', alpha=0.5) + 
    #   geom_hline(yintercept = 0, linetype = 1, colour='dodgerblue', alpha=0.5) + 
    #   theme_bw(14) +
    #   xlab('Pre-Holocene linear slope') +
    #   ylab('Holocene linear slope') +
    #   labs(title = paste0(toupper(type), ' vs ', toupper(varname))) +
    #   coord_fixed()
    # print(p)
    
    p = ggplot(data=BV_CV_slopes_climvar_type)+
      geom_point(aes(x = preHol, y=Hol, colour=sign_case_preHol_Hol)) +#, breaks = seq(-60, 200, by=15)) + 
      # facet_grid(month~type, scales = 'free') +
      # facet_grid(month~type) +
      # facet_wrap(~month) + #, scales='free') +
      geom_vline(xintercept = 0, linetype = 1, colour='dodgerblue', alpha=0.5) + 
      geom_hline(yintercept = 0, linetype = 1, colour='dodgerblue', alpha=0.5) + 
      geom_abline(intercept = 0, slope=1, linetype = 2, colour='dodgerblue', alpha=0.5) + 
      theme_bw(16) +
      xlab('Pre-Holocene linear slope') +
      ylab('Holocene linear slope') +
      labs(title = paste0('outcome: ', toupper(type), '; predictor: ', toupper(varname)), colour='PH / H') +
      # coord_equal()
      coord_obs_pred()
    print(p)
    
    # p = ggplot(data=BV_CV_slopes_climvar)+
    #   geom_point(aes(x = preHol, y=Hol)) +#, breaks = seq(-60, 200, by=15)) + 
    #   # facet_grid(month~type, scales = 'free') +
    #   facet_grid(type~month) +
    #   geom_vline(xintercept = 0, linetype = 1, colour='dodgerblue', alpha=0.5) + 
    #   geom_hline(yintercept = 0, linetype = 1, colour='dodgerblue', alpha=0.5) + 
    #   theme_bw(14) +
    #   xlab('Pre-Holocene linear slope') +
    #   ylab('Holocene linear slope') +
    #   labs(title = varname) +
    #   coord_fixed()
    # print(p)
  }
}
dev.off()

pdf('figures/BV_vs_CV_linear_slope_scatter_preHol_Hol_by_taxon.pdf', width=10, height=8)

for (j in 1:3){
  
  type = types[j]
  
  BV_CV_slopes_type = BV_CV_slopes[which(BV_CV_slopes$type == type),]
  
  for (i in 1:N_climvars){
    
    varname = name_pair$varname[i]
    
    print(paste0('Variable: ', varname))
    
    BV_CV_slopes_climvar_type = subset(BV_CV_slopes_type, clim_var == varname)
    
    
    
    # p = ggplot(data=BV_CV_slopes_climvar_type)+
    #   geom_point(aes(x = preHol, y=Hol)) +#, breaks = seq(-60, 200, by=15)) + 
    #   # facet_grid(month~type, scales = 'free') +
    #   # facet_grid(month~type) +
    #   facet_wrap(~month) +
    #   geom_vline(xintercept = 0, linetype = 1, colour='dodgerblue', alpha=0.5) + 
    #   geom_hline(yintercept = 0, linetype = 1, colour='dodgerblue', alpha=0.5) + 
    #   theme_bw(14) +
    #   xlab('Pre-Holocene linear slope') +
    #   ylab('Holocene linear slope') +
    #   labs(title = paste0(toupper(type), ' vs ', toupper(varname))) +
    #   coord_fixed()
    # print(p)
    
    p = ggplot(data=BV_CV_slopes_climvar_type)+
      geom_point(aes(x = preHol, y=Hol, colour=sign_case_preHol_Hol)) +#, breaks = seq(-60, 200, by=15)) + 
      # geom_point(aes(x = preHol, y=Hol, colour=sign_case_preHol_Hol)) +#, breaks = seq(-60, 200, by=15)) + 
      # facet_grid(month~type, scales = 'free') +
      # facet_grid(month~type) +
      # facet_wrap(~month) + #, scales='free') +
      geom_vline(xintercept = 0, linetype = 1, colour='dodgerblue', alpha=0.5) + 
      geom_hline(yintercept = 0, linetype = 1, colour='dodgerblue', alpha=0.5) + 
      geom_abline(intercept = 0, slope=1, linetype = 2, colour='dodgerblue', alpha=0.5) + 
      theme_bw(14) +
      xlab('Pre-Holocene linear slope') +
      ylab('Holocene linear slope') +
      labs(title = paste0('outcome: ', toupper(type), '; predictor: ', toupper(varname)), colour='PH / H') +
      # coord_equal() +
      coord_obs_pred() +
      facet_wrap(~taxon)
    print(p)
    
    # p = ggplot(data=BV_CV_slopes_climvar)+
    #   geom_point(aes(x = preHol, y=Hol)) +#, breaks = seq(-60, 200, by=15)) + 
    #   # facet_grid(month~type, scales = 'free') +
    #   facet_grid(type~month) +
    #   geom_vline(xintercept = 0, linetype = 1, colour='dodgerblue', alpha=0.5) + 
    #   geom_hline(yintercept = 0, linetype = 1, colour='dodgerblue', alpha=0.5) + 
    #   theme_bw(14) +
    #   xlab('Pre-Holocene linear slope') +
    #   ylab('Holocene linear slope') +
    #   labs(title = varname) +
    #   coord_fixed()
    # print(p)
  }
}
dev.off()

pdf('figures/BV_vs_CV_linear_slope_scatter_preHol_Hol_by_taxon_type.pdf', width=10, height=8)

# for (j in 1:3){
  
  # type = types[j]
  # 
  # BV_CV_slopes_type = BV_CV_slopes[which(BV_CV_slopes$type == type),]
  
  for (i in 1:N_climvars){
    
    varname = name_pair$varname[i]
    
    print(paste0('Variable: ', varname))
    
    BV_CV_slopes_climvar_type = subset(BV_CV_slopes, clim_var == varname)
    
    p = ggplot(data=BV_CV_slopes_climvar_type)+
      geom_point(aes(x = preHol, y=Hol, colour=type)) +#, breaks = seq(-60, 200, by=15)) + 
      # facet_grid(month~type, scales = 'free') +
      # facet_grid(month~type) +
      # facet_wrap(~month) + #, scales='free') +
      geom_vline(xintercept = 0, linetype = 1, colour='dodgerblue', alpha=0.5) + 
      geom_hline(yintercept = 0, linetype = 1, colour='dodgerblue', alpha=0.5) + 
      theme_bw(14) +
      xlab('Pre-Holocene linear slope') +
      ylab('Holocene linear slope') +
      labs(title = paste0('outcome: BV; predictor: ', toupper(varname)), colour='type') +
      # coord_equal() +
      coord_obs_pred() +
      facet_wrap(~taxon)
    print(p)
    
    # p = ggplot(data=BV_CV_slopes_climvar)+
    #   geom_point(aes(x = preHol, y=Hol)) +#, breaks = seq(-60, 200, by=15)) + 
    #   # facet_grid(month~type, scales = 'free') +
    #   facet_grid(type~month) +
    #   geom_vline(xintercept = 0, linetype = 1, colour='dodgerblue', alpha=0.5) + 
    #   geom_hline(yintercept = 0, linetype = 1, colour='dodgerblue', alpha=0.5) + 
    #   theme_bw(14) +
    #   xlab('Pre-Holocene linear slope') +
    #   ylab('Holocene linear slope') +
    #   labs(title = varname) +
    #   coord_fixed()
    # print(p)
  }
# }
dev.off()


# # FIX ME
# pdf('figures/BV_vs_CV_linear_slope_histogram_preHol_Hol_by_type.pdf')#, width=10, height=8)
# 
# for (i in 1:N_climvars){
#   
#   varname = name_pair$varname[i]
#   
#   print(paste0('Variable: ', varname))
#   
#   
#   for (j in 1:3){
#     
#     type = types[j]
#     
#     BV_CV_slopes_climvar_type = BV_CV_slopes_climvar[which(BV_CV_slopes_climvar$type == type),]
#     
#     # p = ggplot(data=BV_CV_slopes_climvar_type)+
#     #   geom_point(aes(x = preHol, y=Hol)) +#, breaks = seq(-60, 200, by=15)) + 
#     #   # facet_grid(month~type, scales = 'free') +
#     #   # facet_grid(month~type) +
#     #   facet_wrap(~month) +
#     #   geom_vline(xintercept = 0, linetype = 1, colour='dodgerblue', alpha=0.5) + 
#     #   geom_hline(yintercept = 0, linetype = 1, colour='dodgerblue', alpha=0.5) + 
#     #   theme_bw(14) +
#     #   xlab('Pre-Holocene linear slope') +
#     #   ylab('Holocene linear slope') +
#     #   labs(title = paste0(toupper(type), ' vs ', toupper(varname))) +
#     #   coord_fixed()
#     # print(p)
#     
#     p = ggplot(data=BV_CV_slopes_climvar)+
#       geom_density(aes(x = preHol, y=Hol, colour=sign_case_preHol_Hol)) +#, breaks = seq(-60, 200, by=15)) + 
#       # facet_grid(month~type, scales = 'free') +
#       # facet_grid(month~type) +
#       # facet_wrap(~month) + #, scales='free') +
#       geom_vline(xintercept = 0, linetype = 1, colour='dodgerblue', alpha=0.5) + 
#       geom_hline(yintercept = 0, linetype = 1, colour='dodgerblue', alpha=0.5) + 
#       theme_bw(16) +
#       xlab('Pre-Holocene linear slope') +
#       ylab('Holocene linear slope') +
#       labs(title = paste0(toupper(type), ' vs ', toupper(varname)), colour='PH / H') +
#       coord_equal()
#     print(p)
#     
#     # p = ggplot(data=BV_CV_slopes_climvar)+
#     #   geom_point(aes(x = preHol, y=Hol)) +#, breaks = seq(-60, 200, by=15)) + 
#     #   # facet_grid(month~type, scales = 'free') +
#     #   facet_grid(type~month) +
#     #   geom_vline(xintercept = 0, linetype = 1, colour='dodgerblue', alpha=0.5) + 
#     #   geom_hline(yintercept = 0, linetype = 1, colour='dodgerblue', alpha=0.5) + 
#     #   theme_bw(14) +
#     #   xlab('Pre-Holocene linear slope') +
#     #   ylab('Holocene linear slope') +
#     #   labs(title = varname) +
#     #   coord_fixed()
#     # print(p)
#   }
# }
# dev.off()


pdf('figures/BV_vs_CV_linear_slope_difference_histogram_all21_Hol.pdf', width=10, height=8)

for (i in 1:N_climvars){
  
  varname = name_pair$varname[i]
  
  print(paste0('Variable: ', varname))
  
  BV_CV_slopes_climvar = subset(BV_CV_slopes, clim_var == varname)
  
  p = ggplot(data=BV_CV_slopes_climvar)+
    geom_histogram(aes(x = diff_all21_Hol, y=after_stat(density))) +#, breaks = seq(-60, 200, by=15)) + 
    facet_grid(month~type, scales = 'free') +
    geom_vline(xintercept = 0, linetype = 1, colour='dodgerblue', alpha=0.5) + 
    theme_bw(14) +
    xlab('Difference in linear slopes (All 21k - Holocene)') +
    labs(title = varname)
  print(p)
}
dev.off()

pdf('figures/BV_vs_CV_linear_slope_difference_histogram_preHol_Hol_by_var.pdf', width=10, height=8)
for (i in 1:N_climvars){
  
  varname = name_pair$varname[i]
  
  print(paste0('Variable: ', varname))
  
  BV_CV_slopes_climvar = subset(BV_CV_slopes, clim_var == varname)
  
  p = ggplot(data=BV_CV_slopes_climvar)+
    geom_histogram(aes(x = diff_preHol_Hol, y=after_stat(density))) +#, breaks = seq(-60, 200, by=15)) + 
    facet_grid(month~type, scales = 'free') +
    geom_vline(xintercept = 0, linetype = 1, colour='dodgerblue', alpha=0.5) + 
    theme_bw(14) +
    xlab('Difference in linear slopes (pre-Hol - Holocene)') +
    labs(title = varname)
  print(p)
}
dev.off()


pdf('figures/BV_vs_CV_linear_slope_difference_histogram_preHol_Hol_by_type.pdf', width=10, height=8)
for (i in 1:3){
  
  varname = name_pair$varname[i]
  
  print(paste0('Variable: ', varname))
  
  BV_CV_slopes_climvar = subset(BV_CV_slopes, type == types[i])
  BV_CV_slopes_climvar$sign = NA
  BV_CV_slopes_climvar$sign = sign(BV_CV_slopes_climvar$diff_preHol_Hol)
  
  p = ggplot(data=BV_CV_slopes_climvar)+
    geom_histogram(aes(x = diff_preHol_Hol, y=after_stat(density)), bins=15) +#, breaks = seq(-60, 200, by=15)) + 
    # geom_histogram(aes(x = diff_preHol_Hol, y=after_stat(density), fill=diff_preHol_Hol>0), bins=15) +#, breaks = seq(-60, 200, by=15)) + 
    facet_grid(month~clim_var, scales = 'free') +
    geom_vline(xintercept = 0, linetype = 1, colour='dodgerblue', alpha=0.5) + 
    theme_bw(14) +
    xlab('Difference in linear slopes (pre-Hol - Holocene)') +
    labs(title = types[i])
  print(p)
}
dev.off()


pdf('figures/BV_vs_CV_linear_slope_difference_density_preHol_Hol.pdf', width=10, height=8)

for (i in 1:N_climvars){
  
  varname = name_pair$varname[i]
  
  print(paste0('Variable: ', varname))
  
  BV_CV_slopes_climvar = subset(BV_CV_slopes, clim_var == varname)
  # BV_CV_slopes_climvar = BV_CV_slopes_climvar[, which(colnames(BV_CV_slopes_climvar) %in% c('taxon', 
  # 'type', 
  #                                                                                           'clim_var',
  #                                                                                           'month', 
  #                                                                                           'diff_all21_Hol',
  #                                                                                           'diff_preHol_Hol'))]
  # 
  # 
  # BV_CV_slopes_climvar_melt = melt(BV_CV_slopes_climvar, 
  #                                  id.vars = c('taxon', 'type', 'clim_var', 'month'))
  # 
  # p = ggplot(data=BV_CV_slopes_climvar)+
  #   # geom_histogram(aes(x = value, y=after_stat(density), colour=variable)) +#, breaks = seq(-60, 200, by=15)) + 
  #   geom_density(aes(x = diff_preHol_Hol, colour=same_sign_preHol_Hol, fill=same_sign_preHol_Hol), alpha=0.5) +#, breaks = seq(-60, 200, by=15)) + 
  #   facet_grid(month~type, scales = 'free') +
  #   geom_vline(xintercept = 0, linetype = 1, colour='grey35', alpha=0.5) + 
  #   theme_bw(14) +
  #   xlab('Difference in linear slopes (pre-Holocene - Holocene)') +
  #   labs(title = varname)
  # print(p)
  
  p = ggplot(data=BV_CV_slopes_climvar)+
    # geom_histogram(aes(x = value, y=after_stat(density), colour=variable)) +#, breaks = seq(-60, 200, by=15)) + 
    geom_density(aes(x = diff_preHol_Hol), colour='dodgerblue', fill='dodgerblue', alpha=0.5) +#, breaks = seq(-60, 200, by=15)) + 
    facet_grid(month~type, scales = 'free') +
    geom_vline(xintercept = 0, linetype = 1, colour='grey35', alpha=0.5) + 
    theme_bw(14) +
    xlab('Difference in linear slopes (pre-Holocene - Holocene)') +
    labs(title = varname)
  print(p)
}
dev.off()

pdf('figures/BV_vs_CV_linear_slope_difference_point_preHol_Hol.pdf', width=10, height=8)

for (i in 1:N_climvars){
  
  varname = name_pair$varname[i]
  
  print(paste0('Variable: ', varname))
  
  BV_CV_slopes_climvar = subset(BV_CV_slopes, clim_var == varname)
  # BV_CV_slopes_climvar = BV_CV_slopes_climvar[, which(colnames(BV_CV_slopes_climvar) %in% c('taxon', 
  # 'type', 
  #                                                                                           'clim_var',
  #                                                                                           'month', 
  #                                                                                           'diff_all21_Hol',
  #                                                                                           'diff_preHol_Hol'))]
  # 
  # 
  # BV_CV_slopes_climvar_melt = melt(BV_CV_slopes_climvar, 
  #                                  id.vars = c('taxon', 'type', 'clim_var', 'month'))
  # 
  p = ggplot(data=BV_CV_slopes_climvar)+
    # geom_histogram(aes(x = value, y=after_stat(density), colour=variable)) +#, breaks = seq(-60, 200, by=15)) + 
    geom_point(aes(x = diff_preHol_Hol, y=month, colour=taxon), alpha=0.6, size=2) +#, breaks = seq(-60, 200, by=15)) +
    # geom_point(aes(x=month, y=diff_preHol_Hol, colour=taxon), alpha=0.5) +#, breaks = seq(-60, 200, by=15)) + 
    facet_grid(type~., scales = 'free') +
    geom_vline(xintercept = 0, linetype = 1, colour='grey35', alpha=0.5) + 
    theme_bw(14) +
    xlab('Difference in linear slopes (pre-Holocene - Holocene)') +
    labs(title = varname)
  print(p)
}
dev.off()


# # new
# pdf('figures/BV_vs_CV_linear_slope_difference_histogram_preHol_Hol.pdf', width=10, height=8)
# 
# for (i in 1:N_climvars){
#   
#   varname = name_pair$varname[i]
#   
#   print(paste0('Variable: ', varname))
#   
#   BV_CV_slopes_climvar = subset(BV_CV_slopes, clim_var == varname)
#   BV_CV_slopes_climvar = BV_CV_slopes_climvar[, which(colnames(BV_CV_slopes_climvar) %in% c('taxon', 
#                                                                                             'type', 
#                                                                                             'clim_var',
#                                                                                             'month', 
#                                                                                             'preHol',
#                                                                                             'Hol',
#                                                                                             'sign_case_preHol_Hol'))]
#   BV_CV_slopes_climvar_melt = melt(BV_CV_slopes_climvar, 
#                                    id.vars = c('taxon', 'type', 'clim_var', 
#                                                'month', 'sign_case_preHol_Hol'))
#   
#   p = ggplot(data=BV_CV_slopes_climvar_melt)+
#     # geom_histogram(aes(x = value, y=after_stat(density), colour=variable)) +#, breaks = seq(-60, 200, by=15)) + 
#     geom_density(aes(x = value, colour=variable, fill=variable), alpha=0.5) +#, breaks = seq(-60, 200, by=15)) + 
#     facet_grid(month~type, scales = 'free') +
#     geom_vline(xintercept = 0, linetype = 1, colour='grey', alpha=0.5) + 
#     theme_bw(14) +
#     xlab('Difference in linear slopes (All 21k - Holocene)') +
#     labs(title = varname)
#   print(p)
#   
#   BV_CV_slopes_climvar_melt_sub = subset(BV_CV_slopes_climvar_melt, type=='bvc')
#   
#   table(BV_CV_slopes_climvar_melt_sub$sign_case_preHol_Hol, BV_CV_slopes_climvar_melt_sub$taxon)
#   
#   p = ggplot(data=BV_CV_slopes_climvar_melt_sub)+
#     # geom_histogram(aes(x = value, y=after_stat(density), colour=variable)) +#, breaks = seq(-60, 200, by=15)) + 
#     geom_density(aes(x = value, colour=variable, fill=variable), alpha=0.5) +#, breaks = seq(-60, 200, by=15)) + 
#     facet_grid(month~sign_case_preHol_Hol, scales = 'free') +
#     geom_vline(xintercept = 0, linetype = 1, colour='grey', alpha=0.5) + 
#     theme_bw(14) +
#     xlab('Linear slopes') +
#     labs(title = varname)
#   print(p)
# }
# dev.off()
# 
# pdf('figures/BV_vs_CV_linear_slope_difference_histogram_all21k_preHol_Hol.pdf', width=10, height=8)
# 
# for (i in 1:N_climvars){
#   
#   varname = name_pair$varname[i]
#   
#   print(paste0('Variable: ', varname))
#   
#   BV_CV_slopes_climvar = subset(BV_CV_slopes, clim_var == varname)
#   BV_CV_slopes_climvar = BV_CV_slopes_climvar[, which(colnames(BV_CV_slopes_climvar) %in% c('taxon', 
#                                                                                             'type', 
#                                                                                             'clim_var',
#                                                                                             'month', 
#                                                                                             'diff_all21_Hol',
#                                                                                             'diff_preHol_Hol'))]
#   
#   
#   BV_CV_slopes_climvar_melt = melt(BV_CV_slopes_climvar, 
#                                    id.vars = c('taxon', 'type', 'clim_var', 'month'))
#   
#   p = ggplot(data=BV_CV_slopes_climvar_melt)+
#     # geom_histogram(aes(x = value, y=after_stat(density), colour=variable)) +#, breaks = seq(-60, 200, by=15)) + 
#     geom_density(aes(x = value, colour=variable, fill=variable), alpha=0.5) +#, breaks = seq(-60, 200, by=15)) + 
#     facet_grid(month~type, scales = 'free') +
#     geom_vline(xintercept = 0, linetype = 1, colour='grey', alpha=0.5) + 
#     theme_bw(14) +
#     xlab('Difference in linear slopes (All 21k - Holocene)') +
#     labs(title = varname)
#   print(p)
# }
# dev.off()


BV_CV_slopes$All_21k.sig = 0
BV_CV_slopes$All_21k.sig[which(BV_CV_slopes$All_21k.p.value<0.05)] = 1

BV_CV_slopes$Hol.sig = 0
BV_CV_slopes$Hol.sig[which(BV_CV_slopes$Hol.p.value<0.05)] = 1

BV_CV_slopes$preHol.sig = 0
BV_CV_slopes$preHol.sig[which(BV_CV_slopes$preHol.p.value<0.05)] = 1

# ggplot(data=BV_CV_slopes) +
#   geom_point(aes(x = All_21k, 
#                  y = taxon, 
#                  shape = factor(All_21k.sig), 
#                  colour=type), 
#              size=2) +
#   facet_grid(clim_var~month, scales = 'free')
# 
# ggplot(data=BV_CV_slopes) +
#   geom_point(aes(x = All_21k, 
#                  y = type, 
#                  shape = factor(All_21k.sig), 
#                  colour=taxon), 
#              size=2,
#              alpha=0.5) +
#   facet_grid(month~clim_var, scales='free')
# 
# ggplot(data=BV_CV_slopes) +
#   geom_histogram(aes(x = All_21k, 
#                      colour = type,
#                      fill = type), 
#                  alpha=0.5)

# for (i in 1:N_climvars){
#   
#   varname = name_pair$varname[i]
#   
#   print(paste0('Variable: ', varname))
#   
#   BV_CV_slopes_climvar = subset(BV_CV_slopes, clim_var == varname)
#   BV_CV_slopes_climvar = BV_CV_slopes_climvar[,c('taxon', 'type', 'clim_var', 'month',
#                                                  'All_21k', 'Hol', 'preHol')]
#   
#   BV_CV_slopes_climvar_melt = melt(BV_CV_slopes_climvar, 
#                                    id.vars = c('taxon', 'type', 'clim_var', 'month'))
#   
#   ggplot(data=BV_CV_slopes_climvar_melt) +
#     geom_histogram(aes(x = variable, y=after_stat(density)), bins = 15, breaks = seq(-0.6, 2, by = 0.1),
#                    alpha=0.5, colour='blue', fill='lightblue') +
#     geom_vline(xintercept = 0, colour = 'black', lty = 2) +
#     facet_grid(type~month) +
#     theme_bw() +
#     theme(text = element_text(size=14))
#   # ggsave('figures/BV_vs_CV_slopes_histogram_LGM.pdf', width=10, height=8)
#   # ggsave('figures/BV_vs_CV_slopes_histogram_LGM.png', width=10, height=8)
#   
#   ggplot(data=BV_CV_slopes) +
#     geom_histogram(aes(x = PostLGM, y=after_stat(density)), bins = 15, breaks = seq(-0.6, 2, by = 0.1),
#                    alpha=0.5, colour='blue', fill='lightblue') +
#     geom_vline(xintercept = 0, colour = 'black', lty = 2) +
#     facet_grid(type~.) +
#     theme_bw() +
#     theme(text = element_text(size=14))
#   # ggsave('figures/BV_vs_CV_slopes_histogram_postLGM.pdf', width=10, height=8)
#   # ggsave('figures/BV_vs_CV_slopes_histogram_postLGM.png', width=10, height=8)
#   
# }

# ggplot(data=BV_CV_slopes) +
#   geom_histogram(aes(x = Holocene, y=after_stat(density)), bins = 15, breaks = seq(-0.6, 2, by = 0.1),
#                  alpha=0.5, colour='blue', fill='lightblue') +
#   geom_vline(xintercept = 0, colour = 'black', lty = 2) +
#   facet_grid(H.sig~type) +
#   theme_bw() +
#   theme(text = element_text(size=14))
# 
# ggplot(data=BV_CV_slopes) +
#   geom_histogram(aes(x = Holocene, y=after_stat(density)), bins = 15, breaks = seq(-0.6, 2, by = 0.1),
#                  alpha=0.5, colour='blue', fill='lightblue') +
#   geom_vline(xintercept = 0, colour = 'black', lty = 2) +
#   facet_grid(type~.) +
#   theme_bw() +
#   theme(text = element_text(size=14))
# ggsave('figures/BV_vs_CV_slopes_histogram_holocene.pdf', width=10, height=8)
# ggsave('figures/BV_vs_CV_slopes_histogram_holocene.png', width=10, height=8)
# 
# ggplot(data=BV_CV_slopes) +
#   geom_histogram(aes(x = diff_all_H, y=after_stat(density)), bins = 15, breaks = seq(-0.50, 2, by = 0.1),
#                  alpha=0.5, colour='blue', fill='lightblue') +
#   geom_vline(xintercept = 0, colour = 'black', lty = 2) +
#   facet_grid(type~.) +
#   theme_bw() +
#   theme(text = element_text(size=14))
# ggsave('figures/BV_vs_CV_slopes_histogram_holocene_LGM_diff.pdf', width=10, height=8)
# ggsave('figures/BV_vs_CV_slopes_histogram_holocene_LGM_diff.png', width=10, height=8)
# 
# ggplot(data=BV_CV_slopes) +
#   geom_histogram(aes(x = diff_LGM_H, y=after_stat(density)), bins = 15, breaks = seq(-0.50, 2, by = 0.1),
#                  alpha=0.5, colour='blue', fill='lightblue') +
#   geom_vline(xintercept = 0, colour = 'black', lty = 2) +
#   facet_grid(type~.) +
#   theme_bw() +
#   theme(text = element_text(size=14))
# ggsave('figures/BV_vs_CV_slopes_histogram_holocene_LGM_diff.pdf', width=10, height=8)
# ggsave('figures/BV_vs_CV_slopes_histogram_holocene_LGM_diff.png', width=10, height=8)

# ggplot(data=slopes) +
#   geom_density(aes(x = estimate), 
#                  alpha=0.5) +
#   facet_grid(type~.)


ggplot(data=BV_CV_slopes) +
  geom_point(aes(x = PostLGM,
                 y = Holocene,
                 colour=type),
             size=2) +
  geom_smooth(aes(x=PostLGM, y= Holocene, colour=type), method=lm) +
  # facet_wrap(~taxon, scales = 'free_y') +
  theme_bw(14) +
  xlab('Post LGM slopes') +
  ylab('Holocene slopes') +
  coord_fixed()


BV_CV_slopes[which(BV_CV_slopes$All_21k < 0),]


BV_CV_slopes_long = BV_CV_slopes[,c('taxon', 'type', 'All_21k', 'PostLGM', 'Holocene')] %>% 
  pivot_longer(cols = c('All_21k', 'PostLGM', 'Holocene'))

ggplot(data=BV_CV_slopes_long) +
  geom_histogram(aes(x=value, y=after_stat(density))) +
  facet_grid(name~type)

ggplot(data=BV_CV_slopes_long) +
  geom_density(aes(x=value, y=after_stat(density)), colour='dodgerblue', fill='dodgerblue') +
  facet_grid(name~type)

ggplot(data=subset(BV_CV_slopes_long, name %in% c('Holocene', 'PostLGM'))) +
  geom_density(aes(x=value, y=after_stat(density), colour=name, fill=name), alpha=0.4) +
  facet_grid(type~.) +
  theme_bw(14) +
  xlab('Regression slopes BV vs CV')
ggsave('figures/BV_vs_CV_slopes_density_split.pdf', width=10, height=8)
ggsave('figures/BV_vs_CV_slopes_density_split.png', width=10, height=8)

ggplot(data=BV_CV_slopes_long) +
  geom_density(aes(x=value, y=after_stat(density), colour=name, fill=name), alpha=0.4) +
  facet_grid(type~.) +
  theme_bw(14) +
  xlab('Regression slopes BV vs CV')
ggsave('figures/BV_vs_CV_slopes_density_all_split.pdf', width=10, height=8)
ggsave('figures/BV_vs_CV_slopes_density_all_split.png', width=10, height=8)


# ######################################################################################################################
# ## climate velocity masked by taxon
# ######################################################################################################################
# 
# BV_colnames = as.vector(sapply(c('bvc', 'bvn', 'bvs'), function(x) paste0(x, c('_mean', '_median', '_lo', '_hi'))))
# pollen_BVs_summary[, BV_colnames] = pollen_BVs_summary[, BV_colnames]/100
# 
# 
# # CVs = readRDS('output/climate_velocity.RDS')
# # CVs = readRDS('output/climate_velocity_taxon.RDS')
# CVs = readRDS('output/TMAX_JUNE_CV_taxon.RDS')
# CVs$voccMag = CVs$voccMag / 100
# CVs$voccMag_pos = abs(CVs$voccMag)
# 
# # CVs_summary = CVs %>% 
# #   group_by(period, taxon) %>% 
# #   dplyr::summarize(cv_mean = mean(voccMag, na.rm=TRUE),
# #                    cv_median = median(voccMag, na.rm=TRUE),
# #                    cv_sd = sd(voccMag, na.rm=TRUE),
# #                    cv_lo = quantile(voccMag, c(0.05), na.rm=TRUE),
# #                    cv_hi = quantile(voccMag, c(0.95), na.rm=TRUE), .groups="keep")
# CVs_summary = CVs %>%
#   group_by(period, taxon) %>%
#   dplyr::summarize(cv_mean = mean(voccMag_pos, na.rm=TRUE),
#                    cv_median = median(voccMag_pos, na.rm=TRUE),
#                    cv_sd = sd(voccMag_pos, na.rm=TRUE),
#                    cv_lo = quantile(voccMag_pos, c(0.05), na.rm=TRUE),
#                    cv_hi = quantile(voccMag_pos, c(0.95), na.rm=TRUE), .groups="keep")
# 
# ggplot(data=CVs_summary) +
#   geom_point(aes(x=period, y=cv_median)) +
#   geom_linerange(aes(x=period, ymin=cv_lo, ymax=cv_hi)) +
#   ylab('Climate velocity (km/decade)') +
#   xlab('Time (YBP)') +
#   theme_bw(14) +
#   theme(axis.text.x = element_text(angle=90,hjust=1)) +
#   facet_wrap(~taxon)
# 
# 
# ggplot(data=CVs_summary) +
#   geom_point(aes(x=cv_lo, y=cv_hi)) +
#   # geom_linerange(aes(x=period, ymin=cv_lo, ymax=cv_hi)) +
#   # ylab('Climate velocity (km/decade)') +
#   # xlab('Time (YBP)') +
#   # theme_bw(14) +
#   # theme(axis.text.x = element_text(angle=90,hjust=1)) +
#   facet_wrap(~period)
# 
# ggplot(data=CVs) +
#   geom_boxplot(aes(x=factor(period), y=voccMag), outlier.shape = NA) +
#   # geom_linerange(aes(x=period, ymin=cv_lo, ymax=cv_hi)) +
#   ylab('Climate velocity (km/decade)') +
#   xlab('Time (YBP)') +
#   # ylim(c(-3,3)) +
#   theme_bw(12) +
#   theme(axis.text.x = element_text(angle=45,hjust=1, size=8)) +
#   geom_hline(yintercept=0, colour='grey', lty=2, lwd=1) +
#   scale_x_discrete(limits=rev) +
#   facet_wrap(~taxon)
# ggsave('figures/CV_vs_time_taxon_boxplot.pdf', width=10, height=6)
# ggsave('figures/CV_vs_time_taxon_boxplot.png', width=10, height=6)
# 
# # ggplot(data=CVs) +
# #   geom_boxplot(aes(x=factor(period), y=voccMag), outlier.shape = NA) +
# #   # geom_linerange(aes(x=period, ymin=cv_lo, ymax=cv_hi)) +
# #   ylab('Climate velocity (km/decade)') +
# #   xlab('Time (YBP)') +
# #   # ylim(c(-3,3)) +
# #   theme_bw(16) +
# #   theme(axis.text.x = element_text(angle=45,hjust=1, size=14)) +
# #   geom_hline(yintercept=0, colour='grey', lty=2, lwd=2) +
# #   scale_x_discrete(limits=rev) +
# #   facet_wrap(~taxon)
# # ggsave('figures/CV_vs_time_boxplot.pdf', width=10, height=6)
# # ggsave('figures/CV_vs_time_boxplot.png', width=10, height=6)
# 
# 
# PC_first_two$period = -(PC_first_two$period)
# 
# PC_CV = merge(PC_first_two, CVs_summary, by.x = 'period', by.y = 'period')
# 
# PC_CV_wide = PC_CV[,colnames(PC_CV) %in% c('period', 'variable', 'type', 'value', 'cv_mean')] %>% 
#   pivot_wider(id_cols = c('period', 'cv_mean'), values_from = 'value', names_from = c('variable', 'type')) 
# 
# PC_CV_cor = cor(PC_CV_wide[,2:ncol(PC_CV_wide)])
# 
# corrplot(BVS_cor)
# 
# 
# BVS_cor_p = rcorr(as.matrix(BVS_taxa_wide[,2:ncol(BVS_taxa_wide)]))
# BVS_r = BVS_cor_p$r
# BVS_p = BVS_cor_p$P
# BVS_p[which(is.na(BVS_p))] = 0
# 
# corrplot(BVS_r, 
#          order = 'FPC', 
#          p.mat = BVS_p, 
#          insig = 'blank',
#          type = 'lower')
# 
# # BVs between two consecutive time bins 990 years apart
# # 
# 
# 
# Vs = merge(pollen_BVs_summary, 
#            CVs_summary, by = c('period', 'taxon'))
# Vs$era_name = NA
# Vs$era_name[which(Vs$period>18000)] = 'near LGM'
# Vs$era_name[which((Vs$period<18000)&(Vs$period>12000))] = 'post LGM'
# Vs$era_name[which((Vs$period<12000)&(Vs$period>2000))] = 'pre-industrial Holocene'
# Vs$era_name[which(Vs$period<2000)] = 'Anthropocene'
# 
# head(Vs)
# 
# # ggplot() +
# #   geom_point(data=Vs, aes(x=cv_mean, y= bv_mean/100, colour=taxon)) +
# #   geom_smooth(data=Vs, aes(x=cv_mean, y= bv_mean/100, colour=taxon), method=lm) +
# #   # facet_wrap(~taxon) +
# #   theme_bw(14) +
# #   xlab('Climate velocity (km/decade)') +
# #   ylab('Biotic velocity (km/decade)')
# # ggsave('figures/BV_vs_CV_mean_scatter.pdf', width=10, height=8)
# # 
# # ggplot() +
# #   geom_point(data=Vs, aes(x=cv_median, y= bv_median/100)) +
# #   geom_smooth(data=Vs, aes(x=cv_median, y= bv_median/100), method=lm) +
# #   facet_wrap(~taxon) +
# #   theme_bw(14) +
# #   xlab('Climate velocity (km/decade)') +
# #   ylab('Biotic velocity (km/decade)')
# # ggsave('figures/BV_vs_CV_median_scatter_taxon.pdf', width=10, height=8)
# 
# # ggplot(data=Vs) +
# #   geom_point(aes(x=cv_mean, y= bvc_mean/100)) +
# #   geom_smooth(aes(x=cv_mean, y= bvc_mean/100), method=lm) +
# #   facet_wrap(~taxon) +
# #   theme_bw(14) +
# #   xlab('Climate velocity (km/decade)') +
# #   ylab('Biotic velocity (km/decade)')
# # ggsave('figures/BVC_vs_CV_mean_scatter_taxon.pdf', width=10, height=8)
# 
# 
# ggplot() +
#   geom_point(data=Vs, aes(x=cv_mean, y= bvc_mean)) +
#   geom_smooth(data=Vs, aes(x=cv_mean, y= bvc_mean), method=lm) +
#   facet_wrap(~taxon, scales = 'free') +
#   theme_bw(14) +
#   xlab('Climate velocity (km/decade)') +
#   ylab('Biotic velocity (km/decade)')
# ggsave('figures/BVC_vs_CV_masked_mean_scatter_taxon.pdf', width=10, height=8)
# ggsave('figures/BVC_vs_CV_masked_mean_scatter_taxon.png', width=10, height=8)
# 
# ggplot(data=Vs) +
#   geom_point(aes(x=cv_mean, y= bvc_mean)) +
#   geom_linerange(aes(x=cv_mean, ymin=bvc_lo, ymax=bvc_hi)) +
#   geom_smooth(aes(x=cv_mean, y= bvc_mean), method=lm) +
#   facet_wrap(~taxon, scales = 'free') +
#   theme_bw(14) +
#   xlab('Climate velocity (km/decade)') +
#   ylab('Biotic velocity (km/decade)')
# ggsave('figures/BVC_vs_CV_mean_scatter_credible_taxon.pdf', width=10, height=8)
# ggsave('figures/BVC_vs_CV_mean_scatter_credible_taxon.png', width=10, height=8)
# 
# ggplot() +
#   geom_point(data=Vs, aes(x=cv_mean, y= bvc_median)) +
#   geom_smooth(data=Vs, aes(x=cv_mean, y= bvc_median), method=lm) +
#   facet_wrap(~taxon, scales = 'free') +
#   theme_bw(14) +
#   xlab('Climate velocity (km/decade)') +
#   ylab('Biotic velocity (km/decade)')
# ggsave('figures/BVC_vs_CV_masked_median_scatter_taxon.pdf', width=10, height=8)
# ggsave('figures/BVC_vs_CV_masked_median_scatter_taxon.png', width=10, height=8)
# 
# ggplot(data=Vs) +
#   geom_point(aes(x=cv_mean, y= bvc_median)) +
#   geom_linerange(aes(x=cv_mean, ymin=bvc_lo, ymax=bvc_hi)) +
#   geom_smooth(aes(x=cv_mean, y= bvc_median), method=lm) +
#   facet_wrap(~taxon, scales = 'free') +
#   theme_bw(14) +
#   xlab('Climate velocity (km/decade)') +
#   ylab('Biotic velocity (km/decade)')
# ggsave('figures/BVC_vs_CV_median_scatter_credible_taxon.pdf', width=10, height=8)
# ggsave('figures/BVC_vs_CV_median_scatter_credible_taxon.png', width=10, height=8)
# 
# ggplot(data=Vs) +
#   geom_point(aes(x=cv_mean, y= bvc_median)) +
#   geom_linerange(aes(x=cv_mean, ymin=bvc_lo, ymax=bvc_hi)) +
#   # geom_abline(intercept = 0, slope = 1, colour = "grey60") +
#   geom_smooth(aes(x=cv_mean, y= bvc_median), method=lm) +
#   facet_wrap(~taxon, scales = 'free') +
#   theme_bw(14) +
#   theme(aspect.ratio=1) +
#   xlab('Climate velocity (km/decade)') +
#   ylab('Biotic velocity (km/decade)') 
# ggsave('figures/BVC_vs_CV_median_scatter_credible_taxon.pdf', width=10, height=8)
# ggsave('figures/BVC_vs_CV_median_scatter_credible_taxon.png', width=10, height=8)
# 
# # 
# # foo = Vs %>%
# #   group_by(taxon) %>%
# #   summarize(max = max(bvc_mean, na.rm=TRUE)) 
# # colnames(foo) = c('taxon', 'bvc_mean')
# # 
# # merge(foo, Vs, by = c('taxon', 'bvc_mean'))
# 
# 
# ggplot(data=Vs) +
#   geom_point(aes(x=cv_mean, y= bvn_mean)) +
#   geom_linerange(aes(x=cv_mean, ymin=bvn_lo, ymax=bvn_hi)) +
#   geom_smooth(aes(x=cv_mean, y= bvn_mean), method=lm) +
#   facet_wrap(~taxon, scales = 'free') +
#   theme_bw(14) +
#   xlab('Climate velocity (km/decade)') +
#   ylab('Biotic velocity (km/decade)')
# ggsave('figures/BVN_vs_CV_mean_scatter_credible_taxon.pdf', width=10, height=8)
# ggsave('figures/BVN_vs_CV_mean_scatter_credible_taxon.png', width=10, height=8)
# 
# 
# ggplot(data=Vs) +
#   geom_point(aes(x=cv_mean, y= bvs_mean)) +
#   geom_linerange(aes(x=cv_mean, ymin=bvs_lo, ymax=bvs_hi)) +
#   geom_smooth(aes(x=cv_mean, y= bvs_mean), method=lm) +
#   facet_wrap(~taxon, scales = 'free') +
#   theme_bw(14) +
#   xlab('Climate velocity (km/decade)') +
#   ylab('Biotic velocity (km/decade)')
# ggsave('figures/BVS_vs_CV_mean_scatter_credible_taxon.pdf', width=10, height=8)
# ggsave('figures/BVS_vs_CV_mean_scatter_credible_taxon.png', width=10, height=8)
# 
# 
# ggplot(data=Vs) +
#   geom_point(aes(x=cv_mean, y= bvn_median)) +
#   geom_linerange(aes(x=cv_mean, ymin=bvn_lo, ymax=bvn_hi)) +
#   geom_smooth(aes(x=cv_mean, y= bvn_median), method=lm) +
#   facet_wrap(~taxon, scales = 'free') +
#   theme_bw(14) +
#   xlab('Climate velocity (km/decade)') +
#   ylab('Biotic velocity (km/decade)')
# ggsave('figures/BVN_vs_CV_median_scatter_credible_taxon.pdf', width=10, height=8)
# ggsave('figures/BVN_vs_CV_median_scatter_credible_taxon.png', width=10, height=8)
# 
# 
# ggplot(data=Vs) +
#   geom_point(aes(x=cv_mean, y= bvs_median)) +
#   geom_linerange(aes(x=cv_mean, ymin=bvs_lo, ymax=bvs_hi)) +
#   geom_smooth(aes(x=cv_mean, y= bvs_median), method=lm) +
#   facet_wrap(~taxon, scales = 'free_y') +
#   theme_bw(14) +
#   xlab('Climate velocity (km/decade)') +
#   ylab('Biotic velocity (km/decade)')
# ggsave('figures/BVS_vs_CV_median_scatter_credible_taxon.pdf', width=10, height=8)
# ggsave('figures/BVS_vs_CV_median_scatter_credible_taxon.png', width=10, height=8)
# 
# 
# ggplot(data=Vs) +
#   geom_point(aes(x=bvc_mean, y= bvn_mean), alpha=0.5) +
#   geom_linerange(aes(x=bvc_mean, ymin=bvn_lo, ymax=bvn_hi), alpha=0.5) +
#   geom_linerange(aes(y=bvn_mean, xmin=bvc_lo, xmax=bvc_hi), alpha=0.5) +
#   geom_smooth(aes(x=bvc_mean, y= bvn_mean), method=lm, fullrange=TRUE) +
#   facet_wrap(~taxon, scales = 'free') +
#   theme_bw(14) +
#   xlab('Centroid biotic velocity (km/decade)') +
#   ylab('Northern Biotic velocity (km/decade)')
# ggsave('figures/BVN_vs_BVC_mean_scatter_credible_taxon.pdf', width=10, height=8)
# ggsave('figures/BVN_vs_BVC_mean_scatter_credible_taxon.png', width=10, height=8)
# 
# ggplot(data=Vs) +
#   geom_point(aes(x=bvc_mean, y= bvs_mean), alpha=0.5) +
#   geom_linerange(aes(x=bvc_mean, ymin=bvs_lo, ymax=bvs_hi), alpha=0.5) +
#   geom_linerange(aes(y=bvs_mean, xmin=bvc_lo, xmax=bvc_hi), alpha=0.5) +
#   geom_smooth(aes(x=bvc_mean, y= bvs_mean), method=lm, fullrange=TRUE) +
#   facet_wrap(~taxon, scales = 'free') +
#   theme_bw(14) +
#   xlab('Centroid biotic velocity (km/decade)') +
#   ylab('Southern Biotic velocity (km/decade)')
# ggsave('figures/BVS_vs_BVC_mean_scatter_credible_taxon.pdf', width=10, height=8)
# ggsave('figures/BVS_vs_BVC_mean_scatter_credible_taxon.png', width=10, height=8)
# 
# 
# ggplot(data=Vs) +
#   geom_point(aes(x=bvs_mean, y= bvn_mean), alpha=0.5) +
#   geom_linerange(aes(x=bvs_mean, ymin=bvn_lo, ymax=bvn_hi), alpha=0.5) +
#   geom_linerange(aes(y=bvn_mean, xmin=bvs_lo, xmax=bvs_hi), alpha=0.5) +
#   geom_smooth(aes(x=bvs_mean, y= bvn_mean), method=lm, fullrange=TRUE) +
#   facet_wrap(~taxon, scales = 'free') +
#   theme_bw(14) +
#   xlab('Southern biotic velocity (km/decade)') +
#   ylab('Northern Biotic velocity (km/decade)') 
# ggsave('figures/BVN_vs_BVS_mean_scatter_credible_taxon.pdf', width=10, height=8)
# ggsave('figures/BVN_vs_BVS_mean_scatter_credible_taxon.png', width=10, height=8)
# 
# 
# Vs_melt_mean = melt(Vs[,c('period', 'era_name', 'taxon', 'bvc_mean', 'bvn_mean', 'bvs_mean', 'cv_mean')], 
#                     id.vars = c('period', 'era_name', 'taxon', 'cv_mean'))
# Vs_melt_median = melt(Vs[,c('period', 'era_name', 'taxon', 'bvc_median', 'bvn_median', 'bvs_median', 'cv_mean')], 
#                     id.vars = c('period', 'era_name', 'taxon', 'cv_mean'))
# Vs_melt_lo = melt(Vs[,c('period', 'era_name','taxon', 'bvc_lo', 'bvn_lo', 'bvs_lo', 'cv_mean')], 
#                   id.vars = c('period', 'era_name', 'taxon', 'cv_mean'))
# Vs_melt_hi = melt(Vs[,c('period', 'era_name', 'taxon', 'bvc_hi', 'bvn_hi', 'bvs_hi', 'cv_mean')], 
#                   id.vars = c('period', 'era_name', 'taxon', 'cv_mean'))
# 
# Vs_melt = data.frame(Vs_melt_mean[,1:5], 
#                      type = substr(Vs_melt_mean$variable, 1, 3),
#                      bv_mean = Vs_melt_mean$value,
#                      bv_median = Vs_melt_median$value,
#                      bv_lo = Vs_melt_lo$value,
#                      bv_hi = Vs_melt_hi$value)
# 
# Vs_melt$split_name = NA
# Vs_melt$split_name[which(Vs_melt$period <= 12000)] = 'Holocene'
# Vs_melt$split_name[which(Vs_melt$period > 12000)] = 'pre_Holocene'
# 
# 
# ggplot(data=Vs_melt) +
#   geom_point(aes(x=cv_mean, y= bv_mean, colour=type), alpha=0.5) +
#   geom_linerange(aes(x=cv_mean, ymin=bv_lo, ymax=bv_hi, colour=type), alpha=0.5) +
#   geom_smooth(aes(x=cv_mean, y= bv_mean, colour=type, group=type), method=lm, fullrange=TRUE, alpha=0.5) +
#   facet_wrap(~taxon, scales = 'free_y') +
#   theme_bw(14) +
#   xlab('Climate velocity (km/decade)') +
#   ylab('Biotic velocity (km/decade)')
# ggsave('figures/BVs_vs_CV_mean_scatter_credible_taxon.pdf', width=10, height=8)
# ggsave('figures/BVs_vs_CV_mean_scatter_credible_taxon.png', width=10, height=8)
# 
# 
# ggplot(data=subset(Vs_melt, type =='bvc')) +
#   geom_point(aes(x=cv_mean, y= bv_mean, colour=era_name)) +
#   geom_linerange(aes(x=cv_mean, ymin=bv_lo, ymax=bv_hi, colour=era_name)) +
#   geom_smooth(aes(x=cv_mean, y= bv_mean), method=lm) +
#   facet_wrap(~taxon, scales = 'free_y') +
#   theme_bw(14) +
#   xlab('Climate velocity (km/decade)') +
#   ylab('Biotic velocity (km/decade)')
# ggsave('figures/BVC_vs_CV_mean_scatter_era_credible_taxon.pdf', width=10, height=8)
# ggsave('figures/BVC_vs_CV_mean_scatter_era_credible_taxon.png', width=10, height=8)
# 
# 
# ggplot(data=subset(Vs_melt, type =='bvc')) +
#   geom_point(aes(x=cv_mean, y= bv_median, colour=era_name)) +
#   geom_linerange(aes(x=cv_mean, ymin=bv_lo, ymax=bv_hi, colour=era_name)) +
#   geom_smooth(aes(x=cv_mean, y = bv_median), method=lm, colour="grey40", fill='grey60', alpha=0.4) +
#   facet_wrap(~taxon, scales = 'free_y') +
#   theme_bw(14) +
#   xlab('Climate velocity (km/decade)') +
#   ylab('Biotic velocity (km/decade)')
# ggsave('figures/BVC_vs_CV_median_scatter_era_credible_taxon.pdf', width=10, height=8)
# ggsave('figures/BVC_vs_CV_median_scatter_era_credible_taxon.png', width=10, height=8)
# 
# 
# 
# ggplot(data=subset(Vs_melt, (period<14000)&(type =='bvc'))) +
#   geom_point(aes(x=cv_mean, y= bv_mean)) +
#   geom_linerange(aes(x=cv_mean, ymin=bv_lo, ymax=bv_hi)) +
#   geom_smooth(aes(x=cv_mean, y= bv_mean), method=lm) +
#   facet_wrap(~taxon, scales = 'free_y') +
#   theme_bw(14) +
#   xlab('Climate velocity (km/decade)') +
#   ylab('Biotic velocity (km/decade)')
# ggsave('figures/BVC_vs_CV_mean_scatter_credible_taxon_14000_mod.pdf', width=10, height=8)
# 
# ggplot(data=subset(Vs_melt, (period<12000)&(type =='bvc'))) +
#   geom_point(aes(x=cv_mean, y= bv_mean)) +
#   geom_linerange(aes(x=cv_mean, ymin=bv_lo, ymax=bv_hi)) +
#   geom_smooth(aes(x=cv_mean, y= bv_mean), method=lm, fullrange=TRUE) +
#   facet_wrap(~taxon, scales = 'free_y') +
#   theme_bw(14) +
#   xlab('Climate velocity (km/decade)') +
#   ylab('Biotic velocity (km/decade)')
# ggsave('figures/BVC_vs_CV_mean_scatter_credible_taxon_holocene.pdf', width=10, height=8)
# ggsave('figures/BVC_vs_CV_mean_scatter_credible_taxon_holocene.png', width=10, height=8)
# 
# ggplot(data=subset(Vs_melt, (period<12000)&(type =='bvc'))) +
#   geom_point(aes(x=cv_mean, y= bv_mean, colour=factor(period))) +
#   geom_linerange(aes(x=cv_mean, ymin=bv_lo, ymax=bv_hi, colour=factor(period))) +
#   geom_smooth(aes(x=cv_mean, y= bv_mean), method=lm) +
#   facet_wrap(~taxon, scales = 'free_y') +
#   theme_bw(14) +
#   xlab('Climate velocity (km/decade)') +
#   ylab('Biotic velocity (km/decade)')
# ggsave('figures/BVC_vs_CV_mean_scatter_period_credible_taxon_holocene.pdf', width=10, height=8)
# ggsave('figures/BVC_vs_CV_mean_scatter_period_credible_taxon_holocene.png', width=10, height=8)
# 
# ggplot(data=subset(Vs_melt, (period<12000)&(type =='bvc'))) +
#   geom_point(aes(x=cv_mean, y= bv_median)) +
#   geom_linerange(aes(x=cv_mean, ymin=bv_lo, ymax=bv_hi)) +
#   geom_smooth(aes(x=cv_mean, y= bv_median), method=lm) +
#   facet_wrap(~taxon, scales = 'free_y') +
#   theme_bw(14) +
#   xlab('Climate velocity (km/decade)') +
#   ylab('Biotic velocity (km/decade)')
# ggsave('figures/BVC_vs_CV_median_scatter_credible_taxon_holocene.pdf', width=10, height=8)
# ggsave('figures/BVC_vs_CV_median_scatter_credible_taxon_holocene.png', width=10, height=8)
# 
# ggplot(data=subset(Vs_melt, (period>14000)&(type =='bvc'))) +
#   geom_point(aes(x=cv_mean, y= bv_median)) +
#   geom_linerange(aes(x=cv_mean, ymin=bv_lo, ymax=bv_hi)) +
#   geom_smooth(aes(x=cv_mean, y= bv_median), method=lm) +
#   facet_wrap(~taxon, scales = 'free_y') +
#   theme_bw(14) +
#   xlab('Climate velocity (km/decade)') +
#   ylab('Biotic velocity (km/decade)')
# ggsave('figures/BVC_vs_CV_median_scatter_credible_taxon_postLGM.pdf', width=10, height=8)
# ggsave('figures/BVC_vs_CV_median_scatter_credible_taxon_postLGM.png', width=10, height=8)
# 
# ggplot(data=subset(Vs_melt, (type =='bvc'))) +
#   geom_point(aes(x=cv_mean, y= bv_median, colour=split_name)) +
#   geom_linerange(aes(x=cv_mean, ymin=bv_lo, ymax=bv_hi, colour=split_name)) +
#   geom_smooth(aes(x=cv_mean, y= bv_median, colour=split_name), method=lm, fullrange=TRUE) +
#   facet_wrap(~taxon, scales = 'free_y') +
#   theme_bw(14) +
#   xlab('Climate velocity (km/decade)') +
#   ylab('Biotic velocity (km/decade)')
# ggsave('figures/BVC_vs_CV_median_scatter_credible_taxon_split.pdf', width=10, height=8)
# ggsave('figures/BVC_vs_CV_median_scatter_credible_taxon_split.png', width=10, height=8)
# 
# ggplot(data=subset(Vs_melt, (type =='bvn'))) +
#   geom_point(aes(x=cv_mean, y= bv_median, colour=split_name)) +
#   geom_linerange(aes(x=cv_mean, ymin=bv_lo, ymax=bv_hi, colour=split_name)) +
#   geom_smooth(aes(x=cv_mean, y= bv_median, colour=split_name), method=lm, fullrange=TRUE) +
#   facet_wrap(~taxon, scales = 'free_y') +
#   theme_bw(14) +
#   xlab('Climate velocity (km/decade)') +
#   ylab('Biotic velocity (km/decade)')
# ggsave('figures/BVN_vs_CV_median_scatter_credible_taxon_split.pdf', width=10, height=8)
# ggsave('figures/BVN_vs_CV_median_scatter_credible_taxon_split.png', width=10, height=8)
# 
# ggplot(data=subset(Vs_melt, (type =='bvs'))) +
#   geom_point(aes(x=cv_mean, y= bv_median, colour=split_name)) +
#   geom_linerange(aes(x=cv_mean, ymin=bv_lo, ymax=bv_hi, colour=split_name)) +
#   geom_smooth(aes(x=cv_mean, y= bv_median, colour=split_name), method=lm, fullrange=TRUE) +
#   facet_wrap(~taxon, scales = 'free_y') +
#   theme_bw(14) +
#   xlab('Climate velocity (km/decade)') +
#   ylab('Biotic velocity (km/decade)')
# ggsave('figures/BVS_vs_CV_median_scatter_credible_taxon_split.pdf', width=10, height=8)
# ggsave('figures/BVS_vs_CV_median_scatter_credible_taxon_split.png', width=10, height=8)
# 
# # ggplot(data=subset(Vs_melt, period<12000)) +
# #   geom_point(aes(x=cv_mean, y= bv_mean/100, colour = type)) +
# #   geom_linerange(aes(x=cv_mean, ymin=bv_lo/100, ymax=bv_hi/100, colour = type)) +
# #   geom_smooth(aes(x=cv_mean, y= bv_mean/100, colour=type, group=type), method=lm) +
# #   facet_wrap(~taxon, scales = 'free_y') +
# #   theme_bw(14) +
# #   xlab('Climate velocity (km/decade)') +
# #   ylab('Biotic velocity (km/decade)')
# # ggsave('figures/BVC_vs_CV_mean_scatter_credible_taxon_holocene.pdf', width=10, height=8)
# # ggsave('figures/BVC_vs_CV_mean_scatter_credible_taxon_holocene.png', width=10, height=8)
# 
# 
# ggplot(subset(Vs_melt,(type =='bvc'))) +
#   geom_point(aes(x=cv_mean, y= bv_mean)) +
#   geom_linerange(aes(x=cv_mean, ymin = bv_lo, ymax=bv_hi)) +
#   geom_smooth(aes(x=cv_mean, y= bv_mean), method=lm) +
#   theme_bw(14) +
#   xlab('Climate velocity (km/decade)') +
#   ylab('Biotic velocity (km/decade)')
# ggsave('figures/BVC_vs_CV_mean_scatter.pdf', width=10, height=8)
# 
# # ggplot() +
# #   geom_point(data=Vs, aes(x=cv_median, y= bv_median/100)) +
# #   geom_smooth(data=Vs, aes(x=cv_median, y= bv_median/100), method=lm) #+
# #   # facet_wrap(~taxon)
# 
# # ggplot(data=Vs_melt) +
# #   geom_point(aes(x=cv_median, y= bv_median)) +
# #   facet_wrap(~period)
# 
# BV_CV_lm = Vs_melt %>%
#   nest_by(taxon, type) %>%
#   mutate(mod = list(lm(bv_mean ~ cv_mean, data = data))) %>%
#   reframe(tidy(mod)) 
# colnames(BV_CV_lm)[which(colnames(BV_CV_lm) == 'estimate')] = 'All_21k'
# colnames(BV_CV_lm)[which(colnames(BV_CV_lm) == 'std.error')] = 'A.std.error'
# colnames(BV_CV_lm)[which(colnames(BV_CV_lm) == 'statistic')] = 'A.statistic'
# colnames(BV_CV_lm)[which(colnames(BV_CV_lm) == 'p.value')] = 'A.p.value'
# 
# BV_CV_lm_post = subset(Vs_melt, period>12000) %>%
#   nest_by(taxon, type) %>%
#   mutate(mod = list(lm(bv_mean ~ cv_mean, data = data))) %>%
#   reframe(tidy(mod)) 
# colnames(BV_CV_lm_post)[which(colnames(BV_CV_lm_post) == 'estimate')] = 'PostLGM'
# colnames(BV_CV_lm_post)[which(colnames(BV_CV_lm_post) == 'std.error')] = 'P.std.error'
# colnames(BV_CV_lm_post)[which(colnames(BV_CV_lm_post) == 'statistic')] = 'P.statistic'
# colnames(BV_CV_lm_post)[which(colnames(BV_CV_lm_post) == 'p.value')] = 'P.p.value'
# 
# BV_CV_lm_hol = subset(Vs_melt, period<12000) %>%
#   nest_by(taxon, type) %>%
#   mutate(mod = list(lm(bv_mean ~ cv_mean, data = data))) %>%
#   reframe(tidy(mod))
# colnames(BV_CV_lm_hol)[which(colnames(BV_CV_lm_hol) == 'estimate')] = 'Holocene'
# colnames(BV_CV_lm_hol)[which(colnames(BV_CV_lm_hol) == 'std.error')] = 'H.std.error'
# colnames(BV_CV_lm_hol)[which(colnames(BV_CV_lm_hol) == 'statistic')] = 'H.statistic'
# colnames(BV_CV_lm_hol)[which(colnames(BV_CV_lm_hol) == 'p.value')] = 'H.p.value'
# 
# 
# BV_CV_slopes = merge(BV_CV_lm, BV_CV_lm_hol, by =c('taxon', 'type', 'term'))
# BV_CV_slopes = merge(BV_CV_slopes, BV_CV_lm_post, by =c('taxon', 'type', 'term'))
# BV_CV_slopes = subset(BV_CV_slopes, term == 'cv_mean')
# BV_CV_slopes$diff_all_H = BV_CV_slopes$All_21k - BV_CV_slopes$Holocene
# BV_CV_slopes$diff_LGM_H = BV_CV_slopes$PostLGM - BV_CV_slopes$Holocene
# 
# 
# ggplot(data=BV_CV_slopes) +
#   geom_histogram(aes(x = diff_all_H)) +#, breaks = seq(-60, 200, by=15)) + 
#   facet_grid(type~.)
# 
# # BV_CV_slopes = rbind(data.frame(subset(BV_CV_lm, term == 'cv_mean'), period = '21k'),
# #                      data.frame(subset(BV_CV_lm_hol, term == 'cv_mean'), period = 'holocene'))
# 
# BV_CV_slopes$A.sig = 0
# BV_CV_slopes$A.sig[which(BV_CV_slopes$A.p.value<0.05)] = 1
# 
# BV_CV_slopes$H.sig = 0
# BV_CV_slopes$H.sig[which(BV_CV_slopes$H.p.value<0.05)] = 1
# 
# 
# ggplot(data=BV_CV_slopes) +
#   geom_point(aes(x = All_21k, 
#                  y = taxon, 
#                  shape = factor(A.sig), 
#                  colour=type), 
#              size=2)
# 
# ggplot(data=BV_CV_slopes) +
#   geom_histogram(aes(x = All_21k, 
#                      colour = type,
#                      fill = type), 
#                  alpha=0.5)
# 
# ggplot(data=BV_CV_slopes) +
#   geom_histogram(aes(x = All_21k, y=after_stat(density)), bins = 15, breaks = seq(-0.6, 2, by = 0.1),
#                  alpha=0.5, colour='blue', fill='lightblue') +
#   geom_vline(xintercept = 0, colour = 'black', lty = 2) +
#   facet_grid(A.sig~type) +
#   theme_bw() +
#   theme(text = element_text(size=14))
# 
# ggplot(data=BV_CV_slopes) +
#   geom_histogram(aes(x = All_21k, y=after_stat(density)), bins = 15, breaks = seq(-0.6, 2, by = 0.1),
#                  alpha=0.5, colour='blue', fill='lightblue') +
#   geom_vline(xintercept = 0, colour = 'black', lty = 2) +
#   facet_grid(type~.) +
#   theme_bw() +
#   theme(text = element_text(size=14))
# ggsave('figures/BV_vs_CV_slopes_histogram_LGM.pdf', width=10, height=8)
# ggsave('figures/BV_vs_CV_slopes_histogram_LGM.png', width=10, height=8)
# 
# ggplot(data=BV_CV_slopes) +
#   geom_histogram(aes(x = PostLGM, y=after_stat(density)), bins = 15, breaks = seq(-0.6, 2, by = 0.1),
#                  alpha=0.5, colour='blue', fill='lightblue') +
#   geom_vline(xintercept = 0, colour = 'black', lty = 2) +
#   facet_grid(type~.) +
#   theme_bw() +
#   theme(text = element_text(size=14))
# ggsave('figures/BV_vs_CV_slopes_histogram_postLGM.pdf', width=10, height=8)
# ggsave('figures/BV_vs_CV_slopes_histogram_postLGM.png', width=10, height=8)
# 
# ggplot(data=BV_CV_slopes) +
#   geom_histogram(aes(x = Holocene, y=after_stat(density)), bins = 15, breaks = seq(-0.6, 2, by = 0.1),
#                  alpha=0.5, colour='blue', fill='lightblue') +
#   geom_vline(xintercept = 0, colour = 'black', lty = 2) +
#   facet_grid(H.sig~type) +
#   theme_bw() +
#   theme(text = element_text(size=14))
# 
# ggplot(data=BV_CV_slopes) +
#   geom_histogram(aes(x = Holocene, y=after_stat(density)), bins = 15, breaks = seq(-0.6, 2, by = 0.1),
#                  alpha=0.5, colour='blue', fill='lightblue') +
#   geom_vline(xintercept = 0, colour = 'black', lty = 2) +
#   facet_grid(type~.) +
#   theme_bw() +
#   theme(text = element_text(size=14))
# ggsave('figures/BV_vs_CV_slopes_histogram_holocene.pdf', width=10, height=8)
# ggsave('figures/BV_vs_CV_slopes_histogram_holocene.png', width=10, height=8)
# 
# ggplot(data=BV_CV_slopes) +
#   geom_histogram(aes(x = diff_all_H, y=after_stat(density)), bins = 15, breaks = seq(-0.50, 2, by = 0.1),
#                  alpha=0.5, colour='blue', fill='lightblue') +
#   geom_vline(xintercept = 0, colour = 'black', lty = 2) +
#   facet_grid(type~.) +
#   theme_bw() +
#   theme(text = element_text(size=14))
# ggsave('figures/BV_vs_CV_slopes_histogram_holocene_LGM_diff.pdf', width=10, height=8)
# ggsave('figures/BV_vs_CV_slopes_histogram_holocene_LGM_diff.png', width=10, height=8)
# 
# ggplot(data=BV_CV_slopes) +
#   geom_histogram(aes(x = diff_LGM_H, y=after_stat(density)), bins = 15, breaks = seq(-0.50, 2, by = 0.1),
#                  alpha=0.5, colour='blue', fill='lightblue') +
#   geom_vline(xintercept = 0, colour = 'black', lty = 2) +
#   facet_grid(type~.) +
#   theme_bw() +
#   theme(text = element_text(size=14))
# ggsave('figures/BV_vs_CV_slopes_histogram_holocene_LGM_diff.pdf', width=10, height=8)
# ggsave('figures/BV_vs_CV_slopes_histogram_holocene_LGM_diff.png', width=10, height=8)
# 
# # ggplot(data=slopes) +
# #   geom_density(aes(x = estimate), 
# #                  alpha=0.5) +
# #   facet_grid(type~.)
# 
# 
# ggplot(data=BV_CV_slopes) +
#   geom_point(aes(x = PostLGM,
#                  y = Holocene,
#                  colour=type),
#              size=2) +
#   geom_smooth(aes(x=PostLGM, y= Holocene, colour=type), method=lm) +
#   # facet_wrap(~taxon, scales = 'free_y') +
#   theme_bw(14) +
#   xlab('Post LGM slopes') +
#   ylab('Holocene slopes') +
#   coord_fixed()
# 
# 
# BV_CV_slopes[which(BV_CV_slopes$All_21k < 0),]
# 
# 
# BV_CV_slopes_long = BV_CV_slopes[,c('taxon', 'type', 'All_21k', 'PostLGM', 'Holocene')] %>% 
#   pivot_longer(cols = c('All_21k', 'PostLGM', 'Holocene'))
# 
# ggplot(data=BV_CV_slopes_long) +
#   geom_histogram(aes(x=value, y=after_stat(density))) +
#   facet_grid(name~type)
# 
# ggplot(data=BV_CV_slopes_long) +
#   geom_density(aes(x=value, y=after_stat(density)), colour='dodgerblue', fill='dodgerblue') +
#   facet_grid(name~type)
# 
# ggplot(data=subset(BV_CV_slopes_long, name %in% c('Holocene', 'PostLGM'))) +
#   geom_density(aes(x=value, y=after_stat(density), colour=name, fill=name), alpha=0.4) +
#   facet_grid(type~.) +
#   theme_bw(14) +
#   xlab('Regression slopes BV vs CV')
# ggsave('figures/BV_vs_CV_slopes_density_split.pdf', width=10, height=8)
# ggsave('figures/BV_vs_CV_slopes_density_split.png', width=10, height=8)
# 
# ggplot(data=BV_CV_slopes_long) +
#   geom_density(aes(x=value, y=after_stat(density), colour=name, fill=name), alpha=0.4) +
#   facet_grid(type~.) +
#   theme_bw(14) +
#   xlab('Regression slopes BV vs CV')
# ggsave('figures/BV_vs_CV_slopes_density_all_split.pdf', width=10, height=8)
# ggsave('figures/BV_vs_CV_slopes_density_all_split.png', width=10, height=8)

######################################################################################################################
## climate velocity
######################################################################################################################

CVs = readRDS('output/climate_velocity.RDS')

CVs_summary = CVs %>% 
  group_by(period) %>% 
  dplyr::summarize(cv_mean = mean(voccMag, na.rm=TRUE),
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
  geom_point(aes(x = Holocene, 
                 y = taxon, 
                 shape = factor(sig), 
                 colour=type), 
             size=2)

ggplot(data=BV_CV_slopes) +
  geom_histogram(aes(x = Holocene, 
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
ggsave('figures/BV_vs_CV_slopes_histogram_LGM.pdf', width=10, height=8)
ggsave('figures/BV_vs_CV_slopes_histogram_LGM.png', width=10, height=8)

ggplot(data=BV_CV_slopes) +
  geom_histogram(aes(x = Holocene), bins = 15, breaks = seq(-60, 200, by = 15),
                 alpha=0.5, colour='blue', fill='lightblue', alpha=0.5) +
  geom_vline(xintercept = 0, colour = 'black', lty = 2) +
  facet_grid(type~.) +
  theme_bw() +
  theme(text = element_text(size=14))
ggsave('figures/BV_vs_CV_slopes_histogram_holocene.pdf', width=10, height=8)
ggsave('figures/BV_vs_CV_slopes_histogram_holocene.png', width=10, height=8)

ggplot(data=BV_CV_slopes) +
  geom_histogram(aes(x = diff), bins = 15, breaks = seq(-50, 160, by = 15),
                 alpha=0.5, colour='blue', fill='lightblue', alpha=0.5) +
  geom_vline(xintercept = 0, colour = 'black', lty = 2) +
  facet_grid(type~.) +
  theme_bw() +
  theme(text = element_text(size=14))
ggsave('figures/BV_vs_CV_slopes_histogram_holocene_LGM_diff.pdf', width=10, height=8)
ggsave('figures/BV_vs_CV_slopes_histogram_holocene_LGM_diff.png', width=10, height=8)

# ggplot(data=slopes) +
#   geom_density(aes(x = estimate), 
#                  alpha=0.5) +
#   facet_grid(type~.)


######################################################################################################################
## NS quant 95 BV
######################################################################################################################

