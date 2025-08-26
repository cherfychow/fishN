
## FishN: comparing fish assemblage abundance surveying methods
# Abundance estimate comparisons

require(dplyr)
require(beepr)
require(glmmTMB)

require(stringr)
require(ggplot2)
require(patchwork)
source('https://gist.githubusercontent.com/cherfychow/e9ae890fd16f4c86730748c067feee2b/raw/b2db138ab8164c237129308ea643de78cd54b252/cherulean.R')

looks <- theme_classic(base_size = 13) + 
  theme(panel.grid = element_blank(), 
  axis.line = element_line(linewidth = .35), 
  axis.ticks = element_line(linewidth = .35))


# PAIRWISE ABUNDANCE COMPARISONS --------------------------------------------
# look at scatterplots first
# 
# # REST - MaxN
# ggplot(dt_allmethods) +
#   geom_jitter(aes(x = REST, y = MaxN, fill = Size_class), width = 0.05, height = 0.05,
#               shape = 21, alpha = 0.6, size = 2) +
#   looks + scale_fill_cherulean(palette = "cheridis", discrete = T, reverse = T)
# 
# # REST - UVC
# ggplot(dt_allmethods) +
#   geom_jitter(aes(x = REST, y = UVC, fill = Size_class), width = 0.05, height = 0.05,
#               shape = 21, alpha = 0.6, size = 2) +
#   looks + scale_fill_cherulean(palette = "cheridis", discrete = T, reverse = T)
# 
# # UVC - MaxN
# ggplot(dt_allmethods) +
#   geom_jitter(aes(x = UVC, y = MaxN, fill = Size_class), width = 0.05, height = 0.05,
#               shape = 21, alpha = 0.6, size = 2) +
#   looks + scale_fill_cherulean(palette = "cheridis", discrete = T, reverse = T)

## Fit GLMs --------------------------------------------

set.seed(240)
# initially started out thinking about fitting Size_class as an effect, 
# but the size spectra makes me think that letting size absorb variation is not the way to go
# dominated by 10-19 cm

# read in traits to add in
traits <- read.csv('../data/traits_group.csv', header = T)
dt_abund <- left_join(dt_allmethods, traits, by = "Taxon")


model_abund <- as.list(rep('', 6))
# MaxN ~ REST
model_abund[[1]] <- glmmTMB(REST ~ log(MaxN + 1) + (1|site_ID), data = dt_abund, 
  family = ziGamma(link = 'log'), ziformula = ~log(MaxN + 1))
# UVC ~ REST
model_abund[[2]] <- glmmTMB(UVC ~ log(REST + 1) + (1|site_ID), data = dt_abund,
  family = truncated_poisson(link = 'log'), ziformula = ~log(REST + 1))
# MaxN ~ UVC
model_abund[[3]] <- glmmTMB(MaxN ~ log(UVC + 1) + (1|site_ID), data = dt_abund,
  family = truncated_poisson(link = 'log'), ziformula = ~log(UVC + 1))

# candidates with aggregation
model_abund[[4]] <- glmmTMB(REST ~ log(MaxN + 1) + Group + (1|site_ID), data = dt_abund, 
  family = ziGamma(link = 'log'), ziformula = ~log(MaxN + 1))
# UVC ~ REST
model_abund[[5]] <- glmmTMB(UVC ~ log(REST + 1) + Group + (1|site_ID), data = dt_abund,
  family = truncated_poisson(link = 'log'), ziformula = ~log(REST + 1))
# MaxN ~ UVC
model_abund[[6]] <- glmmTMB(MaxN ~ log(UVC + 1) + Group+ (1|site_ID), data = dt_abund,
  family = truncated_poisson(link = 'log'), ziformula = ~log(UVC + 1))

require(performance)

sapply(model_abund, AIC)[c(1,4,2,5,3,6)] # look at AICs
# aggregation is consistently lower
# 1, 5, 6 

# model prediction
# make a standard new x that will be used by all models, with aggregation held at a mean
social <- mean(dt_abund$Group)
newx <- data.frame(REST = seq(min(dt_abund$REST), max(dt_abund$REST), length.out = 50), 
                   UVC = seq(min(dt_abund$UVC), max(dt_abund$UVC), length.out = 50), 
                   MaxN = seq(min(dt_abund$MaxN), max(dt_abund$MaxN), length.out = 50), 
                   Group = social)
abund_pred <- as.list(rep('', 3))
selected = c(1,5,6)
for (i in 1:3) { # each model
  abund_pred[[i]] <- predict(model_abund[[selected[i]]], newdata = newx, se.fit = T, 
    type = 'response', re.form = NA)
  abund_pred[[i]] <- with(abund_pred[[i]], data.frame(fitted = fit, lwr = fit - 1.96 * se.fit, 
                                                  upr = fit + 1.96 * se.fit))
}
# add newx to the prediction data frame
abund_pred[[1]] <- abund_pred[[1]] %>% bind_cols(newx %>% select(MaxN))
abund_pred[[2]] <- abund_pred[[2]] %>% bind_cols(newx %>% select(REST, Group))
abund_pred[[3]] <- abund_pred[[3]] %>% bind_cols(newx %>% select(UVC, Group))

# partial regression plots for observed abundances (no aggregation yet)
f_abundpair <- as.list(rep('', 3))
f_abundpair[[1]] <- ggplot(dt_abund) +
  geom_point(aes(x = MaxN , y = REST , fill = Size_class), 
    alpha = 0.6, size = 2, shape = 21) +
  geom_ribbon(data = abund_pred[[1]], 
              aes(x = MaxN , ymin = lwr, ymax = upr), linetype = 'dashed', 
              color = 'black', linewidth= 0.5, fill = 'transparent') +
  geom_line(data = abund_pred[[1]], aes(x = MaxN , y = fitted)) +
  labs(y = 'REST', x = 'MaxN') +
  scale_fill_cherulean(palette = "cheridis", discrete = T, name = "Size class (cm)") +
  scale_x_log10() + scale_y_log10()

# REST - UVC
f_abundpair[[2]] <- ggplot(dt_abund) +
  geom_point(aes(y = UVC , x = REST, fill = Size_class), 
              alpha = 0.6, size = 2, shape =21) +
  geom_ribbon(data = abund_pred[[2]], aes(x = REST, ymin = lwr, ymax = upr), 
              linetype = 'dashed', color = 'black', linewidth= 0.5, fill = 'transparent') +
  geom_line(data = abund_pred[[2]], aes(x = REST, y = fitted)) +
  labs(x = 'REST', y = 'UVC') +
  scale_fill_cherulean(palette = "cheridis", discrete = T, name = "Size class (cm)")+
  scale_x_log10() + scale_y_log10()

# UVC - MaxN
f_abundpair[[3]] <- ggplot(dt_abund) +
  geom_point(aes(y = MaxN , x = UVC , fill = Size_class),  
              alpha = 0.6, size = 2, shape = 21) +
  geom_ribbon(data = abund_pred[[3]], 
              aes(x = UVC, ymin = lwr, ymax = upr), linetype = 'dashed', color = 'black', linewidth= 0.5, fill = 'transparent') +
  geom_line(data = abund_pred[[3]], aes(x = UVC, y = fitted)) +
  labs(x = 'UVC', y = 'MaxN')+
  scale_fill_cherulean(palette = "cheridis", discrete = T, name = "Size class (cm)")+
  scale_x_log10() + scale_y_log10()

# partial regression plots for partial aggregation effects
# generate partial predictions for aggregation, held at abundance mean
meanREST <- mean(dt_abund$REST)
meanPoint <- mean(dt_abund$UVC)
newx <- data.frame(Group = 1:5, REST = meanREST, UVC = meanPoint)
for (i in 2:3) { # each model
  abund_pred[[i+2]] <- predict(model_abund[[selected[i]]], newdata = newx, se.fit = T, re.form = NA, type = 'response')
  abund_pred[[i+2]] <- with(abund_pred[[i+2]], data.frame(fitted = fit, lwr = fit - 1.96 * se.fit, 
                                                  upr = fit + 1.96 * se.fit,
                                                  Group = newx$Group))
}

f_abundsoc <- as.list(rep('', 2))

f_abundsoc[[1]] <- ggplot(dt_abund) +
  geom_point(aes(x = Group, y = UVC), size = 1, alpha = .4, shape = 21, color = 'grey20', 
             fill = 'white', position = position_jitter(width = .1)) +
  geom_errorbar(data = abund_pred[[4]], aes(x = Group + 0.25, ymin = lwr, ymax = upr), 
                color = 'grey20', width = 0.2) +
  geom_point(data = abund_pred[[4]], aes(x = Group + 0.25, y = fitted), 
             color = 'grey20', size = 2.5) +
  labs(x = 'Aggregation', y = 'UVC') + scale_y_log10()

f_abundsoc[[2]] <- ggplot(dt_abund) +
  geom_point(aes(x = Group, y = MaxN), size = 1, alpha = .4, shape = 21, color = 'grey20', 
             fill = 'white', position = position_jitter(width = .1)) +
  geom_errorbar(data = abund_pred[[5]], aes(x = Group + 0.25, ymin = lwr, ymax = upr), 
                color = 'grey20', width = 0.2) +
  geom_point(data = abund_pred[[5]], aes(x = Group + 0.25, y = fitted), 
             color = 'grey20', size = 2.5) +
  labs(x = 'Aggregation', y = 'MaxN') + scale_y_log10()


## Model effects figure ----------------------------------------------------

# method effect
pair_effects <- data.frame(model = c('REST ~ MaxN', 'UVC ~ REST', 'MaxN ~ UVC'), 
                           par = 'method', 
                           estimate = sapply(model_abund[selected], FUN = function (x) {summary(x)$coefficients$cond[2,1]}, simplify = T) %>% as.vector, 
                           se = sapply(model_abund[selected], FUN = function (x) {summary(x)$coefficients$cond[2,2]}, simplify = T) %>% as.vector)

# Aggregation effect
pair_effects <- data.frame(model = c('UVC ~ REST', 'MaxN ~ UVC'), 
                           par = 'aggregation', 
                           estimate = sapply(model_abund[c(5,6)], FUN = function (x) {summary(x)$coefficients$cond[3,1]}, simplify = T) %>% as.vector, 
                           se = sapply(model_abund[c(5,6)], FUN = function (x) {summary(x)$coefficients$cond[3,2]}, simplify = T) %>% as.vector) %>%
                bind_rows(pair_effects, .)

# random intercept                          
pair_effects <- data.frame(model = c('REST ~ MaxN', 'UVC ~ REST', 'MaxN ~ UVC'), 
                           par = rep(c('(1|site)'), 3), 
                           estimate = sapply(model_abund[selected], FUN = function (x) {summary(x)$varcor$cond$site_ID %>% attr(., "stddev")}, simplify = T) %>% as.vector) %>%
                bind_rows(pair_effects, .)
pair_effects$lwr <- with(pair_effects, estimate - 1.96 * se)
pair_effects$upr <- with(pair_effects, estimate + 1.96 * se)

f_effects <- ggplot(data = pair_effects) +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey60') +
  geom_errorbar(aes(xmin = lwr, xmax = upr, y = model, color = par), width = 0.2, linewidth = 0.75) +
  geom_point(aes(x = estimate, y = model, fill = par), size = 3, shape = 21, color = 'grey20') +
  scale_color_cherulean(palette = "cheridis", discrete = T, name = NULL) + looks +
  scale_fill_cherulean(palette = "cheridis", discrete = T, name = NULL) +
  scale_y_discrete(position = "right") +
  # scale_x_continuous(lim = c(0,1.2)) +
  labs(y = NULL, x = 'Effect estimate')


## Size spectra comparison ----------------------------------------------------

# couldn't fit size as a factor so here is the comparison

# make the counts relative by site totals first
totals <- dt_all_long %>% group_by(method, site_ID) %>% summarise(total = sum(n))
dt_all_size <- dt_all_long %>% group_by(site_ID, method, Size_class) %>% summarise(n = sum(n))

for (i in 1:nrow(totals)) {
  dt_all_size$n[with(dt_all_size, site_ID == totals$site_ID[i] & method == totals$method[i])] <- dt_all_size$n[with(dt_all_size, site_ID == totals$site_ID[i] & method == totals$method[i])]/totals$total[i]
}

dt_all_size$midSize <- str_extract(dt_all_size$Size_class, '(?<=\\-|(\\<\\s))[:digit:]+$') %>% 
  as.numeric
dt_all_size$midSize[dt_all_size$midSize > 5] <- dt_all_size$midSize[dt_all_size$midSize > 5] -4
dt_all_size$midSize[dt_all_size$midSize == 5] <- dt_all_size$midSize[dt_all_size$midSize == 5] - 2.5

sizespec <- dt_all_size[rep(1:nrow(dt_all_size), dt_all_size$n*100), c(1,2,5)]
sizespec <- sizespec %>% filter(!method == 'UVC')
temp <- dt_uvc %>% filter(site_ID %in% unique(dt_all_size$site_ID)) %>% select(site_ID, TL_cm, count)
sizespec <- temp[rep(1:nrow(temp), temp$count), 1:2] %>% rename(midSize = TL_cm) %>% mutate(method = 'UVC', .before = midSize) %>% bind_rows(sizespec, .)

f_size <- ggplot(sizespec) +
  geom_density(aes(x = midSize, fill = method, color = method), alpha = 0.4) +
  labs(x = "Observed TL (cm)", y = 'Frequency') +
  facet_grid(cols = vars(site_ID)) +
  looks + scale_fill_cherulean(palette = "gomphosus", discrete = T) + 
  scale_color_cherulean(palette = "gomphosus", discrete = T) + 
  theme(strip.background = element_blank(), strip.text = element_blank())


# Print figures -----------------------------------------------------------

(f_abundpair[[1]] + f_abundpair[[2]] + f_abundpair[[3]]) * looks * theme(legend.position = "none")
ggsave('fig_abundmethod.pdf', width = 20, height = 6.5, units = "cm")
(f_abundsoc[[1]] + f_abundsoc[[2]]) * looks * scale_x_continuous(breaks = c(1:5)) + plot_layout(guides = "collect")
ggsave('fig_aggregation.pdf', width = 13, height = 6.5, units = "cm")
f_effects + theme(legend.position = 'bottom')
ggsave('fig_effects.pdf', width = 8, height = 8, units = "cm")
f_size
ggsave(plot = f_size, '../../MEE/resubmission/figures/fig_size.pdf', width = 18, height = 6.5, units = "cm")
rm(traits, totals, temp, newx)
