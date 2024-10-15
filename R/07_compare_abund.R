
## FishN: comparing fish assemblage abundance surveying methods
# Abundance estimate comparisons

require(dplyr)
require(lubridate)
require(beepr)

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
model_abund[[1]] <- lm(log(MaxN + 1) ~ log(REST + 1), data = dt_abund)
# UVC ~ REST
model_abund[[2]] <- lm(log(UVC + 1) ~ log(REST + 1), data = dt_abund)
# MaxN ~ UVC
model_abund[[3]] <- lm(log(MaxN + 1) ~ log(UVC + 1), data = dt_abund)

# candidates with aggregation
model_abund[[4]] <- lm(log(MaxN + 1) ~ log(REST + 1) + Group, data = dt_abund)
# UVC ~ REST
model_abund[[5]] <- lm(log(UVC + 1) ~ log(REST + 1) + Group, data = dt_abund)
# MaxN ~ UVC
model_abund[[6]] <- lm(log(MaxN + 1) ~ log(UVC + 1) + Group, data = dt_abund)

summary(model_abund[[4]])
summary(model_abund[[5]])
summary(model_abund[[6]])

sapply(model_abund, AIC)[c(1,4,2,5,3,6)] # look at AICs
# aggregation is consistently lower

# model prediction
# make a standard new x that will be used by all models, with aggregation held at a mean
social <- mean(dt_abund$Group)
newx <- data.frame(REST = seq(min(dt_abund$REST), max(dt_abund$REST), length.out = 50), 
                   UVC = seq(min(dt_abund$UVC), max(dt_abund$UVC), length.out = 50), 
                   Group = social)
abund_pred <- as.list(rep('', 3))
for (i in 1:3) { # each model
  abund_pred[[i]] <- predict(model_abund[[i+3]], newdata = newx, se.fit = T)
  abund_pred[[i]] <- with(abund_pred[[i]], data.frame(fitted = fit, lwr = fit - 1.96 * se.fit, 
                                                  upr = fit + 1.96 * se.fit))
}
# add newx to the prediction data frame
abund_pred[[1]] <- abund_pred[[1]] %>% bind_cols(newx %>% select(REST, Group))
abund_pred[[2]] <- abund_pred[[2]] %>% bind_cols(newx %>% select(REST, Group))
abund_pred[[3]] <- abund_pred[[3]] %>% bind_cols(newx %>% select(UVC, Group))

# partial regression plots for observed abundances (no aggregation yet)
f_abundpair <- as.list(rep('', 3))
f_abundpair[[1]] <- ggplot(dt_abund) +
  geom_jitter(aes(y = log(MaxN + 1), x = log(REST + 1), fill = Size_class), width = 0.05, height = 0.05, 
              alpha = 0.6, size = 2, shape = 21) +
  geom_ribbon(data = abund_pred[[1]], 
              aes(x = log(REST + 1), ymin = lwr, ymax = upr), linetype = 'dashed', color = 'black', linewidth= 0.5, fill = 'transparent') +
  geom_line(data = abund_pred[[1]], aes(x = log(REST + 1), y = fitted)) +
  labs(x = 'log(REST)', y = 'log(MaxN)') +
  coord_cartesian(ylim = c(0,4.75)) +
  scale_fill_cherulean(palette = "cheridis", discrete = T, name = "Size class (cm)")

# REST - UVC
f_abundpair[[2]] <- ggplot(dt_allmethods) +
  geom_jitter(aes(y = log(UVC + 1), x = log(REST + 1), fill = Size_class), width = 0.05, height = 0.05, 
              alpha = 0.6, size = 2, shape =21) +
  geom_ribbon(data = abund_pred[[2]], aes(x = log(REST + 1), ymin = lwr, ymax = upr), 
              linetype = 'dashed', color = 'black', linewidth= 0.5, fill = 'transparent') +
  geom_line(data = abund_pred[[2]], aes(x = log(REST + 1), y = fitted)) +
  labs(x = 'log(REST)', y = 'log(UVC)') +
  coord_cartesian(ylim = c(0,4.75)) +
  scale_fill_cherulean(palette = "cheridis", discrete = T, name = "Size class (cm)")

# UVC - MaxN
f_abundpair[[3]] <- ggplot(dt_allmethods) +
  geom_jitter(aes(y = log(MaxN + 1), x = log(UVC + 1), fill = Size_class), width = 0.05, height = 0.05, 
              alpha = 0.6, size = 2, shape = 21) +
  geom_ribbon(data = abund_pred[[3]], 
              aes(x = log(UVC + 1), ymin = lwr, ymax = upr), linetype = 'dashed', color = 'black', linewidth= 0.5, fill = 'transparent') +
  geom_line(data = abund_pred[[3]], aes(x = log(UVC + 1), y = fitted)) +
  labs(x = 'log(UVC)', y = 'log(MaxN)')+
  coord_cartesian(ylim = c(0,4.75)) +
  scale_fill_cherulean(palette = "cheridis", discrete = T, name = "Size class (cm)")

# partial regression plots for partial aggregation effects
# generate partial predictions for aggregation, held at abundance mean
meanREST <- mean(dt_abund$REST)
meanPoint <- mean(dt_abund$UVC)
newx <- data.frame(Group = 1:5, REST = meanREST, UVC = meanPoint)
for (i in 4:6) { # each model
  abund_pred[[i]] <- predict(model_abund[[i]], newdata = newx, se.fit = T)
  abund_pred[[i]] <- with(abund_pred[[i]], data.frame(fitted = fit, lwr = fit - 1.96 * se.fit, 
                                                  upr = fit + 1.96 * se.fit,
                                                  Group = newx$Group))
}

f_abundsoc <- as.list(rep('', 3))

f_abundsoc[[1]] <- ggplot(dt_abund) +
  geom_point(aes(x = Group, y = log(MaxN + 1)), size = 1, alpha = .4, shape = 21, color = 'grey30', 
             fill = 'white', position = position_jitter(width = .1)) +
  geom_errorbar(data = abund_pred[[4]], aes(x = Group + 0.25, ymin = lwr, ymax = upr), 
                color = 'grey20', width = 0.2) +
  geom_point(data = abund_pred[[4]], aes(x = Group + 0.25, y = fitted), 
             color = 'grey20', size = 2.5) +
  labs(x = 'Aggregation', y = 'log(MaxN)')

f_abundsoc[[2]] <- ggplot(dt_abund) +
  geom_point(aes(x = Group, y = log(UVC + 1)), size = 1, alpha = .4, shape = 21, color = 'grey20', 
             fill = 'white', position = position_jitter(width = .1)) +
  geom_errorbar(data = abund_pred[[5]], aes(x = Group + 0.25, ymin = lwr, ymax = upr), 
                color = 'grey20', width = 0.2) +
  geom_point(data = abund_pred[[5]], aes(x = Group + 0.25, y = fitted), 
             color = 'grey20', size = 2.5) +
  labs(x = 'Aggregation', y = 'log(UVC)')

f_abundsoc[[3]] <- ggplot(dt_abund) +
  geom_point(aes(x = Group, y = log(MaxN + 1)), size = 1, alpha = .4, shape = 21, color = 'grey20', 
             fill = 'white', position = position_jitter(width = .1)) +
  geom_errorbar(data = abund_pred[[6]], aes(x = Group + 0.25, ymin = lwr, ymax = upr), 
                color = 'grey20', width = 0.2) +
  geom_point(data = abund_pred[[6]], aes(x = Group + 0.25, y = fitted), 
             color = 'grey20', size = 2.5) +
  labs(x = 'Aggregation', y = 'log(MaxN)')


## Model effects figure ----------------------------------------------------

pair_effects <- data.frame(model = rep(c('MaxN ~ REST', 'UVC ~ REST', 'MaxN ~ UVC'), each = 2), 
                           par = rep(c('method', 'aggregation'), 3), 
                           estimate = sapply(model_abund[c(4:6)], FUN = function (x) {coefficients(x)[2:3]}, simplify = T) %>% as.vector, 
                           se = sapply(model_abund[c(4:6)], FUN = function (x) {summary(x)$coefficients[-1,2]}, simplify = T) %>% as.vector)
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
ggsave('../../MEE/resubmission/figures/fig_abundmethod.pdf', width = 18, height = 6.5, units = "cm")
(f_abundsoc[[1]] + f_abundsoc[[2]] + f_abundsoc[[3]]) * looks * scale_x_continuous(breaks = c(1:5)) * coord_cartesian(xlim = c(0.5,5.7), ylim = c(0,4.75), clip = 'off') + plot_layout(guides = "collect")
ggsave('../../MEE/resubmission/figures/fig_aggregation.pdf', width = 18, height = 6.5, units = "cm")
f_effects + theme(legend.position = 'bottom')
ggsave('../../MEE/resubmission/figures/fig_effects.pdf', width = 8, height = 8, units = "cm")
f_size
ggsave(plot = f_size, '../../MEE/resubmission/figures/fig_size.pdf', width = 18, height = 6.5, units = "cm")
rm(traits, totals, temp, newx)
