

## FishN: comparing fish assemblage abundance surveying methods
# Comparison between REST, MaxN, and UVC

require(dplyr)
require(lubridate)
require(beepr)

require(vegan)
require(stringr)
require(ape)
require(ggplot2)
require(patchwork)
source('https://gist.githubusercontent.com/cherfychow/e9ae890fd16f4c86730748c067feee2b/raw/b2db138ab8164c237129308ea643de78cd54b252/cherulean.R')

# load files generated from previous analyses 
# COMMENT OUT BEFORE SOURCING

# data_maxn <- read.csv('../data/data_MaxN.csv', header = T)
# data_uvc <- read.csv('../data/uvc_himb.csv', header = T)
# output_poisson <- read.csv('../outputs/REST_final.csv', header = T)

# PCOA COMPOSITION COMPARISON -----------------------------------------------------------------

# prep REST to match the assemblage survey data like RUV or MaxN
set.seed(24)

data_uvc$Size_class[data_uvc$TL_cm <= 5] <- "_5"
data_uvc$Size_class[data_uvc$TL_cm > 5 & data_uvc$TL_cm <= 10] <- "5_9"
data_uvc$Size_class[data_uvc$TL_cm > 10 & data_uvc$TL_cm <= 20] <- "10_19"
data_uvc$Size_class[data_uvc$TL_cm > 20 & data_uvc$TL_cm <= 30] <- "20_29"
data_uvc$Size_class[data_uvc$TL_cm > 30 & data_uvc$TL_cm <= 40] <- "30_39"
data_uvc$Size_class[data_uvc$TL_cm > 40 & data_uvc$TL_cm <= 50] <- "40_49"
data_uvc$Size_class[data_uvc$TL_cm > 50 & data_uvc$TL_cm <= 60] <- "50_59"

data_rest <- data_ruv2 %>% distinct(site_ID, Taxon, Size_class, spsize) %>% 
  left_join(., output_poisson[c(1,3,4,12)], by=c("site_ID", "Taxon", "Size_class"))

# cross-method comparison doesn't have to deal with size classes
# construct the community matrix
nmatrix <- data_uvc %>% filter(str_detect(site_ID, 'hinalea|kaku|sunset')) %>% 
  group_by(site_ID, Taxon) %>% 
  summarise(n = sum(count)) %>% mutate(method = 'uvc') %>% ungroup # filter and aggregate UVC
nmatrix <- data_rest %>% group_by(site_ID, Taxon) %>% 
  summarise(n = sum(D, na.rm = T)) %>% mutate(method = 'rest') %>% ungroup %>% filter(!is.na(n), !n == 0) %>% 
  bind_rows(nmatrix, .)
nmatrix <- data_maxn %>% group_by(site_ID, Taxon) %>% summarise(n = sum(MaxN, na.rm = T)) %>% ungroup %>% 
  mutate(method = 'maxN') %>% bind_rows(nmatrix, .)
nmatrix <- nmatrix %>% tidyr::pivot_wider(., names_from = Taxon, values_from = n)
nmatrix <- nmatrix %>% replace(is.na(.), 0) # replace NAs with 0

dis <- vegdist(nmatrix[-(1:2)], method = "euclidean", diag = F, binary = F)
method_pcoa <- ape::pcoa(dis, correction = "none", rn = paste(nmatrix$site_ID, nmatrix$method, sep="_"))
biplot(method_pcoa)

## Plot PCoA with symbology -----------------------------------------------------------------
source('https://raw.githubusercontent.com/cherfychow/FishTraitsCoralRec/main/analysis_code/function_convhull.R')
# convex hull script

# scree plot
ggplot(data=method_pcoa$values[1:7,], aes(x=1:7, y=Relative_eig/sum(method_pcoa$values$Relative_eig))) +
  geom_line() + geom_point(shape=21, fill='white', size=3) + labs(x=NULL, y=NULL) +
  scale_x_continuous(breaks=c(1:10)) +
  theme_bw(base_size = 13) + labs(title = "Scree plot (PCoA)")

pcoadt <- as.data.frame(method_pcoa$vectors) %>% bind_cols(., nmatrix[1:2])

f_pcoa1 <- ggplot(data = pcoadt) +
  geom_polygon(aes(x = Axis.1, y = Axis.2, fill = site_ID, color = site_ID), alpha = 0.5) +
  geom_point(aes(x = Axis.1, y = Axis.2, fill = site_ID, shape = method), size = 3) +
  looks + labs(x = "PCo1", y = "PCo2") +
  scale_color_cherulean(palette = 'cheridis', discrete = T, name = 'Site') +
  scale_fill_cherulean(palette = 'cheridis', discrete = T, name = 'Site') +
  scale_shape_manual(values = 21:23, labels = c('MaxN', 'REST', 'UVC'), name = "Method")

f_pcoa2 <- ggplot(data = pcoadt) +
  geom_polygon(aes(x = Axis.3, y = Axis.4, fill = site_ID)) +
  geom_point(aes(x = Axis.3, y = Axis.4, fill = site_ID, shape = method), size = 3) +
  looks + labs(x = "PCo3", y = "PCo4") +
  scale_fill_cherulean(palette = 'cheridis', discrete = T, name = 'Site')+
  scale_shape_manual(values = 21:23, labels = c('MaxN', 'REST', 'UVC'), name = "Method")

# (f_pcoa1|f_pcoa2) + plot_layout(guides = "collect")
# I wonder if the scale of the abundances are affecting the scalings

nmatrix_rel <- apply(nmatrix[-(1:2)], 1, function(x) x/sum(x)) %>% t %>% as.data.frame
nmatrix_rel <- bind_cols(nmatrix[1:2], nmatrix_rel)
# try again
dis_rel <- vegdist(nmatrix_rel[-(1:2)], method = "euclidean", diag = F, binary = F)
method_pcoa_rel <- ape::pcoa(dis_rel, correction = "none", rn = paste(nmatrix_rel$site_ID, nmatrix_rel$method, sep="_"))

# scree plot
ggplot(data=method_pcoa_rel$values[1:7,], aes(x=1:7, y=Relative_eig/sum(method_pcoa_rel$values$Relative_eig))) +
  geom_line() + geom_point(shape=21, fill='white', size=3) + labs(x=NULL, y=NULL) +
  scale_x_continuous(breaks=c(1:10)) +
  theme_bw(base_size = 13) + labs(title = "Scree plot (PCoA, scaled abundances)")


relpcoadt <- as.data.frame(method_pcoa_rel$vectors) %>% bind_cols(., nmatrix_rel[1:2])
f_relpcoa1 <- ggplot(data = relpcoadt) +
  geom_polygon(aes(x = Axis.1, y = Axis.2, fill = site_ID, color = site_ID), alpha = 0.5) +
  geom_point(aes(x = Axis.1, y = Axis.2, fill = site_ID, shape = method), size = 3) +
  looks + labs(x = "PCo1", y = "PCo2")+
  scale_fill_cherulean(palette = "cheridis", discrete = T, name = "Site") +
  scale_color_cherulean(palette = "cheridis", discrete = T, name = "Site") +
  scale_shape_manual(values = 21:23, labels = c('MaxN', 'REST', 'UVC'), name = "Method")

f_relpcoa2 <- ggplot(data = relpcoadt) +
  geom_polygon(aes(x = Axis.3, y = Axis.4, fill = site_ID, color = site_ID), alpha = 0.5) +
  geom_point(aes(x = Axis.3, y = Axis.4, fill = site_ID, shape = method), size = 3) +
  looks + labs(x = "PCo3", y = "PCo4")+
  scale_fill_cherulean(palette = "cheridis", discrete = T, name = "Site") +
  scale_shape_manual(values = 21:23, labels = c('MaxN', 'REST', 'UVC'), name = "Method")

(f_pcoa1 / f_relpcoa1) + plot_layout(guides = "collect")


# COMPARE ASSEMBLAGES BETWEEN CAMERAS ---------------------------------------
# compare the assemblage composition with PCoA between MaxN cameras

nmatrix_cam <- data_maxn %>% group_by(site_ID, Camera, Taxon) %>% 
  summarise(n = sum(MaxN)) %>% ungroup
nmatrix_cam <- nmatrix_cam %>% tidyr::pivot_wider(names_from = Taxon, values_from = n) %>% 
  replace(is.na(.), 0) # replace NAs with 0

dis_maxn <- vegdist(nmatrix_cam[-(1:2)], method = "euclidean", diag = F, binary = F)
pcoa_maxn <- ape::pcoa(dis_maxn, correction = "none", rn = paste(nmatrix_cam$site_ID, nmatrix_cam$Camera, sep="_"))
biplot(pcoa_maxn)

# check scree
ggplot(data=pcoa_maxn$values[1:7,], aes(x=1:7, y=Cumul_eig)) +
  geom_line() + geom_point(shape=21, fill='white', size=3) + labs(x=NULL, y=NULL) +
  scale_x_continuous(breaks=c(1:10)) +
  looks + labs(title = "Scree plot")

# 4 dimensions is ok
## Plot PCoA results

dt_maxnpcoa <- as.data.frame(pcoa_maxn$vectors) %>% bind_cols(., nmatrix_cam[1:2])
# calculate convex hulls per site
convhull_cam <- list(as.list(rep(NA, 3)), as.list(rep(NA, 3)))
for (i in 1:3) {
  convhull_cam[[1]][[i]] <- convhull.vert(dt_maxnpcoa %>% 
                                            filter(site_ID == unique(site_ID)[i]) %>% 
                                            select(Axis.1, Axis.2)) %>% 
    mutate(site_ID = unique(dt_maxnpcoa$site_ID)[i])
  convhull_cam[[2]][[i]] <- convhull.vert(dt_maxnpcoa %>% 
                                            filter(site_ID == unique(site_ID)[i]) %>% 
                                            select(Axis.3, Axis.4)) %>% 
    mutate(site_ID = unique(dt_maxnpcoa$site_ID)[i])
}
convhull_cam[[1]] <- bind_rows(convhull_cam[[1]])
convhull_cam[[2]] <- bind_rows(convhull_cam[[2]])

f_maxncam <- as.list(rep(NA,3))
f_maxncam[[1]] <- ggplot(data = dt_maxnpcoa) +
  geom_polygon(data = convhull_cam[[1]], 
               aes(x = Axis.1, y = Axis.2, fill = site_ID)) +
  geom_point(aes(x = Axis.1, y = Axis.2, fill = site_ID), size = 3, shape = 21) +
  looks + labs(x = "PCo1", y = "PCo2")+
  scale_fill_cherulean(palette = "cheridis", discrete = T, name = 'Site')

f_maxncam[[2]] <- ggplot(data = dt_maxnpcoa) +
  geom_polygon(data = convhull_cam[[2]],
               aes(x = Axis.3, y = Axis.4, fill = site_ID)) +
  geom_point(aes(x = Axis.3, y = Axis.4, fill = site_ID), size = 3, shape = 21) +
  looks + labs(x = "PCo3", y = "PCo4")+
  scale_fill_cherulean(palette = "cheridis", discrete = T, name = 'Site')

(f_maxncam[[1]]|f_maxncam[[2]]) + plot_layout(guides = "collect")

# PAIRWISE ABUNDANCE COMPARISONS --------------------------------------------
# look at scatterplots first and then maybe fit OLS/GLMs

# make a dataframe of all methods joined
# each row/record an observed species-size class

dt_allmethods <- data_maxn %>% group_by(site_ID, Taxon, Size_class) %>% summarise(MaxN = sum(MaxN)) %>% 
  full_join(., data_rest[-4], by = c('site_ID', 'Taxon', 'Size_class')) %>% 
  rename(., REST = D)
dt_allmethods <- data_uvc %>% group_by(site_ID, Taxon, Size_class) %>% 
  summarise(UVC = sum(count)) %>% ungroup() %>% filter(str_detect(site_ID, 'kaku|hinalea|sunset')) %>% 
  full_join(., dt_allmethods, by = c('site_ID', 'Taxon', 'Size_class'))

# replace NAs with 0s
dt_allmethods$UVC[is.na(dt_allmethods$UVC)] <- 0
dt_allmethods$REST[is.na(dt_allmethods$REST)] <- 0
dt_allmethods$MaxN[is.na(dt_allmethods$MaxN)] <- 0

pretty <- c('^\\_' = '\\< ', '(?<=[:digit:])\\_' = '\\-')
dt_allmethods$Size_class <- str_replace_all(dt_allmethods$Size_class, pretty) %>% 
  factor(., levels = c('< 5', '5-9', '10-19', 
                      '20-29', '30-39', '40-49', '50-59'), ordered = T)

# long format
dt_all_long <- dt_allmethods %>% 
  tidyr::pivot_longer(., cols = UVC:REST, names_to = "method", values_to = "n")

# REST - MaxN
f_pairscatter <- as.list(rep(NA, 3))
f_pairscatter[[1]] <- ggplot(dt_allmethods) +
  geom_jitter(aes(x = REST, y = MaxN, fill = Size_class), width = 0.05, height = 0.05,
              shape = 21, alpha = 0.6, size = 2) +
  looks + scale_fill_cherulean(palette = "cheridis", discrete = T, reverse = T)

# REST - UVC
f_pairscatter[[2]] <- ggplot(dt_allmethods) +
  geom_jitter(aes(x = REST, y = UVC, fill = Size_class), width = 0.05, height = 0.05,
              shape = 21, alpha = 0.6, size = 2) +
  looks + scale_fill_cherulean(palette = "cheridis", discrete = T, reverse = T)

# UVC - MaxN
f_pairscatter[[3]] <- ggplot(dt_allmethods) +
  geom_jitter(aes(x = UVC, y = MaxN, fill = Size_class), width = 0.05, height = 0.05,
              shape = 21, alpha = 0.6, size = 2) +
  looks + scale_fill_cherulean(palette = "cheridis", discrete = T, reverse = T)

(f_pairscatter[[1]] + f_pairscatter[[2]] + f_pairscatter[[3]]) * scale_x_log10() * scale_y_log10() + plot_layout(guides = "collect")
(f_pairscatter[[1]] + f_pairscatter[[2]] + f_pairscatter[[3]]) + plot_layout(guides = "collect")

# size spectra
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

ggplot(sizespec) +
  geom_density(aes(x = midSize, fill = method, color = method), alpha = 0.4) +
  labs(x = "Size class (cm)", y = 'Relative frequency') +
  facet_grid(rows = vars(site_ID)) +
  looks + scale_fill_cherulean(palette = "gomphosus", discrete = T) + 
  scale_color_cherulean(palette = "gomphosus", discrete = T)

pairs(dt_allmethods[,4:6], log = 'xy', pch = 1:7, cex = 1.2)
# hinalea
pairs(dt_allmethods[dt_allmethods$site_ID == 'hale_hinalea', 4:6], log = 'xy', pch = 1:7, cex = 1.2)

## TALLIES --------------------------------------------
# exclusive species check

# calculate the "absences" by method
dt_all_long %>% group_by(method, site_ID, Taxon) %>% summarise(n = sum(n)) %>% # aggregate size classes
  ungroup() %>% group_by(method, site_ID) %>% 
  summarise(absences = length(which(n == 0))) %>% ungroup %>% group_by(method) %>% 
  summarise(meanMiss = mean(absences))

# not splitting by site
temp <- dt_allmethods %>% group_by(Taxon) %>% summarise(UVC = sum(UVC), REST = sum(REST), MaxN = sum(MaxN))
n_distinct(temp$Taxon)

# EXCLUSIVES
temp %>% filter(MaxN == 0 & UVC > 0) %>% pull(Taxon) %>% n_distinct #UVC
temp %>% filter(UVC == 0 & MaxN > 0) %>% pull(Taxon) %>% n_distinct # video
# overlaps, video uvc
temp %>% filter(MaxN > 0, UVC > 0) %>% pull(Taxon) %>% n_distinct
temp %>% filter(MaxN > 0, UVC > 0, REST > 0) %>% pull(Taxon) %>% n_distinct


# overlaps
temp %>% filter(UVC > 0 & MaxN > 0) %>% 
  summarise(overlaps = n()) # overlaps in video and uvc, by site
temp %>% ungroup() %>% filter(UVC > 0 & MaxN > 0) %>% 
 pull(Taxon) %>% n_distinct # overlap total video-uvc

temp %>% filter(MaxN > 0, REST == 0, UVC > 0) %>% View # maxN + UVC
temp %>% ungroup() %>% filter(MaxN > 0, REST > 0, UVC > 0) %>% nrow # all methods overlap

rm(temp)
# Species richness per method
dt_all_long %>% group_by(site_ID, method) %>% filter(n > 0) %>% summarise(S = length(unique(Taxon)))
temp <- dt_allmethods %>% group_by(site_ID, Taxon) %>% summarise(UVC = sum(UVC), REST = sum(REST), MaxN = sum(MaxN))
temp %>% filter(REST == 0, MaxN > 0) %>% summarise(missedsp = n_distinct(Taxon))


## Fit GLMs --------------------------------------------

set.seed(240)
# initially started out thinking about fitting Size_class as an effect, 
# but the size spectra makes me think that letting size absorb variation is not the way to go
# dominated by 10-19 cm

# read in traits to add in
traits <- read.csv('../data/traits_group.csv', header = T)
dt_abund <- left_join(dt_allmethods, traits, by = "Taxon")


model_methodpairs <- as.list(rep('', 6))
# MaxN ~ REST
model_methodpairs[[1]] <- lm(log(MaxN + 1) ~ log(REST + 1), data = dt_abund)
# UVC ~ REST
model_methodpairs[[2]] <- lm(log(UVC + 1) ~ log(REST + 1), data = dt_abund)
# MaxN ~ UVC
model_methodpairs[[3]] <- lm(log(MaxN + 1) ~ log(UVC + 1), data = dt_abund)

# candidates with sociality
model_methodpairs[[4]] <- lm(log(MaxN + 1) ~ log(REST + 1) + Group, data = dt_abund)
# UVC ~ REST
model_methodpairs[[5]] <- lm(log(UVC + 1) ~ log(REST + 1) + Group, data = dt_abund)
# MaxN ~ UVC
model_methodpairs[[6]] <- lm(log(MaxN + 1) ~ log(UVC + 1) + Group, data = dt_abund)

summary(model_methodpairs[[4]])
summary(model_methodpairs[[5]])
summary(model_methodpairs[[6]])

sapply(model_methodpairs, AIC)[c(1,4,2,5,3,6)] # look at AICs
# sociality is consistently lower

# model prediction
# make a standard new x that will be used by all models, with sociality held at a mean
social <- mean(dt_abund$Group)
newx <- data.frame(REST = seq(min(dt_abund$REST), max(dt_abund$REST), length.out = 50), 
                   UVC = seq(min(dt_abund$UVC), max(dt_abund$UVC), length.out = 50), 
                   Group = social)
predicts <- as.list(rep('', 3))
for (i in 1:3) { # each model
    predicts[[i]] <- predict(model_methodpairs[[i+3]], newdata = newx, se.fit = T)
    predicts[[i]] <- with(predicts[[i]], data.frame(fitted = fit, lwr = fit - 1.96 * se.fit, 
                                                    upr = fit + 1.96 * se.fit))
}
# add newx to the prediction data frame
predicts[[1]] <- predicts[[1]] %>% bind_cols(newx %>% select(Group))
predicts[[2]] <- predicts[[2]] %>% bind_cols(newx %>% select(REST, Group))
predicts[[3]] <- predicts[[3]] %>% bind_cols(newx %>% select(UVC, Group))

# partial regression plots for observed abundances (no sociality yet)
f_abundpair <- as.list(rep('', 3))
f_abundpair[[1]] <- ggplot(dt_abund) +
  geom_jitter(aes(y = log(MaxN + 1), x = log(REST + 1), fill = Size_class), width = 0.05, height = 0.05, 
               alpha = 0.6, size = 2, shape = 21) +
  geom_ribbon(data = predicts[[1]], 
              aes(x = log(REST + 1), ymin = lwr, ymax = upr), linetype = 'dashed', color = 'black', linewidth= 0.5, fill = 'transparent') +
  geom_line(data = predicts[[1]], aes(x = log(REST + 1), y = fitted)) +
  labs(x = 'log(REST)', y = 'log(MaxN)') +
  scale_fill_cherulean(palette = "cheridis", discrete = T, name = "Size class (cm)")

# REST - UVC
f_abundpair[[2]] <- ggplot(dt_allmethods) +
  geom_jitter(aes(y = log(UVC + 1), x = log(REST + 1), fill = Size_class), width = 0.05, height = 0.05, 
             alpha = 0.6, size = 2, shape =21) +
  geom_ribbon(data = predicts[[2]], aes(x = log(REST + 1), ymin = lwr, ymax = upr), 
              linetype = 'dashed', color = 'black', linewidth= 0.5, fill = 'transparent') +
  geom_line(data = predicts[[2]], aes(x = log(REST + 1), y = fitted)) +
  labs(x = 'log(REST)', y = 'log(Point count)')+
  scale_fill_cherulean(palette = "cheridis", discrete = T, name = "Size class (cm)")

# UVC - MaxN
f_abundpair[[3]] <- ggplot(dt_allmethods) +
  geom_jitter(aes(y = log(MaxN + 1), x = log(UVC + 1), fill = Size_class), width = 0.05, height = 0.05, 
               alpha = 0.6, size = 2, shape = 21) +
  geom_ribbon(data = predicts[[3]], 
              aes(x = log(UVC + 1), ymin = lwr, ymax = upr), linetype = 'dashed', color = 'black', linewidth= 0.5, fill = 'transparent') +
  geom_line(data = predicts[[3]], aes(x = log(UVC + 1), y = fitted)) +
  labs(x = 'log(Point count)', y = 'log(MaxN)')+
  scale_fill_cherulean(palette = "cheridis", discrete = T, name = "Size class (cm)")

(f_abundpair[[1]] + f_abundpair[[2]] + f_abundpair[[3]]) * looks * coord_cartesian(xlim = c(0,4.75), ylim = c(0,4.75)) + plot_layout(guides = "collect")

# partial regression plots for partial sociality effects
# generate partial predictions for sociality, held at abundance mean
meanREST <- mean(dt_abund$REST)
meanPoint <- mean(dt_abund$UVC)
newx <- data.frame(Group = 1:5, REST = meanREST, UVC = meanPoint)
for (i in 4:6) { # each model
  predicts[[i]] <- predict(model_methodpairs[[i]], newdata = newx, se.fit = T)
  predicts[[i]] <- with(predicts[[i]], data.frame(fitted = fit, lwr = fit - 1.96 * se.fit, 
                                                  upr = fit + 1.96 * se.fit,
                                                  Group = newx$Group))
}

f_abundsoc <- as.list(rep('', 3))

f_abundsoc[[1]] <- ggplot(dt_abund) +
  ggdist::stat_slab(aes(x = Group, y = log(MaxN + 1)), 
                       adjust = .5, justification = -.5, 
                       color = 'grey20', fill = 'grey70', size = 0.5, alpha = 0.4) + 
  geom_point(aes(x = Group, y = log(MaxN + 1)), size = 1, alpha = .4, shape = 21, color = 'grey30', 
             fill = 'white', position = position_jitter(width = .1)) +
  geom_errorbar(data = predicts[[4]], aes(x = Group + 0.25, ymin = lwr, ymax = upr), 
                  color = 'grey20', width = 0.2) +
  geom_point(data = predicts[[4]], aes(x = Group + 0.25, y = fitted), 
                  color = 'grey20', size = 2.5) +
  labs(x = 'Sociality', y = 'log(MaxN)')

f_abundsoc[[2]] <- ggplot(dt_abund) +
  ggdist::stat_slab(aes(x = Group, y = log(UVC + 1)), 
                    adjust = .5, justification = -.5, 
                    color = 'grey20', fill = 'grey70', size = 0.5, alpha = 0.4) + 
  geom_point(aes(x = Group, y = log(UVC + 1)), size = 1, alpha = .4, shape = 21, color = 'grey20', 
             fill = 'white', position = position_jitter(width = .1)) +
  geom_errorbar(data = predicts[[5]], aes(x = Group + 0.25, ymin = lwr, ymax = upr), 
                 color = 'grey20', width = 0.2) +
  geom_point(data = predicts[[5]], aes(x = Group + 0.25, y = fitted), 
             color = 'grey20', size = 2.5) +
  labs(x = 'Sociality', y = 'log(Point count)')

f_abundsoc[[3]] <- ggplot(dt_abund) +
  ggdist::stat_slab(aes(x = Group, y = log(MaxN + 1)), 
                    adjust = .5, justification = -.5, 
                    color = 'grey20', fill = 'grey70', size = 0.5, alpha = 0.4) + 
  geom_point(aes(x = Group, y = log(MaxN + 1)), size = 1, alpha = .4, shape = 21, color = 'grey20', 
             fill = 'white', position = position_jitter(width = .1)) +
  geom_errorbar(data = predicts[[6]], aes(x = Group + 0.25, ymin = lwr, ymax = upr), 
                color = 'grey20', width = 0.2) +
  geom_point(data = predicts[[6]], aes(x = Group + 0.25, y = fitted), 
             color = 'grey20', size = 2.5) +
  labs(x = 'Sociality', y = 'log(MaxN)')

(f_abundsoc[[1]] + f_abundsoc[[2]] + f_abundsoc[[3]]) * looks * scale_x_continuous(breaks = c(1:5)) * coord_cartesian(xlim = c(0.5,5.7), ylim = c(0,4.75), clip = 'off') + plot_layout(guides = "collect")

# SADs --------------------------------------------
# do it on size aggregated data
temp <- dt_all_long %>% group_by(site_ID, method) %>% summarise(totaln = sum(n))
# add relative abundances
dt_all_long <- left_join(dt_all_long, temp, by = c("site_ID", "method")) %>% 
  mutate(reln = n/totaln) %>% select(!totaln)

ggplot(data = dt_all_long %>% filter(reln > 0) %>% group_by(site_ID, method, Taxon) %>% 
         summarise(reln = sum(reln))) +
  geom_density(aes(x = reln, color = method), linewidth = 0.75) +
  labs(x = 'Abundances', y = "Species frequency") + looks +
  scale_color_cherulean(palette = "cheridis", discrete = T, name = 'Method')

SAD <- ggplot(data = dt_all_long %>% filter(reln > 0) %>% group_by(site_ID, method, Taxon) %>% 
         summarise(reln = sum(reln))) +
  geom_density(aes(x = reln, color = method), linewidth = 0.5) +
  labs(x = 'Abundances', y = "Species frequency") + looks +
  scale_color_cherulean(palette = "gomphosus", discrete = T, name = 'Method') +
  facet_wrap(vars(site_ID))

# CALCULATE ASSEMBLAGE METRICS -----------------------------------------------------------------

## ASSEMBLAGE STANDING BIOMASS -----------------------------------------------------------
require(rfishbase) # load fishbase

# extract length weight and length length data
# dt_LW <- length_weight(species_list = sort(unique(dt_all_long$Taxon)))
# dt_LW <- dt_LW %>% select(Species, Type, a, b, Number, Locality)
# write.csv(dt_LW, '../data/data_LW.csv', row.names = F) # select the ones to use manually

dt_LW <- read.csv('../data/data_LWab.csv', header = T)[1:4]

# match length conversions
dt_LL <- length_length(species_list = sort(unique(dt_LW$Species[dt_LW$Type != 'TL'])))
n_distinct(dt_LW$Species[dt_LW$Type != 'TL']) # how many species, 15
# only species that need conversions
dt_LL <- dt_LL %>% select(Species, Length1, Length2, a, b) %>%
  filter(Length2 == 'TL' | Length1 == 'TL')  # TL in length1 could be solved by algebra to inverse
colnames(dt_LL)[4:5] <- c('LLa', 'LLb')

dt_LL$inverse <- 0
dt_LL$inverse[dt_LL$Length1 == 'TL'] <- 1

# need to pick out the LL records that fit what I need
dt_LL$Type <- dt_LL$Length1
dt_LL$Type[dt_LL$inverse == 1] <- dt_LL$Length2[dt_LL$inverse == 1]
# make a column that just represents what is being converted from
dt_LL <- dt_LL[-c(2:3)]
dt_LL <- dt_LL %>% arrange(Species, inverse, Type, LLb) %>% 
  filter(LLb != 1) # anything that is exactly 1 is likely wrong

dt_LL <- dt_LL %>% group_by(Species, inverse, Type) %>% 
  summarise(LLa = mean(LLa), LLb = mean(LLb)) # take the mean if there are multiple

dt_LW <- left_join(dt_LW, dt_LL %>% filter(inverse == 0), by=c("Species", "Type"))
temp_LW <- left_join(dt_LW, dt_LL %>% filter(inverse == 1), by=c("Species", "Type")) %>% 
  filter(is.na(inverse.x), is.na(inverse.y) == F)

# manual join offline
# write.csv(dt_LW, '../data/data_LW_LL.csv', row.names = F)
# clipr::write_clip(temp_LW)
rm(temp_LW)
dt_LW <- read.csv('../data/data_LW_LL.csv', header = T)

# calculate biomass for each species-size
data_biomass <- dt_all_long
data_biomass$midSize <- str_extract(data_biomass$Size_class, '(?<=\\-|(\\<\\s))[:digit:]+$') %>% 
  as.numeric
data_biomass$midSize[data_biomass$midSize > 5] <- data_biomass$midSize[data_biomass$midSize > 5] -4
data_biomass$midSize[data_biomass$midSize == 5] <- data_biomass$midSize[data_biomass$midSize == 5] - 2.5

data_biomass <- left_join(data_biomass, dt_LW[1:7], by=c("Taxon" = "Species"))
data_biomass <- data_biomass %>% filter(!is.na(Type)) # remove uncertain records that can't be matched
# W = a * L^b
# length conversion first, midSize TL -> right length type -> biomass
# dt_LW inverse = 1 means TL = L1
# inverse = 1: L2 = a + b*L1
# inverse = 0:  = (L2 - a) / b 
data_biomass$L <- data_biomass$midSize
data_biomass$L[data_biomass$Type != 'TL' & data_biomass$inverse == 0] <- with(data_biomass[data_biomass$Type != 'TL' & data_biomass$inverse == 0, ], (midSize - LLa) / LLb)
data_biomass$L[data_biomass$Type != 'TL' & data_biomass$inverse == 1] <- with(data_biomass[data_biomass$Type != 'TL' & data_biomass$inverse == 1, ], LLa + LLb * midSize)
data_biomass$biomass <- 0
data_biomass$biomass <- with(data_biomass, (a * L ^ b) * n) # length -> weight * number of individuals
dt_metrics <- data_biomass %>% filter(n > 0) %>% 
  group_by(site_ID, method) %>% summarise(SSB = sum(biomass), S = n_distinct(Taxon))

detach('package:rfishbase')
rm(dt_LL, dt_LW)
dt_metrics$SSB <- dt_metrics$SSB / 1000 # in kilograms

# dt_traits <- read.csv('../data/traits_manual.csv', header = T) %>% select(Species, Schooling, RemarksRefs) %>% 
#   right_join(., dt_LW, by="Species") %>% arrange(Species)


## Richness measures -----------------------------------------------------------

require(SpadeR)

dt_S <- dt_all_long %>% select(!c(Size_class, n)) %>% group_by_at(vars(-reln)) %>% 
  summarise(reln = sum(reln)) %>% 
  mutate(names = paste(site_ID, method, sep="_")) %>% ungroup %>%  # aggregate size classes
  select(!c(site_ID, method)) %>% 
  tidyr::pivot_wider(names_from = names, values_from = reln)
dt_S[is.na(dt_S)] <- 0
dt_S[-1] <- round(dt_S[-1] * 500, 0) # turn proportional abundance into "counts"

expS <- as.list(rep(NA, ncol(dt_S) - 1)) # empty list to populate with results
for (i in 2:ncol(dt_S)) {
  expS[[i-1]] <- ChaoSpecies(dt_S[i][dt_S[i] > 0], datatype = "abundance", k=10, conf=0.95)
} # calculate expected species richness for each site-method without zeroes

dt_metrics$expS <- sapply(expS, function(x) x$Species_table[4])
dt_metrics$expSlwr <- sapply(expS, function(x) x$Species_table[4,3])
dt_metrics$expSupr <- sapply(expS, function(x) x$Species_table[4,4])
# extract the iChao-bc metric, which accounts for detection differences and bias correction
# this only works if the site ID and method order matches exactly with dt_metrics!

require(iNEXT)
# Calculate Hill numbers, effective number of species instead of an evenness index
hill <- as.list(rep(NA, ncol(dt_S) - 1)) # empty list to populate with results
for (i in 2:ncol(dt_S)) {
  hill[[i-1]] <- iNEXT(dt_S[i][dt_S[i] > 0], datatype = "abundance", q = c(1,2))
} # calculate expected species richness for each site-method without zeroes
# q0 = species richness
# q1 = exponential Shannon entropy, rare species have greater influence
# q2 = inverse Simpson index, dominant species have greater influence
dt_metrics$q1 <- sapply(hill, function(x) x$iNextEst$size_based %>% filter(m == max(m), Order.q == 1) %>% pull(qD))
dt_metrics$q2 <- sapply(hill, function(x) x$iNextEst$size_based %>% filter(m == max(m), Order.q == 2) %>% pull(qD))
dt_metrics$q1lwr <- sapply(hill, function(x) x$iNextEst$size_based %>% filter(m == max(m), Order.q == 1) %>% pull(qD.LCL))
dt_metrics$q1upr <- sapply(hill, function(x) x$iNextEst$size_based %>% filter(m == max(m), Order.q == 1) %>% pull(qD.UCL))
dt_metrics$q2lwr <- sapply(hill, function(x) x$iNextEst$size_based %>% filter(m == max(m), Order.q == 2) %>% pull(qD.LCL))
dt_metrics$q2upr <- sapply(hill, function(x) x$iNextEst$size_based %>% filter(m == max(m), Order.q == 2) %>% pull(qD.UCL))

dt_metrics <- dt_metrics[c(1:5,8,9)] %>% tidyr::pivot_longer(cols = -(1:2), names_to = "metrics", values_to = "value") %>% 
  left_join(., dt_metrics[c(1,2,6,7,10:13)] %>% tidyr::pivot_longer(cols = -(1:2), names_pattern = "^(.*)(lwr|upr)$", names_to = c("metrics", ".value")), by = c("site_ID", "method", "metrics"))

f_metrics <- as.list(rep('', 5))
for (i in 1:n_distinct(dt_metrics$metrics)) {
  f_metrics[[i]] <- ggplot(dt_metrics %>% filter(metrics == dt_metrics$metrics[i])) +
  geom_pointrange(aes(x = site_ID, y = value, ymin = lwr, ymax = upr, shape = method, fill = method)) +
  scale_fill_cherulean(palette = 'gomphosus', discrete = T) + looks +
  scale_shape_manual(values = 21:23) +
  scale_x_discrete(labels = c("HH", "HK", "SP")) +
  labs(x = '', y = dt_metrics$metrics[i])
}

f_metrics[[3]] <- f_metrics[[3]] + scale_y_log10()
f_metrics[[1]] <- f_metrics[[1]] + labs(y = 'SSB (kg)', x = '')

f_metrics[[1]] + f_metrics[[2]] + f_metrics[[3]] + f_metrics[[4]] + f_metrics[[5]] + plot_layout(guides = "collect")

## PCoA on metrics -----------------------------------------------------------

set.seed(24)
metrics_dis <- vegdist(dt_metrics[-(1:2)], method = "euclidean", diag = F, binary = F)
metrics_pcoa <- ape::pcoa(metrics_dis, correction = "none")

ggplot(data=metrics_pcoa$values[1:5,], aes(x=1:5, y=Relative_eig/sum(metrics_pcoa$values$Relative_eig))) +
  geom_line() + geom_point(shape=21, fill='white', size=3) + labs(x=NULL, y=NULL) +
  scale_x_continuous(breaks=c(1:10)) +
  theme_bw(base_size = 13) + labs(title = "Scree plot (PCoA)")

biplot(metrics_pcoa, rn = NULL)
metrics_env <- envfit(metrics_pcoa$vectors, env = dt_metrics[-(1:2)], perm = 999)
plot(metrics_env)

dt_metricspcoa <- as.data.frame(metrics_pcoa$vectors) %>% bind_cols(., dt_metrics[1:2])
metrics_vec <- 75 * scores(metrics_env, "vectors") %>% as.data.frame
metrics_vec$metrics <- rownames(metrics_vec)

f_metricpcoa <- ggplot(data = dt_metricspcoa) +
  geom_segment(data = metrics_vec, aes(x = 0, y = 0, xend = Axis.1, yend = Axis.2), 
               arrow = arrow(length = unit(1.5, "mm"), type = "closed"), color = 'grey') +
  geom_text(data = metrics_vec, aes(x = Axis.1 + 5, y = Axis.2 + 5, label = metrics),
                           size = 4, color = 'grey50') +
  geom_polygon(aes(x = Axis.1, y = Axis.2, fill = site_ID)) +
  geom_point(aes(x = Axis.1, y = Axis.2, fill = site_ID, shape = method), size = 3) +
  looks + labs(x = "PCo1", y = "PCo2")+
  scale_fill_cherulean(palette = "cheridis", discrete = T, name = "Site") +
  scale_color_cherulean(palette = "cheridis", discrete = T, name = "Site") +
  scale_shape_manual(values = 21:23, labels = c('MaxN', 'REST', 'UVC'), name = "Method") +
  coord_fixed()

f_metricpcoa

# ACCUMULATION/RAREFACTION CURVE -----------------------------------------------------------
# Model species accumulation curves

temp <- dt_all_long %>% group_by(site_ID, Taxon, method) %>% 
  summarise(n = sum(n), reln=sum(reln)) %>% ungroup %>% filter(n > 0) %>% group_split(site_ID, method)
temp2 <- dt_all_long %>% ungroup() %>% filter(n > 0) %>% group_by(site_ID, method) %>% summarise(total = sum(round(n, 0))) # use for sample size, as in total number of individuals
# temp2 <- dt_all_long %>% ungroup() %>% filter(n > 0) %>% group_by(site_ID, method) %>% summarise(total = length(unique(Taxon))) # use for sample size, as in total number of individuals

# hard to use timestamps for point count species accumulation, so we'll do it by abundance-weighted sampling
dt_SAC <- as.list(rep('', 9))
for(i in 1:9) {
  dt_SAC[[i]] <- data.frame(site_ID = temp[[i]]$site_ID[1], method = temp[[i]]$method[1],
                            sample = 1:temp2$total[i], nsp = NA) # empty dataframe
  tempsp <- data.frame(sample = NA, sp = NA) # dummy to append species accumulation to
  temp3 <- temp[[i]]$Taxon[rep(1:nrow(temp[[i]]), ceiling(temp[[i]]$n))]
  for (j in 1:temp2$total[i]) {
    tempsp <- rbind(tempsp,
      data.frame(sample = j, sp = sample(temp3, size = j, replace = F)))
    # sample the assemblage using relative abundances as probability
    # temporarily store sampled species
    dt_SAC[[i]][j,4] <- tempsp %>% filter(sample <= j, is.na(sp) == F) %>% pull(sp) %>% n_distinct(.) 
    # calculate species accumulative
  } 
  dt_SAC[[i]]$site_ID <- temp[[i]]$site_ID[1]
  dt_SAC[[i]]$method <- temp[[i]]$method[1]
}

ggplot(data = bind_rows(dt_SAC)) +
  geom_line(aes(x = sample, y = nsp, group = interaction(site_ID, method),
                color = site_ID), size = 1) + looks + scale_x_log10()

# now we'll fit accumulation models
fit_SAC <- as.list(rep(NA, 9))
SAC_pred <- as.list(rep(NA, 9))
set.seed(240)
for (i in 1:9) { # beta P model
  optout <- tryCatch(nls(data = na.omit(dt_SAC[[i]]),
                         formula = nsp ~ a - ((a-b)*exp(-c*sample)),
                         start = list(a = 10, b = 1, c = 0.1)),
                     error = function(e) e)
  if(inherits(optout, "error")) next # if an error message gets generated this run, move to the next iteration of the loop
  fit_SAC[[i]] <- optout
}

fit_SAC[[9]] <- nls(data = na.omit(dt_SAC[[9]]),
    formula = nsp ~ a - ((a-b)*exp(-c*sample)),
    start = list(a = 20, b = 1, c = 0.1))

# for (i in 1:9) { # power model
#   optout <- tryCatch(nls(formula = nsp ~ k + s * log(sample), 
#                          data = na.omit(dt_SAC[[i]]), 
#                          start = list(k = 0.1, s = 0.1)),
#                      error = function(e) e)
#   if(inherits(optout, "error")) next # if an error message gets generated this run, move to the next iteration of the loop
#   fit_SAC[[i]] <- optout
# }
require(nlstools)
for (i in 1:9) {
  SAC_pred[[i]] <- nlsBootPredict(
    nlsBoot(fit_SAC[[i]], niter=500), 
    newdata = dt_SAC[[i]][3], interval = "confidence") %>% as.data.frame
  colnames(SAC_pred[[i]])[2:3] <- c('lwr', 'upr')
  SAC_pred[[i]] <- data.frame(site_ID = temp[[i]]$site_ID[1], method = temp[[i]]$method[1], 
                              sample = dt_SAC[[i]][3]) %>% bind_cols(., SAC_pred[[i]])
}

rm(temp, temp2, tempsp)

SAC_rar <- ggplot() +
  geom_point(data = bind_rows(dt_SAC), aes(x = sample, y = nsp, group = interaction(site_ID, method), color = method), 
             alpha = 0.2, shape = 21) +
  geom_ribbon(data = bind_rows(SAC_pred), aes(x = sample, ymin = lwr, ymax = upr, group = interaction(site_ID, method),
                                            fill = method), alpha = 0.4) +
  geom_line(data = bind_rows(SAC_pred), aes(x = sample, y = Median, group = interaction(site_ID, method),
                                            color = method)) +
  scale_color_cherulean(palette = "gomphosus", discrete = T, name = "Method") +
  scale_fill_cherulean(palette = "gomphosus", discrete = T, name = "Method") +
  facet_wrap(vars(site_ID), scales = 'free_x') + looks + scale_x_log10() +
  labs(x = "Number of samples", y = "Number of species")

## SACs for video methods

## SAC data set up and prep
data_ruv2 <- read.csv('../data/ruv_himb_pilot.csv', header = T)
require(lubridate)
data_ruv2$entrytime_c <- ms(data_ruv2$Time_entry) %>% as.duration
data_ruv2$exittime_c <- ms(data_ruv2$Time_exit) %>% as.duration
# calculate staying time duration
data_ruv2$staytime <- with(data_ruv2, exittime_c - entrytime_c) %>% as.numeric
summary(data_ruv2)
data_ruv <- data_ruv %>% filter(staytime > 0) # remove any zero second occurrences
data_ruv2 <- data_ruv2 %>% arrange(site_ID, Camera, VidFile, entrytime_c)
#  timestamps are by video but we want them to be timestamps relative to total observation time, not video time.
cameras <- data_ruv2 %>% distinct(site_ID, Camera, VidFile)
ncam <- cameras %>% group_by(site_ID, Camera) %>% summarise(ncam = n_distinct(VidFile)) %>% pull(ncam)
seq <- ''
  for (i in 1:length(ncam)) {
    seq <- c(seq, 1:ncam[i])
  } # generate a sequence vector to identify the nominal order of video files
seq <- seq[-1] # trim that dummy start

cameras$VidSeq <- seq

data_ruv2 <- left_join(data_ruv2, cameras, by=c('site_ID', 'Camera', 'VidFile'))
data_ruv2$VidSeq <- as.numeric(data_ruv2$VidSeq)

# make the timestamps relative to site and not video file
# all videos duration = 11:48 except for deep_0 = 17:42
for (i in 2:5) {
  # loop for each video file in order(1-4 or 5) add a multiple of the video file length
  # 2:5 because the first files don't need fixing
  # for the 11:48 videos
  data_ruv2$SACentry[data_ruv2$VidSeq == i & data_ruv2$Camera != 'deep_0'] <- data_ruv2$entrytime_c[data_ruv2$VidSeq == i & data_ruv2$Camera != 'deep_0'] + (i-1)*(dminutes(11) + dseconds(48))
  # for deep_0 17:42 videos
  data_ruv2$SACentry[data_ruv2$VidSeq == i & data_ruv2$Camera == 'deep_0'] <- data_ruv2$entrytime_c[data_ruv2$VidSeq == i & data_ruv2$Camera == 'deep_0'] + (i-1)*(dminutes(17) + dseconds(42))
}

# just put the time stamps in for first video files
data_ruv2$SACentry[data_ruv2$VidSeq == 1] <- data_ruv2$entrytime_c[data_ruv2$VidSeq == 1]
data_ruv2$VidSeq <- NULL # get rid of these dummy objects
rm(cameras, ncam, seq)

## Calculate cumulative species richness over time
# make a dataframe to populate species accumulation per timestamp
SAC <- data_ruv2 %>% mutate(site_cam = paste(site_ID, Camera, sep = "_")) %>% 
  select(site_ID, site_cam, SACentry)  # matches the row index of data_ruv

# make the SAC times relative to the time of first observation
for (i in 1:11) {
  time1 <- SAC$SACentry[which(SAC$site_cam == unique(SAC$site_cam)[i])[1]]
  for (j in which(SAC$site_cam == unique(SAC$site_cam)[i])) {
    SAC$SACentry[j] <- SAC$SACentry[j] - time1 # take every time stamp and subtract time1 from it
  }
}
rm(time1)

SAC$spN <- NA
# for every row at every site-camera, calculate the species number accumulation through time
for (i in 1:11) {
  rows = seq(which(SAC$site_cam == unique(SAC$site_cam)[i])[1], max(which(SAC$site_cam == unique(SAC$site_cam)[i])), by = 3)
  for (j in rows) {
    if (j == rows[1]) { # the first row from each site represents the first species record, so these will always start at 1
      SAC$spN[j] <- 1
    }
    else {
      SAC$spN[j] <- data_ruv$Taxon[rows[1]:j] %>% n_distinct # number of unique species from the first row of the site to row j
    }
  }
}
SAC <- SAC %>% filter(!is.na(spN))
SAC$method <- 'MaxN'

# do again, filtering out the species that REST didn't have
SAC2 <- data_ruv2 %>% mutate(site_cam = paste(site_ID, Camera, sep = "_"), site_sp = paste(site_ID, Taxon, sep = "_")) %>% 
  filter(site_sp %in% with(data_rest, paste(site_ID, Taxon, sep="_"))) %>% 
  select(site_ID, site_cam, site_sp, SACentry)  # matches the row index of data_ruv

# make the SAC times relative to the time of first observation
for (i in 1:11) {
  time1 <- SAC2$SACentry[which(SAC2$site_cam == unique(SAC2$site_cam)[i])[1]]
  for (j in which(SAC2$site_cam == unique(SAC2$site_cam)[i])) {
    SAC2$SACentry[j] <- SAC2$SACentry[j] - time1 # take every time stamp and subtract time1 from it
  }
}
rm(time1)

SAC2$spN <- NA
# for every row at every site-camera, calculate the species number accumulation through time
for (i in 1:11) {
  # don't calculate for every row... every 3rd row
  rows = seq(which(SAC2$site_cam == unique(SAC2$site_cam)[i])[1], max(which(SAC2$site_cam == unique(SAC2$site_cam)[i])), by = 3)
  for (j in rows) {
    if (j == rows[1]) { # the first row from each site represents the first species record, so these will always start at 1
      SAC2$spN[j] <- 1
    }
    else {
      SAC2$spN[j] <- SAC2$site_sp[rows[1]:j] %>% n_distinct # number of unique species from the first row of the site to row j
    }
  }
}

# merge + cleanup
SAC2$method <- 'REST'
SAC <- bind_rows(SAC, SAC2)
rm(SAC2, rows)
SAC <- SAC %>% filter(!is.na(spN))

fit_SAC2 <- as.list(rep(NA, 6))
for (i in 1:3) {
  for (j in 1:2) {
    fit_SAC2[[i + (j-1)*3]] <- nls(data = SAC %>% 
                                     filter(site_ID == unique(SAC$site_ID)[i], method == unique(SAC$method)[j]),
                         formula = spN ~ a - ((a-b)*exp(-c*SACentry)),
                         start = list(a = 15, b = 1, c = 0.001))
  }
}

SAC_pred2 <- as.list(rep(NA, 6))
newx <- as_tibble(seq(1, 2700, length.out = 200))
names(newx) <- 'SACentry'
for (i in 1:3) {
  for (j in 1:2) {
    SAC_pred2[[i + (j-1)*3]] <- nlsBootPredict(
      nlsBoot(fit_SAC2[[i + (j-1)*3]], niter=500), 
      newdata = newx, interval = "confidence") %>% as.data.frame
    colnames(SAC_pred2[[i + (j-1)*3]]) <- c('fit', 'lwr', 'upr')
    SAC_pred2[[i + (j-1)*3]] <- data.frame(site_ID = unique(SAC$site_ID)[i], method = unique(SAC$method)[j],
                                 SACentry = newx) %>% bind_cols(., SAC_pred2[[i + (j-1)*3]])
  }
}


SAC_vid <- ggplot() +
  geom_point(data = SAC, aes(x = SACentry, y = spN, group = site_cam, color = method), 
             alpha = 0.1, shape = 21) +
  geom_ribbon(data = bind_rows(SAC_pred2), aes(x = SACentry, ymin = lwr, ymax = upr, group = interaction(site_ID, method),
                                              fill = method), alpha = 0.4) +
  geom_line(data = bind_rows(SAC_pred2), aes(x = SACentry, y = fit, group = interaction(site_ID, method),
                                            color = method)) +
  scale_color_manual(values = cherulean_palettes[['gomphosus']][c(1,4)], name = 'Method') +
  scale_fill_manual(values = cherulean_palettes[['gomphosus']][c(1,4)], name = 'Method') + looks +
  facet_wrap(vars(site_ID)) +
  labs(x = "Time elapsed (s)", y = "Number of species") + guides(color = 'none', fill = 'none')

(SAC_rar / SAC_vid) * theme(legend.position = 'bottom')
