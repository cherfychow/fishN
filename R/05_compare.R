

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

###############################################
# PCOA COMPARISON -----------------------------------------------------------------
sdf
data_uvc <- read.csv('../data/uvc_himb.csv', header = T)
# prep REST to match the assemblage survey data like RUV or MaxN

data_uvc$Size_class[data_uvc$TL_cm <= 5] <- "_5"
data_uvc$Size_class[data_uvc$TL_cm > 5 & data_uvc$TL_cm <= 10] <- "5_9"
data_uvc$Size_class[data_uvc$TL_cm > 10 & data_uvc$TL_cm <= 20] <- "10_19"
data_uvc$Size_class[data_uvc$TL_cm > 20 & data_uvc$TL_cm <= 30] <- "20_29"
data_uvc$Size_class[data_uvc$TL_cm > 30 & data_uvc$TL_cm <= 40] <- "30_39"
data_uvc$Size_class[data_uvc$TL_cm > 40 & data_uvc$TL_cm <= 50] <- "40_49"
data_uvc$Size_class[data_uvc$TL_cm > 50 & data_uvc$TL_cm <= 60] <- "50_59"

data_rest <- data_ruv %>% distinct(site_ID, Taxon, Size_class, spsize) %>% 
  left_join(., output_poisson[c(1,3,4,12)], by=c("site_ID", "Taxon", "Size_class"))

# cross-method comparison doesn't have to deal with size classes
# construct the community matrix
nmatrix <- data_uvc %>% filter(str_detect(site_ID, 'kinalea|kaku|sunset')) %>% 
  group_by(site_ID, Taxon) %>% summarise(n = sum(count)) %>% mutate(method = 'uvc') %>% ungroup # filter and aggregate UVC
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
  geom_polygon(aes(x = Axis.1, y = Axis.2, fill = site_ID), alpha = 0.5) +
  geom_point(aes(x = Axis.1, y = Axis.2, fill = site_ID, shape = method), size = 3) +
  looks + labs(x = "PCo1", y = "PCo2") +
  scale_color_cherulean(palette = 'cheridis', discrete = T, name = 'Site') +
  scale_fill_cherulean(palette = 'cheridis', discrete = T, name = 'Site') +
  scale_shape_manual(values = 21:23, labels = c('MaxN', 'REST', 'UVC'), name = "Method")

f_pcoa2 <- ggplot(data = pcoadt) +
  geom_polygon(aes(x = Axis.3, y = Axis.4, fill = site_ID), alpha = 0.5) +
  geom_point(aes(x = Axis.3, y = Axis.4, fill = site_ID, shape = method), size = 3) +
  looks + labs(x = "PCo3", y = "PCo4") +
  scale_fill_cherulean(palette = 'cheridis', discrete = T, name = 'Site')+
  scale_shape_manual(values = 21:23, labels = c('MaxN', 'REST', 'UVC'), name = "Method")

(f_pcoa1|f_pcoa2) + plot_layout(guides = "collect")
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
  geom_polygon(aes(x = Axis.1, y = Axis.2, fill = site_ID), alpha = 0.5) +
  geom_point(aes(x = Axis.1, y = Axis.2, fill = site_ID, shape = method), size = 3) +
  looks + labs(x = "PCo1", y = "PCo2")+
  scale_fill_cherulean(palette = "cheridis", discrete = T, name = "Site") +
  scale_shape_manual(values = 21:23, labels = c('MaxN', 'REST', 'UVC'), name = "Method")

f_relpcoa2 <- ggplot(data = relpcoadt) +
  geom_polygon(aes(x = Axis.3, y = Axis.4, fill = site_ID), alpha = 0.5) +
  geom_point(aes(x = Axis.3, y = Axis.4, fill = site_ID, shape = method), size = 3) +
  looks + labs(x = "PCo3", y = "PCo4")+
  scale_fill_cherulean(palette = "cheridis", discrete = T, name = "Site") +
  scale_shape_manual(values = 21:23, labels = c('MaxN', 'REST', 'UVC'), name = "Method")

f_relpcoa <- (f_relpcoa1|f_relpcoa2) + plot_layout(guides = "collect")
f_relpcoa

###############################################
# COMPARE ASSEMBLAGES BETWEEN CAMERAS ---------------------------------------
# compare the assemblage composition with PCoA between MaxN cameras

nmatrix_max <- data_maxn %>% group_by(site_ID, Camera, Taxon) %>% summarise(n = sum(MaxN)) %>% ungroup
nmatrix_max <- nmatrix_max %>% tidyr::pivot_wider(names_from = Taxon, values_from = n) %>% 
  replace(is.na(.), 0) # replace NAs with 0

dis_maxn <- vegdist(nmatrix_max[-(1:2)], method = "euclidean", diag = F, binary = F)
pcoa_maxn <- ape::pcoa(dis_maxn, correction = "none", rn = paste(nmatrix_max$site_ID, nmatrix_max$Camera, sep="_"))
biplot(pcoa_maxn)

# check scree
ggplot(data=pcoa_maxn$values[1:7,], aes(x=1:7, y=Cumul_eig)) +
  geom_line() + geom_point(shape=21, fill='white', size=3) + labs(x=NULL, y=NULL) +
  scale_x_continuous(breaks=c(1:10)) +
  looks + labs(title = "Scree plot")

# 4 dimensions is ok
## Plot PCoA results

dt_maxnpcoa <- as.data.frame(pcoa_maxn$vectors) %>% bind_cols(., nmatrix_max[1:2])
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
  geom_polygon(data = convhull_cam[[1]], alpha = 0.3, color = 'grey50', 
               aes(x = Axis.1, y = Axis.2, fill = site_ID)) +
  geom_point(aes(x = Axis.1, y = Axis.2, fill = site_ID), size = 3, shape = 21) +
  looks + labs(x = "PCo1", y = "PCo2")+
  scale_fill_cherulean(palette = "gomphosus", discrete = T, name = 'Site')

f_maxncam[[2]] <- ggplot(data = dt_maxnpcoa) +
  geom_polygon(data = convhull_cam[[2]], alpha = 0.3, color = 'grey50',
               aes(x = Axis.3, y = Axis.4, fill = site_ID)) +
  geom_point(aes(x = Axis.3, y = Axis.4, fill = site_ID), size = 3, shape = 21) +
  looks + labs(x = "PCo3", y = "PCo4")+
  scale_fill_cherulean(palette = "gomphosus", discrete = T, name = 'Site')

(f_maxncam[[1]]|f_maxncam[[2]]) + plot_layout(guides = "collect")

###############################################
# PAIRWISE ABUNDANCE COMPARISONS --------------------------------------------
# look at scatterplots first and then maybe fit OLS/GLMs

# make a dataframe of all methods joined
# each row/record an observed species-size class

dt_allmethods <- data_maxn %>% group_by(site_ID, Taxon, Size_class) %>% summarise(MaxN = sum(MaxN)) %>% 
  full_join(., data_rest[-4], by = c('site_ID', 'Taxon', 'Size_class')) %>% 
  rename(., REST = D)
dt_allmethods <- data_uvc %>% group_by(site_ID, Taxon, Size_class) %>% 
  summarise(UVC = sum(count)) %>% ungroup() %>% filter(str_detect(site_ID, 'kaku|kinalea|sunset')) %>% 
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

dt_all_size$midSize <- str_extract(dt_all_size$Size_class, '(?<=\\_)[:digit:]+$') %>% 
  as.numeric
dt_all_size$midSize[dt_all_size$midSize > 5] <- dt_all_size$midSize[dt_all_size$midSize > 5] -4
dt_all_size$midSize[dt_all_size$midSize == 5] <- dt_all_size$midSize[dt_all_size$midSize == 5] - 2.5

ggplot(dt_all_size) +
  looks + scale_fill_cherulean(palette = "cheridis", discrete = T) +
  geom_col(aes(x = midSize, y = n, fill = method), color = 'transparent', position = 'dodge') +
  labs(x = "Size class (cm)", y = 'Relative frequency') +
  facet_wrap(vars(site_ID))

pairs(dt_allmethods[,4:6], log = 'xy', pch = 1:7, cex = 1.2)
# kinalea
pairs(dt_allmethods[dt_allmethods$site_ID == 'hale_kinalea', 4:6], log = 'xy', pch = 1:7, cex = 1.2)

# exclusive species check
# calculate the "absences" by method
dt_all_long %>% group_by(method, site_ID, Taxon) %>% summarise(n = sum(n)) %>% # aggregate size classes
  ungroup() %>% group_by(method, site_ID) %>% 
  summarise(absences = length(which(n == 0)))
n_distinct(dt_allmethods$Taxon)
# UVC exclusives
temp <- dt_allmethods %>% group_by(site_ID, Taxon) %>% summarise(UVC = sum(UVC), REST = sum(REST), MaxN = sum(MaxN))
temp %>% filter(REST == 0, MaxN == 0, UVC > 0) %>% pull(Taxon) %>% n_distinct
# video/MaxN exclusives
# MaxN + REST presence/absences are the same since they're nested data.
temp %>% filter(UVC == 0) %>% pull(Taxon) %>% n_distinct
# overlaps
temp %>% filter(UVC > 0, (REST > 0 | MaxN > 0)) %>% 
  summarise(overlaps = n()) # video vs UVC
temp %>% ungroup() %>% filter(UVC > 0, MaxN > 0) %>% 
 pull(Taxon) %>% n_distinct # video vs UVC
temp %>% filter(UVC > 0, REST > 0, MaxN == 0) %>% 
  ungroup(Taxon) %>% summarise(overlaps = n()) # REST vs UVC
temp %>% filter(MaxN > 0, REST == 0, UVC > 0) %>% View # maxN + UVC
temp %>% filter(MaxN > 0, REST > 0, UVC > 0) %>% 
  ungroup(Taxon) %>% summarise(overlaps = n()) # all
temp %>% filter(MaxN > 0, REST > 0, UVC > 0) %>% 
  ungroup(Taxon) %>% pull(Taxon) %>% n_distinct # all
rm(temp)
# Species richness per method
dt_all_long %>% group_by(site_ID, method) %>% filter(n > 0) %>% summarise(S = length(unique(Taxon)))


###############################################
## Fit GLMs --------------------------------------------

set.seed(240)
# initially started out thinking about fitting Size_class as an effect, 
# but the size spectra makes me think that letting size absorb variation is not the way to go

model_methodpairs <- as.list(rep('', 3))
# MaxN ~ REST
model_methodpairs[[1]] <- lm(log(MaxN + 1) ~ log(REST + 1), data = dt_allmethods)
# UVC ~ REST
model_methodpairs[[2]] <- lm(log(UVC + 1) ~ log(REST + 1), data = dt_allmethods)
# MaxN ~ UVC
model_methodpairs[[3]] <- lm(log(MaxN + 1) ~ log(UVC + 1), data = dt_allmethods)

summary(model_methodpairs[[1]])
summary(model_methodpairs[[2]])
summary(model_methodpairs[[3]])

# model prediction
# make a standard new x that will be used by all models, including size_class predictor data
newx = data.frame(REST = seq(0, max(dt_allmethods$REST), length.out = 40),
                  MaxN = seq(0, max(dt_allmethods$MaxN), length.out = 40),
                  UVC = seq(0, max(dt_allmethods$UVC), length.out = 40))
predicts <- as.list(rep('', 3))
for (i in 1:3) { # each model
    predicts[[i]] <- predict(model_methodpairs[[i]], newdata = newx, se.fit = T)
    predicts[[i]] <- with(predicts[[i]], data.frame(fitted = fit, lwr = fit - 1.96 * se.fit, 
                                                    upr = fit + 1.96 * se.fit))
}
predicts[[1]]$x = newx$REST
predicts[[2]]$x = newx$REST
predicts[[3]]$x = newx$UVC

f_modelpairs <- as.list(rep('', 3))
# REST - MaxN
f_modelpairs[[1]] <- ggplot(dt_allmethods) +
  geom_jitter(aes(y = log(MaxN + 1), x = log(REST + 1), fill = Size_class), width = 0.05, height = 0.05, 
               alpha = 0.6, size = 2, shape = 21) +
  geom_ribbon(data = predicts[[1]], 
              aes(x = log(x + 1), ymin = lwr, ymax = upr), linetype = 'dashed', color = 'black', linewidth= 0.5, fill = 'transparent') +
  geom_line(data = predicts[[1]], aes(x = log(x + 1), y = fitted)) +
  labs(x = 'log(REST)', y = 'log(MaxN)') +
  scale_fill_cherulean(palette = "cheridis", discrete = T, name = "Size class (cm)")

# REST - UVC
f_modelpairs[[2]] <- ggplot(dt_allmethods) +
  geom_jitter(aes(y = log(UVC + 1), x = log(REST + 1), fill = Size_class), width = 0.05, height = 0.05, 
             alpha = 0.6, size = 2, shape =21) +
  geom_ribbon(data = predicts[[2]], aes(x = log(x + 1), ymin = lwr, ymax = upr), 
              linetype = 'dashed', color = 'black', linewidth= 0.5, fill = 'transparent') +
  geom_line(data = predicts[[2]], aes(x = log(x + 1), y = fitted)) +
  labs(x = 'log(REST)', y = 'log(Point count)')+
  scale_fill_cherulean(palette = "cheridis", discrete = T, name = "Size class (cm)")

# UVC - MaxN
f_modelpairs[[3]] <- ggplot(dt_allmethods) +
  geom_jitter(aes(y = log(MaxN + 1), x = log(UVC + 1), fill = Size_class), width = 0.05, height = 0.05, 
               alpha = 0.6, size = 2, shape = 21) +
  geom_ribbon(data = predicts[[3]], 
              aes(x = log(x + 1), ymin = lwr, ymax = upr), linetype = 'dashed', color = 'black', linewidth= 0.5, fill = 'transparent') +
  geom_line(data = predicts[[3]], aes(x = log(x + 1), y = fitted)) +
  labs(x = 'log(Point count)', y = 'log(MaxN)')+
  scale_fill_cherulean(palette = "cheridis", discrete = T, name = "Size class (cm)")

(f_modelpairs[[1]] + f_modelpairs[[2]] + f_modelpairs[[3]]) * looks + plot_layout(guides = "collect")
