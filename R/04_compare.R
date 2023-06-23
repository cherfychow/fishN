

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
source('https://gist.githubusercontent.com/cherfychow/e9ae890fd16f4c86730748c067feee2b/raw/bb52c82ecea1ebf0e87340b37d9b86aead7debb6/cherulean.R')

###############################################
# PCOA COMPARISON -----------------------------------------------------------------

data_uvc <- read.csv('../data/uvc_himb.csv', header = T)
# prep REST to match the assemblage survey data like RUV or MaxN

data_rest <- data_ruv %>% distinct(site_ID, Taxon, Size_class, spsize) %>% 
  left_join(., output_poisson[c(1,3,4,12)], by=c("site_ID", "Taxon", "Size_class"))

# data_rest$combinedD <- data_rest$D # use MaxN to fill in gaps of species that couldn't be fit with REST
# data_rest$combinedD[is.na(data_rest$combinedD)] <- data_rest$MaxN[is.na(data_rest$combinedD)]

# cross-method comparison doesn't have to deal with size classes
# construct the community matrix
nmatrix <- data_uvc %>% filter(str_detect(site_ID, 'kinalea|kaku|sunset')) %>% 
  group_by(site_ID, Taxon) %>% summarise(n = sum(count)) %>% mutate(method = 'uvc') %>% ungroup # filter and aggregate UVC
nmatrix <- data_rest %>% group_by(site_ID, Taxon) %>% 
  summarise(n = sum(D, na.rm = T)) %>% mutate(method = 'rest') %>% ungroup %>% filter(!is.na(n), !n == 0) %>% 
  bind_rows(nmatrix, .)
# nmatrix <- data_rest %>% group_by(site_ID, Taxon) %>% 
#   summarise(n = sum(combinedD, na.rm = T)) %>% mutate(method = 'restmax') %>% ungroup %>% filter(!is.na(n)) %>% 
#   bind_rows(nmatrix, .)
nmatrix <- data_maxn %>% group_by(site_ID, Taxon) %>% summarise(n = sum(MaxN, na.rm = T)) %>% ungroup %>% 
  mutate(method = 'maxN') %>% bind_rows(nmatrix, .)
nmatrix <- nmatrix %>% tidyr::pivot_wider(., names_from = Taxon, values_from = n)
nmatrix <- nmatrix %>% replace(is.na(.), 0) # replace NAs with 0

dis <- vegdist(nmatrix[-(1:2)], method = "euclidean", diag = F, binary = F)
method_pcoa <- ape::pcoa(dis, correction = "none", rn = paste(nmatrix$site_ID, nmatrix$method, sep="_"))
biplot(method_pcoa)

## Plot PCoA with symbology -----------------------------------------------------------------

# scree plot
ggplot(data=method_pcoa$values[1:7,], aes(x=1:7, y=Relative_eig/sum(method_pcoa$values$Relative_eig))) +
  geom_line() + geom_point(shape=21, fill='white', size=3) + labs(x=NULL, y=NULL) +
  scale_x_continuous(breaks=c(1:10)) +
  theme_bw(base_size = 13) + labs(title = "Scree plot (PCoA)")

pcoadt <- as.data.frame(method_pcoa$vectors) %>% bind_cols(., nmatrix[1:2])

pcoa1 <- ggplot(data = pcoadt) +
  geom_point(aes(x = Axis.1, y = Axis.2, color = site_ID, shape = method), size = 3) +
  theme_bw(base_size=13) + labs(x = "PCo1", y = "PCo2") +
  scale_color_cherulean(palette = "cheridis", discrete = T)

pcoa2 <- ggplot(data = pcoadt) +
  geom_point(aes(x = Axis.3, y = Axis.4, color = site_ID, shape = method), size = 3) +
  theme_bw(base_size=13) + labs(x = "PCo3", y = "PCo4") +
  scale_color_cherulean(palette = "cheridis", discrete = T)

(pcoa1|pcoa2) + plot_layout(guides = "collect")
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
relpcoa1 <- ggplot(data = relpcoadt) +
  geom_point(aes(x = Axis.1, y = Axis.2, color = site_ID, shape = method), size = 3) +
  theme_bw(base_size=13) + labs(x = "PCo1", y = "PCo2", subtitle = "Scaled abundances")+
  scale_color_cherulean(palette = "cheridis", discrete = T)

relpcoa2 <- ggplot(data = relpcoadt) +
  geom_point(aes(x = Axis.3, y = Axis.4, color = site_ID, shape = method), size = 3) +
  theme_bw(base_size=13) + labs(x = "PCo3", y = "PCo4")+
  scale_color_cherulean(palette = "cheridis", discrete = T)

(relpcoa1|relpcoa2) + plot_layout(guides = "collect")

###############################################
# PAIRWISE ABUNDANCE COMPARISONS --------------------------------------------
# look at scatterplots first and then maybe fit OLS/GLMs

# make a dataframe of all methods joined
# each row/record an observed species-size class

dt_allmethods <- data_maxn %>% select(!c(spsize,Camera)) %>% 
  full_join(., data_rest[-4], by = c('site_ID', 'Taxon', 'Size_class')) %>% 
  rename(., REST = D)
dt_allmethods <- data_uvc %>% group_by(site_ID, Taxon, Size_class) %>% 
  summarise(UVC = sum(count)) %>% ungroup() %>% filter(str_detect(site_ID, 'kaku|kinalea')) %>% 
  full_join(., dt_allmethods, by = c('site_ID', 'Taxon', 'Size_class'))

# replace NAs with 0s
dt_allmethods$UVC[is.na(dt_allmethods$UVC)] <- 0
dt_allmethods$REST[is.na(dt_allmethods$REST)] <- 0
dt_allmethods$MaxN[is.na(dt_allmethods$MaxN)] <- 0

# REST - MaxN
r_m <- ggplot(dt_allmethods) +
  geom_point(aes(x = REST, y = MaxN, fill = Size_class, size = Size_class), shape = 21, alpha = 0.6) +
  theme_bw(base_size = 13)

# REST - UVC
r_u <- ggplot(dt_allmethods) +
  geom_point(aes(x = REST, y = UVC, fill = Size_class, size = Size_class), shape = 21, alpha = 0.6) +
  theme_bw(base_size = 13)

# UVC - MaxN
u_m <- ggplot(dt_allmethods) +
  geom_point(aes(x = UVC, y = MaxN, fill = Size_class, size = Size_class), shape = 21, alpha = 0.6) +
  theme_bw(base_size = 13)

(r_m + r_u + u_m) * scale_x_log10() * scale_y_log10() + plot_layout(guides = "collect")

# kaku
pairs(dt_allmethods[,4:6], log = 'xy', pch = 1:7, cex = 1.2)
# kinalea
pairs(dt_allmethods[dt_allmethods$site_ID == 'hale_kinalea', 4:6], log = 'xy', pch = 1:7, cex = 1.2)

# exclusive species check
dt_all_long <- dt_allmethods %>% 
  tidyr::pivot_longer(., cols = UVC:REST, names_to = "method", values_to = "n")
# calculate the "absences" by method
dt_all_long %>% group_by(method, site_ID, Taxon) %>% summarise(n = sum(n)) %>% # aggregate size classes
  ungroup() %>% group_by(method, site_ID) %>% 
  summarise(exclusive = length(which(n == 0)))
n_distinct(dt_allmethods$Taxon)
# UVC exclusives
temp <- dt_allmethods %>% group_by(site_ID, Taxon) %>% summarise(UVC = sum(UVC), REST = sum(REST), MaxN = sum(MaxN))
temp %>% filter(REST == 0, MaxN == 0, UVC > 0) 
# video/MaxN exclusives
# MaxN + REST presence/absences are the same since they're nested data.
temp %>% filter(UVC == 0) %>% ungroup(Taxon) %>% summarise(exclusives = n())
# overlaps
temp %>% filter(UVC > 0 & (REST > 0 | MaxN > 0)) %>% 
  ungroup(Taxon) %>% summarise(overlaps = n()) # video vs UVC
temp %>% filter(UVC > 0, REST > 0, MaxN == 0) %>% 
  ungroup(Taxon) %>% summarise(overlaps = n()) # REST vs UVC
temp %>% filter(MaxN > 0, REST == 0, UVC > 0) %>% View # maxN + UVC
temp %>% filter(MaxN > 0, REST > 0, UVC > 0) %>% 
  ungroup(Taxon) %>% summarise(overlaps = n()) # all
rm(temp)
# Species richness per method
dt_all_long %>% group_by(site_ID, method) %>% filter(n > 0) %>% summarise(S = length(unique(Taxon)))


###############################################
## Fit GLMs --------------------------------------------

set.seed(240)

# REST ~ MaxN
model_RM <- as.list(rep('', 2))
# gamma for REST's continuous data
model_RM[[1]] <- lm(log(MaxN + 1) ~ log(REST + 1), data = dt_allmethods)
model_RM[[2]] <- lm(log(MaxN + 1) ~ log(REST + 1) + Size_class, data = dt_allmethods)

summary(model_RM[[1]])
summary(model_RM[[2]])
c(AIC(model_RM[[1]]), AIC(model_RM[[2]])) # model 2
# plot(model_RM[[1]])
# plot(model_RM[[2]])

# REST ~ UVC
model_RU <- as.list(rep('', 2))
# gamma for REST's continuous data
model_RU[[1]] <- lm(log(UVC + 1) ~ log(REST + 1), data = dt_allmethods)
model_RU[[2]] <- lm(log(UVC + 1) ~ log(REST + 1) + Size_class, data = dt_allmethods)
summary(model_RU[[1]])
summary(model_RU[[2]])
c(AIC(model_RU[[1]]), AIC(model_RU[[2]])) # model 2
# plot(model_RU[[1]])
# plot(model_RU[[2]])

# MaxN ~ UVC
model_MU <- as.list(rep('', 2))
# gamma for REST's continuous data
model_MU[[1]] <- lm(log(MaxN + 1) ~ log(UVC + 1), data = dt_allmethods)
model_MU[[2]] <- lm(log(MaxN + 1) ~ log(UVC + 1) + Size_class, data = dt_allmethods)
summary(model_MU[[1]])
summary(model_MU[[2]])
c(AIC(model_MU[[1]]), AIC(model_MU[[2]])) # model 2
# plot(model_MU[[1]])
# plot(model_MU[[2]])
selectedmodels <- list(model_RM[[2]], model_RU[[1]], model_MU[[2]])
rm(model_RM, model_RU, model_MU)

# model prediction
# make a standard new x that will be used by all models, including size_class predictor data
tempx = data.frame(REST = seq(0, max(dt_allmethods$REST), length.out = 20),
                  MaxN = seq(0, max(dt_allmethods$MaxN), length.out = 20),
                  UVC = seq(0, max(dt_allmethods$UVC), length.out = 20))
# iteratively bind rows with the different size classes, only 5-40
newx <- tempx %>% mutate(Size_class = unique(dt_allmethods$Size_class)[1])
newx <- tempx %>% mutate(Size_class = unique(dt_allmethods$Size_class)[2]) %>% bind_rows(., newx) 
newx <- tempx %>% mutate(Size_class = unique(dt_allmethods$Size_class)[3]) %>% bind_rows(., newx) 
newx <- tempx %>% mutate(Size_class = unique(dt_allmethods$Size_class)[6]) %>% bind_rows(., newx) 
rm(tempx)
predicts <- as.list(rep('', 3))
for (i in 1:3) { # each model
    predicts[[i]] <- predict(selectedmodels[[i]], newdata = newx, se.fit = T)
    predicts[[i]] <- with(predicts[[i]], data.frame(fitted = fit, lwr = fit - 1.96 * se.fit, 
                                                    upr = fit + 1.96 * se.fit, Size_class = newx$Size_class))
}
predicts[[1]]$x = newx$REST
predicts[[2]]$x = newx$REST
predicts[[3]]$x = newx$UVC

  
# REST - MaxN
r_m <- ggplot(dt_allmethods) +
  geom_point(aes(y = log(MaxN + 1), x = log(REST + 1)), alpha = 0.6, size = 2, shape = 21) +
  geom_ribbon(data = predicts[[1]] %>% filter(Size_class == '10_20'), 
              aes(x = log(x + 1), ymin = lwr, ymax = upr), linetype = 'dashed', color = 'black', linewidth= 0.5, fill = 'transparent') +
  geom_line(data = predicts[[1]] %>% filter(Size_class == '10_20'), aes(x = log(x + 1), y = fitted)) +
  labs(x = 'log(REST)', y = 'log(MaxN)')

# REST - UVC
r_u <- ggplot(dt_allmethods) +
  geom_point(aes(y = log(UVC + 1), x = log(REST + 1)), alpha = 0.6, size = 2, shape =21) +
  geom_ribbon(data = predicts[[2]], aes(x = log(x + 1), ymin = lwr, ymax = upr), 
              linetype = 'dashed', color = 'black', linewidth= 0.5, fill = 'transparent') +
  geom_line(data = predicts[[2]], aes(x = log(x + 1), y = fitted)) +
  labs(x = 'log(REST)', y = 'log(Point count)')

# UVC - MaxN
u_m <- ggplot(dt_allmethods) +
  geom_point(aes(y = log(MaxN + 1), x = log(UVC + 1)), alpha = 0.6, size = 2, shape = 21) +
  geom_ribbon(data = predicts[[3]] %>% filter(Size_class == '10_20'), 
              aes(x = log(x + 1), ymin = lwr, ymax = upr), linetype = 'dashed', color = 'black', linewidth= 0.5, fill = 'transparent') +
  geom_line(data = predicts[[3]] %>% filter(Size_class == '10_20'), aes(x = log(x + 1), y = fitted)) +
  labs(x = 'log(Point count)', y = 'log(MaxN)')

(r_m + r_u + u_m) * theme_classic(base_size = 13) * theme(panel.grid = element_blank(), axis.line = element_line(linewidth = 0.5)) + plot_layout(guides = "collect")
ggplot(dt_allmethods) +
  geom_point(aes(y = MaxN, x = REST), alpha = 0.6, size = 2) +
  geom_ribbon(data = predicts[[1]] %>% filter(Size_class == '10_20'), 
              aes(x = x, ymin = exp(lwr) - 1, ymax = exp(upr) - 1), linetype = 'dashed', color = 'black', linewidth= 0.5, fill = 'transparent') +
  geom_line(data = predicts[[1]] %>% filter(Size_class == '10_20'), aes(x = x, y = exp(fitted) - 1)) +
  geom_abline(slope = 1, intercept = 0, linetype = "longdash", color = 'blue') +
  labs(x = 'REST', y = 'MaxN')


###############################################
# Camera compositional differences ------------
# PCoA without pooling cameras

# cross-method comparison doesn't have to deal with size classes
# construct the community matrix
n_cam <- data_uvc %>% filter(str_detect(site_ID, 'kinalea|kaku')) %>% 
  group_by(site_ID, Taxon) %>% summarise(n = sum(count)) %>% mutate(method = 'uvc') %>% ungroup # filter and aggregate UVC
n_cam <- data_maxn %>% group_by(site_ID, Taxon) %>% summarise(n = sum(MaxN, na.rm = T)) %>% ungroup %>% 
  mutate(method = 'maxN') %>% bind_rows(nmatrix, .)
n_cam <- nmatrix %>% tidyr::pivot_wider(., names_from = Taxon, values_from = n)
n_cam <- nmatrix %>% replace(is.na(.), 0) # replace NAs with 0

dis <- vegdist(nmatrix[-(1:2)], method = "euclidean", diag = F, binary = F)
method_pcoa <- ape::pcoa(dis, correction = "none", rn = paste(nmatrix$site_ID, nmatrix$method, sep="_"))
method_nmds <- vegan::metaMDS(dis) # really low stress and can't tell between rests... needs metric/quantitative
plot(method_nmds) 
