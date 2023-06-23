

## FishN: comparing fish assemblage abundance surveying methods
# visualise density densities of REST poisson fishes

# object requirements:
# data_ruv
# output_poisson
# output_nb
# and their bootstrap objects

require(dplyr)
require(ggplot2)
require(patchwork)
require(lubridate)
require(ggthemes)
require(ggdist)
source('https://gist.githubusercontent.com/cherfychow/e9ae890fd16f4c86730748c067feee2b/raw/899dcfc1745421cb4e6ba26826b0bfe55fd8ec14/cherulean.R')

# make sure size class is read as an ordered factor
data_ruv$Size_class <- factor(data_ruv$Size_class, levels = c('_5', '5_9', '10_19', 
                                                              '20_29', '30_39', '40_49', '50_59'), ordered = T)
# some global pars
looks = theme_classic(base_size = 13) + theme(axis.line = element_line(linewidth = 0.3),
                                              axis.ticks = element_line(linewidth = 0.3))
sites = 3

# Raw data ----------------------------------------------------------------


# visualise the distributions of staying time

stay <- as.list(rep('', sites))
for (i in 1:sites) {
  stay[[i]] <- ggplot(data = data_ruv %>% filter(site_ID == unique(data_ruv$site_ID)[i])) +
    geom_density(aes(x = staytime, fill = spsize), color = "white", alpha = 0.4, linewidth = 0.4) +
    labs(x = "Staying time (s)", y = "Frequency", subtitle = unique(data_ruv$site_ID)[i]) +
    guides(fill = F) + scale_x_log10() + looks +
    scale_fill_cherulean(palette = 'cheridis', discrete = T)
}

stay[[1]] + stay[[2]] + stay[[3]]

# distributions of detections

detect <- as.list(rep('', sites))
for (i in 1:sites) {
  detect[[i]] <- ggplot(data = data_ruv %>% filter(site_ID == unique(data_ruv$site_ID)[i]) %>% 
                          group_by(spsize, Camera) %>% summarise(detects = sum(Count))) +
    geom_density(aes(x = detects, fill = spsize), color = "white", alpha = 0.4) +
    labs(x = "Number of detections", y = "Frequency", subtitle = unique(data_ruv$site_ID)[i]) +
    guides(fill = F) + looks +
    scale_fill_cherulean(palette = 'cheridis', discrete = T) + scale_x_log10()
}

detect[[1]] + detect[[2]] + detect[[3]]
(detect[[1]] + detect[[2]] + detect[[3]]) * scale_x_continuous()



# REST estimates ---------------------------------------------------------

# plot the bootstrapped estimates of D
bootD <- data.frame(Taxon = output_poisson$Taxon[1],
                    site_ID = output_poisson$site_ID[1],
                    Size_class = output_poisson$Size_class[1],
                    boots = boot_p[[1]][,3])
for (i in 2:length(boot_p)) {
  # extract the bootstrapped densities
  temp <- data.frame(Taxon = output_poisson$Taxon[i],
                     site_ID = output_poisson$site_ID[i],
                     Size_class = output_poisson$Size_class[i],
                     boots = boot_p[[i]][,3])
  bootD <- rbind(bootD, temp)
} # collapse bootstrap runs for each species-size into a long format dataframe
rm(temp)

densities <- ggplot() +
  # stat_slab(data = bootD %>% filter(site_ID == 'hale_kaku'), 
  #             aes(y = Taxon, x = boots, fill = Size_class), height = 6, alpha = 0.5) +
  geom_pointrange(data = output_poisson, 
                  aes(y = Taxon, x = D, color = Size_class, xmin = lwr, xmax = upr), linewidth = 1) +
  labs(x = "Bootstrapped density", y = "Species") +
  theme(legend.position = "bottom") + theme_linedraw(base_size = 13) + 
  scale_colour_cherulean(palette = 'cheridis', discrete = T) +
  facet_grid(cols = vars(site_ID)) + scale_x_log10(breaks = c(0, 0.1, 1, 10, 100),
                                                   labels = c(0, 0.1, 1, 10, 100))
densities

# raincloud of modelled densities

# Hale Kaku
# D1 <- ggplot(data = output_poisson %>% filter(site_ID == 'hale_kaku')) +
#   geom_pointrange(aes(y = Taxon, x = D, color = Size_class, xmin = lwr, xmax = upr), linewidth = 1) +
#   theme_bw() + labs(x = "Estimated density", y = "Species", subtitle = "Hale Kaku") +
#   theme(legend.position = "bottom")
# 
# D2 <- ggplot(data = output_poisson %>% filter(site_ID == 'hale_kinalea')) +
#   geom_pointrange(aes(y = Taxon, x = D, color = Size_class, xmin = lwr, xmax = upr), linewidth = 1.5) +
#   theme_bw() + labs(x = "Estimated density", y = "Species", subtitle = "Hale Kinalea") + 
#   theme(legend.position = "bottom")


# Compare method abundances -----------------------------------------------

# make a binned version for UVC
summary(data_uvc)
data_uvc$Size_class[data_uvc$TL_cm <= 5] <- "_5"
data_uvc$Size_class[data_uvc$TL_cm > 5 & data_uvc$TL_cm <= 10] <- "5_10"
data_uvc$Size_class[data_uvc$TL_cm > 10 & data_uvc$TL_cm <= 20] <- "10_20"
data_uvc$Size_class[data_uvc$TL_cm > 20 & data_uvc$TL_cm <= 30] <- "20_30"
data_uvc$Size_class[data_uvc$TL_cm > 30 & data_uvc$TL_cm <= 40] <- "30_40"
data_uvc$Size_class[data_uvc$TL_cm > 40 & data_uvc$TL_cm <= 50] <- "40_50"
data_uvc$Size_class[data_uvc$TL_cm > 50 & data_uvc$TL_cm <= 60] <- "50_60"

# construct the community matrix
allN <- data_uvc %>% filter(str_detect(site_ID, 'kinalea|kaku')) %>% 
  group_by(site_ID, Taxon) %>% summarise(n = sum(count)) %>% mutate(method = 'uvc') %>% ungroup # filter and aggregate UVC
allN <- data_rest %>% group_by(site_ID, Taxon) %>% 
  summarise(n = sum(D, na.rm = T)) %>% mutate(method = 'rest') %>% ungroup %>% filter(!is.na(n), !n == 0) %>% 
  bind_rows(allN, .)
allN <- data_rest %>% group_by(site_ID, Taxon) %>% 
  summarise(n = sum(combinedD, na.rm = T)) %>% mutate(method = 'restmax') %>% ungroup %>% filter(!is.na(n)) %>% 
  bind_rows(allN, .)
allN <- data_maxn %>% group_by(site_ID, Taxon) %>% summarise(n = sum(MaxN, na.rm = T)) %>% ungroup %>% 
  mutate(method = 'maxN') %>% bind_rows(allN, .)

# calculate proportional abundances
allN$rel_n <- NA
allN_1 <- allN %>% filter(site_ID == 'hale_kaku')
allN_2 <- allN %>% filter(site_ID == 'hale_kinalea')
for (i in 1:n_distinct(allN$method)) {
  allN_1$rel_n[allN_1$method == unique(allN_1$method)[i]] <- allN_1$n[allN_1$method == unique(allN_1$method)[i]]/sum(allN_1$n[allN_1$method == unique(allN_1$method)[i]])
  allN_2$rel_n[allN_2$method == unique(allN_2$method)[i]] <- allN_2$n[allN_2$method == unique(allN_2$method)[i]]/sum(allN_2$n[allN_2$method == unique(allN_2$method)[i]])
}

allN <- bind_rows(allN_1, allN_2)


N_1 <- ggplot() +
  geom_point(data = allN %>% filter(site_ID == 'hale_kaku'), 
             aes(y = Taxon, x = n, shape = method, color = method), size = 3) +
  theme_bw() + labs(x = "Abundance", y = "Species", subtitle = "Hale Kaku")

N_2 <- ggplot() +
  geom_point(data = allN %>% filter(site_ID == 'hale_kinalea'), 
             aes(y = Taxon, x = n, shape = method, color = method), size = 3) +
  theme_bw() + labs(x = "Abundance", y = "Species", subtitle = "Hale Kinalea")

(N_1 + N_2) + plot_layout(guides = "collect")


rN_1 <- ggplot() +
  geom_point(data = allN %>% filter(site_ID == 'hale_kaku'), 
             aes(y = Taxon, x = rel_n, shape = method, color = method), size = 3) +
  theme_bw() + labs(x = "Proportional abundance", y = "Species", subtitle = "Hale Kaku")

rN_2 <- ggplot() +
  geom_point(data = allN %>% filter(site_ID == 'hale_kinalea'), 
             aes(y = Taxon, x = rel_n, shape = method, color = method), size = 3) +
  theme_bw() + labs(x = "Proportional abundance", y = "Species", subtitle = "Hale Kinalea")

(rN_1 + rN_2) + plot_layout(guides = "collect")
(rN_1 + rN_2) * scale_x_log10() + plot_layout(guides = "collect")
