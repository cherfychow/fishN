

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

data_ruv$spsize <- with(data_ruv, paste(Taxon, Size_class, sep = "_"))


# Raw data ----------------------------------------------------------------

# visualise the distributions of staying time
stay1 <- ggplot(data = data_ruv %>% filter(site_ID == 'hale_kaku')) +
  geom_density(aes(x = staytime, fill = spsize), color = "white", alpha = 0.4) +
  labs(x = "Staying time (s)", y = "Frequency", subtitle = "Hale Kaku") +
  guides(fill = F) + scale_x_log10() + theme_bw()

stay2 <- ggplot(data = data_ruv %>% filter(site_ID == 'hale_kinalea')) +
  geom_density(aes(x = staytime, fill = spsize), color = "white", alpha = 0.4) +
  labs(x = "Staying time (s)", y = "Frequency", subtitle = "Hale Kinalea") +
  guides(fill = F) + scale_x_log10() + theme_bw()

stay1 + stay2

# distributions of detections
detect1 <- ggplot(data = data_ruv %>% filter(site_ID == 'hale_kaku') %>% group_by(spsize, Camera) %>% summarise(detects = sum(Count))) +
  geom_density(aes(x = detects, fill = spsize), color = "white", alpha = 0.4) +
  labs(x = "Number of detections", y = "Frequency", subtitle = "Hale Kaku") +
  guides(fill = F) + theme_bw() # + scale_x_log10()

detect2 <- ggplot(data = data_ruv %>% filter(site_ID == 'hale_kinalea') %>% group_by(spsize, Camera) %>% summarise(detects = sum(Count))) +
  geom_density(aes(x = detects, fill = spsize), color = "white", alpha = 0.4) +
  labs(x = "Number of detections", y = "Frequency", subtitle = "Hale Kinalea") +
  guides(fill = F) + theme_bw() # + scale_x_log10()

detect1 + detect2


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

D1 <- ggplot() +
  stat_slab(data = bootD %>% filter(site_ID == 'hale_kaku'), 
               aes(y = Taxon, x = boots, fill = Size_class), height = 6, alpha = 0.5) +
  geom_pointrange(data = output_poisson %>% filter(site_ID == 'hale_kaku'), 
                  aes(y = Taxon, x = D, color = Size_class, xmin = lwr, xmax = upr), linewidth = 1) +
  theme_bw() + labs(x = "Bootstrapped density", y = "Species", subtitle = "Hale Kaku") +
  theme(legend.position = "bottom")

D2 <- ggplot() +
  stat_slab(data = bootD %>% filter(site_ID == 'hale_kinalea'), 
            aes(y = Taxon, x = boots, fill = Size_class), height = 6, alpha = 0.5) +
  geom_pointrange(data = output_poisson %>% filter(site_ID == 'hale_kinalea'), 
                  aes(y = Taxon, x = D, color = Size_class, xmin = lwr, xmax = upr), linewidth = 1) +
  theme_bw() + labs(x = "Bootstrapped density", y = "Species", subtitle = "Hale Kinalea") +
  theme(legend.position = "bottom")

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

D1 + D2

(D1 + D2) * coord_cartesian(xlim = c(0,5)) # look at the non-outlier distributions

