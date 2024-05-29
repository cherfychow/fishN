

## FishN: comparing fish assemblage abundance surveying methods
# visualise density densities of REST poisson fishes

# object requirements:
# data_ruv
# output_rest
# output_nb
# and their bootstrap objects

require(tidyverse)
require(patchwork)
require(ggdist)
source('https://gist.github.com/cherfychow/e9ae890fd16f4c86730748c067feee2b/raw/b2db138ab8164c237129308ea643de78cd54b252/cherulean.R')

# make sure size class is read as an ordered factor
data_ruv$Size_class <- factor(data_ruv$Size_class, levels = c('_5', '5_9', '10_19', 
                                                              '20_29', '30_39', '40_49', '50_59'), ordered = T)
# some global pars
looks <- theme_classic(base_size = 13) + theme(panel.grid = element_blank(), axis.line = element_line(linewidth = .3), axis.ticks = element_line(linewidth = .3))
sites = 3

# Raw data ----------------------------------------------------------------


# visualise the distributions of Duration
temp <- data_ruv %>% group_by(spsize, uniqueCam) %>% summarise(detects = sum(Count)) %>%
  right_join(., data_ruv, by = c('spsize', 'uniqueCam'))
stay <- as.list(rep('', sites))
for (i in 1:sites) {
  stay[[i]] <- ggplot(data = temp %>% filter(site_ID == unique(data_ruv$site_ID)[i])) +
    geom_density(aes(x = staytime, group = spsize), fill = 'deepskyblue2') +
    labs(x = "Duration (s)", y = "Frequency", subtitle = unique(data_ruv$site_ID)[i]) +
    guides(fill = F, color = F) + scale_x_log10() + looks
}

# distributions of detections

detect <- as.list(rep('', sites))
for (i in 1:sites) {
  detect[[i]] <- ggplot(data = data_ruv %>% filter(site_ID == unique(data_ruv$site_ID)[i]) %>% 
                          group_by(spsize, uniqueCam) %>% summarise(detects = sum(Count))) +
    geom_density(aes(x = detects, group = uniqueCam), fill = 'goldenrod') +
    labs(x = "Number of detections", y = "Frequency", subtitle = unique(data_ruv$site_ID)[i]) +
    guides(fill = F, color = F) + scale_x_log10() + looks
}

(detect[[1]] + detect[[2]] + detect[[3]])/(stay[[1]] + stay[[2]] + stay[[3]])


ggplot(data = temp) +
  geom_point(aes(x = detects, y = staytime), shape = 21) + looks +
  scale_x_log10() + scale_y_log10()


# REST estimates ---------------------------------------------------------

pretty <- c('^\\_' = '\\< ', '(?<=[:digit:])\\_' = '\\-')
output_rest$Size_class <- str_replace_all(output_rest$Size_class, pretty) %>% 
  factor(., levels = c('< 5', '5-9', '10-19', 
                       '20-29', '30-39'), ordered = T)

densities <- ggplot() +
  # stat_slab(data = bootD %>% filter(site_ID == 'hale_kaku'), 
  #             aes(y = Taxon, x = boots, fill = Size_class), height = 6, alpha = 0.5) +
  geom_errorbarh(data = output_rest, 
                  aes(y = Taxon, color = Size_class, xmin = lwr, xmax = upr), linewidth = 0.75, height = 0.3, position = position_dodge2(width = 0.5)) +
  geom_point(data = output_rest, 
                aes(y = Taxon, x = D, color = Size_class), size = 2, position = position_dodge2(width = 0.5)) +
  labs(x = "Bootstrapped density", y = "") +
  scale_colour_cherulean(palette = 'cheridis', discrete = T, name = "Size class (cm)") +
  facet_grid(cols = vars(site_ID)) + scale_x_log10(labels = function(x) format(x, scientific = FALSE)) + 
  theme_classic(base_size = 13) + 
  theme(axis.text.y = element_text(face = 'italic'), 
        panel.grid.major.x = element_line(linewidth = .3, colour = 'grey80'),
        panel.grid.minor.x = element_line(linewidth = .3, colour = 'grey80'),
        axis.line = element_line(linewidth = .3), 
        axis.ticks = element_line(linewidth = .3),
        legend.position = 'bottom')
densities

# REST no size estimates

nosize <- read_csv('../outputs/REST_nosize.csv')

ggplot() +
  # stat_slab(data = bootD %>% filter(site_ID == 'hale_kaku'), 
  #             aes(y = Taxon, x = boots, fill = Size_class), height = 6, alpha = 0.5) +
  geom_errorbarh(data = nosize, 
                 aes(y = Taxon, xmin = lwr, xmax = upr), linewidth = 0.75, height = 0.3, position = position_dodge2(width = 0.5), color = 'grey20') +
  geom_point(data = nosize, 
             aes(y = Taxon, x = D), size = 2, position = position_dodge2(width = 0.5), color = 'grey20') +
  labs(x = "Bootstrapped density", y = "") +
  facet_grid(cols = vars(site_ID)) + scale_x_log10(labels = function(x) format(x, scientific = FALSE)) + 
  theme_classic(base_size = 13) + 
  theme(axis.text.y = element_text(face = 'italic'), 
        panel.grid.major.x = element_line(linewidth = .3, colour = 'grey80'),
        panel.grid.minor.x = element_line(linewidth = .3, colour = 'grey80'),
        axis.line = element_line(linewidth = .3), 
        axis.ticks = element_line(linewidth = .3),
        legend.position = 'bottom')

