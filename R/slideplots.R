## Slide deck version for figures
# replots objects only. relies on saved data in the repo.Rproj

color = '#D74E21'

looks2 = theme_classic(base_size = 13) + theme(panel.grid = element_blank(), 
                                              panel.background = element_rect(fill='transparent'), #transparent
                                              plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
                                              panel.grid.major = element_blank(), #remove major gridlines
                                              panel.grid.minor = element_blank(), #remove minor gridlines
                                              legend.background = element_rect(fill='transparent'), #transparent legend bg
                                              legend.box.background = element_rect(fill='transparent'), #transparent legend panel
                                              line = element_line(color = color),
                                              text = element_text(color = color),
                                              axis.text = element_text(color = color),
                                              axis.ticks = element_line(color = color),
                                              axis.line = element_line(color = color))

# Occurrence minimum
occ <- data_ruv %>% ungroup() %>% group_by(site_ID, spsize) %>% 
  summarise(occurrences = sum(Count)) %>% ungroup

ggplot(data = occ) +
  geom_col(aes(x = reorder(spsize, -occurrences), y = occurrences), fill = color, color = 'transparent') +
  looks2 + labs(x = 'Species-size classes', y = 'Occurrences') +
  geom_hline(yintercept = 5) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank()) +
  scale_y_log10() + facet_wrap(vars(site_ID))

ggsave('occurrences.svg', device = 'svg', width = 250, height = 75, units = 'mm')

# Relative PCoA

ggplot(data = relpcoadt) +
  geom_polygon(aes(x = Axis.1, y = Axis.2, alpha = site_ID), color = color, fill = color) +
  geom_point(aes(x = Axis.1, y = Axis.2, shape = method), size = 3, color = color) +
  looks2 + labs(x = "PCo1", y = "PCo2")+
  scale_shape_manual(values = 21:23, labels = c('MaxN', 'REST', 'UVC'), name = "Method") +
  scale_alpha_manual(values = c(0.2, 0.5, 0.8))
ggsave('pcoa.svg', device = 'svg', width = 150, height = 90, units = 'mm')

# simulated accumulation
ggplot() +
  #geom_point(data = bind_rows(dt_SAC), aes(x = sample, y = nsp, group = interaction(site_ID, method), shape = method), 
  #           alpha = 0.4, size = 3, color = color) +
  # geom_ribbon(data = bind_rows(SAC_pred), aes(x = sample, ymin = lwr, ymax = upr, group = interaction(site_ID, method)), 
  #             alpha = 0.4, fill = color) +
  geom_line(data = bind_rows(SAC_pred), aes(x = sample, y = Median, group = interaction(site_ID, method),
                                            size = method), color = color) +
  facet_wrap(vars(site_ID), scales = 'free_x') + scale_x_log10() + looks2 +
  scale_shape_discrete(solid = F) + scale_size_manual(values = c(0.5,1,1.6)) +
  labs(x = "Number of samples", y = "Number of species")


ggplot() +
  geom_point(data = SAC, aes(x = SACentry, y = spN, group = site_cam, shape = method), 
             alpha = 0.05, color = color) +
  # geom_ribbon(data = bind_rows(SAC_pred2), aes(x = SACentry, ymin = lwr, ymax = upr, group = interaction(site_ID, method),
  #                                              fill = method), alpha = 0.4) +
  geom_line(data = bind_rows(SAC_pred2), aes(x = SACentry, y = fit, group = interaction(site_ID, method),
                                             linetype = method), color = color, size = 1.1) +
  scale_shape_discrete(solid = F) + looks2 + scale_linetype_discrete(name = NULL) +
  facet_wrap(vars(site_ID)) +
  labs(x = "Time elapsed (s)", y = "Number of species") + guides(shape = 'none')

ggsave('sac.svg', device = 'svg', width = 250, height = 75, units = 'mm')


# abundance pairs model

ggplot(data = pair_effects) +
  geom_linerange(aes(xmin = lwr, xmax = upr, y = model), linewidth = 1.1, color = color) +
  geom_point(aes(x = estimate, y = model, shape = par), size = 3, color = color) +
  scale_color_cherulean(palette = "cheridis", discrete = T, name = NULL) + looks2 +
  labs(y = NULL, x = 'Effect estimate')



f_soc <- as.list(rep('', 3))

f_soc[[1]] <- ggplot(dt_abund) +
  ggdist::stat_slab(aes(x = Group, y = log(MaxN + 1)), 
                    adjust = .5, justification = -.5, 
                    color = color, fill = color, size = 0.5, alpha = 0.4) + 
  geom_point(aes(x = Group, y = log(MaxN + 1)), size = 1, alpha = .4, shape = 21, color = color, 
             fill = 'transparent', position = position_jitter(width = .1)) +
  geom_errorbar(data = predicts[[4]], aes(x = Group + 0.25, ymin = lwr, ymax = upr), 
                color = color, width = 0.2) +
  geom_point(data = predicts[[4]], aes(x = Group + 0.25, y = fitted), 
             color = color, size = 2.5) +
  labs(x = 'Sociality', y = 'log(MaxN)')

f_soc[[2]] <- ggplot(dt_abund) +
  ggdist::stat_slab(aes(x = Group, y = log(UVC + 1)), 
                    adjust = .5, justification = -.5, 
                    color = color, fill = color, size = 0.5, alpha = 0.4) + 
  geom_point(aes(x = Group, y = log(UVC + 1)), size = 1, alpha = .4, shape = 21, color = color, 
             fill = 'white', position = position_jitter(width = .1)) +
  geom_errorbar(data = predicts[[5]], aes(x = Group + 0.25, ymin = lwr, ymax = upr), 
                color = color, width = 0.2) +
  geom_point(data = predicts[[5]], aes(x = Group + 0.25, y = fitted), 
             color = color, size = 2.5) +
  labs(x = 'Sociality', y = 'log(UVC)')

f_soc[[3]] <- ggplot(dt_abund) +
  ggdist::stat_slab(aes(x = Group, y = log(MaxN + 1)), 
                    adjust = .5, justification = -.5, 
                    color = color, fill = color, size = 0.5, alpha = 0.4) + 
  geom_point(aes(x = Group, y = log(MaxN + 1)), size = 1, alpha = .4, shape = 21, color = color, 
             fill = 'white', position = position_jitter(width = .1)) +
  geom_errorbar(data = predicts[[6]], aes(x = Group + 0.25, ymin = lwr, ymax = upr), 
                color = color, width = 0.2) +
  geom_point(data = predicts[[6]], aes(x = Group + 0.25, y = fitted), 
             color = color, size = 2.5) +
  labs(x = 'Sociality', y = 'log(MaxN)')

(f_soc[[1]] + f_soc[[2]] + f_soc[[3]]) * looks2 * scale_x_continuous(breaks = c(1:5)) * coord_cartesian(xlim = c(0.5,5.7), ylim = c(0,4.75), clip = 'off') + plot_layout(guides = "collect")
ggsave('social.svg', device = 'svg', width = 250, height = 75, units = 'mm')
