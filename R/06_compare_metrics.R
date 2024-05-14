

## FishN: comparing fish assemblage abundance surveying methods
# Metrics comparisons

require(rfishbase) # load fishbase
require(iNEXT)
require(SpadeR)
require(tidyverse)
require(patchwork)

# CALCULATE ASSEMBLAGE METRICS -----------------------------------------------------------------

## ASSEMBLAGE STANDING BIOMASS -----------------------------------------------------------

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

((f_metrics[[1]] + f_metrics[[2]] + f_metrics[[3]]) / (f_metrics[[4]] + f_metrics[[5]] + f_metrics[[1]])) + plot_layout(guides = "collect")

rm(hill, expS, dt_S, data_biomass)
