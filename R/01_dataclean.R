
## FishN: comparing fish assemblage abundance surveying methods
# data prep: clean and standardise RUV and UVC fish data
# not meant to be run as part of the analysis pipeline

require(dplyr)
require(stringr)

data_ruv <- read.csv('../data/202111_HIMB_ruvpilot.csv', header = T)
data_uvc <- read.csv('../data/202111_HIMB_uvc5.csv', header = T)

# initial checks
# mostly for data types, NA/0 counts
str(data_ruv)
summary(data_ruv)

str(data_uvc)
summary(data_uvc)

# get rid of the uvc records with no new fish
data_uvc <- data_uvc %>% filter(!count == 0, !is.na(count))

# species spelling checks
sort(unique(data_ruv$Taxon))
sort(unique(word(data_ruv$Taxon, 1))) # just look at genus
# replacement object, 'old' = 'new'
replace <- c('perspicilliatus' = 'perspicillatus',
             'Chaeotodon' = 'Chaetodon')
data_ruv$Taxon <- str_replace_all(data_ruv$Taxon, replace)

# check
sort(unique(data_ruv$Taxon))
sort(unique(word(data_ruv$Taxon, 1))) # just look at genus

sort(unique(data_uvc$Taxon))
sort(unique(word(data_uvc$Taxon, 1))) # just look at genus
# replacement object, 'old' = 'new'
replace <- c('histrix' = 'hystrix',
             'varians' = 'varius',
             'veliferum' = 'velifer',
             'Acathurus' = 'Acanthurus')
data_uvc$Taxon <- str_replace_all(data_uvc$Taxon, replace)

# check
sort(unique(data_uvc$Taxon))
sort(unique(word(data_uvc$Taxon, 1))) # just look at genus

# validate with taxize
library(taxize)

allsp <- c(data_ruv$Taxon, data_uvc$Taxon) %>% sort %>% unique
allsp <- allsp[-str_which(allsp, 'Scarini|sp$|sp[:digit:]$|\\/')]
allsp_wormsid <- get_wormsid(allsp, accepted = T, searchtype = 'scientific', ask = T)
allsp_wormsid <- as_tibble(allsp_wormsid)
allsp_wormsid$taxon <- allsp
View(allsp_wormsid %>% filter(match == 'not found')) # which species aren't accepted

# fix Lutjanus but omit Asterropteryx semipunctata because it's cryptic
data_ruv$Taxon <- str_replace_all(data_ruv$Taxon, 'Lutjanus flavus', 'Lutjanus fulvus')
data_uvc$Taxon <- str_replace_all(data_uvc$Taxon, 'Lutjanus flavus', 'Lutjanus fulvus')
data_uvc <- data_uvc %>% filter(!str_detect(Taxon, 'Asterropteryx semipunctatus'))

sort(unique(data_ruv$Taxon))
sort(unique(data_uvc$Taxon))

# for data_ruv, all blank IP_TP for parrotfish observations are IP unless specified
# fill in IP where it's parrotfish and blank
data_ruv[data_ruv$IP_TP == "" & str_detect(data_ruv$Taxon, '^Chlorurus|^Calotomus|^Scarus'), 'IP_TP'] <- 'IP'
data_ruv[data_ruv$IP_TP == "" & str_detect(data_ruv$Taxon, '^Chlorurus|^Calotomus|^Scarus'), ] # check if any got missed
data_ruv[str_detect(data_ruv$Taxon, '^Chlorurus|^Calotomus|^Scarus'), ] %>% View

# consistency check
sort(unique(data_ruv$site_ID))
sort(unique(data_uvc$site_ID))
sort(unique(data_ruv$Size_class))
sort(unique(data_ruv$Camera))
sort(unique(data_ruv$VidFile))
data_ruv[data_ruv$Size_class == '','Size_class'] <- '5_10'

write.csv(data_ruv, '../data/ruv_himb_pilot.csv', row.names = F)
write.csv(data_uvc, '../data/uvc_himb.csv', row.names = F)
rm(allsp_wormsid, allsp, replace)
