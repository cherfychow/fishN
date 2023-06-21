
## FishN: comparing fish assemblage abundance surveying methods
# data prep: clean and standardise RUV and UVC fish data
# not meant to be run as part of the analysis pipeline

require(dplyr)
require(stringr)

data_ruv <- read.csv('../data/202111_HIMB_ruv.csv', header = T)
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
data_ruv %>% filter(Taxon == ''|Taxon == ' ') %>% View
sort(unique(data_ruv$Taxon))
sort(unique(word(data_ruv$Taxon, 1))) # just look at genus
# replacement object, 'old' = 'new'
replace <- c('perspicilliatus' = 'perspicillatus',
             'Chaeotodon' = 'Chaetodon',
             'Abudefduf abdominalis x Abudefduf vaigiensis' = 'Abudefduf abdominalis',
             'Acanthururs' = 'Acanthurus',
             'Chaeotdon' = 'Chaetodon')
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
library(worrms)
n_distinct(data_ruv$Taxon) # 50 unique taxa
allsp <- wm_records_names(sort(unique(data_ruv$Taxon)), marine_only = T, fuzzy = F)
allsp <- bind_rows(allsp)
View(allsp %>% filter(status == 'unaccepted')) # which species aren't accepted
# nothing actually wrong with the spellings.

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
summary(data_ruv)
sort(unique(data_ruv$site_ID))
sort(unique(data_uvc$site_ID))
sort(unique(data_ruv$Size_class))
replace = c('10_20' = '10_19',
            '20_30' = '20_29',
            '30_40' = '30_39',
            '40_50' = '40_49',
            '5_10' = '5_9',
            '50_60' = '50_59')
data_ruv$Size_class <- str_replace_all(data_ruv$Size_class, replace)
data_ruv %>% filter(Size_class == ''|is.na(Size_class)) %>% View

sort(unique(data_ruv$Camera))

data_ruv$Count[is.na(data_ruv$Count)] <- 1

write.csv(data_ruv, '../data/ruv_himb_pilot.csv', row.names = F)
write.csv(data_uvc, '../data/uvc_himb.csv', row.names = F)
rm(allsp, replace)
