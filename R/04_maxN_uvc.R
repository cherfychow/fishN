

## FishN: comparing fish assemblage abundance surveying methods
# Comparison between MaxN and UVC

require(dplyr)
require(lubridate)
require(beepr)

# Calculate MaxN from the occurrence timestamps ---------------------------

data_ruv$interval <- with(data_ruv, interval(start = paste0("2022-11-01 12:", Time_entry) %>% ymd_hms(), 
                                             end = paste0("2022-11-01 12:", Time_exit) %>% ymd_hms()))
# just assign a dummy date and hour

# create a sampling sequence of time at 15 second intervals
time <- seq(as.POSIXct("2022-11-01 12:00:00"), as.POSIXct("2022-11-01 12:45:00"), by=15)
overlap_kaku <- as.list(rep(NA, data_ruv %>% filter(site_ID == 'hale_kaku') %>% distinct(Taxon, Size_class) %>% nrow))
for (i in 1:length(overlap_kaku)) { # per species-size class
  overlap_kaku[[i]] <- data.frame(time = time, overlap = NA)
  
  for (t in 1:length(time)) { # for each sampling time stamp
    temp <- NA # make sure start with a clean slate
    
    for (n in which(data_ruv$spsize == unique(data_ruv$spsize)[i])) { # per record of the species-size class in raw data
      if (time[t] %within% data_ruv$interval[n]) {
        temp <- append(temp, n) 
      } else next
    }
    if (length(temp) == 1) next # no overlaps at all, move on
    temp <- temp[-1]
    overlap_kaku[[i]][t,2] <- temp # store in overlap table for species-size i at row t
  }
  
}
beep()

# now use the generated overlap tables for calculating MaxN
for (i in 1:length(overlap_kaku)) {
  overlap[[i]] <- overlap[[i]] %>% filter(!is.na(overlap1)) # remove time stamps with no occurrences
  
  for (t in overlap_kaku[[i]]$time) {
    overlap[[i]]$MaxN <- sum(data_ruv$Count[overlap[[i]][na.omit]])
  }
  
}