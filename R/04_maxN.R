


## FishN: comparing fish assemblage abundance surveying methods
# Calculate MaxN

data_ruv <- data_ruv %>% select(!c(Observer, VidFile)) %>% mutate(spsize = paste(Taxon, Size_class, sep = "_"))

# turn the timestamps into intervals with a dummy date to make timestamps fake "dates"
data_ruv$interval <- with(data_ruv, interval(start = paste0("2022-11-01 12:", Time_entry) %>% ymd_hms(), 
                                             end = paste0("2022-11-01 12:", Time_exit) %>% ymd_hms()))
# just assign a dummy date and hour

# modify camera names so that they're unique across sites
data_ruv <- data_ruv %>% mutate(uniqueCam = paste(site_ID, Camera, sep = '_'))
# create a sampling sequence of time at 5 second intervals
time <- seq(as.POSIXct("2022-11-01 12:00:00"), as.POSIXct("2022-11-01 12:45:00"), by=5)
overlap = as.list(rep(NA, n_distinct(data_ruv$uniqueCam)))
names(overlap) <- unique(data_ruv$uniqueCam)
# list levels
# site_camera > spsize

# general loop for maxN overlap calculations
for (c in 1:n_distinct(data_ruv$uniqueCam)) { # for each camera at this site s
  
  # make an overlap list object for each site, fetching by the overlap name variable
  tempsp <- data_ruv %>% filter(uniqueCam == unique(data_ruv$uniqueCam)[c]) %>% 
    pull(spsize) %>% unique %>% sort
  overlap[[c]] <- as.list(rep(NA, length(tempsp)))
  
  for (i in 1:length(tempsp)) { # per species-size class in camera c at site s
    overlap[[c]][[i]] <- data.frame(time = time, o1 = NA, o2 = NA, o3 = NA, 
                                    o4 = NA, o5 = NA, o6 = NA, o7 = NA, o8 = NA)
    
    for (t in 1:length(time)) { # for each sampling time stamp
      temp <- NA # make sure start with a clean slate
      
      # per record of the species-size class in raw data
      # in alphabetical order so that the species order is consistent across loops
      for (n in which(data_ruv$uniqueCam == unique(data_ruv$uniqueCam)[c] & 
                      data_ruv$spsize == tempsp[i])) { 
        
        if (time[t] %within% data_ruv$interval[n]) {
          temp <- append(temp, n) 
        } else next
      }
      if (length(temp) == 1) next # no overlaps at all, move on
      temp <- temp[-1]
      # store in overlap table for species-size i at row t
      if (length(temp) == 1) { overlap[[c]][[i]][t,2] <- temp }
      else overlap[[c]][[i]][t,2:(length(temp)+1)] <- temp 
    }
    
  }
  
}

beep() # this loop run takes 40 mins for two sites

# now use the generated overlap tables for calculating MaxN
MaxN <- as.list(rep(NA, n_distinct(data_ruv$uniqueCam)))
for (s in 1:n_distinct(data_ruv$uniqueCam)) {
  
  MaxN[[s]] <- data_ruv %>% filter(uniqueCam == unique(data_ruv$uniqueCam)[s]) %>% 
    distinct(site_ID, Camera, Taxon, Size_class, spsize) %>% arrange(spsize) %>% mutate(MaxN = NA) %>% as.data.frame
  
  for (i in 1:length(overlap[[s]])) {
    overlap[[s]][[i]] <- overlap[[s]][[i]] %>% filter(!is.na(o1)) # remove time stamps with no occurrences
    if (nrow(overlap[[s]][[i]]) < 1) next # move on if the whole overlap matrix is null
    overlap[[s]][[i]]$N <-  NA
    for (t in 1:nrow(overlap[[s]][[i]])) { # add up fish N in overlapping
      temp <- overlap[[s]][[i]][t,-1]
      # variable column numbers depending on overlap counts so just everything but the first column (sample time)
      overlap[[s]][[i]][t,'N'] <- sum(data_ruv$Count[ temp[is.na(temp) == F] ])
    }
    MaxN[[s]]$MaxN[i] <- max(overlap[[s]][[i]]$N) # Max of N
  }
  
}

# the NAs are species that still didn't get sampled by the 5 second intervals
# do those manually

data_maxn <- bind_rows(MaxN)
data_maxn$uniqueCam <- with(data_maxn, paste(site_ID, Camera, sep="_"))
rm(temp)
# rm(overlap) # only if you're sure MaxN loops run properly... 

# RERUN SPECIAL FOR MISSED SPECIES
# just for rare or short occurrence species where they didn't get captured by the 5 second interval
# run an overlap sampling of 1 second intervals

missedsp <- data_maxn %>% filter(is.na(MaxN)) # isolate the species we need to target from each site-camera
time <- seq(as.POSIXct("2022-11-01 12:00:00"), as.POSIXct("2022-11-01 12:45:00"), by=1)
sites = n_distinct(data_ruv$site_ID)
overlap = as.list(rep(NA, nrow(missedsp)))
# each list item is a species from a unique site-camera in the missedsp subset

# general loop for maxN overlap calculations
for (i in 1:nrow(missedsp)) {
  
  overlap[[i]] <- data.frame(time = time, o1 = NA, o2 = NA, o3 = NA, 
                                    o4 = NA, o5 = NA, o6 = NA, o7 = NA, o8 = NA)
    
    for (t in 1:length(time)) { # for each sampling time stamp for missed species i
      temp <- NA # make sure start with a clean slate

      for (n in which(data_ruv$uniqueCam == missedsp$uniqueCam[i] & 
                      data_ruv$spsize == missedsp$spsize[i])) { 
        
        if (time[t] %within% data_ruv$interval[n]) {
          temp <- append(temp, n) 
        } else next
      }
      if (length(temp) == 1) next # no overlaps at all, move on
      temp <- temp[-1]
      # store in overlap table for species-size i at row t
      if (length(temp) == 1) { overlap[[i]][t,2] <- temp }
      else overlap[[i]][t,2:(length(temp)+1)] <- temp 
    }
  
}

beep() # this loop run takes 40 mins for two sites

# now use the generated overlap tables for calculating MaxN

for (s in 1:length(overlap)) {
  for (i in 1:length(overlap[[s]])) {
    overlap[[s]] <- overlap[[s]] %>% filter(!is.na(o1)) # remove time stamps with no occurrences
    if (nrow(overlap[[s]]) < 1) next # move on if the whole overlap matrix is null
    overlap[[s]]$N <-  NA
    for (t in 1:nrow(overlap[[s]])) { # add up fish N in overlapping
      temp <- overlap[[s]][t,-1]
      # variable column numbers depending on overlap counts so just everything but the first column (sample time)
      overlap[[s]][t,'N'] <- sum(data_ruv$Count[ temp[is.na(temp) == F] ])
    }
    missedsp$MaxN[s] <- max(overlap[[s]]$N) # Max of N
  }
}

data_maxn <- data_maxn %>% filter(!is.na(MaxN)) %>% bind_rows(., missedsp) %>% 
  arrange(site_ID, Camera, Taxon, Size_class)

# export data_maxn for manual filling of the NAs
write.csv(data_maxn, '../data/data_MaxN.csv', row.names = F)
# data_maxn <- read.csv('../data/data_maxn.csv', header = T) # read back in
rm(overlap, tempsp, time, boots, i, j, n, s, t, missedsp, temp, MaxN)
