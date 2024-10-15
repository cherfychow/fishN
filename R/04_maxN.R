


## FishN: comparing fish assemblage abundance surveying methods
# Calculate MaxN

dt_ruv <- dt_ruv %>% select(!c(Observer, VidFile)) %>% mutate(spsize = paste(Taxon, Size_class, sep = "_"))

# turn the timestamps into intervals with a dummy date to make timestamps fake "dates"
dt_ruv$interval <- with(dt_ruv, interval(start = paste0("2022-11-01 12:", Time_entry) %>% ymd_hms(), 
                                             end = paste0("2022-11-01 12:", Time_exit) %>% ymd_hms()))
# just assign a dummy date and hour

# modify camera names so that they're unique across sites
dt_ruv <- dt_ruv %>% mutate(uniqueCam = paste(site_ID, Camera, sep = '_'))
# create a sampling sequence of time at 5 second intervals
time <- seq(as.POSIXct("2022-11-01 12:00:00"), as.POSIXct("2022-11-01 12:45:00"), by=1)
overlap = as.list(rep(NA, n_distinct(dt_ruv$uniqueCam)))
names(overlap) <- unique(dt_ruv$uniqueCam)
# list levels
# site_camera > spsize

# general loop for maxN overlap calculations
for (c in 1:n_distinct(dt_ruv$uniqueCam)) { # for each camera at this site s
  
  # make an overlap list object for each site, fetching by the overlap name variable
  tempsp <- dt_ruv %>% filter(uniqueCam == unique(dt_ruv$uniqueCam)[c]) %>% 
    pull(spsize) %>% unique %>% sort
  overlap[[c]] <- as.list(rep(NA, length(tempsp)))
  
  for (i in 1:length(tempsp)) { # per species-size class in camera c at site s
    overlap[[c]][[i]] <- data.frame(time = time, o1 = NA, o2 = NA, o3 = NA, 
                                    o4 = NA, o5 = NA, o6 = NA, o7 = NA, o8 = NA)
    
    for (t in 1:length(time)) { # for each sampling time stamp
      temp <- NA # make sure start with a clean slate
      
      # per record of the species-size class in raw data
      # in alphabetical order so that the species order is consistent across loops
      for (n in which(dt_ruv$uniqueCam == unique(dt_ruv$uniqueCam)[c] & 
                      dt_ruv$spsize == tempsp[i])) { 
        
        if (time[t] %within% dt_ruv$interval[n]) {
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
MaxN <- as.list(rep(NA, n_distinct(dt_ruv$uniqueCam)))
for (s in 1:n_distinct(dt_ruv$uniqueCam)) {
  
  MaxN[[s]] <- dt_ruv %>% filter(uniqueCam == unique(dt_ruv$uniqueCam)[s]) %>% 
    distinct(site_ID, Camera, Taxon, Size_class, spsize) %>% arrange(spsize) %>% mutate(MaxN = NA) %>% as.data.frame
  
  for (i in 1:length(overlap[[s]])) {
    overlap[[s]][[i]] <- overlap[[s]][[i]] %>% filter(!is.na(o1)) # remove time stamps with no occurrences
    if (nrow(overlap[[s]][[i]]) < 1) next # move on if the whole overlap matrix is null
    overlap[[s]][[i]]$N <-  NA
    for (t in 1:nrow(overlap[[s]][[i]])) { # add up fish N in overlapping
      temp <- overlap[[s]][[i]][t,-1]
      # variable column numbers depending on overlap counts so just everything but the first column (sample time)
      overlap[[s]][[i]][t,'N'] <- sum(dt_ruv$Count[ temp[is.na(temp) == F] ])
    }
    MaxN[[s]]$MaxN[i] <- max(overlap[[s]][[i]]$N) # Max of N
  }
  
}

# the NAs are species that still didn't get sampled by the 5 second intervals
# do those manually

dt_maxn <- bind_rows(MaxN)
dt_maxn$uniqueCam <- with(dt_maxn, paste(site_ID, Camera, sep="_"))
rm(temp)
rm(MaxN)
# rm(overlap) # only if you're sure MaxN loops run properly... 

# # RERUN SPECIAL FOR MISSED SPECIES
# # just for rare or short occurrence species where they didn't get captured by the 5 second interval
# # run an overlap sampling of 1 second intervals

# missedsp <- dt_maxn %>% filter(is.na(MaxN)) # isolate the species we need to target from each site-camera
# time <- seq(as.POSIXct("2022-11-01 12:00:00"), as.POSIXct("2022-11-01 12:45:00"), by=1)
# sites = n_distinct(dt_ruv$site_ID)
# overlap = as.list(rep(NA, nrow(missedsp)))
# # each list item is a species from a unique site-camera in the missedsp subset

# # general loop for maxN overlap calculations
# for (i in 1:nrow(missedsp)) {
  
#   overlap[[i]] <- data.frame(time = time, o1 = NA, o2 = NA, o3 = NA, 
#                                     o4 = NA, o5 = NA, o6 = NA, o7 = NA, o8 = NA)
    
#     for (t in 1:length(time)) { # for each sampling time stamp for missed species i
#       temp <- NA # make sure start with a clean slate

#       for (n in which(dt_ruv$uniqueCam == missedsp$uniqueCam[i] & 
#                       dt_ruv$spsize == missedsp$spsize[i])) { 
        
#         if (time[t] %within% dt_ruv$interval[n]) {
#           temp <- append(temp, n) 
#         } else next
#       }
#       if (length(temp) == 1) next # no overlaps at all, move on
#       temp <- temp[-1]
#       # store in overlap table for species-size i at row t
#       if (length(temp) == 1) { overlap[[i]][t,2] <- temp }
#       else overlap[[i]][t,2:(length(temp)+1)] <- temp 
#     }
  
# }

# beep() # this loop run takes 40 mins for two sites

# # now use the generated overlap tables for calculating MaxN

# for (s in 1:length(overlap)) {
#   for (i in 1:length(overlap[[s]])) {
#     overlap[[s]] <- overlap[[s]] %>% filter(!is.na(o1)) # remove time stamps with no occurrences
#     if (nrow(overlap[[s]]) < 1) next # move on if the whole overlap matrix is null
#     overlap[[s]]$N <-  NA
#     for (t in 1:nrow(overlap[[s]])) { # add up fish N in overlapping
#       temp <- overlap[[s]][t,-1]
#       # variable column numbers depending on overlap counts so just everything but the first column (sample time)
#       overlap[[s]][t,'N'] <- sum(dt_ruv$Count[ temp[is.na(temp) == F] ])
#     }
#     missedsp$MaxN[s] <- max(overlap[[s]]$N) # Max of N
#   }
# }

# dt_maxn <- dt_maxn %>% filter(!is.na(MaxN)) %>% bind_rows(., missedsp) %>% 
#   arrange(site_ID, Camera, Taxon, Size_class)

# export dt_maxn for manual filling of the NAs
write.csv(dt_maxn, './data/dt_MaxN_1s.csv', row.names = F)
# dt_maxn <- read.csv('../data/dt_maxn.csv', header = T) # read back in
rm(overlap, tempsp, time, boots, i, j, n, s, t, missedsp, temp, MaxN)
