

## FishN: comparing fish assemblage abundance surveying methods
# Assemblage species capture comparisons

require(dplyr)
require(lubridate)
require(beepr)

require(vegan)
require(stringr)
require(ape)
require(ggplot2)
require(patchwork)
require(nlstools)
source('https://gist.githubusercontent.com/cherfychow/e9ae890fd16f4c86730748c067feee2b/raw/b2db138ab8164c237129308ea643de78cd54b252/cherulean.R')

set.seed(24)
looks <- theme_classic(base_size = 13) + 
  theme(panel.grid = element_blank(), axis.ticks = element_line(linewidth = .75), axis.line = element_line(linewidth = .75))

# load files generated from previous analyses 
# assumes dt_ruv is already in the environment from REST fitting
# dt_rest <- read.csv('../outputs/REST_output.csv', header = T)
dt_maxn <- read.csv('../data/data_maxn.csv', header = T)
dt_uvc <- read.csv('../data/uvc_himb.csv', header = T)
# dt_rest <- rest_output
# 
# rm(rest_output, rest_sp)

# Prep data structure

dt_uvc$Size_class <- ''
dt_uvc$Size_class[dt_uvc$TL_cm <= 5] <- "_5"
dt_uvc$Size_class[dt_uvc$TL_cm > 5 & dt_uvc$TL_cm <= 10] <- "5_9"
dt_uvc$Size_class[dt_uvc$TL_cm > 10 & dt_uvc$TL_cm <= 20] <- "10_19"
dt_uvc$Size_class[dt_uvc$TL_cm > 20 & dt_uvc$TL_cm <= 30] <- "20_29"
dt_uvc$Size_class[dt_uvc$TL_cm > 30 & dt_uvc$TL_cm <= 40] <- "30_39"
dt_uvc$Size_class[dt_uvc$TL_cm > 40 & dt_uvc$TL_cm <= 50] <- "40_49"
dt_uvc$Size_class[dt_uvc$TL_cm > 50 & dt_uvc$TL_cm <= 60] <- "50_59"

dt_rest <- dt_rest %>% select(site_ID, Taxon, Size_class, pred) %>% 
  mutate(spsize = paste(Taxon, Size_class, sep = "_")) %>% rename(D = pred)

# make a dataframe of all methods joined
# each row/record an observed species-size class

dt_allmethods <- dt_maxn %>% group_by(site_ID, Taxon, Size_class) %>% summarise(MaxN = sum(MaxN)) %>% 
  full_join(., dt_rest[-5], by = c('site_ID', 'Taxon', 'Size_class')) %>% 
  rename(., REST = D)
dt_allmethods <- dt_uvc %>% group_by(site_ID, Taxon, Size_class) %>% 
  summarise(UVC = sum(count)) %>% ungroup() %>% filter(str_detect(site_ID, 'kaku|hinalea|sunset')) %>% 
  full_join(., dt_allmethods, by = c('site_ID', 'Taxon', 'Size_class'))

# replace NAs with 0s
dt_allmethods$UVC[is.na(dt_allmethods$UVC)] <- 0
dt_allmethods$REST[is.na(dt_allmethods$REST)] <- 0
dt_allmethods$MaxN[is.na(dt_allmethods$MaxN)] <- 0

pretty <- c('^\\_' = '\\< ', '(?<=[:digit:])\\_' = '\\-')
dt_allmethods$Size_class <- str_replace_all(dt_allmethods$Size_class, pretty) %>% 
  factor(., levels = c('< 5', '5-9', '10-19', 
                       '20-29', '30-39', '40-49', '50-59'), ordered = T)

# long format
dt_all_long <- dt_allmethods %>% 
  tidyr::pivot_longer(., cols = UVC:REST, names_to = "method", values_to = "n")


# PCOA COMPOSITION COMPARISON -----------------------------------------------------------------

# cross-method comparison doesn't have to deal with size classes
# construct the community matrix
pcoa_nmatrix <- dt_uvc %>% filter(str_detect(site_ID, 'hinalea|kaku|sunset')) %>% 
  group_by(site_ID, Taxon) %>% 
  summarise(n = sum(count)) %>% mutate(method = 'uvc') %>% ungroup # filter and aggregate UVC
pcoa_nmatrix <- dt_rest %>% group_by(site_ID, Taxon) %>% 
  summarise(n = sum(D, na.rm = T)) %>% mutate(method = 'rest') %>% ungroup %>% filter(!is.na(n), !n == 0) %>% 
  bind_rows(pcoa_nmatrix, .)
pcoa_nmatrix <- dt_maxn %>% group_by(site_ID, Taxon) %>% summarise(n = sum(MaxN, na.rm = T)) %>% ungroup %>% 
  mutate(method = 'maxN') %>% bind_rows(pcoa_nmatrix, .)
pcoa_nmatrix <- pcoa_nmatrix %>% tidyr::pivot_wider(., names_from = Taxon, values_from = n)
pcoa_nmatrix <- pcoa_nmatrix %>% replace(is.na(.), 0) # replace NAs with 0

pcoa_dis <- vegdist(pcoa_nmatrix[-(1:2)], method = "euclidean", diag = F, binary = F)
method_pcoa <- ape::pcoa(pcoa_dis, correction = "none", rn = paste(pcoa_nmatrix$site_ID, pcoa_nmatrix$method, sep="_"))
biplot(method_pcoa)

## Plot PCoA with symbology -----------------------------------------------------------------
source('https://raw.githubusercontent.com/cherfychow/FishTraitsCoralRec/main/analysis_code/function_convhull.R')
# convex hull script

# scree plot
ggplot(data=method_pcoa$values[1:7,], aes(x=1:7, y=Relative_eig/sum(method_pcoa$values$Relative_eig))) +
  geom_line() + geom_point(shape=21, fill='white', size=3) + labs(x=NULL, y=NULL) +
  scale_x_continuous(breaks=c(1:10)) +
  theme_bw(base_size = 13) + labs(title = "Scree plot (PCoA)")

dt_pcoa <- as.data.frame(method_pcoa$vectors) %>% bind_cols(., pcoa_nmatrix[1:2])

f_pcoa1 <- ggplot(data = dt_pcoa) +
  geom_polygon(aes(x = Axis.1, y = Axis.2, fill = site_ID, color = site_ID), alpha = 0.5) +
  geom_point(aes(x = Axis.1, y = Axis.2, fill = site_ID, shape = method), size = 3) +
  looks + labs(x = "PCo1", y = "PCo2") +
  scale_color_cherulean(palette = 'cheridis', discrete = T, name = 'Site') +
  scale_fill_cherulean(palette = 'cheridis', discrete = T, name = 'Site') +
  scale_shape_manual(values = 21:23, labels = c('MaxN', 'REST', 'UVC'), name = "Method")

# proportional abundances

pcoa_n_rel <- apply(pcoa_nmatrix[-(1:2)], 1, function(x) x/sum(x)) %>% t %>% as.data.frame
pcoa_n_rel <- bind_cols(pcoa_nmatrix[1:2], pcoa_n_rel)
# try again
pcoa_dis_r <- vegdist(pcoa_n_rel[-(1:2)], method = "euclidean", diag = F, binary = F)
method_pcoa_rel <- ape::pcoa(pcoa_dis_r, correction = "none", rn = paste(pcoa_n_rel$site_ID, pcoa_n_rel$method, sep="_"))

# scree plot
ggplot(data=method_pcoa_rel$values[1:7,], aes(x=1:7, y=Relative_eig/sum(method_pcoa_rel$values$Relative_eig))) +
  geom_line() + geom_point(shape=21, fill='white', size=3) + labs(x=NULL, y=NULL) +
  scale_x_continuous(breaks=c(1:10)) +
  theme_bw(base_size = 13) + labs(title = "Scree plot (PCoA, scaled abundances)")


dt_pcoa_r <- as.data.frame(method_pcoa_rel$vectors) %>% bind_cols(., pcoa_n_rel[1:2])
f_relpcoa1 <- ggplot(data = dt_pcoa_r) +
  geom_polygon(aes(x = Axis.1, y = Axis.2, fill = site_ID, color = site_ID), alpha = 0.5) +
  geom_point(aes(x = Axis.1, y = Axis.2, fill = site_ID, shape = method), size = 3) +
  looks + labs(x = "PCo1", y = "PCo2")+
  scale_fill_cherulean(palette = "cheridis", discrete = T, name = "Site") +
  scale_color_cherulean(palette = "cheridis", discrete = T, name = "Site") +
  scale_shape_manual(values = 21:23, labels = c('MaxN', 'REST', 'UVC'), name = "Method")

rm(pcoa_nmatrix, method_pcoa, method_pcoa_rel, pcoa_n_rel, 
   pcoa_dis_r, pcoa_dis)


# COMPARE ASSEMBLAGES BETWEEN CAMERAS ---------------------------------------
# compare the assemblage composition with PCoA between MaxN cameras

pcoa_nmatrix_cam <- dt_maxn %>% group_by(site_ID, Camera, Taxon) %>% 
  summarise(n = sum(MaxN)) %>% ungroup
pcoa_nmatrix_cam <- pcoa_nmatrix_cam %>% tidyr::pivot_wider(names_from = Taxon, values_from = n) %>% 
  replace(is.na(.), 0) # replace NAs with 0

dis_maxn <- vegdist(pcoa_nmatrix_cam[-(1:2)], method = "euclidean", diag = F, binary = F)
pcoa_maxn <- ape::pcoa(dis_maxn, correction = "none", rn = paste(pcoa_nmatrix_cam$site_ID, pcoa_nmatrix_cam$Camera, sep="_"))
biplot(pcoa_maxn)

# check scree
ggplot(data=pcoa_maxn$values[1:7,], aes(x=1:7, y=Cumul_eig)) +
  geom_line() + geom_point(shape=21, fill='white', size=3) + labs(x=NULL, y=NULL) +
  scale_x_continuous(breaks=c(1:10)) +
  looks + labs(title = "Scree plot")

# 4 dimensions is ok
## Plot PCoA results
dt_maxnpcoa <- as.data.frame(pcoa_maxn$vectors) %>% bind_cols(., pcoa_nmatrix_cam[1:2])
# calculate convex hulls per site
convhull_cam <- list(as.list(rep(NA, 3)), as.list(rep(NA, 3)))
for (i in 1:3) {
  convhull_cam[[1]][[i]] <- convhull.vert(dt_maxnpcoa %>% 
                                            filter(site_ID == unique(site_ID)[i]) %>% 
                                            select(Axis.1, Axis.2)) %>% 
    mutate(site_ID = unique(dt_maxnpcoa$site_ID)[i])
  convhull_cam[[2]][[i]] <- convhull.vert(dt_maxnpcoa %>% 
                                            filter(site_ID == unique(site_ID)[i]) %>% 
                                            select(Axis.3, Axis.4)) %>% 
    mutate(site_ID = unique(dt_maxnpcoa$site_ID)[i])
}
convhull_cam[[1]] <- bind_rows(convhull_cam[[1]])
convhull_cam[[2]] <- bind_rows(convhull_cam[[2]])

f_pcoacam <- as.list(rep(NA,2))
f_pcoacam[[1]] <- ggplot(data = dt_maxnpcoa) +
  geom_polygon(data = convhull_cam[[1]], 
               aes(x = Axis.1, y = Axis.2, fill = site_ID)) +
  geom_point(aes(x = Axis.1, y = Axis.2, fill = site_ID), size = 3, shape = 21) +
  looks + labs(x = "PCo1", y = "PCo2")+
  scale_fill_cherulean(palette = "cheridis", discrete = T, name = 'Site')

f_pcoacam[[2]] <- ggplot(data = dt_maxnpcoa) +
  geom_polygon(data = convhull_cam[[2]],
               aes(x = Axis.3, y = Axis.4, fill = site_ID)) +
  geom_point(aes(x = Axis.3, y = Axis.4, fill = site_ID), size = 3, shape = 21) +
  looks + labs(x = "PCo3", y = "PCo4")+
  scale_fill_cherulean(palette = "cheridis", discrete = T, name = 'Site')


rm(pcoa_maxn, pcoa_nmatrix_cam)

# TALLIES --------------------------------------------
# exclusive species check

# calculate the "absences" by method
dt_all_long %>% group_by(method, site_ID, Taxon) %>% summarise(n = sum(n)) %>% # aggregate size classes
  ungroup() %>% group_by(method, site_ID) %>% 
  summarise(absences = length(which(n == 0))) %>% ungroup %>% group_by(method) %>% 
  summarise(meanMiss = mean(absences))

# not splitting by site
temp <- dt_allmethods %>% group_by(Taxon) %>% summarise(UVC = sum(UVC), REST = sum(REST), MaxN = sum(MaxN))
n_distinct(temp$Taxon)

# EXCLUSIVES
temp %>% filter(MaxN == 0 & UVC > 0) %>% pull(Taxon) %>% n_distinct #UVC
temp %>% filter(UVC == 0 & MaxN > 0) %>% pull(Taxon) %>% n_distinct # video
# overlaps, video uvc
temp %>% filter(MaxN > 0, UVC > 0) %>% pull(Taxon) %>% n_distinct
temp %>% filter(MaxN > 0, UVC > 0) %>% pull(Taxon) %>% n_distinct

temp <- dt_allmethods %>% group_by(Taxon, site_ID) %>% summarise(UVC = sum(UVC), REST = sum(REST), MaxN = sum(MaxN))

# overlaps
temp %>% filter(UVC > 0 & MaxN > 0) %>% group_by(site_ID) %>% 
  summarise(overlaps = n()) # overlaps in video and uvc, by site
temp %>% ungroup() %>% filter(UVC > 0 & MaxN > 0) %>% 
 pull(Taxon) %>% n_distinct # overlap total video-uvc

temp %>% filter(MaxN > 0, REST == 0, UVC > 0) %>% View # maxN + UVC
temp %>% ungroup() %>% filter(MaxN > 0, REST > 0, UVC > 0) %>% nrow # all methods overlap

# Species richness per method
dt_all_long %>% group_by(site_ID, method) %>% filter(n > 0) %>% summarise(S = length(unique(Taxon)))
temp <- dt_allmethods %>% group_by(site_ID, Taxon) %>% summarise(UVC = sum(UVC), REST = sum(REST), MaxN = sum(MaxN))
temp %>% filter(REST == 0, MaxN > 0) %>% summarise(missedsp = n_distinct(Taxon))
rm(temp)

# Species abundance distributions --------------------------------------------
# do it on size aggregated data
temp <- dt_all_long %>% group_by(site_ID, method) %>% summarise(totaln = sum(n))
# add relative abundances
dt_all_long <- left_join(dt_all_long, temp, by = c("site_ID", "method")) %>% 
  mutate(reln = n/totaln) %>% select(!totaln)

ggplot(data = dt_all_long %>% filter(reln > 0) %>% group_by(site_ID, method, Taxon) %>% 
         summarise(reln = sum(reln))) +
  geom_density(aes(x = reln, color = method), linewidth = 0.75) +
  labs(x = 'Abundances', y = "Number of species") + looks +
  scale_color_cherulean(palette = "cheridis", discrete = T, name = 'Method')

f_SAD <- ggplot(data = dt_all_long %>% filter(reln > 0) %>% group_by(site_ID, method, Taxon) %>% 
                summarise(reln = sum(reln))) +
  geom_density(aes(x = reln, color = method), linewidth = 0.5) +
  labs(x = 'Abundances', y = "Species frequency") + looks +
  scale_color_cherulean(palette = "gomphosus", discrete = T, name = 'Method') +
  facet_wrap(vars(site_ID))

# Species rarefaction -----------------------------------------------------------
# simulated accumulation

temp <- dt_all_long %>% group_by(site_ID, Taxon, method) %>% 
  summarise(n = sum(n), reln=sum(reln)) %>% ungroup %>% 
  filter(n > 0) %>% group_split(site_ID, method)

temp2 <- dt_all_long %>% ungroup() %>% filter(n > 0) %>% 
  group_by(site_ID, method) %>% summarise(total = sum(round(n, 0))) # use for sample size, as in total number of individuals


# hard to use timestamps for UVC, so we'll do it by abundance-weighted sampling
dt_rar <- as.list(rep('', 9))
for(i in 1:9) {
  dt_rar[[i]] <- data.frame(site_ID = temp[[i]]$site_ID[1], method = temp[[i]]$method[1],
                            sample = 1:temp2$total[i], nsp = NA) # empty dataframe
  tempsp <- data.frame(sample = NA, sp = NA) # dummy to append species accumulation to
  temp3 <- temp[[i]]$Taxon[rep(1:nrow(temp[[i]]), ceiling(temp[[i]]$n))]
  for (j in 1:temp2$total[i]) {
    tempsp <- rbind(tempsp,
      data.frame(sample = j, sp = sample(temp3, size = j, replace = F)))
    # sample the assemblage using relative abundances as probability
    # temporarily store sampled species
    dt_rar[[i]][j,4] <- tempsp %>% filter(sample <= j, is.na(sp) == F) %>% pull(sp) %>% n_distinct(.) 
    # calculate species accumulative
  } 
  dt_rar[[i]]$site_ID <- temp[[i]]$site_ID[1]
  dt_rar[[i]]$method <- temp[[i]]$method[1]
}

ggplot(data = bind_rows(dt_rar)) +
  geom_line(aes(x = sample, y = nsp, group = interaction(site_ID, method),
                color = site_ID), linewidth = 1) + looks + scale_x_log10()

# now we'll fit accumulation models
rar_fit <- as.list(rep(NA, 9))
rar_pred <- as.list(rep(NA, 9))
set.seed(240)
for (i in 1:9) { # beta P model
  optout <- tryCatch(nls(data = na.omit(dt_rar[[i]]),
                         formula = nsp ~ a - ((a-b)*exp(-c*sample)),
                         start = list(a = 10, b = 1, c = 0.1)),
                     error = function(e) e)
  if(inherits(optout, "error")) next # if an error message gets generated this run, move to the next iteration of the loop
  rar_fit[[i]] <- optout
}

rar_fit[[1]] <- nls(data = na.omit(dt_rar[[1]]),
    formula = nsp ~ a - ((a-b)*exp(-c*sample)),
    start = list(a = 30, b = 1, c = 0.1))

require(nlstools)
for (i in 1:9) {
  rar_pred[[i]] <- nlsBootPredict(
    nlsBoot(rar_fit[[i]], niter=500), 
    newdata = dt_rar[[i]][3], interval = "confidence") %>% as.data.frame
  colnames(rar_pred[[i]])[2:3] <- c('lwr', 'upr')
  rar_pred[[i]] <- data.frame(site_ID = temp[[i]]$site_ID[1], method = temp[[i]]$method[1], 
                              sample = dt_rar[[i]][3]) %>% bind_cols(., rar_pred[[i]])
}

rm(temp, temp2, tempsp, optout)

f_rarefaction <- ggplot() +
  geom_point(data = bind_rows(dt_rar), aes(x = sample, y = nsp, group = interaction(site_ID, method), color = method), 
             alpha = 0.2, shape = 21) +
  geom_ribbon(data = bind_rows(rar_pred), aes(x = sample, ymin = lwr, ymax = upr, group = interaction(site_ID, method),
                                            fill = method), alpha = 0.4) +
  geom_line(data = bind_rows(rar_pred), aes(x = sample, y = Median, group = interaction(site_ID, method),
                                            color = method)) +
  scale_color_cherulean(palette = "gomphosus", discrete = T, name = "Method") +
  scale_fill_cherulean(palette = "gomphosus", discrete = T, name = "Method") +
  facet_wrap(vars(site_ID), scales = 'free_x') + looks + scale_x_log10() +
  labs(x = "Number of samples", y = "Number of species") +
  theme(legend.position = 'bottom')

# Species accumulation-----------------------------------------------------------

#  timestamps are by video but we want them to be timestamps relative to total observation time, not video time.
cameras <- dt_ruv %>% distinct(site_ID, Camera, VidFile)
ncam <- cameras %>% group_by(site_ID, Camera) %>% summarise(ncam = n_distinct(VidFile)) %>% pull(ncam)
seq <- ''
  for (i in 1:length(ncam)) {
    seq <- c(seq, 1:ncam[i])
  } # generate a sequence vector to identify the nominal order of video files
seq <- seq[-1] # trim that dummy start

cameras$VidSeq <- seq

dt_ruv <- left_join(dt_ruv, cameras, by=c('site_ID', 'Camera', 'VidFile'))
dt_ruv$VidSeq <- as.numeric(dt_ruv$VidSeq)

# make the timestamps relative to site and not video file
# all videos duration = 11:48 except for deep_0 = 17:42
for (i in 2:5) {
  # loop for each video file in order(1-4 or 5) add a multiple of the video file length
  # 2:5 because the first files don't need fixing
  # for the 11:48 videos
  dt_ruv$SACentry[dt_ruv$VidSeq == i & dt_ruv$Camera != 'deep_0'] <- dt_ruv$entrytime_c[dt_ruv$VidSeq == i & dt_ruv$Camera != 'deep_0'] + (i-1)*(dminutes(11) + dseconds(48))
  # for deep_0 17:42 videos
  dt_ruv$SACentry[dt_ruv$VidSeq == i & dt_ruv$Camera == 'deep_0'] <- dt_ruv$entrytime_c[dt_ruv$VidSeq == i & dt_ruv$Camera == 'deep_0'] + (i-1)*(dminutes(17) + dseconds(42))
}

# just put the time stamps in for first video files
dt_ruv$SACentry[dt_ruv$VidSeq == 1] <- dt_ruv$entrytime_c[dt_ruv$VidSeq == 1]
dt_ruv$VidSeq <- NULL # get rid of these dummy objects
rm(cameras, ncam, seq)

## Calculate cumulative species richness over time
# make a dataframe to populate species accumulation per timestamp
SAC <- dt_ruv %>% mutate(site_cam = paste(site_ID, Camera, sep = "_")) %>% 
  select(site_ID, site_cam, SACentry)  # matches the row index of dt_ruv

# make the SAC times relative to the time of first observation
for (i in 1:11) {
  time1 <- SAC$SACentry[which(SAC$site_cam == unique(SAC$site_cam)[i])[1]]
  for (j in which(SAC$site_cam == unique(SAC$site_cam)[i])) {
    SAC$SACentry[j] <- SAC$SACentry[j] - time1 # take every time stamp and subtract time1 from it
  }
}
rm(time1)

SAC$spN <- NA
# for every row at every site-camera, calculate the species number accumulation through time
for (i in 1:11) {
  rows = seq(which(SAC$site_cam == unique(SAC$site_cam)[i])[1], max(which(SAC$site_cam == unique(SAC$site_cam)[i])), by = 3)
  for (j in rows) {
    if (j == rows[1]) { # the first row from each site represents the first species record, so these will always start at 1
      SAC$spN[j] <- 1
    }
    else {
      SAC$spN[j] <- dt_ruv$Taxon[rows[1]:j] %>% n_distinct # number of unique species from the first row of the site to row j
    }
  }
}
SAC <- SAC %>% filter(!is.na(spN))
SAC$method <- 'MaxN_cameras'

# do again, filtering out the species that REST didn't have
SAC2 <- dt_ruv %>% mutate(site_cam = paste(site_ID, Camera, sep = "_"), site_sp = paste(site_ID, Taxon, sep = "_")) %>% 
  filter(site_sp %in% with(dt_rest %>% filter(!is.na(D)), paste(site_ID, Taxon, sep="_"))) %>% 
  select(site_ID, site_cam, site_sp, SACentry)  # matches the row index of dt_ruv

# make the SAC times relative to the time of first observation
for (i in 1:11) {
  time1 <- SAC2$SACentry[which(SAC2$site_cam == unique(SAC2$site_cam)[i])[1]]
  for (j in which(SAC2$site_cam == unique(SAC2$site_cam)[i])) {
    SAC2$SACentry[j] <- SAC2$SACentry[j] - time1 # take every time stamp and subtract time1 from it
  }
}
rm(time1)

SAC2$spN <- NA
# for every row at every site-camera, calculate the species number accumulation through time
for (i in 1:11) {
  # don't calculate for every row... every 3rd row
  rows = seq(which(SAC2$site_cam == unique(SAC2$site_cam)[i])[1], max(which(SAC2$site_cam == unique(SAC2$site_cam)[i])), by = 3)
  for (j in rows) {
    if (j == rows[1]) { # the first row from each site represents the first species record, so these will always start at 1
      SAC2$spN[j] <- 1
    }
    else {
      SAC2$spN[j] <- SAC2$site_sp[rows[1]:j] %>% n_distinct # number of unique species from the first row of the site to row j
    }
  }
}

# merge + cleanup
SAC2$method <- 'REST_cameras'
SAC <- bind_rows(SAC, SAC2)
rm(SAC2, rows)
SAC <- SAC %>% filter(!is.na(spN))
SAC$site_sp <- NULL
SAC$time <- SAC$SACentry / 60 # to minutes
SAC$SACentry <- NULL

# generate accumulation for all site cameras
SAC2 <- dt_ruv %>% mutate(site_cam = paste(site_ID, Camera, sep = "_"), site_sp = paste(site_ID, Taxon, sep = "_")) %>% 
  filter(site_sp %in% with(dt_rest %>% filter(!is.na(D)), paste(site_ID, Taxon, sep="_"))) %>% 
  select(site_ID, site_cam, site_sp, SACentry)  # matches the row index of dt_ruv

# make the SAC times relative to the time of first observation
for (i in 1:11) {
  time1 <- SAC2$SACentry[which(SAC2$site_cam == unique(SAC2$site_cam)[i])[1]]
  for (j in which(SAC2$site_cam == unique(SAC2$site_cam)[i])) {
    SAC2$SACentry[j] <- SAC2$SACentry[j] - time1 # take every time stamp and subtract time1 from it
  }
}
rm(time1)

SAC2$spN <- NA
SAC2 <- SAC2 %>% select(!site_cam) %>% mutate(time = round(SACentry/60, digits = 1)) %>% 
  arrange(site_ID, time)
SAC2 <- SAC2 %>% distinct(site_ID, time, site_sp, spN)
# for every row at every site-camera, calculate the species number accumulation through time
for (i in 1:3) {
  rows = which(SAC2$site_ID == unique(SAC2$site_ID)[i])
  for (j in rows) {
    if (j == rows[1]) { # the first row from each site represents the first species record, so these will always start at 1
      SAC2$spN[j] <- 1
    }
    else {
      SAC2$spN[j] <- SAC2$site_sp[rows[1]:j] %>% n_distinct # number of unique species from the first row of the site to row j
    }
  }
}

SAC2 <- SAC2 %>% group_by(site_ID, time) %>% summarise(spN = max(spN)) %>% mutate(method = 'REST')
SAC <- bind_rows(SAC, SAC2) # merge in the site-aggregate REST

# now do the same for MaxN
SAC2 <- dt_ruv %>% mutate(site_cam = paste(site_ID, Camera, sep = "_")) %>% # camera level time stamps are still relevant at this point
  select(site_ID, site_cam, SACentry, Taxon)  # matches the row index of dt_ruv

# make the SAC times relative to the time of first observation
for (i in 1:11) {
  time1 <- SAC2$SACentry[which(SAC2$site_cam == unique(SAC2$site_cam)[i])[1]]
  for (j in which(SAC2$site_cam == unique(SAC2$site_cam)[i])) {
    SAC2$SACentry[j] <- SAC2$SACentry[j] - time1 # take every time stamp and subtract time1 from it
  }
}
rm(time1)

SAC2$spN <- NA
SAC2 <- SAC2 %>% select(!site_cam) %>% mutate(time = round(SACentry/60, digits = 1)) %>% 
  arrange(site_ID, time)
SAC2 <- SAC2 %>% distinct(site_ID, time, Taxon, spN)
# for every row at every site-camera, calculate the species number accumulation through time
for (i in 1:3) {
  rows = which(SAC2$site_ID == unique(SAC2$site_ID)[i])
  for (j in rows) {
    if (j == rows[1]) { # the first row from each site represents the first species record, so these will always start at 1
      SAC2$spN[j] <- 1
    }
    else {
      SAC2$spN[j] <- SAC2$Taxon[rows[1]:j] %>% n_distinct # number of unique species from the first row of the site to row j
    }
  }
}

SAC2 <- SAC2 %>% group_by(site_ID, time) %>% summarise(spN = max(spN)) %>% mutate(method = 'MaxN')
SAC <- bind_rows(SAC, SAC2) # merge in

# Add UVC accumulation
temp <- dt_uvc %>% filter(str_detect(site_ID, 'hinalea|kaku|sunset'))  %>% 
  select(site_ID, transect_point, Taxon) %>% arrange(site_ID, transect_point)
temp$spN <- NA # calculate cumulative S
for (i in 1:3) { # for each site
  rows = which(temp$site_ID == unique(temp$site_ID)[i])
  for (j in rows) { # each row for that site
    if (j == rows[1]) { # the first row from each site represents the first species record, so these will always start at 1
      temp$spN[j] <- 1
    }
    else {
      temp$spN[j] <- length(unique(temp$Taxon[rows[1]:j])) # number of unique species from the first row of the site to row j
    }
  }
}
rm(rows)
temp$Taxon <- NULL
SAC <- temp %>% group_by(site_ID, transect_point) %>% summarise(spN = max(spN)) %>% 
  mutate(method = 'UVC') %>% rename(time = transect_point) %>% 
  bind_rows(., SAC)

sac_fit <- as.list(rep(NA, 9))
tempSAC <-  SAC %>% filter(str_detect(method, 'cameras') == F) %>% select(!site_cam)
for (i in 1:3) {
  for (j in 1:3) {
    sac_fit[[i + (j-1)*3]] <- nls(data = tempSAC %>% 
                                     filter(site_ID == unique(tempSAC$site_ID)[i], method == unique(tempSAC$method)[j]),
                         formula = spN ~ a - ((a-b)*exp(-c*time)),
                         start = list(a = 11, b = 1, c = 0.1))
  }
}

sac_pred <- as.list(rep(NA, 9))
newx <- as_tibble(seq(1, 50, length.out = 200))
names(newx) <- 'time'
for (i in 1:3) {
  for (j in 1:3) {
    sac_pred[[i + (j-1)*3]] <- nlsBootPredict(
      nlsBoot(sac_fit[[i + (j-1)*3]], niter=500), 
      newdata = newx, interval = "confidence") %>% as.data.frame
    colnames(sac_pred[[i + (j-1)*3]]) <- c('fit', 'lwr', 'upr')
    sac_pred[[i + (j-1)*3]] <- data.frame(site_ID = unique(tempSAC$site_ID)[i], method = unique(tempSAC$method)[j],
                                 time = newx) %>% bind_cols(., sac_pred[[i + (j-1)*3]])
  }
}



f_SAC <- ggplot() +
  geom_point(data = tempSAC %>% filter(method != 'UVC'), aes(x = time, y = spN, color = method), shape = 21, alpha = 0.1) +
  geom_point(data = tempSAC %>% filter(method == 'UVC'), aes(x = time, y = spN, color = method), shape = 21) +
  geom_ribbon(data = bind_rows(sac_pred), aes(x = time, ymin = lwr, ymax = upr, group = interaction(site_ID, method),
                                              fill = method), alpha = 0.4) +
  geom_line(data = bind_rows(sac_pred), aes(x = time, y = fit, group = interaction(site_ID, method),
                                            color = method)) +
  scale_color_cherulean(discrete = T, palette = 'gomphosus', name = 'Method') +
  scale_fill_cherulean(discrete = T, palette = 'gomphosus', name = 'Method') + looks +
  facet_wrap(vars(site_ID)) +
  labs(x = "Time elapsed (s)", y = "Number of species") + guides(color = 'none', fill = 'none')


rm(sac_fit, rar_fit, newx, temp, dis_maxn, i, j, pretty)


# Print figures -----------------------------------------------------------

(f_relpcoa1 / f_pcoa1) + plot_layout(guides = "collect")
(f_pcoacam[[1]]|f_pcoacam[[2]]) + plot_layout(guides = "collect")
f_SAC * theme(legend.position = 'bottom')
f_rarefaction

((f_SAD + theme(legend.position = 'none')) / f_SAC / f_rarefaction)
ggsave('../outputs/figures/fig2a-c.pdf', width = 140, height = 180, units = 'mm')
