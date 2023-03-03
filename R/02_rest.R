

## FishN: comparing fish assemblage abundance surveying methods
# REST: random encounter staying time modelling (no individual recognition)

# set up working environment
require(dplyr)
require(lubridate)
# custom palettes to be extra
source('https://gist.githubusercontent.com/cherfychow/e9ae890fd16f4c86730748c067feee2b/raw/899dcfc1745421cb4e6ba26826b0bfe55fd8ec14/cherulean.R')

data_ruv <- read.csv('../data/ruv_himb_pilot.csv', header = T)

# convert this to a lubridate duration data type
data_ruv$entrytime_c <- ms(data_ruv$Time_entry) %>% as.duration
data_ruv$exittime_c <- ms(data_ruv$Time_exit) %>% as.duration
# calculate staying time duration
data_ruv$staytime <- with(data_ruv, exittime_c - entrytime_c) %>% as.numeric

# establish the likelihood function
# Poisson for estimating expected detections
# Gamma for estimating expected staying time
# requires inputs
# y = vector of detections per camera i
# t = dataframe or matrix of staying times per detection j for camera i, with staying times as rows, cameras in columns
# n = number of cameras

# Poisson needs lambda = mean = variance
# Gamma needs shape = mean, scale = ratio between variance and mean. Set scale = 1 by default. 

sites = 2 # establish number of sites
cams <- c(4,3) # vector of camera number in the order of unique(data_ruv$site_ID)

for (s in 1:sites) {
  
  # make a data object for staying time per site-camera
  stay <- data_ruv %>% filter(site_ID == unique(data_ruv$site_ID)[s]) %>%
    select(Camera, staytime) %>% tidyr::pivot_wider(names_from = Camera, values_fill = NA)
  
  # vector of detections per camera i at site
  detects <- data_ruv %>% filter(site_ID == unique(data_ruv$site_ID)[s]) %>% 
    group_by(Camera) %>% summarise(detects = sum(Count)) %>% # Counts because some detections are grouped
    pull(detects)
  
  for (n in 1:cams) {
    ### Start MLE ###
    REST_poisson <- function(lambda,shape) {
      
      # do the most nested level first, Tij staying time
      Tij <- rep(0,n)
      for (x in 1:n) {
        Tij[x] <- sum(dgamma(x = stay[n], shape, scale = 1, log = T))
      } 
      return(
        -1 * sum( dpois(x = detects, lambda, log = T) + sum(Tij) )) 
      )
    }
    ### end MLE ###