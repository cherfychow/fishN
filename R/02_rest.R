

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

# site-camera key
sitecam <- distinct(data_ruv, site_ID, Camera)
sites = 2 # establish number of sites
cams <- c(4,3) # vector of camera number in the order of unique(data_ruv$site_ID)

# make a reference vector for species that REST can run on
restsp <- data_ruv %>% group_by(site_ID, Taxon) %>% summarise(occurrences = sum(Count)) %>% ungroup %>% 
  filter(occurrences >= 10) %>% pull(Taxon) %>% sort %>% unique

output_poisson <- data.frame(site_ID = '', Taxon = rep(restsp, each = 2), shapeET = '', lambdaEY = '')

# NEST LEVELS: species-site-camera-detection

for (i in 1:length(restsp)) {
  
  for (s in 1:sites) {
    
    set.seed(240) # replicability
    # vector of detections per camera i at site
    detects <- data_ruv %>% filter(site_ID == unique(data_ruv$site_ID)[s], Taxon == restsp[i]) %>% 
      group_by(Camera) %>% summarise(detects = sum(Count)) %>% # Counts because some detections are grouped
      pull(detects)
    
    # Establish MLE likelihood function
    REST_poisson <- function(par) { # where x is a vector c(gamma shape, poisson lambda)
      Tij <- rep(0,cams[s])
      for (n in 1:cams[s]) {
        # make a data object for staying time per site-camera
        stay <- data_ruv %>% filter(site_ID == unique(data_ruv$site_ID)[s], Taxon == restsp[i],
                                    Camera == sitecam[sitecam$site_ID == unique(data_ruv$site_ID)[s],2][n]) %>% pull(staytime)
        # do the most nested level first, Tij staying time
        for (x in 1:n) {
          Tij[x] <- sum(dgamma(x = stay, shape = par[1], scale = 1, log = T)) # params[1] = shape
        }
      }
      return(
        -1 * sum( dpois(x = detects, lambda = par[2], log = T) + sum(Tij) ) # params[2] = lambda
      )
    }
    # put that through the optimise function
    output_poisson[2*i-s+1,3:4] <- optim(par = c(1,4), REST_poisson)$par
    output_poisson$site_ID[2*i-s+1] <- unique(data_ruv$site_ID)[s]
  }
}

# interpret outputs
# gamma shape = mean = E(T)
# poisson lambda = mean = variance = E(Y)

output_poisson[[3]] <- as.numeric(output_poisson[[3]])
output_poisson[[4]]  <- as.numeric(output_poisson[[4]])
output_poisson$D <- (output_poisson$shapeET*output_poisson$lambdaEY / (2700 * 4))*250 # scale to be density within belt transect
