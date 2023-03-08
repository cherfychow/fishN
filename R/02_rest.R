

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
output_poisson_hess <- as.list(rep(0, length(restsp)*2))

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
          Tij[x] <- dgamma(x = stay, shape = par[1], scale = 1, log = T) %>% sum # params[1] = shape
        }
      }
      return(
        -1 * sum( dpois(x = detects, lambda = par[2], log = T) + sum(Tij) ) # params[2] = lambda
      )
    }
    # put that through the optimise function
    # with error handling
    optout <- tryCatch(optim(par = c(0.5,9), REST_poisson, method = "L-BFGS-B", lower = c(0,0), upper = c(2000,200), hessian = T),
                       error = function(e) e)
    if(inherits(optout, "error")) next # if an error message gets generated this run, move to the next iteration of the loop
    
    output_poisson[2*i-s+1,3:4] <- optout$par
    output_poisson_hess[[2*i-s+1]] <- optout$hessian
    
    # throw in secondary field values for the data output table
    output_poisson$site_ID[2*i-s+1] <- unique(data_ruv$site_ID)[s]
    output_poisson[[3]] <- as.numeric(output_poisson[[3]])
    output_poisson[[4]]  <- as.numeric(output_poisson[[4]])
    output_poisson$D <- (output_poisson$shapeET*output_poisson$lambdaEY / (2700 * 4 * cams[s]))*250 # scale to be density within belt transect
  }
}
# remove temp objects
rm(detects, i, n, s, sp, stay, Tij, x)

notNA <- which(is.na(output_poisson$D) == F) # reference vector for which species with valid Hessian matrices
output_poisson_hess[notNA]
notNA <- notNA[-which(notNA == 6)] # get rid of 6, all zeroes

output_poisson$ET_SE <- NA
output_poisson$EY_SE <- NA
for (i in 1:length(notNA)) {
  se <- tryCatch(sqrt(diag(solve(output_poisson_hess[[notNA[i]]]))), error = function(e) e)
  if(inherits(se, "error")) next
  output_poisson[notNA[i],6:7] <- se
}
purrr::map(output_poisson_hess[notNA], function(x) {sqrt(diag(solve(x)))})

# interpret outputs
# gamma shape = mean = E(T)
# poisson lambda = mean = variance = E(Y)


