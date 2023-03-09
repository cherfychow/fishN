

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

## define some global parameters
# site-camera key
sitecam <- distinct(data_ruv, site_ID, Camera)
sites = 2 # establish number of sites
cams <- c(4,3) # vector of camera number in the order of unique(data_ruv$site_ID)

# make a reference vector for species that REST can run on
# criteria: at least 10 detections total at a site
restsp <- data_ruv %>% group_by(site_ID, Taxon) %>% summarise(occurrences = sum(Count)) %>% ungroup %>% 
  filter(occurrences >= 10) %>% pull(Taxon) %>% sort %>% unique

# establish the likelihood function
# Poisson for estimating expected detections
# Gamma for estimating expected staying time
# requires inputs
# y = vector of detections per camera i
# t = dataframe or matrix of staying times per detection j for camera i, with staying times as rows, cameras in columns
# n = number of cameras

# Poisson needs lambda = mean = variance
# Gamma needs shape = mean, scale = ratio between variance and mean. Set scale = 1 by default. 

output_poisson <- data.frame(site_ID = '', Taxon = rep(restsp, each = 2), shapeT = '', scaleT = '', lambdaY = '')
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
          Tij[x] <- dgamma(x = stay, shape = par[1], scale = par[2], log = T) %>% sum # params[1] = shape
        }
      }
      return(
        -1 * sum( dpois(x = detects, lambda = par[3], log = T) + sum(Tij) ) # params[2] = lambda
      )
    }
    # put that through the optimise function
    # with error handling
    optout <- tryCatch(optim(par = c(0.5,0.5,9), REST_poisson, method = "L-BFGS-B", lower = c(0,0,0), upper = c(2000,9,450), hessian = T),
                       error = function(e) e)
    if(inherits(optout, "error")) next # if an error message gets generated this run, move to the next iteration of the loop
    
    output_poisson[2*i-s+1,3:5] <- optout$par
    output_poisson_hess[[2*i-s+1]] <- optout$hessian
  }
}

# for models that couldn't converge/optimise, set scale = 1

for (i in 1:length(restsp)) {

  for (s in 1:sites) {
    if(is.na(output_poisson$shapeT[2*i-s+1]) == F) next 
    # move on if there's already a model output there
    
    set.seed(240) # replicability
    # vector of detections per camera i at site
    detects <- data_ruv %>% filter(site_ID == unique(data_ruv$site_ID)[s], Taxon == restsp[i]) %>% 
      group_by(Camera) %>% summarise(detects = sum(Count)) %>% # Counts because some detections are grouped
      pull(detects)
    
    # Establish MLE likelihood function
    REST_poisson2 <- function(par) { # where x is a vector c(gamma shape, poisson lambda)
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
    optout <- tryCatch(optim(par = c(0.5,9), REST_poisson2, method = "L-BFGS-B", lower = c(0,0), upper = c(2000,450), hessian = T),
                       error = function(e) e)
    if(inherits(optout, "error")) next # if an error message gets generated this run, move to the next iteration of the loop
    
    # store fitting outputs
    output_poisson[2*i-s+1,c(3,5)] <- optout$par
    output_poisson_hess[[2*i-s+1]] <- optout$hessian
    
  }
}

output_poisson$scaleT[is.na(output_poisson$shapeT) == F & is.na(output_poisson$scaleT) == T] <- 1

# remove temp objects
rm(detects, i, n, s, sp, stay, Tij, x)
output_poisson[[3]] <- as.numeric(output_poisson[[3]])
output_poisson[[4]]  <- as.numeric(output_poisson[[4]])
output_poisson[[5]]  <- as.numeric(output_poisson[[5]])

# calculate D from E(T) and E(Y), D = E(Y)E(T) / (sH)
# s = detection zone, 4 sq m
# H = observation period, 45 minutes = 2700 s, seconds to line up with 
output_poisson$D[output_poisson$site_ID == 'hale_kinalea'] <- (with(output_poisson[output_poisson$site_ID == 'hale_kinalea', ], shapeT * scaleT * lambdaY) / (2700 * 4 * cams[1])) * 250
output_poisson$D[output_poisson$site_ID == 'hale_kaku'] <- (with(output_poisson[output_poisson$site_ID == 'hale_kaku', ], shapeT * scaleT * lambdaY) / (2700 * 4 * cams[2]))*250 # scale to be density within belt transect

# going to need to calculate SEs differently depending on number of parameters fitted

notNA2 <- which(is.na(output_poisson$D) == F & output_poisson$scaleT == 1) # 2 par models
notNA3 <- which(is.na(output_poisson$D) == F & output_poisson$scaleT != 1) # 3 par models
# reference vector for species with valid Hessian matrices
output_poisson_hess[notNA2]
output_poisson_hess[notNA3]
notNA3 <- notNA3[-4] # get rid of 4, all zeroes

output_poisson$shape_var <- NA
output_poisson$scale_var <- NA
output_poisson$lambda_var <- NA
for (i in 1:length(notNA2)) {
  se <- tryCatch(diag(solve(output_poisson_hess[[notNA2[i]]])), error = function(e) e) # variances!
  if(inherits(se, "error")) next
  output_poisson[notNA2[i],c(7,9)] <- se
}
for (i in 1:length(notNA3)) {
  se <- tryCatch(diag(solve(output_poisson_hess[[notNA3[i]]])), error = function(e) e) # variances!
  if(inherits(se, "error")) next
  output_poisson[notNA3[i],7:9] <- se
}

# check how many species have singular standard errors (i.e. unreliable models)
output_poisson %>% 
  filter(is.na(shape_var), is.na(scale_var), is.na(lambda_var), is.na(D) == F) %>% nrow # 4 species
output_poisson %>% filter(is.na(shape_var), is.na(scale_var), is.na(lambda_var), is.na(D) == F) %>% View

# consider running a 2 par model for them

# remove those and rerun
output_poisson[is.na(output_poisson$shape_var) & is.na(output_poisson$scale_var) & is.na(output_poisson$lambda_var), 3:6] <- NA
# also remove the hessian matrices
output_poisson_hess[is.na(output_poisson$shape_var) & is.na(output_poisson$scale_var) & is.na(output_poisson$lambda_var)] <- 0

for (i in 1:length(restsp)) {
  
  for (s in 1:sites) {
    if(is.na(output_poisson$shapeT[2*i-s+1]) == F) next 
    # move on if there's already a model output there
    
    set.seed(240) # replicability
    # vector of detections per camera i at site
    detects <- data_ruv %>% filter(site_ID == unique(data_ruv$site_ID)[s], Taxon == restsp[i]) %>% 
      group_by(Camera) %>% summarise(detects = sum(Count)) %>% # Counts because some detections are grouped
      pull(detects)
    
    # Establish MLE likelihood function
    REST_poisson2 <- function(par) { # where x is a vector c(gamma shape, poisson lambda)
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
    optout <- tryCatch(optim(par = c(0.5,9), REST_poisson2, method = "L-BFGS-B", lower = c(0,0), upper = c(2000,450), hessian = T),
                       error = function(e) e)
    if(inherits(optout, "error")) next # if an error message gets generated this run, move to the next iteration of the loop
    
    # store fitting outputs
    output_poisson[2*i-s+1,c(3,5)] <- optout$par
    output_poisson_hess[[2*i-s+1]] <- optout$hessian
    
  }
}


output_poisson$D[output_poisson$site_ID == 'hale_kinalea' & is.na(output_poisson$D)] <- (with(output_poisson[output_poisson$site_ID == 'hale_kinalea' & is.na(output_poisson$D), ], shapeT * lambdaY) / (2700 * 4 * cams[1])) * 250
output_poisson$D[output_poisson$site_ID == 'hale_kaku' & is.na(output_poisson$D)] <- (with(output_poisson[output_poisson$site_ID == 'hale_kaku' & is.na(output_poisson$D), ], shapeT * lambdaY) / (2700 * 4 * cams[2]))*250 # scale to be density within belt transect
output_poisson$scaleT[is.na(output_poisson$shapeT) == F & is.na(output_poisson$scaleT) == T] <- 1

notNA2 <- which(is.na(output_poisson$D) == F & output_poisson$scaleT == 1) # 2 par models
for (i in 1:length(notNA2)) {
  se <- tryCatch(diag(solve(output_poisson_hess[[notNA2[i]]])), error = function(e) e) # variances!
  if(inherits(se, "error")) next
  output_poisson[notNA2[i],c(7,9)] <- se
}

# no change to these species, filter out. everyone else can use MaxN
output_poisson <- output_poisson %>% filter(is.na(shape_var) == F)

# calculate SE of D

# var(TY) = var(T) * var(Y) + var(T) * EY^2 + var(Y) * ET^2
# SE of D = sqrt(var(TY))
# treat sH as constant scaling factors (units of measurement)
# for staying time, because it's a gamma distribution, E(T) = shape * scale
output_poisson$D_SE <- with(output_poisson, ET_SE*EY_SE + ET_SE*lambdaY^2 + EY_SE*(shapeT*scaleT)^2)
# scale to line up the same "units"
output_poisson$D_SE[output_poisson$site_ID == 'hale_kaku'] <- (output_poisson$D_SE[output_poisson$site_ID == 'hale_kaku'] / (2700 * 4 * 3)) * 250
output_poisson$D_SE[output_poisson$site_ID == 'hale_kinalea'] <- (output_poisson$D_SE[output_poisson$site_ID == 'hale_kinalea'] / (2700 * 4 * 4)) * 250

# interpret outputs
# gamma shape = mean = E(T) expected staying time
# poisson lambda = mean = variance = E(Y) expected number of detections


## Try negative binomial instead of poisson




# export outputs
write.csv(output_poisson, "../outputs/REST_poisson.csv", row.names = F)
