

## FishN: comparing fish assemblage abundance surveying methods
# REST: random encounter staying time modelling (no individual recognition)

# set up working environment
require(dplyr)
require(lubridate)
require(beepr)
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
sites = 2
cams = c(4,3) # in the order of sites shown in unique(data_ruv$site_ID)
# because the model loop runs on data_ruv

# make a reference vector for species that REST can run on
# criteria: at least 10 detections total at a site
output_poisson <- data_ruv %>% group_by(site_ID, Taxon, Size_class) %>% summarise(occurrences = sum(Count)) %>% ungroup %>% 
  filter(occurrences >= 5) %>% select(!occurrences)

data_ruv %>% group_by(site_ID, Taxon, Size_class) %>% summarise(occurrences = sum(Count)) %>% View
# how many species-length classes can't fit
less5 <- data_ruv %>% group_by(site_ID, Taxon, Size_class) %>% summarise(occurrences = sum(Count)) %>% ungroup %>% 
  filter(occurrences < 5) %>% nrow
# total species size class
spsize <- data_ruv %>% group_by(site_ID, Taxon, Size_class) %>% summarise(occurrences = sum(Count)) %>% ungroup %>% nrow
less5/spsize

rm(less5, spsize)

# add a col on the number of cameras for the sites
output_poisson <- distinct(data_ruv, site_ID, Camera) %>% group_by(site_ID) %>% 
  summarise(cam = n_distinct(Camera)) %>% ungroup() %>% 
  full_join(., output_poisson, by="site_ID")

# MLE REST Poisson fitting-------------------------------------------------------------------------

## Fit 1: 3 parameters -----------------------------------------------------

# establish the likelihood function
# Poisson for estimating expected detections
# Gamma for estimating expected staying time
# requires inputs
# y = vector of detections per camera i
# t = dataframe or matrix of staying times per detection j for camera i, with staying times as rows, cameras in columns
# n = number of cameras

# Poisson needs lambda = mean = variance
# Gamma needs shape = mean, scale = ratio between variance and mean. Set scale = 1 by default. 

output_poisson <- mutate(output_poisson, shapeT = NA, scaleT = NA, lambdaY = NA, Lvalue = NA, estPar = NA) %>% as.data.frame() # add blank cols to fill with par estimates
output_poisson_hess <- as.list(rep(0, nrow(output_poisson)))

# NEST LEVELS: species-site-camera-detection
for (i in 1:nrow(output_poisson)) {
  
    set.seed(240) # replicability
    # vector of detections per camera i at site
    detects <- data_ruv %>% 
      filter(site_ID == output_poisson$site_ID[i], Taxon == output_poisson$Taxon[i], Size_class == output_poisson$Size_class[i]) %>% 
      group_by(Camera) %>% summarise(detects = sum(Count)) %>% # Counts because some detections are grouped
      pull(detects)
    
    # Establish MLE likelihood function
    REST_poisson <- function(par) { # where par is a vector c(gamma shape, gamma scale, poisson lambda)
      Tij <- rep(0,output_poisson$cam[i]) # empty vector to fill in with Tij estimates per camera
      for (n in 1:output_poisson$cam[i]) {
        # make a data object for staying time per site-camera
        stay <- data_ruv %>% filter(site_ID == output_poisson$site_ID[i], Taxon == output_poisson$Taxon[i], Size_class == output_poisson$Size_class[i],
                                    Camera == sitecam[sitecam$site_ID == output_poisson$site_ID[i], 2][n]) %>% pull(staytime)
        # do the most nested level first, Tij staying time
        for (x in 1:n) {
          Tij[x] <- dgamma(x = stay, shape = par[1], scale = par[2], log = T) %>% sum # params[1] = shape
        }
      }
      return(
        -1 * ( sum(pois(x = detects, lambda = par[3], log = T)) + sum(Tij) ) # params[2] = lambda
      )
    }
    # put that through the optimise function
    # with error handling
    optout <- tryCatch(optim(par = c(0.5,0.5,9), REST_poisson, method = "L-BFGS-B", lower = c(0,0,0), upper = c(2000,9,450), hessian = T),
                       error = function(e) e)
    if(inherits(optout, "error")) next # if an error message gets generated this run, move to the next iteration of the loop
    
    output_poisson[i,5:7] <- optout$par
    output_poisson$Lvalue[i] <- optout$value
    output_poisson_hess[[i]] <- optout$hessian
    output_poisson$estPar[i] <- 3
}
beep()

## Fit 2: 2 parameter (scale = 1) ------------------------------------------

# for models that couldn't converge/optimise, set scale = 1

for (i in 1:nrow(output_poisson)) {
  
  if(is.na(output_poisson$shapeT[i]) == F) next # skip if previously estimated
  
  set.seed(240) # replicability
  # vector of detections per camera i at site
  detects <- data_ruv %>% 
    filter(site_ID == output_poisson$site_ID[i], Taxon == output_poisson$Taxon[i], Size_class == output_poisson$Size_class[i]) %>% 
    group_by(Camera) %>% summarise(detects = sum(Count)) %>% # Counts because some detections are grouped
    pull(detects)
  
  # Establish MLE likelihood function
  REST_poisson <- function(par) { # where par is a vector c(gamma shape, gamma scale, poisson lambda)
    Tij <- rep(0,output_poisson$cam[i]) # empty vector to fill in with Tij estimates per camera
    for (n in 1:output_poisson$cam[i]) {
      # make a data object for staying time per site-camera
      stay <- data_ruv %>% filter(site_ID == output_poisson$site_ID[i], Taxon == output_poisson$Taxon[i], Size_class == output_poisson$Size_class[i],
                                  Camera == sitecam[sitecam$site_ID == output_poisson$site_ID[i], 2][n]) %>% pull(staytime)
      # do the most nested level first, Tij staying time
      for (x in 1:n) {
        Tij[x] <- dgamma(x = stay, shape = par[1], scale = 1, log = T) %>% sum # params[1] = shape
      }
    }
    return(
      -1 * (sum(dpois(x = detects, lambda = par[2], log = T)) + sum(Tij)) # params[2] = lambda
    )
  }
  # put that through the optimise function
  # with error handling
  optout <- tryCatch(optim(par = c(0.5,9), REST_poisson, method = "L-BFGS-B", lower = c(0,0), upper = c(2000,450), hessian = T),
                     error = function(e) e)
  if(inherits(optout, "error")) next # if an error message gets generated this run, move to the next iteration of the loop
  
  output_poisson[i,c(5,7)] <- optout$par
  output_poisson$Lvalue[i] <- optout$value
  output_poisson_hess[[i]] <- optout$hessian
  output_poisson$estPar[i] <- 2
  output_poisson$scaleT[i] <- 1 # correctly fill in the parameter
}

beep()

##  Calculate variance from Hessian -----------------------------------------------------

# going to need to calculate SEs differently depending on number of parameters fitted

output_poisson_hess[which(sapply(output_poisson_hess, function(e) sum(e) == 0))] # check any that have "improper" fits giving 0 in hessian
output_poisson_hess[which(sapply(output_poisson_hess, function(e) 0 %in% e))] # check any that have "improper" fits giving 0 in hessian

notNA2 <- which(output_poisson$estPar == 2) # 2 par models
# notNA3 <- which(output_poisson$pars == 3) # most were 2 par fits
# reference vector for species with valid Hessian matrices
output_poisson_hess[notNA2]
# output_poisson_hess[notNA3]

output_poisson$shape_se <- NA
output_poisson$scale_se <- NA
output_poisson$lambda_se <- NA
for (i in 1:length(notNA2)) {
  se <- tryCatch(diag(solve(output_poisson_hess[[notNA2[i]]])), error = function(e) e) # variances!
  if(inherits(se, "error")) next
  output_poisson[notNA2[i],c(10,12)] <- se
}
for (i in 1:length(notNA3)) {
  se <- tryCatch(diag(solve(output_poisson_hess[[notNA3[i]]])), error = function(e) e) # variances!
  if(inherits(se, "error")) next
  output_poisson[notNA3[i],10:12] <- se
}

# check how many species have singular standard errors (i.e. unreliable models)
output_poisson %>% 
  filter(is.na(shape_se)| is.na(lambda_se), is.na(D) == F) %>% nrow # 4 species
output_poisson %>% filter(is.na(shape_se) | is.na(lambda_se)) %>% View

output_poisson_hess[which(is.na(output_poisson$shapeT) == F & is.na(output_poisson$shape_se) == T)] # look at the invalid hessians


## Fit 3: 2 parameter fit for species with shit Hessians -------------------

# consider running a 2 par model for them
# otherwise, if they're already 2 par and still shit, keep blank.

# remove those and rerun
output_poisson[is.na(output_poisson$shape_se) & is.na(output_poisson$scale_se) & is.na(output_poisson$lambda_se), 3:8] <- NA

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
    output_poisson$Lvalue[2*i-s+1] <- optout$value
    output_poisson_hess[[2*i-s+1]] <- optout$hessian
    
  }
}

beep()

output_poisson$pars[is.na(output_poisson$shapeT) == F & is.na(output_poisson$scaleT) == T] <- 2
output_poisson$scaleT[is.na(output_poisson$shapeT) == F & is.na(output_poisson$scaleT) == T] <- 1

# recalculate variances for those shit outliers
notNA2 <- which(is.na(output_poisson$D) == F & output_poisson$scaleT == 1) # 2 par models
for (i in 1:length(notNA2)) {
  se <- tryCatch(diag(solve(output_poisson_hess[[notNA2[i]]])), error = function(e) e) # variances!
  if(inherits(se, "error")) next
  output_poisson[notNA2[i],c(9,11)] <- se
}

# ok, if it has not converged in three iterations, remove
output_poisson <- output_poisson %>% filter(!is.na(lambda_se))

# interpret outputs
# gamma shape = mean = E(T) expected staying time
# poisson lambda = mean = variance = E(Y) expected number of detections

# Calculate densities -----------------------------------------------------

# calculate D from E(T) and E(Y), D = E(Y)E(T) / (sH)
# s = detection zone, 4 sq m
# H = observation period, 45 minutes = 2700 s, seconds to line up with

output_poisson$D <- with(output_poisson, (shapeT * scaleT * lambdaY) / (2700 * 4 * cam)) * 250 # scale by 250 to put density values as per transect
# clear up temporary objects from model fitting
rm(s, se, restsp2, notNA2, notNA, notNA3, i, n, sitecam)

# Bootstrap confidence intervals ------------------------------------------

# make empty lists to store bootstrapping per species
bootTY <- as.list(rep(0,nrow(output_poisson)))
# empty dataframe to store the estimated confidence interval bounds
Dpredict <- matrix(nrow = nrow(output_poisson), ncol = 2)
colnames(Dpredict) <- c('pred', 'se')
set.seed(240) # reproducibility :)
boots = 5000 # how many iterations

for (i in 1:nrow(output_poisson)) {

  bootTY[[i]] <- matrix(nrow = boots, ncol = 3) # empty matrix for each bootstrap iteration
  colnames(bootTY[[i]]) <- c('predT', 'predY', 'predD')
  
  # loop for bootstrap iterations per REST model
  # pick out random T and Y from model estimated distributions
  for (j in 1:boots) {
    bootTY[[i]][j,1] <-  rgamma(1, shape = output_poisson$shapeT[i], scale = output_poisson$scaleT[i])
    bootTY[[i]][j,2] <-  rpois(1, lambda = output_poisson$lambdaY[i])
  }
  
  # calculate D from the sampled T and Y parameters
  bootTY[[i]][,3] = (bootTY[[i]][,1] * bootTY[[i]][,2] / (2700 * 4 * output_poisson$cam[i])) * 250
}
beep()

hist(bootTY[[runif(1, 1, nrow(output_poisson))]][,3]) # distribution of a random species' bootstrapped density estimate
# see if it fits log-normal distribution

# will a random species' density distribution look normal if we log it
hist(log(bootTY[[runif(1, 1, nrow(output_poisson))]][,3]))

for (i in 1:nrow(output_poisson)) {
  bootTY[[i]][which(is.infinite(log(bootTY[[i]][,3]))),3] <- NA # replace all infinite values with NA
  Dpredict[i,1] <- mean(log(bootTY[[i]][,3]), na.rm = T)
  Dpredict[i,2] <- sd(log(bootTY[[i]][,3]), na.rm = T)
}

Dpredict <- as_tibble(Dpredict)
# confidence intervals from the SE
Dpredict$lwr <- with(Dpredict, pred - 1.96*se)
Dpredict$upr <- with(Dpredict, pred + 1.96*se)
Dpredict <- exp(Dpredict) # back transform

# merge confidence intervals 
output_poisson <- bind_cols(output_poisson, Dpredict)

# calculate AICc so we can compare with negative binomial models


# export outputs
# write.csv(final_poisson, "../outputs/REST_poisson.csv", row.names = F)
