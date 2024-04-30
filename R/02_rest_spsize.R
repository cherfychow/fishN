

## FishN: comparing fish assemblage abundance surveying methods
# REST: random encounter staying time modelling (no individual recognition)
# Poisson fitting

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
summary(data_ruv)
data_ruv <- data_ruv %>% filter(staytime > 0) # remove any zero second occurrences

## define some global parameters
# site-camera key
sitecam <- distinct(data_ruv, site_ID, Camera)
sites = 3
cams = c(4,3,4) # in the order of sites shown in unique(data_ruv$site_ID)
# because the model loop runs on data_ruv

# make a reference vector for species that REST can run on
# criteria: at least 10 detections total at a site
rest_output <- data_ruv %>% group_by(site_ID, Taxon, Size_class) %>% summarise(occurrences = sum(Count)) %>% ungroup %>% 
  filter(occurrences >= 5) %>% select(!occurrences)

data_ruv %>% group_by(site_ID, Taxon, Size_class) %>% summarise(occurrences = sum(Count)) %>% View
# how many species-length classes can't fit
less5 <- data_ruv %>% group_by(site_ID, Taxon, Size_class) %>% summarise(occurrences = sum(Count)) %>% ungroup %>% 
  filter(occurrences < 5) %>% nrow
# total species size class
spsize <- data_ruv %>% group_by(site_ID, Taxon, Size_class) %>% summarise(occurrences = sum(Count)) %>% ungroup %>% nrow
less5/spsize # 0.4545

rm(less5, spsize)

# add a col on the number of cameras for the sites
rest_output <- distinct(data_ruv, site_ID, Camera) %>% group_by(site_ID) %>% 
  summarise(cam = n_distinct(Camera)) %>% ungroup() %>% 
  full_join(., rest_output, by="site_ID")

# MLE REST Poisson fitting-------------------------------------------------------------------------

# constants

s = 6.75 # estimated GoPro FOV underwater
H = 2700 # 45 minutes in seconds

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

rest_output <- mutate(rest_output, shapeT = NA, scaleT = NA, lambdaY = NA, Lvalue = NA) %>% as.data.frame() # add blank cols to fill with par estimates
rest_output_hess <- as.list(rep(0, nrow(rest_output)))

# NEST LEVELS: species-site-camera-detection
for (i in 1:nrow(rest_output)) {
  
  set.seed(240) # replicability
  # vector of detections per camera i at site
  detects <- data_ruv %>% 
    filter(site_ID == rest_output$site_ID[i], Taxon == rest_output$Taxon[i], Size_class == rest_output$Size_class[i]) %>% 
    group_by(Camera) %>% summarise(detects = sum(Count)) %>% # Counts because some detections are grouped
    pull(detects)
  
  # Establish MLE likelihood function
  REST_poisson <- function(par) { # where par is a vector c(gamma shape, gamma scale, poisson lambda)
    Tij <- rep(0,rest_output$cam[i]) # empty vector to fill in with Tij estimates per camera
    for (n in 1:rest_output$cam[i]) {
      # make a data object for staying time per site-camera
      stay <- data_ruv %>% filter(site_ID == rest_output$site_ID[i], Taxon == rest_output$Taxon[i], Size_class == rest_output$Size_class[i],
                                  Camera == sitecam[sitecam$site_ID == rest_output$site_ID[i], 2][n]) %>% pull(staytime)
      # do the most nested level first, Tij staying time
      for (x in 1:n) {
        Tij[x] <- dgamma(x = stay, shape = par[1], scale = par[2], log = T) %>% sum # params[1] = shape
      }
    }
    return(
      -1 * ( sum(dpois(x = detects, lambda = par[3], log = T)) + sum(Tij) ) # params[2] = lambda
    )
  }
  # put that through the optimise function
  # with error handling
  optout <- tryCatch(optim(par = c(0.5,0.5,9), REST_poisson, method = "L-BFGS-B", lower = c(0,0,0), upper = c(2000,9,450), hessian = T),
                     error = function(e) e)
  if(inherits(optout, "error")) next # if an error message gets generated this run, move to the next iteration of the loop
  
  rest_output[i,5:7] <- optout$par
  rest_output$Lvalue[i] <- optout$value
  rest_output_hess[[i]] <- optout$hessian
}
beep()

rest_output$fit <- NA
rest_output$fit[which(is.na(rest_output$shapeT) == F)] <- 1

rest_output %>% filter(is.na(shapeT)) %>% nrow # how many unconverged
rest_output %>% filter(is.na(shapeT)== F) %>% nrow # 54 models converged


## Fit 2: use the density equation to reduce estimated parameters ------------------------------------------

# for models that couldn't converge/optimise, make lambda = DsH / ktheta

for (i in 1:nrow(rest_output)) {
  
  if(is.na(rest_output$shapeT[i]) == F) next # skip if previously estimated
  
  set.seed(240) # replicability
  # vector of detections per camera i at site
  detects <- data_ruv %>% 
    filter(site_ID == rest_output$site_ID[i], Taxon == rest_output$Taxon[i], Size_class == rest_output$Size_class[i]) %>% 
    group_by(Camera) %>% summarise(detects = sum(Count)) %>% # Counts because some detections are grouped
    pull(detects)
  
  # calculate the cumulative staying time parameter DsH
  DsH <- data_ruv %>% filter(site_ID == rest_output$site_ID[i], 
                             Taxon == rest_output$Taxon[i], 
                             Size_class == rest_output$Size_class[i]) %>% pull(staytime) %>% sum
  
  # Establish MLE likelihood function
  REST_poisson2 <- function(par) { # where par is a vector c(gamma shape, gamma scale, poisson lambda)
    Tij <- rep(0,rest_output$cam[i]) # empty vector to fill in with Tij estimates per camera
    for (n in 1:rest_output$cam[i]) {
      # make a data object for staying time per site-camera
      stay <- data_ruv %>% filter(site_ID == rest_output$site_ID[i], Taxon == rest_output$Taxon[i], Size_class == rest_output$Size_class[i],
                                  Camera == sitecam[sitecam$site_ID == rest_output$site_ID[i], 2][n]) %>% pull(staytime)
      # do the most nested level first, Tij staying time
      for (x in 1:n) {
        Tij[x] <- dgamma(x = stay, shape = par[1], scale = par[2], log = T) %>% sum # params[1] = shape
      }
    }
    return(
      -1 * (sum(dpois(x = detects, lambda = DsH/(par[1] * par[2]), log = T)) + sum(Tij)) # params[2] = lambda
    )
  }
  # put that through the optimise function
  # with error handling
  optout <- tryCatch(optim(par = c(0.5,9), REST_poisson2, method = "L-BFGS-B", lower = c(0,0), upper = c(2000,450), hessian = T),
                     error = function(e) e)
  if(inherits(optout, "error")) next # if an error message gets generated this run,move to the next iteration of the loop
  
  rest_output[i,c(5:6)] <- optout$par # only two parameters this fit
  rest_output[i,7] <- DsH/prod(optout$par) # calculate lambda from shape and scale
  rest_output$Lvalue[i] <- optout$value
  rest_output_hess[[i]] <- optout$hessian
  rest_output$scaleT[i] <- 1 # correctly fill in the parameter
}

beep()

# label these with the type of fitting
rest_output$fit[which(is.na(rest_output$shapeT) == F & is.na(rest_output$fit))] <- 2

##  Calculate variance from Hessian -----------------------------------------------------

# check and calculate variances

rest_output_hess[which(sapply(rest_output_hess, function(e) sum(e) == 0))] # check any that have "improper" fits giving 0 in hessian
rest_output_hess[which(sapply(rest_output_hess, function(e) 0 %in% e))] # check any that have "improper" fits giving 0 in hessian

rest_output$shape_se <- NA
rest_output$scale_se <- NA
rest_output$lambda_se <- NA

# variances for fit 1 (3 parameters)

for (i in 1:length(which(rest_output$fit == 1))) {
  se <- tryCatch(diag(solve(rest_output_hess[[which(rest_output$fit == 1)[i]]])), error = function(e) e) # variances!
  if(inherits(se, "error")) next
  rest_output[which(rest_output$fit == 1)[i],10:12] <- se
}

# for fit 2, 2 parameters

for (i in 1:length(which(rest_output$fit == 2))) {
  se <- tryCatch(diag(solve(rest_output_hess[[which(rest_output$fit == 2)[i]]])), error = function(e) e) # variances!
  if(inherits(se, "error")) next
  rest_output[which(rest_output$fit == 2)[i],c(10,11)] <- se
}


# check how many species have singular standard errors (i.e. unreliable models)
rest_output %>% 
  filter(is.na(shape_se) & is.na(shapeT) == F) %>% nrow # 14 with invalid hessians
rest_output %>% filter(is.na(shape_se) & is.na(shapeT) == F) %>% View # all the same with low likelihoods

rest_output_hess[which(is.na(rest_output$shapeT) == F & is.na(rest_output$shape_se) == T)] # look at the invalid hessians


## Fit 3: 2 parameter fit for species with shit Hessians -------------------

# consider running a 2 par model for them

# remove those and rerun
reruns <- which(is.na(rest_output$shapeT) == F & is.na(rest_output$shape_se) == T)
rest_output[reruns, c(5:7, 10:12)] <- NA

for (i in reruns) {
  
  set.seed(240) # replicability
  # vector of detections per camera i at site
  detects <- data_ruv %>% 
    filter(site_ID == rest_output$site_ID[i], Taxon == rest_output$Taxon[i], Size_class == rest_output$Size_class[i]) %>% 
    group_by(Camera) %>% summarise(detects = sum(Count)) %>% # Counts because some detections are grouped
    pull(detects)
  
  # calculate the cumulative staying time parameter DsH
  DsH <- data_ruv %>% filter(site_ID == rest_output$site_ID[i], 
                             Taxon == rest_output$Taxon[i], 
                             Size_class == rest_output$Size_class[i]) %>% pull(staytime) %>% sum
  
  # Establish MLE likelihood function
  REST_poisson2 <- function(par) { # where par is a vector c(gamma shape, gamma scale, poisson lambda)
    Tij <- rep(0,rest_output$cam[i]) # empty vector to fill in with Tij estimates per camera
    for (n in 1:rest_output$cam[i]) {
      # make a data object for staying time per site-camera
      stay <- data_ruv %>% filter(site_ID == rest_output$site_ID[i], Taxon == rest_output$Taxon[i], Size_class == rest_output$Size_class[i],
                                  Camera == sitecam[sitecam$site_ID == rest_output$site_ID[i], 2][n]) %>% pull(staytime)
      # do the most nested level first, Tij staying time
      for (x in 1:n) {
        Tij[x] <- dgamma(x = stay, shape = par[1], scale = par[2], log = T) %>% sum # params[1] = shape
      }
    }
    return(
      -1 * (sum(dpois(x = detects, lambda = DsH/(par[1] * par[2]), log = T)) + sum(Tij)) # params[2] = lambda
    )
  }
  # put that through the optimise function
  # with error handling
  optout <- tryCatch(optim(par = c(0.5,9), REST_poisson2, method = "L-BFGS-B", lower = c(0,0), upper = c(2000,450), hessian = T),
                     error = function(e) e)
  if(inherits(optout, "error")) next # if an error message gets generated this run,move to the next iteration of the loop
  
  rest_output[i,c(5:6)] <- optout$par # only two parameters this fit
  rest_output[i,7] <- DsH/prod(optout$par) # calculate lambda from shape and scale
  rest_output$Lvalue[i] <- optout$value
  rest_output_hess[[i]] <- optout$hessian
  rest_output$scaleT[i] <- 1 # correctly fill in the parameter
}

beep()

rest_output$fit[which(is.na(rest_output$shapeT) == F & is.na(rest_output$shape_se))] <- 2

# recalculate variances for those shit outliers
var2 <- which(is.na(rest_output$shapeT) == F & is.na(rest_output$shape_se))
for (i in var2) {
  se <- tryCatch(diag(solve(rest_output_hess[[i]])), error = function(e) e) # variances!
  if(inherits(se, "error")) next
  rest_output[i,c(10,11)] <- se
}

# interpret outputs
# gamma shape = mean = E(T) expected staying time
# poisson lambda = mean = variance = E(Y) expected number of detections

# Calculate densities -----------------------------------------------------

# calculate D from E(T) and E(Y), D = E(Y)E(T) / (sH)
# s = detection zone, 4 sq m
# H = observation period, 45 minutes = 2700 s, seconds to line up with

rest_output$D <- with(rest_output, (shapeT * scaleT * lambdaY) / (H * s)) * 250 # scale by 250 to put density values as per transect
# clear up temporary objects from model fitting
rm(se, var2, i, n, sitecam, detects, j, rest_output_hess, optout, cams, reruns, sites, DsH)

# Parametric bootstrapping ------------------------------------------

# make empty lists to store bootstrapping per species
D_boot <- as.list(rep(0,nrow(rest_output)))
# empty dataframe to store the estimated confidence interval bounds
D_pred <- matrix(nrow = nrow(rest_output), ncol = 3)
colnames(D_pred) <- c('pred', 'lwr', 'upr')
set.seed(240) # reproducibility :)
boots = 5000 # how many iterations

for (i in 1:nrow(rest_output)) {
  
  D_boot[[i]] <- matrix(nrow = boots, ncol = 3) # empty matrix for each bootstrap iteration
  colnames(D_boot[[i]]) <- c('predT', 'predY', 'predD')
  
  # loop for bootstrap iterations per REST model
  # pick out random T and Y from model estimated distributions
  for (j in 1:boots) {
    D_boot[[i]][j,1] <-  rgamma(1, shape = rest_output$shapeT[i], scale = rest_output$scaleT[i])
    D_boot[[i]][j,2] <-  rpois(1, lambda = rest_output$lambdaY[i])
  }
  
  # calculate D from the sampled T and Y parameters
  D_boot[[i]][,3] = (D_boot[[i]][,1] * D_boot[[i]][,2] / (H * s)) * 250
}
beep()

hist(D_boot[[runif(1, 1, nrow(rest_output))]][,3]) # distribution of a random species' bootstrapped density estimate
# see if it fits log-normal distribution

# calculate confidence interval using percentile method
for (i in 1:nrow(rest_output)) {
  D_boot[[i]] <- D_boot[[i]][order(D_boot[[i]][,3]),] # order ascending values of D
  D_pred[i,1] <- D_boot[[i]][boots/2,3] # 50th percentile for the predicted D
  D_pred[i,2] <- D_boot[[i]][ round(boots*0.025, 1), 3] # 2.5 percentile for lower
  D_pred[i,3] <- D_boot[[i]][ round(boots*0.975, 1), 3] # 97.5 percentile for lower
}

# merge confidence intervals 
rest_output <- bind_cols(rest_output, D_pred)
rest_output <- rest_output %>% filter(!is.na(shapeT))

# calculate AICc so we can compare with negative binomial models
# 2k - 2log(L) + (2k^2 + 2k) / (n - k - 1)

# first need the number of detections
rest_output <- data_ruv %>% 
  group_by(site_ID, Taxon, Size_class) %>% summarise(detects = sum(Count)) %>% # Counts because some detections are grouped
  left_join(rest_output, ., by=c("site_ID", "Taxon", "Size_class")) %>% ungroup  
rest_output$AICc <- with(rest_output, 2*ifelse(fit == 1, 3, 2) - 2*log(Lvalue) + ((2*ifelse(fit == 1, 3, 2)^2 + 2*ifelse(fit == 1, 3, 2)) / (detects - ifelse(fit == 1, 3, 2) - 1)))

# export outputs
readr::write_csv(rest_output, "../outputs/REST_output.csv")
rm(D_boot, D_pred, boots, H, i, j, s)
