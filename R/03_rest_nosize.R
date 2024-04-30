

## FishN: comparing fish assemblage abundance surveying methods
# REST: random encounter staying time modelling (no individual recognition)
# Poisson fitting

# set up working environment
require(dplyr)
require(lubridate)
require(beepr)
# custom palettes to be extra
source('https://gist.githubusercontent.com/cherfychow/e9ae890fd16f4c86730748c067feee2b/raw/899dcfc1745421cb4e6ba26826b0bfe55fd8ec14/cherulean.R')

## define some global parameters
# site-camera key
sites = 3
cams = c(4,3,4) # in the order of sites shown in unique(data_ruv$site_ID)
# because the model loop runs on data_ruv

# make a reference vector for species that REST can run on
# criteria: at least 5 detections total at a site
output_sp <- data_ruv %>% group_by(site_ID, Taxon) %>% summarise(occurrences = sum(Count)) %>% ungroup %>% 
  filter(occurrences >= 5) %>% select(!occurrences)

data_ruv %>% group_by(site_ID, Taxon) %>% summarise(occurrences = sum(Count)) %>% View
# 31 species records at sites that can't be fit with REST models

# add a col on the number of cameras for the sites
output_sp <- distinct(data_ruv, site_ID, Camera) %>% group_by(site_ID) %>% 
  summarise(cam = n_distinct(Camera)) %>% ungroup() %>% 
  full_join(., output_sp, by="site_ID")

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

output_sp <- mutate(output_sp, shapeT = NA, scaleT = NA, lambdaY = NA, Lvalue = NA, estPar = NA) %>% as.data.frame() # add blank cols to fill with par estimates
output_sp_hess <- as.list(rep(0, nrow(output_sp)))

# NEST LEVELS: species-site-camera-detection
for (i in 1:nrow(output_sp)) {
  
    set.seed(240) # replicability
    # vector of detections per camera i at site
    detects <- data_ruv %>% 
      filter(site_ID == output_sp$site_ID[i], Taxon == output_sp$Taxon[i]) %>% 
      group_by(Camera) %>% summarise(detects = sum(Count)) %>% # Counts because some detections are grouped
      pull(detects)
    
    # Establish MLE likelihood function
    REST_poisson <- function(par) { # where par is a vector c(gamma shape, gamma scale, poisson lambda)
      Tij <- rep(0,output_sp$cam[i]) # empty vector to fill in with Tij estimates per camera
      for (n in 1:output_sp$cam[i]) {
        # make a data object for staying time per site-camera
        stay <- data_ruv %>% filter(site_ID == output_sp$site_ID[i], Taxon == output_sp$Taxon[i],
                                    Camera == sitecam[sitecam$site_ID == output_sp$site_ID[i], 2][n]) %>% pull(staytime)
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
    
    output_sp[i,4:6] <- optout$par
    output_sp$Lvalue[i] <- optout$value
    output_sp_hess[[i]] <- optout$hessian
    output_sp$estPar[i] <- 3
}
beep()

## Fit 2: 2 parameter (scale = 1) ------------------------------------------

# for models that couldn't converge/optimise, set scale = 1

for (i in 1:nrow(output_sp)) {
  
  if(is.na(output_sp$shapeT[i]) == F) next # skip if previously estimated
  
  set.seed(240) # replicability
  # vector of detections per camera i at site
  detects <- data_ruv %>% 
    filter(site_ID == output_sp$site_ID[i], Taxon == output_sp$Taxon[i]) %>% 
    group_by(Camera) %>% summarise(detects = sum(Count)) %>% # Counts because some detections are grouped
    pull(detects)
  
  # Establish MLE likelihood function
  REST_poisson2 <- function(par) { # where par is a vector c(gamma shape, gamma scale, poisson lambda)
    Tij <- rep(0,output_sp$cam[i]) # empty vector to fill in with Tij estimates per camera
    for (n in 1:output_sp$cam[i]) {
      # make a data object for staying time per site-camera
      stay <- data_ruv %>% filter(site_ID == output_sp$site_ID[i], Taxon == output_sp$Taxon[i],
                                  Camera == sitecam[sitecam$site_ID == output_sp$site_ID[i], 2][n]) %>% pull(staytime)
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
  optout <- tryCatch(optim(par = c(0.5,9), REST_poisson2, method = "L-BFGS-B", lower = c(0,0), upper = c(2000,450), hessian = T),
                     error = function(e) e)
  if(inherits(optout, "error")) next # if an error message gets generated this run, move to the next iteration of the loop
  
  output_sp[i,c(4,6)] <- optout$par
  output_sp$Lvalue[i] <- optout$value
  output_sp_hess[[i]] <- optout$hessian
  output_sp$estPar[i] <- 2
  output_sp$scaleT[i] <- 1 # correctly fill in the parameter
}

beep()

##  Calculate variance from Hessian -----------------------------------------------------

# going to need to calculate SEs differently depending on number of parameters fitted

output_sp_hess[which(sapply(output_sp_hess, function(e) sum(e) == 0))] # check any that have "improper" fits giving 0 in hessian
output_sp_hess[which(sapply(output_sp_hess, function(e) 0 %in% e))] # check any that have "improper" fits giving 0 in hessian

notNA2 <- which(output_sp$estPar == 2) # 2 par models
notNA3 <- which(output_sp$estPar == 3) # most were 2 par fits
# reference vector for species with valid Hessian matrices
output_sp_hess[notNA2]
output_sp_hess[notNA3]

output_sp$shape_se <- NA
output_sp$scale_se <- NA
output_sp$lambda_se <- NA
for (i in 1:length(notNA2)) {
  se <- tryCatch(diag(solve(output_sp_hess[[notNA2[i]]])), error = function(e) e) # variances!
  if(inherits(se, "error")) next
  output_sp[notNA2[i],c(9,11)] <- se
}
for (i in 1:length(notNA3)) {
  se <- tryCatch(diag(solve(output_sp_hess[[notNA3[i]]])), error = function(e) e) # variances!
  if(inherits(se, "error")) next
  output_sp[notNA3[i],9:11] <- se
}

# check how many species have singular standard errors (i.e. unreliable models)
output_sp %>% 
  filter(is.na(lambda_se) & is.na(lambdaY) == F) %>% nrow # 6 with invalid hessians
output_sp %>% filter(is.na(lambda_se) & is.na(lambdaY)==F) %>% View # they're all 3 parameter models with low likelihoods

output_sp_hess[which(is.na(output_sp$shapeT) == F & is.na(output_sp$shape_se) == T)] # look at the invalid hessians


## Fit 3: 2 parameter fit for species with shit Hessians -------------------

# consider running a 2 par model for them
# otherwise, if they're already 2 par and still shit, keep blank.

# remove those and rerun
output_sp[is.na(output_sp$lambdaY) == F & is.na(output_sp$lambda_se) & output_sp$estPar == 3, c(4:7, 9:11)] <- NA

if (length(which(is.na(output_sp$shape_se) & is.na(output_sp$estPar) == F)) > 0) {
  for (i in 1:nrow(output_sp)) {
    
    if(is.na(output_sp$shapeT[i]) == F) next # skip if previously estimated
    
    set.seed(240) # replicability
    # vector of detections per camera i at site
    detects <- data_ruv %>% 
      filter(site_ID == output_sp$site_ID[i], Taxon == output_sp$Taxon[i]) %>% 
      group_by(Camera) %>% summarise(detects = sum(Count)) %>% # Counts because some detections are grouped
      pull(detects)
    
    # put that through the optimise function
    # with error handling
    # don't need to specify the 2 par REST function again, just reuse it.
    optout <- tryCatch(optim(par = c(0.5,9), REST_poisson2, method = "L-BFGS-B", lower = c(0,0), upper = c(2000,450), hessian = T),
                       error = function(e) e)
    if(inherits(optout, "error")) next # if an error message gets generated this run, move to the next iteration of the loop
    
    output_sp[i,c(4,6)] <- optout$par
    output_sp$Lvalue[i] <- optout$value
    output_sp_hess[[i]] <- optout$hessian
    output_sp$estPar[i] <- 2
    output_sp$scaleT[i] <- 1 # correctly fill in the parameter
  }
}
beep()

# recalculate variances for those shit outliers
var2 <- which(is.na(output_sp$shapeT) == F & output_sp$estPar == 2)
for (i in var2) {
  se <- tryCatch(diag(solve(output_sp_hess[[i]])), error = function(e) e) # variances!
  if(inherits(se, "error")) next
  output_sp[i,c(9,11)] <- se
}

output_sp_hess[is.na(output_sp$shapeT) == F & is.na(output_sp$shape_se)]

# ok, if it has not converged in three iterations, remove
output_sp <- output_sp %>% filter(!is.na(lambda_se))

# interpret outputs
# gamma shape = mean = E(T) expected staying time
# poisson lambda = mean = variance = E(Y) expected number of detections

# Calculate densities -----------------------------------------------------

# calculate D from E(T) and E(Y), D = E(Y)E(T) / (sH)
# s = detection zone, 4 sq m
# H = observation period, 45 minutes = 2700 s, seconds to line up with

output_sp$D <- with(output_sp, (shapeT * scaleT * lambdaY) / (2700 * 6.75 * cam))
# clear up temporary objects from model fitting
rm(s, se, var2, notNA2, notNA3, i, n, sitecam, detects, j)

# add in number of occurrences
output_sp <- data_ruv %>% group_by(site_ID, Taxon) %>% summarise(occurrences = sum(Count)) %>% full_join(., output_sp, by=c("site_ID", "Taxon"))

# Parametric bootstrapping ------------------------------------------

# make empty lists to store bootstrapping per species
boot_p <- as.list(rep(0,nrow(output_sp)))
# empty dataframe to store the estimated confidence interval bounds
pred_p <- matrix(nrow = nrow(output_sp), ncol = 3)
colnames(pred_p) <- c('pred', 'lwr', 'upr')
set.seed(240) # reproducibility :)
boots = 5000 # how many iterations

for (i in 1:nrow(output_sp)) {

  boot_p[[i]] <- matrix(nrow = boots, ncol = 3) # empty matrix for each bootstrap iteration
  colnames(boot_p[[i]]) <- c('predT', 'predY', 'predD')
  
  # loop for bootstrap iterations per REST model
  # pick out random T and Y from model estimated distributions
  for (j in 1:boots) {
    boot_p[[i]][j,1] <-  rgamma(1, shape = output_sp$shapeT[i], scale = output_sp$scaleT[i])
    boot_p[[i]][j,2] <-  rpois(1, lambda = output_sp$lambdaY[i])
  }
  
  # calculate D from the sampled T and Y parameters
  boot_p[[i]][,3] = (boot_p[[i]][,1] * boot_p[[i]][,2] / (2700 * 6.75 * output_sp$cam[i]))
}
beep()

hist(boot_p[[runif(1, 1, nrow(output_sp))]][,3]) # distribution of a random species' bootstrapped density estimate
# see if it fits log-normal distribution

# calculate confidence interval using percentile method
for (i in 1:nrow(output_sp)) {
  boot_p[[i]] <- boot_p[[i]][order(boot_p[[i]][,3]),] # order ascending values of D
  pred_p[i,1] <- boot_p[[i]][boots/2,3] # 50th percentile for the predicted D
  pred_p[i,2] <- boot_p[[i]][ round(boots*0.025, 1), 3] # 2.5 percentile for lower
  pred_p[i,3] <- boot_p[[i]][ round(boots*0.975, 1), 3] # 97.5 percentile for lower
}

# merge confidence intervals 
output_sp <- bind_cols(output_sp, pred_p)

# calculate AICc so we can compare with negative binomial models
# 2k - 2log(L) + (2k^2 + 2k) / (n - k - 1)

# first need the number of detections

output_sp$AICc <- with(output_sp, 2*estPar - 2*log(Lvalue) + ((2*estPar^2 + 2*estPar) / (occurrences - estPar - 1)))

# export outputs
write.csv(output_sp, "../outputs/REST_species.csv", row.names = F)
