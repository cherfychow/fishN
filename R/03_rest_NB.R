

## FishN: comparing fish assemblage abundance surveying methods
# REST: random encounter staying time modelling (no individual recognition)
# Negative binomial fitting instead of Poisson

# set up working environment
require(dplyr) 
require(lubridate)
require(beepr)
# custom palettes to be extra
source('https://gist.githubusercontent.com/cherfychow/e9ae890fd16f4c86730748c067feee2b/raw/899dcfc1745421cb4e6ba26826b0bfe55fd8ec14/cherulean.R')

# ASSUMES data_ruv is already set up from rest_poisson already!!

## define some global parameters
# site-camera key
sitecam <- distinct(data_ruv, site_ID, Camera)
sites = 2
cams = c(4,3) # in the order of sites shown in unique(data_ruv$site_ID)
# because the model loop runs on data_ruv

# make a reference data frame for species that REST can run on
# criteria: at least 5 detections total at a site
output_nb <- data_ruv %>% group_by(site_ID, Taxon, Size_class) %>% summarise(occurrences = sum(Count)) %>% ungroup %>% 
  filter(occurrences >= 5) %>% select(!occurrences)

# add a col on the number of cameras for the sites
output_nb <- distinct(data_ruv, site_ID, Camera) %>% group_by(site_ID) %>% 
  summarise(cam = n_distinct(Camera)) %>% ungroup() %>% 
  full_join(., output_nb, by="site_ID")

# MLE REST Poisson fitting-------------------------------------------------------------------------

## Fit 1: full 4 parameters -----------------------------------------------------

# establish the likelihood function
# Poisson for estimating expected detections
# Gamma for estimating expected staying time
# requires inputs
# y = vector of detections per camera i
# t = dataframe or matrix of staying times per detection j for camera i, with staying times as rows, cameras in columns
# n = number of cameras

# Neg. binomial needs mu and size. Mu = mmean, and size = dispersion parameter
# Gamma needs shape = mean, scale = ratio between variance and mean. 

output_nb <- mutate(output_nb, shapeT = NA, scaleT = NA, muY = NA, sizeY = NA, Lvalue = NA, estPar = NA) %>% as.data.frame() # add blank cols to fill with par estimates
output_nb_hess <- as.list(rep(0, nrow(output_nb)))

# NEST LEVELS: species-site-camera-detection
for (i in 1:nrow(output_nb)) {
  
  set.seed(240) # replicability
  # vector of detections per camera i at site
  detects <- data_ruv %>% 
    filter(site_ID == output_nb$site_ID[i], Taxon == output_nb$Taxon[i], Size_class == output_nb$Size_class[i]) %>% 
    group_by(Camera) %>% summarise(detects = sum(Count)) %>% # Counts because some detections are grouped
    pull(detects)
  
  # Establish MLE likelihood function
  REST_nb <- function(par) { # where par is a vector c(gamma shape, gamma scale, poisson lambda)
    Tij <- rep(0,output_nb$cam[i]) # empty vector to fill in with Tij estimates per camera
    for (n in 1:output_nb$cam[i]) {
      # make a data object for staying time per site-camera
      stay <- data_ruv %>% filter(site_ID == output_nb$site_ID[i], Taxon == output_nb$Taxon[i], Size_class == output_nb$Size_class[i],
                                  Camera == sitecam[sitecam$site_ID == output_nb$site_ID[i], 2][n]) %>% pull(staytime)
      # do the most nested level first, Tij staying time
      for (x in 1:n) {
        Tij[x] <- dgamma(x = stay, shape = par[1], scale = par[2], log = T) %>% sum # params[1] = shape
      }
    }
    return(
      -1 * ( sum(dnbinom(x = detects, mu = par[3], size = par[4], log = T)) + sum(Tij) ) # params[2] = lambda
    )
  }
  # put that through the optimise function
  # with error handling
  optout <- tryCatch(optim(par = c(0.5,0.5,5,0.5), REST_nb, method = "L-BFGS-B", lower = c(0,0,0,0), upper = c(2000,9,450,20), hessian = T),
                     error = function(e) e)
  if(inherits(optout, "error")) next # if an error message gets generated this run, move to the next iteration of the loop
  
  output_nb[i,5:8] <- optout$par
  output_nb$Lvalue[i] <- optout$value
  output_nb_hess[[i]] <- optout$hessian
  output_nb$estPar[i] <- 4
}
beep()

## Fit 2: 2 parameter (scale = 1) ------------------------------------------

# for models that couldn't converge/optimise, set scale = 1

for (i in 1:nrow(output_nb)) {
  
  if(output_nb$estPar[i] == 4) next # skip if previously estimated
  
  set.seed(240) # replicability
  # vector of detections per camera i at site
  detects <- data_ruv %>% 
    filter(site_ID == output_nb$site_ID[i], Taxon == output_nb$Taxon[i], Size_class == output_nb$Size_class[i]) %>% 
    group_by(Camera) %>% summarise(detects = sum(Count)) %>% # Counts because some detections are grouped
    pull(detects)
  
  # Establish MLE likelihood function
  REST_nb2 <- function(par) { # where par is a vector c(gamma shape, gamma scale, poisson lambda)
    Tij <- rep(0,output_nb$cam[i]) # empty vector to fill in with Tij estimates per camera
    for (n in 1:output_nb$cam[i]) {
      # make a data object for staying time per site-camera
      stay <- data_ruv %>% filter(site_ID == output_nb$site_ID[i], Taxon == output_nb$Taxon[i], Size_class == output_nb$Size_class[i],
                                  Camera == sitecam[sitecam$site_ID == output_nb$site_ID[i], 2][n]) %>% pull(staytime)
      # do the most nested level first, Tij staying time
      for (x in 1:n) {
        Tij[x] <- dgamma(x = stay, shape = par[1], scale = 1, log = T) %>% sum # params[1] = shape
      }
    }
    return(
      -1 * (sum(dnbinom(x = detects, mu = par[2], size = par[3], log = T)) + sum(Tij)) # params[2] = lambda
    )
  }
  # put that through the optimise function
  # with error handling
  optout <- tryCatch(optim(par = c(0.2,5,0.2), REST_nb2, method = "L-BFGS-B", lower = c(0,0,0), upper = c(2000,450,20), hessian = T),
                     error = function(e) e)
  if(inherits(optout, "error")) next # if an error message gets generated this run, move to the next iteration of the loop
  
  output_nb[i,c(5,7,8)] <- optout$par
  output_nb$Lvalue[i] <- optout$value
  output_nb_hess[[i]] <- optout$hessian
  output_nb$estPar[i] <- 3
  output_nb$scaleT[i] <- 1 # correctly fill in the parameter
}

beep()

##  Calculate variance from Hessian -----------------------------------------------------

# going to need to calculate SEs differently depending on number of parameters fitted

output_nb_hess[which(sapply(output_nb_hess, function(e) sum(e) == 0))] # check any that have "improper" fits giving 0 in hessian

output_nb <- mutate(output_nb, shape_se = NA, scale_se = NA, mu_se = NA, size_se = NA)

for (i in 1:nrow(output_nb)) {
  se <- tryCatch(diag(solve(output_nb_hess[[i]])), error = function(e) e) # variances!
  if(inherits(se, "error")) next
  if (output_nb$scaleT[i] == 1) {
    output_nb[i,c(11,13,14)] <- se
  } else {output_nb[i,11:14] <- se}
}

# check how many species have singular standard errors (i.e. unreliable models)
output_nb %>% 
  filter(is.na(shape_se)| is.na(mu_se)) %>% nrow # 24 species
output_nb %>% filter(is.na(shape_se) | is.na(mu_se)) %>% View

output_nb_hess[which(is.na(output_nb$shapeT) == F & is.na(output_nb$shape_se) == T)] # look at the invalid hessians


## Fit 3: 2 parameter fit for species with shit Hessians -------------------

# consider running a 2 par model for them
# otherwise, if they're already 2 par and still shit, keep blank.

# remove those and rerun
# only retry those that were run with 4 parameters
output_nb[is.na(output_nb$shape_se) & output_nb$estPar == 4, 5:8] <- NA

if (length(which(is.na(output_nb$shape_se) & output_nb$estPar == 4)) > 0) {
  for (i in 1:nrow(output_nb)) {
  
  if(is.na(output_nb$shapeT[i]) == F) next # skip if previously estimated
  
  set.seed(240) # replicability
  # vector of detections per camera i at site
  detects <- data_ruv %>% 
    filter(site_ID == output_nb$site_ID[i], Taxon == output_nb$Taxon[i], Size_class == output_nb$Size_class[i]) %>% 
    group_by(Camera) %>% summarise(detects = sum(Count)) %>% # Counts because some detections are grouped
    pull(detects)
  
  # put that through the optimise function
  # with error handling
  optout <- tryCatch(optim(par = c(0.2,5,0.2), REST_nb2, method = "L-BFGS-B", lower = c(0,0,0), upper = c(2000,450,20), hessian = T),
                     error = function(e) e)
  if(inherits(optout, "error")) next # if an error message gets generated this run, move to the next iteration of the loop
  
  output_nb[i,c(5,7,8)] <- optout$par
  output_nb$Lvalue[i] <- optout$value
  output_nb_hess[[i]] <- optout$hessian
  output_nb$estPar[i] <- 3
  output_nb$scaleT[i] <- 1 # correctly fill in the parameter
}

beep()

# recalculate variances for those shit outliers
for (i in which(is.na(output_nb$shape_se) & output_nb$estPar == 3)) {
  se <- tryCatch(diag(solve(output_nb_hess[[i]])), error = function(e) e) # variances!
  if(inherits(se, "error")) next
  output_nb[i,c(11,13,14)] <- se
}}

# ok, if it has not converged in three iterations, remove
output_nb <- output_nb %>% filter(!is.na(mu_se))

# interpret outputs
# gamma shape = mean = E(T) expected staying time
# poisson lambda = mean = variance = E(Y) expected number of detections

# Calculate densities -----------------------------------------------------

# calculate D from E(T) and E(Y), D = E(Y)E(T) / (sH)
# s = detection zone, 4 sq m
# H = observation period, 45 minutes = 2700 s, seconds to line up with

output_nb$D <- with(output_nb, (shapeT * scaleT * muY) / (2700 * 4 * cam)) * 250 # scale by 250 to put density values as per transect
# clear up temporary objects from model fitting
rm(s, se, restsp2, notNA2, notNA, notNA3, i, n, sitecam)

# Bootstrap confidence intervals ------------------------------------------

# make empty lists to store bootstrapping per species
bootTY <- as.list(rep(0,nrow(output_nb)))
# empty dataframe to store the estimated confidence interval bounds
Dpredict <- matrix(nrow = nrow(output_nb), ncol = 2)
colnames(Dpredict) <- c('pred', 'se')
set.seed(240) # reproducibility :)
boots = 5000 # how many iterations

for (i in 1:nrow(output_nb)) {
  
  bootTY[[i]] <- matrix(nrow = boots, ncol = 3) # empty matrix for each bootstrap iteration
  colnames(bootTY[[i]]) <- c('predT', 'predY', 'predD')
  
  # loop for bootstrap iterations per REST model
  # pick out random T and Y from model estimated distributions
  for (j in 1:boots) {
    bootTY[[i]][j,1] <-  rgamma(1, shape = output_nb$shapeT[i], scale = output_nb$scaleT[i])
    bootTY[[i]][j,2] <-  rnbinom(1, mu = output_nb$muY, size = output_nb$sizeY)
  }
  
  # calculate D from the sampled T and Y parameters
  bootTY[[i]][,3] = (bootTY[[i]][,1] * bootTY[[i]][,2] / (2700 * 4 * output_nb$cam[i])) * 250
}
beep()

hist(bootTY[[runif(1, 1, nrow(output_nb))]][,3]) # distribution of a random species' bootstrapped density estimate
# see if it fits log-normal distribution

# will a random species' density distribution look normal if we log it
hist(log(bootTY[[runif(1, 1, nrow(output_nb))]][,3]))

for (i in 1:nrow(output_nb)) {
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
output_nb <- bind_cols(output_nb, Dpredict)

# calculate AICc so we can compare with negative binomial models
# 2k - 2log(L) + (2k^2 + 2k) / (n - k - 1)

# first need the number of detections
output_nb <- data_ruv %>% 
  group_by(site_ID, Taxon, Size_class) %>% summarise(detects = sum(Count)) %>% left_join(output_nb, ., by=c("site_ID", "Taxon", "Size_class")) %>% select(!scale_se) %>% ungroup # Counts because some detections are grouped
output_nb$AICc <- with(output_nb, 2*estPar - 2*log(Lvalue) + ((2*estPar^2 + 2*estPar) / (detects - estPar - 1)))

# export outputs
write.csv(output_nb, "../outputs/REST_nbinom.csv", row.names = F)
