

## FishN: comparing fish assemblage abundance surveying methods
# REST: random encounter staying time modelling (no individual recognition), neg binomial version
# + model comparisons with models fit with Poisson

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
# Neg binomial for estimating expected detections
# Gamma for estimating expected staying time
# requires inputs
# y = vector of detections per camera i
# t = dataframe or matrix of staying times per detection j for camera i, with staying times as rows, cameras in columns
# n = number of cameras

# Neg binomial needs mu = mean, size = dispersion parameter
# Gamma needs shape = mean, scale = ratio between variance and mean. Set scale = 1 by default. 

output_nb <- data.frame(site_ID = rep(c('hale_kaku', 'hale_kinalea'), length(restsp)), Taxon = rep(restsp, each = 2), shapeT = '', scaleT = '', muY = '', sizeY = '')
output_nb_hess <- as.list(rep(0, length(restsp)*2))

# NEST LEVELS: species-site-camera-detection
set.seed(240) # replicability
for (i in 1:length(restsp)) {
  
  for (s in 1:sites) {
    # vector of detections per camera i at site
    detects <- data_ruv %>% filter(site_ID == unique(data_ruv$site_ID)[s], Taxon == restsp[i]) %>% 
      group_by(Camera) %>% summarise(detects = sum(Count)) %>% # Counts because some detections are grouped
      pull(detects)
    
    # Establish MLE likelihood function
    REST_nb <- function(par) { # where x is a vector c(gamma shape, poisson lambda)
      Tij <- rep(0,cams[s])
      for (n in 1:cams[s]) {
        # make a data object for staying time per site-camera
        stay <- data_ruv %>% filter(site_ID == unique(data_ruv$site_ID)[s], Taxon == restsp[i],
                                    Camera == sitecam[sitecam$site_ID == unique(data_ruv$site_ID)[s],2][n]) %>% pull(staytime)
        # do the most nested level first, Tij staying time
        for (x in 1:n) {
          Tij[x] <- dgamma(x = stay, shape = par[1], scale = par[2], log = T) %>% sum
        }
      }
      return(
        -1 * sum( dnbinom(x = detects, mu = par[3], size = par[4], log = T) + sum(Tij) )
      )
    }
    # put that through the optimise function
    # with error handling
    optout <- tryCatch(optim(par = c(0.5,0.5,9,0.1), REST_nb, method = "L-BFGS-B", lower = c(0,0,0,0), upper = c(2000,50,450,50), hessian = T),
                       error = function(e) e)
    if(inherits(optout, "error")) next # if an error message gets generated this run, move to the next iteration of the loop
    
    output_nb[2*i-s+1,3:6] <- optout$par
    output_nb$Lvalue[2*i-s+1] <- optout$value
    output_nb_hess[[2*i-s+1]] <- optout$hessian
  }
}

# for models that couldn't converge/optimise, set scale = 1

for (i in 1:length(restsp)) {

  for (s in 1:sites) {
    if(is.na(output_nb$shapeT[2*i-s+1]) == F) next 
    # move on if there's already a model output there
    
    set.seed(240) # replicability
    # vector of detections per camera i at site
    detects <- data_ruv %>% filter(site_ID == unique(data_ruv$site_ID)[s], Taxon == restsp[i]) %>% 
      group_by(Camera) %>% summarise(detects = sum(Count)) %>% # Counts because some detections are grouped
      pull(detects)
    
    # Establish MLE likelihood function
    REST_nb2 <- function(par) { # where x is a vector c(gamma shape, nb lambda)
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
    optout <- tryCatch(optim(par = c(0.5,9), REST_nb2, method = "L-BFGS-B", lower = c(0,0), upper = c(2000,450), hessian = T),
                       error = function(e) e)
    if(inherits(optout, "error")) next # if an error message gets generated this run, move to the next iteration of the loop
    
    # store fitting outputs
    output_nb[2*i-s+1,c(3,5)] <- optout$par
    output_nb_hess[[2*i-s+1]] <- optout$hessian
    
  }
}

output_nb$scaleT[is.na(output_nb$shapeT) == F & is.na(output_nb$scaleT) == T] <- 1

# remove temp objects
rm(detects, i, n, s, sp, stay, Tij, x)
output_nb[[3]] <- as.numeric(output_nb[[3]])
output_nb[[4]]  <- as.numeric(output_nb[[4]])
output_nb[[5]]  <- as.numeric(output_nb[[5]])

# calculate D from E(T) and E(Y), D = E(Y)E(T) / (sH)
# s = detection zone, 4 sq m
# H = observation period, 45 minutes = 2700 s, seconds to line up with 
output_nb$D[output_nb$site_ID == 'hale_kinalea'] <- (with(output_nb[output_nb$site_ID == 'hale_kinalea', ], shapeT * scaleT * lambdaY) / (2700 * 4 * cams[1])) * 250
output_nb$D[output_nb$site_ID == 'hale_kaku'] <- (with(output_nb[output_nb$site_ID == 'hale_kaku', ], shapeT * scaleT * lambdaY) / (2700 * 4 * cams[2]))*250 # scale to be density within belt transect

# going to need to calculate SEs differently depending on number of parameters fitted

notNA2 <- which(is.na(output_nb$D) == F & output_nb$scaleT == 1) # 2 par models
notNA3 <- which(is.na(output_nb$D) == F & output_nb$scaleT != 1) # 3 par models
# reference vector for species with valid Hessian matrices
output_nb_hess[notNA2]
output_nb_hess[notNA3]
notNA3 <- notNA3[-4] # get rid of 4, all zeroes

output_nb$shape_var <- NA
output_nb$scale_var <- NA
output_nb$lambda_var <- NA
for (i in 1:length(notNA2)) {
  se <- tryCatch(diag(solve(output_nb_hess[[notNA2[i]]])), error = function(e) e) # variances!
  if(inherits(se, "error")) next
  output_nb[notNA2[i],c(7,9)] <- se
}
for (i in 1:length(notNA3)) {
  se <- tryCatch(diag(solve(output_nb_hess[[notNA3[i]]])), error = function(e) e) # variances!
  if(inherits(se, "error")) next
  output_nb[notNA3[i],7:9] <- se
}

# check how many species have singular standard errors (i.e. unreliable models)
output_nb %>% 
  filter(is.na(shape_var), is.na(scale_var), is.na(lambda_var), is.na(D) == F) %>% nrow # 4 species
output_nb %>% filter(is.na(shape_var), is.na(scale_var), is.na(lambda_var), is.na(D) == F) %>% View

# consider running a 2 par model for them

# remove those and rerun
output_nb[is.na(output_nb$shape_var) & is.na(output_nb$scale_var) & is.na(output_nb$lambda_var), 3:6] <- NA
# also remove the hessian matrices
output_nb_hess[is.na(output_nb$shape_var) & is.na(output_nb$scale_var) & is.na(output_nb$lambda_var)] <- 0

for (i in 1:length(restsp)) {
  
  for (s in 1:sites) {
    if(is.na(output_nb$shapeT[2*i-s+1]) == F) next 
    # move on if there's already a model output there
    
    set.seed(240) # replicability
    # vector of detections per camera i at site
    detects <- data_ruv %>% filter(site_ID == unique(data_ruv$site_ID)[s], Taxon == restsp[i]) %>% 
      group_by(Camera) %>% summarise(detects = sum(Count)) %>% # Counts because some detections are grouped
      pull(detects)
    
    # Establish MLE likelihood function
    REST_nb2 <- function(par) { # where x is a vector c(gamma shape, nb lambda)
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
    optout <- tryCatch(optim(par = c(0.5,9), REST_nb2, method = "L-BFGS-B", lower = c(0,0), upper = c(2000,450), hessian = T),
                       error = function(e) e)
    if(inherits(optout, "error")) next # if an error message gets generated this run, move to the next iteration of the loop
    
    # store fitting outputs
    output_nb[2*i-s+1,c(3,5)] <- optout$par
    output_nb$Lvalue[2*i-s+1] <- optout$value
    output_nb_hess[[2*i-s+1]] <- optout$hessian
    
  }
}


output_nb$D[output_nb$site_ID == 'hale_kinalea' & is.na(output_nb$D)] <- (with(output_nb[output_nb$site_ID == 'hale_kinalea' & is.na(output_nb$D), ], shapeT * lambdaY) / (2700 * 4 * cams[1])) * 250
output_nb$D[output_nb$site_ID == 'hale_kaku' & is.na(output_nb$D)] <- (with(output_nb[output_nb$site_ID == 'hale_kaku' & is.na(output_nb$D), ], shapeT * lambdaY) / (2700 * 4 * cams[2]))*250 # scale to be density within belt transect
output_nb$scaleT[is.na(output_nb$shapeT) == F & is.na(output_nb$scaleT) == T] <- 1

notNA2 <- which(is.na(output_nb$D) == F & output_nb$scaleT == 1) # 2 par models
for (i in 1:length(notNA2)) {
  se <- tryCatch(diag(solve(output_nb_hess[[notNA2[i]]])), error = function(e) e) # variances!
  if(inherits(se, "error")) next
  output_nb[notNA2[i],c(7,9)] <- se
}

# interpret outputs
# gamma shape = mean = E(T) expected staying time
# nb lambda = mean = variance = E(Y) expected number of detections

# BOOTSTRAP CONFIDENCE INTERVALS

# make empty lists to store bootstrapping per species
bootTY <- as.list(rep(0,nrow(output_nb)))
# empty dataframe to store the estimated confidence interval bounds
Dpredict <- matrix(nrow = nrow(output_nb), ncol = 3)
colnames(Dpredict) <- c('pred', 'lwr', 'upr')
set.seed(240) # reproducibility :)
for (i in 1:nrow(output_nb)) {
  bootTY[[i]] <- matrix(nrow = 1000, ncol = 3) # empty matrix for each bootstrap iteration
  colnames(bootTY[[i]]) <- c('predT', 'predY', 'predD')
  
  # loop for bootstrap iterations per REST model
  # pick out random T and Y from model estimated distribution parameters
  for (j in 1:1000) {
    bootTY[[i]][j,1] <-  rgamma(1, shape = output_nb$shapeT[i], scale = output_nb$scaleT[i])
    bootTY[[i]][j,2] <-  rpois(1, lambda = output_nb$lambdaY[i])
  }
  
  # calculate D from the sampled T and Y parameters
  if (output_nb$site_ID[i] == 'hale_kaku') { # account for camera number differences
    bootTY[[i]][,3] = (bootTY[[i]][,1] * bootTY[[i]][,2] / (2700 * 4 * 3)) * 250
  } else {bootTY[[i]][,3] = (bootTY[[i]][,1] * bootTY[[i]][,2] / (2700 * 4 * 4)) * 250}
  
  # calculate the predicted values from each bootstrap
  bootTY[[i]] <- bootTY[[i]] %>% as_tibble %>% arrange(predD) %>% as.matrix # arrange predicted values in ascending order
  Dpredict[i,1] <- mean(bootTY[[i]][,3])
  Dpredict[i,2] <- bootTY[[i]][25,3] # because we arranged it, the 25 = the point of 2.5%
  Dpredict[i,3] <- bootTY[[i]][975,3] # 975 = 97.5
}

# produce final output
output_nb <- bind_cols(output_nb, Dpredict)
final_nb <- data_ruv %>% group_by(site_ID, Taxon) %>% summarise(occurrences = sum(Count)) %>% ungroup %>% 
  filter(occurrences >= 10) %>% select(!occurrences) %>% inner_join(., output_nb, by=c("site_ID", "Taxon"))

# export outputs
# write.csv(final_nb, "../outputs/REST_poisson.csv", row.names = F)
