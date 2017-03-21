#########################################################
######   Hierarchical Bayesian Model to understand 
######   how within and between species trait relationships vary
###########################################################



require(rjags)
library(R2jags)
library(lattice)
library(superdiag)
library(mcmcplots)
library(taxize)
library(mvtnorm)
library(lmodel2)

##################################################
####### Load Data ################
##################################################


testspp <- LES$Species[1:5]


t <- tax_name(testspp, get="family", db="ncbi")

which(is.na(t[,3]))

###################################################################
##########      Linear model of w/in species N vs LMA ###########################
###################################################################


##### .Specify Data for model fitting ######################

### get rid of NAs, because
traits.Nclean <- traits.common[-which(is.na(traits.common$log.Nmass)),]

logLMA <- traits.Nclean$log.LMA
logNmass <- traits.Nclean$log.Nmass
species <- as.numeric(traits.Nclean$SP.ID)
# note: there are quite a few sites with a lot of years, but also some sites with very few years

nspecies <- length(unique(species))
n <- length(logNmass)




######### New method for fitting model (from http://www.jkarreth.net/files/bayes-cph_Tutorial-JAGS.pdf) 


##################### .Specify Model ########################
##### first and easiest model: just for well represented species, no Spp mean modeling
mod.NvLMA <- function () {
  
  # Priors (random effects with grand slope and intercept) - not interested in the intercept, mostly just the slope
  for(j in 1:nspecies){  # Slope (b) and intercept (a) for each species
    b.spp[j] ~ dnorm(mu.b, tau.b) #dnorm(mu.b, sigma.b) # slope for each species
    a.spp[j] ~ dnorm(mu.a, tau.a)#dnorm(mu.a, sigma.a) # intercept for each species
  }
  sigma.resid ~ dunif(0,5)  # Residual standard deviation # REML version suggests SD of 0.07477

  # Hyperparameters
  ## -> ML suggests -0.159, w/ spp sd of 0.167
  mu.b ~ dnorm(0, 0.01)  #(sd=10) Mean hyperparameter for random slopes- flat over -2,2 
  sigma.b ~ dunif(0, 5)  # SD hyperparameter for random slopes # flat over great range
  ## -> ML suggests 0.3596 w/ spp sd of 0.33479
  mu.a ~ dnorm(0, 0.01) # Mean hyperparameter for random intercepts
  sigma.a ~ dunif(0,5) # SD hyperparameter for random intercepts
  
  
  # Likelihood
  for(i in 1:n){
    logNmass[i] ~ dnorm(mu.logNmass[i], tau.resid)
    mu.logNmass[i] <- b.spp[species[i]] * logLMA[i] + a.spp[species[i]] # relationship derived by true nirv 
    
  }
  
  # Derived quantities 
  tau.resid <- 1/(sigma.resid*sigma.resid)
  tau.b <- 1/(sigma.b*sigma.b)
  tau.a <- 1/(sigma.a*sigma.a)
  
  # Posterior predictive checks
  for(i in 1:n){
    predictedN[i] <- mu.logNmass[i] # store the fitted value
    resid[i] <- logNmass[i]- mu.logNmass[i] # store the residual
    sq[i] <- pow(resid[i],2)
    # generate perfect data set
    y.new[i] ~ dnorm(mu.logNmass[i], tau.resid)
    sq.new[i] <- pow(y.new[i]- mu.logNmass[i], 2)
  }
  fit <- sum(sq[])
  fit.new <- sum(sq.new[])
  test <- step(fit.new - fit)
  bpvalue <- mean(test)
  

}

# Assemble data into a list
data.list.NvLMA <- list(n=n, logLMA=logLMA, logNmass=logNmass, species=species, nspecies=nspecies)

# define parameters to follow
params.NvLMA <- c("mu.a","mu.b","sigma.a", "sigma.b","fit","fit.new","bpvalue", "sigma.resid")
params.NvLMA2 <- c("mu.a","mu.b","sigma.a", "sigma.b","fit","fit.new","bpvalue", "sigma.resid", "resid",'predictedN',"a.spp","b.spp")


# write an inits function that auto-defines inits for each chain
mod.inits.NvLMA <- function(){
  list("mu.a"=rlnorm(1), "mu.b"=rlnorm(1), "sigma.a" =runif(1,0,4), "sigma.b" =runif(1,0,4), "sigma.resid"=runif(1,0,4))
}

########### ,fit the model ##############
set.seed(531)
mod.fit.NvLMA <- jags(data=data.list.NvLMA, inits = mod.inits.NvLMA
                        , parameters.to.save = params.NvLMA2, n.chains=3
                        , n.iter=5000, n.burnin=2000, n.thin=10
                        , model.file= mod.NvLMA )

# print(mod.fit.NvLMA)
# plot(mod.fit.NvLMA)
# traceplot(mod.fit.NvLMA)

### convert to mcmc object for advanced shit
mod.fit.NvLMA.mcmc <- as.mcmc(mod.fit.NvLMA)


########## .convergence diagnostics ##############
## or super diagnostic all-in-one
superdiag(mod.fit.NvLMA.mcmc, burnin=100) #this won't work for things where you follow fitted and resid.

## or ggplot versions from {mcmcplots}
#denplot(mod.fit.nirverr.mcmc)
#traplot(mod.fit.nirverr.mcmc, parms = c("mu.b", "mu.a"))
mcmcplot(mod.fit.NvLMA.mcmc, parms = c("mu.b","mu.a","sigma.a", "sigma.b", "sigma.resid")) # WHOA SHIT! this makes a html output with ALL the diagnostic plots for each param!
# huh, there is a really long lag for mu.b and mu.a. not sure how important this is...
caterplot(mod.fit.NvLMA.mcmc, parms = c("mu.b", "mu.a", "sigma.a", "sigma.b", "sigma.resid"))
# can also mage plots with {ggmcmc}

###### .pull out model parameters  #############################
mod.params.NvLMA <- summary(mod.fit.NvLMA.mcmc)$statistics[,1]
Predicted <- mod.params.NvLMA[grep('predictedN', names(mod.params.NvLMA))]
Resid <- mod.params.NvLMA[grep('resid', names(mod.params.NvLMA))][1:length(Predicted)] # gotta get rid of 'sigma.resid'
a.est <- mod.params.NvLMA[grep("a.spp", names(mod.params.NvLMA))]
b.est <- mod.params.NvLMA[grep("b.spp", names(mod.params.NvLMA))]
spp.est <- data.frame(SP.ID = levels(traits.Nclean$SP.ID), a=a.est, b=b.est)
global.params <- mod.params.NvLMA[c(grep("mu", names(mod.params.NvLMA)), grep("sigma", names(mod.params.NvLMA)))]
### adding in site data
sitesums <- flux %>% group_by (name) %>% filter (TA_ERA>5) %>% summarise (mnirv = mean(nirv, na.rm=T), mgpp = mean(GPP_NT_VUT_MEAN, na.rm=T), mNDVI= mean(ndvi, na.rm=T), lc = unique(lc), mTA_ERA = mean(TA_ERA, na.rm=T), mmod_gpp = mean(mod_gpp, na.rm=T))
sitesumaries <- merge.data.frame(sitesums, site.est, by.x="name", by.y="site")

##### .model checking ######
plot(Resid~Predicted)
plot(Resid~logLMA)
plot(Resid~traits.Nclean$log.LL); abline(h=0)
plot(Resid~factor(species))
qqp(Resid) # a few small resids aren't normal
qqp(a.est)
qqp(b.est)
plot(a.est~b.est) # super correlated, but probably to be expected. Do I need to model this?


# ##### various plot options
# xyplot(mod.fit.nirverr.mcmc)
# xyplot(mod.fit.nirverr.mcmc, layout=c(2,2), aspect="fill")
# densityplot(mod.fit.nirverr.mcmc)
# autocorr.plot(mod.fit.nirverr.mcmc)
# gelman.plot(mod.fit.nirverr.mcmc)
# geweke.diag(mod.fit.nirverr.mcmc) # don't have a clue what this is showing...
# geweke.plot(mod.fit.nirverr.mcmc)
# ### other diagnostic functions to check out:
# # raftery.diag()
# # heidel.diag()


######################## END: pure w/in spp model ##########################







###################################################################
##########      Linear model of w/in species + spp mean  N vs LMA ###########################
###################################################################


##### .Specify Data for model fitting ######################

### get rid of NAs, because
traits.Nclean <- traits.common[-which(is.na(traits.common$log.Nmass)),]

logLMA <- traits.Nclean$log.LMA
logNmass <- traits.Nclean$log.Nmass
species <- as.numeric(traits.Nclean$SP.ID)
# note: there are quite a few sites with a lot of years, but also some sites with very few years

nspecies <- length(unique(species))
n <- length(logNmass)




######### New method for fitting model (from http://www.jkarreth.net/files/bayes-cph_Tutorial-JAGS.pdf) 


##################### .Specify Model ########################
##### first and easiest model: just for well represented species, no Spp mean modeling
mod.NvLMAmns <- function () {
  
  # Priors (random effects with grand slope and intercept) - not interested in the intercept, mostly just the slope
  for(j in 1:nspecies){  # Slope (b) and intercept (a) for each species
    b.spp[j] ~ dnorm(mu.b, tau.b) #dnorm(mu.b, sigma.b) # slope for each species
    a.spp[j] ~ dnorm(mu.a, tau.a)#dnorm(mu.a, sigma.a) # intercept for each species
  }
  sigma.resid ~ dunif(0,5)  # Residual standard deviation # REML version suggests SD of 0.07477
  
  # Priors for species level means for higher order modeling
  for(k in 1:nspecies){  # mean N and LMA for each species
    muN.spp[k] ~ dnorm(0.2, 0.01) # 0.2 is the species logN grand mean from the LES
    muLMA.spp[k] ~dnorm(2,0.001) # 2ish is the species LMA grand mean from the LES
}
  
  # Hyperparameters
  ## -> ML suggests -0.159, w/ spp sd of 0.167
  mu.b ~ dnorm(0, 0.01)  #(sd=10) Mean hyperparameter for random slopes- flat over -2,2 
  sigma.b ~ dunif(0, 5)  # SD hyperparameter for random slopes # flat over great range
  ## -> ML suggests 0.3596 w/ spp sd of 0.33479
  mu.a ~ dnorm(0, 0.01) # Mean hyperparameter for random intercepts
  sigma.a ~ dunif(0,5) # SD hyperparameter for random intercepts
  ## for the grand slope across spp means
  mu.sppmean ~ dnorm(0,0.01) # grand across species mean relationship
  sigma.sppmean ~ dunif(0,5) #assuming sigma across species means will be pretty reasonable (like <1)
  
  # Likelihood
  for(i in 1:n){
    logNmass[i] ~ dnorm(mu.logNmass[i], tau.resid)
    mu.logNmass[i] <- b.spp[species[i]] * logLMA[i] + a.spp[species[i]] # relationship derived by true nirv 
    logNmass[i] ~
    
  }
  
  # Derived quantities 
  tau.resid <- 1/(sigma.resid*sigma.resid)
  tau.b <- 1/(sigma.b*sigma.b)
  tau.a <- 1/(sigma.a*sigma.a)
  
  # Posterior predictive checks
  for(i in 1:n){
    predictedN[i] <- mu.logNmass[i] # store the fitted value
    resid[i] <- logNmass[i]- mu.logNmass[i] # store the residual
    sq[i] <- pow(resid[i],2)
    # generate perfect data set
    y.new[i] ~ dnorm(mu.logNmass[i], tau.resid)
    sq.new[i] <- pow(y.new[i]- mu.logNmass[i], 2)
  }
  fit <- sum(sq[])
  fit.new <- sum(sq.new[])
  test <- step(fit.new - fit)
  bpvalue <- mean(test)
  
  
}

# Assemble data into a list
data.list.NvLMAmns <- list(n=n, logLMA=logLMA, logNmass=logNmass, species=species, nspecies=nspecies)

# define parameters to follow
params.NvLMAmns <- c("mu.a","mu.b","sigma.a", "sigma.b","fit","fit.new","bpvalue", "sigma.resid")
params.NvLMAmns2 <- c("mu.a","mu.b","sigma.a", "sigma.b","fit","fit.new","bpvalue", "sigma.resid", "resid",'predictedN',"a.spp","b.spp")


# write an inits function that auto-defines inits for each chain
mod.inits.NvLMAmns <- function(){
  list("mu.a"=rlnorm(1), "mu.b"=rlnorm(1), "sigma.a" =runif(1,0,4), "sigma.b" =runif(1,0,4), "sigma.resid"=runif(1,0,4))
}

########### ,fit the model ##############
set.seed(531)
mod.fit.NvLMAmns <- jags(data=data.list.NvLMAmns, inits = mod.inits.NvLMAmns
                      , parameters.to.save = params.NvLMAmns2, n.chains=3
                      , n.iter=5000, n.burnin=2000, n.thin=10
                      , model.file= mod.NvLMAmns )

# print(mod.fit.NvLMAmns)
# plot(mod.fit.NvLMAmns)
# traceplot(mod.fit.NvLMAmns)

### convert to mcmc object for advanced shit
mod.fit.NvLMAmns.mcmc <- as.mcmc(mod.fit.NvLMAmns)


########## .convergence diagnostics ##############
## or super diagnostic all-in-one
superdiag(mod.fit.NvLMAmns.mcmc, burnin=100) #this won't work for things where you follow fitted and resid.

## or ggplot versions from {mcmcplots}
#denplot(mod.fit.nirverr.mcmc)
#traplot(mod.fit.nirverr.mcmc, parms = c("mu.b", "mu.a"))
mcmcplot(mod.fit.NvLMAmns.mcmc, parms = c("mu.b","mu.a","sigma.a", "sigma.b", "sigma.resid")) # WHOA SHIT! this makes a html output with ALL the diagnostic plots for each param!
# huh, there is a really long lag for mu.b and mu.a. not sure how important this is...
caterplot(mod.fit.NvLMAmns.mcmc, parms = c("mu.b", "mu.a", "sigma.a", "sigma.b", "sigma.resid"))
# can also mage plots with {ggmcmc}

###### .pull out model parameters  #############################
mod.params.NvLMAmns <- summary(mod.fit.NvLMAmns.mcmc)$statistics[,1]
Predicted <- mod.params.NvLMAmns[grep('predictedN', names(mod.params.NvLMAmns))]
Resid <- mod.params.NvLMAmns[grep('resid', names(mod.params.NvLMAmns))][1:length(Predicted)] # gotta get rid of 'sigma.resid'
a.est <- mod.params.NvLMAmns[grep("a.spp", names(mod.params.NvLMAmns))]
b.est <- mod.params.NvLMAmns[grep("b.spp", names(mod.params.NvLMAmns))]
spp.est <- data.frame(SP.ID = levels(traits.Nclean$SP.ID), a=a.est, b=b.est)
global.params <- mod.params.NvLMAmns[c(grep("mu", names(mod.params.NvLMAmns)), grep("sigma", names(mod.params.NvLMAmns)))]
### adding in site data
sitesums <- flux %>% group_by (name) %>% filter (TA_ERA>5) %>% summarise (mnirv = mean(nirv, na.rm=T), mgpp = mean(GPP_NT_VUT_MEAN, na.rm=T), mNDVI= mean(ndvi, na.rm=T), lc = unique(lc), mTA_ERA = mean(TA_ERA, na.rm=T), mmod_gpp = mean(mod_gpp, na.rm=T))
sitesumaries <- merge.data.frame(sitesums, site.est, by.x="name", by.y="site")

##### .model checking ######
plot(Resid~Predicted)
plot(Resid~logLMA)
plot(Resid~traits.Nclean$log.LL); abline(h=0)
plot(Resid~factor(species))
qqp(Resid) # a few small resids aren't normal
qqp(a.est)
qqp(b.est)
plot(a.est~b.est) # super correlated, but probably to be expected. Do I need to model this?







###################################################################
###### CORRELATION MULTILEVEL model w/t-distributed error ####################
#------------------------------------------------------------------

#pulled from: https://www.r-bloggers.com/bayesian-estimation-of-correlation-now-robust/

##### .designate Data #################
# lma <- LESspp$log.LMA[-which(is.na(LESspp$log.LMA) | is.na(LESspp$log.Nmass))]
# Nmass <- LESspp$log.Nmass[-which(is.na(LESspp$log.LMA) | is.na(LESspp$log.Nmass))]
lma <- traits.common$log.LMA_PSA[-which(is.na(traits.common$log.LMA_PSA) | is.na(traits.common$log.Nmass))]
Nmass <- traits.common$log.Nmass[-which(is.na(traits.common$log.LMA_PSA) | is.na(traits.common$log.Nmass))]
species <- as.numeric(traits.common$SP.ID[-which(is.na(traits.common$log.LMA_PSA) | is.na(traits.common$log.Nmass))])

###### .define model ##################
mod.robust_corr <- "
  var mu[nspp,2]
  model {

for(i in 1:n) {
  # We've replaced dmnorm with and dmt ...
  x[i,1:2] ~ dmnorm(mu[1,1:2], inverse(cov[(1*2-1):(1*2),1:2]))
      # likelihood function is just a t distribution with a mean, precision, and nu=degrees of freedom 
  }
for (j in 1:nspp){
  #row1 <- (i * 2) -1
  #row2 <- (i * 2) 
  sigma.spp[j,1] ~ dgamma(sh.spplma, ra.spplma) # within-species sd in lma
  sigma.spp[j,2] ~ dgamma(sh.sppNmass, ra.sppNmass) # within-species sd in Nmass
  rho.betascale[j] ~  dbeta(a.spp, b.spp) # within-species correlation between lma and Nmass
  rho[j] <- (rho.betascale[j] * 2) -1
  ## put these together into a covariance matrix
  cov[j*2-1,1] <- sigma.spp[j,1] * sigma.spp[j,1] #lma variance for each species
  cov[j*2,2] <- sigma.spp[j,2] * sigma.spp[j,2] # Nmass variance for each species
  cov[j*2-1,2] <- sigma.spp[j,1] * sigma.spp[j,2] * rho[j] # lma~Nmass covariance for each species
  cov[j*2,1] <- sigma.spp[j,1] * sigma.spp[j,2] * rho[j]
  mu[j,1:2] ~ dmnorm(mu.btwspp[], prec.btwspp[,])
}

# E[gamma(a,b)] = a/b
# Var[gamma(a,b)] = a/b2
# gamma is the conjugate prior for tau.



 #### Hyperpriors for Rho #########

## I think, if I use a beta distribtuion * 2) -1, I should have a random variable that varies between -1, and 1
a.spp <- ((1-mu.5rho)/var.5rho - 1/mu.5rho) * mu.5rho ^2 #parameterized mu.5rho is the global rho *2 -1
b.spp <- a.spp * (1/mu.5rho - 1)
mu.5rho <- (mu.rho + 1 )/ 2 # transformed mean for correlation in beta-space [0,1] rather than correlation space [-1,1]
mu.rho ~ dunif(-1,1) # diffuse prior over mu.rho
var.5rho ~ dunif(0,.25) # diffuse prior over var.5rho (variance of the beta transformed var.rho)
    # NOTE: I think, I'll just calculate the variance of rho post hoc by looking at the variance of rho estimates per species...


 ##### Hyperpriors for LMA and Nmass variances ######
sh.spplma <- pow(global.sigmalma,2) / var.sigmalma # shape parameter from the global variance in lma sd
ra.spplma <- global.sigmalma / var.sigmalma # rate parameter from the global mean sd for lma

sh.sppNmass <- pow(global.sigmaNmass,2) / var.sigmaNmass # shape parameter from the global variance in lma sd
ra.sppNmass <- global.sigmaNmass / var.sigmaNmass # rate parameter from the global mean sd for lma


# priors for the hyperperameters on lma an Nmass sd
global.sigmalma ~ dunif(0,10) # mean global sd for lma
var.sigmalma ~ dunif(0,10) # variance in species level sd of lma

global.sigmaNmass ~ dunif(0,10) # mean global sd for lma
var.sigmaNmass ~ dunif(0,10) # variance in species level sd of lma


######## Hyperparameters on the between species LMA and Nmass correlation #############
# now calculate the between species mean correlations from the species mu[] vector ########

#### precision parameter is the inverse of the covariance matrix between x1 and x2
prec.btwspp[1:2,1:2] <- inverse(cov.btwspp[,])

cov.btwspp[1,1] <- sigma.btw[1] * sigma.btw[1] # variance of btw spp LMA
cov.btwspp[1,2] <- sigma.btw[1] * sigma.btw[2] * rho.btw # covariance between LMA and Nmass btw spp
cov.btwspp[2,1] <- sigma.btw[1] * sigma.btw[2] * rho.btw 
cov.btwspp[2,2] <- sigma.btw[2] * sigma.btw[2] # variance of btw spp Nmass 
#L66
###### HyperPriors on between spp means ####
sigma.btw[1] ~ dunif(0, 10) # these are very large. could probably bring them down to 10 with no problem
sigma.btw[2] ~ dunif(0, 10) # might also throw in a dgamma to help things converge
rho.btw ~ dunif(-1, 1)
mu.btwspp[1] ~ dnorm(0, 0.001)
mu.btwspp[2] ~ dnorm(0, 0.001)

# Generate predicted dataset using random draws from the bivariate t-distr
# x_rand ~ dmt(mu[], prec[ , ], nu)
}
"
#_____________________________________________________________________________________


##########
#### testing out how the beta works for making priors on the w/in species variances

# r.5 <- .75
# var.r.5 <-0.04  # upper bound on the variance is 0.25., if variance is 0.04, gets real skewed. if variance is 0.025, skew begins to decrease.
# a<- ((1-r.5)/var.r.5 - 1/r.5) * r.5 ^2
# b<- a * (1/r.5 - 1)
#   
# test <- rbeta(n = 1000, shape1 = a, shape2=b)
# mean(test*2 -1)
# var(test*2-1)
 # so the mean scales, but the variance does not. Instead I have to take the variance of the scaled values

## gamma distribution for variances
global.sigma <- 10
var.sigma <- 10
sh <- pow(global.sigma,2) / var.sigma # shape parameter from the global variance in lma sd
ra <- global.sigma / var.sigma # rate parameter from the global mean sd for lma
tmp <- rgamma(n=1000, shape= sh, rate=ra)
hist(tmp, breaks=20)
hist(tmp, breaks=100)#, xlim=c(0,1))
hist(spsd$vlm, add=T, col="red")
spsd <- traits.common %>% group_by(SP.ID) %>% summarise(vlm=sd(log.LMA_PSA,na.rm=T), vn =sd(log.Nmass,na.rm=T))
# sp the mean sd appears to be 0.0895 for lma and 0.0797 for Nmass
# the between spp sd in sd is 0.026 and 0.022

## NOTE: uniform priors on global.sigma and var.sigma result in some super non-uniform distributions for individual species sigmas.


########## .Setting data and inits ########

data_list = list(x = cbind(lma, Nmass), n = length(lma), nspp = length(unique(species)), species=species)
# Use robust estimates of the parameters as initial values
inits_list = list(mu.btwspp = c(rnorm(1), rnorm(1)), rho.btw = runif(1,min = -1,1)
                  , sigma.btw = c(runif(1,0,10), runif(1,0,10))
                  , var.sigmaNmass = runif(1,0,1), global.sigmaNmass = runif(1,0,1)
                  , var.sigmalma = runif(1,0,1), global.sigmalma = runif(1,0,1)
                  , mu.rho = runif(1,-1,1) 
                  , var.5rho = runif(1,0,.25)
)

###########
jags_model <- jags.model(textConnection(mod.robust_corr), data = data_list,
                         inits = inits_list, n.adapt = 500, n.chains = 3, quiet = FALSE)



mod.fit.NvLMAmns <- jags(data=data.list.NvLMAmns, inits = mod.inits.NvLMAmns
                         , parameters.to.save = params.NvLMAmns2, n.chains=3
                         , n.iter=5000, n.burnin=2000, n.thin=10
                         , model.file= mod.NvLMAmns )

update(jags_model, 500)
mcmc_samples <- coda.samples(jags_model, c("mu", "rho", "sigma", "nu", "x_rand"),
                             n.iter = 500)


par(mfrow = c(8, 2), mar = rep(2, 4))
plot(mcmc_samples, auto.layout = FALSE)

samples_mat <- as.matrix(mcmc_samples)

plot(Nmass~lma)
dataEllipse(samples_mat[, c("x_rand[1]", "x_rand[2]")], levels = c(0.5, 0.95),
            plot.points = FALSE)



############### END: LMA v Nmass correlation model ###############################################




###################################################################
###### CORRELATION model w/t-distributed error ####################
#------------------------------------------------------------------

#pulled from: https://www.r-bloggers.com/bayesian-estimation-of-correlation-now-robust/

##### .designate Data #################
lma <- LESspp$log.LMA[-which(is.na(LESspp$log.LMA) | is.na(LESspp$log.Nmass))]
Nmass <- LESspp$log.Nmass[-which(is.na(LESspp$log.LMA) | is.na(LESspp$log.Nmass))]

###### .define model ##################
mod.robust_corr <- "
model {
for(i in 1:n) {
# We've replaced dmnorm with and dmt ...
x[i,1:2] ~ dmt(mu[], prec[ , ], nu)
# likelihood function is just a t distribution with a mean, precision, and nu=degrees of freedom 
}

#### precision parameter is the inverse of the covariance matrix between x1 and x2
prec[1:2,1:2] <- inverse(cov[,])

cov[1,1] <- sigma[1] * sigma[1] # variance of X1
cov[1,2] <- sigma[1] * sigma[2] * rho # on-diagonal components
cov[2,1] <- sigma[1] * sigma[2] * rho # could probably make this a random effect
cov[2,2] <- sigma[2] * sigma[2] # variance of X2

###### Priors ####
sigma[1] ~ dunif(0, 10) # these are very large. could probably bring them down to 10 with no problem
sigma[2] ~ dunif(0, 10) # might also throw in a dgamma to help things converge
rho ~ dunif(-1, 1)
mu[1] ~ dnorm(0, 0.001)
mu[2] ~ dnorm(0, 0.001)

# ... and added a prior on the degree of freedom parameter nu.
nu <- nuMinusOne+1
nuMinusOne ~ dexp(1/29)

# Generate predicted dataset using random draws from the bivariate t-distr
x_rand ~ dmt(mu[], prec[ , ], nu)
}
"

########## .Setting data and inits ########

data_list = list(x = cbind(lma, Nmass), n = length(lma))
# Use robust estimates of the parameters as initial values
inits_list = list(mu = c(median(lma), median(Nmass)), rho = cor(lma
                                                                , Nmass, method = "spearman"), sigma = c(mad(lma), mad(Nmass)))
jags_model <- jags.model(textConnection(mod.robust_corr), data = data_list,
                         inits = inits_list, n.adapt = 500, n.chains = 3, quiet = FALSE)
update(jags_model, 500)
mcmc_samples <- coda.samples(jags_model, c("mu", "rho", "sigma", "nu", "x_rand"),
                             n.iter = 500)


par(mfrow = c(8, 2), mar = rep(2, 4))
plot(mcmc_samples, auto.layout = FALSE)

samples_mat <- as.matrix(mcmc_samples)

plot(Nmass~lma)
dataEllipse(samples_mat[, c("x_rand[1]", "x_rand[2]")], levels = c(0.5, 0.95),
            plot.points = FALSE)



############### END: LMA v Nmass correlation model ###############################################










#------------------------------------------------------------------------------
################ CORRELATION 3 way multivariate t correlation ################
#_______________________________________________________________________________


##### .designate Data #################
lma <- LESspp$log.LMA[-which(is.na(LESspp$log.LMA) | is.na(LESspp$log.Nmass) | is.na(LESspp$log.LL))]
Nmass <- LESspp$log.Nmass[-which(is.na(LESspp$log.LMA) | is.na(LESspp$log.Nmass)| is.na(LESspp$log.LL))]
LL <- LESspp$log.LL[-which(is.na(LESspp$log.LMA) | is.na(LESspp$log.Nmass)| is.na(LESspp$log.LL))]

###### .define model ##################
mod.robust_corr_all <- "
model {
for(i in 1:n) {
# We've replaced dmnorm with and dmt ...
  x[i, 1:3] ~ dmnorm(mu[], prec[,])
# x[i,1:3] ~ dmt(mu[], prec[ , ], nu)
# likelihood function is just a t distribution with a mean, precision, and nu=degrees of freedom 
}

#### precision parameter is the inverse of the covariance matrix between x1 and x2
prec[1:3,1:3] <- inverse(cov[,])

cov[1,1] <- sigma[1] * sigma[1] # variance of X1
cov[1,2] <- sigma[1] * sigma[2] * rho12 # off diagonal components
cov[1,3] <- sigma[1] * sigma[3] * rho13
cov[2,1] <- sigma[1] * sigma[2] * rho12 # could probably make this a random effect
cov[2,2] <- sigma[2] * sigma[2] # variance of X2
cov[2,3] <- sigma[2] * sigma[3] * rho23
cov[3,1] <- sigma[1] * sigma[3] * rho13
cov[3,2] <- sigma[2] * sigma[3] * rho23
cov[3,3] <- sigma[3] * sigma[3]

###### Priors ####
sigma[1] ~ dunif(0.0001, 5) # these are very large. could probably bring them down to 10 with no problem
sigma[2] ~ dunif(0.0001, 5) # might also throw in a dgamma to help things converge
sigma[3] ~ dunif(0.0001, 5)
rho12 ~ dunif(-.80, .10)
rho13 ~ dunif(-.10, .80)
rho23 ~ dunif(rho12*rho13-sqrt((1-rho12^2)*(1-rho13^2)), rho12*rho13+sqrt((1-rho12^2)*(1-rho13^2)))
mu[1] ~ dnorm(0, 0.01)
mu[2] ~ dnorm(0, 0.01)
mu[3] ~ dnorm(0, 0.01)

# ... and added a prior on the degree of freedom parameter nu.
# nu <- nuMinusOne+1
# nuMinusOne ~ dexp(1/29)

# Generate predicted dataset using random draws from the bivariate t-distr
#x_rand ~ dmt(mu[], prec[ , ], nu)
x_rand ~ dmnorm(mu[], prec[ , ])
}
"

########## .Setting data and inits ########

data_list = list(x = cbind(lma, Nmass, LL), n = length(lma))
# Use robust estimates of the parameters as initial values
inits_list = list(mu = c(median(lma), median(Nmass), median(LL))
                  , rho12 = cor(lma , Nmass, method = "spearman")
                  , rho13 = cor(lma , LL, method = "spearman")
                  , rho23 = cor(Nmass, LL, method = "spearman")
                  , sigma = c(mad(lma), mad(Nmass), mad(LL)))

inits_list = list(mu = c(median(lma), median(Nmass), median(LL))
                  , rho12 = .1
                  , rho13 = .1
                  , rho23 = .1
                  , sigma = c(.1,.3,.1))

jags_model_all <- jags.model(textConnection(mod.robust_corr_all), data = data_list,
                         inits = inits_list, n.adapt = 500, n.chains = 3, quiet = FALSE)
update(jags_model_all, 500)
mcmc_samples_all <- coda.samples(jags_model_all, c("mu", "rho12", "rho13", "rho23", "sigma", "nu", "x_rand"),
                             n.iter = 500)


par(mfrow = c(8, 2), mar = rep(2, 4))
plot(mcmc_samples, auto.layout = FALSE)

samples_mat_all <- as.matrix(mcmc_samples_all)

plot(Nmass~lma)
dataEllipse(samples_mat[, c("x_rand[1]", "x_rand[2]")], levels = c(0.5, 0.95),
            plot.points = FALSE)






#______________________________________________________________________
############### New 3way correlation model with Wishart prior ##############
#____________________________________________________________________

# example of using a wishart prior on the correlation matrix
# from: https://sourceforge.net/p/mcmc-jags/discussion/610037/thread/eab372de/

mod.corr_all <- "model {
  
  # Priors on the inverse of the correlation matrix
  invW[1:3,1:3] ~ dwish(I[,],4) # df=4 (J+1) to get a uniform prior on correlation
  
  # Scaling matrix V (in term of its inverse)
  invV[1,1] <- 1/lambda[1]
  invV[2,2] <- 1/lambda[2]
  invV[3,3] <- 1/lambda[3]
  invV[1,2] <- 0
  invV[2,1] <- 0
  invV[3,1] <- 0
  invV[1,3] <- 0
  invV[2,3] <- 0
  invV[3,2] <- 0
  
  # Prior on scaling coefficients (using half normal priors here...)
  # I think these are variances? 
  lambda[1] ~ dnorm(0,0.0001)T(0,)
  lambda[2] ~ dnorm(0,0.0001)T(0,)
  lambda[3] ~ dnorm(0,0.0001)T(0,)
  
  # Priors on mu
  mu[1] ~ dnorm(0, 0.01)
  mu[2] ~ dnorm(0, 0.01)
  mu[3] ~ dnorm(0, 0.01)

  # Precision matrix of multivariate node
  invSigma[1:3,1:3] <- invV[,]%*%invW[,]%*%invV[,]
  
  # Covariance parameters # this is the covariance matrix
  Sigma[1:3,1:3] <- inverse(invSigma[,])
  var1 <- Sigma[1,1] # variance 1
  var2 <- Sigma[2,2] # variance 2
  var3 <- Sigma[3,3] # variance 3
  sd[1] <- sqrt(var1) # standard deviation 1 
  sd[2] <- sqrt(var2) # standard deviation 2
  sd[3] <- sqrt(var3) # standard deviation 3
  rho[1] <- Sigma[1,2]/sqrt(var1*var2)
  rho[2] <- Sigma[1,3]/sqrt(var1*var3)
  rho[3] <- Sigma[2,3]/sqrt(var2*var3)
  
  # Multivariate node
  for(i in 1:N){
    x[i,1:3] ~ dmnorm(mu[],invSigma[,])
  }
  
  # Posterior predictive checks
  x_rand ~ dmnorm(mu[], invSigma[ , ])
  # for(i in 1:n){
  #   predicted[i,1:3] <- dmnorm(mu[], invSigma[,]) # store the fitted value
  #   resid[i,1:3] <- predicted[i,1:3]- mu.logNmass[i] # store the residual
  #   sq[i] <- pow(resid[i],2)
  #   # generate perfect data set
  #   y.new[i] ~ dnorm(mu.logNmass[i], tau.resid)
  #   sq.new[i] <- pow(y.new[i]- mu.logNmass[i], 2)
  # }
  # fit <- sum(sq[])
  # fit.new <- sum(sq.new[])
  # test <- step(fit.new - fit)
  # bpvalue <- mean(test)
  
  
}"



# Assemble data into a list
data.corr_all <- list(N=length(lma), I=identityMatrix()[1:3,1:3], x=cbind(lma,Nmass,LL))

# define parameters to follow
params.corr_all <- c("mu", "rho", "sd","x_rand")


# write an inits function that auto-defines inits for each chain
inits.corr_all <- function(){
  list(mu=c(median(lma), median(Nmass),median(LL)), rho=c(0,0,0))
}

########### ,fit the model ##############
set.seed(531)
fit.corr_all <- jags(data=data.corr_all
                         , parameters.to.save = params.corr_all, n.chains=3
                         , n.iter=5000, n.burnin=2000, n.thin=10
                         , model.file= textConnection(mod.corr_all) )

##### This fails because of some sampling thing. Maybe it's in the inits?


#### second try from the Adrian Someone's UNC course
mod.corr_all2 <- "model{
  for(i in 1:n) {
    y[i,1:3]~dmnorm(B.hat[],Tau.B[,])
  } 

    B.hat[1]<-mu.a
    B.hat[2]<-mu.b
    B.hat[3]<-mu.c
  mu.a~dnorm(0,.0001)
  mu.b~dnorm(0,.0001)
  mu.c~dnorm(0,.0001)
  Tau.B[1:3,1:3] ~ dwish(W[,],4)
  sigma.raw[1:3,1:3] <- inverse(Tau.B[,])
  #for comparison with model with separate priors
  sigma.a <- sqrt(sigma.raw[1,1])
  sigma.b <- sqrt(sigma.raw[2,2])
  sigma.c <- sqrt(sigma.raw[3,3])
  rho.ab <- sigma.raw[2,1]/ (sigma.a*sigma.b)
  rho.ac <- sigma.raw[3,1]/ (sigma.a*sigma.c)
  rho.bc <- sigma.raw[2,3]/ (sigma.b*sigma.c)
}"

# Assemble data into a list
data.corr_all <- list(n=length(lma), W=diag(3)*c(mad(lma), mad(Nmass), mad(LL)), y=cbind(lma,Nmass,LL))

# define parameters to follow
params.corr_all <- c("mu.a","mu.b","mu.c", "rho.ab","rho.ac","rho.bc", "sigma.a","sigma.b","sigma.c")

require(MCMCpack)
# write an inits function that auto-defines inits for each chain
inits.corr_all <- function(){
  list(mu.a=rnorm(1), mu.b=rnorm(1), mu.c=rnorm(1), Tau.B=rwish(3,diag(3)))
}

########### ,fit the model ##############
set.seed(531)
fit.corr_all <- jags(data=data.corr_all, inits = inits.corr_all
                     , parameters.to.save = params.corr_all, n.chains=3
                     , n.iter=5000, n.burnin=2000, n.thin=10
                     , model.file= textConnection(mod.corr_all2) )
# It actually initialized!!!

fit.corr_all.mcmc <- as.mcmc(fit.corr_all)

print(fit.corr_all)
## Check Rhat
max(fit.corr_all$BUGSoutput$summary[,"Rhat"]) ## check Rhat
## check neff
min(fit.corr_all$BUGSoutput$summary[,"n.eff"]) ## not too much compression?
## some plots
mcmcplot(fit.corr_all.mcmc)




########### TRYING TO ADD A LEVEL!!!m (incomplete, does not run) #####################
# I can't figure out how to draw the variance-covariance matrix from a global variance-covariance matrix while making sure that all of them are positive definite.
#### second try from the Adrian Someone's UNC course
mod.corr_spp <- "model{
for(j in 1:nspp){
  
}

for(i in 1:n) {
y[i,1:3]~dmnorm(B.hat[species[i], 1:3],Tau.B[,])
} 

B.hat[1]<-mu.a
B.hat[2]<-mu.b
B.hat[3]<-mu.c
mu.a~dnorm(0,.0001)
mu.b~dnorm(0,.0001)
mu.c~dnorm(0,.0001)
Tau.B[1:3,1:3] ~ dwish(W[,],4)
sigma.raw[1:3,1:3] <- inverse(Tau.B[,])
#for comparison with model with separate priors
sigma.a <- sqrt(sigma.raw[1,1])
sigma.b <- sqrt(sigma.raw[2,2])
sigma.c <- sqrt(sigma.raw[3,3])
rho.ab <- sigma.raw[2,1]/ (sigma.a*sigma.b)
rho.ac <- sigma.raw[3,1]/ (sigma.a*sigma.c)
rho.bc <- sigma.raw[2,3]/ (sigma.b*sigma.c)
}"

# Assemble data into a list
data.corr_all <- list(n=length(lma), W=diag(3)*c(mad(lma), mad(Nmass), mad(LL)), y=cbind(lma,Nmass,LL))

# define parameters to follow
params.corr_all <- c("mu.a","mu.b","mu.c", "rho.ab","rho.ac","rho.bc", "sigma.a","sigma.b","sigma.c")

require(MCMCpack)
# write an inits function that auto-defines inits for each chain
inits.corr_all <- function(){
  list(mu.a=rnorm(1), mu.b=rnorm(1), mu.c=rnorm(1), Tau.B=rwish(3,diag(3)))
}

########### ,fit the model ##############
set.seed(531)
fit.corr_all <- jags(data=data.corr_all, inits = inits.corr_all
                     , parameters.to.save = params.corr_all, n.chains=3
                     , n.iter=5000, n.burnin=2000, n.thin=10
                     , model.file= textConnection(mod.corr_all2) )
# It actually initialized!!!

fit.corr_all.mcmc <- as.mcmc(fit.corr_all)

print(fit.corr_all)
## Check Rhat
max(fit.corr_all$BUGSoutput$summary[,"Rhat"]) ## check Rhat
## check neff
min(fit.corr_all$BUGSoutput$summary[,"n.eff"]) ## not too much compression?
## some plots
mcmcplot(fit.corr_all.mcmc)














############################################################
############ **Maximum Likelihood easy models** ###############################################################
############################################################



########### Breaking Data up by Taxonomic Level#####################
levels(LES$Needle.Broad.lf) <- list(B="B",N="N",unknown="")

### common spp and genera and families in GLOPNET
commonspp <- names(which(xtabs(~Species, LES)>=5)) # 6 have >5 records, 14 have >4 records, 36 have >3 records
commongen <- names(which(xtabs(~Genus, LES)>=5))
commonfam <- names(which(xtabs(~Family, LES)>=5))


#### averaging things up to spp level ########
LESspp <- LES %>% group_by(Species, GE.SP, Genus, Family) %>% summarise(GrForm = unique(GF)[1], DecidEver = unique(Decid.E.green)[1], NeedleBroad = unique(Needle.Broad.lf)[1], C3.C4 = unique(C3C4)[1], slog.LL = mean(log.LL, na.rm=T), slog.LMA = mean(log.LMA, na.rm=T), slog.Nmass = mean(log.Nmass, na.rm=T),slog.Narea = mean(log.Narea, na.rm=T)
                                                                        ,rslog.LL = log(mean(10^log.LL, na.rm=T),base=10), rslog.LMA = log(mean(10^log.LMA, na.rm=T),base=10), rslog.Nmass = log(mean(10^log.Nmass, na.rm=T),base=10) )
# have to rename things slog to keep for summarise, but want them just as log
colnames(LESspp) <- gsub("slog", "log", colnames(LESspp))
# xtabs(~Species, LES)[which(xtabs(~Species, LES)>1)]
  # ~500 entries are doubled species. not many of them have more than 5 entries though. 
  # evidently some of those have conflicting binary IDs
### code for cleaning out problem genera that have multiple families
# LESgenprobs <- LESspp %>% group_by(Genus) %>% summarise(nFam= length(unique(na.omit(Family))))
# prob.gen <- LESgenprobs$Genus[which(LESgenprobs$nFam>1)]                                                   
# tmp <- LESspp[which(LESspp$Genus %in% prob.gen),]
traitsspp <- traits %>% group_by(GE.SP, GENUS, Family) %>% summarise(slog.LL = mean(log.LL, na.rm=T), slog.LMA = mean(log.LMA, na.rm=T), slog.Nmass = mean(log.Nmass, na.rm=T),slog.Narea = mean(log.Narea, na.rm=T)
                                                                        ,rslog.LL = log(mean(10^log.LL, na.rm=T),base=10), rslog.LMA = log(mean(10^log.LMA, na.rm=T),base=10), rslog.Nmass = log(mean(10^log.Nmass, na.rm=T),base=10))
traitsspp$Species <- gsub("\\."," ", traitsspp$GE.SP)
traitsspp$GrForm <- rep("T", times=nrow(traitsspp))
shrubs <- c("Rhododendron.macrophyllum", "Ribes.divaricatum","Lithocarpus.densiflorus","Frangula.purshiana","Holodiscus.discolor","Purshia.tridentate" ,
            "Cercocarpus.unknown","Chrysolepis.chrysophylla","Cornus.unknown","Corylus.cornuta","Calocedrus.decurrens","Ceanothus.velutinus","Arbutus.menziesii","Alnus.rubra" )
## need to check on Arbutus
traitsspp$GrForm[which(traitsspp$GE.SP %in% shrubs)] <- "S"
traitsspp$DecidEver <- "E"
# we'll just super quickly call everything in shrubs decid and ignore the rest
traitsspp$DecidEver[which(traitsspp$GE.SP %in% c(shrubs,"Acer.circinatum", "Acer.macrophyllum","Larix.occidentalis" ))] <- "D"
traitsspp$NeedleBroad <- "NA"
traitsspp$C3.C4 <- "C3"
colnames(traitsspp)[2] <- "Genus"
colnames(traitsspp) <- gsub("slog", "log", colnames(traitsspp))

traitsspp <- traitsspp[,match(colnames(LESspp), colnames(traitsspp))]
### combine LES species means and traits species means
allspp <- rbind(data.frame(LESspp), data.frame(traitsspp))
  #*** This has a handful of duplicate spp. But I haven't solved it yet.

#### Genus mean dataframe ####
LESgen <- allspp %>% group_by(Genus)  %>% summarise(Fam = sort(unique(Family))[1], Needle.Broad = unique(na.omit(NeedleBroad))[1], glog.LL = mean(log.LL, na.rm=T), glog.LMA = mean(log.LMA, na.rm=T), glog.Nmass = mean(log.Nmass, na.rm=T),glog.Narea = mean(log.Narea, na.rm=T)
                                                    ,rglog.LL = log(mean(10^rlog.LL, na.rm=T),base=10), rglog.LMA = log(mean(10^rlog.LMA, na.rm=T),base=10), rglog.Nmass = log(mean(10^rlog.Nmass, na.rm=T),base=10), nspp = n() )
  # taking mean of raw values or logged values doesn't matter all that much yet
colnames(LESgen) <- gsub("glog", "log", colnames(LESgen))
colnames(LESgen)[2] <- "Family"


#### Family Mean dataframe ####
LESfam <- LESgen %>% group_by(Family) %>% summarise(flog.LL = mean(log.LL, na.rm=T), flog.LMA = mean(log.LMA, na.rm=T), flog.Nmass = mean(log.Nmass, na.rm=T),flog.Narea = mean(log.Narea, na.rm=T)
                                                  ,rflog.LL = log(mean(10^rlog.LL, na.rm=T),base=10), rflog.LMA = log(mean(10^rlog.LMA, na.rm=T),base=10), rflog.Nmass = log(mean(10^rlog.Nmass, na.rm=T),base=10)
                                                 , tnspp = sum(nspp), ngen=n() )
colnames(LESfam) <- gsub("flog", "log", colnames(LESfam))
  # still doesn't matter too too much. just moves some of the points towards larger values in rlog things

ggplot(LESspp , aes(x=log.Nmass, y=log.LMA)) +
  geom_point(col="grey") +
  geom_point(data=LESgen, size=2) +
  geom_point(data=LESfam, col="darkred", aes(size=log(tnspp)) ,alpha=1/2)


### family on top of genus on top of spp
ggplot(LES, aes(x=log.Nmass, y=log.LMA)) + geom_point(col="grey") + geom_point(data=LESgen, aes(x=rlog.Nmass, y=rlog.LMA), size=3) + geom_point(data=LESfam[which(LESfam$tnspp>2),], aes(x=rlog.Nmass,y=rlog.LMA, size=ngen), col='darkred')

ggplot(data.all, aes(x=log.Nmass, y=log.LMA)) + geom_point(col="grey") + geom_smooth(col="black", method="lm",se=FALSE) +
  geom_point(data=data.all, aes(col=Genus)) + 
  geom_point(data=LESfam[which(LESfam$tnspp>2),], aes(x=rlog.Nmass,y=rlog.LMA, size=ngen), col='darkred')


#### Log.Nmass~log.LMA #################

win.spp.mod <- lmer(log.Nmass~log.LMA + (1|SP.ID) + (0 + log.LMA|SP.ID), traits.Nclean)

win.genus.mod <- lmer(log.Nmass~log.LMA  + (1|Genus) + (0 + log.LMA|Genus), LESspp)








LESgens <- LES[which(LES$Genus %in% names(which(xtabs(~Genus, LES)>4))), ]
### plot the LES + species.means from PACNW.
ggplot(LES, aes(x=log.LMA, y=log.Nmass, col=Needle.Broad.lf)) + geom_point() + geom_smooth(method = "lm") + 
  geom_point(data=species.means, aes(x=mlog.LMA, y=mlog.Nmass, col=NULL))


### looking at species w/>5 points in the GLOPNET dataset
p1 <- ggplot(LES, aes(x=log.LL, y=log.LMA)) + geom_point(color="grey") + geom_point(data=traits.common, aes(y=log.LMA_PSA, colour=SP.ID), alpha=1/5) +
  geom_smooth(method="lm", se=F, colour='darkgrey') + geom_smooth(data=traits.common, aes(y=log.LMA_PSA, colour=SP.ID),method="lm", se=F) +
  geom_point(data=LES[which(LES$Species %in% commonspp),], aes(colour=Species)) + geom_smooth(data=LES[which(LES$Species %in% commonspp),], aes(colour=Species), method="lm",se=F) +
  theme(legend.position="none")


### looking at across genus relationships
p2 <- ggplot(LES, aes(x=log.LL, y=log.LMA)) +
  geom_point(color="grey") + geom_smooth(method="lm", se=F, colour='darkgrey') +
  geom_smooth(data=LES[which(LES$Genus %in% commongna),], aes(colour=Genus), method="lm",se=F) +
  theme(legend.position="none")

geom_point(data=LES[which(LES$Genus %in% commongna),], aes(colour=Genus)) +

multiplot(p1,p2, cols=2)


######### LL vs LMA

ggplot(LES, aes(x=log.LMA, y=log.LL)) + geom_point(color="grey") + geom_point(data=traits.common, aes(x=log.LMA_PSA, colour=SP.ID), alpha=1/5) +
  geom_smooth(method="lm", se=F, colour='darkgrey') + geom_smooth(data=traits.common, aes(x=log.LMA_PSA, colour=SP.ID),method="lm", se=F) +
  geom_point(data=LES[which(LES$Species %in% commonspp),], aes(colour=Species)) + geom_smooth(data=LES[which(LES$Species %in% commonspp),], aes(colour=Species), method="lm",se=F)


### looking at across genus relationships
ggplot(LES, aes(x=log.LMA, y=log.LL)) + geom_point(color="grey") + geom_smooth(method="lm", se=F, colour='darkgrey') +
  geom_point(data=LES[which(LES$Genus %in% commongen),], aes(colour=Genus)) + geom_smooth(data=LES[which(LES$Genus %in% commongen),], aes(colour=Genus), method="lm",se=F)

win.spp.mod <- lmer(log.LL~log.LMA + (1|SP.ID) + (0 + log.LMA|SP.ID), traits.common)
win.genus.mod <- lmer(log.LL~log.LMA  + (1|Genus) + (0 + log.LMA|Genus), LES)
summary(win.spp.mod)
summary(win.genus.mod)
summary(lm(log.LL~log.LMA, LES))

xtabs(~Genus, LES)
commonspp <- names(which(xtabs(~Species, LES)>5))
commongen <- names(which(xtabs(~Genus, LES)>5))
commonfam <- names(which(xtabs(~Family, LES)>5))
                   






#_______________________________________________________________________

############ ***Piecewise hierarchical analysis*** ################
#_______________________________________________________________________

# while I'm trying to get the Bayesian approach up and running, just going to do a first pass with ML/piecewise correlations
  # going to loop over species, then genera, then familes, etc. and calculate all correlations and Major Axis regression slopes
  # then results are the distribution of all of those slopes/rhos at different taxo levels

# going to roll with traits.common5 for the moment, just so I have a larger sample size.

# plot of all repped species w>5records in either LES or traits.common5
ggplot(traits.common5, aes(x=log.Nmass, y=log.LMA, col=SP.ID)) + 
  geom_point() + geom_smooth(method="lm", se = F) + 
  geom_point(data=LES[which(LES$Species %in% commonspp),], aes(col=Species), size=2) + geom_smooth(data=LES[which(LES$Species %in% commonspp),], aes(col=Species), method="lm", se=F)


# Narea vs LMA
ggplot(spp.data, aes(x=log.LMA, y=log.Narea, col=Species)) + geom_point(data=LES, col="grey") + geom_point() + geom_smooth(method="lm", se=F) + geom_smooth(data=LES, aes(col=NULL), method="lm", se=F)
  #in arabidopsis
ggplot(arab, aes(x=log(LMA), y=log(Narea), col=Type)) + geom_point()

###### Data preperation: #########

##### combined common spp from PACNW and LES datasets
commonspp <- names(which(xtabs(~Species, LES)>=5)) # 6 have >5 records, 14 have >4 records, 36 have >3 records
## NOTE: I orginally was using log.LMA_PSA, but as of 3.14.17, I had some NAs in the traits dataset. I had to fill them in.
spp.data1 <- traits.common5 %>% select(FullSpecies,log.Nmass, log.LL, log.LMA_PSA, GENUS, Family, log.Narea)
colnames(spp.data1)[c(1,4,5)]<- c("Species","log.LMA", "Genus")
spp.data2 <- LES %>% filter(Species %in% commonspp) %>% select(Species, log.Nmass, log.LL, log.LMA, Genus, Family, log.Narea)
# NOTE: there are 95 trait measurements where LL is set at 1 year (log.LL = 1.079181). Luckily, I think my variance test in fit.MAR weeds these out...

## created combined dataset with 29 total species that have >5 observations
spp.data <- rbind(spp.data1, spp.data2)
spp.data$Species <- factor(spp.data$Species)
## look to make sure the combination worked 
#ggplot(spp.data, aes(x=log.LMA, y=log.Nmass, col=Species)) + geom_point() + geom_smooth(method="lm", se=F)


##### Genus level data
# currently just working with LES until I combine the PACNW dataset into this.
#gen.data <- LESspp[which(LESspp$Genus %in% names(which(xtabs(~Genus,LESspp)>5))),]
# 634 measurements of 51 genera
# update 03.14.17, added traits spp to LES
gen.data <- allspp[which(allspp$Genus %in% names(which(xtabs(~Genus, allspp)>=5))),]
# 659 measurements of 52 genera. Added Abies, but evidently that was the only additional genus...
gen.data$Genus <- factor(gen.data$Genus) # get rid of unused levels


##### Spp. in Family level data
# currently just working with LES until I combine the PACNW dataset into this.
# ppinfam.data <- LESspp[which(LESspp$Family %in% names(which(xtabs(~Family,LESspp)>5))),]
# 1374 measurements of 62 Families
sppinfam.data <- allspp[which(allspp$Family %in% names(which(xtabs(~Family,allspp)>=5))),]
sppinfam.data$Family <- factor(sppinfam.data$Family) # get rid of unused levels
# 1418 measurements of 63 Families. Added a family!

#### gen in Family level data
# currently just working with LES until I combine the PACNW dataset into this. # good to use LESgen, because it's been updated with allspp above
geninfam.data <- LESgen[which(LESgen$Family %in% names(which(xtabs(~Family,LESgen)>=5))),]
# 556 measurements of 40 Families , added 2 Families from last batch, and probably quite a few obs (old n_LMA.N = 459)
geninfam.data$Family <- factor(geninfam.data$Family) # get rid of unused levels
geninfam.data$Genus <- factor(geninfam.data$Genus)

#### Family Means:
fam.data <- LESfam # 189 families
fam.dataclean <- LESfam[which(LESfam$tnspp>2),] # 97 families


nsp <- spp.data %>% group_by (Species) %>% summarise(nLL = length(which(!is.na(log.LL))), nLMA = length(which(!is.na(log.LMA))), nN = length(which(!is.na(log.Nmass))))
# apply(nsp[,2:4],MARGIN = 2, FUN=mean)
ngen <- gen.data %>% group_by (Genus) %>% summarise(nLL = length(which(!is.na(log.LL))), nLMA = length(which(!is.na(log.LMA))), nN = length(which(!is.na(log.Nmass))))
# apply(ngen[,2:4],MARGIN = 2, FUN=mean)
nfam <- sppinfam.data %>% group_by (Family) %>% summarise(nLL = length(which(!is.na(log.LL))), nLMA = length(which(!is.na(log.LMA))), nN = length(which(!is.na(log.Nmass))))
# apply(nfam[,2:4],MARGIN = 2, FUN=mean)
nfamg <- geninfam.data %>% group_by (Family) %>% summarise(nLL = length(which(!is.na(log.LL))), nLMA = length(which(!is.na(log.LMA))), nN = length(which(!is.na(log.Nmass))))
# apply(nfamg[,2:4],MARGIN = 2, FUN=mean)


#### Function for fitting MARs

fit.MAR <- function(xvar, yvar, data, method="SMA") {
  if(method=="SMA") meth = 3
  if(method=="MA") meth = 2
  if(method=="OLS") meth =1
  if(length(t(as.vector(data[!(is.na(data[,xvar]) | is.na(data[,yvar])), xvar])))<3){
    return(rep(NA, times=7))
  }
  else{
    if(var(data[,yvar], na.rm=T)==0){
      return(rep(NA, times=7))
    }
    else{
      tmp.mod <- suppressMessages(lmodel2(get(yvar)~get(xvar), data))
      intercept <- tmp.mod$regression.results$Intercept[meth]
      slope <- tmp.mod$regression.results$Slope[meth]
      rho <- tmp.mod$r
      r.sq <- tmp.mod$rsquare
      n <- tmp.mod$n
      varX <- var(data[,xvar], na.rm=T)
      varY <- var(data[,yvar], na.rm=T)
      results <- c(intercept, slope,rho, r.sq, n, varX, varY)
      return(results)
    }
  }
}




############# **LMA vs Nmass** ###########
############ .Species level analysis #####################

spp.results <- data.frame(matrix(NA, nrow=length(unique(spp.data$Species)), ncol=8))
colnames(spp.results) <- c("Species", "Int","Slope","Rho","r.sq","n","varLMA","varNmass")
for(i in 1:length(unique(spp.data$Species))){
  species <- levels(spp.data$Species)[i]
  print(species)
  dataz <- spp.data[which(spp.data$Species==species),]
  res <- fit.MAR(xvar='log.LMA',yvar="log.Nmass",data=dataz)
  spp.results[i,1] <- species
  spp.results[i,2:8] <- res
}
  # Seems to be working!!!

############ .Genus level analysis #####################

gen.results <- data.frame(matrix(NA, nrow=length(unique(gen.data$Genus)), ncol=8))
colnames(gen.results) <- c("Genus", "Int","Slope","Rho","r.sq","n","varLMA","varNmass")
for(i in 1:length(unique(gen.data$Genus))){
  genus <- levels(gen.data$Genus)[i]
  print(genus)
  dataz <- gen.data[which(gen.data$Genus==genus),]
  gen.results[i,1] <- genus
  res <- fit.MAR(xvar='log.LMA',yvar="log.Nmass",data=dataz)
  gen.results[i,2:8] <- res 

}


############ .spp w/in Family level analysis #####################

sppinfam.results <- data.frame(matrix(NA, nrow=length(unique(sppinfam.data$Family)), ncol=8))
colnames(sppinfam.results) <- c("Family", "Int","Slope","Rho","r.sq","n","varLMA","varNmass")
for(i in 1:length(unique(sppinfam.data$Family))){
  family <- levels(sppinfam.data$Family)[i]
  print(family)
  dataz <- sppinfam.data[which(sppinfam.data$Family==family),]
  res <- fit.MAR(xvar='log.LMA',yvar="log.Nmass",data=dataz)
  sppinfam.results[i,1] <- family
  sppinfam.results[i,2:8] <- res
}


############ .gen w/in Family level analysis #####################

geninfam.results <- data.frame(matrix(NA, nrow=length(unique(geninfam.data$Family)), ncol=8))
colnames(geninfam.results) <- c("Family", "Int","Slope","Rho","r.sq","n","varLMA","varNmass")
for(i in 1:length(unique(geninfam.data$Family))){
  family <- levels(geninfam.data$Family)[i]
  print(family)
  dataz <- geninfam.data[which(geninfam.data$Family==family),]
  res <- fit.MAR(xvar='log.LMA',yvar="log.Nmass",data=dataz)
  geninfam.results[i,1] <- family
  geninfam.results[i,2:8] <- res
}


############ .Family level analysis #####################

# currently just working with LES until I combine the PACNW dataset into this.
famcorr.all <- cor(x=fam.data[,c("log.LL","log.LMA","log.Nmass")], use = "pairwise.complete.obs")
famcorr.clean <- cor(x=fam.dataclean[,c("log.LL","log.LMA","log.Nmass")], use = "pairwise.complete.obs")

fam.res_LMA.N <- c("fam.all", fit.MAR(xvar='log.LMA',yvar="log.Nmass",data=fam.data),"Fam")
names(fam.res_LMA.N) <- c("Taxo.Unit","Int","Slope","Rho","r.sq", "n","varLMA","varNmass","Type")
fam.resclean_LMA.N <- c("fam.clean", fit.MAR(xvar='log.LMA',yvar="log.Nmass",data=fam.dataclean), "Fam.clean")
names(fam.resclean_LMA.N) <- c("Taxo.Unit","Int","Slope","Rho","r.sq", "n","varLMA","varNmass","Type")


###### .Combining Results into one dataframe #####
# first add a "Type" column
spp.results$Type <- rep("w/inSpp", times=nrow(spp.results))
gen.results$Type <- rep("w/inGen", times=nrow(gen.results))
sppinfam.results$Type <- rep("Sppw/inFam", times=nrow(sppinfam.results))
geninfam.results$Type <- rep("Genw/inFam", times=nrow(geninfam.results))
# now make the column names all match
colnames(spp.results)[1] <- "Taxo.Unit"
colnames(gen.results)[1] <- "Taxo.Unit"
colnames(sppinfam.results)[1] <- "Taxo.Unit"
colnames(geninfam.results)[1] <- "Taxo.Unit"
all.results.LMAN <-rbind(spp.results,gen.results, sppinfam.results, geninfam.results, fam.res_LMA.N, fam.resclean_LMA.N)
all.results.LMAN$Type <- factor(all.results.LMAN$Type)
levels(all.results.LMAN$Type) <- list(w.inSpp = "w/inSpp",   w.inGen= "w/inGen",    Sppw.inFam="Sppw/inFam", Genw.inFam="Genw/inFam", Fam = "Fam", Famclean = "Fam.clean")







####### ***LL and LMA*** #####################

############ .Species level analysis #####################


spp.results <- data.frame(matrix(NA, nrow=length(unique(spp.data$Species)), ncol=8))
colnames(spp.results) <- c("Species", "Int","Slope","Rho","r.sq","n","varLMA","varLL")
for(i in 1:length(unique(spp.data$Species))){
  species <- levels(spp.data$Species)[i]
  print(species)
  dataz <- spp.data[which(spp.data$Species==species),]
  res <- fit.MAR(xvar='log.LMA',yvar="log.LL",data=dataz)
  spp.results[i,1] <- species
  spp.results[i,2:8] <- res
}


############ .Genus level analysis #####################

# currently just working with LES until I combine the PACNW dataset into this.
gen.data <- LESspp[which(LESspp$Genus %in% names(which(xtabs(~Genus,LESspp)>5))),]
# 634 measurements of 51 genera
gen.data$Genus <- factor(gen.data$Genus) # get rid of unused levels

gen.results <- data.frame(matrix(NA, nrow=length(unique(gen.data$Genus)), ncol=8))
colnames(gen.results) <- c("Genus", "Int","Slope","Rho","r.sq","n","varLMA","varLL")
for(i in 1:length(unique(gen.data$Genus))){
  genus <- levels(gen.data$Genus)[i]
  print(genus)
  dataz <- gen.data[which(gen.data$Genus==genus),]
  gen.results[i,1] <- genus
  res <- fit.MAR(xvar='log.LMA',yvar="log.LL",data=dataz)
  gen.results[i,2:8] <- res 
  
}



############ .spp w/in Family level analysis #####################

# currently just working with LES until I combine the PACNW dataset into this.
sppinfam.data <- LESspp[which(LESspp$Family %in% names(which(xtabs(~Family,LESspp)>5))),]
# 1374 measurements of 62 Families
sppinfam.data$Family <- factor(sppinfam.data$Family) # get rid of unused levels

sppinfam.results <- data.frame(matrix(NA, nrow=length(unique(sppinfam.data$Family)), ncol=8))
colnames(sppinfam.results) <- c("Family", "Int","Slope","Rho","r.sq","n","varLMA","varLL")
for(i in 1:length(unique(sppinfam.data$Family))){
  family <- levels(sppinfam.data$Family)[i]
  print(family)
  dataz <- sppinfam.data[which(sppinfam.data$Family==family),]
  res <- fit.MAR(xvar='log.LMA',yvar="log.LL",data=dataz)
  sppinfam.results[i,1] <- family
  sppinfam.results[i,2:8] <- res
}


############ .gen w/in Family level analysis #####################

# currently just working with LES until I combine the PACNW dataset into this.
geninfam.data <- LESgen[which(LESgen$Family %in% names(which(xtabs(~Family,LESgen)>5))),]
# 1374 measurements of 62 Families
geninfam.data$Family <- factor(geninfam.data$Family) # get rid of unused levels

geninfam.results <- data.frame(matrix(NA, nrow=length(unique(geninfam.data$Family)), ncol=8))
colnames(geninfam.results) <- c("Family", "Int","Slope","Rho","r.sq","n","varLMA","varLL")
for(i in 1:length(unique(geninfam.data$Family))){
  family <- levels(geninfam.data$Family)[i]
  print(family)
  dataz <- geninfam.data[which(geninfam.data$Family==family),]
  res <- fit.MAR(xvar='log.LMA',yvar="log.LL",data=dataz)
  geninfam.results[i,1] <- family
  geninfam.results[i,2:8] <- res
}


############ .Family level analysis #####################

# currently just working with LES until I combine the PACNW dataset into this.

famcorr.all <- cor(x=fam.data[,c("log.LL","log.LMA","log.LL")], use = "pairwise.complete.obs")
famcorr.clean <- cor(x=fam.dataclean[,c("log.LL","log.LMA","log.LL")], use = "pairwise.complete.obs")

fam.res_LMA.LL <- c("fam.all", fit.MAR(xvar='log.LMA',yvar="log.LL",data=fam.data),"Fam")
names(fam.res_LMA.LL) <- c("Taxo.Unit","Int","Slope","Rho","r.sq", "n","varLMA","varNmass","Type")
fam.resclean_LMA.LL <- c("fam.clean", fit.MAR(xvar='log.LMA',yvar="log.LL",data=fam.dataclean), "Fam.clean")
names(fam.resclean_LMA.LL) <- c("Taxo.Unit","Int","Slope","Rho","r.sq", "n","varLMA","varNmass","Type")



###### .Combining Results into one dataframe #####

# first add a "Type" column
spp.results$Type <- rep("w/inSpp", times=nrow(spp.results))
gen.results$Type <- rep("w/inGen", times=nrow(gen.results))
sppinfam.results$Type <- rep("Sppw/inFam", times=nrow(sppinfam.results))
geninfam.results$Type <- rep("Genw/inFam", times=nrow(geninfam.results))
# now make the column names all match
colnames(spp.results)[1] <- "Taxo.Unit"
colnames(gen.results)[1] <- "Taxo.Unit"
colnames(sppinfam.results)[1] <- "Taxo.Unit"
colnames(geninfam.results)[1] <- "Taxo.Unit"
all.results.LMALL <-rbind(spp.results,gen.results, sppinfam.results, geninfam.results, fam.res_LMA.LL, fam.resclean_LMA.LL)
all.results.LMALL$Type <- factor(all.results.LMALL$Type)
levels(all.results.LMALL$Type) <- list(w.inSpp = "w/inSpp",   w.inGen= "w/inGen",    Sppw.inFam="Sppw/inFam", Genw.inFam="Genw/inFam", Fam = "Fam", Famclean = "Fam.clean")





####### ***LL and Nmass*** #####################



spp.results <- data.frame(matrix(NA, nrow=length(unique(spp.data$Species)), ncol=8))
colnames(spp.results) <- c("Species", "Int","Slope","Rho","r.sq","n","varNmass","varLL")
for(i in 1:length(unique(spp.data$Species))){
  species <- levels(spp.data$Species)[i]
  print(species)
  dataz <- spp.data[which(spp.data$Species==species),]
  res <- fit.MAR(xvar='log.Nmass',yvar="log.LL",data=dataz)
  spp.results[i,1] <- species
  spp.results[i,2:8] <- res
}


############ .Genus level analysis #####################

# currently just working with LES until I combine the PACNW dataset into this.
gen.data <- LESspp[which(LESspp$Genus %in% names(which(xtabs(~Genus,LESspp)>5))),]
# 634 measurements of 51 genera
gen.data$Genus <- factor(gen.data$Genus) # get rid of unused levels

gen.results <- data.frame(matrix(NA, nrow=length(unique(gen.data$Genus)), ncol=8))
colnames(gen.results) <- c("Genus", "Int","Slope","Rho","r.sq","n","varNmass","varLL")
for(i in 1:length(unique(gen.data$Genus))){
  genus <- levels(gen.data$Genus)[i]
  print(genus)
  dataz <- gen.data[which(gen.data$Genus==genus),]
  gen.results[i,1] <- genus
  res <- fit.MAR(xvar='log.Nmass',yvar="log.LL",data=dataz)
  gen.results[i,2:8] <- res 
  
}



############ .spp w/in Family level analysis #####################

# currently just working with LES until I combine the PACNW dataset into this.
sppinfam.data <- LESspp[which(LESspp$Family %in% names(which(xtabs(~Family,LESspp)>5))),]
# 1374 measurements of 62 Families
sppinfam.data$Family <- factor(sppinfam.data$Family) # get rid of unused levels

sppinfam.results <- data.frame(matrix(NA, nrow=length(unique(sppinfam.data$Family)), ncol=8))
colnames(sppinfam.results) <- c("Family", "Int","Slope","Rho","r.sq","n","varNmass","varLL")
for(i in 1:length(unique(sppinfam.data$Family))){
  family <- levels(sppinfam.data$Family)[i]
  print(family)
  dataz <- sppinfam.data[which(sppinfam.data$Family==family),]
  res <- fit.MAR(xvar='log.Nmass',yvar="log.LL",data=dataz)
  sppinfam.results[i,1] <- family
  sppinfam.results[i,2:8] <- res
}





############ .gen w/in Family level analysis #####################

# currently just working with LES until I combine the PACNW dataset into this.
geninfam.data <- LESgen[which(LESgen$Family %in% names(which(xtabs(~Family,LESgen)>5))),]
# 1374 measurements of 62 Families
geninfam.data$Family <- factor(geninfam.data$Family) # get rid of unused levels

geninfam.results <- data.frame(matrix(NA, nrow=length(unique(geninfam.data$Family)), ncol=8))
colnames(geninfam.results) <- c("Family", "Int","Slope","Rho","r.sq","n","varNmass","varLL")
for(i in 1:length(unique(geninfam.data$Family))){
  family <- levels(geninfam.data$Family)[i]
  print(family)
  dataz <- geninfam.data[which(geninfam.data$Family==family),]
  res <- fit.MAR(xvar='log.Nmass',yvar="log.LL",data=dataz)
  geninfam.results[i,1] <- family
  geninfam.results[i,2:8] <- res
}





############ .Family level analysis #####################

# currently just working with LES until I combine the PACNW dataset into this.

fam.data <- LESfam # 189 families
fam.dataclean <- LESfam[which(LESfam$tnspp>2),] # 97 families

famcorr.all <- cor(x=fam.data[,c("log.LMA","log.Nmass","log.LL")], use = "pairwise.complete.obs")
famcorr.clean <- cor(x=fam.dataclean[,c("log.LMA","log.Nmass","log.LL")], use = "pairwise.complete.obs")

fam.res_N.LL <- c("fam.all", fit.MAR(xvar='log.Nmass',yvar="log.LL",data=fam.data),"Fam")
names(fam.res_N.LL) <- c("Taxo.Unit","Int","Slope","Rho","r.sq", "n","varLMA","varNmass","Type")
fam.resclean_N.LL <- c("fam.clean", fit.MAR(xvar='log.Nmass',yvar="log.LL",data=fam.dataclean), "Fam.clean")
names(fam.res_N.LL) <- c("Taxo.Unit","Int","Slope","Rho","r.sq", "n","varLMA","varNmass","Type")



###### .Combining Results into one dataframe #####

# first add a "Type" column
spp.results$Type <- rep("w/inSpp", times=nrow(spp.results))
gen.results$Type <- rep("w/inGen", times=nrow(gen.results))
sppinfam.results$Type <- rep("Sppw/inFam", times=nrow(sppinfam.results))
geninfam.results$Type <- rep("Genw/inFam", times=nrow(geninfam.results))
# now make the column names all match
colnames(spp.results)[1] <- "Taxo.Unit"
colnames(gen.results)[1] <- "Taxo.Unit"
colnames(sppinfam.results)[1] <- "Taxo.Unit"
colnames(geninfam.results)[1] <- "Taxo.Unit"
all.results.NmassLL <-rbind(spp.results,gen.results, sppinfam.results, geninfam.results, fam.res_N.LL, fam.resclean_N.LL)
all.results.NmassLL$Type <- factor(all.results.NmassLL$Type)
levels(all.results.NmassLL$Type) <- list(w.inSpp = "w/inSpp",   w.inGen= "w/inGen",    Sppw.inFam="Sppw/inFam", Genw.inFam="Genw/inFam", Fam = "Fam", Famclean = "Fam.clean")






####### ***LMA and Narea!!!*** #####################



spp.results <- data.frame(matrix(NA, nrow=length(unique(spp.data$Species)), ncol=8))
colnames(spp.results) <- c("Species", "Int","Slope","Rho","r.sq","n","varLMA","varNarea")
for(i in 1:length(unique(spp.data$Species))){
  species <- levels(spp.data$Species)[i]
  print(species)
  dataz <- spp.data[which(spp.data$Species==species),]
  res <- fit.MAR(xvar='log.LMA',yvar="log.Narea",data=dataz)
  spp.results[i,1] <- species
  spp.results[i,2:8] <- res
}


############ .Genus level analysis #####################

# currently just working with LES until I combine the PACNW dataset into this.
gen.data <- LESspp[which(LESspp$Genus %in% names(which(xtabs(~Genus,LESspp)>5))),]
# 634 measurements of 51 genera
gen.data$Genus <- factor(gen.data$Genus) # get rid of unused levels

gen.results <- data.frame(matrix(NA, nrow=length(unique(gen.data$Genus)), ncol=8))
colnames(gen.results) <- c("Genus", "Int","Slope","Rho","r.sq","n","varLMA","varNarea")
for(i in 1:length(unique(gen.data$Genus))){
  genus <- levels(gen.data$Genus)[i]
  print(genus)
  dataz <- gen.data[which(gen.data$Genus==genus),]
  gen.results[i,1] <- genus
  res <- fit.MAR(xvar='log.LMA',yvar="log.Narea",data=dataz)
  gen.results[i,2:8] <- res 
  
}



############ .spp w/in Family level analysis #####################

# currently just working with LES until I combine the PACNW dataset into this.
sppinfam.data <- LESspp[which(LESspp$Family %in% names(which(xtabs(~Family,LESspp)>5))),]
# 1374 measurements of 62 Families
sppinfam.data$Family <- factor(sppinfam.data$Family) # get rid of unused levels

sppinfam.results <- data.frame(matrix(NA, nrow=length(unique(sppinfam.data$Family)), ncol=8))
colnames(sppinfam.results) <- c("Family", "Int","Slope","Rho","r.sq","n","varLMA","varNarea")
for(i in 1:length(unique(sppinfam.data$Family))){
  family <- levels(sppinfam.data$Family)[i]
  print(family)
  dataz <- sppinfam.data[which(sppinfam.data$Family==family),]
  res <- fit.MAR(xvar='log.LMA',yvar="log.Narea",data=dataz)
  sppinfam.results[i,1] <- family
  sppinfam.results[i,2:8] <- res
}





############ .gen w/in Family level analysis #####################

# currently just working with LES until I combine the PACNW dataset into this.
geninfam.data <- LESgen[which(LESgen$Family %in% names(which(xtabs(~Family,LESgen)>5))),]
# 1374 measurements of 62 Families
geninfam.data$Family <- factor(geninfam.data$Family) # get rid of unused levels

geninfam.results <- data.frame(matrix(NA, nrow=length(unique(geninfam.data$Family)), ncol=8))
colnames(geninfam.results) <- c("Family", "Int","Slope","Rho","r.sq","n","varLMA","varNarea")
for(i in 1:length(unique(geninfam.data$Family))){
  family <- levels(geninfam.data$Family)[i]
  print(family)
  dataz <- geninfam.data[which(geninfam.data$Family==family),]
  res <- fit.MAR(xvar='log.LMA',yvar="log.Narea",data=dataz)
  geninfam.results[i,1] <- family
  geninfam.results[i,2:8] <- res
}





############ .Family level analysis #####################

# currently just working with LES until I combine the PACNW dataset into this.

# fam.data <- LESfam # 189 families
# fam.dataclean <- LESfam[which(LESfam$tnspp>2),] # 97 families

# famcorr.all <- cor(x=fam.data[,c("log.LMA","log.Nmass","log.Narea")], use = "pairwise.complete.obs")
# famcorr.clean <- cor(x=fam.dataclean[,c("log.LMA","log.Nmass","log.Narea")], use = "pairwise.complete.obs")

fam.res_LMA.Narea <- c("fam.all", fit.MAR(xvar='log.LMA',yvar="log.Narea",data=fam.data),"Fam")
names(fam.res_LMA.Narea) <- c("Taxo.Unit","Int","Slope","Rho","r.sq", "n","varLMA","varLMA","Type")
fam.resclean_LMA.Narea <- c("fam.clean", fit.MAR(xvar='log.LMA',yvar="log.Narea",data=fam.dataclean), "Fam.clean")
names(fam.res_LMA.Narea) <- c("Taxo.Unit","Int","Slope","Rho","r.sq", "n","varLMA","varLMA","Type")



###### .Combining Results into one dataframe #####

# first add a "Type" column
spp.results$Type <- rep("w/inSpp", times=nrow(spp.results))
gen.results$Type <- rep("w/inGen", times=nrow(gen.results))
sppinfam.results$Type <- rep("Sppw/inFam", times=nrow(sppinfam.results))
geninfam.results$Type <- rep("Genw/inFam", times=nrow(geninfam.results))
# now make the column names all match
colnames(spp.results)[1] <- "Taxo.Unit"
colnames(gen.results)[1] <- "Taxo.Unit"
colnames(sppinfam.results)[1] <- "Taxo.Unit"
colnames(geninfam.results)[1] <- "Taxo.Unit"
all.results.LMANarea <-rbind(spp.results,gen.results, sppinfam.results, geninfam.results, fam.res_LMA.Narea, fam.resclean_LMA.Narea)
all.results.LMANarea$Type <- factor(all.results.LMANarea$Type)
levels(all.results.LMANarea$Type) <- list(w.inSpp = "w/inSpp",   w.inGen= "w/inGen",    Sppw.inFam="Sppw/inFam", Genw.inFam="Genw/inFam", Fam = "Fam", Famclean = "Fam.clean")





########### ** Combining all ANALYSES into 1 df ** ######

LMALL <- all.results.LMALL
colnames(LMALL)[2:6] <- paste(colnames(all.results.LMALL)[2:6],"LMA.LL", sep="_")
LMAN <- all.results.LMAN
colnames(LMAN)[2:6] <- paste(colnames(all.results.LMAN)[2:6],"LMA.N", sep="_")
NmassLL <- all.results.NmassLL
colnames(NmassLL)[2:6] <- paste(colnames(all.results.NmassLL)[2:6],"N.LL", sep="_")
LMANarea <- all.results.LMANarea
colnames(LMANarea)[2:6] <- paste(colnames(all.results.LMANarea)[2:6],"LMA.Narea", sep="_")

all.results <- cbind(LMALL, LMAN[,-c(1,7,9)], NmassLL[,-c(1,7,8,9)], LMANarea[,-c(1,7,9)]) # drop the duplicate 'var' columns

#write.csv(all.results, "Results_SimpleMAreg_v1_030817.csv")
#write.csv(all.results, "Results_SimpleMAreg_v2_031417.csv")
#write.csv(all.results, "Results_SimpleMAreg_v3_031717.csv")

#all.resultsold <- read.csv("Results_SimpleMAreg_v2_031417.csv")





############# Plotting Slopes of diff taxo levels #################

all.results <- read.csv("Results_SimpleMAreg_v3_031717.csv", row.names = 1)
levels(all.results$Type) <- list(w.inSpp = "w.inSpp", w.inGen = "w.inGen", Sppw.inFam= "Sppw.inFam",Genw.inFam="Genw.inFam", Fam="Fam",Famclean="Famclean")

######## Funnel Plots #######
quartz(width=8,height=5)
par(mar=c(3,3,1,1), mgp=c(2,1,0),mfrow=c(2,3), oma=c(0,0,2,0))
### Rho
plot(Rho~n, all.results.LMALL, pch=16, col=Type, main="LMA ~ LL")
abline(h=0, col="grey", lty=2)
abline(h=all.results.LMALL$Rho[which(all.results.LMALL$Taxo.Unit=="fam.all")])
plot(Rho~n, all.results.NmassLL, pch=16, col=Type, main="Nmass ~ LL")
abline(h=0, col="grey", lty=2)
abline(h=all.results.NmassLL$Rho[which(all.results.NmassLL$Taxo.Unit=="fam.all")])
plot(Rho~n, all.results.LMAN, pch=16, col=Type, main="LMA ~ Nmass")
abline(h=0, col="grey", lty=2)
abline(h=all.results.LMAN$Rho[which(all.results.LMAN$Taxo.Unit=="fam.all")])
legend('topright', legend = levels(all.results.LMAN$Type), pch=16, col=mypal, bty ="n")


### MA slopes
plot(Slope~n, all.results.LMALL, pch=16, col=Type, main="LMA ~ LL")
abline(h=0, col="grey", lty=2)
abline(h=all.results.LMALL$Slope[which(all.results.LMALL$Taxo.Unit=="fam.all")])
plot(Slope~n, all.results.NmassLL, pch=16, col=Type, main="Nmass ~ LL")
abline(h=0, col="grey", lty=2)
abline(h=all.results.NmassLL$Slope[which(all.results.NmassLL$Taxo.Unit=="fam.all")])
plot(Slope~n, all.results.LMAN, pch=16, col=Type, main="LMA ~ Nmass")
abline(h=0, col="grey", lty=2)
abline(h=all.results.LMAN$Slope[which(all.results.LMAN$Taxo.Unit=="fam.all")])
legend('topright', legend = levels(all.results.LMAN$Type), pch=16, col=mypal, bty ="n")


### Slopes as f(varLMA/LL)
plot(Slope~varLMA, all.results.LMALL, pch=16, col=Type)
abline(h=0, col="grey", lty=2)
abline(h=all.results.LMALL$Slope[which(all.results.LMALL$Taxo.Unit=="fam.all")])
plot(Slope~varNmass, all.results.NmassLL, pch=16, col=Type)
abline(h=0, col="grey", lty=2)
abline(h=all.results.NmassLL$Slope[which(all.results.NmassLL$Taxo.Unit=="fam.all")])
plot(Slope~varLMA, all.results.LMAN, pch=16, col=Type)
abline(h=0, col="grey", lty=2)
abline(h=all.results.LMAN$Slope[which(all.results.LMAN$Taxo.Unit=="fam.all")])

### Rho 
plot(Rho~varLMA, all.results.LMALL, pch=16, col=Type)
abline(h=0, col="grey", lty=2)
abline(h=all.results.LMALL$Rho[which(all.results.LMALL$Taxo.Unit=="fam.all")])
plot(Rho~varNmass, all.results.NmassLL, pch=16, col=Type)
abline(h=0, col="grey", lty=2)
abline(h=all.results.NmassLL$Rho[which(all.results.NmassLL$Taxo.Unit=="fam.all")])
plot(Rho~varLMA, all.results.LMAN, pch=16, col=Type)
abline(h=0, col="grey", lty=2)
abline(h=all.results.LMAN$Rho[which(all.results.LMAN$Taxo.Unit=="fam.all")])


### two panel of correlation of Nmass things with varNmass
quartz(width=3,height=4.5)
par(mar=c(4,4,0,1), mfrow=c(2,1), mgp=c(2.5,1,0), oma=c(0,0,2,0))
plot(Rho~varNmass, all.results.LMAN, pch=16, col=Type, xlab="")
#mtext( text="Nmass v LMA", side=3, line=0, font=2)
abline(h=0, col="grey", lty=2)
abline(h=all.results.LMAN$Rho[which(all.results.LMAN$Taxo.Unit=="fam.all")])
legend('topright', legend = levels(all.results.LMAN$Type), pch=16, col=mypal, bty ="n", cex=.7)
plot(Rho~varNmass, all.results.NmassLL, pch=16, col=Type, xlab="Var. in %N")
abline(h=0, col="grey", lty=2)
abline(h=all.results.NmassLL$Rho[which(all.results.NmassLL$Taxo.Unit=="fam.all")])


### 4 panel of Slope/Rho of LMA-LL things with varLL and varLMA
quartz(width=4.5,height=4.5)
par(mar=c(4,4,0,1), mfrow=c(2,2), mgp=c(2.5,1,0), oma=c(0,0,2,0))
plot(Rho_LMA.LL~varLMA, all.results, pch=16, col=Type, xlab="")
#mtext( text="Nmass v LMA", side=3, line=0, font=2)
abline(h=0, col="grey", lty=2)
abline(h=all.results$Rho_LMA.LL[which(all.results$Taxo.Unit=="fam.all")])
legend('topright', legend = levels(all.results.LMAN$Type), pch=16, col=mypal, bty ="n", cex=.7)
plot(Rho~varNmass, all.results.NmassLL, pch=16, col=Type, xlab="Var. in %N")
abline(h=0, col="grey", lty=2)
abline(h=all.results.NmassLL$Rho[which(all.results.NmassLL$Taxo.Unit=="fam.all")])



#_______________________
##### Boxplots 
#_______________________

quartz(width=3, height=4.5)
par(mfrow=c(2,1), mar=c(4,4,2,1), oma=c(2,0,0,0))
boxplot(Slope_LMA.N~Type, all.results, ylim=c(-2,1.5),las=3, main="log(LMA)~log(Nmass", ylab="MA Slope")
points(y=rep(as.numeric(fam.res_LMA.N[3]), times=7),x=seq(4.7,5.3, by=.1), pch=15, cex=.7)
points(y=rep(as.numeric(fam.resclean_LMA.N[3]), times=7),x=seq(5.7,6.3, by=.1), pch=15, cex=.7)
abline(h=0, lty=2)
boxplot(Rho_LMA.N~Type, all.results, ylim=c(-1,1),las=3, ylab="Rho")
points(y=rep(as.numeric(fam.res_LMA.N[4]), times=7),x=seq(4.7,5.3, by=.1), pch=15, cex=.7)
points(y=rep(as.numeric(fam.resclean_LMA.N[4]), times=7),x=seq(5.7,6.3, by=.1), pch=15, cex=.7)
abline(h=0, lty=2)


####LMA v Nmass: Trimming to only well represented taxa (10 or more)
quartz(width=3, height=4.5,title = "LMA v Nmass")
par(mfrow=c(2,1), mar=c(0,4,0,1), oma=c(6,0,2,0))
p <- boxplot(Slope_LMA.N~Type, all.results[which(all.results$n_LMA.N>5),]
        , ylim=c(-2.2,1.5),las=3, ylab="MA Slope" #, main="log(LMA)~log(Nmass) strict"
        , col=paste0(mypal[1:6],"66"), boxcol=paste0(mypal[1:6],"66")
        ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[1:6],"AA")
        , staplelwd=0, outpch=16, outcex=.5, outcol=mypal[1:6]
        , boxwex=.7, xaxt="n")
abline(h=0, lty=2)
#text(y=par()$usr[3]+.2, x=c(1,2,3,4,5,6), labels = p$n)
#text(y=par()$usr[4]-.2,x=.5, labels = "a)")
mtext(text = "b)", side = 3, adj=0, line=.2)
mtext(text= "LMA vs Nmass", side=3, line=.2)
p <- boxplot(Rho_LMA.N~Type, all.results[which(all.results$n_LMA.N>5),]
             , ylim=c(-1,1.1),las=3, ylab="Rho"
             , col=paste0(mypal[1:6],"66"), boxcol=paste0(mypal[1:6],"66")
             ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[1:6],"AA")
             , staplelwd=0, outpch=16, outcex=.5, outcol=mypal[1:6]
             , boxwex=.7)
abline(h=0, lty=2)
text(y=par()$usr[4]-.2, x=c(1,2,3,4,5,6), labels = p$n)
#text(y=par()$usr[4]-.2,x=.5, labels = "b)")




## only things with reasonable correlations 
quartz(width=3, height=4.5)
par(mfrow=c(2,1), mar=c(4,4,2,1), oma=c(2,0,0,0))
boxplot(Slope_LMA.N~Type, all.results[which(abs(all.results$Rho_LMA.N)>.3),], ylim=c(-2,1.5),las=3, main="log(LMA)~log(Nmass) Hcor", ylab="MA Slope")
points(y=rep(as.numeric(fam.res_LMA.N[3]), times=7),x=seq(4.7,5.3, by=.1), pch=15, cex=.7)
points(y=rep(as.numeric(fam.resclean_LMA.N[3]), times=7),x=seq(5.7,6.3, by=.1), pch=15, cex=.7)
abline(h=0, lty=2)
boxplot(Rho_LMA.N~Type, all.results[which(all.results$n_LMA.N>9),], ylim=c(-1,1),las=3, ylab="Rho")
points(y=rep(as.numeric(fam.res_LMA.N[4]), times=7),x=seq(4.7,5.3, by=.1), pch=15, cex=.7)
points(y=rep(as.numeric(fam.resclean_LMA.N[4]), times=7),x=seq(5.7,6.3, by=.1), pch=15, cex=.7)
abline(h=0, lty=2)



####### LMA vs LL ########

quartz(width=3, height=4.5)
par(mfrow=c(2,1), mar=c(4,4,2,1), oma=c(2,0,0,0))
boxplot(Slope_LMA.LL~Type, all.results, ylim=c(-2,2.5),las=3, main="log(LMA)~log(LL)", ylab="MA Slope")
points(y=rep(as.numeric(fam.res_LMA.LL[3]), times=7),x=seq(4.7,5.3, by=.1), pch=15, cex=.7)
points(y=rep(as.numeric(fam.resclean_LMA.LL[3]), times=7),x=seq(5.7,6.3, by=.1), pch=15, cex=.7)
abline(h=0, lty=2)
boxplot(Rho_LMA.LL~Type, all.results, ylim=c(-1,1),las=3, ylab="Rho")
points(y=rep(as.numeric(fam.res_LMA.LL[4]), times=7),x=seq(4.7,5.3, by=.1), pch=15, cex=.7)
points(y=rep(as.numeric(fam.resclean_LMA.LL[4]), times=7),x=seq(5.7,6.3, by=.1), pch=15, cex=.7)
abline(h=0, lty=2)



#### Trimming to only well represented taxa (10 or more)
quartz(width=3, height=4.5)
par(mfrow=c(2,1), mar=c(4,4,2,1), oma=c(2,0,0,0))
boxplot(Slope_LMA.LL~Type, all.results[which(all.results$n_LMA.LL>9),], ylim=c(-2,2.5),las=3, main="log(LMA)~log(LL) strict", ylab="MA Slope")
points(y=rep(as.numeric(fam.res_LMA.LL[3]), times=7),x=seq(4.7,5.3, by=.1), pch=15, cex=.7)
points(y=rep(as.numeric(fam.resclean_LMA.LL[3]), times=7),x=seq(5.7,6.3, by=.1), pch=15, cex=.7)
abline(h=0, lty=2)
boxplot(Rho_LMA.LL~Type, all.results[which(all.results$n_LMA.LL>9),], ylim=c(-1,1),las=3, ylab="Rho")
points(y=rep(as.numeric(fam.res_LMA.LL[4]), times=7),x=seq(4.7,5.3, by=.1), pch=15, cex=.7)
points(y=rep(as.numeric(fam.resclean_LMA.LL[4]), times=7),x=seq(5.7,6.3, by=.1), pch=15, cex=.7)
abline(h=0, lty=2)


## only things with reasonable correlations 
quartz(width=3, height=4.5)
par(mfrow=c(2,1), mar=c(4,4,2,1), oma=c(2,0,0,0))
boxplot(Slope_LMA.LL~Type, all.results[which(abs(all.results$Rho_LMA.LL)>.3),], ylim=c(-2,2.5),las=3, main="log(LMA)~log(LL) Hcor", ylab="MA Slope")
points(y=rep(as.numeric(fam.res_LMA.LL[3]), times=7),x=seq(4.7,5.3, by=.1), pch=15, cex=.7)
points(y=rep(as.numeric(fam.resclean_LMA.LL[3]), times=7),x=seq(5.7,6.3, by=.1), pch=15, cex=.7)
abline(h=0, lty=2)
boxplot(Rho_LMA.LL~Type, all.results[which(all.results$n_LMA.LL>9),], ylim=c(-1,1),las=3, ylab="Rho")
points(y=rep(as.numeric(fam.res_LMA.LL[4]), times=7),x=seq(4.7,5.3, by=.1), pch=15, cex=.7)
points(y=rep(as.numeric(fam.resclean_LMA.LL[4]), times=7),x=seq(5.7,6.3, by=.1), pch=15, cex=.7)
abline(h=0, lty=2)




########### Nmass vs LL #########
## going to throw all the results together so I can boxplot/violin plot them all

quartz(width=3, height=4.5)
par(mfrow=c(2,1), mar=c(4,4,2,1), oma=c(2,0,0,0))
boxplot(Slope_N.LL~Type, all.results, ylim=c(-3,2.5),las=3, main="log(Nmass)~log(LL)", ylab="MA Slope")
points(y=rep(as.numeric(fam.res_N.LL[3]), times=7),x=seq(4.7,5.3, by=.1), pch=15, cex=.7)
points(y=rep(as.numeric(fam.resclean_N.LL[3]), times=7),x=seq(5.7,6.3, by=.1), pch=15, cex=.7)
abline(h=0, lty=2)
boxplot(Rho_N.LL~Type, all.results, ylim=c(-1,1),las=3, ylab="Rho")
points(y=rep(as.numeric(fam.res_N.LL[4]), times=7),x=seq(4.7,5.3, by=.1), pch=15, cex=.7)
points(y=rep(as.numeric(fam.resclean_N.LL[4]), times=7),x=seq(5.7,6.3, by=.1), pch=15, cex=.7)
abline(h=0, lty=2)

# quartz(width=3, height=4.5)
# par(mfrow=c(2,1), mar=c(4,4,2,1), oma=c(2,0,0,0))
# boxplot(Slope_N.LL~Type, all.results[which(all.results$n_N.LL>9),], ylim=c(-3,2.5),las=3, main="log(LL)~log(Nmass) strict", ylab="MA Slope")
# points(y=rep(as.numeric(fam.res_N.LL[3]), times=7),x=seq(4.7,5.3, by=.1), pch=15, cex=.7)
# points(y=rep(as.numeric(fam.resclean_N.LL[3]), times=7),x=seq(5.7,6.3, by=.1), pch=15, cex=.7)
# abline(h=0, lty=2)
# boxplot(Rho_N.LL~Type, all.results[which(all.results$n_N.LL>9),], ylim=c(-1,1),las=3, ylab="Rho")
# points(y=rep(as.numeric(fam.res_N.LL[4]), times=7),x=seq(4.7,5.3, by=.1), pch=15, cex=.7)
# points(y=rep(as.numeric(fam.resclean_N.LL[4]), times=7),x=seq(5.7,6.3, by=.1), pch=15, cex=.7)
# abline(h=0, lty=2)


####LL v Nmass: Trimming to only well represented taxa (10 or more)
quartz(width=3, height=4.5, title="LL v Nmass")
par(mfrow=c(2,1), mar=c(0,4,0,1), oma=c(6,0,2,0))
p <- boxplot(Slope_N.LL~Type, all.results[which(all.results$n_N.LL>5),]
             , ylim=c(-3.5,3),las=3, ylab="MA Slope" #, main="log(LMA)~log(Nmass) strict"
             , col=paste0(mypal[1:6],"66"), boxcol=paste0(mypal[1:6],"66")
             ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[1:6],"AA")
             , staplelwd=0, outpch=16, outcex=.5, outcol=mypal[1:6]
             , boxwex=.7, xaxt="n")
abline(h=0, lty=2)
#text(y=par()$usr[3]+.2, x=c(1,2,3,4,5,6), labels = p$n)
#text(y=par()$usr[4]-.2,x=.5, labels = "b)")
mtext(text = "b)", side = 3, adj=0, line=.2)
mtext(text= "LL vs Nmass", side=3, line=.2)

p <- boxplot(Rho_N.LL~Type, all.results[which(all.results$n_N.LL>5),]
             , ylim=c(-1,1.1),las=3, ylab="Rho"
             , col=paste0(mypal[1:6],"66"), boxcol=paste0(mypal[1:6],"66")
             ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[1:6],"AA")
             , staplelwd=0, outpch=16, outcex=.5, outcol=mypal[1:6]
             , boxwex=.7)
abline(h=0, lty=2)
text(y=par()$usr[4]-.2, x=c(1,2,3,4,5,6), labels = p$n)
#text(y=par()$usr[4]-.2,x=.5, labels = "b)")




## only things with reasonable correlations 
quartz(width=3, height=4.5)
par(mfrow=c(2,1), mar=c(4,4,2,1), oma=c(2,0,0,0))
boxplot(Slope_N.LL~Type, all.results[which(abs(all.results$Rho_N.LL)>.3),], ylim=c(-3,2.5),las=3, main="log(LL)~log(Nmass) Hcor", ylab="MA Slope")
points(y=rep(as.numeric(fam.res_N.LL[3]), times=7),x=seq(4.7,5.3, by=.1), pch=15, cex=.7)
points(y=rep(as.numeric(fam.resclean_N.LL[3]), times=7),x=seq(5.7,6.3, by=.1), pch=15, cex=.7)
abline(h=0, lty=2)
boxplot(Rho_N.LL~Type, all.results[which(all.results$n_N.LL>9),], ylim=c(-1,1),las=3, ylab="Rho")
points(y=rep(as.numeric(fam.res_N.LL[4]), times=7),x=seq(4.7,5.3, by=.1), pch=15, cex=.7)
points(y=rep(as.numeric(fam.resclean_N.LL[4]), times=7),x=seq(5.7,6.3, by=.1), pch=15, cex=.7)
abline(h=0, lty=2)



###### LMA v Narea !! ########
quartz(width=3, height=4.5)
par(mfrow=c(2,1), mar=c(4,4,2,1), oma=c(2,0,0,0))
boxplot(Slope_LMA.Narea~Type, all.results,las=3,ylim=c(-1.5,2), main="log(LMA)~log(Narea)", ylab="MA Slope")
points(y=rep(as.numeric(fam.res_LMA.Narea[3]), times=7),x=seq(4.7,5.3, by=.1), pch=15, cex=.7)
points(y=rep(as.numeric(fam.resclean_LMA.Narea[3]), times=7),x=seq(5.7,6.3, by=.1), pch=15, cex=.7)
abline(h=0, lty=2)
boxplot(Rho_LMA.Narea~Type, all.results, ylim=c(-1,1),las=3, ylab="Rho")
points(y=rep(as.numeric(fam.res_LMA.Narea[4]), times=7),x=seq(4.7,5.3, by=.1), pch=15, cex=.7)
points(y=rep(as.numeric(fam.resclean_LMA.Narea[4]), times=7),x=seq(5.7,6.3, by=.1), pch=15, cex=.7)
abline(h=0, lty=2)


#### Trimming to only well represented taxa (10 or more)
quartz(width=3, height=4.5)
par(mfrow=c(2,1), mar=c(4,4,2,1), oma=c(2,0,0,0))
boxplot(Slope_LMA.Narea~Type, all.results[which(all.results$n_LMA.Narea>9),], ylim=c(-1.5,2),las=3, main="log(LMA)~log(Narea) strict", ylab="MA Slope")
points(y=rep(as.numeric(fam.res_LMA.Narea[3]), times=7),x=seq(4.7,5.3, by=.1), pch=15, cex=.7)
points(y=rep(as.numeric(fam.resclean_LMA.Narea[3]), times=7),x=seq(5.7,6.3, by=.1), pch=15, cex=.7)
abline(h=0, lty=2)
boxplot(Rho_LMA.Narea~Type, all.results[which(all.results$n_LMA.Narea>9),], ylim=c(-1,1),las=3, ylab="Rho")
points(y=rep(as.numeric(fam.res_LMA.Narea[4]), times=7),x=seq(4.7,5.3, by=.1), pch=15, cex=.7)
points(y=rep(as.numeric(fam.resclean_LMA.Narea[4]), times=7),x=seq(5.7,6.3, by=.1), pch=15, cex=.7)
abline(h=0, lty=2)


## only things with reasonable correlations 
quartz(width=3, height=4.5)
par(mfrow=c(2,1), mar=c(4,4,2,1), oma=c(2,0,0,0))
boxplot(Slope_LMA.Narea~Type, all.results[which(abs(all.results$Rho_LMA.Narea)>.3),], ylim=c(-1.5,2),las=3, main="log(LMA)~log(Narea) Hcor", ylab="MA Slope")
points(y=rep(as.numeric(fam.res_LMA.Narea[3]), times=7),x=seq(4.7,5.3, by=.1), pch=15, cex=.7)
points(y=rep(as.numeric(fam.resclean_LMA.Narea[3]), times=7),x=seq(5.7,6.3, by=.1), pch=15, cex=.7)
abline(h=0, lty=2)
boxplot(Rho_LMA.Narea~Type, all.results[which(all.results$n_LMA.Narea>9),], ylim=c(-1,1),las=3, ylab="Rho")
points(y=rep(as.numeric(fam.res_LMA.Narea[4]), times=7),x=seq(4.7,5.3, by=.1), pch=15, cex=.7)
points(y=rep(as.numeric(fam.resclean_LMA.Narea[4]), times=7),x=seq(5.7,6.3, by=.1), pch=15, cex=.7)
abline(h=0, lty=2)



####### Plotting the LMA vs Narea scaling in unit rather than log space ######
plot(Narea~LMA, LES, col="grey", pch=16)
xs <- seq(1,1000, by=10)
## plotting w/in species slopes
for(i in which(all.results$Type=="w.inSpp")){
  ys <- 10^all.results$Int_LMA.Narea[i] * xs^ all.results$Slope_LMA.Narea[i]
  lines(ys~xs, col="darkblue", lwd=.5)
}
## plotting w/in genera slopes
for(i in which(all.results$Type=="w.inGen")){
  ys <- 10^all.results$Int_LMA.Narea[i] * xs^ all.results$Slope_LMA.Narea[i]
  lines(ys~xs, col="darkgreen")
}
## plotting w/in fams
for(i in which(all.results$Type=="Genw.inFam")){
  ys <- 10^all.results$Int_LMA.Narea[i] * xs^ all.results$Slope_LMA.Narea[i]
  lines(ys~xs, col="darkred")
}
## grand across family relationship
ys <- 10^all.results$Int_LMA.Narea[184] * xs^ all.results$Slope_LMA.Narea[184]
lines(ys~xs, lwd=3)

####### some analysis? ##########


###### Analyzing the correlations
LMALL <- lm(Rho_LMA.LL~Type, weights = n_LMA.LL, all.results)
# so Spp are significantly negative
  # and every level above spp is positive and significantly different from w/in.spp
LMAN <- lm(Rho_LMA.N~Type, weights = n_LMA.N, all.results)
  # they're all negative, but becoming increasingly negative at higher Taxo levels
lman <- aov(Rho_LMA.N~Type, all.results,weights=n_LMA.N)
TukeyHSD(lman)
  # everything different from w/in spp. but above spp not different


NmassLL <- lm(Rho_N.LL~Type, weights = n_N.LL, all.results)
  # w.inGen not different from w/.in spp, but everything else is...
nmassll <- aov(Rho_N.LL~Type, all.results,weights=n_N.LL)
TukeyHSD(nmassll)
  # with super conservative TukeyHSD, nothing is significantly different, though w/spp is close



###### Analyzing the MA slopes
LMALL <- lm(Slope_LMA.LL~Type, weights = n_LMA.LL, all.results)
# so Spp are significantly negative
# and every level above spp is positive and significantly different from w/in.spp



#### SLOPE LMA v LL
LMALL <- lm(Slope_LMA.LL~Type, all.results)
LMALLnull <- lm(Slope_LMA.LL~1, all.results)
anova(LMALL, LMALLnull) # p=0.044
LMALL <- lm(Slope_LMA.LL~Type, weights = n_LMA.LL, all.results)
LMALLnull <- lm(Slope_LMA.LL~1, weights = n_LMA.LL, all.results)
anova(LMALL, LMALLnull) # p<0.0001
LMALL <- lm(Slope_LMA.LL~Type, weights = varLL, all.results)
LMALLnull <- lm(Slope_LMA.LL~1, weights = varLL, all.results)
anova(LMALL, LMALLnull) # p=0.002
LMALL <- lm(Slope_LMA.LL~Type, weights = varLMA, all.results)
LMALLnull <- lm(Slope_LMA.LL~1, weights = varLMA, all.results)
anova(LMALL, LMALLnull) # p=0.031
# Rho LMA v Nmass
LMALL <- lm(Rho_LMA.LL~Type, all.results)
LMALLnull <- lm(Rho_LMA.LL~1, all.results)
anova(LMALL, LMALLnull) # p=0.007
LMALL <- lm(Rho_LMA.LL~Type, weights = n_LMA.LL, all.results)
LMALLnull <- lm(Rho_LMA.LL~1, weights = n_LMA.LL, all.results)
anova(LMALL, LMALLnull) # p<0.0001
LMALL <- lm(Rho_LMA.LL~Type, weights = varLL, all.results)
LMALLnull <- lm(Rho_LMA.LL~1, weights = varLL, all.results)
anova(LMALL, LMALLnull) # p<0.0001
LMALL <- lm(Rho_LMA.LL~Type, weights = varLMA, all.results)
LMALLnull <- lm(Rho_LMA.LL~1, weights = varLMA, all.results)
anova(LMALL, LMALLnull) # p=0.0019


#### SLOPE LMA v Nmass
LMAN <- lm(Slope_LMA.N~Type, all.results)
LMANnull <- lm(Slope_LMA.N~1, all.results)
anova(LMAN, LMANnull) # p=0.736
LMAN <- lm(Slope_LMA.N~Type, weights = n_LMA.N, all.results)
LMANnull <- lm(Slope_LMA.N~1, weights = n_LMA.N, all.results)
anova(LMAN, LMANnull) # p=0.96
LMAN <- lm(Slope_LMA.N~Type, weights = varNmass, all.results)
LMANnull <- lm(Slope_LMA.N~1, weights = varNmass, all.results)
anova(LMAN, LMANnull) # p=0.899
LMAN <- lm(Slope_LMA.N~Type, weights = varLMA, all.results)
LMANnull <- lm(Slope_LMA.N~1, weights = varLMA, all.results)
anova(LMAN, LMANnull) # p=0.437
# Rho LMA v Nmass
LMAN <- lm(Rho_LMA.N~Type, all.results)
LMANnull <- lm(Rho_LMA.N~1, all.results)
anova(LMAN, LMANnull) # p=0.0273
LMAN <- lm(Rho_LMA.N~Type, weights = n_LMA.N, all.results)
LMANnull <- lm(Rho_LMA.N~1, weights = n_LMA.N, all.results)
anova(LMAN, LMANnull) # p<0.001
LMAN <- lm(Rho_LMA.N~Type, weights = varNmass, all.results)
LMANnull <- lm(Rho_LMA.N~1, weights = varNmass, all.results)
anova(LMAN, LMANnull) # p=0.186
LMAN <- lm(Rho_LMA.N~Type, weights = varLMA, all.results)
LMANnull <- lm(Rho_LMA.N~1, weights = varLMA, all.results)
anova(LMAN, LMANnull) # p=0.0008



#### SLOPE LL v Nmass
NmassLL <- lm(Slope_N.LL~Type, all.results)
NmassLLnull <- lm(Slope_N.LL~1, all.results)
anova(NmassLL, NmassLLnull) # p=0.897
NmassLL <- lm(Slope_N.LL~Type, weights = n_N.LL, all.results)
NmassLLnull <- lm(Slope_N.LL~1, weights=n_N.LL, all.results)
anova(NmassLL, NmassLLnull) # p=0.0878
NmassLL <- lm(Slope_N.LL~Type, weights=varNmass, all.results)
NmassLLnull <- lm(Slope_N.LL~1, all.results, weights=varNmass)
anova(NmassLL, NmassLLnull) # p=0.61
NmassLL <- lm(Slope_N.LL~Type, weights=varLL, all.results)
NmassLLnull <- lm(Slope_N.LL~1, all.results, weights=varLL)
anova(NmassLL, NmassLLnull) # p=0.31
#### Rho LL v Nmass
NmassLL <- lm(Rho_N.LL~Type, all.results)
NmassLLnull <- lm(Rho_N.LL~1, all.results)
anova(NmassLL, NmassLLnull) # p=0.018
NmassLL <- lm(Rho_N.LL~Type, weights = n_N.LL, all.results)
NmassLLnull <- lm(Rho_N.LL~1, weights=n_N.LL, all.results)
anova(NmassLL, NmassLLnull) # p<0.0001
NmassLL <- lm(Rho_N.LL~Type, weights=varNmass, all.results)
NmassLLnull <- lm(Rho_N.LL~1, all.results, weights=varNmass)
anova(NmassLL, NmassLLnull) # p=0.005
NmassLL <- lm(Rho_N.LL~Type, weights=varLL, all.results)
NmassLLnull <- lm(Rho_N.LL~1, all.results, weights=varLL)
anova(NmassLL, NmassLLnull) # p<0.0001



#### SLOPE LMA v Narea
LMANarea <- lm(Slope_LMA.Narea~Type, all.results)
LMANareanull <- lm(Slope_LMA.Narea~1, all.results)
anova(LMANarea, LMANareanull) # p=0.0632
LMANarea <- lm(Slope_LMA.Narea~Type, weights = n_LMA.Narea, all.results)
LMANareanull <- lm(Slope_LMA.Narea~1, weights=n_LMA.Narea, all.results)
anova(LMANarea, LMANareanull) # p<0.0001
LMANarea <- lm(Slope_LMA.Narea~Type, weights=varNarea, all.results)
LMANareanull <- lm(Slope_LMA.Narea~1, all.results, weights=varNarea)
anova(LMANarea, LMANareanull) # p=0.549
LMANarea <- lm(Slope_LMA.Narea~Type, weights=varLMA, all.results)
LMANareanull <- lm(Slope_LMA.Narea~1, all.results, weights=varLMA)
anova(LMANarea, LMANareanull) # p=0.135
#### Rho LMA v Narea
LMANarea <- lm(Rho_LMA.Narea~Type, all.results)
LMANareanull <- lm(Rho_LMA.Narea~1, all.results)
anova(LMANarea, LMANareanull) # p=0.98
LMANarea <- lm(Rho_LMA.Narea~Type, weights = n_LMA.Narea, all.results)
LMANareanull <- lm(Rho_LMA.Narea~1, weights=n_LMA.Narea, all.results)
anova(LMANarea, LMANareanull) # p<0.45
LMANarea <- lm(Rho_LMA.Narea~Type, weights=varNarea, all.results)
LMANareanull <- lm(Rho_LMA.Narea~1, all.results, weights=varNarea)
anova(LMANarea, LMANareanull) # p=0.609
LMANarea <- lm(Rho_LMA.Narea~Type, weights=varLMA, all.results)
LMANareanull <- lm(Rho_LMA.Narea~1, all.results, weights=varLMA)
anova(LMANarea, LMANareanull) # p<0.20

# w.inGen not different from w/.in spp, but everything else is...
nmassll <- aov(Slope_N.LL~Type, all.results,weights=n_N.LL)
TukeyHSD(nmassll)
# with super conservative TukeyHSD, nothing is significantly different, though w/spp is close

LMA.Narea <- lm(Slope_LMA.Narea~Type, weights = n_LMA.Narea, all.results)
# Everything's different from w.inSpp. and increasingly different at higher Taxon scales...
lmanarea <- aov(Slope_LMA.Narea~Type, all.results)
TukeyHSD(lmanarea)
lmanarea <- aov(Slope_LMA.Narea~Type, all.results,weights=n_N.LL)
TukeyHSD(lmanarea)
lmanarea <- aov(Slope_LMA.Narea~Type, all.results,weights=varLMA)
TukeyHSD(lmanarea)
lmanarea <- aov(Slope_LMA.Narea~Type, all.results,weights=varNarea)
TukeyHSD(lmanarea)

### plotting two things at once...
ggplot(data=NULL,aes(x=all.results.LMAN$Rho, y=all.results.NmassLL$Rho, col=all.results.NmassLL$Type, size=all.results.NmassLL$n)) + geom_point()


ggplot(data=NULL,aes(x=all.results.LMAN$Rho, y=all.results.LMALL$Rho, col=all.results.NmassLL$Type, size=all.results.NmassLL$n)) + geom_point() +
  geom_abline(slope = 0, intercept=0, col="grey")+ geom_vline(xintercept = 0, col="grey")





###### PCA pieces ##########

## PCA on 655 records
pcLES <- prcomp(LES[-which(is.na(LES$log.LMA) | is.na(LES$log.LL) | is.na(LES$log.Nmass)),c("log.LMA","log.LL","log.Nmass")],center = T, scale. = T)
  # PC1 explains 78% of variance
pctraits <- prcomp(traits[-which(is.na(traits$log.LMA)|is.na(traits$log.LL)|is.na(traits$log.Nmass)),c("log.LMA","log.LL","log.Nmass") ],center = T, scale. = T)
  # of the full dataset (1010 trait measurements): PC1 explains 70% of the variance
pcPSME <- prcomp(traits.common.narm[which(traits.common.narm$SP.ID=="PSEMEN" & !is.na(traits.common.narm$log.Nmass)),c("log.LMA","log.LL","log.Nmass") ],center = T, scale. = T)
  # 220 trait measurements: 41%
pcPINPON <- prcomp(traits.common.narm[which(traits.common.narm$SP.ID=="PINPON" & !is.na(traits.common.narm$log.Nmass) & !is.na(traits.common.narm$log.LL)),c("log.LMA","log.LL","log.Nmass") ],center = T, scale. = T)
  # 180 trait measurements: 33%
pcTSUHET <- prcomp(traits.common.narm[which(traits.common.narm$SP.ID=="TSUHET" & !is.na(traits.common.narm$log.Nmass) & !is.na(traits.common.narm$log.LL)),c("log.LMA","log.LL","log.Nmass") ],center = T, scale. = T)
  # 60 trait measurements: 52%