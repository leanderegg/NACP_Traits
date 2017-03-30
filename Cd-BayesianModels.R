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












