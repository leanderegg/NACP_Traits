
#################################################################################################################
#________________________________________________________________________________________________________________
########################## *** +++AVERAGED RAW TRAITS+++ rather than averaged log.traits ****** ##############################
#________________________________________________________________________________________________________________
#################################################################################################################
##### Note: this was deprecated and moved to Cd-Anderegg2018_Manuscript_Analysis_v2.R #######
#****** This is now just a reminder of some of the original decisions made during analysis.


# This script is copied exactly from Cd-TaxonomicAnalysis.R.
# But I'm going to rename previously 'log.TRAIT' columns and switch rlog.TRAIT columns to be log.TRAIT so that we use those instead.

# I've also cut out my running commentary 

## as of 04.01.17 working with the full dataset already combined in Cd-Initial_Analysis for the Variance Decomp
commonspp <-  names(which(xtabs(~Species, data.all)>=5))
spp.data<- data.all %>% filter(Species %in% commonspp)
spp.data$Species <- factor(spp.data$Species)
spp.data$Genus <- factor(spp.data$Genus)
spp.data$Family <- factor(spp.data$Family)


########## Species level trait averages ########

allspp <- data.all %>% group_by(Species, Genus, Family) %>% summarise( lslog.LL = mean(log.LL, na.rm=T), lslog.LMA = mean(log.LMA, na.rm=T), lslog.Nmass = mean(log.Nmass, na.rm=T), lslog.Narea = mean(log.Narea, na.rm=T)
                                                                       ,slog.LL = log(mean(10^log.LL, na.rm=T),base=10), slog.LMA = log(mean(10^log.LMA, na.rm=T),base=10), slog.Nmass = log(mean(10^log.Nmass, na.rm=T),base=10), slog.Narea = log(mean(10^log.Narea, na.rm=T),base=10)
                                                                       ,MAT = mean(MAT, na.rm=T), MAP=mean(MAP, na.rm=T), VPD=mean(VPD, na.rm=T) )
colnames(allspp) <- gsub("slog", "log", colnames(allspp))




#### Genus mean dataframe ####
allgen <- allspp %>% group_by(Genus)  %>% summarise(Fam = sort(unique(Family))[1], lglog.LL = mean(llog.LL, na.rm=T),lglog.LMA = mean(llog.LMA, na.rm=T),lglog.Nmass = mean(llog.Nmass, na.rm=T), lglog.Narea = mean(llog.Narea, na.rm=T)
                                                    ,glog.LL = log(mean(10^log.LL, na.rm=T),base=10), glog.LMA = log(mean(10^log.LMA, na.rm=T),base=10), glog.Nmass = log(mean(10^log.Nmass, na.rm=T),base=10), glog.Narea = log(mean(10^log.Narea, na.rm=T),base=10), nspp = n() 
                                                    ,MAT = mean(MAT, na.rm=T), MAP=mean(MAP, na.rm=T), VPD=mean(VPD, na.rm=T) )
# taking mean of raw values or logged values doesn't matter all that much yet
colnames(allgen) <- gsub("glog", "log", colnames(allgen))
colnames(allgen)[2] <- "Family"
# 939 genera from 211 families
# allgen <- allspp %>% group_by(Genus)  %>% summarise(Fam = sort(unique(Family))[1], Needle.Broad = unique(na.omit(NeedleBroad))[1], glog.LL = mean(log.LL, na.rm=T), glog.LMA = mean(log.LMA, na.rm=T), glog.Nmass = mean(log.Nmass, na.rm=T),glog.Narea = mean(log.Narea, na.rm=T)
#                                                     ,rglog.LL = log(mean(10^rlog.LL, na.rm=T),base=10), rglog.LMA = log(mean(10^rlog.LMA, na.rm=T),base=10), rglog.Nmass = log(mean(10^rlog.Nmass, na.rm=T),base=10), nspp = n() )



#### Family Mean dataframe ####
allfam <- allgen %>% group_by(Family) %>% summarise(lflog.LL = mean(llog.LL, na.rm=T), lflog.LMA = mean(llog.LMA, na.rm=T), lflog.Nmass = mean(llog.Nmass, na.rm=T),lflog.Narea = mean(llog.Narea, na.rm=T)
                                                    ,flog.LL = log(mean(10^log.LL, na.rm=T),base=10), flog.LMA = log(mean(10^log.LMA, na.rm=T),base=10), flog.Nmass = log(mean(10^log.Nmass, na.rm=T),base=10), flog.Narea = log(mean(10^log.Narea, na.rm=T), base=10)
                                                    , tnspp = sum(nspp), ngen=n() )
colnames(allfam) <- gsub("flog", "log", colnames(allfam))
# still doesn't matter too too much. just moves some of the points towards larger values in rlog things


######## Seperating out replicated Genera and Families for taxonomic analysis ########################

##### Genus level data
# currently just working with LES until I combine the PACNW dataset into this.
#gen.data <- LESspp[which(LESspp$Genus %in% names(which(xtabs(~Genus,LESspp)>5))),]
# 634 measurements of 51 genera
# update 03.14.17, added traits spp to LES
# 659 measurements of 52 genera. Added Abies, but evidently that was the only additional genus...
# update 04.01.17, most up to date full dataset
gen.data <- allspp[which(allspp$Genus %in% names(which(xtabs(~Genus, allspp)>=5))),]
# 750 measurements of 73 genera.
gen.data$Genus <- factor(gen.data$Genus) # get rid of unused levels


##### Spp. in Family level data
# currently just working with LES until I combine the PACNW dataset into this.
# ppinfam.data <- LESspp[which(LESspp$Family %in% names(which(xtabs(~Family,LESspp)>5))),]
# 1374 measurements of 62 Families
# 1418 measurements of 63 Families. Added a family!
# 04.01.17 - 1786 measurements of 83 families. my taxonomic cleaning was v. helpful!
sppinfam.data <- allspp[which(allspp$Family %in% names(which(xtabs(~Family,allspp)>=5))),]
sppinfam.data$Family <- factor(sppinfam.data$Family) # get rid of unused levels



#### gen in Family level data
# 556 measurements of 40 Families , added 2 Families from last batch, and probably quite a few obs (old n_LMA.N = 459)
# 04.01.17 - 684 obs of 50 Families. again, big win from taxonomic cleaning!
geninfam.data <- allgen[which(allgen$Family %in% names(which(xtabs(~Family,allgen)>=5))),]
geninfam.data$Family <- factor(geninfam.data$Family) # get rid of unused levels
geninfam.data$Genus <- factor(geninfam.data$Genus)

#### Family Means:
fam.data <- allfam # new: 211 families (up from 189)
fam.dataclean <- allfam[which(allfam$tnspp>2),] # 101 families, up from 97 families


############# End: Dataset Creation #################################################




#_______________________________________________________________________

############ ***Piecewise hierarchical analysis*** ################
#_______________________________________________________________________

# while I'm trying to get the Bayesian approach up and running, just going to do a first pass with ML/piecewise correlations
# going to loop over species, then genera, then familes, etc. and calculate all correlations and Major Axis regression slopes
# then results are the distribution of all of those slopes/rhos at different taxo levels

# going to roll with traits.common5 for the moment, just so I have a larger sample size.

# # plot of all repped species w>5records in either LES or traits.common5
# ggplot(traits.common5, aes(x=log.Nmass, y=log.LMA, col=SP.ID)) + 
#   geom_point() + geom_smooth(method="lm", se = F) + 
#   geom_point(data=LES[which(LES$Species %in% commonspp),], aes(col=Species), size=2) + geom_smooth(data=LES[which(LES$Species %in% commonspp),], aes(col=Species), method="lm", se=F)
# 
# 
# # Narea vs LMA
# ggplot(spp.data, aes(x=log.LMA, y=log.Narea, col=Species)) + geom_point(data=LES, col="grey") + geom_point() + geom_smooth(method="lm", se=F) + geom_smooth(data=LES, aes(col=NULL), method="lm", se=F)
# #in arabidopsis
# ggplot(arab, aes(x=log(LMA), y=log(Narea), col=Type)) + geom_point()



#### Function for fitting MARs

fit.MAR <- function(xvar, yvar, data, method="SMA") {
  if(method=="SMA") meth = 3
  if(method=="MA") meth = 2
  if(method=="OLS") meth =1
  if(length(t(as.vector(data[!(is.na(data[,xvar]) | is.na(data[,yvar])), xvar])))<3){
    return(rep(NA, times=9))
  }
  else{
    if(var(data[,yvar], na.rm=T)==0 | var(data[,xvar],na.rm=T)==0){
      return(rep(NA, times=9))
    }
    else{
      tmp.mod <- suppressMessages(lmodel2(get(yvar)~get(xvar), data))
      intercept <- tmp.mod$regression.results$Intercept[meth]
      slope <- tmp.mod$regression.results$Slope[meth]
      slope.lci <- tmp.mod$confidence.intervals[meth,4]
      slope.uci <- tmp.mod$confidence.intervals[meth,5]
      rho <- tmp.mod$r
      r.sq <- tmp.mod$rsquare
      n <- tmp.mod$n
      varX <- var(data[,xvar], na.rm=T)
      varY <- var(data[,yvar], na.rm=T)
      results <- c(intercept, slope,slope.lci, slope.uci, rho, r.sq, n, varX, varY)
      return(results)
    }
  }
}


# from Cd-NullModel.R
test.sig <- function(x, test){
  if(x<test[1] | x>test[6]) return(0.025)
  else
    if(x<test[2] | x>test[5]) return(0.05)
  else
    if(x<test[3] | x>test[4]) return(0.1)
  else return(1)
}
# from Cd-NullModel.R - fits a null model based on the mean and range of the data
fit.null <- function(xvar, yvar, observed, nulldata, nits){
  # find the trait ranges to sample
  rangeX <- range(observed[,xvar], na.rm=T)
  difX <- rangeX[2]-rangeX[1]
  #meanX <- mean(observed[,xvar], na.rm=T)
  rangeY <- range(observed[,yvar], na.rm=T)
  difY <- rangeY[2]-rangeY[1]
  #meanY <- mean(observed[,yvar], na.rm=T)
  # cut down the null data to just have non NAs
  nulldata <- nulldata[which(!is.na(nulldata[,xvar]) & !is.na(nulldata[,yvar])),]
  #nulldata2 <- nulldata %>% filter(!is.na(xvar) & !is.na(yvar))
  restrictednull <- nulldata[which(nulldata[,xvar]> min(nulldata[,xvar], na.rm=T)+difX/2 & nulldata[,xvar]< max(nulldata[,xvar], na.rm=T)-difX/2 & nulldata[,yvar]> min(nulldata[,yvar], na.rm=T)+difY/2 & nulldata[,yvar]< max(nulldata[,yvar], na.rm=T)-difY/2),] 
  if(nrow(restrictednull)<3){restrictednull <- nulldata}
  nullcor <- c(rep(NA, times=nits))
  for(i in 1:nits){
    center <- restrictednull[sample(nrow(restrictednull),size = 1),]
    #null <- nulldata %>% filter(xvar > center[,xvar]-rangeX/2) #, xvar < center[,xvar]+rangeX/2, yvar > center[,yvar]-rangeY/2, yvar < center[,yvar]-rangeY/2)
    null <- nulldata[which(nulldata[,xvar] > as.numeric(center[,xvar])-difX/2 & nulldata[,xvar] < as.numeric(center[,xvar]+difX/2) & 
                             nulldata[,yvar] > as.numeric(center[,yvar])-difY/2 & nulldata[,yvar] < as.numeric(center[,yvar]+difY/2)),]
    if(nrow(null)<5){# var(null[,xvar])==0 | var(null[,yvar])==0){ # this is crappy as hell, but to get rid of bad centers I'm just going to try repicking them and hope the probability of getting two bad points in a row is low...
      center <- restrictednull[sample(nrow(restrictednull),size = 1),]
      #null <- nulldata %>% filter(xvar > center[,xvar]-rangeX/2) #, xvar < center[,xvar]+rangeX/2, yvar > center[,yvar]-rangeY/2, yvar < center[,yvar]-rangeY/2)
      null <- nulldata[which(nulldata[,xvar] > as.numeric(center[,xvar])-difX/2 & nulldata[,xvar] < as.numeric(center[,xvar]+difX/2) & 
                               nulldata[,yvar] > as.numeric(center[,yvar])-difY/2 & nulldata[,yvar] < as.numeric(center[,yvar]+difY/2)),]
    }
    ndist <- null[sample(nrow(null), size = nrow(observed), replace = TRUE),]
    nullcor[i] <- cor(ndist[,c(xvar,yvar)])[2,1]
  }
  nullquantiles <- quantile (nullcor, probs=c(0.025,0.05,0.1,.9,.95,.975), na.rm=T)
  names(nullquantiles) <- c("lci_2.5", "lci_5","lci_10","uci_10","uci_5","uci_2.5")
  return(nullquantiles)
}

rho.sig <- function(rho, n){
  t <- rho/sqrt((1-rho^2)/(n-2))
  pt(-abs(t), n-2,)*2
}




############# **LMA vs Nmass** ###########
############ .Species level analysis #####################
niters <- 10000

ptm <- proc.time()
set.seed(42)
spp.results <- data.frame(matrix(NA, nrow=length(unique(spp.data$Species)), ncol=17))
colnames(spp.results) <- c("Species", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varLMA","varNmass", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(spp.data$Species))){
  species <- levels(spp.data$Species)[i]
  print(species)
  dataz <- spp.data[which(spp.data$Species==species),]
  res <- fit.MAR(xvar='log.LMA',yvar="log.Nmass",data=dataz)
  spp.results[i,1] <- species
  spp.results[i,2:10] <- res
  if (!is.na(res[1]) & res[7]>4){ # only fit null model if there are >5 data points
    nullbounds <- fit.null(xvar='log.LMA', yvar="log.Nmass", observed = dataz, nulldata = allspp, nits = niters)  
    spp.results[i, 11:16] <- nullbounds
    spp.results[i,17] <- test.sig(x=spp.results$Rho[i], test=nullbounds)
  }
}
proc.time()-ptm





# Seems to be working!!!

############ .Genus level analysis #####################
t0 <- proc.time()
gen.results <- data.frame(matrix(NA, nrow=length(unique(gen.data$Genus)), ncol=17))
colnames(gen.results) <- c("Genus", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varLMA","varNmass", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(gen.data$Genus))){
  genus <- levels(gen.data$Genus)[i]
  print(genus)
  dataz <- gen.data[which(gen.data$Genus==genus),]
  gen.results[i,1] <- genus
  res <- fit.MAR(xvar='log.LMA',yvar="log.Nmass",data=dataz)
  gen.results[i,2:10] <- res 
  if (!is.na(res[1]) & res[7]>4){ # only fit null model if there are >5 data points, and the fit.MAR actually ran
    nullbounds <- fit.null(xvar='log.LMA', yvar="log.Nmass", observed = dataz, nulldata = allgen, nits = niters)  
    gen.results[i, 11:16] <- nullbounds
    gen.results[i,17] <- test.sig(x=gen.results$Rho[i], test=nullbounds)
  }
  
}
proc.time()-t0


############ .spp w/in Family level analysis #####################
t0<- proc.time()
sppinfam.results <- data.frame(matrix(NA, nrow=length(unique(sppinfam.data$Family)), ncol=17))
colnames(sppinfam.results) <- c("Family", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varLMA","varNmass", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(sppinfam.data$Family))){
  family <- levels(sppinfam.data$Family)[i]
  print(family)
  dataz <- sppinfam.data[which(sppinfam.data$Family==family),]
  res <- fit.MAR(xvar='log.LMA',yvar="log.Nmass",data=dataz)
  sppinfam.results[i,1] <- family
  sppinfam.results[i,2:10] <- res
  if (!is.na(res[1]) & res[7]>4){ # only fit null model if there are >5 data points
    nullbounds <- fit.null(xvar='log.LMA', yvar="log.Nmass", observed = dataz, nulldata = allfam, nits = niters)  
    sppinfam.results[i, 11:16] <- nullbounds
    sppinfam.results[i,17] <- test.sig(x=sppinfam.results$Rho[i], test=nullbounds)
  }
  
}
proc.time()-t0

############ .gen w/in Family level analysis #####################
t0<-proc.time()
geninfam.results <- data.frame(matrix(NA, nrow=length(unique(geninfam.data$Family)), ncol=17))
colnames(geninfam.results) <- c("Family", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varLMA","varNmass", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(geninfam.data$Family))){
  family <- levels(geninfam.data$Family)[i]
  print(family)
  dataz <- geninfam.data[which(geninfam.data$Family==family),]
  res <- fit.MAR(xvar='log.LMA',yvar="log.Nmass",data=dataz)
  geninfam.results[i,1] <- family
  geninfam.results[i,2:10] <- res
  if (!is.na(res[1]) & res[7]>4){ # only fit null model if there are >5 data points
    nullbounds <- fit.null(xvar='log.LMA', yvar="log.Nmass", observed = dataz, nulldata = allfam, nits = niters)  
    geninfam.results[i, 11:16] <- nullbounds
    geninfam.results[i,17] <- test.sig(x=geninfam.results$Rho[i], test=nullbounds)
  }
}
proc.time()-t0

############ .Family level analysis #####################

# # currently just working with LES until I combine the PACNW dataset into this.
# famcorr.all <- cor(x=fam.data[,c("log.LL","log.LMA","log.Nmass")], use = "pairwise.complete.obs")
# famcorr.clean <- cor(x=fam.dataclean[,c("log.LL","log.LMA","log.Nmass")], use = "pairwise.complete.obs")

fam.res_LMA.N <- c("fam.all", fit.MAR(xvar='log.LMA',yvar="log.Nmass",data=fam.data),rep(NA, times=7),"Fam") 
names(fam.res_LMA.N) <- c("Taxo.Unit","Int","Slope","Slope.lci","Slope.uci","Rho","r.sq", "n","varLMA","varNmass", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig" ,"Type")
fam.resclean_LMA.N <- c("fam.clean", fit.MAR(xvar='log.LMA',yvar="log.Nmass",data=fam.dataclean),rep(NA, times=7), "Fam.clean")
names(fam.resclean_LMA.N) <- c("Taxo.Unit","Int","Slope","Slope.lci","Slope.uci","Rho","r.sq", "n","varLMA","varNmass", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Type")
global_LMA.N <- c("global", fit.MAR(xvar='log.LMA',yvar="log.Nmass",data=allspp),rep(NA, times=7), "global")
names(global_LMA.N) <- c("Taxo.Unit","Int","Slope","Slope.lci","Slope.uci","Rho","r.sq", "n","varLMA","varNmass", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Type")


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
all.results.LMAN <-rbind(spp.results,gen.results, sppinfam.results, geninfam.results, fam.res_LMA.N, fam.resclean_LMA.N, global_LMA.N)
all.results.LMAN$Type <- factor(all.results.LMAN$Type)
levels(all.results.LMAN$Type) <- list(w.inSpp = "w/inSpp",   w.inGen= "w/inGen",    Sppw.inFam="Sppw/inFam", Genw.inFam="Genw/inFam", Fam = "Fam", Famclean = "Fam.clean", global="global")







####### ***LMA and LL*** #####################

############ .Species level analysis #####################


spp.results <- data.frame(matrix(NA, nrow=length(unique(spp.data$Species)), ncol=17))
colnames(spp.results) <- c("Species", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varLMA","varLL","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(spp.data$Species))){
  species <- levels(spp.data$Species)[i]
  print(species)
  dataz <- spp.data[which(spp.data$Species==species),]
  res <- fit.MAR(xvar='log.LMA',yvar="log.LL",data=dataz)
  spp.results[i,1] <- species
  spp.results[i,2:10] <- res
  if (!is.na(res[1]) &res[7]>4){ # only fit null model if there are >5 data points
    nullbounds <- fit.null(xvar='log.LMA', yvar="log.LL", observed = dataz, nulldata = allspp, nits = niters)  
    spp.results[i, 11:16] <- nullbounds
    spp.results[i,17] <- test.sig(x=spp.results$Rho[i], test=nullbounds)
  }
}


############ .Genus level analysis #####################


gen.results <- data.frame(matrix(NA, nrow=length(unique(gen.data$Genus)), ncol=17))
colnames(gen.results) <- c("Genus", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varLMA","varLL","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(gen.data$Genus))){
  genus <- levels(gen.data$Genus)[i]
  print(genus)
  dataz <- gen.data[which(gen.data$Genus==genus),]
  gen.results[i,1] <- genus
  res <- fit.MAR(xvar='log.LMA',yvar="log.LL",data=dataz)
  gen.results[i,2:10] <- res 
  if (!is.na(res[1]) & res[7]>4){ # only fit null model if there are >5 data points, and the fit.MAR actually ran
    nullbounds <- fit.null(xvar='log.LMA', yvar="log.LL", observed = dataz, nulldata = allgen, nits = niters)  
    gen.results[i, 11:16] <- nullbounds
    gen.results[i,17] <- test.sig(x=gen.results$Rho[i], test=nullbounds)
  }
}



############ .spp w/in Family level analysis #####################


sppinfam.results <- data.frame(matrix(NA, nrow=length(unique(sppinfam.data$Family)), ncol=17))
colnames(sppinfam.results) <- c("Family", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varLMA","varLL","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(sppinfam.data$Family))){
  family <- levels(sppinfam.data$Family)[i]
  print(family)
  dataz <- sppinfam.data[which(sppinfam.data$Family==family),]
  res <- fit.MAR(xvar='log.LMA',yvar="log.LL",data=dataz)
  sppinfam.results[i,1] <- family
  sppinfam.results[i,2:10] <- res
  if (!is.na(res[1]) & res[7]>4){ # only fit null model if there are >5 data points
    nullbounds <- fit.null(xvar='log.LMA', yvar="log.LL", observed = dataz, nulldata = allfam, nits = niters)  
    sppinfam.results[i, 11:16] <- nullbounds
    sppinfam.results[i,17] <- test.sig(x=sppinfam.results$Rho[i], test=nullbounds)
  }
}


############ .gen w/in Family level analysis #####################


geninfam.results <- data.frame(matrix(NA, nrow=length(unique(geninfam.data$Family)), ncol=17))
colnames(geninfam.results) <- c("Family", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varLMA","varLL","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(geninfam.data$Family))){
  family <- levels(geninfam.data$Family)[i]
  print(family)
  dataz <- geninfam.data[which(geninfam.data$Family==family),]
  res <- fit.MAR(xvar='log.LMA',yvar="log.LL",data=dataz)
  geninfam.results[i,1] <- family
  geninfam.results[i,2:10] <- res
  if (!is.na(res[1]) & res[7]>4){ # only fit null model if there are >5 data points
    nullbounds <- fit.null(xvar='log.LMA', yvar="log.LL", observed = dataz, nulldata = allfam, nits = niters)  
    geninfam.results[i, 11:16] <- nullbounds
    geninfam.results[i,17] <- test.sig(x=geninfam.results$Rho[i], test=nullbounds)
  }
}


############ .Family level analysis #####################

fam.res_LMA.LL <- c("fam.all", fit.MAR(xvar='log.LMA',yvar="log.LL",data=fam.data),rep(NA, times=7),"Fam")
names(fam.res_LMA.LL) <- c("Taxo.Unit","Int","Slope","Slope.lci","Slope.uci","Rho","r.sq", "n","varLMA","varLL","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Type")
fam.resclean_LMA.LL <- c("fam.clean", fit.MAR(xvar='log.LMA',yvar="log.LL",data=fam.dataclean),rep(NA, times=7), "Fam.clean")
names(fam.resclean_LMA.LL) <- c("Taxo.Unit","Int","Slope","Slope.lci","Slope.uci","Rho","r.sq", "n","varLMA","varLL","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Type")
global_LMA.LL <- c("global", fit.MAR(xvar='log.LMA',yvar="log.LL",data=allspp),rep(NA, times=7), "global")
names(global_LMA.LL) <- c("Taxo.Unit","Int","Slope","Slope.lci","Slope.uci","Rho","r.sq", "n","varLMA","varLL","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Type")



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
all.results.LMALL <-rbind(spp.results,gen.results, sppinfam.results, geninfam.results, fam.res_LMA.LL, fam.resclean_LMA.LL, global_LMA.LL)
all.results.LMALL$Type <- factor(all.results.LMALL$Type)
levels(all.results.LMALL$Type) <- list(w.inSpp = "w/inSpp",   w.inGen= "w/inGen",    Sppw.inFam="Sppw/inFam", Genw.inFam="Genw/inFam", Fam = "Fam", Famclean = "Fam.clean", global="global")





####### ***LL and Nmass*** #####################



spp.results <- data.frame(matrix(NA, nrow=length(unique(spp.data$Species)), ncol=17))
colnames(spp.results) <- c("Species", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varLL","varNmass","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(spp.data$Species))){
  species <- levels(spp.data$Species)[i]
  print(species)
  dataz <- spp.data[which(spp.data$Species==species),]
  res <- fit.MAR(xvar='log.LL',yvar="log.Nmass",data=dataz)
  spp.results[i,1] <- species
  spp.results[i,2:10] <- res
  if (!is.na(res[1]) &res[7]>4){ # only fit null model if there are >5 data points
    nullbounds <- fit.null(xvar='log.LL', yvar="log.Nmass", observed = dataz, nulldata = allspp, nits = niters)  
    spp.results[i, 11:16] <- nullbounds
    spp.results[i,17] <- test.sig(x=spp.results$Rho[i], test=nullbounds)
  }
}


############ .Genus level analysis #####################

gen.results <- data.frame(matrix(NA, nrow=length(unique(gen.data$Genus)), ncol=17))
colnames(gen.results) <- c("Genus", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varLL","varNmass","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(gen.data$Genus))){
  genus <- levels(gen.data$Genus)[i]
  print(genus)
  dataz <- gen.data[which(gen.data$Genus==genus),]
  gen.results[i,1] <- genus
  res <- fit.MAR(xvar='log.LL',yvar="log.Nmass",data=dataz)
  gen.results[i,2:10] <- res 
  if (!is.na(res[1]) & res[7]>4){ # only fit null model if there are >5 data points, and the fit.MAR actually ran
    nullbounds <- fit.null(xvar='log.LL', yvar="log.Nmass", observed = dataz, nulldata = allgen, nits = niters)  
    gen.results[i, 11:16] <- nullbounds
    gen.results[i,17] <- test.sig(x=gen.results$Rho[i], test=nullbounds)
  }
}



############ .spp w/in Family level analysis #####################

sppinfam.results <- data.frame(matrix(NA, nrow=length(unique(sppinfam.data$Family)), ncol=17))
colnames(sppinfam.results) <- c("Family", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varLL","varNmass","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(sppinfam.data$Family))){
  family <- levels(sppinfam.data$Family)[i]
  print(family)
  dataz <- sppinfam.data[which(sppinfam.data$Family==family),]
  res <- fit.MAR(xvar='log.LL',yvar="log.Nmass",data=dataz)
  sppinfam.results[i,1] <- family
  sppinfam.results[i,2:10] <- res
  if (!is.na(res[1]) & res[7]>4){ # only fit null model if there are >5 data points
    nullbounds <- fit.null(xvar='log.LL', yvar="log.Nmass", observed = dataz, nulldata = allfam, nits = niters)  
    sppinfam.results[i, 11:16] <- nullbounds
    sppinfam.results[i,17] <- test.sig(x=sppinfam.results$Rho[i], test=nullbounds)
  }
}





############ .gen w/in Family level analysis #####################
geninfam.results <- data.frame(matrix(NA, nrow=length(unique(geninfam.data$Family)), ncol=17))
colnames(geninfam.results) <- c("Family", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varLL","varNmass","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(geninfam.data$Family))){
  family <- levels(geninfam.data$Family)[i]
  print(family)
  dataz <- geninfam.data[which(geninfam.data$Family==family),]
  res <- fit.MAR(xvar='log.LL',yvar="log.Nmass",data=dataz)
  geninfam.results[i,1] <- family
  geninfam.results[i,2:10] <- res
  if (!is.na(res[1]) & res[7]>4){ # only fit null model if there are >5 data points
    nullbounds <- fit.null(xvar='log.LL', yvar="log.Nmass", observed = dataz, nulldata = allfam, nits = niters)  
    geninfam.results[i, 11:16] <- nullbounds
    geninfam.results[i,17] <- test.sig(x=geninfam.results$Rho[i], test=nullbounds)
  }
}





############ .Family level analysis #####################


fam.res_LL.N <- c("fam.all", fit.MAR(xvar='log.LL',yvar="log.Nmass",data=fam.data),rep(NA, times=7),"Fam")
names(fam.res_LL.N) <- c("Taxo.Unit","Int","Slope","Slope.lci","Slope.uci","Rho","r.sq", "n","varLL","varNmass","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Type")
fam.resclean_LL.N <- c("fam.clean", fit.MAR(xvar='log.LL',yvar="log.Nmass",data=fam.dataclean),rep(NA, times=7), "Fam.clean")
names(fam.res_LL.N) <- c("Taxo.Unit","Int","Slope","Slope.lci","Slope.uci","Rho","r.sq", "n","varLL","varNmass","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Type")
global_LL.N <- c("global", fit.MAR(xvar='log.LL',yvar="log.Nmass",data=allspp),rep(NA, times=7), "global")
names(global_LL.N) <- c("Taxo.Unit","Int","Slope","Slope.lci","Slope.uci","Rho","r.sq", "n","varLL","varNmass","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Type")





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
all.results.LLNmass <-rbind(spp.results,gen.results, sppinfam.results, geninfam.results, fam.res_LL.N, fam.resclean_LL.N, global_LL.N)
all.results.LLNmass$Type <- factor(all.results.LLNmass$Type)
levels(all.results.LLNmass$Type) <- list(w.inSpp = "w/inSpp",   w.inGen= "w/inGen",    Sppw.inFam="Sppw/inFam", Genw.inFam="Genw/inFam", Fam = "Fam", Famclean = "Fam.clean", global="global")






####### ***LMA and Narea!!!*** #####################



spp.results <- data.frame(matrix(NA, nrow=length(unique(spp.data$Species)), ncol=17))
colnames(spp.results) <- c("Species", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varLMA","varNarea","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(spp.data$Species))){
  species <- levels(spp.data$Species)[i]
  print(species)
  dataz <- spp.data[which(spp.data$Species==species),]
  res <- fit.MAR(xvar='log.LMA',yvar="log.Narea",data=dataz)
  spp.results[i,1] <- species
  spp.results[i,2:10] <- res
  if (!is.na(res[1]) &res[7]>4){ # only fit null model if there are >5 data points
    nullbounds <- fit.null(xvar='log.LMA', yvar="log.Narea", observed = dataz, nulldata = allspp, nits = niters)  
    spp.results[i, 11:16] <- nullbounds
    spp.results[i,17] <- test.sig(x=spp.results$Rho[i], test=nullbounds)
  }
}


############ .Genus level analysis #####################

gen.results <- data.frame(matrix(NA, nrow=length(unique(gen.data$Genus)), ncol=17))
colnames(gen.results) <- c("Genus", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varLMA","varNarea","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(gen.data$Genus))){
  genus <- levels(gen.data$Genus)[i]
  print(genus)
  dataz <- gen.data[which(gen.data$Genus==genus),]
  gen.results[i,1] <- genus
  res <- fit.MAR(xvar='log.LMA',yvar="log.Narea",data=dataz)
  gen.results[i,2:10] <- res 
  if (!is.na(res[1]) & res[7]>4){ # only fit null model if there are >5 data points, and the fit.MAR actually ran
    nullbounds <- fit.null(xvar='log.LMA', yvar="log.Narea", observed = dataz, nulldata = allgen, nits = niters)  
    gen.results[i, 11:16] <- nullbounds
    gen.results[i,17] <- test.sig(x=gen.results$Rho[i], test=nullbounds)
  }
}



############ .spp w/in Family level analysis #####################

sppinfam.results <- data.frame(matrix(NA, nrow=length(unique(sppinfam.data$Family)), ncol=17))
colnames(sppinfam.results) <- c("Family", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varLMA","varNarea","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(sppinfam.data$Family))){
  family <- levels(sppinfam.data$Family)[i]
  print(family)
  dataz <- sppinfam.data[which(sppinfam.data$Family==family),]
  res <- fit.MAR(xvar='log.LMA',yvar="log.Narea",data=dataz)
  sppinfam.results[i,1] <- family
  sppinfam.results[i,2:10] <- res
  if (!is.na(res[1]) & res[7]>4){ # only fit null model if there are >5 data points
    nullbounds <- fit.null(xvar='log.LMA', yvar="log.Narea", observed = dataz, nulldata = allfam, nits = niters)  
    sppinfam.results[i, 11:16] <- nullbounds
    sppinfam.results[i,17] <- test.sig(x=sppinfam.results$Rho[i], test=nullbounds)
  }
}





############ .gen w/in Family level analysis #####################

geninfam.results <- data.frame(matrix(NA, nrow=length(unique(geninfam.data$Family)), ncol=17))
colnames(geninfam.results) <- c("Family", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varLMA","varNarea","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(geninfam.data$Family))){
  family <- levels(geninfam.data$Family)[i]
  print(family)
  dataz <- geninfam.data[which(geninfam.data$Family==family),]
  res <- fit.MAR(xvar='log.LMA',yvar="log.Narea",data=dataz)
  geninfam.results[i,1] <- family
  geninfam.results[i,2:10] <- res
  if (!is.na(res[1]) & res[7]>4){ # only fit null model if there are >5 data points
    nullbounds <- fit.null(xvar='log.LMA', yvar="log.Narea", observed = dataz, nulldata = allfam, nits = niters)  
    geninfam.results[i, 11:16] <- nullbounds
    geninfam.results[i,17] <- test.sig(x=geninfam.results$Rho[i], test=nullbounds)
  }
}





############ .Family level analysis #####################


fam.res_LMA.Narea <- c("fam.all", fit.MAR(xvar='log.LMA',yvar="log.Narea",data=fam.data),rep(NA, times=7),"Fam")
names(fam.res_LMA.Narea) <- c("Taxo.Unit","Int","Slope","Slope.lci","Slope.uci","Rho","r.sq", "n","varLMA","varNarea","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Type")
fam.resclean_LMA.Narea <- c("fam.clean", fit.MAR(xvar='log.LMA',yvar="log.Narea",data=fam.dataclean),rep(NA, times=7), "Fam.clean")
names(fam.res_LMA.Narea) <- c("Taxo.Unit","Int","Slope","Slope.lci","Slope.uci","Rho","r.sq", "n","varLMA","varNarea","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Type")
global_LMA.Narea <- c("global", fit.MAR(xvar='log.LMA',yvar="log.Narea",data=allspp),rep(NA, times=7), "global")
names(global_LMA.Narea) <- c("Taxo.Unit","Int","Slope","Slope.lci","Slope.uci","Rho","r.sq", "n","varLMA","varNarea","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Type")




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
all.results.LMANarea <-rbind(spp.results,gen.results, sppinfam.results, geninfam.results, fam.res_LMA.Narea, fam.resclean_LMA.Narea, global_LMA.Narea)
all.results.LMANarea$Type <- factor(all.results.LMANarea$Type)
levels(all.results.LMANarea$Type) <- list(w.inSpp = "w/inSpp",   w.inGen= "w/inGen",    Sppw.inFam="Sppw/inFam", Genw.inFam="Genw/inFam", Fam = "Fam", Famclean = "Fam.clean", global="global")






####### ***LL and Narea*** #####################



spp.results <- data.frame(matrix(NA, nrow=length(unique(spp.data$Species)), ncol=17))
colnames(spp.results) <- c("Species", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varLL","varNarea","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(spp.data$Species))){
  species <- levels(spp.data$Species)[i]
  print(species)
  dataz <- spp.data[which(spp.data$Species==species),]
  res <- fit.MAR(xvar="log.LL", yvar='log.Narea',data=dataz)
  spp.results[i,1] <- species
  spp.results[i,2:10] <- res
  if (!is.na(res[1]) &res[7]>4){ # only fit null model if there are >5 data points
    nullbounds <- fit.null(xvar='log.LL', yvar="log.Narea", observed = dataz, nulldata = allspp, nits = niters)  
    spp.results[i, 11:16] <- nullbounds
    spp.results[i,17] <- test.sig(x=spp.results$Rho[i], test=nullbounds)
  }
}


############ .Genus level analysis #####################

gen.results <- data.frame(matrix(NA, nrow=length(unique(gen.data$Genus)), ncol=17))
colnames(gen.results) <- c("Genus", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varLL","varNarea","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(gen.data$Genus))){
  genus <- levels(gen.data$Genus)[i]
  print(genus)
  dataz <- gen.data[which(gen.data$Genus==genus),]
  gen.results[i,1] <- genus
  res <- fit.MAR(xvar="log.LL",yvar='log.Narea',data=dataz)
  gen.results[i,2:10] <- res 
  if (!is.na(res[1]) & res[7]>4){ # only fit null model if there are >5 data points, and the fit.MAR actually ran
    nullbounds <- fit.null(xvar='log.LL', yvar="log.Narea", observed = dataz, nulldata = allgen, nits = niters)  
    gen.results[i, 11:16] <- nullbounds
    gen.results[i,17] <- test.sig(x=gen.results$Rho[i], test=nullbounds)
  }
}



############ .spp w/in Family level analysis #####################

sppinfam.results <- data.frame(matrix(NA, nrow=length(unique(sppinfam.data$Family)), ncol=17))
colnames(sppinfam.results) <- c("Family", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varLL","varNarea","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(sppinfam.data$Family))){
  family <- levels(sppinfam.data$Family)[i]
  print(family)
  dataz <- sppinfam.data[which(sppinfam.data$Family==family),]
  res <- fit.MAR(xvar='log.LL',yvar="log.Narea",data=dataz)
  sppinfam.results[i,1] <- family
  sppinfam.results[i,2:10] <- res
  if (!is.na(res[1]) & res[7]>4){ # only fit null model if there are >5 data points
    nullbounds <- fit.null(xvar='log.LL', yvar="log.Narea", observed = dataz, nulldata = allfam, nits = niters)  
    sppinfam.results[i, 11:16] <- nullbounds
    sppinfam.results[i,17] <- test.sig(x=sppinfam.results$Rho[i], test=nullbounds)
  }
}





############ .gen w/in Family level analysis #####################

geninfam.results <- data.frame(matrix(NA, nrow=length(unique(geninfam.data$Family)), ncol=17))
colnames(geninfam.results) <- c("Family", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varLL","varNarea","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(geninfam.data$Family))){
  family <- levels(geninfam.data$Family)[i]
  print(family)
  dataz <- geninfam.data[which(geninfam.data$Family==family),]
  res <- fit.MAR(xvar='log.LL',yvar="log.Narea",data=dataz)
  geninfam.results[i,1] <- family
  geninfam.results[i,2:10] <- res
  if (!is.na(res[1]) & res[7]>4){ # only fit null model if there are >5 data points
    nullbounds <- fit.null(xvar='log.LL', yvar="log.Narea", observed = dataz, nulldata = allfam, nits = niters)  
    geninfam.results[i, 11:16] <- nullbounds
    geninfam.results[i,17] <- test.sig(x=geninfam.results$Rho[i], test=nullbounds)
  }
}





############ .Family level analysis #####################


fam.res_LL.Narea <- c("fam.all", fit.MAR(xvar="log.LL",yvar='log.Narea',data=fam.data),rep(NA, times=7),"Fam")
names(fam.res_LL.Narea) <- c("Taxo.Unit","Int","Slope","Slope.lci","Slope.uci","Rho","r.sq", "n","varLL","varNarea","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Type")
fam.resclean_LL.Narea <- c("fam.clean", fit.MAR(xvar="log.LL",yvar='log.Narea',data=fam.dataclean),rep(NA, times=7), "Fam.clean")
names(fam.res_LL.Narea) <- c("Taxo.Unit","Int","Slope","Slope.lci","Slope.uci","Rho","r.sq", "n","varLL","varNarea","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Type")
global_LL.Narea <- c("global", fit.MAR(xvar="log.LL",yvar='log.Narea',data=allspp),rep(NA, times=7), "global")
names(global_LL.Narea) <- c("Taxo.Unit","Int","Slope","Slope.lci","Slope.uci","Rho","r.sq", "n","varLL","varNarea","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Type")






############ CMW analysis #####################

# note: need to load biomass in Cd-Initial_Analysis and also run the code in Cd-CommunityWeightedMeans


CWM_LMA.LL <- c("CWM", fit.MAR(xvar='log.cw_LMAp_if',yvar="log.cw_LLp_if",data=biomass),rep(NA, times=7), "CWM")
names(CWM_LMA.LL) <- c("Taxo.Unit","Int","Slope","Slope.lci","Slope.uci","Rho","r.sq", "n","varLMA","varLL","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Type")


CWM_LMA.N <- c("CWM", fit.MAR(xvar='log.cw_LMAp_if',yvar="log.cw_Nmassp_if",data=biomass),rep(NA, times=7), "CWM")
names(CWM_LMA.N) <- c("Taxo.Unit","Int","Slope","Slope.lci","Slope.uci","Rho","r.sq", "n","varLMA","varNmass", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Type")

CWM_LL.N <- c("CWM", fit.MAR(xvar='log.cw_LLp_if',yvar="log.cw_Nmassp_if",data=biomass),rep(NA, times=7), "CWM")
names(CWM_LL.N) <- c("Taxo.Unit","Int","Slope","Slope.lci","Slope.uci","Rho","r.sq", "n","varLL","varNmass","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Type")

CWM_LMA.Narea <- c("CWM", fit.MAR(xvar='log.cw_LMAp_if',yvar="log.cw_Nareap_if",data=biomass),rep(NA, times=7), "CWM")
names(CWM_LMA.Narea) <- c("Taxo.Unit","Int","Slope","Slope.lci","Slope.uci","Rho","r.sq", "n","varLMA","varNarea","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Type")

CWM_LL.Narea <- c("CWM", fit.MAR(xvar="log.cw_LLp_if",yvar='log.cw_Nareap_if',data=biomass),rep(NA, times=7), "CWM")
names(CWM_LL.Narea) <- c("Taxo.Unit","Int","Slope","Slope.lci","Slope.uci","Rho","r.sq", "n","varLL","varNarea","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Type")







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
all.results.LLNarea <-rbind(spp.results,gen.results, sppinfam.results, geninfam.results, fam.res_LL.Narea, fam.resclean_LL.Narea, global_LL.Narea)
all.results.LLNarea$Type <- factor(all.results.LLNarea$Type)
levels(all.results.LLNarea$Type) <- list(w.inSpp = "w/inSpp",   w.inGen= "w/inGen",    Sppw.inFam="Sppw/inFam", Genw.inFam="Genw/inFam", Fam = "Fam", Famclean = "Fam.clean", global="global")

















########### ** ==Combining ALL ANALYSES into 1 df == ** ######

LMALL <- all.results.LMALL
colnames(LMALL)[c(2:8, 11:17)] <- paste(colnames(all.results.LMALL)[c(2:8, 11:17)],"LMA.LL", sep="_")
LMAN <- all.results.LMAN
colnames(LMAN)[c(2:8, 11:17)] <- paste(colnames(all.results.LMAN)[c(2:8, 11:17)],"LMA.N", sep="_")
LLNmass <- all.results.LLNmass
colnames(LLNmass)[c(2:8, 11:17)] <- paste(colnames(all.results.LLNmass)[c(2:8, 11:17)],"LL.N", sep="_")
LMANarea <- all.results.LMANarea
colnames(LMANarea)[c(2:8, 11:17)] <- paste(colnames(all.results.LMANarea)[c(2:8, 11:17)],"LMA.Narea", sep="_")
LLNarea <- all.results.LLNarea
colnames(LLNarea)[c(2:8, 11:17)] <- paste(colnames(all.results.LLNarea)[c(2:8, 11:17)],"LL.Narea", sep="_")

all.results <- cbind(LMALL, LMAN[,-c(1,9,18)], LLNmass[,-c(1,9,10,18)], LMANarea[,-c(1,9,18)], LLNarea[,-c(1,9,10,18)]) # drop the duplicate 'var' columns


### add in CWMs to dataframe post hoc
allcwms <- c(CWM_LMA.LL, CWM_LMA.N[-c(1,9,18)], CWM_LL.N[-c(1,9,10,18)], CWM_LMA.Narea[-c(1,9,18)], CWM_LL.Narea[-c(1,9,10,18)])
test <- all.results
test$Taxo.Unit <- as.character(test$Taxo.Unit)
test$Type <- as.character(test$Type)
test[nrow(test)+1,] <- allcwms

all.results.cwm <- test


#---> addition of Mt Rainier species deprecated as of 12.20.17 because code has been updated to include them from the beginning

# ### add in Mt Rainier species as of 12.11.17 
# # takes all.results and adds in 6 rows in w/in spp for Mt Rainier.
# # did this so I didn't have to completely rerun the entire analysis above. 
# # Instead, just make the Rainier.results dfs in Cd-Mt_Rainier_LMA_creation
# 
# empty <-  data.frame(matrix(NA, nrow = 6, ncol=ncol(all.results.cwm)))
# names(empty) <- names(all.results.cwm)
# empty$Taxo.Unit <- c("Aster alpigenus","Castilleja parviflora", "Erythronium montanum", "Lupinus arcticus","Valeriana sitchensis","Veratrum viride")
# empty$Type <- "w.inSpp"
# all.resultstest <- rbind(all.results[which(all.results$Type=="w.inSpp"),], empty, all.results[which(all.results$Type!="w.inSpp"),])
# all.resultstest[which(all.resultstest$Taxo.Unit %in% c("Aster alpigenus","Castilleja parviflora", "Erythronium montanum", "Lupinus arcticus","Valeriana sitchensis","Veratrum viride")), grep("LMA.LL", colnames(all.resultstest))] <- Rain.results.LMALL[,which(!colnames(Rain.results.LMALL) %in% c("Species","varLMA","varLL","rho.sig"))]
# all.resultstest[which(all.resultstest$Taxo.Unit %in% c("Aster alpigenus","Castilleja parviflora", "Erythronium montanum", "Lupinus arcticus","Valeriana sitchensis","Veratrum viride")), c("varLMA","varLL")] <- Rain.results.LMALL[,c("varLMA","varLL")]
# all.resultstest[which(all.resultstest$Taxo.Unit %in% c("Aster alpigenus","Castilleja parviflora", "Erythronium montanum", "Lupinus arcticus","Valeriana sitchensis","Veratrum viride")), grep("LMA.N", colnames(all.resultstest))[1:14]] <- Rain.results.LMAN[,which(!colnames(Rain.results.LMAN) %in% c("Species","varLMA","varNmass","rho.sig"))]
# all.resultstest[which(all.resultstest$Taxo.Unit %in% c("Aster alpigenus","Castilleja parviflora", "Erythronium montanum", "Lupinus arcticus","Valeriana sitchensis","Veratrum viride")), c("varNmass")] <- Rain.results.LMAN[,c("varNmass")]
# all.resultstest[which(all.resultstest$Taxo.Unit %in% c("Aster alpigenus","Castilleja parviflora", "Erythronium montanum", "Lupinus arcticus","Valeriana sitchensis","Veratrum viride")), grep("LMA.Narea", colnames(all.resultstest))] <- Rain.results.LMANarea[,which(!colnames(Rain.results.LMANarea) %in% c("Species","varLMA","varNarea","rho.sig"))]
# all.resultstest[which(all.resultstest$Taxo.Unit %in% c("Aster alpigenus","Castilleja parviflora", "Erythronium montanum", "Lupinus arcticus","Valeriana sitchensis","Veratrum viride")), c("varNarea")] <- Rain.results.LMANarea[,c("varNarea")]
# all.resultstest[which(all.resultstest$Taxo.Unit %in% c("Aster alpigenus","Castilleja parviflora", "Erythronium montanum", "Lupinus arcticus","Valeriana sitchensis","Veratrum viride")), grep("LL.N", colnames(all.resultstest))[1:14]] <- Rain.results.LLN[,which(!colnames(Rain.results.LLN) %in% c("Species","varLL","varNmass","rho.sig"))]
# all.resultstest[which(all.resultstest$Taxo.Unit %in% c("Aster alpigenus","Castilleja parviflora", "Erythronium montanum", "Lupinus arcticus","Valeriana sitchensis","Veratrum viride")), grep("LL.Narea", colnames(all.resultstest))] <- Rain.results.LLNarea[,which(!colnames(Rain.results.LLNarea) %in% c("Species","varLL","varNarea","rho.sig"))]
# 
# 
# ### add into .cwm dataframe
# all.results.cwmtest <- rbind(all.results.cwm[which(all.results.cwm$Type=="w.inSpp"),], empty, all.results.cwm[which(all.results.cwm$Type!="w.inSpp"),])
# all.results.cwmtest[which(all.results.cwmtest$Taxo.Unit %in% c("Aster alpigenus","Castilleja parviflora", "Erythronium montanum", "Lupinus arcticus","Valeriana sitchensis","Veratrum viride")), grep("LMA.LL", colnames(all.results.cwmtest))] <- Rain.results.LMALL[,which(!colnames(Rain.results.LMALL) %in% c("Species","varLMA","varLL","rho.sig"))]
# all.results.cwmtest[which(all.results.cwmtest$Taxo.Unit %in% c("Aster alpigenus","Castilleja parviflora", "Erythronium montanum", "Lupinus arcticus","Valeriana sitchensis","Veratrum viride")), c("varLMA","varLL")] <- Rain.results.LMALL[,c("varLMA","varLL")]
# all.results.cwmtest[which(all.results.cwmtest$Taxo.Unit %in% c("Aster alpigenus","Castilleja parviflora", "Erythronium montanum", "Lupinus arcticus","Valeriana sitchensis","Veratrum viride")), grep("LMA.N", colnames(all.results.cwmtest))[1:14]] <- Rain.results.LMAN[,which(!colnames(Rain.results.LMAN) %in% c("Species","varLMA","varNmass","rho.sig"))]
# all.results.cwmtest[which(all.results.cwmtest$Taxo.Unit %in% c("Aster alpigenus","Castilleja parviflora", "Erythronium montanum", "Lupinus arcticus","Valeriana sitchensis","Veratrum viride")), c("varNmass")] <- Rain.results.LMAN[,c("varNmass")]
# all.results.cwmtest[which(all.results.cwmtest$Taxo.Unit %in% c("Aster alpigenus","Castilleja parviflora", "Erythronium montanum", "Lupinus arcticus","Valeriana sitchensis","Veratrum viride")), grep("LMA.Narea", colnames(all.results.cwmtest))] <- Rain.results.LMANarea[,which(!colnames(Rain.results.LMANarea) %in% c("Species","varLMA","varNarea","rho.sig"))]
# all.results.cwmtest[which(all.results.cwmtest$Taxo.Unit %in% c("Aster alpigenus","Castilleja parviflora", "Erythronium montanum", "Lupinus arcticus","Valeriana sitchensis","Veratrum viride")), c("varNarea")] <- Rain.results.LMANarea[,c("varNarea")]
# all.results.cwmtest[which(all.results.cwmtest$Taxo.Unit %in% c("Aster alpigenus","Castilleja parviflora", "Erythronium montanum", "Lupinus arcticus","Valeriana sitchensis","Veratrum viride")), grep("LL.N", colnames(all.results.cwmtest))[1:14]] <- Rain.results.LLN[,which(!colnames(Rain.results.LLN) %in% c("Species","varLL","varNmass","rho.sig"))]
# all.results.cwmtest[which(all.results.cwmtest$Taxo.Unit %in% c("Aster alpigenus","Castilleja parviflora", "Erythronium montanum", "Lupinus arcticus","Valeriana sitchensis","Veratrum viride")), grep("LL.Narea", colnames(all.results.cwmtest))] <- Rain.results.LLNarea[,which(!colnames(Rain.results.LLNarea) %in% c("Species","varLL","varNarea","rho.sig"))]



#write.csv(all.results, "Results_SimpleMAreg_v1_030817.csv")
#write.csv(all.results, "Results_SimpleMAreg_v2_031417.csv")
#write.csv(all.results, "Results_SimpleMAreg_v3_031717.csv")
#write.csv(all.results, "Results_SimpleMAreg_v4_040117.csv")
#write.csv(all.results, "Results_SimpleMAreg_v5_040217.csv") # updated with NareaLL
#write.csv(all.results, "Results_SimpleMAreg_v6_040717.csv") # switched to LLNarea
#write.csv(all.results, "Results_SimpleMAreg_v7_051717.csv") # switched to LLNmass, and added CIs from null model
#write.csv(all.results, "Results_SimpleMAreg_v8_051717.csv") # fixed bug that screwed up LL.Narea null model (when taxa had more variance than all families)
#write.csv(all.results, "Results_SimpleMAreg_v9rawavgs_20170620.csv")
#write.csv(all.results, "Results_SimpleMAreg_v9rawavgs_20170828_wCWM.csv") # with CWMs now added
# write.csv(all.results, "Results_SimpleMAreg_v10rawavgs_20171120.csv") # added CIs for SMAs.
# write.csv(all.results.cwm, "/Users/leeanderegg/Dropbox/NACP_Traits/NACP_Traits_Rcode/Results_SimpleMAreg_v10rawavgs_20171120_wCWM.csv") # added CIs for SMAs.
# write.csv(all.resultstest, "Results_SimpleMAreg_v11rawavgs_20171211.csv") # added CIs for SMAs.
# write.csv(all.results.cwmtest, "/Users/leeanderegg/Dropbox/NACP_Traits/NACP_Traits_Rcode/Results_SimpleMAreg_v11rawavgs_20171211_wCWM.csv") # added CIs for SMAs.
 write.csv(all.results, "Results_SimpleMAreg_v12rawavgs_20171220.csv") # added CIs for SMAs.
 write.csv(all.results.cwm, "/Users/leeanderegg/Dropbox/NACP_Traits/NACP_Traits_Rcode/Results_SimpleMAreg_v12rawavgs_20171220_wCWM.csv") # added CIs for SMAs.



all.resultsold <- read.csv("Results_SimpleMAreg_v11rawavgs_20171211_wCWM.csv")





##### checking whether new results from 11.20 match old results from 08.28
  # takeaway: only numeric changes are in Genw.inFamily results. not sure why, but it doesn't seem big.

# length(which(as.numeric(all.results$sig_LL.N)==all.resultsold$sig_LL.N)) # 4/74
# length(which(as.numeric(all.results$sig_LL.Narea)==all.resultsold$sig_LL.Narea)) # 7/70
# length(which(as.numeric(all.results$sig_LMA.LL)==all.resultsold$sig_LMA.LL)) # 7/74
# length(which(as.numeric(all.results$sig_LMA.N)==all.resultsold$sig_LMA.N)) # 10/200
# length(which(as.numeric(all.results$sig_LMA.Narea)==all.resultsold$sig_LMA.Narea)) # 19/200
# 
# 
# length(which(as.numeric(all.results$Slope_LMA.LL)==all.resultsold$Slope_LMA.LL)) # 32/136
# plot(x=all.results$Slope_LMA.LL[which(as.numeric(all.results$Slope_LMA.LL)!=all.resultsold$Slope_LMA.LL)], y=all.resultsold$Slope_LMA.LL[which(as.numeric(all.results$Slope_LMA.LL)!=all.resultsold$Slope_LMA.LL)], pch=16, col=factor(all.results$Type[which(as.numeric(all.results$Slope_LMA.LL)!=all.resultsold$Slope_LMA.LL)]))
#   # for LMA v LL, all seems to be gen w/in fams that screwed up
# plot(x=all.results$Rho_LMA.LL[which(as.numeric(all.results$Rho_LMA.LL)!=all.resultsold$Rho_LMA.LL)], y=all.resultsold$Rho_LMA.LL[which(as.numeric(all.results$Slope_LMA.LL)!=all.resultsold$Slope_LMA.LL)], pch=16, col=factor(all.results$Type[which(as.numeric(all.results$Slope_LMA.LL)!=all.resultsold$Slope_LMA.LL)]))
# plot(x=all.results$Slope_LMA.N[which(as.numeric(all.results$Slope_LMA.N)!=all.resultsold$Slope_LMA.N)], y=all.resultsold$Slope_LMA.N[which(as.numeric(all.results$Slope_LMA.N)!=all.resultsold$Slope_LMA.N)], pch=16, col=factor(all.results$Type[which(as.numeric(all.results$Slope_LMA.N)!=all.resultsold$Slope_LMA.N)]))
# 
# 
# 
# xtabs(~all.results$Type[which(as.numeric(all.results$Slope_LMA.N)!=all.resultsold$Slope_LMA.N)])
# xtabs(~all.results$Type[which(as.numeric(all.results$sig_LMA.N)!=all.resultsold$sig_LMA.N)])
# 
# xtabs(~all.results$Type[which(as.numeric(all.results$Slope_LMA.Narea)!=all.resultsold$Slope_LMA.Narea)])
# xtabs(~all.results$Type[which(as.numeric(all.results$sig_LMA.Narea)!=all.resultsold$sig_LMA.Narea)])
# 
# xtabs(~all.results$Type[which(as.numeric(all.results$Slope_LMA.LL)!=all.resultsold$Slope_LMA.LL)])
# xtabs(~all.results$Type[which(as.numeric(all.results$sig_LMA.LL)!=all.resultsold$sig_LMA.LL)])
# 
# xtabs(~all.results$Type[which(as.numeric(all.results$Slope_LL.N)!=all.resultsold$Slope_LL.N)])
# xtabs(~all.results$Type[which(as.numeric(all.results$sig_LL.N)!=all.resultsold$sig_LL.N)])
# 
# plot(all.results$sig_LMA.Narea[which(as.numeric(all.results$sig_LMA.Narea)!=all.resultsold$sig_LMA.Narea)]~all.resultsold$sig_LMA.Narea[which(as.numeric(all.results$sig_LMA.Narea)!=all.resultsold$sig_LMA.Narea)])
# 
# cbind(type=all.results$Type[which(as.numeric(all.results$sig_LMA.Narea)!=all.resultsold$sig_LMA.Narea)],new=all.results$sig_LMA.Narea[which(as.numeric(all.results$sig_LMA.Narea)!=all.resultsold$sig_LMA.Narea)], old=all.resultsold$sig_LMA.Narea[which(as.numeric(all.results$sig_LMA.Narea)!=all.resultsold$sig_LMA.Narea)])
# cbind(type=all.results$Type[which(as.numeric(all.results$sig_LMA.LL)!=all.resultsold$sig_LMA.LL)],new=all.results$sig_LMA.LL[which(as.numeric(all.results$sig_LMA.LL)!=all.resultsold$sig_LMA.LL)], old=all.resultsold$sig_LMA.LL[which(as.numeric(all.results$sig_LMA.LL)!=all.resultsold$sig_LMA.LL)])
# cbind(type=all.results$Type[which(as.numeric(all.results$sig_LMA.N)!=all.resultsold$sig_LMA.N)],new=all.results$sig_LMA.N[which(as.numeric(all.results$sig_LMA.N)!=all.resultsold$sig_LMA.N)], old=all.resultsold$sig_LMA.N[which(as.numeric(all.results$sig_LMA.N)!=all.resultsold$sig_LMA.N)])
# cbind(type=all.results$Type[which(as.numeric(all.results$sig_LL.N)!=all.resultsold$sig_LL.N)],new=all.results$sig_LL.N[which(as.numeric(all.results$sig_LL.N)!=all.resultsold$sig_LL.N)], old=all.resultsold$sig_LL.N[which(as.numeric(all.results$sig_LL.N)!=all.resultsold$sig_LL.N)])
# cbind(type=all.results$Type[which(as.numeric(all.results$sig_LL.Narea)!=all.resultsold$sig_LL.Narea)],new=all.results$sig_LL.Narea[which(as.numeric(all.results$sig_LL.Narea)!=all.resultsold$sig_LL.Narea)], old=all.resultsold$sig_LL.Narea[which(as.numeric(all.results$sig_LL.Narea)!=all.resultsold$sig_LL.Narea)])



#________________________________________________________________________
############# Initial Plotting #################
#________________________________________________________________________


#________________________________________________________________________
########## LOAD RESULTS DATA ##################

# older version without community weighted mean row
# "Results_SimpleMAreg_v9rawavgs_20170620.csv" - version used in First Submission to Ecology Letters
# "Results_SimpleMAreg_v10rawavgs_20171120.csv" version used for revision, pre Mt Rainier addition
all.results <- read.csv("/Users/leeanderegg/Dropbox/NACP_Traits/NACP_Traits_Rcode/Results_SimpleMAreg_v11rawavgs_20171211.csv", row.names = 1)
all.results <- read.csv("/Users/leeanderegg/Desktop/scratch/20171220_results/SMA_Results_All.csv", row.names=1)
levels(all.results$Type) <- list(w.inSpp = "w.inSpp", w.inGen = "w.inGen", Sppw.inFam= "Sppw.inFam",Genw.inFam="Genw.inFam", Fam="Fam",Famclean="Famclean", global="global")
# # add in the Mt Rainier Species
# empty <-  data.frame(matrix(NA, nrow = 6, ncol=ncol(all.results.cwm)))
# names(empty) <- names(all.results.cwm)
# empty$Taxo.Unit <- c("Aster alpigenus","Castilleja parviflora", "Erythronium montanum", "Lupinus arcticus","Valeriana sitchensis","Veratrum viride")
# empty$Type <- "w.inSpp"
# all.resultstest <- rbind(all.results[which(all.results$Type=="w.inSpp"),], empty, all.results[which(all.results$Type!="w.inSpp"),])
# all.resultstest[which(all.resultstest$Taxo.Unit %in% c("Aster alpigenus","Castilleja parviflora", "Erythronium montanum", "Lupinus arcticus","Valeriana sitchensis","Veratrum viride")), grep("LMA.LL", colnames(all.resultstest))] <- Rain.results.LMALL[,which(!colnames(Rain.results.LMALL) %in% c("Species","varLMA","varLL","rho.sig"))]
# all.resultstest[which(all.resultstest$Taxo.Unit %in% c("Aster alpigenus","Castilleja parviflora", "Erythronium montanum", "Lupinus arcticus","Valeriana sitchensis","Veratrum viride")), c("varLMA","varLL")] <- Rain.results.LMALL[,c("varLMA","varLL")]
# all.resultstest[which(all.resultstest$Taxo.Unit %in% c("Aster alpigenus","Castilleja parviflora", "Erythronium montanum", "Lupinus arcticus","Valeriana sitchensis","Veratrum viride")), grep("LMA.N", colnames(all.resultstest))[1:14]] <- Rain.results.LMAN[,which(!colnames(Rain.results.LMAN) %in% c("Species","varLMA","varNmass","rho.sig"))]
# all.resultstest[which(all.resultstest$Taxo.Unit %in% c("Aster alpigenus","Castilleja parviflora", "Erythronium montanum", "Lupinus arcticus","Valeriana sitchensis","Veratrum viride")), c("varNmass")] <- Rain.results.LMAN[,c("varNmass")]
# all.resultstest[which(all.resultstest$Taxo.Unit %in% c("Aster alpigenus","Castilleja parviflora", "Erythronium montanum", "Lupinus arcticus","Valeriana sitchensis","Veratrum viride")), grep("LMA.Narea", colnames(all.resultstest))] <- Rain.results.LMANarea[,which(!colnames(Rain.results.LMANarea) %in% c("Species","varLMA","varNarea","rho.sig"))]
# all.resultstest[which(all.resultstest$Taxo.Unit %in% c("Aster alpigenus","Castilleja parviflora", "Erythronium montanum", "Lupinus arcticus","Valeriana sitchensis","Veratrum viride")), c("varNarea")] <- Rain.results.LMANarea[,c("varNarea")]
# all.resultstest[which(all.resultstest$Taxo.Unit %in% c("Aster alpigenus","Castilleja parviflora", "Erythronium montanum", "Lupinus arcticus","Valeriana sitchensis","Veratrum viride")), grep("LL.N", colnames(all.resultstest))[1:14]] <- Rain.results.LLN[,which(!colnames(Rain.results.LLN) %in% c("Species","varLL","varNmass","rho.sig"))]
# all.resultstest[which(all.resultstest$Taxo.Unit %in% c("Aster alpigenus","Castilleja parviflora", "Erythronium montanum", "Lupinus arcticus","Valeriana sitchensis","Veratrum viride")), grep("LL.Narea", colnames(all.resultstest))] <- Rain.results.LLNarea[,which(!colnames(Rain.results.LLNarea) %in% c("Species","varLL","varNarea","rho.sig"))]




all.results.cl <- all.results %>% filter(Type %in% c("w.inSpp","w.inGen","Genw.inFam","Famclean","global"))
all.results <- read.csv("/Users/leeanderegg/Desktop/scratch/20171220_results/SMA_Results_All.csv", row.names=1)
all.results.cl$Type <- factor(all.results.cl$Type)
all.results.cl$Taxo.Unit <- factor(all.results.cl$Taxo.Unit)
# add columns of whether the correlation is significant (useful for plotting)
all.results.cl$rho.sig_LMA.N <- rho.sig(all.results.cl$Rho_LMA.N, all.results.cl$n_LMA.N)
all.results.cl$rho.sig_LMA.Narea <- rho.sig(all.results.cl$Rho_LMA.Narea, all.results.cl$n_LMA.Narea)
all.results.cl$rho.sig_LMA.LL <- rho.sig(all.results.cl$Rho_LMA.LL, all.results.cl$n_LMA.LL)
all.results.cl$rho.sig_LL.N <- rho.sig(all.results.cl$Rho_LL.N, all.results.cl$n_LL.N)
all.results.cl$rho.sig_LL.Narea <- rho.sig(all.results.cl$Rho_LL.Narea, all.results.cl$n_LL.Narea)

# "Results_SimpleMAreg_v9rawavgs_20170620.csv" - version used in First Submission to Ecology Letters
# "Results_SimpleMAreg_v10rawavgs_20171120_wCWM.csv"
all.results.cwm <- read.csv("/Users/leeanderegg/Dropbox/NACP_Traits/NACP_Traits_Rcode/Results_SimpleMAreg_v11rawavgs_20171211.csv", row.names = 1)
all.results.cwm <- read.csv("/Users/leeanderegg/Desktop/scratch/20171220_results/SMA_Results_All_CWM.csv", row.names=1)
levels(all.results.cwm$Type) <- list(w.inSpp = "w.inSpp", w.inGen = "w.inGen", Sppw.inFam= "Sppw.inFam",Genw.inFam="Genw.inFam", Fam="Fam",Famclean="Famclean", global="global", CWM="CWM")
#   # add in the Mt Rainier Species
# empty <-  data.frame(matrix(NA, nrow = 6, ncol=ncol(all.results.cwm)))
# names(empty) <- names(all.results.cwm)
# empty$Taxo.Unit <- c("Aster alpigenus","Castilleja parviflora", "Erythronium montanum", "Lupinus arcticus","Valeriana sitchensis","Veratrum viride")
# empty$Type <- "w.inSpp"
# all.results.cwmtest <- rbind(all.results.cwm[which(all.results.cwm$Type=="w.inSpp"),], empty, all.results.cwm[which(all.results.cwm$Type!="w.inSpp"),])
# all.results.cwmtest[which(all.results.cwmtest$Taxo.Unit %in% c("Aster alpigenus","Castilleja parviflora", "Erythronium montanum", "Lupinus arcticus","Valeriana sitchensis","Veratrum viride")), grep("LMA.LL", colnames(all.results.cwmtest))] <- spp.results[,which(!colnames(spp.results) %in% c("Species","varLMA","varLL","rho.sig"))]
# all.results.cwmtest[which(all.results.cwmtest$Taxo.Unit %in% c("Aster alpigenus","Castilleja parviflora", "Erythronium montanum", "Lupinus arcticus","Valeriana sitchensis","Veratrum viride")), c("varLMA","varLL")] <- spp.results[,c("varLMA","varLL")]


# filter out the levels I don't care about
all.results.cwm.cl <- all.results.cwm %>% filter(Type %in% c("w.inSpp","w.inGen","Genw.inFam","Famclean","global", "CWM"))
all.results.cwm.cl$Type <- factor(all.results.cwm.cl$Type)
all.results.cwm.cl$Taxo.Unit <- factor(all.results.cwm.cl$Taxo.Unit)
  # add columns of whether the correlation is significant (useful for plotting)
all.results.cwm.cl$rho.sig_LMA.N <- rho.sig(all.results.cwm.cl$Rho_LMA.N, all.results.cwm.cl$n_LMA.N)
all.results.cwm.cl$rho.sig_LMA.Narea <- rho.sig(all.results.cwm.cl$Rho_LMA.Narea, all.results.cwm.cl$n_LMA.Narea)
all.results.cwm.cl$rho.sig_LMA.LL <- rho.sig(all.results.cwm.cl$Rho_LMA.LL, all.results.cwm.cl$n_LMA.LL)
all.results.cwm.cl$rho.sig_LL.N <- rho.sig(all.results.cwm.cl$Rho_LL.N, all.results.cwm.cl$n_LL.N)
all.results.cwm.cl$rho.sig_LL.Narea <- rho.sig(all.results.cwm.cl$Rho_LL.Narea, all.results.cwm.cl$n_LL.Narea)




#________________________________________________________________________
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
##### Boxplots ######
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
boxplot(Slope_LL.N~Type, all.results, ylim=c(-3,2.5),las=3, main="log(Nmass)~log(LL)", ylab="MA Slope")
points(y=rep(as.numeric(fam.res_LL.N[3]), times=7),x=seq(4.7,5.3, by=.1), pch=15, cex=.7)
points(y=rep(as.numeric(fam.resclean_LL.N[3]), times=7),x=seq(5.7,6.3, by=.1), pch=15, cex=.7)
abline(h=0, lty=2)
boxplot(Rho_LL.N~Type, all.results, ylim=c(-1,1),las=3, ylab="Rho")
points(y=rep(as.numeric(fam.res_LL.N[4]), times=7),x=seq(4.7,5.3, by=.1), pch=15, cex=.7)
points(y=rep(as.numeric(fam.resclean_LL.N[4]), times=7),x=seq(5.7,6.3, by=.1), pch=15, cex=.7)
abline(h=0, lty=2)

# quartz(width=3, height=4.5) # not updated for axis flop
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
p <- boxplot(Slope_LL.N~Type, all.results[which(all.results$n_LL.N>5),]
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

p <- boxplot(Rho_LL.N~Type, all.results[which(all.results$n_LL.N>5),]
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
boxplot(Slope_LL.N~Type, all.results[which(abs(all.results$Rho_LL.N)>.3),], ylim=c(-3,2.5),las=3, main="log(LL)~log(Nmass) Hcor", ylab="MA Slope")
points(y=rep(as.numeric(fam.res_LL.N[3]), times=7),x=seq(4.7,5.3, by=.1), pch=15, cex=.7)
points(y=rep(as.numeric(fam.resclean_LL.N[3]), times=7),x=seq(5.7,6.3, by=.1), pch=15, cex=.7)
abline(h=0, lty=2)
boxplot(Rho_LL.N~Type, all.results[which(all.results$n_LL.N>9),], ylim=c(-1,1),las=3, ylab="Rho")
points(y=rep(as.numeric(fam.res_LL.N[4]), times=7),x=seq(4.7,5.3, by=.1), pch=15, cex=.7)
points(y=rep(as.numeric(fam.resclean_LL.N[4]), times=7),x=seq(5.7,6.3, by=.1), pch=15, cex=.7)
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





#________________________________________________________________________
####### **STATISTICAL ANALYSES** ##########
#________________________________________________________________________

###### Analyzing the correlations
LMALL <- lm(Rho_LMA.LL~Type, weights = n_LMA.LL, all.results)
# so Spp are significantly negative
# and every level above spp is positive and significantly different from w/in.spp
LMAN <- lm(Rho_LMA.N~Type, weights = n_LMA.N, all.results)
# they're all negative, but becoming increasingly negative at higher Taxo levels
lman <- aov(Rho_LMA.N~Type, all.results,weights=n_LMA.N)
TukeyHSD(lman)
# everything different from w/in spp. but above spp not different


NmassLL <- lm(Rho_LL.N~Type, weights = n_LL.N, all.results)
# w.inGen not different from w/.in spp, but everything else is...
nmassll <- aov(Rho_LL.N~Type, all.results,weights=n_LL.N)
TukeyHSD(nmassll)
# with super conservative TukeyHSD, nothing is significantly different, though w/spp is close



###### Analyzing the MA slopes
LMALL <- lm(Slope_LMA.LL~Type, weights = n_LMA.LL, all.results)
# so Spp are significantly negative
# and every level above spp is positive and significantly different from w/in.spp



#### SLOPE LMA v LL
  # one outlier (129), but it doesn't seem to really screw the model up too badly. And it goes away with different weightings
# originally misspecified as call all.results but filtering with all.results.cl. only effected family values, but changed some significance.
# on 11.20.2017 updated for Ecol Let manuscript revision.
# results are updated in Table S3, in SI version v8


LMALL <- lm(Slope_LMA.LL~Type, all.results.cl[which(all.results.cl$n_LMA.LL>5),])
LMALLnull <- lm(Slope_LMA.LL~1, all.results.cl[which(all.results.cl$n_LMA.LL>5),])
anova(LMALL, LMALLnull) #newest p=6.511e-05 new average p=0.001475, old p = 0.018 old average of logs:p=0.044
LMALL <- lm(Slope_LMA.LL~Type, weights = n_LMA.LL, all.results.cl[which(all.results.cl$n_LMA.LL>5),])
LMALLnull <- lm(Slope_LMA.LL~1, weights = n_LMA.LL, all.results.cl[which(all.results.cl$n_LMA.LL>5),])
anova(LMALL, LMALLnull) # new p=0.0001136, old p=0.025, old p<0.0001 old
LMALL <- lm(Slope_LMA.LL~Type, weights = varLMA, all.results.cl[which(all.results.cl$n_LMA.LL>5),])
LMALLnull <- lm(Slope_LMA.LL~1, weights = varLMA, all.results.cl[which(all.results.cl$n_LMA.LL>5),])
anova(LMALL, LMALLnull) #new: p= 3.1 e -7, old p=0.037, old p=0.031
LMALL <- lm(Slope_LMA.LL~Type, weights = varLL, all.results.cl[which(all.results.cl$n_LMA.LL>5),])
LMALLnull <- lm(Slope_LMA.LL~1, weights = varLL, all.results.cl[which(all.results.cl$n_LMA.LL>5),])
anova(LMALL, LMALLnull) # new p=1.02e-8, old p=0.00077 old, p=0.002
# Rho LMA v LL
LMALL <- lm(Rho_LMA.LL~Type, all.results.cl[which(all.results.cl$n_LMA.LL>5),])
LMALLnull <- lm(Rho_LMA.LL~1, all.results.cl[which(all.results.cl$n_LMA.LL>5),])
anova(LMALL, LMALLnull) # new & old: p<0.0001, p=0.007
LMALL <- lm(Rho_LMA.LL~Type, weights = n_LMA.LL, all.results.cl[which(all.results.cl$n_LMA.LL>5),])
LMALLnull <- lm(Rho_LMA.LL~1, weights = n_LMA.LL, all.results.cl[which(all.results.cl$n_LMA.LL>5),])
anova(LMALL, LMALLnull) # p<0.0001 new and old and older
LMALL <- lm(Rho_LMA.LL~Type, weights = varLMA, all.results.cl[which(all.results.cl$n_LMA.LL>5),])
LMALLnull <- lm(Rho_LMA.LL~1, weights = varLMA, all.results.cl[which(all.results.cl$n_LMA.LL>5),])
anova(LMALL, LMALLnull) #new & old p<0.0001, old p=0.0019
LMALL <- lm(Rho_LMA.LL~Type, weights = varLL, all.results.cl[which(all.results.cl$n_LMA.LL>5),])
LMALLnull <- lm(Rho_LMA.LL~1, weights = varLL, all.results.cl[which(all.results.cl$n_LMA.LL>5),])
anova(LMALL, LMALLnull) # p<0.0001 new and old and older

#### SLOPE LMA v Nmass
  # Not the world's most normal
  # The weighting scheme introduces pretty large parameter instability...
# originally misspecified as call all.results but filtering with all.results.cl. only effected family values, but changed some significance.
# on 11.20.2017 updated for Ecol Let manuscript revision.
  # results are updated in Table S3, in SI version v8
LMAN <- lm(Slope_LMA.N~Type, all.results.cl[which(all.results.cl$n_LMA.N>5),])
LMANnull <- lm(Slope_LMA.N~1, all.results.cl[which(all.results.cl$n_LMA.N>5),])
anova(LMAN, LMANnull) #new p=0.3292, old raw averaged: p=0.507 p=0.736
LMAN <- lm(Slope_LMA.N~Type, weights = n_LMA.N, all.results.cl[which(all.results.cl$n_LMA.N>5),])
LMANnull <- lm(Slope_LMA.N~1, weights = n_LMA.N, all.results.cl[which(all.results.cl$n_LMA.N>5),])
anova(LMAN, LMANnull) # new p=0.8726, old p=0.208 p=0.96
LMAN <- lm(Slope_LMA.N~Type, weights = varLMA, all.results.cl[which(all.results.cl$n_LMA.N>5),])
LMANnull <- lm(Slope_LMA.N~1, weights = varLMA, all.results.cl[which(all.results.cl$n_LMA.N>5),])
anova(LMAN, LMANnull) #new p=0.07163, old p=0.274, old p=0.437
LMAN <- lm(Slope_LMA.N~Type, weights = varNmass, all.results.cl[which(all.results.cl$n_LMA.N>5),])
LMANnull <- lm(Slope_LMA.N~1, weights = varNmass, all.results.cl[which(all.results.cl$n_LMA.N>5),])
anova(LMAN, LMANnull) #new p=0.8796, old p=0.534, old p=0.899
# Rho LMA v Nmass
LMAN <- lm(Rho_LMA.N~Type, all.results.cl[which(all.results.cl$n_LMA.N>5),])
LMANnull <- lm(Rho_LMA.N~1, all.results.cl[which(all.results.cl$n_LMA.N>5),])
anova(LMAN, LMANnull) #new p=0.000363, old raw avgs p=0.007 , old p=0.0273
LMAN <- lm(Rho_LMA.N~Type, weights = n_LMA.N, all.results.cl[which(all.results.cl$n_LMA.N>5),])
LMANnull <- lm(Rho_LMA.N~1, weights = n_LMA.N, all.results.cl[which(all.results.cl$n_LMA.N>5),])
anova(LMAN, LMANnull) # p<0.001 new and old and older
LMAN <- lm(Rho_LMA.N~Type, weights = varLMA, all.results.cl[which(all.results.cl$n_LMA.N>5),])
LMANnull <- lm(Rho_LMA.N~1, weights = varLMA, all.results.cl[which(all.results.cl$n_LMA.N>5),])
anova(LMAN, LMANnull) #new p<0.001,  p=0.003 old
LMAN <- lm(Rho_LMA.N~Type, weights = varNmass, all.results.cl[which(all.results.cl$n_LMA.N>5),])
LMANnull <- lm(Rho_LMA.N~1, weights = varNmass, all.results.cl[which(all.results.cl$n_LMA.N>5),])
anova(LMAN, LMANnull) # new p=0.02266, old p=0.111, p=0.186


#### SLOPE LL v Nmass
  # looks ok, and parameter uncertainty instability reasonably limited
  # originally misspecified as call all.results but filtering with all.results.cl. only effected family values, but changed some significance.
  # on 11.20.2017 updated for Ecol Let manuscript revision.
LLNmass <- lm(Slope_LL.N~Type, all.results.cl[which(all.results.cl$n_LL.N>5),])
LLNmassnull <- lm(Slope_LL.N~1, all.results.cl[which(all.results.cl$n_LL.N>5),])
anova(LLNmass, LLNmassnull) #with all.results.cl: p=0.2012, old raw avg (w/ all.results but filtered with all.results.cl) p=0.1729, old p=0.897
LLNmass <- lm(Slope_LL.N~Type, weights = n_LL.N, all.results.cl[which(all.results.cl$n_LL.N>5),])
LLNmassnull <- lm(Slope_LL.N~1, weights=n_LL.N, all.results.cl[which(all.results.cl$n_LL.N>5),])
anova(LLNmass, LLNmassnull) #new: p=0.4924, old p=0.687, old p=0.0878
LLNmass <- lm(Slope_LL.N~Type, weights=varLL, all.results.cl[which(all.results.cl$n_LL.N>5),])
LLNmassnull <- lm(Slope_LL.N~1, all.results.cl[which(all.results.cl$n_LL.N>5),], weights=varLL)
anova(LLNmass, LLNmassnull) #new: p=0.1455, old p=0.1953, old p=0.31
LLNmass <- lm(Slope_LL.N~Type, weights=varNmass, all.results.cl[which(all.results.cl$n_LL.N>5),])
LLNmassnull <- lm(Slope_LL.N~1, all.results.cl[which(all.results.cl$n_LL.N>5),], weights=varNmass)
anova(LLNmass, LLNmassnull) #new: p=0.08603,  old: p=0.0051, old p=0.61
#### Rho LL v Nmass
LLNmass <- lm(Rho_LL.N~Type, all.results.cl[which(all.results.cl$n_LL.N>5),])
LLNmassnull <- lm(Rho_LL.N~1, all.results.cl[which(all.results.cl$n_LL.N>5),])
anova(LLNmass, LLNmassnull) #new p=0.001934, with misspecification p=0.066, p=0.018
LLNmass <- lm(Rho_LL.N~Type, weights = n_LL.N, all.results.cl[which(all.results.cl$n_LL.N>5),])
LLNmassnull <- lm(Rho_LL.N~1, weights=n_LL.N, all.results.cl[which(all.results.cl$n_LL.N>5),])
anova(LLNmass, LLNmassnull) #new: p=2.186e-08, with misspecification p=0.243, new and old 
LLNmass <- lm(Rho_LL.N~Type, weights=varLL, all.results.cl[which(all.results.cl$n_LL.N>5),])
LLNmassnull <- lm(Rho_LL.N~1, all.results.cl[which(all.results.cl$n_LL.N>5),], weights=varLL)
anova(LLNmass, LLNmassnull) #new:0.0005164, misspecified: p=0.0379 p<0.0001 old 
LLNmass <- lm(Rho_LL.N~Type, weights=varNmass, all.results.cl[which(all.results.cl$n_LL.N>5),])
LLNmassnull <- lm(Rho_LL.N~1, all.results.cl[which(all.results.cl$n_LL.N>5),], weights=varNmass)
anova(LLNmass, LLNmassnull) #new p=0.001932, misspecified p=0.017, old p=0.005



#### SLOPE LMA v Narea
  # NOTE: Protea repens has a very strange and strong negative relationship, and Abies Alba really screws over the genus Abies, giving it a strong negative trend despite the rest being postive.
  # I have to remove these outliers to keep from violating assumptions of normaility
# updated, already had appropriate all.results.cl, but changed slightly with new version 11.20.17
LMANarea <- lm(Slope_LMA.Narea~Type, all.results.cl[which(all.results.cl$n_LMA.Narea>5 & all.results.cl$Slope_LMA.Narea>-1),])
LMANareanull <- lm(Slope_LMA.Narea~1, all.results.cl[which(all.results.cl$n_LMA.Narea>5 & all.results.cl$Slope_LMA.Narea>-1),])
anova(LMANarea, LMANareanull) #new p=0.0006761, raw avg p=0.0.0009, old: p=0.0632, cl=0.02992
LMANarea <- lm(Slope_LMA.Narea~Type, weights = n_LMA.Narea, all.results.cl[which(all.results.cl$n_LMA.Narea>5 & all.results.cl$Slope_LMA.Narea>-1),])
LMANareanull <- lm(Slope_LMA.Narea~1, weights=n_LMA.Narea, all.results.cl[which(all.results.cl$n_LMA.Narea>5 & all.results.cl$Slope_LMA.Narea>-1),])
anova(LMANarea, LMANareanull) # p<0.0001 &cl new and old  and older
LMANarea <- lm(Slope_LMA.Narea~Type, weights=varLMA, all.results.cl[which(all.results.cl$n_LMA.Narea>5 & all.results.cl$Slope_LMA.Narea>-1),])
LMANareanull <- lm(Slope_LMA.Narea~1, all.results.cl[which(all.results.cl$n_LMA.Narea>5 & all.results.cl$Slope_LMA.Narea>-1),], weights=varLMA)
anova(LMANarea, LMANareanull) #new p=0.003378,  raw avg p=0.00319 p=0.135, cl=0.001735
LMANarea <- lm(Slope_LMA.Narea~Type, weights=varNarea, all.results.cl[which(all.results.cl$n_LMA.Narea>5 & all.results.cl$Slope_LMA.Narea>-1),])
LMANareanull <- lm(Slope_LMA.Narea~1, all.results.cl[which(all.results.cl$n_LMA.Narea>5 & all.results.cl$Slope_LMA.Narea>-1),], weights=varNarea)
anova(LMANarea, LMANareanull) #new p=0.03914, raw avg p=0.05068, old p=0.549, cl=0.4077
#### Rho LMA v Narea
LMANarea <- lm(Rho_LMA.Narea~Type, all.results.cl[which(all.results.cl$n_LMA.Narea>5 & all.results.cl$Slope_LMA.Narea>-1),])
LMANareanull <- lm(Rho_LMA.Narea~1, all.results.cl[which(all.results.cl$n_LMA.Narea>5 & all.results.cl$Slope_LMA.Narea>-1),])
anova(LMANarea, LMANareanull) #new p=0.624, old raw avg p=0.5843, old log avgs p=0.98, cl=0.8787
LMANarea <- lm(Rho_LMA.Narea~Type, weights = n_LMA.Narea, all.results.cl[which(all.results.cl$n_LMA.Narea>5 & all.results.cl$Slope_LMA.Narea>-1),])
LMANareanull <- lm(Rho_LMA.Narea~1, weights=n_LMA.Narea, all.results.cl[which(all.results.cl$n_LMA.Narea>5 & all.results.cl$Slope_LMA.Narea>-1),])
anova(LMANarea, LMANareanull) #new p=0.0039, old raw avg p=0.0036, old log avg p<0.45, cl=0.2666
LMANarea <- lm(Rho_LMA.Narea~Type, weights=varLMA, all.results.cl[which(all.results.cl$n_LMA.Narea>5 & all.results.cl$Slope_LMA.Narea>-1),])
LMANareanull <- lm(Rho_LMA.Narea~1, all.results.cl[which(all.results.cl$n_LMA.Narea>5 & all.results.cl$Slope_LMA.Narea>-1),], weights=varLMA)
anova(LMANarea, LMANareanull) #new p=0.3878, p=0.3736, old p<0.20
LMANarea <- lm(Rho_LMA.Narea~Type, weights=varNarea, all.results.cl[which(all.results.cl$n_LMA.Narea>5 & all.results.cl$Slope_LMA.Narea>-1),])
LMANareanull <- lm(Rho_LMA.Narea~1, all.results.cl[which(all.results.cl$n_LMA.Narea>5 & all.results.cl$Slope_LMA.Narea>-1),], weights=varNarea)
anova(LMANarea, LMANareanull) #new p=0.2859, p=0.209, old  p=0.609


#### SLOPE LL v Narea
  # a little parameter instability, and a bit non-normal in some formulations. But not bad.,
# was using .cl, p values updated slightly with v 11.20.17
LLNarea <- lm(Slope_LL.Narea~Type, all.results.cl[which(all.results.cl$n_LL.Narea>5),])
LLNareanull <- lm(Slope_LL.Narea~1, all.results.cl[which(all.results.cl$n_LL.Narea>5),])
anova(LLNarea, LLNareanull) #new p=0.02021, old raw avg p=0.01799, old  p=0.2545, cl=0.1181
LLNarea <- lm(Slope_LL.Narea~Type, weights = n_LL.Narea, all.results.cl[which(all.results.cl$n_LL.Narea>5),])
LLNareanull <- lm(Slope_LL.Narea~1, weights=n_LL.Narea, all.results.cl[which(all.results.cl$n_LL.Narea>5),])
anova(LLNarea, LLNareanull) # p=0.0007899, old cl p= 1.4 e-6 new and old
LLNarea <- lm(Slope_LL.Narea~Type, weights=varLL, all.results.cl[which(all.results.cl$n_LL.Narea>5),])
LLNareanull <- lm(Slope_LL.Narea~1, all.results.cl[which(all.results.cl$n_LL.Narea>5),], weights=varLL)
anova(LLNarea, LLNareanull) #new p= 0.0008546, old p=0.00074 p=, cl=0.00364
LLNarea <- lm(Slope_LL.Narea~Type, weights=varNarea, all.results.cl[which(all.results.cl$n_LL.Narea>5),])
LLNareanull <- lm(Slope_LL.Narea~1, all.results.cl[which(all.results.cl$n_LL.Narea>5),], weights=varNarea)
anova(LLNarea, LLNareanull) #new p=0.00489, p=0.0041,  old, cl=0.00805 
#### Rho LL v Narea
LLNarea <- lm(Rho_LL.Narea~Type, all.results.cl[which(all.results.cl$n_LL.Narea>5),])
LLNareanull <- lm(Rho_LL.Narea~1, all.results.cl[which(all.results.cl$n_LL.Narea>5),])
anova(LLNarea, LLNareanull) #new & old raw avg p=0.004, old log avg p=, cl=0.22
LLNarea <- lm(Rho_LL.Narea~Type, weights = n_LL.Narea, all.results.cl[which(all.results.cl$n_LL.Narea>5),])
LLNareanull <- lm(Rho_LL.Narea~1, weights=n_LL.Narea, all.results.cl[which(all.results.cl$n_LL.Narea>5),])
anova(LLNarea, LLNareanull) # new p= 3.938e-6, pcl = 1.42 e-8 new and old
LLNarea <- lm(Rho_LL.Narea~Type, weights=varLL, all.results.cl[which(all.results.cl$n_LL.Narea>5),])
LLNareanull <- lm(Rho_LL.Narea~1, all.results.cl[which(all.results.cl$n_LL.Narea>5),], weights=varLL)
anova(LLNarea, LLNareanull) #new= p=0.00056, p=0.00062, old p<0.03126
LLNarea <- lm(Rho_LL.Narea~Type, weights=varNarea, all.results.cl[which(all.results.cl$n_LL.Narea>5),])
LLNareanull <- lm(Rho_LL.Narea~1, all.results.cl[which(all.results.cl$n_LL.Narea>5),], weights=varNarea)
anova(LLNarea, LLNareanull) #new p=0.003433, p=0.0037, old p=0.07745





# w.inGen not different from w/.in spp, but everything else is...
nmassll <- aov(Slope_LL.N~Type, all.results,weights=n_LL.N)
TukeyHSD(nmassll)
# with super conservative TukeyHSD, nothing is significantly different, though w/spp is close

LMA.Narea <- lm(Slope_LMA.Narea~Type, weights = n_LMA.Narea, all.results)
# Everything's different from w.inSpp. and increasingly different at higher Taxon scales...
lmanarea <- aov(Slope_LMA.Narea~Type, all.results)
TukeyHSD(lmanarea)
lmanarea <- aov(Slope_LMA.Narea~Type, all.results,weights=n_LMA.Narea)
TukeyHSD(lmanarea)
lmanarea <- aov(Slope_LMA.Narea~Type, all.results,weights=varLMA)
TukeyHSD(lmanarea)
lmanarea <- aov(Slope_LMA.Narea~Type, all.results,weights=varNarea)
TukeyHSD(lmanarea)










#------------------------------------------------------------------------
##########  Re-analysis only with PNW stands with LAI <4 to seperate the sun/shade issue #####
#------------------------------------------------------------------------

traits.common.sun <- traits.common5 %>% filter(LAI_O < 4)





#traits.common.sun <- traits.common.s[which(traits.common$SP.




####### ***LMA and LL*** #####################

############ .Species level analysis #####################


spp.results <- data.frame(matrix(NA, nrow=length(unique(traits.common.sun$SP.ID)), ncol=15))
colnames(spp.results) <- c("Species", "Int","Slope","Rho","r.sq","n","varLMA","varLL","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(traits.common.sun$SP.ID))){
  species <- levels(traits.common.sun$SP.ID)[i]
  print(species)
  dataz <- traits.common.sun[which(traits.common.sun$SP.ID==species),]
  res <- fit.MAR(xvar='log.LMA',yvar="log.LL",data=dataz)
  spp.results[i,1] <- species
  spp.results[i,2:8] <- res
  if (!is.na(res[1]) &res[5]>4){ # only fit null model if there are >5 data points
    nullbounds <- fit.null(xvar='log.LMA', yvar="log.LL", observed = dataz, nulldata = allspp, nits = 1000)  
    spp.results[i, 9:14] <- nullbounds
    spp.results[i,15] <- test.sig(x=spp.results$Rho[i], test=nullbounds)
  }
}



  


#------------------------------------------------------------------------
##########  Re-analysis only with PNW samples with LAI <median LL to seperate the sun/shade issue #####
#------------------------------------------------------------------------

traits.common.old <- traits.common5 %>% group_by(SP.ID) %>% filter(LLmonths > median(LLmonths, na.rm=T))

traits.common.old <- traits.common5 %>% group_by(SP.ID) %>% filter(LLmonths > quantile(LLmonths,.25, na.rm=T))




#traits.common.sun <- traits.common.s[which(traits.common$SP.




####### ***LMA and LL*** #####################

############ .Species level analysis #####################


spp.results <- data.frame(matrix(NA, nrow=length(unique(traits.common.old$SP.ID)), ncol=15))
colnames(spp.results) <- c("Species", "Int","Slope","Rho","r.sq","n","varLMA","varLL","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(traits.common.old$SP.ID))){
  species <- levels(traits.common.old$SP.ID)[i]
  print(species)
  dataz <- traits.common.old[which(traits.common.old$SP.ID==species),]
  res <- fit.MAR(xvar='log.LMA',yvar="log.LL",data=dataz)
  spp.results.old[i,1] <- species
  spp.results.old[i,2:8] <- res
  if (!is.na(res[1]) &res[5]>4){ # only fit null model if there are >5 data points
    nullbounds <- fit.null(xvar='log.LMA', yvar="log.LL", observed = dataz, nulldata = allspp, nits = 1000)  
    spp.results.old[i, 9:14] <- nullbounds
    spp.results.old[i,15] <- test.sig(x=spp.results.old$Rho[i], test=nullbounds)
  }
}


spp.results.oldmed <- spp.results.old  

                                           