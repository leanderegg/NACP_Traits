#############################################################
###      Analysis of PNW and GLOPNET trait datasets to 
###         assess taxonomic scales of trait variation
###         data downloaded from: http://dx.doi.org/10.3334/ORNLDAAC/1292 
###            on 01/23/16 by LDLA
#############################################################





#________________________________________________________________________
############ **Data preparation** ###############################################################
#________________________________________________________________________

# Requires: LES dataframe (GLOPNET)
#           traits.common (PNW with >10 records)
#           traits.common5 (PNW with 5 or more records)

levels(LES$Needle.Broad.lf) <- list(B="B",N="N",unknown="")

# ### common spp and genera and families in GLOPNE. only used in plotting subsetting
# commonspp <- names(which(xtabs(~Species, LES)>=5)) # 6 have >5 records, 14 have >4 records, 36 have >3 records
# commongen <- names(which(xtabs(~Genus, LES)>=5))
# commonfam <- names(which(xtabs(~Family, LES)>=5))
# 


####### Making complete INTRA-SPECIFIC dataset ##################

# # PNW species with 5 or more records
# spp.data1 <- traits.common5 %>% select(FullSpecies,log.Nmass, log.LL, log.LMA_PSA, GENUS, Family, log.Narea)
# colnames(spp.data1)[c(1,4,5)]<- c("Species","log.LMA", "Genus")
#   # NOTE: there are 95 trait measurements where LL is set at 1 year (log.LL = 1.079181). Luckily, I think my variance test in fit.MAR weeds these out...
#   ## NOTE: I orginally was using log.LMA_PSA, but as of 3.14.17, I had some NAs in the traits dataset. I had to fill them in.
# 
# ## LES species with 5 or more records
# commonspp <- names(which(xtabs(~Species, LES)>=5)) # 9 have >5 records, 16 have >4 records, 39 have >3 records
# spp.data2 <- LES %>% filter(Species %in% commonspp) %>% select(Species, log.Nmass, log.LL, log.LMA, Genus, Family, log.Narea)
# 
# data.supp <- read.csv("/Users/leeanderegg/Dropbox/NACP_Traits/Intra-data/OtherData_Combined_033017.csv", header=T, row.names=1)
# spp.data3 <- data.supp %>% select(Species,log.Nmass,log.LL, log.LMA, Genus, Family,log.Narea)
# 
# ## created combined dataset with 29 total species that have >5 observations
# spp.data <- rbind(spp.data1, spp.data2, spp.data3)
# spp.data$Species <- factor(spp.data$Species)
# ## look to make sure the combination worked 
# #ggplot(spp.data, aes(x=log.LMA, y=log.Nmass, col=Species)) + geom_point() + geom_smooth(method="lm", se=F)

## as of 04.01.17 working with the full dataset already combined in Cd-Initial_Analysis for the Variance Decomp
commonspp <-  names(which(xtabs(~Species, data.all)>=5))
spp.data<- data.all %>% filter(Species %in% commonspp)
spp.data$Species <- factor(spp.data$Species)
spp.data$Genus <- factor(spp.data$Genus)
spp.data$Family <- factor(spp.data$Family)
  # this actually keeps 26 more records than the old version, 44 spp, 26 genera, and 15 families




########## Species level trait averages ########

allspp <- data.all %>% group_by(Species, Genus, Family) %>% summarise( slog.LL = mean(log.LL, na.rm=T), slog.LMA = mean(log.LMA, na.rm=T), slog.Nmass = mean(log.Nmass, na.rm=T),slog.Narea = mean(log.Narea, na.rm=T)
                                                                       ,rslog.LL = log(mean(10^log.LL, na.rm=T),base=10), rslog.LMA = log(mean(10^log.LMA, na.rm=T),base=10), rslog.Nmass = log(mean(10^log.Nmass, na.rm=T),base=10)
                                                                       ,MAT = mean(MAT, na.rm=T), MAP=mean(MAP, na.rm=T), VPD=mean(VPD, na.rm=T) )
colnames(allspp) <- gsub("slog", "log", colnames(allspp))
  # this results in 1993 entries, 67 fewer entries than old (20 that may have been replicates and the 47 w/out family)
## LES species means
# LESspp <- LES %>% group_by(Species, GE.SP, Genus, Family) %>% summarise(GrForm = unique(GF)[1], DecidEver = unique(Decid.E.green)[1], NeedleBroad = unique(Needle.Broad.lf)[1], C3.C4 = unique(C3C4)[1], slog.LL = mean(log.LL, na.rm=T), slog.LMA = mean(log.LMA, na.rm=T), slog.Nmass = mean(log.Nmass, na.rm=T),slog.Narea = mean(log.Narea, na.rm=T)
#                                                                         ,rslog.LL = log(mean(10^log.LL, na.rm=T),base=10), rslog.LMA = log(mean(10^log.LMA, na.rm=T),base=10), rslog.Nmass = log(mean(10^log.Nmass, na.rm=T),base=10) )
# # have to rename things slog to keep for summarise, but want them just as log
# colnames(LESspp) <- gsub("slog", "log", colnames(LESspp))
# # xtabs(~Species, LES)[which(xtabs(~Species, LES)>1)]
# # ~500 entries are doubled species. not many of them have more than 5 entries though. 
# # evidently some of those have conflicting binary IDs
# ### code for cleaning out problem genera that have multiple families
# # LESgenprobs <- LESspp %>% group_by(Genus) %>% summarise(nFam= length(unique(na.omit(Family))))
# # prob.gen <- LESgenprobs$Genus[which(LESgenprobs$nFam>1)]                                                   
# # tmp <- LESspp[which(LESspp$Genus %in% prob.gen),]
# ## PNW species means
# traitsspp <- traits %>% group_by(GE.SP, GENUS, Family) %>% summarise(slog.LL = mean(log.LL, na.rm=T), slog.LMA = mean(log.LMA, na.rm=T), slog.Nmass = mean(log.Nmass, na.rm=T),slog.Narea = mean(log.Narea, na.rm=T)
#                                                                      ,rslog.LL = log(mean(10^log.LL, na.rm=T),base=10), rslog.LMA = log(mean(10^log.LMA, na.rm=T),base=10), rslog.Nmass = log(mean(10^log.Nmass, na.rm=T),base=10))
# # get species names with spaces that can be merged with LES style names
# traitsspp$Species <- gsub("\\."," ", traitsspp$GE.SP)
# traitsspp$GrForm <- rep("T", times=nrow(traitsspp))
# shrubs <- c("Rhododendron.macrophyllum", "Ribes.divaricatum","Frangula.purshiana","Holodiscus.discolor","Purshia.tridentate" ,
#             "Cercocarpus.unknown","Corylus.cornuta","Ceanothus.velutinus")
# ## 03.20.17, removed Arbutus and Alnus, Lithocarpus (actually now Notholighocarpus), Calocedrus, cause they're trees. Chrysolepis.chrysophylla is now a tree too, cause it can be both. Cornus.unknown is also prob a tree.
# traitsspp$GrForm[which(traitsspp$GE.SP %in% shrubs)] <- "S"
# traitsspp$DecidEver <- "E"
# # we'll just super quickly call everything in shrubs decid and ignore the rest ***** NOT CORRECT ********
# traitsspp$DecidEver[which(traitsspp$GE.SP %in% c(shrubs,"Acer.circinatum", "Acer.macrophyllum","Larix.occidentalis"  ))] <- "D"
# traitsspp$NeedleBroad <- "NA"
# traitsspp$C3.C4 <- "C3"
# colnames(traitsspp)[2] <- "Genus"
# colnames(traitsspp) <- gsub("slog", "log", colnames(traitsspp))
# # reorder columns to match LESspp
# traitsspp <- traitsspp[,match(colnames(LESspp), colnames(traitsspp))]
# 
# ### combine LES species means and traits species means
# allspp <- rbind(data.frame(LESspp), data.frame(traitsspp))
# #*** This has a handful of duplicate spp. But I haven't solved it yet.





#### Genus mean dataframe ####
allgen <- allspp %>% group_by(Genus)  %>% summarise(Fam = sort(unique(Family))[1],  glog.LL = mean(log.LL, na.rm=T), glog.LMA = mean(log.LMA, na.rm=T), glog.Nmass = mean(log.Nmass, na.rm=T),glog.Narea = mean(log.Narea, na.rm=T)
                                                    ,rglog.LL = log(mean(10^rlog.LL, na.rm=T),base=10), rglog.LMA = log(mean(10^rlog.LMA, na.rm=T),base=10), rglog.Nmass = log(mean(10^rlog.Nmass, na.rm=T),base=10), nspp = n() 
                                                    ,MAT = mean(MAT, na.rm=T), MAP=mean(MAP, na.rm=T), VPD=mean(VPD, na.rm=T) )
# taking mean of raw values or logged values doesn't matter all that much yet
colnames(allgen) <- gsub("glog", "log", colnames(allgen))
colnames(allgen)[2] <- "Family"
  # 939 genera from 211 families
# allgen <- allspp %>% group_by(Genus)  %>% summarise(Fam = sort(unique(Family))[1], Needle.Broad = unique(na.omit(NeedleBroad))[1], glog.LL = mean(log.LL, na.rm=T), glog.LMA = mean(log.LMA, na.rm=T), glog.Nmass = mean(log.Nmass, na.rm=T),glog.Narea = mean(log.Narea, na.rm=T)
#                                                     ,rglog.LL = log(mean(10^rlog.LL, na.rm=T),base=10), rglog.LMA = log(mean(10^rlog.LMA, na.rm=T),base=10), rglog.Nmass = log(mean(10^rlog.Nmass, na.rm=T),base=10), nspp = n() )



#### Family Mean dataframe ####
allfam <- allgen %>% group_by(Family) %>% summarise(flog.LL = mean(log.LL, na.rm=T), flog.LMA = mean(log.LMA, na.rm=T), flog.Nmass = mean(log.Nmass, na.rm=T),flog.Narea = mean(log.Narea, na.rm=T)
                                                    ,rflog.LL = log(mean(10^rlog.LL, na.rm=T),base=10), rflog.LMA = log(mean(10^rlog.LMA, na.rm=T),base=10), rflog.Nmass = log(mean(10^rlog.Nmass, na.rm=T),base=10)
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





# 
# 
# nsp <- spp.data %>% group_by (Species) %>% summarise(nLL = length(which(!is.na(log.LL))), nLMA = length(which(!is.na(log.LMA))), nN = length(which(!is.na(log.Nmass))))
# # apply(nsp[,2:4],MARGIN = 2, FUN=mean)
# ngen <- gen.data %>% group_by (Genus) %>% summarise(nLL = length(which(!is.na(log.LL))), nLMA = length(which(!is.na(log.LMA))), nN = length(which(!is.na(log.Nmass))))
# # apply(ngen[,2:4],MARGIN = 2, FUN=mean)
# nfam <- sppinfam.data %>% group_by (Family) %>% summarise(nLL = length(which(!is.na(log.LL))), nLMA = length(which(!is.na(log.LMA))), nN = length(which(!is.na(log.Nmass))))
# # apply(nfam[,2:4],MARGIN = 2, FUN=mean)
# nfamg <- geninfam.data %>% group_by (Family) %>% summarise(nLL = length(which(!is.na(log.LL))), nLMA = length(which(!is.na(log.LMA))), nN = length(which(!is.na(log.Nmass))))
# # apply(nfamg[,2:4],MARGIN = 2, FUN=mean)
# 
# 
# 
# # Plotting what happens if you average logged versus unlogged. Turns out it shifts things, but doesn't really change anything....
# ggplot(allspp , aes(x=log.Nmass, y=log.LMA)) +
#   geom_point(col="grey") +
#   geom_point(data=LESgen, size=2) +
#   geom_point(data=LESfam, col="darkred", aes(size=log(tnspp)) ,alpha=1/2)
# 
# 
# ### family on top of genus on top of spp
# ggplot(LES, aes(x=log.Nmass, y=log.LMA)) + geom_point(col="grey") + geom_point(data=LESgen, aes(x=rlog.Nmass, y=rlog.LMA), size=3) + geom_point(data=LESfam[which(LESfam$tnspp>2),], aes(x=rlog.Nmass,y=rlog.LMA, size=ngen), col='darkred')
# 
# ggplot(data.all, aes(x=log.Nmass, y=log.LMA)) + geom_point(col="grey") + geom_smooth(col="black", method="lm",se=FALSE) +
#   geom_point(data=data.all, aes(col=Genus)) + 
#   geom_point(data=LESfam[which(LESfam$tnspp>2),], aes(x=rlog.Nmass,y=rlog.LMA, size=ngen), col='darkred')









# #LESgens <- LES[which(LES$Genus %in% names(which(xtabs(~Genus, LES)>4))), ]
# ### plot the LES + species.means from PACNW.
# ggplot(LES, aes(x=log.LMA, y=log.Nmass, col=Needle.Broad.lf)) + geom_point() + geom_smooth(method = "lm") + 
#  geom_point(data=allspp, aes(x=log.LMA, y=log.Nmass, col=NULL))
# 
# 
# ### looking at species w/>5 points in the GLOPNET dataset
# p1 <- ggplot(LES, aes(x=log.LL, y=log.LMA)) + geom_point(color="grey") + geom_point(data=traits.common, aes(y=log.LMA_PSA, colour=SP.ID), alpha=1/5) +
#   geom_smooth(method="lm", se=F, colour='darkgrey') + geom_smooth(data=traits.common, aes(y=log.LMA_PSA, colour=SP.ID),method="lm", se=F) +
#   geom_point(data=LES[which(LES$Species %in% commonspp),], aes(colour=Species)) + geom_smooth(data=LES[which(LES$Species %in% commonspp),], aes(colour=Species), method="lm",se=F) +
#   theme(legend.position="none")
# 
# 
# ### looking at across genus relationships
# p2 <- ggplot(LES, aes(x=log.LL, y=log.LMA)) +
#   geom_point(color="grey") + geom_smooth(method="lm", se=F, colour='darkgrey') +
#   geom_smooth(data=LES[which(LES$Genus %in% commongna),], aes(colour=Genus), method="lm",se=F) +
#   theme(legend.position="none")
# 
# geom_point(data=LES[which(LES$Genus %in% commongna),], aes(colour=Genus)) +
#   
#   multiplot(p1,p2, cols=2)
# 
# 
# ######### LL vs LMA
# 
# ggplot(LES, aes(x=log.LMA, y=log.LL)) + geom_point(color="grey") + geom_point(data=traits.common, aes(x=log.LMA_PSA, colour=SP.ID), alpha=1/5) +
#   geom_smooth(method="lm", se=F, colour='darkgrey') + geom_smooth(data=traits.common, aes(x=log.LMA_PSA, colour=SP.ID),method="lm", se=F) +
#   geom_point(data=LES[which(LES$Species %in% commonspp),], aes(colour=Species)) + geom_smooth(data=LES[which(LES$Species %in% commonspp),], aes(colour=Species), method="lm",se=F)
# 
# 
# ### looking at across genus relationships
# ggplot(LES, aes(x=log.LMA, y=log.LL)) + geom_point(color="grey") + geom_smooth(method="lm", se=F, colour='darkgrey') +
#   geom_point(data=LES[which(LES$Genus %in% commongen),], aes(colour=Genus)) + geom_smooth(data=LES[which(LES$Genus %in% commongen),], aes(colour=Genus), method="lm",se=F)
# 
# win.spp.mod <- lmer(log.LL~log.LMA + (1|SP.ID) + (0 + log.LMA|SP.ID), traits.common)
# win.genus.mod <- lmer(log.LL~log.LMA  + (1|Genus) + (0 + log.LMA|Genus), LES)
# summary(win.spp.mod)
# summary(win.genus.mod)
# summary(lm(log.LL~log.LMA, LES))
# 








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
    return(rep(NA, times=7))
  }
  else{
    if(var(data[,yvar], na.rm=T)==0 | var(data[,xvar],na.rm=T)==0){
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
ptm <- proc.time()
set.seed(42)
spp.results <- data.frame(matrix(NA, nrow=length(unique(spp.data$Species)), ncol=15))
colnames(spp.results) <- c("Species", "Int","Slope","Rho","r.sq","n","varLMA","varNmass", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(spp.data$Species))){
  species <- levels(spp.data$Species)[i]
  print(species)
  dataz <- spp.data[which(spp.data$Species==species),]
  res <- fit.MAR(xvar='log.LMA',yvar="log.Nmass",data=dataz)
  spp.results[i,1] <- species
  spp.results[i,2:8] <- res
  if (!is.na(res[1]) & res[5]>4){ # only fit null model if there are >5 data points
    nullbounds <- fit.null(xvar='log.LMA', yvar="log.Nmass", observed = dataz, nulldata = allspp, nits = 1000)  
    spp.results[i, 9:14] <- nullbounds
    spp.results[i,15] <- test.sig(x=spp.results$Rho[i], test=nullbounds)
  }
}
proc.time()-ptm





# Seems to be working!!!

############ .Genus level analysis #####################
t0 <- proc.time()
gen.results <- data.frame(matrix(NA, nrow=length(unique(gen.data$Genus)), ncol=15))
colnames(gen.results) <- c("Genus", "Int","Slope","Rho","r.sq","n","varLMA","varNmass", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(gen.data$Genus))){
  genus <- levels(gen.data$Genus)[i]
  print(genus)
  dataz <- gen.data[which(gen.data$Genus==genus),]
  gen.results[i,1] <- genus
  res <- fit.MAR(xvar='log.LMA',yvar="log.Nmass",data=dataz)
  gen.results[i,2:8] <- res 
  if (!is.na(res[1]) & res[5]>4){ # only fit null model if there are >5 data points, and the fit.MAR actually ran
    nullbounds <- fit.null(xvar='log.LMA', yvar="log.Nmass", observed = dataz, nulldata = allgen, nits = 1000)  
    gen.results[i, 9:14] <- nullbounds
    gen.results[i,15] <- test.sig(x=gen.results$Rho[i], test=nullbounds)
  }
  
}
proc.time()-t0


############ .spp w/in Family level analysis #####################
t0<- proc.time()
sppinfam.results <- data.frame(matrix(NA, nrow=length(unique(sppinfam.data$Family)), ncol=15))
colnames(sppinfam.results) <- c("Family", "Int","Slope","Rho","r.sq","n","varLMA","varNmass", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(sppinfam.data$Family))){
  family <- levels(sppinfam.data$Family)[i]
  print(family)
  dataz <- sppinfam.data[which(sppinfam.data$Family==family),]
  res <- fit.MAR(xvar='log.LMA',yvar="log.Nmass",data=dataz)
  sppinfam.results[i,1] <- family
  sppinfam.results[i,2:8] <- res
  if (!is.na(res[1]) & res[5]>4){ # only fit null model if there are >5 data points
    nullbounds <- fit.null(xvar='log.LMA', yvar="log.Nmass", observed = dataz, nulldata = allfam, nits = 1000)  
    sppinfam.results[i, 9:14] <- nullbounds
    sppinfam.results[i,15] <- test.sig(x=sppinfam.results$Rho[i], test=nullbounds)
  }
  
}
proc.time()-t0

############ .gen w/in Family level analysis #####################
t0<-proc.time()
geninfam.results <- data.frame(matrix(NA, nrow=length(unique(geninfam.data$Family)), ncol=15))
colnames(geninfam.results) <- c("Family", "Int","Slope","Rho","r.sq","n","varLMA","varNmass", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(geninfam.data$Family))){
  family <- levels(geninfam.data$Family)[i]
  print(family)
  dataz <- geninfam.data[which(geninfam.data$Family==family),]
  res <- fit.MAR(xvar='log.LMA',yvar="log.Nmass",data=dataz)
  geninfam.results[i,1] <- family
  geninfam.results[i,2:8] <- res
  if (!is.na(res[1]) & res[5]>4){ # only fit null model if there are >5 data points
    nullbounds <- fit.null(xvar='log.LMA', yvar="log.Nmass", observed = dataz, nulldata = allfam, nits = 1000)  
    geninfam.results[i, 9:14] <- nullbounds
    geninfam.results[i,15] <- test.sig(x=geninfam.results$Rho[i], test=nullbounds)
  }
}
proc.time()-t0

############ .Family level analysis #####################

# # currently just working with LES until I combine the PACNW dataset into this.
# famcorr.all <- cor(x=fam.data[,c("log.LL","log.LMA","log.Nmass")], use = "pairwise.complete.obs")
# famcorr.clean <- cor(x=fam.dataclean[,c("log.LL","log.LMA","log.Nmass")], use = "pairwise.complete.obs")

fam.res_LMA.N <- c("fam.all", fit.MAR(xvar='log.LMA',yvar="log.Nmass",data=fam.data),rep(NA, times=7),"Fam") 
names(fam.res_LMA.N) <- c("Taxo.Unit","Int","Slope","Rho","r.sq", "n","varLMA","varNmass", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig" ,"Type")
fam.resclean_LMA.N <- c("fam.clean", fit.MAR(xvar='log.LMA',yvar="log.Nmass",data=fam.dataclean),rep(NA, times=7), "Fam.clean")
names(fam.resclean_LMA.N) <- c("Taxo.Unit","Int","Slope","Rho","r.sq", "n","varLMA","varNmass", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Type")
global_LMA.N <- c("global", fit.MAR(xvar='log.LMA',yvar="log.Nmass",data=allspp),rep(NA, times=7), "global")
names(global_LMA.N) <- c("Taxo.Unit","Int","Slope","Rho","r.sq", "n","varLMA","varNmass", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Type")


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


spp.results <- data.frame(matrix(NA, nrow=length(unique(spp.data$Species)), ncol=15))
colnames(spp.results) <- c("Species", "Int","Slope","Rho","r.sq","n","varLMA","varLL","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(spp.data$Species))){
  species <- levels(spp.data$Species)[i]
  print(species)
  dataz <- spp.data[which(spp.data$Species==species),]
  res <- fit.MAR(xvar='log.LMA',yvar="log.LL",data=dataz)
  spp.results[i,1] <- species
  spp.results[i,2:8] <- res
  if (!is.na(res[1]) &res[5]>4){ # only fit null model if there are >5 data points
    nullbounds <- fit.null(xvar='log.LMA', yvar="log.LL", observed = dataz, nulldata = allspp, nits = 1000)  
    spp.results[i, 9:14] <- nullbounds
    spp.results[i,15] <- test.sig(x=spp.results$Rho[i], test=nullbounds)
  }
}


############ .Genus level analysis #####################


gen.results <- data.frame(matrix(NA, nrow=length(unique(gen.data$Genus)), ncol=15))
colnames(gen.results) <- c("Genus", "Int","Slope","Rho","r.sq","n","varLMA","varLL","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(gen.data$Genus))){
  genus <- levels(gen.data$Genus)[i]
  print(genus)
  dataz <- gen.data[which(gen.data$Genus==genus),]
  gen.results[i,1] <- genus
  res <- fit.MAR(xvar='log.LMA',yvar="log.LL",data=dataz)
  gen.results[i,2:8] <- res 
  if (!is.na(res[1]) & res[5]>4){ # only fit null model if there are >5 data points, and the fit.MAR actually ran
    nullbounds <- fit.null(xvar='log.LMA', yvar="log.LL", observed = dataz, nulldata = allgen, nits = 1000)  
    gen.results[i, 9:14] <- nullbounds
    gen.results[i,15] <- test.sig(x=gen.results$Rho[i], test=nullbounds)
  }
}



############ .spp w/in Family level analysis #####################


sppinfam.results <- data.frame(matrix(NA, nrow=length(unique(sppinfam.data$Family)), ncol=15))
colnames(sppinfam.results) <- c("Family", "Int","Slope","Rho","r.sq","n","varLMA","varLL","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(sppinfam.data$Family))){
  family <- levels(sppinfam.data$Family)[i]
  print(family)
  dataz <- sppinfam.data[which(sppinfam.data$Family==family),]
  res <- fit.MAR(xvar='log.LMA',yvar="log.LL",data=dataz)
  sppinfam.results[i,1] <- family
  sppinfam.results[i,2:8] <- res
  if (!is.na(res[1]) & res[5]>4){ # only fit null model if there are >5 data points
    nullbounds <- fit.null(xvar='log.LMA', yvar="log.LL", observed = dataz, nulldata = allfam, nits = 1000)  
    sppinfam.results[i, 9:14] <- nullbounds
    sppinfam.results[i,15] <- test.sig(x=sppinfam.results$Rho[i], test=nullbounds)
  }
}


############ .gen w/in Family level analysis #####################


geninfam.results <- data.frame(matrix(NA, nrow=length(unique(geninfam.data$Family)), ncol=15))
colnames(geninfam.results) <- c("Family", "Int","Slope","Rho","r.sq","n","varLMA","varLL","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(geninfam.data$Family))){
  family <- levels(geninfam.data$Family)[i]
  print(family)
  dataz <- geninfam.data[which(geninfam.data$Family==family),]
  res <- fit.MAR(xvar='log.LMA',yvar="log.LL",data=dataz)
  geninfam.results[i,1] <- family
  geninfam.results[i,2:8] <- res
  if (!is.na(res[1]) & res[5]>4){ # only fit null model if there are >5 data points
    nullbounds <- fit.null(xvar='log.LMA', yvar="log.LL", observed = dataz, nulldata = allfam, nits = 1000)  
    geninfam.results[i, 9:14] <- nullbounds
    geninfam.results[i,15] <- test.sig(x=geninfam.results$Rho[i], test=nullbounds)
  }
}


############ .Family level analysis #####################

fam.res_LMA.LL <- c("fam.all", fit.MAR(xvar='log.LMA',yvar="log.LL",data=fam.data),rep(NA, times=7),"Fam")
names(fam.res_LMA.LL) <- c("Taxo.Unit","Int","Slope","Rho","r.sq", "n","varLMA","varLL","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Type")
fam.resclean_LMA.LL <- c("fam.clean", fit.MAR(xvar='log.LMA',yvar="log.LL",data=fam.dataclean),rep(NA, times=7), "Fam.clean")
names(fam.resclean_LMA.LL) <- c("Taxo.Unit","Int","Slope","Rho","r.sq", "n","varLMA","varLL","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Type")
global_LMA.LL <- c("global", fit.MAR(xvar='log.LMA',yvar="log.LL",data=allspp),rep(NA, times=7), "global")
names(global_LMA.LL) <- c("Taxo.Unit","Int","Slope","Rho","r.sq", "n","varLMA","varLL","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Type")



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



spp.results <- data.frame(matrix(NA, nrow=length(unique(spp.data$Species)), ncol=15))
colnames(spp.results) <- c("Species", "Int","Slope","Rho","r.sq","n","varLL","varNmass","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(spp.data$Species))){
  species <- levels(spp.data$Species)[i]
  print(species)
  dataz <- spp.data[which(spp.data$Species==species),]
  res <- fit.MAR(xvar='log.LL',yvar="log.Nmass",data=dataz)
  spp.results[i,1] <- species
  spp.results[i,2:8] <- res
  if (!is.na(res[1]) &res[5]>4){ # only fit null model if there are >5 data points
    nullbounds <- fit.null(xvar='log.LL', yvar="log.Nmass", observed = dataz, nulldata = allspp, nits = 1000)  
    spp.results[i, 9:14] <- nullbounds
    spp.results[i,15] <- test.sig(x=spp.results$Rho[i], test=nullbounds)
  }
}


############ .Genus level analysis #####################

gen.results <- data.frame(matrix(NA, nrow=length(unique(gen.data$Genus)), ncol=15))
colnames(gen.results) <- c("Genus", "Int","Slope","Rho","r.sq","n","varLL","varNmass","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(gen.data$Genus))){
  genus <- levels(gen.data$Genus)[i]
  print(genus)
  dataz <- gen.data[which(gen.data$Genus==genus),]
  gen.results[i,1] <- genus
  res <- fit.MAR(xvar='log.LL',yvar="log.Nmass",data=dataz)
  gen.results[i,2:8] <- res 
  if (!is.na(res[1]) & res[5]>4){ # only fit null model if there are >5 data points, and the fit.MAR actually ran
    nullbounds <- fit.null(xvar='log.LL', yvar="log.Nmass", observed = dataz, nulldata = allgen, nits = 1000)  
    gen.results[i, 9:14] <- nullbounds
    gen.results[i,15] <- test.sig(x=gen.results$Rho[i], test=nullbounds)
  }
}



############ .spp w/in Family level analysis #####################

sppinfam.results <- data.frame(matrix(NA, nrow=length(unique(sppinfam.data$Family)), ncol=15))
colnames(sppinfam.results) <- c("Family", "Int","Slope","Rho","r.sq","n","varLL","varNmass","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(sppinfam.data$Family))){
  family <- levels(sppinfam.data$Family)[i]
  print(family)
  dataz <- sppinfam.data[which(sppinfam.data$Family==family),]
  res <- fit.MAR(xvar='log.LL',yvar="log.Nmass",data=dataz)
  sppinfam.results[i,1] <- family
  sppinfam.results[i,2:8] <- res
  if (!is.na(res[1]) & res[5]>4){ # only fit null model if there are >5 data points
    nullbounds <- fit.null(xvar='log.LL', yvar="log.Nmass", observed = dataz, nulldata = allfam, nits = 1000)  
    sppinfam.results[i, 9:14] <- nullbounds
    sppinfam.results[i,15] <- test.sig(x=sppinfam.results$Rho[i], test=nullbounds)
  }
}





############ .gen w/in Family level analysis #####################
geninfam.results <- data.frame(matrix(NA, nrow=length(unique(geninfam.data$Family)), ncol=15))
colnames(geninfam.results) <- c("Family", "Int","Slope","Rho","r.sq","n","varLL","varNmass","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(geninfam.data$Family))){
  family <- levels(geninfam.data$Family)[i]
  print(family)
  dataz <- geninfam.data[which(geninfam.data$Family==family),]
  res <- fit.MAR(xvar='log.LL',yvar="log.Nmass",data=dataz)
  geninfam.results[i,1] <- family
  geninfam.results[i,2:8] <- res
  if (!is.na(res[1]) & res[5]>4){ # only fit null model if there are >5 data points
    nullbounds <- fit.null(xvar='log.LL', yvar="log.Nmass", observed = dataz, nulldata = allfam, nits = 1000)  
    geninfam.results[i, 9:14] <- nullbounds
    geninfam.results[i,15] <- test.sig(x=geninfam.results$Rho[i], test=nullbounds)
  }
}





############ .Family level analysis #####################


fam.res_LL.N <- c("fam.all", fit.MAR(xvar='log.LL',yvar="log.Nmass",data=fam.data),rep(NA, times=7),"Fam")
names(fam.res_LL.N) <- c("Taxo.Unit","Int","Slope","Rho","r.sq", "n","varLL","varNmass","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Type")
fam.resclean_LL.N <- c("fam.clean", fit.MAR(xvar='log.LL',yvar="log.Nmass",data=fam.dataclean),rep(NA, times=7), "Fam.clean")
names(fam.res_LL.N) <- c("Taxo.Unit","Int","Slope","Rho","r.sq", "n","varLL","varNmass","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Type")
global_LL.N <- c("global", fit.MAR(xvar='log.LL',yvar="log.Nmass",data=allspp),rep(NA, times=7), "global")
names(global_LL.N) <- c("Taxo.Unit","Int","Slope","Rho","r.sq", "n","varLL","varNmass","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Type")





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



spp.results <- data.frame(matrix(NA, nrow=length(unique(spp.data$Species)), ncol=15))
colnames(spp.results) <- c("Species", "Int","Slope","Rho","r.sq","n","varLMA","varNarea","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(spp.data$Species))){
  species <- levels(spp.data$Species)[i]
  print(species)
  dataz <- spp.data[which(spp.data$Species==species),]
  res <- fit.MAR(xvar='log.LMA',yvar="log.Narea",data=dataz)
  spp.results[i,1] <- species
  spp.results[i,2:8] <- res
  if (!is.na(res[1]) &res[5]>4){ # only fit null model if there are >5 data points
    nullbounds <- fit.null(xvar='log.LMA', yvar="log.Narea", observed = dataz, nulldata = allspp, nits = 1000)  
    spp.results[i, 9:14] <- nullbounds
    spp.results[i,15] <- test.sig(x=spp.results$Rho[i], test=nullbounds)
  }
}


############ .Genus level analysis #####################

gen.results <- data.frame(matrix(NA, nrow=length(unique(gen.data$Genus)), ncol=15))
colnames(gen.results) <- c("Genus", "Int","Slope","Rho","r.sq","n","varLMA","varNarea","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(gen.data$Genus))){
  genus <- levels(gen.data$Genus)[i]
  print(genus)
  dataz <- gen.data[which(gen.data$Genus==genus),]
  gen.results[i,1] <- genus
  res <- fit.MAR(xvar='log.LMA',yvar="log.Narea",data=dataz)
  gen.results[i,2:8] <- res 
  if (!is.na(res[1]) & res[5]>4){ # only fit null model if there are >5 data points, and the fit.MAR actually ran
    nullbounds <- fit.null(xvar='log.LMA', yvar="log.Narea", observed = dataz, nulldata = allgen, nits = 1000)  
    gen.results[i, 9:14] <- nullbounds
    gen.results[i,15] <- test.sig(x=gen.results$Rho[i], test=nullbounds)
  }
}



############ .spp w/in Family level analysis #####################

sppinfam.results <- data.frame(matrix(NA, nrow=length(unique(sppinfam.data$Family)), ncol=15))
colnames(sppinfam.results) <- c("Family", "Int","Slope","Rho","r.sq","n","varLMA","varNarea","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(sppinfam.data$Family))){
  family <- levels(sppinfam.data$Family)[i]
  print(family)
  dataz <- sppinfam.data[which(sppinfam.data$Family==family),]
  res <- fit.MAR(xvar='log.LMA',yvar="log.Narea",data=dataz)
  sppinfam.results[i,1] <- family
  sppinfam.results[i,2:8] <- res
  if (!is.na(res[1]) & res[5]>4){ # only fit null model if there are >5 data points
    nullbounds <- fit.null(xvar='log.LMA', yvar="log.Narea", observed = dataz, nulldata = allfam, nits = 1000)  
    sppinfam.results[i, 9:14] <- nullbounds
    sppinfam.results[i,15] <- test.sig(x=sppinfam.results$Rho[i], test=nullbounds)
  }
}





############ .gen w/in Family level analysis #####################

geninfam.results <- data.frame(matrix(NA, nrow=length(unique(geninfam.data$Family)), ncol=15))
colnames(geninfam.results) <- c("Family", "Int","Slope","Rho","r.sq","n","varLMA","varNarea","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(geninfam.data$Family))){
  family <- levels(geninfam.data$Family)[i]
  print(family)
  dataz <- geninfam.data[which(geninfam.data$Family==family),]
  res <- fit.MAR(xvar='log.LMA',yvar="log.Narea",data=dataz)
  geninfam.results[i,1] <- family
  geninfam.results[i,2:8] <- res
  if (!is.na(res[1]) & res[5]>4){ # only fit null model if there are >5 data points
    nullbounds <- fit.null(xvar='log.LMA', yvar="log.Narea", observed = dataz, nulldata = allfam, nits = 1000)  
    geninfam.results[i, 9:14] <- nullbounds
    geninfam.results[i,15] <- test.sig(x=geninfam.results$Rho[i], test=nullbounds)
  }
}





############ .Family level analysis #####################


fam.res_LMA.Narea <- c("fam.all", fit.MAR(xvar='log.LMA',yvar="log.Narea",data=fam.data),rep(NA, times=7),"Fam")
names(fam.res_LMA.Narea) <- c("Taxo.Unit","Int","Slope","Rho","r.sq", "n","varLMA","varNarea","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Type")
fam.resclean_LMA.Narea <- c("fam.clean", fit.MAR(xvar='log.LMA',yvar="log.Narea",data=fam.dataclean),rep(NA, times=7), "Fam.clean")
names(fam.res_LMA.Narea) <- c("Taxo.Unit","Int","Slope","Rho","r.sq", "n","varLMA","varNarea","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Type")
global_LMA.Narea <- c("global", fit.MAR(xvar='log.LMA',yvar="log.Narea",data=allspp),rep(NA, times=7), "global")
names(global_LMA.Narea) <- c("Taxo.Unit","Int","Slope","Rho","r.sq", "n","varLMA","varNarea","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Type")




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



spp.results <- data.frame(matrix(NA, nrow=length(unique(spp.data$Species)), ncol=15))
colnames(spp.results) <- c("Species", "Int","Slope","Rho","r.sq","n","varLL","varNarea","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(spp.data$Species))){
  species <- levels(spp.data$Species)[i]
  print(species)
  dataz <- spp.data[which(spp.data$Species==species),]
  res <- fit.MAR(xvar="log.LL", yvar='log.Narea',data=dataz)
  spp.results[i,1] <- species
  spp.results[i,2:8] <- res
  if (!is.na(res[1]) &res[5]>4){ # only fit null model if there are >5 data points
    nullbounds <- fit.null(xvar='log.LL', yvar="log.Narea", observed = dataz, nulldata = allspp, nits = 1000)  
    spp.results[i, 9:14] <- nullbounds
    spp.results[i,15] <- test.sig(x=spp.results$Rho[i], test=nullbounds)
  }
}


############ .Genus level analysis #####################

gen.results <- data.frame(matrix(NA, nrow=length(unique(gen.data$Genus)), ncol=15))
colnames(gen.results) <- c("Genus", "Int","Slope","Rho","r.sq","n","varLL","varNarea","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(gen.data$Genus))){
  genus <- levels(gen.data$Genus)[i]
  print(genus)
  dataz <- gen.data[which(gen.data$Genus==genus),]
  gen.results[i,1] <- genus
  res <- fit.MAR(xvar="log.LL",yvar='log.Narea',data=dataz)
  gen.results[i,2:8] <- res 
  if (!is.na(res[1]) & res[5]>4){ # only fit null model if there are >5 data points, and the fit.MAR actually ran
    nullbounds <- fit.null(xvar='log.LL', yvar="log.Narea", observed = dataz, nulldata = allgen, nits = 1000)  
    gen.results[i, 9:14] <- nullbounds
    gen.results[i,15] <- test.sig(x=gen.results$Rho[i], test=nullbounds)
  }
}



############ .spp w/in Family level analysis #####################

sppinfam.results <- data.frame(matrix(NA, nrow=length(unique(sppinfam.data$Family)), ncol=15))
colnames(sppinfam.results) <- c("Family", "Int","Slope","Rho","r.sq","n","varLL","varNarea","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(sppinfam.data$Family))){
  family <- levels(sppinfam.data$Family)[i]
  print(family)
  dataz <- sppinfam.data[which(sppinfam.data$Family==family),]
  res <- fit.MAR(xvar='log.LL',yvar="log.Narea",data=dataz)
  sppinfam.results[i,1] <- family
  sppinfam.results[i,2:8] <- res
  if (!is.na(res[1]) & res[5]>4){ # only fit null model if there are >5 data points
    nullbounds <- fit.null(xvar='log.LL', yvar="log.Narea", observed = dataz, nulldata = allfam, nits = 1000)  
    sppinfam.results[i, 9:14] <- nullbounds
    sppinfam.results[i,15] <- test.sig(x=sppinfam.results$Rho[i], test=nullbounds)
  }
}





############ .gen w/in Family level analysis #####################

geninfam.results <- data.frame(matrix(NA, nrow=length(unique(geninfam.data$Family)), ncol=15))
colnames(geninfam.results) <- c("Family", "Int","Slope","Rho","r.sq","n","varLL","varNarea","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(geninfam.data$Family))){
  family <- levels(geninfam.data$Family)[i]
  print(family)
  dataz <- geninfam.data[which(geninfam.data$Family==family),]
  res <- fit.MAR(xvar='log.LL',yvar="log.Narea",data=dataz)
  geninfam.results[i,1] <- family
  geninfam.results[i,2:8] <- res
  if (!is.na(res[1]) & res[5]>4){ # only fit null model if there are >5 data points
    nullbounds <- fit.null(xvar='log.LL', yvar="log.Narea", observed = dataz, nulldata = allfam, nits = 1000)  
    geninfam.results[i, 9:14] <- nullbounds
    geninfam.results[i,15] <- test.sig(x=geninfam.results$Rho[i], test=nullbounds)
  }
}





############ .Family level analysis #####################


fam.res_LL.Narea <- c("fam.all", fit.MAR(xvar="log.LL",yvar='log.Narea',data=fam.data),rep(NA, times=7),"Fam")
names(fam.res_LL.Narea) <- c("Taxo.Unit","Int","Slope","Rho","r.sq", "n","varLL","varNarea","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Type")
fam.resclean_LL.Narea <- c("fam.clean", fit.MAR(xvar="log.LL",yvar='log.Narea',data=fam.dataclean),rep(NA, times=7), "Fam.clean")
names(fam.res_LL.Narea) <- c("Taxo.Unit","Int","Slope","Rho","r.sq", "n","varLL","varNarea","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Type")
global_LL.Narea <- c("global", fit.MAR(xvar="log.LL",yvar='log.Narea',data=allspp),rep(NA, times=7), "global")
names(global_LL.Narea) <- c("Taxo.Unit","Int","Slope","Rho","r.sq", "n","varLL","varNarea","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Type")



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
colnames(LMALL)[c(2:6, 9:15)] <- paste(colnames(all.results.LMALL)[c(2:6, 9:15)],"LMA.LL", sep="_")
LMAN <- all.results.LMAN
colnames(LMAN)[c(2:6, 9:15)] <- paste(colnames(all.results.LMAN)[c(2:6, 9:15)],"LMA.N", sep="_")
LLNmass <- all.results.LLNmass
colnames(LLNmass)[c(2:6, 9:15)] <- paste(colnames(all.results.LLNmass)[c(2:6, 9:15)],"LL.N", sep="_")
LMANarea <- all.results.LMANarea
colnames(LMANarea)[c(2:6, 9:15)] <- paste(colnames(all.results.LMANarea)[c(2:6, 9:15)],"LMA.Narea", sep="_")
LLNarea <- all.results.LLNarea
colnames(LLNarea)[c(2:6, 9:15)] <- paste(colnames(all.results.LLNarea)[c(2:6, 9:15)],"LL.Narea", sep="_")

all.results <- cbind(LMALL, LMAN[,-c(1,7,16)], LLNmass[,-c(1,7,8,16)], LMANarea[,-c(1,7,16)], LLNarea[,-c(1,7,8,16)]) # drop the duplicate 'var' columns

#write.csv(all.results, "Results_SimpleMAreg_v1_030817.csv")
#write.csv(all.results, "Results_SimpleMAreg_v2_031417.csv")
#write.csv(all.results, "Results_SimpleMAreg_v3_031717.csv")
#write.csv(all.results, "Results_SimpleMAreg_v4_040117.csv")
#write.csv(all.results, "Results_SimpleMAreg_v5_040217.csv") # updated with NareaLL
#write.csv(all.results, "Results_SimpleMAreg_v6_040717.csv") # switched to LLNarea
#write.csv(all.results, "Results_SimpleMAreg_v7_051717.csv") # switched to LLNmass, and added CIs from null model
write.csv(all.results, "Results_SimpleMAreg_v8_051717.csv") # fixed bug that screwed up LL.Narea null model (when taxa had more variance than all families)
#all.resultsold <- read.csv("Results_SimpleMAreg_v2_031417.csv")








#________________________________________________________________________
############# Initial Plotting #################
#________________________________________________________________________


#________________________________________________________________________
########## LOAD RESULTS DATA ##################

all.results <- read.csv("Results_SimpleMAreg_v8_051717.csv", row.names = 1)
levels(all.results$Type) <- list(w.inSpp = "w.inSpp", w.inGen = "w.inGen", Sppw.inFam= "Sppw.inFam",Genw.inFam="Genw.inFam", Fam="Fam",Famclean="Famclean", global="global")

all.results.cl <- all.results %>% filter(Type %in% c("w.inSpp","w.inGen","Genw.inFam","Famclean","global"))
all.results.cl$Type <- factor(all.results.cl$Type)





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
# Rho LMA v LL
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
LLNmass <- lm(Slope_LL.N~Type, all.results)
LLNmassnull <- lm(Slope_LL.N~1, all.results)
anova(LLNmass, LLNmassnull) # p=0.897
LLNmass <- lm(Slope_LL.N~Type, weights = n_LL.N, all.results)
LLNmassnull <- lm(Slope_LL.N~1, weights=n_LL.N, all.results)
anova(LLNmass, LLNmassnull) # p=0.0878
LLNmass <- lm(Slope_LL.N~Type, weights=varNmass, all.results)
LLNmassnull <- lm(Slope_LL.N~1, all.results, weights=varNmass)
anova(LLNmass, LLNmassnull) # p=0.61
LLNmass <- lm(Slope_LL.N~Type, weights=varLL, all.results)
LLNmassnull <- lm(Slope_LL.N~1, all.results, weights=varLL)
anova(LLNmass, LLNmassnull) # p=0.31
#### Rho LL v Nmass
LLNmass <- lm(Rho_LL.N~Type, all.results)
LLNmassnull <- lm(Rho_LL.N~1, all.results)
anova(LLNmass, LLNmassnull) # p=0.018
LLNmass <- lm(Rho_LL.N~Type, weights = n_LL.N, all.results)
LLNmassnull <- lm(Rho_LL.N~1, weights=n_LL.N, all.results)
anova(LLNmass, LLNmassnull) # p<0.0001
LLNmass <- lm(Rho_LL.N~Type, weights=varNmass, all.results)
LLNmassnull <- lm(Rho_LL.N~1, all.results, weights=varNmass)
anova(LLNmass, LLNmassnull) # p=0.005
LLNmass <- lm(Rho_LL.N~Type, weights=varLL, all.results)
LLNmassnull <- lm(Rho_LL.N~1, all.results, weights=varLL)
anova(LLNmass, LLNmassnull) # p<0.0001



#### SLOPE LMA v Narea
LMANarea <- lm(Slope_LMA.Narea~Type, all.results.cl)
LMANareanull <- lm(Slope_LMA.Narea~1, all.results.cl)
anova(LMANarea, LMANareanull) # p=0.0632, cl=0.02992
LMANarea <- lm(Slope_LMA.Narea~Type, weights = n_LMA.Narea, all.results.cl)
LMANareanull <- lm(Slope_LMA.Narea~1, weights=n_LMA.Narea, all.results.cl)
anova(LMANarea, LMANareanull) # p<0.0001 &cl
LMANarea <- lm(Slope_LMA.Narea~Type, weights=varNarea, all.results.cl)
LMANareanull <- lm(Slope_LMA.Narea~1, all.results.cl, weights=varNarea)
anova(LMANarea, LMANareanull) # p=0.549, cl=0.4077 
LMANarea <- lm(Slope_LMA.Narea~Type, weights=varLMA, all.results.cl)
LMANareanull <- lm(Slope_LMA.Narea~1, all.results.cl, weights=varLMA)
anova(LMANarea, LMANareanull) # p=0.135, cl=0.001735
#### Rho LMA v Narea
LMANarea <- lm(Rho_LMA.Narea~Type, all.results.cl)
LMANareanull <- lm(Rho_LMA.Narea~1, all.results.cl)
anova(LMANarea, LMANareanull) # p=0.98, cl=0.8787
LMANarea <- lm(Rho_LMA.Narea~Type, weights = n_LMA.Narea, all.results.cl)
LMANareanull <- lm(Rho_LMA.Narea~1, weights=n_LMA.Narea, all.results.cl)
anova(LMANarea, LMANareanull) # p<0.45, cl=0.2666
LMANarea <- lm(Rho_LMA.Narea~Type, weights=varNarea, all.results.cl)
LMANareanull <- lm(Rho_LMA.Narea~1, all.results.cl, weights=varNarea)
anova(LMANarea, LMANareanull) # p=0.609
LMANarea <- lm(Rho_LMA.Narea~Type, weights=varLMA, all.results.cl)
LMANareanull <- lm(Rho_LMA.Narea~1, all.results.cl, weights=varLMA)
anova(LMANarea, LMANareanull) # p<0.20



#### SLOPE LL v Narea
LLNarea <- lm(Slope_LL.Narea~Type, all.results.cl)
LLNareanull <- lm(Slope_LL.Narea~1, all.results.cl)
anova(LLNarea, LLNareanull) # p=0.2545, cl=0.1181
LLNarea <- lm(Slope_LL.Narea~Type, weights = n_LL.Narea, all.results.cl)
LLNareanull <- lm(Slope_LL.Narea~1, weights=n_LL.Narea, all.results.cl)
anova(LLNarea, LLNareanull) # cl p= 1.4 e-6
LLNarea <- lm(Slope_LL.Narea~Type, weights=varNarea, all.results.cl)
LLNareanull <- lm(Slope_LL.Narea~1, all.results.cl, weights=varNarea)
anova(LLNarea, LLNareanull) #, cl=0.00805 
LLNarea <- lm(Slope_LL.Narea~Type, weights=varLL, all.results.cl)
LLNareanull <- lm(Slope_LL.Narea~1, all.results.cl, weights=varLL)
anova(LLNarea, LLNareanull) # p=, cl=0.00364
#### Rho LL v Narea
LLNarea <- lm(Rho_LL.Narea~Type, all.results.cl)
LLNareanull <- lm(Rho_LL.Narea~1, all.results.cl)
anova(LLNarea, LLNareanull) # p=, cl=0.22
LLNarea <- lm(Rho_LL.Narea~Type, weights = n_LL.Narea, all.results.cl)
LLNareanull <- lm(Rho_LL.Narea~1, weights=n_LL.Narea, all.results.cl)
anova(LLNarea, LLNareanull) # pcl = 1.42 e-8
LLNarea <- lm(Rho_LL.Narea~Type, weights=varNarea, all.results.cl)
LLNareanull <- lm(Rho_LL.Narea~1, all.results.cl, weights=varNarea)
anova(LLNarea, LLNareanull) # p=0.07745
LLNarea <- lm(Rho_LL.Narea~Type, weights=varLL, all.results.cl)
LLNareanull <- lm(Rho_LL.Narea~1, all.results.cl, weights=varLL)
anova(LLNarea, LLNareanull) # p<0.03126


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
# 
# ### plotting two things at once...
# ggplot(data=NULL,aes(x=all.results.LMAN$Rho, y=all.results.NmassLL$Rho, col=all.results.NmassLL$Type, size=all.results.NmassLL$n)) + geom_point()
# 
# 
# ggplot(data=NULL,aes(x=all.results.LMAN$Rho, y=all.results.LMALL$Rho, col=all.results.NmassLL$Type, size=all.results.NmassLL$n)) + geom_point() +
#   geom_abline(slope = 0, intercept=0, col="grey")+ geom_vline(xintercept = 0, col="grey")








######### LMA v LL, angios vs gymnos #######
tmp <- all.results %>% filter(Type=="w.inSpp" & n_LMA.LL>4)
tmp$Taxo.Unit
lh <- c("g","g","g","g",
        "g","g","g","g",
        "g","g","g","g",
        "g","a","g","a","a","g")
llspecies <- tmp$Taxo.Unit
boxplot(Slope_LMA.LL~lh, tmp)
plot(log.LL~log.LMA, LES, pch=16, col="grey")
points(log.LL~log.LMA, spp.data[which(spp.data$Species %in% llspecies),], col=Species)





#________________________________________________________________________
#________________________________________________________________________
###### Old Code: PCA pieces ##########

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









################### LL summaries ########
tmp <- traits.common5 %>% group_by(SP.ID) %>% summarise(mLL = mean(LEAF_LIFE, na.rm=T), minLL = min(LEAF_LIFE, na.rm=T), maxLL = max(LEAF_LIFE, na.rm=T), sdLL = sd(LEAF_LIFE, na.rm=T), cvLL = sd(LEAF_LIFE, na.rm=T)/mean(LEAF_LIFE, na.rm=T), nLL=length(-which(is.na(LEAF_LIFE))))
plot(sdLL~mLL, tmp) # sd is def a function of mean
plot(cvLL~mLL, tmp) # CV is reasonably unrelated to mean

tmp2 <- spp.data %>% group_by(Species) %>% summarise(mLL = mean(10^log.LL, na.rm=T), sdLL = sd(10^log.LL, na.rm=T), cvLL = sd(10^log.LL, na.rm=T)/mean(10^log.LL, na.rm=T), nLL=length(which(10^log.LL > 0))
                                                     ,mLMA = mean(10^log.LMA, na.rm=T), sdLMA = sd(10^log.LMA, na.rm=T), cvLMA = sd(10^log.LMA, na.rm=T)/mean(10^log.LMA, na.rm=T), nLMA=length(which(10^log.LMA > 0))
                                                     ,mNmass = mean(10^log.Nmass, na.rm=T), sdNmass = sd(10^log.Nmass, na.rm=T), cvNmass = sd(10^log.Nmass, na.rm=T)/mean(10^log.Nmass, na.rm=T), nNmass=length(which(10^log.Nmass > 0)))



##### plotting trait CVs as f(trait mean) ######
    ### Takehome point: LL kinda increases CV with increasing mean, but two major outliers keep it from being significant
                        # LMA is actualy significant, but significantly negative

quartz(width=7, height=3)
par(mfrow=c(1,3), mar=c(3.5,3.5,1,1), mgp=c(2.5,1,0), cex=1)
plot(cvLMA~mLMA, tmp2, ylab="LMA CV", xlab="mean LMA", pch=16)
abline(lm(cvLMA~mLMA, tmp2))
mtext(text = "p=0.054", side=3, line=-1.3, adj=.8)
plot(cvLL~mLL, tmp2[-which(tmp2$nLL<3 | tmp2$cvLL==0),], ylab="Leaf Lifespan CV", xlab="mean Leaf Lifespan (mo)", pch=16)
#summary(lm(cvLL~mLL, tmp2[-which(tmp2$nLL<3 | tmp2$cvLL==0 | tmp2$cvLL>.8 | tmp2$mLL>200),]))
mtext(text= "p=0.578", side=3, line=-1.3, adj=.8)
abline(lm(cvLL~mLL, tmp2[-which(tmp2$nLL<3 | tmp2$cvLL==0 | tmp2$cvLL>.8  | tmp2$mLL>200),]), lty=2)
#mtext(text= "(p=0.001)", side=3, line=-2.3, adj=.8)
plot(cvNmass~mNmass, tmp2, pch=16, xlab="mean Nmass", ylab="Nmass CV")
abline(lm(cvNmass~mNmass, tmp2), lty=2)
mtext(text= "p=0.259", side=3, line=-1.3, adj=.8)
