#######################################################################################
######## Importing and Combining Supplementary within-species trait data ##############
########## Need to make 'WithinSpecies_Trait_Data_040117.csv' for Cd-Anderegg2018_Manuscript_Analysis_v2.R
#######################################################################################




### arabidopsis data from Blonder et al. 2016
arab <- read.csv("/Users/leeanderegg/Dropbox/NACP_Traits/Intra-data/Blonder2015Arabidopsis.csv",header=T)
colnames(arab)<- c("Type","ID","Genotype","LL","LMA","Nmass","Amass","VD","LDMC")
arab$SLA <- 1/arab$LMA
arab$Aarea <- arab$Amass * arab$LMA
arab$Narea <- arab$Nmass/100 * arab$LMA
head(arab)
arab$LLmonths <- arab$LL/30
arab$log.LL <- log(arab$LLmonths, base=10)
arab$log.LMA <- log(arab$LMA, base=10)
arab$log.Narea <- log(arab$Narea, base=10)



resN <- resid(lm(log.Nmass~log.LMA, LES))

plot(resN~log.LL, LES[-which(is.na(LES$log.Nmass) | is.na(LES$log.LMA)),])

ITresN <- resid(lm(log.Nmass~log.LMA, traits.common5))
plot(ITresN~log.LL, traits.common5[-which(is.na(traits.common5$log.Nmass) | is.na(traits.common5$log.LMA)),], col=SP.ID)


#### examining Viburnum dataset:
vib <- read.table("/Users/leeanderegg/Desktop/alltrait_raw.txt", header = T)
vib$LLmonths <- vib$LLS / 4
vib$log.LL <- log(vib$LLmonths, base=10)
vib$log.LMA <- log(vib$LMA, base=10)
vib$log.Nmass <- log(vib$X.N, base=10)


##### Helianthus common garden
hel <- read.csv("/Users/leeanderegg/Desktop/Mason2015_Helianthus/Helianthus_data_Mason2015_commongarden.csv", header=T,na.strings = ".")

hel <- hel[which(hel$Species != "H. divaricatus"),]
hel$Species <- factor(hel$Species)
plot(log(LL_days)~log(LMA), hel, col=Species, pch=as.numeric(Species))
for(i in levels(hel$Species)){
  plot.MAR(xvar = "LMA", yvar = "LL_days",data = hel[which(hel$Species==i),], linecol = "black")
}




################### IMPORT and COMBINE ################
###### combining all the other datasets into something I can add to LES and traits for Variance Decomp and covariance analysis

## arabidposis dataset from Blonder et al. 2015, not going to add because growth chamber

##### Coffee data from Martin et al. 2016
coffee <- read.csv("/Users/leeanderegg/Dropbox/NACP_Traits/Intra-data/Martin2016Coffee.csv",header=T)
coffee$Narea <- coffee$leafn/100 * coffee$lma_g_m2
coffee$log.LMA <- log(coffee$lma_g_m2, base=10)
coffee$log.Narea <- log(coffee$Narea, base=10)
coffee$Species <- "Coffea arabica"
# averaged to individual rather than branch
coffee.ind <- coffee %>% group_by(Species, tree_unique, site) %>% summarise(LMA = mean(lma_g_m2, na.rm=T), Narea = mean(Narea, na.rm=T), Nmass=mean(leafn, na.rm=T))
coffee.ind$log.LMA <- log(coffee.ind$LMA, base=10)
coffee.ind$log.Narea <- log(coffee.ind$Narea, base=10)
coffee.ind$log.Nmass <- log(coffee.ind$Nmass, base=10)
coffee.ind$LL <- NA
coffee.ind$log.LL <- NA
coffee.ind$Genus <- "Coffea"
coffee.ind$Family <- "Rubiacaea"
coffee.ind$Project <- "coffee"
colnames(coffee.ind)[2:3] <- c("Indiv", "Site")


####### LDLA unpublished eucaluptus data from Tasmania
#already averaged to individual
# Bring in Tasmanian Euc data, preliminary for OVAT and OBLI as pf 03.29.17
Euc <- read.csv("/Users/leeanderegg/Dropbox/NACP_Traits/Intra-data/LDLA_TasiEucs_initial032917.csv", header=T, row.names=1)
colnames(Euc)[29] <- "Nmass"
Euc$LMAold <- Euc$LMA
Euc$LMA <- Euc$LMAold * 10000
Euc$SLA <- 1/Euc$LMA
Euc$log.LMA <- log(Euc$LMA, base=10)
Euc$log.Nmass <- log(Euc$Nmass, base=10) #in the GLOPNET dataset, Nmass is in log.% not log g/g
Euc$Narea <- Euc$Nmass/100 * Euc$LMA
Euc$log.Narea <- log(Euc$Narea, base=10)
Euc$SP.ID <- Euc$Species
Euc$Species <- "Eucalyptus obliqua"
Euc$Species[which(Euc$SP.ID=="OVAT")] <- "Eucalyptus ovata"
write.csv(Euc, "Eucalyptus_LeafTrait_Data.csv")

euc.ind <- Euc %>% select(Species, Treetag.x, Site, LMA, Narea, Nmass, log.LMA, log.Narea, log.Nmass)
euc.ind$LL <- NA
euc.ind$log.LL <- NA
euc.ind$Genus <- "Eucalyptus"
euc.ind$Family <- "Myrtacaea"
euc.ind$Project <- "Eucs"
colnames(euc.ind)[2]<- "Indiv"

######## WRL Anderegg unpublished Quercus gambelii data from CO
quga_branch <- read.csv("/Users/leeanderegg/Dropbox/NACP_Traits/Intra-data/QUGA_Dryweights_wLL.csv", header=T)
quga_branch$bLMA <- quga_branch$LEAFMASS/quga_branch$AREA * 1000000
quga <- quga_branch %>% group_by (SPECIES,ELEV,STAND,TREE) %>% summarise(LMA=mean(bLMA),LLmonths=mean(LLmonths)) 
quga$log.LL<- log(quga$LLmonths, base=10)
# quga leaf area is in mm2
quga$SLA <- 1/quga$LMA
quga$log.LMA <- log(quga$LMA, base=10)
quga$Narea <- NA
quga$Nmass <- NA
quga$log.Narea <- NA
quga$log.Nmass <- NA
quga$Species <- "Quercus gambelii"
quga$treetag <- with(quga,paste(SPECIES,ELEV,STAND,TREE, sep="-"))
quga$standtag <- with(quga, paste(SPECIES,ELEV,STAND, sep="-"))
quga.ind <- data.frame(quga) %>% select(Species,treetag,standtag, LMA, Narea, Nmass, log.LMA, log.Narea,log.Nmass, LLmonths, log.LL) 
quga.ind$Genus <- "Quercus"
quga.ind$Family <- "Fagaceae"
quga.ind$Project <- "CO"
colnames(quga.ind)[c(2,3,10)] <- c("Indiv","Site","LL")
write.csv(quga.ind, "Quercus_gambelii_LeafTrait_Data_20180220.csv")

## average to the stand level for trait covariation analysis
mquga <- quga %>% group_by(Species, standtag) %>% summarise(LMA = mean(LMA), Narea = mean(Narea), Nmass=mean(Nmass), LL = mean(LLmonths))


####### LDLA & JHRL 2015 data for Populus Tremuloides in CO
COpotr <- read.csv("/Users/leeanderegg/Dropbox/NACP_Traits/Intra-data/LDLA_JHRL_Traits_all_final_10_25_15.csv", header=T)
COpotrll <- read.csv("/Users/leeanderegg/Dropbox/NACP_Traits/Intra-data/CO_POTR-LL_033017.csv", header=T)
potr <- COpotr %>% filter(Species=="POTR") %>% select(Species, Elev, Plot, Tree, TreeTag,PlotUn, Lat, Long, ElevPlot, SLA)
potr$SP.ID <- potr$Species
potr$Species <- "Populus tremuloides"
## SLA is in cm2 rather than m2
potr$SLAcm2 <- potr$SLA
potr$SLA <- potr$SLAcm2 / 10000
potr$LMA <- 1/potr$SLA
potr$log.LMA <- log(potr$LMA, base=10)
# these are the LL lifespans from the AVHRR pixels.
potr$LLdays <- COpotrll$mean_growing[match(potr$PlotUn, COpotrll$Plot.Name)]
# but average elevation LL might be the best I can really do...
# also, I believe the average AVHRR growing season length of low elevation (195 days) and mid elevation (176 days),
# but I don't believe the H elevation (179 days) because it's probably mostly conifer pixels... So I'll remove them
potr$LLdayselev <- 195 # the average L
potr$LLdayselev[which(potr$Elev=="M")] <- 176
potr$LLdayselev[which(potr$Elev=="H")] <- NA
potr$LLmonths <- potr$LLdays/30
potr$LLmonthselev <- potr$LLdayselev/30
potr$log.LL <- log(potr$LLmonths, base=10)
potr$log.LLelev <- log(potr$LLmonthselev, base=10)
potr$Narea <- NA
potr$Nmass <- NA
potr$log.Narea <- NA
potr$log.Nmass <- NA
potr.ind <- potr %>% select(Species, TreeTag, PlotUn, LMA, Narea, Nmass, log.LMA, log.Narea, log.Nmass, LLmonthselev,log.LLelev)
potr.ind$Genus <- "Populus"
potr.ind$Family <- "Salicaceae"
potr.ind$Project <- "CO"
colnames(potr.ind)[c(2,3,10,11)]<- c("Indiv","Site","LL","log.LL")
write.csv(potr.ind, "Populus_tremuloides_LeafTrait_Data_20180202.csv")
# stand level averages
mpotr <- potr %>% group_by(PlotUn) %>% summarise(LMA=mean(LMA, na.rm=T), LLmonths=mean(LLmonthselev))
mpotr$log.LMA <- log(mpotr$LMA, base=10)
mpotr$log.LL <- log(mpotr$LLmonths, base=10)




######## Combine individual level data for all supplemental datasets ########
suppdata <- rbind(data.frame(coffee.ind), data.frame(euc.ind), data.frame(quga.ind), data.frame(potr.ind))
#write.csv(suppdata,"/Users/leeanderegg/Dropbox/NACP_Traits/Intra-data/OtherData_Combined_033017.csv")
write.csv(suppdata,"/Users/leeanderegg/Dropbox/NACP_Traits/Intra-data/OtherData_Combined_040117.csv")


### Renamed WithinSpecies_Trait_Data_040117.csv ######
