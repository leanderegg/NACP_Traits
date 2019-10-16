#############################################################
###      Analysis of NACP_TERRA trait and stand data
###         to assess within-species trait variation
#############################################################


# load required packages #
require(lme4)
require(lmerTest)
require(reshape)
require(lattice)
require(RColorBrewer)
require(mgcv)
require(dplyr)
require(car)
require(stringr)
require(stringi)
require(MuMIn)
require(lmodel2)

#set color palette
mypal <- brewer.pal(n=9, "Set1")
palette(mypal)


#__________________________________________________________________________________
######## BEGIN: Dataset creation and cleaning: #####################################
#__________________________________________________________________________________

############ . Load GLOPNET data ########
# these data were downloaded from https://www.nature.com/articles/nature02403 on 1/15/17
# species names were then cleaned and family data added using the {taxize} R package, and finalized by hand
# for script detailing this cleaning, please email Leander Anderegg at leanderegg@gmail.com
# 9 measurements were unable to be identified to family and were removed from the analysis
# all measurements ID'd to genus (e.g. "Pelea sp") were treated as a unique species
#note: this results in the averaging of 3 "Rubus sp" and 2 "Nephrolepsis sp", but all other unidentified species are sp1,sp2 etc.


setwd("/home/landeregg/traits")
setwd("/Users/leeanderegg/Dropbox/NACP_Traits/NACP_Traits_Rcode/Data_171220/")
LESclim <- read.csv("Wright2004_sitedata.csv")
#LES <- read.csv("/Users/leeanderegg/Dropbox/NACP_Traits/NACP_Traits_Rcode/LES_taxocleaning_033117_v4.csv", row.names=1)
LES <- read.csv("LES_taxocleaning_121917_v5.csv", row.names=1)
# still have 9 NAs for family that couldn't be taxonomically id'ed


LES$LMA <- 10^LES$log.LMA # create unlogged values for averaging to higher taxonomic levels
LES$LL <- 10^LES$log.LL
LES$Nmass <- 10^LES$log.Nmass
LES$Narea <- 10^LES$log.Narea
LES$Aarea <- 10^LES$log.Aarea
LES$Amass <- 10^LES$log.Amass
LES$SLA <- 1/LES$LMA # in GLOPNET, SLA and LMA are in different units (g/m2 vs cm2/g)

## add in climate
LES$MAT <- LESclim$MAT[match(LES$Dataset, LESclim$Dataset)]
LES$MAP <- LESclim$Rain[match(LES$Dataset, LESclim$Dataset)]
LES$VPD <- LESclim$VPD[match(LES$Dataset, LESclim$Dataset)]
LES$RAD <- LESclim$RAD[match(LES$Dataset, LESclim$Dataset)]
LES$PET <- LESclim$PET[match(LES$Dataset, LESclim$Dataset)]








######## . Load PNW data #########################
#these data were downloaded from  http://dx.doi.org/10.3334/ORNLDAAC/1292 on 1/23/16 by LDLA
# see Berner & Law 2016 for data description


## _____________Species ID lookup table to relate traits and biomass ______________________
speciesID <- read.csv("SpeciesNames_lookup_GOOD_031316.csv", header=T)
# this table relates the species abbreviations used in the NACP_TERRA_PNW forest_biomass and trait measurement datasets.


## ____________________________Soil Data______________________________
soil <- read.csv("NACP_TERRA_PNW_soil_cleaned.csv", header=T, na.strings = "-9999")

# Cleaning some data inconsistencies
soil[which(soil$Layer=="top" & soil$UpperDepth>5),]
# Plot 252 has two top layers. need to switch the second to 'bottom'
# Plot 86 has top and bottom layers flipped
## and something on the order of 9 plots don't have a 'top layer
soil$Layer[which(soil$PLOT_ID==252)] <- c("top", "middle","bottom")
soil$Layer[which(soil$PLOT_ID==86)] <- c("top", "bottom")
probs2 <- names(xtabs(~PLOT_ID, soil)[which(xtabs(~PLOT_ID, soil)==2)])
probs2names <- xtabs(~PLOT_ID, soil[which(soil$PLOT_ID %in% probs2 & soil$Layer != "top"),])
## 8 plots have 'middle' and 'bottom', but middle starts at 0. I'm going to just make all 'middles' into 'tops
soil$Layer[which(soil$PLOT_ID %in% names(probs2names[which(probs2names==2)]) & soil$Layer=="middle")] <- "top" # 8 plots have only 2 layers and one of them is not 'top'
# OK, that seemed to solve it. Everything has a 'top' layer, but below that is unknown. could be 'middle' could be 'bottom'...
# thus, we'll work with only the top layer for now.
soil.top <- soil[which(soil$Layer=="top"), ]



###_____________ Stand characteristics dataset______________________
# on 5.27.16 I deleted the units row in this dataset to make _v2 so that types import correctly
biomass <- read.csv("NACP_TERRA_PNW_forest_biomass_productivity_v2.csv", header= T,na.strings = "-9999" )

# merge soil characteristics and biomass data
biomass$soil_N <- soil.top$soil_N[match(biomass$PLOT_ID, soil.top$PLOT_ID)]
biomass$soil_pH <- soil.top$soil_pH[match(biomass$PLOT_ID, soil.top$PLOT_ID)]

# # calculate relative growth rate
# biomass$RGR <- log(biomass$AG_BIOMASS_TREE_TOTAL_AS_CARBON) - log(biomass$AG_BIOMASS_TREE_TOTAL_AS_CARBON - biomass$AG_PROD_TREE_TOTAL_AS_CARBON)
# # calculate an empirical biomass standardized growth rate
# growthmod <- lm(AG_PROD_TREE_TOTAL_AS_CARBON~AG_BIOMASS_TREE_TOTAL_AS_CARBON, data=biomass)
# growthmod2 <- gam(AG_PROD_TREE_TOTAL_AS_CARBON~s(AG_BIOMASS_TREE_TOTAL_AS_CARBON), data=biomass)
# growthmod3 <- gam(AG_PROD_TREE_TOTAL_AS_CARBON~s(log(AG_BIOMASS_TREE_TOTAL_AS_CARBON)), data=biomass)
# biomass$logAG_BIOMASS <- log(biomass$AG_BIOMASS_TREE_TOTAL_AS_CARBON)
# # there's actually a better relationship between log(biomass) than regular biomass
# biomass$BIOstGROWTH[which(biomass$AG_PROD_TREE_TOTAL_AS_CARBON>-5)] <- resid(growthmod)
# biomass$BIOstGROWTHgam[which(biomass$AG_PROD_TREE_TOTAL_AS_CARBON>-5)] <- resid(growthmod2)
# biomass$logBIOstGROWTHgam[which(biomass$AG_PROD_TREE_TOTAL_AS_CARBON>-5)] <- resid(growthmod3)
# 
# # add leaf mass fraction (leaf mass/total mass) and leaf allocation (leaf production/total production)
# biomass$LeafFrac <- biomass$AG_BIOMASS_TREE_FOLIAGE_AS_CARBON/biomass$AG_BIOMASS_TREE_TOTAL_AS_CARBON
# biomass$LeafAlloc <- biomass$AG_PROD_TREE_FOLIAGE_AS_CARBON/biomass$AG_PROD_TREE_TOTAL_AS_CARBON



### _______________________Traits dataset ____________________________________________
# newest version = NACP_TERRA_PNW_leaf_traits_v1_plusClim.csv
# older versions (on order they became obsolete): NACP_TERRA_PNW_leaf_trait.csv
# this is the traits data as downloaded from the ORNL repo, plus climate variables pulled by L Berner (see main text)
traits <- read.csv("NACP_TERRA_PNW_leaf_trait_v1_plusClim.csv", header=T, na.strings="-9999")
# make sure the correct columns are factors
fac.colstraits <- c(1,5:7,12:16)
for(j in fac.colstraits){ traits[,j] <- factor(traits[,j])}
# make sure the correct columns are numeric
num.colstraits <- c(8:11,17:27)
for (i in num.colstraits){traits[,i] <- as.numeric(as.character(traits[,i]))}

# ### back calculate some missing values:
# # calculate LEAF_CARBON for METOFIRE plots that have LEAF_CARBON_WT, but not elemental analysis
# traits$LEAF_CARBON[which(is.na(traits$LEAF_CARBON))] <- with(traits[which(is.na(traits$LEAF_CARBON)),], LEAF_CARBON_WT/LEAF_DRY_WT * 100)
# # length(which(is.na(traits$LEAF_CARBON))) # this should return 1

# add in Genus.species column, some clim columns and a column for species IDs used in biomass dataset
traits$GE.SP <- paste(traits$GENUS, traits$SPECIES, sep=".")
traits$SP.ID <- speciesID$bio.sp[match(traits$GE.SP, speciesID$traits.sp)]
#traits$MAT <-biomass$MAT_C[match(traits$PLOT_ID, biomass$PLOT_ID)]
#traits$MAP <- biomass$MAP[match(traits$PLOT_ID, biomass$PLOT_ID)] 
#sum soil moisture for total moisture content
traits$soilmoist.all.mm <- apply(traits[, c("soilmoist.lvl1.mm", "soilmoist.lvl2.mm", "soilmoist.lvl3.mm")],MARGIN = 1, FUN=sum)
# Id the measurements that come from a dominant species,,,
traits$FOREST_TYPE <- biomass$SPP_O1_ABBREV[match(traits$PLOT_ID, biomass$PLOT_ID)]
traits$dominance <- biomass$SPP_O1_BASAL_AREA_FRACTION[match(traits$PLOT_ID, biomass$PLOT_ID)]
# note: this only gives the BA fraction of the dominant spp, not of the spp for which the trait is measured. if FOREST_TYPE!=GE.SP, this is meaningless
# Approximate Stand Age
traits$ASA <- biomass$ASA[match(traits$PLOT_ID, biomass$PLOT_ID)]
traits$ELEVATION <- biomass$ELEVATION[match(traits$PLOT_ID, biomass$PLOT_ID)]
# traits$AG_TBIO <- biomass$AG_BIOMASS_TREE_TOTAL_AS_CARBON[match(traits$PLOT_ID, biomass$PLOT_ID)]
traits$AG_TGROWTH <- biomass$AG_PROD_TREE_TOTAL_AS_CARBON[match(traits$PLOT_ID, biomass$PLOT_ID)]
# traits$AG_WBIO <- biomass$AG_BIOMASS_TREE_WOOD_AS_CARBON[match(traits$PLOT_ID, biomass$PLOT_ID)]
# traits$AG_WGROWTH <- biomass$AG_PROD_TREE_WOOD_AS_CARBON[match(traits$PLOT_ID, biomass$PLOT_ID)]
# traits$AG_FBIO <- biomass$AG_BIOMASS_TREE_FOLIAGE_AS_CARBON[match(traits$PLOT_ID, biomass$PLOT_ID)]
# traits$AG_FGROWTH <- biomass$AG_PROD_TREE_FOLIAGE_AS_CARBON[match(traits$PLOT_ID, biomass$PLOT_ID)]
# traits$LeafFrac <- biomass$AG_BIOMASS_TREE_FOLIAGE_AS_CARBON[match(traits$PLOT_ID, biomass$PLOT_ID)]/biomass$AG_BIOMASS_TREE_TOTAL_AS_CARBON[match(traits$PLOT_ID, biomass$PLOT_ID)]

# traits$BIOST_TGROWTH <- biomass$BIOstGROWTH[match(traits$PLOT_ID, biomass$PLOT_ID)]
# traits$BIOST_TGROWTHgam <- biomass$BIOstGROWTHgam[match(traits$PLOT_ID, biomass$PLOT_ID)]
# traits$logBIOST_TGROWTHgam <- biomass$logBIOstGROWTHgam[match(traits$PLOT_ID, biomass$PLOT_ID)]
## relative growth rate (assume TBIO = M2 and TGROWTH=M2-M1), and TGROWTH is per year
# so I need to subtract TGROWTH from TBIO to get M1
# traits$RGR <- log(traits$AG_TBIO) - log(traits$AG_TBIO - traits$AG_TGROWTH)
# note: this also solves the problem of stand 98 where production if higher than biomass
# old way assuming forecasting rather than hindcasting
#traits$RGR <- log(traits$AG_TBIO + traits$AG_TGROWTH) - log(traits$AG_TBIO)


## add in soil characteristics
traits$soil_N <- soil.top$soil_N[match(traits$PLOT_ID, soil.top$PLOT_ID)]
traits$soil_C <- soil.top$soil_C[match(traits$PLOT_ID, soil.top$PLOT_ID)]
traits$soil_pH <- soil.top$soil_pH[match(traits$PLOT_ID, soil.top$PLOT_ID)]
traits$TotalSoilDepth <- soil.top$TotalDepth[match(traits$PLOT_ID, soil.top$PLOT_ID)]
traits$TopSoilDepth <- soil.top$LowerDepth[match(traits$PLOT_ID, soil.top$PLOT_ID)]
traits$Bulk_Density <- soil.top$Bulk_Density[match(traits$PLOT_ID, soil.top$PLOT_ID)]
traits$pSAND <- soil.top$pSAND[match(traits$PLOT_ID, soil.top$PLOT_ID)]
traits$pSILT <- soil.top$pSILT[match(traits$PLOT_ID, soil.top$PLOT_ID)]
traits$pCLAY <- soil.top$pCLAY[match(traits$PLOT_ID, soil.top$PLOT_ID)]


## calculate traits used in LES
traits$LMA <- with(traits, LEAF_DRY_WT/(LEAF_HSA*1/100^2)) # LES works in LMA rather than SLA
# also SLA is in gC rather than LEAF_DRY_WT, and in cm2 rather than m2 as in Wright 2004. so I need to recalculate based on leaf dry weight and hemispheric leaf area
# bunch of EPA plots are missing LEAF_DRY_WT and LEAF_HSA, so have to back calculate from LEAF_CARBON and SLA_HSA
traits$LMA[which(is.na(traits$LMA))] <- with(traits[which(is.na(traits$LMA)),], 1/(SLA_HSA * LEAF_CARBON/100 * 1/100^2))
# Note: recalculating LMA for all plots using this method yeilds LMA values w/ mean difference of 0.00044, so essentially rounding error
traits$SLA_drymass <- 1/traits$LMA # make a non-carbon SLA
traits$log.LMA <- with(traits, log(LMA,base = 10))
traits$LLmonths <- traits$LEAF_LIFE * 12 # LES works in months rather than years
traits$log.LL <- log(traits$LLmonths,base = 10)
traits$log.Nmass <- log(traits$LEAF_NITROGEN,base = 10)
traits$Narea <- traits$LEAF_NITROGEN/100 / traits$SLA_drymass #Narea = Nmass * Ml/Al (LMA)
traits$log.Narea <- log(traits$Narea,base = 10)
# make LMA_PSA but the pines don't have the conversion factor, so I'll just plug in 1.2 (PINPON=1.18, PINCON=1.29)
# traits$PSA_to_HSA[which(traits$SP.ID=="PINPON")] <- 1.18 # fill in missing PINPON values
# traits$PSA_to_HSA[which(traits$SP.ID=="PINCON")] <- 1.29 # fill in missing PINCON values
# traits$PSA_to_HSA[which(traits$GENUS=="Pinus" & is.na(traits$PSA_to_HSA))] <- 1.2
# traits$LMA_PSA <- traits$LMA * traits$PSA_to_HSA # * PSA_to_HSA is same as /HSA_to_PSA. and LMA is mass/area
# traits$log.LMA_PSA <- log(traits$LMA_PSA, base=10)

## make an RGR that only has values when that species makes up >50% of the stand BA
# traits$RGRdom <- traits$RGR
# traits$RGRdom[which(as.character(traits$FOREST_TYPE) != as.character(traits$SP.ID))] <- NA # get rid of all measurements that aren't on the dominant species
# traits$RGRdom[which(as.character(traits$FOREST_TYPE) == as.character(traits$SP.ID) & traits$dominance<50)] <- NA # define dominance as >50% of stand BA
# traits$stGrowthdom <- traits$BIOST_TGROWTHgam
# traits$stGrowthdom[which(as.character(traits$FOREST_TYPE) != as.character(traits$SP.ID))] <- NA
# traits$stGrowthdom[which(as.character(traits$FOREST_TYPE) == as.character(traits$SP.ID) & traits$dominance<50)] <- NA

## make restrictive RGR that only has values when spp makes up >90% of stand BA
# traits$RGRdom90 <- traits$RGR
# traits$RGRdom90[which(as.character(traits$FOREST_TYPE) != as.character(traits$SP.ID))] <- NA# get rid of all measurements that aren't on the dominant species
# traits$RGRdom90[which(as.character(traits$FOREST_TYPE) == as.character(traits$SP.ID) & traits$dominance<90)] <- NA
# traits$stGrowthdom90 <- traits$BIOST_TGROWTHgam
# traits$stGrowthdom90[which(as.character(traits$FOREST_TYPE) != as.character(traits$SP.ID))] <- NA
# traits$stGrowthdom90[which(as.character(traits$FOREST_TYPE) == as.character(traits$SP.ID) & traits$dominance<90)] <- NA
# note: these stGrowth calculations are based on all plots, so they contain the most information possible for the standardization process.

traits$FullSpecies <- paste(traits$GENUS, traits$SPECIES, sep=" ")

# make a PCA of climate variables and add them to the df
climpca.traits <- prcomp(traits[,c(grep("gy", colnames(traits)), which(colnames(traits) %in% c("soilmoist.lvl1.mm","soilmoist.all.mm")))],scale=T)
traits$climPC1 <- climpca.traits$x[,1]
traits$climPC2 <- climpca.traits$x[,2]
traits$climPC3 <- climpca.traits$x[,3]


## ______________ adding in family names to Traits ___________________
PNWfams <- read.csv("PNWfams_031717.csv")
traits$Family <- PNWfams$family[match(traits$FullSpecies, PNWfams$query)]




############### . calculate CWMs for PNW dataset ############

# note on nomenclature:
# TRAIT1 = species mean trait value for species 1
# wTRAIT1 = BA weighted species mean trait value
# TRAITp1 = plot/site mean value for the species species 1 - originally used s, but that was horrific. changed 06.14.17
# wTRAITp1 = BA weighted plot/site mean value for species 1

# total species means
#species.means <- traits %>% group_by(GENUS,SPECIES,SP.ID) %>% summarise(mLMA_HSA = mean(LMA_HSA, na.rm=T), mLLmonths = mean(LLmonths, na.rm=T), mCARBON = mean(LEAF_CARBON, na.rm=T), mNITROGEN = mean(LEAF_NITROGEN, na.rm=T), mCN = mean(LEAF_CN, na.rm=T), N = n(), nPlots = length(unique(PLOT_ID)), nProj = length(unique(PROJECT)), nEcoReg = length(unique(ECOREGION)), mMAP = mean(MAP), mMAT = mean(MAT), mlog.LMA = mean(log.LMA, na.rm=T), mlog.LL =mean(log.LL, na.rm=T), mlog.Nmass=mean(log.Nmass, na.rm=T), mlog.LMA_PSA = mean(log.LMA_PSA, na.rm=T), mRGR = mean(RGRdom, na.rm=T), mstGrowthdom = mean(stGrowthdom, na.rm=T)) 
## site means
# Note: mLMA = LMA_HSA, and mLMA_PSA = LMA_PSA, I've used LMA_HSA in full dataset creation
# also, mlog.Trait = mean of logged traits
# log.Trait = log of mean traits

spp.traits <- traits %>% group_by(SP.ID) %>% summarise(nsample = n(), SLA = mean(SLA_HSA, na.rm=T), nSLA = n()- length(which(is.na(SLA_HSA))), CN = mean(LEAF_CN, na.rm=T), nCN = n()- length(which(is.na(LEAF_CN)))
                                                       , LIFE = mean(LEAF_LIFE, na.rm=T), nLIFE=n()- length(which(is.na(LEAF_LIFE))), nplots =length(unique(PLOT_ID))
                                                       , mLMA = mean(LMA, na.rm=T), mLLmonths= mean(LLmonths, na.rm=T), mNmass = mean(LEAF_NITROGEN, na.rm=T), mNarea=mean(Narea, na.rm=T)
                                                       , mlog.LMA = mean(log.LMA, na.rm=T), mlog.LL = mean(log.LL, na.rm=T), mlog.Narea=mean(log.Narea, na.rm=T), mlog.Nmass=mean(log.Nmass, na.rm=T)
                                                       # , RGR = mean(RGRdom, na.rm=T), stGrowth = mean(stGrowthdom, na.rm=T)
                                                       , climPC1 = mean(climPC1, na.rm=T), climPC2 = mean(climPC2, na.rm=T), climPC3 = mean(climPC3, na.rm=T)
                                                       , soil_N= mean(soil_N, na.rm=T), soil_pH=mean(soil_pH, na.rm=T), ASA = mean(ASA, na.rm=T), LAI_O=mean(LAI_O, na.rm=T), AG_TGROWTH = mean(AG_TGROWTH, na.rm=T))
spp.traits$log.LMA <- log(spp.traits$mLMA, base=10)
spp.traits$log.LL <- log(spp.traits$mLLmonths, base=10)
spp.traits$log.Nmass <- log(spp.traits$mNmass, base=10)
spp.traits$log.Narea <- log(spp.traits$mNarea, base=10)


# calculate species means for each plot in which they occur
spp.plot.traits <- traits %>% group_by(SP.ID, PLOT_ID) %>% summarise(nsample = n(), mSLA = mean(SLA_HSA, na.rm=T), nSLA = n()- length(which(is.na(SLA_HSA))), mCN = mean(LEAF_CN, na.rm=T), nCN = n()- length(which(is.na(LEAF_CN)))
                                                                     , mLIFE = mean(LEAF_LIFE, na.rm=T), nLIFE=n()- length(which(is.na(LEAF_LIFE))), mNmass = mean(LEAF_NITROGEN, na.rm=T), nplots =length(unique(PLOT_ID))
                                                                     , mLMA = mean(LMA, na.rm=T), mLLmonths= mean(LLmonths, na.rm=T), mNarea=mean(Narea, na.rm=T)
                                                                     , mlog.LMA = mean(log.LMA, na.rm=T), mlog.LL = mean(log.LL, na.rm=T), mlog.Narea=mean(log.Narea, na.rm=T), mlog.Nmass=mean(log.Nmass, na.rm=T)
                                                                     #, RGR = unique(RGRdom), stGrowth = mean(stGrowthdom, na.rm=T)
                                                                     , climPC1 = unique(climPC1), climPC2 = unique(climPC2), climPC3 = unique(climPC3)
)
# create unique species-plot tag
spp.plot.traits$SP.PLOT <- paste(spp.plot.traits$SP.ID, spp.plot.traits$PLOT_ID, sep="-")


### make unique identifiers for matching biomass and spp.plot.traits rows
biomass$SP1.PLOT <- paste(biomass$SPP_O1_ABBREV,biomass$PLOT_ID, sep="-")
biomass$SP2.PLOT <- paste(biomass$SPP_O2_ABBREV,biomass$PLOT_ID, sep="-")
biomass$SP3.PLOT <- paste(biomass$SPP_O3_ABBREV,biomass$PLOT_ID, sep="-")
biomass$SP4.PLOT <- paste(biomass$SPP_O4_ABBREV,biomass$PLOT_ID, sep="-")



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
############ ++ CWMs based on species-plot mean trait values #######################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# naming convention:
#   "wTraitp1" = the plot mean trait value for species 1, weighted by its basal area fraction at that plot
#   "wTrait1" = the species mean trait value for species 1, weighted by its basal area fraction (for infilling species with no measurements at a plot)
#   "cw_Traitp" = the sum of basal area fraction-weighted traits, i.e. the community weighted mean (CWM) for a plot based on plot level trait measurements
#    the p indicates trait values were measured for that species, at that plot. no p in variable name indicates that species mean trait value was used.
#   "if" or "_if" indicates values with missing plot-level data that have been infilled with species mean data

#____________________________________ LMA  ______________________________________
biomass$wLMAp1 <- spp.plot.traits$mLMA[match(as.character(biomass$SP1.PLOT), as.character(spp.plot.traits$SP.PLOT))] * biomass$SPP_O1_BASAL_AREA_FRACTION/100
biomass$wLMAp2 <- spp.plot.traits$mLMA[match(biomass$SP2.PLOT, spp.plot.traits$SP.PLOT)] * biomass$SPP_O2_BASAL_AREA_FRACTION/100
biomass$wLMAp3 <- spp.plot.traits$mLMA[match(biomass$SP3.PLOT, spp.plot.traits$SP.PLOT)] * biomass$SPP_O3_BASAL_AREA_FRACTION/100
biomass$wLMAp4 <- spp.plot.traits$mLMA[match(biomass$SP4.PLOT, spp.plot.traits$SP.PLOT)] * biomass$SPP_O4_BASAL_AREA_FRACTION/100

biomass$cw_LMAp <- apply(biomass[,c("wLMAp1", "wLMAp2","wLMAp3","wLMAp4")],MARGIN=1,FUN=sum, na.rm=T)
# now need to remove values that got smoothed over due to na.rm=T in the apply function
biomass$cw_LMAp[which(is.na(biomass$wLMAp1))] <- NA # 86 plots lacking LMA data of SPP01
biomass$cw_LMAp[which(is.na(biomass$wLMAp2) & biomass$SPP_O2_BASAL_AREA_FRACTION>30)] <- NA
#18 sites lacking LMA data for SPP02 where SPP02 makes up >30 % of BA, but only 5 have SPP01
biomass$cw_LMAp[which(biomass$cw_LMAp==0)] <- NA
# 91 plots missing substantial LMA data
# length(which(is.na(biomass$cw_LMAp)))




#____________________________________ Leaf Life  ______________________________________
biomass$wLLp1 <- spp.plot.traits$mLLmonths[match(biomass$SP1.PLOT, spp.plot.traits$SP.PLOT)]* biomass$SPP_O1_BASAL_AREA_FRACTION/100
biomass$wLLp2 <- spp.plot.traits$mLLmonths[match(biomass$SP2.PLOT, spp.plot.traits$SP.PLOT)]* biomass$SPP_O2_BASAL_AREA_FRACTION/100
biomass$wLLp3 <- spp.plot.traits$mLLmonths[match(biomass$SP3.PLOT, spp.plot.traits$SP.PLOT)]* biomass$SPP_O3_BASAL_AREA_FRACTION/100
biomass$wLLp4 <- spp.plot.traits$mLLmonths[match(biomass$SP4.PLOT, spp.plot.traits$SP.PLOT)]* biomass$SPP_O4_BASAL_AREA_FRACTION/100

biomass$cw_LLp <- apply(biomass[,c("wLLp1", "wLLp2","wLLp3","wLLp4")],MARGIN=1,FUN=sum, na.rm=T)
# now need to remove values that got smoothed over due to na.rm=T in the apply function
biomass$cw_LLp[which(is.na(biomass$wLLp1))] <- NA # 94 sites lacking LeafLL of SPP01
biomass$cw_LLp[which(is.na(biomass$wLLp2) & biomass$SPP_O2_BASAL_AREA_FRACTION>30)] <- NA
# 22 sites lack SPP02 leaf life, 9 of them unique
length(which(is.na(biomass$cw_LLp))) # 103 plots sans significant LL

#____________________________________ Nmass  ______________________________________
biomass$wNmassp1 <- spp.plot.traits$mNmass[match(biomass$SP1.PLOT, spp.plot.traits$SP.PLOT)]* biomass$SPP_O1_BASAL_AREA_FRACTION/100
biomass$wNmassp2 <- spp.plot.traits$mNmass[match(biomass$SP2.PLOT, spp.plot.traits$SP.PLOT)]* biomass$SPP_O2_BASAL_AREA_FRACTION/100
biomass$wNmassp3 <- spp.plot.traits$mNmass[match(biomass$SP3.PLOT, spp.plot.traits$SP.PLOT)]* biomass$SPP_O3_BASAL_AREA_FRACTION/100
biomass$wNmassp4 <- spp.plot.traits$mNmass[match(biomass$SP4.PLOT, spp.plot.traits$SP.PLOT)]* biomass$SPP_O4_BASAL_AREA_FRACTION/100

biomass$cw_Nmassp <- apply(biomass[,c("wNmassp1", "wNmassp2","wNmassp3","wNmassp4")],MARGIN=1,FUN=sum, na.rm=T)
# now need to remove values that got smoothed over due to na.rm=T in the apply function
biomass$cw_Nmassp[which(is.na(biomass$wNmassp1))] <- NA # 98 sites lacking LeafNITROGEN of SPP01
biomass$cw_Nmassp[which(is.na(biomass$wNmassp2) & biomass$SPP_O2_BASAL_AREA_FRACTION>30)] <- NA
biomass$cw_Nmassp[which(biomass$cw_Nmassp==0)] <- NA
# 26 sites lack SPP02 leaf NITROGEN, 3 of them unique
# length(which(is.na(biomass$cw_Nmassp))) # 102 sites sans Nmass




#____________________________________ Narea  ______________________________________
biomass$wNareap1 <- spp.plot.traits$mNarea[match(biomass$SP1.PLOT, spp.plot.traits$SP.PLOT)]* biomass$SPP_O1_BASAL_AREA_FRACTION/100
biomass$wNareap2 <- spp.plot.traits$mNarea[match(biomass$SP2.PLOT, spp.plot.traits$SP.PLOT)]* biomass$SPP_O2_BASAL_AREA_FRACTION/100
biomass$wNareap3 <- spp.plot.traits$mNarea[match(biomass$SP3.PLOT, spp.plot.traits$SP.PLOT)]* biomass$SPP_O3_BASAL_AREA_FRACTION/100
biomass$wNareap4 <- spp.plot.traits$mNarea[match(biomass$SP4.PLOT, spp.plot.traits$SP.PLOT)]* biomass$SPP_O4_BASAL_AREA_FRACTION/100

biomass$cw_Nareap <- apply(biomass[,c("wNareap1", "wNareap2","wNareap3","wNareap4")],MARGIN=1,FUN=sum, na.rm=T)
# now need to remove values that got smoothed over due to na.rm=T in the apply function
biomass$cw_Nareap[which(is.na(biomass$wNareap1))] <- NA # 98 sites lacking LeafNarea of SPP01
biomass$cw_Nareap[which(is.na(biomass$wNareap2) & biomass$SPP_O2_BASAL_AREA_FRACTION>30)] <- NA
biomass$cw_Nareap[which(biomass$cw_Nareap==0)] <- NA
# 26 sites lack SPP02 leaf Narea, 3 of them unique



# logging community weighted trait means based on plot level trait measurements
biomass$log.cw_LMAp <- log(biomass$cw_LMAp, base=10)
biomass$log.cw_LLp <- log(biomass$cw_LLp, base=10)
biomass$log.cw_Nmassp <- log(biomass$cw_Nmassp, base=10)
biomass$log.cw_Nareap <- log(biomass$cw_Nareap, base=10)

###### add the climPCs to biomass
biomass$climPC1 <-spp.plot.traits$climPC1[match(biomass$SP1.PLOT, spp.plot.traits$SP.PLOT)]
biomass$climPC2 <-spp.plot.traits$climPC2[match(biomass$SP1.PLOT, spp.plot.traits$SP.PLOT)]



#_____________________________________________________________________________
###### ++ Infilling spp trait values with spp means when they're missing ###########
#_____________________________________________________________________________
### Infilling LMA of spp with <25% of BA with species mean values


###________________________________ Creating BA-weighted LMA values using spp means ________________________________
biomass$wLMA1 <- spp.traits$mLMA[match(biomass$SPP_O1_ABBREV, spp.traits$SP.ID)]*biomass$SPP_O1_BASAL_AREA_FRACTION/100
biomass$wLMA2 <- spp.traits$mLMA[match(biomass$SPP_O2_ABBREV, spp.traits$SP.ID)]*biomass$SPP_O2_BASAL_AREA_FRACTION/100
biomass$wLMA3 <- spp.traits$mLMA[match(biomass$SPP_O3_ABBREV, spp.traits$SP.ID)]*biomass$SPP_O3_BASAL_AREA_FRACTION/100
biomass$wLMA4 <- spp.traits$mLMA[match(biomass$SPP_O4_ABBREV, spp.traits$SP.ID)]*biomass$SPP_O4_BASAL_AREA_FRACTION/100

biomass$cw_LMA <- apply(biomass[,c("wLMA1", "wLMA2","wLMA3","wLMA4")],MARGIN=1,FUN=sum, na.rm=T)
biomass$cw_LMA[which(biomass$cw_LMA==0)] <- NA

wLMAif1 <- biomass$wLMAp1 # pull out the plot-level basal-area weighted values so I can infill missing values with the species mean
wLMAif1[which(is.na(wLMAif1))] <- biomass$wLMA1[which(is.na(wLMAif1))] # infill species means for missing SP1 traits
wLMAif2 <- biomass$wLMAp2
wLMAif2[which(is.na(wLMAif2))] <- biomass$wLMA2[which(is.na(wLMAif2))] # infill species means for missing SP2 traits
wLMAif3 <- biomass$wLMAp3
wLMAif3[which(is.na(wLMAif3))] <- biomass$wLMA3[which(is.na(wLMAif3))] # infill species means for missing SP3 traits
wLMAif4 <- biomass$wLMAp4
wLMAif4[which(is.na(wLMAif4))] <- biomass$wLMA4[which(is.na(wLMAif4))] # infill species means for missing SP4 traits
# calculate CWM, but with all missing trait values infilled with species means. Note: this infills some missing dominant species with species mean traits
cw_LMAp_if  <- apply(data.frame(wLMAif1, wLMAif2, wLMAif3, wLMAif4),MARGIN=1, FUN=sum, na.rm=T)
# get rid of plots that still have no data (na.rm=T means empty plots return 0)
cw_LMAp_if[which(cw_LMAp_if==0)] <- NA
# now remove plots where missing species represent >25% of basal area
cw_LMAp_if[which(is.na(biomass$wLMAp1) & biomass$SPP_O1_BASAL_AREA_FRACTION>25)] <- NA # 80 plots lacking LMA data of SPP01 & SPO1 has >25% of the basal area
cw_LMAp_if[which(is.na(biomass$wLMAp2) & biomass$SPP_O2_BASAL_AREA_FRACTION>25)] <- NA # 26 plots lacking LMA of SPP2 w>25% of BA, 6 plots that weren't killed above
cw_LMAp_if[which(is.na(biomass$wLMAp3) & biomass$SPP_O3_BASAL_AREA_FRACTION>25)] <- NA # 2 plots lacking a prevalent SPP3, only one novel
#93 plots are removed because they are missing trait data from a dominant species (species with >25% of stand basal area)

biomass$cw_LMAp_if <- cw_LMAp_if # add infilled CWM LMA values to biomass




###________________________________ Leaf Lifespan infilling with spp mean values ________________________________
biomass$wLL1 <- spp.traits$mLLmonths[match(biomass$SPP_O1_ABBREV, spp.traits$SP.ID)]*biomass$SPP_O1_BASAL_AREA_FRACTION/100
biomass$wLL2 <- spp.traits$mLLmonths[match(biomass$SPP_O2_ABBREV, spp.traits$SP.ID)]*biomass$SPP_O2_BASAL_AREA_FRACTION/100
biomass$wLL3 <- spp.traits$mLLmonths[match(biomass$SPP_O3_ABBREV, spp.traits$SP.ID)]*biomass$SPP_O3_BASAL_AREA_FRACTION/100
biomass$wLL4 <- spp.traits$mLLmonths[match(biomass$SPP_O4_ABBREV, spp.traits$SP.ID)]*biomass$SPP_O4_BASAL_AREA_FRACTION/100
biomass$cw_LL <- apply(biomass[,c("wLL1", "wLL2","wLL3","wLL4")],MARGIN=1,FUN=sum, na.rm=T)
biomass$cw_LL[which(biomass$cw_LL==0)] <- NA

wLLif1 <- biomass$wLLp1
wLLif1[which(is.na(wLLif1))] <- biomass$wLL1[which(is.na(wLLif1))]
wLLif2 <- biomass$wLLp2
wLLif2[which(is.na(wLLif2))] <- biomass$wLL2[which(is.na(wLLif2))]
wLLif3 <- biomass$wLLp3
wLLif3[which(is.na(wLLif3))] <- biomass$wLL3[which(is.na(wLLif3))]
wLLif4 <- biomass$wLLp4
wLLif4[which(is.na(wLLif4))] <- biomass$wLL4[which(is.na(wLLif4))]
cw_LLp_if  <- apply(data.frame(wLLif1, wLLif2, wLLif3, wLLif4),MARGIN=1, FUN=sum, na.rm=T)
# get rid of plots that still have no data
cw_LLp_if[which(cw_LLp_if==0)] <- NA
cw_LLp_if[which(is.na(biomass$wLLp1) & biomass$SPP_O1_BASAL_AREA_FRACTION>25)] <- NA # 88 plots lacking LL data of SPP01 & SPO1 has >25% of the basal area
cw_LLp_if[which(is.na(biomass$wLLp2) & biomass$SPP_O2_BASAL_AREA_FRACTION>25)] <- NA # 31 plots lacking LL of SPP2 w>25% of BA, 11 plots that weren't killed above
cw_LLp_if[which(is.na(biomass$wLLp3) & biomass$SPP_O3_BASAL_AREA_FRACTION>25)] <- NA # 2 plots lacking a prevalent SPP3, only one novel

biomass$cw_LLp_if <- cw_LLp_if


###________________________________ Nmass infilling with spp mean values ________________________________
biomass$wNmass1 <- spp.traits$mNmass[match(biomass$SPP_O1_ABBREV, spp.traits$SP.ID)]*biomass$SPP_O1_BASAL_AREA_FRACTION/100
biomass$wNmass2 <- spp.traits$mNmass[match(biomass$SPP_O2_ABBREV, spp.traits$SP.ID)]*biomass$SPP_O2_BASAL_AREA_FRACTION/100
biomass$wNmass3 <- spp.traits$mNmass[match(biomass$SPP_O3_ABBREV, spp.traits$SP.ID)]*biomass$SPP_O3_BASAL_AREA_FRACTION/100
biomass$wNmass4 <- spp.traits$mNmass[match(biomass$SPP_O4_ABBREV, spp.traits$SP.ID)]*biomass$SPP_O4_BASAL_AREA_FRACTION/100
biomass$cw_Nmass <- apply(biomass[,c("wNmass1", "wNmass2","wNmass3","wNmass4")],MARGIN=1,FUN=sum, na.rm=T)
biomass$cw_Nmass[which(biomass$cw_Nmass==0)] <- NA

wNmassif1 <- biomass$wNmassp1
wNmassif1[which(is.na(wNmassif1))] <- biomass$wNmass1[which(is.na(wNmassif1))]
wNmassif2 <- biomass$wNmassp2
wNmassif2[which(is.na(wNmassif2))] <- biomass$wNmass2[which(is.na(wNmassif2))]
wNmassif3 <- biomass$wNmassp3
wNmassif3[which(is.na(wNmassif3))] <- biomass$wNmass3[which(is.na(wNmassif3))]
wNmassif4 <- biomass$wNmassp4
wNmassif4[which(is.na(wNmassif4))] <- biomass$wNmass4[which(is.na(wNmassif4))]
cw_Nmassp_if  <- apply(data.frame(wNmassif1, wNmassif2, wNmassif3, wNmassif4),MARGIN=1, FUN=sum, na.rm=T)
# get rid of plots that still have no data
cw_Nmassp_if[which(cw_Nmassp_if==0)] <- NA
cw_Nmassp_if[which(is.na(biomass$wNmassp1) & biomass$SPP_O1_BASAL_AREA_FRACTION>25)] <- NA # 93 plots lacking Nmass data of SPP01 & SPO1 has >25% of the basal area
cw_Nmassp_if[which(is.na(biomass$wNmassp2) & biomass$SPP_O2_BASAL_AREA_FRACTION>25)] <- NA # 36 plots lacking Nmass of SPP2 w>25% of BA, 4 plots that weren't killed above
cw_Nmassp_if[which(is.na(biomass$wNmassp3) & biomass$SPP_O3_BASAL_AREA_FRACTION>25)] <- NA # 2 plots lacking a prevalent SPP3, only one novel

biomass$cw_Nmassp_if <- cw_Nmassp_if



##________________________________ Narea infilling with spp mean values ________________________________
biomass$wNarea1 <- spp.traits$mNarea[match(biomass$SPP_O1_ABBREV, spp.traits$SP.ID)]*biomass$SPP_O1_BASAL_AREA_FRACTION/100
biomass$wNarea2 <- spp.traits$mNarea[match(biomass$SPP_O2_ABBREV, spp.traits$SP.ID)]*biomass$SPP_O2_BASAL_AREA_FRACTION/100
biomass$wNarea3 <- spp.traits$mNarea[match(biomass$SPP_O3_ABBREV, spp.traits$SP.ID)]*biomass$SPP_O3_BASAL_AREA_FRACTION/100
biomass$wNarea4 <- spp.traits$mNarea[match(biomass$SPP_O4_ABBREV, spp.traits$SP.ID)]*biomass$SPP_O4_BASAL_AREA_FRACTION/100
biomass$cw_Narea <- apply(biomass[,c("wNarea1", "wNarea2","wNarea3","wNarea4")],MARGIN=1,FUN=sum, na.rm=T)
biomass$cw_Narea[which(biomass$cw_Narea==0)] <- NA

wNareaif1 <- biomass$wNareap1
wNareaif1[which(is.na(wNareaif1))] <- biomass$wNarea1[which(is.na(wNareaif1))]
wNareaif2 <- biomass$wNareap2
wNareaif2[which(is.na(wNareaif2))] <- biomass$wNarea2[which(is.na(wNareaif2))]
wNareaif3 <- biomass$wNareap3
wNareaif3[which(is.na(wNareaif3))] <- biomass$wNarea3[which(is.na(wNareaif3))]
wNareaif4 <- biomass$wNareap4
wNareaif4[which(is.na(wNareaif4))] <- biomass$wNarea4[which(is.na(wNareaif4))]
cw_Nareap_if  <- apply(data.frame(wNareaif1, wNareaif2, wNareaif3, wNareaif4),MARGIN=1, FUN=sum, na.rm=T)
# get rid of plots that still have no data
cw_Nareap_if[which(cw_Nareap_if==0)] <- NA
cw_Nareap_if[which(is.na(biomass$wNareap1) & biomass$SPP_O1_BASAL_AREA_FRACTION>25)] <- NA # 93 plots lacking Narea data of SPP01 & SPO1 has >25% of the basal area
cw_Nareap_if[which(is.na(biomass$wNareap2) & biomass$SPP_O2_BASAL_AREA_FRACTION>25)] <- NA # 36 plots lacking Narea of SPP2 w>25% of BA, 4 plots that weren't killed above
cw_Nareap_if[which(is.na(biomass$wNareap3) & biomass$SPP_O3_BASAL_AREA_FRACTION>25)] <- NA # 2 plots lacking a prevalent SPP3, only one novel

biomass$cw_Nareap_if <- cw_Nareap_if


biomass$log.cw_LMAp_if <- log(biomass$cw_LMAp_if, base=10)
biomass$log.cw_LLp_if <- log(biomass$cw_LLp_if, base=10)
biomass$log.cw_Nmassp_if <- log(biomass$cw_Nmassp_if, base=10)
biomass$log.cw_Nareap_if <- log(biomass$cw_Nareap_if, base=10)




#________________ End CWM trait calculation for PNW dataset _________________________







############### . Load supplemental w/in spp data ############
# This comes from Martin et al. 2015 on Coffea arabica
# plus a number of unpublished datasets collected primarily by LDL Anderegg
# these were combined into one dataframe. please email LDL Anderegg at leanderegg@gmail.com for additional metadata.


data.supp <- read.csv("OtherData_Combined_040117.csv", header=T, row.names=1)
# make climate columns that align with PNW and LES
data.supp$MAT <- NA
data.supp$MAP <- NA
data.supp$VPD <- NA




############### . Load additional w/in spp data from Mt Rainier############
# this includes data from 5 annual or perennial wildflowers, collected by ML Sethi in 2017
Rainier <- read.csv("Leaf Economic Spectrum_alltraits_v1_171211.csv", row.names=1)
Rainier$MAT <- NA
Rainier$MAP <- NA
Rainier$VPD <- NA
Rainier$Project <- "Rainier"






############### . Making Master Dataset ###################


###### combine PNW and glopnet dataset. ###
data1 <- traits %>% select(FullSpecies,log.LMA, log.LL, log.Nmass, log.Narea, GENUS, Family, tmean.gy.c,ppt.gy.mm,vpd.gy.max)
colnames(data1)[c(1,2,6,8:10)]<- c("Species", "log.LMA","Genus","MAT","MAP","VPD")
data1$Project <- rep("PACNW", times=nrow(data1))
data1$log.LL[which(data1$log.LL<1.2)] <- NA # all the deciduous species have '1yr' lifespan, but really that's wrong
data2 <- LES %>% filter(!is.na(Family)) %>% select(Species, log.LMA, log.LL, log.Nmass, log.Narea,Genus, Family, MAT, MAP, VPD)
# currently removes 9 records as of 12.19.17
data2$Project <- rep("GLOPNET", times=nrow(data2))

# select relevant columns from supplemental data
data3 <- data.supp %>% select(Species,log.LMA,log.LL, log.Nmass,log.Narea,Genus,Family,MAT, MAP, VPD, Project)

# select relevant columns from Mt. Rainier data
data4 <- Rainier %>% select(Species=Species_full, log.LMA, log.LL, log.Nmass, log.Narea, Genus, Family, MAT, MAP, VPD, Project)

## make combined dataset with all the taxonomically resolved species
data.all <- rbind(data1, data2, data3, data4) # new w/ Rainier: 4267 measurements total
data.all$Species <- factor(data.all$Species)
data.all$Genus <- factor(data.all$Genus)
data.all$Family <- factor(data.all$Family)
data.all$Project <- factor(data.all$Project)

#__________________________________________________________________________________
################ END: Dataset creation and cleaning #####################
#__________________________________________________________________________________





#__________________________________________________________________________________
################ Begin: Variance Decomposition and Figure 1 #####################
#__________________________________________________________________________________

###### . Full dataset, global variation #######################


## with all spp used for the hierarchical analysis
#logLMAvar <- lmer(log.LMA~ Project + (1|Family/Genus/Species), data.all)
logLMAvar <- lmer(log.LMA~ Project + (1|Family) + (1|Genus) + (1|Species), data.all)
logLLvar <- lmer(log.LL~ Project + (1|Family) + (1|Genus) + (1|Species), data.all[-which(data.all$Project=="CO"),])
logNmassvar <- lmer(log.Nmass~ Project + (1|Family) + (1|Genus) + (1|Species), data.all)
logNareavar <- lmer(log.Narea~ Project + (1|Family) + (1|Genus) + (1|Species), data.all)

# LMAvar <- lmer(10^log.LMA~ Project + (1|Family) + (1|Genus) + (1|Species), data.all)
# LLvar <- lmer(10^log.LL~ Project + (1|Family) + (1|Genus) + (1|Species), data.all[-which(data.all$Project=="CO"),])
# Nmassvar <- lmer(10^log.Nmass~ Project + (1|Family) + (1|Genus) + (1|Species), data.all)
# Nareavar <- lmer(10^log.Narea~ Project + (1|Family) + (1|Genus) + (1|Species), data.all)
# #logNmassvar <- lmer(log.Nmass~ Project + (1|Family) + (1|Genus) + (1|Species), data.all)
# # in absolute scale (months), the variation appears to be HUGE w/ind spp, then w/in genera, then between Families
# # all levels of this analysis are hugely non-normal, except Families aren't too bad
LMAvariance <- data.frame(VarCorr(logLMAvar))
LLvariance <- data.frame(VarCorr(logLLvar))
Nmassvariance <- data.frame(VarCorr(logNmassvar))
Nareavariance <- data.frame(VarCorr(logNareavar))
# raw traits
rLMAvariance <- data.frame(VarCorr(LMAvar))
rLLvariance <- data.frame(VarCorr(LLvar))
rNmassvariance <- data.frame(VarCorr(Nmassvar))
rNareavariance <- data.frame(VarCorr(Nareavar))

traitvars <- data.frame(LMAvariance[,4], LLvariance[,4], Nmassvariance[,4], Nareavariance[,4])
colnames(traitvars) <- c("logLMA", "logLL", "logNmass", "logNarea")
rownames(traitvars) <- c("BtwSpecies", "BtwGenera", "BtwFamilies", "WtinSpecies")

traitvars_scaled <- traitvars  
for(i in 1:ncol(traitvars)){
  traitvars_scaled[,i] <- traitvars[,i]/sum(traitvars[,i])
}


# rtraitvars <- data.frame(rLMAvariance[,4], rLLvariance[,4], rNmassvariance[,4], rNareavariance[,4])
# colnames(rtraitvars) <- c("LMA", "LL", "Nmass", "Narea")
# rownames(rtraitvars) <- c("BtwSpecies", "BtwGenera", "BtwFamilies", "WtinSpecies")

# rtraitvars_scaled <- rtraitvars  
# for(i in 1:ncol(rtraitvars)){
#   rtraitvars_scaled[,i] <- rtraitvars[,i]/sum(rtraitvars[,i])
# }
# rtraitvars_scaled2 <- rtraitvars_scaled[c(3,2,1,4),]
traitvars_scaled2 <- traitvars_scaled[c(3,2,1,4),]
# 








########## END: Variance Decomposition and Figure 1 ################







#__________________________________________________________________________________
################ Begin: Trait-trait scale dependence Analysis #####################
#__________________________________________________________________________________


########### . Create w/in spp, species means, genus means, and family means datasets ###############
# then creating a dataset of genera w/ >5 species, and families w/ >5 genera
# select species w/ > 5 measurements
commonspp <-  names(which(xtabs(~Species, data.all)>=5))
spp.data<- data.all %>% filter(Species %in% commonspp)
spp.data$Species <- factor(spp.data$Species)
spp.data$Genus <- factor(spp.data$Genus)
spp.data$Family <- factor(spp.data$Family)


########## Species level trait averages ______________________________________________

allspp <- data.all %>% group_by(Species, Genus, Family) %>% summarise( lslog.LL = mean(log.LL, na.rm=T), lslog.LMA = mean(log.LMA, na.rm=T), lslog.Nmass = mean(log.Nmass, na.rm=T), lslog.Narea = mean(log.Narea, na.rm=T)
                                                                       ,slog.LL = log(mean(10^log.LL, na.rm=T),base=10), slog.LMA = log(mean(10^log.LMA, na.rm=T),base=10), slog.Nmass = log(mean(10^log.Nmass, na.rm=T),base=10), slog.Narea = log(mean(10^log.Narea, na.rm=T),base=10)
                                                                       ,MAT = mean(MAT, na.rm=T), MAP=mean(MAP, na.rm=T), VPD=mean(VPD, na.rm=T) )
colnames(allspp) <- gsub("slog", "log", colnames(allspp))




#### Genus means ___________________________________________________________________
allgen <- allspp %>% group_by(Genus)  %>% summarise(Fam = sort(unique(Family))[1], lglog.LL = mean(llog.LL, na.rm=T),lglog.LMA = mean(llog.LMA, na.rm=T),lglog.Nmass = mean(llog.Nmass, na.rm=T), lglog.Narea = mean(llog.Narea, na.rm=T)
                                                    ,glog.LL = log(mean(10^log.LL, na.rm=T),base=10), glog.LMA = log(mean(10^log.LMA, na.rm=T),base=10), glog.Nmass = log(mean(10^log.Nmass, na.rm=T),base=10), glog.Narea = log(mean(10^log.Narea, na.rm=T),base=10), nspp = n() 
                                                    ,MAT = mean(MAT, na.rm=T), MAP=mean(MAP, na.rm=T), VPD=mean(VPD, na.rm=T) )
# taking mean of raw values or logged values doesn't matter all that much yet
colnames(allgen) <- gsub("glog", "log", colnames(allgen))
colnames(allgen)[2] <- "Family"
# 939 genera from 211 families
# allgen <- allspp %>% group_by(Genus)  %>% summarise(Fam = sort(unique(Family))[1], Needle.Broad = unique(na.omit(NeedleBroad))[1], glog.LL = mean(log.LL, na.rm=T), glog.LMA = mean(log.LMA, na.rm=T), glog.Nmass = mean(log.Nmass, na.rm=T),glog.Narea = mean(log.Narea, na.rm=T)
#                                                     ,rglog.LL = log(mean(10^rlog.LL, na.rm=T),base=10), rglog.LMA = log(mean(10^rlog.LMA, na.rm=T),base=10), rglog.Nmass = log(mean(10^rlog.Nmass, na.rm=T),base=10), nspp = n() )



#### Family Means _______________________________________________________________________
allfam <- allgen %>% group_by(Family) %>% summarise(lflog.LL = mean(llog.LL, na.rm=T), lflog.LMA = mean(llog.LMA, na.rm=T), lflog.Nmass = mean(llog.Nmass, na.rm=T),lflog.Narea = mean(llog.Narea, na.rm=T)
                                                    ,flog.LL = log(mean(10^log.LL, na.rm=T),base=10), flog.LMA = log(mean(10^log.LMA, na.rm=T),base=10), flog.Nmass = log(mean(10^log.Nmass, na.rm=T),base=10), flog.Narea = log(mean(10^log.Narea, na.rm=T), base=10)
                                                    , tnspp = sum(nspp), ngen=n() )
colnames(allfam) <- gsub("flog", "log", colnames(allfam))

#### Family Means for fams w/ > 3 species
fam.data <- allfam # new: 211 families (up from 189)
fam.dataclean <- allfam[which(allfam$tnspp>2),] # 101 families, up from 97 families




##### dataset of only genera w/ >5 species _____________________________________________
# for w/in genus analysis
gen.data <- allspp[which(allspp$Genus %in% names(which(xtabs(~Genus, allspp)>=5))),]
# 750 measurements of 73 genera.
gen.data$Genus <- factor(gen.data$Genus) # get rid of unused levels


##### Spp. in Family level data, for families w/ >5 species in them ________________________
sppinfam.data <- allspp[which(allspp$Family %in% names(which(xtabs(~Family,allspp)>=5))),]
sppinfam.data$Family <- factor(sppinfam.data$Family) # get rid of unused levels



#### gen in Family level data, for families w/ >5 genera in them _________________________
# 556 measurements of 40 Families , added 2 Families from last batch, and probably quite a few obs (old n_LMA.N = 459)
# 04.01.17 - 684 obs of 50 Families. again, big win from taxonomic cleaning!
geninfam.data <- allgen[which(allgen$Family %in% names(which(xtabs(~Family,allgen)>=5))),]
geninfam.data$Family <- factor(geninfam.data$Family) # get rid of unused levels
geninfam.data$Genus <- factor(geninfam.data$Genus)




#_______________________________________________________________________

############ . Functions for fitting SMA regressions and null models ################
#_______________________________________________________________________



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
  pt(-abs(t), n-2)*2
}




############ . Fitting SMA slopes and null models ##############################
niters <- 10000 # number of null model iterations to use
niters.spinfam <- 10 # number of null model iterations to use for species-in-family analysis (not presented in ms) 

dirname <- "20171220_results"
dir.create(dirname)

############# + LMA vs Nmass ###########

###_______________ Species level analysis _______________________________________________


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

###_______________ Species w/in Genus level analysis _______________________________________________
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


###_______________ Species w./in Fam level analysis  (not in ms)_______________________________________________
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
    nullbounds <- fit.null(xvar='log.LMA', yvar="log.Nmass", observed = dataz, nulldata = allfam, nits = niters.spinfam)  
    sppinfam.results[i, 11:16] <- nullbounds
    sppinfam.results[i,17] <- test.sig(x=sppinfam.results$Rho[i], test=nullbounds)
  }
  
}
proc.time()-t0

#_______________ Genus w/in Family level analysis _______________________________________________
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

###_______________ Family level analysis _______________________________________________

# # currently just working with LES until I combine the PACNW dataset into this.
# famcorr.all <- cor(x=fam.data[,c("log.LL","log.LMA","log.Nmass")], use = "pairwise.complete.obs")
# famcorr.clean <- cor(x=fam.dataclean[,c("log.LL","log.LMA","log.Nmass")], use = "pairwise.complete.obs")

fam.res_LMA.N <- c("fam.all", fit.MAR(xvar='log.LMA',yvar="log.Nmass",data=fam.data),rep(NA, times=7),"Fam") 
names(fam.res_LMA.N) <- c("Taxo.Unit","Int","Slope","Slope.lci","Slope.uci","Rho","r.sq", "n","varLMA","varNmass", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig" ,"Type")
fam.resclean_LMA.N <- c("fam.clean", fit.MAR(xvar='log.LMA',yvar="log.Nmass",data=fam.dataclean),rep(NA, times=7), "Fam.clean")
names(fam.resclean_LMA.N) <- c("Taxo.Unit","Int","Slope","Slope.lci","Slope.uci","Rho","r.sq", "n","varLMA","varNmass", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Type")
global_LMA.N <- c("global", fit.MAR(xvar='log.LMA',yvar="log.Nmass",data=allspp),rep(NA, times=7), "global")
names(global_LMA.N) <- c("Taxo.Unit","Int","Slope","Slope.lci","Slope.uci","Rho","r.sq", "n","varLMA","varNmass", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Type")


###_______________ Combining into 1 df _______________________________________________
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

write.csv(all.results.LMAN, file = paste0("./",dirname,"/SMA_Results_LMAvNmass.csv"))




# still to run
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


