#############################################################
###      Analysis of NACP_TERRA trait and stand data
###         to assess within-species trait variation
###
###         data downloaded from: http://dx.doi.org/10.3334/ORNLDAAC/1292 
###            on 01/23/16 by LDLA
#############################################################

# load required packages #
require(lme4)
require(lmerTest)
require(reshape)
require(lattice)
require(RColorBrewer)
require(mgcv)
require(dplyr)
require(rgl) # the plot3d function
require(scatterplot3d) # scatterplot3d() function
require(car)
require(ggplot2)
require(stringr)
require(stringi)
# optimized pairs function
source("/Users/leeanderegg/Desktop/ZuurMixedModelling/AllRCode/HighstatLibV6.R")
# ggplot funcitons (sets default w/ no background, also greats multiplot() function)
source("/Users/leeanderegg/Desktop/R functions, general/ggplot_helpers.R")
# error_bars function
source("/Users/leeanderegg/Desktop/R functions, general/error_bars.R")
# standard error function
sterr <- function(x,...){
  sd(x, na.rm=T)/ length(na.omit(x))
}

# set color palette
mypal <- brewer.pal(n=12, "Set3")
mypal <- brewer.pal(n=9, "Set1")

palette(mypal)



#write.csv(biomass, "PACNW_Biomass_plus_traits_111916.csv")


####### LOAD DATA ########

## public TRY traits 
# TRY <- read.csv("/Users/leeanderegg/Dropbox/NACP_Traits/Try_public_LES traits/TRY_publicLES_011717.csv")
# TRY$DataName <- gsub(TRY$DataName, pattern=" ",replacement = ".")
# #TRYnew <- TRY[-which(is.na(TRY$StdValue)),]
# TRYtest <- TRY %>% select(ObservationID, DatasetID, AccSpeciesName, DataName, OrigValueStr) %>% filter(!duplicated(.))
#     # there are 102 rows that are somehow duplicated but have different DataIDs... grrrr
#     #length(which(duplicated(TRYtest)))
# # probs <- stri_trans_nfc(TRYtest$OrigValueStr)
#   # some weird characters in lat,longs in OrigValueStr
# TRYtest$OrigValueStr <- stri_trans_nfc(TRYtest$OrigValueStr)
# # TRYtest$OrigValueStr[403688:403689]
# 
# TRYwide <- spread(data = TRYtest, value = OrigValueStr, key =  DataName, convert=T)
# cols <- c(grep('Specific', colnames(TRYwide)), grep('Leaf.lifespan', colnames(TRYwide)), grep('Leaf.nitrogen', colnames(TRYwide)), grep('Specific', colnames(TRYwide)), grep('Leaf.area.index', colnames(TRYwide)),grep('Sun.vers', colnames(TRYwide)),grep('Leaf.area.index', colnames(TRYwide)), grep('Family', colnames(TRYwide), grep("Genus", colnames(TRYwide))) )
# 
# TRYred <- TRYwide[,c(1:3, cols)]
# colnames(TRYred) <- c( "ObservationID","DatasetID","AccSpeciesName", "SLA","SLA_HSA","Leaf_Life","Narea","Nmass","SLA1","SLA_HSA1","standLAI","SunShade","LAI1","Family","Family.APG" )
# TRYred$LMA <- 1/(TRYred$SLA)
# TRYred$log.LL <- log(as.numeric(TRYred$Leaf_Life), base=10)
# TRYred$log.Nmass <- log(as.numeric(TRYred$Nmass), base=10)
# TRYred$log.LMA <- log(as.numeric(TRYred$LMA), base=10)


# vars <- unique(TRY$DataName)[c(1:7,9:11,13,14,17,19,20,22,217:220,233,239,252)]
# TRYred <- TRY %>% filter(DataName %in% vars) %>% select(DatasetID, Dataset, AccSpeciesID, AccSpeciesName,ObservationID,TraitName,DataID,DataName,OrigValueStr,Replicates,StdValue ) %>% arrange(ObservationID, DataID)
# TRYlong <- TRYred %>% select(DatasetID, Dataset, AccSpeciesID, AccSpeciesName,ObservationID,DataName,OrigValueStr)
# TRYtest <- TRYlong %>% reshape(idvar=1:5, timevar="DataName", direction="wide")
# TRYtest <- spread(TRYlong, key=DataName, value=OrigValueStr)

# TRYsla <- TRY %>% filter(TraitID=="11") %>% arrange(AccSpeciesName, DatasetID, ObservationID)
# TRYnmass <- TRY %>% filter(TraitID=="14") %>% arrange(AccSpeciesName, DatasetID, ObservationID)
# TRYll <- TRY %>% filter(TraitID=="12") %>% arrange(AccSpeciesName, DatasetID, ObservationID)

#### NOTE: 
TRY <- read.csv("/Users/leeanderegg/Dropbox/NACP_Traits/Try_public_LES traits/try_traits_new_031417.csv", header=T)

tmp <- TRY[which(TRY$LL>0, TRY$SLA>0),]
tmp1 <- tmp[which(tmp$SpeciesName %in% names(which(xtabs(~SpeciesName, tmp)>2))),]
tmp1$SpeciesName <- factor(tmp1$SpeciesName)

## Wright et al 2004 LES data
    # initial import with taxonomic association. still have 500 family NAs
# LES <- read.csv("Writght2004_LESdata.csv", header=T)#[,-c(22:25)]
# LES$GE.SP <- str_replace(LES$Species, " ",".") # get a genus.species column to compare with traits
# LES$Genus <- matrix(unlist(str_split(string=LES$Species, pattern=" ", n=2)), byrow=T, ncol=2)[,1]
#   # note, things are log10, not ln in this dataset!
#   # also note: LMA is g/m2, and SLA in traits is cm2/gC
#   # -> so to convert I have to recalculate SLA w/ leaf dry wt, then multiply by 1/100^2, then take 1/SLAnew
#   # also turns out that EPA doesn't report LEAF_DRY_WT or LEAF_HSA. so I'll have to back calculate it from SLA_HSA and %C
# 
# ## getting family 
# family.itis1 <- tax_name(LES$Species[1:500],get="family",db="itis")
# newfams1 <- tax_name(LES$Species[which(is.na(family.itis1$family))], get="family", db="ncbi")
# families500 <- family.itis1$family
# families500[which(is.na(family.itis1$family))] <- newfams1$family
#   #  so 99 out of the first 500 spp don't have a family in either database...
# 
# # middle 1000
# family.itis2 <- tax_name(LES$Species[501:1500],get="family",db="itis")
# newfams2 <- tax_name(LES$Species[501:1500][which(is.na(family.itis2$family))], get="family", db="ncbi")
# families1500 <- family.itis2$family
# families1500[which(is.na(family.itis2$family))] <- newfams2$family
# 
# 
# # l1501 to end
# family.itis3 <- tax_name(LES$Species[1501:2548],get="family",db="itis")
# newfams3 <- tax_name(LES$Species[1501:2548][which(is.na(family.itis3$family))], get="family", db="ncbi")
# families2500 <- family.itis3$family
# families2500[which(is.na(family.itis3$family))] <- newfams3$family
# 
# families <- c(families500, families1500, families2500)
# LES$family <- families
# 
# write.csv(LES, "Writght2004_LESdata_taxonomy012517.csv")


###### updated dataframe w/ taxonomy #####
LES <- read.csv("Writght2004_LESdata_taxonomy012517.csv", header=T)[,-1]
#### LES site climate
LESclim <- read.csv("/Users/leeanderegg/Dropbox/NACP_Traits/Wright2004_sitedata.csv")

##### Code for fixing problem spp that result in problem genera (genera w/ multiple families)
# some Acers classified in old family - Aceraceae
LES$Family[which(LES$Genus=="Acer" & LES$Family != "Sapindaceae")] <- "Sapindaceae"                                                   
# one Atriplex in old family - Chenopodiaceae
LES$Family[which(LES$Genus=="Atriplex")] <- "Amaranthaceae"
# misnamed Betula and Lycaenidae
LES[which(LES$Species=="Betula pumilla"), c("Species","GE.SP")] <- t(gsub(pattern = "ll", replacement = "l", x = t(LES[which(LES$Species=="Betula pumilla"), c("Species","GE.SP")]))) 
LES$Family[which(LES$Genus=="Betula")] <- "Betulaceae"
# misnamed Carex Asteraceae (Carex dimorphotheca also a synonymn for Carex stenophylla subsp. stenophylloides)
LES$Family[which(LES$Genus=="Carex")] <- "Cyperaceae"
# Diospyros mislabeld Fabaceae
LES$Family[which(LES$Genus=="Diospyros")] <- "Ebenaceae"
# Eugenia mislabeld a lot
LES$Family[which(LES$Genus=="Eugenia")] <- "Myrtaceae"
# Ficus mislabeled Agridae
LES$Family[which(LES$Genus=="Ficus")] <- "Moraceae"
# Grevillea often mislabeled
LES$Family[which(LES$Genus=="Grevillea")] <- "Proteaceae"
# Hakea often mislabeled
LES$Family[which(LES$Genus=="Hakea")] <- "Proteaceae"
# Hibbertia usually mislabeled
LES$Family[which(LES$Genus=="Hibbertia")] <- "Dilleniaceae"
# Ilex
LES$Family[which(LES$Genus=="Ilex")] <- "Aquifoliaceae"
# Luzula
LES$Family[which(LES$Genus=="Luzula")] <- "Juncaceae"
# Miconia
LES$Family[which(LES$Genus=="Miconia")] <- "Melastomataceae"
#Olea Oleaceae -> note: a couple spp have subspp that don't get pooled.
LES$Family[which(LES$Genus=="Olea")] <- "Oleaceae"
# Salix --> SALIX IS CLEAN
LES$Species[which(LES$Species=="Salix alba fraxinus")] <- "Salix alba"
LES$GE.SP[which(LES$Species=="Salix.alba fraxinus")] <- "Salix.alba"
LES$Species[which(LES$Species=="Salix dasyclados phylicifolia")] <- "Salix dasyclados"
LES$GE.SP[which(LES$Species=="Salix.dasyclados phylicifolia")] <- "Salix.dasyclados"
LES$Species[which(LES$Species=="Salix fragilus")] <- "Salix fragilis"
LES$GE.SP[which(LES$Species=="Salix.fragilus")] <- "Salix.fragilis"

LES$Family[which(LES$Genus=="Salix")] <- "Salicaceae"
LES$LMA <- 10^LES$log.LMA
LES$LL <- 10^LES$log.LL
LES$Nmass <- 10^LES$log.Nmass
LES$Narea <- 10^LES$log.Narea
LES$Aarea <- 10^LES$log.Aarea
LES$Amass <- 10^LES$log.Amass
LES$SLA <- 1/LES$LMA
LES$Acats <- rep("H", times=nrow(LES)) # call everything >25 H
LES$Acats[which(LES$Aarea<25 & LES$Aarea>11.5)] <- "MH" # everything median (11.5) to 25 = MH
LES$Acats[which(LES$Aarea<11.5 & LES$Aarea>7)] <- "ML" # everything from median to 1st qu (7) = ML
LES$Acats[which(LES$Aarea<7)] <- "L"
LES$Acats <- factor(LES$Acats)
LES$Acats[which(is.na(LES$Aarea))] <- NA
levels(LES$Acats) <- list(H="H",MH="MH",ML="ML",L="L")
## add in climate
LES$MAT <- LESclim$MAT[match(LES$Dataset, LESclim$Dataset)]
LES$MAP <- LESclim$Rain[match(LES$Dataset, LESclim$Dataset)]
LES$VPD <- LESclim$VPD[match(LES$Dataset, LESclim$Dataset)]
LES$RAD <- LESclim$RAD[match(LES$Dataset, LESclim$Dataset)]
LES$PET <- LESclim$PET[match(LES$Dataset, LESclim$Dataset)]

# ggplot(LES, aes(x=log.Nmass, y=log.LMA, col=families )) + geom_point()+ geom_smooth(method="lm", se=F) + theme(legend.position="none")
# length(which(is.na(families)))
# family.itis1$family[grep(pattern =  '^c', family.itis1$family)] <- newfams1



## Species ID lookup table to relate traits and biomass
speciesID <- read.csv("SpeciesNames_lookup_GOOD_031316.csv", header=T)


## Soil Data
soil <- read.csv("NACP_TERRA_PNW_soil_cleaned.csv", header=T, na.strings = "-9999")

#### Having to clean some shit up ...
soil[which(soil$Layer=="top" & soil$UpperDepth>5),]
# Plot 252 has two top layers. need to switch the second to 'bottom'
# Plot 86 has top and bottom layers flipped
## and something on the order of 9 plots don't have a 'top layer
soil$Layer[which(soil$PLOT_ID==252)] <- c("top", "middle","bottom")
soil$Layer[which(soil$PLOT_ID==86)] <- c("top", "bottom")
#probs1 <- names(xtabs(~PLOT_ID, soil)[which(xtabs(~PLOT_ID, soil)<2)])
#soil[which(soil$PLOT_ID %in% probs & soil$Layer!="top"),] # all plots with only 1 layer have 'top'
probs2 <- names(xtabs(~PLOT_ID, soil)[which(xtabs(~PLOT_ID, soil)==2)])
probs2names <- xtabs(~PLOT_ID, soil[which(soil$PLOT_ID %in% probs2 & soil$Layer != "top"),])
## 8 plots have 'middle' and 'bottom', but middle starts at 0. I'm going to just make all 'middles' into 'tops
soil$Layer[which(soil$PLOT_ID %in% names(probs2names[which(probs2names==2)]) & soil$Layer=="middle")] <- "top" # 8 plots have only 2 layers and one of them is not 'top'
# OK, that seemed to solve it, I think... Everything has a 'top' layer, but below that is unknown. could be 'middle' could be 'bottom'...
# which(xtabs(~PLOT_ID, soil[which(soil$Layer=='top'),])>1)
# soiltest <- soil %>% melt(id.vars = 1:16)
# soiltest$variable <- paste(soiltest$Layer, soiltest$variable,sep=".")
# soilwide <- soiltest %>% spread(variable, value)


## tried to spread things out, but we'll try just working with the top layer for shits and giggles...
soil.top <- soil[which(soil$Layer=="top"), ]
# plot(soil_N~MAP, soil.top) # N generally increases with Precip
# plot(soil_N~MAT, soil.top) # N generally increases with T, but probably only if there's precip?
# ggplot(soil.top, aes(x=MAT, y=soil_N, col=MAP)) + geom_point()
    # there seems to be an interaction between MAP and MAT, MAT increases soil_N only at high MAPS.
# ggplot(soil.top, aes(x=LowerDepth, y=soil_N, col=MAP)) + geom_point()
    # high N is only achieved with deep soils, and then probably only when MAP is high...


### Stand characteristics dataset ###
# on 5.27.16 I deleted the units row in this dataset to make _v2 so that types import correctly.
biomass <- read.csv("NACP_TERRA_PNW_forest_biomass_productivity_v2.csv", header= T,na.strings = "-9999" )
biomass$soil_N <- soil.top$soil_N[match(biomass$PLOT_ID, soil.top$PLOT_ID)]
biomass$RGR <- log(biomass$AG_BIOMASS_TREE_TOTAL_AS_CARBON) - log(biomass$AG_BIOMASS_TREE_TOTAL_AS_CARBON - biomass$AG_PROD_TREE_TOTAL_AS_CARBON)
## control growth for existing biomass
growthmod <- lm(AG_PROD_TREE_TOTAL_AS_CARBON~AG_BIOMASS_TREE_TOTAL_AS_CARBON, data=biomass)
growthmod2 <- gam(AG_PROD_TREE_TOTAL_AS_CARBON~s(AG_BIOMASS_TREE_TOTAL_AS_CARBON), data=biomass)
growthmod3 <- gam(AG_PROD_TREE_TOTAL_AS_CARBON~s(log(AG_BIOMASS_TREE_TOTAL_AS_CARBON)), data=biomass)
biomass$logAG_BIOMASS <- log(biomass$AG_BIOMASS_TREE_TOTAL_AS_CARBON)
# growthmod4 <- gls(AG_PROD_TREE_TOTAL_AS_CARBON~logAG_BIOMASS, data=biomass, weights = varExp(~logAG_BIOMASS))
# there's actually a better relationship between log(biomass) than regular biomass
biomass$BIOstGROWTH[which(biomass$AG_PROD_TREE_TOTAL_AS_CARBON>-5)] <- resid(growthmod)
biomass$BIOstGROWTHgam[which(biomass$AG_PROD_TREE_TOTAL_AS_CARBON>-5)] <- resid(growthmod2)
biomass$logBIOstGROWTHgam[which(biomass$AG_PROD_TREE_TOTAL_AS_CARBON>-5)] <- resid(growthmod3)


#biomass <- biomass.raw[-1,] # drop second row of csv with measurement units
#biomass$PROJECT <- factor(biomass$PROJECT) # reformat columns after units row dropped
#biomass$PLOT_ID_AMERIFLUX <- factor(biomass$PLOT_ID_AMERIFLUX)
#biomass$ECOREGION <- factor(biomass$ECOREGION)
#biomass$PLOT_ID <- factor(biomass$PLOT_ID)
#num.cols <- c(10:14,16,18,20,22,24:32)
#for(i in num.cols){ biomass[,i] <- as.numeric(as.character(biomass[,i]))}

### Traits dataset ####
# newest version = NACP_TERRA_PNW_leaf_traits_v1_plusClim.csv
  # older versions (on order they became obsolete): NACP_TERRA_PNW_leaf_trait.csv

traits <- read.csv("NACP_TERRA_PNW_leaf_trait_v1_plusClim.csv", header=T, na.strings="-9999")
# traits <- traits.raw[-1,]
fac.colstraits <- c(1,5:7,12:16)
for(j in fac.colstraits){ traits[,j] <- factor(traits[,j])}
num.colstraits <- c(8:11,17:27)
for (i in num.colstraits){traits[,i] <- as.numeric(as.character(traits[,i]))}

### back calculate some missing values:
# calculate LEAF_CARBON for METOFIRE plots that have LEAF_CARBON_WT, but not elemental analysis
traits$LEAF_CARBON[which(is.na(traits$LEAF_CARBON))] <- with(traits[which(is.na(traits$LEAF_CARBON)),], LEAF_CARBON_WT/LEAF_DRY_WT * 100)
# length(which(is.na(traits$LEAF_CARBON))) # this should return 1

# add in Genus.species column, some clim columns and a column for species IDs used in biomass dataset
traits$GE.SP <- paste(traits$GENUS, traits$SPECIES, sep=".")
traits$SP.ID <- speciesID$bio.sp[match(traits$GE.SP, speciesID$traits.sp)]
traits$MAT <-biomass$MAT_C[match(traits$PLOT_ID, biomass$PLOT_ID)]
traits$MAP <- biomass$MAP[match(traits$PLOT_ID, biomass$PLOT_ID)] 
traits$soilmoist.all.mm <- apply(traits[, c("soilmoist.lvl1.mm", "soilmoist.lvl2.mm", "soilmoist.lvl3.mm")],MARGIN = 1, FUN=sum)
traits$FOREST_TYPE <- biomass$SPP_O1_ABBREV[match(traits$PLOT_ID, biomass$PLOT_ID)]
traits$dominance <- biomass$SPP_O1_BASAL_AREA_FRACTION[match(traits$PLOT_ID, biomass$PLOT_ID)]
traits$ASA <- biomass$ASA[match(traits$PLOT_ID, biomass$PLOT_ID)]
traits$ELEVATION <- biomass$ELEVATION[match(traits$PLOT_ID, biomass$PLOT_ID)]
traits$AG_TBIO <- biomass$AG_BIOMASS_TREE_TOTAL_AS_CARBON[match(traits$PLOT_ID, biomass$PLOT_ID)]
traits$AG_TGROWTH <- biomass$AG_PROD_TREE_TOTAL_AS_CARBON[match(traits$PLOT_ID, biomass$PLOT_ID)]
traits$AG_WBIO <- biomass$AG_BIOMASS_TREE_WOOD_AS_CARBON[match(traits$PLOT_ID, biomass$PLOT_ID)]
traits$AG_WGROWTH <- biomass$AG_PROD_TREE_WOOD_AS_CARBON[match(traits$PLOT_ID, biomass$PLOT_ID)]
traits$AG_FBIO <- biomass$AG_BIOMASS_TREE_FOLIAGE_AS_CARBON[match(traits$PLOT_ID, biomass$PLOT_ID)]
traits$AG_FGROWTH <- biomass$AG_PROD_TREE_FOLIAGE_AS_CARBON[match(traits$PLOT_ID, biomass$PLOT_ID)]

traits$BIOST_TGROWTH <- biomass$BIOstGROWTH[match(traits$PLOT_ID, biomass$PLOT_ID)]
traits$BIOST_TGROWTHgam <- biomass$BIOstGROWTHgam[match(traits$PLOT_ID, biomass$PLOT_ID)]
traits$logBIOST_TGROWTHgam <- biomass$logBIOstGROWTHgam[match(traits$PLOT_ID, biomass$PLOT_ID)]
## relative growth rate (assume TBIO = M2 and TGROWTH=M2-M1), and TGROWTH is per year
  # so I need to subtract TGROWTH from TBIO to get M1
traits$RGR <- log(traits$AG_TBIO) - log(traits$AG_TBIO - traits$AG_TGROWTH)
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

## LES traits

traits$LMA <- with(traits, LEAF_DRY_WT/(LEAF_HSA*1/100^2)) # LES works in LMA rather than SLA
  # also SLA is in gC rather than LEAF_DRY_WT, and in cm2 rather than m2 as in Wright 2004 
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
traits$PSA_to_HSA[which(traits$SP.ID=="PINPON")] <- 1.18 # fill in missing PINPON values
traits$PSA_to_HSA[which(traits$SP.ID=="PINCON")] <- 1.29 # fill in missing PINCON values
traits$PSA_to_HSA[which(traits$GENUS=="Pinus" & is.na(traits$PSA_to_HSA))] <- 1.2
traits$LMA_PSA <- traits$LMA * traits$PSA_to_HSA # * PSA_to_HSA is same as /HSA_to_PSA. and LMA is mass/area
traits$log.LMA_PSA <- log(traits$LMA_PSA, base=10)
## make an RGR that only has values when that species makes up >50% of the stand BA
traits$RGRdom <- traits$RGR
traits$RGRdom[which(as.character(traits$FOREST_TYPE) != as.character(traits$SP.ID))] <- NA
traits$RGRdom[which(as.character(traits$FOREST_TYPE) == as.character(traits$SP.ID) & traits$dominance<50)] <- NA
traits$stGrowthdom <- traits$BIOST_TGROWTHgam
traits$stGrowthdom[which(as.character(traits$FOREST_TYPE) != as.character(traits$SP.ID))] <- NA
traits$stGrowthdom[which(as.character(traits$FOREST_TYPE) == as.character(traits$SP.ID) & traits$dominance<50)] <- NA

########### END Data import and format ##############


traits$FullSpecies <- paste(traits$GENUS, traits$SPECIES, sep=" ")
# PNWspecies <- sort(unique((traits$FullSpecies)))
# PNWfams <- tax_name(PNWspecies, get="family", db="itis")
# PNWfams$family[grep("Cercocarpus", PNWfams$query)] <- "Rosaceae"
# PNWfams$family[grep("Cornus", PNWfams$query)] <- "Cornaceae"
# PNWfams$family[grep("Frangula", PNWfams$query)] <- "Rhamnaceae" # actually genus Rhamnus
# PNWfams$family[grep("Purshia tridentate", PNWfams$query)] <- "Rosaceae" # actually tridentata
# 
# write.csv(PNWfams, "PNWfams_031717.csv")

PNWfams <- read.csv("PNWfams_031717.csv")
traits$Family <- PNWfams$family[match(traits$FullSpecies, PNWfams$query)]







#----------------------------------------------------------------------------------
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
###### **Analysis of intra-specific trait variation in common species** ###########
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#----------------------------------------------------------------------------------

# definition of common species = >20 records




######### Creating traits.common (and summaries) ################ 

### Look at species coverage
speciescoverage <- traits %>% group_by(GENUS) %>% summarize(nspecies=length(unique(SPECIES)))
speciescoverage <- speciescoverage[order(speciescoverage$nspecies, decreasing=T),]
  ## So it looks like we've got 24 genera, 6 of which have more than one species. only 3 have >2 spp.

### Look at coverage within genera
genera <- traits%>%group_by(GENUS)%>% summarise(nspp=length(unique(SPECIES)))
common.genera <- genera$GENUS[which(genera$nspp>1)]


### pull out common species (>10 records)
common.species10 <- names(which(xtabs(~GE.SP, traits)>10)) # species with >20 records (12), >50 (7), >10 (16), 
common.species5 <- names(which(xtabs(~GE.SP, traits)>=5)) # 23 spp
#lesscommon.species <- names(which(xtabs(~GE.SP, traits)>10)) # species with >20 records (12), >50 (7), >10 (16), 
#traits.lesscommon <- subset(traits, subset=GE.SP %in% lesscommon.species)

### make traits.common dataset
traits.common <- subset(traits, subset=GE.SP %in% common.species10)
  # but Purshia tridendata has a bunch of measurements from only 2 sites. need to remove I think.
traits.common <- traits.common[-which(traits.common$GE.SP == "Purshia.tridentate"),]
traits.common$PLOT_ID <- factor(traits.common$PLOT_ID) # have 3 plots with no common species, so dropping them
traits.common$GE.SP <- factor(traits.common$GE.SP)
traits.common$GENUS <- factor(traits.common$GENUS)
traits.common$SPECIES <- factor(traits.common$SPECIES)
traits.common$SLA_HSAsc <- c(scale(traits.common$SLA_HSA))
traits.common$FOREST_TYPE <- factor(traits.common$FOREST_TYPE)
traits.common$SP.ID <- factor(traits.common$SP.ID)
# make a PCA of climate variables and add them to the df
climpca <- prcomp(traits.common[,c(grep("gy", colnames(traits.common)), which(colnames(traits.common) %in% c("soilmoist.lvl1.mm","soilmoist.all.mm")))],scale=T)
traits.common$climPC1 <- climpca$x[,1]
traits.common$climPC2 <- climpca$x[,2]
traits.common$climPC3 <- climpca$x[,3]

##### more inclusive traits.common5, for traits analysis w/ as many species as possible
traits.common5 <- subset(traits, subset=GE.SP %in% common.species5)
traits.common5$PLOT_ID <- factor(traits.common5$PLOT_ID) # have 3 plots with no common5 species, so dropping them
traits.common5$GE.SP <- factor(traits.common5$GE.SP)
traits.common5$GENUS <- factor(traits.common5$GENUS)
traits.common5$SPECIES <- factor(traits.common5$SPECIES)
traits.common5$SLA_HSAsc <- c(scale(traits.common5$SLA_HSA))
traits.common5$FOREST_TYPE <- factor(traits.common5$FOREST_TYPE)
traits.common5$SP.ID <- factor(traits.common5$SP.ID)
climpca5 <- prcomp(traits.common5[,c(grep("gy", colnames(traits.common5)), which(colnames(traits.common5) %in% c("soilmoist.lvl1.mm","soilmoist.all.mm")))],scale=T)
traits.common5$climPC1 <- climpca5$x[,1]
traits.common5$climPC2 <- climpca5$x[,2]
traits.common5$climPC3 <- climpca5$x[,3]
  # climate pca on sites rather than measurements. I suppose the above PCA has replicated climate data for plots with multiple measurements
# plotclims <- traits %>% group_by(PLOT_ID) %>% summarise(MAT = unique(tmean.gy.c), MAP=unique(ppt.gy.mm),
#                                                         CMI = unique(cmi.gy.mm), VPD = unique(vpd.gy.max),
#                                                         soil1 = unique(soilmoist.lvl1.mm), soilall=unique(soilmoist.all.mm))
# climpcanew <- prcomp(plotclims[,-1],scale=T)
  # results are essentially identical. It flips the sign of PC2 so it's 'coldness', but loadings are identical and
  # proportion explained by PC1 and PC2 are essentially identical






### also making one with the common trait spp and common genera
traits.commongenera <- subset(traits, subset=GE.SP %in% common.species | GENUS %in% common.genera)
traits.commongenera <- traits.commongenera[-which(traits.commongenera$GE.SP == "Purshia.tridentate"),]
traits.commongenera$PLOT_ID <- factor(traits.commongenera$PLOT_ID) # have 3 plots with no commongenera species, so dropping them
traits.commongenera$GE.SP <- factor(traits.commongenera$GE.SP)
traits.commongenera$GENUS <- factor(traits.commongenera$GENUS)
traits.commongenera$SPECIES <- factor(traits.commongenera$SPECIES)
traits.commongenera$SLA_HSAsc <- scale(traits.commongenera$SLA_HSA)
traits.commongenera$FOREST_TYPE <- factor(traits.commongenera$FOREST_TYPE)
traits.commongenera$SP.ID <- factor(traits.commongenera$SP.ID)

# removing all NAs drops out 148 rows.
traits.common.narm <- traits.common[-which(is.na(traits.common$SLA_HSAsc)|is.na(traits.common$LAI_O)|is.na(traits.common$HEIGHTC_m)|is.na(traits.common$soilmoist.lvl1.mm)),]



### assess missing trait values by species, and number of plots for each species
traits.common %>% group_by(GE.SP) %>% summarize(nplots = length(unique(PLOT_ID)), nsamples = n(), naLAI =length(which(is.na(LAI_O))), naHeight = length(which(is.na(HEIGHTC_m))), naSLA_PSA =length(which(is.na(SLA_PSA))), anSLA_HSA =length(which(is.na(SLA_HSA))),  naCN = length(which(is.na(LEAF_CN))), naLIFE = length(which(is.na(LEAF_LIFE))))   # decent plot replication for many species
  # take home: use SLA_HSA (because it has fewer missing values), Junocc doesn't have LEAF_LIFE (scales rather than leaves), other species missing values all just reduce replication rather than eliminate species
traits.common %>% group_by(GENUS) %>% summarize(length(unique(SPECIES)))# so not great species replication within genera
plotrep <- traits.common[!is.na(traits.common$SLA_PSA),] %>% group_by(GENUS, SPECIES, PLOT_ID) %>% summarize(nsamples = n())
  # most plots only have 2-4 samples per plot. 
  #   but 19 plots have >5 samples
  #   and 10 plots have >10 samples





#-----------------------------------------------------------------------------------

########### Initial Data Exploration (visualization, colinearity, etc.) ###########

#-----------------------------------------------------------------------------------



######## Are projects comparable? ##################
# the main takehome from this is that I probably need a random intercept for PROJECT in most analyses. 
# It only REALLY shows up in % C, where EPA & COHO are pretty consistently higher than CADIS and ORCA

## all species
p1 <- ggplot(traits.common, aes(y=LEAF_CARBON, x=LEAF_NITROGEN, colour=PROJECT)) + geom_point()
p2 <- ggplot(traits.common, aes(y=LEAF_CARBON, x=AG_FGROWTH, colour=PROJECT)) + geom_point()
p3 <- ggplot(traits.common, aes(y=LEAF_CARBON, x=MAP, colour=PROJECT)) + geom_point()
p4 <- ggplot(traits.common[which(traits.common$SP.ID == "PINPON"),], aes(x=LEAF_NITROGEN, y=LEAF_CARBON, colour=PROJECT)) + geom_point() + geom_smooth(method='lm') + ggtitle("PIPO")
p5 <- ggplot(traits.common[which(traits.common$SP.ID == "PSEMEN"),], aes(x=LEAF_NITROGEN, y=LEAF_CARBON, colour=PROJECT)) + geom_point() + geom_smooth(method='lm') + ggtitle("PSME")
p6 <- ggplot(traits.common[which(traits.common$SP.ID == "TSUHET"),], aes(x=LEAF_NITROGEN, y=LEAF_CARBON, colour=PROJECT)) + geom_point() + geom_smooth(method='lm') + ggtitle("TSHE")
p7 <- ggplot(traits.common, aes(x=SLA_HSA, y=LEAF_LIFE, colour=PROJECT)) + geom_point() + geom_smooth(method='lm')
p8 <- ggplot(traits.common, aes(y=AG_FGROWTH, x=LAI_O, colour=PROJECT)) + geom_point() + geom_smooth(method='lm')
p9 <- ggplot(traits.common, aes(y=ASA, x=SP.ID)) + geom_boxplot() + geom_point(aes(colour=PROJECT))


# +++++++ 3 Summary Figures  +++++++++++++++
  # exploring how traits differ by project and what not
quartz(width=4, height=6.5) # LeafCarbon_ProjectComparison.pdf
multiplot(plotlist=list(p1,p2,p3)) # EPA and COHO pretty consistently higher than CADIS and ORCA
quartz(width=4, height=6.5) # LeafCarbon_BySpecies_Comparison.pdf
multiplot(plotlist=list(p4,p5,p6)) # diff spp sampled in diff projects (eg. only EPA for TSHE), but Carbon thing holds
quartz(width=7, height=5) # ProjectOverviews_traitsgrowthASA.pdf
multiplot(plotlist=list(p7,p8,p9), layout=matrix(c(1,2,3,3), nrow = 2, byrow = T))




#################### Spatial distribution of plots #######################

######RgoogleMaps
library(RgoogleMaps)

mapcols <- brewer.pal(n=6, "Set1")
palette(mapcols)
# ## getting the appropriate google map
siteTerrain <- GetMap(center=c(lat=42.02,lon=-120.511),
                      destfile="RgoogleMap2.png",
                      maptype="hybrid",zoom=6)
# 
siteTerrain <- GetMap(center=c(lat=mean(traits$LATITUDE),lon=mean(traits$LONGITUDE)),
                      destfile="RgoogleMap3.png",
                      maptype="hybrid",zoom=9)
# ## adding missing trees as pts
tmp <- PlotOnStaticMap(siteTerrain, lat = traits$LATITUDE,
                       lon = traits$LONGITUDE, 
                       destfile = "RgoogleMap3.png", cex=1.5,pch=20,
                       col=traits$PROJECT,FUN=points,add=F)
legend("right", legend=levels(traits$PROJECT), pch=16, col=mapcols[1:length(unique(traits$PROJECT))], bg="white")


######
genus.summaries <- traits %>% group_by(GENUS) %>% summarize(nspp = length(unique(SP.ID)), ntraits=n())





#### INITAL VISUALIZATION ####################
pairs(traits.common[,c(10,11,23,24,25,26,27,29,30)],col=traits.common$SPECIES, pch=16, cex=.8)
Mypairs(traits.common[,c(10,11,23,24,25,26,27,29,30)])
## based on this, it looks like traits that are worth looking at are:
# SLA_PSA (SLA from projected leaf area)
# LEAF_CARBON, LEAF_NITROGEN, LEAF_CN (probably just leaf CN)
# NEW OBSERVATION 3/22/16: leaf C (%) AND leaf N (%) both seem to decrease with LEAF_LIFE
#     this is the strongest corr in the whole shebang. I probably need to analyze these seperately!!!
# LEAF_LIFE



### Visualize SLA data
xyplot(SLA_HSA~soilmoist.lvl1.mm|GE.SP,traits.common)
xyplot(SLA_HSA~vpd.gy.max|GE.SP,traits.common)
xyplot(SLA_HSA~cmi.gy.mm|GE.SP,traits.common)
xyplot(SLA_HSA~ppt.gy.mm|GE.SP,traits.common) # fills in a few points that aren't in plot data? so better than MAP
# xyplot(SLA_HSA~MAP|GE.SP,traits.common)
xyplot(SLA_HSA~tmean.gy.c|GE.SP,traits.common)
# xyplot(SLA_HSA~MAT|GE.SP,traits.common) # tmean.gy.c better because fills in missing plots?
xyplot(SLA_HSA~FOREST_TYPE|GE.SP, traits.common)

# Take home: no extremly large/obvious climate related changes in SLA that just jump out


### Visualize LEAF_CN data
xyplot(LEAF_CN~soilmoist.lvl1.mm|GE.SP,traits.common, type=c('p','smooth',"r"), col.line= "red3")
xyplot(LEAF_CN~vpd.gy.max|GE.SP,traits.common, type=c('p','smooth',"r"), col.line= "red3")
xyplot(LEAF_CN~cmi.gy.mm|GE.SP,traits.common, type=c('p','smooth',"r"), col.line= "red3")
xyplot(LEAF_CN~ppt.gy.mm|GE.SP,traits.common, type=c('p','smooth',"r"), col.line= "red3") # fills in a few points that aren't in plot data? so better than MAP
# xyplot(LEAF_CN~MAP|GE.SP,traits.common)
xyplot(LEAF_CN~tmean.gy.c|GE.SP,traits.common)
# xyplot(LEAF_CN~MAT|GE.SP,traits.common) # tmean.gy.c better because fills in missing plots?

# Take home: a few species show trends in response to one variable, but rarely the same variable, and typically not the super common species


### Visualize LEAF_LIFE data
xyplot(LEAF_LIFE~soilmoist.lvl1.mm|GE.SP,traits.common, type=c('p','smooth',"r"), col.line= "red3")
xyplot(LEAF_LIFE~vpd.gy.max|GE.SP,traits.common, type=c('p','smooth',"r"), col.line= "red3")
xyplot(LEAF_LIFE~cmi.gy.mm|GE.SP,traits.common, type=c('p','smooth',"r"), col.line= "red3")
xyplot(LEAF_LIFE~ppt.gy.mm|GE.SP,traits.common,  type=c('p','smooth',"r"), col.line= "red3") # fills in a few points that aren't in plot data? so better than MAP
# xyplot(LEAF_LIFE~MAP|GE.SP,traits.common)
xyplot(LEAF_LIFE~tmean.gy.c|GE.SP,traits.common,  type=c('p','smooth',"r"), col.line= "red3")
# xyplot(LEAF_LIFE~MAT|GE.SP,traits.common) # tmean.gy.c better because fills in missing plots?

# Take home: Some STRONG temp relationships in Doug fir, but not always consistent with other species.


xyplot(LEAF_LIFE~LEAF_CN|GE.SP,traits.common, type=c('p','smooth',"r"), col.line= "red3") # TSHE and ABGR show visible patterns, but other than that nothing
xyplot(LEAF_LIFE~SLA_HSA|GE.SP,traits.common, type=c('p',"r"), col.line= "red3") # leaf lifespan and SLA may be related in true firs and doug fir
xyplot(SLA_HSA~LEAF_LIFE|GE.SP,traits.common, type=c('p',"r"), col.line= "red3") # leaf lifespan and SLA may be related in true firs and doug fir
xyplot(SLA_HSA~LEAF_CN|GE.SP,traits.common,  type=c('p','smooth',"r"), col.line= "red3") # nothing obvious except maybe a negative relationship in TSHE
xyplot(LEAF_NITROGEN~LEAF_LIFE|GE.SP, traits.common, type=c('p','smooth',"r"), col.line= "red3")
xyplot(LEAF_CARBON~LEAF_LIFE|GE.SP, traits.common, type=c('p',"r"), col.line= "red3")
xyplot(LEAF_NITROGEN~SLA_HSA|GE.SP, traits.common, type=c('p','smooth',"r"), col.line= "red3") # not much going on except maybe TSHE
xyplot(LEAF_CARBON~SLA_HSA|GE.SP, traits.common, type=c('p','smooth',"r"), col.line= "red3") # not much going on


xyplot(LEAF_LIFE~tmean.gy.c|GE.SP, traits.common)
plot(LEAF_LIFE~tmean.gy.c, traits.common, subset=GENUS=="Abies", col=SPECIES, pch=16, cex=1)
for(i in unique(traits.common$GE.SP[which(traits.common$GENUS=="Abies")])){
  abline(lm(LEAF_LIFE~tmean.gy.c, traits.common, subset=GE.SP==i), lwd=2)
}

plot(LEAF_LIFE~tmean.gy.c, traits.common, subset=GENUS=="Pinus", col=SPECIES, pch=16, cex=1)
for(i in unique(traits.common$GE.SP[which(traits.common$GENUS=="Pinus")])){
  abline(lm(LEAF_LIFE~tmean.gy.c, traits.common, subset=GE.SP==i), lwd=2)
}

plot(LEAF_LIFE~tmean.gy.c, traits.common, col=GE.SP, pch=16, cex=1)
for(i in levels(traits.common$GE.SP)[-c(6,7,8)]){
  sig <- summary(lm(LEAF_LIFE~tmean.gy.c, traits.common, subset=GE.SP==i))
  if(sig$coefficients[2,4]>0.05)
  abline(lm(LEAF_LIFE~tmean.gy.c, traits.common, subset=GE.SP==i), lwd=1,lty=2, col=mypal[which(levels(traits.common$GE.SP)==i)])
  else
    abline(lm(LEAF_LIFE~tmean.gy.c, traits.common, subset=GE.SP==i), lwd=3,lty=1, col=mypal[which(levels(traits.common$GE.SP)==i)])
  
  }







#-----------------------------------------------------------------------------
########### do things relate to site productivity? ######
#-----------------------------------------------------------------------------


plot(LEAF_LIFE~AG_TGROWTH, traits.common, col=SP.ID)
for(i in levels(traits.common$SP.ID)){
  tmp <- traits.common[which(traits.common$SP.ID==i),]
  abline(lm(LEAF_LIFE~AG_TGROWTH, tmp),col=mypal[which(levels(traits.common$GE.SP)==i)])
}

plot(LEAF_NITROGEN~AG_TGROWTH, traits.common, col=SP.ID)
for(i in levels(traits.common$SP.ID)){
  tmp <- traits.common[which(traits.common$SP.ID==i & !is.na(traits.common$LEAF_NITROGEN) & !is.na(traits.common$AG_TGROWTH)),]
  mod <- lm(LEAF_NITROGEN~AG_TGROWTH, tmp)
  points (fitted(mod)~AG_TGROWTH,data=tmp, col=SP.ID, type="l")
}
  ## Leaf N does seem to increase a bit with growth...

plot(AG_TGROWTH~vpd.gy.max, traits.common, col=SP.ID)


####### 3d plotting options:

## plot3d makes rotatable 3d pop out window
with(traits.common, plot3d(tmean.gy.c, ppt.gy.mm, SLA_HSA, box=F, col=as.numeric(GE.SP)))

with(traits.common[which(traits.common$SP.ID=="PSEMEN"),], plot3d(SLA_HSA, LEAF_LIFE, LEAF_CARBON, box=F, col=as.numeric(GE.SP)))

with(traits.common[which(traits.common$GENUS=="Pinus"),], plot3d(tmean.gy.c, ppt.gy.mm, SLA_HSA, box=F, col=as.numeric(GE.SP)))
with(traits.common[which(traits.common$GENUS=="Abies"),], plot3d(tmean.gy.c, ppt.gy.mm, SLA_HSA, box=F, col=as.numeric(GE.SP)))
with(traits.common[which(traits.common$GE.SP=="Pseudotsuga.menziesii"),], plot3d(tmean.gy.c, ppt.gy.mm, SLA_HSA, box=F, col=as.numeric(GE.SP)))
    # for doug fir, there's a strong effect of PPT, but only at higher temperatures!!!
# with(traits.common[which(traits.common$GE.SP=="Pseudotsuga.menziesii"),], plot3d(vpd.gy.max, ppt.gy.mm, SLA_HSA, box=F, col=as.numeric(GE.SP)))
with(traits.common[which(traits.common$GE.SP=="Pseudotsuga.menziesii"),], plot3d(tmean.gy.c, ppt.gy.mm, LEAF_LIFE, box=F, col=as.numeric(GE.SP)))

with(traits.common[which(traits.common$GE.SP=="Pinus.ponderosa"),], plot3d(tmean.gy.c, ppt.gy.mm, SLA_HSA, box=F, col=as.numeric(GE.SP)))


with(traits.common[which(traits.common$SP.ID=="PSEMEN"),], plot3d(vpd.gy.max, soilmoist.all.mm, AG_TGROWTH, box=F, col=as.numeric(GE.SP)))

### not yet implemented
# ## scatterplot3d makes static 3d plot in plotting device, but can add tons of stuff.
# s3d <- with(climates, scatterplot3d(PPT_wt, DD5_sm, CMD_sm, pch=pchs, highlight.3d = F, type="h",  color=colsdark,angle = 102, xlab="Winter Precip", ylab="Summer GDD",zlab="Summer CMD", box = T, grid=T))
# s3d$points3d(x = climates$PPT_wt, y=rep(1400, times=nrow(climates)), z=climates$CMD_sm,pch=climates$pchs+14, col=colslight,cex=0.8)
# s3d$points3d(x = climates$PPT_wt, y=climates$DD5_sm, z=rep(50, times=26), pch=climates$pchs+14, col=colslight, cex=0.8)
# s3d$points3d(x = rep(1000, times=26), y=climates$DD5_sm, z=climates$CMD_sm, pch=climates$pchs+14, col=colslight, cex=0.8)







#-----------------------------------------------------------------------
######## Climate colinearity and associated infrences: ########
#-----------------------------------------------------------------------------


Mypairs(traits.common[,c(10,11,grep("gy", colnames(traits.common)), grep("soilmoist", colnames(traits.common)))])
  ## based on this, it looks like:
  # ppt and cmi are VERY similar (vpd has less of an effect and is also colinear with ppt)
  # the three soil moistures are pretty colinear
  # ALL the soil moisture is in lvl3, deep soil.
  # there is apparent no relationship between climate, LAI, and canopy height when everything is pooled.

pairs(traits.common[,c(10,11,grep("gy", colnames(traits.common)), grep("soilmoist", colnames(traits.common)))], bg=traits.common$GE.SP, pch=16, cex=.6, diag.panel = panel.hist, lower.panel = panel.smooth)


##### from dicking around with 3d plot, noticed that soil moisture is a pretty solid function of tmean and ppt
# more detailed observations
plot(tmean.gy.c~ppt.gy.mm, traits.common, cex=2+scale(soilmoist.all.mm))
plot(tmean.gy.c~ppt.gy.mm, traits.common, cex=2+scale(soilmoist.all.mm), subset=GE.SP=="Pseudotsuga.menziesii")
summary(lm(soilmoist.all.mm~scale(tmean.gy.c) + scale(ppt.gy.mm), traits.common))
  # so ppt drives most of the variation, but tmean is still a really significant part. together r2=0.78
summary(lm(soilmoist.all.mm~scale(tmean.gy.c) + scale(ppt.gy.mm), traits.common, subset=GE.SP=="Pseudotsuga.menziesii"))
  # same is true for just PSME.
summary(lm(soilmoist.all.mm~scale(tmean.gy.c) + scale(ppt.gy.mm) + scale(vpd.gy.max), traits.common))
  # also worth noting that VPD is highly colinear with ppt and explains no additional soil moisture variation

### additional observations:
  # 1) soilmoist.lvl1 is driven by ppt. BUT IT HOLDS ~ no water!!
plot(tmean.gy.c~ppt.gy.mm, traits.common, cex=2+scale(soilmoist.lvl1.mm), subset=GE.SP=="Pseudotsuga.menziesii")
summary(lm(soilmoist.lvl1.mm~scale(tmean.gy.c) + scale(ppt.gy.mm), traits.common))
  # 2) soilmoist.lvl3 is driven by both (and looks like soilmoist.all, becuase it's the largest contributer!)
plot(tmean.gy.c~ppt.gy.mm, traits.common, cex=2+scale(soilmoist.lvl3.mm), subset=GE.SP=="Pseudotsuga.menziesii")
summary(lm(soilmoist.lvl3.mm~scale(tmean.gy.c) + scale(ppt.gy.mm), traits.common))

plot(SLA_HSA~ppt.gy.mm, traits.common, cex=2+scale(tmean.gy.c), col=FOREST_TYPE)


# if just looking for non-colinear combinations of variables: tmean, vpd, and soilmoist.lvl1 are it.
corvif(traits.common[,c(grep("gy", colnames(traits.common)), which(colnames(traits.common) %in% c("soilmoist.lvl1.mm","soilmoist.all.mm")))])
  # highest vif is cmi.gy.mm (168.67)
corvif(traits.common[,c("tmean.gy.c", "ppt.gy.mm", "vpd.gy.max", "soilmoist.lvl1.mm", "soilmoist.all.mm")])
  # next is soilmoist.all.mm
corvif(traits.common[,c("tmean.gy.c", "ppt.gy.mm", "vpd.gy.max", "soilmoist.lvl1.mm")])
  # next is ppt.gy.mm
corvif(traits.common[,c("tmean.gy.c", "vpd.gy.max", "soilmoist.lvl1.mm")])
  # now all vifs <2 (a vif of 2.5 = r2 of .6 with other predictors)

## NOTE: the climate PCA has been moved up top to create all the necessary columns for traits.common
# climpca <- prcomp(traits.common[,c(grep("gy", colnames(traits.common)), which(colnames(traits.common) %in% c("soilmoist.lvl1.mm","soilmoist.all.mm")))],scale=T)
#                       PC1          PC2         PC3         PC4         PC5          PC6
# tmean.gy.c         0.08084123  0.942945191  0.09699542  0.22594441  0.20933746  0.006156666
# ppt.gy.mm          0.46537136  0.136683503  0.13610779 -0.53502678 -0.29889495  0.608793844
# cmi.gy.mm          0.47421045  0.055304613  0.16471261 -0.36378724 -0.09306079 -0.777131995
# vpd.gy.max        -0.40060054  0.258417605 -0.68980547 -0.51693791 -0.12837637 -0.114902898
# soilmoist.lvl1.mm  0.43885823  0.009095837 -0.53418600  0.50928174 -0.51197628 -0.021873057
# soilmoist.all.mm   0.44602836 -0.149167668 -0.42866076 -0.06133211  0.76130094  0.108244664
# from loadings it looks like PC1 = water availability
#                             PC2 = Temp
#                             PC3 = evap demand (VPD + surface soil.)
# Importance of components:
#                          PC1    PC2     PC3     PC4     PC5     PC6
# Standard deviation     2.0734 1.0385 0.65791 0.34004 0.26509 0.06161
# Proportion of Variance 0.7165 0.1797 0.07214 0.01927 0.01171 0.00063
# Cumulative Proportion  0.7165 0.8962 0.96838 0.98766 0.99937 1.00000
biplot(climpca) #yup. same interpretation

climpcaPSME <- prcomp(traits.common[which(traits.common$SP.ID=="PSEMEN"),c(grep("gy", colnames(traits.common)), which(colnames(traits.common) %in% c("soilmoist.lvl1.mm","soilmoist.all.mm")))],scale=T)
  # note: no real problems with skewness that might screw things up with a PCA
  # almost identical results to the above, except PC3 is pretty much all soilmoist.lvl1 and doesn't have quite as much of the variance




## Quick visualization of PCs vs traits in PSME
pairs(traits.common[which(traits.common$SP.ID=="PSEMEN"),c("AG_TBIO","AG_TGROWTH",'SLA_HSA', "LEAF_CARBON","LEAF_NITROGEN","LEAF_CN","LEAF_LIFE","climPC1","climPC2","climPC3")], diag.panel = panel.hist, upper.panel=panel.smooth, lower.panel = function(x, y, digits=2, prefix="", cex.cor = 7) {
  panel.cor(x, y, digits, prefix, cex.cor)})

## PIPO
pairs(traits.common[which(traits.common$SP.ID=="PINPON"),c("AG_TBIO","AG_TGROWTH",'SLA_HSA', "LEAF_CARBON","LEAF_NITROGEN","LEAF_CN","LEAF_LIFE","climPC1","climPC2","climPC3")], diag.panel = panel.hist, upper.panel=panel.smooth, lower.panel = function(x, y, digits=2, prefix="", cex.cor = 7) {
  panel.cor(x, y, digits, prefix, cex.cor)})

# TSHE
pairs(traits.common[which(traits.common$SP.ID=="TSUHET"),c("AG_TBIO","AG_TGROWTH",'SLA_HSA', "LEAF_CARBON","LEAF_NITROGEN","LEAF_CN","LEAF_LIFE","climPC1","climPC2","climPC3")], diag.panel = panel.hist, upper.panel=panel.smooth, lower.panel = function(x, y, digits=2, prefix="", cex.cor = 7) {
  panel.cor(x, y, digits, prefix, cex.cor)})

# ABICON
pairs(traits.common[which(traits.common$SP.ID=="ABICON"),c("AG_TBIO","AG_TGROWTH",'SLA_HSA', "LEAF_CARBON","LEAF_NITROGEN","LEAF_CN","LEAF_LIFE","climPC1","climPC2","climPC3")], diag.panel = panel.hist, upper.panel=panel.smooth, lower.panel = function(x, y, digits=2, prefix="", cex.cor = 7) {
  panel.cor(x, y, digits, prefix, cex.cor)})

# ABIGRAN
pairs(traits.common[which(traits.common$SP.ID=="ABIGRA"),c("AG_TBIO","AG_TGROWTH",'SLA_HSA', "LEAF_CARBON","LEAF_NITROGEN","LEAF_CN","LEAF_LIFE","climPC1","climPC2","climPC3")], diag.panel = panel.hist, upper.panel=panel.smooth, lower.panel = function(x, y, digits=2, prefix="", cex.cor = 7) {
  panel.cor(x, y, digits, prefix, cex.cor)})




## trying a quick variance decomposition ##


### The first, over the top modeling attempt, that is definitely not the best solution.
SLA1 <- lmer(SLA_HSAsc ~ LAI_O + HEIGHTC_m + tmean.gy.c  + vpd.gy.max + soilmoist.lvl1.mm + (1|GENUS) + (LAI_O|GE.SP) + (HEIGHTC_m|GE.SP) + (tmean.gy.c|GE.SP) + (vpd.gy.max|GE.SP) + (soilmoist.lvl1.mm|GE.SP), traits.common)
  # this model often fails to converge. the current output is stored in the environment (as of 3/17/16)
  
nlmeControl(maxIter=150, msMaxIter = 100, pnlsMaxIter = 10, msVerbose = TRUE)
SLA1nlme <- lme(SLA_HSAsc ~ LAI_O + HEIGHTC_m + tmean.gy.c  + vpd.gy.max , random= ~ LAI_O + HEIGHTC_m + tmean.gy.c  + vpd.gy.max | GE.SP, data=traits.common.narm, control=nlmeControl(maxIter=150, msMaxIter = 100, pnlsMaxIter = 10, msVerbose = TRUE))




#-----------------------------------------------------------------------------
################# VARIANCE DECOMPOSITION by fitting only random effects ####################
#-----------------------------------------------------------------------------



SLAvar <- lmer(SLA_HSA~ (1|PLOT_ID) + (1|GENUS) + (1|GE.SP), traits.common)
SLAvar2 <- lmer(SLA_HSA~ (1|GENUS/SPECIES/PLOT_ID), traits.common)
  ## These two are slightly different formulations??? but they give roughly the same answers.

LeafLifevar <- lmer(LEAF_LIFE~ (1|PLOT_ID) + (1|GENUS) + (1|GE.SP), traits.common)
LeafLifevar2 <- lmer(LEAF_LIFE~ (1|GENUS/SPECIES/PLOT_ID), traits.common)

## probably need better models for these.
LeafCNvar <- lmer(LEAF_CN~ (1|PLOT_ID) + (1|GENUS) + (1|GE.SP), traits.common)
LeafCNvar2 <- lmer(LEAF_CN~ (1|GENUS/SPECIES/PLOT_ID), traits.common)

LeafCarbvar <- lmer(LEAF_CARBON~ (1|PLOT_ID) + (1|GENUS) + (1|GE.SP), traits.common)
LeafCarbvar2 <- lmer(LEAF_CARBON~ (1|GENUS/SPECIES/PLOT_ID), traits.common)

LeafNitvar <- lmer(LEAF_NITROGEN~ (1|PLOT_ID) + (1|GENUS) + (1|GE.SP), traits.common)
LeafNitvar2 <- lmer(LEAF_NITROGEN~ (1|GENUS/SPECIES/PLOT_ID), traits.common)

SLAvariance <- data.frame(VarCorr(SLAvar2))
LeafLifevariance <- data.frame(VarCorr(LeafLifevar2))
LeafCNvariance <- data.frame(VarCorr(LeafCNvar2))
LeafCarvariance <- data.frame(VarCorr(LeafCarbvar2))
LeafNitvariance <- data.frame(VarCorr(LeafNitvar2))

traitvars <- data.frame(SLAvariance[,4], LeafLifevariance[,4], LeafCNvariance[,4], LeafCarvariance[,4], LeafNitvariance[,4])
colnames(traitvars) <- c("SLA", "LEAF_Life", "LEAF_CN", "LEAF_CARBON", "LEAF_NITROGEN")
rownames(traitvars) <- c("BtwPlot", "BtwSpecies", "BtwGenus", "WtinPlot")

traitvars_scaled <- traitvars  
for(i in 1:ncol(traitvars)){
traitvars_scaled[,i] <- traitvars[,i]/sum(traitvars[,i])
}

## making stacked bar plots of variance decomposition:
quartz(width=5, height=4)
cols <- brewer.pal(11, "RdBu")[c(1,11,9, 6)]
barplot(as.matrix(traitvars_scaled),beside=F,legend.text = F,xpd = T, names.arg = c("SLA", "Life","C/N","%C","%N"),args.legend = list(x=4, y=1.3, ncol=2), col = paste0(cols,"CC"), ylab="Proportion of total Variance")
legend(xpd=T, x = 1, y=1.3, legend=c("W/in Plot", "Btw Plots", "Btw Spp", "Btw Genera"), fill=paste0(cols,"CC")[c(4,1,2,3)], ncol=2, bty="n",  cex=1.2)



######## making stacked bar plots with intra-specific variation highlighted
traitvars_scaled2 <- traitvars_scaled[c(2,3,1,4),]
## making stacked bar plots of variance decomposition:
quartz(width=5, height=4)
cols <- brewer.pal(11, "RdBu")[c(11,9,1, 6)]
barplot(as.matrix(traitvars_scaled2),beside=F,legend.text = F,xpd = T, names.arg = c("SLA", "Leaf\nLife","C/N","%C","%N"),args.legend = list(x=4, y=1.3, ncol=2), col = paste0(cols,"CC"), ylab="Proportion of total Variance", xlab="Traits.common")
legend(xpd=T, x = 1, y=1.3, legend=c("W/in Plot", "Btw Plots", "Btw Spp", "Btw Genera"), fill=paste0(cols,"CC")[c(4,3,1,2)], ncol=2, bty="n",  cex=1.2)

## dropping C and N solo
# traitvars_scaled3 <- traitvars_scaled[c(3,2,1,4),]
# quartz(width=4.5, height=3.5)
# par(mar=c(3,4,1,6.5), mgp=c(2.5,1,0))
# cols <- brewer.pal(11, "RdBu")[c(11,9,1, 6)]
# barplot(as.matrix(traitvars_scaled3[,1:3]),beside=F,legend.text = F,xpd = T, names.arg = c("SLA", "Leaf\nLife","C/N"),args.legend = list(x=4, y=1.3, ncol=2), col = paste0(cols,"CC"), ylab="Proportion of total Variance")
# legend(xpd=T, x = 3.6, y=.7, legend=c("W/in Plot", "Btw Plots", "Btw Spp", "Btw Genera"), fill=paste0(cols,"CC")[c(4,3,2,1)], ncol=1, bty="n",  cex=1)






#-----------------------------------------------------------------------------
################ VARIANCE DECOMPOSITION of only common needle leafed conifers ####################
#-----------------------------------------------------------------------------

conifers <- names(which(xtabs(~SPP_O1_ABBREV, biomass)>1))
traits.conifers <- subset(traits, subset=SP.ID %in% conifers) # 1099 trait measurements
conifers <- c("Abies","Picea","Pinus","Pseudotsuga","Thuja","Tsuga") # 1065 trait measurements
traits.conifers <- subset(traits, subset=GENUS %in% conifers)

traits.conifers$SP.ID <- factor(traits.conifers$SP.ID)

dominants <- levels(biomass$SPP_O1_ABBREV)
traits.dominants <- subset(traits, subset=SP.ID %in% dominants)



variance.decomp <- function(dataz, lab){
  
  #### logged trait values
  SLAvar2 <- lmer(SLA_HSA~ (1|GENUS/SPECIES/PLOT_ID), dataz)
  LMAvar2 <- lmer(log.LMA~ (1|GENUS/SPECIES/PLOT_ID), dataz)
  LeafLifevar2 <- lmer(log.LL~ (1|GENUS/SPECIES/PLOT_ID), dataz)
  LeafCNvar2 <- lmer(LEAF_CN~ (1|GENUS/SPECIES/PLOT_ID), dataz)
  LeafCarbvar2 <- lmer(LEAF_CARBON~ (1|GENUS/SPECIES/PLOT_ID), dataz)
  LeafNitvar2 <- lmer(log.Nmass~ (1|GENUS/SPECIES/PLOT_ID), dataz)
  Nareavar2 <- lmer(log.Narea~ (1|GENUS/SPECIES/PLOT_ID), dataz)
  
  SLAvariance <- data.frame(VarCorr(SLAvar2))
  LMAvariance <- data.frame(VarCorr(LMAvar2))
  LeafLifevariance <- data.frame(VarCorr(LeafLifevar2))
  LeafCNvariance <- data.frame(VarCorr(LeafCNvar2))
  LeafCarvariance <- data.frame(VarCorr(LeafCarbvar2))
  LeafNitvariance <- data.frame(VarCorr(LeafNitvar2))
  Nareavariance <- data.frame(VarCorr(Nareavar2))
  
  traitvars <- data.frame(LMAvariance[,4], LeafLifevariance[,4], LeafNitvariance[,4], Nareavariance[,4], LeafCarvariance[,4], LeafCNvariance[,4])
  colnames(traitvars) <- c("log.LMA", "log.LL","log.Nmass","log.Narea","Cmass", "LEAF_CN")
  rownames(traitvars) <- c("BtwPlot", "BtwSpecies", "BtwGenus", "WtinPlot")
  traitvars_scaled <- traitvars  
  for(i in 1:ncol(traitvars)){
    traitvars_scaled[,i] <- traitvars[,i]/sum(traitvars[,i])
  }
  # reorder appropriately so that it plots with w.inPlots on top and btw Genera on bottom.
  traitvars_scaled1 <- traitvars_scaled[c(3,2,1,4),]
  
  
  quartz(width=5, height=4, "logged")
  cols <- brewer.pal(11, "PRGn")[c(1,3,9,11)]
  barplot(as.matrix(traitvars_scaled1),beside=F,legend.text = F,xpd = T, names.arg = c("log\n(LMA)", "log\n(LL)","log\n(Nm)","log\n(Na)", "%C","C/N"),args.legend = list(x=4, y=1.3, ncol=2), col = paste0(cols,"CC"), ylab="Proportion of total Variance", xlab=lab)
  # almost works: expression(paste("log\n","(",N[mass],")", sep="")) # except it offsets the upper and lower lines
  legend(xpd=T, x = 1, y=1.3, legend=c("W/in Plot", "Btw Plots", "Btw Spp", "Btw Genera"), fill=paste0(cols,"CC")[c(4,3,2,1)], ncol=2, bty="n",  cex=1.2)
  
  
  #### raw values
  rLMAvar2 <- lmer(LMA~ (1|GENUS/SPECIES/PLOT_ID), dataz)
  rLeafLifevar2 <- lmer(LEAF_LIFE~ (1|GENUS/SPECIES/PLOT_ID), dataz)
  rLeafCNvar2 <- lmer(LEAF_CN~ (1|GENUS/SPECIES/PLOT_ID), dataz)
  rLeafCarbvar2 <- lmer(LEAF_CARBON~ (1|GENUS/SPECIES/PLOT_ID), dataz)
  rLeafNitvar2 <- lmer(LEAF_NITROGEN~ (1|GENUS/SPECIES/PLOT_ID), dataz)
  rNareavar2 <- lmer(Narea~ (1|GENUS/SPECIES/PLOT_ID), dataz)
  
  #rSLAvariance <- data.frame(VarCorr(rSLAvar2))
  rLMAvariance <- data.frame(VarCorr(rLMAvar2))
  rLeafLifevariance <- data.frame(VarCorr(rLeafLifevar2))
  rLeafCNvariance <- data.frame(VarCorr(rLeafCNvar2))
  rLeafCarvariance <- data.frame(VarCorr(rLeafCarbvar2))
  rLeafNitvariance <- data.frame(VarCorr(rLeafNitvar2))
  rNareavariance <- data.frame(VarCorr(rNareavar2))
  
  rtraitvars <- data.frame(rLMAvariance[,4], rLeafLifevariance[,4], rLeafNitvariance[,4], rNareavariance[,4], rLeafCarvariance[,4], rLeafCNvariance[,4])
  colnames(rtraitvars) <- c("LMA", "LEAF_Life", "Nmass","Narea","Cmass","LEAF_CN")
  rownames(rtraitvars) <- c("BtwPlot", "BtwSpecies", "BtwGenus", "WtinPlot")
  rtraitvars_scaled <- rtraitvars  
  for(i in 1:ncol(rtraitvars)){
    rtraitvars_scaled[,i] <- rtraitvars[,i]/sum(rtraitvars[,i])
  }
  rtraitvars_scaled1 <- rtraitvars_scaled[c(3,2,1,4),]
  
  quartz(width=5, height=4, "raw")
  cols <- brewer.pal(11, "PRGn")[c(1,3,9,11)]
  barplot(as.matrix(rtraitvars_scaled1),beside=F,legend.text = F,xpd = T, names.arg = c("LMA", "Leaf\nLife",expression(N[mass]),expression(N[area]), "%C","C/N"),args.legend = list(x=4, y=1.3, ncol=2), col = paste0(cols,"CC"), ylab="Proportion of total Variance", xlab=lab)
  legend(xpd=T, x = 1, y=1.3, legend=c("W/in Plot", "Btw Plots", "Btw Spp", "Btw Genera"), fill=paste0(cols,"CC")[c(4,3,2,1)], ncol=2, bty="n",  cex=1.2)
  
  return(list(raw = rtraitvars_scaled1, logged = traitvars_scaled1))
  
}

variance.decomp(dataz=traits.common, lab="traits.common")
variance.decomp(dataz=traits.dominants, lab="traits.dominants")
variance.decomp(dataz=traits.conifers, lab="traits.conifers")
traits.domcon <- traits.dominants[-which(traits.dominants$SP.ID=="QUECHR"),]
domconvars <- variance.decomp(dataz=traits.domcon, lab="traits.domcon")

SLAvar <- lmer(SLA_HSA~ (1|PLOT_ID) + (1|GENUS) + (1|GE.SP), traits.conifers)
SLAvar2 <- lmer(SLA_HSA~ (1|GENUS/SPECIES/PLOT_ID), traits.conifers)
## These two are slightly different formulations??? but they give roughly the same answers.
LMAvar <- lmer(LMA~ (1|PLOT_ID) + (1|GENUS) + (1|GE.SP), traits.conifers) # uses LMA_HSA rather than PSA
LMAvar2 <- lmer(LMA~ (1|GENUS/SPECIES/PLOT_ID), traits.conifers)


LeafLifevar <- lmer(LEAF_LIFE~ (1|PLOT_ID) + (1|GENUS) + (1|GE.SP), traits.conifers)
LeafLifevar2 <- lmer(LEAF_LIFE~ (1|GENUS/SPECIES/PLOT_ID), traits.conifers)

## probably need better models for these.
LeafCNvar <- lmer(LEAF_CN~ (1|PLOT_ID) + (1|GENUS) + (1|GE.SP), traits.conifers)
LeafCNvar2 <- lmer(LEAF_CN~ (1|GENUS/SPECIES/PLOT_ID), traits.conifers)

LeafCarbvar <- lmer(LEAF_CARBON~ (1|PLOT_ID) + (1|GENUS) + (1|GE.SP), traits.conifers)
LeafCarbvar2 <- lmer(LEAF_CARBON~ (1|GENUS/SPECIES/PLOT_ID), traits.conifers)

LeafNitvar <- lmer(LEAF_NITROGEN~ (1|PLOT_ID) + (1|GENUS) + (1|GE.SP), traits.conifers)
LeafNitvar2 <- lmer(LEAF_NITROGEN~ (1|GENUS/SPECIES/PLOT_ID), traits.conifers)

SLAvariance <- data.frame(VarCorr(SLAvar2))
LMAvariance <- data.frame(VarCorr(LMAvar2))
LeafLifevariance <- data.frame(VarCorr(LeafLifevar2))
LeafCNvariance <- data.frame(VarCorr(LeafCNvar2))
LeafCarvariance <- data.frame(VarCorr(LeafCarbvar2))
LeafNitvariance <- data.frame(VarCorr(LeafNitvar2))


traitvars <- data.frame(SLAvariance[,4], LeafLifevariance[,4], LeafCNvariance[,4], LeafCarvariance[,4], LeafNitvariance[,4])
colnames(traitvars) <- c("SLA", "LEAF_Life", "LEAF_CN", "LEAF_CARBON", "LEAF_NITROGEN")
rownames(traitvars) <- c("BtwPlot", "BtwSpecies", "BtwGenus", "WtinPlot")

## dataframe with LMA
traitvars2 <- data.frame(LMAvariance[,4], LeafLifevariance[,4], LeafCNvariance[,4], LeafCarvariance[,4], LeafNitvariance[,4])
colnames(traitvars2) <- c("LMA", "LEAF_Life", "LEAF_CN", "LEAF_CARBON", "LEAF_NITROGEN")
rownames(traitvars2) <- c("BtwPlot", "BtwSpecies", "BtwGenus", "WtinPlot")


traitvars_scaled <- traitvars  
for(i in 1:ncol(traitvars)){
  traitvars_scaled[,i] <- traitvars[,i]/sum(traitvars[,i])
}

traitvars_scaled2 <- traitvars2  
for(i in 1:ncol(traitvars2)){
  traitvars_scaled2[,i] <- traitvars2[,i]/sum(traitvars2[,i])
}

######## making stacked bar plots with intra-specific variation highlighted
traitvars_scaled <- traitvars_scaled[c(3,2,1,4),]
traitvars_scaled2 <- traitvars_scaled2[c(3,2,1,4),]

## making stacked bar plots of variance decomposition:
quartz(width=5, height=4)
cols <- brewer.pal(11, "RdBu")[c(11,9,1, 6)]
barplot(as.matrix(traitvars_scaled),beside=F,legend.text = F,xpd = T, names.arg = c("SLA", "Leaf\nLife","C/N","%C","%N"),args.legend = list(x=4, y=1.3, ncol=2), col = paste0(cols,"CC"), ylab="Proportion of total Variance", xlab="dominants")
legend(xpd=T, x = 1, y=1.3, legend=c("W/in Plot", "Btw Plots", "Btw Spp", "Btw Genera"), fill=paste0(cols,"CC")[c(4,3,1,2)], ncol=2, bty="n",  cex=1.2)


quartz(width=5, height=4, )
cols <- brewer.pal(11, "PRGn")[c(1,3,9,11)]
barplot(as.matrix(traitvars_scaled2),beside=F,legend.text = F,xpd = T, names.arg = c("LMA", "Leaf\nLife","C/N","%C","%N"),args.legend = list(x=4, y=1.3, ncol=2), col = paste0(cols,"CC"), ylab="Proportion of total Variance", xlab="dominants")
legend(xpd=T, x = 1, y=1.3, legend=c("W/in Plot", "Btw Plots", "Btw Spp", "Btw Genera"), fill=paste0(cols,"CC")[c(4,3,2,1)], ncol=2, bty="n",  cex=1.2)


###### showing how much of trait space these conifers still cover
quartz(width=4, height=4)
par(mar=c(3,3,1,1))
plot(log.LL~log.LMA, LES, col="grey", pch=16, ylab="log Leaf Lifespan (months)", xlab="log LMA (g/m^2)")
points(log.LL~log.LMA, traits.conifers, col=SP.ID)

quartz(width=2, height=2)
par(mar=c(3,3,1,1))

plot(LL~LMA, LES, col="grey", pch=16, ylab="log Leaf Lifespan (months)", xlab="log LMA (g/m^2)", xlim=c(0,800))
points(LLmonths~LMA, traits.conifers, col=SP.ID)






#-----------------------------------------------------------------------------
############### VARIANCE DECOMPOSITION of only dominant species in the biomass dataset ####################
#-----------------------------------------------------------------------------

### 
dominants <- levels(biomass$SPP_O1_ABBREV)
traits.dominants <- subset(traits, subset=SP.ID %in% dominants)

SLAvar <- lmer(SLA_HSA~ (1|PLOT_ID) + (1|GENUS) + (1|GE.SP), traits.dominants)
SLAvar2 <- lmer(SLA_HSA~ (1|GENUS/SPECIES/PLOT_ID), traits.dominants)
## These two are slightly different formulations??? but they give roughly the same answers.
rLMAvar <- lmer(LMA~ (1|PLOT_ID) + (1|GENUS) + (1|GE.SP), traits.dominants) # uses LMA_HSA rather than PSA
rLMAvar2 <- lmer(LMA~ (1|GENUS/SPECIES/PLOT_ID), traits.dominants)
LMAvar <- lmer(LMA~ (1|PLOT_ID) + (1|GENUS) + (1|GE.SP), traits.dominants) # uses LMA_HSA rather than PSA
LMAvar2 <- lmer(LMA~ (1|GENUS/SPECIES/PLOT_ID), traits.dominants)


rLeafLifevar <- lmer(LEAF_LIFE~ (1|PLOT_ID) + (1|GENUS) + (1|GE.SP), traits.dominants)
rLeafLifevar2 <- lmer(LEAF_LIFE~ (1|GENUS/SPECIES/PLOT_ID), traits.dominants)

## probably need better models for these.
rLeafCNvar <- lmer(LEAF_CN~ (1|PLOT_ID) + (1|GENUS) + (1|GE.SP), traits.dominants)
rLeafCNvar2 <- lmer(LEAF_CN~ (1|GENUS/SPECIES/PLOT_ID), traits.dominants)

rLeafCarbvar <- lmer(LEAF_CARBON~ (1|PLOT_ID) + (1|GENUS) + (1|GE.SP), traits.dominants)
rLeafCarbvar2 <- lmer(LEAF_CARBON~ (1|GENUS/SPECIES/PLOT_ID), traits.dominants)

rLeafNitvar <- lmer(LEAF_NITROGEN~ (1|PLOT_ID) + (1|GENUS) + (1|GE.SP), traits.dominants)
rLeafNitvar2 <- lmer(LEAF_NITROGEN~ (1|GENUS/SPECIES/PLOT_ID), traits.dominants)

rNareavar2 <- lmer(Narea~ (1|GENUS/SPECIES/PLOT_ID), traits.dominants)

rSLAvariance <- data.frame(VarCorr(rSLAvar2))
rLMAvariance <- data.frame(VarCorr(rLMAvar2))
rLeafLifevariance <- data.frame(VarCorr(rLeafLifevar2))
rLeafCNvariance <- data.frame(VarCorr(rLeafCNvar2))
rLeafCarvariance <- data.frame(VarCorr(rLeafCarbvar2))
rLeafNitvariance <- data.frame(VarCorr(rLeafNitvar2))
rNareavariance <- data.frame(VarCorr(rNareavar2))


## dataframe with SLA
rtraitvars <- data.frame(rSLAvariance[,4], rLeafLifevariance[,4], rLeafCNvariance[,4], rLeafCarvariance[,4], rLeafNitvariance[,4])
colnames(rtraitvars) <- c("SLA", "LEAF_Life", "LEAF_CN", "LEAF_CARBON", "LEAF_NITROGEN")
rownames(rtraitvars) <- c("BtwPlot", "BtwSpecies", "BtwGenus", "WtinPlot")
## dataframe with LMA
rtraitvars2 <- data.frame(rLMAvariance[,4], rLeafLifevariance[,4], rLeafCNvariance[,4], rLeafCarvariance[,4], rLeafNitvariance[,4], rNareavariance[,4])
colnames(rtraitvars2) <- c("LMA", "LEAF_Life", "LEAF_CN", "LEAF_CARBON", "LEAF_NITROGEN", "Narea")
rownames(rtraitvars2) <- c("BtwPlot", "BtwSpecies", "BtwGenus", "WtinPlot")


rtraitvars_scaled <- rtraitvars  
for(i in 1:ncol(rtraitvars)){
  rtraitvars_scaled[,i] <- rtraitvars[,i]/sum(rtraitvars[,i])
}

rtraitvars_scaled2 <- rtraitvars2  
for(i in 1:ncol(rtraitvars2)){
  rtraitvars_scaled2[,i] <- rtraitvars2[,i]/sum(rtraitvars2[,i])
}

######## making stacked bar plots with intra-specific variation highlighted
rtraitvars_scaled <- traitvars_scaled[c(3,2,1,4),]
rtraitvars_scaled2 <- traitvars_scaled2[c(3,2,1,4),]


####### Logged trait variance decomp of dominant conifers ####
SLAvar <- lmer(SLA_HSA~ (1|PLOT_ID) + (1|GENUS) + (1|GE.SP), traits.dominants)
SLAvar2 <- lmer(SLA_HSA~ (1|GENUS/SPECIES/PLOT_ID), traits.dominants)
## These two are slightly different formulations??? but they give roughly the same answers.
LMAvar <- lmer(log.LMA~ (1|PLOT_ID) + (1|GENUS) + (1|GE.SP), traits.dominants) # uses LMA_HSA rather than PSA
LMAvar2 <- lmer(log.LMA~ (1|GENUS/SPECIES/PLOT_ID), traits.dominants)


LeafLifevar <- lmer(log.LL~ (1|PLOT_ID) + (1|GENUS) + (1|GE.SP), traits.dominants)
LeafLifevar2 <- lmer(log.LL~ (1|GENUS/SPECIES/PLOT_ID), traits.dominants)

## probably need better models for these.
LeafCNvar <- lmer(LEAF_CN~ (1|PLOT_ID) + (1|GENUS) + (1|GE.SP), traits.dominants)
LeafCNvar2 <- lmer(LEAF_CN~ (1|GENUS/SPECIES/PLOT_ID), traits.dominants)

LeafCarbvar <- lmer(LEAF_CARBON~ (1|PLOT_ID) + (1|GENUS) + (1|GE.SP), traits.dominants)
LeafCarbvar2 <- lmer(LEAF_CARBON~ (1|GENUS/SPECIES/PLOT_ID), traits.dominants)

LeafNitvar <- lmer(log.Nmass~ (1|PLOT_ID) + (1|GENUS) + (1|GE.SP), traits.dominants)
LeafNitvar2 <- lmer(log.Nmass~ (1|GENUS/SPECIES/PLOT_ID), traits.dominants)

Nareavar <- lmer(log.Narea~ (1|PLOT_ID) + (1|GENUS) + (1|GE.SP), traits.dominants)
Nareavar2 <- lmer(log.Narea~ (1|GENUS/SPECIES/PLOT_ID), traits.dominants)


SLAvariance <- data.frame(VarCorr(SLAvar2))
LMAvariance <- data.frame(VarCorr(LMAvar2))
LeafLifevariance <- data.frame(VarCorr(LeafLifevar2))
LeafCNvariance <- data.frame(VarCorr(LeafCNvar2))
LeafCarvariance <- data.frame(VarCorr(LeafCarbvar2))
LeafNitvariance <- data.frame(VarCorr(LeafNitvar2))
Nareavariance <- data.frame(VarCorr(Nareavar2))

## dataframe with SLA
traitvars <- data.frame(SLAvariance[,4], LeafLifevariance[,4], LeafCNvariance[,4], LeafCarvariance[,4], LeafNitvariance[,4])
colnames(traitvars) <- c("SLA", "LEAF_Life", "LEAF_CN", "LEAF_CARBON", "LEAF_NITROGEN")
rownames(traitvars) <- c("BtwPlot", "BtwSpecies", "BtwGenus", "WtinPlot")
## dataframe with LMA
traitvars2 <- data.frame(LMAvariance[,4], LeafLifevariance[,4], LeafCNvariance[,4], LeafCarvariance[,4], LeafNitvariance[,4], Nareavariance[,4])
colnames(traitvars2) <- c("LMA", "LEAF_Life", "LEAF_CN", "LEAF_CARBON", "LEAF_NITROGEN","Narea")
rownames(traitvars2) <- c("BtwPlot", "BtwSpecies", "BtwGenus", "WtinPlot")


traitvars_scaled <- traitvars  
for(i in 1:ncol(traitvars)){
  traitvars_scaled[,i] <- traitvars[,i]/sum(traitvars[,i])
}

traitvars_scaled2 <- traitvars2  
for(i in 1:ncol(traitvars2)){
  traitvars_scaled2[,i] <- traitvars2[,i]/sum(traitvars2[,i])
}

######## making stacked bar plots with intra-specific variation highlighted
traitvars_scaled <- traitvars_scaled[c(3,2,1,4),]
traitvars_scaled2 <- traitvars_scaled2[c(3,2,1,4),]

#..........



## RAW TRAITS making stacked bar plots of variance decomposition:

quartz(width=5, height=4)
cols <- brewer.pal(11, "RdBu")[c(11,9,1, 6)]
barplot(as.matrix(rtraitvars_scaled),beside=F,legend.text = F,xpd = T, names.arg = c("SLA", "Leaf\nLife","C/N","%C","%N"),args.legend = list(x=4, y=1.3, ncol=2), col = paste0(cols,"CC"), ylab="Proportion of total Variance", xlab="dominants")
legend(xpd=T, x = 1, y=1.3, legend=c("W/in Plot", "Btw Plots", "Btw Spp", "Btw Genera"), fill=paste0(cols,"CC")[c(4,3,1,2)], ncol=2, bty="n",  cex=1.2)


quartz(width=5, height=4, title = "Raw traits")
cols <- brewer.pal(11, "PRGn")[c(1,3,9,11)]
barplot(as.matrix(rtraitvars_scaled2),beside=F,legend.text = F,xpd = T, names.arg = c("LMA", "Leaf\nLife","C/N","%C","%N","Narea"),args.legend = list(x=4, y=1.3, ncol=2), col = paste0(cols,"CC"), ylab="Proportion of total Variance", xlab="dominants")
legend(xpd=T, x = 1, y=1.3, legend=c("W/in Plot", "Btw Plots", "Btw Spp", "Btw Genera"), fill=paste0(cols,"CC")[c(4,3,2,1)], ncol=2, bty="n",  cex=1.2)



## Logged traits making stacked bar plots of variance decomposition:
quartz(width=5, height=4)
cols <- brewer.pal(11, "RdBu")[c(11,9,1, 6)]
barplot(as.matrix(traitvars_scaled),beside=F,legend.text = F,xpd = T, names.arg = c("SLA", "Leaf\nLife","C/N","%C","%N"),args.legend = list(x=4, y=1.3, ncol=2), col = paste0(cols,"CC"), ylab="Proportion of total Variance", xlab="dominants")
legend(xpd=T, x = 1, y=1.3, legend=c("W/in Plot", "Btw Plots", "Btw Spp", "Btw Genera"), fill=paste0(cols,"CC")[c(4,3,1,2)], ncol=2, bty="n",  cex=1.2)


quartz(width=5, height=4, title="logged traits")
cols <- brewer.pal(11, "PRGn")[c(1,3,9,11)]
barplot(as.matrix(traitvars_scaled2),beside=F,legend.text = F,xpd = T, names.arg = c("LMA", "Leaf\nLife","C/N","%C","%N","Narea"),args.legend = list(x=4, y=1.3, ncol=2), col = paste0(cols,"CC"), ylab="Proportion of total Variance", xlab="dominants")
legend(xpd=T, x = 1, y=1.3, legend=c("W/in Plot", "Btw Plots", "Btw Spp", "Btw Genera"), fill=paste0(cols,"CC")[c(4,3,2,1)], ncol=2, bty="n",  cex=1.2)



#-----------------------------------------------------------------------------
##################### VARIANCE DECOMPOSITION of entire dataset ####################
#-----------------------------------------------------------------------------


SLAvar <- lmer(SLA_HSA~ (1|PLOT_ID) + (1|GENUS) + (1|GE.SP), traits)
SLAvar2 <- lmer(SLA_HSA~ (1|GENUS/SPECIES/PLOT_ID), traits)
## These two are slightly different formulations??? but they give roughly the same answers.

LeafLifevar <- lmer(LEAF_LIFE~ (1|PLOT_ID) + (1|GENUS) + (1|GE.SP), traits)
LeafLifevar2 <- lmer(LEAF_LIFE~ (1|GENUS/SPECIES/PLOT_ID), traits)

## probably need better models for these.
LeafCNvar <- lmer(LEAF_CN~ (1|PLOT_ID) + (1|GENUS) + (1|GE.SP), traits)
LeafCNvar2 <- lmer(LEAF_CN~ (1|GENUS/SPECIES/PLOT_ID), traits)

LeafCarbvar <- lmer(LEAF_CARBON~ (1|PLOT_ID) + (1|GENUS) + (1|GE.SP), traits)
LeafCarbvar2 <- lmer(LEAF_CARBON~ (1|GENUS/SPECIES/PLOT_ID), traits)

LeafNitvar <- lmer(LEAF_NITROGEN~ (1|PLOT_ID) + (1|GENUS) + (1|GE.SP), traits)
LeafNitvar2 <- lmer(LEAF_NITROGEN~ (1|GENUS/SPECIES/PLOT_ID), traits)

SLAvariance <- data.frame(VarCorr(SLAvar2))
LeafLifevariance <- data.frame(VarCorr(LeafLifevar2))
LeafCNvariance <- data.frame(VarCorr(LeafCNvar2))
LeafCarvariance <- data.frame(VarCorr(LeafCarbvar2))
LeafNitvariance <- data.frame(VarCorr(LeafNitvar2))

traitvars <- data.frame(SLAvariance[,4], LeafLifevariance[,4], LeafCNvariance[,4], LeafCarvariance[,4], LeafNitvariance[,4])
colnames(traitvars) <- c("SLA", "LEAF_Life", "LEAF_CN", "LEAF_CARBON", "LEAF_NITROGEN")
rownames(traitvars) <- c("BtwPlot", "BtwSpecies", "BtwGenus", "WtinPlot")

traitvars_scaled <- traitvars  
for(i in 1:ncol(traitvars)){
  traitvars_scaled[,i] <- traitvars[,i]/sum(traitvars[,i])
}

## making stacked bar plots of variance decomposition:
traitvars_scaled2 <- traitvars_scaled[c(2,3,1,4),]
## making stacked bar plots of variance decomposition:
quartz(width=5, height=4)
cols <- brewer.pal(11, "RdBu")[c(11,9,1, 6)]
barplot(as.matrix(traitvars_scaled2),beside=F,legend.text = F,xpd = T, names.arg = c("SLA", "Leaf\nLife","C/N","%C","%N"),args.legend = list(x=4, y=1.3, ncol=2), col = paste0(cols,"CC"), ylab="Proportion of total Variance", xlab="Entire Dataset")
legend(xpd=T, x = 1, y=1.3, legend=c("W/in Plot", "Btw Plots", "Btw Spp", "Btw Genera"), fill=paste0(cols,"CC")[c(4,3,1,2)], ncol=2, bty="n",  cex=1.2)



### VARIANCE DECOMPOSITION of only common species and genera ####################
SLAvar <- lmer(SLA_HSA~ (1|PLOT_ID) + (1|GENUS) + (1|GE.SP), traits.commongenera)
SLAvar2 <- lmer(SLA_HSA~ (1|GENUS/SPECIES/PLOT_ID), traits.commongenera)
## These two are slightly different formulations??? but they give roughly the same answers.

LeafLifevar <- lmer(LEAF_LIFE~ (1|PLOT_ID) + (1|GENUS) + (1|GE.SP), traits.commongenera)
LeafLifevar2 <- lmer(LEAF_LIFE~ (1|GENUS/SPECIES/PLOT_ID), traits.commongenera)

## probably need better models for these.
LeafCNvar <- lmer(LEAF_CN~ (1|PLOT_ID) + (1|GENUS) + (1|GE.SP), traits.commongenera)
LeafCNvar2 <- lmer(LEAF_CN~ (1|GENUS/SPECIES/PLOT_ID), traits.commongenera)

LeafCarbvar <- lmer(LEAF_CARBON~ (1|PLOT_ID) + (1|GENUS) + (1|GE.SP), traits.commongenera)
LeafCarbvar2 <- lmer(LEAF_CARBON~ (1|GENUS/SPECIES/PLOT_ID), traits.commongenera)

LeafNitvar <- lmer(LEAF_NITROGEN~ (1|PLOT_ID) + (1|GENUS) + (1|GE.SP), traits.commongenera)
LeafNitvar2 <- lmer(LEAF_NITROGEN~ (1|GENUS/SPECIES/PLOT_ID), traits.commongenera)

SLAvariance <- data.frame(VarCorr(SLAvar2))
LeafLifevariance <- data.frame(VarCorr(LeafLifevar2))
LeafCNvariance <- data.frame(VarCorr(LeafCNvar2))
LeafCarvariance <- data.frame(VarCorr(LeafCarbvar2))
LeafNitvariance <- data.frame(VarCorr(LeafNitvar2))

traitvars <- data.frame(SLAvariance[,4], LeafLifevariance[,4], LeafCNvariance[,4], LeafCarvariance[,4], LeafNitvariance[,4])
colnames(traitvars) <- c("SLA", "LEAF_Life", "LEAF_CN", "LEAF_CARBON", "LEAF_NITROGEN")
rownames(traitvars) <- c("BtwPlot", "BtwSpecies", "BtwGenus", "WtinPlot")

traitvars_scaled <- traitvars  
for(i in 1:ncol(traitvars)){
  traitvars_scaled[,i] <- traitvars[,i]/sum(traitvars[,i])
}


traitvars_scaled2 <- traitvars_scaled[c(2,3,1,4),]
## making stacked bar plots of variance decomposition:
quartz(width=5, height=4)
cols <- brewer.pal(11, "RdBu")[c(11,9,1, 6)]
barplot(as.matrix(traitvars_scaled2),beside=F,legend.text = F,xpd = T, names.arg = c("SLA", "Leaf\nLife","C/N","%C","%N"),args.legend = list(x=4, y=1.3, ncol=2), col = paste0(cols,"CC"), ylab="Proportion of total Variance", xlab="Genera w/ replication")
legend(xpd=T, x = 1, y=1.3, legend=c("W/in Plot", "Btw Plots", "Btw Spp", "Btw Genera"), fill=paste0(cols,"CC")[c(4,3,1,2)], ncol=2, bty="n",  cex=1.2)







#-----------------------------------------------------------------------------
##################### VARIANCE DECOMPOSITION log(LES) traits ####################
#-----------------------------------------------------------------------------

# ## with the quick traits.common dataset. This looks similar, but not identical to the unlogged trait space decomp
# logLMAvar <- lmer(log.LMA~ (1|GENUS/SPECIES/PLOT_ID), traits.common)
# logLLvar <- lmer(log.LL~ (1|GENUS/SPECIES/PLOT_ID), traits.common)
# logNmassvar <- lmer(log.Nmass~ (1|GENUS/SPECIES/PLOT_ID), traits.common)
# # logLMAvar <- lmer(log.LMA~ (1|GENUS) + (1|SP.ID) + (1|PLOT_ID), traits.common)
# # logLLvar <- lmer(log.LL~  (1|GENUS) + (1|SP.ID) + (1|PLOT_ID), traits.common)
# # logNmassvar <- lmer(log.Nmass~  (1|GENUS) + (1|SP.ID) + (1|PLOT_ID), traits.common)
#   # note: this unnested formulation is probably getting plot ID wrong, because some spp have the same plots. But the results are pretty similar, just with the w/in plot bars increasing for LL and Nmass
# 
# LMAvariance <- data.frame(VarCorr(logLMAvar))
# LLvariance <- data.frame(VarCorr(logLLvar))
# Nmassvariance <- data.frame(VarCorr(logNmassvar))
# 
# traitvars <- data.frame(LMAvariance[,4], LLvariance[,4], Nmassvariance[,4])
# colnames(traitvars) <- c("logLMA", "logLL", "logNmass")
# rownames(traitvars) <- c("BtwPlot", "BtwSpecies", "BtwGenus", "WtinPlot")
# 
# traitvars_scaled <- traitvars  
# for(i in 1:ncol(traitvars)){
#   traitvars_scaled[,i] <- traitvars[,i]/sum(traitvars[,i])
# }
# 
# traitvars_scaled2 <- traitvars_scaled[c(2,3,1,4),]
# quartz(width=5, height=4)
# cols <- brewer.pal(11, "RdBu")[c(11,9,1, 6)]
# barplot(as.matrix(traitvars_scaled2),beside=F,legend.text = F,xpd = T, names.arg = c("logLMA", "Leaf\nLife","logN"),args.legend = list(x=4, y=1.3, ncol=2), col = paste0(cols,"CC"), ylab="Proportion of total Variance", xlab="Traits.common")
# legend(xpd=T, x = 1, y=1.3, legend=c("W/in Plot", "Btw Plots", "Btw Spp", "Btw Genera"), fill=paste0(cols,"CC")[c(4,3,1,2)], ncol=2, bty="n",  cex=1.2)


############## **Making Master Dataset for Variance Decomposition** ###################
###### trying with combined PACNW and glopnet dataset. ###
data1 <- traits %>% select(FullSpecies,log.LMA, log.LL, log.Nmass, log.Narea, GENUS, Family)
colnames(data1)[c(1,2,6)]<- c("Species", "log.LMA","Genus")
data1$Project <- rep("PACNW", times=nrow(data1))
data1$log.LL[which(data1$log.LL<1.2)] <- NA # all the deciduous species have '1yr' lifespan, but really that's wrong
data2 <- LES %>% filter(!is.na(Family)) %>% select(Species, log.LMA, log.LL, log.Nmass, log.Narea,Genus, Family)
data2$Project <- rep("GLOPNET", times=nrow(data2))
data.supp <- read.csv("/Users/leeanderegg/Dropbox/NACP_Traits/Intra-data/OtherData_Combined_033017.csv", header=T, row.names=1)
data3 <- data.supp %>% select(Species,log.LMA,log.LL, log.Nmass,log.Narea,Genus,Family)

## make combined dataset with all the taxonomically resolved species
data.all <- rbind(data1, data2, data3)
data.all$Species <- factor(data.all$Species)
data.all$Genus <- factor(data.all$Genus)
data.all$Family <- factor(data.all$Family)
data.all$Project <- factor(data.all$Project)

## with all spp used for the hierarchical analysis
#logLMAvar <- lmer(log.LMA~ Project + (1|Family/Genus/Species), data.all)
logLMAvar <- lmer(log.LMA~ Project + (1|Family) + (1|Genus) + (1|Species), data.all)
logLLvar <- lmer(log.LL~ Project + (1|Family) + (1|Genus) + (1|Species), data.all)
logNmassvar <- lmer(log.Nmass~ Project + (1|Family) + (1|Genus) + (1|Species), data.all)
logNareavar <- lmer(log.Narea~ Project + (1|Family) + (1|Genus) + (1|Species), data.all)

LMAvar <- lmer(10^log.LMA~ Project + (1|Family) + (1|Genus) + (1|Species), data.all)
LLvar <- lmer(10^log.LL~ Project + (1|Family) + (1|Genus) + (1|Species), data.all)
Nmassvar <- lmer(10^log.Nmass~ Project + (1|Family) + (1|Genus) + (1|Species), data.all)
Nareavar <- lmer(10^log.Narea~ Project + (1|Family) + (1|Genus) + (1|Species), data.all)
#logNmassvar <- lmer(log.Nmass~ Project + (1|Family) + (1|Genus) + (1|Species), data.all)
    # in absolute scale (months), the variation appears to be HUGE w/ind spp, then w/in genera, then between Families
    # all levels of this analysis are hugely non-normal, except Families aren't too bad
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


rtraitvars <- data.frame(rLMAvariance[,4], rLLvariance[,4], rNmassvariance[,4], rNareavariance[,4])
colnames(rtraitvars) <- c("LMA", "LL", "Nmass", "Narea")
rownames(rtraitvars) <- c("BtwSpecies", "BtwGenera", "BtwFamilies", "WtinSpecies")

rtraitvars_scaled <- rtraitvars  
for(i in 1:ncol(rtraitvars)){
  rtraitvars_scaled[,i] <- rtraitvars[,i]/sum(rtraitvars[,i])
}
rtraitvars_scaled2 <- rtraitvars_scaled[c(3,2,1,4),]
traitvars_scaled2 <- traitvars_scaled[c(3,2,1,4),]

quartz(width=5, height=4)
cols <- brewer.pal(11, "RdBu")[c(11,9,1, 6)]
barplot(as.matrix(traitvars_scaled2),beside=F,legend.text = F,xpd = T, names.arg = c("logLMA", "logLL","logNmass","logNarea"),args.legend = list(x=4, y=1.3, ncol=2), col = paste0(cols,"CC"), ylab="Proportion of total Variance\n(log traits)", xlab="", mgp=c(2,1,0))
legend(xpd=T, x = 0, y=1.3, legend=rownames(traitvars_scaled2), fill=paste0(cols,"CC"), ncol=2, bty="n",  cex=1.2)


quartz(width=5, height=4)
cols <- brewer.pal(11, "RdBu")[c(11,9,1, 6)]
barplot(as.matrix(rtraitvars_scaled2),beside=F,legend.text = F,xpd = T, names.arg = c("LMA", "Leaf\nLife","Nmass","Narea"),args.legend = list(x=4, y=1.3, ncol=2), col = paste0(cols,"CC"), ylab="Proportion of total Variance\n(raw traits)", mgp=c(2,1,0))
mtext(text="28 spp w/ replication, 1000+ spp,\n500+ genera, 150+ families",side = 1,line = 3.3)
legend(xpd=T, x = 0, y=1.3, legend=rownames(rtraitvars_scaled2), fill=paste0(cols,"CC"), ncol=2, bty="n",  cex=1.2)


#### Check that these normal distributions are a decent way to do this #######
qqp(ranef(logLMAvar)$Family[,1]) # 179 families
qqp(ranef(logLMAvar)$Genus[,1]) # 752 genera
qqp(ranef(logLMAvar)$Species[,1]) # 1538 species
qqp(resid(logLMAvar)) # 3248 trait measurements
  # has some small and large outliers, particularly for Genus and Species, but looks pretty symmetrical
  # resids have some serious high and low outliers. They're symmetrical, but...

qqp(ranef(logLLvar)$Family[,1]) # 109 families
qqp(ranef(logLLvar)$Genus[,1]) # 331 genera
qqp(ranef(logLLvar)$Species[,1]) # 548 species
qqp(resid(logLLvar)) #1596 trait measurements
  # Family and Genus REAL normal. Species has some outliers both high and low.


qqp(ranef(logNmassvar)$Family[,1]) # 175 families
qqp(ranef(logNmassvar)$Genus[,1]) # 695 genera
qqp(ranef(logNmassvar)$Species[,1]) # 1342 species
qqp(resid(logNmassvar)) #2890 trait measurments
  #Family and Species pretty normal. Genus has some assymetrically low outliers...


######### plotting Random Effects estimates of Family and Genus values #####
logLLfam <- ranef(logLLvar)$Family
logNmassfam <- ranef(logNmassvar)$Family
logLMAfam <- ranef(logLMAvar)$Family

logLMAfam$log.Nmass <- logNmassfam[match(row.names(logLMAfam),row.names(logNmassfam)),1] # 8 families w/Nmass don't have LMA, and 12 with LMA don't have Nmass
logLMAfam$log.LL <- logLLfam[match(row.names(logLMAfam),row.names(logLLfam)),1] # three fams w/ LL don't have LMA & 73 families with LMA don't have LL
colnames(logLMAfam)[1] <-"log.LMA"
plot(log.LMA~log.LL, logLMAfam)
plot(log.LMA~log.LL, LESfam)

      # so the emergent covariation between N, LL, and LMA at the family level looks pretty similar w/ family averages vs modeled family values




#-----------------------------------------------------------------------------

################# LEAF ECONOMIC SPECTRUM ANALYSIS ###################

#-----------------------------------------------------------------------------



xyplot(SLA_HSA~LEAF_LIFE|GE.SP,traits.common[-which(traits.common$SP.ID %in% c("CALDEC" ,"JUNOCC", "ARBMEN")),], type=c('p',"r"), col.line= "red3") # leaf lifespan and SLA may be related in true firs and doug fir

colslight <- paste0(c(mapcols,mypal), "66")
palette(colslight)
plot(SLA_HSA~LEAF_LIFE, traits.common, pch=16, col=SP.ID)
for(i in levels(traits.common$SP.ID)[-c(6,7,8)]){
  tmp <- subset(traits.common, subset=SP.ID==i)
  mod <- lm(SLA_HSA~LEAF_LIFE, tmp)
  lines(fitted(mod)~tmp$LEAF_LIFE[which(tmp$LEAF_LIFE>0 & tmp$SLA_HSA>0)], col=c(mapcols,mypal)[which(levels(traits.common$SP.ID)==i)])
}
LESmod1 <- lmer(SLA_HSA~LEAF_LIFE + (LEAF_LIFE|SP.ID), traits.common[which(traits.common$LEAF_LIFE>0),])
LESmod1null <- lmer(SLA_HSA~1 + (LEAF_LIFE|SP.ID), traits.common[which(traits.common$LEAF_LIFE>0),])
LESmod1null1 <- lmer(SLA_HSA~1 + (1|SP.ID), traits.common[which(traits.common$LEAF_LIFE>0),])
anova(LESmod1null, LESmod1) # marginally significant
anova(LESmod1null, LESmod1null1) # significnat
sjp.lmer(LESmod1)
sjp.lmer(LESmod1, type="rs.ri")

plot(LEAF_NITROGEN~LEAF_LIFE, traits.common, pch=16, col=SP.ID)
for(i in levels(traits.common$SP.ID)[-c(6,7,8)]){
  tmp <- subset(traits.common, subset=SP.ID==i)
  mod <- lm(LEAF_NITROGEN~LEAF_LIFE, tmp)
  lines(fitted(mod)~tmp$LEAF_LIFE[which(tmp$LEAF_LIFE>0 & tmp$LEAF_NITROGEN>0)], col=c(mapcols,mypal)[which(levels(traits.common$SP.ID)==i)])
}
  ## yeah, on average it sort of does. but really ns
LESmod2 <- lmer(LEAF_NITROGEN~LEAF_LIFE + (LEAF_LIFE|SP.ID), traits.common[which(traits.common$LEAF_LIFE>0),])
LESmod2null <- lmer(LEAF_NITROGEN~1 + (LEAF_LIFE|SP.ID), traits.common[which(traits.common$LEAF_LIFE>0),])
LESmod2null1 <- lmer(LEAF_NITROGEN~1 + (1|SP.ID), traits.common[which(traits.common$LEAF_LIFE>0),])
anova(LESmod2null, LESmod2) 
anova(LESmod2null, LESmod2null1)
    # so the random slopes are super significnat (lots of within/species variation)
    # but it's not significantly in the direction of the LES
sjp.lmer(LESmod2)
sjp.lmer(LESmod2, type="rs.ri")



plot(log(LEAF_NITROGEN)~log(SLA_HSA), traits.common, pch=16, col=SP.ID)
for(i in levels(traits.common$SP.ID)){
  tmp <- subset(traits.common, subset=SP.ID==i)
  mod <- lm(log(LEAF_NITROGEN)~log(SLA_HSA), tmp)
  lines(fitted(mod)~log(tmp$SLA_HSA[which(tmp$SLA_HSA>0 & tmp$LEAF_NITROGEN>0)]), col=c(mapcols,mypal)[which(levels(traits.common$SP.ID)==i)])
}
## yeah, on average it sort of does. but really ns
LESmod3 <- lmer(log(LEAF_NITROGEN)~log(SLA_HSA) + (log(SLA_HSA)|SP.ID), traits.common[which(traits.common$SLA_HSA>0),])
LESmod3null <- lmer(log(LEAF_NITROGEN)~1 + (log(SLA_HSA)|SP.ID), traits.common[which(traits.common$SLA_HSA>0),])
LESmod3null1 <- lmer(LEAF_NITROGEN~1 + (1|SP.ID), traits.common[which(traits.common$SLA_HSA>0),])
anova(LESmod3null, LESmod3) 
anova(LESmod3null, LESmod3null1)
# The relationship is actually significnat with either Log Log or not values. p=0.001 unlog, 0.01769 logged
sjp.lmer(LESmod3)
sjp.lmer(LESmod3, type="rs.ri")



######### .Comparison to Wright et al. 2004 Nature ###########

# overlapping species in both datasets
LESoverlap <- data.frame(LES[which(LES$GE.SP %in% speciesID$traits.sp),] %>% arrange(GE.SP)) [-c(3,4), ] # get rid of the second 2 Arbutus menziesii
traitsoverlap <- traits[which(traits$GE.SP %in% LES$GE.SP),]
traitsoverlap.means <- data.frame(traitsoverlap %>% group_by(GE.SP, SP.ID) %>% summarise(mlog.LL = mean(log.LL, na.rm=T), mlog.LMA=mean(log.LMA, na.rm=T), mlog.Nmass = mean(log.Nmass, na.rm=T), mlog.Narea = mean(log.Narea, na.rm=T), mlog.LMA_PSA = mean(log.LMA_PSA, na.rm=T)
                                                                                  , selog.LL = sterr(log.LL, na.rm=T), selog.LMA=sterr(log.LMA, na.rm=T), selog.Nmass = sterr(log.Nmass, na.rm=T), selog.Narea = sterr(log.Narea, na.rm=T), selog.LMA_PSA = sterr(log.LMA_PSA, na.rm=T)
                                                                                  , sdlog.LL = sd(log.LL, na.rm=T), sdlog.LMA=sd(log.LMA, na.rm=T), sdlog.Nmass = sd(log.Nmass, na.rm=T), sdlog.Narea = sd(log.Narea, na.rm=T), sdlog.LMA_PSA = sd(log.LMA_PSA, na.rm=T)))
##--- comparing trait means for species present in both datasets ---###
      # takehome, not that much trait overlap
      # somethings wring with my Narea
      # traits are medium similar. though not bang on.
quartz(width=5, height=5)
par(mfrow=c(2,2))
for (i in c("log.LL", "log.LMA", "log.Nmass", "log.LMA_PSA")){
#  plot(LESoverlap[,i]~traitsoverlap.means[,paste0("m",i)], pch=16, col=factor(LESoverlap$GE.SP), main=i, ylab="LES", xlab="PACNW traits")
  plot(traitsoverlap.means[,paste0("m",i)]~LESoverlap[,i], pch=16, col=factor(LESoverlap$GE.SP), main=i, xlab="LES", ylab="PACNW traits")
  arrows(x0=LESoverlap[,i], y0=traitsoverlap.means[,paste0("m",i)], y1=traitsoverlap.means[,paste0("m",i)]+traitsoverlap.means[,paste0("sd",i)], length=.01, angle=90)
  arrows(x0=LESoverlap[,i], y0=traitsoverlap.means[,paste0("m",i)], y1=traitsoverlap.means[,paste0("m",i)]-traitsoverlap.means[,paste0("sd",i)], length=.01, angle=90)
  abline(a=0, b=1)
  if(i == "log.LL"){legend('topleft', legend=traitsoverlap.means$SP.ID, pch=16, col=mypal[1:nrow(LESoverlap)], cex=.6, bty="n")}
}

## comparing LMA (with HSA) and recalculated LMA with PSA
i <- "log.LMA"
plot(traitsoverlap.means[,paste0("m",i)]~LESoverlap[,i], pch=16, col=factor(LESoverlap$GE.SP), main=i, xlab="LES", ylab="PACNW traits", ylim=c(1.2,2.9))
arrows(x0=LESoverlap[,i], y0=traitsoverlap.means[,paste0("m",i)], y1=traitsoverlap.means[,paste0("m",i)]+traitsoverlap.means[,paste0("sd",i)], length=.01, angle=90)
arrows(x0=LESoverlap[,i], y0=traitsoverlap.means[,paste0("m",i)], y1=traitsoverlap.means[,paste0("m",i)]-traitsoverlap.means[,paste0("sd",i)], length=.01, angle=90)
abline(a=0,b=1)
points(traitsoverlap.means$mlog.LMA_PSA~LESoverlap$log.LMA, col=factor(LESoverlap$GE.SP))

## using PSA based LMA looks a little more like Wright, but actually makes PSME worse.
plot(I(10^traitsoverlap.means$mlog.LMA) ~ I(10^LESoverlap$log.LMA), pch=16, ylim=c(45, 300), ylab="PACNW LMA", xlab="Wright et al LMA")
abline(a=0,b=1)
points(I(10^traitsoverlap.means$mlog.LMA_PSA) ~ I(10^LESoverlap$log.LMA))
legend('bottomright', legend=c("HSA","PSA"), pch=c(16, 1), bty='n')



##----- other visuals of how they compare -----###

# N vs LL
ggplot(traits, aes(x=log.Nmass, y=log.LL, colour=SP.ID)) + geom_point(aes(alpha=1)) + geom_point(data=LESoverlap, aes(colour=NULL))
  # so the decid species all stack up at 1 yr, but should be lower

# N vs LMA of overlapping spp
ggplot(traitsoverlap, aes(x=log.Nmass, y=log.LMA, colour=GE.SP)) + geom_point(aes(alpha=1)) + geom_point(data=LESoverlap, aes(colour=GE.SP, size=3))

# LMA vs LL of overlapping spp
ggplot(traitsoverlap, aes(x=log.LL, y=log.LMA, colour=GE.SP)) + geom_point(aes(alpha=1)) + geom_point(data=LESoverlap, aes(colour=GE.SP, size=3))


####### ..Trait Space Overlap ###############
ggplot(LES, aes(x=log.Nmass, y=log.LL, colour=Needle.Broad.lf, alpha =1)) + geom_point() + geom_smooth() + geom_point(data=traits, aes(x=log.Nmass, y=log.LL, colour=NA)) + geom_smooth(data=traits, aes(x=log.Nmass, y=log.LL, colour=NA))
    # take home: PACNW traits are all way on the long lived side, but cover a medium portion of the N spectrum
ggplot(LES, aes(x=log.Nmass, y=log.LMA, colour=Needle.Broad.lf, alpha =1)) + geom_point() + geom_smooth() + geom_point(data=traits, aes(x=log.Nmass, y=log.LMA_PSA, colour=NA)) + geom_smooth(data=traits, aes(x=log.Nmass, y=log.LMA_PSA, colour=NA))

ggplot(LES, aes(x=log.Nmass, y=10^log.LMA, colour=Needle.Broad.lf, alpha =1)) + geom_point() + geom_smooth() + geom_point(data=traits, aes(x=log.Nmass, y=LMA_PSA, colour=NA)) + geom_smooth(data=traits, aes(x=log.Nmass, y=log.LMA, colour=NA))


# unlogged LL and Nmass
ggplot(LES, aes(x=I(10^log.Nmass), y=10^log.LL, colour=Needle.Broad.lf, alpha =1)) + geom_point() + geom_smooth() + geom_point(data=traits, aes(x=10^log.Nmass, y=10^log.LL, colour=NA)) + geom_smooth(data=traits, aes(x=10^log.Nmass, y=10^log.LL, colour=NA))
# unlogged LMA and Nmass
ggplot(LES, aes(x=I(10^log.Nmass), y=10^log.LMA, colour=Needle.Broad.lf, alpha =1)) + geom_point() + geom_smooth() + geom_point(data=traits, aes(x=10^log.Nmass, y=10^log.LMA, colour=NA)) + geom_smooth(data=traits, aes(x=10^log.Nmass, y=10^log.LMA, colour=NA))
# unlogged LMA_PSA and Nmass


p <- ggplot(LES, aes(x=log.Nmass, y=log.LL, colour=Needle.Broad.lf))
p + geom_point() + scale_color_grey()  +  theme_classic() +theme(
  axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))  # theoretically throws down a nice looking classic plot. Mine is missing axis lines
  # also (start=0.8, end=0.2) reverses the greyscale rampe
p + geom_point() + scale_color_hue(l=65, c=100)  # bigger l= lighter, bigger c = brighter 
p + geom_point(aes(colour=log.LMA)) + scale_color_gradient(low=mypal[1], high=mypal[2])#low="blue", high="darkred")



###### Nmass ~ LL : PACNW vs LES trends ##
quartz(width=5, height=4)
ggplot(LES, aes(x=log.LL, y=log.Nmass)) + geom_point(color="grey") + geom_point(data=traits.common, aes(colour=SP.ID), alpha=1/5) +
  geom_smooth(method="lm", se=F, colour='darkgrey') + geom_smooth(data=traits.common, aes(colour=SP.ID),method="lm", se=F)

###### LMA ~ LL : PACNW vs LES trends ##
ggplot(LES, aes(x=log.LMA, y=log.LL)) + geom_point(color="grey") + geom_point(data=traits.common, aes(x=log.LMA_PSA, colour=SP.ID), alpha=1/5) +
  geom_smooth(method="lm", se=F, colour='darkgrey') + geom_smooth(data=traits.common, aes(x=log.LMA_PSA, colour=SP.ID),method="lm", se=F)

ggplot(LES, aes(x=log.LMA, y=log.LL)) + geom_point(color="grey") + geom_point(data=traits.common, aes(colour=SP.ID), alpha=1/5) +
  geom_smooth(method="lm", se=F, colour='darkgrey') + geom_smooth(data=traits.common, aes(colour=SP.ID),method="lm", se=F)

###### LMA ~ N : PACNW vs LES trends ##
ggplot(LES, aes(x=log.LMA, y=log.Nmass)) + geom_point(color="grey") + geom_point(data=traits.common, aes(x=log.LMA_PSA, colour=SP.ID), alpha=1/5) +
  geom_smooth(method="lm", se=F, colour='darkgrey') + geom_smooth(data=traits.common, aes(x=log.LMA_PSA, colour=SP.ID),method="lm", se=F)

ggplot(LES, aes(x=log.LMA, y=log.Nmass)) + geom_point(color="grey") + geom_point(data=traits.common, aes( colour=SP.ID), alpha=1/5) +
  geom_smooth(method="lm", se=F, colour='darkgrey') + geom_smooth(data=traits.common, aes( colour=SP.ID),method="lm", se=F)




###### statistics ##########





#-----------------------------------------------------------------------------
################# SLA ANALYSIS ###################
#-----------------------------------------------------------------------------

## .Douglas Fir : a P~T interaction ###########
## Colored continuously by tmean
colors <- rev(heat.colors(101)) 
colors <- rev(brewer.pal(n=11,"RdYlBu"))
colramp <- colorRampPalette(colors[c(1,6,11)])
z <- traits.common$tmean.gy.c[which(traits.common$GE.SP=="Pseudotsuga.menziesii")]
zcolor <- colramp(101)[(z - min(z))/diff(range(z))*100 + 1] 
zcolorlight <- paste0(zcolor, "BB")
plot(SLA_HSA~ppt.gy.mm, traits.common, subset=GE.SP=="Pseudotsuga.menziesii", col=zcolorlight, pch=16, main="Doug fir SLA ~f(clim)")
abline(lm(SLA_HSA~ppt.gy.mm, traits.common, subset=GE.SP=="Pseudotsuga.menziesii" & tmean.gy.c > 10), col=colors[11], lwd=2)
abline(lm(SLA_HSA~ppt.gy.mm, traits.common, subset=GE.SP=="Pseudotsuga.menziesii" & tmean.gy.c <= 10 & tmean.gy.c > 8), col="grey", lwd=2)
abline(lm(SLA_HSA~ppt.gy.mm, traits.common, subset=GE.SP=="Pseudotsuga.menziesii" & tmean.gy.c <= 8), col=colors[1], lwd=2)
legend(xpd=T, "topleft", pch=16, col=colors[c(11,1)], legend=c("Tmean > 10C", "Tmean < 8C"), bty="n")

plot(SLA_HSA~ppt.gy.mm, traits.common, subset=GE.SP=="Pseudotsuga.menziesii" & tmean.gy.c > 10, col=colors[11], pch=16)
points(SLA_HSA~ppt.gy.mm, traits.common, subset=GE.SP=="Pseudotsuga.menziesii" & tmean.gy.c <= 10 & tmean.gy.c > 8, col=colors[6], pch=16)
points(SLA_HSA~ppt.gy.mm, traits.common, subset=GE.SP=="Pseudotsuga.menziesii" & tmean.gy.c <= 8, col=colors[1], pch=16)


### or we'll try the added variable plots/partial plots (using the {car} package)
testSLA <- lm(SLA_HSA ~ ppt.gy.mm * tmean.gy.c, traits.common, subset=GE.SP=="Pseudotsuga.menziesii")
avPlots(testSLA)
avPlot(testSLA, variable='ppt.gy.mm')
  # NOTE: because there's and INTERACTION, I'm not entirely sure how to plot this...

tmpdata <- traits.common[which(traits.common$SP.ID=="PSEMEN"),]
tmpdata$ppt.gy.mmsc <- scale(tmpdata$ppt.gy.mm)
tmpdata$tmean.gy.csc <- scale(tmpdata$tmean.gy.c)
SLApsme <- lmer(scale(SLA_HSA)~ ppt.gy.mmsc * tmean.gy.csc + (1|PLOT_ID), tmpdata, REML=F)
SLApsmenull <- lmer(scale(SLA_HSA)~ 1 + (1|PLOT_ID), traits.common, subset=SP.ID=="PSEMEN")
reduced <- step(SLApsme)

## trying out the same model but with climate PCs instead
SLApsme1 <- lmer(scale(SLA_HSA)~ climPC1 * climPC2 + climPC3 + (1|PLOT_ID), tmpdata, REML=F)
SLApsme1 <- lmer(scale(SLA_HSA)~ climPC1 * climPC2 + climPC3 + (1|PLOT_ID), tmpdata, REML=F)


### for shits, let's just look at the whole common traits dataset
tmp <- lm(SLA_HSA ~ scale(ppt.gy.mm) * scale(tmean.gy.c), traits.common[-which(traits.common$SP.ID == "PSEMEN"),])
z <- traits.common$tmean.gy.c
zcolor <- colramp(101)[(z - min(z))/diff(range(z))*100 + 1] 
zcolorlight <- paste0(zcolor, "BB")
plot(SLA_HSA~ppt.gy.mm, traits.common, subset= SP.ID != "PSEMEN", col=zcolorlight, pch=16)
abline(lm(SLA_HSA~ppt.gy.mm, traits.common, subset= tmean.gy.c > 9), col=colors[11], lwd=2)
abline(lm(SLA_HSA~ppt.gy.mm, traits.common, subset= tmean.gy.c <= 9), col=colors[1], lwd=2)
  # this appears to hold in the whole common traits dataset, with or without PSME. However, doesn't hold in the entire dataset quite as strongly. But that's probably due to HUGE outliers.



####### .visualizing SLA trends in other Species #####
# Ponderosa Pine
z <- traits.common$tmean.gy.c[which(traits.common$GE.SP=="Pinus.ponderosa")]
zcolor <- colramp(101)[(z - min(z))/diff(range(z))*100 + 1] 
zcolorlight <- paste0(zcolor, "BB")
plot(SLA_HSA~ppt.gy.mm, traits.common, subset=GE.SP=="Pinus.ponderosa", col=zcolorlight, pch=16)
abline(lm(SLA_HSA~ppt.gy.mm, traits.common, subset=GE.SP=="Pinus.ponderosa" & tmean.gy.c > 10), col=colors[11], lwd=2)
abline(lm(SLA_HSA~ppt.gy.mm, traits.common, subset=GE.SP=="Pinus.ponderosa" & tmean.gy.c <= 10), col=colors[1], lwd=2)

testSLA <- lm(SLA_HSA ~ ppt.gy.mm * tmean.gy.c, traits.common, subset=SP.ID=="PINPON")
## Shit. generally no relationships with ANYTHING. No climate variables whatsoever!


# Other species:
sp <- "ABICON"

z <- traits.common$tmean.gy.c[which(traits.common$SP.ID==sp)]
zcolor <- colramp(101)[(z - min(z))/diff(range(z))*100 + 1] 
zcolorlight <- paste0(zcolor, "FF")
plot(SLA_HSA~ppt.gy.mm, traits.common, subset=SP.ID==sp, col=zcolorlight, pch=16)
abline(lm(SLA_HSA~ppt.gy.mm, traits.common, subset=SP.ID==sp & tmean.gy.c > 10), col=colors[11], lwd=2)
abline(lm(SLA_HSA~ppt.gy.mm, traits.common, subset=SP.ID==sp & tmean.gy.c <= 10), col=colors[1], lwd=2)

testSLA <- lm(SLA_HSA ~ ppt.gy.mm * tmean.gy.c, traits.common, subset=SP.ID==sp)
## Hmmm. Generally some evidence for a relationship with Temp (less so but possibly still there with PPT)

#### .3d SLA plotting ######
open3d()
with(traits.common[which(traits.common$GE.SP=="Pseudotsuga.menziesii"),], plot3d(tmean.gy.c, ppt.gy.mm, SLA_HSA, box=F, col=as.numeric(GE.SP)))
with(traits.common[which(traits.common$GE.SP=="Pseudotsuga.menziesii"),], plot3d(climPC1, climPC2, LEAF_LIFE, box=F, col=as.numeric(GE.SP)))

my.lm1 <- lm(SLA_HSA~tmean.gy.c*ppt.gy.mm, traits.common, subset= GE.SP=="Pseudotsuga.menziesii" )
my.grd <- expand.grid(tmean.gy.c=seq(4,12,by=.5), ppt.gy.mm=seq(450,3000, by=50))
my.grd$pred <- predict(my.lm1, newdata = my.grd)
predictions <- matrix(my.grd[[3]], length(unique(my.grd[[1]])), length(unique(my.grd[[2]])))
persp3d(x=unique(my.grd[[1]]), y=unique(my.grd[[2]])
        , z=predictions, add=T,color="#FF2222",alpha=0.5)



### All species 

open3d()
with(traits.common, plot3d(tmean.gy.c, ppt.gy.mm, SLA_HSA, box=F, col=as.numeric(GE.SP)))

my.lm1 <- lm(SLA_HSA~tmean.gy.c*ppt.gy.mm, traits.common, subset= GE.SP=="Pseudotsuga.menziesii" )
my.grd <- expand.grid(tmean.gy.c=seq(4,12,by=.5), ppt.gy.mm=seq(450,3000, by=50))
my.grd$pred <- predict(my.lm1, newdata = my.grd)
predictions <- matrix(my.grd[[3]], length(unique(my.grd[[1]])), length(unique(my.grd[[2]])))
persp3d(x=unique(my.grd[[1]]), y=unique(my.grd[[2]])
        , z=predictions, add=T,color="#FF2222",alpha=0.5)



### NOTE in the gamm land, I don't know what I'm doing with gam, much less gamm (from {mgcv}) or gamm4 {gamm4}
# THUS: I need to spend a lot more time figuring out exactly what this is doing.
SLAgam <- gamm4(I(log(SLA_HSA))~s(MAT) + DIVISION, random= ~ (1|PLOT_ID) + (1|GENUS) + (1|COMMON_NAME), data=traits.common, subset =MAT>-10 & SLA_HSA>0)
SLAgam1 <- gamm4(I(log(SLA_HSA))~s(MAT), random= ~ (1|PLOT_ID) + (1|GENUS) + (1|COMMON_NAME), data=traits.common, subset =MAT>-10 & SLA_HSA>0)
SLAgam2 <- gamm4(I(log(SLA_HSA))~1, random= ~ (1|PLOT_ID) + (1|GENUS) + (1|COMMON_NAME), data=traits.common, subset =MAT>-10 & SLA_HSA>0)
anova(SLAgam2$mer, SLAgam1$mer, SLAgam$mer) # well, that would suggest that MAT isn't important (if same smoother applied to all species)






#-----------------------------------------------------------------------------
################# LEAF_LIFE ANALYSIS ###################
#-----------------------------------------------------------------------------


## .Douglas Fir : a P~T interaction ###########
## Colored continuously by tmean
colors <- rev(heat.colors(101)) 
colors <- rev(brewer.pal(n=11,"RdYlBu"))
colramp <- colorRampPalette(colors[c(1,6,11)])
z <- traits.common$tmean.gy.c[which(traits.common$GE.SP=="Pseudotsuga.menziesii")]
zcolor <- colramp(101)[(z - min(z))/diff(range(z))*100 + 1] 
zcolorlight <- paste0(zcolor, "BB")
plot(LEAF_LIFE~ppt.gy.mm, traits.common, subset=GE.SP=="Pseudotsuga.menziesii", col=zcolorlight, pch=16, main="Doug fir SLA ~f(clim)")
abline(lm(LEAF_LIFE~ppt.gy.mm, traits.common, subset=GE.SP=="Pseudotsuga.menziesii" & tmean.gy.c > 10), col=colors[11], lwd=2)
abline(lm(LEAF_LIFE~ppt.gy.mm, traits.common, subset=GE.SP=="Pseudotsuga.menziesii" & tmean.gy.c <= 10 & tmean.gy.c > 8), col="grey", lwd=2)
abline(lm(LEAF_LIFE~ppt.gy.mm, traits.common, subset=GE.SP=="Pseudotsuga.menziesii" & tmean.gy.c <= 8), col=colors[1], lwd=2)
legend(xpd=T, "topright", pch=16, col=colors[c(11,1)], legend=c("Tmean > 10C", "Tmean < 8C"), bty="n")

plot(SLA_HSA~ppt.gy.mm, traits.common, subset=GE.SP=="Pseudotsuga.menziesii" & tmean.gy.c > 10, col=colors[11], pch=16)
points(SLA_HSA~ppt.gy.mm, traits.common, subset=GE.SP=="Pseudotsuga.menziesii" & tmean.gy.c <= 10 & tmean.gy.c > 8, col=colors[6], pch=16)
points(SLA_HSA~ppt.gy.mm, traits.common, subset=GE.SP=="Pseudotsuga.menziesii" & tmean.gy.c <= 8, col=colors[1], pch=16)

LL_test0 <- lmer(LEAF_LIFE~ 1 + (1|PLOT_ID), traits.common, subset=GE.SP=="Pseudotsuga.menziesii",REML=F)
LL_test1 <- lmer(LEAF_LIFE~ ppt.gy.mm + (1|PLOT_ID), traits.common, subset=GE.SP=="Pseudotsuga.menziesii",REML=F)
LL_test2 <- lmer(LEAF_LIFE~ tmean.gy.c + (1|PLOT_ID), traits.common, subset=GE.SP=="Pseudotsuga.menziesii",REML=F)
LL_test3 <- lmer(LEAF_LIFE~ ppt.gy.mm + tmean.gy.c + (1|PLOT_ID), traits.common, subset=GE.SP=="Pseudotsuga.menziesii",REML=F)
LL_test4 <- lmer(LEAF_LIFE~ ppt.gy.mm * tmean.gy.c + (1|PLOT_ID), traits.common, subset=GE.SP=="Pseudotsuga.menziesii",REML=F)
LL_testnew <- lmer(LEAF_LIFE~ climPC1 * climPC2  +  climPC3 + (1|PLOT_ID), traits.common, subset=GE.SP=="Pseudotsuga.menziesii",REML=F)
step(LL_testnew) #looks like climPC1 * climPC2 + climPC3 is probs the best model
AIC(LL_test0, LL_test1, LL_test2, LL_test3, LL_test4, LL_testnew)
anova(LL_test3, LL_test4)
  # LL decreases w/ PPT and w/ T

#### .3d LEAF_LIFE plotting ######

## PSEMEN
with(traits.common[which(traits.common$GE.SP=="Pseudotsuga.menziesii"),], plot3d(tmean.gy.c, ppt.gy.mm, LEAF_LIFE, box=F, col=as.numeric(GE.SP)))

my.lm2 <- lm(LEAF_LIFE~tmean.gy.c+ppt.gy.mm, traits.common, subset= GE.SP=="Pseudotsuga.menziesii" )
my.grd2 <- expand.grid(tmean.gy.c=seq(4,12,by=.5), ppt.gy.mm=seq(450,3000, by=50))
my.grd2$pred1 <- predict(my.lm2, newdata = my.grd)
predictions <- matrix(my.grd2[[3]], length(unique(my.grd2[[1]])), length(unique(my.grd2[[2]])))
persp3d(x=unique(my.grd2[[1]]), y=unique(my.grd2[[2]])
        , z=predictions, add=T,color="#FF2222",alpha=0.5, front="lines")

## all species
with(traits.common, plot3d(tmean.gy.c, ppt.gy.mm, LEAF_LIFE, box=F, col=as.numeric(GE.SP)))

my.lm2 <- lm(LEAF_LIFE~tmean.gy.c*ppt.gy.mm, traits.common, subset= GE.SP=="Pseudotsuga.menziesii" )
my.grd2 <- expand.grid(tmean.gy.c=seq(4,12,by=.5), ppt.gy.mm=seq(450,3000, by=50))
my.grd2$pred1 <- predict(my.lm2, newdata = my.grd)
predictions <- matrix(my.grd2[[3]], length(unique(my.grd2[[1]])), length(unique(my.grd2[[2]])))
persp3d(x=unique(my.grd2[[1]]), y=unique(my.grd2[[2]])
        , z=predictions, add=T,color="#FF2222",alpha=0.5, front="lines")




###### ***Full LEAF_LIFE modeling***: #######

## Illustrate the interactions for Leaf Life
avPlots(model=lm(LEAF_LIFE~climPC1*climPC2 + climPC3, data=traits.common[which(traits.common$SP.ID=="PSEMEN"),]))
avPlots(model=lm(SLA_HSA~climPC1*climPC2 + BIOST_TGROWTHgam, data=traits.common[which(traits.common$SP.ID=="PSEMEN"),]))

avPlots(model=lm(LEAF_NITROGEN~climPC1*climPC2 + climPC3 + BIOST_TGROWTH, data=traits.common[which(traits.common$SP.ID=="PSEMEN"),]))

## model Leaf Life linearly
LLmodclim <- lmer(LEAF_LIFE~climPC1*climPC2 + climPC3 + (climPC1 + climPC2 + climPC3|SP.ID) + (1|PLOT_ID), traits.common[which(traits.common$BIOST_TGROWTH > -500),], REML=F)
LLmodG <- lmer(LEAF_LIFE~scale(AG_TGROWTH) + (scale(AG_TGROWTH)|SP.ID), traits.common[which(traits.common$BIOST_TGROWTH > -500),], REML=F)
LLmodAll <- lmer(LEAF_LIFE~climPC1*climPC2 + climPC3 + scale(AG_TGROWTH) + (scale(AG_TGROWTH)|SP.ID)+ (climPC1 + climPC2 + climPC3|SP.ID) + (1|PLOT_ID), traits.common[which(traits.common$BIOST_TGROWTH > -500),], REML=F)
AIC(LLmodclim, LLmodG, LLmodAll)

## model Leaf Life as a poisson
LLpois <- glm(LEAF_LIFE~climPC1 * climPC2 + climPC3 + GENUS , traits.common, family = poisson)
LLpois <- glmer(LEAF_LIFE~climPC1 * climPC2 + climPC3 + (1|GENUS/SP.ID/PLOT_ID ), traits.common, family = poisson)
LLpois <- glmer(LEAF_LIFE~scale(AG_TGROWTH) + (scale(AG_TGROWTH)|SP.ID), traits.common, family=poisson)

LLpoisG <- glmer(LEAF_LIFE~scale(BIOST_TGROWTH) + (scale(BIOST_TGROWTH)|SP.ID), traits.common[which(traits.common$BIOST_TGROWTH > -500),], family=poisson)
LLpoisC <- glmer(LEAF_LIFE~climPC1 * climPC2 + climPC3 + (climPC1|SP.ID) + (climPC2|SP.ID) , traits.common[which(traits.common$BIOST_TGROWTH> -500),], family = poisson)
LLpoisC1 <- glmer(LEAF_LIFE~climPC1 + climPC2 + climPC3 + (climPC1|SP.ID) , traits.common[which(traits.common$BIOST_TGROWTH> -500),], family = poisson)
LLpoisAll <- glmer(LEAF_LIFE~climPC1 + climPC2 + climPC3 + scale(BIOST_TGROWTH) + (scale(BIOST_TGROWTH)|SP.ID)  + (climPC1|SP.ID) , traits.common[which(traits.common$BIOST_TGROWTH> -500),], family = poisson)
AIC(LLpoisG, LLpoisC1, LLpoisAll)
## TAKEHOME: LEAF_LIFE is WAAAY more strongly related to climate than growth


## model SLA linearly
SLAmodclim <- lmer(SLA_HSA~climPC1*climPC2  + (climPC1 + climPC2 |SP.ID) + (1|PLOT_ID), traits.common[which(traits.common$BIOST_TGROWTH > -500),], REML=F)
SLAmodG <- lmer(SLA_HSA~scale(AG_TGROWTH) + (scale(AG_TGROWTH)|SP.ID), traits.common[which(traits.common$BIOST_TGROWTH > -500),], REML=F)
SLAmodAll <- lmer(SLA_HSA~climPC1*climPC2  + scale(AG_TGROWTH) + (scale(AG_TGROWTH)|SP.ID)+ (climPC1 + climPC2|SP.ID) + (1|PLOT_ID), traits.common[which(traits.common$BIOST_TGROWTH > -500),], REML=F)
AIC(SLAmodclim, SLAmodG, SLAmodAll)
  ### SLA is best described by BOTH..., but really the random slopes for GROWTH


## model LEAF_N linearly
LNmodclim1 <- lmer(LEAF_NITROGEN~climPC1  + (climPC1|SP.ID) + (1|PLOT_ID), traits.common[which(traits.common$BIOST_TGROWTH > -500),], REML=F)
LNmodclim2 <- lmer(LEAF_NITROGEN~climPC1 + climPC2+ (climPC1 + climPC2 |SP.ID) + (1|PLOT_ID), traits.common[which(traits.common$BIOST_TGROWTH > -500),], REML=F)
LNmodclim3 <- lmer(LEAF_NITROGEN~climPC1 + climPC2+ climPC3  + (climPC1 + climPC2 + climPC3|SP.ID) + (1|PLOT_ID), traits.common[which(traits.common$BIOST_TGROWTH > -500),], REML=F)
LNmodclim4 <- lmer(LEAF_NITROGEN~ climPC2+ climPC3  + (climPC3 + climPC2 |SP.ID) + (1|PLOT_ID), traits.common[which(traits.common$BIOST_TGROWTH > -500),], REML=F)
LNmodclim5 <- lmer(LEAF_NITROGEN~ climPC3  + (climPC3|SP.ID) + (1|PLOT_ID), traits.common[which(traits.common$BIOST_TGROWTH > -500),], REML=F)
LNmodclim6 <- lmer(LEAF_NITROGEN~ climPC3 +climPC1*climPC2 + (climPC3|SP.ID) + (1|PLOT_ID), traits.common[which(traits.common$BIOST_TGROWTH > -500),], REML=F)
AIC(LNmodclim1, LNmodclim2, LNmodclim3, LNmodclim4, LNmodclim5, LNmodclim6)

LNmodG <- lmer(LEAF_NITROGEN~scale(AG_TGROWTH) + (scale(AG_TGROWTH)|SP.ID), traits.common[which(traits.common$BIOST_TGROWTH > -500),], REML=F)
LNmodAll <- lmer(LEAF_NITROGEN~climPC1*climPC2+climPC3  + scale(AG_TGROWTH) + (scale(AG_TGROWTH)|SP.ID)+ (climPC3|SP.ID) + (1|PLOT_ID), traits.common[which(traits.common$BIOST_TGROWTH > -500),], REML=F)
AIC(LNmodclim6, LNmodG, LNmodAll)
### LN is best described by BOTH...
  # an interaction between PC1 and PC2, plus a direct effect of PC3 and the random slopes of GROWTH







#-----------------------------------------------------------------------------
####### GROWTH relationships, Leaf Biomass : Stem Biomass and SLA ###############
#-----------------------------------------------------------------------------

### visualization of how Growth relates to Biomass and how I calculated relative growth rate:
ggplot(biomass, aes(x=log(AG_BIOMASS_TREE_TOTAL_AS_CARBON), y=AG_PROD_TREE_TOTAL_AS_CARBON, col=SPP_O1_ABBREV)) + 
  geom_point() + xlab("log(Stand Biomass)") + ylab("Biomass Growth") + 
  geom_smooth(aes(colour=NULL))



# observation: AG_FBIO is strongly related to LAI_O across species and sites. BUT, AG_WBIO is less so
# -> at the whole stand scale, there might be a relationship between the SLA (investment in area vs mass) and Leaf:Stem investment
# *** This would probably best be done at the Community-weighted-trait level?

## Leaf mass to stem mass ratio for these plots ###
traits.common$LtoS <- with(traits.common, AG_FBIO/AG_WBIO)
traits.common$RGR.old <- with(traits.common, AG_TGROWTH/AG_TBIO)
traits.common$lnRGR.old <- with(traits.common, log(AG_TGROWTH)/AG_TBIO)
traits.common$RGR[which(traits.common$RGR>2)] <- NA

plot(AG_TGROWTH~HEIGHTC_m, traits.common, col=FOREST_TYPE)
for(i in unique(traits.common$FOREST_TYPE)[-5]){
  mod <- lm(AG_TGROWTH~AG_TBIO, traits.common[which(traits.common$FOREST_TYPE==i),])
  lines(fitted(mod)~na.omit(traits.common$AG_TBIO[which(traits.common$FOREST_TYPE==i)]))
}
abline(lm(AG_TGROWTH~HEIGHTC_m, traits.common), lwd=2)
  ## NOTE: Height explains 0.079 of GROWTH
  #         TBIO explains .181 of GROWTH (or 24% in biomass)
  #         FBIO explains .36 of GROWTH
  #         LAI_O explains .48 of GROWTH
# I'm not entirely sure how to swing this...

Mypairs(traits.common[,c("AG_TGROWTH", "BIOST_TGROWTH", "BIOST_TGROWTHgam", "SLA_HSA","LEAF_LIFE","LEAF_CARBON","LEAF_NITROGEN")])



######## First pass with regular Biomass and non-LES (logged) traits #######
palette(mypal)
quartz(width=5, height=4)
par(mfrow=c(2,2), mar=c(2,4,0,0), oma=c(2,0,1,1), mgp=c(2.5,1,0))
## SLA vs Growth
plot(SLA_HSA~BIOST_TGROWTHgam, traits.common, pch=16, col=SP.ID, cex=.5)
for(i in levels(traits.common$SP.ID)){
  tmp <- subset(traits.common, subset=SP.ID==i)
  mod <- lm(SLA_HSA~BIOST_TGROWTHgam, tmp)
  lines(fitted(mod)~tmp$BIOST_TGROWTHgam[which(tmp$BIOST_TGROWTHgam>-400 & tmp$SLA_HSA>0)], col=c(mypal, mypal[1])[which(levels(traits.common$SP.ID)==i)], lwd=2)
}
mtext(3, text = "p=0.4/0.186")
# either with only dominants
tmp <- traits.common[which(traits.common$BIOST_TGROWTHgam> -500 & as.character(traits.common$FOREST_TYPE)==traits.common$SP.ID),]
# or with all spp and plots
tmp <- traits.common[which(traits.common$BIOST_TGROWTHgam> -500),]
STHmod1 <- lmer(SLA_HSA~scale(BIOST_TGROWTHgam) + (scale(BIOST_TGROWTHgam)|SP.ID), tmp)
STHmod1null <- lmer(SLA_HSA~1 + (scale(BIOST_TGROWTHgam)|SP.ID), tmp)
STHmod1null1 <- lmer(SLA_HSA~1 + (1|SP.ID), tmp)
anova(STHmod1null, STHmod1) # ns p=0.4, just for dominant spp p=0.186
anova(STHmod1null, STHmod1null1) # significnat
#sjp.lmer(STHmod1)
#sjp.lmer(STHmod1, type="rs.ri")

## LEAF LIFE vs growth
plot(log(LEAF_LIFE)~logBIOST_TGROWTHgam, traits.common, pch=16, col=SP.ID, cex=.5)
for(i in levels(traits.common$SP.ID)[-c(6,7,8)]){
  tmp <- subset(traits.common, subset=SP.ID==i)
  mod <- lm(LEAF_LIFE~BIOST_TGROWTHgam, tmp)
  lines(fitted(mod)~tmp$BIOST_TGROWTHgam[which(tmp$BIOST_TGROWTHgam>-400 & tmp$LEAF_LIFE>0)], col=c(mypal, mypal[1])[which(levels(traits.common$SP.ID)==i)], lwd=2)
}
mtext(3, text = "p=0.04/0.024")
# either with only dominants
tmp <- traits.common[which(traits.common$BIOST_TGROWTHgam> -500 & as.character(traits.common$FOREST_TYPE)==traits.common$SP.ID),]
# or with all spp and plots
tmp <- traits.common[which(traits.common$BIOST_TGROWTHgam> -500),]

STHmod2 <- lmer(LEAF_LIFE~scale(BIOST_TGROWTHgam) + (scale(BIOST_TGROWTHgam)|SP.ID), tmp)
STHmod2null <- lmer(LEAF_LIFE~1 + (scale(BIOST_TGROWTHgam)|SP.ID), tmp)
STHmod2null1 <- lmer(LEAF_LIFE~1 + (1|SP.ID), tmp)
anova(STHmod2null, STHmod2) # marginally significant p=0.047, only dominants (p=0.0236)
anova(STHmod2null, STHmod2null1) # significnat
#sjp.lmer(STHmod2)
#sjp.lmer(STHmod2, type="rs.ri")


## LEAF NITROGEN vs growth
plot(LEAF_NITROGEN~BIOST_TGROWTHgam, traits.common, pch=16, col=SP.ID, cex=0.5)
for(i in levels(traits.common$SP.ID)){
  tmp <- subset(traits.common, subset=SP.ID==i)
  mod <- lm(LEAF_NITROGEN~BIOST_TGROWTHgam, tmp)
  lines(fitted(mod)~tmp$BIOST_TGROWTHgam[which(tmp$BIOST_TGROWTHgam>-400 & tmp$LEAF_NITROGEN>0)], col=c(mypal, mypal[1])[which(levels(traits.common$SP.ID)==i)], lwd=2)
}
mtext(3, text = "p<0.0001/0.005",line = -1)
# either with only dominants
tmp <- traits.common[which(traits.common$BIOST_TGROWTHgam> -500 & as.character(traits.common$FOREST_TYPE)==traits.common$SP.ID),]
# or with all spp and plots
tmp <- traits.common[which(traits.common$BIOST_TGROWTHgam> -500),]

STHmod3 <- lmer(LEAF_NITROGEN~scale(BIOST_TGROWTHgam) + (scale(BIOST_TGROWTHgam)|SP.ID), tmp)
STHmod3null <- lmer(LEAF_NITROGEN~1 + (scale(BIOST_TGROWTHgam)|SP.ID), tmp)
STHmod3null1 <- lmer(LEAF_NITROGEN~1 + (1|SP.ID), tmp)
anova(STHmod3null, STHmod3) # significant p=0.0001, only dominants p=0.005
anova(STHmod3null, STHmod3null1) # significnat
#sjp.lmer(STHmod3)
#sjp.lmer(STHmod3, type="rs.ri")



## LEAF CARBON vs growth
plot(LEAF_CARBON~BIOST_TGROWTHgam, traits.common, pch=16, col=SP.ID, cex=0.5)
for(i in levels(traits.common$SP.ID)){
  tmp <- subset(traits.common, subset=SP.ID==i)
  mod <- lm(LEAF_CARBON~BIOST_TGROWTHgam, tmp)
  lines(fitted(mod)~tmp$BIOST_TGROWTHgam[which(tmp$BIOST_TGROWTHgam>-400 & tmp$LEAF_CARBON>0)], col=c(mypal, mypal[1])[which(levels(traits.common$SP.ID)==i)], lwd=2)
}
mtext(3, text = "p=0.05/0.061",line = -1)
# either with only dominants
tmp <- traits.common[which(traits.common$BIOST_TGROWTHgam> -500 & as.character(traits.common$FOREST_TYPE)==traits.common$SP.ID),]
# or with all spp and plots
tmp <- traits.common[which(traits.common$BIOST_TGROWTHgam> -500),]

STHmod4 <- lmer(LEAF_CARBON~scale(BIOST_TGROWTHgam) + (scale(BIOST_TGROWTHgam)|SP.ID), tmp)
STHmod4null <- lmer(LEAF_CARBON~1 + (scale(BIOST_TGROWTHgam)|SP.ID), tmp)
STHmod4null1 <- lmer(LEAF_CARBON~1 + (1|SP.ID), tmp)
anova(STHmod4null, STHmod4) # marginally significant p=0.047, dominants p=0.061
anova(STHmod4null, STHmod4null1) # significnat
#sjp.lmer(STHmod4)
#sjp.lmer(STHmod4, type="rs.ri")




############## Second Pass with log(Biomass) and LES traits
quartz(width=5, height=4)
par(mfrow=c(2,2), mar=c(2,4,0,0), oma=c(2,0,1,1), mgp=c(2.5,1,0))
## SLA vs Growth
plot(SLA_HSA~logBIOST_TGROWTHgam, traits.common, pch=16, col=SP.ID, cex=.5)
for(i in levels(traits.common$SP.ID)){
  tmp <- subset(traits.common, subset=SP.ID==i)
  mod <- lm(SLA_HSA~logBIOST_TGROWTHgam, tmp)
  lines(fitted(mod)~tmp$logBIOST_TGROWTHgam[which(tmp$logBIOST_TGROWTHgam>-400 & tmp$SLA_HSA>0)], col=c(mypal, mypal[1])[which(levels(traits.common$SP.ID)==i)], lwd=2)
}
mtext(3, text = "p=0.4/0.186")
# either with only dominants
tmp <- traits.common[which(traits.common$logBIOST_TGROWTHgam> -500 & as.character(traits.common$FOREST_TYPE)==traits.common$SP.ID),]
# or with all spp and plots
tmp <- traits.common[which(traits.common$logBIOST_TGROWTHgam> -500),]
STHmod1 <- lmer(SLA_HSA~scale(logBIOST_TGROWTHgam) + (scale(logBIOST_TGROWTHgam)|SP.ID), tmp)
STHmod1null <- lmer(SLA_HSA~1 + (scale(logBIOST_TGROWTHgam)|SP.ID), tmp)
STHmod1null1 <- lmer(SLA_HSA~1 + (1|SP.ID), tmp)
anova(STHmod1null, STHmod1) # ns p=0.4, just for dominant spp p=0.186
anova(STHmod1null, STHmod1null1) # significnat
#sjp.lmer(STHmod1)
#sjp.lmer(STHmod1, type="rs.ri")

## LEAF LIFE vs growth
plot(log(LEAF_LIFE)~logBIOST_TGROWTHgam, traits.common, pch=16, col=SP.ID, cex=.5)
for(i in levels(traits.common$SP.ID)[-c(6,7,8)]){
  tmp <- subset(traits.common, subset=SP.ID==i)
  mod <- lm(log(LEAF_LIFE)~logBIOST_TGROWTHgam, tmp)
  lines(fitted(mod)~tmp$logBIOST_TGROWTHgam[which(tmp$logBIOST_TGROWTHgam>-400 & tmp$LEAF_LIFE>0)], col=c(mypal, mypal[1])[which(levels(traits.common$SP.ID)==i)], lwd=2)
}
mtext(3, text = "p=0.04/0.024")
# either with only dominants
tmp <- traits.common[which(traits.common$logBIOST_TGROWTHgam> -500 & as.character(traits.common$FOREST_TYPE)==traits.common$SP.ID),]
# or with all spp and plots
tmp <- traits.common[which(traits.common$logBIOST_TGROWTHgam> -500),]

STHmod2 <- lmer(LEAF_LIFE~scale(logBIOST_TGROWTHgam) + (scale(logBIOST_TGROWTHgam)|SP.ID), tmp)
STHmod2null <- lmer(LEAF_LIFE~1 + (scale(logBIOST_TGROWTHgam)|SP.ID), tmp)
STHmod2null1 <- lmer(LEAF_LIFE~1 + (1|SP.ID), tmp)
anova(STHmod2null, STHmod2) # marginally significant p=0.047, only dominants (p=0.0236)
anova(STHmod2null, STHmod2null1) # significnat
#sjp.lmer(STHmod2)
#sjp.lmer(STHmod2, type="rs.ri")


## LEAF NITROGEN vs growth
plot(LEAF_NITROGEN~logBIOST_TGROWTHgam, traits.common, pch=16, col=SP.ID, cex=0.5)
for(i in levels(traits.common$SP.ID)){
  tmp <- subset(traits.common, subset=SP.ID==i)
  mod <- lm(LEAF_NITROGEN~logBIOST_TGROWTHgam, tmp)
  lines(fitted(mod)~tmp$logBIOST_TGROWTHgam[which(tmp$logBIOST_TGROWTHgam>-400 & tmp$LEAF_NITROGEN>0)], col=c(mypal, mypal[1])[which(levels(traits.common$SP.ID)==i)], lwd=2)
}
mtext(3, text = "p<0.0001/0.005",line = -1)
# either with only dominants
tmp <- traits.common[which(traits.common$logBIOST_TGROWTHgam> -500 & as.character(traits.common$FOREST_TYPE)==traits.common$SP.ID),]
# or with all spp and plots
tmp <- traits.common[which(traits.common$logBIOST_TGROWTHgam> -500),]

STHmod3 <- lmer(LEAF_NITROGEN~scale(logBIOST_TGROWTHgam) + (scale(logBIOST_TGROWTHgam)|SP.ID), tmp)
STHmod3null <- lmer(LEAF_NITROGEN~1 + (scale(logBIOST_TGROWTHgam)|SP.ID), tmp)
STHmod3null1 <- lmer(LEAF_NITROGEN~1 + (1|SP.ID), tmp)
anova(STHmod3null, STHmod3) # significant p=0.0001, only dominants p=0.005
anova(STHmod3null, STHmod3null1) # significnat
#sjp.lmer(STHmod3)
#sjp.lmer(STHmod3, type="rs.ri")



## LEAF CARBON vs growth
plot(LEAF_CARBON~logBIOST_TGROWTHgam, traits.common, pch=16, col=SP.ID, cex=0.5)
for(i in levels(traits.common$SP.ID)){
  tmp <- subset(traits.common, subset=SP.ID==i)
  mod <- lm(LEAF_CARBON~logBIOST_TGROWTHgam, tmp)
  lines(fitted(mod)~tmp$logBIOST_TGROWTHgam[which(tmp$logBIOST_TGROWTHgam>-400 & tmp$LEAF_CARBON>0)], col=c(mypal, mypal[1])[which(levels(traits.common$SP.ID)==i)], lwd=2)
}
mtext(3, text = "p=0.05/0.061",line = -1)
# either with only dominants
tmp <- traits.common[which(traits.common$logBIOST_TGROWTHgam> -500 & as.character(traits.common$FOREST_TYPE)==traits.common$SP.ID),]
# or with all spp and plots
tmp <- traits.common[which(traits.common$logBIOST_TGROWTHgam> -500),]

STHmod4 <- lmer(LEAF_CARBON~scale(logBIOST_TGROWTHgam) + (scale(logBIOST_TGROWTHgam)|SP.ID), tmp)
STHmod4null <- lmer(LEAF_CARBON~1 + (scale(logBIOST_TGROWTHgam)|SP.ID), tmp)
STHmod4null1 <- lmer(LEAF_CARBON~1 + (1|SP.ID), tmp)
anova(STHmod4null, STHmod4) # marginally significant p=0.047, dominants p=0.061
anova(STHmod4null, STHmod4null1) # significnat
#sjp.lmer(STHmod4)
#sjp.lmer(STHmod4, type="rs.ri")








plot(SLA_HSA~LtoS, traits.common[which(traits.common$dominance>90),], col=SP.ID, xlim=c(0,1)) # note: there are a couple weird plots to be removed
for(i in unique(traits.common$SP.ID)){
  abline(lm(SLA_HSA~LtoS, traits.common, subset=SP.ID==i&dominance>90), col=traits.common$SP.ID[which(traits.common$SP.ID==i)])
}
abline(lm(SLA_HSA~LtoS, traits.common, subset=dominace>90), lwd=3)
    ## Hmm. less of a pattern there than I would have thought...

## SLA vs RGR and SLA vs GROWTH
plot(SLA_HSA~RGR, traits.common[which(traits.common$dominance>90),], col=SP.ID, xlim=c(0,.03))
for(i in unique(traits.common$SP.ID)){
  abline(lm(SLA_HSA~RGR, traits.common, subset=SP.ID==i&dominance>90), col=traits.common$SP.ID[which(traits.common$SP.ID==i)])
}
abline(lm(SLA_HSA~RGR, traits.common), lwd=3)
# growth
plot(SLA_HSA~AG_TGROWTH, traits.common, col=SP.ID)
for(i in unique(traits.common$SP.ID)){
  abline(lm(SLA_HSA~AG_TGROWTH, traits.common, subset=SP.ID==i), col=traits.common$SP.ID[which(traits.common$SP.ID==i)])
}
abline(lm(SLA_HSA~AG_TGROWTH, traits.common), lwd=3)
 



####### Third pass with RGR (log growth) ################
#NOTE: I'm not entirely sure how to structure the random effects for these models, because each stand only has one y variable, but multiple x's
# this means I can't do a random intercept per stand. I can do a random slope per stand, but this seems...wrong. things are way less significant with this model

ggplot(traits.common, aes(x=log(LEAF_NITROGEN), y=RGR, col=SP.ID)) + geom_point() + geom_smooth(method="lm", se=F)

ggplot(traits.common, aes(x=log(LMA), y=RGR, col=SP.ID)) + geom_point() + geom_smooth(method="lm", se=F)

ggplot(traits.common, aes(x=log(LLmonths), y=RGR, col=SP.ID)) + geom_point() + geom_smooth(method="lm", se=F)
  ## LL is really the only thing that follows the expected LES traits ('fast' traits == fast RGR)
ggplot(traits.common, aes(x=log(LLmonths), y=AG_TGROWTH, col=SP.ID)) + geom_point() + geom_smooth(method="lm", se=F)
ggplot(traits.common, aes(x=log(LLmonths), y=BIOST_TGROWTH, col=SP.ID)) + geom_point() + geom_smooth(method="lm", se=F)
ggplot(traits.common, aes(x=log(LLmonths), y=BIOST_TGROWTHgam, col=SP.ID)) + geom_point() + geom_smooth(method="lm", se=F)
  ## LL is related to RGR, and BIOST_TGROWTHgam, but not to raw growth...


############## Growth v traits, only in dominants ########################################
### Let's look only at the traits from the dominant species in each plot
traits.mono <- traits[which(as.character(traits$SP.ID)==as.character(traits$FOREST_TYPE)),]
  # lose about half of traits n=756
  # let's also get rid of anything without 3 or more stands
reps <- traits.mono %>% group_by(SP.ID) %>% summarise(nSites = length(unique(PLOT_ID)))
traits.mono <- traits.mono[which(traits.mono$SP.ID %in% reps$SP.ID[which(reps$nSites>2)] & traits.mono$dominance>75),]
  # drops down to 719, 625 have >50% dominance, 578 have >60% dominance, 497 have >75% dominance
traits.mono$SP.ID <- factor(traits.mono$SP.ID) #drop unused species levels


### not really a relationship between RGR and log LMA
options(na.action="na.omit")
modLMAvRGR <- lmer(RGR~log.LMA + (log.LMA|SP.ID), traits.mono, REML=F)
modLMAvRGRnull <- lmer(RGR~1 + (log.LMA|SP.ID), traits.mono, REML=F)
anova(modLMAvRGRnull, modLMAvRGR) # LRT p=0.3557 with full dataset, p=.21 with .mono, p=0.346 with restricted traits.mono, p=0.52 .mono75%, p=0.586 logRGR .mono75%
#summary( lmer(log.LMA~RGR + (RGR|SP.ID), traits.mono))

### nearly marginally significant relationship between RGR and log Nmass
# modNmassvRGR <- lmer(RGR~log.Nmass + (log.Nmass|SP.ID), traits.common)
# modNmassvRGRnull <-lmer(RGR~1 + (log.Nmass|SP.ID), traits.common)
# anova(modNmassvRGRnull, modNmassvRGR) # p = 0.1306
# summary( lmer(log.Nmass~RGR + (RGR|SP.ID), traits.common))
modNareavRGR <- lmer(log(RGR)~log.Narea + (log.Narea|SP.ID), traits.mono, REML=F)
modNareavRGRnull <-lmer(log(RGR)~1 + (log.Narea|SP.ID), traits.mono, REML=F)
anova(modNareavRGRnull, modNareavRGR) # p = 0.054, # p=0.8 with traits.mono, p=0.543 with restricted traits.mono, p=0.29 .mono75%, p=0.046 logRGR .mono75%
summary( lmer(log.Narea~RGR + (RGR|SP.ID), traits.mono))


### significant relationship between log(RGR) and log LL
modLLvRGR <- lmer(log(RGR)~log.LL + (log.LL|SP.ID), traits.mono)
modLLvRGRnull <- lmer(log(RGR)~1 + (log.LL|SP.ID), traits.mono)
anova(modLLvRGRnull, modLLvRGR) # p = 0.0058 all traits.common, p=0.032 with traits.mono ,p= 0.047 with restricted traits.mono, p=0.005 .mono75%
summary( lmer(log.LL~RGR + (RGR|SP.ID), traits.mono))
modLLvRGR <- glmer(RGR~log.LL + (log.LL|SP.ID), traits.mono,family=(gaussian(link='log')))
modLLvRGRnull <- glmer(RGR~1 + (log.LL|SP.ID), traits.mono,family=(gaussian(link='log')))
anova(modLLvRGRnull, modLLvRGR) # p = 0.0058 all traits.common, p=0.032 with traits.mono ,p= 0.047 with restricted traits.mono, p=0.005 .mono75%



## with different biomass...
modLMAvRGR <- lmer(BIOST_TGROWTHgam~log.LMA + (log.LMA|SP.ID), traits.mono)
modLMAvRGRnull <- lmer(BIOST_TGROWTHgam~1 + (log.LMA|SP.ID), traits.mono)
anova(modLMAvRGRnull, modLMAvRGR) # LRT p=0.78 .common, p=0.944 .mono, p=0.75 with only spp w/ >=3 plots, p=0.1427 .mono75%

modNareavRGR <- lmer(BIOST_TGROWTHgam~log.Narea + (log.Narea|SP.ID), traits.mono)
modNareavRGRnull <-lmer(BIOST_TGROWTHgam~1 + (log.Narea|SP.ID), traits.mono)
anova(modNareavRGRnull, modNareavRGR) # p=0.0478 .common,  p = 0.044 # p=0.105 with only spp w >=3 plots, p=0.098 .mono75%

modLLvRGR <- lmer(BIOST_TGROWTHgam~log.LL + (log.LL|SP.ID), traits.mono)
modLLvRGRnull <- lmer(BIOST_TGROWTHgam~1 + (log.LL|SP.ID), traits.mono)
anova(modLLvRGRnull, modLLvRGR) # p=0.0146 .common,  p=0.0121 # p=0.0752 with only spp w >=3 plots, p=0.05224 .mono 75%




#_____________________________________________________________________
#### FIGURE: RGR versus log TRAITS #################################
#_____________________________________________________________________

## note: need to run models up above


quartz(width=7.1, height=3)
par(mar=c(4,2,1.2,1), oma=c(0,2,0,0),mfrow=c(1,3), cex=1.1,mgp=c(2.5,1,0))

### plotting Leaf Lifespan v RGR
f0LL <- predict(modLLvRGR, re.form=NA, type="response")
f1LL <- predict(modLLvRGR, type="response")
# sort lma values so lines draw right
I <- order(traits.mono$log.LL[-which(is.na(traits.mono$RGR) | is.na(traits.mono$log.LL))])
LLs <- sort(traits.mono$log.LL[-which(is.na(traits.mono$RGR) | is.na(traits.mono$log.LL))])
spid <- traits.mono$SP.ID[-which(is.na(traits.mono$RGR) | is.na(traits.mono$log.LL))][I]
plot(RGR~log.LL, traits.mono, col="gray", pch=16, cex=.4, xlab="log(Leaf Lifespan)")
mtext("Relative Growth Rate", side=2, line=2.5, cex=1.1)
for(i in levels(spid)){
  lines(LLs[which(spid==i)], f1LL[I][which(spid==i)],  col=spid[which(spid==i)], lwd=2)
}
lines(LLs, f0LL[I], lwd=4, lty=1)
mtext(text = "p=0.015", side=3, adj = .8, line=-1, font=2 )
mtext("a)", side=3, adj=0)

### plotting LMA v RGR
f0lma <- predict(modLMAvRGR, re.form=NA)
f1lma <- fitted(modLMAvRGR)
# sort lma values so lines draw right
I <- order(traits.mono$log.LMA[-which(is.na(traits.mono$RGR) | is.na(traits.mono$log.LMA))])
LMAs <- sort(traits.mono$log.LMA[-which(is.na(traits.mono$RGR) | is.na(traits.mono$log.LMA))])
spid <- traits.mono$SP.ID[-which(is.na(traits.mono$RGR) | is.na(traits.mono$log.LMA))][I]
plot(RGR~log.LMA, traits.mono, col="gray", pch=16, cex=.4, xlab="log(LMA)")
for(i in levels(spid)){
  lines(LMAs[which(spid==i)], f1lma[I][which(spid==i)],  col=spid[which(spid==i)], lwd=2)
}
lines(c(min(LMAs), max(LMAs)), c(min(f0lma), max(f0lma)), lwd=4, lty=3)
mtext(text = "p=0.35", side=3, adj = .2, line=-1 )
mtext("b)", side=3, adj=0)

### plotting Narea v RGR
f0Narea <- predict(modNareavRGR, re.form=NA)
f1Narea <- fitted(modNareavRGR)
# sort lma values so lines draw right
I <- order(traits.mono$log.Narea[-which(is.na(traits.mono$RGR) | is.na(traits.mono$log.Narea))])
Nareas <- sort(traits.mono$log.Narea[-which(is.na(traits.mono$RGR) | is.na(traits.mono$log.Narea))])
spid <- traits.mono$SP.ID[-which(is.na(traits.mono$RGR) | is.na(traits.mono$log.Narea))][I]
plot(RGR~log.Narea, traits.mono, col="gray", pch=16, cex=.4, xlab="log(Narea)")
for(i in levels(spid)){
  lines(Nareas[which(spid==i)], f1Narea[I][which(spid==i)],  col=spid[which(spid==i)], lwd=2)
}
lines(Nareas, f0Narea[I], lwd=4, lty=1)
mtext(text = "p=0.054", side=3, adj = .2, line=-1 )
mtext("c)", side=3, adj=0)




####### Figure reploted with traits.mono and biost_growth instead of RGR.
quartz(width=7.1, height=3)
par(mar=c(4,2,1.2,1), oma=c(0,2,0,0),mfrow=c(1,3), cex=1.1,mgp=c(2.5,1,0))

### plotting Leaf Lifespan v RGR
f0LL <- predict(modLLvRGR, re.form=NA)
f1LL <- fitted(modLLvRGR)
# sort lma values so lines draw right
I <- order(traits.mono$log.LL[-which(is.na(traits.mono$BIOST_TGROWTHgam) | is.na(traits.mono$log.LL))])
LLs <- sort(traits.mono$log.LL[-which(is.na(traits.mono$BIOST_TGROWTHgam) | is.na(traits.mono$log.LL))])
spid <- traits.mono$SP.ID[-which(is.na(traits.mono$BIOST_TGROWTHgam) | is.na(traits.mono$log.LL))][I]
plot(BIOST_TGROWTHgam~log.LL, traits.mono, col="gray", pch=16, cex=.4, xlab="log(Leaf Lifespan)")
mtext("Relative Growth Rate", side=2, line=2.5, cex=1.1)
for(i in levels(spid)){
  lines(LLs[which(spid==i)], f1LL[I][which(spid==i)],  col=spid[which(spid==i)], lwd=2)
}
lines(LLs, f0LL[I], lwd=4, lty=1)
mtext(text = "p=0.015", side=3, adj = .8, line=-1, font=2 )
mtext("a)", side=3, adj=0)

### plotting LMA v BIOST_TGROWTHgam
f0lma <- predict(modLMAvRGR, re.form=NA)
f1lma <- fitted(modLMAvRGR)
# sort lma values so lines draw right
I <- order(traits.mono$log.LMA[-which(is.na(traits.mono$BIOST_TGROWTHgam) | is.na(traits.mono$log.LMA))])
LMAs <- sort(traits.mono$log.LMA[-which(is.na(traits.mono$BIOST_TGROWTHgam) | is.na(traits.mono$log.LMA))])
spid <- traits.mono$SP.ID[-which(is.na(traits.mono$BIOST_TGROWTHgam) | is.na(traits.mono$log.LMA))][I]
plot(BIOST_TGROWTHgam~log.LMA, traits.mono, col="gray", pch=16, cex=.4, xlab="log(LMA)")
for(i in levels(spid)){
  lines(LMAs[which(spid==i)], f1lma[I][which(spid==i)],  col=spid[which(spid==i)], lwd=2)
}
lines(c(min(LMAs), max(LMAs)), c(min(f0lma), max(f0lma)), lwd=4, lty=3)
mtext(text = "p=0.35", side=3, adj = .2, line=-1 )
mtext("b)", side=3, adj=0)

### plotting Narea v BIOST_TGROWTHgam
f0Narea <- predict(modNareavRGR, re.form=NA)
f1Narea <- fitted(modNareavRGR)
# sort lma values so lines draw right
I <- order(traits.mono$log.Narea[-which(is.na(traits.mono$BIOST_TGROWTHgam) | is.na(traits.mono$log.Narea))])
Nareas <- sort(traits.mono$log.Narea[-which(is.na(traits.mono$BIOST_TGROWTHgam) | is.na(traits.mono$log.Narea))])
spid <- traits.mono$SP.ID[-which(is.na(traits.mono$BIOST_TGROWTHgam) | is.na(traits.mono$log.Narea))][I]
plot(BIOST_TGROWTHgam~log.Narea, traits.mono, col="gray", pch=16, cex=.4, xlab="log(Narea)")
for(i in levels(spid)){
  lines(Nareas[which(spid==i)], f1Narea[I][which(spid==i)],  col=spid[which(spid==i)], lwd=2)
}
lines(Nareas, f0Narea[I], lwd=4, lty=1)
mtext(text = "p=0.054", side=3, adj = .2, line=-1 )
mtext("c)", side=3, adj=0)


#_______________________________________________________________________
############ ****** Trait v RGR (final with plot means) *** ##################
#_______________________________________________________________________
# ok, I realized that all of the above are taking plot level averages for RGR and I really should be taking plot level averages for traits to.

spp.plot.traits <- traits.common %>% group_by(SP.ID, PLOT_ID) %>% summarise(nsample = n(), mSLA = mean(SLA_HSA, na.rm=T), nSLA = n()- length(which(is.na(SLA_HSA))), mCN = mean(LEAF_CN, na.rm=T), nCN = n()- length(which(is.na(LEAF_CN)))
                                                                     , mLIFE = mean(LEAF_LIFE, na.rm=T), nLIFE=n()- length(which(is.na(LEAF_LIFE))), mNITROGEN = mean(LEAF_NITROGEN, na.rm=T), nplots =length(unique(PLOT_ID))
                                                                     , mLMA = mean(LMA_PSA, na.rm=T), mLLmonths= mean(LLmonths, na.rm=T), mNarea=mean(Narea, na.rm=T)
                                                                     , mlog.LMA = mean(log.LMA, na.rm=T), mlog.LL = mean(log.LL, na.rm=T), mlog.Narea=mean(log.Narea, na.rm=T), mlog.Nmass=mean(log.Nmass, na.rm=T), RGR = unique(RGRdom), stGrowth = mean(stGrowthdom, na.rm=T)
                                                                     , climPC1 = unique(climPC1), climPC2 = unique(climPC2), climPC3=unique(climPC3), soil_N = unique(soil_N), soil_pH = unique(soil_pH),ASA=unique(ASA) )
# Note: I used LMA_PSA in other trait analysis in TaxonomicAnalysis. So I've switched to it here
# let's just pull out the common dominant species in an easy to name df
goodspp <- spp.plot.traits %>% filter(RGR>0) %>% group_by(SP.ID) %>% summarise(nplots= length(unique(PLOT_ID)))
plotavs <- spp.plot.traits%>% subset(SP.ID %in% goodspp$SP.ID[which(goodspp$nplots>2)] & RGR>0)
plotavs$SP.ID <- factor(plotavs$SP.ID)
# just to make it look better, I'm going to move the common species to the front of the queue
levels(plotavs$SP.ID) <- list(PSEMEN="PSEMEN", PINPON="PINPON",ABICON="ABICON",ABIGRA="ABIGRA",PINJEF="PINJEF",PINCON="PINCON", JUNOCC="JUNOCC", TSUHET="TSUHET", ABIPRO='ABIPRO')


### significant relationship between mlog(RGR) and mlog LL
  # haven't rerun logs since updated plotavs
modLLvRGR <- lmer(RGR~mlog.LL + (mlog.LL|SP.ID), plotavs)
modLLvRGRnull <- lmer(RGR~1 + (mlog.LL|SP.ID), plotavs)
anova(modLLvRGRnull, modLLvRGR) #p= 0.0007 / log = 0.01146      IGNORE -> p = 0.0058 all traits.common, p=0.032 with plotavs ,p= 0.047 with restricted plotavs, p=0.005 .mono75%
# using an exponential link but a gaussian error distribution
#modLLvRGR <- glmer(RGR~scale(mlog.LL) + (scale(mlog.LL)|SP.ID), plotavs,family=(gaussian(link='log')))
#modLLvRGRnull <- glmer(RGR~1 + (mlog.LL|SP.ID), plotavs,family=(gaussian(link='log')))
#anova(modLLvRGRnull, modLLvRGR) # 0.0116     IGNORE ->p = 0.0058 all traits.common, p=0.032 with plotavs ,p= 0.047 with restricted plotavs, p=0.005 .mono75%


modLMAvRGR <- lmer(RGR~mlog.LMA + (mlog.LMA|SP.ID), plotavs, REML=F)
modLMAvRGRnull <- lmer(RGR~1 + (mlog.LMA|SP.ID), plotavs, REML=F)
anova(modLMAvRGRnull, modLMAvRGR) #p=0.554/ log = 0.488      IGNORE -> LRT p=0.3557 with full dataset, p=.21 with .mono, p=0.346 with restricted plotavs, p=0.52 .mono75%, p=0.586 mlogRGR .mono75%

modNareavRGR <- lmer(RGR~mlog.Narea + (mlog.Narea|SP.ID), plotavs, REML=F)
modNareavRGRnull <-lmer(RGR~1 + (scale(mlog.Narea)|SP.ID), plotavs, REML=F)
anova(modNareavRGRnull, modNareavRGR) #p=0.18   / log=0.24   IGNORE -> p = 0.054, # p=0.8 with plotavs, p=0.543 with restricted plotavs, p=0.29 .mono75%, p=0.046 mlogRGR .mono75%




## with different biomass...
modLLvGrowth <- lmer(stGrowth~mlog.LL + (mlog.LL|SP.ID), plotavs, REML=F)
modLLvGrowthnull <- lmer(stGrowth~1 + (mlog.LL|SP.ID), plotavs, REML=F)
anova(modLLvGrowthnull, modLLvGrowth) # p=0.09677     IGNORE -> p=0.0146 .common,  p=0.0121 # p=0.0752 with only spp w >=3 plots, p=0.05224 .mono 75%

modLMAvGrowth <- lmer(stGrowth~mlog.LMA + (mlog.LMA|SP.ID), plotavs, REML=F)
modLMAvGrowthnull <- lmer(stGrowth~1 + (mlog.LMA|SP.ID), plotavs, REML=F)
anova(modLMAvGrowthnull, modLMAvGrowth) #p=0.747      IGNORE ->LRT p=0.78 .common, p=0.944 .mono, p=0.75 with only spp w/ >=3 plots, p=0.1427 .mono75%

modNareavGrowth <- lmer(stGrowth~mlog.Narea + (mlog.Narea|SP.ID), plotavs, REML=F)
modNareavGrowthnull <-lmer(stGrowth~1 + (mlog.Narea|SP.ID), plotavs, REML=F)
anova(modNareavGrowthnull, modNareavGrowth) #p=0.1884     IGNORE ->p=0.0478 .common,  p = 0.044 # p=0.105 with only spp w >=3 plots, p=0.098 .mono75%






####### Figure 6: reploted with traits.mono and biost_growth instead of RGR.
quartz(width=7.1, height=6)
pdf(file = "Traits-v-Growth_v1.pdf",width = 7.1, height=6)
par(mar=c(2.5,2,1.2,1), oma=c(2,2,0,0),mfrow=c(2,3), cex=1.1,mgp=c(2.5,1,0))

#### top three panels: RGR
f0LL <- predict(modLLvRGR, re.form=NA)
f1LL <- fitted(modLLvRGR)
# sort lma values so lines draw right
I <- order(plotavs$mlog.LL[-which(is.na(plotavs$RGR) | is.na(plotavs$mlog.LL))])
LLs <- sort(plotavs$mlog.LL[-which(is.na(plotavs$RGR) | is.na(plotavs$mlog.LL))])
spid <- plotavs$SP.ID[-which(is.na(plotavs$RGR) | is.na(plotavs$mlog.LL))][I]
plot(RGR~mlog.LL, plotavs, col=SP.ID, pch=16, cex=.4, xlab="")#log(Leaf Lifespan)")
mtext("Relative Growth Rate", side=2, line=2.5, cex=1.1)
for(i in levels(spid)){
  lines(LLs[which(spid==i)], f1LL[I][which(spid==i)],  col=spid[which(spid==i)], lwd=2)
}
lines(LLs, f0LL[I], lwd=4, lty=1)
mtext(text = "p=0.0007", side=3, adj = .8, line=-1, font=1 )
mtext("a)", side=3, adj=0)

### plotting LMA v RGR
f0lma <- predict(modLMAvRGR, re.form=NA)
f1lma <- fitted(modLMAvRGR)
# sort lma values so lines draw right # no NA values in LMA
I <- order(plotavs$mlog.LMA)
LMAs <- sort(plotavs$mlog.LMA)
spid <- plotavs$SP.ID[I]
plot(RGR~mlog.LMA, plotavs, col=SP.ID, pch=16, cex=.4, xlab="")#log(LMA)")
for(i in levels(spid)){
  lines(LMAs[which(spid==i)], f1lma[I][which(spid==i)],  col=spid[which(spid==i)], lwd=2)
}
lines(c(min(LMAs), max(LMAs)), c(min(f0lma), max(f0lma)), lwd=4, lty=3)
mtext(text = "p=0.55", side=3, adj = .2, line=-1 )
mtext("b)", side=3, adj=0)

### plotting Narea v RGR
f0Narea <- predict(modNareavRGR, re.form=NA)
f1Narea <- fitted(modNareavRGR)
# sort lma values so lines draw right
I <- order(plotavs$mlog.Narea[-which(is.na(plotavs$RGR) | is.na(plotavs$mlog.Narea))])
Nareas <- sort(plotavs$mlog.Narea[-which(is.na(plotavs$RGR) | is.na(plotavs$mlog.Narea))])
spid <- plotavs$SP.ID[-which(is.na(plotavs$RGR) | is.na(plotavs$mlog.Narea))][I]
plot(RGR~mlog.Narea, plotavs, col=SP.ID, pch=16, cex=.4, xlab="")#log(Narea)")
for(i in levels(spid)){
  lines(Nareas[which(spid==i)], f1Narea[I][which(spid==i)],  col=spid[which(spid==i)], lwd=2)
}
lines(Nareas, f0Narea[I], lwd=4, lty=3)
mtext(text = "p=0.18", side=3, adj = .2, line=-1 )
mtext("c)", side=3, adj=0)




### plotting Leaf Lifespan v RGR
f0LL <- predict(modLLvGrowth, re.form=NA)
f1LL <- fitted(modLLvGrowth)
# sort lma values so lines draw right
I <- order(plotavs$mlog.LL[-which(is.na(plotavs$stGrowth) | is.na(plotavs$mlog.LL))])
LLs <- sort(plotavs$mlog.LL[-which(is.na(plotavs$stGrowth) | is.na(plotavs$mlog.LL))])
spid <- plotavs$SP.ID[-which(is.na(plotavs$stGrowth) | is.na(plotavs$mlog.LL))][I]
plot(stGrowth~mlog.LL, plotavs, col=SP.ID, pch=16, cex=.4, xlab="log(Leaf Lifespan)")
mtext("Standardized Growth", side=2, line=2.5, cex=1.1)

for(i in levels(spid)){
  lines(LLs[which(spid==i)], f1LL[I][which(spid==i)],  col=spid[which(spid==i)], lwd=2)
}
lines(LLs, f0LL[I], lwd=4, lty=1)
mtext(text = "p=0.097", side=3, adj = .8, line=-1, font=1 )
mtext("d)", side=3, adj=0)
mtext("log(Leaf Lifespan)", side=1, line=2.5, cex=1.1)

### plotting LMA v stGrowth
f0lma <- predict(modLMAvGrowth, re.form=NA)
f1lma <- fitted(modLMAvGrowth)
# sort lma values so lines draw right # no NA values in LMA
I <- order(plotavs$mlog.LMA)
LMAs <- sort(plotavs$mlog.LMA)
spid <- plotavs$SP.ID[I]
plot(stGrowth~mlog.LMA, plotavs, col=SP.ID, pch=16, cex=.4, xlab="log(LMA)")
for(i in levels(spid)){
  lines(LMAs[which(spid==i)], f1lma[I][which(spid==i)],  col=spid[which(spid==i)], lwd=2)
}
lines(c(min(LMAs), max(LMAs)), c(min(f0lma), max(f0lma)), lwd=4, lty=3)
mtext(text = "p=0.74", side=3, adj = .2, line=-1 )
mtext("e)", side=3, adj=0)
mtext("log(LMA)", side=1, line=2.5, cex=1.1)

### plotting Narea v stGrowth
f0Narea <- predict(modNareavGrowth, re.form=NA)
f1Narea <- fitted(modNareavGrowth)
# sort lma values so lines draw right
I <- order(plotavs$mlog.Narea[-which(is.na(plotavs$stGrowth) | is.na(plotavs$mlog.Narea))])
Nareas <- sort(plotavs$mlog.Narea[-which(is.na(plotavs$stGrowth) | is.na(plotavs$mlog.Narea))])
spid <- plotavs$SP.ID[-which(is.na(plotavs$stGrowth) | is.na(plotavs$mlog.Narea))][I]
plot(stGrowth~mlog.Narea, plotavs, col=SP.ID, pch=16, cex=.4, xlab="log(Narea)")
for(i in levels(spid)){
  lines(Nareas[which(spid==i)], f1Narea[I][which(spid==i)],  col=spid[which(spid==i)], lwd=2)
}
lines(Nareas, f0Narea[I], lwd=4, lty=3)
mtext(text = "p=0.19", side=3, adj = .2, line=-1 )
mtext("f)", side=3, adj=0)
mtext("log(Narea)", side=1, line=2.5, cex=1.1)













#_______________________________________________________________________
############### ****Trait - environment relationships**** ##############
#_______________________________________________________________________

#### LMA vs LL in PSME vs precip #####
  # take home seems to be: 
#             - relationship may become more negative at higher aridities.
#             - LL is higher at the same LMA in drier places (direct intercept offset)
#             - dry values differ a lot from MVNorm, so SMA probably the best.

ggplot(traits[which(traits$SP.ID=="PSEMEN"),], aes(x=log.LMA, y=log.LL, col=soilmoist.lvl1.mm, size=soilmoist.lvl1.mm, alpha=1/2)) + geom_point() +
  geom_smooth(data=traits[which(traits$SP.ID=="PSEMEN" & traits$soilmoist.lvl1.mm<15),], method='lm', col="black") +
  geom_smooth(data=traits[which(traits$SP.ID=="PSEMEN" & traits$soilmoist.lvl1.mm>15),], method='lm')

pPSME <- ggplot(traits[which(traits$SP.ID=="PSEMEN"),], aes(x=log.LMA, y=log.LL, col=MAP)) + geom_point() +
  geom_smooth(data=traits[which(traits$SP.ID=="PSEMEN" & traits$MAP<1500),], method='lm', col="black") +
  geom_smooth(data=traits[which(traits$SP.ID=="PSEMEN" & traits$MAP>1500),], method='lm')+ ggtitle("PSEMEN")



plot(log.LL~log.LMA, traits[which(traits$SP.ID=="PSEMEN" & traits$MAP<1500),], pch=16, cex=MAP/max(MAP, na.rm=T))
  ## standardized major axis regressions (note: MA loogs FUNKY)
MAwet <- lmodel2(log.LL~log.LMA, traits[which(traits$SP.ID=="PSEMEN" & traits$MAP>1500),]) 
lines.lmodel2(MAwet, method = "SMA", confidence = F )
MAdry <- lmodel2(log.LL~log.LMA, traits[which(traits$SP.ID=="PSEMEN" & traits$MAP<1500),]) 
lines.lmodel2(MAdry, method = "SMA", confidence=F)
## OLS regressions
abline(lm(log.LL~log.LMA,traits[which(traits$SP.ID=="PSEMEN" & traits$MAP>1500),]), lwd=2)
abline(lm(log.LL~log.LMA,traits[which(traits$SP.ID=="PSEMEN" & traits$MAP<1500),]), lwd=1)

plot(LLmonths~LMA, traits[which(traits$SP.ID=="PSEMEN"),], pch=16, cex=MAP/max(MAP, na.rm=T))
abline(lm(LLmonths~LMA,traits[which(traits$SP.ID=="PSEMEN" & traits$MAP<1500),]), lwd=1)
abline(lm(LLmonths~LMA,traits[which(traits$SP.ID=="PSEMEN" & traits$MAP>1500),]), lwd=2)




##### LES LLvLMA relationship as f(MAP) ######


### Look at PSME LL as a f(MAP, MAT) -> could this be a surface?
# for LL, there is definitely a MAP,MAT interaction
plot(MAT~MAP, traits[which(traits$SP.ID=="PSEMEN"),], cex=LMA/max(LMA, na.rm = T)*2)
# a subtle pattern in MAT
plot(LLmonths~MAT, traits[which(traits$SP.ID=="PSEMEN"),])
# a much stronger MAP pattern
plot(LLmonths~MAP, traits[which(traits$SP.ID=="PSEMEN"),])
# decent pattern with VPD
plot(LLmonths~vpd.gy.max, traits[which(traits$SP.ID=="PSEMEN"),])
# not really an stronger pattern with CMI than PPT
plot(LLmonths~cmi.gy.mm, traits[which(traits$SP.ID=="PSEMEN"),])
# 
plot(MAT~MAP, traits[which(traits$SP.ID=="PSEMEN"),], cex=LMA/max(LMA, na.rm = T)*3)

# a quick model selection to verify
mod1 <- lm(log.LL~1, traits[which(traits$SP.ID=="PSEMEN" & traits$MAP>0 & traits$log.LMA >0),]) # there are 5 NAs in MAP
mod2 <- lm(log.LL~log.LMA, traits[which(traits$SP.ID=="PSEMEN"& traits$MAP>0 & traits$log.LMA >0),])
mod3 <- lm(log.LL~MAP, traits[which(traits$SP.ID=="PSEMEN"& traits$MAP>0 & traits$log.LMA >0),])
mod4 <- lm(log.LL~log.LMA + MAP, traits[which(traits$SP.ID=="PSEMEN"& traits$MAP>0 & traits$log.LMA >0),])
mod5 <- lm(log.LL~log.LMA * MAP, traits[which(traits$SP.ID=="PSEMEN"& traits$MAP>0 & traits$log.LMA >0),])
# but only the intercept shift is significant...
AIC(mod1, mod2, mod3, mod4, mod5) # best model only involves MAP
anova(mod4, mod5) # both main effects are obviously significant, but interaction is only marginally significant (p=0.056)


### plotting PSME LL vs LMA by MAP interaction on top of LES LLvLMA with MAP interaction

  # first make predictions for fitting high rain and low rain lines
xs <- seq(.8,3.3, by=.1)
llmod <- lm(log.LL~log.LMA * MAP, LES)
predswet <- predict(llmod, newdata = list(log.LMA=xs, MAP=rep(3000, length(xs)))) 
predsdry <- predict(llmod, newdata = list(log.LMA=xs, MAP=rep(500, length(xs)))) 

llmodPSME <- lm(log.LL~log.LMA_PSA * MAP, traits[which(traits$SP.ID=="PSEMEN"),])
predswetPM <- predict(llmodPSME, newdata = list(log.LMA_PSA=xs, MAP=rep(3000, length(xs)))) 
predsdryPM <- predict(llmodPSME, newdata = list(log.LMA_PSA=xs, MAP=rep(500, length(xs)))) 

  # now plot the data and predictions
plot(log.LL~log.LMA, LES, pch=16, col="grey")
lines(predswet~xs, col="darkblue")
lines(predsdry~xs, col="darkred")

points(log.LL~log.LMA_PSA, traits[which(traits$SP.ID=="PSEMEN"),], pch=16, cex=.9, col=paste0(mypal[3],"66"))
lines(predswetPM~xs, col="darkblue", lwd=3, lty=2)
lines(predsdryPM~xs, col="darkred", lwd=3, lty=2)




#LL for species other than PSME ######

#Ponderosa - way less obvious pattern than PSME
plot(LLmonths~LMA, traits[which(traits$SP.ID=="PINPON"),], cex=MAP/max(MAP, na.rm=T)+1)
# a subtle pattern in MAT
plot(LLmonths~MAT, traits[which(traits$SP.ID=="PINPON"),]); abline(lm(LLmonths~MAT, traits[which(traits$SP.ID=="PINPON"),]))
# a subtle pattern with MAP
plot(LLmonths~MAP, traits[which(traits$SP.ID=="PINPON"),]); abline(lm(LLmonths~MAP, traits[which(traits$SP.ID=="PINPON"),]))
# no pattern with vpd.
plot(LLmonths~vpd.gy.max, traits[which(traits$SP.ID=="PINPON"),]); abline(lm(LLmonths~vpd.gy.max, traits[which(traits$SP.ID=="PINPON"),]))
# not really an stronger pattern with CMI than PPT
plot(LLmonths~cmi.gy.mm, traits[which(traits$SP.ID=="PINPON"),]); abline(lm(LLmonths~cmi.gy.mm, traits[which(traits$SP.ID=="PINPON"),]))
# similarly subtle soilmoist pattern
plot(LLmonths~soilmoist.lvl1.mm, traits[which(traits$SP.ID=="PINPON"),]); abline(lm(LLmonths~soilmoist.lvl1.mm, traits[which(traits$SP.ID=="PINPON"),]))

(pPINPON <- ggplot(traits[which(traits$SP.ID=="PINPON"),], aes(x=log.LMA, y=log.LL, col=MAP)) + geom_point() +
  geom_smooth(data=traits[which(traits$SP.ID=="PINPON" & traits$MAP<1000),], method='lm', col="black") +
  geom_smooth(data=traits[which(traits$SP.ID=="PINPON" & traits$MAP>1000),], method='lm')+ ggtitle("PINPON"))

# the same shifting up of LL for a given LMA at low Precip. But the slope goes from postive to no slope rather than 
# no slope to negative slope...
  # a quick model selection to verify
mod1 <- lm(log.LL~1, traits[which(traits$SP.ID=="PINPON" & traits$MAP>0),]) # there are 5 NAs in MAP
mod2 <- lm(log.LL~log.LMA, traits[which(traits$SP.ID=="PINPON"& traits$MAP>0),])
mod3 <- lm(log.LL~MAP, traits[which(traits$SP.ID=="PINPON"& traits$MAP>0),])
mod4 <- lm(log.LL~log.LMA + MAP, traits[which(traits$SP.ID=="PINPON"& traits$MAP>0),])
mod5 <- lm(log.LL~log.LMA * MAP, traits[which(traits$SP.ID=="PINPON"& traits$MAP>0),])
  # but only the intercept shift is significant...
AIC(mod1, mod2, mod3, mod4, mod5) # best model only involves MAP
anova(mod1, mod3) # MAP has a significant negative effect on LL, p=0.0169


#ABIGRA - Similar pattern to PSME, but not significant...
plot(LLmonths~LMA, traits[which(traits$SP.ID=="ABIGRA"),], cex=MAP/max(MAP, na.rm=T)+1)
# a strange pattern with MAT, would be v-strong without very low T plot
plot(LLmonths~MAT, traits[which(traits$SP.ID=="ABIGRA"),]); abline(lm(LLmonths~MAT, traits[which(traits$SP.ID=="ABIGRA"),]))
# a decent pattern with MAP
plot(LLmonths~MAP, traits[which(traits$SP.ID=="ABIGRA"),]); abline(lm(LLmonths~MAP, traits[which(traits$SP.ID=="ABIGRA"),]))
# decent pattern with vpd
plot(LLmonths~vpd.gy.max, traits[which(traits$SP.ID=="ABIGRA"),]); abline(lm(LLmonths~vpd.gy.max, traits[which(traits$SP.ID=="ABIGRA"),]))
# not really a pattern with cmi
plot(LLmonths~cmi.gy.mm, traits[which(traits$SP.ID=="ABIGRA"),]); abline(lm(LLmonths~cmi.gy.mm, traits[which(traits$SP.ID=="ABIGRA"),]))
# similarly subtle soilmoist pattern
plot(LLmonths~soilmoist.lvl1.mm, traits[which(traits$SP.ID=="ABIGRA"),]); abline(lm(LLmonths~soilmoist.lvl1.mm, traits[which(traits$SP.ID=="ABIGRA"),]))

(pABIGRA <- ggplot(traits[which(traits$SP.ID=="ABIGRA"),], aes(x=log.LMA, y=log.LL, col=MAP)) + geom_point() +
  geom_smooth(data=traits[which(traits$SP.ID=="ABIGRA" & traits$MAP<1000),], method='lm', col="black") +
  geom_smooth(data=traits[which(traits$SP.ID=="ABIGRA" & traits$MAP>1000),], method='lm') + ggtitle("ABIGRA"))

# the same shifting up of LL for a given LMA at low Precip. Similar slope pattern to PSEMEN
# a quick model selection to verify
mod1 <- lm(log.LL~1, traits[which(traits$SP.ID=="ABIGRA" & traits$MAP>0 & traits$log.LMA >0),]) # there are 5 NAs in MAP
mod2 <- lm(log.LL~log.LMA, traits[which(traits$SP.ID=="ABIGRA"& traits$MAP>0 & traits$log.LMA >0),])
mod3 <- lm(log.LL~MAP, traits[which(traits$SP.ID=="ABIGRA"& traits$MAP>0 & traits$log.LMA >0),])
mod4 <- lm(log.LL~log.LMA + MAP, traits[which(traits$SP.ID=="ABIGRA"& traits$MAP>0 & traits$log.LMA >0),])
mod5 <- lm(log.LL~log.LMA * MAP, traits[which(traits$SP.ID=="ABIGRA"& traits$MAP>0 & traits$log.LMA >0),])
# but only the intercept shift is significant...
AIC(mod1, mod2, mod3, mod4, mod5) # best model only involves MAP
anova(mod1, mod2) # actually, only a marginally significant log.LMA effect.
  # This effect becomes very significant with an SMA from lmodel2




#ABICON
plot(LLmonths~LMA, traits[which(traits$SP.ID=="ABICON"),], cex=MAP/max(MAP, na.rm=T)+1)
Mypairs( traits[which(traits$SP.ID=="ABICON"), c("LLmonths","LMA","MAT", "MAP","vpd.gy.max","cmi.gy.mm","soilmoist.lvl1.mm")])
  # best pattern really seems to be with vpd. not really ANY LMA pattern right off the bat... 
  #also vpd looks like a pattern I don't necessarily believe...

ggplot(traits[which(traits$SP.ID=="ABICON"),], aes(x=log.LMA, y=log.LL, col=vpd.gy.max)) + geom_point() +
  geom_smooth(data=traits[which(traits$SP.ID=="ABICON" & traits$vpd.gy.max>1.6),], method='lm', col="black") +
  geom_smooth(data=traits[which(traits$SP.ID=="ABICON" & traits$vpd.gy.max<1.6),], method='lm')  + ggtitle("ABICON")


(pABICON <- ggplot(traits[which(traits$SP.ID=="ABICON"),], aes(x=log.LMA, y=log.LL, col=MAP)) + geom_point() +
  geom_smooth(data=traits[which(traits$SP.ID=="ABICON" & traits$MAP<800),], method='lm', col="black") +
  geom_smooth(data=traits[which(traits$SP.ID=="ABICON" & traits$MAP>800),], method='lm')  + ggtitle("ABICON"))


# the same shifting up of LL for a given LMA at low Precip. Similar slope pattern to PSEMEN
# a quick model selection to verify
mod1 <- lm(log.LL~1, traits[which(traits$SP.ID=="ABICON" & traits$MAP>0 & traits$log.LMA >0),]) # there are 5 NAs in MAP
mod2 <- lm(log.LL~log.LMA, traits[which(traits$SP.ID=="ABICON"& traits$MAP>0 & traits$log.LMA >0),])
mod3 <- lm(log.LL~vpd.gy.max, traits[which(traits$SP.ID=="ABICON"& traits$MAP>0 & traits$log.LMA >0),])
mod4 <- lm(log.LL~log.LMA + vpd.gy.max, traits[which(traits$SP.ID=="ABICON"& traits$MAP>0 & traits$log.LMA >0),])
mod5 <- lm(log.LL~log.LMA * vpd.gy.max, traits[which(traits$SP.ID=="ABICON"& traits$MAP>0 & traits$log.LMA >0),])
# but only the intercept shift is significant...
AIC(mod1, mod2, mod3, mod4, mod5) 
  # as expected, no Precip signal or log.LMA signal.
  # there is a positive vpd signal, though. Higher vpds = higher leaf life spans..???
lmodel2(log.LL~log.LMA, traits[which(traits$SP.ID=="ABICON"),])
  # the SMA major axis slope seems to be pretty solidly negative. But I'm not sure if I buy the MA slope...




#TSUHET
plot(LLmonths~LMA, traits[which(traits$SP.ID=="TSUHET"),], cex=MAP/max(MAP, na.rm=T)+1)
Mypairs( traits[which(traits$SP.ID=="TSUHET"), c("LLmonths","LMA","MAT", "MAP","vpd.gy.max","cmi.gy.mm","soilmoist.lvl1.mm")])
# no pattern with LMA, but strongly INCREASES with vpd. Much more reasonable than ABICON. Maybe dgreases with MAT and soil moisture?
# ALSO, LMA seems to strongly decrease with MAT...

ggplot(traits[which(traits$SP.ID=="TSUHET"),], aes(x=log.LMA, y=log.LL, col=vpd.gy.max)) + geom_point() +
  geom_smooth(data=traits[which(traits$SP.ID=="TSUHET" & traits$vpd.gy.max>1.2),], method='lm', col="black") +
  geom_smooth(data=traits[which(traits$SP.ID=="TSUHET" & traits$vpd.gy.max<1.2),], method='lm')

pTSUHET <- ggplot(traits[which(traits$SP.ID=="TSUHET"),], aes(x=log.LMA, y=log.LL, col=MAP)) + geom_point() +
  geom_smooth(data=traits[which(traits$SP.ID=="TSUHET" & traits$MAP<2200),], method='lm', col="black") +
  geom_smooth(data=traits[which(traits$SP.ID=="TSUHET" & traits$MAP>2200),], method='lm') + ggtitle("TSUHET")


# the same shifting up of LL for a given LMA at low Precip. Similar slope pattern to PSEMEN
# a quick model selection to verify
mod1 <- lm(log.LL~1, traits[which(traits$SP.ID=="TSUHET" & traits$MAP>0 & traits$log.LMA >0),]) # there are 5 NAs in MAP
mod2 <- lm(log.LL~log.LMA, traits[which(traits$SP.ID=="TSUHET"& traits$MAP>0 & traits$log.LMA >0),])
mod3 <- lm(log.LL~vpd.gy.max, traits[which(traits$SP.ID=="TSUHET"& traits$MAP>0 & traits$log.LMA >0),])
mod4 <- lm(log.LL~log.LMA + vpd.gy.max, traits[which(traits$SP.ID=="TSUHET"& traits$MAP>0 & traits$log.LMA >0),])
mod5 <- lm(log.LL~log.LMA * vpd.gy.max, traits[which(traits$SP.ID=="TSUHET"& traits$MAP>0 & traits$log.LMA >0),])
# but only the intercept shift is significant...
AIC(mod1, mod2, mod3, mod4, mod5) 
anova(mod1, mod3) # seems to be a very significant increase of LL at higher VPDs...
# again, not entirely sure how to interpret this...
xtabs(~SP.ID, traits.common)



# PINCON
pPINCON <- ggplot(traits[which(traits$SP.ID=="PINCON"),], aes(x=log.LMA, y=log.LL, col=MAP)) + geom_point() +
  geom_smooth(data=traits[which(traits$SP.ID=="PINCON" & traits$MAP<1000),], method='lm', col="black") +
  geom_smooth(data=traits[which(traits$SP.ID=="PINCON" & traits$MAP>1000),], method='lm') + ggtitle("PINCON")

pPINJEF <- ggplot(traits[which(traits$SP.ID=="PINJEF"),], aes(x=log.LMA, y=log.LL, col=MAP)) + geom_point() +
  geom_smooth(data=traits[which(traits$SP.ID=="PINJEF" & traits$MAP<1000),], method='lm', col="black") +
  geom_smooth(data=traits[which(traits$SP.ID=="PINJEF" & traits$MAP>1000),], method='lm') + ggtitle("PINJEF")

pPICSIT <- ggplot(traits[which(traits$SP.ID=="PICSIT"),], aes(x=log.LMA, y=log.LL, col=MAP)) + geom_point() +
  geom_smooth(data=traits[which(traits$SP.ID=="PICSIT" & traits$MAP<2500),], method='lm', col="black") +
  geom_smooth(data=traits[which(traits$SP.ID=="PICSIT" & traits$MAP>2500),], method='lm') + ggtitle("PICSIT")

pABIAMA <- ggplot(traits[which(traits$SP.ID=="ABIAMA"),], aes(x=log.LMA, y=log.LL, col=MAP)) + geom_point() +
  geom_smooth(data=traits[which(traits$SP.ID=="ABIAMA" & traits$MAP<2400),], method='lm', col="black") +
  geom_smooth(data=traits[which(traits$SP.ID=="ABIAMA" & traits$MAP>2400),], method='lm') + ggtitle("ABIAMA")

######## All well represented species together
multiplot(pPSME, pPINPON,pABIGRA,pABICON,pTSUHET,pPINCON, pPINJEF, pPICSIT,pABIAMA,cols=3)

# JUNOCC (no LL, just looking at lma)
Mypairs( traits[which(traits$SP.ID=="JUNOCC"), c("LMA","MAT", "MAP","vpd.gy.max","cmi.gy.mm","soilmoist.lvl1.mm")])
  # seems to be one anomolous stand.. very high P and low T, but somehow low soil moisture??\
  # 


########## Full mixed-effect analysis of LL ~f(LMA * MAP) effects ############
xtabs(~SP.ID, traits.common[which(traits.common$LEAF_LIFE>1),])
  # 13 of the 16 species in the full analysis are PNW
sp <- names(xtabs(~Species, LES[-which(is.na(spp.data$log.LL) | is.na(spp.data$log.LMA)),])[which(xtabs(~Species, LES[-which(is.na(spp.data$log.LL)| is.na(spp.data$log.LMA)),])>4)])
tmp <- LES[which(LES$Species %in% sp),]
# Acer rubrum              Acer saccharum          Crataegus monogyna 
# 8                           5                           5 
# Leucadendron salignum (m&f)         Populus tremuloides           Protea laurifolia 
# 8                           5                           5 
# Protea neriifolia               Protea repens               Quercus rubra 
# 5                           6                           8 
# Vaccinium myrtillus 
# 6 
### Nope. something's wrong here, because only Quercus chrysolepsis and thithocarpus densiflora make it through.
  # turns out none of them actually have good LL data. grrrr...


# so with traits.common5 we go!
fullmod <- lmer(log.LL~log.LMA * scale(MAP) + (log.LMA|SP.ID) + (0 + log.LMA:scale(MAP)|SP.ID), traits.common5)

######### What are the climate relationships??? ##########
Mypairs(traits.common[which(traits.common$SP.ID=="PSEMEN"), c("log.LL","log.LMA","log.Narea","climPC1","climPC2","soil_N","TotalSoilDepth","ASA",
                                                            "LAI_O","RGR","dominance","AG_TGROWTH")])
# LL is the only thing related to climate or stand environment
# PSMEN: LL related to Precip and not Temp! decreases with soil_N, increases with ASA, and decreases with RGR

  # one take-away:all the moisture variables are pretty well related, 
# and all except vpd are equally negatively related to LL for PSEMEN.
# This suggests that I can probably just stick with climPC1 & 2 for a Wetness and a Temp clim var.
Mypairs(traits.common[which(traits.common$SP.ID=="PINPON"), c("log.LL","log.LMA","log.Narea","climPC1","climPC2","soil_N","TotalSoilDepth","ASA",
                                                              "LAI_O","RGR","dominance")])
# PINPON: stand age and heat affect LL, but not much effect.
# - again, no environmental relationships with LMA or Narea

Mypairs(traits.common[which(traits.common$SP.ID=="ABIGRA" & traits.common$ASA<400), c("log.LL","log.LMA","log.Narea","climPC1","climPC2","soil_N","TotalSoilDepth","ASA",
                                                              "LAI_O","RGR","dominance")])

# ABIGRA: there's one old old stand. things look different with ASA if you pull that stand out.
# - LL related to T, soilN, soil depth??, and RGR. 
# - also first time there's a relationship between Narea and LL. 

Mypairs(traits.common[which(traits.common$SP.ID=="ABICON"), c("log.LL","log.LMA","log.Narea","climPC1","climPC2","soil_N","ASA",
                                                                                      "LAI_O","RGR","dominance")])

# ABICON: LL a little related to soil N, and increases with LAI_O. But not really a strong relationship with RGR or Narea
# -LMA and Narea actually positively related to growth rate...

Mypairs(traits.common[which(traits.common$SP.ID=="TSUHET"), c("log.LL","log.LMA","log.Narea","climPC1","climPC2","soil_N","TotalSoilDepth","ASA",
                                                                                      "LAI_O","RGR","dominance")])
# TSUHET : LL decreases with Temp, decreases with soil_N, increases with ASA, 
# Narea actually partly related to LL.
# -LMA strongly related to temp, soil_N, and stand age.

Mypairs(traits.common[which(traits.common$SP.ID=="PINCON"), c("log.LL","log.LMA","log.Narea","climPC1","climPC2","soil_N","ASA",
                                                                                      "LAI_O","RGR","dominance")])
# PINCON : LL related to temp, soil_N, stand age, and RGR.
# - lma not related to anything, Narea only related to soil_N

Mypairs(traits.common[which(traits.common$SP.ID=="PINJEF"), c("log.LL","log.LMA","log.Narea","climPC1","climPC2","soil_N","ASA",
                                                                                      "LAI_O","RGR","dominance")])
# PINJEF : LL strongly related to temp, DECREASE with stand age, and INCREASE with LAI_O. no relationship to RGR
# - lma and Narea actually negatively related to LL.

Mypairs(traits.common[which(traits.common$SP.ID=="ABIAMA"), c("log.LL","log.LMA","log.Narea","climPC1","climPC2","soil_N","TotalSoilDepth","ASA",
                                                                                      "LAI_O","RGR","dominance")])

# ABIAMA: LL kinda related to precip, NOT to temp, increases with stand age and decreases with LAI
#- lma and Narea actually related to LL! and ASA and LAI.



#### Modeling the trait-climate relationships 

require(MuMIn)
options(na.action = "na.fail")


trait.mods <- function(traitdata =traits.common, species, trait="log.LL", modcrit=F ){
  dataz <- data.frame(traitdata) %>% filter(SP.ID==species) %>% select(PLOT_ID, get(trait), climPC1, climPC2, soil_N, ASA, LAI_O, AG_TGROWTH) %>% filter(complete.cases(.))
   print(nrow(dataz))
  dataz$log.ASA <- log(dataz$ASA)
  datazsc <- scale(dataz[,-1])
  colnames(datazsc) <- paste0(colnames(dataz[,-1]),'sc')
  datazall <- data.frame(dataz, datazsc)
  traitmod <- lmer(get(trait)~climPC1sc + climPC2sc + soil_Nsc + log.ASAsc + LAI_Osc + AG_TGROWTHsc + (1|PLOT_ID), datazall, REML=F)
  traitdredge <- dredge(traitmod, extra=list(r.squaredGLMM))
  traitmodcalls <- dredge(traitmod, evaluate = F)
#  traitmodavg <- model.avg(traitdredge,subset=cumsum(weight)<=.90)
  traitmodavg <- model.avg(traitdredge,subset=delta<=4)
  
  bestmodobject <- eval(traitmodcalls[[rownames(traitdredge)[1]]])
  if(modcrit==T){
  scatter.smooth(resid(bestmodobject)~fitted(bestmodobject)); abline(h=0)
  plot(dataz[,trait]~predict(bestmodobject, re.form=NA));abline(a=0,b=1)
  qqp(resid(bestmodobject), main="residuals")
  qqp(ranef(bestmodobject)[[1]][,1], main='Random Effects')
  }
  results <- list(best = traitdredge[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc","r.squaredGLMM.R2m","r.squaredGLMM.R2c",'AICc')]
                  , avg = traitmodavg$coefficients[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc")]
                  , imp = as.vector(traitmodavg$importance)[match(c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc"),names(traitmodavg$importance))]
                  , fulldredge = traitdredge
                  , fullmodavg = traitmodavg
                  , modnum = rownames(traitdredge)[1]
                  , bestmodcall = traitmodcalls[[rownames(traitdredge)[1]]]
                  , bestmodobject = bestmodobject
                  , n=nrow(datazall)
                  , deltaNULL = traitdredge[which(rownames(traitdredge)==1),"delta"])
  return(results)
}

#PSEMENllraw <- trait.mods(traitdata = traits.common, species = "PSEMEN",trait = "LLmonths", modcrit = T)
  # raw values almost identical, variable rankings very close
PSEMENll <- trait.mods(traitdata = traits.common, species = "PSEMEN",trait = "log.LL", modcrit=T)
# PSEMENlmaraw <- trait.mods(traitdata = traits.common, species = "PSEMEN",trait = "LMA", modcrit=T)
  # raw values less predictable, and resids less normal, variable rankings v similar
PSEMENlma <- trait.mods(traitdata = traits.common, species = "PSEMEN",trait = "log.LMA", modcrit=T)
# PSEMENnarearaw <- trait.mods(traitdata = traits.common, species = "PSEMEN",trait = "Narea", modcrit=T)
  # raw values pretty similar R2, variable rankings nearly identical
PSEMENnarea <- trait.mods(traitdata = traits.common, species = "PSEMEN",trait = "log.Narea", modcrit=T)
  # n=220

PINPONll <- trait.mods(traitdata = traits.common, species = "PINPON",trait = "log.LL")
PINPONlma <- trait.mods(traitdata = traits.common, species = "PINPON",trait = "log.LMA")
PINPONnarea <- trait.mods(traitdata = traits.common, species = "PINPON",trait = "log.Narea")
  # n= 97

ABICONll <- trait.mods(traitdata = traits.common, species = "ABICON",trait = "log.LL")
ABICONlma <- trait.mods(traitdata = traits.common, species = "ABICON",trait = "log.LMA")
ABICONnarea <- trait.mods(traitdata = traits.common, species = "ABICON",trait = "log.Narea")

TSUHETll <- trait.mods(traitdata = traits.common, species = "TSUHET",trait = "log.LL")
TSUHETlma <- trait.mods(traitdata = traits.common, species = "TSUHET",trait = "log.LMA")
TSUHETnarea <- trait.mods(traitdata = traits.common, species = "TSUHET",trait = "log.Narea") # might be one outlier that's leveraging things?

PINCONll <- trait.mods(traitdata = traits.common, species = "PINCON",trait = "log.LL")
PINCONlma <- trait.mods(traitdata = traits.common, species = "PINCON",trait = "log.LMA")
PINCONnarea <- trait.mods(traitdata = traits.common, species = "PINCON",trait = "log.Narea")

PINJEFll <- trait.mods(traitdata = traits.common, species = "PINJEF",trait = "log.LL")
PINJEFlma <- trait.mods(traitdata = traits.common, species = "PINJEF",trait = "log.LMA") # two outliers
PINJEFnarea <- trait.mods(traitdata = traits.common, species = "PINJEF",trait = "log.Narea")


####### Combining all environmental trait models ######
llbestmods <- rbind(PSEMENll$best, PINPONll$best,PINCONll$best,PINJEFll$best, ABICONll$best, TSUHETll$best)
llbestmods$n <- rbind(PSEMENll$n, PINPONll$n,PINCONll$n,PINJEFll$n, ABICONll$n, TSUHETll$n)
llbestmods$deltas <-rbind(PSEMENll$deltaNULL, PINPONll$deltaNULL,PINCONll$deltaNULL,PINJEFll$deltaNULL, ABICONll$deltaNULL, TSUHETll$deltaNULL)
llbestmods$SP.ID <-  c("PSEMEN","PINPON","PINCON","PINJEF","ABICON","TSUHET")

lmabestmods <- rbind(PSEMENlma$best, PINPONlma$best,PINCONlma$best,PINJEFlma$best, ABICONlma$best, TSUHETlma$best)
lmabestmods$n <- rbind(PSEMENlma$n, PINPONlma$n,PINCONlma$n,PINJEFlma$n, ABICONlma$n, TSUHETlma$n)
lmabestmods$deltas <-rbind(PSEMENlma$deltaNULL, PINPONlma$deltaNULL,PINCONlma$deltaNULL,PINJEFlma$deltaNULL, ABICONlma$deltaNULL, TSUHETlma$deltaNULL)
lmabestmods$SP.ID <-  c("PSEMEN","PINPON","PINCON","PINJEF","ABICON","TSUHET")


nareabestmods <- rbind(PSEMENnarea$best, PINPONnarea$best,PINCONnarea$best,PINJEFnarea$best, ABICONnarea$best, TSUHETnarea$best)
nareabestmods$n <- rbind(PSEMENnarea$n, PINPONnarea$n,PINCONnarea$n,PINJEFnarea$n, ABICONnarea$n, TSUHETnarea$n)
nareabestmods$deltas <-rbind(PSEMENnarea$deltaNULL, PINPONnarea$deltaNULL,PINCONnarea$deltaNULL,PINJEFnarea$deltaNULL, ABICONnarea$deltaNULL, TSUHETnarea$deltaNULL)
nareabestmods$SP.ID <-  c("PSEMEN","PINPON","PINCON","PINJEF","ABICON","TSUHET")
nareabestmods$Species <-  c("Pseudotsuga menziesii","Pinus ponderosa","Pinus contorta","Pinus jeffreyii","Abies concolor","Tsuga heterophylla")

write.csv(llbestmods, "Trait_models/LeafLife_bestmodels_031617.csv")
write.csv(lmabestmods, "Trait_models/LMA_bestmodels_031617.csv")
write.csv(nareabestmods, "Trait_models/Narea_bestmodels_031617.csv")


#### Variable Importances ######
llimps <- data.frame(rbind(PSEMENll$imp, PINPONll$imp,PINCONll$imp,PINJEFll$imp, ABICONll$imp, TSUHETll$imp)[,-1])
colnames(llimps) <- c("climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc")
llimps$SP.ID <- c("PSEMEN","PINPON","PINCON","PINJEF","ABICON","TSUHET")
llimps$trait <- rep("LeafLife", times=nrow(llimps))

lmaimps <- data.frame(rbind(PSEMENlma$imp, PINPONlma$imp,PINCONlma$imp,PINJEFlma$imp, ABICONlma$imp, TSUHETlma$imp)[,-1])
colnames(lmaimps) <- c("climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc")
lmaimps$SP.ID <- c("PSEMEN","PINPON","PINCON","PINJEF","ABICON","TSUHET")
lmaimps$trait <- rep("LMA", times=nrow(lmaimps))

nareaimps <- data.frame(rbind(PSEMENnarea$imp, PINPONnarea$imp,PINCONnarea$imp,PINJEFnarea$imp, ABICONnarea$imp, TSUHETnarea$imp)[,-1])
colnames(nareaimps) <- c("climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc")
nareaimps$SP.ID <- c("PSEMEN","PINPON","PINCON","PINJEF","ABICON","TSUHET")
nareaimps$trait <- rep("Narea", times=nrow(nareaimps))

importances <- rbind(llimps, lmaimps, nareaimps)
colnames(importances) <- c("climPC1", "climPC2","soil_N","Stand_Age","LAI","Growth", "SP.ID", "trait")
# write.csv(importances, "Trait_models/VariableImportances_wide_041617aiccuttoff.csv")


importances <- read.csv("Trait_models/VariableImportances_wide_041617aiccuttoff.csv", row.names = 1)
impslong <- melt(importances, id.vars = c("SP.ID","trait"))

p1 <- ggplot(impslong, aes(x=variable, y=value, col=variable)) + geom_boxplot() + 
  geom_abline(slope = 0, intercept = 0) + facet_grid(trait~.) +
  theme(axis.text.x=element_text(angle=90), legend.position="none") + ylab("Variable Importnce")

##### Model Averaged Effects Sizes ######
llavgmods <- data.frame(rbind(PSEMENll$avg, PINPONll$avg,PINCONll$avg,PINJEFll$avg, ABICONll$avg, TSUHETll$avg))
lmaavgmods <- data.frame(rbind(PSEMENlma$avg, PINPONlma$avg,PINCONlma$avg,PINJEFlma$avg, ABICONlma$avg, TSUHETlma$avg))
nareaavgmods <- data.frame(rbind(PSEMENnarea$avg, PINPONnarea$avg,PINCONnarea$avg,PINJEFnarea$avg, ABICONnarea$avg, TSUHETnarea$avg))
llavgmods$SP.ID <- c("PSEMEN","PINPON","PINCON","PINJEF","ABICON","TSUHET")
llavgmods$trait <- rep("LeafLife", times=nrow(llavgmods))
lmaavgmods$SP.ID <- c("PSEMEN","PINPON","PINCON","PINJEF","ABICON","TSUHET")
lmaavgmods$trait <- rep("LMA", times=nrow(lmaavgmods))
nareaavgmods$SP.ID <- c("PSEMEN","PINPON","PINCON","PINJEF","ABICON","TSUHET")
nareaavgmods$trait <- rep("Narea", times=nrow(nareaavgmods))

avgmods <- rbind(llavgmods, lmaavgmods, nareaavgmods)[,-1]
colnames(avgmods) <- c("climPC1", "climPC2","soil_N","Stand_Age","LAI","Growth", "SP.ID", "trait")
#write.csv(avgmods, "Trait_models/Model_Averages_wide_031617aiccuttoff.csv")

avgmods <- read.csv("Trait_models/Model_Averages_wide_031617aiccuttoff.csv", row.names = 1)
avglong <- melt(avgmods, id.vars = c("SP.ID","trait"))

(p2 <- ggplot(avglong, aes(x=variable, y=value, col=variable)) + geom_boxplot() + 
  geom_abline(slope = 0, intercept = 0) + facet_grid(trait~.) +
  theme(axis.text.x=element_text(angle=90), legend.position="None") + ylab("Model-averaged Coefficients"))


quartz(width=6, height=5)
multiplot(p2,p1, cols = 2)

### PICSIT and ABIAMA are right on the edge of useful ###
  #### yeeee, not sure I trust these.
PICSITll <- trait.mods(traitdata = traits.common, species = "PICSIT",trait = "log.LL")
PICSITlma <- trait.mods(traitdata = traits.common, species = "PICSIT",trait = "log.LMA")
PICSITnarea <- trait.mods(traitdata = traits.common, species = "PICSIT",trait = "log.Narea")


ABIAMAll <- trait.mods(traitdata = traits.common, species = "ABIAMA",trait = "log.LL")
ABIAMAlma <- trait.mods(traitdata = traits.common, species = "ABIAMA",trait = "log.LMA")
ABIAMAnarea <- trait.mods(traitdata = traits.common, species = "ABIAMA",trait = "log.Narea")



ABIGRAll <- trait.mods(traitdata = traits.common, species = "ABIGRA",trait = "log.LL")
ABIGRAlma <- trait.mods(traitdata = traits.common, species = "ABIGRA",trait = "log.LMA")
ABIGRAnarea <- trait.mods(traitdata = traits.common, species = "ABIGRA",trait = "log.Narea")
# n= 12... and the main troubles are soil_N, ASA, log.Narea
# PLOT_ID     log.LL    climPC1    climPC2     soil_N        ASA      LAI_O AG_TGROWTH 
#     0          5          0          0         51         34          8          8 

# JUNOCCll <- trait.mods(traitdata = traits.common, species = "JUNOCC",trait = "log.LL")
JUNOCClma <- trait.mods(traitdata = traits.common, species = "JUNOCC",trait = "log.LMA")
JUNOCCnarea <- trait.mods(traitdata = traits.common, species = "JUNOCC",trait = "log.Narea")
# damn. missing lots of rows. NAs:
# PLOT_ID     log.LL    climPC1    climPC2     soil_N        ASA      LAI_O AG_TGROWTH 
# 0         68          0          0         55         10          0          0 


# ll <- lmer(log.LL~climPC1sc + climPC2sc + soil_Nsc + log.ASAsc + LAI_Osc + AG_TGROWTHsc + (1|PLOT_ID), datazall, REML=F)
# # lllm <- lm(log.LL~climPC1sc + climPC2sc + soil_Nsc + scale(log(ASA)) + LAI_Osc + AG_TGROWTHsc , datazall)
# # AIC(ll1,ll2,ll3,ll4)
# lldredge <- dredge(ll, extra=list(r.squaredGLMM))
# llmods <- dredge(ll, evaluate = F)
# llmodavg <- model.avg(lldredge,subset=cumsum(weight)<=.90)
# lma <- lmer(log.LMA~climPC1 + climPC2 + soil_N + scale(log(ASA)) + LAI_O + AG_TGROWTH + (1|PLOT_ID), dataz, REML=F)
# lmadredge <- dredge(lma,extra=list(r.squaredGLMM) )
# lmamodavg <- model.avg(lmadredge,subset=cumsum(weight)<=.90)
# Narea1 <- lmer(log.Narea~climPC1 + climPC2 + soil_N + scale(log(ASA)) + LAI_O + AG_TGROWTH + (1|PLOT_ID), dataz, REML=F)
# Nareadredge <- dredge(Narea1, extra=list(r.squaredGLMM))
# Nareamodavg <- model.avg(Nareadredge, subset=cumsum(weight)<=.90)








#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------

############ ***Community-weighted trait calculations *** ####################

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# note on nomenclature:
# TRAIT1 = species mean trait value for species 1
# wTRAIT1 = BA weighted species mean trait value
# TRAITs1 = site mean value for the species species 1
# wTRAITs1 = BA weighted site mean value for species 1



# total species means
species.means <- traits %>% group_by(GENUS,SPECIES,SP.ID) %>% summarise(mSLA = mean(SLA_HSA, na.rm=T), mLEAF_LIFE = mean(LEAF_LIFE, na.rm=T), mCARBON = mean(LEAF_CARBON, na.rm=T), mNITROGEN = mean(LEAF_NITROGEN, na.rm=T), mCN = mean(LEAF_CN, na.rm=T), N = n(), nPlots = length(unique(PLOT_ID)), nProj = length(unique(PROJECT)), nEcoReg = length(unique(ECOREGION)), mMAP = mean(MAP), mMAT = mean(MAT), mlog.LMA = mean(log.LMA, na.rm=T), mlog.LL =mean(log.LL, na.rm=T), mlog.Nmass=mean(log.Nmass, na.rm=T), mlog.LMA_PSA = mean(log.LMA_PSA, na.rm=T), mRGR = mean(RGRdom, na.rm=T), mstGrowthdom = mean(stGrowthdom, na.rm=T)) 
## site means
pairs(species.means[,c(1,4,5,6,7,8,9,13,14)], upper.panel = panel.smooth)

common.genera <- c("Tsuga","Picea","Acer","Quercus","Pinus","Abies")
plot(mSLA~mMAT, species.means[which(species.means$GENUS %in% common.genera),], col=GENUS, pch=16, ylim=c(0,200))


####### Plot Trait values based on species mean values ######
biomass$SLA1 <- species.means$mSLA[match(biomass$SPP_O1_ABBREV, species.means$SP.ID)]
biomass$SLA2 <- species.means$mSLA[match(biomass$SPP_O2_ABBREV, species.means$SP.ID)]
biomass$SLA3 <- species.means$mSLA[match(biomass$SPP_O3_ABBREV, species.means$SP.ID)]
biomass$SLA4 <- species.means$mSLA[match(biomass$SPP_O4_ABBREV, species.means$SP.ID)]

biomass$wSLA1 <- biomass$SLA1 * biomass$SPP_O1_BASAL_AREA_FRACTION/100
biomass$wSLA2 <- biomass$SLA2 * biomass$SPP_O2_BASAL_AREA_FRACTION/100
biomass$wSLA3 <- biomass$SLA3 * biomass$SPP_O3_BASAL_AREA_FRACTION/100
biomass$wSLA4 <- biomass$SLA4 * biomass$SPP_O4_BASAL_AREA_FRACTION/100
biomass$ComM_SLA <- apply(biomass[,c("wSLA1", "wSLA2","wSLA3","wSLA4")],MARGIN=1,FUN=sum, na.rm=T)


plot(ComM_SLA~MAP, biomass, col=SPP_O1_ABBREV, pch=14 + as.numeric(biomass$PROJECT))
plot(ComM_SLA~MAT_C, biomass, col=SPP_O1_ABBREV, pch=14 + as.numeric(biomass$PROJECT))
abline(lm(ComM_SLA~MAT_C, biomass))
summary(lm(ComM_SLA~MAP + MAT_C + SPP_O1_ABBREV, biomass))

## let's look at how many plots are dominated primarily by one species:
length(biomass[which(biomass$SPP_O1_BASAL_AREA_FRACTION==100),1])
  #55 plots are monoculture
length(biomass[which(biomass$SPP_O1_BASAL_AREA_FRACTION>90),1])
  #116/265 plots are >90% 1 spp
length(biomass[which(biomass$SPP_O1_BASAL_AREA_FRACTION>75),1])
  #150 plots are >75% 1 species.
length(biomass[which(biomass$SPP_O1_BASAL_AREA_FRACTION<50),1])
  #only 39 of 265 plots are bonafide co-dominant...




####### **Comparing Plot Level Values based on SPP averages, vs Plot values*** ######
spp.traits <- traits %>% group_by(GE.SP) %>% summarise(nsample = n(), SLA = mean(SLA_HSA, na.rm=T), nSLA = n()- length(which(is.na(SLA_HSA))), CN = mean(LEAF_CN, na.rm=T), nCN = n()- length(which(is.na(LEAF_CN))), LIFE = mean(LEAF_LIFE, na.rm=T), nLIFE=n()- length(which(is.na(LEAF_LIFE))), nplots =length(unique(PLOT_ID)) )
## An idea: leaf lifespan seems to vary quite a bit with MAT. What if I calculated plot averaged leaf lifespan based on spp averages vs based on the actual data from the plot? could say something about using a single value for a species?
spp.plot.traits <- traits %>% group_by(SP.ID, PLOT_ID) %>% summarise(nsample = n(), mSLA = mean(SLA_HSA, na.rm=T), nSLA = n()- length(which(is.na(SLA_HSA))), mCN = mean(LEAF_CN, na.rm=T), nCN = n()- length(which(is.na(LEAF_CN)))
                                                                     , mLIFE = mean(LEAF_LIFE, na.rm=T), nLIFE=n()- length(which(is.na(LEAF_LIFE))), mNITROGEN = mean(LEAF_NITROGEN, na.rm=T), nplots =length(unique(PLOT_ID))
                                                                     , mLMA = mean(LMA_PSA, na.rm=T), mLLmonths= mean(LLmonths, na.rm=T), mNarea=mean(Narea, na.rm=T)
                                                                     , mlog.LMA = mean(log.LMA, na.rm=T), mlog.LL = mean(log.LL, na.rm=T), mlog.Narea=mean(log.Narea, na.rm=T), mlog.Nmass=mean(log.Nmass, na.rm=T), RGR = unique(RGRdom), stGrowth = mean(stGrowthdom, na.rm=T))
  # Note: I used LMA_PSA in other trait analysis in TaxonomicAnalysis. So I've switched to it here
spp.plot.traits$SP.PLOT <- paste(spp.plot.traits$SP.ID, spp.plot.traits$PLOT_ID, sep="-")
# let's just pull out the common dominant species in an easy to name df
plotavs <- spp.plot.traits%>% subset(spp.plot.traits$SP.ID %in% names(which(xtabs(~SPP_O1_ABBREV, biomass)>2)))

### make unique identifiers for matching biomass and spp.plot.traits rows
biomass$SP1.PLOT <- paste(biomass$SPP_O1_ABBREV,biomass$PLOT_ID, sep="-")
biomass$SP2.PLOT <- paste(biomass$SPP_O2_ABBREV,biomass$PLOT_ID, sep="-")
biomass$SP3.PLOT <- paste(biomass$SPP_O3_ABBREV,biomass$PLOT_ID, sep="-")
biomass$SP4.PLOT <- paste(biomass$SPP_O4_ABBREV,biomass$PLOT_ID, sep="-")


#### Plot average trait values for each species ############################
biomass$SLAs1 <- spp.plot.traits$mSLA[match(as.character(biomass$SP1.PLOT), as.character(spp.plot.traits$SP.PLOT))]
biomass$SLAs2 <- spp.plot.traits$mSLA[match(biomass$SP2.PLOT, spp.plot.traits$SP.PLOT)]
biomass$SLAs3 <- spp.plot.traits$mSLA[match(biomass$SP3.PLOT, spp.plot.traits$SP.PLOT)]
biomass$SLAs4 <- spp.plot.traits$mSLA[match(biomass$SP4.PLOT, spp.plot.traits$SP.PLOT)]

biomass$wSLAs1 <- biomass$SLAs1 * biomass$SPP_O1_BASAL_AREA_FRACTION/100
biomass$wSLAs2 <- biomass$SLAs2 * biomass$SPP_O2_BASAL_AREA_FRACTION/100
biomass$wSLAs3 <- biomass$SLAs3 * biomass$SPP_O3_BASAL_AREA_FRACTION/100
biomass$wSLAs4 <- biomass$SLAs4 * biomass$SPP_O4_BASAL_AREA_FRACTION/100
biomass$ComM_SLAs <- apply(biomass[,c("wSLAs1", "wSLAs2","wSLAs3","wSLAs4")],MARGIN=1,FUN=sum, na.rm=T)
# now need to remove values that got smoothed over due to na.rm=T in the apply function
biomass$ComM_SLAs[which(is.na(biomass$SLAs1))] <- NA # 85 plots lacking SLA data of SPP01
biomass$ComM_SLAs[which(is.na(biomass$SLAs2) & biomass$SPP_O2_BASAL_AREA_FRACTION>30)] <- NA
  #18 sites lacking SLA data for SPP02 where SPP02 makes up >30 % of BA, but only 5 have SPP01
biomass$ComM_SLAs[which(biomass$ComM_SLAs==0)] <- NA


#### Community weighted LEAf_LIFE based on plot values
biomass$LIFEs1 <- spp.plot.traits$mLIFE[match(biomass$SP1.PLOT, spp.plot.traits$SP.PLOT)]
biomass$LIFEs2 <- spp.plot.traits$mLIFE[match(biomass$SP2.PLOT, spp.plot.traits$SP.PLOT)]
biomass$LIFEs3 <- spp.plot.traits$mLIFE[match(biomass$SP3.PLOT, spp.plot.traits$SP.PLOT)]
biomass$LIFEs4 <- spp.plot.traits$mLIFE[match(biomass$SP4.PLOT, spp.plot.traits$SP.PLOT)]

biomass$wLIFEs1 <- biomass$LIFEs1 * biomass$SPP_O1_BASAL_AREA_FRACTION/100
biomass$wLIFEs2 <- biomass$LIFEs2 * biomass$SPP_O2_BASAL_AREA_FRACTION/100
biomass$wLIFEs3 <- biomass$LIFEs3 * biomass$SPP_O3_BASAL_AREA_FRACTION/100
biomass$wLIFEs4 <- biomass$LIFEs4 * biomass$SPP_O4_BASAL_AREA_FRACTION/100
biomass$ComM_LIFEs <- apply(biomass[,c("wLIFEs1", "wLIFEs2","wLIFEs3","wLIFEs4")],MARGIN=1,FUN=sum, na.rm=T)
# now need to remove values that got smoothed over due to na.rm=T in the apply function
biomass$ComM_LIFEs[which(is.na(biomass$LIFEs1))] <- NA # 94 sites lacking LeafLife of SPP01
biomass$ComM_LIFEs[which(is.na(biomass$LIFEs2) & biomass$SPP_O2_BASAL_AREA_FRACTION>30)] <- NA
  # 22 sites lack SPP02 leaf life, 9 of them unique


#### Community weighted LEAf_NITROGEN based on plot values
biomass$NITROGENs1 <- spp.plot.traits$mNITROGEN[match(biomass$SP1.PLOT, spp.plot.traits$SP.PLOT)]
biomass$NITROGENs2 <- spp.plot.traits$mNITROGEN[match(biomass$SP2.PLOT, spp.plot.traits$SP.PLOT)]
biomass$NITROGENs3 <- spp.plot.traits$mNITROGEN[match(biomass$SP3.PLOT, spp.plot.traits$SP.PLOT)]
biomass$NITROGENs4 <- spp.plot.traits$mNITROGEN[match(biomass$SP4.PLOT, spp.plot.traits$SP.PLOT)]

biomass$wNITROGENs1 <- biomass$NITROGENs1 * biomass$SPP_O1_BASAL_AREA_FRACTION/100
biomass$wNITROGENs2 <- biomass$NITROGENs2 * biomass$SPP_O2_BASAL_AREA_FRACTION/100
biomass$wNITROGENs3 <- biomass$NITROGENs3 * biomass$SPP_O3_BASAL_AREA_FRACTION/100
biomass$wNITROGENs4 <- biomass$NITROGENs4 * biomass$SPP_O4_BASAL_AREA_FRACTION/100
biomass$ComM_NITROGENs <- apply(biomass[,c("wNITROGENs1", "wNITROGENs2","wNITROGENs3","wNITROGENs4")],MARGIN=1,FUN=sum, na.rm=T)
# now need to remove values that got smoothed over due to na.rm=T in the apply function
biomass$ComM_NITROGENs[which(is.na(biomass$NITROGENs1))] <- NA # 98 sites lacking LeafNITROGEN of SPP01
biomass$ComM_NITROGENs[which(is.na(biomass$NITROGENs2) & biomass$SPP_O2_BASAL_AREA_FRACTION>30)] <- NA
biomass$ComM_NITROGENs[which(biomass$ComM_NITROGENs==0)] <- NA
# 26 sites lack SPP02 leaf NITROGEN, 3 of them unique





#### Community weighted Narea based on plot values
biomass$Nareas1 <- spp.plot.traits$mNarea[match(biomass$SP1.PLOT, spp.plot.traits$SP.PLOT)]
biomass$Nareas2 <- spp.plot.traits$mNarea[match(biomass$SP2.PLOT, spp.plot.traits$SP.PLOT)]
biomass$Nareas3 <- spp.plot.traits$mNarea[match(biomass$SP3.PLOT, spp.plot.traits$SP.PLOT)]
biomass$Nareas4 <- spp.plot.traits$mNarea[match(biomass$SP4.PLOT, spp.plot.traits$SP.PLOT)]

biomass$wNareas1 <- biomass$Nareas1 * biomass$SPP_O1_BASAL_AREA_FRACTION/100
biomass$wNareas2 <- biomass$Nareas2 * biomass$SPP_O2_BASAL_AREA_FRACTION/100
biomass$wNareas3 <- biomass$Nareas3 * biomass$SPP_O3_BASAL_AREA_FRACTION/100
biomass$wNareas4 <- biomass$Nareas4 * biomass$SPP_O4_BASAL_AREA_FRACTION/100
biomass$ComM_Nareas <- apply(biomass[,c("wNareas1", "wNareas2","wNareas3","wNareas4")],MARGIN=1,FUN=sum, na.rm=T)
# now need to remove values that got smoothed over due to na.rm=T in the apply function
biomass$ComM_Nareas[which(is.na(biomass$Nareas1))] <- NA # 98 sites lacking LeafNarea of SPP01
biomass$ComM_Nareas[which(is.na(biomass$Nareas2) & biomass$SPP_O2_BASAL_AREA_FRACTION>30)] <- NA
biomass$ComM_Nareas[which(biomass$ComM_Nareas==0)] <- NA
# 26 sites lack SPP02 leaf Narea, 3 of them unique







### making a 'canopy complexity' type metric #######
biomass$Ccomplex <- biomass$AG_BIOMASS_TREE_WOOD_AS_CARBON/biomass$LAI_O


#-----------------------
###### Infilling spp trait values with spp means when they're missing
#------------------------


### Infilling SLA with species mean values
wSLAif1 <- biomass$wSLAs1
wSLAif1[which(is.na(wSLAif1))] <- biomass$wSLA1[which(is.na(wSLAif1))]
wSLAif2 <- biomass$wSLAs2
wSLAif2[which(is.na(wSLAif2))] <- biomass$wSLA2[which(is.na(wSLAif2))]
wSLAif3 <- biomass$wSLAs3
wSLAif3[which(is.na(wSLAif3))] <- biomass$wSLA3[which(is.na(wSLAif3))]
wSLAif4 <- biomass$wSLAs4
wSLAif4[which(is.na(wSLAif4))] <- biomass$wSLA4[which(is.na(wSLAif4))]
ComM_SLAs_if  <- apply(data.frame(wSLAif1, wSLAif2, wSLAif3, wSLAif4),MARGIN=1, FUN=sum, na.rm=T)
ComM_SLAs_if[which(ComM_SLAs_if==0)] <- NA
# infilled 85 plots

##### Infilling LEAF_LIFE with species mean values ######
# note, need to make species average values first.
LIFE1 <- species.means$mLEAF_LIFE[match(biomass$SPP_O1_ABBREV, species.means$SP.ID)]
LIFE2 <- species.means$mLEAF_LIFE[match(biomass$SPP_O2_ABBREV, species.means$SP.ID)]
LIFE3 <- species.means$mLEAF_LIFE[match(biomass$SPP_O3_ABBREV, species.means$SP.ID)]
LIFE4 <- species.means$mLEAF_LIFE[match(biomass$SPP_O4_ABBREV, species.means$SP.ID)]
wLIFE1 <- LIFE1 * biomass$SPP_O1_BASAL_AREA_FRACTION/100 # 13 total plots either lack LIFE for the dominant spp or don't have SPP data
wLIFE2 <- LIFE2 * biomass$SPP_O2_BASAL_AREA_FRACTION/100
wLIFE3 <- LIFE3 * biomass$SPP_O3_BASAL_AREA_FRACTION/100
wLIFE4 <- LIFE4 * biomass$SPP_O4_BASAL_AREA_FRACTION/100

### now infill missing values with species means
wLIFEif1 <- biomass$wLIFEs1
wLIFEif1[which(is.na(wLIFEif1))] <- wLIFE1[which(is.na(wLIFEif1))]
wLIFEif2 <- biomass$wLIFEs2
wLIFEif2[which(is.na(wLIFEif2))] <- wLIFE2[which(is.na(wLIFEif2))]
wLIFEif3 <- biomass$wLIFEs3
wLIFEif3[which(is.na(wLIFEif3))] <- wLIFE3[which(is.na(wLIFEif3))]
wLIFEif4 <- biomass$wLIFEs4
wLIFEif4[which(is.na(wLIFEif4))] <- wLIFE4[which(is.na(wLIFEif4))]
ComM_LIFEs_if  <- apply(data.frame(wLIFEif1, wLIFEif2, wLIFEif3, wLIFEif4),MARGIN=1, FUN=sum, na.rm=T)
ComM_LIFEs_if[which(is.na(wLIFEif1))] <- NA
which(is.na(biomass$ComM_LIFEs) & ComM_LIFEs_if>0)
# infilled 90 plots

##### Infilling LEAF_NITROGEN with species mean values ######
# note, need to make species average values first.

NITROGEN1 <- species.means$mNITROGEN[match(biomass$SPP_O1_ABBREV, species.means$SP.ID)]
NITROGEN2 <- species.means$mNITROGEN[match(biomass$SPP_O2_ABBREV, species.means$SP.ID)]
NITROGEN3 <- species.means$mNITROGEN[match(biomass$SPP_O3_ABBREV, species.means$SP.ID)]
NITROGEN4 <- species.means$mNITROGEN[match(biomass$SPP_O4_ABBREV, species.means$SP.ID)]
wNITROGEN1 <- NITROGEN1 * biomass$SPP_O1_BASAL_AREA_FRACTION/100 # 13 total plots either lack NITROGEN for the dominant spp or don't have SPP data
wNITROGEN2 <- NITROGEN2 * biomass$SPP_O2_BASAL_AREA_FRACTION/100
wNITROGEN3 <- NITROGEN3 * biomass$SPP_O3_BASAL_AREA_FRACTION/100
wNITROGEN4 <- NITROGEN4 * biomass$SPP_O4_BASAL_AREA_FRACTION/100

### now infill missing values with species means
wNITROGENif1 <- biomass$wNITROGENs1
wNITROGENif1[which(is.na(wNITROGENif1))] <- wNITROGEN1[which(is.na(wNITROGENif1))]
wNITROGENif2 <- biomass$wNITROGENs2
wNITROGENif2[which(is.na(wNITROGENif2))] <- wNITROGEN2[which(is.na(wNITROGENif2))]
wNITROGENif3 <- biomass$wNITROGENs3
wNITROGENif3[which(is.na(wNITROGENif3))] <- wNITROGEN3[which(is.na(wNITROGENif3))]
wNITROGENif4 <- biomass$wNITROGENs4
wNITROGENif4[which(is.na(wNITROGENif4))] <- wNITROGEN4[which(is.na(wNITROGENif4))]
ComM_NITROGENs_if  <- apply(data.frame(wNITROGENif1, wNITROGENif2, wNITROGENif3, wNITROGENif4),MARGIN=1, FUN=sum, na.rm=T)
ComM_NITROGENs_if[which(is.na(wNITROGENif1))] <- NA
length(which(is.na(biomass$ComM_NITROGENs) & ComM_NITROGENs_if>0))
# infilled 96 plots

# plot(AG_PROD_TREE_TOTAL_AS_CARBON~ComM_NITROGENs_if, biomass, pch=16, cex=1.1)
# points(AG_PROD_TREE_TOTAL_AS_CARBON~ComM_NITROGENs, biomass, pch=16, cex=.9, col="red")

length(biomass$PLOT_ID[which(biomass$LAI_O>0 & biomass$AG_PROD_TREE_TOTAL_AS_CARBON>0 & ComM_NITROGENs_if>0 & ComM_LIFEs_if>0 & ComM_SLAs_if>0)]) 
length(biomass$PLOT_ID[which(biomass$LAI_O>0 & biomass$AG_PROD_TREE_TOTAL_AS_CARBON>0 & biomass$ComM_NITROGENs>0 & biomass$ComM_LIFEs>0 & biomass$ComM_SLAs>0)]) 

##### Export best dataset ## NOTE: THIS IS HARD CODED AND MIGHT NOT WORK IF COLUMNS ARE CHANGED
  # combine things together and remove the columns we don't actually need
biomassbest <- data.frame(biomass[,-c(35:43,48:55,57:64,66:73)],ComM_SLAs_if,ComM_LIFEs_if,ComM_NITROGENs_if)
colnames(biomassbest)[grep("ComM", colnames(biomassbest))] <- c("plotSLA", "plotLIFE","plotNITROGEN","plotSLA_if","plotLIFE_if","plotNITROGEN_if")
write.csv(biomassbest, "PACNW_Biomass_plus_traits_120716.csv")

 
biomass <- read.csv("PACNW_Biomass_plus_traits_120716.csv")

###### Plotting Community Weighted Means ##########
length(biomass$plotLIFE[which(biomass$plotLIFE>0 & biomass$plotSLA>0 & biomass$plotNITROGEN>0)])
plot(plotLIFE~plotSLA, biomass, col=SPP_O1_ABBREV) # wow. not much covariation between SLA and leaf Life at plot scale
plot(plotLIFE~plotNITROGEN, biomass, col=SPP_O1_ABBREV) # also no variation between leafLife and Nitrogen at the plot scale
plot(plotNITROGEN~plotSLA, biomass, col=SPP_O1_ABBREV) # positive relationship between leafN and SLA at plot level, not all of it is cover type driven

plot(Ccomplex~ASA, biomass, col=SPP_O1_ABBREV, pch=16, ylim=c(0,20000))

ggplot(biomass[which(biomass$SPP_O1_ABBREV %in% c("ABICON", "ABIGRA","JUNOCC", "LAROCC","PINCON","PINPON","PSEMEN")),], aes(x=ASA, y=Ccomplex, col=SPP_O1_ABBREV)) + geom_point()

#### community weighted Traits with LAI
plot(plotSLA~LAI_O, biomass)
plot(plotLIFE~LAI_O, biomass)
plot(plotNITROGEN~LAI_O, biomass)
  # 
plot(plotSLA~AG_PROD_TREE_TOTAL_AS_CARBON, biomass)
plot(plotLIFE~AG_PROD_TREE_TOTAL_AS_CARBON, biomass)
plot(plotNITROGEN~AG_PROD_TREE_TOTAL_AS_CARBON, biomass)


### figuring out how the names correspond between the abbreviations and the traits.
    ### Used this to great SP.ID column for cross referencing....
bio.sp <- levels(as.factor(c(levels(biomass$SPP_O1_ABBREV), levels(biomass$SPP_O2_ABBREV), levels(biomass$SPP_O3_ABBREV), levels(biomass$SPP_O4_ABBREV))))
trait.sp <- levels(as.factor(traits$GE.SP))
sp.names.raw <- cbind(c(bio.sp,NA,NA,NA,NA),trait.sp)
write.csv(sp.names.raw, "SpeciesNames_lookup_raw_031316.csv")





###### 







#-----------------------------------------------------------------------------
##### Traits available in TRY #########
#-----------------------------------------------------------------------------

SLA_TRY <- read.table("/Users/leeanderegg/Dropbox/Trade-offs project/TRY_SLA_SpeciesTable_03222016.txt",skip = 3,sep="\t", header=T )
# 516 species w/ >50 SLA measurements
# 142 species w/ >100 SLA measurements
# 10 species w/ >500 SLA measurements
  # of the 137000 obs, 107500 are georefed, only 49500 are public

LeafCN_TRY <- read.csv("/Users/leeanderegg/Dropbox/Trade-offs project/TRY_LeafCN_SpeciesTable_03222016imported.csv",header=T )
# 7 species w/ >50 C/N measurements
# 2 species w/ >100 C/N
# 1 species w/ >500 C/N
  # of the 8560 obs, 7020 are goerefed, only 2313 are public

LeafLife_TRY <- read.table("/Users/leeanderegg/Dropbox/Trade-offs project/TRY_LeafLife_SpeciesTable_03222016.txt",skip = 3,sep="\t", header=T )
# 5 species w/ >50 LeafLife
# 3 species w/ >100 LeafLife
# 0 species w/ >500 LeafLife
  # of the 4795 obs, 3097 are georefed, only 1114 are public






