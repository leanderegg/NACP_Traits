

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


spp.traits <- traits %>% group_by(GE.SP) %>% summarise(nsample = n(), SLA = mean(SLA_HSA, na.rm=T), nSLA = n()- length(which(is.na(SLA_HSA))), CN = mean(LEAF_CN, na.rm=T), nCN = n()- length(which(is.na(LEAF_CN)))
                                                       , LIFE = mean(LEAF_LIFE, na.rm=T), nLIFE=n()- length(which(is.na(LEAF_LIFE))), nplots =length(unique(PLOT_ID))
                                                       , mLMA = mean(LMA_PSA, na.rm=T), mLLmonths= mean(LLmonths, na.rm=T), mNarea=mean(Narea, na.rm=T)
                                                       , mlog.LMA = mean(log.LMA, na.rm=T), mlog.LL = mean(log.LL, na.rm=T), mlog.Narea=mean(log.Narea, na.rm=T), mlog.Nmass=mean(log.Nmass, na.rm=T)
                                                       , RGR = unique(RGRdom), stGrowth = mean(stGrowthdom, na.rm=T) )
## An idea: leaf lifespan seems to vary quite a bit with MAT. What if I calculated plot averaged leaf lifespan based on spp averages vs based on the actual data from the plot? could say something about using a single value for a species?
spp.plot.traits <- traits %>% group_by(SP.ID, PLOT_ID) %>% summarise(nsample = n(), mSLA = mean(SLA_HSA, na.rm=T), nSLA = n()- length(which(is.na(SLA_HSA))), mCN = mean(LEAF_CN, na.rm=T), nCN = n()- length(which(is.na(LEAF_CN)))
                                                                     , mLIFE = mean(LEAF_LIFE, na.rm=T), nLIFE=n()- length(which(is.na(LEAF_LIFE))), mNITROGEN = mean(LEAF_NITROGEN, na.rm=T), nplots =length(unique(PLOT_ID))
                                                                     , mLMA = mean(LMA_PSA, na.rm=T), mLLmonths= mean(LLmonths, na.rm=T), mNarea=mean(Narea, na.rm=T)
                                                                     , mlog.LMA = mean(log.LMA, na.rm=T), mlog.LL = mean(log.LL, na.rm=T), mlog.Narea=mean(log.Narea, na.rm=T), mlog.Nmass=mean(log.Nmass, na.rm=T)
                                                                     , RGR = unique(RGRdom), stGrowth = mean(stGrowthdom, na.rm=T))
# Note: I used LMA_PSA in other trait analysis in TaxonomicAnalysis. So I've switched to it here
spp.plot.traits$SP.PLOT <- paste(spp.plot.traits$SP.ID, spp.plot.traits$PLOT_ID, sep="-")
# let's just pull out the common dominant species in an easy to name df
plotavs <- spp.plot.traits%>% subset(spp.plot.traits$SP.ID %in% names(which(xtabs(~SPP_O1_ABBREV, biomass)>2)))

### make unique identifiers for matching biomass and spp.plot.traits rows
biomass$SP1.PLOT <- paste(biomass$SPP_O1_ABBREV,biomass$PLOT_ID, sep="-")
biomass$SP2.PLOT <- paste(biomass$SPP_O2_ABBREV,biomass$PLOT_ID, sep="-")
biomass$SP3.PLOT <- paste(biomass$SPP_O3_ABBREV,biomass$PLOT_ID, sep="-")
biomass$SP4.PLOT <- paste(biomass$SPP_O4_ABBREV,biomass$PLOT_ID, sep="-")




#pairs(species.means[,c(1,4,5,6,7,8,9,13,14)], upper.panel = panel.smooth)

#common.genera <- c("Tsuga","Picea","Acer","Quercus","Pinus","Abies")
#plot(mSLA~mMAT, species.means[which(species.means$GENUS %in% common.genera),], col=GENUS, pch=16, ylim=c(0,200))




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


biomass.com <- read.csv("PACNW_Biomass_plus_traits_120716.csv")

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