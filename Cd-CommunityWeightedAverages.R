

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------

############ ***Community-weighted trait calculations *** ####################

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# note on nomenclature:
# TRAIT1 = species mean trait value for species 1
# wTRAIT1 = BA weighted species mean trait value
# TRAITp1 = plot/site mean value for the species species 1 - originally used s, but that was horrific. changed 06.14.17
# wTRAITp1 = BA weighted plot/site mean value for species 1



# total species means
#species.means <- traits %>% group_by(GENUS,SPECIES,SP.ID) %>% summarise(mLMA_HSA = mean(LMA_HSA, na.rm=T), mLLmonths = mean(LLmonths, na.rm=T), mCARBON = mean(LEAF_CARBON, na.rm=T), mNITROGEN = mean(LEAF_NITROGEN, na.rm=T), mCN = mean(LEAF_CN, na.rm=T), N = n(), nPlots = length(unique(PLOT_ID)), nProj = length(unique(PROJECT)), nEcoReg = length(unique(ECOREGION)), mMAP = mean(MAP), mMAT = mean(MAT), mlog.LMA = mean(log.LMA, na.rm=T), mlog.LL =mean(log.LL, na.rm=T), mlog.Nmass=mean(log.Nmass, na.rm=T), mlog.LMA_PSA = mean(log.LMA_PSA, na.rm=T), mRGR = mean(RGRdom, na.rm=T), mstGrowthdom = mean(stGrowthdom, na.rm=T)) 
## site means
  # Note: mLMA = LMA_PSA, and mLMA_HSA = LMA_HSA
  # also, mlog.Trait = mean of logged traits
  # log.Trait = log of mean traits

spp.traits <- traits %>% group_by(SP.ID) %>% summarise(nsample = n(), SLA = mean(SLA_HSA, na.rm=T), nSLA = n()- length(which(is.na(SLA_HSA))), CN = mean(LEAF_CN, na.rm=T), nCN = n()- length(which(is.na(LEAF_CN)))
                                                       , LIFE = mean(LEAF_LIFE, na.rm=T), nLIFE=n()- length(which(is.na(LEAF_LIFE))), nplots =length(unique(PLOT_ID))
                                                       , mLMA = mean(LMA_PSA, na.rm=T), mLLmonths= mean(LLmonths, na.rm=T), mNmass = mean(LEAF_NITROGEN, na.rm=T), mNarea=mean(Narea, na.rm=T)
                                                       , mlog.LMA = mean(log.LMA, na.rm=T), mlog.LL = mean(log.LL, na.rm=T), mlog.Narea=mean(log.Narea, na.rm=T), mlog.Nmass=mean(log.Nmass, na.rm=T)
                                                       , RGR = mean(RGRdom, na.rm=T), stGrowth = mean(stGrowthdom, na.rm=T)
                                                       , climPC1 = mean(climPC1, na.rm=T), climPC2 = mean(climPC2, na.rm=T), climPC3 = mean(climPC3, na.rm=T)
                                                       , soil_N= mean(soil_N, na.rm=T), soil_pH=mean(soil_pH, na.rm=T), ASA = mean(ASA, na.rm=T), LAI_O=mean(LAI_O, na.rm=T), AG_TGROWTH = mean(AG_TGROWTH, na.rm=T))
spp.traits$log.LMA <- log(spp.traits$mLMA, base=10)
spp.traits$log.LL <- log(spp.traits$mLLmonths, base=10)
spp.traits$log.Nmass <- log(spp.traits$mNmass, base=10)
spp.traits$log.Narea <- log(spp.traits$mNarea, base=10)

## An idea: leaf lifespan seems to vary quite a bit with MAT. What if I calculated plot averaged leaf lifespan based on spp averages vs based on the actual data from the plot? could say something about using a single value for a species?
spp.plot.traits <- traits %>% group_by(SP.ID, PLOT_ID) %>% summarise(nsample = n(), mSLA = mean(SLA_HSA, na.rm=T), nSLA = n()- length(which(is.na(SLA_HSA))), mCN = mean(LEAF_CN, na.rm=T), nCN = n()- length(which(is.na(LEAF_CN)))
                                                                     , mLIFE = mean(LEAF_LIFE, na.rm=T), nLIFE=n()- length(which(is.na(LEAF_LIFE))), mNmass = mean(LEAF_NITROGEN, na.rm=T), nplots =length(unique(PLOT_ID))
                                                                     , mLMA = mean(LMA_PSA, na.rm=T), mLLmonths= mean(LLmonths, na.rm=T), mNarea=mean(Narea, na.rm=T)
                                                                     , mlog.LMA = mean(log.LMA, na.rm=T), mlog.LL = mean(log.LL, na.rm=T), mlog.Narea=mean(log.Narea, na.rm=T), mlog.Nmass=mean(log.Nmass, na.rm=T)
                                                                     , RGR = unique(RGRdom), stGrowth = mean(stGrowthdom, na.rm=T)
                                                                     , climPC1 = unique(climPC1), climPC2 = unique(climPC2), climPC3 = unique(climPC3)
                                                                     )

# Note: I used LMA_PSA in other trait analysis in TaxonomicAnalysis. So I've switched to it here
spp.plot.traits$SP.PLOT <- paste(spp.plot.traits$SP.ID, spp.plot.traits$PLOT_ID, sep="-")


### make unique identifiers for matching biomass and spp.plot.traits rows
biomass$SP1.PLOT <- paste(biomass$SPP_O1_ABBREV,biomass$PLOT_ID, sep="-")
biomass$SP2.PLOT <- paste(biomass$SPP_O2_ABBREV,biomass$PLOT_ID, sep="-")
biomass$SP3.PLOT <- paste(biomass$SPP_O3_ABBREV,biomass$PLOT_ID, sep="-")
biomass$SP4.PLOT <- paste(biomass$SPP_O4_ABBREV,biomass$PLOT_ID, sep="-")


# let's just pull out the common dominant species in an easy to name df
# plotavs <- spp.plot.traits%>% subset(spp.plot.traits$SP.ID %in% names(which(xtabs(~SPP_O1_ABBREV, biomass)>2)))






#pairs(species.means[,c(1,4,5,6,7,8,9,13,14)], upper.panel = panel.smooth)

#common.genera <- c("Tsuga","Picea","Acer","Quercus","Pinus","Abies")
#plot(mSLA~mMAT, species.means[which(species.means$GENUS %in% common.genera),], col=GENUS, pch=16, ylim=c(0,200))




#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
############ CWMs based on species mean trait values #######################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 
# ####### Plot Trait values based on spp mean values ######
biomass$wLMA1 <- spp.traits$mLMA[match(biomass$SPP_O1_ABBREV, spp.traits$SP.ID)]*biomass$SPP_O1_BASAL_AREA_FRACTION/100
biomass$wLMA2 <- spp.traits$mLMA[match(biomass$SPP_O2_ABBREV, spp.traits$SP.ID)]*biomass$SPP_O2_BASAL_AREA_FRACTION/100
biomass$wLMA3 <- spp.traits$mLMA[match(biomass$SPP_O3_ABBREV, spp.traits$SP.ID)]*biomass$SPP_O3_BASAL_AREA_FRACTION/100
biomass$wLMA4 <- spp.traits$mLMA[match(biomass$SPP_O4_ABBREV, spp.traits$SP.ID)]*biomass$SPP_O4_BASAL_AREA_FRACTION/100
# # # cut out this step to limit the number of columns I create.
# # biomass$wSLA1 <- biomass$SLA1 * biomass$SPP_O1_BASAL_AREA_FRACTION/100
# # biomass$wSLA2 <- biomass$SLA2 * biomass$SPP_O2_BASAL_AREA_FRACTION/100
# # biomass$wSLA3 <- biomass$SLA3 * biomass$SPP_O3_BASAL_AREA_FRACTION/100
# # biomass$wSLA4 <- biomass$SLA4 * biomass$SPP_O4_BASAL_AREA_FRACTION/100
biomass$cw_LMA <- apply(biomass[,c("wLMA1", "wLMA2","wLMA3","wLMA4")],MARGIN=1,FUN=sum, na.rm=T) 
biomass$cw_LMA[which(biomass$cw_LMA==0)] <- NA
# 


biomass$wLL1 <- spp.traits$mLLmonths[match(biomass$SPP_O1_ABBREV, spp.traits$SP.ID)]*biomass$SPP_O1_BASAL_AREA_FRACTION/100
biomass$wLL2 <- spp.traits$mLLmonths[match(biomass$SPP_O2_ABBREV, spp.traits$SP.ID)]*biomass$SPP_O2_BASAL_AREA_FRACTION/100
biomass$wLL3 <- spp.traits$mLLmonths[match(biomass$SPP_O3_ABBREV, spp.traits$SP.ID)]*biomass$SPP_O3_BASAL_AREA_FRACTION/100
biomass$wLL4 <- spp.traits$mLLmonths[match(biomass$SPP_O4_ABBREV, spp.traits$SP.ID)]*biomass$SPP_O4_BASAL_AREA_FRACTION/100

biomass$cw_LL <- apply(biomass[,c("wLL1", "wLL2","wLL3","wLL4")],MARGIN=1,FUN=sum, na.rm=T) 
biomass$cw_LL[which(biomass$cw_LL==0)] <- NA


biomass$wNmass1 <- spp.traits$mNmass[match(biomass$SPP_O1_ABBREV, spp.traits$SP.ID)]*biomass$SPP_O1_BASAL_AREA_FRACTION/100
biomass$wNmass2 <- spp.traits$mNmass[match(biomass$SPP_O2_ABBREV, spp.traits$SP.ID)]*biomass$SPP_O2_BASAL_AREA_FRACTION/100
biomass$wNmass3 <- spp.traits$mNmass[match(biomass$SPP_O3_ABBREV, spp.traits$SP.ID)]*biomass$SPP_O3_BASAL_AREA_FRACTION/100
biomass$wNmass4 <- spp.traits$mNmass[match(biomass$SPP_O4_ABBREV, spp.traits$SP.ID)]*biomass$SPP_O4_BASAL_AREA_FRACTION/100

biomass$cw_Nmass <- apply(biomass[,c("wNmass1", "wNmass2","wNmass3","wNmass4")],MARGIN=1,FUN=sum, na.rm=T) 
biomass$cw_Nmass[which(biomass$cw_Nmass==0)] <- NA



biomass$wNarea1 <- spp.traits$mNarea[match(biomass$SPP_O1_ABBREV, spp.traits$SP.ID)]*biomass$SPP_O1_BASAL_AREA_FRACTION/100
biomass$wNarea2 <- spp.traits$mNarea[match(biomass$SPP_O2_ABBREV, spp.traits$SP.ID)]*biomass$SPP_O2_BASAL_AREA_FRACTION/100
biomass$wNarea3 <- spp.traits$mNarea[match(biomass$SPP_O3_ABBREV, spp.traits$SP.ID)]*biomass$SPP_O3_BASAL_AREA_FRACTION/100
biomass$wNarea4 <- spp.traits$mNarea[match(biomass$SPP_O4_ABBREV, spp.traits$SP.ID)]*biomass$SPP_O4_BASAL_AREA_FRACTION/100

biomass$cw_Narea <- apply(biomass[,c("wNarea1", "wNarea2","wNarea3","wNarea4")],MARGIN=1,FUN=sum, na.rm=T) 
biomass$cw_Narea[which(biomass$cw_Narea==0)] <- NA


# quartz(width=5, height=5)
# par(mfrow=c(2,2), mar=c(3.5,3.5,2,2), oma=c(0,0,0,0), mgp=c(2.5,1,0))
# 
# plot(log.cw_LMAp_if~log(cw_LMA, base=10), biomass, ylab="CW LMA from plot data"
#      , xlab="CW LMA from spp data"
#      , ylim=c(2.1,2.7), xlim=c(2.1,2.7))
# abline(a=0,b=1)
# mtext(text = paste("R2=", round(summary(lm(log.cw_LMAp_if~log(cw_LMA, base=10), biomass))$r.squared,3)))
# 
# plot(log.cw_LLp_if~log(cw_LL, base=10), biomass, ylab="CW LL from plot data", xlab="CW LL from spp data"
#      , ylim=c(1.4,2.4), xlim=c(1.4,2.4))
# abline(a=0,b=1)
# mtext(text = paste("R2=", round(summary(lm(log.cw_LLp_if~log(cw_LL, base=10), biomass))$r.squared,3)))
# 
# plot(log.cw_Nmassp_if~log(cw_Nmass, base=10), biomass, ylab="CW Nmass from plot data", xlab="CW Nmass from spp data"
#      , ylim=c(-0.15,0.2), xlim=c(-0.15,0.2))
# abline(a=0,b=1)
# mtext(text = paste("R2=", round(summary(lm(log.cw_Nmassp_if~log(cw_Nmass, base=10), biomass))$r.squared,3)))
# 
# plot(log.cw_Nareap_if~log(cw_Narea, base=10), biomass, ylab="CW Narea from plot data", xlab="CW Narea from spp data"
#      , ylim=c(0.05,0.6), xlim=c(0.05,0.6))
# abline(a=0,b=1)
# mtext(text = paste("R2=", round(summary(lm(log.cw_Nareap_if~log(cw_Narea, base=10), biomass))$r.squared,3)))

# plot(cw_SLA~MAP, biomass, col=SPP_O1_ABBREV, pch=14 + as.numeric(biomass$PROJECT))
#   # surprisingly linear increase in SLA w/ MAP, despite using spp means
# plot(cw_SLA~MAT_C, biomass, col=SPP_O1_ABBREV, pch=14 + as.numeric(biomass$PROJECT))
# abline(lm(cw_SLA~MAT_C, biomass))
# summary(lm(cw_SLA~MAP + MAT_C + SPP_O1_ABBREV, biomass))
# 
# # ## let's look at how many plots are dominated primarily by one species:
# # length(biomass[which(biomass$SPP_O1_BASAL_AREA_FRACTION==100),1])
# # #55 plots are monoculture
# # length(biomass[which(biomass$SPP_O1_BASAL_AREA_FRACTION>90),1])
# # #116/265 plots are >90% 1 spp
# # length(biomass[which(biomass$SPP_O1_BASAL_AREA_FRACTION>75),1])
# # #150 plots are >75% 1 species.
# # length(biomass[which(biomass$SPP_O1_BASAL_AREA_FRACTION<50),1])
# # #only 39 of 265 plots are bonafide co-dominant...
# 
# 
# 
# 
# ####### **Comparing Plot Level Values based on SPP averages, vs Plot values*** ######
# 
# 
# #### Plot average trait values for each species ############################
# biomass$wSLAp1 <- spp.plot.traits$mSLA[match(as.character(biomass$SP1.PLOT), as.character(spp.plot.traits$SP.PLOT))] * biomass$SPP_O1_BASAL_AREA_FRACTION/100
# biomass$wSLAp2 <- spp.plot.traits$mSLA[match(biomass$SP2.PLOT, spp.plot.traits$SP.PLOT)] * biomass$SPP_O2_BASAL_AREA_FRACTION/100
# biomass$wSLAp3 <- spp.plot.traits$mSLA[match(biomass$SP3.PLOT, spp.plot.traits$SP.PLOT)] * biomass$SPP_O3_BASAL_AREA_FRACTION/100
# biomass$wSLAp4 <- spp.plot.traits$mSLA[match(biomass$SP4.PLOT, spp.plot.traits$SP.PLOT)] * biomass$SPP_O4_BASAL_AREA_FRACTION/100
# 
# # biomass$wSLAp1 <- biomass$SLAp1 * biomass$SPP_O1_BASAL_AREA_FRACTION/100
# # biomass$wSLAp2 <- biomass$SLAp2 * biomass$SPP_O2_BASAL_AREA_FRACTION/100
# # biomass$wSLAp3 <- biomass$SLAp3 * biomass$SPP_O3_BASAL_AREA_FRACTION/100
# # biomass$wSLAp4 <- biomass$SLAp4 * biomass$SPP_O4_BASAL_AREA_FRACTION/100
# biomass$cw_SLAp <- apply(biomass[,c("wSLAp1", "wSLAp2","wSLAp3","wSLAp4")],MARGIN=1,FUN=sum, na.rm=T)
# # now need to remove values that got smoothed over due to na.rm=T in the apply function
# biomass$cw_SLAp[which(is.na(biomass$wSLAp1))] <- NA # 86 plots lacking SLA data of SPP01
# biomass$cw_SLAp[which(is.na(biomass$wSLAp2) & biomass$SPP_O2_BASAL_AREA_FRACTION>30)] <- NA
# #18 sites lacking SLA data for SPP02 where SPP02 makes up >30 % of BA, but only 5 have SPP01
# biomass$cw_SLAp[which(biomass$cw_SLAp==0)] <- NA
#   # 91 plots missing substantial SLA data
# # length(which(is.na(biomass$cw_SLAp)))
#   
# ##### Comparing SLA w/ SLAp ######
# plot(cw_SLAp~cw_SLA, biomass, col=SPP_O1_ABBREV)
# ggplot(biomass, aes(x=cw_SLA, y=cw_SLAp, col=SPP_O1_ABBREV)) + geom_point() + geom_abline(slope=1)
# abline(a=0,b=1)
# summary(lm(cw_SLAp~cw_SLA, biomass))
#   # actually does a pretty shitty job getting cwunity-weighted traits from means rather than the species in the plots
#   # R2= .57
# quartz(width=3.5, height=6)
# par(mfrow=c(2,1), mar=c(3.5,3.5, 1,1), mgp=c(2.5,1,0))
# plot(cw_SLAp~MAP, biomass, col=SPP_O1_ABBREV, pch=16, ylim=c(40,250), ylab="CWM SLA, per site")
# abline(lm(cw_SLAp~MAP, biomass))
# mtext(text = paste("R2=", round(summary(lm(cw_SLAp~MAP, biomass))$r.squared, 3)), side = 3, line=-1.5, adj=.8)
# plot(cw_SLA~MAP, biomass[which(biomass$cw_SLAp>0),], col=SPP_O1_ABBREV, pch=16, ylim=c(40,250), ylab="CWM SLA, spp avg")
# abline(lm(cw_SLA~MAP, biomass[which(biomass$cw_SLAp>0),]))
# mtext(text = paste("R2=", round(summary(lm(cw_SLA~MAP, biomass[which(biomass$cw_SLAp>0),]))$r.squared, 3)), side = 3, line=-1.5, adj=.8)
# 



######### LMA rather than SLA ##############
# can't just calculate this from SLA, because I'm using LMA_PSA rather than HSA
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




#### cwunity weighted new: LLmonths  old: LEAf_LIFE based on plot values
biomass$wLLp1 <- spp.plot.traits$mLLmonths[match(biomass$SP1.PLOT, spp.plot.traits$SP.PLOT)]* biomass$SPP_O1_BASAL_AREA_FRACTION/100
biomass$wLLp2 <- spp.plot.traits$mLLmonths[match(biomass$SP2.PLOT, spp.plot.traits$SP.PLOT)]* biomass$SPP_O2_BASAL_AREA_FRACTION/100
biomass$wLLp3 <- spp.plot.traits$mLLmonths[match(biomass$SP3.PLOT, spp.plot.traits$SP.PLOT)]* biomass$SPP_O3_BASAL_AREA_FRACTION/100
biomass$wLLp4 <- spp.plot.traits$mLLmonths[match(biomass$SP4.PLOT, spp.plot.traits$SP.PLOT)]* biomass$SPP_O4_BASAL_AREA_FRACTION/100

# biomass$wLLp1 <- biomass$LLp1 * biomass$SPP_O1_BASAL_AREA_FRACTION/100
# biomass$wLLp2 <- biomass$LLp2 * biomass$SPP_O2_BASAL_AREA_FRACTION/100
# biomass$wLLp3 <- biomass$LLp3 * biomass$SPP_O3_BASAL_AREA_FRACTION/100
# biomass$wLLp4 <- biomass$LLp4 * biomass$SPP_O4_BASAL_AREA_FRACTION/100
biomass$cw_LLp <- apply(biomass[,c("wLLp1", "wLLp2","wLLp3","wLLp4")],MARGIN=1,FUN=sum, na.rm=T)
# now need to remove values that got smoothed over due to na.rm=T in the apply function
biomass$cw_LLp[which(is.na(biomass$wLLp1))] <- NA # 94 sites lacking LeafLL of SPP01
biomass$cw_LLp[which(is.na(biomass$wLLp2) & biomass$SPP_O2_BASAL_AREA_FRACTION>30)] <- NA
# 22 sites lack SPP02 leaf life, 9 of them unique
length(which(is.na(biomass$cw_LLp))) # 103 plots sans significant LL

#### cwunity weighted LEAf_NITROGEN based on plot values
biomass$wNmassp1 <- spp.plot.traits$mNmass[match(biomass$SP1.PLOT, spp.plot.traits$SP.PLOT)]* biomass$SPP_O1_BASAL_AREA_FRACTION/100
biomass$wNmassp2 <- spp.plot.traits$mNmass[match(biomass$SP2.PLOT, spp.plot.traits$SP.PLOT)]* biomass$SPP_O2_BASAL_AREA_FRACTION/100
biomass$wNmassp3 <- spp.plot.traits$mNmass[match(biomass$SP3.PLOT, spp.plot.traits$SP.PLOT)]* biomass$SPP_O3_BASAL_AREA_FRACTION/100
biomass$wNmassp4 <- spp.plot.traits$mNmass[match(biomass$SP4.PLOT, spp.plot.traits$SP.PLOT)]* biomass$SPP_O4_BASAL_AREA_FRACTION/100

# biomass$wNmassp1 <- biomass$Nmassp1 * biomass$SPP_O1_BASAL_AREA_FRACTION/100
# biomass$wNmassp2 <- biomass$Nmassp2 * biomass$SPP_O2_BASAL_AREA_FRACTION/100
# biomass$wNmassp3 <- biomass$Nmassp3 * biomass$SPP_O3_BASAL_AREA_FRACTION/100
# biomass$wNmassp4 <- biomass$Nmassp4 * biomass$SPP_O4_BASAL_AREA_FRACTION/100
biomass$cw_Nmassp <- apply(biomass[,c("wNmassp1", "wNmassp2","wNmassp3","wNmassp4")],MARGIN=1,FUN=sum, na.rm=T)
# now need to remove values that got smoothed over due to na.rm=T in the apply function
biomass$cw_Nmassp[which(is.na(biomass$wNmassp1))] <- NA # 98 sites lacking LeafNITROGEN of SPP01
biomass$cw_Nmassp[which(is.na(biomass$wNmassp2) & biomass$SPP_O2_BASAL_AREA_FRACTION>30)] <- NA
biomass$cw_Nmassp[which(biomass$cw_Nmassp==0)] <- NA
# 26 sites lack SPP02 leaf NITROGEN, 3 of them unique
# length(which(is.na(biomass$cw_Nmassp))) # 102 sites sans Nmass




#### cwunity weighted Narea based on plot values
biomass$wNareap1 <- spp.plot.traits$mNarea[match(biomass$SP1.PLOT, spp.plot.traits$SP.PLOT)]* biomass$SPP_O1_BASAL_AREA_FRACTION/100
biomass$wNareap2 <- spp.plot.traits$mNarea[match(biomass$SP2.PLOT, spp.plot.traits$SP.PLOT)]* biomass$SPP_O2_BASAL_AREA_FRACTION/100
biomass$wNareap3 <- spp.plot.traits$mNarea[match(biomass$SP3.PLOT, spp.plot.traits$SP.PLOT)]* biomass$SPP_O3_BASAL_AREA_FRACTION/100
biomass$wNareap4 <- spp.plot.traits$mNarea[match(biomass$SP4.PLOT, spp.plot.traits$SP.PLOT)]* biomass$SPP_O4_BASAL_AREA_FRACTION/100
# 
# biomass$wNareap1 <- biomass$Nareap1 * biomass$SPP_O1_BASAL_AREA_FRACTION/100
# biomass$wNareap2 <- biomass$Nareap2 * biomass$SPP_O2_BASAL_AREA_FRACTION/100
# biomass$wNareap3 <- biomass$Nareap3 * biomass$SPP_O3_BASAL_AREA_FRACTION/100
# biomass$wNareap4 <- biomass$Nareap4 * biomass$SPP_O4_BASAL_AREA_FRACTION/100
biomass$cw_Nareap <- apply(biomass[,c("wNareap1", "wNareap2","wNareap3","wNareap4")],MARGIN=1,FUN=sum, na.rm=T)
# now need to remove values that got smoothed over due to na.rm=T in the apply function
biomass$cw_Nareap[which(is.na(biomass$wNareap1))] <- NA # 98 sites lacking LeafNarea of SPP01
biomass$cw_Nareap[which(is.na(biomass$wNareap2) & biomass$SPP_O2_BASAL_AREA_FRACTION>30)] <- NA
biomass$cw_Nareap[which(biomass$cw_Nareap==0)] <- NA
# 26 sites lack SPP02 leaf Narea, 3 of them unique




biomass$log.cw_LMAp <- log(biomass$cw_LMAp, base=10)
biomass$log.cw_LLp <- log(biomass$cw_LLp, base=10)
biomass$log.cw_Nmassp <- log(biomass$cw_Nmassp, base=10)
biomass$log.cw_Nareap <- log(biomass$cw_Nareap, base=10)

###### add the climPCs to biomass########
biomass$climPC1 <-spp.plot.traits$climPC1[match(biomass$SP1.PLOT, spp.plot.traits$SP.PLOT)]
biomass$climPC2 <-spp.plot.traits$climPC2[match(biomass$SP1.PLOT, spp.plot.traits$SP.PLOT)]



#_____________________________________________________________________________
###### Infilling spp trait values with spp means when they're missing ###########
#_____________________________________________________________________________
### Infilling LMA of spp with <30% of BA with species mean values

####### Creating BA-weighted LMA values using spp means
biomass$wLMA1 <- spp.traits$mLMA[match(biomass$SPP_O1_ABBREV, spp.traits$SP.ID)]*biomass$SPP_O1_BASAL_AREA_FRACTION/100
biomass$wLMA2 <- spp.traits$mLMA[match(biomass$SPP_O2_ABBREV, spp.traits$SP.ID)]*biomass$SPP_O2_BASAL_AREA_FRACTION/100
biomass$wLMA3 <- spp.traits$mLMA[match(biomass$SPP_O3_ABBREV, spp.traits$SP.ID)]*biomass$SPP_O3_BASAL_AREA_FRACTION/100
biomass$wLMA4 <- spp.traits$mLMA[match(biomass$SPP_O4_ABBREV, spp.traits$SP.ID)]*biomass$SPP_O4_BASAL_AREA_FRACTION/100
# # cut out this step to limit the number of columns I create.
# biomass$wLMA1 <- biomass$LMA1 * biomass$SPP_O1_BASAL_AREA_FRACTION/100
# biomass$wLMA2 <- biomass$LMA2 * biomass$SPP_O2_BASAL_AREA_FRACTION/100
# biomass$wLMA3 <- biomass$LMA3 * biomass$SPP_O3_BASAL_AREA_FRACTION/100
# biomass$wLMA4 <- biomass$LMA4 * biomass$SPP_O4_BASAL_AREA_FRACTION/100
biomass$cw_LMA <- apply(biomass[,c("wLMA1", "wLMA2","wLMA3","wLMA4")],MARGIN=1,FUN=sum, na.rm=T)
biomass$cw_LMA[which(biomass$cw_LMA==0)] <- NA
# summary(lm(cw_LMAp~cw_LMA, biomass))
# cw based on spp averages captures 81% of the variation in actual CW based on plot means (though this has some noise from the points I'm trying to infill...)

wLMAif1 <- biomass$wLMAp1
wLMAif1[which(is.na(wLMAif1))] <- biomass$wLMA1[which(is.na(wLMAif1))]
wLMAif2 <- biomass$wLMAp2
wLMAif2[which(is.na(wLMAif2))] <- biomass$wLMA2[which(is.na(wLMAif2))]
wLMAif3 <- biomass$wLMAp3
wLMAif3[which(is.na(wLMAif3))] <- biomass$wLMA3[which(is.na(wLMAif3))]
wLMAif4 <- biomass$wLMAp4
wLMAif4[which(is.na(wLMAif4))] <- biomass$wLMA4[which(is.na(wLMAif4))]
cw_LMAp_if  <- apply(data.frame(wLMAif1, wLMAif2, wLMAif3, wLMAif4),MARGIN=1, FUN=sum, na.rm=T)
# get rid of plots that still have no data
cw_LMAp_if[which(cw_LMAp_if==0)] <- NA
cw_LMAp_if[which(is.na(biomass$wLMAp1) & biomass$SPP_O1_BASAL_AREA_FRACTION>25)] <- NA # 80 plots lacking LMA data of SPP01 & SPO1 has >25% of the basal area
cw_LMAp_if[which(is.na(biomass$wLMAp2) & biomass$SPP_O2_BASAL_AREA_FRACTION>25)] <- NA # 26 plots lacking LMA of SPP2 w>25% of BA, 6 plots that weren't killed above
cw_LMAp_if[which(is.na(biomass$wLMAp3) & biomass$SPP_O3_BASAL_AREA_FRACTION>25)] <- NA # 2 plots lacking a prevalent SPP3, only one novel
#cw_LMAp_if[which(is.na(biomass$wLMAp4) & biomass$SPP_O4_BASAL_AREA_FRACTION>25)] <- NA
# Yup, there were a fair # of plots that I biased cw_LMAp down because of missing spp that were pretty well represented

biomass$cw_LMAp_if <- cw_LMAp_if
# summary(lm(cw_LMAp_if~cw_LMA, biomass))
# but still miss 18% of the plot to plot variation with species means
# length(which(biomass$cw_LMAp_if != biomass$cw_LMAp))
# corrected 78 plots



##### Leaf Lifespan infilling
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
#cw_LLp_if[which(is.na(biomass$wLLp4) & biomass$SPP_O4_BASAL_AREA_FRACTION>25)] <- NA
# Yup, there were a fair # of plots that I biased cw_LLp down because of missing spp that were pretty well represented
#plot(cw_LLp_if~cw_LLp, biomass)

biomass$cw_LLp_if <- cw_LLp_if
# summary(lm(cw_LLp_if~cw_LL, biomass))
# but still miss 36% !!! of the plot to plot variation with species means
# length(which(biomass$cw_LLp_if != biomass$cw_LLp))
# corrected 74 plots


##### Nmass infilling
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
#cw_Nmassp_if[which(is.na(biomass$wNmassp4) & biomass$SPP_O4_BASAL_AREA_FRACTION>25)] <- NA
# Yup, there were a fair # of plots that I biased cw_Nmassp down because of missing spp that were pretty well represented
#plot(cw_Nmassp_if~cw_Nmassp, biomass)
# Shit! this was actually pretty damn important!

biomass$cw_Nmassp_if <- cw_Nmassp_if
# summary(lm(cw_Nmassp_if~cw_Nmass, biomass))
# but still miss 55% !!! of the plot to plot variation with species means
# length(which(biomass$cw_Nmassp_if != biomass$cw_Nmassp))
# corrected 68 plots



##### Narea infilling
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
#cw_Nareap_if[which(is.na(biomass$wNareap4) & biomass$SPP_O4_BASAL_AREA_FRACTION>25)] <- NA
# Yup, there were a fair # of plots that I biased cw_Nareap down because of missing spp that were pretty well represented
# plot(cw_Nareap_if~cw_Nareap, biomass)


biomass$cw_Nareap_if <- cw_Nareap_if
# summary(lm(cw_Nareap_if~cw_Narea, biomass))
# but still miss 45% !!! of the plot to plot variation with species means
# length(which(biomass$cw_Nareap_if != biomass$cw_Nareap))
# corrected 68 plots




biomass$log.cw_LMAp_if <- log(biomass$cw_LMAp_if, base=10)
biomass$log.cw_LLp_if <- log(biomass$cw_LLp_if, base=10)
biomass$log.cw_Nmassp_if <- log(biomass$cw_Nmassp_if, base=10)
biomass$log.cw_Nareap_if <- log(biomass$cw_Nareap_if, base=10)













############ OLD CODE: ###############################


# #______________________________________________________
# ######### Calculating with mean log(traits) rather than mean traits that are then logged ##########
#   #---> I realized this is a stupid way to calculate this. I definetly want to log AFTER I cw things
# ######### log.LMA 
# # can't just calculate this from SLA, because I'm using LMA_PSA rather than HSA
# biomass$wlog.LMAp1 <- spp.plot.traits$mlog.LMA[match(as.character(biomass$SP1.PLOT), as.character(spp.plot.traits$SP.PLOT))] * biomass$SPP_O1_BASAL_AREA_FRACTION/100
# biomass$wlog.LMAp2 <- spp.plot.traits$mlog.LMA[match(biomass$SP2.PLOT, spp.plot.traits$SP.PLOT)] * biomass$SPP_O2_BASAL_AREA_FRACTION/100
# biomass$wlog.LMAp3 <- spp.plot.traits$mlog.LMA[match(biomass$SP3.PLOT, spp.plot.traits$SP.PLOT)] * biomass$SPP_O3_BASAL_AREA_FRACTION/100
# biomass$wlog.LMAp4 <- spp.plot.traits$mlog.LMA[match(biomass$SP4.PLOT, spp.plot.traits$SP.PLOT)] * biomass$SPP_O4_BASAL_AREA_FRACTION/100
# 
# biomass$cw_log.LMAp <- apply(biomass[,c("wlog.LMAp1", "wlog.LMAp2","wlog.LMAp3","wlog.LMAp4")],MARGIN=1,FUN=sum, na.rm=T)
# # now need to remove values that got smoothed over due to na.rm=T in the apply function
# biomass$cw_log.LMAp[which(is.na(biomass$wlog.LMAp1))] <- NA # 86 plots lacking log.LMA data of SPP01
# biomass$cw_log.LMAp[which(is.na(biomass$wlog.LMAp2) & biomass$SPP_O2_BASAL_AREA_FRACTION>30)] <- NA
# #18 sites lacking log.LMA data for SPP02 where SPP02 makes up >30 % of BA, but only 5 have SPP01
# biomass$cw_log.LMAp[which(biomass$cw_log.LMAp==0)] <- NA
# # 91 plots missing substantial LMA data
# # length(which(is.na(biomass$cw_LMAp)))
# 
# 
# ###### log.LL
# biomass$wlog.LLp1 <- spp.plot.traits$mlog.LL[match(biomass$SP1.PLOT, spp.plot.traits$SP.PLOT)]* biomass$SPP_O1_BASAL_AREA_FRACTION/100
# biomass$wlog.LLp2 <- spp.plot.traits$mlog.LL[match(biomass$SP2.PLOT, spp.plot.traits$SP.PLOT)]* biomass$SPP_O2_BASAL_AREA_FRACTION/100
# biomass$wlog.LLp3 <- spp.plot.traits$mlog.LL[match(biomass$SP3.PLOT, spp.plot.traits$SP.PLOT)]* biomass$SPP_O3_BASAL_AREA_FRACTION/100
# biomass$wlog.LLp4 <- spp.plot.traits$mlog.LL[match(biomass$SP4.PLOT, spp.plot.traits$SP.PLOT)]* biomass$SPP_O4_BASAL_AREA_FRACTION/100
# 
# biomass$cw_log.LLp <- apply(biomass[,c("wlog.LLp1", "wlog.LLp2","wlog.LLp3","wlog.LLp4")],MARGIN=1,FUN=sum, na.rm=T)
# # now need to remove values that got smoothed over due to na.rm=T in the apply function
# biomass$cw_log.LLp[which(is.na(biomass$wlog.LLp1))] <- NA # 94 sites lacking Leaflog.LL of SPP01
# biomass$cw_log.LLp[which(is.na(biomass$wlog.LLp2) & biomass$SPP_O2_BASAL_AREA_FRACTION>30)] <- NA
# # 22 sites lack SPP02 leaf life, 9 of them unique
# length(which(is.na(biomass$cw_log.LLp))) # 103 plots sans significant log.LL
# 
# ###### log.Nmass
# biomass$wlog.Nmassp1 <- spp.plot.traits$mlog.Nmass[match(biomass$SP1.PLOT, spp.plot.traits$SP.PLOT)]* biomass$SPP_O1_BASAL_AREA_FRACTION/100
# biomass$wlog.Nmassp2 <- spp.plot.traits$mlog.Nmass[match(biomass$SP2.PLOT, spp.plot.traits$SP.PLOT)]* biomass$SPP_O2_BASAL_AREA_FRACTION/100
# biomass$wlog.Nmassp3 <- spp.plot.traits$mlog.Nmass[match(biomass$SP3.PLOT, spp.plot.traits$SP.PLOT)]* biomass$SPP_O3_BASAL_AREA_FRACTION/100
# biomass$wlog.Nmassp4 <- spp.plot.traits$mlog.Nmass[match(biomass$SP4.PLOT, spp.plot.traits$SP.PLOT)]* biomass$SPP_O4_BASAL_AREA_FRACTION/100
# 
# biomass$cw_log.Nmassp <- apply(biomass[,c("wlog.Nmassp1", "wlog.Nmassp2","wlog.Nmassp3","wlog.Nmassp4")],MARGIN=1,FUN=sum, na.rm=T)
# # now need to remove values that got smoothed over due to na.rm=T in the apply function
# biomass$cw_log.Nmassp[which(is.na(biomass$wlog.Nmassp1))] <- NA # 98 sites lacking LeafNITROGEN of SPP01
# biomass$cw_log.Nmassp[which(is.na(biomass$wlog.Nmassp2) & biomass$SPP_O2_BASAL_AREA_FRACTION>30)] <- NA
# biomass$cw_log.Nmassp[which(biomass$cw_log.Nmassp==0)] <- NA
# # 26 sites lack SPP02 leaf NITROGEN, 3 of them unique
# # length(which(is.na(biomass$cw_log.Nmassp))) # 102 sites sans log.Nmass
# 
# ##### log.Narea
# biomass$wlog.Nareap1 <- spp.plot.traits$mlog.Narea[match(biomass$SP1.PLOT, spp.plot.traits$SP.PLOT)]* biomass$SPP_O1_BASAL_AREA_FRACTION/100
# biomass$wlog.Nareap2 <- spp.plot.traits$mlog.Narea[match(biomass$SP2.PLOT, spp.plot.traits$SP.PLOT)]* biomass$SPP_O2_BASAL_AREA_FRACTION/100
# biomass$wlog.Nareap3 <- spp.plot.traits$mlog.Narea[match(biomass$SP3.PLOT, spp.plot.traits$SP.PLOT)]* biomass$SPP_O3_BASAL_AREA_FRACTION/100
# biomass$wlog.Nareap4 <- spp.plot.traits$mlog.Narea[match(biomass$SP4.PLOT, spp.plot.traits$SP.PLOT)]* biomass$SPP_O4_BASAL_AREA_FRACTION/100
# #
# biomass$cw_log.Nareap <- apply(biomass[,c("wlog.Nareap1", "wlog.Nareap2","wlog.Nareap3","wlog.Nareap4")],MARGIN=1,FUN=sum, na.rm=T)
# # now need to remove values that got smoothed over due to na.rm=T in the apply function
# biomass$cw_log.Nareap[which(is.na(biomass$wlog.Nareap1))] <- NA # 98 sites lacking Leaflog.Narea of SPP01
# biomass$cw_log.Nareap[which(is.na(biomass$wlog.Nareap2) & biomass$SPP_O2_BASAL_AREA_FRACTION>30)] <- NA
# biomass$cw_log.Nareap[which(biomass$cw_log.Nareap==0)] <- NA
# # 26 sites lack SPP02 leaf log.Narea, 3 of them unique
# 
# # lost roughly 100 sites per trait, leaving ~160 sites. not bad, but infilling might be better?




### Plot what happens if I log before or after averaging ##
# 
# quartz(width=7, height=6)
# par(mfrow=c(2,2),mar=c(3.5,3.5,1,1), mgp=c(2.5,1,0))
# plot(log.cw_LMAp~cw_log.LMAp, biomass)
# plot(log.cw_LLp~cw_log.LLp, biomass)
# plot(log.cw_Nareap~cw_log.Nareap, biomass)
# plot(log.cw_Nmassp~cw_log.Nmassp, biomass)
# 








###### looking for trade-offs at the CW mean level ####
myspecies <- c("ABICON","ABIGRA","PINPON","PSEMEN","PINJEF")
plot(log.cw_LLp_if~log.cw_LMAp_if, biomass[which(!is.na(biomass$log.cw_LLp_if)),], col=factor(SPP_O1_ABBREV), pch=16)
for(i in myspecies){
  abline(lm(log.cw_LLp_if~log.cw_LMAp_if, biomass[which(!is.na(biomass$log.cw_LLp_if) & biomass$SPP_O1_ABBREV==i),]))
  
}
#### woof. doesn't work at all for CW LL vs LMA
plot(log.cw_Nareap_if~log.cw_LMAp_if, biomass[which(!is.na(biomass$log.cw_LLp_if)),], col=factor(SPP_O1_ABBREV), pch=16)
for(i in myspecies){
  abline(lm(log.cw_Nareap_if~log.cw_LMAp_if, biomass[which(!is.na(biomass$log.cw_LLp_if) & biomass$SPP_O1_ABBREV==i),]))
  
}
#### works super slick for Narea and LMA, which says N is mass-distributed, I think

plot(log.cw_Nmassp_if~log.cw_LMAp_if, biomass[which(!is.na(biomass$log.cw_LLp_if)),], col=factor(SPP_O1_ABBREV), pch=16)
for(i in myspecies){
  abline(lm(log.cw_Nmassp_if~log.cw_LMAp_if, biomass[which(!is.na(biomass$log.cw_LLp_if) & biomass$SPP_O1_ABBREV==i),]))
  
}
# doesn't work at all for Nmass

plot(log.cw_Nareap_if~log.cw_LLp_if, biomass[which(!is.na(biomass$log.cw_LLp_if)),], col=factor(SPP_O1_ABBREV), pch=16)
for(i in myspecies){
  abline(lm(log.cw_Nareap_if~log.cw_LLp_if, biomass[which(!is.na(biomass$log.cw_LLp_if) & biomass$SPP_O1_ABBREV==i),]))
  
}

# nada for Narea vs LL












### making a 'canopy complexity' type metric #######
biomass$Ccomplex <- biomass$AG_BIOMASS_TREE_WOOD_AS_CARBON/biomass$LAI_O






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
wLIFEif1 <- biomass$wLIFEp1
wLIFEif1[which(is.na(wLIFEif1))] <- wLIFE1[which(is.na(wLIFEif1))]
wLIFEif2 <- biomass$wLIFEp2
wLIFEif2[which(is.na(wLIFEif2))] <- wLIFE2[which(is.na(wLIFEif2))]
wLIFEif3 <- biomass$wLIFEp3
wLIFEif3[which(is.na(wLIFEif3))] <- wLIFE3[which(is.na(wLIFEif3))]
wLIFEif4 <- biomass$wLIFEp4
wLIFEif4[which(is.na(wLIFEif4))] <- wLIFE4[which(is.na(wLIFEif4))]
cw_LIFEp_if  <- apply(data.frame(wLIFEif1, wLIFEif2, wLIFEif3, wLIFEif4),MARGIN=1, FUN=sum, na.rm=T)
cw_LIFEp_if[which(is.na(wLIFEif1))] <- NA
which(is.na(biomass$cw_LIFEp) & cw_LIFEp_if>0)
# infilled 90 plots

##### Infilling LEAF_NITROGEN with species mean values ######
# note, need to make species average values first.

Nmass1 <- species.means$mNmass[match(biomass$SPP_O1_ABBREV, species.means$SP.ID)]
Nmass2 <- species.means$mNmass[match(biomass$SPP_O2_ABBREV, species.means$SP.ID)]
Nmass3 <- species.means$mNmass[match(biomass$SPP_O3_ABBREV, species.means$SP.ID)]
Nmass4 <- species.means$mNmass[match(biomass$SPP_O4_ABBREV, species.means$SP.ID)]
wNmass1 <- Nmass1 * biomass$SPP_O1_BASAL_AREA_FRACTION/100 # 13 total plots either lack Nmass for the dominant spp or don't have SPP data
wNmass2 <- Nmass2 * biomass$SPP_O2_BASAL_AREA_FRACTION/100
wNmass3 <- Nmass3 * biomass$SPP_O3_BASAL_AREA_FRACTION/100
wNmass4 <- Nmass4 * biomass$SPP_O4_BASAL_AREA_FRACTION/100

### now infill missing values with species means
wNmassif1 <- biomass$wNmassp1
wNmassif1[which(is.na(wNmassif1))] <- wNmass1[which(is.na(wNmassif1))]
wNmassif2 <- biomass$wNmassp2
wNmassif2[which(is.na(wNmassif2))] <- wNmass2[which(is.na(wNmassif2))]
wNmassif3 <- biomass$wNmassp3
wNmassif3[which(is.na(wNmassif3))] <- wNmass3[which(is.na(wNmassif3))]
wNmassif4 <- biomass$wNmassp4
wNmassif4[which(is.na(wNmassif4))] <- wNmass4[which(is.na(wNmassif4))]
cw_Nmassp_if  <- apply(data.frame(wNmassif1, wNmassif2, wNmassif3, wNmassif4),MARGIN=1, FUN=sum, na.rm=T)
cw_Nmassp_if[which(is.na(wNmassif1))] <- NA
length(which(is.na(biomass$cw_Nmassp) & cw_Nmassp_if>0))
# infilled 96 plots

# plot(AG_PROD_TREE_TOTAL_AS_CARBON~cw_Nmassp_if, biomass, pch=16, cex=1.1)
# points(AG_PROD_TREE_TOTAL_AS_CARBON~cw_Nmassp, biomass, pch=16, cex=.9, col="red")

length(biomass$PLOT_ID[which(biomass$LAI_O>0 & biomass$AG_PROD_TREE_TOTAL_AS_CARBON>0 & cw_Nmassp_if>0 & cw_LIFEp_if>0 & cw_SLAs_if>0)]) 
length(biomass$PLOT_ID[which(biomass$LAI_O>0 & biomass$AG_PROD_TREE_TOTAL_AS_CARBON>0 & biomass$cw_Nmassp>0 & biomass$cw_LIFEp>0 & biomass$cw_SLAs>0)]) 

##### Export best dataset ## NOTE: THIS IS HARD CODED AND MIGHT NOT WORK IF COLUMNS ARE CHANGED
# combine things together and remove the columns we don't actually need
biomassbest <- data.frame(biomass[,-c(35:43,48:55,57:64,66:73)],cw_SLAs_if,cw_LIFEp_if,cw_Nmassp_if)
colnames(biomassbest)[grep("cw", colnames(biomassbest))] <- c("plotSLA", "plotLIFE","plotNmass","plotSLA_if","plotLIFE_if","plotNmass_if")
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