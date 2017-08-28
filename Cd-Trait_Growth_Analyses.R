

#-----------------------------------------------------------------------------
####### Trait-Growth Analyses ###############
#-----------------------------------------------------------------------------






#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
######### OLD ANALYSES AND PLOTS #########################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#________________________________________________________________________________
####### .OLD relationships, Leaf Biomass : Stem Biomass and SLA ###############
#________________________________________________________________________________


### visualization of how Growth relates to Biomass and how I calculated relative growth rate:
ggplot(biomass, aes(x=log(AG_BIOMASS_TREE_TOTAL_AS_CARBON), y=AG_PROD_TREE_TOTAL_AS_CARBON, col=SPP_O1_ABBREV)) + 
  geom_point() + xlab("log(Stand Biomass)") + ylab("Biomass Growth") + 
  geom_smooth(aes(colour=NULL))


# observation: AG_FBIO is strongly related to LAI_O across species and sites. BUT, AG_WBIO is less so
# -> at the whole stand scale, there might be a relationship between the SLA (investment in area vs mass) and Leaf:Stem investment
# *** This would probably best be done at the Community-weighted-trait level?


## visualizing with traitst.common instead
plot(AG_TGROWTH~log(AG_TBIO), traits.common, col=FOREST_TYPE)
for(i in unique(traits.common$FOREST_TYPE)[-5]){
  mod <- lm(AG_TGROWTH~log(AG_TBIO), traits.common[which(traits.common$FOREST_TYPE==i),])
  lines(fitted(mod)~log(na.omit(traits.common$AG_TBIO[which(traits.common$FOREST_TYPE==i)])))
}
abline(lm(AG_TGROWTH~log(HEIGHTC_m), traits.common), lwd=2)

ggplot(traits.common, aes(x=log(AG_TBIO), y=AG_TGROWTH, col=FOREST_TYPE)) + geom_point() + geom_smooth(se=F, method = "lm") 
  # wow, within cover types, the relationship between growth and raw above ground total biomass is pretty flat but changes a lot between cover types
  # but relationships are reasonably linear for different forest types when plotted against log(Biomass)
# summary(lm(AG_PROD_TREE_TOTAL_AS_CARBON~log(AG_BIOMASS_TREE_TOTAL_AS_CARBON), biomass))
## NOTE: Height explains 0.079 of GROWTH
#         TBIO explains .181 of GROWTH (or 24% in biomass), actually 24% (raw TBIO) or 29% (log(TBIO)) if run on Biomass
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



#_________________________
##### .WRONG Old plots, using all traits rather than plot means ###########
#_________________________


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





#--- old RGR calculations ---
# ## Leaf mass to stem mass ratio for these plots ###
 traits.common$LtoS <- with(traits.common, AG_FBIO/AG_WBIO)
# traits.common$RGR.old <- with(traits.common, AG_TGROWTH/AG_TBIO)
# traits.common$lnRGR.old <- with(traits.common, log(AG_TGROWTH)/AG_TBIO)
# traits.common$RGR[which(traits.common$RGR>2)] <- NA



######### Playing around with Leaf to Shoot Ratio 
plot(SLA_HSA~LtoS, traits.common[which(traits.common$dominance>90),], col=SP.ID, xlim=c(0,1)) # note: there are a couple weird plots to be removed
for(i in unique(traits.common$SP.ID)){
  abline(lm(SLA_HSA~LtoS, traits.common, subset=SP.ID==i&dominance>90), col=traits.common$SP.ID[which(traits.common$SP.ID==i)])
}
abline(lm(SLA_HSA~LtoS, traits.common, subset=dominace>90), lwd=3)
## Hmm. less of a pattern there than I would have thought...


#______________________________
################ .More OLD PLOTS on indiv trait measurements ##################
#______________________________

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








####### .Third pass with RGR (log growth) ################
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






############## .Growth v traits, only in dominants ########################################
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
#### .OLD FIGURE: RGR versus log TRAITS #################################
#_____________________________________________________________________
  # this was the original analysis, before I realized I should just run it on stand average trait values
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
quartz(width=6.8, height=3) # Eco Let width is 173mm, 110mm or 82mm
par(mar=c(4,2,1.2,1), oma=c(0,2,0,0),mfrow=c(1,3), cex=1,mgp=c(2.5,1,0))

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



####################### END: old analyses and Figures with individual triat measurements ############################








#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
############ ****** Trait v RGR (final with plot means) *** ##################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# ok, I realized that all of the above are taking plot level averages for RGR and I really should be taking plot level averages for traits to.

spp.plot.traits <- traits.common5 %>% group_by(SP.ID, PLOT_ID) %>% summarise(nsample = n(), mSLA = mean(SLA_HSA, na.rm=T), nSLA = n()- length(which(is.na(SLA_HSA))), mCN = mean(LEAF_CN, na.rm=T), nCN = n()- length(which(is.na(LEAF_CN)))
                                                                            , mLIFE = mean(LEAF_LIFE, na.rm=T), nLIFE=n()- length(which(is.na(LEAF_LIFE))), mNITROGEN = mean(LEAF_NITROGEN, na.rm=T), nplots =length(unique(PLOT_ID))
                                                                            , mLMA = mean(LMA_PSA, na.rm=T), mLLmonths= mean(LLmonths, na.rm=T), mNarea=mean(Narea, na.rm=T)
                                                                            , mlog.LMA = mean(log.LMA, na.rm=T), mlog.LL = mean(log.LL, na.rm=T), mlog.Narea=mean(log.Narea, na.rm=T), mlog.Nmass=mean(log.Nmass, na.rm=T)
                                                                            , RGR = unique(RGRdom), stGrowth = mean(stGrowthdom, na.rm=T)
                                                                            , RGR90 = unique(RGRdom90), stGrowth90 = mean(stGrowthdom90, na.rm=T)
                                                                            , climPC1 = unique(climPC1), climPC2 = unique(climPC2), climPC3=unique(climPC3), soil_N = unique(soil_N), soil_pH = unique(soil_pH),ASA=unique(ASA)
                                                                            , NPP = unique(AG_TGROWTH), Biomass = unique(AG_TBIO), LeafFrac = unique(LeafFrac), LAI = unique(LAI_O), Height = unique(HEIGHTC_m))
# Note: I used LMA_PSA in other trait analysis in TaxonomicAnalysis. So I've switched to it here
# Second Note: from here on out RGR is only of dominant spp because I used RGRdom, which was only calculated for things
  # trait measurements of the FOREST_TYPE species that also made up >50% of the stand Basal Area.
  # anything that's not a 'dominant' species ends up being NA.
#--> I just need to update RGRdom in traits df in order to rerun this with more restrictive criteria


##### 50% dominance threshold
# let's just pull out the common dominant species in an easy to name df
goodspp <- spp.plot.traits %>% filter(RGR>0) %>% group_by(SP.ID) %>% summarise(nplots= length(unique(PLOT_ID)))
# remove species with fewer than 3 plots
plotavs <- spp.plot.traits%>% subset(SP.ID %in% goodspp$SP.ID[which(goodspp$nplots>2)] & RGR>0)
plotavs$SP.ID <- factor(plotavs$SP.ID)
# final sample size.
# ABICON ABIGRA ABIPRO JUNOCC PINCON PINJEF PINPON PSEMEN TSUHET 
# 16     12      3      6      6      7     36     60      5 

# just to make it look better, I'm going to move the common species to the front of the queue
levels(plotavs$SP.ID) <- list(PSEMEN="PSEMEN", PINPON="PINPON",ABICON="ABICON",ABIGRA="ABIGRA",PINJEF="PINJEF",PINCON="PINCON", JUNOCC="JUNOCC", TSUHET="TSUHET", ABIPRO='ABIPRO')


#### 90% dominance threshold
# let's just pull out the common dominant species in an easy to name df
goodspp90 <- spp.plot.traits %>% filter(RGR90>0) %>% group_by(SP.ID) %>% summarise(nplots= length(unique(PLOT_ID)))
# remove species with fewer than 3 plots, as of 06.22.17 bumped up to 5 plots
plotavs90 <- spp.plot.traits %>% subset(SP.ID %in% goodspp90$SP.ID[which(goodspp90$nplots>5)] & RGR90>0)
plotavs90$SP.ID <- factor(plotavs90$SP.ID)
# final sample size.
  # if dominance >80% stand BA
# ABIAMA ABICON ABIGRA ABIMAG ABIPRO ALNRUB JUNOCC PINCON PINJEF PINPON PSEMEN QUECHR TSUHET 
# 1     11      5      1      3      1      6      6      4     31     35      1      1 
  # if dominance >85% stand BA, lose ABIGRA, ABICON
# ABICON ABIGRA ABIMAG ABIPRO ALNRUB JUNOCC PINCON PINJEF PINPON PSEMEN QUECHR TSUHET 
# 10      2      1      3      1      6      6      4     29     33      1      1 
  # if dominance >90% stand BA, lose ABIGRA, ABICON.
# ABICON ABIGRA ABIPRO ALNRUB JUNOCC PINCON PINJEF PINPON PSEMEN QUECHR TSUHET 
# 9      2      3      1      6      6      4     29     31      1      1 
  #--> Don't actually lose that much from 85% to 90%, and 90% sounds like a pretty strict cuttoff

# just to make it look better, I'm going to move the common species to the front of the queue
levels(plotavs90$SP.ID) <- list(PSEMEN="PSEMEN", PINPON="PINPON",ABICON="ABICON", JUNOCC="JUNOCC",PINCON="PINCON")#,PINJEF="PINJEF", ABIPRO='ABIPRO')


domSPPavs <- spp.plot.traits %>% subset(SP.ID %in% goodspp90$SP.ID & RGR90 > 0 ) %>% summarise(RGR=mean(RGR90, na.rm=T), stGrowth=mean(stGrowth90, na.rm=T)
                     ,  LLmonths=mean(mLLmonths, na.rm=T), LMA = mean(mLMA, na.rm=T), Nmass = mean(mNITROGEN, na.rm=T), Narea= mean(mNarea, na.rm=T)
                     , climPC1 = mean(climPC1), climPC2 = mean(climPC2), soil_N = mean(soil_N,na.rm=T), soil_pH = mean(soil_pH), ASA= mean(ASA, na.rm=T))
domSPPavs$log.LL <- log(domSPPavs$LLmonths, base=10)
domSPPavs$log.LMA <- log(domSPPavs$LMA, base=10)
domSPPavs$log.Nmass <- log(domSPPavs$Nmass, base=10)
domSPPavs$log.Narea <- log(domSPPavs$Narea, base=10)

#### Analyses with 50% dominance ########################################


### significant relationship between mlog(RGR) and mlog LL
# haven't rerun logs since updated plotavs
modLLvRGR <- lmer(RGR~mlog.LL + (mlog.LL|SP.ID), plotavs)
modLLvRGRnull <- lmer(RGR~1 + (mlog.LL|SP.ID), plotavs)
anova(modLLvRGRnull, modLLvRGR) #p= 0.0007 / log = 0.01146      IGNORE -> p = 0.0058 all traits.common, p=0.032 with plotavs ,p= 0.047 with restricted plotavs, p=0.005 .mono75%
r.squaredGLMM(modLLvRGR) # marginal R2= 0.26, conditional = 0.4457
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

modNmassvRGR <- lmer(RGR~mlog.Nmass + (mlog.Nmass|SP.ID), plotavs, REML=F)
modNmassvRGRnull <-lmer(RGR~1 + (scale(mlog.Nmass)|SP.ID), plotavs, REML=F)
anova(modNmassvRGRnull, modNmassvRGR) #p=0.03327   / log=0.02596   
r.squaredGLMM(modNmassvRGR) # marginal R2=0.058, conditional = 0.4017


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

modNmassvGrowth <- lmer(stGrowth~mlog.Nmass + (mlog.Nmass|SP.ID), plotavs, REML=F)
modNmassvGrowthnull <-lmer(stGrowth~1 + (mlog.Nmass|SP.ID), plotavs, REML=F)
anova(modNmassvGrowthnull, modNmassvGrowth) #p=0.1009 


### Test whether traits or env better predict RGR ###
tmplotavs <- plotavs[with(plotavs, which(!is.na(mlog.Nmass) & !is.na(mlog.Narea) & !is.na(mlog.LL) & !is.na(mlog.LMA) & !is.na(soil_N) & !is.na(ASA))),]
traitsmod <- lmer(RGR ~ mlog.Nmass + mlog.LL + mlog.LMA + mlog.Narea + (1|SP.ID), tmplotavs)
enviromod <- lmer(RGR ~ climPC1 + climPC2 + soil_N + log(ASA) + (1|SP.ID), tmplotavs)
# things are even worse with stGrowth
traitsmod <- lmer(stGrowth ~ mlog.Nmass + mlog.LL + mlog.LMA + mlog.Narea + (1|SP.ID), tmplotavs)
# R2m = 0.173
enviromod <- lmer(stGrowth ~ climPC1 + climPC2 + soil_N + log(ASA) + (1|SP.ID), tmplotavs)
# R2m = 0.36

traitsmod <- lmer(log(RGR) ~ mlog.Nmass + mlog.LL + mlog.LMA + mlog.Narea + (1|SP.ID), tmplotavs)
enviromod <- lmer(log(RGR) ~ climPC1 + climPC2 + soil_N + log(ASA) + (1|SP.ID), tmplotavs)
r.squaredGLMM(traitsmod) # raw R2m = .418, loged r2m=41
r.squaredGLMM(enviromod) # raw R2m = .653, logged R2m =.73

####### Figure 6: reploted with traits.mono and biost_growth instead of RGR.
# for 6 panel figure without Nmass, but with rGR top row and biostGRWOTH bottom
quartz(width=7.1, height=6) # Eco Let width is 173mm, 110mm or 82mm
#pdf(file = "Traits-v-Growth_v2.pdf",width = 7.1, height=6)
par(mar=c(2.5,2,1.2,1), oma=c(2,2,0,0),mfrow=c(2,4), cex=1.1,mgp=c(2.5,1,0))

# for 4 panel figure with only RGR, but Nmass added
quartz(width=4.3, height=5)
pdf(file="Traits-v-Growth_v2_RGRonly.pdf", width=4.3, height=5)
par(mar=c(2.5,2,1.2,1), oma=c(2,2,0,0),mfrow=c(2,2), cex=1,mgp=c(2.5,.5,0))

#### top three panels: RGR
f0LL <- predict(modLLvRGR, re.form=NA)
f1LL <- fitted(modLLvRGR)
# sort lma values so lines draw right
I <- order(plotavs$mlog.LL[-which(is.na(plotavs$RGR) | is.na(plotavs$mlog.LL))])
LLs <- sort(plotavs$mlog.LL[-which(is.na(plotavs$RGR) | is.na(plotavs$mlog.LL))])
spid <- plotavs$SP.ID[-which(is.na(plotavs$RGR) | is.na(plotavs$mlog.LL))][I]
plot(RGR~mlog.LL, plotavs, col=SP.ID, pch=16, cex=.6, xlab="log(Leaf Lifespan)")
mtext("Relative Growth Rate", side=2, line=2, cex=1.1)
mtext("log(Leaf Life)", side=1, line=1.5, cex=1.1)

for(i in levels(spid)){
  lines(LLs[which(spid==i)], f1LL[I][which(spid==i)],  col=spid[which(spid==i)], lwd=2)
}
lines(LLs, f0LL[I], lwd=4, lty=1)
mtext(text = "p=0.0007", side=3, adj = .8, line=-1, font=1 )
mtext(text = expression(paste(R[m]^2,"=0.26", sep="")), side=3, adj = .8, line=-2.2 )
mtext("a)", side=3, adj=0)

### plotting LMA v RGR
f0lma <- predict(modLMAvRGR, re.form=NA)
f1lma <- fitted(modLMAvRGR)
# sort lma values so lines draw right # no NA values in LMA
I <- order(plotavs$mlog.LMA)
LMAs <- sort(plotavs$mlog.LMA)
spid <- plotavs$SP.ID[I]
plot(RGR~mlog.LMA, plotavs, col=SP.ID, pch=16, cex=.6, xlab="log(LMA)")#log(LMA)")
mtext("log(LMA)", side=1, line=1.5, cex=1.1)

for(i in levels(spid)){
  lines(LMAs[which(spid==i)], f1lma[I][which(spid==i)],  col=spid[which(spid==i)], lwd=2)
}
lines(c(min(LMAs), max(LMAs)), c(min(f0lma), max(f0lma)), lwd=4, lty=3)
mtext(text = "p=0.55", side=3, adj = .9, line=-1 )
mtext("b)", side=3, adj=0)

### plotting Narea v RGR
f0Narea <- predict(modNareavRGR, re.form=NA)
f1Narea <- fitted(modNareavRGR)
# sort lma values so lines draw right
I <- order(plotavs$mlog.Narea[-which(is.na(plotavs$RGR) | is.na(plotavs$mlog.Narea))])
Nareas <- sort(plotavs$mlog.Narea[-which(is.na(plotavs$RGR) | is.na(plotavs$mlog.Narea))])
spid <- plotavs$SP.ID[-which(is.na(plotavs$RGR) | is.na(plotavs$mlog.Narea))][I]
plot(RGR~mlog.Narea, plotavs, col=SP.ID, pch=16, cex=.6, xlab="")#log(Narea)")
mtext("Relative Growth Rate", side=2, line=2, cex=1.1)
mtext("log(Narea)", side=1, line=1.5, cex=1.1)

for(i in levels(spid)){
  lines(Nareas[which(spid==i)], f1Narea[I][which(spid==i)],  col=spid[which(spid==i)], lwd=2)
}
lines(Nareas, f0Narea[I], lwd=4, lty=3)
mtext(text = "p=0.18", side=3, adj = .1, line=-1 )
mtext("c)", side=3, adj=0)

### plotting Nmass v RGR
f0Nmass <- predict(modNmassvRGR, re.form=NA)
f1Nmass <- fitted(modNmassvRGR)
# sort lma values so lines draw right
I <- order(plotavs$mlog.Nmass[-which(is.na(plotavs$RGR) | is.na(plotavs$mlog.Nmass))])
Nmasss <- sort(plotavs$mlog.Nmass[-which(is.na(plotavs$RGR) | is.na(plotavs$mlog.Nmass))])
spid <- plotavs$SP.ID[-which(is.na(plotavs$RGR) | is.na(plotavs$mlog.Nmass))][I]
plot(RGR~mlog.Nmass, plotavs, col=SP.ID, pch=16, cex=.6, xlab="", ylim=c(0,.37), xlim=c(-.22,.25))#log(Nmass)")
mtext("log(Nmass)", side=1, line=1.5, cex=1.1)

for(i in levels(spid)){
  lines(Nmasss[which(spid==i)], f1Nmass[I][which(spid==i)],  col=spid[which(spid==i)], lwd=2)
}
lines(Nmasss, f0Nmass[I], lwd=4, lty=1)
mtext(text = "p=0.033", side=3, adj = .01, line=-1 )
mtext(text = expression(paste(R[m]^2,"=0.06", sep="")), side=3, adj = .01, line=-2.2 )
mtext("d)", side=3, adj=0)
# 

###### bio ST Growth ######
### plotting Leaf Lifespan v biomass Standardized growth
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
mtext("e)", side=3, adj=0)
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
mtext("f)", side=3, adj=0)
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
mtext(text = "p=0.10", side=3, adj = .2, line=-1 )
mtext("g)", side=3, adj=0)
mtext("log(Narea)", side=1, line=2.5, cex=1.1)


### plotting Nmass v stGrowth
f0Nmass <- predict(modNmassvGrowth, re.form=NA)
f1Nmass <- fitted(modNmassvGrowth)
# sort lma values so lines draw right
I <- order(plotavs$mlog.Nmass[-which(is.na(plotavs$stGrowth) | is.na(plotavs$mlog.Nmass))])
Nmasss <- sort(plotavs$mlog.Nmass[-which(is.na(plotavs$stGrowth) | is.na(plotavs$mlog.Nmass))])
spid <- plotavs$SP.ID[-which(is.na(plotavs$stGrowth) | is.na(plotavs$mlog.Nmass))][I]
plot(stGrowth~mlog.Nmass, plotavs, col=SP.ID, pch=16, cex=.4, xlab="log(Nmass)")
for(i in levels(spid)){
  lines(Nmasss[which(spid==i)], f1Nmass[I][which(spid==i)],  col=spid[which(spid==i)], lwd=2)
}
lines(Nmasss, f0Nmass[I], lwd=4, lty=3)
mtext(text = "p=0.19", side=3, adj = .2, line=-1 )
mtext("h)", side=3, adj=0)
mtext("log(Nmass)", side=1, line=2.5, cex=1.1)







############# ***Analyses with 90% Dominance #####################


### significant relationship between mlog(RGR, base=10) and mlog LL
# haven't rerun logs since updated plotavs909
options(na.action= na.omit)
modLLvRGR <- lmer(log(RGR, base=10)~mlog.LL + (mlog.LL|SP.ID), plotavs90)
modLLvRGRnull <- lmer(log(RGR, base=10)~1 + (mlog.LL|SP.ID), plotavs90)
anova(modLLvRGRnull, modLLvRGR) #p= 0.09138 / log = 0.0124   
r.squaredGLMM(modLLvRGR) # marginal R2= 0.298, conditional = 0.390 with logged RGR
# using an exponential link but a gaussian error distribution
# modLLvRGR <- glmer(RGR~scale(mlog.LL) + (scale(mlog.LL)|SP.ID), plotavs90,family=(gaussian(link='log')))
# modLLvRGRnull <- glmer(RGR~1 + (mlog.LL|SP.ID), plotavs90,family=(gaussian(link='log')))
# anova(modLLvRGRnull, modLLvRGR) # 0.081     


modLMAvRGR <- lmer(log(RGR, base=10)~mlog.LMA + (mlog.LMA|SP.ID), plotavs90, REML=F)
modLMAvRGRnull <- lmer(log(RGR, base=10)~1 + (mlog.LMA|SP.ID), plotavs90, REML=F)
anova(modLMAvRGRnull, modLMAvRGR) #p=0.2616/ log = 0.098    
r.squaredGLMM(modLMAvRGR) #m 0.1095632 c 0.1375287 

modNareavRGR <- lmer(log(RGR, base=10)~mlog.Narea + (mlog.Narea|SP.ID), plotavs90, REML=F)
modNareavRGRnull <-lmer(log(RGR, base=10)~1 + (scale(mlog.Narea)|SP.ID), plotavs90, REML=F)
anova(modNareavRGRnull, modNareavRGR) #p=0.6595   / log = 0.6838   

modNmassvRGR <- lmer(log(RGR, base=10)~mlog.Nmass + (mlog.Nmass|SP.ID), plotavs90, REML=F)
modNmassvRGRnull <-lmer(log(RGR, base=10)~1 + (scale(mlog.Nmass)|SP.ID), plotavs90, REML=F)
anova(modNmassvRGRnull, modNmassvRGR) #p=0.1278   / log=0.0351  
r.squaredGLMM(modNmassvRGR) # marginal R2=0.09963567, conditional = 0.21170685


## with different biomass...
modLLvGrowth <- lmer(stGrowth~mlog.LL + (mlog.LL|SP.ID), plotavs90, REML=F)
modLLvGrowthnull <- lmer(stGrowth~1 + (mlog.LL|SP.ID), plotavs90, REML=F)
anova(modLLvGrowthnull, modLLvGrowth) # p=0.156    

modLMAvGrowth <- lmer(stGrowth~mlog.LMA + (mlog.LMA|SP.ID), plotavs90, REML=F)
modLMAvGrowthnull <- lmer(stGrowth~1 + (mlog.LMA|SP.ID), plotavs90, REML=F)
anova(modLMAvGrowthnull, modLMAvGrowth) #p=0.1144

modNareavGrowth <- lmer(stGrowth~mlog.Narea + (mlog.Narea|SP.ID), plotavs90, REML=F)
modNareavGrowthnull <-lmer(stGrowth~1 + (mlog.Narea|SP.ID), plotavs90, REML=F)
anova(modNareavGrowthnull, modNareavGrowth) #p=0.8852

modNmassvGrowth <- lmer(stGrowth~mlog.Nmass + (mlog.Nmass|SP.ID), plotavs90, REML=F)
modNmassvGrowthnull <-lmer(stGrowth~1 + (mlog.Nmass|SP.ID), plotavs90, REML=F)
anova(modNmassvGrowthnull, modNmassvGrowth) #p=0.2458


### Test whether traits or env better predict RGR ###
  # make a dataset without any NAs for model selection
tmplotavs90 <- plotavs90 %>% filter(complete.cases(mlog.Nmass, mlog.Narea, mlog.LL, mlog.LMA, soil_N, ASA))#[with(plotavs90, which(!is.na(mlog.Nmass) & !is.na(mlog.Narea) & !is.na(mlog.LL) & !is.na(mlog.LMA) & !is.na(soil_N) & !is.na(ASA) & which(plotavs90$SP.ID %in% c("PINJEF","ABIPRO")),])),]
traitsmod <- lmer(RGR ~  mlog.LL + mlog.LMA + mlog.Nmass + (1|SP.ID), tmplotavs90) # note: mlog.Narea is a linear combo of LMA and Nmass, so it actually gets dropped no matter what
r.squaredGLMM(traitsmod) # R2m = .404, R2c = 0.484
enviromod <- lmer(RGR ~ climPC1 + climPC2 + soil_N + log(ASA) + (1|SP.ID), tmplotavs90)
r.squaredGLMM(enviromod) #R2m =0.737, R2c= 0.743

# things are even worse with stGrowth
traitsmod <- lmer(stGrowth ~ mlog.Nmass + mlog.LL + mlog.LMA + mlog.Narea + (1|SP.ID), tmplotavs90)
# R2m = 0.169
enviromod <- lmer(stGrowth ~ climPC1 + climPC2 + soil_N + log(ASA) + (1|SP.ID), tmplotavs90)
# R2m = 0.538

traitsmod <- lmer(log(RGR, base=10) ~ mlog.Nmass + mlog.LL + mlog.LMA + mlog.Narea + (1|SP.ID), tmplotavs90)
enviromod <- lmer(log(RGR, base=10) ~ climPC1 + climPC2 + soil_N + log(ASA) + (1|SP.ID), tmplotavs90)
r.squaredGLMM(traitsmod) # raw R2m = .404, loged r2m=0.422
r.squaredGLMM(enviromod) # raw R2m = 0.7425527, logged R2m =0.782
  



r.squaredGLMM(lmer(log(RGR, base=10)~log(ASA) + log(LeafFrac) + (1|SP.ID), tmplotavs90[-which(tmplotavs90$PLOT_ID %in% c(43,55)),]))
r.squaredGLMM(lmer(log(RGR, base=10)~log(LeafFrac) + (1|SP.ID), tmplotavs90[-which(tmplotavs90$PLOT_ID %in% c(43,55)),]))
r.squaredGLMM(lmer(log(RGR, base=10)~LAI + (1|SP.ID), tmplotavs90))
r.squaredGLMM(lmer(log(RGR, base=10)~I(NPP * LeafFrac/mLMA) + (1|SP.ID), tmplotavs90))
r.squaredGLMM(lmer(log(RGR, base=10)~log(ASA) + log(Height) + (1|SP.ID), tmplotavs90))


####### Figure 6: reploted with traits.mono and biost_growth instead of RGR.
# for 6 panel figure without Nmass, but with rGR top row and biostGRWOTH bottom
quartz(width=6.8, height=6) # Eco Let width is 173mm, 110mm or 82mm
#pdf(file = "Traits-v-Growth_v2.pdf",width = 7.1, height=6)
par(mar=c(0,0,0,0), oma=c(4.5,4,1.2,1),mfrow=c(2,4), cex=1.1,mgp=c(2.5,1,0))

# for 4 panel figure with only RGR, but Nmass added
quartz(width=4.3, height=5)
#pdf(file="Traits-v-Growth-90Dom_v1_RGRonly.pdf", width=4.3, height=5)
par(mar=c(2.5,2,1.2,1), oma=c(2,2,0,0),mfrow=c(2,2), cex=1,mgp=c(2.5,.5,0))

#### top three panels: RGR
f0LL <- predict(modLLvRGR, re.form=NA)
f1LL <- fitted(modLLvRGR)
# sort lma values so lines draw right
I <- order(plotavs90$mlog.LL[-which(is.na(plotavs90$RGR) | is.na(plotavs90$mlog.LL))])
LLs <- sort(plotavs90$mlog.LL[-which(is.na(plotavs90$RGR) | is.na(plotavs90$mlog.LL))])
spid <- plotavs90$SP.ID[-which(is.na(plotavs90$RGR) | is.na(plotavs90$mlog.LL))][I]
plot(log(RGR, base=10)~mlog.LL, plotavs90, col=SP.ID, pch=16, cex=.6,xaxt="n", xlab="") #xlab="log(Leaf Lifespan)")
mtext("Relative Growth Rate", side=2, line=2, cex=1.1)
#mtext("log(Leaf Life)", side=1, line=1.5, cex=1.1)

for(i in levels(spid)){
  lines(LLs[which(spid==i)], f1LL[I][which(spid==i)],  col=spid[which(spid==i)], lwd=2)
}
lines(LLs, f0LL[I], lwd=4, lty=1)
mtext(text = paste0("p=",round(anova(modLLvRGR, modLLvRGRnull)$`Pr(>Chisq)`[2],3)), side=3, adj = .9, line=-1, font=1 )
r2tmp <- round(r.squaredGLMM(modLLvRGR)[1], 2)
mtext(text = bquote(~R[m]^2==.(r2tmp)), side=3, adj = .9, line=-2.2 )
mtext("a)", side=3, adj=0)

### plotting LMA v RGR
f0lma <- predict(modLMAvRGR, re.form=NA)
f1lma <- fitted(modLMAvRGR)
# sort lma values so lines draw right # no NA values in LMA
I <- order(plotavs90$mlog.LMA)
LMAs <- sort(plotavs90$mlog.LMA)
spid <- plotavs90$SP.ID[I]
plot(log(RGR, base=10)~mlog.LMA, plotavs90, col=SP.ID, pch=16, cex=5*LeafFrac,xaxt="n", yaxt='n') #xlab="log(LMA)")#, xlim=c(2,2.6))#log(LMA)")
#mtext("log(LMA)", side=1, line=1.5, cex=1.1)

for(i in levels(spid)){
  lines(LMAs[which(spid==i)], f1lma[I][which(spid==i)],  col=spid[which(spid==i)], lwd=2)
}
lines(c(min(LMAs), max(LMAs)), c(max(f0lma), min(f0lma)), lwd=4, lty=3)
mtext(text = paste0("p=",round(anova(modLMAvRGR, modLMAvRGRnull)$`Pr(>Chisq)`[2],3)), side=3, adj = .05, line=-1 )
mtext("b)", side=3, adj=0)

### plotting Nmass v RGR
f0Nmass <- predict(modNmassvRGR, re.form=NA)
f1Nmass <- fitted(modNmassvRGR)
# sort lma values so lines draw right
I <- order(plotavs90$mlog.Nmass[which(!is.na(plotavs90$RGR) & !is.na(plotavs90$mlog.Nmass))])
Nmasss <- sort(plotavs90$mlog.Nmass[which(!is.na(plotavs90$RGR) & !is.na(plotavs90$mlog.Nmass))])
spid <- plotavs90$SP.ID[which(!is.na(plotavs90$RGR) & !is.na(plotavs90$mlog.Nmass))][I]
plot(log(RGR, base=10)~mlog.Nmass, plotavs90, col=SP.ID, pch=16, cex=.6, xlab="", xlim=c(-.22,.25), xaxt="n", yaxt="n")#log(Nmass)")
#mtext("log(Nmass)", side=1, line=1.5, cex=1.1)

for(i in levels(spid)){
  lines(Nmasss[which(spid==i)], f1Nmass[I][which(spid==i)],  col=spid[which(spid==i)], lwd=2)
}
lines(Nmasss, f0Nmass[I], lwd=4, lty=1)
mtext(text = paste0("p=",round(anova(modNmassvRGR, modNmassvRGRnull)$`Pr(>Chisq)`[2],3)), side=3, adj = .01, line=-1 )
r2tmp <- round(r.squaredGLMM(modNmassvRGR)[1], 2)
mtext(text = bquote(~R[m]^2==.(r2tmp)), side=3, adj = .01, line=-2.2 )
mtext("c)", side=3, adj=0)
# 

### plotting Narea v RGR
f0Narea <- predict(modNareavRGR, re.form=NA)
f1Narea <- fitted(modNareavRGR)
# sort lma values so lines draw right
I <- order(plotavs90$mlog.Narea[which(!is.na(plotavs90$RGR) & !is.na(plotavs90$mlog.Narea))])
Nareas <- sort(plotavs90$mlog.Narea[which(!is.na(plotavs90$RGR) & !is.na(plotavs90$mlog.Narea))])
spid <- plotavs90$SP.ID[which(!is.na(plotavs90$RGR) & !is.na(plotavs90$mlog.Narea))][I]
plot(log(RGR, base=10)~mlog.Narea, plotavs90, col=SP.ID, pch=16, cex=.6, xlab="", xaxt="n", yaxt="n")#log(Narea)")
#mtext("Relative Growth Rate", side=2, line=2, cex=1.1)
#mtext("log(Narea)", side=1, line=1.5, cex=1.1)

for(i in levels(spid)){
  lines(Nareas[which(spid==i)], f1Narea[I][which(spid==i)],  col=spid[which(spid==i)], lwd=2)
}
lines(Nareas, f0Narea[I], lwd=4, lty=3)
mtext(text = paste0("p=",round(anova(modNareavRGR, modNareavRGRnull)$`Pr(>Chisq)`[2],3)), side=3, adj = .05, line=-1 )
mtext("d)", side=3, adj=0)





############# ****CWtrait Growth~Trait analysis with infilling #############
# as of 06.21, I replaced TRAITp with TRAITp_if because infilling's important
# old p values are without infilling, if p values are the best in class


LLmod <- lm(RGR~log.cw_LLp_if, biomass)
LLmod2 <- lm(log(RGR, base=10)~log.cw_LLp_if, biomass)
#LLmod2 <- lm(RGR~cw_log.LLp_if, biomass)
summary(LLmod)
summary(LLmod2) #if p=4.144e-08 p = 8.5e-7

LMAmod <- lm(RGR~log.cw_LMAp_if, biomass)
LMAmod2 <- lm(log(RGR, base=10)~log.cw_LMAp_if, biomass)
#LMAmod2 <- lm(RGR~cw_log.LMAp_if, biomass)
summary(LMAmod)
summary(LMAmod2) #if p=0.3699, p=0.84


Nmassmod <- lm(RGR~log.cw_Nmassp_if, biomass)
Nmassmod2 <- lm(log(RGR, base=10)~log.cw_Nmassp_if, biomass)
#Nmassmod2 <- lm(RGR~cw_log.Nmassp_if, biomass)
summary(Nmassmod)
summary(Nmassmod2) # if p=0.002344, old p=0.019

Nareamod <- lm(RGR~log.cw_Nareap_if, biomass)
Nareamod2 <- lm(log(RGR, base=10)~log.cw_Nareap_if, biomass)
#Nareamod2 <- lm(RGR~cw_log.Nareap_if, biomass)
summary(Nareamod)
summary(Nareamod2) #if p= p-value: 0.3675, p=0.96

## on second thought, I think I should use average and then log, so I don't strangely weight things
# so I should use all the log.cw variables
# also definitely need to log(RGR, base=10) in order to normalize it.




############# ****SPPtrait Growth~Trait analysis #############
# run this with the mean trait and RGR values from only those species that ever achieve 90% dominance (and only the stands where they are 90% dominant)

LLmodSPP <- lm(RGR~log.LL, domSPPavs)
LLmodSPP2 <- lm(log(RGR, base=10)~log.LL, domSPPavs)
#LLmodSPP2 <- lm(RGR~cw_log.LLp_if, biomass)
summary(LLmodSPP) #p= 0.145
summary(LLmodSPP2) #p=0.03942

LMAmodSPP <- lm(RGR~log.LMA, domSPPavs)
LMAmodSPP2 <- lm(log(RGR, base=10)~log.LMA, domSPPavs)
#LMAmodSPP2 <- lm(RGR~cw_log.LMAp_if, biomass)
summary(LMAmodSPP) # p=0.48
summary(LMAmodSPP2) #if p=0.2948


NmassmodSPP <- lm(RGR~log.Nmass, domSPPavs)
NmassmodSPP2 <- lm(log(RGR, base=10)~log.Nmass, domSPPavs)
#NmassmodSPP2 <- lm(RGR~cw_log.Nmassp_if, domSPPavs)
summary(NmassmodSPP) # p=0.2193
summary(NmassmodSPP2) # if p=0.1121

NareamodSPP <- lm(RGR~log.Narea, domSPPavs)
NareamodSPP2 <- lm(log(RGR, base=10)~log.Narea, domSPPavs)
#NareamodSPP2 <- lm(RGR~cw_log.Nareap_if, domSPPavs)
summary(NareamodSPP) #p=0.5303
summary(NareamodSPP2) #if p= p-value: 0.6825



### SPP: Test whether traits or env better predict RGR ###
# make a dataset without any NAs for model selection
tmdomSPPavs <- domSPPavs %>% filter(complete.cases(log.Nmass, log.LL, log.LMA, soil_N, ASA))#[with(domSPPavs, which(!is.na(mlog.Nmass) & !is.na(mlog.Narea) & !is.na(mlog.LL) & !is.na(mlog.LMA) & !is.na(soil_N) & !is.na(ASA) & which(domSPPavs$SP.ID %in% c("PINJEF","ABIPRO")),])),]
traitsmod <- lm(RGR ~  log.LL + log.LMA + log.Nmass , tmdomSPPavs) # note: mlog.Narea is a linear combo of LMA and Nmass, so it actually gets dropped no matter what
traitsmod <- lm(RGR ~  log.LL  , tmdomSPPavs) # note: mlog.Narea is a linear combo of LMA and Nmass, so it actually gets dropped no matter what
summary(traitsmod) # R2 = .24
enviromod <- lm(RGR ~ climPC1 + climPC2 + soil_N + log(ASA), tmdomSPPavs)
enviromod <- lm(RGR ~  log(ASA), tmdomSPPavs)
summary(enviromod) # r2=0.77

# things are even worse with stGrowth
traitsmod <- lmer(stGrowth ~ mlog.Nmass + mlog.LL + mlog.LMA + mlog.Narea + (1|SP.ID), tmdomSPPavs)
# R2m = 0.169
enviromod <- lmer(stGrowth ~ climPC1 + climPC2 + soil_N + log(ASA) + (1|SP.ID), tmdomSPPavs)
# R2m = 0.538

traitsmod <- lmer(log(RGR, base=10) ~ mlog.Nmass + mlog.LL + mlog.LMA + mlog.Narea + (1|SP.ID), tmdomSPPavs)
enviromod <- lmer(log(RGR, base=10) ~ climPC1 + climPC2 + soil_N + log(ASA) + (1|SP.ID), tmdomSPPavs)
r.squaredGLMM(traitsmod) # raw R2m = .404, loged r2m=0.422
r.squaredGLMM(enviromod) # raw R2m = 0.7425527, logged R2m =0.782










############## New Figure ##########

quartz(width=4.3, height=5)
#pdf(file="Traits-v-Growth-90Dom_v1_RGRonly.pdf", width=4.3, height=5)
par(mar=c(2.5,2,1.2,1), oma=c(2,2,0,0),mfrow=c(2,2), cex=1,mgp=c(2.5,.5,0))

plot(log(RGR, base=10)~log.cw_LLp_if, biomass, pch=16)
abline(lm(log(RGR, base=10)~log.cw_LLp_if, biomass), lwd=2)
mtext(text = paste0("p=", round(anova(LLmod2)$`Pr(>F)`[1],8), ","), side=1, adj=.1, line=-1, cex=.9)
mtext(text = bquote(~R^2==.(round(summary(LLmod2)$r.squared,2))), side=1, adj = .95, line=-1 , cex=.9)
mtext("log(RGR, base=10)", side=2, line=2, cex=1.1)
mtext("log(cwm Leaf Life)", side=1, line=1.5, cex=1.1)
mtext("a)", side=3, adj=0)

plot(log(RGR, base=10)~log.cw_LMAp_if, biomass, pch=16)
abline(lm(log(RGR, base=10)~log.cw_LMAp_if, biomass), lwd=2, lty=2)
mtext(text = paste0("p=", round(anova(LMAmod2)$`Pr(>F)`[1],2)), side=1, adj=.05, line=-1, cex=.9)
mtext("log(cwm LMA)", side=1, line=1.5, cex=1.1)
mtext("b)", side=3, adj=0)

plot(log(RGR, base=10)~log.cw_Nmassp_if, biomass, pch=16)
abline(lm(log(RGR, base=10)~log.cw_Nmassp_if, biomass), lwd=2, lty=1)
mtext(text = paste0("p=", round(anova(Nmassmod2)$`Pr(>F)`[1],3),","), side=1, adj=.05, line=-1, cex=.9)
mtext(text = bquote(~R^2==.(round(summary(Nmassmod2)$r.squared,2))), side=1, adj = .95, line=-1, cex=.9 )
mtext("log(RGR, base=10)", side=2, line=2, cex=1.1)
mtext("log(cwm Nmass)", side=1, line=1.5, cex=1.1)
mtext("c)", side=3, adj=0)

plot(log(RGR, base=10)~log.cw_Nareap_if, biomass, pch=16)#, yaxt="n")
#axis(2,at = log(c(.01,.02,.03,0.04,.1,0.2,0.3)), labels =c(.01,.02,.03,0.04,.1,0.2,0.3) )
abline(lm(log(RGR, base=10)~log.cw_Nareap_if, biomass), lwd=2, lty=2)
mtext(text = paste0("p=", round(anova(Nareamod2)$`Pr(>F)`[1],2)), side=1, adj=.05, line=-1)
mtext("log(cwm Narea)", side=1, line=1.5, cex=1.1)
mtext("d)", side=3, adj=0)





########### Analysis with stGrowth ##########


LLmod.g <- lm(BIOstGROWTHgam~log.cw_LLp_if, biomass)
#LLmod2.g <- lm(log(BIOstGROWTHgam)~log.cw_LLp_if, biomass)
#LLmod2 <- lm(BIOstGROWTHgam~cw_log.LLp_if, biomass)
summary(LLmod.g) #if p=2.041e-05, old p = 7.073e-05, r2 = 0.09
#summary(LLmod2.g)

LMAmod.g <- lm(BIOstGROWTHgam~log.cw_LMAp_if, biomass)
#LMAmod2 <- lm(log(BIOstGROWTHgam)~log.cw_LMAp_if, biomass)
#LMAmod2 <- lm(BIOstGROWTHgam~cw_log.LMAp_if, biomass)
summary(LMAmod.g) # if p=p-value: 0.0001929, old p=0.00131, R2 = 0.054
#summary(LMAmod2)


Nmassmod.g <- lm(BIOstGROWTHgam~log.cw_Nmassp_if, biomass)
# Nmassmod2 <- lm(log(BIOstGROWTHgam)~log.cw_Nmassp_if, biomass)
#Nmassmod2 <- lm(BIOstGROWTHgam~cw_log.Nmassp_if, biomass)
summary(Nmassmod.g)#if p-value: 0.0002197, p=0.007553, r2 = 0.03898
#summary(Nmassmod2) 

Nareamod.g <- lm(BIOstGROWTHgam~log.cw_Nareap_if, biomass)
# Nareamod2 <- lm(log(BIOstGROWTHgam)~log.cw_Nareap_if, biomass)
#Nareamod2 <- lm(BIOstGROWTHgam~cw_log.Nareap_if, biomass)
summary(Nareamod.g) #if p-value: 0.06943, old p=0.04692, r2 - 0.02523
#summary(Nareamod2)






quartz(width=4.3, height=5)
#pdf(file="Traits-v-Growth-90Dom_v1_RGRonly.pdf", width=4.3, height=5)
par(mar=c(2.5,2,1.2,1), oma=c(2,2,0,0),mfrow=c(2,2), cex=1,mgp=c(2.5,.5,0))

plot(BIOstGROWTHgam~log.cw_LLp_if, biomass, pch=16, ylim=c(-400,850))
abline(lm(BIOstGROWTHgam~log.cw_LLp_if, biomass), lwd=2)
mtext(text = paste0("p<", round(anova(LLmod.g)$`Pr(>F)`[1],6), ","), side=3, adj=.1, line=-1, cex=.9)
mtext(text = bquote(~R^2==.(round(summary(LLmod.g)$r.squared,2))), side=3, adj = .95, line=-1 , cex=.9)
mtext("Emp st. Growth", side=2, line=2, cex=1.1)
mtext("log(cwm Leaf Life)", side=1, line=1.5, cex=1.1)
mtext("a)", side=3, adj=0)

plot(BIOstGROWTHgam~log.cw_LMAp_if, biomass, pch=16, ylim=c(-400,850))
abline(lm(BIOstGROWTHgam~log.cw_LMAp_if, biomass), lwd=2, lty=1)
mtext(text = paste0("p=", round(anova(LMAmod.g)$`Pr(>F)`[1],4),","), side=3, adj=.1, line=-1, cex=.9)
mtext(text = bquote(~R^2==.(round(summary(LMAmod.g)$r.squared,2))), side=3, adj = .95, line=-1 , cex=.9)
mtext("log(cwm LMA)", side=1, line=1.5, cex=1.1)
mtext("b)", side=3, adj=0)

plot(BIOstGROWTHgam~log.cw_Nmassp_if, biomass, pch=16, ylim=c(-400,850))
abline(lm(BIOstGROWTHgam~log.cw_Nmassp_if, biomass), lwd=2, lty=1)
mtext(text = paste0("p=", round(anova(Nmassmod2)$`Pr(>F)`[1],3),","), side=3, adj=.05, line=-1, cex=.9)
mtext(text = bquote(~R^2==.(round(summary(Nmassmod.g)$r.squared,2))), side=3, adj = .95, line=-1, cex=.9 )
mtext("Emp st. Growth", side=2, line=2, cex=1.1)
mtext("log(cwm Nmass)", side=1, line=1.5, cex=1.1)
mtext("c)", side=3, adj=0)

plot(BIOstGROWTHgam~log.cw_Nareap_if, biomass, pch=16, ylim=c(-400,850))#, yaxt="n")
#axis(2,at = log(c(.01,.02,.03,0.04,.1,0.2,0.3)), labels =c(.01,.02,.03,0.04,.1,0.2,0.3) )
abline(lm(BIOstGROWTHgam~log.cw_Nareap_if, biomass), lwd=2, lty=2)
mtext(text = paste0("p=", round(anova(Nareamod.g)$`Pr(>F)`[1],2)), side=3, adj=.05, line=-1)
#mtext(text = bquote(~R^2==.(round(summary(Nareamod.g)$r.squared,2))), side=3, adj = .95, line=-1, cex=.9 )
mtext("log(cwm Narea)", side=1, line=1.5, cex=1.1)
mtext("d)", side=3, adj=0)











#____________________________________________________________________________
########### Test whether traits or env better predict RGR ###################
#____________________________________________________________________________


# make a dataset without any NAs for model selection
tmpbiomass <- biomass[with(biomass, which(!is.na(log.cw_Nmassp_if) & !is.na(log.cw_Nareap_if) & !is.na(log.cw_LLp_if) & !is.na(log.cw_LMAp_if) & !is.na(soil_N) & !is.na(ASA))),]
traitsmod <- lm(RGR ~ log.cw_Nmassp_if + log.cw_LLp_if + log.cw_LMAp_if + log.cw_Nareap_if, tmpbiomass)
summary(traitsmod) # R2= 0.1567, LL is all that really matters
enviromod <- lm(RGR ~ climPC1 + climPC2 + soil_N + log(ASA), tmpbiomass)
summary(enviromod) #R2m =0.5527, climPC1, climPC2, and log(ASA) are important

# things are even worse with stGrowth
traitsmod <- lm(BIOstGROWTHgam ~ log.cw_Nmassp_if + log.cw_LLp_if + log.cw_LMAp_if + log.cw_Nareap_if , tmpbiomass)
# R2m = 0.2087 # also only LL matters
enviromod <- lm(BIOstGROWTHgam ~ climPC1 + climPC2 + soil_N + log(ASA), tmpbiomass)
# R2m = 0.4848, climPC2 and log(ASA) matter

traitsmod <- lm(log(RGR, base=10) ~ log.cw_Nmassp_if + log.cw_LLp_if + log.cw_LMAp_if + log.cw_Nareap_if, tmpbiomass)
enviromod <- lm(log(RGR, base=10) ~ climPC1 + climPC2 + soil_N + log(LeafFrac), tmpbiomass)
summary(traitsmod) # log R2m = .1923 , only LL matters
summary(enviromod) # log R2m = 0.6424, climPC1, climPC2, and ASA matter.



summary(lm(log(RGR, base=10)~log.cw_LMAp_if, tmpbiomass[-which(tmpbiomass$PLOT_ID %in% c(43,55)),]))
summary(lm(log(RGR, base=10)~log(LeafFrac), tmpbiomass[-which(tmpbiomass$PLOT_ID %in% c(43,55)),]))

summary(lm(log(RGR, base=10)~log(ASA), tmpbiomass))















#______________________________________
###### 4 panel 2x2 for NPS poster #######
#______________________________________
quartz(width=4.3, height=5)
#pdf(file="Traits-v-Growth-90Dom_v1_RGRonly.pdf", width=4.3, height=5)
par(mar=c(2.5,0,1.2,0), oma=c(2,4,1,1),mfrow=c(2,2), cex=1,mgp=c(2.5,.5,0))

mypal <- brewer.pal(n=9, "Set1")
palette(mypal)

ylims <- c(-5.5,0)
cex.text <- .8
cex.pts <- .9
#### top three panels: RGR
f0LL <- predict(modLLvRGR, re.form=NA)
f1LL <- fitted(modLLvRGR)
# sort lma values so lines draw right
I <- order(plotavs90$mlog.LL[-which(is.na(plotavs90$RGR) | is.na(plotavs90$mlog.LL))])
LLs <- sort(plotavs90$mlog.LL[-which(is.na(plotavs90$RGR) | is.na(plotavs90$mlog.LL))])
spid <- plotavs90$SP.ID[-which(is.na(plotavs90$RGR) | is.na(plotavs90$mlog.LL))][I]
plot(log(RGR, base=10)~mlog.LL, plotavs90, col=SP.ID, pch=16, cex=cex.pts, xlab="", ylim=ylims) #xlab="log(Leaf Lifespan)")
#mtext("log(Relative Growth Rate)", side=2, line=2, cex=1.1)
mtext("log(Leaf Life)", side=1, line=1.5, cex=1.1)
#lines(LLs, f0LL[I], lwd=4, lty=1)
for(i in levels(spid)){
  lines(LLs[which(spid==i)], f1LL[I][which(spid==i)],  col=spid[which(spid==i)], lwd=3)
}
mtext(text = paste0("p=",round(anova(modLLvRGR, modLLvRGRnull)$`Pr(>Chisq)`[2],3)), side=3, adj = .9, line=-1, font=1, cex=cex.text )
r2tmp <- round(r.squaredGLMM(modLLvRGR)[1], 2)
mtext(text = bquote(~R[m]^2==.(r2tmp)), side=3, adj = .9, line=-2.2 ,cex=cex.text)
mtext("a)", side=3, adj=0)

### plotting LMA v RGR
f0lma <- predict(modLMAvRGR, re.form=NA)
f1lma <- fitted(modLMAvRGR)
# sort lma values so lines draw right # no NA values in LMA
I <- order(plotavs90$mlog.LMA)
LMAs <- sort(plotavs90$mlog.LMA)
spid <- plotavs90$SP.ID[I]
plot(log(RGR, base=10)~mlog.LMA, plotavs90, col=SP.ID, pch=16, cex=cex.pts, yaxt='n',ylim=ylims) #xlab="log(LMA)")#, xlim=c(2,2.6))#log(LMA)")
mtext("log(LMA)", side=1, line=1.5, cex=1.1)
#lines(c(min(LMAs), max(LMAs)), c(max(f0lma), min(f0lma)), lwd=4, lty=3)
for(i in levels(spid)){
  lines(LMAs[which(spid==i)], f1lma[I][which(spid==i)],  col=spid[which(spid==i)], lwd=3, lty=2)
}
mtext(text = paste0("p=",round(anova(modLMAvRGR, modLMAvRGRnull)$`Pr(>Chisq)`[2],3)), side=3, adj = .05, line=-1 ,cex=cex.text)
mtext("b)", side=3, adj=0)

### plotting Nmass v RGR
f0Nmass <- predict(modNmassvRGR, re.form=NA)
f1Nmass <- fitted(modNmassvRGR)
# sort lma values so lines draw right
I <- order(plotavs90$mlog.Nmass[which(!is.na(plotavs90$RGR) & !is.na(plotavs90$mlog.Nmass))])
Nmasss <- sort(plotavs90$mlog.Nmass[which(!is.na(plotavs90$RGR) & !is.na(plotavs90$mlog.Nmass))])
spid <- plotavs90$SP.ID[which(!is.na(plotavs90$RGR) & !is.na(plotavs90$mlog.Nmass))][I]
plot(log(RGR, base=10)~mlog.Nmass, plotavs90, col=SP.ID, pch=16, cex=cex.pts, xlab="", xlim=c(-.22,.25),ylim=ylims)#log(Nmass)")
mtext("log(Nmass)", side=1, line=1.5, cex=1.1)
#lines(Nmasss, f0Nmass[I], lwd=4, lty=1)
mtext("log(Relative Growth Rate)", side=2, line=2, cex=1.1, at=1)

for(i in levels(spid)){
  lines(Nmasss[which(spid==i)], f1Nmass[I][which(spid==i)],  col=spid[which(spid==i)], lwd=3)
}
mtext(text = paste0("p=",round(anova(modNmassvRGR, modNmassvRGRnull)$`Pr(>Chisq)`[2],3)), side=3, adj = .01, line=-1 ,cex=cex.text)
r2tmp <- round(r.squaredGLMM(modNmassvRGR)[1], 2)
mtext(text = bquote(~R[m]^2==.(r2tmp)), side=3, adj = .01, line=-2.2 ,cex=cex.text)
mtext("c)", side=3, adj=0)
# 

### plotting Narea v RGR
f0Narea <- predict(modNareavRGR, re.form=NA)
f1Narea <- fitted(modNareavRGR)
# sort lma values so lines draw right
I <- order(plotavs90$mlog.Narea[which(!is.na(plotavs90$RGR) & !is.na(plotavs90$mlog.Narea))])
Nareas <- sort(plotavs90$mlog.Narea[which(!is.na(plotavs90$RGR) & !is.na(plotavs90$mlog.Narea))])
spid <- plotavs90$SP.ID[which(!is.na(plotavs90$RGR) & !is.na(plotavs90$mlog.Narea))][I]
plot(log(RGR, base=10)~mlog.Narea, plotavs90, col=SP.ID, pch=16, cex=cex.pts, xlab="", yaxt="n",ylim=ylims)#log(Narea)")
#mtext("Relative Growth Rate", side=2, line=2, cex=1.1)
mtext("log(Narea)", side=1, line=1.5, cex=1.1)
#lines(Nareas, f0Narea[I], lwd=4, lty=3)
for(i in levels(spid)){
  lines(Nareas[which(spid==i)], f1Narea[I][which(spid==i)],  col=spid[which(spid==i)], lwd=3, lty=2)
}

mtext(text = paste0("p=",round(anova(modNareavRGR, modNareavRGRnull)$`Pr(>Chisq)`[2],3)), side=3, adj = .05, line=-1 ,cex=cex.text)
mtext("d)", side=3, adj=0)
#mtext("w/in spp", side=4, line=0)






#______________________________________
###### 4 panel 1x4 for NPS poster #######
#______________________________________
quartz(width=6.8, height=2.75)
#pdf(file="Traits-v-Growth-90Dom_v1_RGRonly.pdf", width=4.3, height=5)
par(mar=c(2.5,0,1.2,0), oma=c(2,4,1,1),mfrow=c(1,4), cex=1,mgp=c(2.5,.5,0))

mypal <- brewer.pal(n=9, "Set1")
palette(mypal)

ylims <- c(-5.5,0)
cex.text <- .8
cex.pts <- .9
#### top three panels: RGR
f0LL <- predict(modLLvRGR, re.form=NA)
f1LL <- fitted(modLLvRGR)
# sort lma values so lines draw right
I <- order(plotavs90$mlog.LL[-which(is.na(plotavs90$RGR) | is.na(plotavs90$mlog.LL))])
LLs <- sort(plotavs90$mlog.LL[-which(is.na(plotavs90$RGR) | is.na(plotavs90$mlog.LL))])
spid <- plotavs90$SP.ID[-which(is.na(plotavs90$RGR) | is.na(plotavs90$mlog.LL))][I]
plot(log(RGR, base=10)~mlog.LL, plotavs90, col=SP.ID, pch=16, cex=cex.pts, xlab="", ylim=ylims) #xlab="log(Leaf Lifespan)")
mtext("log(Relative Growth Rate)", side=2, line=2, cex=1.1)
mtext("log(Leaf Life)", side=1, line=1.5, cex=1.1)
#lines(LLs, f0LL[I], lwd=4, lty=1)
for(i in levels(spid)){
  lines(LLs[which(spid==i)], f1LL[I][which(spid==i)],  col=spid[which(spid==i)], lwd=3)
}
mtext(text = paste0("p=",round(anova(modLLvRGR, modLLvRGRnull)$`Pr(>Chisq)`[2],3)), side=3, adj = .9, line=-1, font=1, cex=cex.text )
r2tmp <- round(r.squaredGLMM(modLLvRGR)[1], 2)
mtext(text = bquote(~R[m]^2==.(r2tmp)), side=3, adj = .9, line=-2.2 ,cex=cex.text)
mtext("a)", side=3, adj=0)

### plotting LMA v RGR
f0lma <- predict(modLMAvRGR, re.form=NA)
f1lma <- fitted(modLMAvRGR)
# sort lma values so lines draw right # no NA values in LMA
I <- order(plotavs90$mlog.LMA)
LMAs <- sort(plotavs90$mlog.LMA)
spid <- plotavs90$SP.ID[I]
plot(log(RGR, base=10)~mlog.LMA, plotavs90, col=SP.ID, pch=16, cex=cex.pts, yaxt='n',ylim=ylims) #xlab="log(LMA)")#, xlim=c(2,2.6))#log(LMA)")
mtext("log(LMA)", side=1, line=1.5, cex=1.1)
#lines(c(min(LMAs), max(LMAs)), c(max(f0lma), min(f0lma)), lwd=4, lty=3)
for(i in levels(spid)){
  lines(LMAs[which(spid==i)], f1lma[I][which(spid==i)],  col=spid[which(spid==i)], lwd=3, lty=2)
}
mtext(text = paste0("p=",round(anova(modLMAvRGR, modLMAvRGRnull)$`Pr(>Chisq)`[2],3)), side=3, adj = .05, line=-1 ,cex=cex.text)
mtext("b)", side=3, adj=0)

### plotting Nmass v RGR
f0Nmass <- predict(modNmassvRGR, re.form=NA)
f1Nmass <- fitted(modNmassvRGR)
# sort lma values so lines draw right
I <- order(plotavs90$mlog.Nmass[which(!is.na(plotavs90$RGR) & !is.na(plotavs90$mlog.Nmass))])
Nmasss <- sort(plotavs90$mlog.Nmass[which(!is.na(plotavs90$RGR) & !is.na(plotavs90$mlog.Nmass))])
spid <- plotavs90$SP.ID[which(!is.na(plotavs90$RGR) & !is.na(plotavs90$mlog.Nmass))][I]
plot(log(RGR, base=10)~mlog.Nmass, plotavs90, col=SP.ID, pch=16, cex=cex.pts, xlab="", yaxt='n', xlim=c(-.22,.25),ylim=ylims)#log(Nmass)")
mtext("log(Nmass)", side=1, line=1.5, cex=1.1)
#lines(Nmasss, f0Nmass[I], lwd=4, lty=1)
# mtext("log(Relative Growth Rate)", side=2, line=2, cex=1.1, at=1)

for(i in levels(spid)){
  lines(Nmasss[which(spid==i)], f1Nmass[I][which(spid==i)],  col=spid[which(spid==i)], lwd=3)
}
mtext(text = paste0("p=",round(anova(modNmassvRGR, modNmassvRGRnull)$`Pr(>Chisq)`[2],3)), side=3, adj = .01, line=-1 ,cex=cex.text)
r2tmp <- round(r.squaredGLMM(modNmassvRGR)[1], 2)
mtext(text = bquote(~R[m]^2==.(r2tmp)), side=3, adj = .01, line=-2.2 ,cex=cex.text)
mtext("c)", side=3, adj=0)
# 

### plotting Narea v RGR
f0Narea <- predict(modNareavRGR, re.form=NA)
f1Narea <- fitted(modNareavRGR)
# sort lma values so lines draw right
I <- order(plotavs90$mlog.Narea[which(!is.na(plotavs90$RGR) & !is.na(plotavs90$mlog.Narea))])
Nareas <- sort(plotavs90$mlog.Narea[which(!is.na(plotavs90$RGR) & !is.na(plotavs90$mlog.Narea))])
spid <- plotavs90$SP.ID[which(!is.na(plotavs90$RGR) & !is.na(plotavs90$mlog.Narea))][I]
plot(log(RGR, base=10)~mlog.Narea, plotavs90, col=SP.ID, pch=16, cex=cex.pts, xlab="", yaxt="n",ylim=ylims)#log(Narea)")
#mtext("Relative Growth Rate", side=2, line=2, cex=1.1)
mtext("log(Narea)", side=1, line=1.5, cex=1.1)
#lines(Nareas, f0Narea[I], lwd=4, lty=3)
for(i in levels(spid)){
  lines(Nareas[which(spid==i)], f1Narea[I][which(spid==i)],  col=spid[which(spid==i)], lwd=3, lty=2)
}

mtext(text = paste0("p=",round(anova(modNareavRGR, modNareavRGRnull)$`Pr(>Chisq)`[2],3)), side=3, adj = .05, line=-1 ,cex=cex.text)
mtext("d)", side=3, adj=0)
#mtext("w/in spp", side=4, line=0)





###### bio ST Growth ######
quartz(width=4.3, height=5)
#pdf(file="Traits-v-Growth-90Dom_v1_RGRonly.pdf", width=4.3, height=5)
par(mar=c(2.5,2,1.2,1), oma=c(2,2,0,0),mfrow=c(2,2), cex=1,mgp=c(2.5,.5,0))

### plotting Leaf Lifespan v biomass Standardized growth
f0LL <- predict(modLLvGrowth, re.form=NA)
f1LL <- fitted(modLLvGrowth)
# sort lma values so lines draw right
I <- order(plotavs90$mlog.LL[-which(is.na(plotavs90$stGrowth) | is.na(plotavs90$mlog.LL))])
LLs <- sort(plotavs90$mlog.LL[-which(is.na(plotavs90$stGrowth) | is.na(plotavs90$mlog.LL))])
spid <- plotavs90$SP.ID[-which(is.na(plotavs90$stGrowth) | is.na(plotavs90$mlog.LL))][I]
plot(stGrowth~mlog.LL, plotavs90, col=SP.ID, pch=16, cex=.4, xlab="log(Leaf Lifespan)")
mtext("Standardized Growth", side=2, line=2.5, cex=1.1)

for(i in levels(spid)){
  lines(LLs[which(spid==i)], f1LL[I][which(spid==i)],  col=spid[which(spid==i)], lwd=2)
}
lines(LLs, f0LL[I], lwd=4, lty=3)
mtext(text = paste0("p=",round(anova(modLLvGrowth, modLLvGrowthnull)$`Pr(>Chisq)`[2],3)), side=3, adj = .8, line=-1, font=1 )
mtext("a)", side=3, adj=0)
mtext("log(Leaf Lifespan)", side=1, line=1.5, cex=1.1)

### plotting LMA v stGrowth
f0lma <- predict(modLMAvGrowth, re.form=NA)
f1lma <- fitted(modLMAvGrowth)
# sort lma values so lines draw right # no NA values in LMA
I <- order(plotavs90$mlog.LMA)
LMAs <- sort(plotavs90$mlog.LMA)
spid <- plotavs90$SP.ID[I]
plot(stGrowth~mlog.LMA, plotavs90, col=SP.ID, pch=16, cex=.4, xlab="log(LMA)")
for(i in levels(spid)){
  lines(LMAs[which(spid==i)], f1lma[I][which(spid==i)],  col=spid[which(spid==i)], lwd=2)
}
lines(c(min(LMAs), max(LMAs)), c(min(f0lma), max(f0lma)), lwd=4, lty=3)
mtext(text = paste0("p=",round(anova(modLMAvGrowth, modLMAvGrowthnull)$`Pr(>Chisq)`[2],3)), side=3, adj = .2, line=-1 )
mtext("b)", side=3, adj=0)
mtext("log(LMA)", side=1, line=1.5, cex=1.1)

### plotting Narea v stGrowth
f0Narea <- predict(modNareavGrowth, re.form=NA)
f1Narea <- fitted(modNareavGrowth)
# sort lma values so lines draw right
I <- order(plotavs90$mlog.Narea[which(!is.na(plotavs90$stGrowth) & !is.na(plotavs90$mlog.Narea))])
Nareas <- sort(plotavs90$mlog.Narea[which(!is.na(plotavs90$stGrowth) & !is.na(plotavs90$mlog.Narea))])
spid <- plotavs90$SP.ID[which(!is.na(plotavs90$stGrowth) & !is.na(plotavs90$mlog.Narea))][I]
plot(stGrowth~mlog.Narea, plotavs90, col=SP.ID, pch=16, cex=.4, xlab="log(Narea)")
for(i in levels(spid)){
  lines(Nareas[which(spid==i)], f1Narea[I][which(spid==i)],  col=spid[which(spid==i)], lwd=2)
}
lines(Nareas, f0Narea[I], lwd=4, lty=3)
mtext(text = paste0("p=",round(anova(modNareavGrowth, modNareavGrowthnull)$`Pr(>Chisq)`[2],3)), side=3, adj = .95, line=-1 )
mtext("c)", side=3, adj=0)
mtext("log(Narea)", side=1, line=1.5, cex=1.1)


### plotting Nmass v stGrowth
f0Nmass <- predict(modNmassvGrowth, re.form=NA)
f1Nmass <- fitted(modNmassvGrowth)
# sort lma values so lines draw right
I <- order(plotavs90$mlog.Nmass[which(!is.na(plotavs90$stGrowth) & !is.na(plotavs90$mlog.Nmass))])
Nmasss <- sort(plotavs90$mlog.Nmass[which(!is.na(plotavs90$stGrowth) & !is.na(plotavs90$mlog.Nmass))])
spid <- plotavs90$SP.ID[which(!is.na(plotavs90$stGrowth) & !is.na(plotavs90$mlog.Nmass))][I]
plot(stGrowth~mlog.Nmass, plotavs90, col=SP.ID, pch=16, cex=.4, xlab="log(Nmass)")
for(i in levels(spid)){
  lines(Nmasss[which(spid==i)], f1Nmass[I][which(spid==i)],  col=spid[which(spid==i)], lwd=2)
}
lines(Nmasss, f0Nmass[I], lwd=4, lty=3)
mtext(text = paste0("p=",round(anova(modNmassvGrowth, modNmassvGrowthnull)$`Pr(>Chisq)`[2],3)), side=3, adj = .2, line=-1 )
mtext("d)", side=3, adj=0)
mtext("log(Nmass)", side=1, line=1.5, cex=1.1)








############### Something Allometric #############

# strong positive relationship between RGR and Leaf Fraction
plot(log(RGR, base=10)~log(LeafFrac), biomass, col=SPP_O1_ABBREV, pch=16,cex=(ASA+200)/200)
plot(log(RGR, base=10)~log(ASA), biomass, col=SPP_O1_ABBREV, pch=16,cex=(ASA+200)/200)

plot(log(RGR, base=10)~log(I(LeafFrac*mSLA)), plotavs90, col=SP.ID, pch=16, cex=mSLA/50)
  # if you multiply LeafFrac*mSLA, things kinda linearize and seperate, but look pretty similar to just LeafFrac?
plot((RGR)~(LAI), plotavs90, col=SP.ID, pch=16, cex=mSLA/50)
  # NO relationship between RGR and LAI (or Logged values)!!!
plot(log(LeafFrac)~log(ASA), plotavs90,col=SP.ID, pch=16, cex=mSLA/50)
legend("bottomleft", legend=levels(plotavs90$SP.ID), col=mypal[1:5], pch=16)
# Strong negative relationship between log(LeafFrac) and log(ASA), UP UNTIL ABOUT ASA=90 YRS
  # then things get messy in older stands and PINCON, ABICON. (PSEMEN and PINPON are pretty dam good)
  # also, stands 43, 55 are bad again.
plot(LAI~LeafFrac, plotavs90,col=SP.ID, pch=16, cex=mSLA/50)
# not much pattern (or a hump shaped pattern) between LAI and LeafFrac

plot(LAI~ASA, plotavs90,col=SP.ID, pch=16, cex=mSLA/50)
# No relationship between LAI and ASA

plot(NPP~LeafFrac, plotavs90, col=SP.ID, pch=16, cex=(ASA+200)/200)
plot(NPP~ASA, plotavs90, col=SP.ID, pch=16, cex=(ASA+200)/200)
plot(NPP~LAI, plotavs90, col=SP.ID, pch=16, cex=(ASA+200)/200)

plot((AG_PROD_TREE_TOTAL_AS_CARBON)~LeafFrac , biomass, col=SPP_O1_ABBREV, pch=16,cex=(ASA+200)/200)
plot((AG_PROD_TREE_TOTAL_AS_CARBON)~ASA , biomass, col=SPP_O1_ABBREV, pch=16,cex=(ASA+200)/200)
plot((AG_PROD_TREE_TOTAL_AS_CARBON)~LeafAlloc , biomass, col=SPP_O1_ABBREV, pch=16,cex=(ASA+200)/200)
plot((AG_PROD_TREE_TOTAL_AS_CARBON)~LAI_O , biomass, col=SPP_O1_ABBREV, pch=16,cex=(ASA+200)/200)


## Note, the CV of SLA is .22 and LeafFrac is .96. so It looks like, within-species, LeafFrac is just WAY more variable.

plot(log(LeafFrac)~log(ASA), plotavs90, col=SP.ID, pch=16)
  #strong negative relationship between log(ASA) and log(LeafFrac), leaf frac decreases in older stands, 
  # though spp other than PSEMEN and PINPON don't really show this as much

plot(log(LeafFrac)~log(ASA), biomass, col=SPP_O1_ABBREV, pch=16)
  # This relationship breaks down as you get more old, mixed species stands

plot(log(RGR, base=10)~log(ASA), biomass, col=SPP_O1_ABBREV, pch=16)
plot(log(RGR, base=10)~log(LeafFrac), biomass, col=SPP_O1_ABBREV, pch=16)
  