
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



#_______________________________________________________________________
####################### Model Selection/Ensambling ####################
#### Modeling the trait-climate relationships 
#________________________________________________________________________

### Step 1: add in Ecoregion as a predictor
### Step 2: figure out whether I should include other species.
trait <- "log.LL"
sLL <- samplsizes <- traits.common %>% select(GE.SP, PLOT_ID, get(trait), climPC1, climPC2, soil_N, ASA, LAI_O, AG_TGROWTH,ECOREGION) %>% filter(complete.cases(.)) %>% group_by(GE.SP) %>% summarise( n= n(), nplots = length(unique(PLOT_ID)), nEcos = length(unique(ECOREGION))) %>% filter(n>15 & nplots>5) %>% arrange(n)
trait <- "log.LMA"
sLMA <- samplsizes <- traits.common %>% select(GE.SP, PLOT_ID, get(trait), climPC1, climPC2, soil_N, ASA, LAI_O, AG_TGROWTH,ECOREGION) %>% filter(complete.cases(.)) %>% group_by(GE.SP) %>% summarise( n= n(), nplots = length(unique(PLOT_ID)), nEcos = length(unique(ECOREGION))) %>% filter(n>15 & nplots>5)%>% arrange(n)
trait <- "log.Narea"
sNarea <- samplsizes <- traits.common %>% select(GE.SP, PLOT_ID, get(trait), climPC1, climPC2, soil_N, ASA, LAI_O, AG_TGROWTH,ECOREGION) %>% filter(complete.cases(.)) %>% group_by(GE.SP) %>% summarise( n= n(), nplots = length(unique(PLOT_ID)), nEcos = length(unique(ECOREGION))) %>% filter(n>15 & nplots>5)%>% arrange(n)
trait <- "log.Nmass"
sNmass <- samplsizes <- traits.common %>% select(GE.SP, PLOT_ID, get(trait), climPC1, climPC2, soil_N, ASA, LAI_O, AG_TGROWTH,ECOREGION) %>% filter(complete.cases(.)) %>% group_by(GE.SP) %>% summarise( n= n(), nplots = length(unique(PLOT_ID)), nEcos = length(unique(ECOREGION))) %>% filter(n>15 & nplots>5)%>% arrange(n)



require(MuMIn)
options(na.action = "na.fail")


trait.mods <- function(traitdata =traits.common, species, trait="log.LL", modcrit=F ){
  dataz <- data.frame(traitdata) %>% filter(SP.ID==species) %>% select(PLOT_ID, ECOREGION, get(trait), climPC1, climPC2, soil_N, ASA, LAI_O, AG_TGROWTH) %>% filter(complete.cases(.))
  print(nrow(dataz))
  dataz$log.ASA <- log(dataz$ASA)
  datazsc <- scale(dataz[,-c(1:2)]) # originally, this just scaled the predictors, but I think I also want to scale the traits.
  colnames(datazsc) <- paste0(colnames(dataz[,-c(1:2)]),'sc')
  datazall <- data.frame(dataz, datazsc)
  traitsc <- paste0(trait, "sc")
  if(length(unique(dataz$ECOREGION))==1){
    traitmod <- lmer(get(traitsc)~climPC1sc + climPC2sc + soil_Nsc + log.ASAsc + LAI_Osc + AG_TGROWTHsc + (1|PLOT_ID), datazall, REML=F)
  }
  else{
  traitmod <- lmer(get(traitsc)~climPC1sc + climPC2sc + soil_Nsc + log.ASAsc + LAI_Osc + AG_TGROWTHsc + ECOREGION + (1|PLOT_ID), datazall, REML=F)
  }
  traitdredge <- dredge(traitmod, extra=list(r.squaredGLMM))
  traitmodcalls <- dredge(traitmod, evaluate = F)
  #  traitmodavg <- model.avg(traitdredge,subset=cumsum(weight)<=.90)
  traitmodavg <- model.avg(traitdredge,subset=delta<=4)
  
  bestmodobject <- eval(traitmodcalls[[rownames(traitdredge)[1]]])
  if(modcrit==T){
    scatter.smooth(resid(bestmodobject)~fitted(bestmodobject)); abline(h=0)
    plot(datazall[,traitsc]~predict(bestmodobject, re.form=NA));abline(a=0,b=1)
    qqp(resid(bestmodobject), main="residuals")
    qqp(ranef(bestmodobject)[[1]][,1], main='Random Effects')
  }
  if(length(unique(dataz$ECOREGION))==1){
    results <- list(best = traitdredge[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc","r.squaredGLMM.R2m","r.squaredGLMM.R2c",'AICc')]
                    , avg = traitmodavg$coefficients[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc")] # note: If I try to get ecoregion effects, this can blow up. I think I need to leave ECOREGION out of this
                    , imp = as.vector(traitmodavg$importance)[match(c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc","ECOREGION"),names(traitmodavg$importance))]
                    , fulldredge = traitdredge
                    , fullmodavg = traitmodavg
                    , modnum = rownames(traitdredge)[1]
                    , bestmodcall = traitmodcalls[[rownames(traitdredge)[1]]]
                    , bestmodobject = bestmodobject
                    , n=nrow(datazall)
                    , nplots = length(ranef(bestmodobject)[[1]][,1])
                    , necos = length(unique(datazall$ECOREGION))
                    , deltaNULL = traitdredge[which(rownames(traitdredge)==1),"delta"])
  }
  else{
  results <- list(best = traitdredge[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc","ECOREGION","r.squaredGLMM.R2m","r.squaredGLMM.R2c",'AICc')]
                  , avg = traitmodavg$coefficients[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc")] # note: If I try to get ecoregion effects, this can blow up. I think I need to leave ECOREGION out of this
                  , imp = as.vector(traitmodavg$importance)[match(c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc","ECOREGION"),names(traitmodavg$importance))]
                  , fulldredge = traitdredge
                  , fullmodavg = traitmodavg
                  , modnum = rownames(traitdredge)[1]
                  , bestmodcall = traitmodcalls[[rownames(traitdredge)[1]]]
                  , bestmodobject = bestmodobject
                  , n=nrow(datazall)
                  , nplots = length(ranef(bestmodobject)[[1]][,1])
                  , necos = length(unique(datazall$ECOREGION))
                  , deltaNULL = traitdredge[which(rownames(traitdredge)==1),"delta"])
  }
  return(results)
}


####### Function for fitting models without ECOREGION ######
trait.mods.ne <- function(traitdata =traits.common, species, trait="log.LL", modcrit=F ){
  dataz <- data.frame(traitdata) %>% filter(SP.ID==species) %>% select(PLOT_ID, ECOREGION, get(trait), climPC1, climPC2, soil_N, ASA, LAI_O, AG_TGROWTH) %>% filter(complete.cases(.))
  print(nrow(dataz))
  dataz$log.ASA <- log(dataz$ASA)
  datazsc <- scale(dataz[,-c(1:2)]) # originally, this just scaled the predictors, but I think I also want to scale the traits.
  colnames(datazsc) <- paste0(colnames(dataz[,-c(1:2)]),'sc')
  datazall <- data.frame(dataz, datazsc)
  traitsc <- paste0(trait, "sc")
  traitmod <- lmer(get(traitsc)~climPC1sc + climPC2sc + soil_Nsc + log.ASAsc + LAI_Osc + AG_TGROWTHsc + (1|PLOT_ID), datazall, REML=F)
  traitdredge <- dredge(traitmod, extra=list(r.squaredGLMM))
  traitmodcalls <- dredge(traitmod, evaluate = F)
  #  traitmodavg <- model.avg(traitdredge,subset=cumsum(weight)<=.90)
  traitmodavg <- model.avg(traitdredge,subset=delta<=4)
  
  bestmodobject <- eval(traitmodcalls[[rownames(traitdredge)[1]]])
  if(modcrit==T){
    scatter.smooth(resid(bestmodobject)~fitted(bestmodobject)); abline(h=0)
    plot(datazall[,traitsc]~predict(bestmodobject, re.form=NA));abline(a=0,b=1)
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
                    , nplots = length(ranef(bestmodobject)[[1]][,1])
                    , necos = length(unique(datazall$ECOREGION))
                    , deltaNULL = traitdredge[which(rownames(traitdredge)==1),"delta"])

  return(results)
}



### run for all traits for species with 
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
PSEMENnmass <- trait.mods(traitdata = traits.common, species = "PSEMEN",trait = "log.Nmass", modcrit=T)
# n=220



PINPONll <- trait.mods(traitdata = traits.common, species = "PINPON",trait = "log.LL")
PINPONlma <- trait.mods(traitdata = traits.common, species = "PINPON",trait = "log.LMA")
PINPONnarea <- trait.mods(traitdata = traits.common, species = "PINPON",trait = "log.Narea")
# n= 97
PINPONnmass <- trait.mods(traitdata = traits.common, species = "PINPON",trait = "log.Nmass")

ABICONll <- trait.mods(traitdata = traits.common, species = "ABICON",trait = "log.LL")
ABICONlma <- trait.mods(traitdata = traits.common, species = "ABICON",trait = "log.LMA")
ABICONnarea <- trait.mods(traitdata = traits.common, species = "ABICON",trait = "log.Narea")
ABICONnmass <- trait.mods(traitdata = traits.common, species = "ABICON",trait = "log.Nmass")

TSUHETll <- trait.mods(traitdata = traits.common, species = "TSUHET",trait = "log.LL")
TSUHETlma <- trait.mods(traitdata = traits.common, species = "TSUHET",trait = "log.LMA")
TSUHETnarea <- trait.mods(traitdata = traits.common, species = "TSUHET",trait = "log.Narea") # might be one outlier that's leveraging things?
TSUHETnmass <- trait.mods(traitdata = traits.common, species = "TSUHET",trait = "log.Nmass") # might be one outlier that's leveraging things?

PINCONll <- trait.mods(traitdata = traits.common, species = "PINCON",trait = "log.LL")
PINCONlma <- trait.mods(traitdata = traits.common, species = "PINCON",trait = "log.LMA")
PINCONnarea <- trait.mods(traitdata = traits.common, species = "PINCON",trait = "log.Narea")
PINCONnmass <- trait.mods(traitdata = traits.common, species = "PINCON",trait = "log.Nmass")


PINJEFll <- trait.mods(traitdata = traits.common, species = "PINJEF",trait = "log.LL")
PINJEFlma <- trait.mods(traitdata = traits.common, species = "PINJEF",trait = "log.LMA") # two outliers
PINJEFnarea <- trait.mods(traitdata = traits.common, species = "PINJEF",trait = "log.Narea")
PINJEFnmass <- trait.mods(traitdata = traits.common, species = "PINJEF",trait = "log.Nmass")

#ARBMENll <- trait.mods(traitdata = traits.common, species = "ARBMEN",trait = "log.LL")
ARBMENlma <- trait.mods(traitdata = traits.common, species = "ARBMEN",trait = "log.LMA") 
ARBMENnarea <- trait.mods(traitdata = traits.common, species = "ARBMEN",trait = "log.Narea")
ARBMENnmass <- trait.mods(traitdata = traits.common, species = "ARBMEN",trait = "log.Nmass")


###### Apply to CWmeans  ###########
###applying to the cw traits ####
# just copied the function code to taylor to biomass

cw.trait.mods <- function(traitdata =biomass, trait="log.cw_LLp", modcrit=F ){
  dataz <- data.frame(biomass)  %>% select(PLOT_ID, ECOREGION, get(trait), climPC1, climPC2, soil_N, ASA, LAI_O, AG_PROD_TREE_TOTAL_AS_CARBON) %>% filter(complete.cases(.))
  print(nrow(dataz))
  dataz$log.ASA <- log(dataz$ASA)
  datazsc <- scale(dataz[,-c(1:2)]) # originally, this just scaled the predictors, but I think I also want to scale the traits.
  colnames(datazsc) <- paste0(colnames(dataz[,-c(1:2)]),'sc')
  colnames(datazsc)[grep("AG_", colnames(datazsc))] <- "AG_TGROWTHsc" # change this column name for compatibility
  datazall <- data.frame(dataz, datazsc)
  tn <- trait
  traitsc <- paste(tn, "sc", sep="")
  traitmod <- lm(get(traitsc)~climPC1sc + climPC2sc + soil_Nsc + log.ASAsc + LAI_Osc + AG_TGROWTHsc + ECOREGION , datazall)
  traitdredge <- dredge(traitmod, extra=list(r.squaredGLMM, "R^2"))
  traitmodcalls <- dredge(traitmod, evaluate = F)
  #  traitmodavg <- model.avg(traitdredge,subset=cumsum(weight)<=.90)
  traitmodavg <- model.avg(traitdredge,subset=delta<=4)
  
  bestmodobject <- eval(traitmodcalls[[rownames(traitdredge)[1]]])
  if(modcrit==T){
    scatter.smooth(resid(bestmodobject)~fitted(bestmodobject)); abline(h=0)
    plot(datazall[,traitsc]~predict(bestmodobject, re.form=NA));abline(a=0,b=1)
    qqp(resid(bestmodobject), main="residuals")
    # qqp(ranef(bestmodobject)[[1]][,1], main='Random Effects')
  }
  
  results <- list(best = traitdredge[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc","ECOREGION","r.squaredGLMM.R2m","R^2",'AICc')]
                  , avg = traitmodavg$coefficients[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc")] # note: If I try to get ecoregion effects, this can blow up. I think I need to leave ECOREGION out of this
                  , imp = as.vector(traitmodavg$importance)[match(c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc","ECOREGION"),names(traitmodavg$importance))]
                  , fulldredge = traitdredge
                  , fullmodavg = traitmodavg
                  , modnum = rownames(traitdredge)[1]
                  , bestmodcall = traitmodcalls[[rownames(traitdredge)[1]]]
                  , bestmodobject = bestmodobject
                  , n=nrow(datazall)
                  , nplots = nrow(datazall)
                  , necos = length(unique(datazall$ECOREGION))
                  , deltaNULL = traitdredge[which(rownames(traitdredge)==1),"delta"])
  
  return(results)
}


####### .CW Leaf Lifespan ######
trait <- "log.cw_LLp"
dataz <- data.frame(biomass)  %>% select(PLOT_ID, ECOREGION, get(trait), climPC1, climPC2, soil_N, ASA, LAI_O, AG_PROD_TREE_TOTAL_AS_CARBON) %>% filter(complete.cases(.))
print(nrow(dataz))
dataz$log.ASA <- log(dataz$ASA)
datazsc <- scale(dataz[,-c(1:2)]) # originally, this just scaled the predictors, but I think I also want to scale the traits.
colnames(datazsc) <- paste0(colnames(dataz[,-c(1:2)]),'sc')
colnames(datazsc)[grep("AG_", colnames(datazsc))] <- "AG_TGROWTHsc" # change this column name for compatibility
datazall <- data.frame(dataz, datazsc)
tn <- trait
traitsc <- paste(tn, "sc", sep="")
traitmod <- lm(get(traitsc)~climPC1sc + climPC2sc + soil_Nsc + log.ASAsc + LAI_Osc + AG_TGROWTHsc + ECOREGION , datazall)
traitdredge <- dredge(traitmod, extra=list(r.squaredGLMM, "R^2"))
traitmodcalls <- dredge(traitmod, evaluate = F)
traitmodavg <- model.avg(traitdredge,subset=delta<=4)
bestmodobject <- eval(traitmodcalls[[rownames(traitdredge)[1]]])
  scatter.smooth(resid(bestmodobject)~fitted(bestmodobject)); abline(h=0)
  plot(datazall[,traitsc]~predict(bestmodobject, re.form=NA));abline(a=0,b=1)
  qqp(resid(bestmodobject), main="residuals")
  # qqp(ranef(bestmodobject)[[1]][,1], main='Random Effects')

CWll <- list(best = traitdredge[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc","ECOREGION","r.squaredGLMM.R2m","R^2",'AICc')]
                , avg = traitmodavg$coefficients[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc")] # note: If I try to get ecoregion effects, this can blow up. I think I need to leave ECOREGION out of this
                , imp = as.vector(traitmodavg$importance)[match(c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc","ECOREGION"),names(traitmodavg$importance))]
                , fulldredge = traitdredge
                , fullmodavg = traitmodavg
                , modnum = rownames(traitdredge)[1]
                , bestmodcall = traitmodcalls[[rownames(traitdredge)[1]]]
                , bestmodobject = bestmodobject
                , n=nrow(datazall)
                , nplots = nrow(datazall)
                , necos = length(unique(datazall$ECOREGION))
                , deltaNULL = traitdredge[which(rownames(traitdredge)==1),"delta"])


####### .CW LMA ####
## NOTE: This is non-replicable, and I don't know why. You actually have to run all the code above for some strange assignment error reason
trait="log.cw_LMAp"
dataz <- data.frame(biomass)  %>% select(PLOT_ID, ECOREGION, get(trait), climPC1, climPC2, soil_N, ASA, LAI_O, AG_PROD_TREE_TOTAL_AS_CARBON) %>% filter(complete.cases(.))
print(nrow(dataz))
dataz$log.ASA <- log(dataz$ASA)
datazsc <- scale(dataz[,-c(1:2)]) # originally, this just scaled the predictors, but I think I also want to scale the traits.
colnames(datazsc) <- paste0(colnames(dataz[,-c(1:2)]),'sc')
colnames(datazsc)[grep("AG_", colnames(datazsc))] <- "AG_TGROWTHsc" # change this column name for compatibility
datazall <- data.frame(dataz, datazsc)
tn <- trait
traitsc <- paste(tn, "sc", sep="")
traitmod <- lm(get(traitsc)~climPC1sc + climPC2sc + soil_Nsc + log.ASAsc + LAI_Osc + AG_TGROWTHsc + ECOREGION , datazall)
traitdredge <- dredge(traitmod, extra=list(r.squaredGLMM, "R^2"))
traitmodcalls <- dredge(traitmod, evaluate = F)
#  traitmodavg <- model.avg(traitdredge,subset=cumsum(weight)<=.90)
traitmodavg <- model.avg(traitdredge,subset=delta<=4)

bestmodobject <- eval(traitmodcalls[[rownames(traitdredge)[1]]])

scatter.smooth(resid(bestmodobject)~fitted(bestmodobject)); abline(h=0)
plot(datazall[,traitsc]~predict(bestmodobject, re.form=NA));abline(a=0,b=1)
qqp(resid(bestmodobject), main="residuals")
# qqp(ranef(bestmodobject)[[1]][,1], main='Random Effects')


CWlma <- list(best = traitdredge[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc","ECOREGION","r.squaredGLMM.R2m","R^2",'AICc')]
              , avg = traitmodavg$coefficients[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc")] # note: If I try to get ecoregion effects, this can blow up. I think I need to leave ECOREGION out of this
              , imp = as.vector(traitmodavg$importance)[match(c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc","ECOREGION"),names(traitmodavg$importance))]
              , fulldredge = traitdredge
              , fullmodavg = traitmodavg
              , modnum = rownames(traitdredge)[1]
              , bestmodcall = traitmodcalls[[rownames(traitdredge)[1]]]
              , bestmodobject = bestmodobject
              , n=nrow(datazall)
              , nplots = nrow(datazall)
              , necos = length(unique(datazall$ECOREGION))
              , deltaNULL = traitdredge[which(rownames(traitdredge)==1),"delta"])


#### Nmass ####
trait="log.cw_Nmassp"
dataz <- data.frame(biomass)  %>% select(PLOT_ID, ECOREGION, get(trait), climPC1, climPC2, soil_N, ASA, LAI_O, AG_PROD_TREE_TOTAL_AS_CARBON) %>% filter(complete.cases(.))
print(nrow(dataz))
dataz$log.ASA <- log(dataz$ASA)
datazsc <- scale(dataz[,-c(1:2)]) # originally, this just scaled the predictors, but I think I also want to scale the traits.
colnames(datazsc) <- paste0(colnames(dataz[,-c(1:2)]),'sc')
colnames(datazsc)[grep("AG_", colnames(datazsc))] <- "AG_TGROWTHsc" # change this column name for compatibility
datazall <- data.frame(dataz, datazsc)
tn <- trait
traitsc <- paste(tn, "sc", sep="")
traitmod <- lm(get(traitsc)~climPC1sc + climPC2sc + soil_Nsc + log.ASAsc + LAI_Osc + AG_TGROWTHsc + ECOREGION , datazall)
traitdredge <- dredge(traitmod, extra=list(r.squaredGLMM, "R^2"))
traitmodcalls <- dredge(traitmod, evaluate = F)
#  traitmodavg <- model.avg(traitdredge,subset=cumsum(weight)<=.90)
traitmodavg <- model.avg(traitdredge,subset=delta<=4)

bestmodobject <- eval(traitmodcalls[[rownames(traitdredge)[1]]])

scatter.smooth(resid(bestmodobject)~fitted(bestmodobject)); abline(h=0)
plot(datazall[,traitsc]~predict(bestmodobject, re.form=NA));abline(a=0,b=1)
qqp(resid(bestmodobject), main="residuals")
# qqp(ranef(bestmodobject)[[1]][,1], main='Random Effects')


CWnmass <- list(best = traitdredge[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc","ECOREGION","r.squaredGLMM.R2m","R^2",'AICc')]
                , avg = traitmodavg$coefficients[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc")] # note: If I try to get ecoregion effects, this can blow up. I think I need to leave ECOREGION out of this
                , imp = as.vector(traitmodavg$importance)[match(c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc","ECOREGION"),names(traitmodavg$importance))]
                , fulldredge = traitdredge
                , fullmodavg = traitmodavg
                , modnum = rownames(traitdredge)[1]
                , bestmodcall = traitmodcalls[[rownames(traitdredge)[1]]]
                , bestmodobject = bestmodobject
                , n=nrow(datazall)
                , nplots = nrow(datazall)
                , necos = length(unique(datazall$ECOREGION))
                , deltaNULL = traitdredge[which(rownames(traitdredge)==1),"delta"])


###### Narea ###########
trait <- "log.cw_Nareap"
dataz <- data.frame(biomass)  %>% select(PLOT_ID, ECOREGION, get(trait), climPC1, climPC2, soil_N, ASA, LAI_O, AG_PROD_TREE_TOTAL_AS_CARBON) %>% filter(complete.cases(.))
print(nrow(dataz))
dataz$log.ASA <- log(dataz$ASA)
datazsc <- scale(dataz[,-c(1:2)]) # originally, this just scaled the predictors, but I think I also want to scale the traits.
colnames(datazsc) <- paste0(colnames(dataz[,-c(1:2)]),'sc')
colnames(datazsc)[grep("AG_", colnames(datazsc))] <- "AG_TGROWTHsc" # change this column name for compatibility
datazall <- data.frame(dataz, datazsc)
tn <- trait
traitsc <- paste(tn, "sc", sep="")
traitmod <- lm(get(traitsc)~climPC1sc + climPC2sc + soil_Nsc + log.ASAsc + LAI_Osc + AG_TGROWTHsc + ECOREGION , datazall)
traitdredge <- dredge(traitmod, extra=list(r.squaredGLMM, "R^2"))
traitmodcalls <- dredge(traitmod, evaluate = F)
#  traitmodavg <- model.avg(traitdredge,subset=cumsum(weight)<=.90)
traitmodavg <- model.avg(traitdredge,subset=delta<=4)

bestmodobject <- eval(traitmodcalls[[rownames(traitdredge)[1]]])

scatter.smooth(resid(bestmodobject)~fitted(bestmodobject)); abline(h=0)
plot(datazall[,traitsc]~predict(bestmodobject, re.form=NA));abline(a=0,b=1)
qqp(resid(bestmodobject), main="residuals")
# qqp(ranef(bestmodobject)[[1]][,1], main='Random Effects')


CWnarea <- list(best = traitdredge[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc","ECOREGION","r.squaredGLMM.R2m","R^2",'AICc')]
                , avg = traitmodavg$coefficients[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc")] # note: If I try to get ecoregion effects, this can blow up. I think I need to leave ECOREGION out of this
                , imp = as.vector(traitmodavg$importance)[match(c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc","ECOREGION"),names(traitmodavg$importance))]
                , fulldredge = traitdredge
                , fullmodavg = traitmodavg
                , modnum = rownames(traitdredge)[1]
                , bestmodcall = traitmodcalls[[rownames(traitdredge)[1]]]
                , bestmodobject = bestmodobject
                , n=nrow(datazall)
                , nplots = nrow(datazall)
                , necos = length(unique(datazall$ECOREGION))
                , deltaNULL = traitdredge[which(rownames(traitdredge)==1),"delta"])





colnames(CWll$best)[10] <- "r.squaredGLMM.R2c"
colnames(CWlma$best)[10] <- "r.squaredGLMM.R2c"
colnames(CWnmass$best)[10] <- "r.squaredGLMM.R2c"
colnames(CWnarea$best)[10] <- "r.squaredGLMM.R2c"



####### Combining all environmental trait models ######
llbestmods <- rbind(PSEMENll$best, PINPONll$best,PINCONll$best,PINJEFll$best, ABICONll$best, TSUHETll$best, CWll$best)
llbestmods$n <- rbind(PSEMENll$n, PINPONll$n,PINCONll$n,PINJEFll$n, ABICONll$n, TSUHETll$n, CWll$n)
llbestmods$nplots <- rbind(PSEMENll$nplots, PINPONll$nplots,PINCONll$nplots,PINJEFll$nplots, ABICONll$nplots, TSUHETll$nplots, CWll$nplots)
llbestmods$necos <- rbind(PSEMENll$necos, PINPONll$necos,PINCONll$necos,PINJEFll$necos, ABICONll$necos, TSUHETll$necos, CWll$necos)
llbestmods$deltas <-rbind(PSEMENll$deltaNULL, PINPONll$deltaNULL,PINCONll$deltaNULL,PINJEFll$deltaNULL, ABICONll$deltaNULL, TSUHETll$deltaNULL,CWll$deltaNULL)
llbestmods$SP.ID <-  c("PSEMEN","PINPON","PINCON","PINJEF","ABICON","TSUHET","CWmean")

lmabestmods <- rbind(PSEMENlma$best, PINPONlma$best,PINCONlma$best,PINJEFlma$best, ABICONlma$best, TSUHETlma$best, CWlma$best)
lmabestmods$n <- rbind(PSEMENlma$n, PINPONlma$n,PINCONlma$n,PINJEFlma$n, ABICONlma$n, TSUHETlma$n, CWlma$n)
lmabestmods$nplots <- rbind(PSEMENlma$nplots, PINPONlma$nplots,PINCONlma$nplots,PINJEFlma$nplots, ABICONlma$nplots, TSUHETlma$nplots,CWlma$n)
lmabestmods$necos <- rbind(PSEMENlma$necos, PINPONlma$necos,PINCONlma$necos,PINJEFlma$necos, ABICONlma$necos, TSUHETlma$necos,CWlma$necos)
lmabestmods$deltas <-rbind(PSEMENlma$deltaNULL, PINPONlma$deltaNULL,PINCONlma$deltaNULL,PINJEFlma$deltaNULL, ABICONlma$deltaNULL, TSUHETlma$deltaNULL,CWlma$deltaNULL)
lmabestmods$SP.ID <-  c("PSEMEN","PINPON","PINCON","PINJEF","ABICON","TSUHET", "CWmean")


nareabestmods <- rbind(PSEMENnarea$best, PINPONnarea$best,PINCONnarea$best,PINJEFnarea$best, ABICONnarea$best, TSUHETnarea$best, CWnarea$best)
nareabestmods$n <- rbind(PSEMENnarea$n, PINPONnarea$n,PINCONnarea$n,PINJEFnarea$n, ABICONnarea$n, TSUHETnarea$n, CWnarea$n)
nareabestmods$nplots <- rbind(PSEMENnarea$nplots, PINPONnarea$nplots,PINCONnarea$nplots,PINJEFnarea$nplots, ABICONnarea$nplots, TSUHETnarea$nplots, CWnarea$nplots)
nareabestmods$necos <- rbind(PSEMENnarea$necos, PINPONnarea$necos,PINCONnarea$necos,PINJEFnarea$necos, ABICONnarea$necos, TSUHETnarea$necos, CWnarea$necos)
nareabestmods$deltas <-rbind(PSEMENnarea$deltaNULL, PINPONnarea$deltaNULL,PINCONnarea$deltaNULL,PINJEFnarea$deltaNULL, ABICONnarea$deltaNULL, TSUHETnarea$deltaNULL, CWnarea$deltaNULL)
nareabestmods$SP.ID <-  c("PSEMEN","PINPON","PINCON","PINJEF","ABICON","TSUHET", "CWmean")
nareabestmods$Species <-  c("Pseudotsuga menziesii","Pinus ponderosa","Pinus contorta","Pinus jeffreyii","Abies concolor","Tsuga heterophylla", "Community-Weighted Mean")


nmassbestmods <- rbind(PSEMENnmass$best, PINPONnmass$best,PINCONnmass$best,PINJEFnmass$best, ABICONnmass$best, TSUHETnmass$best, CWnmass$best)
nmassbestmods$n <- rbind(PSEMENnmass$n, PINPONnmass$n,PINCONnmass$n,PINJEFnmass$n, ABICONnmass$n, TSUHETnmass$n, CWnmass$n)
nmassbestmods$nplots <- rbind(PSEMENnmass$nplots, PINPONnmass$nplots,PINCONnmass$nplots,PINJEFnmass$nplots, ABICONnmass$nplots, TSUHETnmass$nplots, CWnmass$nplots)
nmassbestmods$necos <- rbind(PSEMENnmass$necos, PINPONnmass$necos,PINCONnmass$necos,PINJEFnmass$necos, ABICONnmass$necos, TSUHETnmass$necos, CWnmass$necos)
nmassbestmods$deltas <-rbind(PSEMENnmass$deltaNULL, PINPONnmass$deltaNULL,PINCONnmass$deltaNULL,PINJEFnmass$deltaNULL, ABICONnmass$deltaNULL, TSUHETnmass$deltaNULL, CWnmass$deltaNULL)
nmassbestmods$SP.ID <-  c("PSEMEN","PINPON","PINCON","PINJEF","ABICON","TSUHET", "CWmean")
nmassbestmods$Species <-  c("Pseudotsuga menziesii","Pinus ponderosa","Pinus contorta","Pinus jeffreyii","Abies concolor","Tsuga heterophylla", "Community-Weighted Mean")




#___________________________________________________________________________________
############ Fitting all models without ecoregions #############
#___________________________________________________________________________________


PSEMENll.ne <- trait.mods.ne(traitdata = traits.common, species = "PSEMEN",trait = "log.LL", modcrit=T)
PSEMENlma.ne <- trait.mods.ne(traitdata = traits.common, species = "PSEMEN",trait = "log.LMA", modcrit=T)
PSEMENnarea.ne <- trait.mods.ne(traitdata = traits.common, species = "PSEMEN",trait = "log.Narea", modcrit=T)
PSEMENnmass.ne <- trait.mods.ne(traitdata = traits.common, species = "PSEMEN",trait = "log.Nmass", modcrit=T)

PINPONll.ne <- trait.mods.ne(traitdata = traits.common, species = "PINPON",trait = "log.LL")
PINPONlma.ne <- trait.mods.ne(traitdata = traits.common, species = "PINPON",trait = "log.LMA")
PINPONnarea.ne <- trait.mods.ne(traitdata = traits.common, species = "PINPON",trait = "log.Narea")
PINPONnmass.ne <- trait.mods.ne(traitdata = traits.common, species = "PINPON",trait = "log.Nmass")

ABICONll.ne <- trait.mods.ne(traitdata = traits.common, species = "ABICON",trait = "log.LL")
ABICONlma.ne <- trait.mods.ne(traitdata = traits.common, species = "ABICON",trait = "log.LMA")
ABICONnarea.ne <- trait.mods.ne(traitdata = traits.common, species = "ABICON",trait = "log.Narea")
ABICONnmass.ne <- trait.mods.ne(traitdata = traits.common, species = "ABICON",trait = "log.Nmass")

TSUHETll.ne <- trait.mods.ne(traitdata = traits.common, species = "TSUHET",trait = "log.LL")
TSUHETlma.ne <- trait.mods.ne(traitdata = traits.common, species = "TSUHET",trait = "log.LMA")
TSUHETnarea.ne <- trait.mods.ne(traitdata = traits.common, species = "TSUHET",trait = "log.Narea") # might be one outlier that's leveraging things?
TSUHETnmass.ne <- trait.mods.ne(traitdata = traits.common, species = "TSUHET",trait = "log.Nmass") # might be one outlier that's leveraging things?

PINCONll.ne <- trait.mods.ne(traitdata = traits.common, species = "PINCON",trait = "log.LL")
PINCONlma.ne <- trait.mods.ne(traitdata = traits.common, species = "PINCON",trait = "log.LMA")
PINCONnarea.ne <- trait.mods.ne(traitdata = traits.common, species = "PINCON",trait = "log.Narea")
PINCONnmass.ne <- trait.mods.ne(traitdata = traits.common, species = "PINCON",trait = "log.Nmass")

PINJEFll.ne <- trait.mods.ne(traitdata = traits.common, species = "PINJEF",trait = "log.LL")
PINJEFlma.ne <- trait.mods.ne(traitdata = traits.common, species = "PINJEF",trait = "log.LMA") # two outliers
PINJEFnarea.ne <- trait.mods.ne(traitdata = traits.common, species = "PINJEF",trait = "log.Narea")
PINJEFnmass.ne <- trait.mods.ne(traitdata = traits.common, species = "PINJEF",trait = "log.Nmass")

# #ARBMENll.ne <- trait.mods.ne(traitdata = traits.common, species = "ARBMEN",trait = "log.LL")
# ARBMENlma.ne <- trait.mods.ne(traitdata = traits.common, species = "ARBMEN",trait = "log.LMA") 
# ARBMENnarea.ne <- trait.mods.ne(traitdata = traits.common, species = "ARBMEN",trait = "log.Narea")
# ARBMENnmass.ne <- trait.mods.ne(traitdata = traits.common, species = "ARBMEN",trait = "log.Nmass")



####### Combining all environmental trait models without ECOREGION######
llbestmods.ne <- rbind(PSEMENll.ne$best, PINPONll.ne$best,PINCONll.ne$best,PINJEFll.ne$best, ABICONll.ne$best, TSUHETll.ne$best)
llbestmods.ne$n <- rbind(PSEMENll.ne$n, PINPONll.ne$n,PINCONll.ne$n,PINJEFll.ne$n, ABICONll.ne$n, TSUHETll.ne$n)
llbestmods.ne$nplots <- rbind(PSEMENll.ne$nplots, PINPONll.ne$nplots,PINCONll.ne$nplots,PINJEFll.ne$nplots, ABICONll.ne$nplots, TSUHETll.ne$nplots)
llbestmods.ne$necos <- rbind(PSEMENll.ne$necos, PINPONll.ne$necos,PINCONll.ne$necos,PINJEFll.ne$necos, ABICONll.ne$necos, TSUHETll.ne$necos)
llbestmods.ne$deltas <-rbind(PSEMENll.ne$deltaNULL, PINPONll.ne$deltaNULL,PINCONll.ne$deltaNULL,PINJEFll.ne$deltaNULL, ABICONll.ne$deltaNULL, TSUHETll.ne$deltaNULL)
llbestmods.ne$SP.ID <-  c("PSEMEN","PINPON","PINCON","PINJEF","ABICON","TSUHET")

lmabestmods.ne <- rbind(PSEMENlma.ne$best, PINPONlma.ne$best,PINCONlma.ne$best,PINJEFlma.ne$best, ABICONlma.ne$best, TSUHETlma.ne$best)
lmabestmods.ne$n <- rbind(PSEMENlma.ne$n, PINPONlma.ne$n,PINCONlma.ne$n,PINJEFlma.ne$n, ABICONlma.ne$n, TSUHETlma.ne$n)
lmabestmods.ne$nplots <- rbind(PSEMENlma.ne$nplots, PINPONlma.ne$nplots,PINCONlma.ne$nplots,PINJEFlma.ne$nplots, ABICONlma.ne$nplots, TSUHETlma.ne$nplots)
lmabestmods.ne$necos <- rbind(PSEMENlma.ne$necos, PINPONlma.ne$necos,PINCONlma.ne$necos,PINJEFlma.ne$necos, ABICONlma.ne$necos, TSUHETlma.ne$necos)
lmabestmods.ne$deltas <-rbind(PSEMENlma.ne$deltaNULL, PINPONlma.ne$deltaNULL,PINCONlma.ne$deltaNULL,PINJEFlma.ne$deltaNULL, ABICONlma.ne$deltaNULL, TSUHETlma.ne$deltaNULL)
lmabestmods.ne$SP.ID <-  c("PSEMEN","PINPON","PINCON","PINJEF","ABICON","TSUHET")


nareabestmods.ne <- rbind(PSEMENnarea.ne$best, PINPONnarea.ne$best,PINCONnarea.ne$best,PINJEFnarea.ne$best, ABICONnarea.ne$best, TSUHETnarea.ne$best)
nareabestmods.ne$n <- rbind(PSEMENnarea.ne$n, PINPONnarea.ne$n,PINCONnarea.ne$n,PINJEFnarea.ne$n, ABICONnarea.ne$n, TSUHETnarea.ne$n)
nareabestmods.ne$nplots <- rbind(PSEMENnarea.ne$nplots, PINPONnarea.ne$nplots,PINCONnarea.ne$nplots,PINJEFnarea.ne$nplots, ABICONnarea.ne$nplots, TSUHETnarea.ne$nplots)
nareabestmods.ne$necos <- rbind(PSEMENnarea.ne$necos, PINPONnarea.ne$necos,PINCONnarea.ne$necos,PINJEFnarea.ne$necos, ABICONnarea.ne$necos, TSUHETnarea.ne$necos)
nareabestmods.ne$deltas <-rbind(PSEMENnarea.ne$deltaNULL, PINPONnarea.ne$deltaNULL,PINCONnarea.ne$deltaNULL,PINJEFnarea.ne$deltaNULL, ABICONnarea.ne$deltaNULL, TSUHETnarea.ne$deltaNULL)
nareabestmods.ne$SP.ID <-  c("PSEMEN","PINPON","PINCON","PINJEF","ABICON","TSUHET")
nareabestmods.ne$Species <-  c("Pseudotsuga menziesii","Pinus ponderosa","Pinus contorta","Pinus jeffreyii","Abies concolor","Tsuga heterophylla")


nmassbestmods.ne <- rbind(PSEMENnmass.ne$best, PINPONnmass.ne$best,PINCONnmass.ne$best,PINJEFnmass.ne$best, ABICONnmass.ne$best, TSUHETnmass.ne$best)
nmassbestmods.ne$n <- rbind(PSEMENnmass.ne$n, PINPONnmass.ne$n,PINCONnmass.ne$n,PINJEFnmass.ne$n, ABICONnmass.ne$n, TSUHETnmass.ne$n)
nmassbestmods.ne$nplots <- rbind(PSEMENnmass.ne$nplots, PINPONnmass.ne$nplots,PINCONnmass.ne$nplots,PINJEFnmass.ne$nplots, ABICONnmass.ne$nplots, TSUHETnmass.ne$nplots)
nmassbestmods.ne$necos <- rbind(PSEMENnmass.ne$necos, PINPONnmass.ne$necos,PINCONnmass.ne$necos,PINJEFnmass.ne$necos, ABICONnmass.ne$necos, TSUHETnmass.ne$necos)
nmassbestmods.ne$deltas <-rbind(PSEMENnmass.ne$deltaNULL, PINPONnmass.ne$deltaNULL,PINCONnmass.ne$deltaNULL,PINJEFnmass.ne$deltaNULL, ABICONnmass.ne$deltaNULL, TSUHETnmass.ne$deltaNULL)
nmassbestmods.ne$SP.ID <-  c("PSEMEN","PINPON","PINCON","PINJEF","ABICON","TSUHET")
nmassbestmods.ne$Species <-  c("Pseudotsuga menziesii","Pinus ponderosa","Pinus contorta","Pinus jeffreyii","Abies concolor","Tsuga heterophylla")







# write.csv(llbestmods, "Trait_models/LeafLife_bestmodels_031617.csv")
# write.csv(lmabestmods, "Trait_models/LMA_bestmodels_031617.csv")
# write.csv(nareabestmods, "Trait_models/Narea_bestmodels_031617.csv")
# write.csv(nmassbestmods, "Trait_models/Nmass_bestmodels_041517.csv")
# with ECOREGION in the models
# write.csv(llbestmods, "Trait_models/LeafLife_bestmodels_wECOREG_061617.csv")
# write.csv(lmabestmods, "Trait_models/LMA_bestmodels_wECOREG_061617.csv")
# write.csv(nareabestmods, "Trait_models/Narea_bestmodels_wECOREG_061617.csv")
# write.csv(nmassbestmods, "Trait_models/Nmass_bestmodels_wECOREG_061617.csv")

#### Rerun for quality control with and without Ecoregions #### -> also switched my date convention to to go year month day
# write.csv(llbestmods, "Trait_models/LeafLife_bestmodels_20170619.csv")
# write.csv(lmabestmods, "Trait_models/LMA_bestmodels_20170619.csv")
# write.csv(nareabestmods, "Trait_models/Narea_bestmodels_20170619.csv")
# write.csv(nmassbestmods, "Trait_models/Nmass_bestmodels_20170619.csv")
write.csv(llbestmods.ne, "Trait_models/LeafLife_bestmodels_noECOREG_20170619.csv")
write.csv(lmabestmods.ne, "Trait_models/LMA_bestmodels_noECOREG_20170619.csv")
write.csv(nareabestmods.ne, "Trait_models/Narea_bestmodels_noECOREG_20170619.csv")
write.csv(nmassbestmods.ne, "Trait_models/Nmass_bestmodels_noECOREG_20170619.csv")
# _v2 has CWmeans added on
write.csv(llbestmods, "Trait_models/LeafLife_bestmodels_20170619_v2.csv")
write.csv(lmabestmods, "Trait_models/LMA_bestmodels_20170619_v2.csv")
write.csv(nareabestmods, "Trait_models/Narea_bestmodels_20170619_v2.csv")
write.csv(nmassbestmods, "Trait_models/Nmass_bestmodels_20170619_v2.csv")



#### Variable Importances ######
llimps <- data.frame(rbind(PSEMENll$imp, PINPONll$imp,PINCONll$imp,PINJEFll$imp, ABICONll$imp, TSUHETll$imp, CWll$imp)[,-1])
colnames(llimps) <- c("climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc","ECOREGION")
llimps$SP.ID <- c("PSEMEN","PINPON","PINCON","PINJEF","ABICON","TSUHET", "CWmean")
llimps$trait <- rep("LeafLife", times=nrow(llimps))

lmaimps <- data.frame(rbind(PSEMENlma$imp, PINPONlma$imp,PINCONlma$imp,PINJEFlma$imp, ABICONlma$imp, TSUHETlma$imp, CWlma$imp)[,-1])
colnames(lmaimps) <- c("climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc","ECOREGION")
lmaimps$SP.ID <- c("PSEMEN","PINPON","PINCON","PINJEF","ABICON","TSUHET", "CWmean")
lmaimps$trait <- rep("LMA", times=nrow(lmaimps))

nareaimps <- data.frame(rbind(PSEMENnarea$imp, PINPONnarea$imp,PINCONnarea$imp,PINJEFnarea$imp, ABICONnarea$imp, TSUHETnarea$imp, CWnarea$imp)[,-1])
colnames(nareaimps) <- c("climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc","ECOREGION")
nareaimps$SP.ID <- c("PSEMEN","PINPON","PINCON","PINJEF","ABICON","TSUHET", "CWmean")
nareaimps$trait <- rep("Narea", times=nrow(nareaimps))


nmassimps <- data.frame(rbind(PSEMENnmass$imp, PINPONnmass$imp,PINCONnmass$imp,PINJEFnmass$imp, ABICONnmass$imp, TSUHETnmass$imp, CWnmass$imp)[,-1])
colnames(nmassimps) <- c("climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc","ECOREGION")
nmassimps$SP.ID <- c("PSEMEN","PINPON","PINCON","PINJEF","ABICON","TSUHET", "CWmean")
nmassimps$trait <- rep("Nmass", times=nrow(nmassimps))


importances <- rbind(llimps, lmaimps, nareaimps, nmassimps)
colnames(importances) <- c("climPC1", "climPC2","soil_N","Stand_Age","LAI","Growth", "ECOREGION", "SP.ID", "trait")
importances$ECOREGION[which(is.na(importances$ECOREGION))] <- 0 # change NAs to 0s (i.e. when no model contained ECOREGION, it should have an importance of 0)
# write.csv(importances, "Trait_models/VariableImportances_wide_031617aiccuttoff.csv")
# write.csv(importances, "Trait_models/VariableImportances_wide_041517aiccuttoff.csv")
# write.csv(importances, "Trait_models/VariableImportances_wide_wECOREGION_061617aiccuttoff.csv")
# write.csv(importances, "Trait_models/VariableImportances_wide_20170619aiccuttoff.csv")
write.csv(importances, "Trait_models/VariableImportances_wide_20170619_v2aiccuttoff.csv") # same models but with CWmeans added



##### Model Averaged Effects Sizes ######
llavgmods <- data.frame(rbind(PSEMENll$avg, PINPONll$avg,PINCONll$avg,PINJEFll$avg, ABICONll$avg, TSUHETll$avg, CWll$avg))
llavgmods$SP.ID <- c("PSEMEN","PINPON","PINCON","PINJEF","ABICON","TSUHET","CWmean")
llavgmods$trait <- rep("LeafLife", times=nrow(llavgmods))

lmaavgmods <- data.frame(rbind(PSEMENlma$avg, PINPONlma$avg,PINCONlma$avg,PINJEFlma$avg, ABICONlma$avg, TSUHETlma$avg, CWlma$avg))
lmaavgmods$SP.ID <- c("PSEMEN","PINPON","PINCON","PINJEF","ABICON","TSUHET", "CWmean")
lmaavgmods$trait <- rep("LMA", times=nrow(lmaavgmods))

nareaavgmods <- data.frame(rbind(PSEMENnarea$avg, PINPONnarea$avg,PINCONnarea$avg,PINJEFnarea$avg, ABICONnarea$avg, TSUHETnarea$avg, CWnarea$avg))
nareaavgmods$SP.ID <- c("PSEMEN","PINPON","PINCON","PINJEF","ABICON","TSUHET","CWmean")
nareaavgmods$trait <- rep("Narea", times=nrow(nareaavgmods))

nmassavgmods <- data.frame(rbind(PSEMENnmass$avg, PINPONnmass$avg,PINCONnmass$avg,PINJEFnmass$avg, ABICONnmass$avg, TSUHETnmass$avg, CWnmass$avg))
nmassavgmods$SP.ID <- c("PSEMEN","PINPON","PINCON","PINJEF","ABICON","TSUHET","CWmean")
nmassavgmods$trait <- rep("Nmass", times=nrow(nmassavgmods))

avgmods <- rbind(llavgmods, lmaavgmods, nareaavgmods, nmassavgmods)[,-1]
colnames(avgmods) <- c("climPC1", "climPC2","soil_N","Stand_Age","LAI","Growth", "SP.ID", "trait")
#write.csv(avgmods, "Trait_models/Model_Averages_wide_031617aiccuttoff.csv")
#write.csv(avgmods, "Trait_models/Model_Averages_wide_041517aiccuttoff.csv")
#write.csv(avgmods, "Trait_models/Model_Averages_wide_061617aiccutoff.csv")
#write.csv(avgmods, "Trait_models/Model_Averages_wide_20170619aiccutoff.csv")
write.csv(avgmods, "Trait_models/Model_Averages_wide_20170619_v2aiccutoff.csv")


#_______________________________________________________________________________________________
########### Plotting Model Ensamble Output ########################
### Currently no code for replotting the importances w/out ECOREGION (ie .ne results)
#_______________________________________________________________________________________________

importances <- read.csv("Trait_models/VariableImportances_wide_20170619_v2aiccuttoff.csv", row.names = 1)
#importances.we <- read.csv("Trait_models/VariableImportances_wide_wECOREGION_061617aiccuttoff.csv", row.names = 1)
#importances.we$ECOREGION[which(is.na(importances.we$ECOREGION))] <- 0
impslong <- melt(importances, id.vars = c("SP.ID","trait"))
#impslong.we <- melt(importances.we, id.vars= c("SP.ID", "trait"))

#avgmods <- read.csv("Trait_models/Model_Averages_wide_041517aiccuttoff.csv", row.names = 1)
avgmods <- read.csv("Trait_models/Model_Averages_wide_20170619_v2aiccutoff.csv", row.names = 1)
avglong <- melt(avgmods, id.vars = c("SP.ID","trait"))

avglong$xvals <- as.numeric(avglong$variable)
impslong$xvals <- as.numeric(impslong$variable)



# old code for adding ECOREGION in to boxplots. Not needed with dotplots
#levels(avglong$variable) <- list(climPC1="climPC1",climPC2= "climPC2", soil_N="soil_N",Stand_Age="Stand_Age",LAI="LAI",Growth="Growth", Ecoregion="ECOREGION")
#avglong$variable <- as.character(avglong$variable)
#tmp <- data.frame(SP.ID="PSEMEN", trait=c("LeafLife","LMA","Narea","Nmass"), vairable = "Ecoregion", value=NA)
#avglong[145,] <- data.frame(SP.ID=NA, trait=NA, vairable = "Ecoregion", value=NA)
#avglong$variable <- as.factor(avglong$variable)



# ##### OLD: Boxplots #########
# p1 <- ggplot(impslong, aes(x=variable, y=value, col=variable)) + geom_boxplot() + 
#   geom_abline(slope = 0, intercept = 0) + facet_grid(trait~.) +
#   theme(axis.text.x=element_text(angle=90), legend.position="none") + ylab("Variable Importnce")
# 
# (p2 <- ggplot(avglong, aes(x=variable, y=value, col=variable)) + geom_boxplot() + 
#     geom_abline(slope = 0, intercept = 0) + facet_grid(trait~.) +
#     theme(axis.text.x=element_text(angle=90), legend.position="None") + ylab("Model-averaged Standardized Effect Sizes"))
# 
# 
# quartz(width=6, height=5)
# multiplot(p2,p1, cols = 2)
# 
# #__________________________________________________




########## NEW: Dotplots with CW traits added on top #########



quartz(width=4.33, height=6) # full page = 6.81
par(mfcol=c(4,2), mar=c(0,3.5,0,0), oma=c(5,0,2,1), mgp=c(2.5,.75,0), cex=.9)
panlab.cex <- .9
panlab.ln <- -1.2
CWbg <- "white"

### Effect Sizes #####
boxplot(value~xvals, avglong[which(avglong$trait=="LeafLife" & avglong$SP.ID!="CWmean"),], at=c(1,2,3,4,5,6)
        , ylim=c(-.7,1), xaxt="n", xlim=c(0.5,7.5)
        ,outpch=1, border="grey", col="grey", outcol="grey"
        , staplewex=0, notch=F, medcol="black" #, medlwd=0
        , whisklty=1, whisklwd=2, boxwex=.7, whiskcol="grey" #, boxwex=1, medlwd=3
        ,las=2)
#plot(value~xvals, avglong[which(avglong$trait=="LeafLife" & avglong$SP.ID!="CWmean"),], pch=16, xaxt="n", ylab="", ylim=c(-0.7,1), xlab="", xlim=c(0.6,6.4))
abline(h=0)
points(value~xvals, avglong[which(avglong$trait=="LeafLife" & avglong$SP.ID=="CWmean"),], pch=24, bg=CWbg)
mtext(text = "B)",side = 3,line = panlab.ln,adj=0.02,cex=panlab.cex )
mtext(text = "log(LL)",side = 3,line = panlab.ln,adj=0.9,cex=panlab.cex)
text(x=7, y=0, pos=3, labels = "NA",cex=.7)

boxplot(value~xvals, avglong[which(avglong$trait=="LMA" & avglong$SP.ID!="CWmean"),], at=c(1,2,3,4,5,6)
        , ylim=c(-.7,1), xaxt="n", xlim=c(0.5,7.5)
        ,outpch=1, border="grey", col="grey", outcol="grey"
        , staplewex=0, notch=F, medcol="black" #, medlwd=0
        , whisklty=1, whisklwd=2, boxwex=.7, whiskcol="grey" #, boxwex=1, medlwd=3
        ,las=2)
#plot(value~xvals, avglong[which(avglong$trait=="LMA" & avglong$SP.ID!="CWmean"),], pch=16, xaxt="n", ylab="", ylim=c(-0.7,1), xlab="", xlim=c(0.6,6.4))
abline(h=0)
points(value~xvals, avglong[which(avglong$trait=="LMA" & avglong$SP.ID=="CWmean"),], pch=24, bg=CWbg)
mtext(text = "D)",side = 3,line = panlab.ln,adj=0.02,cex=panlab.cex )
mtext(text = "log(LMA)",side = 3,line = panlab.ln,adj=0.9,cex=panlab.cex)
text(x=7, y=0, pos=3, labels = "NA",cex=.7)

boxplot(value~xvals, avglong[which(avglong$trait=="Nmass" & avglong$SP.ID!="CWmean"),], at=c(1,2,3,4,5,6)
        , ylim=c(-.7,1), xaxt="n", xlim=c(0.5,7.5)
        ,outpch=1, border="grey", col="grey", outcol="grey"
        , staplewex=0, notch=F, medcol="black" #, medlwd=0
        , whisklty=1, whisklwd=2, boxwex=.7, whiskcol="grey" #, boxwex=1, medlwd=3
        ,las=2)
# plot(value~xvals, avglong[which(avglong$trait=="Nmass" & avglong$SP.ID!="CWmean"),], pch=16, xaxt="n", ylab="", ylim=c(-0.7,1), xlab="", xlim=c(0.6,6.4))
abline(h=0)
points(value~xvals, avglong[which(avglong$trait=="Nmass" & avglong$SP.ID=="CWmean"),], pch=24, bg=CWbg)
mtext(text = "F)",side = 3,line = panlab.ln,adj=0.02,cex=panlab.cex )
mtext(text = "log(Nmass)",side = 3,line = panlab.ln,adj=0.9,cex=panlab.cex)
mtext(text="Standardized Effect Sizes",side = 2, adj=-.1, line=2.1)
text(x=7, y=0, pos=3, labels = "NA",cex=.7)

boxplot(value~xvals, avglong[which(avglong$trait=="Narea" & avglong$SP.ID!="CWmean"),], at=c(1,2,3,4,5,6)
        , ylim=c(-.7,1), xaxt="n", xlim=c(0.5,7.5)
        ,outpch=1, border="grey", col="grey", outcol="grey"
        , staplewex=0, notch=F, medcol="black" #, medlwd=0
        , whisklty=1, whisklwd=2, boxwex=.7, whiskcol="grey" #, boxwex=1, medlwd=3
        ,las=2)
#plot(value~xvals, avglong[which(avglong$trait=="Narea" & avglong$SP.ID!="CWmean"),], pch=16, xaxt="n", ylab="", ylim=c(-0.7,1), xlab="", xlim=c(0.6,6.4))
abline(h=0)
points(value~xvals, avglong[which(avglong$trait=="Narea" & avglong$SP.ID=="CWmean"),], pch=24, bg=CWbg)
mtext(text = "H)",side = 3,line = panlab.ln,adj=0.02,cex=panlab.cex )
mtext(text = "log(Narea)",side = 3,line = panlab.ln,adj=0.9,cex=panlab.cex)
text(x=7, y=0, pos=3, labels = "NA",cex=.7)

axis(side = 1,at = c(1,2,3,4,5,6,7), labels = c("Wetness", "Warmth","Soil N", "log(Age)","LAI","Gr Rate", "Ecoreg"), xpd=NA, las=3, srt=50, xlim=c(0.6,6.4))


### Variable Importances #####

boxplot(value~xvals, impslong[which(impslong$trait=="LeafLife" & impslong$SP.ID!="CWmean"),], at=c(1,2,3,4,5,6,7)
        , ylim=c(0,1.3), yaxt="n", xaxt="n"
        ,outpch=1, border="grey", col="grey", outcol="grey"
        , staplewex=0, notch=F, medcol="black" #, medlwd=0
        , whisklty=1, whisklwd=2, boxwex=.7, whiskcol="grey" #, boxwex=1, medlwd=3
        ,las=2)
points(value~xvals, impslong[which(impslong$trait=="LeafLife" & impslong$SP.ID=="CWmean"),], pch=24, bg=CWbg)
axis(2,at=c(0,.25,.5,.75,1), labels = c(0,.25,.5,.75,1))
mtext(text = "C)",side = 3,line = panlab.ln,adj=0.02,cex=panlab.cex )

boxplot(value~xvals, impslong[which(impslong$trait=="LMA" & impslong$SP.ID!="CWmean"),], at=c(1,2,3,4,5,6,7)
        , ylim=c(0,1.3), yaxt="n", xaxt="n"
        ,outpch=1, border="grey", col="grey", outcol="grey"
        , staplewex=0, notch=F, medcol="black" #, medlwd=0
        , whisklty=1, whisklwd=2, boxwex=.7, whiskcol="grey" #, boxwex=1, medlwd=3
        ,las=2)
points(value~xvals, impslong[which(impslong$trait=="LMA" & impslong$SP.ID=="CWmean"),], pch=24, bg=CWbg)
axis(2,at=c(0,.25,.5,.75,1), labels = c(0,.25,.5,.75,1))
mtext(text = "E)",side = 3,line = panlab.ln,adj=0.02,cex=panlab.cex )

boxplot(value~xvals, impslong[which(impslong$trait=="Nmass" & impslong$SP.ID!="CWmean"),], at=c(1,2,3,4,5,6,7)
        , ylim=c(0,1.3), yaxt="n", xaxt="n"
        ,outpch=1, border="grey", col="grey", outcol="grey"
        , staplewex=0, notch=F, medcol="black" #, medlwd=0
        , whisklty=1, whisklwd=2, boxwex=.7, whiskcol="grey" #, boxwex=1, medlwd=3
        ,las=2)
points(value~xvals, impslong[which(impslong$trait=="Nmass" & impslong$SP.ID=="CWmean"),], pch=24, bg=CWbg)
axis(2,at=c(0,.25,.5,.75,1), labels = c(0,.25,.5,.75,1))
mtext(text = "G)",side = 3,line = panlab.ln,adj=0.02,cex=panlab.cex )
mtext(text="Vairable Importances",side = 2, adj=-.2, line=2.1)

boxplot(value~xvals, impslong[which(impslong$trait=="Narea" & impslong$SP.ID!="CWmean"),], at=c(1,2,3,4,5,6,7)
        , ylim=c(0,1.3), yaxt="n", xaxt="n"
        ,outpch=1, border="grey", col="grey", outcol="grey"
        , staplewex=0, notch=F, medcol="black" #, medlwd=0
        , whisklty=1, whisklwd=2, boxwex=.7, whiskcol="grey" #, boxwex=1, medlwd=3
        ,las=2)
points(value~xvals, impslong[which(impslong$trait=="Narea" & impslong$SP.ID=="CWmean"),], pch=24, bg=CWbg)
axis(2,at=c(0,.25,.5,.75,1), labels = c(0,.25,.5,.75,1))
mtext(text = "I)",side = 3,line = panlab.ln,adj=0.02,cex=panlab.cex )
axis(side = 1,at = c(1,2,3,4,5,6,7), labels = c("Wetness", "Warmth","Soil N", "log(Age)","LAI","Gr Rate","Ecoreg"), xpd=NA, las=3, srt=50, xlim=c(0.6,6.4))



###### Old version with dotplot
# plot(value~xvals, impslong[which(impslong$trait=="LeafLife"),], pch=16, xaxt="n", ylab="", ylim=c(0,1.3), xlab="", xlim=c(0.6,7.4))
# #abline(h=0)
# mtext(text = "E)",side = 3,line = panlab.ln,adj=0.02,cex=panlab.cex )
# #mtext(text = "log(LL)",side = 3,line = panlab.ln,adj=0.9,cex=panlab.cex)
# plot(value~xvals, impslong[which(impslong$trait=="LMA"),], pch=16, xaxt="n", ylab="", ylim=c(0,1.3), xlab="", xlim=c(0.6,7.4))
# #abline(h=0)
# mtext(text = "F)",side = 3,line = panlab.ln,adj=0.02,cex=panlab.cex )
# #mtext(text = "log(LMA)",side = 3,line = panlab.ln,adj=0.9,cex=panlab.cex)
# plot(value~xvals, impslong[which(impslong$trait=="Nmass"),], pch=16, xaxt="n", ylab="", ylim=c(0,1.3), xlab="", xlim=c(0.6,7.4))
# #abline(h=0)
# mtext(text = "G)",side = 3,line = panlab.ln,adj=0.02,cex=panlab.cex )
# #mtext(text = "log(Nmass)",side = 3,line = panlab.ln,adj=0.9,cex=panlab.cex)
# mtext(text="Standardized Effect Sizes",side = 2, adj=-.1, line=2.1)
# plot(value~xvals, impslong[which(impslong$trait=="Narea"),], pch=16, xaxt="n", ylab="", ylim=c(0,1.3), xlab="", xlim=c(0.6,7.4))
# #abline(h=0)
# mtext(text = "H)",side = 3,line = panlab.ln,adj=0.02,cex=panlab.cex )
# #mtext(text = "log(Narea)",side = 3,line = panlab.ln,adj=0.9,cex=panlab.cex)






#________________________________________________________
#### figure of R2s for the different traits



quartz(width=4.33, height=2)
par(mar=c(.2,5,3,3), mgp=c(2,.75,0))
boxplot(llbestmods$r.squaredGLMM.R2m[-which(llbestmods$SP.ID=="CWmean")]~rep(1, times=6), horizontal=T, at=4, xlim=c(0,5), ylim=c(0,1)
        , yaxt="n", xaxt="n"
        ,outpch=1, border="grey", col="grey", outcol="grey"
        , staplewex=0, notch=F, medcol="black" #, medlwd=0
        , whisklty=1, whisklwd=2, boxwex=.7, whiskcol="grey" )
boxplot(lmabestmods$r.squaredGLMM.R2m[-which(lmabestmods$SP.ID=="CWmean")]~rep(1, times=6), horizontal=T, at=3, xlim=c(0,5), ylim=c(0,1)
        , yaxt="n", xaxt="n"
        ,outpch=1, border="grey", col="grey", outcol="grey"
        , staplewex=0, notch=F, medcol="black" #, medlwd=0
        , whisklty=1, whisklwd=2, boxwex=.7, whiskcol="grey", add=T )
boxplot(nmassbestmods$r.squaredGLMM.R2m[-which(nmassbestmods$SP.ID=="CWmean")]~rep(1, times=6), horizontal=T, at=2, xlim=c(0,5), ylim=c(0,1)
        , yaxt="n", xaxt="n"
        ,outpch=1, border="grey", col="grey", outcol="grey"
        , staplewex=0, notch=F, medcol="black" #, medlwd=0
        , whisklty=1, whisklwd=2, boxwex=.7, whiskcol="grey", add=T )
boxplot(llbestmods$r.squaredGLMM.R2m[-which(llbestmods$SP.ID=="CWmean")]~rep(1, times=6), horizontal=T, at=1, xlim=c(0,5), ylim=c(0,1)
        , yaxt="n", xaxt="n"
        ,outpch=1, border="grey", col="grey", outcol="grey"
        , staplewex=0, notch=F, medcol="black" #, medlwd=0
        , whisklty=1, whisklwd=2, boxwex=.7, whiskcol="grey", add=T )


#plot(llbestmods$r.squaredGLMM.R2m[-which(llbestmods$SP.ID=="CWmean")]~rep(1, times=6), ylim=c(0,1), xlim=c(0,5), xaxt="n", ylab=expression(paste(R^2," or marginal ",R^2)), pch=16)
#points(lmabestmods$r.squaredGLMM.R2m[-which(lmabestmods$SP.ID=="CWmean")]~rep(2, times=6), pch=16)
#points(nmassbestmods$r.squaredGLMM.R2m[-which(nmassbestmods$SP.ID=="CWmean")]~rep(3, times=6), pch=16)
#points(nareabestmods$r.squaredGLMM.R2m[-which(nareabestmods$SP.ID=="CWmean")]~rep(4, times=6), pch=16)
points(c(4)~llbestmods$r.squaredGLMM.R2m[which(llbestmods$SP.ID=="CWmean")], pch=24,bg="white")
points(c(3)~lmabestmods$r.squaredGLMM.R2m[which(lmabestmods$SP.ID=="CWmean")], pch=24,bg="white")
points(c(2)~nmassbestmods$r.squaredGLMM.R2m[which(nmassbestmods$SP.ID=="CWmean")], pch=24,bg="white")
points(c(1)~nareabestmods$r.squaredGLMM.R2m[which(nareabestmods$SP.ID=="CWmean")], pch=24,bg="white")
axis(3)
axis(2, labels = c("Narea","Nmass","LMA","Leaf Life"), at=c(1,2,3,4), las=2)
mtext(side=3,text = expression(paste(R^2," or marginal ", R^2)), line=1.7)
legend(x=.675, y=1, legend = "Ind. Spp", fill="gray", bty="n",col = "grey",border = "grey", cex=.8)
legend(x=.7, y=1.5, legend = "CW mean", bg="white", bty="n",pch=24, cex=.8)
mtext("A)", cex=panlab.cex, side=3, adj=0.02, line=panlab.ln)



########### Other SPP without large enough sample sizes #######################
# ### PICSIT and ABIAMA are right on the edge of useful ###
# #### yeeee, not sure I trust these.
# PICSITll <- trait.mods(traitdata = traits.common, species = "PICSIT",trait = "log.LL")
# PICSITlma <- trait.mods(traitdata = traits.common, species = "PICSIT",trait = "log.LMA")
# PICSITnarea <- trait.mods(traitdata = traits.common, species = "PICSIT",trait = "log.Narea")
# 
# 
# ABIAMAll <- trait.mods(traitdata = traits.common, species = "ABIAMA",trait = "log.LL")
# ABIAMAlma <- trait.mods(traitdata = traits.common, species = "ABIAMA",trait = "log.LMA")
# ABIAMAnarea <- trait.mods(traitdata = traits.common, species = "ABIAMA",trait = "log.Narea")
# 
# 
# 
# ABIGRAll <- trait.mods(traitdata = traits.common, species = "ABIGRA",trait = "log.LL")
# ABIGRAlma <- trait.mods(traitdata = traits.common, species = "ABIGRA",trait = "log.LMA")
# ABIGRAnarea <- trait.mods(traitdata = traits.common, species = "ABIGRA",trait = "log.Narea")
# # n= 12... and the main troubles are soil_N, ASA, log.Narea
# # PLOT_ID     log.LL    climPC1    climPC2     soil_N        ASA      LAI_O AG_TGROWTH 
# #     0          5          0          0         51         34          8          8 
# 
# # JUNOCCll <- trait.mods(traitdata = traits.common, species = "JUNOCC",trait = "log.LL")
# JUNOCClma <- trait.mods(traitdata = traits.common, species = "JUNOCC",trait = "log.LMA")
# JUNOCCnarea <- trait.mods(traitdata = traits.common, species = "JUNOCC",trait = "log.Narea")
# # damn. missing lots of rows. NAs:
# # PLOT_ID     log.LL    climPC1    climPC2     soil_N        ASA      LAI_O AG_TGROWTH 
# # 0         68          0          0         55         10          0          0 
# 
# 
# # ll <- lmer(log.LL~climPC1sc + climPC2sc + soil_Nsc + log.ASAsc + LAI_Osc + AG_TGROWTHsc + (1|PLOT_ID), datazall, REML=F)
# # # lllm <- lm(log.LL~climPC1sc + climPC2sc + soil_Nsc + scale(log(ASA)) + LAI_Osc + AG_TGROWTHsc , datazall)
# # # AIC(ll1,ll2,ll3,ll4)
# # lldredge <- dredge(ll, extra=list(r.squaredGLMM))
# # llmods <- dredge(ll, evaluate = F)
# # llmodavg <- model.avg(lldredge,subset=cumsum(weight)<=.90)
# # lma <- lmer(log.LMA~climPC1 + climPC2 + soil_N + scale(log(ASA)) + LAI_O + AG_TGROWTH + (1|PLOT_ID), dataz, REML=F)
# # lmadredge <- dredge(lma,extra=list(r.squaredGLMM) )
# # lmamodavg <- model.avg(lmadredge,subset=cumsum(weight)<=.90)
# # Narea1 <- lmer(log.Narea~climPC1 + climPC2 + soil_N + scale(log(ASA)) + LAI_O + AG_TGROWTH + (1|PLOT_ID), dataz, REML=F)
# # Nareadredge <- dredge(Narea1, extra=list(r.squaredGLMM))
# # Nareamodavg <- model.avg(Nareadredge, subset=cumsum(weight)<=.90)
# 
# 
# 




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
######### Ecoregion Analysis ###############
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++


ggplot(traits.common[which(traits.common$GE.SP %in% c("Pseudotsuga.menziesii","Pinus.ponderosa")),], aes(x=ECOREGION, y=log.LL, col=ECOREGION)) + geom_boxplot() + facet_grid(GE.SP~.)
# it actually looks like ecoregion might explain a fair amount of LL variation, and maybe Nmass
# though the effect is not consistent between ponderosa and doug fir, maybe data driven param but not mechanistic?
plot(log.LMA~log.Nmass, traits.common[which(traits.common$GE.SP %in% c("Pseudotsuga.menziesii")),], col=ECOREGION, pch=16)
# also, ecoregions cluster rather suprisingly in log.LMA vs log.Nmass space...
# good for doug fir, less good for ponderosa




