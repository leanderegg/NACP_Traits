
#_______________________________________________________________________
############### ****Trait - environment relationships**** ##############
#_______________________________________________________________________

# function for calculating to coefficient of variation for data with NAs
CV <- function(x){ return(sd(as.numeric(x), na.rm=T)/mean(as.numeric(x),na.rm=T))}

#____________________
####### Old Data Exploration and Visualization #############
#____________________

#### .LMA vs LL in PSME vs precip #####
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




##### .LES LLvLMA relationship as f(MAP) ######


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




### .LL for species other than PSME ######

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


########## .Full mixed-effect analysis of LL ~f(LMA * MAP) effects ############
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




######### .Pairplots of traits~Env Variables ##########
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

################## END: Old Exploratory Code ############################################
#__________________________________________________________________________________________






#_______________________________________________________________________
####################### ******Model Selection/Ensambling ********** ####################
#### Modeling the trait-climate relationships 
#________________________________________________________________________

### Summarize the number of branches, plots and ecoregions available for common species
  # PINCON, PINJEF, TSUHET, ABICON, PINPON,PSEMEN generally have >30 branches from >15 plots
trait <- "log.LL"
sLL <- traits.common %>% select(GE.SP, PLOT_ID, get(trait), climPC1, climPC2, soil_N, ASA, LAI_O, AG_TGROWTH,ECOREGION) %>% filter(complete.cases(.)) %>% group_by(GE.SP) %>% summarise( n= n(), nplots = length(unique(PLOT_ID)), nEcos = length(unique(ECOREGION))) %>% filter(n>15 & nplots>5) %>% arrange(n)
trait <- "log.LMA"
sLMA  <- traits.common %>% select(GE.SP, PLOT_ID, get(trait), climPC1, climPC2, soil_N, ASA, LAI_O, AG_TGROWTH,ECOREGION) %>% filter(complete.cases(.)) %>% group_by(GE.SP) %>% summarise( n= n(), nplots = length(unique(PLOT_ID)), nEcos = length(unique(ECOREGION))) %>% filter(n>15 & nplots>5)%>% arrange(n)
trait <- "log.Narea"
sNarea  <- traits.common %>% select(GE.SP, PLOT_ID, get(trait), climPC1, climPC2, soil_N, ASA, LAI_O, AG_TGROWTH,ECOREGION) %>% filter(complete.cases(.)) %>% group_by(GE.SP) %>% summarise( n= n(), nplots = length(unique(PLOT_ID)), nEcos = length(unique(ECOREGION))) %>% filter(n>15 & nplots>5)%>% arrange(n)
trait <- "log.Nmass"
sNmass  <- traits.common %>% select(GE.SP, PLOT_ID, get(trait), climPC1, climPC2, soil_N, ASA, LAI_O, AG_TGROWTH,ECOREGION) %>% filter(complete.cases(.)) %>% group_by(GE.SP) %>% summarise( n= n(), nplots = length(unique(PLOT_ID)), nEcos = length(unique(ECOREGION))) %>% filter(n>15 & nplots>5)%>% arrange(n)

# including pH massively decreases my sample size. Completely lose ABICON, PINCON, and PINJEF. Only TSUHET and PSEMEN have >30 branches from >15 plots.
sLLph <- traits.common %>% select(GE.SP, PLOT_ID, get(trait), climPC1, climPC2, soil_N, soil_pH, ASA, LAI_O, AG_TGROWTH,ECOREGION) %>% filter(complete.cases(.)) %>% group_by(GE.SP) %>% summarise( n= n(), nplots = length(unique(PLOT_ID)), nEcos = length(unique(ECOREGION))) %>% filter(n>15 & nplots>5) %>% arrange(n)

### Define functions for automated model ensembling
require(MuMIn)
options(na.action = "na.fail")

## analysis with ecoregion
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
 # if(length(unique(dataz$ECOREGION))==1){
 #    ecoreg <- lmer(get(traitsc)~climPC1sc + climPC2sc + soil_Nsc + log.ASAsc + LAI_Osc + AG_TGROWTHsc + (1|PLOT_ID), datazall, REML=F)
  #}
  # note: this bit of code only commented out on 8.28. may have fucked all this shit up previously...
  # else{
  #   traitmod <- lmer(get(traitsc)~climPC1sc + climPC2sc + soil_Nsc + log.ASAsc + LAI_Osc + AG_TGROWTHsc + ECOREGION + (1|PLOT_ID), datazall, REML=F)
  # }
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


####### Function for fitting models without ECOREGION + soil pH ######
trait.mods.ne.ph <- function(traitdata =traits.common, species, trait="log.LL", modcrit=F ){
  dataz <- data.frame(traitdata) %>% filter(SP.ID==species) %>% select(PLOT_ID, ECOREGION, get(trait), climPC1, climPC2, soil_N, soil_pH, ASA, LAI_O, AG_TGROWTH) %>% filter(complete.cases(.))
  print(nrow(dataz))
  dataz$log.ASA <- log(dataz$ASA)
  datazsc <- scale(dataz[,-c(1:2)]) # originally, this just scaled the predictors, but I think I also want to scale the traits.
  colnames(datazsc) <- paste0(colnames(dataz[,-c(1:2)]),'sc')
  datazall <- data.frame(dataz, datazsc)
  traitsc <- paste0(trait, "sc")
  traitmod <- lmer(get(traitsc)~climPC1sc + climPC2sc + soil_Nsc + soil_pHsc + log.ASAsc + LAI_Osc + AG_TGROWTHsc + (1|PLOT_ID), datazall, REML=F)
  traitdredge <- dredge(traitmod, extra=list(r.squaredGLMM))
  traitmodcalls <- dredge(traitmod, evaluate = F)
  #  traitmodavg <- model.avg(traitdredge,subset=cumsum(weight)<=.90)
  traitmodavg <- model.avg(traitdredge,subset=delta<=4)
  # if(length(unique(dataz$ECOREGION))==1){
  #   ecoreg <- lmer(get(traitsc)~climPC1sc + climPC2sc + soil_Nsc + log.ASAsc + LAI_Osc + AG_TGROWTHsc + (1|PLOT_ID), datazall, REML=F)
  # }
  # else{
  #   traitmod <- lmer(get(traitsc)~climPC1sc + climPC2sc + soil_Nsc + log.ASAsc + LAI_Osc + AG_TGROWTHsc + ECOREGION + (1|PLOT_ID), datazall, REML=F)
  # }
  bestmodobject <- eval(traitmodcalls[[rownames(traitdredge)[1]]])
  if(modcrit==T){
    scatter.smooth(resid(bestmodobject)~fitted(bestmodobject)); abline(h=0)
    plot(datazall[,traitsc]~predict(bestmodobject, re.form=NA));abline(a=0,b=1)
    qqp(resid(bestmodobject), main="residuals")
    qqp(ranef(bestmodobject)[[1]][,1], main='Random Effects')
  }
  results <- list(best = traitdredge[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","soil_pHsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc","r.squaredGLMM.R2m","r.squaredGLMM.R2c",'AICc')]
                  , avg = traitmodavg$coefficients[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","soil_pHsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc")]
                  , imp = as.vector(traitmodavg$importance)[match(c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","soil_pHsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc"),names(traitmodavg$importance))]
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







################ Fitting models with Ecoregion ###########################

##### .W/in Spp models for 6 conifers #####
PSEMENll <- trait.mods(traitdata = traits.common, species = "PSEMEN",trait = "log.LL", modcrit=T)
PSEMENlma <- trait.mods(traitdata = traits.common, species = "PSEMEN",trait = "log.LMA", modcrit=T)
PSEMENnarea <- trait.mods(traitdata = traits.common, species = "PSEMEN",trait = "log.Narea", modcrit=T)
PSEMENnmass <- trait.mods(traitdata = traits.common, species = "PSEMEN",trait = "log.Nmass", modcrit=T)

# --- testing out whether results are similar with Raw trait values ---
#PSEMENllraw <- trait.mods(traitdata = traits.common, species = "PSEMEN",trait = "LLmonths", modcrit = T)
# raw values almost identical, variable rankings very close
# PSEMENlmaraw <- trait.mods(traitdata = traits.common, species = "PSEMEN",trait = "LMA", modcrit=T)
# raw values less predictable, and resids less normal, variable rankings v similar
# PSEMENnarearaw <- trait.mods(traitdata = traits.common, species = "PSEMEN",trait = "Narea", modcrit=T)
# raw values pretty similar R2, variable rankings nearly identical
  # so results aren't super sensitive to transformation

PINPONll <- trait.mods(traitdata = traits.common, species = "PINPON",trait = "log.LL")
PINPONlma <- trait.mods(traitdata = traits.common, species = "PINPON",trait = "log.LMA")
PINPONnarea <- trait.mods(traitdata = traits.common, species = "PINPON",trait = "log.Narea")
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

# Except for LL, Arbutus is almost worth including (16 branches from 11 sites)
#ARBMENll <- trait.mods(traitdata = traits.common, species = "ARBMEN",trait = "log.LL")
ARBMENlma <- trait.mods(traitdata = traits.common, species = "ARBMEN",trait = "log.LMA") 
ARBMENnarea <- trait.mods(traitdata = traits.common, species = "ARBMEN",trait = "log.Narea")
ARBMENnmass <- trait.mods(traitdata = traits.common, species = "ARBMEN",trait = "log.Nmass")


######## .Double checking w/in spp with only LAI <4 stands ###############
  # This restriction knocks out three species, but remaining species look pretty similar
PSEMENll.ne.sun <- trait.mods.ne(traitdata = traits.common.sun, species = "PSEMEN",trait = "log.LL", modcrit=T)
PINPONll.ne.sun <- trait.mods.ne(traitdata = traits.common.sun, species = "PINPON",trait = "log.LL")
ABICONll.ne.sun <- trait.mods.ne(traitdata = traits.common.sun, species = "ABICON",trait = "log.LL")
#TSUHETll.ne.sun <- trait.mods(traitdata = traits.common.sun, species = "TSUHET",trait = "log.LL")
#PINCONll.ne.sun <- trait.mods(traitdata = traits.common.sun, species = "PINCON",trait = "log.LL")
#PINJEFll.ne.sun <- trait.mods(traitdata = traits.common.sun, species = "PINJEF",trait = "log.LL")




###### .Apply to CWmeans  ###########

## deleted the function (moved to old code at end), because something never quite worked and I got tired of futzing

####### .CW Leaf Lifespan ######
trait <- "log.cw_LLp_if"
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
trait="log.cw_LMAp_if"
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


#### .CWM Nmass ####
trait="log.cw_Nmassp_if"
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


###### .CWM Narea ###########
trait <- "log.cw_Nareap_if"
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




# rename the R2 variable I care about to align with w/in spp models for later concatenation
  # note: I don't know why the R2 from r.squaredGLMM differs from the R2 calculated by dredge for non mixed-models
  # they're generally very close but not identical. I will use th R2 from dredge because I trust it more.
colnames(CWll$best)[10] <- "r.squaredGLMM.R2c"
colnames(CWlma$best)[10] <- "r.squaredGLMM.R2c"
colnames(CWnmass$best)[10] <- "r.squaredGLMM.R2c"
colnames(CWnarea$best)[10] <- "r.squaredGLMM.R2c"











#___________________________________________________________________________________
############ Fitting all models without ecoregions #############
#___________________________________________________________________________________

#######. W/in Spp models, NE ##################
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




####### .CW traits without ECOREGION ##########

####### .CW Leaf Lifespan, NE ######
trait <- "log.cw_LLp_if"
dataz <- data.frame(biomass)  %>% select(PLOT_ID, ECOREGION, get(trait), climPC1, climPC2, soil_N, ASA, LAI_O, AG_PROD_TREE_TOTAL_AS_CARBON) %>% filter(complete.cases(.))
print(nrow(dataz))
dataz$log.ASA <- log(dataz$ASA)
datazsc <- scale(dataz[,-c(1:2)]) # originally, this just scaled the predictors, but I think I also want to scale the traits.
colnames(datazsc) <- paste0(colnames(dataz[,-c(1:2)]),'sc')
colnames(datazsc)[grep("AG_", colnames(datazsc))] <- "AG_TGROWTHsc" # change this column name for compatibility
datazall <- data.frame(dataz, datazsc)
tn <- trait
traitsc <- paste(tn, "sc", sep="")
traitmod <- lm(get(traitsc)~climPC1sc + climPC2sc + soil_Nsc + log.ASAsc + LAI_Osc + AG_TGROWTHsc , datazall)
traitdredge <- dredge(traitmod, extra=list(r.squaredGLMM, "R^2"))
traitmodcalls <- dredge(traitmod, evaluate = F)
traitmodavg <- model.avg(traitdredge,subset=delta<=4)
bestmodobject <- eval(traitmodcalls[[rownames(traitdredge)[1]]])
scatter.smooth(resid(bestmodobject)~fitted(bestmodobject)); abline(h=0)
plot(datazall[,traitsc]~predict(bestmodobject, re.form=NA));abline(a=0,b=1)
qqp(resid(bestmodobject), main="residuals")
# qqp(ranef(bestmodobject)[[1]][,1], main='Random Effects')

CWll.ne <- list(best = traitdredge[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc","r.squaredGLMM.R2m","R^2",'AICc')]
             , avg = traitmodavg$coefficients[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc")] # note: If I try to get ecoregion effects, this can blow up. I think I need to leave ECOREGION out of this
             , imp = as.vector(traitmodavg$importance)[match(c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc"),names(traitmodavg$importance))]
             , fulldredge = traitdredge
             , fullmodavg = traitmodavg
             , modnum = rownames(traitdredge)[1]
             , bestmodcall = traitmodcalls[[rownames(traitdredge)[1]]]
             , bestmodobject = bestmodobject
             , n=nrow(datazall)
             , nplots = nrow(datazall)
             , necos = length(unique(datazall$ECOREGION))
             , deltaNULL = traitdredge[which(rownames(traitdredge)==1),"delta"])


####### .CW LMA NE ####
trait="log.cw_LMAp_if"
dataz <- data.frame(biomass)  %>% select(PLOT_ID, ECOREGION, get(trait), climPC1, climPC2, soil_N, ASA, LAI_O, AG_PROD_TREE_TOTAL_AS_CARBON) %>% filter(complete.cases(.))
print(nrow(dataz))
dataz$log.ASA <- log(dataz$ASA)
datazsc <- scale(dataz[,-c(1:2)]) # originally, this just scaled the predictors, but I think I also want to scale the traits.
colnames(datazsc) <- paste0(colnames(dataz[,-c(1:2)]),'sc')
colnames(datazsc)[grep("AG_", colnames(datazsc))] <- "AG_TGROWTHsc" # change this column name for compatibility
datazall <- data.frame(dataz, datazsc)
tn <- trait
traitsc <- paste(tn, "sc", sep="")
traitmod <- lm(get(traitsc)~climPC1sc + climPC2sc + soil_Nsc + log.ASAsc + LAI_Osc + AG_TGROWTHsc , datazall)
traitdredge <- dredge(traitmod, extra=list(r.squaredGLMM, "R^2"))
traitmodcalls <- dredge(traitmod, evaluate = F)
#  traitmodavg <- model.avg(traitdredge,subset=cumsum(weight)<=.90)
traitmodavg <- model.avg(traitdredge,subset=delta<=4)

bestmodobject <- eval(traitmodcalls[[rownames(traitdredge)[1]]])

scatter.smooth(resid(bestmodobject)~fitted(bestmodobject)); abline(h=0)
plot(datazall[,traitsc]~predict(bestmodobject, re.form=NA));abline(a=0,b=1)
qqp(resid(bestmodobject), main="residuals")
# qqp(ranef(bestmodobject)[[1]][,1], main='Random Effects')


CWlma.ne <- list(best = traitdredge[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc","r.squaredGLMM.R2m","R^2",'AICc')]
              , avg = traitmodavg$coefficients[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc")] # note: If I try to get ecoregion effects, this can blow up. I think I need to leave ECOREGION out of this
              , imp = as.vector(traitmodavg$importance)[match(c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc"),names(traitmodavg$importance))]
              , fulldredge = traitdredge
              , fullmodavg = traitmodavg
              , modnum = rownames(traitdredge)[1]
              , bestmodcall = traitmodcalls[[rownames(traitdredge)[1]]]
              , bestmodobject = bestmodobject
              , n=nrow(datazall)
              , nplots = nrow(datazall)
              , necos = length(unique(datazall$ECOREGION))
              , deltaNULL = traitdredge[which(rownames(traitdredge)==1),"delta"])


#### .Nmass  NE####
trait="log.cw_Nmassp_if"
dataz <- data.frame(biomass)  %>% select(PLOT_ID, ECOREGION, get(trait), climPC1, climPC2, soil_N, ASA, LAI_O, AG_PROD_TREE_TOTAL_AS_CARBON) %>% filter(complete.cases(.))
print(nrow(dataz))
dataz$log.ASA <- log(dataz$ASA)
datazsc <- scale(dataz[,-c(1:2)]) # originally, this just scaled the predictors, but I think I also want to scale the traits.
colnames(datazsc) <- paste0(colnames(dataz[,-c(1:2)]),'sc')
colnames(datazsc)[grep("AG_", colnames(datazsc))] <- "AG_TGROWTHsc" # change this column name for compatibility
datazall <- data.frame(dataz, datazsc)
tn <- trait
traitsc <- paste(tn, "sc", sep="")
traitmod <- lm(get(traitsc)~climPC1sc + climPC2sc + soil_Nsc + log.ASAsc + LAI_Osc + AG_TGROWTHsc , datazall)
traitdredge <- dredge(traitmod, extra=list(r.squaredGLMM, "R^2"))
traitmodcalls <- dredge(traitmod, evaluate = F)
#  traitmodavg <- model.avg(traitdredge,subset=cumsum(weight)<=.90)
traitmodavg <- model.avg(traitdredge,subset=delta<=4)

bestmodobject <- eval(traitmodcalls[[rownames(traitdredge)[1]]])

scatter.smooth(resid(bestmodobject)~fitted(bestmodobject)); abline(h=0)
plot(datazall[,traitsc]~predict(bestmodobject, re.form=NA));abline(a=0,b=1)
qqp(resid(bestmodobject), main="residuals")
# qqp(ranef(bestmodobject)[[1]][,1], main='Random Effects')


CWnmass.ne <- list(best = traitdredge[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc","r.squaredGLMM.R2m","R^2",'AICc')]
                , avg = traitmodavg$coefficients[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc")] # note: If I try to get ecoregion effects, this can blow up. I think I need to leave ECOREGION out of this
                , imp = as.vector(traitmodavg$importance)[match(c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc"),names(traitmodavg$importance))]
                , fulldredge = traitdredge
                , fullmodavg = traitmodavg
                , modnum = rownames(traitdredge)[1]
                , bestmodcall = traitmodcalls[[rownames(traitdredge)[1]]]
                , bestmodobject = bestmodobject
                , n=nrow(datazall)
                , nplots = nrow(datazall)
                , necos = length(unique(datazall$ECOREGION))
                , deltaNULL = traitdredge[which(rownames(traitdredge)==1),"delta"])


###### .Narea NE###########
trait <- "log.cw_Nareap_if"
dataz <- data.frame(biomass)  %>% select(PLOT_ID, ECOREGION, get(trait), climPC1, climPC2, soil_N, ASA, LAI_O, AG_PROD_TREE_TOTAL_AS_CARBON) %>% filter(complete.cases(.))
print(nrow(dataz))
dataz$log.ASA <- log(dataz$ASA)
datazsc <- scale(dataz[,-c(1:2)]) # originally, this just scaled the predictors, but I think I also want to scale the traits.
colnames(datazsc) <- paste0(colnames(dataz[,-c(1:2)]),'sc')
colnames(datazsc)[grep("AG_", colnames(datazsc))] <- "AG_TGROWTHsc" # change this column name for compatibility
datazall <- data.frame(dataz, datazsc)
tn <- trait
traitsc <- paste(tn, "sc", sep="")
traitmod <- lm(get(traitsc)~climPC1sc + climPC2sc + soil_Nsc + log.ASAsc + LAI_Osc + AG_TGROWTHsc , datazall)
traitdredge <- dredge(traitmod, extra=list(r.squaredGLMM, "R^2"))
traitmodcalls <- dredge(traitmod, evaluate = F)
#  traitmodavg <- model.avg(traitdredge,subset=cumsum(weight)<=.90)
traitmodavg <- model.avg(traitdredge,subset=delta<=4)

bestmodobject <- eval(traitmodcalls[[rownames(traitdredge)[1]]])

scatter.smooth(resid(bestmodobject)~fitted(bestmodobject)); abline(h=0)
plot(datazall[,traitsc]~predict(bestmodobject, re.form=NA));abline(a=0,b=1)
qqp(resid(bestmodobject), main="residuals")
# qqp(ranef(bestmodobject)[[1]][,1], main='Random Effects')


CWnarea.ne <- list(best = traitdredge[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc","r.squaredGLMM.R2m","R^2",'AICc')]
                , avg = traitmodavg$coefficients[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc")] # note: If I try to get ecoregion effects, this can blow up. I think I need to leave ECOREGION out of this
                , imp = as.vector(traitmodavg$importance)[match(c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc"),names(traitmodavg$importance))]
                , fulldredge = traitdredge
                , fullmodavg = traitmodavg
                , modnum = rownames(traitdredge)[1]
                , bestmodcall = traitmodcalls[[rownames(traitdredge)[1]]]
                , bestmodobject = bestmodobject
                , n=nrow(datazall)
                , nplots = nrow(datazall)
                , necos = length(unique(datazall$ECOREGION))
                , deltaNULL = traitdredge[which(rownames(traitdredge)==1),"delta"])




### rename R2 column calculated by dredge so it aligns with w/in species conditional R2. Will use this R2 rather than the marginal R2 from r.squaredGLMM
colnames(CWll.ne$best)[9] <- "r.squaredGLMM.R2c"
colnames(CWlma.ne$best)[9] <- "r.squaredGLMM.R2c"
colnames(CWnmass.ne$best)[9] <- "r.squaredGLMM.R2c"
colnames(CWnarea.ne$best)[9] <- "r.squaredGLMM.R2c"






###### .Apply to SPPmeans  ###########
  # uses object spp.traits calculated in Cd-CommunityWeightedAverages.R 

####### .SPP Leaf Lifespan NE ######
trait <- "mlog.LL"
dataz <- data.frame(spp.traits)  %>% select( get(trait), climPC1, climPC2, soil_N, ASA, LAI_O, AG_TGROWTH) %>% filter(complete.cases(.))
print(nrow(dataz))
dataz$log.ASA <- log(dataz$ASA)
datazsc <- scale(dataz) # originally, this just scaled the predictors, but I think I also want to scale the traits.
colnames(datazsc) <- paste0(colnames(dataz),'sc')
#colnames(datazsc)[grep("AG_", colnames(datazsc))] <- "AG_TGROWTHsc" # change this column name for compatibility
datazall <- data.frame(dataz, datazsc)
tn <- trait
traitsc <- paste(tn, "sc", sep="")
traitmod <- lm(get(traitsc)~climPC1sc + climPC2sc + soil_Nsc + log.ASAsc + LAI_Osc + AG_TGROWTHsc , datazall)
traitdredge <- dredge(traitmod, extra=list(r.squaredGLMM, "R^2"))
traitmodcalls <- dredge(traitmod, evaluate = F)
traitmodavg <- model.avg(traitdredge,subset=delta<=4)
bestmodobject <- eval(traitmodcalls[[rownames(traitdredge)[1]]])
scatter.smooth(resid(bestmodobject)~fitted(bestmodobject)); abline(h=0)
plot(datazall[,traitsc]~predict(bestmodobject, re.form=NA));abline(a=0,b=1)
qqp(resid(bestmodobject), main="residuals")
# qqp(ranef(bestmodobject)[[1]][,1], main='Random Effects')

SPPll.ne <- list(best = traitdredge[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc","r.squaredGLMM.R2m","R^2",'AICc')]
              , avg = traitmodavg$coefficients[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc")] # note: If I try to get ecoregion effects, this can blow up. I think I need to leave ECOREGION out of this
              , imp = as.vector(traitmodavg$importance)[match(c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc"),names(traitmodavg$importance))]
              , fulldredge = traitdredge
              , fullmodavg = traitmodavg
              , modnum = rownames(traitdredge)[1]
              , bestmodcall = traitmodcalls[[rownames(traitdredge)[1]]]
              , bestmodobject = bestmodobject
              , n=nrow(datazall)
              #             , nplots = nrow(datazall)
              #             , necos = length(unique(datazall$ECOREGION))
              , deltaNULL = traitdredge[which(rownames(traitdredge)==1),"delta"])


####### .SPP LMA NE####
trait="mlog.LMA"
dataz <- data.frame(spp.traits)  %>% select( get(trait), climPC1, climPC2, soil_N, ASA, LAI_O, AG_TGROWTH) %>% filter(complete.cases(.))
print(nrow(dataz))
dataz$log.ASA <- log(dataz$ASA)
datazsc <- scale(dataz) # originally, this just scaled the predictors, but I think I also want to scale the traits.
colnames(datazsc) <- paste0(colnames(dataz),'sc')
#colnames(datazsc)[grep("AG_", colnames(datazsc))] <- "AG_TGROWTHsc" # change this column name for compatibility
datazall <- data.frame(dataz, datazsc)
tn <- trait
traitsc <- paste(tn, "sc", sep="")
traitmod <- lm(get(traitsc)~climPC1sc + climPC2sc + soil_Nsc + log.ASAsc + LAI_Osc + AG_TGROWTHsc , datazall)
traitdredge <- dredge(traitmod, extra=list(r.squaredGLMM, "R^2"))
traitmodcalls <- dredge(traitmod, evaluate = F)
traitmodavg <- model.avg(traitdredge,subset=delta<=4)
bestmodobject <- eval(traitmodcalls[[rownames(traitdredge)[1]]])
scatter.smooth(resid(bestmodobject)~fitted(bestmodobject)); abline(h=0)
plot(datazall[,traitsc]~predict(bestmodobject, re.form=NA));abline(a=0,b=1)
qqp(resid(bestmodobject), main="residuals")
# qqp(ranef(bestmodobject)[[1]][,1], main='Random Effects')

SPPlma.ne <- list(best = traitdredge[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc","r.squaredGLMM.R2m","R^2",'AICc')]
               , avg = traitmodavg$coefficients[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc")] # note: If I try to get ecoregion effects, this can blow up. I think I need to leave ECOREGION out of this
               , imp = as.vector(traitmodavg$importance)[match(c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc"),names(traitmodavg$importance))]
               , fulldredge = traitdredge
               , fullmodavg = traitmodavg
               , modnum = rownames(traitdredge)[1]
               , bestmodcall = traitmodcalls[[rownames(traitdredge)[1]]]
               , bestmodobject = bestmodobject
               , n=nrow(datazall)
               #             , nplots = nrow(datazall)
               #             , necos = length(unique(datazall$ECOREGION))
               , deltaNULL = traitdredge[which(rownames(traitdredge)==1),"delta"])

#### .SPP Nmass NE####
trait="mlog.Nmass"
dataz <- data.frame(spp.traits)  %>% select( get(trait), climPC1, climPC2, soil_N, ASA, LAI_O, AG_TGROWTH) %>% filter(complete.cases(.))
print(nrow(dataz))
dataz$log.ASA <- log(dataz$ASA)
datazsc <- scale(dataz) # originally, this just scaled the predictors, but I think I also want to scale the traits.
colnames(datazsc) <- paste0(colnames(dataz),'sc')
#colnames(datazsc)[grep("AG_", colnames(datazsc))] <- "AG_TGROWTHsc" # change this column name for compatibility
datazall <- data.frame(dataz, datazsc)
tn <- trait
traitsc <- paste(tn, "sc", sep="")
traitmod <- lm(get(traitsc)~climPC1sc + climPC2sc + soil_Nsc + log.ASAsc + LAI_Osc + AG_TGROWTHsc , datazall)
traitdredge <- dredge(traitmod, extra=list(r.squaredGLMM, "R^2"))
traitmodcalls <- dredge(traitmod, evaluate = F)
traitmodavg <- model.avg(traitdredge,subset=delta<=4)
bestmodobject <- eval(traitmodcalls[[rownames(traitdredge)[1]]])
scatter.smooth(resid(bestmodobject)~fitted(bestmodobject)); abline(h=0)
plot(datazall[,traitsc]~predict(bestmodobject, re.form=NA));abline(a=0,b=1)
qqp(resid(bestmodobject), main="residuals")
# qqp(ranef(bestmodobject)[[1]][,1], main='Random Effects')

SPPnmass.ne <- list(best = traitdredge[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc","r.squaredGLMM.R2m","R^2",'AICc')]
                 , avg = traitmodavg$coefficients[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc")] # note: If I try to get ecoregion effects, this can blow up. I think I need to leave ECOREGION out of this
                 , imp = as.vector(traitmodavg$importance)[match(c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc"),names(traitmodavg$importance))]
                 , fulldredge = traitdredge
                 , fullmodavg = traitmodavg
                 , modnum = rownames(traitdredge)[1]
                 , bestmodcall = traitmodcalls[[rownames(traitdredge)[1]]]
                 , bestmodobject = bestmodobject
                 , n=nrow(datazall)
                 #             , nplots = nrow(datazall)
                 #             , necos = length(unique(datazall$ECOREGION))
                 , deltaNULL = traitdredge[which(rownames(traitdredge)==1),"delta"])

###### .SPP Narea NE ###########
trait <- "mlog.Narea"
dataz <- data.frame(spp.traits)  %>% select( get(trait), climPC1, climPC2, soil_N, ASA, LAI_O, AG_TGROWTH) %>% filter(complete.cases(.))
print(nrow(dataz))
dataz$log.ASA <- log(dataz$ASA)
datazsc <- scale(dataz) # originally, this just scaled the predictors, but I think I also want to scale the traits.
colnames(datazsc) <- paste0(colnames(dataz),'sc')
#colnames(datazsc)[grep("AG_", colnames(datazsc))] <- "AG_TGROWTHsc" # change this column name for compatibility
datazall <- data.frame(dataz, datazsc)
tn <- trait
traitsc <- paste(tn, "sc", sep="")
traitmod <- lm(get(traitsc)~climPC1sc + climPC2sc + soil_Nsc + log.ASAsc + LAI_Osc + AG_TGROWTHsc , datazall)
traitdredge <- dredge(traitmod, extra=list(r.squaredGLMM, "R^2"))
traitmodcalls <- dredge(traitmod, evaluate = F)
traitmodavg <- model.avg(traitdredge,subset=delta<=4)
bestmodobject <- eval(traitmodcalls[[rownames(traitdredge)[1]]])
scatter.smooth(resid(bestmodobject)~fitted(bestmodobject)); abline(h=0)
plot(datazall[,traitsc]~predict(bestmodobject, re.form=NA));abline(a=0,b=1)
qqp(resid(bestmodobject), main="residuals")
# qqp(ranef(bestmodobject)[[1]][,1], main='Random Effects')

SPPnarea.ne <- list(best = traitdredge[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc","r.squaredGLMM.R2m","R^2",'AICc')]
                 , avg = traitmodavg$coefficients[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc")] # note: If I try to get ecoregion effects, this can blow up. I think I need to leave ECOREGION out of this
                 , imp = as.vector(traitmodavg$importance)[match(c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc"),names(traitmodavg$importance))]
                 , fulldredge = traitdredge
                 , fullmodavg = traitmodavg
                 , modnum = rownames(traitdredge)[1]
                 , bestmodcall = traitmodcalls[[rownames(traitdredge)[1]]]
                 , bestmodobject = bestmodobject
                 , n=nrow(datazall)
                 #             , nplots = nrow(datazall)
                 #             , necos = length(unique(datazall$ECOREGION))
                 , deltaNULL = traitdredge[which(rownames(traitdredge)==1),"delta"])



###### rename R2 column for later concatenation with w/in spp models
colnames(SPPll.ne$best)[9] <- "r.squaredGLMM.R2c"
colnames(SPPlma.ne$best)[9] <- "r.squaredGLMM.R2c"
colnames(SPPnmass.ne$best)[9] <- "r.squaredGLMM.R2c"
colnames(SPPnarea.ne$best)[9] <- "r.squaredGLMM.R2c"








#___________________________________________________________________________________
############ **Fitting all models without ecoregions + soil_pH #############
#___________________________________________________________________________________


PSEMENll.ne.ph <- trait.mods.ne.ph(traitdata = traits.common, species = "PSEMEN",trait = "log.LL", modcrit=T)
PSEMENlma.ne.ph <- trait.mods.ne.ph(traitdata = traits.common, species = "PSEMEN",trait = "log.LMA", modcrit=T)
PSEMENnarea.ne.ph <- trait.mods.ne.ph(traitdata = traits.common, species = "PSEMEN",trait = "log.Narea", modcrit=T)
PSEMENnmass.ne.ph <- trait.mods.ne.ph(traitdata = traits.common, species = "PSEMEN",trait = "log.Nmass", modcrit=T)

PINPONll.ne.ph <- trait.mods.ne.ph(traitdata = traits.common, species = "PINPON",trait = "log.LL")
PINPONlma.ne.ph <- trait.mods.ne.ph(traitdata = traits.common, species = "PINPON",trait = "log.LMA")
PINPONnarea.ne.ph <- trait.mods.ne.ph(traitdata = traits.common, species = "PINPON",trait = "log.Narea")
PINPONnmass.ne.ph <- trait.mods.ne.ph(traitdata = traits.common, species = "PINPON",trait = "log.Nmass")

ABICONll.ne.ph <- trait.mods.ne.ph(traitdata = traits.common, species = "ABICON",trait = "log.LL")
ABICONlma.ne.ph <- trait.mods.ne.ph(traitdata = traits.common, species = "ABICON",trait = "log.LMA")
ABICONnarea.ne.ph <- trait.mods.ne.ph(traitdata = traits.common, species = "ABICON",trait = "log.Narea")
ABICONnmass.ne.ph <- trait.mods.ne.ph(traitdata = traits.common, species = "ABICON",trait = "log.Nmass")

TSUHETll.ne.ph <- trait.mods.ne.ph(traitdata = traits.common, species = "TSUHET",trait = "log.LL")
TSUHETlma.ne.ph <- trait.mods.ne.ph(traitdata = traits.common, species = "TSUHET",trait = "log.LMA")
TSUHETnarea.ne.ph <- trait.mods.ne.ph(traitdata = traits.common, species = "TSUHET",trait = "log.Narea") # might be one outlier that's leveraging things?
TSUHETnmass.ne.ph <- trait.mods.ne.ph(traitdata = traits.common, species = "TSUHET",trait = "log.Nmass") # might be one outlier that's leveraging things?

PINCONll.ne.ph <- trait.mods.ne.ph(traitdata = traits.common, species = "PINCON",trait = "log.LL")
PINCONlma.ne.ph <- trait.mods.ne.ph(traitdata = traits.common, species = "PINCON",trait = "log.LMA")
PINCONnarea.ne.ph <- trait.mods.ne.ph(traitdata = traits.common, species = "PINCON",trait = "log.Narea")
PINCONnmass.ne.ph <- trait.mods.ne.ph(traitdata = traits.common, species = "PINCON",trait = "log.Nmass")

PINJEFll.ne.ph <- trait.mods.ne.ph(traitdata = traits.common, species = "PINJEF",trait = "log.LL")
PINJEFlma.ne.ph <- trait.mods.ne.ph(traitdata = traits.common, species = "PINJEF",trait = "log.LMA") # two outliers
PINJEFnarea.ne.ph <- trait.mods.ne.ph(traitdata = traits.common, species = "PINJEF",trait = "log.Narea")
PINJEFnmass.ne.ph <- trait.mods.ne.ph(traitdata = traits.common, species = "PINJEF",trait = "log.Nmass")

# #ARBMENll.ne.ph <- trait.mods.ne.ph(traitdata = traits.common, species = "ARBMEN",trait = "log.LL")
# ARBMENlma.ne.ph <- trait.mods.ne.ph(traitdata = traits.common, species = "ARBMEN",trait = "log.LMA") 
# ARBMENnarea.ne.ph <- trait.mods.ne.ph(traitdata = traits.common, species = "ARBMEN",trait = "log.Narea")
# ARBMENnmass.ne.ph <- trait.mods.ne.ph(traitdata = traits.common, species = "ARBMEN",trait = "log.Nmass")




####### CW traits without ECOREGION + pH

####### .CW Leaf Lifespan, NE ph ######
trait <- "log.cw_LLp_if"
dataz <- data.frame(biomass)  %>% select(PLOT_ID, ECOREGION, get(trait), climPC1, climPC2, soil_N, soil_pH, ASA, LAI_O, AG_PROD_TREE_TOTAL_AS_CARBON) %>% filter(complete.cases(.))
print(nrow(dataz))
dataz$log.ASA <- log(dataz$ASA)
datazsc <- scale(dataz[,-c(1:2)]) # originally, this just scaled the predictors, but I think I also want to scale the traits.
colnames(datazsc) <- paste0(colnames(dataz[,-c(1:2)]),'sc')
colnames(datazsc)[grep("AG_", colnames(datazsc))] <- "AG_TGROWTHsc" # change this column name for compatibility
datazall <- data.frame(dataz, datazsc)
tn <- trait
traitsc <- paste(tn, "sc", sep="")
traitmod <- lm(get(traitsc)~climPC1sc + climPC2sc + soil_Nsc + soil_pHsc + log.ASAsc + LAI_Osc + AG_TGROWTHsc , datazall)
traitdredge <- dredge(traitmod, extra=list(r.squaredGLMM, "R^2"))
traitmodcalls <- dredge(traitmod, evaluate = F)
traitmodavg <- model.avg(traitdredge,subset=delta<=4)
bestmodobject <- eval(traitmodcalls[[rownames(traitdredge)[1]]])
scatter.smooth(resid(bestmodobject)~fitted(bestmodobject)); abline(h=0)
plot(datazall[,traitsc]~predict(bestmodobject, re.form=NA));abline(a=0,b=1)
qqp(resid(bestmodobject), main="residuals")
# qqp(ranef(bestmodobject)[[1]][,1], main='Random Effects')

CWll.ne.ph <- list(best = traitdredge[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","soil_pHsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc","r.squaredGLMM.R2m","R^2",'AICc')]
                , avg = traitmodavg$coefficients[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","soil_pHsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc")] # note: If I try to get ecoregion effects, this can blow up. I think I need to leave ECOREGION out of this
                , imp = as.vector(traitmodavg$importance)[match(c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","soil_pHsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc"),names(traitmodavg$importance))]
                , fulldredge = traitdredge
                , fullmodavg = traitmodavg
                , modnum = rownames(traitdredge)[1]
                , bestmodcall = traitmodcalls[[rownames(traitdredge)[1]]]
                , bestmodobject = bestmodobject
                , n=nrow(datazall)
                , nplots = nrow(datazall)
                , necos = length(unique(datazall$ECOREGION))
                , deltaNULL = traitdredge[which(rownames(traitdredge)==1),"delta"])


####### .CW LMA NE ph ####
## NOTE: This is non-replicable, and I don't know why. You actually have to run all the code above for some strange assignment error reason
trait="log.cw_LMAp_if"
dataz <- data.frame(biomass)  %>% select(PLOT_ID, ECOREGION, get(trait), climPC1, climPC2, soil_N, soil_pH, ASA, LAI_O, AG_PROD_TREE_TOTAL_AS_CARBON) %>% filter(complete.cases(.))
print(nrow(dataz))
dataz$log.ASA <- log(dataz$ASA)
datazsc <- scale(dataz[,-c(1:2)]) # originally, this just scaled the predictors, but I think I also want to scale the traits.
colnames(datazsc) <- paste0(colnames(dataz[,-c(1:2)]),'sc')
colnames(datazsc)[grep("AG_", colnames(datazsc))] <- "AG_TGROWTHsc" # change this column name for compatibility
datazall <- data.frame(dataz, datazsc)
tn <- trait
traitsc <- paste(tn, "sc", sep="")
traitmod <- lm(get(traitsc)~climPC1sc + climPC2sc + soil_Nsc + soil_pHsc + log.ASAsc + LAI_Osc + AG_TGROWTHsc , datazall)
traitdredge <- dredge(traitmod, extra=list(r.squaredGLMM, "R^2"))
traitmodcalls <- dredge(traitmod, evaluate = F)
traitmodavg <- model.avg(traitdredge,subset=delta<=4)
bestmodobject <- eval(traitmodcalls[[rownames(traitdredge)[1]]])
scatter.smooth(resid(bestmodobject)~fitted(bestmodobject)); abline(h=0)
plot(datazall[,traitsc]~predict(bestmodobject, re.form=NA));abline(a=0,b=1)
qqp(resid(bestmodobject), main="residuals")
# qqp(ranef(bestmodobject)[[1]][,1], main='Random Effects')

CWlma.ne.ph <- list(best = traitdredge[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","soil_pHsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc","r.squaredGLMM.R2m","R^2",'AICc')]
                   , avg = traitmodavg$coefficients[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","soil_pHsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc")] # note: If I try to get ecoregion effects, this can blow up. I think I need to leave ECOREGION out of this
                   , imp = as.vector(traitmodavg$importance)[match(c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","soil_pHsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc"),names(traitmodavg$importance))]
                   , fulldredge = traitdredge
                 , fullmodavg = traitmodavg
                 , modnum = rownames(traitdredge)[1]
                 , bestmodcall = traitmodcalls[[rownames(traitdredge)[1]]]
                 , bestmodobject = bestmodobject
                 , n=nrow(datazall)
                 , nplots = nrow(datazall)
                 , necos = length(unique(datazall$ECOREGION))
                 , deltaNULL = traitdredge[which(rownames(traitdredge)==1),"delta"])


#### .Nmass  NE ph ####
trait="log.cw_Nmassp_if"
dataz <- data.frame(biomass)  %>% select(PLOT_ID, ECOREGION, get(trait), climPC1, climPC2, soil_N, soil_pH, ASA, LAI_O, AG_PROD_TREE_TOTAL_AS_CARBON) %>% filter(complete.cases(.))
print(nrow(dataz))
dataz$log.ASA <- log(dataz$ASA)
datazsc <- scale(dataz[,-c(1:2)]) # originally, this just scaled the predictors, but I think I also want to scale the traits.
colnames(datazsc) <- paste0(colnames(dataz[,-c(1:2)]),'sc')
colnames(datazsc)[grep("AG_", colnames(datazsc))] <- "AG_TGROWTHsc" # change this column name for compatibility
datazall <- data.frame(dataz, datazsc)
tn <- trait
traitsc <- paste(tn, "sc", sep="")
traitmod <- lm(get(traitsc)~climPC1sc + climPC2sc + soil_Nsc + soil_pHsc + log.ASAsc + LAI_Osc + AG_TGROWTHsc , datazall)
traitdredge <- dredge(traitmod, extra=list(r.squaredGLMM, "R^2"))
traitmodcalls <- dredge(traitmod, evaluate = F)
traitmodavg <- model.avg(traitdredge,subset=delta<=4)
bestmodobject <- eval(traitmodcalls[[rownames(traitdredge)[1]]])
scatter.smooth(resid(bestmodobject)~fitted(bestmodobject)); abline(h=0)
plot(datazall[,traitsc]~predict(bestmodobject, re.form=NA));abline(a=0,b=1)
qqp(resid(bestmodobject), main="residuals")
# qqp(ranef(bestmodobject)[[1]][,1], main='Random Effects')

CWnmass.ne.ph <- list(best = traitdredge[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","soil_pHsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc","r.squaredGLMM.R2m","R^2",'AICc')]
                   , avg = traitmodavg$coefficients[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","soil_pHsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc")] # note: If I try to get ecoregion effects, this can blow up. I think I need to leave ECOREGION out of this
                   , imp = as.vector(traitmodavg$importance)[match(c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","soil_pHsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc"),names(traitmodavg$importance))]
                   , fulldredge = traitdredge
                   , fullmodavg = traitmodavg
                   , modnum = rownames(traitdredge)[1]
                   , bestmodcall = traitmodcalls[[rownames(traitdredge)[1]]]
                   , bestmodobject = bestmodobject
                   , n=nrow(datazall)
                   , nplots = nrow(datazall)
                   , necos = length(unique(datazall$ECOREGION))
                   , deltaNULL = traitdredge[which(rownames(traitdredge)==1),"delta"])


###### .Narea NE ph ###########
trait <- "log.cw_Nareap_if"
dataz <- data.frame(biomass)  %>% select(PLOT_ID, ECOREGION, get(trait), climPC1, climPC2, soil_N, soil_pH, ASA, LAI_O, AG_PROD_TREE_TOTAL_AS_CARBON) %>% filter(complete.cases(.))
print(nrow(dataz))
dataz$log.ASA <- log(dataz$ASA)
datazsc <- scale(dataz[,-c(1:2)]) # originally, this just scaled the predictors, but I think I also want to scale the traits.
colnames(datazsc) <- paste0(colnames(dataz[,-c(1:2)]),'sc')
colnames(datazsc)[grep("AG_", colnames(datazsc))] <- "AG_TGROWTHsc" # change this column name for compatibility
datazall <- data.frame(dataz, datazsc)
tn <- trait
traitsc <- paste(tn, "sc", sep="")
traitmod <- lm(get(traitsc)~climPC1sc + climPC2sc + soil_Nsc + soil_pHsc + log.ASAsc + LAI_Osc + AG_TGROWTHsc , datazall)
traitdredge <- dredge(traitmod, extra=list(r.squaredGLMM, "R^2"))
traitmodcalls <- dredge(traitmod, evaluate = F)
traitmodavg <- model.avg(traitdredge,subset=delta<=4)
bestmodobject <- eval(traitmodcalls[[rownames(traitdredge)[1]]])
scatter.smooth(resid(bestmodobject)~fitted(bestmodobject)); abline(h=0)
plot(datazall[,traitsc]~predict(bestmodobject, re.form=NA));abline(a=0,b=1)
qqp(resid(bestmodobject), main="residuals")
# qqp(ranef(bestmodobject)[[1]][,1], main='Random Effects')

CWnarea.ne.ph <- list(best = traitdredge[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","soil_pHsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc","r.squaredGLMM.R2m","R^2",'AICc')]
                   , avg = traitmodavg$coefficients[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","soil_pHsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc")] # note: If I try to get ecoregion effects, this can blow up. I think I need to leave ECOREGION out of this
                   , imp = as.vector(traitmodavg$importance)[match(c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","soil_pHsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc"),names(traitmodavg$importance))]
                   , fulldredge = traitdredge
                   , fullmodavg = traitmodavg
                   , modnum = rownames(traitdredge)[1]
                   , bestmodcall = traitmodcalls[[rownames(traitdredge)[1]]]
                   , bestmodobject = bestmodobject
                   , n=nrow(datazall)
                   , nplots = nrow(datazall)
                   , necos = length(unique(datazall$ECOREGION))
                   , deltaNULL = traitdredge[which(rownames(traitdredge)==1),"delta"])




# NOTE: have to switch from 10 to 9 without ecoregion
colnames(CWll.ne.ph$best)[9] <- "r.squaredGLMM.R2c"
colnames(CWlma.ne.ph$best)[9] <- "r.squaredGLMM.R2c"
colnames(CWnmass.ne.ph$best)[9] <- "r.squaredGLMM.R2c"
colnames(CWnarea.ne.ph$best)[9] <- "r.squaredGLMM.R2c"






###### Apply to SPPmeans NE ph ###########
###applying to the cw traits ####
# just copied the function code to taylor to biomass

## deleted the function (moved to old code at end), because something never quite worked and I got tired of futzing

####### .SPP Leaf Lifespan ne ph ######
trait <- "log.LL"
dataz <- data.frame(spp.traits)  %>% select( get(trait), climPC1, climPC2, soil_N, soil_pH, ASA, LAI_O, AG_TGROWTH) %>% filter(complete.cases(.))
print(nrow(dataz))
dataz$log.ASA <- log(dataz$ASA)
datazsc <- scale(dataz) # originally, this just scaled the predictors, but I think I also want to scale the traits.
colnames(datazsc) <- paste0(colnames(dataz),'sc')
#colnames(datazsc)[grep("AG_", colnames(datazsc))] <- "AG_TGROWTHsc" # change this column name for compatibility
datazall <- data.frame(dataz, datazsc)
tn <- trait
traitsc <- paste(tn, "sc", sep="")
traitmod <- lm(get(traitsc)~climPC1sc + climPC2sc + soil_Nsc + soil_pHsc + log.ASAsc + LAI_Osc + AG_TGROWTHsc , datazall)
traitdredge <- dredge(traitmod, extra=list(r.squaredGLMM, "R^2"))
traitmodcalls <- dredge(traitmod, evaluate = F)
traitmodavg <- model.avg(traitdredge,subset=delta<=4)
bestmodobject <- eval(traitmodcalls[[rownames(traitdredge)[1]]])
scatter.smooth(resid(bestmodobject)~fitted(bestmodobject)); abline(h=0)
plot(datazall[,traitsc]~predict(bestmodobject, re.form=NA));abline(a=0,b=1)
qqp(resid(bestmodobject), main="residuals")
# qqp(ranef(bestmodobject)[[1]][,1], main='Random Effects')

SPPll.ne.ph <- list(best = traitdredge[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","soil_pHsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc","r.squaredGLMM.R2m","R^2",'AICc')]
                 , avg = traitmodavg$coefficients[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","soil_pHsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc")] # note: If I try to get ecoregion effects, this can blow up. I think I need to leave ECOREGION out of this
                 , imp = as.vector(traitmodavg$importance)[match(c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","soil_pHsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc"),names(traitmodavg$importance))]
                 , fulldredge = traitdredge
                 , fullmodavg = traitmodavg
                 , modnum = rownames(traitdredge)[1]
                 , bestmodcall = traitmodcalls[[rownames(traitdredge)[1]]]
                 , bestmodobject = bestmodobject
                 , n=nrow(datazall)
                 #             , nplots = nrow(datazall)
                 #             , necos = length(unique(datazall$ECOREGION))
                 , deltaNULL = traitdredge[which(rownames(traitdredge)==1),"delta"])


####### .SPP LMA NE ph####
## NOTE: This is non-replicable, and I don't know why. You actually have to run all the code above for some strange assignment error reason
trait="log.LMA"
dataz <- data.frame(spp.traits)  %>% select( get(trait), climPC1, climPC2, soil_N, ASA, LAI_O, AG_TGROWTH) %>% filter(complete.cases(.))
print(nrow(dataz))
dataz$log.ASA <- log(dataz$ASA)
datazsc <- scale(dataz) # originally, this just scaled the predictors, but I think I also want to scale the traits.
colnames(datazsc) <- paste0(colnames(dataz),'sc')
#colnames(datazsc)[grep("AG_", colnames(datazsc))] <- "AG_TGROWTHsc" # change this column name for compatibility
datazall <- data.frame(dataz, datazsc)
tn <- trait
traitsc <- paste(tn, "sc", sep="")
traitmod <- lm(get(traitsc)~climPC1sc + climPC2sc + soil_Nsc + log.ASAsc + LAI_Osc + AG_TGROWTHsc , datazall)
traitdredge <- dredge(traitmod, extra=list(r.squaredGLMM, "R^2"))
traitmodcalls <- dredge(traitmod, evaluate = F)
traitmodavg <- model.avg(traitdredge,subset=delta<=4)
bestmodobject <- eval(traitmodcalls[[rownames(traitdredge)[1]]])
scatter.smooth(resid(bestmodobject)~fitted(bestmodobject)); abline(h=0)
plot(datazall[,traitsc]~predict(bestmodobject, re.form=NA));abline(a=0,b=1)
qqp(resid(bestmodobject), main="residuals")
# qqp(ranef(bestmodobject)[[1]][,1], main='Random Effects')

SPPlma.ne.ph <- list(best = traitdredge[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc","r.squaredGLMM.R2m","R^2",'AICc')]
                  , avg = traitmodavg$coefficients[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc")] # note: If I try to get ecoregion effects, this can blow up. I think I need to leave ECOREGION out of this
                  , imp = as.vector(traitmodavg$importance)[match(c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc"),names(traitmodavg$importance))]
                  , fulldredge = traitdredge
                  , fullmodavg = traitmodavg
                  , modnum = rownames(traitdredge)[1]
                  , bestmodcall = traitmodcalls[[rownames(traitdredge)[1]]]
                  , bestmodobject = bestmodobject
                  , n=nrow(datazall)
                  #             , nplots = nrow(datazall)
                  #             , necos = length(unique(datazall$ECOREGION))
                  , deltaNULL = traitdredge[which(rownames(traitdredge)==1),"delta"])

#### .SPP Nmass NE ph ####
trait="log.Nmass"
dataz <- data.frame(spp.traits)  %>% select( get(trait), climPC1, climPC2, soil_N, ASA, LAI_O, AG_TGROWTH) %>% filter(complete.cases(.))
print(nrow(dataz))
dataz$log.ASA <- log(dataz$ASA)
datazsc <- scale(dataz) # originally, this just scaled the predictors, but I think I also want to scale the traits.
colnames(datazsc) <- paste0(colnames(dataz),'sc')
#colnames(datazsc)[grep("AG_", colnames(datazsc))] <- "AG_TGROWTHsc" # change this column name for compatibility
datazall <- data.frame(dataz, datazsc)
tn <- trait
traitsc <- paste(tn, "sc", sep="")
traitmod <- lm(get(traitsc)~climPC1sc + climPC2sc + soil_Nsc + log.ASAsc + LAI_Osc + AG_TGROWTHsc , datazall)
traitdredge <- dredge(traitmod, extra=list(r.squaredGLMM, "R^2"))
traitmodcalls <- dredge(traitmod, evaluate = F)
traitmodavg <- model.avg(traitdredge,subset=delta<=4)
bestmodobject <- eval(traitmodcalls[[rownames(traitdredge)[1]]])
scatter.smooth(resid(bestmodobject)~fitted(bestmodobject)); abline(h=0)
plot(datazall[,traitsc]~predict(bestmodobject, re.form=NA));abline(a=0,b=1)
qqp(resid(bestmodobject), main="residuals")
# qqp(ranef(bestmodobject)[[1]][,1], main='Random Effects')

SPPnmass.ne.ph <- list(best = traitdredge[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc","r.squaredGLMM.R2m","R^2",'AICc')]
                    , avg = traitmodavg$coefficients[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc")] # note: If I try to get ecoregion effects, this can blow up. I think I need to leave ECOREGION out of this
                    , imp = as.vector(traitmodavg$importance)[match(c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc"),names(traitmodavg$importance))]
                    , fulldredge = traitdredge
                    , fullmodavg = traitmodavg
                    , modnum = rownames(traitdredge)[1]
                    , bestmodcall = traitmodcalls[[rownames(traitdredge)[1]]]
                    , bestmodobject = bestmodobject
                    , n=nrow(datazall)
                    #             , nplots = nrow(datazall)
                    #             , necos = length(unique(datazall$ECOREGION))
                    , deltaNULL = traitdredge[which(rownames(traitdredge)==1),"delta"])

###### .SPP Narea NE ph ###########
trait <- "log.Narea"
dataz <- data.frame(spp.traits)  %>% select( get(trait), climPC1, climPC2, soil_N, ASA, LAI_O, AG_TGROWTH) %>% filter(complete.cases(.))
print(nrow(dataz))
dataz$log.ASA <- log(dataz$ASA)
datazsc <- scale(dataz) # originally, this just scaled the predictors, but I think I also want to scale the traits.
colnames(datazsc) <- paste0(colnames(dataz),'sc')
#colnames(datazsc)[grep("AG_", colnames(datazsc))] <- "AG_TGROWTHsc" # change this column name for compatibility
datazall <- data.frame(dataz, datazsc)
tn <- trait
traitsc <- paste(tn, "sc", sep="")
traitmod <- lm(get(traitsc)~climPC1sc + climPC2sc + soil_Nsc + log.ASAsc + LAI_Osc + AG_TGROWTHsc , datazall)
traitdredge <- dredge(traitmod, extra=list(r.squaredGLMM, "R^2"))
traitmodcalls <- dredge(traitmod, evaluate = F)
traitmodavg <- model.avg(traitdredge,subset=delta<=4)
bestmodobject <- eval(traitmodcalls[[rownames(traitdredge)[1]]])
scatter.smooth(resid(bestmodobject)~fitted(bestmodobject)); abline(h=0)
plot(datazall[,traitsc]~predict(bestmodobject, re.form=NA));abline(a=0,b=1)
qqp(resid(bestmodobject), main="residuals")
# qqp(ranef(bestmodobject)[[1]][,1], main='Random Effects')

SPPnarea.ne.ph <- list(best = traitdredge[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc","r.squaredGLMM.R2m","R^2",'AICc')]
                    , avg = traitmodavg$coefficients[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc")] # note: If I try to get ecoregion effects, this can blow up. I think I need to leave ECOREGION out of this
                    , imp = as.vector(traitmodavg$importance)[match(c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc"),names(traitmodavg$importance))]
                    , fulldredge = traitdredge
                    , fullmodavg = traitmodavg
                    , modnum = rownames(traitdredge)[1]
                    , bestmodcall = traitmodcalls[[rownames(traitdredge)[1]]]
                    , bestmodobject = bestmodobject
                    , n=nrow(datazall)
                    #             , nplots = nrow(datazall)
                    #             , necos = length(unique(datazall$ECOREGION))
                    , deltaNULL = traitdredge[which(rownames(traitdredge)==1),"delta"])




colnames(SPPll.ne.ph$best)[9] <- "r.squaredGLMM.R2c"
colnames(SPPlma.ne.ph$best)[9] <- "r.squaredGLMM.R2c"
colnames(SPPnmass.ne.ph$best)[9] <- "r.squaredGLMM.R2c"
colnames(SPPnarea.ne.ph$best)[9] <- "r.squaredGLMM.R2c"














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







####### Combining all environmental trait models NE (Table S5)######
  # combine all info on the best models for creating Table S5

llbestmods.ne <- rbind(PSEMENll.ne$best, PINPONll.ne$best,PINCONll.ne$best,PINJEFll.ne$best, ABICONll.ne$best, TSUHETll.ne$best,SPPll.ne$best, CWll.ne$best)
llbestmods.ne$n <- rbind(PSEMENll.ne$n, PINPONll.ne$n,PINCONll.ne$n,PINJEFll.ne$n, ABICONll.ne$n, TSUHETll.ne$n,SPPll.ne$n,CWll.ne$n)
llbestmods.ne$nplots <- rbind(PSEMENll.ne$nplots, PINPONll.ne$nplots,PINCONll.ne$nplots,PINJEFll.ne$nplots, ABICONll.ne$nplots, TSUHETll.ne$nplots, NA, CWll.ne$nplots)
llbestmods.ne$necos <- rbind(PSEMENll.ne$necos, PINPONll.ne$necos,PINCONll.ne$necos,PINJEFll.ne$necos, ABICONll.ne$necos, TSUHETll.ne$necos,NA, CWll.ne$necos)
llbestmods.ne$deltas <-rbind(PSEMENll.ne$deltaNULL, PINPONll.ne$deltaNULL,PINCONll.ne$deltaNULL,PINJEFll.ne$deltaNULL, ABICONll.ne$deltaNULL, TSUHETll.ne$deltaNULL,SPPll.ne$deltaNULL,CWll.ne$deltaNULL)
llbestmods.ne$SP.ID <-  c("PSEMEN","PINPON","PINCON","PINJEF","ABICON","TSUHET", "SPPmean","CWmean")

lmabestmods.ne <- rbind(PSEMENlma.ne$best, PINPONlma.ne$best,PINCONlma.ne$best,PINJEFlma.ne$best, ABICONlma.ne$best, TSUHETlma.ne$best, SPPlma.ne$best, CWlma.ne$best)
lmabestmods.ne$n <- rbind(PSEMENlma.ne$n, PINPONlma.ne$n,PINCONlma.ne$n,PINJEFlma.ne$n, ABICONlma.ne$n, TSUHETlma.ne$n, SPPlma.ne$n, CWlma.ne$n)
lmabestmods.ne$nplots <- rbind(PSEMENlma.ne$nplots, PINPONlma.ne$nplots,PINCONlma.ne$nplots,PINJEFlma.ne$nplots, ABICONlma.ne$nplots, TSUHETlma.ne$nplots, NA, CWlma.ne$nplots)
lmabestmods.ne$necos <- rbind(PSEMENlma.ne$necos, PINPONlma.ne$necos,PINCONlma.ne$necos,PINJEFlma.ne$necos, ABICONlma.ne$necos, TSUHETlma.ne$necos, NA, CWlma.ne$necos)
lmabestmods.ne$deltas <-rbind(PSEMENlma.ne$deltaNULL, PINPONlma.ne$deltaNULL,PINCONlma.ne$deltaNULL,PINJEFlma.ne$deltaNULL, ABICONlma.ne$deltaNULL, TSUHETlma.ne$deltaNULL, SPPlma.ne$deltaNULL, CWlma.ne$deltaNULL)
lmabestmods.ne$SP.ID <-  c("PSEMEN","PINPON","PINCON","PINJEF","ABICON","TSUHET","SPPmean", "CWmean")


nareabestmods.ne <- rbind(PSEMENnarea.ne$best, PINPONnarea.ne$best,PINCONnarea.ne$best,PINJEFnarea.ne$best, ABICONnarea.ne$best, TSUHETnarea.ne$best, SPPnarea.ne$best, CWnarea.ne$best)
nareabestmods.ne$n <- rbind(PSEMENnarea.ne$n, PINPONnarea.ne$n,PINCONnarea.ne$n,PINJEFnarea.ne$n, ABICONnarea.ne$n, TSUHETnarea.ne$n, SPPnarea.ne$n, CWnarea.ne$n)
nareabestmods.ne$nplots <- rbind(PSEMENnarea.ne$nplots, PINPONnarea.ne$nplots,PINCONnarea.ne$nplots,PINJEFnarea.ne$nplots, ABICONnarea.ne$nplots, TSUHETnarea.ne$nplots, NA, CWnarea.ne$nplots)
nareabestmods.ne$necos <- rbind(PSEMENnarea.ne$necos, PINPONnarea.ne$necos,PINCONnarea.ne$necos,PINJEFnarea.ne$necos, ABICONnarea.ne$necos, TSUHETnarea.ne$necos, NA, CWnarea.ne$necos)
nareabestmods.ne$deltas <-rbind(PSEMENnarea.ne$deltaNULL, PINPONnarea.ne$deltaNULL,PINCONnarea.ne$deltaNULL,PINJEFnarea.ne$deltaNULL, ABICONnarea.ne$deltaNULL, TSUHETnarea.ne$deltaNULL,  SPPnarea.ne$deltaNULL,CWnarea.ne$deltaNULL)
nareabestmods.ne$SP.ID <-  c("PSEMEN","PINPON","PINCON","PINJEF","ABICON","TSUHET", "SPPmean", "CWmean")
#nareabestmods.ne$Species <-  c("Pseudotsuga menziesii","Pinus ponderosa","Pinus contorta","Pinus jeffreyii","Abies concolor","Tsuga heterophylla")


nmassbestmods.ne <- rbind(PSEMENnmass.ne$best, PINPONnmass.ne$best,PINCONnmass.ne$best,PINJEFnmass.ne$best, ABICONnmass.ne$best, TSUHETnmass.ne$best,  SPPnmass.ne$best,CWnmass.ne$best)
nmassbestmods.ne$n <- rbind(PSEMENnmass.ne$n, PINPONnmass.ne$n,PINCONnmass.ne$n,PINJEFnmass.ne$n, ABICONnmass.ne$n, TSUHETnmass.ne$n,  SPPnmass.ne$n,CWnmass.ne$n)
nmassbestmods.ne$nplots <- rbind(PSEMENnmass.ne$nplots, PINPONnmass.ne$nplots,PINCONnmass.ne$nplots,PINJEFnmass.ne$nplots, ABICONnmass.ne$nplots, TSUHETnmass.ne$nplots, NA,CWnmass.ne$nplots)
nmassbestmods.ne$necos <- rbind(PSEMENnmass.ne$necos, PINPONnmass.ne$necos,PINCONnmass.ne$necos,PINJEFnmass.ne$necos, ABICONnmass.ne$necos, TSUHETnmass.ne$necos, NA,CWnmass.ne$necos)
nmassbestmods.ne$deltas <-rbind(PSEMENnmass.ne$deltaNULL, PINPONnmass.ne$deltaNULL,PINCONnmass.ne$deltaNULL,PINJEFnmass.ne$deltaNULL, ABICONnmass.ne$deltaNULL, TSUHETnmass.ne$deltaNULL,  SPPnmass.ne$deltaNULL,CWnmass.ne$deltaNULL)
nmassbestmods.ne$SP.ID <-  c("PSEMEN","PINPON","PINCON","PINJEF","ABICON","TSUHET", "SPPmean","CWmean")
#nmassbestmods.ne$Species <-  c("Pseudotsuga menziesii","Pinus ponderosa","Pinus contorta","Pinus jeffreyii","Abies concolor","Tsuga heterophylla")






####### Combining all environmental trait models NE ph ######
llbestmods.ne.ph <- rbind(PSEMENll.ne.ph$best, PINPONll.ne.ph$best, TSUHETll.ne.ph$best,SPPll.ne.ph$best, CWll.ne.ph$best)
llbestmods.ne.ph$n <- rbind(PSEMENll.ne.ph$n, PINPONll.ne.ph$n, ABICONll.ne.ph$n, TSUHETll.ne.ph$n,SPPll.ne.ph$n,CWll.ne.ph$n)
llbestmods.ne.ph$nplots <- rbind(PSEMENll.ne.ph$nplots, PINPONll.ne.ph$nplots, ABICONll.ne.ph$nplots, TSUHETll.ne.ph$nplots, NA, CWll.ne.ph$nplots)
llbestmods.ne.ph$necos <- rbind(PSEMENll.ne.ph$necos, PINPONll.ne.ph$necos, ABICONll.ne.ph$necos, TSUHETll.ne.ph$necos,NA, CWll.ne.ph$necos)
llbestmods.ne.ph$deltas <-rbind(PSEMENll.ne.ph$deltaNULL, PINPONll.ne.ph$deltaNULL, ABICONll.ne.ph$deltaNULL, TSUHETll.ne.ph$deltaNULL,SPPll.ne.ph$deltaNULL,CWll.ne.ph$deltaNULL)
llbestmods.ne.ph$SP.ID <-  c("PSEMEN","PINPON","ABICON","TSUHET", "SPPmean","CWmean")

lmabestmods.ne.ph <- rbind(PSEMENlma.ne.ph$best, PINPONlma.ne.ph$best, ABICONlma.ne.ph$best, TSUHETlma.ne.ph$best, SPPlma.ne.ph$best, CWlma.ne.ph$best)
lmabestmods.ne.ph$n <- rbind(PSEMENlma.ne.ph$n, PINPONlma.ne.ph$n, ABICONlma.ne.ph$n, TSUHETlma.ne.ph$n, SPPlma.ne.ph$n, CWlma.ne.ph$n)
lmabestmods.ne.ph$nplots <- rbind(PSEMENlma.ne.ph$nplots, PINPONlma.ne.ph$nplots, ABICONlma.ne.ph$nplots, TSUHETlma.ne.ph$nplots, NA, CWlma.ne.ph$nplots)
lmabestmods.ne.ph$necos <- rbind(PSEMENlma.ne.ph$necos, PINPONlma.ne.ph$necos, ABICONlma.ne.ph$necos, TSUHETlma.ne.ph$necos, NA, CWlma.ne.ph$necos)
lmabestmods.ne.ph$deltas <-rbind(PSEMENlma.ne.ph$deltaNULL, PINPONlma.ne.ph$deltaNULL, ABICONlma.ne.ph$deltaNULL, TSUHETlma.ne.ph$deltaNULL, SPPlma.ne.ph$deltaNULL, CWlma.ne.ph$deltaNULL)
lmabestmods.ne.ph$SP.ID <-  c("PSEMEN","PINPON","ABICON","TSUHET","SPPmean", "CWmean")


nareabestmods.ne.ph <- rbind(PSEMENnarea.ne.ph$best, PINPONnarea.ne.ph$best, ABICONnarea.ne.ph$best, TSUHETnarea.ne.ph$best, SPPnarea.ne.ph$best, CWnarea.ne.ph$best)
nareabestmods.ne.ph$n <- rbind(PSEMENnarea.ne.ph$n, PINPONnarea.ne.ph$n, ABICONnarea.ne.ph$n, TSUHETnarea.ne.ph$n, SPPnarea.ne.ph$n, CWnarea.ne.ph$n)
nareabestmods.ne.ph$nplots <- rbind(PSEMENnarea.ne.ph$nplots, PINPONnarea.ne.ph$nplots, ABICONnarea.ne.ph$nplots, TSUHETnarea.ne.ph$nplots, NA, CWnarea.ne.ph$nplots)
nareabestmods.ne.ph$necos <- rbind(PSEMENnarea.ne.ph$necos, PINPONnarea.ne.ph$necos,ABICONnarea.ne.ph$necos, TSUHETnarea.ne.ph$necos, NA, CWnarea.ne.ph$necos)
nareabestmods.ne.ph$deltas <-rbind(PSEMENnarea.ne.ph$deltaNULL, PINPONnarea.ne.ph$deltaNULL, ABICONnarea.ne.ph$deltaNULL, TSUHETnarea.ne.ph$deltaNULL,  SPPnarea.ne.ph$deltaNULL,CWnarea.ne.ph$deltaNULL)
nareabestmods.ne.ph$SP.ID <-  c("PSEMEN","PINPON","ABICON","TSUHET", "SPPmean", "CWmean")
#nareabestmods.ne.ph$Species <-  c("Pseudotsuga menziesii","Pinus ponderosa","Abies concolor","Tsuga heterophylla")


nmassbestmods.ne.ph <- rbind(PSEMENnmass.ne.ph$best, PINPONnmass.ne.ph$best, ABICONnmass.ne.ph$best, TSUHETnmass.ne.ph$best,  SPPnmass.ne.ph$best,CWnmass.ne.ph$best)
nmassbestmods.ne.ph$n <- rbind(PSEMENnmass.ne.ph$n, PINPONnmass.ne.ph$n, ABICONnmass.ne.ph$n, TSUHETnmass.ne.ph$n,  SPPnmass.ne.ph$n,CWnmass.ne.ph$n)
nmassbestmods.ne.ph$nplots <- rbind(PSEMENnmass.ne.ph$nplots, PINPONnmass.ne.ph$nplots, ABICONnmass.ne.ph$nplots, TSUHETnmass.ne.ph$nplots, NA,CWnmass.ne.ph$nplots)
nmassbestmods.ne.ph$necos <- rbind(PSEMENnmass.ne.ph$necos, PINPONnmass.ne.ph$necos, ABICONnmass.ne.ph$necos, TSUHETnmass.ne.ph$necos, NA,CWnmass.ne.ph$necos)
nmassbestmods.ne.ph$deltas <-rbind(PSEMENnmass.ne.ph$deltaNULL, PINPONnmass.ne.ph$deltaNULL, ABICONnmass.ne.ph$deltaNULL, TSUHETnmass.ne.ph$deltaNULL,  SPPnmass.ne.ph$deltaNULL,CWnmass.ne.ph$deltaNULL)
nmassbestmods.ne.ph$SP.ID <-  c("PSEMEN","PINPON","ABICON","TSUHET", "SPPmean","CWmean")
#nmassbestmods.ne.ph$Species <-  c("Pseudotsuga menziesii","Pinus ponderosa","Pinus contorta","Pinus jeffreyii","Abies concolor","Tsuga heterophylla")







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
# write.csv(llbestmods.ne, "Trait_models/LeafLife_bestmodels_noECOREG_20170619.csv") # shit. accidentally overwrote this with 0621 version
# write.csv(lmabestmods.ne, "Trait_models/LMA_bestmodels_noECOREG_20170619.csv")
# write.csv(nareabestmods.ne, "Trait_models/Narea_bestmodels_noECOREG_20170619.csv")
# write.csv(nmassbestmods.ne, "Trait_models/Nmass_bestmodels_noECOREG_20170619.csv")

# _v2 has CWmeans added on
# write.csv(llbestmods, "Trait_models/LeafLife_bestmodels_20170619_v2.csv")
# write.csv(lmabestmods, "Trait_models/LMA_bestmodels_20170621_v2.csv")
# write.csv(nareabestmods, "Trait_models/Narea_bestmodels_20170619_v2.csv")
# write.csv(nmassbestmods, "Trait_models/Nmass_bestmodels_20170619_v2.csv")

# New versions with _if CW means (and CW means at all for NE versions)
# write.csv(llbestmods, "Trait_models/LeafLife_bestmodels_20170621.csv")
# write.csv(lmabestmods, "Trait_models/LMA_bestmodels_20170621.csv")
# write.csv(nareabestmods, "Trait_models/Narea_bestmodels_20170621.csv")
# write.csv(nmassbestmods, "Trait_models/Nmass_bestmodels_20170621.csv")
#   # no ecoregions
# write.csv(llbestmods.ne, "Trait_models/LeafLife_bestmodels_noECOREG_20170621.csv")
# write.csv(lmabestmods.ne, "Trait_models/LMA_bestmodels_noECOREG_20170621.csv")
# write.csv(nareabestmods.ne, "Trait_models/Narea_bestmodels_noECOREG_20170621.csv")
# write.csv(nmassbestmods.ne, "Trait_models/Nmass_bestmodels_noECOREG_20170621.csv")


### now with SPPmeans
write.csv(llbestmods.ne, "Trait_models/LeafLife_bestmodels_noECOREG_20170726.csv")
write.csv(lmabestmods.ne, "Trait_models/LMA_bestmodels_noECOREG_20170726.csv")
write.csv(nareabestmods.ne, "Trait_models/Narea_bestmodels_noECOREG_20170726.csv")
write.csv(nmassbestmods.ne, "Trait_models/Nmass_bestmodels_noECOREG_20170726.csv")


# New versions 
write.csv(llbestmods, "Trait_models/LeafLife_bestmodels_201706")
write.csv(lmabestmods, "Trait_models/LMA_bestmodels_201706")
write.csv(nareabestmods, "Trait_models/Narea_bestmodels_201706")
write.csv(nmassbestmods, "Trait_models/Nmass_bestmodels_201706")
  # no ecoregions
# write.csv(llbestmods.ne, "Trait_models/LeafLife_bestmodels_noECOREG_20170726.csv")
# write.csv(lmabestmods.ne, "Trait_models/LMA_bestmodels_noECOREG_20170726.csv")
# write.csv(nareabestmods.ne, "Trait_models/Narea_bestmodels_noECOREG_20170726.csv")
# write.csv(nmassbestmods.ne, "Trait_models/Nmass_bestmodels_noECOREG_201707.csv")




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
# write.csv(importances, "Trait_models/VariableImportances_wide_20170619_v2aiccuttoff.csv") # same models but with CWmeans added
write.csv(importances, "Trait_models/VariableImportances_wide_20170621aiccuttoff.csv") # CW models with infilled values





#### Variable Importances NE ######
llimps.ne <- data.frame(rbind(PSEMENll.ne$imp, PINPONll.ne$imp,PINCONll.ne$imp,PINJEFll.ne$imp, ABICONll.ne$imp, TSUHETll.ne$imp, SPPll.ne$imp, CWll.ne$imp)[,-1])
colnames(llimps.ne) <- c("climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc")
llimps.ne$SP.ID <- c("PSEMEN","PINPON","PINCON","PINJEF","ABICON","TSUHET","SPPmean", "CWmean")
llimps.ne$trait <- rep("LeafLife", times=nrow(llimps.ne))

lmaimps.ne <- data.frame(rbind(PSEMENlma.ne$imp, PINPONlma.ne$imp,PINCONlma.ne$imp,PINJEFlma.ne$imp, ABICONlma.ne$imp, TSUHETlma.ne$imp, SPPlma.ne$imp, CWlma.ne$imp)[,-1])
colnames(lmaimps.ne) <- c("climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc")
lmaimps.ne$SP.ID <- c("PSEMEN","PINPON","PINCON","PINJEF","ABICON","TSUHET", "SPPmean","CWmean")
lmaimps.ne$trait <- rep("LMA", times=nrow(lmaimps.ne))

nareaimps.ne <- data.frame(rbind(PSEMENnarea.ne$imp, PINPONnarea.ne$imp,PINCONnarea.ne$imp,PINJEFnarea.ne$imp, ABICONnarea.ne$imp, TSUHETnarea.ne$imp, SPPnarea.ne$imp, CWnarea.ne$imp)[,-1])
colnames(nareaimps.ne) <- c("climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc")
nareaimps.ne$SP.ID <- c("PSEMEN","PINPON","PINCON","PINJEF","ABICON","TSUHET", "SPPmean","CWmean")
nareaimps.ne$trait <- rep("Narea", times=nrow(nareaimps.ne))


nmassimps.ne <- data.frame(rbind(PSEMENnmass.ne$imp, PINPONnmass.ne$imp,PINCONnmass.ne$imp,PINJEFnmass.ne$imp, ABICONnmass.ne$imp, TSUHETnmass.ne$imp, SPPnmass.ne$imp, CWnmass.ne$imp)[,-1])
colnames(nmassimps.ne) <- c("climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc")
nmassimps.ne$SP.ID <- c("PSEMEN","PINPON","PINCON","PINJEF","ABICON","TSUHET", "SPPmean","CWmean")
nmassimps.ne$trait <- rep("Nmass", times=nrow(nmassimps.ne))


importances.ne1 <- rbind(llimps.ne, lmaimps.ne, nareaimps.ne, nmassimps.ne)
colnames(importances.ne1) <- c("climPC1", "climPC2","soil_N","Stand_Age","LAI","Growth", "SP.ID", "trait")
#importances.ne$ECOREGION[which(is.na(importances.ne$ECOREGION))] <- 0 # change NAs to 0s (i.e. when no model contained ECOREGION, it should have an importance of 0)
# write.csv(importances, "Trait_models/VariableImportances_wide_031617aiccuttoff.csv")
# write.csv(importances, "Trait_models/VariableImportances_wide_041517aiccuttoff.csv")
# write.csv(importances, "Trait_models/VariableImportances_wide_wECOREGION_061617aiccuttoff.csv")
# write.csv(importances, "Trait_models/VariableImportances_wide_20170619aiccuttoff.csv")
# write.csv(importances, "Trait_models/VariableImportances_wide_20170619_v2aiccuttoff.csv") # same models but with CWmeans added
# write.csv(importances.ne, "Trait_models/VariableImportances_wide_noECOREG_20170621aiccuttoff.csv") # CW models with infilled values
# write.csv(importances.ne, "Trait_models/VariableImportances_wide_noECOREG_20170726.csv") # CW models with infilled values
write.csv(importances.ne, "Trait_models/VariableImportances_wide_noECOREG_") # CW models with infilled values



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
#write.csv(avgmods, "Trait_models/Model_Averages_wide_20170619_v2aiccutoff.csv")
#write.csv(avgmods, "Trait_models/Model_Averages_wide_20170621aiccutoff.csv") # with _if CWmean values





##### Model Averaged Effects Sizes NE ######
llavgmods.ne <- data.frame(rbind(PSEMENll.ne$avg, PINPONll.ne$avg,PINCONll.ne$avg,PINJEFll.ne$avg, ABICONll.ne$avg, TSUHETll.ne$avg,SPPll.ne$avg, CWll.ne$avg))
llavgmods.ne$SP.ID <- c("PSEMEN","PINPON","PINCON","PINJEF","ABICON","TSUHET","SPPmean", "CWmean")
llavgmods.ne$trait <- rep("LeafLife", times=nrow(llavgmods.ne))

lmaavgmods.ne <- data.frame(rbind(PSEMENlma.ne$avg, PINPONlma.ne$avg,PINCONlma.ne$avg,PINJEFlma.ne$avg, ABICONlma.ne$avg, TSUHETlma.ne$avg,SPPlma.ne$avg, CWlma.ne$avg))
lmaavgmods.ne$SP.ID <- c("PSEMEN","PINPON","PINCON","PINJEF","ABICON","TSUHET","SPPmean","CWmean")
lmaavgmods.ne$trait <- rep("LMA", times=nrow(lmaavgmods.ne))

nareaavgmods.ne <- data.frame(rbind(PSEMENnarea.ne$avg, PINPONnarea.ne$avg,PINCONnarea.ne$avg,PINJEFnarea.ne$avg, ABICONnarea.ne$avg, TSUHETnarea.ne$avg, SPPnarea.ne$avg, CWnarea.ne$avg))
nareaavgmods.ne$SP.ID <- c("PSEMEN","PINPON","PINCON","PINJEF","ABICON","TSUHET","SPPmean","CWmean")
nareaavgmods.ne$trait <- rep("Narea", times=nrow(nareaavgmods.ne))

nmassavgmods.ne <- data.frame(rbind(PSEMENnmass.ne$avg, PINPONnmass.ne$avg,PINCONnmass.ne$avg,PINJEFnmass.ne$avg, ABICONnmass.ne$avg, TSUHETnmass.ne$avg, SPPnmass.ne$avg, CWnmass.ne$avg))
nmassavgmods.ne$SP.ID <- c("PSEMEN","PINPON","PINCON","PINJEF","ABICON","TSUHET","SPPmean","CWmean")
nmassavgmods.ne$trait <- rep("Nmass", times=nrow(nmassavgmods.ne))

avgmods.ne <- rbind(llavgmods.ne, lmaavgmods.ne, nareaavgmods.ne, nmassavgmods.ne)[,-1]
colnames(avgmods.ne) <- c("climPC1", "climPC2","soil_N","Stand_Age","LAI","Growth", "SP.ID", "trait")
#write.csv(avgmods, "Trait_models/Model_Averages_wide_031617aiccuttoff.csv")
#write.csv(avgmods, "Trait_models/Model_Averages_wide_041517aiccuttoff.csv")
#write.csv(avgmods, "Trait_models/Model_Averages_wide_061617aiccutoff.csv")
#write.csv(avgmods, "Trait_models/Model_Averages_wide_20170619aiccutoff.csv")
#write.csv(avgmods, "Trait_models/Model_Averages_wide_20170619_v2aiccutoff.csv")
#write.csv(avgmods.ne, "Trait_models/Model_Averages_wide_noECOREG_20170621aiccutoff.csv") # with _if CWmean values
write.csv(avgmods.ne, "Trait_models/Model_Averages_wide_noECOREG_20170726aiccutoff.csv") # with _if CWmean values












#_______________________________________________________
########## Modeling with Ecoregion Alone ################
#_______________________________________________________

spp <- "PSEMEN"
trait <- "log.LL"
traitmod <- lmer(get(trait)~ECOREGION + (1|PLOT_ID), traits.common[which(traits.common$SP.ID==spp & !is.na(traits.common[,trait])),], REML=T)
r.squaredGLMM(traitmod)


ecoregmods <- data.frame(SP.ID = c("PSEMEN","PINPON","PINCON","PINJEF","ABICON","TSUHET")
                        , LLr2 = rep(NA, times=6), LMAr2 = rep(NA, times=6), Nmassr2= rep(NA, times=6), Narear2= rep(NA, times=6), nEcos= rep(NA, times=6))

for(i in 1:6){ # loop over species
  spp <- c("PSEMEN","PINPON","PINCON","PINJEF","ABICON","TSUHET")[i]
  for(j in 1:4) {# loop over traits
    trait <- c("log.LL", "log.LMA","log.Nmass","log.Narea")[j]
    tmpmod <- lmer(get(trait)~ECOREGION + (1|PLOT_ID), traits.common[which(traits.common$SP.ID==spp & !is.na(traits.common[,trait])),], REML=T)
    ecoregmods[i, j+1] <- r.squaredGLMM(tmpmod)[1]
  }
  ecoregmods[i, 6] <- length(unique(traits.common$ECOREGION[which(traits.common$SP.ID==spp & !is.na(traits.common[,trait]))]))
}


CWecoregmods <- data.frame(LLr2 = summary(lm(log.cw_LLp_if~ECOREGION, biomass[which(!is.na(biomass$log.cw_LLp_if)),]))$r.squared
                           , LMAr2 = summary(lm(log.cw_LMAp_if~ECOREGION, biomass[which(!is.na(biomass$log.cw_LMAp_if)),]))$r.squared
                           , Nmassr2 = summary(lm(log.cw_Nmassp_if~ECOREGION, biomass[which(!is.na(biomass$log.cw_Nmassp_if)),]))$r.squared
                           , Narear2 = summary(lm(log.cw_Nareap_if~ECOREGION, biomass[which(!is.na(biomass$log.cw_Nareap_if)),]))$r.squared)
CWecoregmods <- unlist(CWecoregmods)



#_______________________________________________________________________________________________
########### **Fig5: Plotting Model Ensamble Output, w/ECOREG ########################
###
#_______________________________________________________________________________________________
#importances.we <- read.csv("Trait_models/VariableImportances_wide_wECOREGION_061617aiccuttoff.csv", row.names = 1)
# importances <- read.csv("Trait_models/VariableImportances_wide_20170612_v2aiccuttoff.csv", row.names = 1)
importances <- read.csv("Trait_models/VariableImportances_wide_20170621aiccuttoff.csv", row.names = 1)
importances.ne <-read.csv("Trait_models/VariableImportances_wide_noECOREG_20170726.csv", row.names = 1)
#importances.we$ECOREGION[which(is.na(importances.we$ECOREGION))] <- 0
impslong <- melt(importances, id.vars = c("SP.ID","trait"))
impslong.ne <- melt(importances.ne, id.vars=c("SP.ID","trait"))
#impslong.we <- melt(importances.we, id.vars= c("SP.ID", "trait"))

#avgmods <- read.csv("Trait_models/Model_Averages_wide_041517aiccuttoff.csv", row.names = 1)
#avgmods <- read.csv("Trait_models/Model_Averages_wide_20170619_v2aiccutoff.csv", row.names = 1)
avgmods <- read.csv("Trait_models/Model_Averages_wide_20170621aiccutoff.csv", row.names = 1)
avgmods.ne <- read.csv("Trait_models/Model_Averages_wide_noECOREG_20170726aiccutoff.csv", row.names = 1)
avglong <- melt(avgmods, id.vars = c("SP.ID","trait"))
avglong.ne <- melt(avgmods.ne, id.vars = c("SP.ID","trait"))

avglong$xvals <- as.numeric(avglong$variable)
impslong$xvals <- as.numeric(impslong$variable)
avglong.ne$xvals <- as.numeric(avglong.ne$variable)
impslong.ne$xvals <- as.numeric(impslong.ne$variable)



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




### New boxplots, with ECOREG

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
boxplot(nareabestmods$r.squaredGLMM.R2m[-which(nareabestmods$SP.ID=="CWmean")]~rep(1, times=6), horizontal=T, at=1, xlim=c(0,5), ylim=c(0,1)
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





#_______________________________________________________________________________
################## **R2 comparison Ecoreg vs all Env Vars** ######################
#_______________________________________________________________________________

panlab.cex <- .9
panlab.ln <- -1.2
CWbg <- "white"#mypal[6]#"white"
SPPbg <- "black"#mypal[colchoices[2]]#"black"#mypal[2]
boxcol <- "grey"#paste0(mypal[colchoices[1]],"66")#"grey"
ptcex=1.1

quartz(width=4.33, height=2.5)
par(mar=c(.2,5,3,3), mgp=c(2,.75,0))
boxplot(llbestmods.ne$r.squaredGLMM.R2m[-which(llbestmods.ne$SP.ID %in% c("SPPmean","CWmean"))]~rep(1, times=6), horizontal=T, at=4, xlim=c(0,5), ylim=c(0,1)
        , yaxt="n", xaxt="n"
        ,outpch=1, border=boxcol, col=boxcol, outcol=boxcol
        , staplewex=0, notch=F, medcol="black" #, medlwd=0
        , whisklty=1, whisklwd=2, boxwex=.5, whiskcol=boxcol
        , range=0)
boxplot(lmabestmods.ne$r.squaredGLMM.R2m[-which(lmabestmods.ne$SP.ID%in% c("SPPmean","CWmean"))]~rep(1, times=6), horizontal=T, at=3, xlim=c(0,5), ylim=c(0,1)
        , yaxt="n", xaxt="n"
        ,outpch=1, border=boxcol, col=boxcol, outcol=boxcol
        , staplewex=0, notch=F, medcol="black" #, medlwd=0
        , whisklty=1, whisklwd=2, boxwex=.5, whiskcol=boxcol, add=T , range=0)
boxplot(nmassbestmods.ne$r.squaredGLMM.R2m[-which(nmassbestmods.ne$SP.ID%in% c("SPPmean","CWmean"))]~rep(1, times=6), horizontal=T, at=2, xlim=c(0,5), ylim=c(0,1)
        , yaxt="n", xaxt="n"
        ,outpch=1, border=boxcol, col=boxcol, outcol=boxcol
        , staplewex=0, notch=F, medcol="black" #, medlwd=0
        , whisklty=1, whisklwd=2, boxwex=.5, whiskcol=boxcol, add=T , range=0)
boxplot(nareabestmods.ne$r.squaredGLMM.R2m[-which(nareabestmods.ne$SP.ID%in% c("SPPmean","CWmean"))]~rep(1, times=6), horizontal=T, at=1, xlim=c(0,5), ylim=c(0,1)
        , yaxt="n", xaxt="n"
        ,outpch=1, border=boxcol, col=boxcol, outcol=boxcol
        , staplewex=0, notch=F, medcol="black" #, medlwd=0
        , whisklty=1, whisklwd=2, boxwex=.5, whiskcol=boxcol, add=T , range=0)


#plot(llbestmods.ne$r.squaredGLMM.R2m[-which(llbestmods.ne$SP.ID=="CWmean")]~rep(1, times=6), ylim=c(0,1), xlim=c(0,5), xaxt="n", ylab=expression(paste(R^2," or marginal ",R^2)), pch=16)
#points(lmabestmods.ne$r.squaredGLMM.R2m[-which(lmabestmods.ne$SP.ID=="CWmean")]~rep(2, times=6), pch=16)
#points(nmassbestmods.ne$r.squaredGLMM.R2m[-which(nmassbestmods.ne$SP.ID=="CWmean")]~rep(3, times=6), pch=16)
#points(nareabestmods.ne$r.squaredGLMM.R2m[-which(nareabestmods.ne$SP.ID=="CWmean")]~rep(4, times=6), pch=16)
points(c(4)~llbestmods.ne$r.squaredGLMM.R2c[which(llbestmods.ne$SP.ID=="CWmean")], pch=24,bg=CWbg) # the simple R^2 is stored in R2c for CWmean
points(c(3)~lmabestmods.ne$r.squaredGLMM.R2c[which(lmabestmods.ne$SP.ID=="CWmean")], pch=24,bg=CWbg)
points(c(2)~nmassbestmods.ne$r.squaredGLMM.R2c[which(nmassbestmods.ne$SP.ID=="CWmean")], pch=24,bg=CWbg)
points(c(1)~nareabestmods.ne$r.squaredGLMM.R2c[which(nareabestmods.ne$SP.ID=="CWmean")], pch=24,bg=CWbg)
#SPP means
#points(c(4)~llbestmods.ne$r.squaredGLMM.R2c[which(llbestmods.ne$SP.ID=="SPPmean")], pch=25, bg=SPPbg) # the simple R^2 is stored in R2c for CWmean
#points(c(3)~lmabestmods.ne$r.squaredGLMM.R2c[which(lmabestmods.ne$SP.ID=="SPPmean")], pch=25, bg=SPPbg)
#points(c(2)~nmassbestmods.ne$r.squaredGLMM.R2c[which(nmassbestmods.ne$SP.ID=="SPPmean")], pch=25, bg=SPPbg)
#points(c(1)~nareabestmods.ne$r.squaredGLMM.R2c[which(nareabestmods.ne$SP.ID=="SPPmean")], pch=25, bg=SPPbg)




axis(3)
axis(2, labels = c("Narea","Nmass","LMA","Leaf Life"), at=c(1,2,3,4), las=2)
mtext(side=3,text = expression(paste(R^2," or marginal ", R^2)), line=1.7)
#legend(x=.675, y=1, legend = "Ind. Spp", fill=boxcol, bty="n",col = boxcol,border = boxcol, cex=.8)
#legend(x=.7, y=2.3, legend = c("CW mean","Spp mean"), pt.bg=c(CWbg, SPPbg), bty="n",pch=c(24,25), cex=.8)
# legend(x=.7-.1, y=3+2.5, legend = c("CWM Ecoreg", "Spp Ecoreg"), bty="n", pch=c(23,4), cex=.8)

 legend(x=.73, y=2.5, legend = c("CWM-eco","CWM-full","Spp-eco"), pt.bg=c(CWbg, mypal[1], mypal[1]), bty="n",pch=c(24,24,22), cex=.8)
 legend(x=.71, y=1.1, legend = "Spp-full", fill="gray", bty="n",col = boxcol,border = boxcol, cex=.8)
mtext("A)", cex=panlab.cex, side=3, adj=0.02, line=panlab.ln)

#### Add in ecoregion only models
boxplot(ecoregmods$LLr2~rep(1.2, times=6), horizontal=T, at=3.7, xlim=c(0,5), ylim=c(0,1)
        , yaxt="n", xaxt="n", add=T
        ,outpch=1, border=mypal[1], col=mypal[1], outcol=mypal[1]
        , staplewex=0, notch=F, medcol="black" #, medlwd=0
        , whisklty=1, whisklwd=2, boxwex=.5, whiskcol=mypal[1]
        , range=0)

boxplot(ecoregmods$LMAr2~rep(1.2, times=6), horizontal=T, at=2.7, xlim=c(0,5), ylim=c(0,1)
        , yaxt="n", xaxt="n", add=T
        ,outpch=1, border=mypal[1], col=mypal[1], outcol=mypal[1]
        , staplewex=0, notch=F, medcol="black" #, medlwd=0
        , whisklty=1, whisklwd=2, boxwex=.5, whiskcol=mypal[1]
        , range=0)

boxplot(ecoregmods$Nmassr2~rep(1.2, times=6), horizontal=T, at=1.7, xlim=c(0,5), ylim=c(0,1)
        , yaxt="n", xaxt="n", add=T
        ,outpch=1, border=mypal[1], col=mypal[1], outcol=mypal[1]
        , staplewex=0, notch=F, medcol="black" #, medlwd=0
        , whisklty=1, whisklwd=2, boxwex=.5, whiskcol=mypal[1]
        , range=0)

boxplot(ecoregmods$Narear2~rep(1.2, times=6), horizontal=T, at=0.7, xlim=c(0,5), ylim=c(0,1)
        , yaxt="n", xaxt="n", add=T
        ,outpch=1, border=mypal[1], col=mypal[1], outcol=mypal[1]
        , staplewex=0, notch=F, medcol="black" #, medlwd=0
        , whisklty=1, whisklwd=2, boxwex=.5, whiskcol=mypal[1]
        , range=0)

# add in ecoregion only rsquareds
#points(c(4,3,2,1)-0.3~colMeans(ecoregmods[,-1]), pch=4,bg=mypal[1]) # the simple R^2 is stored in R2c for CWmean
points(c(4,3,2,1)-0.3~CWecoregmods, pch=24,bg=mypal[1]) # the simple R^2 is stored in R2c for CWmean





#__________________________________________________________________________________
########## **Fig5: Variable Imps, Ensemble Effecs, R2 NO ECOREG (.ne) #########
#_________________________________________________________________________________


quartz(width=4.33, height=6) # full page = 6.81
par(mfcol=c(4,2), mar=c(0,3.5,0,0), oma=c(5,0,2,1), mgp=c(2.5,.75,0), cex=.9)
panlab.cex <- .9
panlab.ln <- -1.2
CWbg <- "white"#mypal[6]#"white"
SPPbg <- "black"#mypal[colchoices[2]]#"black"#mypal[2]
boxcol <- "grey"#paste0(mypal[colchoices[1]],"66")#"grey"
ptcex=1.1
### Effect Sizes #####
boxplot(value~xvals, avglong.ne[which(avglong.ne$trait=="LeafLife" & avglong.ne$SP.ID!="CWmean" & avglong.ne$SP.ID!="SPPmean" ),], at=c(1,2,3,4,5,6)
        , ylim=c(-.7,1), xaxt="n", xlim=c(0.5,6.5)
        ,outpch=1, border=boxcol, col=boxcol, outcol=boxcol
        , staplewex=0, notch=F, medcol="black" #, medlwd=0
        , whisklty=1, whisklwd=2, boxwex=.7, whiskcol=boxcol #, boxwex=1, medlwd=3
        ,las=2, range=0)
#plot(value~xvals, avglong.ne[which(avglong.ne$trait=="LeafLife" & avglong.ne$SP.ID!="CWmean"),], pch=16, xaxt="n", ylab="", ylim=c(-0.7,1), xlab="", xlim=c(0.6,6.4))
abline(h=0)
points(value~xvals, avglong.ne[which(avglong.ne$trait=="LeafLife" & avglong.ne$SP.ID=="SPPmean"),], pch=25,cex=ptcex, bg=SPPbg)
points(value~xvals, avglong.ne[which(avglong.ne$trait=="LeafLife" & avglong.ne$SP.ID=="CWmean"),], pch=24, bg=CWbg, cex=ptcex)
mtext(text = "B)",side = 3,line = panlab.ln,adj=0.02,cex=panlab.cex )
mtext(text = "log(LL)",side = 3,line = panlab.ln,adj=0.9,cex=panlab.cex)
#text(x=7, y=0, pos=3, labels = "NA",cex=.7)

boxplot(value~xvals, avglong.ne[which(avglong.ne$trait=="LMA" & avglong.ne$SP.ID!="CWmean"),], at=c(1,2,3,4,5,6)
        , ylim=c(-.7,1), xaxt="n", xlim=c(0.5,6.5)
        ,outpch=1, border=boxcol, col=boxcol, outcol=boxcol
        , staplewex=0, notch=F, medcol="black" #, medlwd=0
        , whisklty=1, whisklwd=2, boxwex=.7, whiskcol=boxcol #, boxwex=1, medlwd=3
        ,las=2, range=0)
#plot(value~xvals, avglong.ne[which(avglong.ne$trait=="LMA" & avglong.ne$SP.ID!="CWmean"),], pch=16, xaxt="n", ylab="", ylim=c(-0.7,1), xlab="", xlim=c(0.6,6.4))
abline(h=0)
points(value~xvals, avglong.ne[which(avglong.ne$trait=="LMA" & avglong.ne$SP.ID=="SPPmean"),], pch=25, cex=ptcex,bg=SPPbg)
points(value~xvals, avglong.ne[which(avglong.ne$trait=="LMA" & avglong.ne$SP.ID=="CWmean"),], pch=24, bg=CWbg, cex=ptcex)
mtext(text = "D)",side = 3,line = panlab.ln,adj=0.02,cex=panlab.cex )
mtext(text = "log(LMA)",side = 3,line = panlab.ln,adj=0.9,cex=panlab.cex)
#text(x=7, y=0, pos=3, labels = "NA",cex=.7)

boxplot(value~xvals, avglong.ne[which(avglong.ne$trait=="Nmass" & avglong.ne$SP.ID!="CWmean"),], at=c(1,2,3,4,5,6)
        , ylim=c(-.7,1), xaxt="n", xlim=c(0.5,6.5)
        ,outpch=1, border=boxcol, col=boxcol, outcol=boxcol
        , staplewex=0, notch=F, medcol="black" #, medlwd=0
        , whisklty=1, whisklwd=2, boxwex=.7, whiskcol=boxcol #, boxwex=1, medlwd=3
        ,las=2, range=0)
# plot(value~xvals, avglong.ne[which(avglong.ne$trait=="Nmass" & avglong.ne$SP.ID!="CWmean"),], pch=16, xaxt="n", ylab="", ylim=c(-0.7,1), xlab="", xlim=c(0.6,6.4))
abline(h=0)
points(value~xvals, avglong.ne[which(avglong.ne$trait=="Nmass" & avglong.ne$SP.ID=="SPPmean"),],cex=ptcex, pch=25, bg=SPPbg)
points(value~xvals, avglong.ne[which(avglong.ne$trait=="Nmass" & avglong.ne$SP.ID=="CWmean"),], pch=24, bg=CWbg, cex=ptcex)
mtext(text = "F)",side = 3,line = panlab.ln,adj=0.02,cex=panlab.cex )
mtext(text = "log(Nmass)",side = 3,line = panlab.ln,adj=0.9,cex=panlab.cex)
mtext(text="Standardized Effect Sizes",side = 2, adj=-.1, line=2.1)
#text(x=7, y=0, pos=3, labels = "NA",cex=.7)

boxplot(value~xvals, avglong.ne[which(avglong.ne$trait=="Narea" & avglong.ne$SP.ID!="CWmean"),], at=c(1,2,3,4,5,6)
        , ylim=c(-.7,1), xaxt="n", xlim=c(0.5,6.5)
        ,outpch=1, border=boxcol, col=boxcol, outcol=boxcol
        , staplewex=0, notch=F, medcol="black" #, medlwd=0
        , whisklty=1, whisklwd=2, boxwex=.7, whiskcol=boxcol #, boxwex=1, medlwd=3
        ,las=2, range=0)
#plot(value~xvals, avglong.ne[which(avglong.ne$trait=="Narea" & avglong.ne$SP.ID!="CWmean"),], pch=16, xaxt="n", ylab="", ylim=c(-0.7,1), xlab="", xlim=c(0.6,6.4))
abline(h=0)
points(value~xvals, avglong.ne[which(avglong.ne$trait=="Narea" & avglong.ne$SP.ID=="SPPmean"),], pch=25, bg=SPPbg, cex=ptcex)
points(value~xvals, avglong.ne[which(avglong.ne$trait=="Narea" & avglong.ne$SP.ID=="CWmean"),], pch=24, bg=CWbg, cex=ptcex)
mtext(text = "H)",side = 3,line = panlab.ln,adj=0.02,cex=panlab.cex )
mtext(text = "log(Narea)",side = 3,line = panlab.ln,adj=0.9,cex=panlab.cex)
#text(x=7, y=0, pos=3, labels = "NA",cex=.7)

axis(side = 1,at = c(1,2,3,4,5,6), labels = c("Wetness", "Warmth","Soil N", "log(Age)","LAI","Gr Rate"), xpd=NA, las=3, srt=50, xlim=c(0.6,6.4))


### Variable Importances #####

boxplot(value~xvals, impslong.ne[which(impslong.ne$trait=="LeafLife" & impslong.ne$SP.ID!="CWmean"),], at=c(1,2,3,4,5,6)
        , ylim=c(0,1.3), yaxt="n", xaxt="n"
        ,outpch=1, border=boxcol, col=boxcol, outcol=boxcol
        , staplewex=0, notch=F, medcol="black" #, medlwd=0
        , whisklty=1, whisklwd=2, boxwex=.7, whiskcol=boxcol #, boxwex=1, medlwd=3
        ,las=2, range=0)
points(value~xvals, impslong.ne[which(impslong.ne$trait=="LeafLife" & impslong.ne$SP.ID=="SPPmean"),], pch=25,cex=ptcex, bg=SPPbg)
points(value~xvals, impslong.ne[which(impslong.ne$trait=="LeafLife" & impslong.ne$SP.ID=="CWmean"),], pch=24, bg=CWbg, cex=ptcex)
axis(2,at=c(0,.25,.5,.75,1), labels = c(0,.25,.5,.75,1))
mtext(text = "C)",side = 3,line = panlab.ln,adj=0.02,cex=panlab.cex )

boxplot(value~xvals, impslong.ne[which(impslong.ne$trait=="LMA" & impslong.ne$SP.ID!="CWmean"),], at=c(1,2,3,4,5,6)
        , ylim=c(0,1.3), yaxt="n", xaxt="n"
        ,outpch=1, border=boxcol, col=boxcol, outcol=boxcol
        , staplewex=0, notch=F, medcol="black" #, medlwd=0
        , whisklty=1, whisklwd=2, boxwex=.7, whiskcol=boxcol #, boxwex=1, medlwd=3
        ,las=2, range=0)
points(value~xvals, impslong.ne[which(impslong.ne$trait=="LMA" & impslong.ne$SP.ID=="SPPmean"),], pch=25,cex=ptcex, bg=SPPbg)
points(value~xvals, impslong.ne[which(impslong.ne$trait=="LMA" & impslong.ne$SP.ID=="CWmean"),], pch=24, bg=CWbg, cex=ptcex)
axis(2,at=c(0,.25,.5,.75,1), labels = c(0,.25,.5,.75,1))
mtext(text = "E)",side = 3,line = panlab.ln,adj=0.02,cex=panlab.cex )

boxplot(value~xvals, impslong.ne[which(impslong.ne$trait=="Nmass" & impslong.ne$SP.ID!="CWmean"),], at=c(1,2,3,4,5,6)
        , ylim=c(0,1.3), yaxt="n", xaxt="n"
        ,outpch=1, border=boxcol, col=boxcol, outcol=boxcol
        , staplewex=0, notch=F, medcol="black" #, medlwd=0
        , whisklty=1, whisklwd=2, boxwex=.7, whiskcol=boxcol #, boxwex=1, medlwd=3
        ,las=2, range=0)
points(value~xvals, impslong.ne[which(impslong.ne$trait=="Nmass" & impslong.ne$SP.ID=="SPPmean"),], pch=25,cex=ptcex, bg=SPPbg)
points(value~xvals, impslong.ne[which(impslong.ne$trait=="Nmass" & impslong.ne$SP.ID=="CWmean"),], pch=24, cex=ptcex,bg=CWbg)
axis(2,at=c(0,.25,.5,.75,1), labels = c(0,.25,.5,.75,1))
mtext(text = "G)",side = 3,line = panlab.ln,adj=0.02,cex=panlab.cex )
mtext(text="Vairable Importances",side = 2, adj=-.2, line=2.1)

boxplot(value~xvals, impslong.ne[which(impslong.ne$trait=="Narea" & impslong.ne$SP.ID!="CWmean"),], at=c(1,2,3,4,5,6)
        , ylim=c(0,1.3), yaxt="n", xaxt="n"
        ,outpch=1, border=boxcol, col=boxcol, outcol=boxcol
        , staplewex=0, notch=F, medcol="black" #, medlwd=0
        , whisklty=1, whisklwd=2, boxwex=.7, whiskcol=boxcol #, boxwex=1, medlwd=3
        ,las=2, range=0)
points(value~xvals, impslong.ne[which(impslong.ne$trait=="Narea" & impslong.ne$SP.ID=="SPPmean"),], pch=25,cex=ptcex, bg=SPPbg)
points(value~xvals, impslong.ne[which(impslong.ne$trait=="Narea" & impslong.ne$SP.ID=="CWmean"),], pch=24, cex=ptcex,bg=CWbg)
axis(2,at=c(0,.25,.5,.75,1), labels = c(0,.25,.5,.75,1))
mtext(text = "I)",side = 3,line = panlab.ln,adj=0.02,cex=panlab.cex )
axis(side = 1,at = c(1,2,3,4,5,6), labels = c("Wetness", "Warmth","Soil N", "log(Age)","LAI","Gr Rate"), xpd=NA, las=3, srt=50, xlim=c(0.6,6.4))



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
boxplot(llbestmods.ne$r.squaredGLMM.R2m[-which(llbestmods.ne$SP.ID %in% c("SPPmean","CWmean"))]~rep(1, times=6), horizontal=T, at=4, xlim=c(0,5), ylim=c(0,1)
        , yaxt="n", xaxt="n"
        ,outpch=1, border=boxcol, col=boxcol, outcol=boxcol
        , staplewex=0, notch=F, medcol="black" #, medlwd=0
        , whisklty=1, whisklwd=2, boxwex=.7, whiskcol=boxcol
        , range=0)
boxplot(lmabestmods.ne$r.squaredGLMM.R2m[-which(lmabestmods.ne$SP.ID%in% c("SPPmean","CWmean"))]~rep(1, times=6), horizontal=T, at=3, xlim=c(0,5), ylim=c(0,1)
        , yaxt="n", xaxt="n"
        ,outpch=1, border=boxcol, col=boxcol, outcol=boxcol
        , staplewex=0, notch=F, medcol="black" #, medlwd=0
        , whisklty=1, whisklwd=2, boxwex=.7, whiskcol=boxcol, add=T , range=0)
boxplot(nmassbestmods.ne$r.squaredGLMM.R2m[-which(nmassbestmods.ne$SP.ID%in% c("SPPmean","CWmean"))]~rep(1, times=6), horizontal=T, at=2, xlim=c(0,5), ylim=c(0,1)
        , yaxt="n", xaxt="n"
        ,outpch=1, border=boxcol, col=boxcol, outcol=boxcol
        , staplewex=0, notch=F, medcol="black" #, medlwd=0
        , whisklty=1, whisklwd=2, boxwex=.7, whiskcol=boxcol, add=T , range=0)
boxplot(nareabestmods.ne$r.squaredGLMM.R2m[-which(nareabestmods.ne$SP.ID%in% c("SPPmean","CWmean"))]~rep(1, times=6), horizontal=T, at=1, xlim=c(0,5), ylim=c(0,1)
        , yaxt="n", xaxt="n"
        ,outpch=1, border=boxcol, col=boxcol, outcol=boxcol
        , staplewex=0, notch=F, medcol="black" #, medlwd=0
        , whisklty=1, whisklwd=2, boxwex=.7, whiskcol=boxcol, add=T , range=0)


#plot(llbestmods.ne$r.squaredGLMM.R2m[-which(llbestmods.ne$SP.ID=="CWmean")]~rep(1, times=6), ylim=c(0,1), xlim=c(0,5), xaxt="n", ylab=expression(paste(R^2," or marginal ",R^2)), pch=16)
#points(lmabestmods.ne$r.squaredGLMM.R2m[-which(lmabestmods.ne$SP.ID=="CWmean")]~rep(2, times=6), pch=16)
#points(nmassbestmods.ne$r.squaredGLMM.R2m[-which(nmassbestmods.ne$SP.ID=="CWmean")]~rep(3, times=6), pch=16)
#points(nareabestmods.ne$r.squaredGLMM.R2m[-which(nareabestmods.ne$SP.ID=="CWmean")]~rep(4, times=6), pch=16)
points(c(4)~llbestmods.ne$r.squaredGLMM.R2c[which(llbestmods.ne$SP.ID=="CWmean")], pch=24,bg=CWbg) # the simple R^2 is stored in R2c for CWmean
points(c(3)~lmabestmods.ne$r.squaredGLMM.R2c[which(lmabestmods.ne$SP.ID=="CWmean")], pch=24,bg=CWbg)
points(c(2)~nmassbestmods.ne$r.squaredGLMM.R2c[which(nmassbestmods.ne$SP.ID=="CWmean")], pch=24,bg=CWbg)
points(c(1)~nareabestmods.ne$r.squaredGLMM.R2c[which(nareabestmods.ne$SP.ID=="CWmean")], pch=24,bg=CWbg)
  #SPP means
points(c(4)~llbestmods.ne$r.squaredGLMM.R2c[which(llbestmods.ne$SP.ID=="SPPmean")], pch=25, bg=SPPbg) # the simple R^2 is stored in R2c for CWmean
points(c(3)~lmabestmods.ne$r.squaredGLMM.R2c[which(lmabestmods.ne$SP.ID=="SPPmean")], pch=25, bg=SPPbg)
points(c(2)~nmassbestmods.ne$r.squaredGLMM.R2c[which(nmassbestmods.ne$SP.ID=="SPPmean")], pch=25, bg=SPPbg)
points(c(1)~nareabestmods.ne$r.squaredGLMM.R2c[which(nareabestmods.ne$SP.ID=="SPPmean")], pch=25, bg=SPPbg)

# add in ecoregion only rsquareds
# points(c(4,3,2,1)~colMeans(ecoregmods[,-1]), pch=4,bg=CWbg) # the simple R^2 is stored in R2c for CWmean
# points(c(4,3,2,1)~CWecoregmods, pch=23,bg=CWbg) # the simple R^2 is stored in R2c for CWmean



axis(3)
axis(2, labels = c("Narea","Nmass","LMA","Leaf Life"), at=c(1,2,3,4), las=2)
mtext(side=3,text = expression(paste(R^2," or marginal ", R^2)), line=1.7)
legend(x=.675, y=1, legend = "Ind. Spp", fill=boxcol, bty="n",col = boxcol,border = boxcol, cex=.8)
legend(x=.7, y=2.3, legend = c("CW mean","Spp mean"), pt.bg=c(CWbg, SPPbg), bty="n",pch=c(24,25), cex=.8)
# legend(x=.7-.1, y=3+2.5, legend = c("CWM Ecoreg", "Spp Ecoreg"), bty="n", pch=c(23,4), cex=.8)

# legend(x=.75, y=2.5, legend = c("CWM-eco","CWM-full","Spp-eco"), bg=CWbg, bty="n",pch=c(23,24,4), cex=.7)
# legend(x=.73, y=.9, legend = "Spp-env", fill="gray", bty="n",col = boxcol,border = boxcol, cex=.7)
mtext("A)", cex=panlab.cex, side=3, adj=0.02, line=panlab.ln)























#_________________________________________________________________
########## Effect Sizes for Poster Presentation #########
#_________________________________________________________________


quartz(width=4.33, height=4) # full page = 6.81
par(mfcol=c(2,2), mar=c(0,0,0,0), oma=c(5,3.5,2,1), mgp=c(2.5,.75,0), cex=.9)
panlab.cex <- .9
panlab.ln <- -1.2
CWbg <- mypal[6]#"white"
ylims <- c(-.6,.7)
### Effect Sizes #####
boxplot(value~xvals, avglong.ne[which(avglong.ne$trait=="LeafLife" & avglong.ne$SP.ID!="CWmean"),], at=c(1,2,3,4,5,6)
        , ylim=ylims, xaxt="n", xlim=c(0.5,6.5)
        ,outpch=1, border=paste0(mypal[colchoices[1]],"66"), col=paste0(mypal[colchoices[1]],"66"), outcol=paste0(mypal[colchoices[1]],"66")
        , staplewex=0, notch=F, medcol="black" #, medlwd=0
        , whisklty=1, whisklwd=2, boxwex=.7, whiskcol=paste0(mypal[colchoices[1]],"66") #, boxwex=1, medlwd=3
        ,las=2)
#plot(value~xvals, avglong.ne[which(avglong.ne$trait=="LeafLife" & avglong.ne$SP.ID!="CWmean"),], pch=16, xaxt="n", ylab="", ylim=c(-0.7,1), xlab="", xlim=c(0.6,6.4))
abline(h=0)
points(value~xvals, avglong.ne[which(avglong.ne$trait=="LeafLife" & avglong.ne$SP.ID=="CWmean"),], pch=24, bg=CWbg)
#mtext(text = "B)",side = 3,line = panlab.ln,adj=0.02,cex=panlab.cex )
mtext(text = "log(LL)",side = 3,line = panlab.ln,adj=0.9,cex=panlab.cex)
#text(x=7, y=0, pos=3, labels = "NA",cex=.7)
legend(x=0, y=1.15, legend = "w/in Spp variation", fill=paste0(mypal[colchoices[1]],"66"), bty="n",col = paste0(mypal[colchoices[1]],"66"),border = paste0(mypal[colchoices[1]],"66"), cex=1.1, xpd=NA)
legend(x=7, y=1.15, legend = "var in CW mean", bty="n",pch=24, cex=1.1, xpd=NA, pt.bg=CWbg)

boxplot(value~xvals, avglong.ne[which(avglong.ne$trait=="LMA" & avglong.ne$SP.ID!="CWmean"),], at=c(1,2,3,4,5,6)
        , ylim=ylims, xaxt="n", xlim=c(0.5,6.5)
        ,outpch=1, border=paste0(mypal[colchoices[1]],"66"), col=paste0(mypal[colchoices[1]],"66"), outcol=paste0(mypal[colchoices[1]],"66")
        , staplewex=0, notch=F, medcol="black" #, medlwd=0
        , whisklty=1, whisklwd=2, boxwex=.7, whiskcol=paste0(mypal[colchoices[1]],"66") #, boxwex=1, medlwd=3
        ,las=0)
#plot(value~xvals, avglong.ne[which(avglong.ne$trait=="LMA" & avglong.ne$SP.ID!="CWmean"),], pch=16, xaxt="n", ylab="", ylim=c(-0.7,1), xlab="", xlim=c(0.6,6.4))
abline(h=0)
points(value~xvals, avglong.ne[which(avglong.ne$trait=="LMA" & avglong.ne$SP.ID=="CWmean"),], pch=24, bg=CWbg)
#mtext(text = "D)",side = 3,line = panlab.ln,adj=0.02,cex=panlab.cex )
mtext(text = "log(LMA)",side = 3,line = panlab.ln,adj=0.9,cex=panlab.cex)
#text(x=7, y=0, pos=3, labels = "NA",cex=.7)
axis(side = 1,at = c(1,2,3,4,5,6), labels = c("Wetness", "Warmth","Soil N", "log(Age)","LAI","NPP"), xpd=NA, las=3, srt=50, xlim=c(0.6,6.4))
mtext(text="Standardized Effect Sizes",side = 2, adj=-.7, line=2.1)

boxplot(value~xvals, avglong.ne[which(avglong.ne$trait=="Nmass" & avglong.ne$SP.ID!="CWmean"),], at=c(1,2,3,4,5,6)
        , ylim=ylims, xaxt="n", xlim=c(0.5,6.5), yaxt="n"
        ,outpch=1, border=paste0(mypal[colchoices[1]],"66"), col=paste0(mypal[colchoices[1]],"66"), outcol=paste0(mypal[colchoices[1]],"66")
        , staplewex=0, notch=F, medcol="black" #, medlwd=0
        , whisklty=1, whisklwd=2, boxwex=.7, whiskcol=paste0(mypal[colchoices[1]],"66") #, boxwex=1, medlwd=3
        ,las=0)
# plot(value~xvals, avglong.ne[which(avglong.ne$trait=="Nmass" & avglong.ne$SP.ID!="CWmean"),], pch=16, xaxt="n", ylab="", ylim=c(-0.7,1), xlab="", xlim=c(0.6,6.4))
abline(h=0)
points(value~xvals, avglong.ne[which(avglong.ne$trait=="Nmass" & avglong.ne$SP.ID=="CWmean"),], pch=24, bg=CWbg)
#mtext(text = "F)",side = 3,line = panlab.ln,adj=0.02,cex=panlab.cex )
mtext(text = "log(Nmass)",side = 3,line = panlab.ln,adj=0.9,cex=panlab.cex)
#mtext(text="Standardized Effect Sizes",side = 2, adj=-.1, line=2.1)
#text(x=7, y=0, pos=3, labels = "NA",cex=.7)

boxplot(value~xvals, avglong.ne[which(avglong.ne$trait=="Narea" & avglong.ne$SP.ID!="CWmean"),], at=c(1,2,3,4,5,6)
        , ylim=ylims, xaxt="n", xlim=c(0.5,6.5), yaxt="n"
        ,outpch=1, border=paste0(mypal[colchoices[1]],"66"), col=paste0(mypal[colchoices[1]],"66"), outcol=paste0(mypal[colchoices[1]],"66")
        , staplewex=0, notch=F, medcol="black" #, medlwd=0
        , whisklty=1, whisklwd=2, boxwex=.7, whiskcol=paste0(mypal[colchoices[1]],"66") #, boxwex=1, medlwd=3
        ,las=0)
#plot(value~xvals, avglong.ne[which(avglong.ne$trait=="Narea" & avglong.ne$SP.ID!="CWmean"),], pch=16, xaxt="n", ylab="", ylim=c(-0.7,1), xlab="", xlim=c(0.6,6.4))
abline(h=0)
points(value~xvals, avglong.ne[which(avglong.ne$trait=="Narea" & avglong.ne$SP.ID=="CWmean"),], pch=24, bg=CWbg)
#mtext(text = "H)",side = 3,line = panlab.ln,adj=0.02,cex=panlab.cex )
mtext(text = "log(Narea)",side = 3,line = panlab.ln,adj=0.9,cex=panlab.cex)
#text(x=7, y=0, pos=3, labels = "NA",cex=.7)

axis(side = 1,at = c(1,2,3,4,5,6), labels = c("Wetness", "Warmth","Soil N", "log(Age)","LAI","NPP"), xpd=NA, las=3, srt=50, xlim=c(0.6,6.4))






####### R2 for Poster
quartz(width=4.33, height=2)
par(mar=c(.2,5,3,3), mgp=c(2,.75,0))
boxwidth = 1
boxplot(llbestmods.ne$r.squaredGLMM.R2m[-which(llbestmods.ne$SP.ID=="CWmean")]~rep(1, times=6), horizontal=T, at=4, xlim=c(0.5,4.75), ylim=c(0,1)
        , yaxt="n", xaxt="n"
        ,outpch=1, border=paste0(mypal[colchoices[1]],"66"), col=paste0(mypal[colchoices[1]],"66"), outcol=paste0(mypal[colchoices[1]],"66")
        , staplewex=0, notch=F, medcol="black" #, medlwd=0
        , whisklty=1, whisklwd=2, boxwex=boxwidth, whiskcol=paste0(mypal[colchoices[1]],"66") )
boxplot(lmabestmods.ne$r.squaredGLMM.R2m[-which(lmabestmods.ne$SP.ID=="CWmean")]~rep(1, times=6), horizontal=T, at=3, xlim=c(0,5), ylim=c(0,1)
        , yaxt="n", xaxt="n"
        ,outpch=1, border=paste0(mypal[colchoices[1]],"66"), col=paste0(mypal[colchoices[1]],"66"), outcol=paste0(mypal[colchoices[1]],"66")
        , staplewex=0, notch=F, medcol="black" #, medlwd=0
        , whisklty=1, whisklwd=2, boxwex=boxwidth, whiskcol=paste0(mypal[colchoices[1]],"66"), add=T )
boxplot(nmassbestmods.ne$r.squaredGLMM.R2m[-which(nmassbestmods.ne$SP.ID=="CWmean")]~rep(1, times=6), horizontal=T, at=2, xlim=c(0,5), ylim=c(0,1)
        , yaxt="n", xaxt="n"
        ,outpch=1, border=paste0(mypal[colchoices[1]],"66"), col=paste0(mypal[colchoices[1]],"66"), outcol=paste0(mypal[colchoices[1]],"66")
        , staplewex=0, notch=F, medcol="black" #, medlwd=0
        , whisklty=1, whisklwd=2, boxwex=boxwidth, whiskcol=paste0(mypal[colchoices[1]],"66"), add=T )
boxplot(nareabestmods.ne$r.squaredGLMM.R2m[-which(llbestmods.ne$SP.ID=="CWmean")]~rep(1, times=6), horizontal=T, at=1, xlim=c(0,5), ylim=c(0,1)
        , yaxt="n", xaxt="n"
        ,outpch=1, border=paste0(mypal[colchoices[1]],"66"), col=paste0(mypal[colchoices[1]],"66"), outcol=paste0(mypal[colchoices[1]],"66")
        , staplewex=0, notch=F, medcol="black" #, medlwd=0
        , whisklty=1, whisklwd=2, boxwex=boxwidth, whiskcol=paste0(mypal[colchoices[1]],"66"), add=T )


#plot(llbestmods.ne$r.squaredGLMM.R2m[-which(llbestmods.ne$SP.ID=="CWmean")]~rep(1, times=6), ylim=c(0,1), xlim=c(0,5), xaxt="n", ylab=expression(paste(R^2," or marginal ",R^2)), pch=16)
#points(lmabestmods.ne$r.squaredGLMM.R2m[-which(lmabestmods.ne$SP.ID=="CWmean")]~rep(2, times=6), pch=16)
#points(nmassbestmods.ne$r.squaredGLMM.R2m[-which(nmassbestmods.ne$SP.ID=="CWmean")]~rep(3, times=6), pch=16)
#points(nareabestmods.ne$r.squaredGLMM.R2m[-which(nareabestmods.ne$SP.ID=="CWmean")]~rep(4, times=6), pch=16)
points(c(4)~llbestmods.ne$r.squaredGLMM.R2c[which(llbestmods.ne$SP.ID=="CWmean")], pch=24,bg=CWbg) # the simple R^2 is stored in R2c for CWmean
points(c(3)~lmabestmods.ne$r.squaredGLMM.R2c[which(lmabestmods.ne$SP.ID=="CWmean")], pch=24,bg=CWbg)
points(c(2)~nmassbestmods.ne$r.squaredGLMM.R2c[which(nmassbestmods.ne$SP.ID=="CWmean")], pch=24,bg=CWbg)
points(c(1)~nareabestmods.ne$r.squaredGLMM.R2c[which(nareabestmods.ne$SP.ID=="CWmean")], pch=24,bg=CWbg)
# add in ecoregion only rsquareds
# points(c(4,3,2,1)~colMeans(ecoregmods[,-1]), pch=4,bg="white") # the simple R^2 is stored in R2c for CWmean
# points(c(4,3,2,1)~CWecoregmods, pch=23,bg="white") # the simple R^2 is stored in R2c for CWmean



axis(3, mgp=c(2.5,.2,0),tck=.05)
axis(2, labels = c("Narea","Nmass","LMA","Leaf Life"), at=c(1,2,3,4), las=2, cex.axis=1.3)
mtext(side=3,text = expression(paste(R^2, "of conifer Trait-Env relationships")), line=1.1, cex=1.3, at=0.3)
#legend(x=.675, y=1, legend = "Ind. Spp", fill="gray", bty="n",col = paste0(mypal[colchoices[1]],"66"),border = paste0(mypal[colchoices[1]],"66"), cex=.8)
#legend(x=.7, y=1.5, legend = "CW mean", bg="white", bty="n",pch=24, cex=.8)
# legend(x=.7-.1, y=3+2.5, legend = c("CWM Ecoreg", "Spp Ecoreg"), bty="n", pch=c(23,4), cex=.8)

# legend(x=.75, y=2.5, legend = c("CWM-eco","CWM-full","Spp-eco"), bg="white", bty="n",pch=c(23,24,4), cex=.7)
# legend(x=.73, y=.9, legend = "Spp-env", fill="gray", bty="n",col = paste0(mypal[colchoices[1]],"66"),border = paste0(mypal[colchoices[1]],"66"), cex=.7)
#mtext("A)", cex=panlab.cex, side=3, adj=0.02, line=panlab.ln)













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


fit.ecoreg <- function( dataz, Species){
  datazall <- dataz[which(dataz$SP.ID==Species)]
}

lleco <- lmer(get(traitsc)~climPC1sc + climPC2sc + soil_Nsc + log.ASAsc + LAI_Osc + AG_TGROWTHsc + (1|PLOT_ID), datazall, REML=F)




ggplot(traits.common[which(traits.common$GE.SP %in% c("Pseudotsuga.menziesii","Pinus.ponderosa")),], aes(x=ECOREGION, y=log.LL, col=ECOREGION)) + geom_boxplot() + facet_grid(GE.SP~.)
# it actually looks like ecoregion might explain a fair amount of LL variation, and maybe Nmass
# though the effect is not consistent between ponderosa and doug fir, maybe data driven param but not mechanistic?
plot(log.LMA~log.Nmass, traits.common[which(traits.common$GE.SP %in% c("Pseudotsuga.menziesii")),], col=ECOREGION, pch=16)
# also, ecoregions cluster rather suprisingly in log.LMA vs log.Nmass space...
# good for doug fir, less good for ponderosa




##### old function for cw traits model ensambling. Something never worked here

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





