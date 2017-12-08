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
# stand level averages
mpotr <- potr %>% group_by(PlotUn) %>% summarise(LMA=mean(LMA, na.rm=T), LLmonths=mean(LLmonthselev))
mpotr$log.LMA <- log(mpotr$LMA, base=10)
mpotr$log.LL <- log(mpotr$LLmonths, base=10)




######## Combine individual level data for all supplemental datasets ########
suppdata <- rbind(data.frame(coffee.ind), data.frame(euc.ind), data.frame(quga.ind), data.frame(potr.ind))
#write.csv(suppdata,"/Users/leeanderegg/Dropbox/NACP_Traits/Intra-data/OtherData_Combined_033017.csv")
write.csv(suppdata,"/Users/leeanderegg/Dropbox/NACP_Traits/Intra-data/OtherData_Combined_040117.csv")



### bill's weird potr values
# bll <- c(244.375,209.1818182,176.6666667,179.6153846,168.6153846,191.9230769,182.3076923,180.1538462,192.8461538)
# bsla <- c(0.009281022,0.007095145,0.007528245,0.006679047,0.002422823,0.009146608,0.005097827,0.002286325,0.006019923)
# blma <- 1/bsla
# bllmonths <- bll/30
# blog.LMA <- log(blma, base=10)
# blog.LL <- log(bllmonths, base=10)
# 
# plot(log.LL~log.LMA, LES, pch=16, col='grey')
# points(log.LLelev~log.LMA, potr, pch=16)
# points(log.LL~log.LMA, LES[which(LES$Species=="Populus tremuloides"),], pch=16, col='blue')
# points(blog.LL~blog.LMA)
# legend("topleft", pch=c(16,16,1), col=c('blue','black','black'), legend = c("GLOPNET","Lee","Bill"))




plot(SLA~Nmass, LES, col=Acats, pch=16)
 abline(lm(SLA~Nmass, LES[which(LES$Acats=="H"),]), col=brewer.pal(n=4, "RdYlBu")[1], lwd=2)
 abline(lm(SLA~Nmass, LES[which(LES$Acats=="MH"),]), col=brewer.pal(n=4, "RdYlBu")[2], lwd=2)
 abline(lm(SLA~Nmass, LES[which(LES$Acats=="ML"),]), col=brewer.pal(n=4, "RdYlBu")[3], lwd=2)
 abline(lm(SLA~Nmass, LES[which(LES$Acats=="L"),]), col=brewer.pal(n=4, "RdYlBu")[4], lwd=2)
 legend(x=0, y=.095,xpd=NA,title = "Amax/area",legend = c("H","MH","ML","L"), pch=16, col=brewer.pal(n=4, "RdYlBu"), lwd=2, bty="n",ncol=4)
 
 points(SLA_drymass~LEAF_NITROGEN, traits[which(traits$SP.ID=="PSEMEN"),], col="darkblue")
 points(SLA_drymass~LEAF_NITROGEN, traits[which(traits$SP.ID=="PINPON"),], col="darkred")
 points(SLA_drymass~LEAF_NITROGEN, traits[which(traits$SP.ID=="TSUHET"),], col="darkgreen")
 
 ### plot Jen's way of looking at things
ggplot(LES, aes(x=Nmass, y=SLA)) + geom_point(col="grey") + 
  geom_point(data=LES[which(LES$Species%in%commonspp),], aes(col=Species)) + geom_smooth(data=LES[which(LES$Species%in%commonspp),], aes(col=Species), method="lm",se=F) +
  geom_point(data=traits[which(traits$SP.ID %in% names(which(xtabs(~SP.ID, traits)>5))),], aes(x=LEAF_NITROGEN, y=SLA_drymass, col=SP.ID)) + 
  geom_smooth(data=traits[which(traits$SP.ID %in% names(which(xtabs(~SP.ID, traits)>5))),], aes(x=LEAF_NITROGEN, y=SLA_drymass, col=SP.ID), method="lm", se=F) +
  theme(legend.position="none")
  
### Simple model to explain SLA
mod1 <- lm(SLA~Nmass * Aarea * LL, LES[-which(is.na(LES$SLA) | is.na(LES$Nmass) | is.na(LES$LL) | is.na(LES$Aarea)),])

library(MuMIn)
options(na.action="na.fail")
tmp <- dredge(mod1)
  # whoa shit. it looks like all of these interactions are needed in the best model...

ndAarea <- data.frame(Nmass=rep(seq(.5,6, by=.5), times=3), Aarea=c(rep(28, times=12), rep(13, times=12), rep(5, times=12)), LL=rep(13, times=12*3))
predsAarea <- predict(mod1,newdata = ndAarea)
 
  #                 Species           Family    BIOME C3C4
  #     Veronica chamaedrys   Plantaginaceae   ALPINE   C3
  # Brachypodium distachyon          Poaceae  GRASS/M   C3
  #           Crepis sancta       Asteraceae  GRASS/M   C3
  #    Drosera rotundifolia      Droseraceae   TUNDRA   C3
  # Gymnocarpium dryopteris      Woodsiaceae   TUNDRA   C3
  #   Melampyrum sylvaticum    Orobanchaceae   TUNDRA   C3
  #     Pinguicula vulgaris Lentibulariaceae   TUNDRA   C3
  #      Erodium cicutarium      Geraniaceae   BOREAL     
  # Sanguisorba officinalis         Rosaceae   BOREAL     
  #     Veronica chamaedrys   Plantaginaceae   BOREAL   C3
  #         Viola mirabilis        Violaceae   BOREAL     
  #        Linnaea borealis      Linnaeaceae   BOREAL     
  #      Acer pensylvanicum      Sapindaceae TEMP_FOR   C3


####### Mesophyll conductance ########

# mesophyll conductance data from Muir 2016
mesdat <- read.csv("/Users/leeanderegg/Dropbox/NACP_Traits/gmes/Muir2016_LMAvGmes.csv")

# match with any overlap from GLOPNET
# mesdat$Aarea <- LES$Aarea[match(mesdat$species,LES$Species)]
# mesdat$Amass <- LES$Amass[match(mesdat$species,LES$Species)] 
# actually more overlap from Maire et al. 2015
Maire <- read.csv("/Users/leeanderegg/Dropbox/NACP_Traits/globamax_data_160609.csv", header=T)
mesdat$Aarea <- Maire$Aarea[match(mesdat$species,Maire$Genus.spp)]
mesdat$Amass <- Maire$Amass[match(mesdat$species,Maire$Genus.spp)]
mesdat$SLA <- 1/mesdat$lma

#Maire$Genus.spp[match(mesdat$species,Maire$Genus.spp)]
#LES$Species[match(mesdat$species,LES$Species)] 

# plot of lma-specific gm versus Am
ggplot(mesdat, aes(x=Aarea, y=gm/lma, col=Amass)) + geom_point()
# or just as a function of A per unit mass
ggplot(mesdat, aes(x=Amass, y=gm, col=Aarea)) + geom_point()
# gm vs sla w/ Am as color
quartz(width=4,height=3.5)
ggplot(mesdat, aes(x=SLA, y=gm, col=Aarea)) + geom_point() +
  geom_line(data=predictions[which(predictions$Aarea==5),], col="black") +
  geom_line(data=predictions[which(predictions$Aarea==12),], col="darkblue") +
  geom_line(data=predictions[which(predictions$Aarea==20),], col="lightblue") 

  
ggplot(mesdat, aes(x=lma, y=gm)) + geom_point() 


mestest <- lm(gm ~ SLA * Aarea, mesdat)
slas <- seq(0.00,0.05, by=.01)
preds <- predict(mestest,newdata = list(SLA=c(slas, slas, slas), Aarea = c(rep(5, times=length(slas)), rep(12, times=length(slas)), rep(20, times=length(slas)))))
predictions <- data.frame(gm = preds, SLA = c(slas, slas, slas),Aarea = c(rep(5, times=length(slas)), rep(12, times=length(slas)), rep(20, times=length(slas))))

####### Random test of LMA~Aarea relationship in Brachychitons

brach <- read.csv("/Users/leeanderegg/Dropbox/Trade-offs project/Tori's Stuff/Brachychiton_Traits_Data_final031816.csv")
brach$LMA <- 1/brach$SLA
brach$Amass <- brach$A * brach$SLA # A/cm2? * cm2/g = A/g

alec <- read.csv("/Users/leeanderegg/Dropbox/Aspen_Leaf_Manuscript/Manuscript DataFiles/MasterDataFile_ASB.csv")
alec$Amass <- alec$amax / alec$lma # umol m-2 s-1 / g/cm2 = A/g
levels(alec$treatment) <- list(C="C",LW="LW",LL="LL")

plot(Amass~LMA, brach, col=Species, ylim=c(0,300), xlim=c(0.01,.15), pch=as.numeric(Treatment.x))
points(Amass/10~I(lma*10), alec, pch=as.numeric(treatment))
points(amass*1000~I(lma_g_m2/1000), coffee, col="green3")


plot(Amass~LMA, brach, col=Species, pch=as.numeric(Treatment.x), ylim=c(0,15), xlim=c(0.01,0.15))
points(amax~I(lma*10), alec, pch=as.numeric(treatment))
points(amax~I(lma_g_m2/1000), coffee, pch=4, col="blue4")

plot(A~LMA, brach, col=Species, pch=as.numeric(Treatment.x), ylim=c(0,15), xlim=c(0.01,0.15))
points(amax~I(lma*10), alec, pch=as.numeric(treatment))
points(amax~I(lma_g_m2/1000), coffee, pch=4, col="blue4")

ggplot(brach, aes(x=LMA, y=Amass, col=Species, pch=Treatment.x)) + geom_point() + geom_smooth(method="lm")



##### Comparison of my results with Blonder's results ######

require(lmodel2)
plot.MAR <- function(xvar, yvar, data, method="SMA", linecol, lwd=1) {
  if(method=="SMA") meth = 3
  if(method=="MA") meth = 2
  if(method=="OLS") meth =1
  if(length(t(as.vector(data[!(is.na(data[,xvar]) | is.na(data[,yvar])), xvar])))<3){
    #return(rep(NA, times=7))
    break()
  }
  else{
    if(var(data[,yvar], na.rm=T)==0){
      break()
    }
    else{
      tmp.mod <- suppressMessages(lmodel2(get(yvar)~get(xvar), data))
      intercept <- tmp.mod$regression.results$Intercept[meth]
      slope <- tmp.mod$regression.results$Slope[meth]
      yhat <- tmp.mod$x * slope + intercept
      lines(yhat~tmp.mod$x, col=linecol, lwd=lwd)
      
    }
  }
}



plot(log.LL~log.LMA, allspp, col="grey", pch=16, ylab="log(LL)", xlab="log(LMA)")
abline(a=all.results.cl$Int_LMA.LL[nrow(all.results.cl)], b=all.results.cl$Slope_LMA.LL[nrow(all.results.cl)], lwd=3, col="#666666", lty=3)
#abline(a=all.results.cl$Int_LMA.LL[nrow(all.results.cl)-1], b=all.results.cl$Slope_LMA.LL[nrow(all.results.cl)-1], lwd=3, col="black")
palette(paste0(genpal,"88"))
points(log.LL~log.LMA, spp.data[which(spp.data$Species %in% all.results.cl$Taxo.Unit[which(all.results.cl$n_LMA.LL>5 & all.results.cl$Type=="w.inSpp")]),], pch=16)#, col=factor(Species))
palette(genpal)
tax <- "w.inSpp"
for (j in 1:length(as.character(all.results.cl$Taxo.Unit[which(all.results.cl$Type==tax & all.results.cl$n_LMA.LL>5)]))){
  i <- as.character(all.results.cl$Taxo.Unit[which(all.results.cl$Type==tax & all.results.cl$n_LMA.LL>5)])[j]
  #plot.MAR(xvar = "log.LMA", yvar = "log.LL",data= spp.data[which(spp.data$Species==i),], linecol = genpal[j], lwd=3)
  plot.MAR(xvar = "log.LMA", yvar = "log.LL",data= spp.data[which(spp.data$Species==i),], linecol = mypal[colchoices[1]], lwd=3)
  
}
# arabidopsis line
points(log.LL~log.LMA, arab, pch=1)
plot.MAR(xvar="log.LMA", yvar="log.LL", data=arab, linecol=mypal[colchoices[2]],lwd = 3)
legend("bottomright", legend = c("GLOPNET","w/in SPP geographic", "Blonder 2015"), lwd=c(3,3,3), lty=c(3,1,1), col = c("#666666",mypal[colchoices[1]],mypal[colchoices[2]]), cex=.6,pt.cex = 1, bty="n")

