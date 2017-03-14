arab <- read.csv("/Users/leeanderegg/Dropbox/NACP_Traits/Blonder2015Arabidopsis.csv",header=T)
colnames(arab)<- c("Type","ID","Genotype","LL","LMA","Nmass","Amass","VD","LDMC")
arab$SLA <- 1/arab$LMA
arab$Aarea <- arab$Amass * arab$LMA
arab$Narea <- arab$Nmass * arab$LMA

LES$LMA <- 10^LES$log.LMA
LES$Nmass <- 10^LES$log.Nmass
LES$Amass <- 10^LES$log.Amass
LES$Aarea <- 10^LES$log.Aarea
traits$SLA_drymass <- 1/traits$LMA


LES$SLA <- 1/LES$LMA
plot(SLA~Nmass, LES)


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
mesdat$Aarea <- LES$Aarea[match(mesdat$species,LES$Species)]
mesdat$Amass <- LES$Amass[match(mesdat$species,LES$Species)] 

# plot of lma-specific gm versus Am
ggplot(mesdat, aes(x=Aarea, y=gm*lma, col=Aarea)) + geom_point()
# or just as a function of A per unit mass
ggplot(mesdat, aes(x=Amass, y=gm, col=Aarea)) + geom_point()
# gm vs sla w/ Am as color
ggplot(mesdat, aes(x=1/lma, y=gm, col=Aarea)) + geom_point()

