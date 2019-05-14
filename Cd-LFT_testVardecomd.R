require(tidyr)
require(lme4)
require(lmerTest)
require(stringr)
require(RColorBrewer)

##### Variance Decomp of LES dataset for LFT analysis ######


mypal <- brewer.pal(n=9, "Set1")
palette(mypal)
# set the colors that wil denote within-species, within-genus, within family and across CWMs
colchoices <- c(1,2,4,3,6)


## load LES from Cd-Initial_Analysis.R
higher.taxa <- function(fam.names){
  output <- data.frame(Order=rep(NA, length(fam.names)),Subclass=rep(NA, length(fam.names)),Class=rep(NA, length(fam.names)), Subdivision=rep(NA, length(fam.names)),Division=rep(NA, length(fam.names)))
  for(i in 1:length(fam.names)){
    tmp <- taxize::classification(fam.names[i], db="itis")[[1]]
    output[i,] <- cbind(tmp$name[which(tmp$rank=="order")], tmp$name[which(tmp$rank=="subclass" | tmp$rank=="superorder")], tmp$name[which(tmp$rank=="class")], tmp$name[which(tmp$rank=="subdivision")], tmp$name[which(tmp$rank=="division")])
  }
  return(output)
}
LEStaxa <- higher.taxa(fam.names=c("Pinaceae","Fagaceae"))


### Thinking about log-transformations and variance decomps:
# most traits are log-normal at global, family and genus scales. However, they look more normally destributed w/in species
  # certainly true for LL
  # defo true for Nmass   
  # LMA is maybe more log-normal within some species (but not all)
  # Al.As seems maybe somewhat log-normal within species


#logLMAvar <- lmer(log.LMA~1 + (1|Family) + (1|Genus), LES)
logLMAvar <- lmer(log.LMA~1 + (1|Family/Genus), LES)
logLLvar <- lmer(log.LL~1 +  (1|Family/Genus), LES)
logNmassvar <- lmer(log.Nmass~ 1 + (1|Family/Genus) ,LES)
logNareavar <- lmer(log.Narea~ 1 + (1|Family/Genus), LES)
logAmassvar <- lmer(log.Amass~1 + (1|Family/Genus), LES)
logAareavar <- lmer(log.Aarea~1 + (1|Family/Genus), LES)
logRdmassvar <- lmer(log.Rdmass~1 + (1|Family/Genus), LES)
logRdareavar <- lmer(log.Rdarea~1 + (1|Family/Genus), LES)
cacivar <-  lmer(`Ca...Ci`~1 + (1|Family/Genus), LES)

LMAvariance <- data.frame(VarCorr(logLMAvar))
LLvariance <- data.frame(VarCorr(logLLvar))
Nmassvariance <- data.frame(VarCorr(logNmassvar))
Nareavariance <- data.frame(VarCorr(logNareavar))
Amassvariance <- data.frame(VarCorr(logAmassvar))
Aareavariance <- data.frame(VarCorr(logAareavar))
Rdmassvariance <- data.frame(VarCorr(logRdmassvar))
Rdareavariance <- data.frame(VarCorr(logRdareavar))
cacivariance <- data.frame(VarCorr(cacivar))

traitvars <- data.frame(LMAvariance[,4], LLvariance[,4], Nmassvariance[,4], Nareavariance[,4], Amassvariance[,4],Aareavariance[,4],Rdmassvariance[,4],Rdareavariance[,4],cacivariance[,4] )
colnames(traitvars) <- c("logLMA", "logLL", "logNmass", "logNarea", "logAmass", "logAarea","logRdmass","logRdarea","caci")
rownames(traitvars) <- c("BtwGenera", "BtwFamilies", "BtwSpecies")

traitvars_scaled <- traitvars  
for(i in 1:ncol(traitvars)){
  traitvars_scaled[,i] <- traitvars[,i]/sum(traitvars[,i])
}


traitvars_scaled2 <- traitvars_scaled[c(2,1,3),]

traitvars_scaled3 <- traitvars_scaled2[,order(traitvars_scaled2[1,], decreasing =T)]

#quartz(width=3.75, height=4)
quartz(width=6.81, height=3.5)
par(mgp=c(3,.7,0), cex.lab=1.3, cex.axis=1.1, mfrow=c(1,1), mar=c(5.5,2,3,6), oma=c(0,3,0,0))
cols <- rev(c(mypal[colchoices[c(1,2,3)]],"#CCCCCC"))#brewer.pal(11, "RdBu")[c(11,9,1, 6)]
# With old colors, alpha=CC, with new colors, alpha=99
bp <-barplot(as.matrix(traitvars_scaled3),beside=F,legend.text = F,xpd = T, names.arg = colnames(traitvars_scaled3),las=2,args.legend = list(x=4, y=1.3, ncol=2), col = paste0(cols,"99")
             , ylab= "Proportion of total Variance" #ylab="Proportion of total Variance\n(log traits)"
             , xlab="", mgp=c(2.5,.8,0))
#legend(xpd=T, x = 0, y=1.2, legend=rev(c("btw Fams","w/in Fam","w/in Gen","w/in Spp")), fill=rev(paste0(cols,"99")), ncol=2, bty="n",  cex=1.2)
#text(x = bp, y= par("usr")[3]-.05,labels =  c(expression(paste(log[10](LMA))), expression(paste(log[10](LL))),expression(paste(log[10](N[mass]))),expression(paste(log[10](N[area])))), srt=40, adj=1,xpd=NA, cex=1.4, font=1)
#mtext("a)", side=3, adj=-.1, line=1.3)
mtext("Proportion of Total Variance", side=2, line=2.8)
mtext("GLOPNET data", side=3, line=0.3)
#quartz(width=5, height=4)
#cols <- brewer.pal(11, "RdBu")[c(11,9,1, 6)]
legend(xpd=T, x = 11, y=0.7, legend=rev(c("btw Fams","btw Gen","btw Spp")), fill=rev(paste0(cols[1:3],"99")), ncol=1, bty="n",  cex=1.2)
# mtext(text="28 spp w/ replication, 1000+ spp,\n500+ genera, 150+ families",side = 1,line = 3.3)
# legend(xpd=T, x = 0, y=1.3, legend=rownames(rtraitvars_scaled2), fill=paste0(cols,"CC"), ncol=2, bty="n",  cex=1.2)





#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
###### load BAAD data
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

baad.data::baad_data_version_current()
  # currently 1.0.1 on 10/02/17
baad.all <- baad.data::baad_data()
baad.all$dictionary[,c("variable", "label", "units")]

baad <- baad.all$data
baad$genus <- sapply(X=baad$speciesMatched, FUN = function(x){strsplit(x, split = " ")[[1]][1]})
baad$studyName <- factor(baad$studyName)
baad$genus <- factor(baad$genus)
baad$family <- factor(baad$family)
baad$Al.As <- baad$a.lf/baad$a.ssbh
baad$Ab.Bl <- baad$m.so/baad$m.rt
baad$RSR <- baad$m.rt/baad$m.so # making root:shoot ratio for matching with Ledo 2017 below
baad$lmf <- baad$m.lf / (baad$m.lf + baad$m.st)
baad$taper <- baad$d.ba/baad$h.t
baad$d.bh[which(is.na(baad$d.bh))] <- 4/5 * baad$d.ba[which(is.na(baad$d.bh))]



#colnames(baad)[which(colnames(baad) %in% c("family"))] <- c("Family", "a.lf")

# a.lf = leaf area
# a.ssbh = sapwood area at breast height
# h.t = height
# m.lf = leaf mass
# m.st = stem mass
# m.so = aboveground mass
# m.rf/c = mass of fine/coarse roots
# m.rt = root mass
# r.st = wood density
# r.ss = sapwood density
# r.sh = heartwood desnity

# ### Al.As
# xtabs(~factor(family), baad[which(!is.na(baad$Al.As)),])
#   # 10 families
# xtabs(~factor(genus), baad[which(!is.na(baad$Al.As)),])
#   #13 genera (might be perfectly nested)
# xtabs(~factor(speciesMatched), baad[which(!is.na(baad$Al.As)),])
#   #20 spp, with multiple pines, oaks and 2 eucs
# 
# ### wood density
# length(xtabs(~factor(family), baad[which(!is.na(baad$r.st)),]))
# # 66 families
# length(xtabs(~factor(genus), baad[which(!is.na(baad$r.st)),]))
# # 121 genera (might be perfectly nested)
# length(xtabs(~factor(speciesMatched), baad[which(!is.na(baad$r.st)),]))
# # 179 spp, so some congeners
# 
# ### aboveground:belowground (Ab.Bl)
# length(xtabs(~factor(family), baad[which(!is.na(baad$Ab.Bl)),]))
# # 87 families
# length(xtabs(~factor(genus), baad[which(!is.na(baad$Ab.Bl)),]))
# #174 genera (might be perfectly nested)
# length(xtabs(~factor(speciesMatched), baad[which(!is.na(baad$Ab.Bl)),]))
# #261 spp
# 
# ### leaf mass fraction (lmf)
# length(xtabs(~factor(family), baad[which(!is.na(baad$lmf)),]))
# # 113 families
# length(xtabs(~factor(genus), baad[which(!is.na(baad$lmf)),]))
# # 314 genera (might be perfectly nested)
# length(xtabs(~factor(speciesMatched), baad[which(!is.na(baad$lmf)),]))
# #535 spp
# 
# ### height (h.t)
# length(xtabs(~factor(family), baad[which(!is.na(baad$h.t)),]))
# # 105 families
# length(xtabs(~factor(genus), baad[which(!is.na(baad$h.t)),]))
# # 309 genera (might be perfectly nested)
# length(xtabs(~factor(speciesMatched), baad[which(!is.na(baad$h.t)),]))
# # 584 spp
# 
# ## see if things are normal or lognormal
# hist(log(baad$Al.As)) #lognormal
# hist(log(baad$r.st)) # kinda normal, kinda lognormal
# hist(log(baad$Ab.Bl)) # very lognormal
# hist((baad$lmf)) # uniform. not sure how to deal with this...
# hist((baad$h.t)) # lognormal

baad$log.Al.As <- log(baad$Al.As, base=10)
baad$log.r.st <- log(baad$r.st, base=10)
baad$log.Ab.Bl <- log(baad$Ab.Bl, base=10)
baad$log.RSR <- log(baad$RSR, base=10)
baad$log.h.t <- log(baad$h.t, base=10)
baad$lm.dbh <- with(baad, m.lf/(3.1415 *(d.bh/2)^2))
baad$log.lm.dbh <- log(baad$lm.dbh, base=10)

#### plots of attributes with size:
quartz(width=6, height=3)
par(mfrow=c(1,4), mar=c(3,3,1,0), oma=c(0,0,0,1), mgp=c(2,1,0), cex=1)
plot(m.lf~d.bh, baad)
plot(log(m.lf)~log(d.bh), baad)
#abline(a=0,b=1, col="red")
plot(lmf~log(d.bh), baad)
#plot(log(Al.As)~log(d.bh), baad, xlim=c())
plot(Al.As~d.bh, baad, xlim=c(0,1), col=studyName)
#abline(a=0,b=1, col="red")
plot(log.lm.dbh~map, baad)


#### Variance Decomp for BAAD
WDvar <- lmer(r.st~ 1 + (1|family/genus/speciesMatched), baad)
logWDvar <- lmer(log.r.st~ 1 + (1|family/genus/speciesMatched), baad)
logAbBlvar <- lmer(log.Ab.Bl~ 1 + (1|family/genus/speciesMatched), baad)
logHeightvar <- lmer(h.t~ 1 + (1|family/genus/speciesMatched), baad)
LMFvar <- lmer(lmf~ 1 + (1|family/genus/speciesMatched), baad)
lm.dbhvar <- lmer(log.lm.dbh ~ 1 + (1|family/genus/speciesMatched), baad)

lm.dbhvar <- lmer(log.lm.dbh ~ 1 + (1|family) + (1|genus) + (1|speciesMatched), baad)


rawWDvariance <- data.frame(VarCorr(WDvar))
WDvariance <- data.frame(VarCorr(logWDvar))
AbBlvariance <- data.frame(VarCorr(logAbBlvar))
Heightvariance <- data.frame(VarCorr(logHeightvar))
LMFvariance <- data.frame(VarCorr(LMFvar))


write.csv(LMFvariance, "/Users/leeanderegg/Dropbox/NACP_Traits/NACP_Traits_Rcode/LMF_VarianceDecomp_053118.csv")



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
####### Global Wood Density Database ######
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

GWD <- read.csv("/Users/leeanderegg/Dropbox/WD project/GlobalWoodDensityDatabase.csv", header=T)
GWD$Genus <- factor(sapply(X=as.character(GWD$Binomial), FUN = function(x){strsplit(x, split = " ")[[1]][1]}))
# bunch of replicated genera, and 100+ families.

colnames(GWD)[4] <- "WD"

# #WDvar2 <- lmer(WD~ 1 + (1|Genus/Binomial), GWD)
# WDvar2 <- lmer(WD~ 1 + (1|Family) + (1|Genus) + (1|Binomial), GWD)
# logWDvar2 <- lmer(log(WD, base=10)~ 1 + (1|Family) + (1|Genus) + (1|Binomial), GWD)
# 
# rawWDvariance2 <- data.frame(VarCorr(WDvar2))
# WDvariance2 <- data.frame(VarCorr(logWDvar2))




##### ####### Al:As & WD variance decomp for Anna's New Phyt Review #############
colnames(GWD)
colnames(baad)
wdbaad <- baad %>% filter(r.st>0) %>% select(WD=r.st, Binomial=speciesMatched, Family=family, Genus=genus, Study=studyName) %>% mutate(WD=WD/1000)
wdGWD <- GWD %>% select(WD, Binomial, Family, Genus,Study=Reference.Number)
wdbaad$Study <- as.character(wdbaad$Study)
wdGWD$Study <- as.character(wdGWD$Study)
allwd <- rbind(wdbaad, wdGWD)

WDvar3 <- lmer(WD~ 1 + (1|Family) + (1|Genus) + (1|Binomial), allwd)
logWDvar3 <- lmer(log(WD, base=10)~ 1 + (1|Family) + (1|Genus) + (1|Binomial), allwd)

rawWDvariance3 <- data.frame(VarCorr(WDvar3))
rawWDvariance3$scaledVar <- rawWDvariance3$vcov/ sum(rawWDvariance3$vcov)
WDvariance3 <- data.frame(VarCorr(logWDvar3))



## make a flag for how Al.As was calcualted
baad$Al.As_type <- NA
baad$Al.As_type[which(!is.na(baad$Al.As))] <- "direct.dbh"
baad$Al.As.others <- NA
baad$Al.As.others[which(is.na(baad$a.ssbh) & !is.na(baad$a.ssba))] <- with(baad[which(is.na(baad$a.ssbh) & !is.na(baad$a.ssba)),], a.lf/a.ssba)
baad$Al.As_type[which(is.na(baad$a.ssbh) & !is.na(baad$a.ssba))] <- "direct.basal"
baad$Al.As_type[which(is.na(baad$a.ssbh) & !is.na(baad$a.stbh) & !is.na(baad$a.lf) & is.na(baad$Al.As.others))] <- "nosw.dbh"
baad$Al.As.others[which(is.na(baad$a.ssbh) & !is.na(baad$a.stbh) & !is.na(baad$a.lf) & is.na(baad$Al.As.others))] <- with(baad[which(is.na(baad$a.ssbh) & !is.na(baad$a.stbh) & !is.na(baad$a.lf) & is.na(baad$Al.As.others)),], a.lf/(a.stbh))
baad$Al.As.all <- baad$Al.As
baad$Al.As.all[which(is.na(baad$Al.As))] <- baad$Al.As.others[which(is.na(baad$Al.As))]

# this new Al.As.all has 3798 measurements, with 234 species, 68 studies, 136 genera, 62 families
AlAsvar <- lmer(log(Al.As.all)~ Al.As_type + (1|family/genus/speciesMatched), baad)
AlAsvariance <- data.frame(VarCorr(AlAsvar))
AlAsvariance$scaledVar <- AlAsvariance$vcov/ sum(AlAsvariance$vcov)


AlAsvar1 <- lmer(log(Al.As.all)~ Al.As_type + (1|family/genus/speciesMatched), baad[which(!is.na(baad$a.stbh)),])
AlAsvar2 <- lmer(log(Al.As.all)~ log(a.stbh) + Al.As_type + (1|family/genus/speciesMatched), baad[which(!is.na(baad$a.stbh)),])
AlAsvar3 <- lmer(log(Al.As.all)~ log(a.stbh) + Al.As_type + (log(a.stbh)|family/genus/speciesMatched), baad[which(!is.na(baad$a.stbh)),])


new.vardecomps.forann <- data.frame(rawWDvariance3[c(3,2,1,4),c("grp","scaledVar")], AlAsvariance$scaledVar[c(3,2,1,4)])
# note: as of 4/29/19 i just copied the new Al_As column to the old spreadsheet

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
######## R:S database (2017) ###########
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


RtS <- read.csv("/Users/leeanderegg/Desktop/LFT analysis/nph14863-sup-0002-NotesS1.csv", header=T)
RtS$vegetationTYPE[which(RtS$vegetationTYPE=="TropMF ")] <- "TropMF"
RtS$vegetationTYPE[which(RtS$vegetationTYPE=="PlantTem ")] <- "PlantTem"
RtS$vegetationTYPE <- factor(RtS$vegetationTYPE)
RtS$log.AGB <- log(RtS$AGB_kg, base=10)
RtS$BGB_kg[which(RtS$BGB_kg==0)] <- NA # there are 94 measurements with 0 for BGB
RtS$log.BGB <- log(RtS$BGB_kg, base=10)
RtS$RS_RATIO[which(RtS$RS_RATIO==0)] <- NA # there are 94 measurements with 0 for BGB
RtS$log.RSR <- log(RtS$RS_RATIO, base=10)
RtS$taper <- RtS$DBH_cm/RtS$H_m
# 
# RtSvar <- lmer(RS_RATIO~ 1 + (1|FAMILY/GENUS/SPECIES), RtS)
# logRtSvar <- lmer(log.RSR ~ 1 + (1|FAMILY/GENUS/SPECIES), RtS[which(!is.na(RtS$SPECIES)& !is.na(RtS$GENUS) & !is.na(RtS$FAMILY)),])
# RtSvariance <- data.frame(VarCorr(logRtSvar))
# 


pal <- brewer.pal(n=9, "Set1")
palette(pal)

plot(log.AGB~log.BGB, RtS, col=vegetationTYPE)
for(i in 1:length(levels(RtS$vegetationTYPE))){
  plot.MAR(xvar = "log.BGB", yvar="log.AGB", data = RtS[which(RtS$vegetationTYPE==levels(RtS$vegetationTYPE)[i]),], linecol = mypal[i], lwd=2 )
}

ggplot(RtS, aes(x=log.BGB, y=log.AGB, col= vegetationTYPE)) + geom_point() + geom_smooth(method="lm") +
  facet_wrap(facets= ~vegetationTYPE)


xtabs(~FAMILY, RtS)






######## Xylem Functional Traits Database XFT #################

require(data.table)
TRYdata <- fread("/Users/leeanderegg/Desktop/LFT analysis/XFT_database/6275.txt", header = T, sep = "\t", dec = ".", fill=T, quote = "", data.table = T)

# strip out all of the info besides traits and a few metadata I want:
TRY <- TRYdata[which(!is.na(TRYdata$TraitID) | TRYdata$DataID %in% c(1861,	413,	2539,	2541,	193,	59,	60)),]
  # note: there should only by 15353 unique obs in the final dataset

TRYmeta <- TRY[which(is.na(TRY$TraitID)),]
TRYmeta$DataName <- factor(TRYmeta$DataName)
levels(TRYmeta$DataName) <- list(Organ = "Plant organ measured" , Dev_stage= "Plant developmental status / plant age / maturity / plant life stage", Curve_Type= "Curve Type", P50_Method="P50 method", Biome="Vegetation type / Biome", Lat= "Latitude", Lon = "Longitude" )
TRYmeta.wide <- spread(TRYmeta %>% select(ObservationID, DataName,OrigValueStr), key = DataName, value = OrigValueStr)

TRYtraits <- TRY[which(!is.na(TRY$TraitID)),]
TRYtraits.trim <- TRYtraits %>% select(1,7,8,10,11,15,16,25)
TRYtraits.long <- TRYtraits.trim %>% gather(key="Type","value",OrigValueStr:ErrorRisk)




##### Using the formatted dataset Bill sent me (might be older version) ####
  # note: exported from Xylem functional traits database Master 13 March 2015.xls but had to change a bunch of column names
xft <- read.csv("/Users/leeanderegg/Desktop/LFT analysis/XFT_traits_13_March_2015.csv", header=T,na.strings = "")
  # turns out there are 300+ uncleaned species
xft$Cleaned.family <- as.character(xft$Cleaned.family)
xft$Cleaned.family[which(is.na(xft$Cleaned.binomial))] <- xft$Family[which(is.na(xft$Cleaned.binomial))]
xft$Cleaned.genus <- as.character(xft$Cleaned.genus)
xft$Cleaned.genus[which(is.na(xft$Cleaned.binomial))] <- xft$Genus[which(is.na(xft$Cleaned.binomial))]
xft$Species <- as.character(xft$Species)
xft$Species[which(xft$Species==" limon")] <- "limon"
xft$Species.cleaning <- xft$Species
xft$Species.cleaning[grep(" x ", xft$Species)] <- "hybrid"
xft$Species.cleaning <- sstr_replace(xft$Species.cleaning, " \\(shade\\)", "")
xft$Species.cleaning[grep("�", xft$Species.cleaning)] <- "tremuloides"
xft$Species.cleaning <- str_replace(xft$Species.cleaning, " ","")
xft$Cleaned.binomial <- as.character(xft$Cleaned.binomial)
xft$Cleaned.binomial[which(is.na(xft$Cleaned.binomial))] <- paste(xft$Genus[which(is.na(xft$Cleaned.binomial))], xft$Species.cleaning[which(is.na(xft$Cleaned.binomial))])
# replacing 's' with 'S'
xft$Plant.organ[which(xft$Plant.organ=="s")] <- "S"


xft$Cleaned.family <- factor(xft$Cleaned.family)
xft$Cleaned.genus <- factor(xft$Cleaned.genus)
xft$Cleaned.binomial <- factor(xft$Cleaned.binomial)
xft$P50 <- str_replace(xft$P50, "_","-")
xft$P50 <- as.numeric(xft$P50)

## turns out there are 300+ rows that don't have 'Cleaned' family, genus, spp

# 
# 
# p50 <- lmer(log(-1*P50)~ Plant.organ + (1|Cleaned.family) + (1|Cleaned.genus) + (1|Cleaned.binomial), xft)
# Ks <- lmer(log(Ks)~ Plant.organ + (1|Cleaned.family) + (1|Cleaned.genus) + (1|Cleaned.binomial), xft)
# KL <- lmer(log(KL)~ Plant.organ + (1|Cleaned.family) + (1|Cleaned.genus) + (1|Cleaned.binomial), xft)
# 
# 





#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#### Combine all the variance decomps together
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 
# # original decomps for easy vis
# baadtraitvars <- data.frame(rawWDvariance2[,4], WDvariance2[,4], rawWDvariance[,4], WDvariance[,4], AbBlvariance[,4], Heightvariance[,4], LMFvariance[,4], RtSvariance[,4])
# colnames(baadtraitvars) <- c("GWDrawWD","GWDlogWD","BAADrawWD", "BAADlogWD", "logAGB_BGB", "logHeight", "LMF", "RtoS")
# rownames(baadtraitvars) <- c("BtwSpecies", "BtwGenera" , "BtwFamilies", "WtinSpp")
# 
# baadtraitvars_scaled <- baadtraitvars  
# for(i in 1:ncol(baadtraitvars)){
#   baadtraitvars_scaled[,i] <- baadtraitvars[,i]/sum(baadtraitvars[,i])
# }
# 
# 
# baadtraitvars_scaled2 <- baadtraitvars_scaled[c(3,2,1,4),]
# 
# baadtraitvars_scaled3 <- baadtraitvars_scaled2[,order(baadtraitvars_scaled2[1,], decreasing =T)]



########. NEW DECOMS for LFT paper #########

# the plan: includ - LMA, Nmass, LL, WD, RSR, P50, Ks

##### LMA, Nmass & LL ###########
data.all <- read.csv("/Users/leeanderegg/Dropbox/NACP_Traits/NACP_Traits_Rcode/FinalExample/DerivedData/AllTraitData_for_Anderegg_etal_2018_EcolLet.csv", header=T, row.names = 1)

## with all spp used for the hierarchical analysis
logLMAvar <- lmer(log.LMA~ 1 + (1|Project) + (1|Family) + (1|Genus) + (1|Species), data.all)
logLLvar <- lmer(log.LL~ 1 + (1|Project)  + (1|Family) + (1|Genus) + (1|Species), data.all[-which(data.all$Project=="CO"),])
logNmassvar <- lmer(log.Nmass~ 1 + (1|Project)  + (1|Family) + (1|Genus) + (1|Species), data.all)
logNareavar <- lmer(log.Narea~ 1 + (1|Project)  + (1|Family) + (1|Genus) + (1|Species), data.all)
# technically, Genus and Species random effects should be nested within Family, but the nested model never converged.
# tests on subsets yeild extremely similar variance estimates for nested and non-nested effects.

# create dataframes with variance parameter estimates
LMAvariance <- data.frame(VarCorr(logLMAvar))
LLvariance <- data.frame(VarCorr(logLLvar))
Nmassvariance <- data.frame(VarCorr(logNmassvar))
Nareavariance <- data.frame(VarCorr(logNareavar))



####### WD (duplicated from above for Anna's NP review) ######3
wdbaad <- baad %>% filter(r.st>0) %>% select(WD=r.st, Binomial=speciesMatched, Family=family, Genus=genus, Study=studyName) %>% mutate(WD=WD/1000)
wdGWD <- GWD %>% select(WD, Binomial, Family, Genus,Study=Reference.Number)
wdbaad$Study <- as.character(wdbaad$Study)
wdGWD$Study <- as.character(wdGWD$Study)
allwd <- rbind(wdbaad, wdGWD)

WDvar3 <- lmer(WD~ 1 + (1|Family) + (1|Genus) + (1|Binomial), allwd)
logWDvar3 <- lmer(log(WD, base=10)~ 1 + (1|Family) + (1|Genus) + (1|Binomial), allwd)

rawWDvariance3 <- data.frame(VarCorr(WDvar3))
#rawWDvariance3$scaledVar <- rawWDvariance3$vcov/ sum(rawWDvariance3$vcov)
#WDvariance3 <- data.frame(VarCorr(logWDvar3))



######## P50, Ks from XFT ###################

logp50var <- lmer(log(-1*P50)~ Plant.organ + (1|Cleaned.family) + (1|Cleaned.genus) + (1|Cleaned.binomial), xft)
logKsvar <- lmer(log(Ks)~ Plant.organ + (1|Cleaned.family) + (1|Cleaned.genus) + (1|Cleaned.binomial), xft)

# create dataframes with variance parameter estimates
P50variance <- data.frame(VarCorr(logp50var))
Ksvariance <- data.frame(VarCorr(logKsvar))



### R:S from baad (supplemented with Ledo 2017) ######

RSbaad <- baad[which(baad$RSR>0),] %>% select(latitude, longitude, vegetationTYPE=vegetation, Binomial=speciesMatched, Family=family, Genus=genus, RSR, log.RSR)
baadspp <- str_replace(unique(RSbaad$speciesMatched),pattern = " ",replacement = "_")
RtSunique <- RtS[which(!RtS$SPECIES %in% baadspp),] %>% select(latitude, longitude, vegetationTYPE, Binomial=SPECIES, Family=FAMILY, Genus=GENUS, RSR=RS_RATIO, log.RSR)
allRSR <- rbind(RSbaad, RtSunique)
allRSR$Binomial <- str_replace(allRSR$Binomial, "_"," ")

logRSRvar <- lmer(log.RSR ~ 1 + (1|Family) + (1|Genus) + (1|Binomial), allRSR)

logRSRvariance <- data.frame(VarCorr(logRSRvar))


# combine all variance estimates, leaving out "Project"
traitvars <- data.frame(LMAvariance[which(LMAvariance$grp !="Project"),4], LLvariance[which(LLvariance$grp !="Project"),4], Nmassvariance[which(Nmassvariance$grp !="Project"),4], Nareavariance[which(Nareavariance$grp !="Project"),4],
                        rawWDvariance3[,4], P50variance[,4],Ksvariance[,4], logRSRvariance[,4])
colnames(traitvars) <- c("logLMA", "logLL", "logNmass", "logNarea", "WD","logP50","logKs","logRSR")
rownames(traitvars) <- c("BtwSpecies", "BtwGenera", "BtwFamilies", "WtinSpecies")


traitvars_scaled <- traitvars  
for(i in 1:ncol(traitvars)){
  traitvars_scaled[,i] <- traitvars[,i]/sum(traitvars[,i])
}
# re-order rows
traitvars_scaled2 <- traitvars_scaled[c(3,2,1,4),]







###### Plotting them #####
#quartz(width=3.75, height=4)
quartz(width=6.81, height=3.5)
par(mgp=c(3,.7,0), cex.lab=1.3, cex.axis=1.1, mfrow=c(1,1), mar=c(6,2,3,6), oma=c(0,3,0,0))
cols <- rev(c(mypal[colchoices[c(1,2,3)]],"#CCCCCC"))#brewer.pal(11, "RdBu")[c(11,9,1, 6)]
# With old colors, alpha=CC, with new colors, alpha=99
bp <-barplot(as.matrix(traitvars_scaled2[,c("logNmass","logLL","logLMA","WD","logP50","logKs","logNarea","logRSR")]),beside=F,legend.text = F,xpd = T,las=2,args.legend = list(x=4, y=1.3, ncol=2), col = paste0(cols,"99")
             , names.arg = c(expression(paste(log(N[mass]))),"log(LeafLife)","log(LMA)", "WD",expression(paste(log(P[50]))),expression(paste(log(K[s]))), expression(paste(log(N[area]))), "log(R:S)")#log(Height)","AboveGrnd:\nBelowGrnd\nbiomass")
             , ylab= "Proportion of total Variance" #ylab="Proportion of total Variance\n(log traits)"
             , xlab="", mgp=c(2.5,.8,0))
#legend(xpd=T, x = 0, y=1.2, legend=rev(c("btw Fams","w/in Fam","w/in Gen","w/in Spp")), fill=rev(paste0(cols,"99")), ncol=2, bty="n",  cex=1.2)
#text(x = bp, y= par("usr")[3]-.05,labels =  c(expression(paste(log[10](LMA))), expression(paste(log[10](LL))),expression(paste(log[10](N[mass]))),expression(paste(log[10](N[area])))), srt=40, adj=1,xpd=NA, cex=1.4, font=1)
#mtext("a)", side=3, adj=-.1, line=1.3)
mtext("Proportion of Total Variance", side=2, line=2.8)
#mtext("BAAD Allometry data", side=3, line=0.3)
#quartz(width=5, height=4)
#cols <- brewer.pal(11, "RdBu")[c(11,9,1, 6)]
#legend(xpd=T, x = 11, y=0.7, legend=rev(c("btw Fams","btw Gen","btw Spp")), fill=rev(paste0(cols[1:3],"99")), ncol=1, bty="n",  cex=1.2)
# mtext(text="28 spp w/ replication, 1000+ spp,\n500+ genera, 150+ families",side = 1,line = 3.3)
# legend(xpd=T, x = 0, y=1.3, legend=rownames(rtraitvars_scaled2), fill=paste0(cols,"CC"), ncol=2, bty="n",  cex=1.2)
legend(xpd=NA, x = 10, y=0.7, legend=rev(c("btw Fams","btw Gen","btw Spp","w/in Spp")), fill=rev(paste0(cols,"99")), ncol=1, bty="n",  cex=1)






#### Plotting Subset for Presentation ########
## making everything the same size as the Ecol Let Vardecmp FIG 1
quartz(width=4.33, height=4.73)
par(mfrow=c(2,2), mgp=c(2,.7,0), cex.lab=1.1, cex.axis=1.1, mar=c(4.5,2,1.5,2),oma=c(0,2,3.8,0))

bp <-barplot(as.matrix(baadtraitvars_scaled3[,c("BAADrawWD","LMF","logHeight","logAGB_BGB"),]),beside=F,legend.text = F,xpd = T,las=2,args.legend = list(x=4, y=1.3, ncol=2), col = paste0(cols,"99")
             , names.arg = c("","","","")#c("Wood Density","Leaf Mass Fraction","log(Height)","AboveGrnd:BelowGrnd\nbiomass") #c("log(WD)","raw WD", "Leaf Mass\nFraction", "log(Height)","AboveGrnd:\nBelowGrnd\nbiomass")
             , ylab= "Proportion of total Variance" #ylab="Proportion of total Variance\n(log traits)"
             , xlab="", mgp=c(2.5,.8,0))
#legend(xpd=T, x = 0, y=1.2, legend=rev(c("btw Fams","w/in Fam","w/in Gen","w/in Spp")), fill=rev(paste0(cols,"99")), ncol=2, bty="n",  cex=1.2)
#text(x = bp, y= par("usr")[3]-.05,labels =  c(expression(paste(log[10](LMA))), expression(paste(log[10](LL))),expression(paste(log[10](N[mass]))),expression(paste(log[10](N[area])))), srt=40, adj=1,xpd=NA, cex=1.4, font=1)
#mtext("a)", side=3, adj=-.1, line=1.3)
mtext("Proportion of Total Variance", side=2, line=2.8)
text(x = bp, y= par("usr")[3]-.05,labels =  c("Wood Density", "Leaf Mass Fraction",expression(paste(log[10](Height))),expression(paste(M[AbvGr]:M[BlwGr]))), srt=40, adj=1,xpd=NA, cex=1.1, font=1)















########## attempt to add in higher classifications than Family ########


genera <- unique(c(as.character(data.all$Genus), as.character(baad$genus), as.character(GWD$Genus), as.character(RtS$GENUS), as.character(xft$Cleaned.genus)))

output <- data.frame(Order=rep(NA, length(genera)),Subclass=rep(NA, length(genera)),Class=rep(NA, length(genera)))

for (i in genera[1:4]){
    tmp <- taxize::classification(i, db="itis")[[1]]
    output[i,] <- cbind(tmp$name[which(tmp$rank=="order")], tmp$name[which(tmp$rank=="subclass" | tmp$rank=="superorder")], tmp$name[which(tmp$rank=="class")])
  }

