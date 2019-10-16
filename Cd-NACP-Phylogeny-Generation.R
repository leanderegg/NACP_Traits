###########################################3
###  Making a phylo tree for PNW trees    ###
############################################

require(stringr)
require(taxize)
require(brranching)
require(ape)
require(dplyr)

### step one: Make a Species List ###############3

biomass <- read.csv("/Users/leeanderegg/Dropbox/NACP_Traits/NACP_Traits_Rcode/FinalExample/DerivedData/PNW_Biomass_data_for_Anderegg_etal_2018_EcolLet.csv", header=T, row.names=1)
traits <- read.csv("/Users/leeanderegg/Dropbox/NACP_Traits/NACP_Traits_Rcode/FinalExample/DerivedData/PNW_Trait_data_for_Anderegg_etal_2018_EcolLet.csv", header=T, row.names=1)


spp.codes <- traits %>% group_by(GE.SP, Family) %>% summarise(Code=unique(SP.ID), fullspp = unique(FullSpecies), n_traits=n())
spplist_raw <- unique(traits$GE.SP)
spplist <- str_replace(spplist_raw[-grep("unknown", spplist_raw)], "\\.", " ")

spplist1 <- spplist[order(spplist)]

sppordered <- spplist1[c(10,33,34,7,8,19,30,31,32,13,14,9,12,15,16,29,18,28,20,21,22,27,26,24,25,23,37,36,4,2,1,6,3,5,11,35,17)]
LeafLife <- c("E","E","D","D","D","E","E","D","D","E","D","D","D","D","D","D","D",rep("E", times=20))

Spp.descriptions <- data.frame(sppordered, LeafLife)
Spp.descriptions$SP.ID <- spp.codes$Code[match(sppordered,spp.codes$fullspp)]
Spp.descriptions$n_traits <- spp.codes$n_traits[match(sppordered,spp.codes$fullspp)]
Spp.descriptions$Family <- spp.codes$Family[match(sppordered,spp.codes$fullspp)]

## clean the names with taxize:
#taxize::use_entrez()
#usethis::edit_r_environ()
spplist_phylo <- phylomatic_names(spplist)
spp_phylo <- phylomatic(spplist_phylo)


pinus1 <- read.tree("/Users/leeanderegg/Desktop/LFT analysis/subtree-node-mrcaott8444ott8454-Pinus--Cathaya.tre")

### summarize species statistics to plot on tips ####

biomass$bio1 <- with(biomass, SPP_O1_BASAL_AREA_FRACTION * AG_BIOMASS_TREE_TOTAL_AS_CARBON)
biomass$bio2 <- with(biomass, SPP_O2_BASAL_AREA_FRACTION * AG_BIOMASS_TREE_TOTAL_AS_CARBON)
biomass$bio3 <- with(biomass, SPP_O3_BASAL_AREA_FRACTION * AG_BIOMASS_TREE_TOTAL_AS_CARBON)
biomass$bio4 <- with(biomass, SPP_O4_BASAL_AREA_FRACTION * AG_BIOMASS_TREE_TOTAL_AS_CARBON)

sp1 <- data.frame(Species=biomass$SPP_O1_ABBREV, biomass=biomass$bio1, perc_BA = biomass$SPP_O1_BASAL_AREA_FRACTION)
sp2 <- data.frame(Species=biomass$SPP_O2_ABBREV, biomass=biomass$bio2, perc_BA = biomass$SPP_O2_BASAL_AREA_FRACTION)
sp3 <- data.frame(Species=biomass$SPP_O3_ABBREV, biomass=biomass$bio3, perc_BA = biomass$SPP_O3_BASAL_AREA_FRACTION)
sp4 <- data.frame(Species=biomass$SPP_O4_ABBREV, biomass=biomass$bio4, perc_BA = biomass$SPP_O4_BASAL_AREA_FRACTION)

biomasslong <- rbind(sp1[-which(is.na(sp1$biomass)|is.na(sp1$Species)),], sp2[-which(is.na(sp2$biomass) | sp2$Species=="none"),],sp3[-which(is.na(sp3$biomass) | sp2$Species=="none"),],sp4[-which(is.na(sp4$biomass) | sp2$Species=="none"),])
tot_bio <- biomasslong %>% group_by(Species) %>% summarise(n_plots=n(), total_biomass=sum(biomass, na.rm=T), tot_dominance = sum(perc_BA, na.rm=T), max_dominance = max(perc_BA, na.rm=T), mean_dominance=mean(perc_BA, na.rm=T), median_dominance=median(perc_BA, na.rm=T))
tot_bio_sorted <- tot_bio[match(Spp.descriptions$SP.ID,tot_bio$Species),]
  # note: this drops CERLED and PRUEMA, both of which don't have any trait measurements. but they're only tiny (<5%) fraction of 1&2 plots, respectively

nacp.sums <- left_join(Spp.descriptions, tot_bio, by=c("SP.ID"="Species"))


xtabs(~GE.SP, traits)
write.csv(nacp.sums, "SpeciesSummaries_forLFTanalysis_V1.csv")

# I added in some columns about dominance and also the first few LFTs
nacp.sums <- read.csv("SpeciesSummaries_forLFTanalysis_V2.csv")




#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
################# trying some variance decompositions with LFT concepts #############
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


traits$PFT <- nacp.sums$PFT[match(traits$SP.ID, nacp.sums$SP.ID)]
traits$Deep.LFT <- nacp.sums$Deep[match(traits$SP.ID, nacp.sums$SP.ID)]
traits$Mid.LFT <- nacp.sums$Mid[match(traits$SP.ID, nacp.sums$SP.ID)]
traits$Shallow.LFT <- nacp.sums$Shallow[match(traits$SP.ID, nacp.sums$SP.ID)]

# with Larix broken out
traits$PFT.larix <- nacp.sums$PFT.larix[match(traits$SP.ID, nacp.sums$SP.ID)]
traits$Deep.LFT.larix <- nacp.sums$Deep.larix[match(traits$SP.ID, nacp.sums$SP.ID)]
traits$Mid.LFT.larix <- nacp.sums$Mid.larix[match(traits$SP.ID, nacp.sums$SP.ID)]
traits$Shallow.LFT.larix <- nacp.sums$Shallow.larix[match(traits$SP.ID, nacp.sums$SP.ID)]



######### Univariate explained variance ##############
Nmass.PFTmod <- lm(log.Nmass~PFT, traits)
Nmass.Deep.LFTmod <- lm(log.Nmass~Deep.LFT, traits)
Nmass.Mid.LFTmod <- lm(log.Nmass~Mid.LFT, traits)
Nmass.Shallow.LFTmod <- lm(log.Nmass~Shallow.LFT, traits)

LMA.PFTmod <- lm(log.LMA~PFT, traits)
LMA.Deep.LFTmod <- lm(log.LMA~Deep.LFT, traits)
LMA.Mid.LFTmod <- lm(log.LMA~Mid.LFT, traits)
LMA.Shallow.LFTmod <- lm(log.LMA~Shallow.LFT, traits)
LMA.clust1mod <- lm(log.LMA~clust1, traits)


LL.PFTmod <- lm(log.LL~PFT, traits)
LL.Deep.LFTmod <- lm(log.LL~Deep.LFT, traits)
LL.Mid.LFTmod <- lm(log.LL~Mid.LFT, traits)
LL.Shallow.LFTmod <- lm(log.LL~Shallow.LFT, traits)
LL.clust1mod <- lm(log.LL~clust1, traits)



## a cluster analysis, unsupervised

clust1 <- hclust(d=dist(traits[which(complete.cases(traits[,c("log.LMA","log.Nmass","log.LL")])==T),c("log.LMA","log.Nmass","log.LL")]),method = "ward.D", members = NULL)
clust1.group8 <-  rect.hclust(clust1, k = 8)
clust1.groups <- rep(NA, nrow(traits[which(complete.cases(traits[,c("log.LMA","log.Nmass","log.LL")])==T),]))
for(i in 1:length(clust1.groups)){
  for(j in 1:8){
    if(i %in% clust1.group8[[j]]) clust1.groups[i] <- j
  }
  
}
traits$clust1 <- NA
traits$clust1[which(complete.cases(traits[,c("log.LMA","log.Nmass","log.LL")])==T)] <- clust1.groups


clust.kmeans <- kmeans(x=traits[which(complete.cases(traits[,c("log.LMA","log.Nmass","log.LL")])==T),c("log.LMA","log.Nmass","log.LL")],centers = 8,nstart = 100)
traits$kmeans <- NA
traits$kmeans[which(complete.cases(traits[,c("log.LMA","log.Nmass","log.LL")])==T)] <- clust.kmeans$cluster


####### . Analysis of kmeans consistency with different n ########
clust.kmeans.5 <- kmeans(x=traits[which(complete.cases(traits[,c("log.LMA","log.Nmass","log.LL")])==T),c("log.LMA","log.Nmass","log.LL")],centers = 5,nstart = 100)
traits$kmeans.5 <- NA
traits$kmeans.5[which(complete.cases(traits[,c("log.LMA","log.Nmass","log.LL")])==T)] <- clust.kmeans.5$cluster

clust.kmeans.6 <- kmeans(x=traits[which(complete.cases(traits[,c("log.LMA","log.Nmass","log.LL")])==T),c("log.LMA","log.Nmass","log.LL")],centers = 6,nstart = 100)
traits$kmeans.6 <- NA
traits$kmeans.6[which(complete.cases(traits[,c("log.LMA","log.Nmass","log.LL")])==T)] <- clust.kmeans.6$cluster

clust.kmeans.10 <- kmeans(x=traits[which(complete.cases(traits[,c("log.LMA","log.Nmass","log.LL")])==T),c("log.LMA","log.Nmass","log.LL")],centers = 10,nstart = 100)
traits$kmeans.10 <- NA
traits$kmeans.10[which(complete.cases(traits[,c("log.LMA","log.Nmass","log.LL")])==T)] <- clust.kmeans.10$cluster



######### Multivariate variance explained ###############
data.clean <- traits[which(complete.cases(traits[,c("log.LMA","log.Nmass","log.LL","PFT")])==T ),c("log.LMA","log.Nmass","log.LL")]
predictors.clean <- traits[which(complete.cases(traits[,c("log.LMA","log.Nmass","log.LL","PFT")])==T ),c("SP.ID","PFT","Deep.LFT","Mid.LFT","Shallow.LFT","clust1","kmeans","kmeans.5","kmeans.6","kmeans.10", "PFT.larix","Deep.LFT.larix","Mid.LFT.larix","Shallow.LFT.larix")]
# & traits$SP.ID!="LAROCC" & traits$SP.ID!="PURTRI"

nacp.sums.abund <- nacp.sums[which(nacp.sums$n_plots>0),]
data.abund <- traits[which(complete.cases(traits[,c("log.LMA","log.Nmass","log.LL","PFT")])==T & traits$SP.ID %in% nacp.sums.abund$SP.ID),c("log.LMA","log.Nmass","log.LL")]
predictors.abund <- traits[which(complete.cases(traits[,c("log.LMA","log.Nmass","log.LL","PFT")])==T & traits$SP.ID %in% nacp.sums.abund$SP.ID),c("SP.ID","PFT","Deep.LFT","Mid.LFT","Shallow.LFT","clust1","kmeans","kmeans.5", "PFT.larix","Deep.LFT.larix","Mid.LFT.larix","Shallow.LFT.larix")]

# removing Larix all together
data.clean.nl <- traits[which(complete.cases(traits[,c("log.LMA","log.Nmass","log.LL","PFT")])==T & traits$SP.ID!="LAROCC"),c("log.LMA","log.Nmass","log.LL")]
predictors.clean.nl <- traits[which(complete.cases(traits[,c("log.LMA","log.Nmass","log.LL","PFT")])==T & traits$SP.ID!="LAROCC"),c("SP.ID","PFT","Deep.LFT","Mid.LFT","Shallow.LFT","clust1","kmeans", "PFT.larix","Deep.LFT.larix","Mid.LFT.larix","Shallow.LFT.larix")]
data.abund.nl <- traits[which(complete.cases(traits[,c("log.LMA","log.Nmass","log.LL","PFT")])==T & traits$SP.ID %in% nacp.sums.abund$SP.ID & traits$SP.ID!="LAROCC"),c("log.LMA","log.Nmass","log.LL")]
predictors.abund.nl <- traits[which(complete.cases(traits[,c("log.LMA","log.Nmass","log.LL","PFT")])==T & traits$SP.ID %in% nacp.sums.abund$SP.ID& traits$SP.ID!="LAROCC"),c("SP.ID","PFT","Deep.LFT","Mid.LFT","Shallow.LFT","clust1","kmeans", "PFT.larix","Deep.LFT.larix","Mid.LFT.larix","Shallow.LFT.larix")]


### maybe I can summarize in 3d traits space?
require(vegan)


##### .everything, no splitting Larix ####
PFTmod <- adonis(data.clean~factor(predictors.clean$PFT), method="euc")
Deep.LFTmod <- adonis(data.clean~factor(predictors.clean$Deep.LFT), method="euc")
Mid.LFTmod <- adonis(data.clean~factor(predictors.clean$Mid.LFT), method="euc")
Shallow.LFTmod <- adonis(data.clean~factor(predictors.clean$Shallow.LFT), method="euc")
clust1mod <- adonis(data.clean~factor(predictors.clean$clust1), method="euc")
kmeansmod <- adonis(data.clean~factor(predictors.clean$kmeans), method="euc")
kmeans5mod <- adonis(data.clean~factor(predictors.clean$kmeans.5), method="euc")
kmeans6mod <- adonis(data.clean~factor(predictors.clean$kmeans.6), method="euc")
kmeans10mod <- adonis(data.clean~factor(predictors.clean$kmeans.10), method="euc")

PFTmod # R2 = 47.8                      # 51.3 no larix # 35.6 no larix or PURTRI # 26.0 abund only # 34.9 abund w/ larix
Deep.LFTmod # R2 = 45.0                 # 48.3 no larix # 30.1 no larix or PURTRI # 21.1 abund only # 30.1
Mid.LFTmod  # R2 = 49.1  #54.2 w/ larix # 52.1 no larix # 35.6 no larix or PURTRI # 27.5 abund only # 35.5
Shallow.LFTmod #R2= 57.6 #62.2 w/ larix # 60.5 no larix # 47.5 no larix or PURTRI # 47.9 abund only # 47.9
clust1mod # R2 = 83.7                   # 83 no larix   # no larix or PURTRI
kmeans5mod
kmeans6mod
kmeansmod # R2 = 85.3                   # 84.7 no larix # no larix or PURTRI
kmeans10mod

### In MS used PFT, Deep, Mid and Shallow.larix


##### .everything, splitting Larix ####
PFTmod.larix <- adonis(data.clean~factor(predictors.clean$PFT.larix), method="euc")
Deep.LFTmod.larix <- adonis(data.clean~factor(predictors.clean$Deep.LFT.larix), method="euc")
Mid.LFTmod.larix <- adonis(data.clean~factor(predictors.clean$Mid.LFT.larix), method="euc")
Shallow.LFTmod.larix <- adonis(data.clean~factor(predictors.clean$Shallow.LFT.larix), method="euc")
#clust1mod <- adonis(data.clean~factor(predictors.clean$clust1), method="euc")
#kmeansmod <- adonis(data.clean~factor(predictors.clean$kmeans), method="euc")


PFTmod.larix
Deep.LFTmod.larix 
Mid.LFTmod.larix  
Shallow.LFTmod.larix 



##### .Abund, no splitting Larix ####
PFTmod.abund <- adonis(data.abund~factor(predictors.abund$PFT), method="euc")
Deep.LFTmod.abund <- adonis(data.abund~factor(predictors.abund$Deep.LFT), method="euc")
Mid.LFTmod.abund <- adonis(data.abund~factor(predictors.abund$Mid.LFT), method="euc")
Shallow.LFTmod.abund <- adonis(data.abund~factor(predictors.abund$Shallow.LFT), method="euc")
clust1mod.abund <- adonis(data.abund~factor(predictors.abund$clust1), method="euc")
kmeansmod.abund <- adonis(data.abund~factor(predictors.abund$kmeans), method="euc")


PFTmod.abund #
Deep.LFTmod.abund # 
Mid.LFTmod.abund  # 
Shallow.LFTmod.abund
clust1mod.abund 
kmeansmod.abund 


##### .Abund, splitting Larix ####
PFTmod.abund.larix <- adonis(data.abund~factor(predictors.abund$PFT.larix), method="euc")
Deep.LFTmod.abund.larix <- adonis(data.abund~factor(predictors.abund$Deep.LFT.larix), method="euc")
Mid.LFTmod.abund.larix <- adonis(data.abund~factor(predictors.abund$Mid.LFT.larix), method="euc")
Shallow.LFTmod.abund.larix <- adonis(data.abund~factor(predictors.abund$Shallow.LFT.larix), method="euc")
#clust1mod <- adonis(data.abund~factor(predictors.abund$clust1), method="euc")
#kmeansmod <- adonis(data.abund~factor(predictors.abund$kmeans), method="euc")


PFTmod.abund.larix
Deep.LFTmod.abund.larix 
Mid.LFTmod.abund.larix  
Shallow.LFTmod.abund.larix 




##### .everything, no Larix ####
PFTmod.nl <- adonis(data.clean.nl~factor(predictors.clean.nl$PFT), method="euc")
Deep.LFTmod.nl <- adonis(data.clean.nl~factor(predictors.clean.nl$Deep.LFT), method="euc")
Mid.LFTmod.nl <- adonis(data.clean.nl~factor(predictors.clean.nl$Mid.LFT), method="euc")
Shallow.LFTmod.nl <- adonis(data.clean.nl~factor(predictors.clean.nl$Shallow.LFT), method="euc")
clust1mod.nl <- adonis(data.clean.nl~factor(predictors.clean.nl$clust1), method="euc")
kmeansmod.nl <- adonis(data.clean.nl~factor(predictors.clean.nl$kmeans), method="euc")


PFTmod.nl # R2 = 47.8                      # 51.3 no larix # 35.6 no larix or PURTRI # 26.0 abund only # 34.9 abund w/ larix
Deep.LFTmod.nl # R2 = 45.0                 # 48.3 no larix # 30.1 no larix or PURTRI # 21.1 abund only # 30.1
Mid.LFTmod.nl  # R2 = 49.1  #54.2 w/ larix # 52.1 no larix # 35.6 no larix or PURTRI # 27.5 abund only # 35.5
Shallow.LFTmod.nl #R2= 57.6 #62.2 w/ larix # 60.5 no larix # 47.5 no larix or PURTRI # 47.9 abund only # 47.9
clust1mod.nl # R2 = 83.7                   # 83 no larix   # no larix or PURTRI
kmeansmod.nl # R2 = 85.3                   # 84.7 no larix # no larix or PURTRI




##### .Abund, no Larix ####
PFTmod.abund.nl <- adonis(data.abund.nl~factor(predictors.abund.nl$PFT), method="euc")
Deep.LFTmod.abund.nl <- adonis(data.abund.nl~factor(predictors.abund.nl$Deep.LFT), method="euc")
Mid.LFTmod.abund.nl <- adonis(data.abund.nl~factor(predictors.abund.nl$Mid.LFT), method="euc")
Shallow.LFTmod.abund.nl <- adonis(data.abund.nl~factor(predictors.abund.nl$Shallow.LFT), method="euc")
#clust1mod <- adonis(data.abund.nl~factor(predictors.abund$clust1), method="euc")
#kmeansmod <- adonis(data.abund.nl~factor(predictors.abund$kmeans), method="euc")


PFTmod.abund.nl
Deep.LFTmod.abund.nl 
Mid.LFTmod.abund.nl  
Shallow.LFTmod.abund.nl 




###### boy was this wrong....
# WCSS <- function(dataz, classification){
#   class <- factor(classification)
#   classes <- levels(class)
#   nclasses <- length(classes)
#   wcss <- rep(NA, times=nclasses)
#   for(i in 1:nclasses){
#     u <- rowMeans(dataz[which(classification==classes[i] & complete.cases(dataz)==T),])
#     wcss[i] <- sum(dist(dataz[which(classification==classes[i] & complete.cases(dataz)==T),])^2)
#   }
#   return(wcss)
# }
# 
# 
# 
# PFTss <- WCSS(dataz = traits[,c("log.Nmass","log.LMA","log.LL")], classification = traits$PFT)
# Deep.LFTss <- WCSS(dataz = traits[,c("log.Nmass","log.LMA","log.LL")], classification = traits$Deep.LFT)
# Mid.LFTss <- WCSS(dataz = traits[,c("log.Nmass","log.LMA","log.LL")], classification = traits$Mid.LFT)
# Shallow.LFTss <- WCSS(dataz = traits[,c("log.Nmass","log.LMA","log.LL")], classification = traits$Shallow.LFT)
# clust1.ss <- WCSS(dataz = traits[,c("log.Nmass","log.LMA","log.LL")], classification = traits$clust1)
# kmeans.ss <- WCSS(dataz = traits[,c("log.Nmass","log.LMA","log.LL")], classification = traits$kmeans)
# 
# 
# 1-sum(PFTss)/sum(dist(dataz[which(complete.cases(dataz)==T),])^2)
# 1-sum(Deep.LFTss)/sum(dist(dataz[which(complete.cases(dataz)==T),])^2)
# 1-sum(Mid.LFTss)/sum(dist(dataz[which(complete.cases(dataz)==T),])^2)
# 1-sum(Shallow.LFTss)/sum(dist(dataz[which(complete.cases(dataz)==T),])^2)
# 1-sum(clust1.ss)/sum(dist(dataz[which(complete.cases(dataz)==T),])^2)
# 1-sum(kmeans.ss)/sum(dist(dataz[which(complete.cases(dataz)==T),])^2)

test <- adonis(dataz~clust.kmeans$cluster, method="euc")

### plotting boxplot of biomass to show hyperdominance
quartz(width=2.5, height=5)
par(mar=c(3.5,1,1,1), mgp=c(2.5,1,0))
barplot(height = rev(nacp.sums$total_biomass/1000), names.arg = NA, horiz = T, las=1, xlab="total biomass (kg C)")






######## Plotting example with LMA and Nmass #########3
mypal <- paste0(brewer.pal(n=5, "Set1")[c(5,2,4,3,1)],"66")
mypal <- brewer.pal(n=5, "Set1")[c(5,2,4,3,1)]

palette(mypal)


quartz(width=2.5, height=2.5)
par(mar=c(3,3,1,1), mgp=c(2,.5,0))
plot(log.LMA~log.Nmass, traits, col=Mid.LFT, cex=.8, pch=16)
points(log.LMA~log.Nmass, traits[which(traits$Mid.LFT=="Asterid"),], col=Mid.LFT, cex=.8, pch=16)

quartz(width=2.5, height=5)
par(mar=c(3,3,1,1), mgp=c(2,.5,0), mfrow=c(2,1))
plot(log.LMA~log.Nmass, traits, col=Mid.LFT, cex=.8, ylab=expression(paste(log[10](LMA))), xlab=expression(paste(log[10](N[mass]))))
points(log.LMA~log.Nmass, traits[which(traits$Mid.LFT=="Asterid"),], col=Mid.LFT, cex=.8)

plot(log.LMA~log.Nmass, traits, col=kmeans.5, cex=.8, ylab=expression(paste(log[10](LMA))), xlab=expression(paste(log[10](N[mass]))))






#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
########### Variance of CWMs explained by LFTs ############################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# for .v1 and .v2 used PFT, Deep, Mid and Shallow.larix

biomass$PFT1 <- nacp.sums$PFT[match(biomass$SPP_O1_ABBREV, nacp.sums$SP.ID)]
biomass$Deep.LFT1 <- nacp.sums$Deep[match(biomass$SPP_O1_ABBREV, nacp.sums$SP.ID)]
biomass$Mid.LFT1 <- nacp.sums$Mid[match(biomass$SPP_O1_ABBREV, nacp.sums$SP.ID)]
biomass$Shallow.LFT1 <- nacp.sums$Shallow.larix[match(biomass$SPP_O1_ABBREV, nacp.sums$SP.ID)]

biomass$PFT2 <- nacp.sums$PFT[match(biomass$SPP_O2_ABBREV, nacp.sums$SP.ID)]
biomass$Deep.LFT2 <- nacp.sums$Deep[match(biomass$SPP_O2_ABBREV, nacp.sums$SP.ID)]
biomass$Mid.LFT2 <- nacp.sums$Mid[match(biomass$SPP_O2_ABBREV, nacp.sums$SP.ID)]
biomass$Shallow.LFT2 <- nacp.sums$Shallow.larix[match(biomass$SPP_O2_ABBREV, nacp.sums$SP.ID)]


biomass$PFT3 <- nacp.sums$PFT[match(biomass$SPP_O3_ABBREV, nacp.sums$SP.ID)]
biomass$Deep.LFT3 <- nacp.sums$Deep[match(biomass$SPP_O3_ABBREV, nacp.sums$SP.ID)]
biomass$Mid.LFT3 <- nacp.sums$Mid[match(biomass$SPP_O3_ABBREV, nacp.sums$SP.ID)]
biomass$Shallow.LFT3 <- nacp.sums$Shallow.larix[match(biomass$SPP_O3_ABBREV, nacp.sums$SP.ID)]


biomass$PFT4 <- nacp.sums$PFT[match(biomass$SPP_O4_ABBREV, nacp.sums$SP.ID)]
biomass$Deep.LFT4 <- nacp.sums$Deep[match(biomass$SPP_O4_ABBREV, nacp.sums$SP.ID)]
biomass$Mid.LFT4 <- nacp.sums$Mid[match(biomass$SPP_O4_ABBREV, nacp.sums$SP.ID)]
biomass$Shallow.LFT4 <- nacp.sums$Shallow.larix[match(biomass$SPP_O4_ABBREV, nacp.sums$SP.ID)]


PFT.traitmeans <- traits %>% group_by(PFT) %>% summarise(LL=mean(LLmonths, na.rm=T), LMA = mean(LMA, na.rm=T), Narea=mean(Narea, na.rm=T), Nmass=mean(LEAF_NITROGEN, na.rm=T)) %>% mutate(log.LL = log(LL, base=10), log.LMA=log(LMA, base=10), log.Narea = log(Narea, base=10), log.Nmass=log(Nmass, base=10))
Deep.LFT.traitmeans <- traits %>% group_by(Deep.LFT) %>% summarise(LL=mean(LLmonths, na.rm=T), LMA = mean(LMA, na.rm=T), Narea=mean(Narea, na.rm=T), Nmass=mean(LEAF_NITROGEN, na.rm=T)) %>% mutate(log.LL = log(LL, base=10), log.LMA=log(LMA, base=10), log.Narea = log(Narea, base=10), log.Nmass=log(Nmass, base=10))
Mid.LFT.traitmeans <- traits %>% group_by(Mid.LFT) %>% summarise(LL=mean(LLmonths, na.rm=T), LMA = mean(LMA, na.rm=T), Narea=mean(Narea, na.rm=T), Nmass=mean(LEAF_NITROGEN, na.rm=T)) %>% mutate(log.LL = log(LL, base=10), log.LMA=log(LMA, base=10), log.Narea = log(Narea, base=10), log.Nmass=log(Nmass, base=10))
Shallow.LFT.traitmeans <- traits %>% group_by(Shallow.LFT.larix) %>% summarise(LL=mean(LLmonths, na.rm=T), LMA = mean(LMA, na.rm=T), Narea=mean(Narea, na.rm=T), Nmass=mean(LEAF_NITROGEN, na.rm=T)) %>% mutate(log.LL = log(LL, base=10), log.LMA=log(LMA, base=10), log.Narea = log(Narea, base=10), log.Nmass=log(Nmass, base=10))




#######. LMA CWM by disaggregation method #################3
# calcualte community weighted LMA by PFT
PFT.LMA1 <- PFT.traitmeans$LMA[match(as.character(biomass$PFT1), as.character(PFT.traitmeans$PFT))] * biomass$SPP_O1_BASAL_AREA_FRACTION/100
PFT.LMA2 <- PFT.traitmeans$LMA[match(biomass$PFT2, PFT.traitmeans$PFT)] * biomass$SPP_O2_BASAL_AREA_FRACTION/100
PFT.LMA3 <- PFT.traitmeans$LMA[match(biomass$PFT3, PFT.traitmeans$PFT)] * biomass$SPP_O3_BASAL_AREA_FRACTION/100
PFT.LMA4 <- PFT.traitmeans$LMA[match(biomass$PFT4, PFT.traitmeans$PFT)] * biomass$SPP_O4_BASAL_AREA_FRACTION/100

biomass$PFT.LMA <- apply(data.frame(PFT.LMA1,PFT.LMA2,PFT.LMA3,PFT.LMA4),MARGIN=1,FUN=sum, na.rm=T)

# calcualte community weighted LMA by Deep LFT
Deep.LFT.LMA1 <- Deep.LFT.traitmeans$LMA[match(as.character(biomass$Deep.LFT1), as.character(Deep.LFT.traitmeans$Deep.LFT))] * biomass$SPP_O1_BASAL_AREA_FRACTION/100
Deep.LFT.LMA2 <- Deep.LFT.traitmeans$LMA[match(biomass$Deep.LFT2, Deep.LFT.traitmeans$Deep.LFT)] * biomass$SPP_O2_BASAL_AREA_FRACTION/100
Deep.LFT.LMA3 <- Deep.LFT.traitmeans$LMA[match(biomass$Deep.LFT3, Deep.LFT.traitmeans$Deep.LFT)] * biomass$SPP_O3_BASAL_AREA_FRACTION/100
Deep.LFT.LMA4 <- Deep.LFT.traitmeans$LMA[match(biomass$Deep.LFT4, Deep.LFT.traitmeans$Deep.LFT)] * biomass$SPP_O4_BASAL_AREA_FRACTION/100

biomass$Deep.LFT.LMA <- apply(data.frame(Deep.LFT.LMA1,Deep.LFT.LMA2,Deep.LFT.LMA3,Deep.LFT.LMA4),MARGIN=1,FUN=sum, na.rm=T)


# calcualte community weighted LMA by Mid LFT
Mid.LFT.LMA1 <- Mid.LFT.traitmeans$LMA[match(as.character(biomass$Mid.LFT1), as.character(Mid.LFT.traitmeans$Mid.LFT))] * biomass$SPP_O1_BASAL_AREA_FRACTION/100
Mid.LFT.LMA2 <- Mid.LFT.traitmeans$LMA[match(biomass$Mid.LFT2, Mid.LFT.traitmeans$Mid.LFT)] * biomass$SPP_O2_BASAL_AREA_FRACTION/100
Mid.LFT.LMA3 <- Mid.LFT.traitmeans$LMA[match(biomass$Mid.LFT3, Mid.LFT.traitmeans$Mid.LFT)] * biomass$SPP_O3_BASAL_AREA_FRACTION/100
Mid.LFT.LMA4 <- Mid.LFT.traitmeans$LMA[match(biomass$Mid.LFT4, Mid.LFT.traitmeans$Mid.LFT)] * biomass$SPP_O4_BASAL_AREA_FRACTION/100

biomass$Mid.LFT.LMA <- apply(data.frame(Mid.LFT.LMA1,Mid.LFT.LMA2,Mid.LFT.LMA3,Mid.LFT.LMA4),MARGIN=1,FUN=sum, na.rm=T)

# calcualte community weighted LMA by Shallow LFT
Shallow.LFT.LMA1 <- Shallow.LFT.traitmeans$LMA[match(as.character(biomass$Shallow.LFT1), as.character(Shallow.LFT.traitmeans$Shallow.LFT.larix))] * biomass$SPP_O1_BASAL_AREA_FRACTION/100
Shallow.LFT.LMA2 <- Shallow.LFT.traitmeans$LMA[match(biomass$Shallow.LFT2, Shallow.LFT.traitmeans$Shallow.LFT.larix)] * biomass$SPP_O2_BASAL_AREA_FRACTION/100
Shallow.LFT.LMA3 <- Shallow.LFT.traitmeans$LMA[match(biomass$Shallow.LFT3, Shallow.LFT.traitmeans$Shallow.LFT.larix)] * biomass$SPP_O3_BASAL_AREA_FRACTION/100
Shallow.LFT.LMA4 <- Shallow.LFT.traitmeans$LMA[match(biomass$Shallow.LFT4, Shallow.LFT.traitmeans$Shallow.LFT.larix)] * biomass$SPP_O4_BASAL_AREA_FRACTION/100

biomass$Shallow.LFT.LMA <- apply(data.frame(Shallow.LFT.LMA1,Shallow.LFT.LMA2,Shallow.LFT.LMA3,Shallow.LFT.LMA4),MARGIN=1,FUN=sum, na.rm=T)



#######. LL CWM by disaggregation method #################3
# calcualte community weighted LL by PFT
PFT.LL1 <- PFT.traitmeans$LL[match(as.character(biomass$PFT1), as.character(PFT.traitmeans$PFT))] * biomass$SPP_O1_BASAL_AREA_FRACTION/100
PFT.LL2 <- PFT.traitmeans$LL[match(biomass$PFT2, PFT.traitmeans$PFT)] * biomass$SPP_O2_BASAL_AREA_FRACTION/100
PFT.LL3 <- PFT.traitmeans$LL[match(biomass$PFT3, PFT.traitmeans$PFT)] * biomass$SPP_O3_BASAL_AREA_FRACTION/100
PFT.LL4 <- PFT.traitmeans$LL[match(biomass$PFT4, PFT.traitmeans$PFT)] * biomass$SPP_O4_BASAL_AREA_FRACTION/100

biomass$PFT.LL <- apply(data.frame(PFT.LL1,PFT.LL2,PFT.LL3,PFT.LL4),MARGIN=1,FUN=sum, na.rm=T)

# calcualte community weighted LL by Deep LFT
Deep.LFT.LL1 <- Deep.LFT.traitmeans$LL[match(as.character(biomass$Deep.LFT1), as.character(Deep.LFT.traitmeans$Deep.LFT))] * biomass$SPP_O1_BASAL_AREA_FRACTION/100
Deep.LFT.LL2 <- Deep.LFT.traitmeans$LL[match(biomass$Deep.LFT2, Deep.LFT.traitmeans$Deep.LFT)] * biomass$SPP_O2_BASAL_AREA_FRACTION/100
Deep.LFT.LL3 <- Deep.LFT.traitmeans$LL[match(biomass$Deep.LFT3, Deep.LFT.traitmeans$Deep.LFT)] * biomass$SPP_O3_BASAL_AREA_FRACTION/100
Deep.LFT.LL4 <- Deep.LFT.traitmeans$LL[match(biomass$Deep.LFT4, Deep.LFT.traitmeans$Deep.LFT)] * biomass$SPP_O4_BASAL_AREA_FRACTION/100

biomass$Deep.LFT.LL <- apply(data.frame(Deep.LFT.LL1,Deep.LFT.LL2,Deep.LFT.LL3,Deep.LFT.LL4),MARGIN=1,FUN=sum, na.rm=T)


# calcualte community weighted LL by Mid LFT
Mid.LFT.LL1 <- Mid.LFT.traitmeans$LL[match(as.character(biomass$Mid.LFT1), as.character(Mid.LFT.traitmeans$Mid.LFT))] * biomass$SPP_O1_BASAL_AREA_FRACTION/100
Mid.LFT.LL2 <- Mid.LFT.traitmeans$LL[match(biomass$Mid.LFT2, Mid.LFT.traitmeans$Mid.LFT)] * biomass$SPP_O2_BASAL_AREA_FRACTION/100
Mid.LFT.LL3 <- Mid.LFT.traitmeans$LL[match(biomass$Mid.LFT3, Mid.LFT.traitmeans$Mid.LFT)] * biomass$SPP_O3_BASAL_AREA_FRACTION/100
Mid.LFT.LL4 <- Mid.LFT.traitmeans$LL[match(biomass$Mid.LFT4, Mid.LFT.traitmeans$Mid.LFT)] * biomass$SPP_O4_BASAL_AREA_FRACTION/100

biomass$Mid.LFT.LL <- apply(data.frame(Mid.LFT.LL1,Mid.LFT.LL2,Mid.LFT.LL3,Mid.LFT.LL4),MARGIN=1,FUN=sum, na.rm=T)


# calcualte community weighted LL by Shallow LFT
Shallow.LFT.LL1 <- Shallow.LFT.traitmeans$LL[match(as.character(biomass$Shallow.LFT1), as.character(Shallow.LFT.traitmeans$Shallow.LFT.larix))] * biomass$SPP_O1_BASAL_AREA_FRACTION/100
Shallow.LFT.LL2 <- Shallow.LFT.traitmeans$LL[match(biomass$Shallow.LFT2, Shallow.LFT.traitmeans$Shallow.LFT.larix)] * biomass$SPP_O2_BASAL_AREA_FRACTION/100
Shallow.LFT.LL3 <- Shallow.LFT.traitmeans$LL[match(biomass$Shallow.LFT3, Shallow.LFT.traitmeans$Shallow.LFT.larix)] * biomass$SPP_O3_BASAL_AREA_FRACTION/100
Shallow.LFT.LL4 <- Shallow.LFT.traitmeans$LL[match(biomass$Shallow.LFT4, Shallow.LFT.traitmeans$Shallow.LFT.larix)] * biomass$SPP_O4_BASAL_AREA_FRACTION/100

biomass$Shallow.LFT.LL <- apply(data.frame(Shallow.LFT.LL1,Shallow.LFT.LL2,Shallow.LFT.LL3,Shallow.LFT.LL4),MARGIN=1,FUN=sum, na.rm=T)







#######. Nmass CWM by disaggregation method #################3
# calcualte community weighted Nmass by PFT
PFT.Nmass1 <- PFT.traitmeans$Nmass[match(as.character(biomass$PFT1), as.character(PFT.traitmeans$PFT))] * biomass$SPP_O1_BASAL_AREA_FRACTION/100
PFT.Nmass2 <- PFT.traitmeans$Nmass[match(biomass$PFT2, PFT.traitmeans$PFT)] * biomass$SPP_O2_BASAL_AREA_FRACTION/100
PFT.Nmass3 <- PFT.traitmeans$Nmass[match(biomass$PFT3, PFT.traitmeans$PFT)] * biomass$SPP_O3_BASAL_AREA_FRACTION/100
PFT.Nmass4 <- PFT.traitmeans$Nmass[match(biomass$PFT4, PFT.traitmeans$PFT)] * biomass$SPP_O4_BASAL_AREA_FRACTION/100

biomass$PFT.Nmass <- apply(data.frame(PFT.Nmass1,PFT.Nmass2,PFT.Nmass3,PFT.Nmass4),MARGIN=1,FUN=sum, na.rm=T)

# calcualte community weighted Nmass by Deep LFT
Deep.LFT.Nmass1 <- Deep.LFT.traitmeans$Nmass[match(as.character(biomass$Deep.LFT1), as.character(Deep.LFT.traitmeans$Deep.LFT))] * biomass$SPP_O1_BASAL_AREA_FRACTION/100
Deep.LFT.Nmass2 <- Deep.LFT.traitmeans$Nmass[match(biomass$Deep.LFT2, Deep.LFT.traitmeans$Deep.LFT)] * biomass$SPP_O2_BASAL_AREA_FRACTION/100
Deep.LFT.Nmass3 <- Deep.LFT.traitmeans$Nmass[match(biomass$Deep.LFT3, Deep.LFT.traitmeans$Deep.LFT)] * biomass$SPP_O3_BASAL_AREA_FRACTION/100
Deep.LFT.Nmass4 <- Deep.LFT.traitmeans$Nmass[match(biomass$Deep.LFT4, Deep.LFT.traitmeans$Deep.LFT)] * biomass$SPP_O4_BASAL_AREA_FRACTION/100

biomass$Deep.LFT.Nmass <- apply(data.frame(Deep.LFT.Nmass1,Deep.LFT.Nmass2,Deep.LFT.Nmass3,Deep.LFT.Nmass4),MARGIN=1,FUN=sum, na.rm=T)


# calcualte community weighted Nmass by Mid LFT
Mid.LFT.Nmass1 <- Mid.LFT.traitmeans$Nmass[match(as.character(biomass$Mid.LFT1), as.character(Mid.LFT.traitmeans$Mid.LFT))] * biomass$SPP_O1_BASAL_AREA_FRACTION/100
Mid.LFT.Nmass2 <- Mid.LFT.traitmeans$Nmass[match(biomass$Mid.LFT2, Mid.LFT.traitmeans$Mid.LFT)] * biomass$SPP_O2_BASAL_AREA_FRACTION/100
Mid.LFT.Nmass3 <- Mid.LFT.traitmeans$Nmass[match(biomass$Mid.LFT3, Mid.LFT.traitmeans$Mid.LFT)] * biomass$SPP_O3_BASAL_AREA_FRACTION/100
Mid.LFT.Nmass4 <- Mid.LFT.traitmeans$Nmass[match(biomass$Mid.LFT4, Mid.LFT.traitmeans$Mid.LFT)] * biomass$SPP_O4_BASAL_AREA_FRACTION/100

biomass$Mid.LFT.Nmass <- apply(data.frame(Mid.LFT.Nmass1,Mid.LFT.Nmass2,Mid.LFT.Nmass3,Mid.LFT.Nmass4),MARGIN=1,FUN=sum, na.rm=T)



# calcualte community weighted Nmass by Shallow LFT
Shallow.LFT.Nmass1 <- Shallow.LFT.traitmeans$Nmass[match(as.character(biomass$Shallow.LFT1), as.character(Shallow.LFT.traitmeans$Shallow.LFT.larix))] * biomass$SPP_O1_BASAL_AREA_FRACTION/100
Shallow.LFT.Nmass2 <- Shallow.LFT.traitmeans$Nmass[match(biomass$Shallow.LFT2, Shallow.LFT.traitmeans$Shallow.LFT.larix)] * biomass$SPP_O2_BASAL_AREA_FRACTION/100
Shallow.LFT.Nmass3 <- Shallow.LFT.traitmeans$Nmass[match(biomass$Shallow.LFT3, Shallow.LFT.traitmeans$Shallow.LFT.larix)] * biomass$SPP_O3_BASAL_AREA_FRACTION/100
Shallow.LFT.Nmass4 <- Shallow.LFT.traitmeans$Nmass[match(biomass$Shallow.LFT4, Shallow.LFT.traitmeans$Shallow.LFT.larix)] * biomass$SPP_O4_BASAL_AREA_FRACTION/100

biomass$Shallow.LFT.Nmass <- apply(data.frame(Shallow.LFT.Nmass1,Shallow.LFT.Nmass2,Shallow.LFT.Nmass3,Shallow.LFT.Nmass4),MARGIN=1,FUN=sum, na.rm=T)






#######. Narea CWM by disaggregation method #################3
# calcualte community weighted Narea by PFT
PFT.Narea1 <- PFT.traitmeans$Narea[match(as.character(biomass$PFT1), as.character(PFT.traitmeans$PFT))] * biomass$SPP_O1_BASAL_AREA_FRACTION/100
PFT.Narea2 <- PFT.traitmeans$Narea[match(biomass$PFT2, PFT.traitmeans$PFT)] * biomass$SPP_O2_BASAL_AREA_FRACTION/100
PFT.Narea3 <- PFT.traitmeans$Narea[match(biomass$PFT3, PFT.traitmeans$PFT)] * biomass$SPP_O3_BASAL_AREA_FRACTION/100
PFT.Narea4 <- PFT.traitmeans$Narea[match(biomass$PFT4, PFT.traitmeans$PFT)] * biomass$SPP_O4_BASAL_AREA_FRACTION/100

biomass$PFT.Narea <- apply(data.frame(PFT.Narea1,PFT.Narea2,PFT.Narea3,PFT.Narea4),MARGIN=1,FUN=sum, na.rm=T)

# calcualte community weighted Narea by Deep LFT
Deep.LFT.Narea1 <- Deep.LFT.traitmeans$Narea[match(as.character(biomass$Deep.LFT1), as.character(Deep.LFT.traitmeans$Deep.LFT))] * biomass$SPP_O1_BASAL_AREA_FRACTION/100
Deep.LFT.Narea2 <- Deep.LFT.traitmeans$Narea[match(biomass$Deep.LFT2, Deep.LFT.traitmeans$Deep.LFT)] * biomass$SPP_O2_BASAL_AREA_FRACTION/100
Deep.LFT.Narea3 <- Deep.LFT.traitmeans$Narea[match(biomass$Deep.LFT3, Deep.LFT.traitmeans$Deep.LFT)] * biomass$SPP_O3_BASAL_AREA_FRACTION/100
Deep.LFT.Narea4 <- Deep.LFT.traitmeans$Narea[match(biomass$Deep.LFT4, Deep.LFT.traitmeans$Deep.LFT)] * biomass$SPP_O4_BASAL_AREA_FRACTION/100

biomass$Deep.LFT.Narea <- apply(data.frame(Deep.LFT.Narea1,Deep.LFT.Narea2,Deep.LFT.Narea3,Deep.LFT.Narea4),MARGIN=1,FUN=sum, na.rm=T)


# calcualte community weighted Narea by Mid LFT
Mid.LFT.Narea1 <- Mid.LFT.traitmeans$Narea[match(as.character(biomass$Mid.LFT1), as.character(Mid.LFT.traitmeans$Mid.LFT))] * biomass$SPP_O1_BASAL_AREA_FRACTION/100
Mid.LFT.Narea2 <- Mid.LFT.traitmeans$Narea[match(biomass$Mid.LFT2, Mid.LFT.traitmeans$Mid.LFT)] * biomass$SPP_O2_BASAL_AREA_FRACTION/100
Mid.LFT.Narea3 <- Mid.LFT.traitmeans$Narea[match(biomass$Mid.LFT3, Mid.LFT.traitmeans$Mid.LFT)] * biomass$SPP_O3_BASAL_AREA_FRACTION/100
Mid.LFT.Narea4 <- Mid.LFT.traitmeans$Narea[match(biomass$Mid.LFT4, Mid.LFT.traitmeans$Mid.LFT)] * biomass$SPP_O4_BASAL_AREA_FRACTION/100

biomass$Mid.LFT.Narea <- apply(data.frame(Mid.LFT.Narea1,Mid.LFT.Narea2,Mid.LFT.Narea3,Mid.LFT.Narea4),MARGIN=1,FUN=sum, na.rm=T)



# calcualte community weighted Narea by Shallow LFT
Shallow.LFT.Narea1 <- Shallow.LFT.traitmeans$Narea[match(as.character(biomass$Shallow.LFT1), as.character(Shallow.LFT.traitmeans$Shallow.LFT.larix))] * biomass$SPP_O1_BASAL_AREA_FRACTION/100
Shallow.LFT.Narea2 <- Shallow.LFT.traitmeans$Narea[match(biomass$Shallow.LFT2, Shallow.LFT.traitmeans$Shallow.LFT.larix)] * biomass$SPP_O2_BASAL_AREA_FRACTION/100
Shallow.LFT.Narea3 <- Shallow.LFT.traitmeans$Narea[match(biomass$Shallow.LFT3, Shallow.LFT.traitmeans$Shallow.LFT.larix)] * biomass$SPP_O3_BASAL_AREA_FRACTION/100
Shallow.LFT.Narea4 <- Shallow.LFT.traitmeans$Narea[match(biomass$Shallow.LFT4, Shallow.LFT.traitmeans$Shallow.LFT.larix)] * biomass$SPP_O4_BASAL_AREA_FRACTION/100

biomass$Shallow.LFT.Narea <- apply(data.frame(Shallow.LFT.Narea1,Shallow.LFT.Narea2,Shallow.LFT.Narea3,Shallow.LFT.Narea4),MARGIN=1,FUN=sum, na.rm=T)


biomass$log.PFT.LMA <- log(biomass$PFT.LMA, base=10)
biomass$log.Deep.LFT.LMA <- log(biomass$Deep.LFT.LMA, base=10)
biomass$log.Mid.LFT.LMA <- log(biomass$Mid.LFT.LMA, base=10)
biomass$log.Shallow.LFT.LMA <- log(biomass$Shallow.LFT.LMA, base=10)

biomass$log.PFT.LL <- log(biomass$PFT.LL, base=10)
biomass$log.Deep.LFT.LL <- log(biomass$Deep.LFT.LL, base=10)
biomass$log.Mid.LFT.LL <- log(biomass$Mid.LFT.LL, base=10)
biomass$log.Shallow.LFT.LL <- log(biomass$Shallow.LFT.LL, base=10)

biomass$log.PFT.Nmass <- log(biomass$PFT.Nmass, base=10)
biomass$log.Deep.LFT.Nmass <- log(biomass$Deep.LFT.Nmass, base=10)
biomass$log.Mid.LFT.Nmass <- log(biomass$Mid.LFT.Nmass, base=10)
biomass$log.Shallow.LFT.Nmass <- log(biomass$Shallow.LFT.Nmass, base=10)




################ . CW variance explained by each predictor ###########

log.CWexplans <- data.frame(Type = c("PFT","Deep LFT", "Mid LFT", "Shallow LFT"), 
                        log.LMA = c(summary(lm(log.cw_LMAp_if~log.PFT.LMA, biomass))$r.squared,
                                    summary(lm(log.cw_LMAp_if~log.Deep.LFT.LMA, biomass))$r.squared,
                                    summary(lm(log.cw_LMAp_if~log.Mid.LFT.LMA, biomass))$r.squared,
                                    summary(lm(log.cw_LMAp_if~log.Shallow.LFT.LMA, biomass))$r.squared),
                        log.LL = c(summary(lm(log.cw_LLp_if~log.PFT.LL, biomass))$r.squared,
                                    summary(lm(log.cw_LLp_if~log.Deep.LFT.LL, biomass))$r.squared,
                                    summary(lm(log.cw_LLp_if~log.Mid.LFT.LL, biomass))$r.squared,
                                    summary(lm(log.cw_LLp_if~log.Shallow.LFT.LL, biomass))$r.squared),
                        log.Nmass = c(summary(lm(log.cw_Nmassp_if~log.PFT.Nmass, biomass))$r.squared,
                                    summary(lm(log.cw_Nmassp_if~log.Deep.LFT.Nmass, biomass))$r.squared,
                                    summary(lm(log.cw_Nmassp_if~log.Mid.LFT.Nmass, biomass))$r.squared,
                                    summary(lm(log.cw_Nmassp_if~log.Shallow.LFT.Nmass, biomass))$r.squared))



CWexplans <- data.frame(Type = c("PFT","Deep LFT", "Mid LFT", "Shallow LFT"), 
                        LMA = c(summary(lm(cw_LMAp_if~PFT.LMA, biomass))$r.squared,
                                    summary(lm(cw_LMAp_if~Deep.LFT.LMA, biomass))$r.squared,
                                    summary(lm(cw_LMAp_if~Mid.LFT.LMA, biomass))$r.squared,
                                    summary(lm(cw_LMAp_if~Shallow.LFT.LMA, biomass))$r.squared),
                        LL = c(summary(lm(cw_LLp_if~PFT.LL, biomass))$r.squared,
                                   summary(lm(cw_LLp_if~Deep.LFT.LL, biomass))$r.squared,
                                   summary(lm(cw_LLp_if~Mid.LFT.LL, biomass))$r.squared,
                                   summary(lm(cw_LLp_if~Shallow.LFT.LL, biomass))$r.squared),
                        Nmass = c(summary(lm(cw_Nmassp_if~PFT.Nmass, biomass))$r.squared,
                                      summary(lm(cw_Nmassp_if~Deep.LFT.Nmass, biomass))$r.squared,
                                      summary(lm(cw_Nmassp_if~Mid.LFT.Nmass, biomass))$r.squared,
                                      summary(lm(cw_Nmassp_if~Shallow.LFT.Nmass, biomass))$r.squared))

rowMeans(CWexplans[,-1])
