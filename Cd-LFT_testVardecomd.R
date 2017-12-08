##### Variance Decomp of LES dataset for LFT analysis ######


## load LES from Cd-Initial_Analysis.R
higher.taxa <- function(fam.names){
  output <- data.frame(Order=rep(NA, length(fam.names)),Subclass=rep(NA, length(fam.names)),Class=rep(NA, length(fam.names)), Subdivision=rep(NA, length(fam.names)),Division=rep(NA, length(fam.names)))
  for(i in 1:length(fam.names)){
    tmp <- classification(fam.names[i], db="itis")[[1]]
    output[i,] <- cbind(tmp$name[which(tmp$rank=="order")], tmp$name[which(tmp$rank=="subclass" | tmp$rank=="superorder")], tmp$name[which(tmp$rank=="class")], tmp$name[which(tmp$rank=="subdivision")], tmp$name[which(tmp$rank=="division")])
  }
  return(output)
}
LEStaxa <- higher.taxa(fam.names=c("Pinaceae","Fagaceae"))

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
baad.all$dictionary[,c("variable", "label")]

baad <- baad.all$data
baad$genus <- sapply(X=baad$speciesMatched, FUN = function(x){strsplit(x, split = " ")[[1]][1]})
baad$studyName <- factor(baad$studyName)
baad$genus <- factor(baad$genus)
baad$family <- factor(baad$family)
baad$Al.As <- baad$a.lf/baad$a.ssbh
baad$Ab.Bl <- baad$m.so/baad$m.rt
baad$lmf <- baad$m.lf / (baad$m.lf + baad$m.st)
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
baad$log.h.t <- log(baad$h.t, base=10)




#### Variance Decomp for BAAD
WDvar <- lmer(r.st~ 1 + (1|family/genus/speciesMatched), baad)
logWDvar <- lmer(log.r.st~ 1 + (1|family/genus/speciesMatched), baad)
logAbBlvar <- lmer(log.Ab.Bl~ 1 + (1|family/genus/speciesMatched), baad)
logHeightvar <- lmer(h.t~ 1 + (1|family/genus/speciesMatched), baad)
LMFvar <- lmer(lmf~ 1 + (1|family/genus/speciesMatched), baad)

rawWDvariance <- data.frame(VarCorr(WDvar))
WDvariance <- data.frame(VarCorr(logWDvar))
AbBlvariance <- data.frame(VarCorr(logAbBlvar))
Heightvariance <- data.frame(VarCorr(logHeightvar))
LMFvariance <- data.frame(VarCorr(LMFvar))




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
####### Global Wood Density Database ######
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

GWD <- read.csv("/Users/leeanderegg/Dropbox/WD project/GlobalWoodDensityDatabase.csv", header=T)
GWD$Genus <- factor(sapply(X=as.character(GWD$Binomial), FUN = function(x){strsplit(x, split = " ")[[1]][1]}))
# bunch of replicated genera, and 100+ families.

colnames(GWD)[4] <- "WD"

#WDvar2 <- lmer(WD~ 1 + (1|Genus/Binomial), GWD)
WDvar2 <- lmer(WD~ 1 + (1|Family) + (1|Genus) + (1|Binomial), GWD)
logWDvar2 <- lmer(log(WD, base=10)~ 1 + (1|Family) + (1|Genus) + (1|Binomial), GWD)

rawWDvariance2 <- data.frame(VarCorr(WDvar2))
WDvariance2 <- data.frame(VarCorr(logWDvar2))



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
######## R:S database (2017) ###########
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


RtS <- read.csv("/Users/leeanderegg/Desktop/LFT analysis/nph14863-sup-0002-NotesS1.csv", header=T)
RtS$vegetationTYPE[which(RtS$vegetationTYPE=="TropMF ")] <- "TropMF"
RtS$vegetationTYPE <- factor(RtS$vegetationTYPE)
RtS$log.AGB <- log(RtS$AGB_kg, base=10)
RtS$log.BGB <- log(RtS$BGB_kg, base=10)
RtS$RS_RATIO[which(RtS$RS_RATIO==0)] <- NA
RtS$log.RSR <- log(RtS$RS_RATIO, base=10)

RtSvar <- lmer(RS_RATIO~ 1 + (1|FAMILY/GENUS/SPECIES), RtS)
logRtSvar <- lmer(log.RSR~ 1 + (1|FAMILY/GENUS/SPECIES), RtS[which(!is.na(RtS$SPECIES)& !is.na(RtS$GENUS) & !is.na(RtS$FAMILY)),])
RtSvariance <- data.frame(VarCorr(logRtSvar))


pal <- brewer.pal(n=9, "Set1")
palette(pal)

plot(log.AGB~log.BGB, RtS, col=vegetationTYPE)
for(i in 1:length(levels(RtS$vegetationTYPE))){
  plot.MAR(xvar = "log.BGB", yvar="log.AGB", data = RtS[which(RtS$vegetationTYPE==levels(RtS$vegetationTYPE)[i]),], linecol = mypal[i], lwd=2 )
}

ggplot(RtS, aes(x=log.BGB, y=log.AGB, col= vegetationTYPE)) + geom_point() + geom_smooth(method="lm") +
  facet_wrap(facets= ~vegetationTYPE)


xtabs(~FAMILY, RtS)


#### Combine all the variance decomps together
baadtraitvars <- data.frame(rawWDvariance2[,4], WDvariance2[,4], rawWDvariance[,4], WDvariance[,4], AbBlvariance[,4], Heightvariance[,4], LMFvariance[,4], RtSvariance[,4])
colnames(baadtraitvars) <- c("GWDrawWD","GWDlogWD","BAADrawWD", "BAADlogWD", "logAGB_BGB", "logHeight", "LMF", "RtoS")
rownames(baadtraitvars) <- c("BtwSpecies", "BtwGenera" , "BtwFamilies", "WtinSpp")

baadtraitvars_scaled <- baadtraitvars  
for(i in 1:ncol(baadtraitvars)){
  baadtraitvars_scaled[,i] <- baadtraitvars[,i]/sum(baadtraitvars[,i])
}


baadtraitvars_scaled2 <- baadtraitvars_scaled[c(3,2,1,4),]

baadtraitvars_scaled3 <- baadtraitvars_scaled2[,order(baadtraitvars_scaled2[1,], decreasing =T)]



###### Plotting them #####
#quartz(width=3.75, height=4)
quartz(width=6.81, height=3.5)
par(mgp=c(3,.7,0), cex.lab=1.3, cex.axis=1.1, mfrow=c(1,1), mar=c(6,2,3,6), oma=c(0,3,0,0))
cols <- rev(c(mypal[colchoices[c(1,2,3)]],"#CCCCCC"))#brewer.pal(11, "RdBu")[c(11,9,1, 6)]
# With old colors, alpha=CC, with new colors, alpha=99
bp <-barplot(as.matrix(baadtraitvars_scaled3),beside=F,legend.text = F,xpd = T,las=2,args.legend = list(x=4, y=1.3, ncol=2), col = paste0(cols,"99")
            , names.arg = colnames(baadtraitvars_scaled3) #c("log(WD)","raw WD", "Leaf Mass\nFraction", "log(Height)","AboveGrnd:\nBelowGrnd\nbiomass")
             , ylab= "Proportion of total Variance" #ylab="Proportion of total Variance\n(log traits)"
             , xlab="", mgp=c(2.5,.8,0))
#legend(xpd=T, x = 0, y=1.2, legend=rev(c("btw Fams","w/in Fam","w/in Gen","w/in Spp")), fill=rev(paste0(cols,"99")), ncol=2, bty="n",  cex=1.2)
#text(x = bp, y= par("usr")[3]-.05,labels =  c(expression(paste(log[10](LMA))), expression(paste(log[10](LL))),expression(paste(log[10](N[mass]))),expression(paste(log[10](N[area])))), srt=40, adj=1,xpd=NA, cex=1.4, font=1)
#mtext("a)", side=3, adj=-.1, line=1.3)
mtext("Proportion of Total Variance", side=2, line=2.8)
mtext("BAAD Allometry data", side=3, line=0.3)
#quartz(width=5, height=4)
#cols <- brewer.pal(11, "RdBu")[c(11,9,1, 6)]
#legend(xpd=T, x = 11, y=0.7, legend=rev(c("btw Fams","btw Gen","btw Spp")), fill=rev(paste0(cols[1:3],"99")), ncol=1, bty="n",  cex=1.2)
# mtext(text="28 spp w/ replication, 1000+ spp,\n500+ genera, 150+ families",side = 1,line = 3.3)
# legend(xpd=T, x = 0, y=1.3, legend=rownames(rtraitvars_scaled2), fill=paste0(cols,"CC"), ncol=2, bty="n",  cex=1.2)
legend(xpd=NA, x = 10, y=0.7, legend=rev(c("btw Fams","btw Gen","btw Spp","w/in Spp")), fill=rev(paste0(cols,"99")), ncol=1, bty="n",  cex=1)


