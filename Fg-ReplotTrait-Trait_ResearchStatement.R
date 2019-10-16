 
#set color palette
mypal <- brewer.pal(n=9, "Set1")
palette(mypal)
# set the colors that wil denote within-species, within-genus, within family and across CWMs
colchoices <- c(1,2,4,3,6)

options(na.action = na.omit)


biomass <- read.csv("DerivedData/PNW_Biomass_data_for_Anderegg_etal_2018_EcolLet.csv", header=T, row.names = 1)
traits <- read.csv("DerivedData/PNW_Trait_data_for_Anderegg_etal_2018_EcolLet.csv", header=T, row.names=1)
data.all <- read.csv("DerivedData/AllTraitData_for_Anderegg_etal_2018_EcolLet.csv", header=T, row.names = 1)
# create a dataset of species mean traits just for PNW dataset:
spp.traits <- traits %>% group_by(SP.ID) %>% summarise(nsample = n(), SLA = mean(SLA_HSA, na.rm=T), nSLA = n()- length(which(is.na(SLA_HSA))), CN = mean(LEAF_CN, na.rm=T), nCN = n()- length(which(is.na(LEAF_CN)))
                                                       , LIFE = mean(LEAF_LIFE, na.rm=T), nLIFE=n()- length(which(is.na(LEAF_LIFE))), nplots =length(unique(PLOT_ID))
                                                       , mLMA = mean(LMA, na.rm=T), mLLmonths= mean(LLmonths, na.rm=T), mNmass = mean(LEAF_NITROGEN, na.rm=T), mNarea=mean(Narea, na.rm=T)
                                                       , mlog.LMA = mean(log.LMA, na.rm=T), mlog.LL = mean(log.LL, na.rm=T), mlog.Narea=mean(log.Narea, na.rm=T), mlog.Nmass=mean(log.Nmass, na.rm=T)
                                                       , climPC1 = mean(climPC1, na.rm=T), climPC2 = mean(climPC2, na.rm=T), climPC3 = mean(climPC3, na.rm=T)
                                                       , soil_N= mean(soil_N, na.rm=T), soil_pH=mean(soil_pH, na.rm=T), ASA = mean(ASA, na.rm=T), LAI_O=mean(LAI_O, na.rm=T), AG_TGROWTH = mean(AG_TGROWTH, na.rm=T))
# log species mean traits for between species analysis
spp.traits$log.LMA <- log(spp.traits$mLMA, base=10)
spp.traits$log.LL <- log(spp.traits$mLLmonths, base=10)
spp.traits$log.Nmass <- log(spp.traits$mNmass, base=10)
spp.traits$log.Narea <- log(spp.traits$mNarea, base=10)



########### . Create w/in spp, species means, genus means, and family means datasets ###############
# then creating a dataset of genera w/ >5 species, and families w/ >5 genera
# select species w/ > 5 measurements
commonspp <-  names(which(xtabs(~Species, data.all)>=5))
spp.data<- data.all %>% filter(Species %in% commonspp)
spp.data$Species <- factor(spp.data$Species)
spp.data$Genus <- factor(spp.data$Genus)
spp.data$Family <- factor(spp.data$Family)


########## Species level trait averages ______________________________________________

allspp <- data.all %>% group_by(Species, Genus, Family) %>% summarise( lslog.LL = mean(log.LL, na.rm=T), lslog.LMA = mean(log.LMA, na.rm=T), lslog.Nmass = mean(log.Nmass, na.rm=T), lslog.Narea = mean(log.Narea, na.rm=T)
                                                                       ,slog.LL = log(mean(10^log.LL, na.rm=T),base=10), slog.LMA = log(mean(10^log.LMA, na.rm=T),base=10), slog.Nmass = log(mean(10^log.Nmass, na.rm=T),base=10), slog.Narea = log(mean(10^log.Narea, na.rm=T),base=10)
                                                                       ,MAT = mean(MAT, na.rm=T), MAP=mean(MAP, na.rm=T), VPD=mean(VPD, na.rm=T) )
colnames(allspp) <- gsub("slog", "log", colnames(allspp))




#### Genus means ___________________________________________________________________
allgen <- allspp %>% group_by(Genus)  %>% summarise(Fam = sort(unique(Family))[1], lglog.LL = mean(llog.LL, na.rm=T),lglog.LMA = mean(llog.LMA, na.rm=T),lglog.Nmass = mean(llog.Nmass, na.rm=T), lglog.Narea = mean(llog.Narea, na.rm=T)
                                                    ,glog.LL = log(mean(10^log.LL, na.rm=T),base=10), glog.LMA = log(mean(10^log.LMA, na.rm=T),base=10), glog.Nmass = log(mean(10^log.Nmass, na.rm=T),base=10), glog.Narea = log(mean(10^log.Narea, na.rm=T),base=10), nspp = n() 
                                                    ,MAT = mean(MAT, na.rm=T), MAP=mean(MAP, na.rm=T), VPD=mean(VPD, na.rm=T) )
# taking mean of raw values or logged values doesn't matter all that much yet
colnames(allgen) <- gsub("glog", "log", colnames(allgen))
colnames(allgen)[2] <- "Family"
# 939 genera from 211 families
# allgen <- allspp %>% group_by(Genus)  %>% summarise(Fam = sort(unique(Family))[1], Needle.Broad = unique(na.omit(NeedleBroad))[1], glog.LL = mean(log.LL, na.rm=T), glog.LMA = mean(log.LMA, na.rm=T), glog.Nmass = mean(log.Nmass, na.rm=T),glog.Narea = mean(log.Narea, na.rm=T)
#                                                     ,rglog.LL = log(mean(10^rlog.LL, na.rm=T),base=10), rglog.LMA = log(mean(10^rlog.LMA, na.rm=T),base=10), rglog.Nmass = log(mean(10^rlog.Nmass, na.rm=T),base=10), nspp = n() )



#### Family Means _______________________________________________________________________
allfam <- allgen %>% group_by(Family) %>% summarise(lflog.LL = mean(llog.LL, na.rm=T), lflog.LMA = mean(llog.LMA, na.rm=T), lflog.Nmass = mean(llog.Nmass, na.rm=T),lflog.Narea = mean(llog.Narea, na.rm=T)
                                                    ,flog.LL = log(mean(10^log.LL, na.rm=T),base=10), flog.LMA = log(mean(10^log.LMA, na.rm=T),base=10), flog.Nmass = log(mean(10^log.Nmass, na.rm=T),base=10), flog.Narea = log(mean(10^log.Narea, na.rm=T), base=10)
                                                    , tnspp = sum(nspp), ngen=n() )
colnames(allfam) <- gsub("flog", "log", colnames(allfam))

#### Family Means for fams w/ > 3 species
fam.data <- allfam # new: 211 families (up from 189)
fam.dataclean <- allfam[which(allfam$tnspp>2),] # 101 families, up from 97 families




##### dataset of only genera w/ >5 species _____________________________________________
# for w/in genus analysis
gen.data <- allspp[which(allspp$Genus %in% names(which(xtabs(~Genus, allspp)>=5))),]
# 750 measurements of 73 genera.
gen.data$Genus <- factor(gen.data$Genus) # get rid of unused levels


##### Spp. in Family level data, for families w/ >5 species in them ________________________
sppinfam.data <- allspp[which(allspp$Family %in% names(which(xtabs(~Family,allspp)>=5))),]
sppinfam.data$Family <- factor(sppinfam.data$Family) # get rid of unused levels



#### gen in Family level data, for families w/ >5 genera in them _________________________
# 556 measurements of 40 Families , added 2 Families from last batch, and probably quite a few obs (old n_LMA.N = 459)
# 04.01.17 - 684 obs of 50 Families. again, big win from taxonomic cleaning!
geninfam.data <- allgen[which(allgen$Family %in% names(which(xtabs(~Family,allgen)>=5))),]
geninfam.data$Family <- factor(geninfam.data$Family) # get rid of unused levels
geninfam.data$Genus <- factor(geninfam.data$Genus)




all.results.cwm.cl <- read.csv("DerivedData/SlopesAndCorrelations_Alltaxa_clean.csv", header=T, row.name=1)


# calculate significance of a pearson correlation coefficient
rho.sig <- function(rho, n){
  t <- rho/sqrt((1-rho^2)/(n-2))
  pt(-abs(t), n-2)*2
}





#_____________________________________________________________________
######## BEGIN: Figure 2-4 Plotting ##########
#_____________________________________________________________________

# function that plots the major axis regression line for a dataset
plot.MAR <- function(xvar, yvar, data, method="SMA", linecol, lwd=1, lty=1) {
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
      xvals <- sort(tmp.mod$x)
      yhat <- xvals * slope + intercept
      lines(yhat~xvals, col=linecol, lwd=lwd, lty=lty)
      
    }
  }
}


###### .FIG 3:  LMA vs LL ###################
colchoices <- c(1,2,4,3,6)
palette(mypal[colchoices])
cex.sig <- 1.1
cex.ns <- .9
crit <- 0.05

mat <- matrix(c(1,3,
                2,3), nrow=2, byrow = T)

#jpeg(width=7.008, height=4, units = "in", res=600,filename = paste0("./",results_dirname,"/Fig3_LMA_LL.jpeg"))
#pdf(width=7.008, height=4, file = paste0("./",results_dirname,"/Fig3_LMA_LL3.pdf"))
quartz(width=7.008, height=4)
layout(mat, heights=c(1,1.5))
par(mar=c(0,4,1.5,1))
p <- boxplot(Slope_LMA.LL~Type, all.results.cwm.cl[which(all.results.cwm.cl$n_LMA.LL>5& !all.results.cwm.cl$Type %in% c("global","Famclean")),], plot=F
             , xlim=c(.5,6.5), ylim=c(-2.5,4))
boxplot(Slope_LMA.LL~Type, all.results.cwm.cl[which(all.results.cwm.cl$n_LMA.LL>5 & !is.na(all.results.cwm.cl$Taxo.Unit)),]
        , ylim=c(-2.5,4),las=3, ylab="SMA Slope" #, main="log(LMA)~log(Nmass) strict"
        , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
        ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
        , staplelwd=0, outpch=NA, outcex=.7, outcol=mypal[colchoices]
        , boxwex=.7, xaxt="n"
        , xlim=c(.5,6.5))
abline(h=0, lty=2)
m <- lmodel2(log.LL~log.LMA, allspp)
points(y=m$regression.results$Slope[3],x=5, pch=16, cex=1.3, col="darkgrey")
arrows(x0=5,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`2.5%-Slope`[3], length = 0,lwd=3, col="darkgrey")
arrows(x0=5,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`97.5%-Slope`[3], length = 0,lwd=3, col="darkgrey")
m <- lmodel2(log.LL~log.LMA, fam.dataclean)
points(y=m$regression.results$Slope[3],x=4, pch=16, cex=1.3, col="black")
arrows(x0=4,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`2.5%-Slope`[3], length = 0,lwd=3)
arrows(x0=4,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`97.5%-Slope`[3], length = 0,lwd=3)
m <- lmodel2(log.cw_LLp_if~log.cw_LMAp_if, biomass)
arrows(x0=6,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`2.5%-Slope`[3], length = 0,lwd=3, col=mypal[5])
arrows(x0=6,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`97.5%-Slope`[3], length = 0,lwd=3, col=mypal[5])
#points(y=m$regression.results$Slope[3],x=6, pch=24, cex=1.3, col="black", bg=mypal[5])
points(y=m$regression.results$Slope[3],x=6, pch=17, cex=1.3, col=mypal[5])
# # only evergreens
# m <- lmodel2(log.LL~log.LMA, LES[which(LES$Decid.E.green=="E"),])
# points(y=m$regression.results$Slope[3],x=5, pch=1, cex=1.3, col="darkgrey")
# arrows(x0=5,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`2.5%-Slope`[3], length = 0,lwd=3, col="darkgrey")
# arrows(x0=5,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`97.5%-Slope`[3], length = 0,lwd=3, col="darkgrey")

tmp <- all.results.cwm.cl[which(all.results.cwm.cl$n_LMA.LL>5 & !is.na(all.results.cwm.cl$sig_LMA.LL)),]
tmp$rho.sig <- rho.sig(rho=tmp$Rho_LMA.LL, n= tmp$n_LMA.LL)
tmp.ns <- tmp[which(tmp$rho.sig>0.1),]
tmp.sig <- tmp[which(tmp$rho.sig<=0.1),]
#palette(paste0(mypal[colchoices], "55"))
set.seed(42)
points(Slope_LMA.LL~jitter(as.numeric(Type)), tmp.ns, pch=1, col=Type, cex=cex.ns)
#palette(mypal[colchoices])
set.seed(42)
points(Slope_LMA.LL~jitter(as.numeric(Type)), tmp.sig, pch=16, col=Type,cex=cex.sig)
#set.seed(42)
#arrows(x0 = jitter(as.numeric(tmp.sig$Type)),y0=tmp.sig$Slope.lci_LMA.LL, y1=tmp.sig$Slope.uci_LMA.LL, length = 0, col=tmp.sig$Type)

#points(Slope~jitter(rep(1, times=nrow(vib.results[which(vib.results$rho.sig<=0.1),])),amount=.3), vib.results[which(vib.results$rho.sig<=0.1),], pch=16, col=mypal[colchoices[1]])
#points(Slope~jitter(rep(1, times=nrow(vib.results[which(vib.results$rho.sig>0.1),])), amount=.3), vib.results[which(vib.results$rho.sig>0.1),], pch=1, col=mypal[colchoices[1]])

#text(y=par()$usr[3]+.2, x=c(1,2,3,4,5,6), labels = p$n)
#text(y=par()$usr[4]-.2,x=.5, labels = "a)")
mtext(text = "a)", side = 3, adj=0, line=.2)
mtext(text= "LL vs LMA", side=3, line=.2)
par(mar=c(6,4,0,1))
#p <- boxplot(Rho_LMA.LL~Type, all.results.cwm.cl[which(all.results.cwm.cl$n_LMA.LL>5& !all.results.cwm.cl$Type %in% c("global","Famclean")),]
boxplot(Rho_LMA.LL~Type, all.results.cwm.cl[which(all.results.cwm.cl$n_LMA.LL>5 &  is.na(all.results.cwm.cl$Taxo.Unit)),]
        , ylim=c(-1,1.4),las=3, ylab="Correlation"
        , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
        ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
        , staplelwd=0, outpch=NA, outcex=.7, outcol=mypal[colchoices]
        , boxwex=.7, xaxt="n"
        , xlim=c(.5,6.5))
axis(1,labels = c("w/in Spp","w/in Gen", "w/in Fam","btw Fam","Global","PNW CWM"), at=c(1,2,3,4,5,6), las=3)
abline(h=0, lty=2)
text(y=par()$usr[4]-.2, x=c(1,2,3), labels = p$n[1:3])
text(y=par()$usr[4]-.2, x=c(4,5,6), labels=paste0("(",all.results.cwm.cl$n_LMA.LL[c((nrow(all.results.cwm.cl)-2),(nrow(all.results.cwm.cl)-1),nrow(all.results.cwm.cl))],")"),cex=.9)
#polygon(x = c(1,2,3,3,2,1), y=na.omit(c(meancis$m10_LMA.LL, rev(meancis$m90_LMA.LL))), border="#EEEEEE",col="#EEEEEE")#lightgrey",col = "lightgrey")
m <- cor.test(x=allspp$log.LMA, y=allspp$log.LL)
points(y=m$estimate,x=5, pch=16, cex=1.3, col="darkgrey")
arrows(x0=5,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3, col="darkgrey")
m <- cor.test(x=fam.dataclean$log.LMA, y=fam.dataclean$log.LL)
points(y=m$estimate,x=4, pch=16, cex=1.3, col="black")
arrows(x0=4,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3)
m <- cor.test(x=biomass$log.cw_LMAp_if, y=biomass$log.cw_LLp_if)
arrows(x0=6,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3, col=mypal[5])
#points(y=m$estimate,x=6, pch=24, cex=1.3, col="black", bg=mypal[5])
points(y=m$estimate,x=6, pch=17, cex=1.3, col=mypal[5])
# # only evergreens
# m <- cor.test(x=LES$log.LMA[which(LES$Decid.E.green=="E")], y=LES$log.LL[which(LES$Decid.E.green=="E")])
# points(y=m$estimate,x=5, pch=1, cex=1.3, col="darkgrey")
# arrows(x0=5,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3, col="darkgrey")

# layering points on top of null model
set.seed(62)
points(Rho_LMA.LL~jitter(as.numeric(Type)), all.results.cwm.cl[which(all.results.cwm.cl$n_LMA.LL>5 & all.results.cwm.cl$sig_LMA.LL<=crit),], pch=16, col=Type, cex=cex.sig)
#palette(paste0(mypal[colchoices], "55"))
set.seed(42)
points(Rho_LMA.LL~jitter(as.numeric(Type)), all.results.cwm.cl[which(all.results.cwm.cl$n_LMA.LL>5 & all.results.cwm.cl$sig_LMA.LL>crit),], pch=1, col=Type ,cex=cex.ns)
#palette(mypal[colchoices])
#axis(1,labels = c("w/in Spp","w/in Gen", "w/in Fam","btw Fam","global"), at=c(1,2,3,4,5), las=3)

#points(Rho~jitter(rep(1, times=nrow(vib.results[which(vib.results$rho.sig<=0.1),])),amount = .5), vib.results[which(vib.results$rho.sig<=0.1),], pch=16, col=mypal[colchoices[1]])
#points(Rho~jitter(rep(1, times=nrow(vib.results[which(vib.results$rho.sig>0.1),]))), vib.results[which(vib.results$rho.sig>0.1),], pch=1, col=paste0(mypal[colchoices[1]], "55"))


## scatterplot LMA LL _________________________
par(mar=c(4,4,1.5,1.5), mgp=c(2.5,1,0), cex.lab=1.3)
ltysig <- 1
ltyns <-  2
lwdsig <-  1.8
lwdns <- 1

plot(log.LL~log.LMA,allspp, col="grey", pch=16, ylab=expression(paste(log[10],"(Leaf Lifespan)")), xlab=expression(paste(log[10](LMA))))

tax <- "w.inGen"
for (i in as.character(all.results.cwm.cl$Taxo.Unit[which(all.results.cwm.cl$Type==tax & all.results.cwm.cl$n_LMA.LL>5)])){
  linety <- ltyns
  linewd <- lwdns
  if(all.results.cwm.cl$rho.sig_LMA.LL[which(all.results.cwm.cl$Taxo.Unit==i)]<0.05) {linety <- ltysig; linewd <- lwdsig}
  plot.MAR(xvar = "log.LMA", yvar = "log.LL",data= gen.data[which(gen.data$Genus==i),], linecol = mypal[colchoices[2]], lty=linety, lwd=linewd)
}

tax <- "Genw.inFam"
for (i in as.character(all.results.cwm.cl$Taxo.Unit[which(all.results.cwm.cl$Type==tax & all.results.cwm.cl$n_LMA.LL>5)])){
  linety <- ltyns
  linewd <- lwdns
  if(all.results.cwm.cl$rho.sig_LMA.LL[which(all.results.cwm.cl$Taxo.Unit==i)]<0.05) {linety <-ltysig ; linewd=lwdsig}
  plot.MAR(xvar = "log.LMA", yvar = "log.LL",data= geninfam.data[which(geninfam.data$Family==i),], linecol = mypal[colchoices[3]], lty=linety, lwd=linewd)
}


#abline(a=all.results.cwm.cl$Int_LMA.LL[which(all.results.cwm.cl$Type=="Famclean")], b=all.results.cwm.cl$Slope_LMA.LL[which(all.results.cwm.cl$Type=="Famclean")], lwd=3, col="black")
abline(a=all.results.cwm.cl$Int_LMA.LL[which(all.results.cwm.cl$Type=="global")], b=all.results.cwm.cl$Slope_LMA.LL[which(all.results.cwm.cl$Type=="global")], lwd=3, col="black", lty=1)
tax <- "w.inSpp"
for (i in as.character(all.results.cwm.cl$Taxo.Unit[which(all.results.cwm.cl$Type==tax & all.results.cwm.cl$n_LMA.LL>5)])){
  linety <- ltyns
  linewd <- lwdns
  if(all.results.cwm.cl$rho.sig_LMA.LL[which(all.results.cwm.cl$Taxo.Unit==i)]<0.05) {linety <- ltysig; linewd <- lwdsig}
  plot.MAR(xvar = "log.LMA", yvar = "log.LL",data= spp.data[which(spp.data$Species==i),], linecol = mypal[colchoices[1]], lty =linety, lwd=linewd)
}
# # add in Mt. Rainier
# for (i in c("Aster alpigenus","Castilleja parviflora", "Erythronium montanum", "Lupinus arcticus","Valeriana sitchensis","Veratrum viride")){
#   linety <- ltyns
#   linewd <- lwdns
#   if(all.results.cwm.cl$rho.sig_LMA.LL[which(all.results.cwm.cl$Taxo.Unit==i)]<0.05) {linety <- ltysig; linewd <- lwdsig}
#   plot.MAR(xvar = "log.LMA", yvar = "log.LL",data= Rainier[which(Rainier$Species_full==i),], linecol = mypal[colchoices[1]], lty =linety, lwd=linewd)
# }

#points(biomass$log.cw_LLp_if~biomass$log.cw_LMAp_if, pch=24, bg=mypal[5])
#all.results.cwm.cl$rho.sig_LMA.LL[170] # p of rho=0.08
#plot.MAR(xvar = "log.cw_LMAp_if", yvar="log.cw_LLp_if", data=biomass, linecol = mypal[5], lwd = 3, lty=ltyns)

legend('bottomright',legend = c("Spp means","Global", "Gen w/in Fam","Spp w/in Gen","Ind w/in Spp"),col = c("grey","black",mypal[colchoices[c(3,2,1)]]), pch = c(16,NA,NA,NA,NA), lty=c(NA,1,1,1,1), lwd=c(NA,3,2,2,2), cex=.9, bg="white")
      