#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
###### Plotting the LMA-LL relationship
####### broken apart by scale for
####### for presentations
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#created: 12/8/2018 for UCSB job talk

# updated object calls to use dataframes created in Cd-Anderegg2018_Manuscript_Analysis_v2.R

require(RColorBrewer)
mypal <- brewer.pal(n=9, "Set1")
palette(mypal)

plotbars <- F
plotpoints <- T

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
if(plotbars==T){
  boxplot(Slope_LMA.LL~Type, all.results.cwm.cl[which(all.results.cwm.cl$n_LMA.LL>5 & !all.results.cwm.cl$Type %in% c("global","Famclean", "CWM")),]#!is.na(all.results.cwm.cl$Taxo.Unit)),]
          , ylim=c(-2.5,4),las=3, ylab="Slope" #, main="log(LMA)~log(Nmass) strict"
          , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
          ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
          , staplelwd=0, outpch=NA, outcex=.7, outcol=mypal[colchoices]
          , boxwex=.7, xaxt="n"
          , xlim=c(.5,6.5))
} else{
  boxplot(Slope_LMA.LL~Type, all.results.cwm.cl[which(all.results.cwm.cl$n_LMA.LL>5 & is.na(all.results.cwm.cl$Taxo.Unit)),]
        , ylim=c(-2.5,4),las=3, ylab="SMA Slope" #, main="log(LMA)~log(Nmass) strict"
        , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
        ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
        , staplelwd=0, outpch=NA, outcex=.7, outcol=mypal[colchoices]
        , boxwex=.7, xaxt="n"
        , xlim=c(.5,6.5))
}
abline(h=0, lty=1, col='darkgrey')
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
if(plotpoints==T) {
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
}
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
if(plotbars == T){
  boxplot(Rho_LMA.LL~Type, all.results.cwm.cl[which(all.results.cwm.cl$n_LMA.LL>5 & !all.results.cwm.cl$Type %in% c("global","Famclean", "CWM")),]
          , ylim=c(-1,1.4),las=3, ylab="Correlation"
          , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
          ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
          , staplelwd=0, outpch=NA, outcex=.7, outcol=mypal[colchoices]
          , boxwex=.7, xaxt="n"
          , xlim=c(.5,6.5))
} else{
  boxplot(Rho_LMA.LL~Type, all.results.cwm.cl[which(all.results.cwm.cl$n_LMA.LL>5 &  is.na(all.results.cwm.cl$Taxo.Unit)),]
        , ylim=c(-1,1.4),las=3, ylab="Correlation"
        , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
        ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
        , staplelwd=0, outpch=NA, outcex=.7, outcol=mypal[colchoices]
        , boxwex=.7, xaxt="n"
        , xlim=c(.5,6.5))
}
axis(1,labels = c("w/in Spp","w/in Gen", "w/in Fam","btw Fam","Global","PNW CWM"), at=c(1,2,3,4,5,6), las=3)
abline(h=0, lty=1, col='darkgrey')
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

if(plotpoints==T){
  
# layering points on top of null model
set.seed(62)
points(Rho_LMA.LL~jitter(as.numeric(Type)), all.results.cwm.cl[which(all.results.cwm.cl$n_LMA.LL>5 & all.results.cwm.cl$sig_LMA.LL<=crit),], pch=16, col=Type, cex=cex.sig)
#palette(paste0(mypal[colchoices], "55"))
set.seed(42)
points(Rho_LMA.LL~jitter(as.numeric(Type)), all.results.cwm.cl[which(all.results.cwm.cl$n_LMA.LL>5 & all.results.cwm.cl$sig_LMA.LL>crit),], pch=1, col=Type ,cex=cex.ns)
#palette(mypal[colchoices])
#axis(1,labels = c("w/in Spp","w/in Gen", "w/in Fam","btw Fam","global"), at=c(1,2,3,4,5), las=3)
}
#points(Rho~jitter(rep(1, times=nrow(vib.results[which(vib.results$rho.sig<=0.1),])),amount = .5), vib.results[which(vib.results$rho.sig<=0.1),], pch=16, col=mypal[colchoices[1]])
#points(Rho~jitter(rep(1, times=nrow(vib.results[which(vib.results$rho.sig>0.1),]))), vib.results[which(vib.results$rho.sig>0.1),], pch=1, col=paste0(mypal[colchoices[1]], "55"))


## scatterplot LMA LL _________________________
par(mar=c(4,4,1.5,1.5), mgp=c(2.5,1,0), cex.lab=1.3)
ltysig <- 1
ltyns <-  2
lwdsig <-  1.8
lwdns <- 1

plot(log.LL~log.LMA,allspp, col="grey", pch=16, ylab="Leaf Lifespan (months)", xlab="Leaf Mass per Area (g/m2)", #ylab=expression(paste(log[10],"(Leaf Lifespan)")), xlab=expression(paste(log[10](LMA))),
     xaxt="n", yaxt="n")

axis(2, at=log(c(1,2,3,4,5,10,20,30,40,50,100,200,300), base=10), labels = F)#c(1,"","","",5,"","","","",50,"",200,""),)
axis(2, at=log(c(1,5,50,300), base=10), labels = c(1,5,50,300))
axis(1, at=log(c(10,20,30,40,50,100,200,300,400,500,1000), base=10), labels = F)# c(10,"","","",50,"","","","","","1000"))
axis(1, at=log(c(20,50,500), base=10), labels =  c(20,50,500))

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





#________________________________________________________________
###### FIG 3 FOR PRESENTATIONS: LMA v LL Relations at different scales ######
#________________________________________________________________
# 2 panel figure with boxplots and 
# a-d: boxplots of MA slope and Rho for NmassvLMA and NmassvLL
# e&f: funnel plots 
# width = 1.5 columns -> 11.4 cm,4.89
#       = 2 columns -> 17.8 cm,7.008
# I think the easiest thing will be to make two quartez: top with boxplots and bottom with funnel

mat <- matrix(c(1,3,
                2,3), nrow=2, byrow = T)
colchoices <- c(1,2,4,3,6,5)
palette(mypal[colchoices])

quartz(width=7.008, height=4)
layout(mat)
par(mar=c(0,4,6,1))
p <- boxplot(Slope_LMA.LL~Type, all.results.cwm.cl[which(all.results.cwm.cl$n_LMA.LL>5 & !all.results.cwm.cl$Type %in% c("Famclean","global","CWM")),]
             , ylim=c(-2.5,4),las=3, ylab="Slope" #, main="log(LMA)~log(Nmass) strict"
             , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
             ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
             , staplelwd=0, outpch=16, outcex=.5, outcol=mypal[colchoices]
             , boxwex=.7, xaxt="n")
m <- lmodel2(log.LL~log.LMA, allspp)
points(y=m$regression.results$Slope[3],x=5, pch=16, cex=1.3, col="darkgrey")
arrows(x0=5,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`2.5%-Slope`[3], length = 0,lwd=3, col="darkgrey")
arrows(x0=5,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`97.5%-Slope`[3], length = 0,lwd=3, col="darkgrey")
m <- lmodel2(log.LL~log.LMA, fam.dataclean)
points(y=m$regression.results$Slope[3],x=4, pch=17, cex=1.3, col="black")
arrows(x0=4,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`2.5%-Slope`[3], length = 0,lwd=3)
arrows(x0=4,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`97.5%-Slope`[3], length = 0,lwd=3)


abline(h=0, lty=2)
#text(y=par()$usr[3]+.2, x=c(1,2,3,4,5,6), labels = p$n)
#text(y=par()$usr[4]-.2,x=.5, labels = "a)")
mtext(text = "a)", side = 3, adj=0, line=.2)
mtext(text= "LMA vs LL", side=3, line=.2)
par(mar=c(6,4,0,1))
p <- boxplot(Rho_LMA.LL~Type, all.results.cwm.cl[which(all.results.cwm.cl$n_LMA.LL>5& !all.results.cwm.cl$Type %in% c("Famclean","global", "CWM")),]
             , ylim=c(-1,1.4),las=3, ylab="Correlation"
             , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
             ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
             , staplelwd=0, outpch=16, outcex=.5, outcol=mypal[colchoices]
             , boxwex=.7, xaxt="n")
m <- cor.test(x=allspp$log.LL, y=allspp$log.LMA)
points(y=m$estimate,x=5, pch=16, cex=1.3, col="darkgrey")
arrows(x0=5,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3, col="darkgrey")
m <- cor.test(x=fam.dataclean$log.LL, y=fam.dataclean$log.LMA)
points(y=m$estimate,x=4, pch=17, cex=1.3, col="black")
arrows(x0=4,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3)


axis(1,labels = c("w/in Spp","w/in Gen", "w/in Fam","btw Fam","global"), at=c(1,2,3,4,5), las=3)
abline(h=0, lty=2)
text(y=par()$usr[4]-.2, x=c(1,2,3), labels = p$n[1:3])
text(y=par()$usr[4]-.2, x=c(4,5), labels=paste0("(",all.results.cwm.cl$n_LMA.LL[c((nrow(all.results.cwm.cl)-1),nrow(all.results.cwm.cl))],")"),cex=.8)

par(mar=c(4,4,1.5,1.5), mgp=c(2.5,1,0))
plot(log.LL~log.LMA, allspp, col="grey", pch=16, ylab="Leaf Lifespan (months)", xlab="Leaf Mass per Area (g/m2)", xaxt="n", yaxt="n")
abline(a=all.results.cwm.cl$Int_LMA.LL[nrow(all.results.cwm.cl)], b=all.results.cwm.cl$Slope_LMA.LL[nrow(all.results.cwm.cl)], lwd=3, col="black", lty=3)
#points(log.LL~log.LMA, allspp[which(allspp$Species=="Solanum straminifolia"),])
#points(log.LL~log.LMA, allspp[which(allspp$Species=="Juniperus monosperma"),])
#points(log.LL~log.LMA, allspp[which(allspp$Species=="Araucaria araucana"),])
#points(log.LL~log.LMA, allspp[which(allspp$Species=="Abies magnifica"),])
axis(2, at=log(c(1,2,3,4,5,10,20,30,40,50,100,200,300), base=10), labels = F)#c(1,"","","",5,"","","","",50,"",200,""),)
axis(2, at=log(c(1,5,50,300), base=10), labels = c(1,5,50,300))
axis(1, at=log(c(10,20,30,40,50,100,200,300,400,500,1000), base=10), labels = F)# c(10,"","","",50,"","","","","","1000"))
axis(1, at=log(c(20,50,500), base=10), labels =  c(20,50,500))



#### Global Spp means
quartz(width=7.008/2, height=4)
par(mar=c(4,4,1.5,1.5))
## scatterplot
####### Plotting the LMA vs Narea scaling in unit rather than log space ######
plot(log.LL~log.LMA, allspp, col="grey", pch=16, ylab="Leaf Lifespan (months)", xlab="Leaf Mass per Area (g/m2)", xaxt="n", yaxt="n")
abline(a=all.results.cwm.cl$Int_LMA.LL[nrow(all.results.cwm.cl)], b=all.results.cwm.cl$Slope_LMA.LL[nrow(all.results.cwm.cl)], lwd=3, col="black", lty=3)
points(log.LL~log.LMA, allspp[which(allspp$Species=="Solanum straminifolia"),])
#points(log.LL~log.LMA, allspp[which(allspp$Species=="Juniperus monosperma"),])
#points(log.LL~log.LMA, allspp[which(allspp$Species=="Araucaria araucana"),])
points(log.LL~log.LMA, allspp[which(allspp$Species=="Abies magnifica"),])
axis(2, at=log(c(1,2,3,4,5,10,20,30,40,50,100,200,300), base=10), labels = F)#c(1,"","","",5,"","","","",50,"",200,""),)
axis(2, at=log(c(1,5,50,300), base=10), labels = c(1,5,50,300),)
axis(1, at=log(c(10,20,30,40,50,100,200,300,400,500,1000), base=10), labels = F)# c(10,"","","",50,"","","","","","1000"))
axis(1, at=log(c(20,50,500), base=10), labels =  c(20,50,500))

##### Across Families ########
quartz(width=7.008/2, height=4)
par(mar=c(4,4,1.5,1.5))
## scatterplot
####### Plotting the LMA vs Narea scaling in unit rather than log space ######
plot(log.LL~log.LMA, allspp, col="grey", pch=16, ylab="log(LL)", xlab="log(LMA)")
points(log.LL~log.LMA, fam.dataclean, pch=16)
abline(a=all.results.cwm.cl$Int_LMA.LL[nrow(all.results.cwm.cl)], b=all.results.cwm.cl$Slope_LMA.LL[nrow(all.results.cwm.cl)], lwd=3, col="black", lty=3)
abline(a=all.results.cwm.cl$Int_LMA.LL[nrow(all.results.cwm.cl)-1], b=all.results.cwm.cl$Slope_LMA.LL[nrow(all.results.cwm.cl)-1], lwd=3, col="black")


##### Within Families ########
quartz(width=7.008/2, height=4)
par(mar=c(4,4,1.5,1.5))
## scatterplot
# need 14 colors, cause I've got 14 genera
genpal <- c(brewer.pal(12,"Set3"), brewer.pal(3, "Set1"))
genpal <- c(rev(brewer.pal(9,"Set1")), brewer.pal(3, "Set2"))
#genpal <- c(brewer.pal(8,"Dark2"), brewer.pal(3, "Set1"))
plot(log.LL~log.LMA, allspp, col="grey", pch=16, ylab="log(LL)", xlab="log(LMA)")
abline(a=all.results.cwm.cl$Int_LMA.LL[nrow(all.results.cwm.cl)], b=all.results.cwm.cl$Slope_LMA.LL[nrow(all.results.cwm.cl)], lwd=3, col="black", lty=3)
abline(a=all.results.cwm.cl$Int_LMA.LL[nrow(all.results.cwm.cl)-1], b=all.results.cwm.cl$Slope_LMA.LL[nrow(all.results.cwm.cl)-1], lwd=3, col="black")
palette(paste0(genpal,"88"))
points(log.LL~log.LMA, geninfam.data[which(geninfam.data$Family %in% all.results.cwm.cl$Taxo.Unit[which(all.results.cwm.cl$n_LMA.LL>5 & all.results.cwm.cl$Type=="Genw.inFam")]),], pch=16, col=factor(Family))
palette(genpal)
tax <- "Genw.inFam"
for (j in 1:length(as.character(all.results.cwm.cl$Taxo.Unit[which(all.results.cwm.cl$Type==tax & all.results.cwm.cl$n_LMA.LL>5)]))){
  i <- as.character(all.results.cwm.cl$Taxo.Unit[which(all.results.cwm.cl$Type==tax & all.results.cwm.cl$n_LMA.LL>5)])[j]
  # different colored lines
  # plot.MAR(xvar = "log.LMA", yvar = "log.LL",data= geninfam.data[which(geninfam.data$Family==i),], linecol = genpal[j], lwd=3)
  # same colored lines
  plot.MAR(xvar = "log.LMA", yvar = "log.LL",data= geninfam.data[which(geninfam.data$Family==i),], linecol = mypal[colchoices[3]], lwd=1)
}



##### Within Genera ########
quartz(width=7.008/2, height=4)
par(mar=c(4,4,1.5,1.5))
## scatterplot
# need 14 colors, cause I've got 14 genera
genpal <- c(brewer.pal(12,"Set3"), brewer.pal(3, "Set1"))
genpal <- c(rev(brewer.pal(9,"Set1")), brewer.pal(3, "Set2"))
genpal <- c(brewer.pal(8,"Dark2"), brewer.pal(3, "Set1"))
plot(log.LL~log.LMA, allspp, col="grey", pch=16, ylab="log(LL)", xlab="log(LMA)")
abline(a=all.results.cwm.cl$Int_LMA.LL[nrow(all.results.cwm.cl)], b=all.results.cwm.cl$Slope_LMA.LL[nrow(all.results.cwm.cl)], lwd=3, col="black", lty=3)
abline(a=all.results.cwm.cl$Int_LMA.LL[nrow(all.results.cwm.cl)-1], b=all.results.cwm.cl$Slope_LMA.LL[nrow(all.results.cwm.cl)-1], lwd=3, col="black")
palette(paste0(genpal,"88"))
points(log.LL~log.LMA, gen.data[which(gen.data$Genus %in% all.results.cwm.cl$Taxo.Unit[which(all.results.cwm.cl$n_LMA.LL>5 & all.results.cwm.cl$Type=="w.inGen")]),], pch=16)#, col=factor(Genus))
palette(genpal)
tax <- "w.inGen"
for (j in 1:length(as.character(all.results.cwm.cl$Taxo.Unit[which(all.results.cwm.cl$Type==tax & all.results.cwm.cl$n_LMA.LL>5)]))){
  i <- as.character(all.results.cwm.cl$Taxo.Unit[which(all.results.cwm.cl$Type==tax & all.results.cwm.cl$n_LMA.LL>5)])[j]
  # plot.MAR(xvar = "log.LMA", yvar = "log.LL",data= gen.data[which(gen.data$Genus==i),], linecol = genpal[j], lwd=3)
  plot.MAR(xvar = "log.LMA", yvar = "log.LL",data= gen.data[which(gen.data$Genus==i),], linecol = mypal[colchoices[2]], lwd=1)
  
}







##### Within Species ########
quartz(width=7.008/2, height=4)
par(mar=c(4,4,1.5,1.5))
## scatterplot
# need 18 colors, cause I've got 18 Spp
genpal <- c(brewer.pal(12,"Set3"), brewer.pal(3, "Set1"))
genpal <- c(rev(brewer.pal(9,"Set1")), brewer.pal(3, "Set2"))
genpal <- c(brewer.pal(8,"Dark2"), brewer.pal(8, "Set1"), brewer.pal(7,"Dark2")[c(1,7)])
plot(log.LL~log.LMA, allspp, col="grey", pch=16, ylab="log(LL)", xlab="log(LMA)")
abline(a=all.results.cwm.cl$Int_LMA.LL[nrow(all.results.cwm.cl)], b=all.results.cwm.cl$Slope_LMA.LL[nrow(all.results.cwm.cl)], lwd=3, col="black", lty=3)
abline(a=all.results.cwm.cl$Int_LMA.LL[nrow(all.results.cwm.cl)-1], b=all.results.cwm.cl$Slope_LMA.LL[nrow(all.results.cwm.cl)-1], lwd=3, col="black")
palette(paste0(genpal,"88"))
points(log.LL~log.LMA, spp.data[which(spp.data$Species %in% all.results.cwm.cl$Taxo.Unit[which(all.results.cwm.cl$n_LMA.LL>5 & all.results.cwm.cl$Type=="w.inSpp")]),], pch=16)#, col=factor(Species))
palette(genpal)
tax <- "w.inSpp"
for (j in 1:length(as.character(all.results.cwm.cl$Taxo.Unit[which(all.results.cwm.cl$Type==tax & all.results.cwm.cl$n_LMA.LL>5)]))){
  i <- as.character(all.results.cwm.cl$Taxo.Unit[which(all.results.cwm.cl$Type==tax & all.results.cwm.cl$n_LMA.LL>5)])[j]
  plot.MAR(xvar = "log.LMA", yvar = "log.LL",data= spp.data[which(spp.data$Species==i),], linecol = genpal[j], lwd=3)
  #plot.MAR(xvar = "log.LMA", yvar = "log.LL",data= spp.data[which(spp.data$Species==i),], linecol = mypal[colchoices[1]], lwd=3)
  
}










#________________________________________________________________
###### Full Global only ######
#________________________________________________________________

mat <- matrix(c(1,3,
                2,3), nrow=2, byrow = T)
colchoices <- c(1,2,4,3,6)
palette(mypal[colchoices])

quartz(width=7.008, height=4)
layout(mat)
par(mar=c(0,4,6,1))
p <- boxplot(Slope_LMA.LL~Type, all.results.cwm.cl[which(all.results.cwm.cl$n_LMA.LL>5 & !all.results.cwm.cl$Type %in% c("Famclean","global","w.inSpp","w.inGen","Genw.inFam")),]
             , ylim=c(-2.5,4),las=3, ylab="Slope" #, main="log(LMA)~log(Nmass) strict"
             , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
             ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
             , staplelwd=0, outpch=16, outcex=.5, outcol=mypal[colchoices]
             , boxwex=.7, xaxt="n")
m <- lmodel2(log.LL~log.LMA, allspp)
points(y=m$regression.results$Slope[3],x=5, pch=16, cex=1.3, col="darkgrey")
arrows(x0=5,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`2.5%-Slope`[3], length = 0,lwd=3, col="darkgrey")
arrows(x0=5,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`97.5%-Slope`[3], length = 0,lwd=3, col="darkgrey")
# m <- lmodel2(log.LL~log.LMA, fam.dataclean)
# points(y=m$regression.results$Slope[3],x=4, pch=17, cex=1.3, col="black")
# arrows(x0=4,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`2.5%-Slope`[3], length = 0,lwd=3)
# arrows(x0=4,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`97.5%-Slope`[3], length = 0,lwd=3)


abline(h=0, lty=2)
#text(y=par()$usr[3]+.2, x=c(1,2,3,4,5,6), labels = p$n)
#text(y=par()$usr[4]-.2,x=.5, labels = "a)")
mtext(text = "a)", side = 3, adj=0, line=.2)
mtext(text= "LMA vs LL", side=3, line=.2)
par(mar=c(6,4,0,1))
p <- boxplot(Rho_LMA.LL~Type, all.results.cwm.cl[which(all.results.cwm.cl$n_LMA.LL>5& !all.results.cwm.cl$Type %in% c("Famclean","global","w.inSpp","w.inGen","Genw.inFam")),]
             , ylim=c(-1,1.4),las=3, ylab="Correlation"
             , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
             ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
             , staplelwd=0, outpch=16, outcex=.5, outcol=mypal[colchoices]
             , boxwex=.7, xaxt="n")
m <- cor.test(x=allspp$log.LL, y=allspp$log.LMA)
points(y=m$estimate,x=5, pch=16, cex=1.3, col="darkgrey")
arrows(x0=5,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3, col="darkgrey")
# m <- cor.test(x=fam.dataclean$log.LL, y=fam.dataclean$log.LMA)
# points(y=m$estimate,x=4, pch=17, cex=1.3, col="black")
# arrows(x0=4,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3)


axis(1,labels = c("w/in Spp","w/in Gen", "w/in Fam","btw Fam","global"), at=c(1,2,3,4,5), las=3)
abline(h=0, lty=2)
# text(y=par()$usr[4]-.2, x=c(1,2,3), labels = p$n[1:3])
text(y=par()$usr[4]-.2, x=c(5), labels=paste0("(",all.results.cwm.cl$n_LMA.LL[nrow(all.results.cwm.cl)],")"),cex=.8)

par(mar=c(4,4,1.5,1.5), mgp=c(2.5,1,0))
plot(log.LL~log.LMA, allspp, col="grey", pch=16, ylab="Leaf Lifespan (months)", xlab="Leaf Mass per Area (g/m2)", xaxt="n", yaxt="n")
abline(a=all.results.cwm.cl$Int_LMA.LL[nrow(all.results.cwm.cl)], b=all.results.cwm.cl$Slope_LMA.LL[nrow(all.results.cwm.cl)], lwd=3, col="black", lty=3)
#points(log.LL~log.LMA, allspp[which(allspp$Species=="Solanum straminifolia"),])
#points(log.LL~log.LMA, allspp[which(allspp$Species=="Juniperus monosperma"),])
#points(log.LL~log.LMA, allspp[which(allspp$Species=="Araucaria araucana"),])
#points(log.LL~log.LMA, allspp[which(allspp$Species=="Abies magnifica"),])
axis(2, at=log(c(1,2,3,4,5,10,20,30,40,50,100,200,300), base=10), labels = F)#c(1,"","","",5,"","","","",50,"",200,""),)
axis(2, at=log(c(1,5,50,300), base=10), labels = c(1,5,50,300))
axis(1, at=log(c(10,20,30,40,50,100,200,300,400,500,1000), base=10), labels = F)# c(10,"","","",50,"","","","","","1000"))
axis(1, at=log(c(20,50,500), base=10), labels =  c(20,50,500))





#________________________________________________________________
###### Full Across Families only ######
#________________________________________________________________

tmpdata <- all.results.cwm.cl[-which(all.results.cwm.cl$Type=="global"),]
tmpdata$Type <- factor(tmpdata$Type) # remove global level for boxplots
mat <- matrix(c(1,3,
                2,3), nrow=2, byrow = T)
colchoices <- c(1,2,4,3,6)
palette(mypal[colchoices])

quartz(width=7.008, height=4)
layout(mat, heights=c(1,1.5))
par(mar=c(0,4,1.5,1), mgp=c(2.6,1,0))
p <- boxplot(Slope_LMA.LL~Type, tmpdata[which(tmpdata$n_LMA.LL>5 & !tmpdata$Type %in% c("Famclean","global","w.inSpp","w.inGen","Genw.inFam","CWM")),]
             , ylim=c(-2.5,4),las=3, ylab="Slope" #, main="log(LMA)~log(Nmass) strict"
             , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
             ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
             , staplelwd=0, outpch=16, outcex=.5, outcol=mypal[colchoices]
             , boxwex=.7, xaxt="n", cex.lab=1.3)
# m <- lmodel2(log.LL~log.LMA, allspp)
# points(y=m$regression.results$Slope[3],x=5, pch=16, cex=1.3, col="darkgrey")
# arrows(x0=5,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`2.5%-Slope`[3], length = 0,lwd=3, col="darkgrey")
# arrows(x0=5,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`97.5%-Slope`[3], length = 0,lwd=3, col="darkgrey")
m <- lmodel2(log.LL~log.LMA, fam.dataclean)
points(y=m$regression.results$Slope[3],x=4, pch=16, cex=1.3, col="black")
arrows(x0=4,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`2.5%-Slope`[3], length = 0,lwd=3)
arrows(x0=4,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`97.5%-Slope`[3], length = 0,lwd=3)


abline(h=0, lty=2)
#text(y=par()$usr[3]+.2, x=c(1,2,3,4,5,6), labels = p$n)
#text(y=par()$usr[4]-.2,x=.5, labels = "a)")
#mtext(text = "a)", side = 3, adj=0, line=.2)
mtext(text= "LMA vs LL", side=3, line=.2)
par(mar=c(6,4,0,1))
p <- boxplot(Rho_LMA.LL~Type, tmpdata[which(tmpdata$n_LMA.LL>5& !tmpdata$Type %in% c("Famclean","global","w.inSpp","w.inGen","Genw.inFam","CWM")),]
             , ylim=c(-1,1.4),las=3, ylab="Correlation"
             , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
             ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
             , staplelwd=0, outpch=16, outcex=.5, outcol=mypal[colchoices]
             , boxwex=.7, xaxt="n", cex.lab=1.3)
# m <- cor.test(x=allspp$log.LL, y=allspp$log.LMA)
# points(y=m$estimate,x=5, pch=16, cex=1.3, col="darkgrey")
# arrows(x0=5,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3, col="darkgrey")
m <- cor.test(x=fam.dataclean$log.LL, y=fam.dataclean$log.LMA)
points(y=m$estimate,x=4, pch=16, cex=1.3, col="black")
arrows(x0=4,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3)


axis(1,labels = c("w/in Spp","w/in Gen", "w/in Fam","btw Fam","PNW CWMs"), at=c(1,2,3,4,5), las=3)
abline(h=0, lty=2)
# text(y=par()$usr[4]-.2, x=c(1,2,3), labels = p$n[1:3])
#text(y=par()$usr[4]-.2, x=c(4,5), labels=paste0("(",all.results.cwm.cl$n_LMA.LL[c((nrow(all.results.cwm.cl)-1),nrow(all.results.cwm.cl))],")"),cex=.8)
#text(y=par()$usr[4]-.2, x=c(1,2,3), labels = p$n[1:3])
text(y=par()$usr[4]-.2, x=c(4), labels=paste0("(",all.results.cwm.cl$n_LMA.LL[c((nrow(all.results.cwm.cl)-2))],")"),cex=.9)

par(mar=c(4,4,1.5,1.5), mgp=c(2.5,1,0))
plot(log.LL~log.LMA, allspp,type="n", col="grey", pch=16, ylab="Leaf Lifespan (months)", xlab="Leaf Mass per Area (g/m2)", xaxt="n", yaxt="n")
#abline(a=all.results.cwm.cl$Int_LMA.LL[nrow(all.results.cwm.cl)], b=all.results.cwm.cl$Slope_LMA.LL[nrow(all.results.cwm.cl)], lwd=3, col="black", lty=3)

#points(log.LL~log.LMA, allspp[which(allspp$Species=="Solanum straminifolia"),])
#points(log.LL~log.LMA, allspp[which(allspp$Species=="Juniperus monosperma"),])
#points(log.LL~log.LMA, allspp[which(allspp$Species=="Araucaria araucana"),])
#points(log.LL~log.LMA, allspp[which(allspp$Species=="Abies magnifica"),])
axis(2, at=log(c(1,2,3,4,5,10,20,30,40,50,100,200,300), base=10), labels = F)#c(1,"","","",5,"","","","",50,"",200,""),)
axis(2, at=log(c(1,5,50,300), base=10), labels = c(1,5,50,300))
axis(1, at=log(c(10,20,30,40,50,100,200,300,400,500,1000), base=10), labels = F)# c(10,"","","",50,"","","","","","1000"))
axis(1, at=log(c(20,50,500), base=10), labels =  c(20,50,500))
points(log.LL~log.LMA, fam.dataclean, pch=16)
abline(a=all.results.cwm.cl$Int_LMA.LL[nrow(all.results.cwm.cl)-1], b=all.results.cwm.cl$Slope_LMA.LL[nrow(all.results.cwm.cl)-1], lwd=3, col="black")









#________________________________________________________________
###### Full Within Families only ######
#________________________________________________________________

mat <- matrix(c(1,3,
                2,3), nrow=2, byrow = T)
colchoices <- c(1,2,4,3,6)
palette(mypal[colchoices])

quartz(width=7.008, height=4)
layout(mat, heights=c(1,1.5))
par(mar=c(0,4,1.5,1), mgp=c(2.6,1,0))
p <- boxplot(Slope_LMA.LL~Type, tmpdata[which(tmpdata$n_LMA.LL>5 & !tmpdata$Type %in% c("Famclean","global","w.inSpp","w.inGen","CWM")),]
             , ylim=c(-2.5,4),las=3, ylab="Slope" #, main="log(LMA)~log(Nmass) strict"
             , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
             ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
             , staplelwd=0, outpch=16, outcex=.5, outcol=mypal[colchoices]
             , boxwex=.7, xaxt="n", cex.lab=1.3)
# m <- lmodel2(log.LL~log.LMA, allspp)
# points(y=m$regression.results$Slope[3],x=5, pch=16, cex=1.3, col="darkgrey")
# arrows(x0=5,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`2.5%-Slope`[3], length = 0,lwd=3, col="darkgrey")
# arrows(x0=5,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`97.5%-Slope`[3], length = 0,lwd=3, col="darkgrey")
m <- lmodel2(log.LL~log.LMA, fam.dataclean)
points(y=m$regression.results$Slope[3],x=4, pch=16, cex=1.3, col="black")
arrows(x0=4,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`2.5%-Slope`[3], length = 0,lwd=3)
arrows(x0=4,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`97.5%-Slope`[3], length = 0,lwd=3)


# tmp <- all.results.cwm.cl[which(all.results.cwm.cl$n_LMA.LL>5 & !is.na(all.results.cwm.cl$sig_LMA.LL) & !all.results.cwm.cl$Type %in% c("Famclean","global","w.inSpp","w.inGen","CWM")),]
# tmp$rho.sig <- rho.sig(rho=tmp$Rho_LMA.LL, n= tmp$n_LMA.LL)
# tmp.ns <- tmp[which(tmp$rho.sig>0.1),]
# tmp.sig <- tmp[which(tmp$rho.sig<=0.1),]
# #palette(paste0(mypal[colchoices], "55"))
# set.seed(42)
# points(Slope_LMA.LL~jitter(as.numeric(Type), factor=3), tmp.ns, pch=1, col=Type, cex=cex.ns)
# #palette(mypal[colchoices])
# set.seed(42)
# points(Slope_LMA.LL~jitter(as.numeric(Type), factor=3), tmp.sig, pch=16, col=Type,cex=cex.sig)

abline(h=0, lty=2)
#text(y=par()$usr[3]+.2, x=c(1,2,3,4,5,6), labels = p$n)
#text(y=par()$usr[4]-.2,x=.5, labels = "a)")
#mtext(text = "a)", side = 3, adj=0, line=.2)
mtext(text= "LMA vs LL", side=3, line=.2)
par(mar=c(6,4,0,1))
p <- boxplot(Rho_LMA.LL~Type, tmpdata[which(tmpdata$n_LMA.LL>5& !tmpdata$Type %in% c("Famclean","global","w.inSpp","w.inGen","CWM")),]
             , ylim=c(-1,1.4),las=3, ylab="Correlation"
             , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
             ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
             , staplelwd=0, outpch=16, outcex=.5, outcol=mypal[colchoices]
             , boxwex=.7, xaxt="n", cex.lab=1.3)
# m <- cor.test(x=allspp$log.LL, y=allspp$log.LMA)
# points(y=m$estimate,x=5, pch=16, cex=1.3, col="darkgrey")
# arrows(x0=5,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3, col="darkgrey")
m <- cor.test(x=fam.dataclean$log.LL, y=fam.dataclean$log.LMA)
points(y=m$estimate,x=4, pch=16, cex=1.3, col="black")
arrows(x0=4,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3)

# 
# set.seed(62)
# points(Rho_LMA.LL~jitter(as.numeric(Type)), all.results.cwm.cl[which(all.results.cwm.cl$n_LMA.LL>5 & all.results.cwm.cl$sig_LMA.LL<=crit),], pch=16, col=Type, cex=cex.sig)
# #palette(paste0(mypal[colchoices], "55"))
# set.seed(42)
# points(Rho_LMA.LL~jitter(as.numeric(Type)), all.results.cwm.cl[which(all.results.cwm.cl$n_LMA.LL>5 & all.results.cwm.cl$sig_LMA.LL>crit),], pch=1, col=Type ,cex=cex.ns)
# 

axis(1,labels = c("w/in Spp","w/in Gen", "w/in Fam","btw Fam","PNW CWMs"), at=c(1,2,3,4,5), las=3)
abline(h=0, lty=2)
# text(y=par()$usr[4]-.2, x=c(1,2,3), labels = p$n[1:3])
#text(y=par()$usr[4]-.2, x=c(4,5), labels=paste0("(",all.results.cwm.cl$n_LMA.LL[c((nrow(all.results.cwm.cl)-1),nrow(all.results.cwm.cl))],")"),cex=.8)
text(y=par()$usr[4]-.2, x=c(3), labels = p$n[3])
text(y=par()$usr[4]-.2, x=c(4), labels=paste0("(",all.results.cwm.cl$n_LMA.LL[c((nrow(all.results.cwm.cl)-2))],")"),cex=.9)

par(mar=c(4,4,1.5,1.5), mgp=c(2.5,1,0))
plot(log.LL~log.LMA, allspp, typ="n",col="grey", pch=16, ylab="Leaf Lifespan (months)", xlab="Leaf Mass per Area (g/m2)", xaxt="n", yaxt="n")
# abline(a=tmpdata$Int_LMA.LL[nrow(tmpdata)], b=tmpdata$Slope_LMA.LL[nrow(tmpdata)], lwd=3, col="black", lty=3)

axis(2, at=log(c(1,2,3,4,5,10,20,30,40,50,100,200,300), base=10), labels = F)#c(1,"","","",5,"","","","",50,"",200,""),)
axis(2, at=log(c(1,5,50,300), base=10), labels = c(1,5,50,300))
axis(1, at=log(c(10,20,30,40,50,100,200,300,400,500,1000), base=10), labels = F)# c(10,"","","",50,"","","","","","1000"))
axis(1, at=log(c(20,50,500), base=10), labels =  c(20,50,500))
points(log.LL~log.LMA, geninfam.data, pch=16, col="grey")
# abline(a=tmpdata$Int_LMA.LL[nrow(tmpdata)-1], b=tmpdata$Slope_LMA.LL[nrow(tmpdata)-1], lwd=3, col="black")
abline(a=tmpdata$Int_LMA.LL[nrow(tmpdata)], b=tmpdata$Slope_LMA.LL[nrow(tmpdata)], lwd=3, col="black")
tax <- "Genw.inFam"
for (j in 1:length(as.character(tmpdata$Taxo.Unit[which(tmpdata$Type==tax & tmpdata$n_LMA.LL>5)]))){
  i <- as.character(tmpdata$Taxo.Unit[which(tmpdata$Type==tax & tmpdata$n_LMA.LL>5)])[j]
  # different colored lines
  # plot.MAR(xvar = "log.LMA", yvar = "log.LL",data= geninfam.data[which(geninfam.data$Family==i),], linecol = genpal[j], lwd=3)
  # same colored lines
  plot.MAR(xvar = "log.LMA", yvar = "log.LL",data= geninfam.data[which(geninfam.data$Family==i),], linecol = mypal[colchoices[3]], lwd=3)
}




#________________________________________________________________
###### Full Within Genera only ######
#________________________________________________________________

mat <- matrix(c(1,3,
                2,3), nrow=2, byrow = T)
colchoices <- c(1,2,4,3,6)
palette(mypal[colchoices])

quartz(width=7.008, height=4)
layout(mat, heights=c(1,1.5))
par(mar=c(0,4,1.5,1), mgp=c(2.6,1,0))
p <- boxplot(Slope_LMA.LL~Type, tmpdata[which(tmpdata$n_LMA.LL>5 & !tmpdata$Type %in% c("Famclean","global","w.inSpp","CWM")),]
             , ylim=c(-2.5,4),las=3, ylab="Slope" #, main="log(LMA)~log(Nmass) strict"
             , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
             ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
             , staplelwd=0, outpch=16, outcex=.5, outcol=mypal[colchoices]
             , boxwex=.7, xaxt="n", cex.lab=1.3)
# m <- lmodel2(log.LL~log.LMA, allspp)
# points(y=m$regression.results$Slope[3],x=5, pch=16, cex=1.3, col="darkgrey")
# arrows(x0=5,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`2.5%-Slope`[3], length = 0,lwd=3, col="darkgrey")
# arrows(x0=5,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`97.5%-Slope`[3], length = 0,lwd=3, col="darkgrey")
m <- lmodel2(log.LL~log.LMA, fam.dataclean)
points(y=m$regression.results$Slope[3],x=4, pch=16, cex=1.3, col="black")
arrows(x0=4,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`2.5%-Slope`[3], length = 0,lwd=3)
arrows(x0=4,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`97.5%-Slope`[3], length = 0,lwd=3)


abline(h=0, lty=2)
#text(y=par()$usr[3]+.2, x=c(1,2,3,4,5,6), labels = p$n)
#text(y=par()$usr[4]-.2,x=.5, labels = "a)")
#mtext(text = "a)", side = 3, adj=0, line=.2)
mtext(text= "LMA vs LL", side=3, line=.2)
par(mar=c(6,4,0,1))
p <- boxplot(Rho_LMA.LL~Type, tmpdata[which(tmpdata$n_LMA.LL>5& !tmpdata$Type %in% c("Famclean","global","w.inSpp","CWM")),]
             , ylim=c(-1,1.4),las=3, ylab="Correlation"
             , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
             ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
             , staplelwd=0, outpch=16, outcex=.5, outcol=mypal[colchoices]
             , boxwex=.7, xaxt="n", cex.lab=1.3)
# m <- cor.test(x=allspp$log.LL, y=allspp$log.LMA)
# points(y=m$estimate,x=5, pch=16, cex=1.3, col="darkgrey")
# arrows(x0=5,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3, col="darkgrey")
m <- cor.test(x=fam.dataclean$log.LL, y=fam.dataclean$log.LMA)
points(y=m$estimate,x=4, pch=16, cex=1.3, col="black")
arrows(x0=4,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3)


axis(1,labels = c("w/in Spp","w/in Gen", "w/in Fam","btw Fam","PNW CWMs"), at=c(1,2,3,4,5), las=3)
abline(h=0, lty=2)
# text(y=par()$usr[4]-.2, x=c(1,2,3), labels = p$n[1:3])
#text(y=par()$usr[4]-.2, x=c(4,5), labels=paste0("(",all.results.cwm.cl$n_LMA.LL[c((nrow(all.results.cwm.cl)-1),nrow(all.results.cwm.cl))],")"),cex=.8)
text(y=par()$usr[4]-.2, x=c(2,3), labels = p$n[c(2,3)])
text(y=par()$usr[4]-.2, x=c(4), labels=paste0("(",all.results.cwm.cl$n_LMA.LL[c((nrow(all.results.cwm.cl)-2))],")"),cex=.9)



par(mar=c(4,4,1.5,1.5), mgp=c(2.5,1,0))
plot(log.LL~log.LMA, allspp, typ="n",col="grey", pch=16, ylab="Leaf Lifespan (months)", xlab="Leaf Mass per Area (g/m2)", xaxt="n", yaxt="n")
abline(a=tmpdata$Int_LMA.LL[nrow(tmpdata)], b=tmpdata$Slope_LMA.LL[nrow(tmpdata)], lwd=3, col="black", lty=3)
#points(log.LL~log.LMA, allspp[which(allspp$Species=="Solanum straminifolia"),])
#points(log.LL~log.LMA, allspp[which(allspp$Species=="Juniperus monosperma"),])
#points(log.LL~log.LMA, allspp[which(allspp$Species=="Araucaria araucana"),])
#points(log.LL~log.LMA, allspp[which(allspp$Species=="Abies magnifica"),])
axis(2, at=log(c(1,2,3,4,5,10,20,30,40,50,100,200,300), base=10), labels = F)#c(1,"","","",5,"","","","",50,"",200,""),)
axis(2, at=log(c(1,5,50,300), base=10), labels = c(1,5,50,300))
axis(1, at=log(c(10,20,30,40,50,100,200,300,400,500,1000), base=10), labels = F)# c(10,"","","",50,"","","","","","1000"))
axis(1, at=log(c(20,50,500), base=10), labels =  c(20,50,500))
points(log.LL~log.LMA, gen.data, pch=16, col="grey")
# abline(a=tmpdata$Int_LMA.LL[nrow(tmpdata)-1], b=tmpdata$Slope_LMA.LL[nrow(tmpdata)-1], lwd=3, col="black")
abline(a=tmpdata$Int_LMA.LL[nrow(tmpdata)], b=tmpdata$Slope_LMA.LL[nrow(tmpdata)], lwd=3, col="black")
tax <- "w.inGen"
for (j in 1:length(as.character(tmpdata$Taxo.Unit[which(tmpdata$Type==tax & tmpdata$n_LMA.LL>5)]))){
  i <- as.character(tmpdata$Taxo.Unit[which(tmpdata$Type==tax & tmpdata$n_LMA.LL>5)])[j]
  # different colored lines
  # plot.MAR(xvar = "log.LMA", yvar = "log.LL",data= geninfam.data[which(geninfam.data$Family==i),], linecol = genpal[j], lwd=3)
  # same colored lines
  plot.MAR(xvar = "log.LMA", yvar = "log.LL",data= gen.data[which(gen.data$Genus==i),], linecol = mypal[colchoices[2]], lwd=3)
}




#________________________________________________________________
###### Full Within Species only ######
#________________________________________________________________

mat <- matrix(c(1,3,
                2,3), nrow=2, byrow = T)
colchoices <- c(1,2,4,3,6)
palette(mypal[colchoices])


quartz(width=7.008, height=4)
layout(mat, heights=c(1,1.5))
par(mar=c(0,4,1.5,1), mgp=c(2.6,1,0))
p <- boxplot(Slope_LMA.LL~Type, tmpdata[which(tmpdata$n_LMA.LL>5 & !tmpdata$Type %in% c("Famclean","global","CWM")),]
             , ylim=c(-2.5,4),las=3, ylab="Slope" #, main="log(LMA)~log(Nmass) strict"
             , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
             ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
             , staplelwd=0, outpch=16, outcex=.5, outcol=mypal[colchoices]
             , boxwex=.7, xaxt="n", cex.lab=1.3)
# m <- lmodel2(log.LL~log.LMA, allspp)
# points(y=m$regression.results$Slope[3],x=5, pch=16, cex=1.3, col="darkgrey")
# arrows(x0=5,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`2.5%-Slope`[3], length = 0,lwd=3, col="darkgrey")
# arrows(x0=5,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`97.5%-Slope`[3], length = 0,lwd=3, col="darkgrey")
m <- lmodel2(log.LL~log.LMA, fam.dataclean)
points(y=m$regression.results$Slope[3],x=4, pch=16, cex=1.3, col="black")
arrows(x0=4,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`2.5%-Slope`[3], length = 0,lwd=3)
arrows(x0=4,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`97.5%-Slope`[3], length = 0,lwd=3)


abline(h=0, lty=2)
#text(y=par()$usr[3]+.2, x=c(1,2,3,4,5,6), labels = p$n)
#text(y=par()$usr[4]-.2,x=.5, labels = "a)")
#mtext(text = "a)", side = 3, adj=0, line=.2)
mtext(text= "LMA vs LL", side=3, line=.2)
par(mar=c(6,4,0,1))
p <- boxplot(Rho_LMA.LL~Type, tmpdata[which(tmpdata$n_LMA.LL>5& !tmpdata$Type %in% c("Famclean","global","CWM")),]
             , ylim=c(-1,1.4),las=3, ylab="Correlation"
             , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
             ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
             , staplelwd=0, outpch=16, outcex=.5, outcol=mypal[colchoices]
             , boxwex=.7, xaxt="n", cex.lab=1.3)
# m <- cor.test(x=allspp$log.LL, y=allspp$log.LMA)
# points(y=m$estimate,x=5, pch=16, cex=1.3, col="darkgrey")
# arrows(x0=5,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3, col="darkgrey")
m <- cor.test(x=fam.dataclean$log.LL, y=fam.dataclean$log.LMA)
points(y=m$estimate,x=4, pch=16, cex=1.3, col="black")
arrows(x0=4,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3)


axis(1,labels = c("w/in Spp","w/in Gen", "w/in Fam","btw Fam","PNW CWMs"), at=c(1,2,3,4,5), las=3)
abline(h=0, lty=2)
# text(y=par()$usr[4]-.2, x=c(1,2,3), labels = p$n[1:3])
#text(y=par()$usr[4]-.2, x=c(4,5), labels=paste0("(",all.results.cwm.cl$n_LMA.LL[c((nrow(all.results.cwm.cl)-1),nrow(all.results.cwm.cl))],")"),cex=.8)
text(y=par()$usr[4]-.2, x=c(1,2,3), labels = p$n[c(1,2,3)])
text(y=par()$usr[4]-.2, x=c(4), labels=paste0("(",all.results.cwm.cl$n_LMA.LL[c((nrow(all.results.cwm.cl)-2))],")"),cex=.9)



par(mar=c(4,4,1.5,1.5), mgp=c(2.5,1,0))
plot(log.LL~log.LMA, allspp, typ="n",col="grey", pch=16, ylab="Leaf Lifespan (months)", xlab="Leaf Mass per Area (g/m2)", xaxt="n", yaxt="n")
abline(a=tmpdata$Int_LMA.LL[nrow(tmpdata)], b=tmpdata$Slope_LMA.LL[nrow(tmpdata)], lwd=3, col="black", lty=3)
#points(log.LL~log.LMA, allspp[which(allspp$Species=="Solanum straminifolia"),])
#points(log.LL~log.LMA, allspp[which(allspp$Species=="Juniperus monosperma"),])
#points(log.LL~log.LMA, allspp[which(allspp$Species=="Araucaria araucana"),])
#points(log.LL~log.LMA, allspp[which(allspp$Species=="Abies magnifica"),])
axis(2, at=log(c(1,2,3,4,5,10,20,30,40,50,100,200,300), base=10), labels = F)#c(1,"","","",5,"","","","",50,"",200,""),)
axis(2, at=log(c(1,5,50,300), base=10), labels = c(1,5,50,300))
axis(1, at=log(c(10,20,30,40,50,100,200,300,400,500,1000), base=10), labels = F)# c(10,"","","",50,"","","","","","1000"))
axis(1, at=log(c(20,50,500), base=10), labels =  c(20,50,500))
points(log.LL~log.LMA, spp.data, pch=16, col="grey")
# abline(a=tmpdata$Int_LMA.LL[nrow(tmpdata)-1], b=tmpdata$Slope_LMA.LL[nrow(tmpdata)-1], lwd=3, col="black")
abline(a=tmpdata$Int_LMA.LL[nrow(tmpdata)], b=tmpdata$Slope_LMA.LL[nrow(tmpdata)], lwd=3, col="black")

# tax <- "Genw.inFam"
# for (j in 1:length(as.character(tmpdata$Taxo.Unit[which(tmpdata$Type==tax & tmpdata$n_LMA.LL>5)]))){
#   i <- as.character(tmpdata$Taxo.Unit[which(tmpdata$Type==tax & tmpdata$n_LMA.LL>5)])[j]
#   # different colored lines
#   # plot.MAR(xvar = "log.LMA", yvar = "log.LL",data= geninfam.data[which(geninfam.data$Family==i),], linecol = genpal[j], lwd=3)
#   # same colored lines
#   plot.MAR(xvar = "log.LMA", yvar = "log.LL",data= geninfam.data[which(geninfam.data$Family==i),], linecol = mypal[colchoices[3]], lwd=1)
# }

# tax <- "w.inGen"
# for (j in 1:length(as.character(tmpdata$Taxo.Unit[which(tmpdata$Type==tax & tmpdata$n_LMA.LL>5)]))){
#   i <- as.character(tmpdata$Taxo.Unit[which(tmpdata$Type==tax & tmpdata$n_LMA.LL>5)])[j]
#   # different colored lines
#   # plot.MAR(xvar = "log.LMA", yvar = "log.LL",data= geninfam.data[which(geninfam.data$Family==i),], linecol = genpal[j], lwd=3)
#   # same colored lines
#   plot.MAR(xvar = "log.LMA", yvar = "log.LL",data= gen.data[which(gen.data$Genus==i),], linecol = mypal[colchoices[2]], lwd=1)
# }




tax <- "w.inSpp"
for (j in 1:length(as.character(tmpdata$Taxo.Unit[which(tmpdata$Type==tax & tmpdata$n_LMA.LL>5)]))){
  i <- as.character(tmpdata$Taxo.Unit[which(tmpdata$Type==tax & tmpdata$n_LMA.LL>5)])[j]
  #plot.MAR(xvar = "log.LMA", yvar = "log.LL",data= spp.data[which(spp.data$Species==i),], linecol = genpal[j], lwd=3)
  plot.MAR(xvar = "log.LMA", yvar = "log.LL",data= spp.data[which(spp.data$Species==i),], linecol = mypal[colchoices[1]], lwd=3)
  
}




#________________________________________________________________
###### Full Btw Communities only ######
#________________________________________________________________

mat <- matrix(c(1,3,
                2,3), nrow=2, byrow = T)
colchoices <- c(1,2,4,3,6)
palette(mypal[colchoices])


quartz(width=7.008, height=4)
layout(mat, heights=c(1,1.5))
par(mar=c(0,4,1.5,1), mgp=c(2.6,1,0))
p <- boxplot(Slope_LMA.LL~Type, tmpdata[which(tmpdata$n_LMA.LL>5 & !tmpdata$Type %in% c("Famclean","global","CWM")),]
             , ylim=c(-2.5,4),las=3, ylab="Slope" #, main="log(LMA)~log(Nmass) strict"
             , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
             ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
             , staplelwd=0, outpch=16, outcex=.5, outcol=mypal[colchoices]
             , boxwex=.7, xaxt="n", cex.lab=1.3)
# m <- lmodel2(log.LL~log.LMA, allspp)
# points(y=m$regression.results$Slope[3],x=5, pch=16, cex=1.3, col="darkgrey")
# arrows(x0=5,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`2.5%-Slope`[3], length = 0,lwd=3, col="darkgrey")
# arrows(x0=5,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`97.5%-Slope`[3], length = 0,lwd=3, col="darkgrey")
m <- lmodel2(log.LL~log.LMA, fam.dataclean)
points(y=m$regression.results$Slope[3],x=4, pch=16, cex=1.3, col="black")
arrows(x0=4,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`2.5%-Slope`[3], length = 0,lwd=3)
arrows(x0=4,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`97.5%-Slope`[3], length = 0,lwd=3)
m <- lmodel2(log.cw_LLp_if~log.cw_LMAp_if, biomass)
arrows(x0=5,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`2.5%-Slope`[3], length = 0,lwd=3, col=mypal[5])
arrows(x0=5,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`97.5%-Slope`[3], length = 0,lwd=3, col=mypal[5])
points(y=m$regression.results$Slope[3],x=5, pch=17, cex=1.3, col=mypal[5])


abline(h=0, lty=2)
#text(y=par()$usr[3]+.2, x=c(1,2,3,4,5,6), labels = p$n)
#text(y=par()$usr[4]-.2,x=.5, labels = "a)")
#mtext(text = "a)", side = 3, adj=0, line=.2)
mtext(text= "LMA vs LL", side=3, line=.2)
par(mar=c(6,4,0,1))
p <- boxplot(Rho_LMA.LL~Type, tmpdata[which(tmpdata$n_LMA.LL>5& !tmpdata$Type %in% c("Famclean","global","CWM")),]
             , ylim=c(-1,1.4),las=3, ylab="Correlation"
             , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
             ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
             , staplelwd=0, outpch=16, outcex=.5, outcol=mypal[colchoices]
             , boxwex=.7, xaxt="n", cex.lab=1.3)
# m <- cor.test(x=allspp$log.LL, y=allspp$log.LMA)
# points(y=m$estimate,x=5, pch=16, cex=1.3, col="darkgrey")
# arrows(x0=5,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3, col="darkgrey")
m <- cor.test(x=fam.dataclean$log.LL, y=fam.dataclean$log.LMA)
points(y=m$estimate,x=4, pch=16, cex=1.3, col="black")
arrows(x0=4,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3)
m <- cor.test(x=biomass$log.cw_LMAp_if, y=biomass$log.cw_LLp_if)
arrows(x0=5,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3, col=mypal[5])
points(y=m$estimate,x=5, pch=17, cex=1.3, col=mypal[5])

axis(1,labels = c("w/in Spp","w/in Gen", "w/in Fam","btw Fam","PNW CWMs"), at=c(1,2,3,4,5), las=3)
abline(h=0, lty=2)
# text(y=par()$usr[4]-.2, x=c(1,2,3), labels = p$n[1:3])
#text(y=par()$usr[4]-.2, x=c(4,5), labels=paste0("(",all.results.cwm.cl$n_LMA.LL[c((nrow(all.results.cwm.cl)-1),nrow(all.results.cwm.cl))],")"),cex=.8)
text(y=par()$usr[4]-.2, x=c(1,2,3), labels = p$n[c(1,2,3)])
text(y=par()$usr[4]-.2, x=c(4,5), labels=paste0("(",all.results.cwm.cl$n_LMA.LL[c((nrow(all.results.cwm.cl)-2),nrow(all.results.cwm.cl))],")"),cex=.9)




par(mar=c(4,4,1.5,1.5), mgp=c(2.5,1,0))
plot(log.LL~log.LMA, allspp, typ="n",col="grey", pch=16, ylab="Leaf Lifespan (months)", xlab="Leaf Mass per Area (g/m2)", xaxt="n", yaxt="n")
abline(a=tmpdata$Int_LMA.LL[nrow(tmpdata)], b=tmpdata$Slope_LMA.LL[nrow(tmpdata)], lwd=3, col="black", lty=3)
#points(log.LL~log.LMA, allspp[which(allspp$Species=="Solanum straminifolia"),])
#points(log.LL~log.LMA, allspp[which(allspp$Species=="Juniperus monosperma"),])
#points(log.LL~log.LMA, allspp[which(allspp$Species=="Araucaria araucana"),])
#points(log.LL~log.LMA, allspp[which(allspp$Species=="Abies magnifica"),])
axis(2, at=log(c(1,2,3,4,5,10,20,30,40,50,100,200,300), base=10), labels = F)#c(1,"","","",5,"","","","",50,"",200,""),)
axis(2, at=log(c(1,5,50,300), base=10), labels = c(1,5,50,300))
axis(1, at=log(c(10,20,30,40,50,100,200,300,400,500,1000), base=10), labels = F)# c(10,"","","",50,"","","","","","1000"))
axis(1, at=log(c(20,50,500), base=10), labels =  c(20,50,500))
points(log.cw_LLp_if~log.cw_LMAp_if, biomass, pch=16, col="grey")
# abline(a=tmpdata$Int_LMA.LL[nrow(tmpdata)-1], b=tmpdata$Slope_LMA.LL[nrow(tmpdata)-1], lwd=3, col="black")
abline(a=tmpdata$Int_LMA.LL[nrow(tmpdata)], b=tmpdata$Slope_LMA.LL[nrow(tmpdata)], lwd=3, col="black")
#plot.MAR(xvar = "log.cw_LMAp_if", yvar="log.cw_LLp_if", data=biomass, linecol = mypal[5], lwd = 4, lty=ltyns)


