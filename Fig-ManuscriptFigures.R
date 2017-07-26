################################################################
###      FIGURES for of NACP_TERRA trait and stand data
###         to assess within-species trait variation
###
###         data downloaded from: http://dx.doi.org/10.3334/ORNLDAAC/1292 
###            on 01/23/16 by LDLA
###############################################################
#### Ecology Letters widths: single column (82 mm), two-thirds page width (110 mm) or full page width (173 mm)
# 3.23', 4.33', or 6.81'


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





#________________________________________________________________
###### FIG 2: Nmass Relations ######
#________________________________________________________________
# 6 panel figure with
# a-d: boxplots of MA slope and Rho for NmassvLMA and NmassvLL
# e&f: funnel plots 
# I think the easiest thing will be to make two quartez: top with boxplots and bottom with funnel

# 
# #________________________________________________________________________________________________
# ########## .Old Version with boxplots and funnel plots and not examples ########
# colchoices <- c(1,2,4,3,6)
# ##### Uppper Panel: boxplots ............................................................
# quartz(width=5.5, height=4.5,title = "LMA v Nmass")
# par(mfrow=c(2,2), mar=c(0,4,0,1), oma=c(6,0,2,0), cex=1)
# palette(mypal[colchoices])
# ### Slope boxplots (a & b)
#   ## LMA v Nmass
# p <- boxplot(Slope_LMA.N~Type, all.results.cl[which(all.results.cl$n_LMA.N>5& !all.results.cl$Type %in% c("global","Famclean") ),]
#              , ylim=c(-2.2,1.5),las=3, ylab="MA Slope" #, main="log(LMA)~log(Nmass) strict"
#              , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
#              ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
#              , staplelwd=0, outpch=16, outcex=.5, outcol=mypal[colchoices]
#              , boxwex=.7, xaxt="n")
# m <- lmodel2(log.Nmass~log.LMA, allspp)
# points(y=m$regression.results$Slope[3],x=5, pch=16, cex=1.3, col="darkgrey")
# arrows(x0=5,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`2.5%-Slope`[3], length = 0,lwd=3, col="darkgrey")
# arrows(x0=5,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`97.5%-Slope`[3], length = 0,lwd=3, col="darkgrey")
# m <- lmodel2(log.Mmass~log.LMA, fam.dataclean)
# points(y=m$regression.results$Slope[3],x=4, pch=16, cex=1.3, col="black")
# arrows(x0=4,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`2.5%-Slope`[3], length = 0,lwd=3)
# arrows(x0=4,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`97.5%-Slope`[3], length = 0,lwd=3)
# 
# abline(h=0, lty=2)
# #text(y=par()$usr[3]+.2, x=c(1,2,3,4,5,6), labels = p$n)
# #text(y=par()$usr[4]-.2,x=.5, labels = "a)")
# mtext(text = "a)", side = 3, adj=0, line=.2)
# mtext(text= "LMA vs Nmass", side=3, line=.2)
#   # LL v Nmass
# p <- boxplot(Slope_LL.N~Type, all.results.cl[which(all.results.cl$n_LL.N>5& !all.results.cl$Type %in% c("global","Famclean") ),]
#              , ylim=c(-1.2,1.2),las=3, ylab="MA Slope" #, main="log(LMA)~log(Nmass) strict"
#              , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
#              ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
#              , staplelwd=0, outpch=16, outcex=.5, outcol=mypal[colchoices]
#              , boxwex=.7, xaxt="n")
# m <- lmodel2(log.Nmass~log.LL, allspp)
# points(y=m$regression.results$Slope[3],x=5, pch=16, cex=1.3, col="darkgrey")
# arrows(x0=5,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`2.5%-Slope`[3], length = 0,lwd=3, col="darkgrey")
# arrows(x0=5,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`97.5%-Slope`[3], length = 0,lwd=3, col="darkgrey")
# m <- lmodel2(log.Mmass~log.LL, fam.dataclean)
# points(y=m$regression.results$Slope[3],x=4, pch=16, cex=1.3, col="black")
# arrows(x0=4,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`2.5%-Slope`[3], length = 0,lwd=3)
# arrows(x0=4,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`97.5%-Slope`[3], length = 0,lwd=3)
# 
# 
# abline(h=0, lty=2)
# #text(y=par()$usr[3]+.2, x=c(1,2,3,4,5,6), labels = p$n)
# #text(y=par()$usr[4]-.2,x=.5, labels = "b)")
# mtext(text = "b)", side = 3, adj=0, line=.2)
# mtext(text= "LL vs Nmass", side=3, line=.2)
# 
# ### Lower boxplots with Rho
#   # LMA v Nmass
# # p <- boxplot(Rho_LMA.N~Type, all.results.cl[which(all.results.cl$n_LMA.N>5 & is.na(all.results.cl$Taxo.Unit)),]
# #              , ylim=c(-1,1.1),las=3, ylab="Rho"
# #              , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
# #              ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
# #              , staplelwd=0, outpch=16, outcex=.5, outcol=mypal[colchoices]
# #              , boxwex=.7, xaxt="n")
# # #lines(meancis$m10_LMA.N~c(1,2,3,4,5))
# # #lines(meancis$m90_LMA.N~c(1,2,3,4,5))
# #polygon(x = c(1,2,3,3,2,1), y=na.omit(c(meancis$m10_LMA.N, rev(meancis$m90_LMA.N))), border="#EEEEEE",col="#EEEEEE")#lightgrey",col = "lightgrey")
# p <- boxplot(Rho_LMA.N~Type, all.results.cl[which(all.results.cl$n_LMA.N>5 & !all.results.cl$Type %in% c("global","Famclean") ),]
#              , ylim=c(-1,1.1),las=3, ylab="Rho"
#              , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
#              ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
#              , staplelwd=0, outpch=16, outcex=.5, outcol=mypal[colchoices]
#              , boxwex=.7, xaxt="n")
# 
# 
#  m <- cor.test(x=allspp$log.LMA, y=allspp$log.Nmass)
#  points(y=m$estimate,x=5, pch=16, cex=1.3, col="darkgrey")
#  arrows(x0=5,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3, col="darkgrey")
# m <- cor.test(x=fam.dataclean$log.LMA, y=fam.dataclean$log.Nmass)
# points(y=m$estimate,x=4, pch=16, cex=1.3, col="black")
# arrows(x0=4,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3)
# 
# 
# 
# axis(1,labels = c("w/in Spp","w/in Gen", "w/in Fam","btw Fam","global"), at=c(1,2,3,4,5), las=3)
# abline(h=0, lty=2)
# text(y=par()$usr[4]-.2, x=c(1,2,3), labels = p$n[1:3])
# text(y=par()$usr[4]-.2,x=c(4,5), labels=paste0("(",all.results.cl$n_LMA.N[c((nrow(all.results.cl)-1),nrow(all.results.cl))],")"), cex=0.7)
# #text(y=par()$usr[4]-.2,x=.5, labels = "b)")
#   # LL an Nmass
# p <- boxplot(Rho_LL.N~Type, all.results.cl[which(all.results.cl$n_LL.N>5& !all.results.cl$Type %in% c("global","Famclean") ),]
#              , ylim=c(-1,1.1),las=3, ylab="Rho"
#              , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
#              ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
#              , staplelwd=0, outpch=16, outcex=.5, outcol=mypal[colchoices]
#              , boxwex=.7, xaxt="n")
# axis(1,labels = c("w/in Spp","w/in Gen", "w/in Fam","btw Fam","global"), at=c(1,2,3,4,5), las=3)
# abline(h=0, lty=2)
# text(y=par()$usr[4]-.2, x=c(1,2,3), labels = p$n[1:3])
# text(y=par()$usr[4]-.2, x=c(4,5), labels=paste0("(",all.results.cl$n_LL.N[c((nrow(all.results.cl)-1),nrow(all.results.cl))],")"),cex=.7)
# m <- cor.test(x=allspp$log.LL, y=allspp$log.Nmass)
# points(y=m$estimate,x=5, pch=16, cex=1.3, col="darkgrey")
# arrows(x0=5,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3, col="darkgrey")
# m <- cor.test(x=fam.dataclean$log.LL, y=fam.dataclean$log.Nmass)
# points(y=m$estimate,x=4, pch=16, cex=1.3, col="black")
# arrows(x0=4,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3)
# 
# #text(y=par()$usr[4]-.2,x=.5, labels = "b)")
# 
# 
# 
# 





#________________________________________________________________________________________________
########## .New Version similar to all other covariation plots ########
#_________________________________________________________________________________________________

## LMA v Nmass boxplots
mat <- matrix(c(1,3,
                2,3), nrow=2, byrow = T)
colchoices <- c(1,2,4,3,6)
palette(mypal[colchoices])

quartz(width=7.008, height=4)
layout(mat)
par(mar=c(0,4,5,1))

# Slope boxplots
p <- boxplot(Slope_LMA.N~Type, all.results.cl[which(all.results.cl$n_LMA.N>5& !all.results.cl$Type %in% c("global","Famclean") ),]
             , ylim=c(-2.2,1.5),las=3, ylab="MA Slope" #, main="log(LMA)~log(Nmass) strict"
             , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
             ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
             , staplelwd=0, outpch=16, outcex=.7, outcol=mypal[colchoices]
             , boxwex=.7, xaxt="n")
  # add global and bewteen family points and error bars
m <- lmodel2(log.Nmass~log.LMA, allspp)
points(y=m$regression.results$Slope[3],x=5, pch=16, cex=1.3, col="darkgrey")
arrows(x0=5,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`2.5%-Slope`[3], length = 0,lwd=3, col="darkgrey")
arrows(x0=5,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`97.5%-Slope`[3], length = 0,lwd=3, col="darkgrey")
m <- lmodel2(log.Mmass~log.LMA, fam.dataclean)
points(y=m$regression.results$Slope[3],x=4, pch=16, cex=1.3, col="black")
arrows(x0=4,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`2.5%-Slope`[3], length = 0,lwd=3)
arrows(x0=4,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`97.5%-Slope`[3], length = 0,lwd=3)
abline(h=0, lty=2)
mtext(text = "a)", side = 3, adj=0, line=.2)
mtext(text= expression(paste("LMA vs ", N[mass], sep=" ")), side=3, line=.2)

# Rho boxplots
par(mar=c(5,4,0,1))
p <- boxplot(Rho_LMA.N~Type, all.results.cl[which(all.results.cl$n_LMA.N>5 & !all.results.cl$Type %in% c("global","Famclean") ),]
             , ylim=c(-1,1.1),las=3, ylab="Rho"
             , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
             ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
             , staplelwd=0, outpch=16, outcex=.7, outcol=mypal[colchoices]
             , boxwex=.7, xaxt="n")
  # add global and between family points
m <- cor.test(x=allspp$log.LMA, y=allspp$log.Nmass)
points(y=m$estimate,x=5, pch=16, cex=1.3, col="darkgrey")
arrows(x0=5,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3, col="darkgrey")
m <- cor.test(x=fam.dataclean$log.LMA, y=fam.dataclean$log.Nmass)
points(y=m$estimate,x=4, pch=16, cex=1.3, col="black")
arrows(x0=4,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3)
  # add axes and sample sizes
axis(1,labels = c("w/in Spp","w/in Gen", "w/in Fam","btw Fam","global"), at=c(1,2,3,4,5), las=3)
abline(h=0, lty=2)
text(y=par()$usr[4]-.2, x=c(1,2,3), labels = p$n[1:3])
text(y=par()$usr[4]-.2,x=c(4,5), labels=paste0("(",all.results.cl$n_LMA.N[c((nrow(all.results.cl)-1),nrow(all.results.cl))],")"), cex=0.7)



# scatterplot w/ lines
par(mar=c(4,4,1.5,1.5), cex.lab=1.2, mgp=c(2.5,.7,0))
plot(log.Nmass~log.LMA, LES, col="grey", pch=16, xlab=expression(paste(log[10],"(LMA)")), ylab=expression(paste(log[10],(N[mass]), sep=" ")))
tax <- "w.inGen"
for (i in as.character(all.results.cl$Taxo.Unit[which(all.results.cl$Type==tax & all.results.cl$n_LMA.N>5)])){
  plot.MAR(xvar = "log.LMA", yvar = "log.Nmass",data= gen.data[which(gen.data$Genus==i),], linecol = mypal[colchoices[2]])
}
tax <- "Genw.inFam"
for (i in as.character(all.results.cl$Taxo.Unit[which(all.results.cl$Type==tax & all.results.cl$n_LMA.N>5)])){
  plot.MAR(xvar = "log.LMA", yvar = "log.Nmass",data= geninfam.data[which(geninfam.data$Family==i),], linecol = mypal[colchoices[3]])
}
abline(a=all.results.cl$Int_LMA.N[nrow(all.results.cl)-1], b=all.results.cl$Slope_LMA.N[nrow(all.results.cl)-1], lwd=3, col="black")
abline(a=all.results.cl$Int_LMA.N[nrow(all.results.cl)], b=all.results.cl$Slope_LMA.N[nrow(all.results.cl)], lwd=3, col="black", lty=3)
tax <- "w.inSpp"
for (i in as.character(all.results.cl$Taxo.Unit[which(all.results.cl$Type==tax & all.results.cl$n_LMA.N>5)])){
  plot.MAR(xvar = "log.LMA", yvar = "log.Nmass",data= spp.data[which(spp.data$Species==i),], linecol = mypal[colchoices[1]])
}
mtext(text = "b)", side = 3, adj=0, line=.2)





# LL v Nmass ________________________
mat <- matrix(c(1,3,
                2,3), nrow=2, byrow = T)
colchoices <- c(1,2,4,3,6)
palette(mypal[colchoices])

quartz(width=7.008, height=4)
layout(mat)
par(mar=c(0,4,5,1))

p <- boxplot(Slope_LL.N~Type, all.results.cl[which(all.results.cl$n_LL.N>5& !all.results.cl$Type %in% c("global","Famclean") ),]
             , ylim=c(-1.2,1.2),las=3, ylab="MA Slope" #, main="log(LMA)~log(Nmass) strict"
             , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
             ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
             , staplelwd=0, outpch=16, outcex=.7, outcol=mypal[colchoices]
             , boxwex=.7, xaxt="n")
m <- lmodel2(log.Nmass~log.LL, allspp)
points(y=m$regression.results$Slope[3],x=5, pch=16, cex=1.3, col="darkgrey")
arrows(x0=5,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`2.5%-Slope`[3], length = 0,lwd=3, col="darkgrey")
arrows(x0=5,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`97.5%-Slope`[3], length = 0,lwd=3, col="darkgrey")
m <- lmodel2(log.Mmass~log.LL, fam.dataclean)
points(y=m$regression.results$Slope[3],x=4, pch=16, cex=1.3, col="black")
arrows(x0=4,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`2.5%-Slope`[3], length = 0,lwd=3)
arrows(x0=4,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`97.5%-Slope`[3], length = 0,lwd=3)
abline(h=0, lty=2)
mtext(text = "c)", side = 3, adj=0, line=.2)
mtext(text= expression(paste("LL vs ",N[mass])), side=3, line=.2)


# Rho boxplot
par(mar=c(5,4,0,1))
p <- boxplot(Rho_LL.N~Type, all.results.cl[which(all.results.cl$n_LL.N>5& !all.results.cl$Type %in% c("global","Famclean") ),]
             , ylim=c(-1,1.1),las=3, ylab="Rho"
             , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
             ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
             , staplelwd=0, outpch=16, outcex=.7, outcol=mypal[colchoices]
             , boxwex=.7, xaxt="n")
axis(1,labels = c("w/in Spp","w/in Gen", "w/in Fam","btw Fam","global"), at=c(1,2,3,4,5), las=3)
abline(h=0, lty=2)
text(y=par()$usr[4]-.2, x=c(1,2,3), labels = p$n[1:3])
text(y=par()$usr[4]-.2, x=c(4,5), labels=paste0("(",all.results.cl$n_LL.N[c((nrow(all.results.cl)-1),nrow(all.results.cl))],")"),cex=.7)
m <- cor.test(x=allspp$log.LL, y=allspp$log.Nmass)
points(y=m$estimate,x=5, pch=16, cex=1.3, col="darkgrey")
arrows(x0=5,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3, col="darkgrey")
m <- cor.test(x=fam.dataclean$log.LL, y=fam.dataclean$log.Nmass)
points(y=m$estimate,x=4, pch=16, cex=1.3, col="black")
arrows(x0=4,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3)

# scatterplot w/ lines
par(mar=c(4,4,1.5,1.5), cex.lab=1.2, mgp=c(2.5,.7,0))
plot(log.Nmass~log.LL, LES, col="grey", pch=16, xlab=expression(paste(log[10](LL))), ylab=expression(paste(log[10],(N[mass]))))
tax <- "w.inGen"
for (i in as.character(all.results.cl$Taxo.Unit[which(all.results.cl$Type==tax & all.results.cl$n_LL.N>5)])){
  plot.MAR(xvar = "log.LL", yvar = "log.Nmass",data= gen.data[which(gen.data$Genus==i),], linecol = mypal[colchoices[2]])
}
tax <- "Genw.inFam"
for (i in as.character(all.results.cl$Taxo.Unit[which(all.results.cl$Type==tax & all.results.cl$n_LL.N>5)])){
  plot.MAR(xvar = "log.LL", yvar = "log.Nmass",data= geninfam.data[which(geninfam.data$Family==i),], linecol = mypal[colchoices[3]])
}
abline(a=all.results.cl$Int_LL.N[nrow(all.results.cl)-1], b=all.results.cl$Slope_LL.N[nrow(all.results.cl)-1], lwd=3, col="black")
abline(a=all.results.cl$Int_LL.N[nrow(all.results.cl)], b=all.results.cl$Slope_LL.N[nrow(all.results.cl)], lwd=3, col="black", lty=3)
tax <- "w.inSpp"
for (i in as.character(all.results.cl$Taxo.Unit[which(all.results.cl$Type==tax & all.results.cl$n_LL.N>5)])){
  plot.MAR(xvar = "log.LL", yvar = "log.Nmass",data= spp.data[which(spp.data$Species==i),], linecol = mypal[colchoices[1]])
}
mtext(text = "d)", side = 3, adj=0, line=.2)



#_________________________________________________________________________________________________











#________________________________________________________________
###### FIG 3: LMA v LL Relations ######
#________________________________________________________________
# 2 panel figure with boxplots and 
# a-d: boxplots of MA slope and Rho for NmassvLMA and NmassvLL
# e&f: funnel plots 
# width = 1.5 columns -> 11.4 cm,4.89
#       = 2 columns -> 17.8 cm,7.008
# I think the easiest thing will be to make two quartez: top with boxplots and bottom with funnel

mat <- matrix(c(1,3,
                2,3), nrow=2, byrow = T)
colchoices <- c(1,2,4,3,6)
palette(mypal[colchoices])

quartz(width=7.008, height=4)
layout(mat)
par(mar=c(0,4,6,1))
p <- boxplot(Slope_LMA.LL~Type, all.results.cl[which(all.results.cl$n_LMA.LL>5& !all.results.cl$Type %in% c("global","Famclean")),]
             , ylim=c(-2.5,4),las=3, ylab="MA Slope" #, main="log(LMA)~log(Nmass) strict"
             , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
             ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
             , staplelwd=0, outpch=16, outcex=.7, outcol=mypal[colchoices]
             , boxwex=.7, xaxt="n")
abline(h=0, lty=2)
m <- lmodel2(log.LL~log.LMA, allspp)
points(y=m$regression.results$Slope[3],x=5, pch=16, cex=1.3, col="darkgrey")
arrows(x0=5,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`2.5%-Slope`[3], length = 0,lwd=3, col="darkgrey")
arrows(x0=5,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`97.5%-Slope`[3], length = 0,lwd=3, col="darkgrey")
m <- lmodel2(log.LL~log.LMA, fam.dataclean)
points(y=m$regression.results$Slope[3],x=4, pch=16, cex=1.3, col="black")
arrows(x0=4,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`2.5%-Slope`[3], length = 0,lwd=3)
arrows(x0=4,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`97.5%-Slope`[3], length = 0,lwd=3)

#text(y=par()$usr[3]+.2, x=c(1,2,3,4,5,6), labels = p$n)
#text(y=par()$usr[4]-.2,x=.5, labels = "a)")
mtext(text = "a)", side = 3, adj=0, line=.2)
mtext(text= "LMA vs LL", side=3, line=.2)
par(mar=c(6,4,0,1))
p <- boxplot(Rho_LMA.LL~Type, all.results.cl[which(all.results.cl$n_LMA.LL>5& !all.results.cl$Type %in% c("global","Famclean")),]
             , ylim=c(-1,1.4),las=3, ylab="Rho"
             , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
             ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
             , staplelwd=0, outpch=16, outcex=.7, outcol=mypal[colchoices]
             , boxwex=.7, xaxt="n")
axis(1,labels = c("w/in Spp","w/in Gen", "w/in Fam","btw Fam","global"), at=c(1,2,3,4,5), las=3)
abline(h=0, lty=2)
text(y=par()$usr[4]-.2, x=c(1,2,3), labels = p$n[1:3])
text(y=par()$usr[4]-.2, x=c(4,5), labels=paste0("(",all.results.cl$n_LMA.LL[c((nrow(all.results.cl)-1),nrow(all.results.cl))],")"),cex=.8)
#polygon(x = c(1,2,3,3,2,1), y=na.omit(c(meancis$m10_LMA.LL, rev(meancis$m90_LMA.LL))), border="#EEEEEE",col="#EEEEEE")#lightgrey",col = "lightgrey")
m <- cor.test(x=allspp$log.LMA, y=allspp$log.LL)
points(y=m$estimate,x=5, pch=16, cex=1.3, col="darkgrey")
arrows(x0=5,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3, col="darkgrey")
m <- cor.test(x=fam.dataclean$log.LMA, y=fam.dataclean$log.LL)
points(y=m$estimate,x=4, pch=16, cex=1.3, col="black")
arrows(x0=4,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3)


par(mar=c(4,4,1.5,1.5))
## scatterplot
####### Plotting the LMA vs Narea scaling in unit rather than log space ######
plot(log.LL~log.LMA, LES, col="grey", pch=16, ylab=expression(paste(log[10](LL))), xlab=expression(paste(log[10](LMA))))

tax <- "w.inGen"
for (i in as.character(all.results.cl$Taxo.Unit[which(all.results.cl$Type==tax & all.results.cl$n_LMA.LL>5)])){
  plot.MAR(xvar = "log.LMA", yvar = "log.LL",data= gen.data[which(gen.data$Genus==i),], linecol = mypal[colchoices[2]])
}

tax <- "Genw.inFam"
for (i in as.character(all.results.cl$Taxo.Unit[which(all.results.cl$Type==tax & all.results.cl$n_LMA.LL>5)])){
  plot.MAR(xvar = "log.LMA", yvar = "log.LL",data= geninfam.data[which(geninfam.data$Family==i),], linecol = mypal[colchoices[3]])
}


abline(a=all.results.cl$Int_LMA.LL[nrow(all.results.cl)-1], b=all.results.cl$Slope_LMA.LL[nrow(all.results.cl)-1], lwd=3, col="black")
abline(a=all.results.cl$Int_LMA.LL[nrow(all.results.cl)], b=all.results.cl$Slope_LMA.LL[nrow(all.results.cl)], lwd=3, col="black", lty=3)
tax <- "w.inSpp"
for (i in as.character(all.results.cl$Taxo.Unit[which(all.results.cl$Type==tax & all.results.cl$n_LMA.LL>5)])){
  plot.MAR(xvar = "log.LMA", yvar = "log.LL",data= spp.data[which(spp.data$Species==i),], linecol = mypal[colchoices[1]])
}

#plot.MAR(xvar="log.LMA", yvar="log.LL", data= fam.data, linecol = mypal[colchoices[4]], lwd=2)

#plot.MAR(xvar="log.LMA", yvar="log.LL", data= LES, linecol = "black", lwd=2)
mtext(text = "b)", side = 3, adj=0, line=.2)










#________________________________________________________________
###### FIG 4: LMA v Narea Relations ######
#________________________________________________________________
# 2 panel figure with boxplots and scatterplot
# width = 1.5 columns -> 11.4 cm,4.89
#       = 2 columns -> 17.8 cm,7.008
# I think the easiest thing will be to make two quartez: top with boxplots and bottom with funnel
mat <- matrix(c(1,3,
                2,3), nrow=2, byrow = T)
colchoices <- c(1,2,4,3,6)
palette(mypal[colchoices])

quartz(width=7.008, height=4)
layout(mat)
par(mar=c(0,4,6,1))
p <- boxplot(Slope_LMA.Narea~Type, all.results.cl[which(all.results.cl$n_LMA.Narea>5& !all.results.cl$Type %in% c("global","Famclean")& all.results.cl$Slope_LMA.Narea>-1),]
             , ylim=c(-0.5,2.5),las=3, ylab="MA Slope" #, main="log(LMA)~log(Nmass) strict"
             , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
             ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
             , staplelwd=0, outpch=16, outcex=.7, outcol=mypal[colchoices]
             , boxwex=.7, xaxt="n")
  # had to remove the two crazy outliers (Protea repens and genus Abies)
abline(h=0, lty=2)
abline(h=1, col="grey")
m <- lmodel2(log.Narea~log.LMA, allspp)
points(y=m$regression.results$Slope[3],x=5, pch=16, cex=1.3, col="darkgrey")
arrows(x0=5,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`2.5%-Slope`[3], length = 0,lwd=3, col="darkgrey")
arrows(x0=5,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`97.5%-Slope`[3], length = 0,lwd=3, col="darkgrey")
m <- lmodel2(log.Marea~log.LMA, fam.dataclean)
points(y=m$regression.results$Slope[3],x=4, pch=16, cex=1.3, col="black")
arrows(x0=4,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`2.5%-Slope`[3], length = 0,lwd=3)
arrows(x0=4,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`97.5%-Slope`[3], length = 0,lwd=3)

#text(y=par()$usr[3]+.2, x=c(1,2,3,4,5,6), labels = p$n)
#text(y=par()$usr[4]-.2,x=.5, labels = "a)")
mtext(text = "a)", side = 3, adj=0, line=.2)
mtext(text= "LMA vs Narea", side=3, line=.2)
par(mar=c(6,4,0,1))
p <- boxplot(Rho_LMA.Narea~Type, all.results.cl[which(all.results.cl$n_LMA.Narea>5& !all.results.cl$Type %in% c("global","Famclean")& all.results.cl$Slope_LMA.Narea>-1),]
             , ylim=c(-1,1.4),las=3, ylab="Rho"
             , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
             ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
             , staplelwd=0, outpch=16, outcex=.7, outcol=mypal[colchoices]
             , boxwex=.7, xaxt="n")
axis(1,labels = c("w/in Spp","w/in Gen", "w/in Fam","btw Fam","global"), at=c(1,2,3,4,5), las=3)
abline(h=0, lty=2)
text(y=par()$usr[4]-.2, x=c(1,2,3), labels = p$n[1:3])
text(y=par()$usr[4]-.2, x=c(4,5), labels=paste0("(",all.results.cl$n_LMA.Narea[c((nrow(all.results.cl)-1),nrow(all.results.cl))],")"),cex=.8)
m <- cor.test(x=allspp$log.LMA, y=allspp$log.Narea)
points(y=m$estimate,x=5, pch=16, cex=1.3, col="darkgrey")
arrows(x0=5,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3, col="darkgrey")
m <- cor.test(x=fam.dataclean$log.LMA, y=fam.dataclean$log.Narea)
points(y=m$estimate,x=4, pch=16, cex=1.3, col="black")
arrows(x0=4,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3)

par(mar=c(4,4,1.5,1.5))
## scatterplot
####### Plotting the LMA vs Narea scaling in unit rather than log space ######
plot(log.Narea~log.LMA, LES, col="grey", pch=16, ylab=expression(paste(log[10](N[area]))), xlab=expression(paste(log[10](LMA))))

tax <- "w.inGen"
for (i in as.character(all.results.cl$Taxo.Unit[which(all.results.cl$Type==tax & all.results.cl$n_LMA.Narea>5)& all.results.cl$Slope_LMA.Narea>-1])){
  plot.MAR(xvar = "log.LMA", yvar = "log.Narea",data= gen.data[which(gen.data$Genus==i),], linecol = mypal[2])
}

tax <- "Genw.inFam"
for (i in as.character(all.results.cl$Taxo.Unit[which(all.results.cl$Type==tax & all.results.cl$n_LMA.Narea>5)])){
  plot.MAR(xvar = "log.LMA", yvar = "log.Narea",data= geninfam.data[which(geninfam.data$Family==i),], linecol = mypal[4])
}

abline(a=all.results.cl$Int_LMA.Narea[nrow(all.results.cl)-1], b=all.results.cl$Slope_LMA.Narea[nrow(all.results.cl)-1], lwd=3, col="black")
abline(a=all.results.cl$Int_LMA.Narea[nrow(all.results.cl)], b=all.results.cl$Slope_LMA.Narea[nrow(all.results.cl)], lwd=3, col="black", lty=3)

# plot.MAR(xvar="log.LMA", yvar="log.Narea", data= fam.data, linecol = mypal[5], lwd=2)
# 
# plot.MAR(xvar="log.LMA", yvar="log.Narea", data= LES, linecol = "black", lwd=2)
tax <- "w.inSpp"
for (i in as.character(all.results.cl$Taxo.Unit[which(all.results.cl$Type==tax & all.results.cl$n_LMA.Narea>5& all.results.cl$Slope_LMA.Narea>-1)])){
  plot.MAR(xvar = "log.LMA", yvar = "log.Narea",data= spp.data[which(spp.data$Species==i),], linecol = mypal[1])
}
mtext(text = "b)", side = 3, adj=0, line=.2)






#________________________________________________________________
###### FIG 5cd:LL v Narea Relations ######
#________________________________________________________________
# 2 panel figure with boxplots and scatterplot
# width = 1.5 columns -> 11.4 cm,4.89
#       = 2 columns -> 17.8 cm,7.008
# I think the easiest thing will be to make two quartez: top with boxplots and bottom with funnel
mat <- matrix(c(1,3,
                2,3), nrow=2, byrow = T)
colchoices <- c(1,2,4,3,6)
palette(mypal[colchoices])

quartz(width=7.008, height=4)
layout(mat)
par(mar=c(0,4,6,1))
p <- boxplot(Slope_LL.Narea~Type, all.results.cl[which(all.results.cl$n_LL.Narea>5),]
             , ylim=c(-2,2),las=3, ylab="MA Slope" #, main="log(LMA)~log(Nmass) strict"
             , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
             ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
             , staplelwd=0, outpch=16, outcex=.5, outcol=mypal[colchoices]
             , boxwex=.7, xaxt="n")

abline(h=0, lty=2)
#text(y=par()$usr[3]+.2, x=c(1,2,3,4,5,6), labels = p$n)
#text(y=par()$usr[4]-.2,x=.5, labels = "a)")
mtext(text = "c)", side = 3, adj=0, line=.2)
mtext(text= "LL vs Narea", side=3, line=.2)
par(mar=c(6,4,0,1))
p <- boxplot(Rho_LL.Narea~Type, all.results.cl[which(all.results.cl$n_LL.Narea>5),]
             , ylim=c(-1,1.4),las=3, ylab="Rho"
             , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
             ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
             , staplelwd=0, outpch=16, outcex=.5, outcol=mypal[colchoices]
             , boxwex=.7, xaxt="n")
abline(h=0, lty=2)
axis(1,labels = c("w/in Spp","w/in Gen", "w/in Fam","btw Fam","global"), at=c(1,2,3,4,5), las=3)
text(y=par()$usr[4]-.2, x=c(1,2,3), labels = p$n[1:3])
text(y=par()$usr[4]-.2, x=c(4,5), labels=paste0("(",all.results.cl$n_LL.Narea[c((nrow(all.results.cl)-1),nrow(all.results.cl))],")"),cex=1)

par(mar=c(4,4,1.5,1.5), mgp=c(2.5,1,0))
## scatterplot

###### testing out leaf lifespan vs Narea relations, I think they'll be ns ########

plot(log.Narea~log.LL, LES, col="grey", pch=16, xlab="log(LL)", ylab="log(Narea)")

tax <- "w.inSpp"
for (i in as.character(all.results.cl$Taxo.Unit[which(all.results.cl$Type==tax & all.results.cl$n_LL.Narea>5)])){
  plot.MAR(yvar = "log.Narea", xvar = "log.LL",data= spp.data[which(spp.data$Species==i),], linecol = mypal[1])
}

tax <- "w.inGen"
for (i in as.character(all.results.cl$Taxo.Unit[which(all.results.cl$Type==tax & all.results.cl$n_LL.Narea>5)])){
  plot.MAR(yvar = "log.Narea", xvar = "log.LL",data= gen.data[which(gen.data$Genus==i),], linecol = mypal[2])
}

tax <- "Genw.inFam"
for (i in as.character(all.results.cl$Taxo.Unit[which(all.results.cl$Type==tax & all.results.cl$n_LL.Narea>5)])){
  plot.MAR(yvar = "log.Narea", xvar = "log.LL",data= geninfam.data[which(geninfam.data$Family==i),], linecol = mypal[4])
}

abline(a=all.results.cl$Int_LL.Narea[nrow(all.results.cl)-1], b=all.results.cl$Slope_LL.Narea[nrow(all.results.cl)-1], lwd=3, col="black")
abline(a=all.results.cl$Int_LL.Narea[nrow(all.results.cl)], b=all.results.cl$Slope_LL.Narea[nrow(all.results.cl)], lwd=3, col="black", lty=3)

# plot.MAR(xvar="log.Narea", yvar="log.LL", data= fam.dataclean, linecol = mypal[5], lwd=2)
# 
# plot.MAR(xvar="log.Narea", yvar="log.LL", data= LES, linecol = "black", lwd=2)
mtext(text = "d)", side = 3, adj=0, line=.2)
















#________________________________________________________________
###### FIG 6: Growth~Trait relationships ######
#________________________________________________________________


####### Figure 6: reploted with traits.mono and biost_growth instead of RGR.
# for 6 panel figure without Nmass, but with rGR top row and biostGRWOTH bottom
quartz(width=6.8, height=4) # Eco Let width is 173mm, 110mm or 82mm
#pdf(file = "Traits-v-Growth_v2.pdf",width = 7.1, height=6)
par(mar=c(0,0,0,0), oma=c(3.5,3.5,1.2,2),mfrow=c(2,4), cex=1.1,mgp=c(2.5,.5,0))

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
plot(log(RGR)~mlog.LL, plotavs90, col=SP.ID, pch=16, cex=cex.pts,xaxt="n", xlab="", ylim=ylims) #xlab="log(Leaf Lifespan)")
#mtext("Relative Growth Rate", side=2, line=2, cex=1.1)
#mtext("log(Leaf Life)", side=1, line=1.5, cex=1.1)
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
plot(log(RGR)~mlog.LMA, plotavs90, col=SP.ID, pch=16, cex=cex.pts,xaxt="n", yaxt='n',ylim=ylims) #xlab="log(LMA)")#, xlim=c(2,2.6))#log(LMA)")
#mtext("log(LMA)", side=1, line=1.5, cex=1.1)
#lines(c(min(LMAs), max(LMAs)), c(max(f0lma), min(f0lma)), lwd=4, lty=3)
for(i in levels(spid)){
  lines(LMAs[which(spid==i)], f1lma[I][which(spid==i)],  col=spid[which(spid==i)], lwd=3)
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
plot(log(RGR)~mlog.Nmass, plotavs90, col=SP.ID, pch=16, cex=cex.pts, xlab="", xlim=c(-.22,.25), xaxt="n", yaxt="n",ylim=ylims)#log(Nmass)")
#mtext("log(Nmass)", side=1, line=1.5, cex=1.1)
#lines(Nmasss, f0Nmass[I], lwd=4, lty=1)

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
plot(log(RGR)~mlog.Narea, plotavs90, col=SP.ID, pch=16, cex=cex.pts, xlab="", xaxt="n", yaxt="n",ylim=ylims)#log(Narea)")
#mtext("Relative Growth Rate", side=2, line=2, cex=1.1)
#mtext("log(Narea)", side=1, line=1.5, cex=1.1)
#lines(Nareas, f0Narea[I], lwd=4, lty=3)
for(i in levels(spid)){
  lines(Nareas[which(spid==i)], f1Narea[I][which(spid==i)],  col=spid[which(spid==i)], lwd=3)
}

mtext(text = paste0("p=",round(anova(modNareavRGR, modNareavRGRnull)$`Pr(>Chisq)`[2],3)), side=3, adj = .05, line=-1 ,cex=cex.text)
mtext("d)", side=3, adj=0)
mtext("w/in spp", side=4, line=0)

###### Community-weighted traits

plot(log(RGR)~log.cw_LLp_if, biomass, pch=16, ylim=ylims, cex=cex.pts)
abline(lm(log(RGR)~log.cw_LLp_if, biomass), lwd=3)
mtext(text = paste0("p=", round(anova(LLmod2)$`Pr(>F)`[1],8), ","), side=3, adj=.1, line=-1, cex=cex.text)
mtext(text = bquote(~R^2==.(round(summary(LLmod2)$r.squared,2))), side=3, adj = .95, line=-1 , cex=cex.text)
mtext("log(RGR)", side=2, line=2, cex=1.1, adj=1.5)
mtext("log(Leaf Life)", side=1, line=2, cex=1.1)
#mtext("a)", side=3, adj=0)

plot(log(RGR)~log.cw_LMAp_if, biomass, pch=16, yaxt="n", ylim=ylims, cex=cex.pts)
abline(lm(log(RGR)~log.cw_LMAp_if, biomass), lwd=3, lty=2)
mtext(text = paste0("p=", round(anova(LMAmod2)$`Pr(>F)`[1],2)), side=3, adj=.05, line=-1, cex=cex.text)
mtext("log(LMA)", side=1, line=2, cex=1.1)
#mtext("b)", side=3, adj=0)

plot(log(RGR)~log.cw_Nmassp_if, biomass, pch=16, ylim=ylims, yaxt="n", cex=cex.pts)
abline(lm(log(RGR)~log.cw_Nmassp_if, biomass), lwd=3, lty=1)
mtext(text = paste0("p=", round(anova(Nmassmod2)$`Pr(>F)`[1],3),","), side=3, adj=.05, line=-1, cex=cex.text)
mtext(text = bquote(~R^2==.(round(summary(Nmassmod2)$r.squared,2))), side=3, adj = .95, line=-1, cex=cex.text)
#mtext("log(RGR)", side=2, line=2, cex=1.1)
mtext("log(Nmass)", side=1, line=2, cex=1.1)
#mtext("c)", side=3, adj=0)

plot(log(RGR)~log.cw_Nareap_if, biomass, pch=16, ylim=ylims, yaxt="n", cex=cex.pts)#, yaxt="n")
#axis(2,at = log(c(.01,.02,.03,0.04,.1,0.2,0.3)), labels =c(.01,.02,.03,0.04,.1,0.2,0.3) )
abline(lm(log(RGR)~log.cw_Nareap_if, biomass), lwd=3, lty=2)
mtext(text = paste0("p=", round(anova(Nareamod2)$`Pr(>F)`[1],2)), side=3, adj=.05, line=-1, cex=cex.text)
mtext("log(Narea)", side=1, line=2, cex=1.1)
#mtext("d)", side=3, adj=0)
mtext("CW means", side=4, line=0)










#________________________________________________________________
###### FIG S1: Funnel Plots for LMA vs LL ######
#________________________________________________________________


quartz(width=5.5,height=4.5)
par(mar=c(4,4,0,1), mfrow=c(2,2), mgp=c(2.5,1,0), oma=c(0,0,2,0))
plot(Slope_LMA.LL~varLMA, all.results.cl, pch=16, col=Type, xlab="Var. in LMA", ylab="Slope   (LMA v LL)")
#mtext( text="Nmass v LMA", side=3, line=0, font=2)
abline(h=0, col="grey", lty=2)
abline(h=all.results$Slope_LMA.LL[which(all.results$Taxo.Unit=="fam.all")])
mtext(text = "a)", side = 3, adj=0, line=.2)
points(Slope_LMA.LL~varLMA, all.results.cl[c(nrow(all.results.cl)-1,nrow(all.results.cl)),], pch=c(24,25), cex=1.4, bg=mypal[colchoices[4:5]])
abline(h=mean(all.results$Slope_LMA.LL[which(all.results$Type=="w.inSpp")], na.rm=T), col=mypal[1])


plot(Rho_LMA.LL~varLMA, all.results.cl, pch=16, col=Type, xlab="Var. in LMA", ylab="Rho   (LMA v LL)")
#mtext( text="Nmass v LMA", side=3, line=0, font=2)
abline(h=0, col="grey", lty=2)
abline(h=all.results.cl$Rho_LMA.LL[which(all.results.cl$Taxo.Unit=="fam.clean")])
mtext(text = "b)", side = 3, adj=0, line=.2)
points(Rho_LMA.LL~varLMA, all.results.cl[c(nrow(all.results.cl)-1,nrow(all.results.cl)),], pch=c(24,25), cex=1.4, bg=mypal[colchoices[4:5]])
abline(h=mean(all.results.cl$Rho_LMA.LL[which(all.results.cl$Type=="w.inSpp")], na.rm=T), col=mypal[1])


plot(Slope_LMA.LL~varLL, all.results.cl, pch=16, col=Type, xlab="Var. in LL", ylab="Slope   (LMA v LL)")
abline(h=0, col="grey", lty=2)
abline(h=all.results.cl$Slope_LMA.LL[which(all.results.cl$Taxo.Unit=="fam.clean")])
mtext(text = "c)", side = 3, adj=0, line=.2)
points(Slope_LMA.LL~varLL, all.results.cl[c(nrow(all.results.cl)-1,nrow(all.results.cl)),], pch=c(24,25), cex=1.4, bg=mypal[5:6])
abline(h=mean(all.results.cl$Slope_LMA.LL[which(all.results.cl$Type=="w.inSpp")], na.rm=T), col=mypal[1])

plot(Rho_LMA.LL~varLL, all.results.cl, pch=16, col=Type, xlab="Var. in LL", ylab="Rho   (LMA v LL)")
abline(h=0, col="grey", lty=2)
abline(h=all.results.cl$Rho_LMA.LL[which(all.results.cl$Taxo.Unit=="fam.clean")])
#legend('bottomright', legend = levels(all.results.cl$Type), pch=c(16,16,16,16,24,25), ncol=2, col=c(mypal[1:4], "black","black"), pt.bg= mypal[1:6], bty ="n", cex=.7)
mtext(text = "d)", side = 3, adj=0, line=.2)
points(Rho_LMA.LL~varLL, all.results.cl[c(nrow(all.results.cl)-1,nrow(all.results.cl)),], pch=c(24,25), cex=1.4, bg=mypal[colchoices[4:5]])
abline(h=mean(all.results.cl$Rho_LMA.LL[which(all.results.cl$Type=="w.inSpp")], na.rm=T), col=mypal[1])
legend('bottomright', legend = levels(all.results.cl$Type), pch=c(16,16,16,24,25), ncol=2, col=c(mypal[colchoices[1:3]], "black","black"), pt.bg= mypal[colchoices[1:5]], bty ="n", cex=.7)




# ##### Lower Panel: funnel plots ............................................................

quartz(width=5.5,height=3)
par(mar=c(4,4,0,1), mfrow=c(1,2), mgp=c(2.5,1,0), oma=c(0,0,2,0))
palette(mypal[colchoices])
plot(Rho_LMA.N~varNmass, all.results.cl, pch=16, col=Type, xlab="Var. in %N", ylab="Rho   (LMA v Nmass)")
#mtext( text="Nmass v LMA", side=3, line=0, font=2)
abline(h=0, col="grey", lty=2)
abline(h=all.results.cl$Rho_LMA.N[which(all.results.cl$Taxo.Unit=="global")])
mtext(text = "c)", side = 3, adj=0, line=.2)
points(Rho_LMA.N~varNmass, all.results.cl[c(nrow(all.results.cl)-1,nrow(all.results.cl)),], pch=c(24,25), cex=1.4, bg=mypal[colchoices[c(length(colchoices)-1,length(colchoices))]])

plot(Rho_LL.N~varNmass, all.results.cl, pch=16, col=Type, xlab="Var. in %N", ylab="Rho   (LL v Nmass)")
abline(h=0, col="grey", lty=2)
abline(h=all.results.cl$Rho_LL.N[which(all.results.cl$Taxo.Unit=="global")])
legend('topright', legend = levels(all.results.cl$Type), pch=c(16,16,16,24,25), col=c(mypal[colchoices[1:3]], "black","black"), pt.bg= mypal[colchoices], bty ="n", cex=.7)
mtext(text = "d)", side = 3, adj=0, line=.2)
points(Rho_LL.N~varNmass, all.results.cl[c(nrow(all.results.cl)-1,nrow(all.results.cl)),], pch=c(24,25), cex=1.4, bg=mypal[colchoices[c(length(colchoices)-1,length(colchoices))]])

#

quartz(width=5.5,height=3)
par(mar=c(4,4,0,1), mfrow=c(1,2), mgp=c(2.5,1,0), oma=c(0,0,2,0))
palette(mypal[colchoices])
plot(Rho_LMA.Narea~varNarea, all.results.cl, pch=16, col=Type, xlab="Var. in %N", ylab="Rho   (LMA v Narea)")
#mtext( text="Narea v LMA", side=3, line=0, font=2)
abline(h=0, col="grey", lty=2)
abline(h=all.results.cl$Rho_LMA.Narea[which(all.results.cl$Taxo.Unit=="global")])
mtext(text = "c)", side = 3, adj=0, line=.2)
points(Rho_LMA.Narea~varNarea, all.results.cl[c(nrow(all.results.cl)-1,nrow(all.results.cl)),], pch=c(24,25), cex=1.4, bg=mypal[colchoices[c(length(colchoices)-1,length(colchoices))]])

plot(Rho_LL.Narea~varNarea, all.results.cl, pch=16, col=Type, xlab="Var. in %N", ylab="Rho   (LL v Narea)")
abline(h=0, col="grey", lty=2)
abline(h=all.results.cl$Rho_LL.Narea[which(all.results.cl$Taxo.Unit=="global")])
legend('bottomright', legend = levels(all.results.cl$Type), pch=c(16,16,16,24,25), col=c(mypal[colchoices[1:3]], "black","black"), pt.bg= mypal[colchoices], bty ="n", cex=.7)
mtext(text = "d)", side = 3, adj=0, line=.2)
points(Rho_LL.Narea~varNarea, all.results.cl[c(nrow(all.results.cl)-1,nrow(all.results.cl)),], pch=c(24,25), cex=1.4, bg=mypal[colchoices[c(length(colchoices)-1,length(colchoices))]])






############ Table of correlations ################
lma.n <- all.results.cl %>% group_by(Type) %>% filter(n_LMA.N>5) %>% summarise(cor.m = mean(Rho_LMA.N, na.rm=T), cor.se = sterr(Rho_LMA.N), slope.m = mean(Slope_LMA.N, na.rm=T), slope.se = sterr(Slope_LMA.N))
lma.ll <- all.results.cl %>% group_by(Type) %>% filter(n_LMA.LL>5) %>% summarise(cor.m = mean(Rho_LMA.LL, na.rm=T), cor.se = sterr(Rho_LMA.LL), slope.m = mean(Slope_LMA.LL, na.rm=T), slope.se = sterr(Slope_LMA.LL))
lma.narea <- all.results.cl %>% group_by(Type) %>% filter(n_LMA.Narea>5) %>% summarise(cor.m = mean(Rho_LMA.Narea, na.rm=T), cor.se = sterr(Rho_LMA.Narea), slope.m = mean(Slope_LMA.Narea, na.rm=T), slope.se = sterr(Slope_LMA.Narea))
ll.narea <- all.results.cl %>% group_by(Type) %>% filter(n_LL.Narea>5) %>% summarise(cor.m = mean(Rho_LL.Narea, na.rm=T), cor.se = sterr(Rho_LL.Narea), slope.m = mean(Slope_LL.Narea, na.rm=T), slope.se = sterr(Slope_LL.Narea))
n.ll <- all.results.cl %>% group_by(Type) %>% filter(n_N.LL>5) %>% summarise(cor.m = mean(Rho_N.LL, na.rm=T), cor.se = sterr(Rho_N.LL), slope.m = mean(Slope_N.LL, na.rm=T), slope.se = sterr(Slope_N.LL))

extract.values <- function (dataz, taxo, variable){
  if(is.na(dataz[taxo,paste(variable,"se", sep=".")])){
    value <- as.character(round(dataz[taxo, paste(variable,"m", sep=".")],digits = 2))
  }
  else{
  value <- paste(round(dataz[taxo,paste(variable, "m", sep=".")], digits=2), round(dataz[taxo, paste(variable, "se", sep=".")], digits=2), sep="±")
  }
  return(value)
}

lmas <- list(lma.ll, lma.n, lma.narea)
lls <- list(n.ll, ll.narea)


###### Correlations ########
c1 <- unlist(lapply(X=lmas, FUN = function(X){extract.values(dataz=X, taxo=1, variable="cor")}))
c2 <- unlist(lapply(X=lls, FUN = function(X){extract.values(dataz=X, taxo=1, variable="cor")}))
table.cor.spp <- data.frame(LMA = c1, LL = c("-", c2), row.names=c("LL","Nmass","Narea"))
c1 <- unlist(lapply(X=lmas, FUN = function(X){extract.values(dataz=X, taxo=4, variable="cor")}))
c2 <- unlist(lapply(X=lls, FUN = function(X){extract.values(dataz=X, taxo=4, variable="cor")}))
table.cor.btwfam <- data.frame(LMA = c1, LL = c("-", c2), row.names=c("LL","Nmass","Narea"))
c1 <- unlist(lapply(X=lmas, FUN = function(X){extract.values(dataz=X, taxo=5, variable="cor")}))
c2 <- unlist(lapply(X=lls, FUN = function(X){extract.values(dataz=X, taxo=5, variable="cor")}))
table.cor.global <- data.frame(LMA = c1, LL = c("-", c2), row.names=c("LL","Nmass","Narea"))

##### SMA Slopes #############
c1 <- unlist(lapply(X=lmas, FUN = function(X){extract.values(dataz=X, taxo=1, variable="slope")}))
c2 <- unlist(lapply(X=lls, FUN = function(X){extract.values(dataz=X, taxo=1, variable="slope")}))
table.slope.spp <- data.frame(LMA = c1, LL = c("-", c2), row.names=c("LL","Nmass","Narea"))
c1 <- unlist(lapply(X=lmas, FUN = function(X){extract.values(dataz=X, taxo=4, variable="slope")}))
c2 <- unlist(lapply(X=lls, FUN = function(X){extract.values(dataz=X, taxo=4, variable="slope")}))
table.slope.btwfam <- data.frame(LMA = c1, LL = c("-", c2), row.names=c("LL","Nmass","Narea"))
c1 <- unlist(lapply(X=lmas, FUN = function(X){extract.values(dataz=X, taxo=5, variable="slope")}))
c2 <- unlist(lapply(X=lls, FUN = function(X){extract.values(dataz=X, taxo=5, variable="slope")}))
table.slope.global <- data.frame(LMA = c1, LL = c("-", c2), row.names=c("LL","Nmass","Narea"))

write.csv(table.cor.spp, "./correlation_tables/Correlations_spp.csv")
write.csv(table.cor.btwfam, "./correlation_tables/Correlations_btwfam.csv")
write.csv(table.cor.global, "./correlation_tables/Correlations_global.csv")
write.csv(table.slope.spp, "./correlation_tables/Slopes_spp.csv")
write.csv(table.slope.btwfam, "./correlation_tables/Slopes_btwfam.csv")
write.csv(table.slope.global, "./correlation_tables/Slopes_global.csv")











######## test figure showing relaitonships at different levels

quartz(width=9, height= 3.5)
par(mar=c(1,1,1,1),oma=c(3,3,0,0), mgp=c(2,1,0), mfrow=c(1,3), cex=1.2)
plot(log.LMA~log.Narea, allspp, pch=16, cex=.8, col="grey")
points(log.LMA~log.Narea, spp.data[which(spp.data$Species %in% names(which(xtabs(~Species, spp.data)>20))),]
       , col=Species, pch=16, cex=.8)
xvar <- "log.Narea";yvar="log.LMA";meth=3; lwd=1
for(i in names(which(xtabs(~Species, spp.data)>20)) ){
 # plot.MAR(xvar = "log.Narea", yvar = "log.LMA",data= spp.data[which(spp.data$Species==i),], linecol = eval("data$Species") )

  tmp.mod <- suppressMessages(lmodel2(get(yvar)~get(xvar), spp.data[which(spp.data$Species==i),]))
  intercept <- tmp.mod$regression.results$Intercept[meth]
  slope <- tmp.mod$regression.results$Slope[meth]
  yhat <- tmp.mod$x * slope + intercept
  lines(yhat~tmp.mod$x, col="black")#col=spp.data$Species[which(spp.data$Species==i)][1], lwd=lwd)
  
}
mtext("log(LMA)",side=2, line=2.5)
abline(lmodel2(log.LMA~log.Narea, fam.data)$regression.results$Intercept[3], lmodel2(log.LMA~log.Narea, fam.data)$regression.results$Slope[3], lwd=3, lty=2 )


plot(log.LMA~log.Narea, allspp, pch=16, cex=.8, col="grey")
points(log.LMA~log.Narea, gen.data[which(gen.data$Genus %in% names(which(xtabs(~Genus, gen.data)>10))),]
       , col=Genus, pch=16, cex=.8)
xvar <- "log.Narea";yvar="log.LMA";meth=3; lwd=1
for(i in names(which(xtabs(~Genus, gen.data)>10)) ){
  # plot.MAR(xvar = "log.Narea", yvar = "log.LMA",data= gen.data[which(gen.data$Genus==i),], linecol = eval("data$Genus") )
  
  tmp.mod <- suppressMessages(lmodel2(get(yvar)~get(xvar), gen.data[which(gen.data$Genus==i),]))
  intercept <- tmp.mod$regression.results$Intercept[meth]
  slope <- tmp.mod$regression.results$Slope[meth]
  yhat <- tmp.mod$x * slope + intercept
  lines(yhat~tmp.mod$x, col="black")#col=gen.data$Genus[which(gen.data$Genus==i)][1], lwd=lwd)
  
}
mtext("log(Narea)",side=1, line=2.5)
abline(lmodel2(log.LMA~log.Narea, fam.data)$regression.results$Intercept[3], lmodel2(log.LMA~log.Narea, fam.data)$regression.results$Slope[3], lwd=3, lty=2 )


plot(log.LMA~log.Narea, allspp, pch=16, cex=.8, col="grey")
points(log.LMA~log.Narea, geninfam.data[which(geninfam.data$Family %in% names(which(xtabs(~Family, geninfam.data)>10))),]
       , col=Family, pch=16, cex=.8)
xvar <- "log.Narea";yvar="log.LMA";meth=3; lwd=1
for(i in names(which(xtabs(~Family, geninfam.data)>10)) ){
  # plot.MAR(xvar = "log.Narea", yvar = "log.LMA",data= geninfam.data[which(geninfam.data$Family==i),], linecol = eval("data$Family") )
  
  tmp.mod <- suppressMessages(lmodel2(get(yvar)~get(xvar), geninfam.data[which(geninfam.data$Family==i),]))
  intercept <- tmp.mod$regression.results$Intercept[meth]
  slope <- tmp.mod$regression.results$Slope[meth]
  yhat <- tmp.mod$x * slope + intercept
  lines(yhat~tmp.mod$x, col="black")#col=geninfam.data$Family[which(geninfam.data$Family==i)][1], lwd=lwd)
  
}

#plot.MAR(xvar="log.Narea", yvar="log.LMA", data= fam.data, linecol = "black", lwd=2)
abline(lmodel2(log.LMA~log.Narea, fam.data)$regression.results$Intercept[3], lmodel2(log.LMA~log.Narea, fam.data)$regression.results$Slope[3], lwd=3, lty=2 )
















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
p <- boxplot(Slope_LMA.LL~Type, all.results.cl[which(all.results.cl$n_LMA.LL>5 & !all.results.cl$Type %in% c("Famclean","global")),]
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
p <- boxplot(Rho_LMA.LL~Type, all.results.cl[which(all.results.cl$n_LMA.LL>5& !all.results.cl$Type %in% c("Famclean","global")),]
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
text(y=par()$usr[4]-.2, x=c(4,5), labels=paste0("(",all.results.cl$n_LMA.LL[c((nrow(all.results.cl)-1),nrow(all.results.cl))],")"),cex=.8)

par(mar=c(4,4,1.5,1.5), mgp=c(2.5,1,0))
plot(log.LL~log.LMA, allspp, col="grey", pch=16, ylab="Leaf Lifespan (months)", xlab="Leaf Mass per Area (g/m2)", xaxt="n", yaxt="n")
abline(a=all.results.cl$Int_LMA.LL[nrow(all.results.cl)], b=all.results.cl$Slope_LMA.LL[nrow(all.results.cl)], lwd=3, col="black", lty=3)
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
abline(a=all.results.cl$Int_LMA.LL[nrow(all.results.cl)], b=all.results.cl$Slope_LMA.LL[nrow(all.results.cl)], lwd=3, col="black", lty=3)
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
abline(a=all.results.cl$Int_LMA.LL[nrow(all.results.cl)], b=all.results.cl$Slope_LMA.LL[nrow(all.results.cl)], lwd=3, col="black", lty=3)
abline(a=all.results.cl$Int_LMA.LL[nrow(all.results.cl)-1], b=all.results.cl$Slope_LMA.LL[nrow(all.results.cl)-1], lwd=3, col="black")


##### Within Families ########
quartz(width=7.008/2, height=4)
par(mar=c(4,4,1.5,1.5))
## scatterplot
  # need 14 colors, cause I've got 14 genera
genpal <- c(brewer.pal(12,"Set3"), brewer.pal(3, "Set1"))
genpal <- c(rev(brewer.pal(9,"Set1")), brewer.pal(3, "Set2"))
#genpal <- c(brewer.pal(8,"Dark2"), brewer.pal(3, "Set1"))
plot(log.LL~log.LMA, allspp, col="grey", pch=16, ylab="log(LL)", xlab="log(LMA)")
abline(a=all.results.cl$Int_LMA.LL[nrow(all.results.cl)], b=all.results.cl$Slope_LMA.LL[nrow(all.results.cl)], lwd=3, col="black", lty=3)
abline(a=all.results.cl$Int_LMA.LL[nrow(all.results.cl)-1], b=all.results.cl$Slope_LMA.LL[nrow(all.results.cl)-1], lwd=3, col="black")
palette(paste0(genpal,"88"))
points(log.LL~log.LMA, geninfam.data[which(geninfam.data$Family %in% all.results.cl$Taxo.Unit[which(all.results.cl$n_LMA.LL>5 & all.results.cl$Type=="Genw.inFam")]),], pch=16, col=factor(Family))
palette(genpal)
tax <- "Genw.inFam"
for (j in 1:length(as.character(all.results.cl$Taxo.Unit[which(all.results.cl$Type==tax & all.results.cl$n_LMA.LL>5)]))){
  i <- as.character(all.results.cl$Taxo.Unit[which(all.results.cl$Type==tax & all.results.cl$n_LMA.LL>5)])[j]
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
abline(a=all.results.cl$Int_LMA.LL[nrow(all.results.cl)], b=all.results.cl$Slope_LMA.LL[nrow(all.results.cl)], lwd=3, col="black", lty=3)
abline(a=all.results.cl$Int_LMA.LL[nrow(all.results.cl)-1], b=all.results.cl$Slope_LMA.LL[nrow(all.results.cl)-1], lwd=3, col="black")
palette(paste0(genpal,"88"))
points(log.LL~log.LMA, gen.data[which(gen.data$Genus %in% all.results.cl$Taxo.Unit[which(all.results.cl$n_LMA.LL>5 & all.results.cl$Type=="w.inGen")]),], pch=16)#, col=factor(Genus))
palette(genpal)
tax <- "w.inGen"
for (j in 1:length(as.character(all.results.cl$Taxo.Unit[which(all.results.cl$Type==tax & all.results.cl$n_LMA.LL>5)]))){
  i <- as.character(all.results.cl$Taxo.Unit[which(all.results.cl$Type==tax & all.results.cl$n_LMA.LL>5)])[j]
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
abline(a=all.results.cl$Int_LMA.LL[nrow(all.results.cl)], b=all.results.cl$Slope_LMA.LL[nrow(all.results.cl)], lwd=3, col="black", lty=3)
abline(a=all.results.cl$Int_LMA.LL[nrow(all.results.cl)-1], b=all.results.cl$Slope_LMA.LL[nrow(all.results.cl)-1], lwd=3, col="black")
palette(paste0(genpal,"88"))
points(log.LL~log.LMA, spp.data[which(spp.data$Species %in% all.results.cl$Taxo.Unit[which(all.results.cl$n_LMA.LL>5 & all.results.cl$Type=="w.inSpp")]),], pch=16)#, col=factor(Species))
palette(genpal)
tax <- "w.inSpp"
for (j in 1:length(as.character(all.results.cl$Taxo.Unit[which(all.results.cl$Type==tax & all.results.cl$n_LMA.LL>5)]))){
  i <- as.character(all.results.cl$Taxo.Unit[which(all.results.cl$Type==tax & all.results.cl$n_LMA.LL>5)])[j]
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
p <- boxplot(Slope_LMA.LL~Type, all.results.cl[which(all.results.cl$n_LMA.LL>5 & !all.results.cl$Type %in% c("Famclean","global","w.inSpp","w.inGen","Genw.inFam")),]
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
p <- boxplot(Rho_LMA.LL~Type, all.results.cl[which(all.results.cl$n_LMA.LL>5& !all.results.cl$Type %in% c("Famclean","global","w.inSpp","w.inGen","Genw.inFam")),]
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
text(y=par()$usr[4]-.2, x=c(5), labels=paste0("(",all.results.cl$n_LMA.LL[nrow(all.results.cl)],")"),cex=.8)

par(mar=c(4,4,1.5,1.5), mgp=c(2.5,1,0))
plot(log.LL~log.LMA, allspp, col="grey", pch=16, ylab="Leaf Lifespan (months)", xlab="Leaf Mass per Area (g/m2)", xaxt="n", yaxt="n")
abline(a=all.results.cl$Int_LMA.LL[nrow(all.results.cl)], b=all.results.cl$Slope_LMA.LL[nrow(all.results.cl)], lwd=3, col="black", lty=3)
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

tmpdata <- all.results.cl[-which(all.results.cl$Type=="global"),]
tmpdata$Type <- factor(tmpdata$Type) # remove global level for boxplots
mat <- matrix(c(1,3,
                2,3), nrow=2, byrow = T)
colchoices <- c(1,2,4,3,6)
palette(mypal[colchoices])

quartz(width=7.008, height=4)
layout(mat)
par(mar=c(0,4,6,1))
p <- boxplot(Slope_LMA.LL~Type, tmpdata[which(tmpdata$n_LMA.LL>5 & !tmpdata$Type %in% c("Famclean","global","w.inSpp","w.inGen","Genw.inFam")),]
             , ylim=c(-2.5,4),las=3, ylab="Slope" #, main="log(LMA)~log(Nmass) strict"
             , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
             ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
             , staplelwd=0, outpch=16, outcex=.5, outcol=mypal[colchoices]
             , boxwex=.7, xaxt="n")
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
mtext(text = "a)", side = 3, adj=0, line=.2)
mtext(text= "LMA vs LL", side=3, line=.2)
par(mar=c(6,4,0,1))
p <- boxplot(Rho_LMA.LL~Type, tmpdata[which(tmpdata$n_LMA.LL>5& !tmpdata$Type %in% c("Famclean","global","w.inSpp","w.inGen","Genw.inFam")),]
             , ylim=c(-1,1.4),las=3, ylab="Correlation"
             , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
             ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
             , staplelwd=0, outpch=16, outcex=.5, outcol=mypal[colchoices]
             , boxwex=.7, xaxt="n")
# m <- cor.test(x=allspp$log.LL, y=allspp$log.LMA)
# points(y=m$estimate,x=5, pch=16, cex=1.3, col="darkgrey")
# arrows(x0=5,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3, col="darkgrey")
m <- cor.test(x=fam.dataclean$log.LL, y=fam.dataclean$log.LMA)
points(y=m$estimate,x=4, pch=16, cex=1.3, col="black")
arrows(x0=4,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3)


axis(1,labels = c("w/in Spp","w/in Gen", "w/in Fam","btw Fam","global"), at=c(1,2,3,4,5), las=3)
abline(h=0, lty=2)
# text(y=par()$usr[4]-.2, x=c(1,2,3), labels = p$n[1:3])
text(y=par()$usr[4]-.2, x=c(4,5), labels=paste0("(",all.results.cl$n_LMA.LL[c((nrow(all.results.cl)-1),nrow(all.results.cl))],")"),cex=.8)

par(mar=c(4,4,1.5,1.5), mgp=c(2.5,1,0))
plot(log.LL~log.LMA, allspp,type="n", col="grey", pch=16, ylab="Leaf Lifespan (months)", xlab="Leaf Mass per Area (g/m2)", xaxt="n", yaxt="n")
#abline(a=all.results.cl$Int_LMA.LL[nrow(all.results.cl)], b=all.results.cl$Slope_LMA.LL[nrow(all.results.cl)], lwd=3, col="black", lty=3)

#points(log.LL~log.LMA, allspp[which(allspp$Species=="Solanum straminifolia"),])
#points(log.LL~log.LMA, allspp[which(allspp$Species=="Juniperus monosperma"),])
#points(log.LL~log.LMA, allspp[which(allspp$Species=="Araucaria araucana"),])
#points(log.LL~log.LMA, allspp[which(allspp$Species=="Abies magnifica"),])
axis(2, at=log(c(1,2,3,4,5,10,20,30,40,50,100,200,300), base=10), labels = F)#c(1,"","","",5,"","","","",50,"",200,""),)
axis(2, at=log(c(1,5,50,300), base=10), labels = c(1,5,50,300))
axis(1, at=log(c(10,20,30,40,50,100,200,300,400,500,1000), base=10), labels = F)# c(10,"","","",50,"","","","","","1000"))
axis(1, at=log(c(20,50,500), base=10), labels =  c(20,50,500))
points(log.LL~log.LMA, fam.dataclean, pch=16)
abline(a=all.results.cl$Int_LMA.LL[nrow(all.results.cl)-1], b=all.results.cl$Slope_LMA.LL[nrow(all.results.cl)-1], lwd=3, col="black")









#________________________________________________________________
###### Full Within Families only ######
#________________________________________________________________

mat <- matrix(c(1,3,
                2,3), nrow=2, byrow = T)
colchoices <- c(1,2,4,3,6)
palette(mypal[colchoices])

quartz(width=7.008, height=4)
layout(mat)
par(mar=c(0,4,6,1))
p <- boxplot(Slope_LMA.LL~Type, tmpdata[which(tmpdata$n_LMA.LL>5 & !tmpdata$Type %in% c("Famclean","global","w.inSpp","w.inGen")),]
             , ylim=c(-2.5,4),las=3, ylab="Slope" #, main="log(LMA)~log(Nmass) strict"
             , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
             ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
             , staplelwd=0, outpch=16, outcex=.5, outcol=mypal[colchoices]
             , boxwex=.7, xaxt="n")
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
mtext(text = "a)", side = 3, adj=0, line=.2)
mtext(text= "LMA vs LL", side=3, line=.2)
par(mar=c(6,4,0,1))
p <- boxplot(Rho_LMA.LL~Type, tmpdata[which(tmpdata$n_LMA.LL>5& !tmpdata$Type %in% c("Famclean","global","w.inSpp","w.inGen")),]
             , ylim=c(-1,1.4),las=3, ylab="Correlation"
             , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
             ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
             , staplelwd=0, outpch=16, outcex=.5, outcol=mypal[colchoices]
             , boxwex=.7, xaxt="n")
# m <- cor.test(x=allspp$log.LL, y=allspp$log.LMA)
# points(y=m$estimate,x=5, pch=16, cex=1.3, col="darkgrey")
# arrows(x0=5,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3, col="darkgrey")
m <- cor.test(x=fam.dataclean$log.LL, y=fam.dataclean$log.LMA)
points(y=m$estimate,x=4, pch=16, cex=1.3, col="black")
arrows(x0=4,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3)


axis(1,labels = c("w/in Spp","w/in Gen", "w/in Fam","btw Fam","global"), at=c(1,2,3,4,5), las=3)
abline(h=0, lty=2)
text(y=par()$usr[4]-.2, x=c(3), labels = p$n[3])
# text(y=par()$usr[4]-.2, x=c(4,5), labels=paste0("(",tmpdata$n_LMA.LL[c((nrow(tmpdata)-1),nrow(tmpdata))],")"),cex=.8)
text(y=par()$usr[4]-.2, x=c(4,5), labels=paste0("(",tmpdata$n_LMA.LL[nrow(tmpdata)],")"),cex=.9)

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
layout(mat)
par(mar=c(0,4,6,1))
p <- boxplot(Slope_LMA.LL~Type, tmpdata[which(tmpdata$n_LMA.LL>5 & !tmpdata$Type %in% c("Famclean","global","w.inSpp")),]
             , ylim=c(-2.5,4),las=3, ylab="Slope" #, main="log(LMA)~log(Nmass) strict"
             , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
             ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
             , staplelwd=0, outpch=16, outcex=.5, outcol=mypal[colchoices]
             , boxwex=.7, xaxt="n")
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
mtext(text = "a)", side = 3, adj=0, line=.2)
mtext(text= "LMA vs LL", side=3, line=.2)
par(mar=c(6,4,0,1))
p <- boxplot(Rho_LMA.LL~Type, tmpdata[which(tmpdata$n_LMA.LL>5& !tmpdata$Type %in% c("Famclean","global","w.inSpp")),]
             , ylim=c(-1,1.4),las=3, ylab="Correlation"
             , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
             ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
             , staplelwd=0, outpch=16, outcex=.5, outcol=mypal[colchoices]
             , boxwex=.7, xaxt="n")
# m <- cor.test(x=allspp$log.LL, y=allspp$log.LMA)
# points(y=m$estimate,x=5, pch=16, cex=1.3, col="darkgrey")
# arrows(x0=5,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3, col="darkgrey")
m <- cor.test(x=fam.dataclean$log.LL, y=fam.dataclean$log.LMA)
points(y=m$estimate,x=4, pch=16, cex=1.3, col="black")
arrows(x0=4,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3)


axis(1,labels = c("w/in Spp","w/in Gen", "w/in Fam","btw Fam","global"), at=c(1,2,3,4,5), las=3)
abline(h=0, lty=2)
text(y=par()$usr[4]-.2, x=c(2,3), labels = p$n[2:3])
# text(y=par()$usr[4]-.2, x=c(4,5), labels=paste0("(",tmpdata$n_LMA.LL[c((nrow(tmpdata)-1),nrow(tmpdata))],")"),cex=.8)
text(y=par()$usr[4]-.2, x=c(4,5), labels=paste0("(",tmpdata$n_LMA.LL[nrow(tmpdata)],")"),cex=.8)

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
layout(mat)
par(mar=c(0,4,6,1))
p <- boxplot(Slope_LMA.LL~Type, tmpdata[which(tmpdata$n_LMA.LL>5 & !tmpdata$Type %in% c("Famclean","global")),]
             , ylim=c(-2.5,4),las=3, ylab="Slope" #, main="log(LMA)~log(Nmass) strict"
             , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
             ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
             , staplelwd=0, outpch=16, outcex=.5, outcol=mypal[colchoices]
             , boxwex=.7, xaxt="n")
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
mtext(text = "a)", side = 3, adj=0, line=.2)
mtext(text= "LMA vs LL", side=3, line=.2)
par(mar=c(6,4,0,1))
p <- boxplot(Rho_LMA.LL~Type, tmpdata[which(tmpdata$n_LMA.LL>5& !tmpdata$Type %in% c("Famclean","global")),]
             , ylim=c(-1,1.4),las=3, ylab="Correlation"
             , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
             ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
             , staplelwd=0, outpch=16, outcex=.5, outcol=mypal[colchoices]
             , boxwex=.7, xaxt="n")
# m <- cor.test(x=allspp$log.LL, y=allspp$log.LMA)
# points(y=m$estimate,x=5, pch=16, cex=1.3, col="darkgrey")
# arrows(x0=5,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3, col="darkgrey")
m <- cor.test(x=fam.dataclean$log.LL, y=fam.dataclean$log.LMA)
points(y=m$estimate,x=4, pch=16, cex=1.3, col="black")
arrows(x0=4,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3)


axis(1,labels = c("w/in Spp","w/in Gen", "w/in Fam","btw Fam","global"), at=c(1,2,3,4,5), las=3)
abline(h=0, lty=2)
text(y=par()$usr[4]-.2, x=c(1,2,3), labels = p$n[1:3])
# text(y=par()$usr[4]-.2, x=c(4,5), labels=paste0("(",tmpdata$n_LMA.LL[c((nrow(tmpdata)-1),nrow(tmpdata))],")"),cex=.8)
text(y=par()$usr[4]-.2, x=c(4,5), labels=paste0("(",tmpdata$n_LMA.LL[nrow(tmpdata)],")"),cex=.8)

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

tax <- "Genw.inFam"
for (j in 1:length(as.character(tmpdata$Taxo.Unit[which(tmpdata$Type==tax & tmpdata$n_LMA.LL>5)]))){
  i <- as.character(tmpdata$Taxo.Unit[which(tmpdata$Type==tax & tmpdata$n_LMA.LL>5)])[j]
  # different colored lines
  # plot.MAR(xvar = "log.LMA", yvar = "log.LL",data= geninfam.data[which(geninfam.data$Family==i),], linecol = genpal[j], lwd=3)
  # same colored lines
  plot.MAR(xvar = "log.LMA", yvar = "log.LL",data= geninfam.data[which(geninfam.data$Family==i),], linecol = mypal[colchoices[3]], lwd=1)
}

tax <- "w.inGen"
for (j in 1:length(as.character(tmpdata$Taxo.Unit[which(tmpdata$Type==tax & tmpdata$n_LMA.LL>5)]))){
  i <- as.character(tmpdata$Taxo.Unit[which(tmpdata$Type==tax & tmpdata$n_LMA.LL>5)])[j]
  # different colored lines
  # plot.MAR(xvar = "log.LMA", yvar = "log.LL",data= geninfam.data[which(geninfam.data$Family==i),], linecol = genpal[j], lwd=3)
  # same colored lines
  plot.MAR(xvar = "log.LMA", yvar = "log.LL",data= gen.data[which(gen.data$Genus==i),], linecol = mypal[colchoices[2]], lwd=1)
}




tax <- "w.inSpp"
for (j in 1:length(as.character(tmpdata$Taxo.Unit[which(tmpdata$Type==tax & tmpdata$n_LMA.LL>5)]))){
  i <- as.character(tmpdata$Taxo.Unit[which(tmpdata$Type==tax & tmpdata$n_LMA.LL>5)])[j]
  #plot.MAR(xvar = "log.LMA", yvar = "log.LL",data= spp.data[which(spp.data$Species==i),], linecol = genpal[j], lwd=3)
  plot.MAR(xvar = "log.LMA", yvar = "log.LL",data= spp.data[which(spp.data$Species==i),], linecol = mypal[colchoices[1]], lwd=3)
  
}




#________________________________________________________________
###### LEAF LIFE ~ ENVIRONMENT example plot ######
#________________________________________________________________

palette(paste(brewer.pal(n=6, "Set1"), "77",sep=""))
tmpdata <- traits.common[which(traits.common$SP.ID %in% c("PSEMEN","PINPON","PINJEF","PINCON","TSUHET","ABICON")),]
tmpdata$SP.ID <- factor(tmpdata$SP.ID)
plot(log.LL~climPC1,tmpdata, col=SP.ID, pch=16)

mod <- lmer(log.LL~climPC2 + (climPC2|SP.ID),tmpdata , REML=T)

f0LL <- predict(mod, re.form=NA)
f1LL <- fitted(mod)
# sort lma values so lines draw right
I <- order(tmpdata$climPC2[-which(is.na(tmpdata$log.LL) | is.na(tmpdata$climPC2))])
LLs <- sort(tmpdata$climPC2[-which(is.na(tmpdata$log.LL) | is.na(tmpdata$climPC2))])
spid <- tmpdata$SP.ID[-which(is.na(tmpdata$log.LL) | is.na(tmpdata$climPC2))][I]
plot(log.LL~climPC2, tmpdata, col=SP.ID, pch=16, cex=.9, xlab="")#log(Leaf Lifespan)")
mtext("Relative Growth Rate", side=2, line=2.5, cex=1.1)
palette(brewer.pal(n=6, "Set1"))
for(i in levels(spid)){
  lines(LLs[which(spid==i)], f1LL[I][which(spid==i)],  col=spid[which(spid==i)], lwd=2)
}

for(i in levels(spid)){
  modsimp <- lm(log.LL~climPC2, tmpdata[which(tmpdata$SP.ID==i),])
  lines(LLs[which(spid==i)], f1LL[I][which(spid==i)],  col=spid[which(spid==i)], lwd=2)
}







#________________________________________________________________
###### Trait Covariation for Presentations ######
#________________________________________________________________


##### LMA v Nmass
quartz(width=3.5, height=3.5)

par(mar=c(4,4,1.5,1.5), mgp=c(2.2,.7,0), cex.lab=1.3, cex.axis=1.1)
plot(log.Nmass~log.LMA, LES, col="grey", pch=16, xlab=expression(paste(log[10],"(LMA)")), ylab=expression(paste(log[10],(N[mass]), sep=" ")))
tax <- "w.inGen"
for (i in as.character(all.results.cl$Taxo.Unit[which(all.results.cl$Type==tax & all.results.cl$n_LMA.N>5)])){
  plot.MAR(xvar = "log.LMA", yvar = "log.Nmass",data= gen.data[which(gen.data$Genus==i),], linecol = mypal[colchoices[2]])
}
tax <- "Genw.inFam"
for (i in as.character(all.results.cl$Taxo.Unit[which(all.results.cl$Type==tax & all.results.cl$n_LMA.N>5)])){
  plot.MAR(xvar = "log.LMA", yvar = "log.Nmass",data= geninfam.data[which(geninfam.data$Family==i),], linecol = mypal[colchoices[3]])
}
abline(a=all.results.cl$Int_LMA.N[nrow(all.results.cl)-1], b=all.results.cl$Slope_LMA.N[nrow(all.results.cl)-1], lwd=3, col="black")
abline(a=all.results.cl$Int_LMA.N[nrow(all.results.cl)], b=all.results.cl$Slope_LMA.N[nrow(all.results.cl)], lwd=3, col="black", lty=3)
tax <- "w.inSpp"
for (i in as.character(all.results.cl$Taxo.Unit[which(all.results.cl$Type==tax & all.results.cl$n_LMA.N>5)])){
  plot.MAR(xvar = "log.LMA", yvar = "log.Nmass",data= spp.data[which(spp.data$Species==i),], linecol = mypal[colchoices[1]])
}
mtext(text = expression(paste("A) Scale Independent: LMA-",N[mass])), side = 3, adj=0.9, line=.2, cex=1.3)
legend ("topright", bty="n", legend=c("w/in Spp", "w/in Gen", "w/in Fam"), lwd=c(2,2,2), lty=c(1,1,1), col=c(mypal[colchoices[c(1,2,3)]]), cex=.9)
legend("bottomleft", bty="n", legend=c("Btw Fams","Global"), lwd=3, lty=c(1,3), cex=1)




##### LMA v Narea
quartz(width=3.5, height=3.5)

par(mar=c(4,4,1.5,1.5), mgp=c(2.2,.7,0), cex.lab=1.3, cex.axis=1.1)
plot(log.Narea~log.LMA, LES, col="grey", pch=16, ylab=expression(paste(log[10],(N[area]))), xlab=expression(paste(log[10],"(LMA)")), ylim=c(-.6,1.1))

tax <- "w.inGen"
for (i in as.character(all.results.cl$Taxo.Unit[which(all.results.cl$Type==tax & all.results.cl$n_LMA.Narea>5 & all.results.cl$Slope_LMA.Narea>-1)])){
  plot.MAR(xvar = "log.LMA", yvar = "log.Narea",data= gen.data[which(gen.data$Genus==i),], linecol = mypal[2])
}

tax <- "Genw.inFam"
for (i in as.character(all.results.cl$Taxo.Unit[which(all.results.cl$Type==tax & all.results.cl$n_LMA.Narea>5)])){
  plot.MAR(xvar = "log.LMA", yvar = "log.Narea",data= geninfam.data[which(geninfam.data$Family==i),], linecol = mypal[4])
}

abline(a=all.results.cl$Int_LMA.Narea[nrow(all.results.cl)-1], b=all.results.cl$Slope_LMA.Narea[nrow(all.results.cl)-1], lwd=3, col="black")
abline(a=all.results.cl$Int_LMA.Narea[nrow(all.results.cl)], b=all.results.cl$Slope_LMA.Narea[nrow(all.results.cl)], lwd=3, col="black", lty=3)

# plot.MAR(xvar="log.LMA", yvar="log.Narea", data= fam.data, linecol = mypal[5], lwd=2)
# 
# plot.MAR(xvar="log.LMA", yvar="log.Narea", data= LES, linecol = "black", lwd=2)
tax <- "w.inSpp"
for (i in as.character(all.results.cl$Taxo.Unit[which(all.results.cl$Type==tax & all.results.cl$n_LMA.Narea>5& all.results.cl$Slope_LMA.Narea>-1)])){
  plot.MAR(xvar = "log.LMA", yvar = "log.Narea",data= spp.data[which(spp.data$Species==i),], linecol = mypal[1])
}
mtext(text = expression(paste("B) Stronger w/in spp: LMA-",N[area])), side = 3, adj=0.7, line=.2, cex=1.3)
#legend ("topright", bty="n", legend=c("w/in Spp", "w/in Gen", "w/in Fam"), lwd=c(2,2,2), lty=c(1,1,1), col=c(mypal[colchoices[c(1,2,3)]]), cex=.9)
#legend("bottomleft", bty="n", legend=c("Btw Fams","Global"), lwd=3, lty=c(1,3), cex=.9)





##### LMA v LL
quartz(width=3.5, height=3.5)

par(mar=c(4,4,1.5,1.5), mgp=c(2.2,.7,0), cex.lab=1.3, cex.axis=1.1)
plot(log.LL~log.LMA, LES, col="grey", pch=16, ylab=expression(paste(log[10],"(LL)")), xlab=expression(paste(log[10],"(LMA)")))

tax <- "w.inGen"
for (i in as.character(all.results.cl$Taxo.Unit[which(all.results.cl$Type==tax & all.results.cl$n_LMA.LL>5)])){
  plot.MAR(xvar = "log.LMA", yvar = "log.LL",data= gen.data[which(gen.data$Genus==i),], linecol = mypal[colchoices[2]])
}

tax <- "Genw.inFam"
for (i in as.character(all.results.cl$Taxo.Unit[which(all.results.cl$Type==tax & all.results.cl$n_LMA.LL>5)])){
  plot.MAR(xvar = "log.LMA", yvar = "log.LL",data= geninfam.data[which(geninfam.data$Family==i),], linecol = mypal[colchoices[3]])
}


abline(a=all.results.cl$Int_LMA.LL[nrow(all.results.cl)-1], b=all.results.cl$Slope_LMA.LL[nrow(all.results.cl)-1], lwd=3, col="black")
abline(a=all.results.cl$Int_LMA.LL[nrow(all.results.cl)], b=all.results.cl$Slope_LMA.LL[nrow(all.results.cl)], lwd=3, col="black", lty=3)
tax <- "w.inSpp"
for (i in as.character(all.results.cl$Taxo.Unit[which(all.results.cl$Type==tax & all.results.cl$n_LMA.LL>5)])){
  plot.MAR(xvar = "log.LMA", yvar = "log.LL",data= spp.data[which(spp.data$Species==i),], linecol = mypal[colchoices[1]])
}

mtext(text = expression(paste("C) Switch sign w/in spp: LMA-LL")), side = 3, adj=0.7, line=.2, cex=1.3)
#legend ("topright", bty="n", legend=c("w/in Spp", "w/in Gen", "w/in Fam"), lwd=c(2,2,2), lty=c(1,1,1), col=c(mypal[colchoices[c(1,2,3)]]), cex=.9)
#legend("bottomleft", bty="n", legend=c("Btw Fams","Global"), lwd=3, lty=c(1,3), cex=.9)



###### LMA vs LL boxplots
quartz(width=2.75, height=3.5)
par(mfrow=c(2,1), mar=c(0,4,0,1), oma=c(4,0,2,0),mgp=c(2.2,.7,0), cex.lab=1.4, cex.axis=1.2)
p <- boxplot(Slope_LMA.LL~Type, all.results.cl[which(all.results.cl$n_LMA.LL>5& !all.results.cl$Type %in% c("global","Famclean")),]
             , ylim=c(-2.5,4),las=3, ylab="MA Slope" #, main="log(LMA)~log(Nmass) strict"
             , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
             ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
             , staplelwd=0, outpch=16, outcex=.7, outcol=mypal[colchoices]
             , boxwex=.7, xaxt="n")
abline(h=0, lty=2)
m <- lmodel2(log.LL~log.LMA, allspp)
points(y=m$regression.results$Slope[3],x=5, pch=16, cex=1.3, col="darkgrey")
arrows(x0=5,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`2.5%-Slope`[3], length = 0,lwd=3, col="darkgrey")
arrows(x0=5,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`97.5%-Slope`[3], length = 0,lwd=3, col="darkgrey")
m <- lmodel2(log.LL~log.LMA, fam.dataclean)
points(y=m$regression.results$Slope[3],x=4, pch=16, cex=1.3, col="black")
arrows(x0=4,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`2.5%-Slope`[3], length = 0,lwd=3)
arrows(x0=4,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`97.5%-Slope`[3], length = 0,lwd=3)

#text(y=par()$usr[3]+.2, x=c(1,2,3,4,5,6), labels = p$n)
#text(y=par()$usr[4]-.2,x=.5, labels = "a)")
mtext(text = "D)", side = 3, adj=0, line=.2, cex=1.4)
mtext(text= "LMA vs LL", side=3, line=.2, cex=1.4)
# par(mar=c(6,4,0,1))
p <- boxplot(Rho_LMA.LL~Type, all.results.cl[which(all.results.cl$n_LMA.LL>5& !all.results.cl$Type %in% c("global","Famclean")),]
             , ylim=c(-1,1.4),las=3, ylab="Rho"
             , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
             ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
             , staplelwd=0, outpch=16, outcex=.7, outcol=mypal[colchoices]
             , boxwex=.7, xaxt="n")
axis(1,labels = c("","", "","",""), at=c(1,2,3,4,5), las=3,tck=0.075)
#mtext(text=c("w/in Spp","w/in Gen", "w/in Fam","btw Fam","global"), side = 1,at = c(1,2,3,4,5), pos=2, srt=50)
text(x = c(1,2,3,4,5), y= par("usr")[3]-.2,labels = c("w/in Spp","w/in Gen", "w/in Fam","btw Fam","global"), srt=45, adj=1,xpd=NA, cex=1.4, font=1)

abline(h=0, lty=2)
text(y=par()$usr[4]-.2, x=c(1,2,3), labels = p$n[1:3])
text(y=par()$usr[4]-.25, x=c(4,5), labels=paste0("(",all.results.cl$n_LMA.LL[c((nrow(all.results.cl)-1),nrow(all.results.cl))],")"),cex=1)
#polygon(x = c(1,2,3,3,2,1), y=na.omit(c(meancis$m10_LMA.LL, rev(meancis$m90_LMA.LL))), border="#EEEEEE",col="#EEEEEE")#lightgrey",col = "lightgrey")
m <- cor.test(x=allspp$log.LMA, y=allspp$log.LL)
points(y=m$estimate,x=5, pch=16, cex=1.3, col="darkgrey")
arrows(x0=5,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3, col="darkgrey")
m <- cor.test(x=fam.dataclean$log.LMA, y=fam.dataclean$log.LL)
points(y=m$estimate,x=4, pch=16, cex=1.3, col="black")
arrows(x0=4,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3)








#________________________________________________________________
###### Trait Env Examples for Presentation ######
#________________________________________________________________


quartz(width=3.5, height=3.5)
par(mar=c(3.5,3.5,1,1), mgp=c(2,.7,0))
plot(log.LL~climPC2 , traits.common[which(traits.common$SP.ID %in% c("PSEMEN","PINPON","PINCON","PINJEF","TSUHET","ABICON")),], pch=16, col="grey", xlim=c(-3.6,3.6), ylim=c(1.38,2.35), xlab="Plot Warmth (climate PC2)", ylab=expression(paste(log[10](LL))))
for (i in c("PSEMEN","PINPON","PINCON","PINJEF","TSUHET","ABICON")){
  tmpmod <- lm(log.LL~climPC2, traits.common[which(traits.common$SP.ID==i),])
  #lines(fitted(tmpmod)~traits.common$climPC2[which(traits.common$SP.ID==i & !is.na(traits.common$log.LL))], lwd=2, col=mypal[colchoices[1]])
  lty <- ifelse(summary(tmpmod)$coefficients[2,4]>0.05,3,1)
  lines(fitted(tmpmod)[order(fitted(tmpmod))]~tmpmod$model$climPC2[order(fitted(tmpmod))], lwd=2, lty=lty, col=mypal[colchoices[1]])
}
legend(x=.2, y=2.4, legend="w/in spp", lwd=2, col=mypal[colchoices[1]], bty="n")
legend(x=.7, y=2.4, legend="", pch=16, col="grey", bty="n")

points(log.cw_LLp_if~climPC2, biomass, pch=17, col="black")
abline(lm(log.cw_LLp_if~climPC2, biomass), lwd=3)

legend(x=.2, y=2.33, legend="CWM", lwd=2, col="black", bty="n")
legend(x=.7, y=2.33, legend="", pch=17, col="black", bty="n")



quartz(width=3.5, height=3.5)
par(mar=c(3.5,3.5,1,1), mgp=c(2,.7,0))

plot(log.LMA~climPC2 , traits.common[which(traits.common$SP.ID %in% c("PSEMEN","PINPON","PINCON","PINJEF","TSUHET","ABICON")),], pch=16, col="grey", xlim=c(-3.5,3.5), ylim=c(1.8,2.8))
for (i in c("PSEMEN","PINPON","PINCON","PINJEF","TSUHET","ABICON")){
  tmpmod <- lm(log.LMA~climPC2, traits.common[which(traits.common$SP.ID==i),])
  lty <- ifelse(summary(tmpmod)$coefficients[2,4]>0.05,3,1)
lines(fitted(tmpmod)[order(fitted(tmpmod))]~tmpmod$model$climPC2[order(fitted(tmpmod))], lwd=2, lty=lty, col=mypal[colchoices[1]])
#  print(summary(tmpmod))
  }
#legend(x=.4, y=3, legend="w/in spp", lwd=2, col=mypal[colchoices[1]], bty="n")
#legend(x=.8, y=2.6, legend="", pch=16, col="grey", bty="n")


points(log.cw_LMAp_if~climPC2, biomass, pch=17, col="black")
abline(lm(log.cw_LMAp_if~climPC2, biomass), lwd=3)
#legend(x=.4, y=2.5, legend="CW mean", lwd=2, col="black", bty="n")
#legend(x=.8, y=2.5, legend="", pch=17, col="black", bty="n")






#________________________________________________________________
###### Null Model Plots - for supplement ######
#________________________________________________________________



meancis <- all.results.cl %>% group_by(Type) %>% summarise(m2.5_LMA.LL = mean(lci_2.5_LMA.LL, na.rm=T), m5_LMA.LL = mean(lci_5_LMA.LL, na.rm=T), m10_LMA.LL = mean(lci_10_LMA.LL, na.rm=T),
                                                           m97.5_LMA.LL = mean(uci_2.5_LMA.LL, na.rm=T), m95_LMA.LL = mean(uci_5_LMA.LL, na.rm=T), m90_LMA.LL = mean(uci_10_LMA.LL, na.rm=T),
                                                           m2.5_LMA.N = mean(lci_2.5_LMA.N, na.rm=T), m5_LMA.N = mean(lci_5_LMA.N, na.rm=T), m10_LMA.N = mean(lci_10_LMA.N, na.rm=T),
                                                           m97.5_LMA.N = mean(uci_2.5_LMA.N, na.rm=T), m95_LMA.N = mean(uci_5_LMA.N, na.rm=T), m90_LMA.N = mean(uci_10_LMA.N, na.rm=T),
                                                           m2.5_LL.N = mean(lci_2.5_LL.N, na.rm=T), m5_LL.N = mean(lci_5_LL.N, na.rm=T), m10_LL.N = mean(lci_10_LL.N, na.rm=T),
                                                           m97.5_LL.N = mean(uci_2.5_LL.N, na.rm=T), m95_LL.N = mean(uci_5_LL.N, na.rm=T), m90_LL.N = mean(uci_10_LL.N, na.rm=T),
                                                           m2.5_LMA.Narea = mean(lci_2.5_LMA.Narea, na.rm=T), m5_LMA.Narea = mean(lci_5_LMA.Narea, na.rm=T), m10_LMA.Narea = mean(lci_10_LMA.Narea, na.rm=T),
                                                           m97.5_LMA.Narea = mean(uci_2.5_LMA.Narea, na.rm=T), m95_LMA.Narea = mean(uci_5_LMA.Narea, na.rm=T), m90_LMA.Narea = mean(uci_10_LMA.Narea, na.rm=T),
                                                           m2.5_LL.Narea = mean(lci_2.5_LL.Narea, na.rm=T), m5_LL.Narea = mean(lci_5_LL.Narea, na.rm=T), m10_LL.Narea = mean(lci_10_LL.Narea, na.rm=T),
                                                           m97.5_LL.Narea = mean(uci_2.5_LL.Narea, na.rm=T), m95_LL.Narea = mean(uci_5_LL.Narea, na.rm=T), m90_LL.Narea = mean(uci_10_LL.Narea, na.rm=T))
# I think the easiest thing will be to make two quartez: top with boxplots and bottom with funnel
colchoices <- c(1,2,4,3,6)
palette(mypal[colchoices])

xvals <- as.numeric(all.results.cl$Type)
pchs <- rep(1, times=nrow(all.results.cl))
pchs[which(all.results.cl$sig_LMA.N>0)] <- 16
crit <-.1
### Lower boxplots with Rho

quartz(width=5, height=6)
par(mfrow=c(3,2), mar=c(5,3,1,1), mgp=c(2,1,0), cex=.8)
# LMA v Nmass
p <- boxplot(Rho_LMA.N~Type, all.results.cl[which(all.results.cl$n_LMA.N>5 & is.na(all.results.cl$Taxo.Unit)),]
             , ylim=c(-1,1.3),las=3, ylab="Rho"
             , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
             ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
             , staplelwd=0, outpch=16, outcex=.5, outcol=mypal[colchoices]
             , boxwex=.7, xaxt="n")
#lines(meancis$m10_LMA.N~c(1,2,3,4,5))
#lines(meancis$m90_LMA.N~c(1,2,3,4,5))
polygon(x = c(.8,2,3.2,3.2,2,.8), y=na.omit(c(meancis$m10_LMA.N, rev(meancis$m90_LMA.N))), border="#DDDDDD",col="#DDDDDD")#lightgrey",col = "lightgrey")
m <- cor.test(x=allspp$log.LMA, y=allspp$log.Nmass)
points(y=m$estimate,x=5, pch=15, cex=1.3, col="darkgrey")
arrows(x0=5,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3, col="darkgrey")
m <- cor.test(x=fam.dataclean$log.LMA, y=fam.dataclean$log.Nmass)
points(y=m$estimate,x=4, pch=17, cex=1.3, col="black")
arrows(x0=4,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3)
# layering points on top of null model
points(Rho_LMA.N~jitter(as.numeric(Type)), all.results.cl[which(all.results.cl$n_LMA.N>5 & all.results.cl$sig_LMA.N<crit),], pch=16)
points(Rho_LMA.N~jitter(as.numeric(Type)), all.results.cl[which(all.results.cl$n_LMA.N>5 & all.results.cl$sig_LMA.N>=crit),], pch=1)
axis(1,labels = c("w/in Spp","w/in Gen", "w/in Fam","btw Fam","global"), at=c(1,2,3,4,5), las=3)
abline(h=0, lty=2)
text(y=par()$usr[4]-.2,x=3, labels = "a) LMA vs Nmass")

# LL v Nmass
p <- boxplot(Rho_LL.N~Type, all.results.cl[which(all.results.cl$n_LL.N>5 & is.na(all.results.cl$Taxo.Unit)),]
             , ylim=c(-1,1.3),las=3, ylab="Rho"
             , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
             ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
             , staplelwd=0, outpch=16, outcex=.5, outcol=mypal[colchoices]
             , boxwex=.7, xaxt="n")
#lines(meancis$m10_LL.N~c(1,2,3,4,5))
#lines(meancis$m90_LL.N~c(1,2,3,4,5))
polygon(x = c(.8,2,3.2,3.2,2,.8), y=na.omit(c(meancis$m10_LL.N, rev(meancis$m90_LL.N))), border="#DDDDDD",col="#DDDDDD")#lightgrey",col = "lightgrey")
m <- cor.test(x=allspp$log.LL, y=allspp$log.Nmass)
points(y=m$estimate,x=5, pch=15, cex=1.3, col="darkgrey")
arrows(x0=5,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3, col="darkgrey")
m <- cor.test(x=fam.dataclean$log.LL, y=fam.dataclean$log.Nmass)
points(y=m$estimate,x=4, pch=17, cex=1.3, col="black")
arrows(x0=4,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3)
# layering points on top of null model
points(Rho_LL.N~jitter(as.numeric(Type)), all.results.cl[which(all.results.cl$n_LL.N>5 & all.results.cl$sig_LL.N<crit),], pch=16)
points(Rho_LL.N~jitter(as.numeric(Type)), all.results.cl[which(all.results.cl$n_LL.N>5 & all.results.cl$sig_LL.N>=crit),], pch=1)
axis(1,labels = c("w/in Spp","w/in Gen", "w/in Fam","btw Fam","global"), at=c(1,2,3,4,5), las=3)
abline(h=0, lty=2)
#text(y=par()$usr[4]-.2,x=.5, labels = "b)")
text(y=par()$usr[4]-.2,x=3, labels = "b) Leaf Life vs Nmass")


# LMA v Narea
p <- boxplot(Rho_LMA.Narea~Type, all.results.cl[which(all.results.cl$n_LMA.Narea>5 & is.na(all.results.cl$Taxo.Unit)),]
             , ylim=c(-1,1.3),las=3, ylab="Rho"
             , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
             ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
             , staplelwd=0, outpch=16, outcex=.5, outcol=mypal[colchoices]
             , boxwex=.7, xaxt="n")
#lines(meancis$m10_LMA.Narea~c(1,2,3,4,5))
#lines(meancis$m90_LMA.Narea~c(1,2,3,4,5))
polygon(x = c(.8,2,3.2,3.2,2,.8), y=na.omit(c(meancis$m10_LMA.Narea, rev(meancis$m90_LMA.Narea))), border="#DDDDDD",col="#DDDDDD")#lightgrey",col = "lightgrey")
m <- cor.test(x=allspp$log.LMA, y=allspp$log.Narea)
points(y=m$estimate,x=5, pch=15, cex=1.3, col="darkgrey")
arrows(x0=5,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3, col="darkgrey")
m <- cor.test(x=fam.dataclean$log.LMA, y=fam.dataclean$log.Narea)
points(y=m$estimate,x=4, pch=17, cex=1.3, col="black")
arrows(x0=4,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3)
# layering points on top of null model
points(Rho_LMA.Narea~jitter(as.numeric(Type)), all.results.cl[which(all.results.cl$n_LMA.Narea>5 & all.results.cl$sig_LMA.Narea<crit),], pch=16)
points(Rho_LMA.Narea~jitter(as.numeric(Type)), all.results.cl[which(all.results.cl$n_LMA.Narea>5 & all.results.cl$sig_LMA.Narea>=crit),], pch=1)
axis(1,labels = c("w/in Spp","w/in Gen", "w/in Fam","btw Fam","global"), at=c(1,2,3,4,5), las=3)
abline(h=0, lty=2)
#text(y=par()$usr[4]-.2,x=.5, labels = "b)")
text(y=par()$usr[4]-.2,x=3, labels = "c) LMA vs Narea")




# LL v Narea
p <- boxplot(Rho_LL.Narea~Type, all.results.cl[which(all.results.cl$n_LL.Narea>5 & is.na(all.results.cl$Taxo.Unit)),]
             , ylim=c(-1,1.3),las=3, ylab="Rho"
             , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
             ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
             , staplelwd=0, outpch=16, outcex=.5, outcol=mypal[colchoices]
             , boxwex=.7, xaxt="n")
#lines(meancis$m10_LL.Narea~c(1,2,3,4,5))
#lines(meancis$m90_LL.Narea~c(1,2,3,4,5))
polygon(x = c(.8,2,3.2,3.2,2,.8), y=na.omit(c(meancis$m10_LL.Narea, rev(meancis$m90_LL.Narea))), border="#DDDDDD",col="#DDDDDD")#lightgrey",col = "lightgrey")
m <- cor.test(x=allspp$log.LL, y=allspp$log.Narea)
points(y=m$estimate,x=5, pch=15, cex=1.3, col="darkgrey")
arrows(x0=5,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3, col="darkgrey")
m <- cor.test(x=fam.dataclean$log.LL, y=fam.dataclean$log.Narea)
points(y=m$estimate,x=4, pch=17, cex=1.3, col="black")
arrows(x0=4,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3)
# layering points on top of null model
points(Rho_LL.Narea~jitter(as.numeric(Type)), all.results.cl[which(all.results.cl$n_LL.Narea>5 & all.results.cl$sig_LL.Narea<crit),], pch=16)
points(Rho_LL.Narea~jitter(as.numeric(Type)), all.results.cl[which(all.results.cl$n_LL.Narea>5 & all.results.cl$sig_LL.Narea>=crit),], pch=1)
axis(1,labels = c("w/in Spp","w/in Gen", "w/in Fam","btw Fam","global"), at=c(1,2,3,4,5), las=3)
abline(h=0, lty=2)
#text(y=par()$usr[4]-.2,x=.5, labels = "b)")
text(y=par()$usr[4]-.2,x=3, labels = "d) Leaf Life vs Narea")




# LMA v LL
p <- boxplot(Rho_LMA.LL~Type, all.results.cl[which(all.results.cl$n_LMA.LL>5 & is.na(all.results.cl$Taxo.Unit)),]
             , ylim=c(-1,1.3),las=3, ylab="Rho"
             , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
             ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
             , staplelwd=0, outpch=16, outcex=.5, outcol=mypal[colchoices]
             , boxwex=.7, xaxt="n")
#lines(meancis$m10_LMA.LL~c(1,2,3,4,5))
#lines(meancis$m90_LMA.LL~c(1,2,3,4,5))
polygon(x = c(.8,2,3.2,3.2,2,.8), y=na.omit(c(meancis$m10_LMA.LL, rev(meancis$m90_LMA.LL))), border="#DDDDDD",col="#DDDDDD")#lightgrey",col = "lightgrey")
m <- cor.test(x=allspp$log.LMA, y=allspp$log.LL)
points(y=m$estimate,x=5, pch=15, cex=1.3, col="darkgrey")
arrows(x0=5,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3, col="darkgrey")
m <- cor.test(x=fam.dataclean$log.LMA, y=fam.dataclean$log.LL)
points(y=m$estimate,x=4, pch=17, cex=1.3, col="black")
arrows(x0=4,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3)
# layering points on top of null model
points(Rho_LMA.LL~jitter(as.numeric(Type)), all.results.cl[which(all.results.cl$n_LMA.LL>5 & all.results.cl$sig_LMA.LL<crit),], pch=16)
points(Rho_LMA.LL~jitter(as.numeric(Type)), all.results.cl[which(all.results.cl$n_LMA.LL>5 & all.results.cl$sig_LMA.LL>=crit),], pch=1)
axis(1,labels = c("w/in Spp","w/in Gen", "w/in Fam","btw Fam","global"), at=c(1,2,3,4,5), las=3)
abline(h=0, lty=2)
#text(y=par()$usr[4]-.2,x=.5, labels = "b)")
text(y=par()$usr[4]-.2,x=3, labels = "e) Leaf Life vs LMA")

plot(1~1, type="n", bty="n", xlab="", ylab="", xaxt="n", yaxt="n")
legend("center",legend = c("not different from null","different from null"), pch = c(1,16),bty="n")
legend("bottomleft", bty="n", legend = "mean of 10%-90%\nnull distributions", fill = "#DDDDDD")




#### sample sizes
xtabs(~sig_LMA.N+Type, all.results.cl[which(!is.na(all.results.cl$sig_LMA.N) & all.results.cl$n_LMA.N>5),])

xtabs(~sig_LL.N+Type, all.results.cl[which(!is.na(all.results.cl$sig_LL.N) & all.results.cl$n_LL.N>5),])
xtabs(~sig_LMA.LL+Type, all.results.cl[which(!is.na(all.results.cl$sig_LMA.LL) & all.results.cl$n_LMA.LL>5),])
xtabs(~sig_LMA.Narea+Type, all.results.cl[which(!is.na(all.results.cl$sig_LMA.Narea) & all.results.cl$n_LMA.Narea>5),])

xtabs(~sig_LL.Narea+Type, all.results.cl[which(!is.na(all.results.cl$sig_LL.Narea) & all.results.cl$n_LL.Narea>5),])






######## Figure S6: Leaf Fraction Figure #################


quartz(width=5, height=3)
par(mar=c(4,4,2,1), mfrow=c(1,2), mgp=c(2.5,1,0))

plot(RGR~LeafFrac, plotavs90[-which(plotavs90$PLOT_ID %in% c(43,55)),], pch=16, col=SP.ID, log="xy", main="a) W/in Species", xlab="Leaf Fraction", cex.main=.9)
legend("topleft", legend=levels(plotavs90$SP.ID), col=mypal[1:5], pch=16, bty="n", cex=.7)
#mtext(side=3,adj=-.1, "a)")as variation in leaf resource use strategy

plot(RGR~LeafFrac, biomass[-which(biomass$PLOT_ID %in% c(43,55)),], pch=16, log="xy", xlab="Leaf Fraction", main="b) Across Communities", cex.main=.9)
#mtext(side=3,adj=-.1, "b)")