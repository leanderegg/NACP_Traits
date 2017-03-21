################################################################
###      FIGURES for of NACP_TERRA trait and stand data
###         to assess within-species trait variation
###
###         data downloaded from: http://dx.doi.org/10.3334/ORNLDAAC/1292 
###            on 01/23/16 by LDLA
###############################################################


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


##### Uppper Panel: boxplots ............................................................
quartz(width=5.5, height=4.5,title = "LMA v Nmass")
par(mfrow=c(2,2), mar=c(0,4,0,1), oma=c(6,0,2,0), cex=1)

### Slope boxplots (a & b)
  ## LMA v Nmass
p <- boxplot(Slope_LMA.N~Type, all.results[which(all.results$n_LMA.N>5),]
             , ylim=c(-2.2,1.5),las=3, ylab="MA Slope" #, main="log(LMA)~log(Nmass) strict"
             , col=paste0(mypal[1:6],"66"), boxcol=paste0(mypal[1:6],"66")
             ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[1:6],"AA")
             , staplelwd=0, outpch=16, outcex=.5, outcol=mypal[1:6]
             , boxwex=.7, xaxt="n")
abline(h=0, lty=2)
#text(y=par()$usr[3]+.2, x=c(1,2,3,4,5,6), labels = p$n)
#text(y=par()$usr[4]-.2,x=.5, labels = "a)")
mtext(text = "a)", side = 3, adj=0, line=.2)
mtext(text= "LMA vs Nmass", side=3, line=.2)
  # LL v Nmass
p <- boxplot(Slope_N.LL~Type, all.results[which(all.results$n_N.LL>5),]
             , ylim=c(-3.5,3),las=3, ylab="MA Slope" #, main="log(LMA)~log(Nmass) strict"
             , col=paste0(mypal[1:6],"66"), boxcol=paste0(mypal[1:6],"66")
             ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[1:6],"AA")
             , staplelwd=0, outpch=16, outcex=.5, outcol=mypal[1:6]
             , boxwex=.7, xaxt="n")
abline(h=0, lty=2)
#text(y=par()$usr[3]+.2, x=c(1,2,3,4,5,6), labels = p$n)
#text(y=par()$usr[4]-.2,x=.5, labels = "b)")
mtext(text = "b)", side = 3, adj=0, line=.2)
mtext(text= "LL vs Nmass", side=3, line=.2)

### Lower boxplots with Rho
  # LMA v Nmass
p <- boxplot(Rho_LMA.N~Type, all.results[which(all.results$n_LMA.N>5),]
             , ylim=c(-1,1.1),las=3, ylab="Rho"
             , col=paste0(mypal[1:6],"66"), boxcol=paste0(mypal[1:6],"66")
             ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[1:6],"AA")
             , staplelwd=0, outpch=16, outcex=.5, outcol=mypal[1:6]
             , boxwex=.7)
abline(h=0, lty=2)
text(y=par()$usr[4]-.2, x=c(1,2,3,4,5,6), labels = p$n)
#text(y=par()$usr[4]-.2,x=.5, labels = "b)")
  # LL an Nmass
p <- boxplot(Rho_N.LL~Type, all.results[which(all.results$n_N.LL>5),]
             , ylim=c(-1,1.1),las=3, ylab="Rho"
             , col=paste0(mypal[1:6],"66"), boxcol=paste0(mypal[1:6],"66")
             ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[1:6],"AA")
             , staplelwd=0, outpch=16, outcex=.5, outcol=mypal[1:6]
             , boxwex=.7)
abline(h=0, lty=2)
text(y=par()$usr[4]-.2, x=c(1,2,3,4,5,6), labels = p$n)
#text(y=par()$usr[4]-.2,x=.5, labels = "b)")


##### Lower Panel: funnel plots ............................................................

quartz(width=5.5,height=3)
par(mar=c(4,4,0,1), mfrow=c(1,2), mgp=c(2.5,1,0), oma=c(0,0,2,0))
plot(Rho_LMA.N~varNmass, all.results, pch=16, col=Type, xlab="Var. in %N", ylab="Rho   (LMA v Nmass)")
#mtext( text="Nmass v LMA", side=3, line=0, font=2)
abline(h=0, col="grey", lty=2)
abline(h=all.results$Rho_LMA.N[which(all.results$Taxo.Unit=="fam.all")])
mtext(text = "c)", side = 3, adj=0, line=.2)
points(Rho_LMA.N~varNmass, all.results[c(nrow(all.results)-1,nrow(all.results)),], pch=c(24,25), cex=1.4, bg=mypal[5:6])

plot(Rho_N.LL~varNmass, all.results, pch=16, col=Type, xlab="Var. in %N", ylab="Rho   (LL v Nmass)")
abline(h=0, col="grey", lty=2)
abline(h=all.results$Rho_N.LL[which(all.results$Taxo.Unit=="fam.all")])
legend('topright', legend = levels(all.results.LMAN$Type), pch=c(16,16,16,16,24,25), col=c(mypal[1:4], "black","black"), pt.bg= mypal[1:6], bty ="n", cex=.7)
mtext(text = "d)", side = 3, adj=0, line=.2)
points(Rho_N.LL~varNmass, all.results[c(nrow(all.results)-1,nrow(all.results)),], pch=c(24,25), cex=1.4, bg=mypal[5:6])




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

quartz(width=7.008, height=4)
layout(mat)
par(mar=c(0,4,6,1))
p <- boxplot(Slope_LMA.LL~Type, all.results[which(all.results$n_LMA.LL>5),]
             , ylim=c(-2.5,4),las=3, ylab="MA Slope" #, main="log(LMA)~log(Nmass) strict"
             , col=paste0(mypal[1:6],"66"), boxcol=paste0(mypal[1:6],"66")
             ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[1:6],"AA")
             , staplelwd=0, outpch=16, outcex=.5, outcol=mypal[1:6]
             , boxwex=.7, xaxt="n")
abline(h=0, lty=2)
#text(y=par()$usr[3]+.2, x=c(1,2,3,4,5,6), labels = p$n)
#text(y=par()$usr[4]-.2,x=.5, labels = "a)")
mtext(text = "a)", side = 3, adj=0, line=.2)
mtext(text= "LMA vs LL", side=3, line=.2)
par(mar=c(6,4,0,1))
p <- boxplot(Rho_LMA.LL~Type, all.results[which(all.results$n_LMA.LL>5),]
             , ylim=c(-1,1.4),las=3, ylab="Rho"
             , col=paste0(mypal[1:6],"66"), boxcol=paste0(mypal[1:6],"66")
             ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[1:6],"AA")
             , staplelwd=0, outpch=16, outcex=.5, outcol=mypal[1:6]
             , boxwex=.7)
abline(h=0, lty=2)
text(y=par()$usr[4]-.2, x=c(1,2,3,4,5,6), labels = p$n)

par(mar=c(4,4,1.5,1.5))
## scatterplot
####### Plotting the LMA vs Narea scaling in unit rather than log space ######
plot(log.LL~log.LMA, LES, col="grey", pch=16, ylab="log(LL)", xlab="log(LMA)")

tax <- "w.inSpp"
for (i in as.character(all.results$Taxo.Unit[which(all.results$Type==tax & all.results$n_LMA.LL>5)])){
  plot.MAR(xvar = "log.LMA", yvar = "log.LL",data= spp.data[which(spp.data$Species==i),], linecol = mypal[1])
}

tax <- "w.inGen"
for (i in as.character(all.results$Taxo.Unit[which(all.results$Type==tax & all.results$n_LMA.LL>5)])){
  plot.MAR(xvar = "log.LMA", yvar = "log.LL",data= gen.data[which(gen.data$Genus==i),], linecol = mypal[2])
}

tax <- "Genw.inFam"
for (i in as.character(all.results$Taxo.Unit[which(all.results$Type==tax & all.results$n_LMA.LL>5)])){
  plot.MAR(xvar = "log.LMA", yvar = "log.LL",data= geninfam.data[which(geninfam.data$Family==i),], linecol = mypal[4])
}


plot.MAR(xvar="log.LMA", yvar="log.LL", data= fam.data, linecol = mypal[5], lwd=2)

plot.MAR(xvar="log.LMA", yvar="log.LL", data= LES, linecol = "black", lwd=2)
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

quartz(width=7.008, height=4)
layout(mat)
par(mar=c(0,4,6,1))
p <- boxplot(Slope_LMA.Narea~Type, all.results[which(all.results$n_LMA.LL>5),]
             , ylim=c(-0.5,2.5),las=3, ylab="MA Slope" #, main="log(LMA)~log(Nmass) strict"
             , col=paste0(mypal[1:6],"66"), boxcol=paste0(mypal[1:6],"66")
             ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[1:6],"AA")
             , staplelwd=0, outpch=16, outcex=.5, outcol=mypal[1:6]
             , boxwex=.7, xaxt="n")
abline(h=0, lty=2)
#text(y=par()$usr[3]+.2, x=c(1,2,3,4,5,6), labels = p$n)
#text(y=par()$usr[4]-.2,x=.5, labels = "a)")
mtext(text = "a)", side = 3, adj=0, line=.2)
mtext(text= "LMA vs Narea", side=3, line=.2)
par(mar=c(6,4,0,1))
p <- boxplot(Rho_LMA.Narea~Type, all.results[which(all.results$n_LMA.Narea>5),]
             , ylim=c(-1,1.4),las=3, ylab="Rho"
             , col=paste0(mypal[1:6],"66"), boxcol=paste0(mypal[1:6],"66")
             ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[1:6],"AA")
             , staplelwd=0, outpch=16, outcex=.5, outcol=mypal[1:6]
             , boxwex=.7)
abline(h=0, lty=2)
text(y=par()$usr[4]-.2, x=c(1,2,3,4,5,6), labels = p$n)

par(mar=c(4,4,1.5,1.5))
## scatterplot
####### Plotting the LMA vs Narea scaling in unit rather than log space ######
plot(log.Narea~log.LMA, LES, col="grey", pch=16, ylab="log(Narea)", xlab="log(LMA)")

tax <- "w.inGen"
for (i in as.character(all.results$Taxo.Unit[which(all.results$Type==tax & all.results$n_LMA.Narea>5)])){
  plot.MAR(xvar = "log.LMA", yvar = "log.Narea",data= gen.data[which(gen.data$Genus==i),], linecol = mypal[2])
}

tax <- "Genw.inFam"
for (i in as.character(all.results$Taxo.Unit[which(all.results$Type==tax & all.results$n_LMA.Narea>5)])){
  plot.MAR(xvar = "log.LMA", yvar = "log.Narea",data= geninfam.data[which(geninfam.data$Family==i),], linecol = mypal[4])
}


plot.MAR(xvar="log.LMA", yvar="log.Narea", data= fam.data, linecol = mypal[5], lwd=2)

plot.MAR(xvar="log.LMA", yvar="log.Narea", data= LES, linecol = "black", lwd=2)
tax <- "w.inSpp"
for (i in as.character(all.results$Taxo.Unit[which(all.results$Type==tax & all.results$n_LMA.Narea>5)])){
  plot.MAR(xvar = "log.LMA", yvar = "log.Narea",data= spp.data[which(spp.data$Species==i),], linecol = mypal[1])
}
mtext(text = "b)", side = 3, adj=0, line=.2)







#________________________________________________________________
###### FIG S1: Funnel Plots for LMA vs LL ######
#________________________________________________________________


quartz(width=5.5,height=4.5)
par(mar=c(4,4,0,1), mfrow=c(2,2), mgp=c(2.5,1,0), oma=c(0,0,2,0))
plot(Slope_LMA.LL~varLMA, all.results, pch=16, col=Type, xlab="Var. in LMA", ylab="Slope   (LMA v LL)")
#mtext( text="Nmass v LMA", side=3, line=0, font=2)
abline(h=0, col="grey", lty=2)
abline(h=all.results$Slope_LMA.LL[which(all.results$Taxo.Unit=="fam.all")])
mtext(text = "a)", side = 3, adj=0, line=.2)
points(Slope_LMA.LL~varLMA, all.results[c(nrow(all.results)-1,nrow(all.results)),], pch=c(24,25), cex=1.4, bg=mypal[5:6])
abline(h=mean(all.results$Slope_LMA.LL[which(all.results$Type=="w.inSpp")], na.rm=T), col=mypal[1])


plot(Rho_LMA.LL~varLMA, all.results, pch=16, col=Type, xlab="Var. in LMA", ylab="Rho   (LMA v LL)")
#mtext( text="Nmass v LMA", side=3, line=0, font=2)
abline(h=0, col="grey", lty=2)
abline(h=all.results$Rho_LMA.LL[which(all.results$Taxo.Unit=="fam.all")])
mtext(text = "b)", side = 3, adj=0, line=.2)
points(Rho_LMA.LL~varLMA, all.results[c(nrow(all.results)-1,nrow(all.results)),], pch=c(24,25), cex=1.4, bg=mypal[5:6])
abline(h=mean(all.results$Rho_LMA.LL[which(all.results$Type=="w.inSpp")], na.rm=T), col=mypal[1])


plot(Slope_LMA.LL~n_LMA.LL, all.results, pch=16, col=Type, xlab="Var. in LL", ylab="Slope   (LMA v LL)")
abline(h=0, col="grey", lty=2)
abline(h=all.results$Slope_LMA.LL[which(all.results$Taxo.Unit=="fam.all")])
mtext(text = "c)", side = 3, adj=0, line=.2)
points(Slope_LMA.LL~n_LMA.LL, all.results[c(nrow(all.results)-1,nrow(all.results)),], pch=c(24,25), cex=1.4, bg=mypal[5:6])
abline(h=mean(all.results$Slope_LMA.LL[which(all.results$Type=="w.inSpp")], na.rm=T), col=mypal[1])

plot(Rho_LMA.LL~n_LMA.LL, all.results, pch=16, col=Type, xlab="Var. in LL", ylab="Rho   (LMA v LL)")
abline(h=0, col="grey", lty=2)
abline(h=all.results$Rho_LMA.LL[which(all.results$Taxo.Unit=="fam.all")])
legend('bottomright', legend = levels(all.results$Type), pch=c(16,16,16,16,24,25), ncol=2, col=c(mypal[1:4], "black","black"), pt.bg= mypal[1:6], bty ="n", cex=.7)
mtext(text = "d)", side = 3, adj=0, line=.2)
points(Rho_LMA.LL~n_LMA.LL, all.results[c(nrow(all.results)-1,nrow(all.results)),], pch=c(24,25), cex=1.4, bg=mypal[5:6])
abline(h=mean(all.results$Rho_LMA.LL[which(all.results$Type=="w.inSpp")], na.rm=T), col=mypal[1])

