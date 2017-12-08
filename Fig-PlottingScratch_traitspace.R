


quartz(width=4.33, height=4.73)
par(mfrow=c(2,2), mgp=c(2,.7,0), cex.lab=1.1, cex.axis=1.1, mar=c(4.5,2,1.5,2),oma=c(0,2,3.8,0))


plot(log.Narea~log.LMA, allspp, col="grey", pch=16, xlab=expression(paste(log[10](LMA))))
mtext(expression(paste(log[10](N[area]))), side=2, line=2)
points(log.Narea~log.LMA, traits)
#points(log.Narea~log.LMA, spp.data, pch=3, col="darkblue")

points(log.Narea~log.LMA, traits.domcon1, pch=16, cex=.8, col=mypal[3])
points(biomass$log.cw_Nareap_if~biomass$log.cw_LMAp_if, bg=mypal[5], pch=24)

mtext("a)", side=3, adj=-.1, line=.3)
legend(x=1.1, y=2.5, xpd=NA,legend = c("Global","PNW woody plants", "Evgrn Needle PFT", "PNW CWMs"), pch=c(16,1,16,24), col=c("grey","black",mypal[3],"black"), pt.bg = mypal[5], bty="n")



quartz(width=4.33, height=4.73)
par(mfrow=c(2,2), mgp=c(2,.7,0), cex.lab=1.1, cex.axis=1.1, mar=c(4.5,2,1.5,2),oma=c(0,2,3.8,0))

plot(I(10^log.Narea)~I(10^log.LMA), allspp, col="grey", pch=16, xlab="raw LMA")
mtext("raw Narea", side=2, line=2)
#points(I(10^log.Narea)~I(10^log.LMA), spp.data, pch=3, col="darkblue")
points(I(10^log.Narea)~I(10^log.LMA), traits)
points(I(10^log.Narea)~I(10^log.LMA), traits.domcon1, pch=16, cex=.8, col=mypal[3])
points(I(10^biomass$log.cw_Nareap_if)~I(10^biomass$log.cw_LMAp_if), bg=mypal[5], pch=24)
mtext("a)", side=3, adj=-.1, line=.3)
legend(x=1.1, y=2.5, xpd=NA,legend = c("Global","PNW woody plants", "Evgrn Needle PFT", "PNW CWMs"), pch=c(16,1,16,24), col=c("grey","black",mypal[3],"black"), pt.bg = mypal[5], bty="n")


plot(log.Nmass~log.LL, allspp, col="grey", pch=16, xlab=expression(paste(log[10](LL))))
mtext(expression(paste(log[10](N[mass]))), side=2, line=2)
points(log.Nmass~log.LL, spp.data, pch=3, col="darkblue")
points(log.Nmass~log.LL, traits)
points(log.Nmass~log.LL, traits.domcon1, pch=16, cex=.8, col=mypal[3])
points(biomass$log.cw_Nmassp_if~biomass$log.cw_LLp_if, bg=mypal[5], pch=24)
mtext("a)", side=3, adj=-.1, line=.3)
legend(x=1.1, y=2.5, xpd=NA,legend = c("Global","PNW woody plants", "Evgrn Needle PFT", "PNW CWMs"), pch=c(16,1,16,24), col=c("grey","black",mypal[3],"black"), pt.bg = mypal[5], bty="n")


plot(I(10^log.Nmass)~I(10^log.LL), allspp, col="grey", pch=16, xlab="Leaf Life (months)")
mtext("raw Nmass", side=2, line=2)
points(I(10^log.Nmass)~I(10^log.LL), spp.data, pch=3, col="darkblue")

points(I(10^log.Nmass)~I(10^log.LL), traits)
points(I(10^log.Nmass)~I(10^log.LL), traits.domcon1, pch=16, cex=.8, col=mypal[3])
points(I(10^biomass$log.cw_Nmassp_if)~I(10^biomass$log.cw_LLp_if), bg=mypal[5], pch=24)
mtext("a)", side=3, adj=-.1, line=.3)
legend(x=1.1, y=2.5, xpd=NA,legend = c("Global","PNW woody plants", "Evgrn Needle PFT", "PNW CWMs"), pch=c(16,1,16,24), col=c("grey","black",mypal[3],"black"), pt.bg = mypal[5], bty="n")





plot(log.Nmass~log.LMA, allspp, col="grey", pch=16, xlab=expression(paste(log[10](LMA))))
mtext(expression(paste(log[10](N[area]))), side=2, line=2)
points(log.Nmass~log.LMA, traits)
points(log.Nmass~log.LMA, traits.domcon1, pch=16, cex=.8, col=mypal[3])
points(biomass$log.cw_Nmassp_if~biomass$log.cw_LMAp_if, bg=mypal[5], pch=24)
mtext("a)", side=3, adj=-.1, line=.3)
legend(x=1.1, y=2.5, xpd=NA,legend = c("Global","PNW woody plants", "Evgrn Needle PFT", "PNW CWMs"), pch=c(16,1,16,24), col=c("grey","black",mypal[3],"black"), pt.bg = mypal[5], bty="n")


plot(I(10^log.Nmass)~I(10^log.LMA), allspp, col="grey", pch=16, xlab=expression(paste(log[10](LMA))))
mtext(expression(paste(log[10](N[area]))), side=2, line=2)
points(I(10^log.Nmass)~I(10^log.LMA), traits)
points(I(10^log.Nmass)~I(10^log.LMA), traits.domcon1, pch=16, cex=.8, col=mypal[3])
points(I(10^biomass$log.cw_Nmassp_if)~I(10^biomass$log.cw_LMAp_if), bg=mypal[5], pch=24)
mtext("a)", side=3, adj=-.1, line=.3)
legend(x=1.1, y=2.5, xpd=NA,legend = c("Global","PNW woody plants", "Evgrn Needle PFT", "PNW CWMs"), pch=c(16,1,16,24), col=c("grey","black",mypal[3],"black"), pt.bg = mypal[5], bty="n")



plot(log.LMA~log.LL, allspp, col="grey", pch=16, xlab=expression(paste(log[10](LL))))
mtext(expression(paste(log[10](N[area]))), side=2, line=2)
points(log.LMA~log.LL, traits)
points(log.LMA~log.LL, traits.domcon1, pch=16, cex=.8, col=mypal[3])
points(biomass$log.cw_LMAp_if~biomass$log.cw_LLp_if, bg=mypal[5], pch=24)
mtext("a)", side=3, adj=-.1, line=.3)
#legend(x=1.1, y=2.5, xpd=NA,legend = c("Global","PNW woody plants", "Evgrn Needle PFT", "PNW CWMs"), pch=c(16,1,16,24), col=c("grey","black",mypal[3],"black"), pt.bg = mypal[5], bty="n")
points(log.LMA~log.LL, vib, pch=4, col=species)
points(log.LMA~log.LL, arab, pch=5)


for(i in levels(vib$species)){
  plot.MAR(xvar = "log.LMA", yvar = "log.LL",data = vib[which(vib$species==i),], linecol = "black")
}

plot(I(10^log.LMA)~I(10^log.LL), allspp, col="grey", pch=16, xlab=expression(paste(log[10](LL))))
mtext(expression(paste(log[10](N[area]))), side=2, line=2)
points(I(10^log.LMA)~I(10^log.LL), traits)
points(I(10^log.LMA)~I(10^log.LL), traits.domcon1, pch=16, cex=.8, col=mypal[3])
points(I(10^biomass$log.cw_LMAp_if)~I(10^biomass$log.cw_LLp_if), bg=mypal[5], pch=24)
mtext("a)", side=3, adj=-.1, line=.3)
legend(x=1.1, y=2.5, xpd=NA,legend = c("Global","PNW woody plants", "Evgrn Needle PFT", "PNW CWMs"), pch=c(16,1,16,24), col=c("grey","black",mypal[3],"black"), pt.bg = mypal[5], bty="n")



vib.data <- vib
colnames(vib.data)[1] <- "Species"
vib.results <- data.frame(matrix(NA, nrow=length(unique(vib.data$Species)), ncol=17))
colnames(vib.results) <- c("Species", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varLMA","varLL","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(vib.data$Species))){
  species <- levels(vib.data$Species)[i]
  print(species)
  dataz <- vib.data[which(vib.data$Species==species),]
  res <- fit.MAR(xvar='log.LMA',yvar="log.LL",data=dataz)
  vib.results[i,1] <- species
  vib.results[i,2:10] <- res
  if (!is.na(res[1]) &res[7]>4){ # only fit null model if there are >5 data points
    nullbounds <- fit.null(xvar='log.LMA', yvar="log.LL", observed = dataz, nulldata = allspp, nits = 1000)  
    vib.results[i, 11:16] <- nullbounds
    vib.results[i,17] <- test.sig(x=vib.results$Rho[i], test=nullbounds)
  }
}
vib.results$rho.sig <- rho.sig(vib.results$Rho, vib.results$n)


hel.data <- hel
hel.data$log.LL <- log(hel.data$LL_days/30, base=10)
hel.data$log.LMA <- log(hel.data$LMA, base=10)
hel.results <- data.frame(matrix(NA, nrow=length(unique(hel.data$Species)), ncol=17))
colnames(hel.results) <- c("Species", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varLMA","varLL","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(hel.data$Species))){
  species <- levels(hel.data$Species)[i]
  print(species)
  dataz <- hel.data[which(hel.data$Species==species),]
  res <- fit.MAR(xvar='log.LMA',yvar="log.LL",data=dataz)
  hel.results[i,1] <- species
  hel.results[i,2:10] <- res
  if (!is.na(res[1]) &res[7]>4){ # only fit null model if there are >5 data points
    nullbounds <- fit.null(xvar='log.LMA', yvar="log.LL", observed = dataz, nulldata = allspp, nits = 1000)  
    hel.results[i, 11:16] <- nullbounds
    hel.results[i,17] <- test.sig(x=hel.results$Rho[i], test=nullbounds)
  }
}
hel.results$rho.sig <- rho.sig(hel.results$Rho, hel.results$n)



plot(log.LL~log.LMA, allspp, col="grey", pch=16, xlab=expression(paste(log[10](LL))))
points(log.LL~log.LMA, hel.data, col=Species, pch=as.numeric(Species))



for(i in levels(hel.data$Species)){
  lty <- 2
  lwd <- 1
  tmp <- hel.data[which(hel.data$Species==i),]
  if(hel.results$rho.sig[which(hel.results$Species==i)]<0.2){ lty<-1; lwd<-2}
  plot.MAR(xvar = "log.LMA", yvar = "log.LL",data = tmp, linecol = "black", lty=lty, lwd=lwd)
}
