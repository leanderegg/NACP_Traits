######## attempt at defining a null model for correlations #####################


colnames(allspp)

colnames(spp.data)

levels(spp.data$Species)

tax <- "w.inSpp"
for (i in as.character(all.results.cl$Taxo.Unit[which(all.results.cl$Type==tax & all.results.cl$n_LMA.N>5)])){
  plot.MAR(xvar = "log.LMA", yvar = "log.Nmass",data= spp.data[which(spp.data$Species==i),], linecol = mypal[colchoices[1]])
}
xvar <- "log.Nmass"
yvar <- "log.LMA"
observed <- spp.data[which(spp.data$Species=="Abies amabilis"),]
nulldata <- allspp
nits <- 100


test.sig <- function(x){
  if(x[3]>x[2] | x[3]<x[1]) return(1)
  else return(0)
}

fit.null <- function(xvar, yvar, observed, nulldata, nits){
  # find the trait ranges to sample
  rangeX <- range(observed[,xvar], na.rm=T)
  difX <- rangeX[2]-rangeX[1]
  #meanX <- mean(observed[,xvar], na.rm=T)
  rangeY <- range(observed[,yvar], na.rm=T)
  difY <- rangeY[2]-rangeY[1]
  #meanY <- mean(observed[,yvar], na.rm=T)
  # cut down the null data to just have non NAs
  nulldata <- nulldata[which(!is.na(nulldata[,xvar]) & !is.na(nulldata[,yvar])),]
  #nulldata2 <- nulldata %>% filter(!is.na(xvar) & !is.na(yvar))
  restrictednull <- nulldata[which(nulldata[,xvar]> min(nulldata[,xvar], na.rm=T)+difX/2 & nulldata[,xvar]< max(nulldata[,xvar], na.rm=T)-difX/2 & nulldata[,yvar]> min(nulldata[,yvar], na.rm=T)+difY/2 & nulldata[,yvar]< max(nulldata[,yvar], na.rm=T)-difY/2),] 
  nullcor <- c(rep(NA, times=nits))
  for(i in 1:nits){
    center <- restrictednull[sample(nrow(restrictednull),size = 1),]
    #null <- nulldata %>% filter(xvar > center[,xvar]-rangeX/2) #, xvar < center[,xvar]+rangeX/2, yvar > center[,yvar]-rangeY/2, yvar < center[,yvar]-rangeY/2)
    null <- nulldata[which(nulldata[,xvar] > as.numeric(center[,xvar])-difX/2 & nulldata[,xvar] < as.numeric(center[,xvar]+difX/2) & 
                             nulldata[,yvar] > as.numeric(center[,yvar])-difY/2 & nulldata[,yvar] < as.numeric(center[,yvar]+difY/2)),]
    ndist <- null[sample(nrow(null), size = nrow(observed), replace = TRUE),]
    nullcor[i] <- cor(ndist[,c(xvar,yvar)])[2,1]
  }
  return(nullcor)
}

tax <- "w.inSpp"
taxlist <- as.character(all.results.cl$Taxo.Unit[which(all.results.cl$Type==tax & all.results.cl$n_LMA.N>5)])
LMA.Nsppnulls <- data.frame(Species=taxlist, lci = rep(NA, times=length(taxlist)), uci = rep(NA, times=length(taxlist)), Rho.obs=all.results.cl$Rho_LMA.N[which(all.results.cl$Type==tax & all.results.cl$n_LMA.N>5)])
for(i in 1:length(taxlist)){
  species <- taxlist[i]
  nullcorrelations <- fit.null(xvar="log.Nmass", yvar="log.LMA", observed=spp.data[which(spp.data$Species==species),], nulldata= allspp, nits=100)
  LMA.Nsppnulls[i,2:3] <- quantile(nullcorrelations,probs = c(0.05, 0.95), na.rm = T)
}

LMA.Nsppnulls$sig <- apply(X=LMA.Nsppnulls[,c("lci","uci","Rho.obs")],MARGIN = 1,FUN = test.sig)





tax <- "w.inSpp"
taxlist <- as.character(all.results.cl$Taxo.Unit[which(all.results.cl$Type==tax & all.results.cl$n_LMA.LL>5)])
LMA.LLsppnulls <- data.frame(Species=taxlist, lci = rep(NA, times=length(taxlist)), uci = rep(NA, times=length(taxlist)), Rho.obs=all.results.cl$Rho_LMA.LL[which(all.results.cl$Type==tax & all.results.cl$n_LMA.LL>5)])
for(i in 1:length(taxlist)){
  species <- taxlist[i]
  nullcorrelations <- fit.null(xvar="log.LMA", yvar="log.LL", observed=spp.data[which(spp.data$Species==species),], nulldata= allspp, nits=100)
  LMA.LLsppnulls[i,2:3] <- quantile(nullcorrelations,probs = c(0.05, 0.95), na.rm = T)
}

LMA.LLsppnulls$sig <- apply(X=LMA.LLsppnulls[,c("lci","uci","Rho.obs")],MARGIN = 1,FUN = test.sig)

LMA.LLsppnulls$n <- all.results.cl$n_LMA.LL[which(all.results.cl$Type==tax & all.results.cl$n_LMA.LL>5)]
LMA.LLsppnulls_10$n <- all.results.cl$n_LMA.LL[which(all.results.cl$Type==tax & all.results.cl$n_LMA.LL>5)]





tax <- "w.inSpp"
taxlist <- as.character(all.results.cl$Taxo.Unit[which(all.results.cl$Type==tax & all.results.cl$n_LMA.Narea>5)])
LMA.Nareasppnulls <- data.frame(Species=taxlist, lci = rep(NA, times=length(taxlist)), uci = rep(NA, times=length(taxlist)), Rho.obs=all.results.cl$Rho_LMA.Narea[which(all.results.cl$Type==tax & all.results.cl$n_LMA.Narea>5)])
for(i in 1:length(taxlist)){
  species <- taxlist[i]
  nullcorrelations <- fit.null(xvar="log.LMA", yvar="log.Narea", observed=spp.data[which(spp.data$Species==species),], nulldata= allspp, nits=100)
  LMA.Nareasppnulls[i,2:3] <- quantile(nullcorrelations,probs = c(0.05, 0.95), na.rm = T)
}

LMA.Nareasppnulls$sig <- apply(X=LMA.Nareasppnulls[,c("lci","uci","Rho.obs")],MARGIN = 1,FUN = test.sig)
