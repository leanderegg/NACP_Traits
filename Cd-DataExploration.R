

quartz(width=4.4, height=3.8)
# focusing on the Pines
ggplot(allspp, aes(x=MAT,y=log.LL)) + geom_point(col="grey") +
  geom_point(data=traits[which(traits$GENUS=="Pinus"),], aes(col=SP.ID)) + 
  geom_smooth(data=traits[which(traits$GENUS=="Pinus"),], aes(col=SP.ID), method="lm", se=F) +
  geom_point(data=allspp[which(allspp$Genus=="Pinus"),]) + geom_smooth(data=allspp[which(allspp$Genus=="Pinus"),], method="lm",se=F, col="black") +
  geom_point(data=allgen[which(allgen$Family=="Pinaceae" ),], shape=2) + geom_smooth(data=allgen[which(allgen$Family=="Pinaceae")[-2],], method="lm", se=F, col="black", lwd=1.4)
  
#focusing on the Pinaceae
p1 <- ggplot(allspp, aes(x=MAT,y=log.LL)) + geom_point(col="grey") +
  geom_point(data=traits[which(traits$Family=="Pinaceae"),], aes(col=SP.ID), alpha=0.06) + 
  geom_smooth(data=traits[which(traits$Family=="Pinaceae"),], aes(col=SP.ID), method="lm", se=F) + theme(legend.position = "none") + ggtitle("Within Species")

p2 <- ggplot(allspp, aes(x=MAT,y=log.LL)) + geom_point(col="grey") +
  geom_point(data=allspp[which(allspp$Family=="Pinaceae"),], aes(col=Genus)) + geom_smooth(data=allspp[which(allspp$Family=="Pinaceae"),],aes(col=Genus), method="lm",se=F) + theme(legend.position = "none") + ggtitle("Within Genera")

p3 <- ggplot(allspp, aes(x=MAT,y=log.LL)) + geom_point(col="grey") +
  geom_point(data=allgen[which(allgen$Family=="Pinaceae" ),], aes(col=Genus)) + geom_smooth(data=allgen[which(allgen$Family=="Pinaceae")[-2],], method="lm", se=F, lwd=1)+ theme(legend.position = "none") + ggtitle("Pinaceae")

quartz(width=7, height=3)
multiplot(p1,p2,p3, cols=3)



(p1 <- ggplot(allspp, aes(x=MAT,y=log.LL)) + geom_point(col="grey") +
  geom_point(data=allspp[which(allspp$Genus %in% c("Quercus", "Salix","Betula","Populus")),], aes(col=Genus)) + geom_smooth(data=allspp[which(allspp$Genus %in% c("Quercus", "Salix","Betula","Populus")),], method="lm",se=F, aes(col=Genus))+
    ggtitle("Within Genera")
) # "Acacia","Piper","Eucalyptus",
(p2 <- ggplot(allspp, aes(x=MAT,y=log.LL)) + geom_point(col="grey") +
    geom_point(data=allgen[which(allgen$Family %in% c("Fabaceae", "Asteraceae", "Ericaceae","Rosaceae","Myrtaceae")),], aes(col=Family)) + 
    geom_smooth(data=allgen[which(allgen$Family %in% c("Fabaceae", "Asteraceae", "Ericaceae","Rosaceae","Myrtaceae")),], method="lm", se=F,aes(col=Family)) +
     ggtitle("Within Angiosperm Families")
) # "Acacia","Piper","Eucalyptus",

quartz(width=7, height=3)
multiplot(p1,p2, cols=2)

  geom_point(data=allgen[which(allgen$Family=="Fabaceae" ),], shape=2) + geom_smooth(data=allgen[which(allgen$Family=="Fabaceae"),], method="lm", se=F, col="black", lwd=2) +
  geom_point(data=allgen[which(allgen$Family=="Asteraceae" ),], shape=3) + geom_smooth(data=allgen[which(allgen$Family=="Asteraceae"),], method="lm", se=F, col="green", lty=2)
