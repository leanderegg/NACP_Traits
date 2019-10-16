####### Playing around with traits for CA project with Phil ###########

require(dplyr)


biomass <- read.csv("./FinalExample/DerivedData/PNW_Biomass_data_for_Anderegg_etal_2018_EcolLet.csv", header=T, row.names = 1)
traits <- read.csv("./FinalExample/DerivedData/PNW_Trait_data_for_Anderegg_etal_2018_EcolLet.csv", header=T, row.names=1)


colnames(traits)

ggplot(traits[which(traits$FullSpecies %in% c("Abies amabilis","Abies concolor","Abies grandis","Arbutus menziesii","Calocedrus decurrens","Juniperus occidentalis","Pinus contorta","Pinus Jeffreyi","Pinus ponderosa","Pseudotsuga menziesii","Tsuga heterophylla")),], aes(x=cmi.gy.mm, y=LMA)) + geom_point()+geom_smooth(method="lm") + facet_wrap(facets = ~FullSpecies) + theme(legend.position = "none")

traits %>% filter(FullSpecies %in% c("Abies amabilis","Abies concolor","Abies grandis","Arbutus menziesii","Calocedrus decurrens","Juniperus occidentalis","Pinus contorta","Pinus Jeffreyi","Pinus ponderosa","Pseudotsuga menziesii","Tsuga heterophylla", "Quercus chrysolepis","Quercus garryana","Quercus kelloggii")) %>%
  group_by (FullSpecies) %>% summarize(q05 = quantile(LMA,probs = c(0.05), na.rm=T),
            q25 = quantile(LMA,probs =.25, na.rm=T),
            q50 = quantile(LMA,probs = .5, na.rm=T),
            q75 = quantile(LMA,probs = .75, na.rm=T),
            q95 = quantile(LMA,probs = .95, na.rm=T))
