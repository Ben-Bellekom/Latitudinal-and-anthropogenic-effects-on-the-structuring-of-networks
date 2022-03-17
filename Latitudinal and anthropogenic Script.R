#Selection of packages
library("bipartite")
library("readr")
library("dplyr")
library("ggplot2")
library("skimr")
library("tidyr")
library("stringr")
library("HiveR")
library(igraph)
library(HiveR)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(gridGraphics)
library(ggplotify)
library(plotrix)


# Change file names to actual file name downloaded from GitHub

ecologicaldatafull<- read.csv(file= "Latitudinal and anthropogenic full data.csv", header=T, stringsAsFactors = FALSE)
siteinfo<-read.csv(file= "Site information.csv", header=T, stringsAsFactors = F)


da1<-ecologicaldatafull %>%
  mutate(host =coalesce(host_species,host_order, host_class), vector=coalesce(diptera_species, diptera_genus))

#make webs
allweb<- frame2webs(da1, varnames = c("host", "vector", "webID", "frequency"))

#calculate metrics
metrics<-as.data.frame(lapply(allweb, networklevel, index="interaction evenness"))
metrics1<-as.data.frame(lapply(allweb, networklevel, index="H2"))


m1<-pivot_longer(metrics, cols=starts_with ("P"))
m2<-as.data.frame(m1)

#adding the row names using the web name column and then removing webid column
row.names(m2)<-m2$name
m2[1]<-NULL

# add site info
evenness1<- m2 %>%
  mutate(habitat= siteinfo$Site, 
         lat=siteinfo$latitude,
         websize=siteinfo$websize,
         numhost=siteinfo$numhost,
         numinsect=siteinfo$numinsect,
         site=siteinfo$name,
         matrixsize=siteinfo$matrixsize,
         richness=siteinfo$sprichness
  )

evenness1<- rename(evenness1, Evenness=value)


h2prime<- m2 %>%
  mutate(habitat= siteinfo$Site, 
         lat=siteinfo$latitude,
         websize=siteinfo$websize,
         numhost=siteinfo$numhost,
         numinsect=siteinfo$numinsect,
         site=siteinfo$name,
         matrixsize=siteinfo$matrixsize,
         richness=siteinfo$sprichness
  )
h2prime<- rename(h2prime, H2prime=value)


aggregate(h2prime$value, list(h2prime$habitat), FUN=mean)


# glm for evenness/h2
model1<-glm(Evenness~habitat+lat+richness+log(matrixsize), family = gaussian , data=evenness1)
model2<-glm(Evenness~lat+richness+log(matrixsize), family = gaussian , data=evenness1)

model1<-glm(H2prime~habitat+lat+richness+log(matrixsize), family = gaussian , data=h2prime)
model2<-glm(H2prime~lat+richness+log(matrixsize), family = gaussian , data=h2prime)


#compare models
anova(model2, model1, test = "Chisq")


#posthoc
result<-glm(Evenness~habitat, data=evenness1)
result1<-aov(result)
TukeyHSD(result1, "habitat")



# number of hosts/insects with latitude controlling for matrix size
lme1<-lm(numinsect~lat+matrixsize, data=evenness1)
lme1<-lm(numhost~lat+matrixsize, data=evenness1)



### Creating figures

evennessbox<-ggplot(evenness1, aes(x=habitat, y=Evenness))+
  geom_boxplot(outlier.color = "black")+
  labs(x="Habitat type", y= "Evenness value")+
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


evennesslat<-ggplot(evenness1, aes(x=lat, y=Evenness, shape= habitat))+
  geom_point(size=3)+
  labs(x="Latitude", y= "Evenness value", shape="Habitat type\n")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

primebox<-ggplot(h2prime, aes(x=habitat, y=value))+
  geom_boxplot(outlier.color = "black")+
  theme_bw() +
  labs(x="Habitat type", y= "H2' value")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

primerlat<-ggplot(h2prime, aes(x=lat, y=value, shape= habitat))+
  geom_point(size=3)+
  labs(x="Latitude", y= "H2' value", shape="Habitat type\n")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

combinedfigure<- ggarrange(evennessbox,evennesslat,primebox, primerlat,
                           labels=c("A", "B", "C", "D"),
                           common.legend = T, legend="right",
                           ncol=2, nrow=2)



# Plot networks by habitat type and a full network- data subsets by habitat and removes webID to allow for full network to be plotted


farmlandalldata<-read.csv(file= "farmland data.csv", header=T, stringsAsFactors = FALSE)

seminaturalalldata<-read.csv(file= "seminatural data.csv", header=T, stringsAsFactors = FALSE)

urbanalldata<-read.csv(file= "urban data.csv", header=T, stringsAsFactors = FALSE)


farmlandfull<-farmlandalldata %>%
  mutate(host =coalesce(host_species,host_order, host_class), vector=coalesce(diptera_species, diptera_genus))

farmlandweb<- frame2webs(farmlandfull, varnames = c("host", "diptera_genus", "webID", "frequency"))


seminaturalfull<-seminaturalalldata %>%
  mutate(host =coalesce(host_species,host_order, host_class), vector=coalesce(diptera_species, diptera_genus))

seminaturalweb<- frame2webs(seminaturalfull, varnames = c("host", "diptera_genus", "webID", "frequency"))


urbanfull<-urbanalldata %>%
  mutate(host =coalesce(host_species,host_order, host_class), vector=coalesce(diptera_species, diptera_genus))

urbanweb<- frame2webs(urbanfull, varnames = c("host", "diptera_genus", "webID", "frequency"))


plotweb(entireweb$Pall,text.rot = 90, high.lablength = 0, low.lablength = 0)

plotweb(farmlandweb$Pall, text.rot = 90, high.lablength = 0, low.lablength = 0 )

plotweb(seminaturalweb$Pall, text.rot = 90,
        x.lim = c(0.22,6.5), high.lablength = 0, low.lablength = 0)
plotweb(urbanweb$Pall, text.rot = 90,high.lablength = 0, low.lablength = 0 )



# Accumulation figure

colnames(allhabitatalldata)<-c("Urban host-biting interactions", "Urban biting Diptera","Urban hosts",
                               "Near-natural host-biting interactions","Near-natural biting Diptera",
                               "Near-natural hosts", "Agricultural host-biting interactions",
                               "Agricultural biting Diptera", "Agricultural hosts")

allhabitatalldata <- iNEXT(allhabitatalldata, q=0, datatype="incidence_freq", se=F)



allhabitatdipterainext<-ggiNEXT(allhabitatalldata, type=1, color.var="site") + 
  theme_bw() + 
  xlab("Number of blood meals")+
  ylab("Number of species or interactions")+ 
  guides(linetype=guide_legend(title="Method"),
         colour=guide_legend(title=""), 
         fill=guide_legend(title=""), 
         shape=guide_legend(title="")) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks.length = unit(0.25, "cm"),
        axis.line.x = element_line(colour="black", size = 0.1),
        axis.line.y = element_line(colour="black", size = 0.1),
        axis.text = element_text(colour = "black", size = 12),
        axis.title = element_text(colour ="black", size = 12))



allhabitatdipterainext<-allhabitatdipterainext + scale_colour_manual(values=c("cornflowerblue", "cornflowerblue", 
                                                                              "cornflowerblue","sienna4",
                                                                              "sienna4","sienna4", "gray29",
                                                                              "gray29","gray29"))

allhabitatdipterainext<-allhabitatdipterainext+scale_shape_manual(values=c(16,15,17,16,15,17,16,15,17))



# Domestic animal loops, removing random host species equal to number of domestic species in the network

# Agricultural loop- data is subset by habitat type

farmlandallspecies<-read.csv(file="farmlandnodomtest.csv", header=T,stringsAsFactors = F )


farmlandallspecies<-farmlandallspecies %>%
  mutate(host =coalesce(host_species,host_order, host_class), vector=coalesce(diptera_species, diptera_genus))

specieslist<-c("H.sapiens", "C.a. hircus", "O.aries", "C.l. familiaris", "G.gallus", "B.taurus", "E.ferus", "F.catus" )


splitfarmland<-split( farmlandallspecies , f = farmlandallspecies$webID )


farmdata<-c()
spnb<-list()
matchspm<-list()
matchlen<-list()
samplesremove<-list()
farmdatawithremoved<-data.frame()
farmlist=list()
farmweb=c()
farmweb1=c()
evetest<-c()
evetest2<-data.frame()



for (Dm in 1:length(splitfarmland)) {
  
  farmdata<-splitfarmland[[Dm]]
  # number of unique host species
  spnb<-unique(farmdata$host_species)
  # match the host species to domestic list
  matchspm<-match(spnb, specieslist)
  #turn nas to 0 in the match list
  matchspm[is.na(matchspm)] <- 0
  #find the number of domestic animals in the network
  matchlen<-length(matchspm[matchspm>0])
  #select species to remove by number of domestic animals
  samplesremove<-sample (spnb, matchlen, replace=F)
  #samplesremove<-replicate(100, sample (spnb, matchlen, replace=F))
  
  #remove the rows with that species name
  # farmdatawithremoved<-farmdata[ ! farmdata$host_species %in% samplesremove, ]
  farmdatawithremoved<-data.frame(farmdata[ ! farmdata$host_species %in% samplesremove, ])
  farmweb<-frame2webs(farmdatawithremoved, varnames =c("host_species", "diptera_species", "webID", "frequency"))
  farmlist[[Dm]]<-farmdatawithremoved
  #farmweb1[[Dm]]<-farmweb
  evetest[Dm]<-as.data.frame(lapply(farmweb, networklevel, index="interaction evenness"))
  #evetest1<-as.data.frame(t(evetest))
  evetest2<-rbind.data.frame(evetest2, evetest)
  
  
}




#Calculating z-scores

farmnrandomnodommeans<-colMeans(evetest2[sapply(evetest2, is.numeric)])


zscoretest<-farmnrandomnodommeans
targetzscore<-farmlandmetric


farmzscores<-(targetzscore-mean(zscoretest))/sd(zscoretest)






# Urban- data is subset by habitat type

rbanallspecies<-read.csv(file="urbannodomtest.csv", header=T,stringsAsFactors = F )


urbanallspecies<-urbanallspecies %>%
  mutate(host =coalesce(host_species,host_order, host_class), vector=coalesce(diptera_species, diptera_genus))


spliturban<-split( urbanallspecies , f = urbanallspecies$webID )


reps<-100
urbandata<-c()
urbanspnb<-list()
urbanmatchspm<-list()
urbanmatchlen<-list()
urbansamplesremove<-list()
urbandatawithremoved<-data.frame()
urbanist=list()
urbanweb=c()
urbanweb1=c()
urbanevetest<-c()
urbanevetest1<-as.data.frame(matrix(0, ncol = 2, nrow = 120))


urbanevetest2<-data.frame()


for (Dm in 1:length(spliturban)) {
  
  urbandata<-spliturban[[Dm]]
  # number of unique host species
  urbanspnb<-unique(urbandata$host_species)
  # match the host species to domestic list
  urbanmatchspm<-match(urbanspnb, specieslist)
  #turn nas to 0 in the match list
  urbanmatchspm[is.na(urbanmatchspm)] <- 0
  #find the number of domestic animals in the network
  urbanmatchlen<-length(urbanmatchspm[urbanmatchspm>0])
  #select species to remove by number of domestic animals
  urbansamplesremove<-sample (urbanspnb, urbanmatchlen, replace=F)
  
  #remove the rows with that species name
  urbandatawithremoved<-data.frame(urbandata[ ! urbandata$host_species %in% urbansamplesremove, ])
  urbanweb<-frame2webs(urbandatawithremoved, varnames =c("host_species", "diptera_species", "webID", "frequency"))
  urbanist[[Dm]]<-urbandatawithremoved
  urbanevetest[Dm]<-as.data.frame(lapply(urbanweb, networklevel, index="interaction evenness"))
  urbanevetest2<-rbind.data.frame(urbanevetest2, urbanevetest)
  
  
}




#Calculating  z-scores

urbanrandomdommeans<-colMeans(urbanevetest2[sapply(urbanevetest2, is.numeric)])


urbantargetzscore<-urbannodomesticmetric
urbanzscoretest<-urbanrandomdommeans


urbanzscores<-(urbantargetzscore-mean(urbanzscoretest))/sd(urbanzscoretest)



# Near natural - data is subset by habitat type

semiallspecies<-read.csv(file="seminodomtest.csv", header=T,stringsAsFactors = F )


semiallspecies<-semiallspecies %>%
  mutate(host =coalesce(host_species,host_order, host_class), vector=coalesce(diptera_species, diptera_genus))


splitsemi<-split( semiallspecies , f = semiallspecies$webID )


semidata<-c()
semispnb<-list()
semimatchspm<-list()
semimatchlen<-list()
semisamplesremove<-list()
semidatawithremoved<-data.frame()
semilist=list()
semiweb=c()
semiweb1=c()
semievetest<-c()
semievetest1<-as.data.frame(matrix(0, ncol = 2, nrow = 120))

semievetest2<-data.frame()


for (Dm in 1:length(splitsemi)) {
  
  semidata<-splitsemi[[Dm]]
  # number of unique host species
  semispnb<-unique(semidata$host_species)
  # match the host species to domestic list
  semimatchspm<-match(semispnb, specieslist)
  #turn nas to 0 in the match list
  semimatchspm[is.na(semimatchspm)] <- 0
  #find the number of domestic animals in the network
  semimatchlen<-length(semimatchspm[semimatchspm>0])
  #select species to remove by number of domestic animals
  semisamplesremove<-sample (semispnb, semimatchlen, replace=F)
  
  #remove the rows with that species name
  semidatawithremoved<-data.frame(semidata[ ! semidata$host_species %in% semisamplesremove, ])
  semiweb<-frame2webs(semidatawithremoved, varnames =c("host_species", "diptera_species", "webID", "frequency"))
  semilist[[Dm]]<-semidatawithremoved
  semievetest[Dm]<-as.data.frame(lapply(semiweb, networklevel, index="interaction evenness"))
  semievetest2<-rbind.data.frame(semievetest2, semievetest)
  
}





#Calculating z-scores

semirandomdommeans<-colMeans(semievetest2[sapply(semievetest2, is.numeric)])



semizscoretest<-semirandomdommeans
semitargetzscore<-seminodommetric


semizscores<-(semitargetzscore-mean(semizscoretest))/sd(semizscoretest)




# plotting null model data with empirical values


par(oma = c(4, 4, 0.2, 0.2), mar=c(0.1, 0.1, 1, 1), mfrow=c(3,1)) 



evetest3<-evetest2

boxplot(evetest3, xaxt="n", ylim=c(0, 0.9))
points(farmlandmetric, col="red", cex=1.5, pch = 15)
title(xlab = "Networks", ylab = "Evenness value")

urbanevetest3<-urbanevetest2

boxplot(urbanevetest3, xaxt="n")
points(urbannodomesticmetric, col="red", cex=1.5, pch = 15)


semievetest3<-semievetest2

boxplot(semievetest3, xaxt="n")
points(seminodommetric, col="red", cex=1.5, pch = 15)


mtext("Network",side=1,line=2,outer=TRUE,cex=1.3)
mtext("Evenness value",side=2,line=2.2,outer=TRUE,cex=1.3,las=0)



# Calculating IE for networks with targeted removal of domestic species

#Agricultural- data is subset to remove domestic hosts

farmnodomestic<-read.csv(file= "farmland no domestic data 3.csv", header=T, stringsAsFactors = F)


farmnodomestic<-farmnodomestic %>%
  mutate(host =coalesce(host_species,host_order, host_class), vector=coalesce(diptera_species, diptera_genus))

farmnodomesticweb<- frame2webs(farmnodomestic, varnames = c("host", "vector", "webID", "frequency"))

farmlandmetric<-as.data.frame(lapply(farmnodomesticweb, networklevel, index="interaction evenness"))


#Near natural- data is subset to remove domestic hosts

seminaturalnodomestic<-read.csv(file= "semi natural no domestic 1.csv", header=T, stringsAsFactors = F)


seminaturalnodomestic<-seminaturalnodomestic %>%
  mutate(host =coalesce(host_species,host_order, host_class), vector=coalesce(diptera_species, diptera_genus))

seminaturalnodomesticweb<- frame2webs(seminaturalnodomestic, varnames = c("host", "vector", "webID", "frequency"))

seminodommetric<-as.data.frame(lapply(seminaturalnodomesticweb, networklevel, index="interaction evenness"))


#Urban- data is subset to remove domestic hosts
urbannodomestic<-read.csv(file= "urban no domestic 1.csv", header=T, stringsAsFactors = F)


urbannodomestic<-urbannodomestic %>%
  mutate(host =coalesce(host_species,host_order, host_class), vector=coalesce(diptera_species, diptera_genus))

urbannodomesticweb<- frame2webs(urbannodomestic, varnames = c("host", "vector", "webID", "frequency"))

urbannodomesticmetric<-as.data.frame(lapply(urbannodomesticweb, networklevel, index="interaction evenness"))










