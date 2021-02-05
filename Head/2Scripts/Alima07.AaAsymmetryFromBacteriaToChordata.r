rm(list=ls(all=TRUE)) 

######################################################################
####### PARSING OF VALERIA FILES: make datasets AlphaProteoBacteria, Fungi, Invertebrates
######################################################################

######## parsing of AlphaProteoBacteria: all genes except ATP8, ND6 & ND4L (ATP6 COX1 COX2 COX3 CYTB  ND1  ND2  ND3  ND4  ND5)
data = read.table("../../Body/1Raw/BestiariiFromValeriya/AlphaproteobacteriaAll.txt", sep = ';', header = TRUE, quote = '')
agg = aggregate(data[,c(5:24)], by = list(data$Species,data$Gene), FUN = mean)
names(agg)=c('Species','Gene',names(data)[c(5:24)])
species = data.frame(table(agg$Species));
VecOfGenes =  data.frame(table(agg$Gene));
VecOfSpeciesWithTenGenes = species[species$Freq == 10,]$Var1; length(VecOfSpeciesWithTenGenes)
agg = agg[agg$Species %in% VecOfSpeciesWithTenGenes,]
agg = agg[order(agg$Species,agg$Gene),]
nrow(agg)
length(unique(agg$Species)) 
table(agg$Gene)
AlphaProteoBacteria = agg
dim(AlphaProteoBacteria)

######## parsing of Fungi
data = read.table("../../Body/1Raw/BestiariiFromValeriya/FungiAll.txt", sep = ';', header = TRUE, quote = '')
agg = aggregate(data[,c(5:24)], by = list(data$Species,data$Gene), FUN = mean)
names(agg)=c('Species','Gene',names(data)[c(5:24)])
agg$Species = gsub('\\[','',agg$Species)
species = data.frame(table(agg$Species));
VecOfGenes =  data.frame(table(agg$Gene));
VecOfSpeciesWithTenGenes = species[species$Freq == 10,]$Var1; length(VecOfSpeciesWithTenGenes)
agg = agg[agg$Species %in% VecOfSpeciesWithTenGenes,]
agg = agg[order(agg$Species,agg$Gene),]
nrow(agg)
length(unique(agg$Species)) 
table(agg$Gene)
Fungi = agg
dim(Fungi)

######## parsing of Invertebrate
data = read.table("../../Body/1Raw/BestiariiFromValeriya/InvertebrateAll.txt", sep = ';', header = TRUE, quote = '')
agg = aggregate(data[,c(5:24)], by = list(data$Species,data$Gene), FUN = mean)
names(agg)=c('Species','Gene',names(data)[c(5:24)])
agg$Species = gsub('\\[','',agg$Species)
VecOfGenes =  data.frame(table(agg$Gene)); # VecOfGenes$Var1
#VecOfGenes = c('ATP6','COX1','COX2','COX3','CYTB','ND1','ND2','ND3','ND4','ND4L','ND5','ND6')
#VecOfGenes = c('ATP6','COX1','COX2','COX3','CYTB','ND1','ND2','ND3','ND4','ND4L','ND5')
VecOfGenes = c('ATP6','COX1','COX2','COX3','CYTB','ND1','ND2','ND3','ND4','ND5')
agg = agg[agg$Gene %in% VecOfGenes,]
species = data.frame(table(agg$Species));
VecOfSpecies = species[species$Freq == 10,]$Var1; length(VecOfSpecies)
agg = agg[agg$Species %in% VecOfSpecies,]
agg = agg[order(agg$Species,agg$Gene),]
nrow(agg) # 1199/11 gens = 109 species
length(unique(agg$Species)) 
table(agg$Gene)
Invertebrates = agg
dim(Invertebrates)

######## parsing of Plants
data = read.table("../../Body/1Raw/BestiariiFromValeriya/PlantsAll.txt", sep = ';', header = TRUE, quote = '')
agg = aggregate(data[,c(5:24)], by = list(data$Species,data$Gene), FUN = mean)
names(agg)=c('Species','Gene',names(data)[c(5:24)])
VecOfGenes =  data.frame(table(agg$Gene)); # VecOfGenes$Var1
#VecOfGenes = c('ATP6','COX1','COX2','COX3','CYTB','ND1','ND2','ND3','ND4','ND4L','ND5','ND6')
#VecOfGenes = c('ATP6','COX1','COX2','COX3','CYTB','ND1','ND2','ND3','ND4','ND4L','ND5')
VecOfGenes = c('ATP6','COX1','COX2','COX3','CYTB','ND1','ND2','ND3','ND4','ND5')
agg = agg[agg$Gene %in% VecOfGenes,]
species = data.frame(table(agg$Species));
VecOfSpecies = species[species$Freq == 10,]$Var1; length(VecOfSpecies)
agg = agg[agg$Species %in% VecOfSpecies,]
agg = agg[order(agg$Species,agg$Gene),]
nrow(agg) # 1199/11 gens = 109 species
length(unique(agg$Species)) 
table(agg$Gene)
Plants = agg
dim(Plants)

######################################################################
####### ANALYSIS (we compare ALL species using only 10genes)
######################################################################
VecOf10UniversalGenes = c('ATP6','COX1','COX2','COX3','CYTB','ND1','ND2','ND3','ND4','ND5')
## HERE SUBSTITUTE IT BY AlphaProteoBacteria or Fungi or Invertebrates or Plants
# data = AlphaProteoBacteria  
# data = Fungi
data = Invertebrates  
# data = Plants
data = data[data$Gene %in% VecOf10UniversalGenes,]
data$Gainers = data$Pro + data$Thr + data$His + data$Gln + data$Asn + data$Lys # 16 codons 
data$Losers = data$Phe + data$Val + data$Gly + data$Cys + data$Trp  # 12 codons
data$All = apply(as.matrix(data[,c(3:22)]),1,FUN = sum)
data$Intermediate = data$All - data$Gainers - data$Losers # 64 - 16 - 12 - 4 (Stops) = 32

agg = aggregate(list(data$Gainers,data$Losers, data$All), by = list(data$Species), FUN = sum)
names(agg)=c('Species','Gainers','Losers','All')
agg$FrOfGainers = agg$Gainers / agg$All
agg$FrOfLosers = agg$Losers / agg$All

boxplot(agg$FrOfGainers,agg$FrOfLosers)
median(c(agg$FrOfGainers,agg$FrOfLosers)) # 0.2337901

######################################################################
####### PLOTTING BY ALINA
######################################################################

####### plotting

library(ggplot2)
library(ggbeeswarm) 
library(cowplot)
library(tidyr)
library(dplyr)

pdf("../../Body/4Figures/Alima07.AaAsymmetryBacteriaBlastP.r.pdf")
data = read.table("../../Body/1Raw/AminoAcidsFromBacteria.csv", sep = ';', header = TRUE, quote = '')
descr = data.frame(table(data$GeneName))

data$Gainers = data$Pro + data$Thr + data$His + data$Gln + data$Asn + data$Lys # 16 codons 
data$Loosers = data$Phe + data$Val + data$Gly + data$Cys + data$Trp  # 12 codons
names(data)
data$All = apply(as.matrix(data[,c(4:23)]),1,FUN = sum)
data$Intermediate = data$All - data$Gainers - data$Loosers # 64 - 16 - 12 - 4 (Stops) = 32

data$FrOfGainers = data$Gainers / data$All
data$FrOfLoosers = data$Loosers / data$All

Agg1 = aggregate(list(data$FrOfGainers,data$FrOfLoosers), by = list(data$GeneName), FUN = mean)
names(Agg1)=c('GeneName','FrOfGainers','FrOfLoosers')
summary(Agg1$FrOfGainers)
summary(Agg1$FrOfLoosers)

boxplot(Agg1$FrOfGainers,Agg1$FrOfLoosers, col = c(rgb(1,0.1,0.1,0.5),rgb(0.1,0.1,0.1,0.5)), names = c('gainers','losers'), main = 'alphaproteobacteria') # ALINA!!!
dev.off()

AggToPlot = Agg1 %>%
  gather(key = 'GainersOrLoosers', value = 'Fraction', FrOfGainers:FrOfLoosers)

bacteria = ggplot(AggToPlot, aes(as.factor(GainersOrLoosers), Fraction, colour = as.factor(GainersOrLoosers))) +
  geom_quasirandom(shape = 1, cex = 0.5) + 
  stat_summary(aes(group = factor(GainersOrLoosers)), 
               fun.y = 'median', geom = 'point', shape = 8, size = 2, col = 'midnightblue') + # to show a median
  theme_minimal() + # поменять серый фон на белый
  theme(axis.text.x = element_blank()) +
  scale_color_manual(name = '', labels = c('Gainers', 'Losers'), 
                     values = c(rgb(1,0.1,0.1,0.5), rgb(0.1,0.1,0.1,0.5))) +
  # scale_fill_discrete(labels = c('Fraction of gainers', 'Fraction of losers'))
  labs(x = '', y = '', title = 'alphaproteobacteria')

save_plot('../../Body/4Figures/Alima07.AaAsymmetryBacteriaBlastP01.r.pdf', bacteria)
