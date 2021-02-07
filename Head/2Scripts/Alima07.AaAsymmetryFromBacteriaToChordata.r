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

# reading table with vertebrates
Vertebrates = read.table("../../Body/1Raw/AminoAcidFreqsChordata.txt", sep = ',', header = TRUE, quote = '')
colnames(Vertebrates)
Vertebrates = Vertebrates[,-c(23:33)]
names(Vertebrates)[23] = 'Taxon'

# adding taxons to all other tables
AlphaProteoBacteria[, 'Taxon'] = 'AlphaProteoBacteria'
Fungi[, 'Taxon'] = 'Fungi'
Invertebrates[, 'Taxon'] = 'Invertebrates'
Plants[, 'Taxon'] = 'Plants'

bestiari = bind_rows(AlphaProteoBacteria, Fungi, Invertebrates, Plants, Vertebrates) %>%
  filter(Gene %in% VecOfGenes) %>%
  mutate(
    Gainers = Pro + Thr + His + Gln + Asn + Lys,
    Losers = Phe + Val + Gly + Cys + Trp,
    All = rowSums(.[3:22])) %>%
  select(Species, Gainers, Losers, All, Taxon) %>%
  group_by(Species) %>% 
  summarise(GainersSum = sum(Gainers), LosersSum = sum(Losers), AllSum = sum(All), 
            Taxon = Taxon) %>%
  mutate(FrOfGainers = GainersSum/AllSum,
         FrOfLosers = LosersSum/AllSum) %>%
  distinct()

table(bestiari$Taxon)

# Actinopterygii AlphaProteoBacteria            Amphibia         AncientFish                Aves 
# 1770                 770                 205                 126                 432 
# Fungi       Invertebrates            Mammalia              Plants            Reptilia 
# 167                 104                 788                  84                 269 

# which taxa have more losers ? 
moreLosers = bestiari[bestiari$FrOfGainers < bestiari$FrOfLosers,]

table(moreLosers$Taxon)
# Actinopterygii AlphaProteoBacteria               Fungi       Invertebrates              Plants 
# 6                 770                 167                  90                  82 

# All bacteria, fungi, most of the invertebrates and plants, few fish have more losers

DataToPlot = bestiari %>%
  filter(Taxon != 'AncientFish') %>%
  select(Species, FrOfGainers, FrOfLosers, Taxon) %>%
  gather(key = 'GainersOrLosers', value = 'Fraction', FrOfGainers:FrOfLosers)

alltaxa = ggplot(DataToPlot, aes(as.factor(Taxon), Fraction, colour = as.factor(GainersOrLosers))) +
  geom_quasirandom(shape = 1, cex = 1) + 
  stat_summary(aes(group = factor(GainersOrLosers)), 
               fun.y = 'median', geom = 'point', shape = 8, size = 2, col = 'midnightblue') + # to show a median
  theme_minimal() + # поменять серый фон на белый
  theme(axis.text.x = element_text(color = "black", size = 12), 
        axis.text.y = element_text(color = "black", size = 12)) +
  scale_color_manual(name = '', labels = c('Gainers', 'Losers'), 
                     values = c(rgb(1,0.1,0.1,0.5), rgb(0.1,0.1,0.1,0.5))) +
  scale_x_discrete(limits = c('AlphaProteoBacteria', 'Fungi', 'Plants', 'Invertebrates',
                              'Actinopterygii', 'Amphibia', 'Reptilia',
                              'Mammalia', 'Aves')) +
  # scale_fill_discrete(labels = c('Fraction of gainers', 'Fraction of losers'))
  labs(x = '', y = '')


save_plot('../../Body/4Figures/Alima07.AaAsymmetryFromBacteriaToChordata.r.pdf', alltaxa,
          base_height = 8)
