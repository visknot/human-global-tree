rm(list=ls(all=TRUE)) 

library(ggplot2)
library(ggbeeswarm) 
library(cowplot)
library(tidyr)
library(dplyr)

pdf("../../Body/4Figures/Alima06.AaAsymmetryChordataRefSeqs.r.pdf")
data = read.table("../../Body/1Raw/AminoAcidFreqsChordata.txt", sep = ',', header = TRUE, quote = '')
colnames(data)
data = data[,-c(23:33)]
data$Gainers = data$Pro + data$Thr + data$His + data$Gln + data$Asn + data$Lys # 16 codons 
data$Loosers = data$Phe + data$Val + data$Gly + data$Cys + data$Trp  # 14 codons (not 16, because I don't discriminate between LeuCT and LeuTT and use both of them as Leu - intermediate)
data$All = apply(as.matrix(data[,c(3:22)]),1,FUN = sum)
data$Intermediate = data$All - data$Gainers - data$Loosers # 64 - 16 - 14 - 4 (Stops) = 30

nrow(data) # 46670
length(unique(data$Species)) # 3590 species


###################### TASK 1: FRACTION OF GAINERS VS FRACTION OF LOOSERS
###################
#### 12 genes
###################
data1 = data[data$Gene != 'ND6',]
Agg1 = aggregate(list(data1$Gainers,data1$Loosers,data1$Intermediate,data1$All), by = list(data1$Species,data1$Class), FUN = sum)
names(Agg1) = c('Species','Class','Gainers','Loosers','Intermediate','All')
Agg1$GainersToLosers = Agg1$Gainers - Agg1$Loosers
Agg1$FrOfGainers = Agg1$Gainers/Agg1$All
Agg1$FrOfLoosers = Agg1$Loosers/Agg1$All

## 1 FractionOfGainers and FractionOfLoosers in general and in different classes

boxplot(Agg1$FrOfGainers,Agg1$FrOfLoosers, col = c(rgb(1,0.1,0.1,0.5),rgb(0.1,0.1,0.1,0.5)), notch = TRUE, outline = FALSE, main = '12 genes, all chordata') # dev.off() 
# boxplot(Agg1$Gainers/16,Agg1$Loosers/12,Agg1$Intermediate/32) 
# table(Agg1$Class) dev.off()
Cold = c('Actinopterygii','Amphibia','Reptilia')
Warm = c('Aves','Mammalia')

par(mfrow=c(1,2))
boxplot(Agg1[Agg1$Class %in% Cold,]$FrOfGainers,Agg1[Agg1$Class %in% Cold,]$FrOfLoosers, ylim = c(0.16,0.31), names = c('gainers','loosers'), main = 'Coldblood', ylab = 'NumberOfAaPer12Genes', notch = TRUE, outline = FALSE, col = c(rgb(1,0.1,0.1,0.5),rgb(0.1,0.1,0.1,0.5)))
boxplot(Agg1[Agg1$Class %in% Warm,]$FrOfGainers,Agg1[Agg1$Class %in% Warm,]$FrOfLoosers, ylim = c(0.16,0.31), names = c('gainers','loosers'), main = 'Warmblood', notch = TRUE, outline = FALSE, col = c(rgb(1,0.1,0.1,0.5),rgb(0.1,0.1,0.1,0.5)))

par(mfrow=c(1,3))
boxplot(Agg1[Agg1$Class == 'Actinopterygii',]$FrOfGainers,Agg1[Agg1$Class == 'Actinopterygii',]$FrOfLoosers, ylim = c(0.16,0.31), names = c('gainers','loosers'), main = 'Actinopterygii', ylab = 'NumberOfAaPer12Genes', notch = TRUE, outline = FALSE,col = c(rgb(1,0.1,0.1,0.5),rgb(0.1,0.1,0.1,0.5)))
boxplot(Agg1[Agg1$Class == 'Amphibia',]$FrOfGainers,Agg1[Agg1$Class == 'Amphibia',]$FrOfLoosers, ylim = c(0.16,0.31), names = c('gainers','loosers'), main = 'Amphibia', notch = TRUE, outline = FALSE, col = c(rgb(1,0.1,0.1,0.5),rgb(0.1,0.1,0.1,0.5)))
boxplot(Agg1[Agg1$Class == 'Reptilia',]$FrOfGainers,Agg1[Agg1$Class == 'Reptilia',]$FrOfLoosers, ylim = c(0.16,0.31), names = c('gainers','loosers'), main = 'Reptilia', notch = TRUE, outline = FALSE, col = c(rgb(1,0.1,0.1,0.5),rgb(0.1,0.1,0.1,0.5)))

# "ALINA - 12 GENES"
par(mfrow=c(1,5))
boxplot(Agg1[Agg1$Class == 'Actinopterygii',]$FrOfGainers,Agg1[Agg1$Class == 'Actinopterygii',]$FrOfLoosers, ylim = c(0.16,0.31), names = c('gainers','loosers'), main = 'Actinopterygii', ylab = 'NumberOfAaPer12Genes', notch = TRUE, outline = FALSE,col = c(rgb(1,0.1,0.1,0.5),rgb(0.1,0.1,0.1,0.5)))
boxplot(Agg1[Agg1$Class == 'Amphibia',]$FrOfGainers,Agg1[Agg1$Class == 'Amphibia',]$FrOfLoosers, ylim = c(0.16,0.31), names = c('gainers','loosers'), main = 'Amphibia', notch = TRUE, outline = FALSE, col = c(rgb(1,0.1,0.1,0.5),rgb(0.1,0.1,0.1,0.5)))
boxplot(Agg1[Agg1$Class == 'Reptilia',]$FrOfGainers,Agg1[Agg1$Class == 'Reptilia',]$FrOfLoosers, ylim = c(0.16,0.31), names = c('gainers','loosers'), main = 'Reptilia', notch = TRUE, outline = FALSE, col = c(rgb(1,0.1,0.1,0.5),rgb(0.1,0.1,0.1,0.5)))
boxplot(Agg1[Agg1$Class == 'Mammalia',]$FrOfGainers,Agg1[Agg1$Class == 'Mammalia',]$FrOfLoosers, ylim = c(0.16,0.31), names = c('gainers','loosers'), main = 'Mammalia', notch = TRUE, outline = FALSE, col = c(rgb(1,0.1,0.1,0.5),rgb(0.1,0.1,0.1,0.5)))
boxplot(Agg1[Agg1$Class == 'Aves',]$FrOfGainers,Agg1[Agg1$Class == 'Aves',]$FrOfLoosers, ylim = c(0.16,0.31), names = c('gainers','loosers'), main = 'Aves', notch = TRUE, outline = FALSE, col = c(rgb(1,0.1,0.1,0.5),rgb(0.1,0.1,0.1,0.5)))


###

AggToPlot = Agg1 %>%
  select(Species, FrOfGainers, FrOfLoosers, Class) %>%
  filter(Class != 'AncientFish') %>%
  gather(key = 'GainersOrLoosers', value = 'Fraction', FrOfGainers:FrOfLoosers) %>%
  mutate(
    WarmOrCold = as.factor(case_when(.$Class %in% c('Mammalia', 'Aves') ~ 1,
                                     .$Class %in% c('Actinopterygii', 'Amphibia', 'Reptilia') ~ 0))
  )

all12Genes = ggplot(AggToPlot, aes(as.factor(Class), Fraction, colour = as.factor(GainersOrLoosers))) +
  geom_quasirandom(shape = 1, cex = 0.5) + 
  stat_summary(aes(group = factor(GainersOrLoosers)), 
               fun.y = 'median', geom = 'point', shape = 8, size = 2, col = 'midnightblue') + # to show a median
  theme_minimal() + # поменять серый фон на белый
  theme(axis.text.x = element_text(color = "black", size = 12), 
        axis.text.y = element_text(color = "black", size = 12)) +
  scale_color_manual(name = '', labels = c('Gainers', 'Losers'), 
                     values = c(rgb(1,0.1,0.1,0.5), rgb(0.1,0.1,0.1,0.5))) +
  scale_x_discrete(labels = c('Actinopterygii', 'Amphibia', 'Reptilia',
                              'Mammalia', 'Aves')) +
  # scale_fill_discrete(labels = c('Fraction of gainers', 'Fraction of losers'))
  labs(x = '', y = '12 genes')

### genes separately

data1$FrOfGainers = data1$Gainers/data1$All
data1 = data1[data1$Class != 'AncientFish',]

gene_order = ggplot(data1, aes(as.factor(Gene), FrOfGainers)) + 
   geom_quasirandom(shape = 1, cex = 0.5, color = rgb(1,0.1,0.1,0.5)) +
   scale_x_discrete(limits = c('COX1', 'COX2', 'ATP8', 'ATP6', 'COX3', 'ND3', 'ND4L', 
                               'ND4', 'ND5', 'CytB', 'ND1', 'ND2')) +
  facet_grid(as.factor(Class) ~ .)

###################
#### ND6
###################
data1 = data[data$Gene == 'ND6',]
Agg1 = aggregate(list(data1$Gainers,data1$Loosers,data1$Intermediate,data1$All), by = list(data1$Species,data1$Class), FUN = sum)
names(Agg1) = c('Species','Class','Gainers','Loosers','Intermediate','All')
Agg1$GainersToLosers = Agg1$Gainers - Agg1$Loosers
Agg1$FrOfGainers = Agg1$Gainers/Agg1$All
Agg1$FrOfLoosers = Agg1$Loosers/Agg1$All

## 1 FractionOfGainers and FractionOfLoosers in general and in different classes

par(mfrow=c(1,1))
boxplot(Agg1$FrOfGainers,Agg1$FrOfLoosers, col = c(rgb(1,0.1,0.1,0.5),rgb(0.1,0.1,0.1,0.5), ylim = c(0,0.55)), notch = TRUE, outline = FALSE) # dev.off() 
# boxplot(Agg1$Gainers/16,Agg1$Loosers/12,Agg1$Intermediate/32) 
# table(Agg1$Class)
Cold = c('Actinopterygii','Amphibia','Reptilia')
Warm = c('Aves','Mammalia')

par(mfrow=c(1,2))
boxplot(Agg1[Agg1$Class %in% Cold,]$FrOfGainers,Agg1[Agg1$Class %in% Cold,]$FrOfLoosers, ylim = c(0,0.55), names = c('gainers','loosers'), main = 'Coldblood', ylab = 'NumberOfAaPer12Genes', notch = TRUE, outline = FALSE, col = c(rgb(1,0.1,0.1,0.5),rgb(0.1,0.1,0.1,0.5)))
boxplot(Agg1[Agg1$Class %in% Warm,]$FrOfGainers,Agg1[Agg1$Class %in% Warm,]$FrOfLoosers, ylim = c(0,0.55), names = c('gainers','loosers'), main = 'Warmblood', notch = TRUE, outline = FALSE, col = c(rgb(1,0.1,0.1,0.5),rgb(0.1,0.1,0.1,0.5)))

par(mfrow=c(1,3))
boxplot(Agg1[Agg1$Class == 'Actinopterygii',]$FrOfGainers,Agg1[Agg1$Class == 'Actinopterygii',]$FrOfLoosers, ylim = c(0,0.55), names = c('gainers','loosers'), main = 'Actinopterygii', ylab = 'NumberOfAaPer12Genes', notch = TRUE, outline = FALSE,col = c(rgb(1,0.1,0.1,0.5),rgb(0.1,0.1,0.1,0.5)))
boxplot(Agg1[Agg1$Class == 'Amphibia',]$FrOfGainers,Agg1[Agg1$Class == 'Amphibia',]$FrOfLoosers, ylim = c(0,0.55), names = c('gainers','loosers'), main = 'Amphibia', notch = TRUE, outline = FALSE, col = c(rgb(1,0.1,0.1,0.5),rgb(0.1,0.1,0.1,0.5)))
boxplot(Agg1[Agg1$Class == 'Reptilia',]$FrOfGainers,Agg1[Agg1$Class == 'Reptilia',]$FrOfLoosers, ylim = c(0,0.55), names = c('gainers','loosers'), main = 'Reptilia', notch = TRUE, outline = FALSE, col = c(rgb(1,0.1,0.1,0.5),rgb(0.1,0.1,0.1,0.5)))

par(mfrow=c(1,2))
boxplot(Agg1[Agg1$Class == 'Mammalia',]$FrOfGainers,Agg1[Agg1$Class == 'Mammalia',]$FrOfLoosers, ylim = c(0,0.55), names = c('gainers','loosers'), main = 'Mammalia', notch = TRUE, outline = FALSE, col = c(rgb(1,0.1,0.1,0.5),rgb(0.1,0.1,0.1,0.5)))
boxplot(Agg1[Agg1$Class == 'Aves',]$FrOfGainers,Agg1[Agg1$Class == 'Aves',]$FrOfLoosers, ylim = c(0,0.55), names = c('gainers','loosers'), main = 'Aves', notch = TRUE, outline = FALSE, col = c(rgb(1,0.1,0.1,0.5),rgb(0.1,0.1,0.1,0.5)))

# "ALINA - ND6"
par(mfrow=c(1,5))
boxplot(Agg1[Agg1$Class == 'Actinopterygii',]$FrOfGainers,Agg1[Agg1$Class == 'Actinopterygii',]$FrOfLoosers, ylim = c(0,0.55), names = c('gainers','loosers'), main = 'Actinopterygii', ylab = 'NumberOfAaPer12Genes', notch = TRUE, outline = FALSE,col = c(rgb(1,0.1,0.1,0.5),rgb(0.1,0.1,0.1,0.5)))
boxplot(Agg1[Agg1$Class == 'Amphibia',]$FrOfGainers,Agg1[Agg1$Class == 'Amphibia',]$FrOfLoosers, ylim = c(0,0.55), names = c('gainers','loosers'), main = 'Amphibia', notch = TRUE, outline = FALSE, col = c(rgb(1,0.1,0.1,0.5),rgb(0.1,0.1,0.1,0.5)))
boxplot(Agg1[Agg1$Class == 'Reptilia',]$FrOfGainers,Agg1[Agg1$Class == 'Reptilia',]$FrOfLoosers, ylim = c(0,0.55), names = c('gainers','loosers'), main = 'Reptilia', notch = TRUE, outline = FALSE, col = c(rgb(1,0.1,0.1,0.5),rgb(0.1,0.1,0.1,0.5)))
boxplot(Agg1[Agg1$Class == 'Mammalia',]$FrOfGainers,Agg1[Agg1$Class == 'Mammalia',]$FrOfLoosers, ylim = c(0,0.55), names = c('gainers','loosers'), main = 'Mammalia', notch = TRUE, outline = FALSE, col = c(rgb(1,0.1,0.1,0.5),rgb(0.1,0.1,0.1,0.5)))
boxplot(Agg1[Agg1$Class == 'Aves',]$FrOfGainers,Agg1[Agg1$Class == 'Aves',]$FrOfLoosers, ylim = c(0,0.55), names = c('gainers','loosers'), main = 'Aves', notch = TRUE, outline = FALSE, col = c(rgb(1,0.1,0.1,0.5),rgb(0.1,0.1,0.1,0.5)))

AggToPlot = Agg1 %>%
  select(Species, FrOfGainers, FrOfLoosers, Class) %>%
  filter(Class != 'AncientFish') %>%
  gather(key = 'GainersOrLoosers', value = 'Fraction', FrOfGainers:FrOfLoosers) %>%
  mutate(
    WarmOrCold = as.factor(case_when(.$Class %in% c('Mammalia', 'Aves') ~ 1,
                                     .$Class %in% c('Actinopterygii', 'Amphibia', 'Reptilia') ~ 0))
  )

nd6 = ggplot(AggToPlot, aes(as.factor(Class), Fraction, colour = as.factor(GainersOrLoosers))) +
  geom_quasirandom(shape = 1, cex = 0.5) + 
  stat_summary(aes(group = factor(GainersOrLoosers)), 
               fun.y = 'median', geom = 'point', shape = 8, size = 2, col = 'midnightblue') + # to show a median
  theme_minimal() + # поменять серый фон на белый
  theme(axis.text.x = element_text(color = "black", size = 12), 
        axis.text.y = element_text(color = "black", size = 12)) +
  scale_color_manual(name = '', labels = c('Gainers', 'Losers'), 
                     values = c(rgb(1,0.1,0.1,0.5), rgb(0.1,0.1,0.1,0.5))) +
  scale_x_discrete(labels = c('Actinopterygii', 'Amphibia', 'Reptilia',
                              'Mammalia', 'Aves')) +
  # scale_fill_discrete(labels = c('Fraction of gainers', 'Fraction of losers'))
  labs(x = '', y = 'ND6')

plots = plot_grid(all12Genes, nd6, nrow = 2)

save_plot('../../Body/4Figures/Alima06.AaAsymmetryChordataRefSeqs01.r.pdf', plots, nrow = 2,
          base_height = 5)

save_plot('../../Body/4Figures/Alima06.AaAsymmetryChordataRefSeqs_GeneOrder.r.pdf', gene_order,
          base_height = 5)

# old 
#par(mfrow=c(1,1))
#wilcox.test(Agg1[Agg1$Class == 'Mammalia',]$FrOfGainers,Agg1[Agg1$Class == 'Aves',]$FrOfGainers)
#hist(Agg1[Agg1$Class == 'Mammalia',]$FrOfGainers, xlim = c(0.2,0.35), ylim = c(0,250), col = rgb(0.2,0.2,0.2,0.5), main = 'fraction of gainers in warmblooded', xlab = '') # dev.off()
#par(new = TRUE)
#hist(Agg1[Agg1$Class == 'Aves',]$FrOfGainers, xlim = c(0.2,0.35), ylim = c(0,250), col =  rgb(1,0.1,0.1,0.5), main = '', xlab = '') # dev.off()
#hist(Agg1[Agg1$Class == 'Actinopterygii',]$FrOfGainers, xlim = c(0.2,0.35), ylim = c(0,250), col = rgb(0.1,0.1,1,0.5), main = 'fraction of gainers in coldblooded', xlab = '')
#par(new = TRUE)
#hist(Agg1[Agg1$Class == 'Amphibia',]$FrOfGainers, xlim = c(0.2,0.35), ylim = c(0,250), col =  rgb(0.1,1,0.1,0.5), main = '', xlab = '')
#par(new = TRUE)
#hist(Agg1[Agg1$Class == 'Reptilia',]$FrOfGainers, xlim = c(0.2,0.35), ylim = c(0,250), col =  rgb(1,0.1,0.1,0.5), main = '', xlab = '')


######################## TASK 2: Fr of gainers to loosers is higher in longlived mammals (aged oocytes with strong asymmetric damage!)
###################
#### 12 genes only (bigger statistics)
###################
data1 = data[data$Gene != 'ND6',]
Agg1 = aggregate(list(data1$Gainers,data1$Loosers,data1$Intermediate,data1$All), by = list(data1$Species,data1$Class), FUN = sum)
names(Agg1) = c('Species','Class','Gainers','Loosers','Intermediate','All')
Agg1$GainersToLosers = Agg1$Gainers - Agg1$Loosers
Agg1$FrOfGainers = Agg1$Gainers/Agg1$All
Agg1$FrOfLoosers = Agg1$Loosers/Agg1$All
Agg1$FrOfGainersMinusFrOfLoosers = Agg1$FrOfGainers - Agg1$FrOfLoosers

GL = read.table("../../Body/1Raw/GenerationLenghtforMammals.xlsx.txt", header = TRUE, sep = '\t')
GL$Species = gsub(' ','_',GL$Scientific_name)

par(mfrow=c(2,1))
Mammals = merge(Agg1,GL, by = 'Species')
cor.test(Mammals$FrOfGainers,Mammals$GenerationLength_d, method = 'spearman') #  "ALINA PICs" (MIN)
plot(log2(Mammals$GenerationLength_d),Mammals$FrOfGainers) # dev.off()

cor.test(Mammals$FrOfLoosers,Mammals$GenerationLength_d, method = 'spearman')  # "ALINA PICs"
plot(log2(Mammals$GenerationLength_d),Mammals$FrOfLoosers) # dev.off()


### PICs 

library(ape)
library(geiger)
library(caper)

tree <- read.tree("../../Body/1Raw/mtalign.aln.treefile.rooted")

Mammals = Mammals[-404, ] #remove duplicates
row.names(Mammals) = Mammals$Species

tree_w = treedata(tree, Mammals[, c('Species', 'FrOfGainers', 'FrOfLoosers', 'GenerationLength_d')], 
                  sort=T, warnings=T)$phy

data_tree = as.data.frame(treedata(tree_w, Mammals[, c('Species', 'FrOfGainers', 'FrOfLoosers', 'GenerationLength_d')], 
                             sort=T, warnings=T)$data)

data_tree$FrOfGainers = as.numeric(as.character(data_tree$FrOfGainers))
data_tree$GenerationLength_d = as.numeric(as.character(data_tree$GenerationLength_d))
data_tree$FrOfLoosers = as.numeric(as.character(data_tree$FrOfLoosers))

cor.test(pic(data_tree$FrOfGainers, tree_w), pic(data_tree$GenerationLength_d, tree_w), method = 'spearman')
# rho = 0.07723414, p-value = 0.04939

cor.test(pic(data_tree$FrOfGainers, tree_w), pic(log2(data_tree$GenerationLength_d), tree_w), method = 'spearman')
# rho = 0.08615036, p-value = 0.02832

plot(pic(data_tree$FrOfGainers, tree_w), pic(log2(data_tree$GenerationLength_d), tree_w))

cor.test(pic(data_tree$FrOfLoosers, tree_w), pic(log2(data_tree$GenerationLength_d), tree_w), method = 'spearman')
# rho = -0.05329016, p-value = 0.1755

###

MammComp = comparative.data(tree_w, data_tree, Species, vcv=TRUE)

model1 = pgls(log2(GenerationLength_d) ~ FrOfGainers, MammComp, lambda="ML")
summary(model1)

# lambda [ ML]  : 0.952

# Coefficients:
#   Estimate Std. Error t value  Pr(>|t|)    
# (Intercept)   6.3979     1.4141  4.5242 7.212e-06 ***
#   FrOfGainers  18.0385     5.3341  3.3818 0.0007636 ***
# Multiple R-squared: 0.01737,	Adjusted R-squared: 0.01585 


############################### TASK 3 fishes and temperature
###################
#### 12 genes only (bigger statistics)
###################
data1 = data[data$Gene != 'ND6',]
Agg1 = aggregate(list(data1$Gainers,data1$Loosers,data1$Intermediate,data1$All), by = list(data1$Species,data1$Class), FUN = sum)
names(Agg1) = c('Species','Class','Gainers','Loosers','Intermediate','All')
Agg1$GainersToLosers = Agg1$Gainers - Agg1$Loosers
Agg1$FrOfGainers = Agg1$Gainers/Agg1$All
Agg1$FrOfLoosers = Agg1$Loosers/Agg1$All
Agg1$FrOfGainersMinusFrOfLoosers = Agg1$FrOfGainers - Agg1$FrOfLoosers

Agg1 = Agg1[Agg1$Class == 'Actinopterygii',]
Temp = read.table("../../Body/1Raw/FishBaseTemperature.txt", header = TRUE)
Temp = aggregate(Temp$Temperature, by = list(Temp$Species), FUN = median)
names(Temp) = c('Species','Temperature')
nrow(Temp)

Fish = merge(Agg1,Temp, by = 'Species') # 305

par(mfrow=c(2,1))
cor.test(Fish$FrOfGainers,Fish$Temperature, method = 'spearman') # "ALINA PICs" (MIN)
plot(Fish$Temperature, Fish$FrOfGainers) # ALINA PICs

cor.test(Fish$FrOfLoosers,Fish$Temperature, method = 'spearman') # "ALINA PICs"
plot(Fish$Temperature, Fish$FrOfLoosers) # ALINA PICs


### PICs

library(ape)
library(geiger)
library(caper)

tree <- read.tree("../../Body/1Raw/mtalign.aln.treefile.rooted")

row.names(Fish) = Fish$Species

tree_w = treedata(tree, Fish, 
                  sort=T, warnings=T)$phy

data_tree = as.data.frame(treedata(tree_w, Fish, 
                                   sort=T, warnings=T)$data)

data_tree$FrOfGainers = as.numeric(as.character(data_tree$FrOfGainers))
data_tree$Temperature = as.numeric(as.character(data_tree$Temperature))
data_tree$FrOfLoosers = as.numeric(as.character(data_tree$FrOfLoosers))

cor.test(pic(data_tree$FrOfGainers, tree_w), pic(data_tree$Temperature, tree_w), method = 'spearman')
# rho = 0.1870207, p-value = 0.001052

cor.test(pic(data_tree$FrOfLoosers, tree_w), pic(data_tree$Temperature, tree_w), method = 'spearman')
# rho = -0.1908986, p-value = 0.0008212

plot(pic(data_tree$FrOfGainers, tree_w), pic(data_tree$Temperature, tree_w))
plot(pic(data_tree$FrOfLoosers, tree_w), pic(data_tree$Temperature, tree_w))

# PGLS

FishComp = comparative.data(tree_w, data_tree, Species, vcv=TRUE)

model1 = pgls(scale(Temperature) ~ scale(FrOfGainers), FishComp, lambda="ML")
summary(model1)

# lambda [ ML]  : 0.838
# Coefficients:
#   Estimate Std. Error t value  Pr(>|t|)    
# (Intercept)        -0.657014   0.402888 -1.6308 0.1039797    
# scale(FrOfGainers)  0.249319   0.066749  3.7352 0.0002241 ***
# Multiple R-squared: 0.04402,	Adjusted R-squared: 0.04086 

################### 5. gainers / loosers in second (only Ah>Gh) and fourth (only Ch>Th) quartiles in the relationship with Generation Length of Mammals (it seems that second quartile is stronger)
###############################
######### second quartile (works good)
###############################
data$Gainers = data$Pro
data$Loosers = data$Phe
data$All = data$Pro + data$Phe + data$Leu + data$Ser
data$Intermediate = data$All - data$Gainers - data$Loosers # 64 - 16 - 14 - 4 (Stops) = 30

#### 12 genes
data1 = data[data$Gene != 'ND6',]
Agg1 = aggregate(list(data1$Gainers,data1$Loosers,data1$Intermediate,data1$All), by = list(data1$Species,data1$Class), FUN = sum)
names(Agg1) = c('Species','Class','Gainers','Loosers','Intermediate','All')
Agg1$GainersToLosers = Agg1$Gainers - Agg1$Loosers
Agg1$FrOfGainers = Agg1$Gainers/Agg1$All
Agg1$FrOfLoosers = Agg1$Loosers/Agg1$All

GL = read.table("../../Body/1Raw/GenerationLenghtforMammals.xlsx.txt", header = TRUE, sep = '\t')
GL$Species = gsub(' ','_',GL$Scientific_name)

par(mfrow=c(2,1))
Mammals = merge(Agg1,GL, by = 'Species')
cor.test(Mammals$FrOfGainers,Mammals$GenerationLength_d, method = 'spearman') #  "ALINA PICs" (MIN 2)
plot(log2(Mammals$GenerationLength_d),Mammals$FrOfGainers) # 0.4297347

cor.test(Mammals$FrOfLoosers,Mammals$GenerationLength_d, method = 'spearman')  # "ALINA PICs"
plot(log2(Mammals$GenerationLength_d),Mammals$FrOfLoosers) # -0.3426591 

###############################
######### fourth quartile (works good)
###############################
data$Gainers = data$Asn + data$Lys
data$Loosers = data$Gly
data$All = data$Asn + data$Lys + data$Ser + data$Asp + data$Glu
data$Intermediate = data$All - data$Gainers - data$Loosers # 64 - 16 - 14 - 4 (Stops) = 30

#### 12 genes
data1 = data[data$Gene != 'ND6',]
Agg1 = aggregate(list(data1$Gainers,data1$Loosers,data1$Intermediate,data1$All), by = list(data1$Species,data1$Class), FUN = sum)
names(Agg1) = c('Species','Class','Gainers','Loosers','Intermediate','All')
Agg1$GainersToLosers = Agg1$Gainers - Agg1$Loosers
Agg1$FrOfGainers = Agg1$Gainers/Agg1$All
Agg1$FrOfLoosers = Agg1$Loosers/Agg1$All

GL = read.table("../../Body/1Raw/GenerationLenghtforMammals.xlsx.txt", header = TRUE, sep = '\t')
GL$Species = gsub(' ','_',GL$Scientific_name)

par(mfrow=c(2,1))
Mammals = merge(Agg1,GL, by = 'Species')
cor.test(Mammals$FrOfGainers,Mammals$GenerationLength_d, method = 'spearman') #  "ALINA PICs"
plot(log2(Mammals$GenerationLength_d),Mammals$FrOfGainers) # 0.2972999

cor.test(Mammals$FrOfLoosers,Mammals$GenerationLength_d, method = 'spearman')  # "ALINA PICs"
plot(log2(Mammals$GenerationLength_d),Mammals$FrOfLoosers) # -0.3426591 - a bit positive

################### 6. gainers / loosers in second (only Ah>Gh) and fourth (only Ch>Th) quartiles in the relationship with temperature of fishes (it seems that fourth quartile is stronger)
###############################
######### second quartile (works good)
###############################
data$Gainers = data$Pro
data$Loosers = data$Phe
data$All = data$Pro + data$Phe + data$Leu + data$Ser
data$Intermediate = data$All - data$Gainers - data$Loosers # 64 - 16 - 14 - 4 (Stops) = 30

#### 12 genes
data1 = data[data$Gene != 'ND6',]
Agg1 = aggregate(list(data1$Gainers,data1$Loosers,data1$Intermediate,data1$All), by = list(data1$Species,data1$Class), FUN = sum)
names(Agg1) = c('Species','Class','Gainers','Loosers','Intermediate','All')
Agg1$GainersToLosers = Agg1$Gainers - Agg1$Loosers
Agg1$FrOfGainers = Agg1$Gainers/Agg1$All
Agg1$FrOfLoosers = Agg1$Loosers/Agg1$All

Agg1 = Agg1[Agg1$Class == 'Actinopterygii',]
Temp = read.table("../../Body/1Raw/FishBaseTemperature.txt", header = TRUE)
Temp = aggregate(Temp$Temperature, by = list(Temp$Species), FUN = median)
names(Temp) = c('Species','Temperature')
nrow(Temp)

Fish = merge(Agg1,Temp, by = 'Species') # 305
nrow(Fish)

par(mfrow=c(2,1))
cor.test(Fish$FrOfGainers,Fish$Temperature, method = 'spearman') # "ALINA PICs"
plot(Fish$Temperature, Fish$FrOfGainers) # ALINA PICs

cor.test(Fish$FrOfLoosers,Fish$Temperature, method = 'spearman') # "ALINA PICs"
plot(Fish$Temperature, Fish$FrOfLoosers) # ALINA PICs

###############################
######### fourth quartile (works good)
###############################
data$Gainers = data$Asn + data$Lys
data$Loosers = data$Gly
data$All = data$Asn + data$Lys + data$Ser + data$Asp + data$Glu
data$Intermediate = data$All - data$Gainers - data$Loosers # 64 - 16 - 14 - 4 (Stops) = 30

#### 12 genes
data1 = data[data$Gene != 'ND6',]
Agg1 = aggregate(list(data1$Gainers,data1$Loosers,data1$Intermediate,data1$All), by = list(data1$Species,data1$Class), FUN = sum)
names(Agg1) = c('Species','Class','Gainers','Loosers','Intermediate','All')
Agg1$GainersToLosers = Agg1$Gainers - Agg1$Loosers
Agg1$FrOfGainers = Agg1$Gainers/Agg1$All
Agg1$FrOfLoosers = Agg1$Loosers/Agg1$All

Agg1 = Agg1[Agg1$Class == 'Actinopterygii',]
Temp = read.table("../../Body/1Raw/FishBaseTemperature.txt", header = TRUE)
Temp = aggregate(Temp$Temperature, by = list(Temp$Species), FUN = median)
names(Temp) = c('Species','Temperature')
nrow(Temp)

Fish = merge(Agg1,Temp, by = 'Species') # 305
nrow(Fish)

par(mfrow=c(2,1))
cor.test(Fish$FrOfGainers,Fish$Temperature, method = 'spearman') # "ALINA PICs"
plot(Fish$Temperature, Fish$FrOfGainers) # ALINA PICs

cor.test(Fish$FrOfLoosers,Fish$Temperature, method = 'spearman') # "ALINA PICs"
plot(Fish$Temperature, Fish$FrOfLoosers) # ALINA PICs

dev.off()

################### 4 analysis of all paired trajectories within each quartet:
# first quartet:
#Agg = aggregate(list(data$Cys,data$Trp,data$Tyr,data$His,data$Gln,data$Arg), by = list(data$Species,data$Class), FUN = sum)
#names(Agg) = c('Species','Class','Cys','Trp','Tyr','His','Gln','Arg')
#boxplot(Agg$Cys,Agg$Tyr,Agg$His,Agg$Cys,Agg$Trp,Agg$Arg,Agg$His,Agg$Gln, names = c('Cys','Tyr','His','Cys','Trp','Arg','His','Gln'), outline = FALSE, notch = TRUE, main = 'first quartet')

# second quartet:
#Agg = aggregate(list(data$Phe,data$Leu,data$Ser,data$Pro), by = list(data$Species,data$Class), FUN = sum)
#names(Agg) = c('Species','Class','Phe','Leu','Ser','Pro')
#boxplot(Agg$Phe,Agg$Leu,Agg$Pro,Agg$Phe,Agg$Ser,Agg$Pro, names = c('Phe','Leu','Pro','Phe','Ser','Pro'), outline = FALSE, notch = TRUE, main = 'second quartet')

# third quartet:
#Agg = aggregate(list(data$Val,data$Ala,data$Thr,data$Met,data$Ile), by = list(data$Species,data$Class), FUN = sum)
#names(Agg) = c('Species','Class','Val','Ala','Thr','Met','Ile')
#boxplot(Agg$Val,Agg$Met,Agg$Ile,Agg$Thr,Agg$Val,Agg$Ala,Agg$Thr, names = c('Val','Met','Ile','Thr','Val','Ala','Thr'), outline = FALSE, notch = TRUE, main = 'third quartet')
# dev.off()
