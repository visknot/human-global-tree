rm(list=ls(all=TRUE)) 

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
