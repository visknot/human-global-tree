rm(list=ls(all=TRUE))
library('dplyr')

pdf("../../Body/4Figures/05.MutSpecGermVsSoma.R.pdf", width = 22, height = 14)
par(mfrow=c(3,1))

##### 12:
Soma = read.table("../../Body/3Results/MutSpec12ObsToExpCancerNeutralClassesOfMut.txt", header = TRUE)
Soma = select(Soma,Subst,ObsToExpFreq); names(Soma)[2] = c('ObsToExpFreqSoma')
Soma$ObsToExpFreqSoma = Soma$ObsToExpFreqSoma/sum(Soma$ObsToExpFreqSoma)
Germ = read.table("../../Body/3Results/MutSpec12ObsToExpSynNoNd6.txt", header = TRUE)
Germ = select(Germ,Subst,ObsToExpFreq); names(Germ)[2] = c('ObsToExpFreqGerm')
Germ$ObsToExpFreqGerm = Germ$ObsToExpFreqGerm/sum(Germ$ObsToExpFreqGerm)
GermSoma = merge(Germ,Soma)
GermSoma$SomaToGermPercent = (GermSoma$ObsToExpFreqSoma - GermSoma$ObsToExpFreqGerm) / (GermSoma$ObsToExpFreqGerm)
ColVec=c('gray',rgb(0,1,0,0.2),'gray','gray','gray',rgb(1,0,0,0.2),rgb(1,0,0,1),'gray','gray','gray',rgb(0,1,0,1),'gray'); length(ColVec)

barplot(GermSoma$ObsToExpFreqGerm, names = GermSoma$Subst, col = ColVec, main = paste("MutSpec in germa"), cex.names = 2, cex.axis = 2)
barplot(GermSoma$ObsToExpFreqSoma, names = GermSoma$Subst, col = ColVec, main = paste("MutSpec in soma"), cex.names = 2, cex.axis = 2)
barplot(GermSoma$SomaToGermPercent, names = GermSoma$Subst, col = ColVec, main = paste("Differences: soma-germ/germ"), cex.names = 2, cex.axis = 2)

## asymmetry in transitions:
names(GermSoma)
GermSoma[GermSoma$Subst == 'GA',]$ObsToExpFreqGerm / GermSoma[GermSoma$Subst == 'CT',]$ObsToExpFreqGerm  # 10.485  
GermSoma[GermSoma$Subst == 'GA',]$ObsToExpFreqSoma / GermSoma[GermSoma$Subst == 'CT',]$ObsToExpFreqSoma  # 11.267

GermSoma[GermSoma$Subst == 'TC',]$ObsToExpFreqGerm / GermSoma[GermSoma$Subst == 'AG',]$ObsToExpFreqGerm  # 1.801  
GermSoma[GermSoma$Subst == 'TC',]$ObsToExpFreqSoma / GermSoma[GermSoma$Subst == 'AG',]$ObsToExpFreqSoma  # 5.602

##### 192:
Soma = read.table("../../Body/3Results/MutSpec192ObsToExpCancerNeutralClassesOfMut.txt", header = TRUE, sep = '\t')
Soma = select(Soma,names,ObsToExpFreq); names(Soma)[2] = c('ObsToExpFreqSoma')
Soma$ObsToExpFreqSoma = Soma$ObsToExpFreqSoma/sum(Soma$ObsToExpFreqSoma)
Germ = read.table("../../Body/3Results/MutSpec192ObsToExpSynNoNd6.txt", header = TRUE, sep = '\t')
Germ = select(Germ,names,ObsToExpFreq); names(Germ)[2] = c('ObsToExpFreqGerm')
Germ[is.na(Germ)] <- 0
Germ$ObsToExpFreqGerm = Germ$ObsToExpFreqGerm/sum(Germ$ObsToExpFreqGerm)
GermSoma = merge(Germ,Soma, by = 'names')
GermSoma$SomaToGermPercent = (GermSoma$ObsToExpFreqSoma - GermSoma$ObsToExpFreqGerm) / (GermSoma$ObsToExpFreqGerm)
GermSoma$SomaToGermPercent[is.na(GermSoma$SomaToGermPercent)] <- 0
GermSoma$SomaToGermPercent[GermSoma$SomaToGermPercent == Inf] <- 0

ColVec=c(rep('gray',16),rep(rgb(0,1,0,0.2),16),rep('gray',16),rep('gray',16),rep('gray',16),rep(rgb(1,0,0,0.2),16),rep(rgb(1,0,0,1),16),rep('gray',16),rep('gray',16),rep('gray',16),rep(rgb(0,1,0,1),16),rep('gray',16)); length(ColVec)
barplot(GermSoma$ObsToExpFreqGerm, names = GermSoma$names, col = ColVec, main = paste("MutSpec in germa"), cex.names = 1, cex.axis = 2, las = 2)
barplot(GermSoma$ObsToExpFreqSoma, names = GermSoma$names, col = ColVec, main = paste("MutSpec in soma"), cex.names = 1, cex.axis = 2, las = 2)
barplot(GermSoma$SomaToGermPercent, names = GermSoma$names, col = ColVec, main = paste("Differences: soma-germ/germ"), cex.names = 1, cex.axis = 2, las = 2)
dev.off()