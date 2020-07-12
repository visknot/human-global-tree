# compare Germ vs Soma MutSpec

rm(list=ls(all=TRUE))
library('dplyr')

pdf("../../Body/4Figures/05.MutSpecGermVsSoma.R.pdf", width = 22, height = 14)
par(mfrow=c(3,1))

######### 1 COMPARE GERM MUTSPEC12 INTERNAL VERSUS EXTERNAL
Germ = read.table("../../Body/3Results/03A.MutSpecGermObsToExp.R.MutSpec12.txt", header = TRUE)
GermInt = Germ[Germ$Branches == 'InternalBranches',]
GermExt = Germ[Germ$Branches == 'TerminalBranches',]
GermAll = Germ[Germ$Branches == 'AllBranches',]

## contrasts of two MutSpecs:
GermInt$ObsToExpFreqNorm = GermInt$ObsToExpFreq/sum(GermInt$ObsToExpFreq); GermInt = GermInt[order(GermInt$SubstHeavy),]
GermExt$ObsToExpFreqNorm = GermExt$ObsToExpFreq/sum(GermExt$ObsToExpFreq); GermExt = GermExt[order(GermExt$SubstHeavy),]

Germ.ExtToInt.MutSpec12.Heavy = (GermExt$ObsToExpFreqNorm - GermInt$ObsToExpFreqNorm) * ((GermExt$ObsToExpFreqNorm + GermInt$ObsToExpFreqNorm) / 2)
ColVec=c('gray',rgb(0,1,0,1),'gray','gray','gray',rgb(1,0,0,1),rgb(1,0,0,0.2),'gray','gray','gray',rgb(0,1,0,0.2),'gray'); length(ColVec)

barplot(GermExt$ObsToExpFreqNorm, names = GermExt$SubstHeavy, col = ColVec, main = paste("GERM ObsToExpNormalized, internal, Heavy Chain"), cex.names = 2, cex.axis = 2)
barplot(GermInt$ObsToExpFreqNorm, names = GermInt$SubstHeavy, col = ColVec, main = paste("GERM ObsToExpNormalized, external, Heavy Chain"), cex.names = 2, cex.axis = 2)
barplot(Germ.ExtToInt.MutSpec12.Heavy, names = GermInt$SubstHeavy, main = paste("GERM ObsToExpNormalized, internal, Heavy Chain"), cex.names = 2, cex.axis = 2)

# plot bar to bar - close to each other in order to compare by eye - first is external!!!!!!!!
par(mfrow=c(1,3))

GroupedBarplot = data.frame(c(GermExt$ObsToExpFreqNorm,GermInt$ObsToExpFreqNorm),c(GermExt$SubstHeavy,GermInt$SubstHeavy),c(GermExt$Branches,GermInt$Branches))
names(GroupedBarplot)=c('ObsToExpFreqNorm','SubstHeavy','Branches')
row.names(GroupedBarplot) = paste(GroupedBarplot$SubstHeavy,GroupedBarplot$Branches, sep = '_')
GroupedBarplot = GroupedBarplot[order(GroupedBarplot$SubstHeavy,-as.numeric(as.factor(GroupedBarplot$Branches))),] # -as.numeric(as.factor(

ColVec=c('gray','gray',rgb(0,1,0,1),rgb(0,1,0,1),'gray','gray','gray','gray','gray','gray',rgb(1,0,0,1),rgb(1,0,0,1),rgb(1,0,0,0.2),rgb(1,0,0,0.2),'gray','gray','gray','gray','gray','gray',rgb(0,1,0,0.2),rgb(0,1,0,0.2),'gray','gray');
BorderVec = rep(c(rgb(0,0,0,1),rgb(1,1,1,0)),12) 
barplot(GroupedBarplot$ObsToExpFreqNorm, names = GroupedBarplot$SubstHeavy, main = paste("Int vs Ext GERM MutSpecs, Heavy Chain"), col = ColVec, border = BorderVec, cex.names = 2, cex.axis = 2, las = 2)

##### check asymmetry between Ts and fraction of Tv
TsVec = c('AG','TC','CT','GA')
### external branches - a bit higher asymmetry between key transitions - chemical damage!
GermExt[GermExt$SubstHeavy == 'CT',]$ObsToExpFreqNorm/GermExt[GermExt$SubstHeavy == 'GA',]$ObsToExpFreqNorm # 13.52577
GermInt[GermInt$SubstHeavy == 'CT',]$ObsToExpFreqNorm/GermInt[GermInt$SubstHeavy == 'GA',]$ObsToExpFreqNorm # 12.67414

GermExt[GermExt$SubstHeavy == 'AG',]$ObsToExpFreqNorm/GermExt[GermExt$SubstHeavy == 'TC',]$ObsToExpFreqNorm # 1.962073
GermInt[GermInt$SubstHeavy == 'AG',]$ObsToExpFreqNorm/GermInt[GermInt$SubstHeavy == 'TC',]$ObsToExpFreqNorm # 1.833251

### external branches - a bit less Tv (more Ts - which might be chemical induced)
AllTvExt = sum(GermExt[!GermExt$SubstHeavy %in% TsVec,]$ObsToExpFreqNorm) # 0.03020045
AllTvInt = sum(GermInt[!GermInt$SubstHeavy %in% TsVec,]$ObsToExpFreqNorm) # 0.03295581
barplot(c(AllTvExt,AllTvInt), names = c('TvExt','TvInt'), main = "GERM: Fraction of Tv in Ext vs Int Branches", cex.names = 2, cex.axis = 2, las = 2)
plot.new()

######### 2 COMPARE GERM MUTSPEC12 ALLBRANCHES WITH CANCER
par(mfrow=c(3,1))

Soma = read.table("../../Body/3Results/04A.MutSpecCancerObsToExp.MutSpec12.txt", header = TRUE)
SomaAllExceptStop = Soma[Soma$SubstType == 'AllExceptStopGains',]
SomaAllExceptStop$ObsToExpFreqNorm = SomaAllExceptStop$ObsToExpFreq/sum(SomaAllExceptStop$ObsToExpFreq); 
SomaAllExceptStop = SomaAllExceptStop[order(SomaAllExceptStop$SubstHeavy),]

GermAll$ObsToExpFreqNorm = GermAll$ObsToExpFreq/sum(GermAll$ObsToExpFreq); 
GermAll = GermAll[order(GermAll$SubstHeavy),]

ColVec=c('gray',rgb(0,1,0,1),'gray','gray','gray',rgb(1,0,0,1),rgb(1,0,0,0.2),'gray','gray','gray',rgb(0,1,0,0.2),'gray'); length(ColVec)

barplot(SomaAllExceptStop$ObsToExpFreqNorm, names = SomaAllExceptStop$SubstHeavy, col = ColVec, main = paste("SOMA ObsToExpNormalized,Heavy Chain"), cex.names = 2, cex.axis = 2)
barplot(GermAll$ObsToExpFreqNorm, names = GermAll$SubstHeavy, col = ColVec, main = paste("GERM ObsToExpNormalized, AllBranches, Heavy Chain"), cex.names = 2, cex.axis = 2)
plot.new()

##### check asymmetry between Ts and fraction of Tv
TsVec = c('AG','TC','CT','GA')
### external branches - a bit higher asymmetry between key transitions - chemical damage!
SomaAllExceptStop[SomaAllExceptStop$SubstHeavy == 'CT',]$ObsToExpFreqNorm/SomaAllExceptStop[SomaAllExceptStop$SubstHeavy == 'GA',]$ObsToExpFreqNorm # 14.30442
GermAll[GermAll$SubstHeavy == 'CT',]$ObsToExpFreqNorm/GermAll[GermAll$SubstHeavy == 'GA',]$ObsToExpFreqNorm # 13.18345

SomaAllExceptStop[SomaAllExceptStop$SubstHeavy == 'AG',]$ObsToExpFreqNorm/SomaAllExceptStop[SomaAllExceptStop$SubstHeavy == 'TC',]$ObsToExpFreqNorm # 5.905732
GermAll[GermAll$SubstHeavy == 'AG',]$ObsToExpFreqNorm/GermAll[GermAll$SubstHeavy == 'TC',]$ObsToExpFreqNorm # 1.910779

### external branches - a bit less Tv (more Ts - which might be chemical induced)
AllTvExt = sum(SomaAllExceptStop[!SomaAllExceptStop$SubstHeavy %in% TsVec,]$ObsToExpFreqNorm) # 0.06013331
AllTvInt = sum(GermAll[!GermAll$SubstHeavy %in% TsVec,]$ObsToExpFreqNorm) # 0.03127012

# plot bar to bar - close to each other in order to compare by eye - first is external!!!!!!!!
par(mfrow=c(1,3))

GroupedBarplot = data.frame(c(GermAll$ObsToExpFreqNorm,SomaAllExceptStop$ObsToExpFreqNorm),c(GermAll$SubstHeavy,SomaAllExceptStop$SubstHeavy),c(rep('GermAll',12),rep('SomaAllExceptStop',12)))
names(GroupedBarplot)=c('ObsToExpFreqNorm','SubstHeavy','GermSoma')
GroupedBarplot = GroupedBarplot[order(GroupedBarplot$SubstHeavy,-as.numeric(as.factor(GroupedBarplot$GermSoma))),] # -as.numeric(as.factor(

ColVec=c('gray','gray',rgb(0,1,0,1),rgb(0,1,0,1),'gray','gray','gray','gray','gray','gray',rgb(1,0,0,1),rgb(1,0,0,1),rgb(1,0,0,0.2),rgb(1,0,0,0.2),'gray','gray','gray','gray','gray','gray',rgb(0,1,0,0.2),rgb(0,1,0,0.2),'gray','gray');
BorderVec = rep(c(rgb(0,0,0,1),rgb(1,1,1,0)),12) 
barplot(GroupedBarplot$ObsToExpFreqNorm, names = GroupedBarplot$SubstHeavy, main = paste("SOMA vs GERM MutSpecs, Heavy Chain"), col = ColVec, border = BorderVec, cex.names = 2, cex.axis = 2, las = 2)
plot.new()
plot.new()

###### 3 GERM MUTSPEC 192 INTERNAL VS EXTERNAL
par(mfrow=c(3,1))
Germ = read.table("../../Body/3Results/03A.MutSpecGermObsToExp.R.MutSpec192.txt", header = TRUE)
GermInt = Germ[Germ$Branches == 'InternalBranches',]
GermExt = Germ[Germ$Branches == 'TerminalBranches',]
GermAll = Germ[Germ$Branches == 'AllBranches',]

GermInt$ObsToExpFreqNorm = GermInt$ObsToExpFreq/sum(GermInt$ObsToExpFreq); GermInt = GermInt[order(GermInt$NamesHeavy),]
GermExt$ObsToExpFreqNorm = GermExt$ObsToExpFreq/sum(GermExt$ObsToExpFreq); GermExt = GermExt[order(GermExt$NamesHeavy),]

ColVec=c(rep('gray',16),rep(rgb(0,1,0,1),16),rep('gray',16),rep('gray',16),rep('gray',16),rep(rgb(1,0,0,1),16),rep(rgb(1,0,0,0.2),16),rep('gray',16),rep('gray',16),rep('gray',16),rep(rgb(0,1,0,0.2),16),rep('gray',16)); length(ColVec)
BorderVec = rep(c(rgb(0.1,0.1,0.1,0.1),rgb(0,0,0,1),rgb(0.1,0.1,0.1,0.1),rgb(0.1,0.1,0.1,0.1)),48) # 192/4 The second colum is with border

barplot(GermExt$ObsToExpFreqNorm, names = GermExt$NamesHeavy, col = ColVec, border = BorderVec, main = paste("GERM ObsToExpNormalized, external, Heavy Chain"), cex.names = 1, cex.axis = 1, las = 2, ylim = c(0,0.08))
barplot(GermInt$ObsToExpFreqNorm, names = GermInt$NamesHeavy, col = ColVec, border = BorderVec, main = paste("GERM ObsToExpNormalized, internal, Heavy Chain"), cex.names = 1, cex.axis = 1, las = 2, ylim = c(0,0.08))
plot.new()

######### 4 MUTSPEC192 GERMALL VS SOMA

GermAll$ObsToExpFreqNorm = GermAll$ObsToExpFreq/sum(GermAll$ObsToExpFreq); GermAll = GermAll[order(GermInt$NamesHeavy),]
GermAll = GermAll[order(GermAll$NamesHeavy),]

Soma = read.table("../../Body/3Results/04A.MutSpecCancerObsToExp.MutSpec192.txt", header = TRUE)
SomaAllExceptStop = Soma[Soma$SubstType == 'AllExceptStopGains',]
SomaAllExceptStop$ObsToExpFreqNorm = SomaAllExceptStop$ObsToExpFreq/sum(SomaAllExceptStop$ObsToExpFreq); 
SomaAllExceptStop = SomaAllExceptStop[order(SomaAllExceptStop$NamesHeavy),]

ColVec=c(rep('gray',16),rep(rgb(0,1,0,1),16),rep('gray',16),rep('gray',16),rep('gray',16),rep(rgb(1,0,0,1),16),rep(rgb(1,0,0,0.2),16),rep('gray',16),rep('gray',16),rep('gray',16),rep(rgb(0,1,0,0.2),16),rep('gray',16)); length(ColVec)
BorderVec = rep(c(rgb(0.1,0.1,0.1,0.1),rgb(0,0,0,1),rgb(0.1,0.1,0.1,0.1),rgb(0.1,0.1,0.1,0.1)),48) # 192/4 The second colum is with border

barplot(GermAll$ObsToExpFreqNorm, names = GermAll$NamesHeavy, col = ColVec, border = BorderVec, main = paste("GERM All ObsToExpNormalized, Heavy Chain"), cex.names = 1, cex.axis = 1, las = 2, ylim = c(0,0.08))
barplot(SomaAllExceptStop$ObsToExpFreqNorm, names = SomaAllExceptStop$NamesHeavy, col = ColVec, border = BorderVec, main = paste("SOMA ObsToExpNormalized, Heavy Chain"), cex.names = 1, cex.axis = 1, las = 2, ylim = c(0,0.08))
plot.new()

#### plot only Ts
par(mfrow=c(2,1))
GermAllTs = GermAll[grepl('CT:',GermAll$NamesHeavy) | grepl('GA:',GermAll$NamesHeavy) | grepl('TC:',GermAll$NamesHeavy) | grepl('AG:',GermAll$NamesHeavy),]
SomaTs = SomaAllExceptStop[grepl('CT:',SomaAllExceptStop$NamesHeavy) | grepl('GA:',SomaAllExceptStop$NamesHeavy) | grepl('TC:',SomaAllExceptStop$NamesHeavy) | grepl('AG:',SomaAllExceptStop$NamesHeavy),]
ColVec=c(rep(rgb(0,1,0,1),16),rep(rgb(1,0,0,1),16),rep(rgb(1,0,0,0.2),16),rep(rgb(0,1,0,0.2),16)); length(ColVec)
BorderVec = rep(c(rgb(0.1,0.1,0.1,0.1),rgb(0,0,0,1),rgb(0.1,0.1,0.1,0.1),rgb(0.1,0.1,0.1,0.1)),16) # 192/4 The second colum is with border

barplot(GermAllTs$ObsToExpFreqNorm, names = GermAllTs$NamesHeavy, col = ColVec, border = BorderVec, main = paste("GERM All ObsToExpNormalized, Heavy Chain"), cex.names = 1, cex.axis = 1, las = 2, ylim = c(0,0.08))
barplot(SomaTs$ObsToExpFreqNorm, names = SomaTs$NamesHeavy, col = ColVec, border = BorderVec, main = paste("SOMA ObsToExpNormalized, Heavy Chain"), cex.names = 1, cex.axis = 1, las = 2, ylim = c(0,0.08))
plot.new()

dev.off()

