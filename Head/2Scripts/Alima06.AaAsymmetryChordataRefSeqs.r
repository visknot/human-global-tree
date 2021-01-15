rm(list=ls(all=TRUE)) 

pdf("../../Body/4Figures/Alima06.AaAsymmetryChordataRefSeqs.r.pdf")
data = read.table("../../Body/1Raw/AminoAcidFreqsChordata.txt", sep = ',', header = TRUE, quote = '')
colnames(data)
data = data[,-c(23:33)]
data$Gainers = data$Pro + data$Thr + data$His + data$Gln + data$Asn + data$Lys # 16 codons 
data$Loosers = data$Phe + data$Val + data$Gly + data$Cys  # 12 codons
data$All = apply(as.matrix(data[,c(3:22)]),1,FUN = sum)
data$Intermediate = data$All - data$Gainers - data$Loosers # 64 - 16 - 12 - 4 (Stops) = 32
Agg1 = aggregate(list(data$Gainers,data$Loosers,data$Intermediate,data$All), by = list(data$Species,data$Class), FUN = sum)
names(Agg1) = c('Species','Class','Gainers','Loosers','Intermediate','All')
Agg1$GainersToLosers = Agg1$Gainers - Agg1$Loosers
Agg1$FrOfGainers = Agg1$Gainers/Agg1$All
Agg1$FrOfLoosers = Agg1$Loosers/Agg1$All

## 1 gainers to loosers between classes
boxplot(Agg1$Gainers,Agg1$Loosers,Agg1$Intermediate) 
boxplot(Agg1$Gainers/16,Agg1$Loosers/12,Agg1$Intermediate/32) 
table(Agg1$Class)
Cold = c('Actinopterygii','Amphibia','Reptilia')
Warm = c('Aves','Mammalia')

boxplot(Agg1$Gainers,Agg1$Loosers)
par(mfrow=c(1,2))
boxplot(Agg1[Agg1$Class %in% Cold,]$Gainers,Agg1[Agg1$Class %in% Cold,]$Loosers, ylim = c(500,1200), names = c('gainers','loosers'), main = 'Coldblood', ylab = 'NumberOfAaPer12Genes')
boxplot(Agg1[Agg1$Class %in% Warm,]$Gainers,Agg1[Agg1$Class %in% Warm,]$Loosers, ylim = c(500,1200), names = c('gainers','loosers'), main = 'Warmblood')

par(mfrow=c(1,3))
boxplot(Agg1[Agg1$Class == 'Actinopterygii',]$Gainers,Agg1[Agg1$Class == 'Actinopterygii',]$Loosers, ylim = c(500,1200), names = c('gainers','loosers'), main = 'Actinopterygii', ylab = 'NumberOfAaPer12Genes')
boxplot(Agg1[Agg1$Class == 'Amphibia',]$Gainers,Agg1[Agg1$Class == 'Amphibia',]$Loosers, ylim = c(500,1200), names = c('gainers','loosers'), main = 'Amphibia')
boxplot(Agg1[Agg1$Class == 'Reptilia',]$Gainers,Agg1[Agg1$Class == 'Reptilia',]$Loosers, ylim = c(500,1200), names = c('gainers','loosers'), main = 'Reptilia')

par(mfrow=c(1,2))
boxplot(Agg1[Agg1$Class == 'Mammalia',]$Gainers,Agg1[Agg1$Class == 'Mammalia',]$Loosers, ylim = c(500,1200), names = c('gainers','loosers'), main = 'Mammalia')
boxplot(Agg1[Agg1$Class == 'Aves',]$Gainers,Agg1[Agg1$Class == 'Aves',]$Loosers, ylim = c(500,1200), names = c('gainers','loosers'), main = 'Aves')

par(mfrow=c(1,1))
wilcox.test(Agg1[Agg1$Class == 'Mammalia',]$FrOfGainers,Agg1[Agg1$Class == 'Aves',]$FrOfGainers)
hist(Agg1[Agg1$Class == 'Mammalia',]$FrOfGainers, xlim = c(0.2,0.35), ylim = c(0,250), col = rgb(0.2,0.2,0.2,0.5), main = 'fraction of gainers in warmblooded', xlab = '') # dev.off()
par(new = TRUE)
hist(Agg1[Agg1$Class == 'Aves',]$FrOfGainers, xlim = c(0.2,0.35), ylim = c(0,250), col =  rgb(1,0.1,0.1,0.5), main = '', xlab = '') # dev.off()

hist(Agg1[Agg1$Class == 'Actinopterygii',]$FrOfGainers, xlim = c(0.2,0.35), ylim = c(0,250), col = rgb(0.1,0.1,1,0.5), main = 'fraction of gainers in coldblooded', xlab = '')
par(new = TRUE)
hist(Agg1[Agg1$Class == 'Amphibia',]$FrOfGainers, xlim = c(0.2,0.35), ylim = c(0,250), col =  rgb(0.1,1,0.1,0.5), main = '', xlab = '')
par(new = TRUE)
hist(Agg1[Agg1$Class == 'Reptilia',]$FrOfGainers, xlim = c(0.2,0.35), ylim = c(0,250), col =  rgb(1,0.1,0.1,0.5), main = '', xlab = '')


## 2 gainers to loosers is higher in longlived mammals (aged oocytes with strong asymmetric damage!)

GL = read.table("../../Body/1Raw/GenerationLenghtforMammals.xlsx.txt", header = TRUE, sep = '\t')
GL$Species = gsub(' ','_',GL$Scientific_name)

par(mfrow=c(1,1))
Mammals = merge(Agg1,GL, by = 'Species')
cor.test(Mammals$FrOfGainers,Mammals$GenerationLength_d, method = 'spearman')
plot(log2(Mammals$GenerationLength_d),Mammals$FrOfGainers) # dev.off()

cor.test(Mammals$GainersToLosers,Mammals$GenerationLength_d, method = 'spearman')
plot(log2(Mammals$GenerationLength_d),Mammals$GainersToLosers)

## 3 analysis of all paired trajectories within each quartet:

# first quartet:
Agg = aggregate(list(data$Cys,data$Trp,data$Tyr,data$His,data$Gln,data$Arg), by = list(data$Species,data$Class), FUN = sum)
names(Agg) = c('Species','Class','Cys','Trp','Tyr','His','Gln','Arg')
boxplot(Agg$Cys,Agg$Tyr,Agg$His,Agg$Cys,Agg$Trp,Agg$Arg,Agg$His,Agg$Gln, names = c('Cys','Tyr','His','Cys','Trp','Arg','His','Gln'), outline = FALSE, notch = TRUE, main = 'first quartet')

# second quartet:
Agg = aggregate(list(data$Phe,data$Leu,data$Ser,data$Pro), by = list(data$Species,data$Class), FUN = sum)
names(Agg) = c('Species','Class','Phe','Leu','Ser','Pro')
boxplot(Agg$Phe,Agg$Leu,Agg$Pro,Agg$Phe,Agg$Ser,Agg$Pro, names = c('Phe','Leu','Pro','Phe','Ser','Pro'), outline = FALSE, notch = TRUE, main = 'second quartet')

# third quartet:
Agg = aggregate(list(data$Val,data$Ala,data$Thr,data$Met,data$Ile), by = list(data$Species,data$Class), FUN = sum)
names(Agg) = c('Species','Class','Val','Ala','Thr','Met','Ile')
boxplot(Agg$Val,Agg$Met,Agg$Ile,Agg$Thr,Agg$Val,Agg$Ala,Agg$Thr, names = c('Val','Met','Ile','Thr','Val','Ala','Thr'), outline = FALSE, notch = TRUE, main = 'third quartet')


dev.off()
