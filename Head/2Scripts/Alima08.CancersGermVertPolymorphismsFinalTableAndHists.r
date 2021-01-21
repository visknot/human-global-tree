###################################
##### prepare final tables, stats and hists from all analyses
###################################

rm(list=ls(all=TRUE)) 

pdf("../../Body/4Figures/Alima08.CancersGermVertPolymorphismsFinalTableAndHists.r.pdf")

SuperFinal = data.frame()

####################
## CANCERS
####################  

AllExceptNd6 = read.table("../../Body/3Results/Alima04.AsymmetryInCancers.12genes.txt", header = TRUE)
names(AllExceptNd6)
AllExceptNd6 = AllExceptNd6[,c(2:8)]
AllExceptNd6 = AllExceptNd6[]
Nd6 = read.table("../../Body/3Results/Alima04.AsymmetryInCancers.ND6.txt", header = TRUE)
names(Nd6)
Nd6 = Nd6[,c(4,6,7,8)]
names(Nd6)[2:4]=c('ExpectedMoreThanOne.Nd6','NumberOfExpectedAaSubst.Nd6','NumberOfUnexpectedAaSubst.Nd6')
All = merge(AllExceptNd6,Nd6, all = TRUE)
All[is.na(All$NumberOfExpectedAaSubst.Nd6),]$NumberOfExpectedAaSubst.Nd6<-0
All[is.na(All$NumberOfUnexpectedAaSubst.Nd6),]$NumberOfUnexpectedAaSubst.Nd6<-0
# All = All[order(-All$L),] - к черту
names(All)
Final = All[,c(1,6,7,5,9,10,8)]
Final$ExpectedMoreThanOne = round(log2(Final$ExpectedMoreThanOne),2); summary(Final$ExpectedMoreThanOne)
Final$ExpectedMoreThanOne.Nd6 = round(log2(Final$ExpectedMoreThanOne.Nd6),2); summary(Final$ExpectedMoreThanOne.Nd6)

Final = Final[order(Final$ExpectedMoreThanOne),]

hist(Final$ExpectedMoreThanOne, xlim = c(-4,8), ylim = c(0,9), col = rgb(1,0.1,0.1,0.3), xlab = '', main = ''); par(new = TRUE)
hist(Final$ExpectedMoreThanOne.Nd6, xlim = c(-4,8),ylim = c(0,9), col = rgb(0.1,0.1,0.1,0.3), xlab = 'log2(expected/unexpected)', main = 'cancer');
abline(v=0, lt = 2, lwd = 2, col = 'black')

wilcox.test(Final$ExpectedMoreThanOne, mu = 0)     # 4.939e-05
wilcox.test(Final$ExpectedMoreThanOne.Nd6, mu = 0) # 0.02225

Final$DataSet = 'Cancer'
SuperFinal=rbind(SuperFinal,Final)

#########
## GERM - 12 genes and Nd6 - all positions
########

AllExceptNd6 = read.table("../../Body/3Results/Alima02.AaAsymmetry.HumanTree.12Genes.ExpectedVsObservedAaChanges.txt", header = TRUE)
names(AllExceptNd6)
AllExceptNd6 = AllExceptNd6[,c(2:8)]
# AllExceptNd6 = AllExceptNd6[]
Nd6 = read.table("../../Body/3Results/Alima02.AaAsymmetry.HumanTree.ND6.ExpectedVsObservedAaChanges.txt", header = TRUE)
names(Nd6)
Nd6 = Nd6[,c(4,6,7,8)]
names(Nd6)[2:4]=c('ExpectedMoreThanOne.Nd6','NumberOfExpectedAaSubst.Nd6','NumberOfUnexpectedAaSubst.Nd6')
All = merge(AllExceptNd6,Nd6, all = TRUE)
All[is.na(All$NumberOfExpectedAaSubst.Nd6),]$NumberOfExpectedAaSubst.Nd6<-0
All[is.na(All$NumberOfUnexpectedAaSubst.Nd6),]$NumberOfUnexpectedAaSubst.Nd6<-0
# All = All[order(-All$L),] - к черту
names(All)
Final = All[,c(1,6,7,5,9,10,8)]
Final$ExpectedMoreThanOne = round(log2(Final$ExpectedMoreThanOne),2); summary(Final$ExpectedMoreThanOne)
Final$ExpectedMoreThanOne.Nd6 = round(log2(Final$ExpectedMoreThanOne.Nd6),2); summary(Final$ExpectedMoreThanOne.Nd6)

Final = Final[order(Final$ExpectedMoreThanOne),]
breaks = seq(-1,1,0.05)
hist(Final$ExpectedMoreThanOne, xlim = c(-0.5,0.6), ylim = c(0,9), col = rgb(1,0.1,0.1,0.3), xlab = '', main = '', breaks = breaks); par(new = TRUE)
hist(Final$ExpectedMoreThanOne.Nd6, xlim = c(-0.5,0.6),ylim = c(0,9), col = rgb(0.1,0.1,0.1,0.3), xlab = 'log2(expected/unexpected)', main = 'germ', breaks = breaks);
abline(v=0, lt = 2, lwd = 2, col = 'black')

wilcox.test(Final$ExpectedMoreThanOne, mu = 0)     # 4.939e-05
wilcox.test(Final$ExpectedMoreThanOne.Nd6, mu = 0) # 0.0107

Final$DataSet = 'Germ'
SuperFinal=rbind(SuperFinal,Final)

#########
## GermConst - 12 genes: different constraints:
########

LowConst = read.table("../../Body/3Results/Alima02.AaAsymmetry.HumanTree.12GenesLowConst.ExpectedVsObservedAaChanges.txt", header = TRUE)
names(LowConst)
LowConst = LowConst[,c(2:8)]

HighConst = read.table("../../Body/3Results/Alima02.AaAsymmetry.HumanTree.12GenesHighConst.ExpectedVsObservedAaChanges.txt", header = TRUE)
names(HighConst)
HighConst = HighConst[,c(4,6,7,8)]
names(HighConst)[2:4]=c('ExpectedMoreThanOne.HighConst','NumberOfExpectedAaSubst.HighConst','NumberOfUnexpectedAaSubst.HighConst')
All = merge(LowConst,HighConst, all = TRUE)
#All[is.na(All$NumberOfExpectedAaSubst.HighPh),]$NumberOfExpectedAaSubst.HighPh<-0
#All[is.na(All$NumberOfUnexpectedAaSubst.HighPh),]$NumberOfUnexpectedAaSubst.Nd6<-0
# All = All[order(-All$L),] - к черту
#names(All)

All$Ratio = All$ExpectedMoreThanOne / All$ExpectedMoreThanOne.HighConst
summary(All$Ratio) # a bit more than one (because LocConst to HighConst)
wilcox.test(All$Ratio, mu = 1)
wilcox.test(All$ExpectedMoreThanOne,All$ExpectedMoreThanOne.HighConst, paired = TRUE)
# cor.test(All$Ratio,All$GranthamDistance, method = 'spearman')

#Final = All[,c(1,6,7,5,9,10,8)]
#Final$ExpectedMoreThanOne = round(log2(Final$ExpectedMoreThanOne),2); summary(Final$ExpectedMoreThanOne)
#Final$ExpectedMoreThanOne.HighPh = round(log2(Final$ExpectedMoreThanOne.HighPh),2); summary(Final$ExpectedMoreThanOne.HighPh)

#Final = Final[order(Final$ExpectedMoreThanOne),] 
# dev.off()
#breaks = seq(-1,1,0.05)
#hist(Final$ExpectedMoreThanOne, xlim = c(-0.3,0.6), ylim = c(0,9), col = rgb(1,0.1,0.1,0.3), xlab = '', main = '', breaks = breaks); par(new = TRUE)
#hist(Final$ExpectedMoreThanOne.HighPh, xlim = c(-0.3,0.6),ylim = c(0,9), col = rgb(0.1,0.1,0.1,0.3), xlab = 'log2(expected/unexpected)', main = '', breaks = breaks);
#abline(v=0, lt = 2, lwd = 2, col = 'black')

#wilcox.test(Final$ExpectedMoreThanOne, mu = 0)     # 9.199e-05
#wilcox.test(Final$ExpectedMoreThanOne.HighPh, mu = 0) # 0.0007255

GermConstr = Final

########
## VertebratePolym (one type of files - all genes without ND6)
########

for (set in 1:6)
{ # set = 1
  if (set == 1) {input = '../../Body/3Results/Alima05.AaAsymmetryVertebratePolymorphisms.r.AllClasses.txt'; class = 'AllClasses'}
  if (set == 2) {input = '../../Body/3Results/Alima05.AaAsymmetryVertebratePolymorphisms.r.Actinopterygii.txt'; class = 'Actinopterygii'}
  if (set == 3) {input = '../../Body/3Results/Alima05.AaAsymmetryVertebratePolymorphisms.r.Amphibia.txt'; class = 'Amphibia'}
  if (set == 4) {input = '../../Body/3Results/Alima05.AaAsymmetryVertebratePolymorphisms.r.Reptilia.txt'; class = 'Reptilia'}
  if (set == 5) {input = '../../Body/3Results/Alima05.AaAsymmetryVertebratePolymorphisms.r.Mammalia.txt'; class = 'Mammalia'}
  if (set == 6) {input = '../../Body/3Results/Alima05.AaAsymmetryVertebratePolymorphisms.r.Aves.txt'; class = 'Aves'}    
  
AllExceptNd6 = read.table(input, header = TRUE)
names(AllExceptNd6)
AllExceptNd6 = AllExceptNd6[,c(2:8)]
AllExceptNd6 = AllExceptNd6[]

names(AllExceptNd6)
Final = AllExceptNd6[,c(3,6,7,5)]
Final$ExpectedMoreThanOne = round(log2(Final$ExpectedMoreThanOne),2); summary(Final$ExpectedMoreThanOne)
Final = Final[order(Final$ExpectedMoreThanOne),]

wilcox.test(Final$ExpectedMoreThanOne, mu = 0)     # 6.333e-05
breaks = seq(-1,3,0.25)
hist(Final$ExpectedMoreThanOne, xlim = c(-1,3), ylim = c(0,9), col = rgb(1,0.1,0.1,0.3), xlab = '', main = class, breaks = breaks); # par(new = TRUE)
abline(v=0, lt = 2, lwd = 2, col = 'black')

names(SuperFinal)
Final$NumberOfExpectedAaSubst.Nd6 = '';
Final$NumberOfUnexpectedAaSubst.Nd6 = '';
Final$ExpectedMoreThanOne.Nd6 = '';
Final$DataSet = class  # rbind in the loop!!!

SuperFinal=rbind(SuperFinal,Final)
}

###### PRINT OUR EVERYTHING

write.table(SuperFinal,'../../Body/3Results/Alima08.CancersGermVertPolymorphismsFinalTableAndHists.r.txt') 
# copy and paste into LibreCalc => copy and special paste (rft)

dev.off()
