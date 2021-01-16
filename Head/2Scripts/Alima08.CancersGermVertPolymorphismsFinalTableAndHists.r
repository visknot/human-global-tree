###################################
##### prepare final cancer table and hist - make it a separate script, collecting cancer, germ data and vertebrate polymorphisms
###################################

rm(list=ls(all=TRUE)) 

pdf("../../Body/4Figures/Alima08.CancersGermVertPolymorphismsFinalTableAndHists.r.pdf")

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
All = All[order(-All$GranthamDistance),]
names(All)
FinalCancer = All[,c(1,6,7,5,9,10,8)]
FinalCancer$ExpectedMoreThanOne = round(log2(FinalCancer$ExpectedMoreThanOne),2); summary(FinalCancer$ExpectedMoreThanOne)
FinalCancer$ExpectedMoreThanOne.Nd6 = round(log2(FinalCancer$ExpectedMoreThanOne.Nd6),2); summary(FinalCancer$ExpectedMoreThanOne.Nd6)

hist(FinalCancer$ExpectedMoreThanOne)
hist(FinalCancer$ExpectedMoreThanOne.Nd6)

hist(FinalCancer$ExpectedMoreThanOne, xlim = c(-4,8), ylim = c(0,9), col = rgb(1,0.1,0.1,0.3), xlab = '', main = ''); par(new = TRUE)
hist(FinalCancer$ExpectedMoreThanOne.Nd6, xlim = c(-4,8),ylim = c(0,9), col = rgb(0.1,0.1,0.1,0.3), xlab = 'log2(expected/unexpected)', main = '');
abline(v=0, lt = 2, lwd = 2, col = 'black')

wilcox.test(FinalCancer$ExpectedMoreThanOne, mu = 0)     # 4.939e-05
wilcox.test(FinalCancer$ExpectedMoreThanOne.Nd6, mu = 0) # 0.02225

######
## GERM
#####






write.table(FinalCancer,'../../Body/3Results/Alima08.CancersGermVertPolymorphismsFinalTableAndHists.r.txt') 



dev.off()
