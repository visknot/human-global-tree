#### READ EXPECTED, OBSERVED, GET A RATIO, WRITE AND PLOTG

rm(list=ls(all=TRUE))
library('dplyr')
library('gtools')
library(seqinr)

pdf("../../Body/4Figures/03A.MutSpecGermObsToExp.R.pdf", width = 22, height = 14)
par(mfrow=c(3,1))

###### MutSpec12
ColVec=c('gray',rgb(0,1,0,1),'gray','gray','gray',rgb(1,0,0,1),rgb(1,0,0,0.2),'gray','gray','gray',rgb(0,1,0,0.2),'gray'); length(ColVec)
Obs = read.table("../../Body/2Derived/02A.GlobalHumanTree.MutSpecsObserved192&12.MutSpecObs12.txt")
Exp = read.table("../../Body/2Derived/03.MutSpecMtDnaExpected.MutSpecExp12.txt")
Exp = Exp[Exp$SubstType == 'synonymous',]
ObsExp = merge(Obs,Exp, by = c('Subst','Genes','SubstType'))
VecOfBranches = unique(ObsExp$Branches)
MS12 = data.frame()
for (i in 1:length(VecOfBranches))
{ # i = 1
  temp = ObsExp[ObsExp$Branches == VecOfBranches[i],]
  temp$FreqObs = temp$NumbObs/sum(temp$NumbObs);
  temp$FreqExp = temp$NumbExp/sum(temp$NumbExp);
  temp$ObsToExpFreq = temp$FreqObs/temp$FreqExp
  MS12 = rbind(MS12,temp)
}

# add heavy chain Subst 
Comp <- function(x) {y = paste(comp(s2c(as.character(x))),collapse = ''); return(toupper(y))}
MS12$SubstHeavy = apply(as.matrix(MS12$Subst), 1, FUN = Comp)

write.table(MS12,"../../Body/3Results/03A.MutSpecGermObsToExp.R.MutSpec12.txt")

# plot MutSpec in light chain notation:
barplot(MS12[MS12$Branches == 'AllBranches',]$ObsToExpFreq, names = MS12[MS12$Branches == 'AllBranches',]$Subst, col = ColVec, main = paste("ObsFreq/ExpFreq SynSubst, AllBranches, Light Chain"), cex.names = 2, cex.axis = 2)
barplot(MS12[MS12$Branches == 'TerminalBranches',]$ObsToExpFreq, names = MS12[MS12$Branches == 'AllBranches',]$Subst, col = ColVec, main = paste("ObsFreq/ExpFreq SynSubst, TerminalBranches, Light Chain"), cex.names = 2, cex.axis = 2)
barplot(MS12[MS12$Branches == 'InternalBranches',]$ObsToExpFreq, names = MS12[MS12$Branches == 'AllBranches',]$Subst, col = ColVec, main = paste("ObsFreq/ExpFreq SynSubst, InternalBranches, Light Chain"), cex.names = 2, cex.axis = 2)

MS12 = MS12[order(MS12$Branches,MS12$SubstHeavy),]
# plot MutSpec in heavy chain notation:
barplot(MS12[MS12$Branches == 'AllBranches',]$ObsToExpFreq, names = MS12[MS12$Branches == 'AllBranches',]$SubstHeavy, col = ColVec, main = paste("ObsFreq/ExpFreq SynSubst, AllBranches, Heavy Chain"), cex.names = 2, cex.axis = 2)
barplot(MS12[MS12$Branches == 'TerminalBranches',]$ObsToExpFreq, names = MS12[MS12$Branches == 'AllBranches',]$SubstHeavy, col = ColVec, main = paste("ObsFreq/ExpFreq SynSubst, TerminalBranches, Heavy Chain"), cex.names = 2, cex.axis = 2)
barplot(MS12[MS12$Branches == 'InternalBranches',]$ObsToExpFreq, names = MS12[MS12$Branches == 'AllBranches',]$SubstHeavy, col = ColVec, main = paste("ObsFreq/ExpFreq SynSubst, InternalBranches, Heavy Chain"), cex.names = 2, cex.axis = 2)

####### MutSpec192
ColVec=c(rep('gray',16),rep(rgb(0,1,0,1),16),rep('gray',16),rep('gray',16),rep('gray',16),rep(rgb(1,0,0,1),16),rep(rgb(1,0,0,0.2),16),rep('gray',16),rep('gray',16),rep('gray',16),rep(rgb(0,1,0,0.2),16),rep('gray',16)); length(ColVec)
BorderVec = rep(c(rgb(0.1,0.1,0.1,0.1),rgb(0.1,0.1,0.1,0.1),rgb(0.1,0.1,0.1,0.1),rgb(0.1,0.1,0.1,0.1),
                  rgb(0.1,0.1,0.1,0.1),rgb(0.1,0.1,0.1,0.1),rgb(0.1,0.1,0.1,0.1),rgb(0.1,0.1,0.1,0.1),
                rgb( 0,0,0,1),rgb( 0,0,0,1),rgb( 0,0,0,1),rgb( 0,0,0,1),
                rgb(0.1,0.1,0.1,0.1),rgb(0.1,0.1,0.1,0.1),rgb(0.1,0.1,0.1,0.1),rgb(0.1,0.1,0.1,0.1)),12) 
length(BorderVec) # the third tetrade is with border (G is the first)

Obs = read.table("../../Body/2Derived/02A.GlobalHumanTree.MutSpecsObserved192&12.MutSpecObs192.txt")
Exp = read.table("../../Body/2Derived/03.MutSpecMtDnaExpected.MutSpecExp192.txt")
Exp = Exp[Exp$SubstType == 'synonymous',]
ObsExp = merge(Obs,Exp, by = c('names','Genes','SubstType'))
VecOfBranches = unique(ObsExp$Branches)

MS192 = data.frame()
for (i in 1:length(VecOfBranches))
{ # i = 3
  temp = ObsExp[ObsExp$Branches == VecOfBranches[i],]
  temp$FreqObs = temp$NumbObs/sum(temp$NumbObs);
  temp$FreqExp = temp$NumbExp/sum(temp$NumbExp);
  temp$ObsToExpFreq = temp$FreqObs/temp$FreqExp
  MS192 = rbind(MS192,temp)
}

MS192[MS192$ObsToExpFreq == Inf | is.na(MS192$ObsToExpFreq),]$ObsToExpFreq <- 0 # WHY THERE ARE SOME NA???? I AHVE A BIT OF IBSERVED BUT THERE IS NO EXPECTED!!!!!!!????
# MS192$names[1] = AC: AAA
SwitchNamesToHeavy <- function(x) {
  subst = unlist(strsplit(x,': '))[1];  subst = toupper(paste(comp(s2c(as.character(subst))),collapse = ''));
  context = unlist(strsplit(x,': '))[2];  context = toupper(paste(rev(comp(s2c(as.character(context)))),collapse = ''));
  paste(subst,context, sep = ': ')  }
MS192$NamesHeavy = apply(as.matrix(MS192$names), 1, FUN = SwitchNamesToHeavy)

write.table(MS192,"../../Body/3Results/03A.MutSpecGermObsToExp.R.MutSpec192.txt")

# plot in light chain notation 
barplot(MS192[MS192$Branches == 'AllBranches',]$ObsToExpFreq, names = MS192[MS192$Branches == 'AllBranches',]$names, col = ColVec, border = BorderVec, main = paste("ObsFreq/ExpFreq SynSubst, AllBranches, Light Chain"), cex.names = 1, cex.axis = 1, las = 2)
barplot(MS192[MS192$Branches == 'TerminalBranches',]$ObsToExpFreq, names = MS192[MS192$Branches == 'AllBranches',]$names, col = ColVec, border = BorderVec, main = paste("ObsFreq/ExpFreq SynSubst, TerminalBranches, Light Chain"), cex.names = 1, cex.axis = 1, las = 2)
barplot(MS192[MS192$Branches == 'InternalBranches',]$ObsToExpFreq, names = MS192[MS192$Branches == 'AllBranches',]$names, col = ColVec, border = BorderVec, main = paste("ObsFreq/ExpFreq SynSubst, InternalBranches, Light Chain"), cex.names = 1, cex.axis = 1, las = 2)

MS192 = MS192[order(MS192$NamesHeavy),]
BorderVec = rep(c(rgb(0.1,0.1,0.1,0.1),rgb(0,0,0,1),rgb(0.1,0.1,0.1,0.1),rgb(0.1,0.1,0.1,0.1)),48) # 192/4 The second colum is with border
barplot(MS192[MS192$Branches == 'AllBranches',]$ObsToExpFreq, names = MS192[MS192$Branches == 'AllBranches',]$NamesHeavy, col = ColVec, border = BorderVec, main = paste("ObsFreq/ExpFreq SynSubst, AllBranches, Heavy Chain"), cex.names = 1, cex.axis = 1, las = 2)
barplot(MS192[MS192$Branches == 'TerminalBranches',]$ObsToExpFreq, names = MS192[MS192$Branches == 'AllBranches',]$NamesHeavy, col = ColVec, border = BorderVec, main = paste("ObsFreq/ExpFreq SynSubst, TerminalBranches, Heavy Chain"), cex.names = 1, cex.axis = 1, las = 2)
barplot(MS192[MS192$Branches == 'InternalBranches',]$ObsToExpFreq, names = MS192[MS192$Branches == 'AllBranches',]$NamesHeavy, col = ColVec, border = BorderVec, main = paste("ObsFreq/ExpFreq SynSubst, InternalBranches, Heavy Chain"), cex.names = 1, cex.axis = 1, las = 2)

dev.off()
