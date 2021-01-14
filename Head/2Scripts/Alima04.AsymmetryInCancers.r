rm(list=ls(all=TRUE)) 

pdf("../../Body/4Figures/Alima04.AsymmetryInVertebratePolymorphisms.r.pdf")

data = read.table("../../Body/1Raw/VertebratePolymorphisms.MutSpecData.OnlyFourFoldDegCytb.txt", header = TRUE)

FourNonZeroTransitions = data[data$T_C > 0 & data$C_T > 0 & data$A_G > 0 & data$G_A > 0,]
FourNonZeroTransitions$Ah_GhDivideByTh_Ch = FourNonZeroTransitions$T_C / FourNonZeroTransitions$A_G
FourNonZeroTransitions$Ch_ThDivideByGh_Ah = FourNonZeroTransitions$G_A / FourNonZeroTransitions$C_T
summary(FourNonZeroTransitions$Ah_GhDivideByTh_Ch)
summary(FourNonZeroTransitions$Ch_ThDivideByGh_Ah)
boxplot(FourNonZeroTransitions$Ah_GhDivideByTh_Ch,FourNonZeroTransitions$Ch_ThDivideByGh_Ah, names = c('Ah_GhDivideByTh_Ch','Ch_ThDivideByGh_Ah'), notch = TRUE, outline = FALSE, main = 'all chordata', ylab = 'asymmetry')

GL = read.table("../../Body/1Raw/GenerationLenghtforMammals.xlsx.txt", header = TRUE, sep = '\t')
GL$Species = gsub(' ','_',GL$Scientific_name)

FourNonZeroTransitionsMammals = merge(FourNonZeroTransitions,GL, by = 'Species')
summary(FourNonZeroTransitionsMammals$Ah_GhDivideByTh_Ch)
summary(FourNonZeroTransitionsMammals$Ch_ThDivideByGh_Ah)
boxplot(FourNonZeroTransitionsMammals$Ah_GhDivideByTh_Ch,FourNonZeroTransitionsMammals$Ch_ThDivideByGh_Ah, names = c('Ah_GhDivideByTh_Ch','Ch_ThDivideByGh_Ah'), notch = TRUE, outline = FALSE, main = 'mammals', ylab = 'asymmetry')

cor.test(FourNonZeroTransitionsMammals$GenerationLength_d,FourNonZeroTransitionsMammals$Ah_GhDivideByTh_Ch, method = 'spearman') # nothing
cor.test(FourNonZeroTransitionsMammals$GenerationLength_d,FourNonZeroTransitionsMammals$Ch_ThDivideByGh_Ah, method = 'spearman') # a bit positive

median(FourNonZeroTransitionsMammals$GenerationLength_d) 
FourNonZeroTransitionsMammalsShortLived = FourNonZeroTransitionsMammals[FourNonZeroTransitionsMammals$GenerationLength_d < median(FourNonZeroTransitionsMammals$GenerationLength_d),]
FourNonZeroTransitionsMammalsLongLived = FourNonZeroTransitionsMammals[FourNonZeroTransitionsMammals$GenerationLength_d >= median(FourNonZeroTransitionsMammals$GenerationLength_d),]

summary(FourNonZeroTransitionsMammalsShortLived$Ah_GhDivideByTh_Ch)
summary(FourNonZeroTransitionsMammalsShortLived$Ch_ThDivideByGh_Ah)

summary(FourNonZeroTransitionsMammalsLongLived$Ah_GhDivideByTh_Ch)
summary(FourNonZeroTransitionsMammalsLongLived$Ch_ThDivideByGh_Ah)

dev.off()
