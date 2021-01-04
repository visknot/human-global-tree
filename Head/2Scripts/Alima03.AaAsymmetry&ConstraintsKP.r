rm(list=ls(all=TRUE)) 

pdf("../../Body/4Figures/Alima03.AaAsymmetry&ConstraintsKP.pdf",)

data = read.csv("../../Body/2Derived/fulltreeCodons.csv", header = TRUE, sep = ";")

#### FILTER OUT: 
data <- subset(data, synonymous == "non-synonymous" & derived_aa != "Ambiguous" & derived_aa != "Asn/Asp" & derived_aa != "Gln/Glu" & derived_aa != "Leu/Ile" & gene_info != "mRNA_ND6")  #убираю лишнее
table(data$derived_aa)
data$AaSubst = paste(data$ancestral_aa,data$derived_aa, sep = '_')

## frequencies of changes in a given position versus AaSymmetry 

#temp = data[(data$ancestral_aa == 'Leu' & data$derived_aa == 'Pro' ) | (data$ancestral_aa == 'Pro' & data$derived_aa == 'Leu') ,]
#temp$AaSubst = paste(temp$ancestral_aa,temp$derived_aa, sep = '_'); nrow(temp)

EvolConstOfPositions = data.frame(table(data$position)); 
names(EvolConstOfPositions)=c('position','HowManySubstitutions')
EvolConstOfPositions = EvolConstOfPositions[order(-EvolConstOfPositions$HowManySubstitutions),]
EvolConstOfPositionsVec1 = EvolConstOfPositions[EvolConstOfPositions$HowManySubstitutions <= quantile(EvolConstOfPositions$HowManySubstitutions,0.95),]$position
EvolConstOfPositionsVec2 = EvolConstOfPositions[EvolConstOfPositions$HowManySubstitutions > quantile(EvolConstOfPositions$HowManySubstitutions,0.95),]$position

length(EvolConstOfPositionsVec1)
length(EvolConstOfPositionsVec2)

temp1 = data[data$position %in% EvolConstOfPositionsVec1,]; nrow(temp1) # 158
temp2 = data[data$position %in% EvolConstOfPositionsVec2,]; nrow(temp2) # 2781

data1 = data.frame(table(temp1$AaSubst)); 
data2 = data.frame(table(temp2$AaSubst)); 

#### load expectations: 
exp = read.table("../../Body/2Derived/ExpectedAaSubstitutionDirection.txt", header = TRUE, sep = " ")
exp$ExpectedAminoAcidSubstBias = gsub('>','_',exp$ExpectedAminoAcidSubstBias)
for (i in 1:nrow(exp)) 
{ # i = 1
exp$ContrExpectedAminoAcidSubstBias[i] = paste(unlist(strsplit(exp$ExpectedAminoAcidSubstBias[i],"_"))[2],unlist(strsplit(exp$ExpectedAminoAcidSubstBias[i],"_"))[1],sep = '_')
}

#### derive numbers for expected and contrexpected:
exp$ExpectedAminoAcidSubstBiasNumberRare = 0
exp$ExpectedAminoAcidSubstBiasNumberCommon = 0
exp$ContrExpectedAminoAcidSubstBiasNumberRare = 0
exp$ContrExpectedAminoAcidSubstBiasNumberCommon = 0

for (i in 1:nrow(exp))
{ # i = 1
  if (nrow(data1[data1$Var1 == exp$ExpectedAminoAcidSubstBias[i],])>0) {exp$ExpectedAminoAcidSubstBiasNumberRare[i] = data1[data1$Var1 == exp$ExpectedAminoAcidSubstBias[i],]$Freq}
  if (nrow(data2[data2$Var1 == exp$ExpectedAminoAcidSubstBias[i],])>0) {exp$ExpectedAminoAcidSubstBiasNumberCommon[i] = data2[data2$Var1 == exp$ExpectedAminoAcidSubstBias[i],]$Freq}
  if (nrow(data1[data1$Var1 == exp$ContrExpectedAminoAcidSubstBias[i],])>0) {exp$ContrExpectedAminoAcidSubstBiasNumberRare[i] = data1[data1$Var1 == exp$ContrExpectedAminoAcidSubstBias[i],]$Freq}
  if (nrow(data2[data2$Var1 == exp$ContrExpectedAminoAcidSubstBias[i],])>0) {exp$ContrExpectedAminoAcidSubstBiasNumberCommon[i] = data2[data2$Var1 == exp$ContrExpectedAminoAcidSubstBias[i],]$Freq}
}

exp$Odds = 0
for (i in 1:nrow(exp))
{ # i = 1
  X = rbind(c(exp$ExpectedAminoAcidSubstBiasNumberRare[i],exp$ContrExpectedAminoAcidSubstBiasNumberRare[i]),c(exp$ExpectedAminoAcidSubstBiasNumberCommon[i],exp$ContrExpectedAminoAcidSubstBiasNumberCommon[i]))
  FT = fisher.test(X, alternative = 'less')
  exp$Odds[i] = FT$estimate
}

# altogether:
X = rbind(c(sum(exp$ExpectedAminoAcidSubstBiasNumberRare),sum(exp$ContrExpectedAminoAcidSubstBiasNumberRare)),c(sum(exp$ExpectedAminoAcidSubstBiasNumberCommon),sum(exp$ContrExpectedAminoAcidSubstBiasNumberCommon)))
FT = fisher.test(X, alternative = 'less')
FT

# trend is correct, but fisher is very weakly significant => cor.test? log regr(!!!!!!)? add grantham?

exp = exp[order(exp$GranthamDistance),]
names(exp)
ExpFinal <- exp[c("NuclPosInCodon", "NuclSubstLightChainNotation", "ExpectedAminoAcidSubstBias", "Odds", "GranthamDistance","ExpectedAminoAcidSubstBiasNumberRare","ContrExpectedAminoAcidSubstBiasNumberRare","ExpectedAminoAcidSubstBiasNumberCommon","ContrExpectedAminoAcidSubstBiasNumberCommon")]

cor.test(ExpFinal[ExpFinal$Odds >0,]$Odds,ExpFinal[ExpFinal$Odds >0,]$GranthamDistance, method = 'spearman') # POSITIVE!!!!!
plot(ExpFinal[ExpFinal$Odds >0,]$Odds,ExpFinal[ExpFinal$Odds >0,]$GranthamDistance)

wilcox.test(ExpFinal[ExpFinal$Odds >0 & ExpFinal$NuclSubstLightChainNotation == 'T>C',]$Odds,ExpFinal[ExpFinal$Odds >0 & ExpFinal$NuclSubstLightChainNotation == 'G>A',]$Odds)
boxplot(ExpFinal[ExpFinal$Odds >0 & ExpFinal$NuclSubstLightChainNotation == 'T>C',]$Odds,ExpFinal[ExpFinal$Odds >0 & ExpFinal$NuclSubstLightChainNotation == 'G>A',]$Odds, names = c('T>C','G>A'))

write.table(ExpFinal, "../../Body/3Results/Alima03.AaAsymmetry&ConstraintsKP.ExpFinal.txt")
dev.off()