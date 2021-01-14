rm(list=ls(all=TRUE)) 

pdf("../../Body/4Figures/Alima03.AaAsymmetry&ConstraintsKP.pdf",)

data = read.csv("../../Body/2Derived/fulltreeCodons.csv", header = TRUE, sep = ";")
data = data[data$note == 'normal',] # filter out everything except protein-coding mutations:
table(data$ancestral_aa)
VecOfNormalAa = unique(data$ancestral_aa); length(VecOfNormalAa)
table(data$derived_aa) # Ambiguous, Asn/Asp, Gln/Glu, Leu/Ile => why all noisy AA only among derived? All of them are on external branches?
data = data[data$derived_aa %in% VecOfNormalAa,]

#### mutations are only a t g c;  Derive Subst and Context (two before and two after)
# Subst
ExtractThird<-function(x) {unlist(strsplit(x,''))[3]}
data$temp1 = apply(as.matrix(data$ancestor),1,FUN = ExtractThird); data = data[data$temp1 %in% c('A','T','G','C'),]
data$temp2 = apply(as.matrix(data$descendant),1,FUN = ExtractThird); data = data[data$temp2 %in% c('A','T','G','C'),]
data$Subst = paste(data$temp1,data$temp2,sep='')
table(data$Subst)

# FILTER FOR THE SAME BACKGROUND (Pos1Anc == Pos1Der; Pos5Anc == Pos5Der; ):
# very first (first) and very last (fifth) should be the same (important for codons - we will garantie, that in the codon there is only one substitution)
nrow(data) # 292532
ExtractFirst<-function(x) {unlist(strsplit(x,''))[1]}
data$Pos1Anc = apply(as.matrix(data$ancestor),1,FUN = ExtractFirst);   data = data[data$Pos1Anc %in% c('a','t','g','c'),]
data$Pos1Der = apply(as.matrix(data$descendant),1,FUN = ExtractFirst); data = data[data$Pos1Der %in% c('a','t','g','c'),]
data=data[data$Pos1Anc == data$Pos1Der,]; nrow(data)

ExtractFifth<-function(x) {unlist(strsplit(x,''))[5]}
data$Pos5Anc = apply(as.matrix(data$ancestor),1,FUN = ExtractFifth);   data = data[data$Pos5Anc %in% c('a','t','g','c'),]
data$Pos5Der = apply(as.matrix(data$descendant),1,FUN = ExtractFifth); data = data[data$Pos5Der %in% c('a','t','g','c'),]
data=data[data$Pos5Anc == data$Pos5Der,]; nrow(data)

# Context
ExtractSecond<-function(x) {unlist(strsplit(x,''))[2]}
data$Pos2Anc = apply(as.matrix(data$ancestor),1,FUN = ExtractSecond); data = data[data$Pos2Anc %in% c('a','t','g','c'),]
data$Pos2Der = apply(as.matrix(data$descendant),1,FUN = ExtractSecond); data = data[data$Pos2Der %in% c('a','t','g','c'),]
data=data[data$Pos2Anc == data$Pos2Der,]

ExtractFourth<-function(x) {unlist(strsplit(x,''))[4]}
data$Pos4Anc = apply(as.matrix(data$ancestor),1,FUN = ExtractFourth); data = data[data$Pos4Anc %in% c('a','t','g','c'),]
data$Pos4Der = apply(as.matrix(data$descendant),1,FUN = ExtractFourth); data = data[data$Pos4Der %in% c('a','t','g','c'),]
data=data[data$Pos4Anc == data$Pos4Der,]

data$Context = paste(toupper(data$Pos2Anc),data$temp1,toupper(data$Pos4Anc),sep='')

#### filter out synonymous
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