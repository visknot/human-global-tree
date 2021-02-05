rm(list=ls(all=TRUE)) 

data = read.table("../../Body/1Raw/VertebratePolymorphisms.MutSpecAminoAcidSubsWithRefSeqData.txt", header = TRUE)
data$TypesOfAASub = gsub('_','>',data$TypesOfAASub)
data$ancestral_aa = gsub(">(.*)",'',data$TypesOfAASub)
data$derived_aa = gsub("(.*)>",'',data$TypesOfAASub)

?????????????



################ OLD
pdf("../../Body/4Figures/Alima05.AaAsymmetryVertebratePolymorphisms.r.pdf")
data = read.table("../../Body/1Raw/VertebratePolymorphisms.MutSpecAminoAcidSubs.txt", header = TRUE)
data$TypesOfAASub = gsub('_','>',data$TypesOfAASub)
data$ancestral_aa = gsub(">(.*)",'',data$TypesOfAASub)
data$derived_aa = gsub("(.*)>",'',data$TypesOfAASub)

table(data$Gene) 
#  ATP6  ATP8  COX1  COX2  COX3  CytB   ND1   ND2   ND3   ND4  ND4L 
#  2436   706   996   839   845 17200  2194  7746   795  2916   548 
# there is no ND6 and ND5 and the most common is CYTB.
# => can we just count amino acid frequencies from all the rest 11 genes from RefSeqs and normalize on them? 
# in this way we don't take into account gene-specific AA content, but we take the maximal statistics (if we want - we can rerun the same only for CYTB later)

table(data$ancestral_aa)
table(data$derived_aa)

nrow(data) # 37221 reconstructed amino acid substitutions
length(unique(data$Species)) # 2020 species
table(data$Class) # 
#Actinopterygii           Amphibia               Aves   Branchiostomidae     Chondrichthyes          Cladistia           Mammalia Petromyzontiformes 
#10344               2106               3344                303                215                 58              13792                 67 
#Reptilia 
#5858 

### here choose classes!!!
for (set in 1:6)
{ # set = 1
if (set == 1) {data1 = data; output = '../../Body/3Results/Alima05.AaAsymmetryVertebratePolymorphisms.r.AllClasses.txt'}
if (set == 2) {data1 = data[data$Class == 'Actinopterygii',]; output = '../../Body/3Results/Alima05.AaAsymmetryVertebratePolymorphisms.r.Actinopterygii.txt'}
if (set == 3) {data1 = data[data$Class == 'Amphibia',]; output = '../../Body/3Results/Alima05.AaAsymmetryVertebratePolymorphisms.r.Amphibia.txt'}
if (set == 4) {data1 = data[data$Class == 'Reptilia',]; output = '../../Body/3Results/Alima05.AaAsymmetryVertebratePolymorphisms.r.Reptilia.txt'}
if (set == 5) {data1 = data[data$Class == 'Mammalia',]; output = '../../Body/3Results/Alima05.AaAsymmetryVertebratePolymorphisms.r.Mammalia.txt'}  
if (set == 6) {data1 = data[data$Class == 'Aves',]; output = '../../Body/3Results/Alima05.AaAsymmetryVertebratePolymorphisms.r.Aves.txt'}    

FromTo = aggregate(data1$FreqOfSub, by = list(data1$TypesOfAASub), FUN = sum)
names(FromTo) = c('FromAncestralToDerived','NumberOfEvents')

##### ADD DUMMY MATRIX WITH ZEROES 
AllAa1 = data.frame(unique(data1$ancestral_aa)); nrow(AllAa1); names(AllAa1) = c('AA1')
AllAa2 = data.frame(unique(data1$ancestral_aa)); nrow(AllAa2); names(AllAa2) = c('AA2')
DummyZeroes = merge(AllAa1,AllAa2)
DummyZeroes$FromAncestralToDerived=paste(DummyZeroes$AA1,DummyZeroes$AA2,sep='>')
DummyZeroes$NumberOfEvents = 0
DummyZeroes = DummyZeroes[c(3,4)]

##### rbind and aggregate FromTo and DummyZeroes
FromTo1 = rbind(FromTo,DummyZeroes)
FromTo2 = aggregate(FromTo1$NumberOfEvents, by = list(FromTo1$FromAncestralToDerived), FUN = sum)
names(FromTo2)=c('FromAncestralToDerived','NumberOfEvents')
FromTo2$From = gsub(">.*",'',FromTo2$FromAncestralToDerived)
FromTo2$To = gsub(".*>",'',FromTo2$FromAncestralToDerived)
FromTo2 = FromTo2[FromTo2$From != FromTo2$To,]; nrow(FromTo2) # filter out synonymous changes
# for (i in 1:nrow(FromTo2)) {FromTo2$AaPairId[i] = paste(sort(unlist(strsplit(FromTo2$ExpectedAminoAcidSubstBias[i],'>'))),collapse = '>') }

##### DERIVE UNIQUE IDENTIFIER FOR EACH PAIR OF AMINOACIDS BASED ON ALPHABET:
for (i in 1:nrow(FromTo2)) {FromTo2$AaPairId[i] = paste(sort(unlist(strsplit(FromTo2$FromAncestralToDerived[i],'>'))),collapse = '>') }

##### from 420 to 210 rows and NumberOfEvents1To2
A = FromTo2[FromTo2$FromAncestralToDerived == FromTo2$AaPairId,]; A= A[c(1,2,5)]; names(A)=c('FromAncestralToDerived1','NumberOfEvents1','AaPairId')
B = FromTo2[FromTo2$FromAncestralToDerived != FromTo2$AaPairId,]; B= B[c(1,2,5)]; names(B)=c('FromAncestralToDerived2','NumberOfEvents2','AaPairId')
FromTo3 = merge(A,B)
FromTo3$NumberOfEvents1To2 = FromTo3$NumberOfEvents1 / FromTo3$NumberOfEvents2 

##### read expectations, add AaPairId and merge
exp = read.table("../../Body/2Derived/ExpectedAaSubstitutionDirection22AA.txt", header = TRUE, sep = " ")
for (i in 1:nrow(exp)) {exp$AaPairId[i] = paste(sort(unlist(strsplit(exp$ExpectedAminoAcidSubstBias[i],'>'))),collapse = '>') }
# exp=exp[c(3,4)]
FromTo4 = merge(FromTo3,exp, by = 'AaPairId', all.x = TRUE)

##### 
FromTo4$ExpectedMoreThanOne=NA
for (i in 1:nrow(FromTo4))
{ # i = 17
  if (!is.na(FromTo4$ExpectedAminoAcidSubstBias[i]) & FromTo4$ExpectedAminoAcidSubstBias[i] == FromTo4$FromAncestralToDerived1[i]) {FromTo4$ExpectedMoreThanOne[i] = FromTo4$NumberOfEvents1[i] / FromTo4$NumberOfEvents2[i]}
  if (!is.na(FromTo4$ExpectedAminoAcidSubstBias[i]) & FromTo4$ExpectedAminoAcidSubstBias[i] == FromTo4$FromAncestralToDerived2[i]) {FromTo4$ExpectedMoreThanOne[i] = FromTo4$NumberOfEvents2[i] / FromTo4$NumberOfEvents1[i]}
}

FromTo4 = FromTo4[order(FromTo4$ExpectedMoreThanOne),]

# delete variants from stop or to stop (too rare and probably deleterious => mistakes):
FromTo5 = FromTo4[!grepl('Stop',FromTo4$AaPairId),] 

summary(FromTo5$ExpectedMoreThanOne)

# derive 'NumberOfExpectedAaSubst' 'NumberOfUnexpectedAaSubst'
FromTo5$NumberOfExpectedAaSubst = 0
FromTo5$NumberOfUnexpectedAaSubst = 0
for (i in 1:nrow(FromTo5[!is.na(FromTo5$ExpectedMoreThanOne),]))
{ # i = 1
  if (FromTo5$FromAncestralToDerived1[i] == FromTo5$ExpectedAminoAcidSubstBias[i]) {FromTo5$NumberOfExpectedAaSubst[i] = FromTo5$NumberOfEvents1[i]; FromTo5$NumberOfUnexpectedAaSubst[i] = FromTo5$NumberOfEvents2[i];}
  if (FromTo5$FromAncestralToDerived2[i] == FromTo5$ExpectedAminoAcidSubstBias[i]) {FromTo5$NumberOfExpectedAaSubst[i] = FromTo5$NumberOfEvents2[i]; FromTo5$NumberOfUnexpectedAaSubst[i] = FromTo5$NumberOfEvents1[i];}
}

hist(FromTo5$ExpectedMoreThanOne,breaks = 100)
wilcox.test(FromTo5$ExpectedMoreThanOne, mu = 1) # PAPER !!!! p-value = 3.338e-06

FromTo5$TotalSubst = FromTo5$NumberOfExpectedAaSubst + FromTo5$NumberOfUnexpectedAaSubst
Short = FromTo5[!is.na(FromTo5$ExpectedMoreThanOne),]
Short$DummyAhGh = 0
for (i in 1:nrow(Short))
{
  if (Short$NuclSubstLightChainNotation[i] == 'T>C') {Short$DummyAhGh[i] = 1}
  if (Short$NuclSubstLightChainNotation[i] == 'G>A') {Short$DummyAhGh[i] = 0}
}

cor.test(Short$ExpectedMoreThanOne,Short$GranthamDistance, method = 'spearman') # negative!!! 
cor.test(Short$TotalSubst,Short$GranthamDistance, method = 'spearman') # a bit negative
cor.test(Short$ExpectedMoreThanOne,Short$TotalSubst, method = 'spearman') 
summary(lm(Short$ExpectedMoreThanOne ~ Short$TotalSubst*Short$GranthamDistance))
summary(lm(Short[Short$NuclSubstLightChainNotation == 'T>C',]$ExpectedMoreThanOne ~ Short[Short$NuclSubstLightChainNotation == 'T>C',]$TotalSubst + Short[Short$NuclSubstLightChainNotation == 'T>C',]$GranthamDistance))

## bias is higher in Ah>Gh than Ch>Th!!! why? Ah>Gh is more asymmetric on average???!!! may be yes!? check human global tree and mammalian average piechart
wilcox.test(Short[Short$NuclSubstLightChainNotation == 'G>A',]$ExpectedMoreThanOne,Short[Short$NuclSubstLightChainNotation == 'T>C',]$ExpectedMoreThanOne)
boxplot(Short[Short$NuclSubstLightChainNotation == 'G>A',]$ExpectedMoreThanOne,Short[Short$NuclSubstLightChainNotation == 'T>C',]$ExpectedMoreThanOne, names = c('Ch>Th','Ah>Gh'), ylab = 'expected shift')

summary(lm(Short$ExpectedMoreThanOne ~ Short$DummyAhGh)) # significant

Short = Short[,-c(2:6)]
write.table(Short,output)
}
dev.off()
