rm(list=ls(all=TRUE))

pdf("../../Body/4Figures/Alima04.AsymmetryInCancers.r.pdf")

data = read.table("../../Body/1Raw/CancerDataFromCampbell/mtDNA_snv_Oct2016.txt", header = TRUE, sep = '\t')
for (i in 1:nrow(data)) {data$AaChanges[i] = unlist(strsplit(data$Annot[i],','))[5]}
data$ancestral_aa = gsub("\\d(.*)",'',data$AaChanges)
data$derived_aa = gsub("(.*)\\d",'',data$AaChanges)

### read mtDNA annotation by Victor:
#mt = read.table("../../Body/1Raw/mtNucAnnotation_MergeRazerV6TurboMaxPro.csv", header = TRUE, sep = ';')
### read cancers, keep only nsSNPs and derive ancestral_aa and derived_aa
#data=data[grepl('nsSNP', data$Annot),]
#data$NucSub = paste(data$ref,data$var, sep = '>')

### merge data with mt and rename Leu and Ser according to their codons
#data = merge(data,mt,by.x = 'position', by.y = 'pos')
#data$AaSub = paste(data$ancestral_aa,data$derived_aa, sep = '>')
#for (i in 1:nrow(data))
#{ # i = 1
#  if (data$ancestral_aa[i] == 'L' & paste(unlist(strsplit(data$codon[i],''))[1:2],collapse = '') == 'TT')
#  {data$ancestral_aa[i] = 'LeuTT'}
#  if (data$ancestral_aa[i] == 'L' & paste(unlist(strsplit(data$codon[i],''))[1:2],collapse = '') == 'CT')
#  {data$ancestral_aa[i] = 'LeuCT'}
#  if (data$ancestral_aa[i] == 'S' & paste(unlist(strsplit(data$codon[i],''))[1:2],collapse = '') == 'TC')
#  {data$ancestral_aa[i] = 'SerTC'}
#  if (data$ancestral_aa[i] == 'S' & paste(unlist(strsplit(data$codon[i],''))[1:2],collapse = '') == 'AG')
#  {data$ancestral_aa[i] = 'SerAG'}
#}

#table(data$ancestral_aa)  # till now there are some 'L' and 'S' => ND6 - parasha!!! victor!!!!
#temp = data[data$ancestral_aa == 'L' | data$ancestral_aa == 'S',]
#data$AaSub = paste(data$ancestral_aa,data$derived_aa, sep = '>')

# in case of derived aminoacids it is more difficult to derive them:
#1) I can assume that if from Gly to Ser it means: from Gly to SerAG etc..
#temp = data[data$derived_aa == 'S',]
#table(temp$AaSub)
# A>S     C>S     F>S     G>S LeuTT>S     N>S     P>S     R>S     T>S     W>S     Y>S 
# 3       3      76     168      11      24      25       2       4       4       2 

#2) if from Phe to Leu and NucSub is 'T>C' it means from Phe to LeuCT, if NucSub is 'C>A' => LeuTT
#temp = data[data$ancestral_aa == 'F' & data$derived_aa == 'L',]
#names(temp)
#temp$NucSub = paste(temp$ref,temp$var, sep = '>')
#table(temp$NucSub)

#for (i in 1:nrow(data))
#{
#  if (temp$AaSub[i] == 'Y>S') {temp$derived_aa[i] = 'SerTC'}
#  if (temp$AaSub[i] == 'Y>S') {temp$derived_aa[i] = 'SerTC'}
#}


####################################
### 2A: for all genes except ND6
####################################

data1=data[!grepl('MT-ND6', data$Annot),]

# from one letter to three letter code AA:
A = c("G","A","L","M","F","W","K","Q","E","S","P","V","I","C","Y","H","R","N","D","T","*")
AAA = c("Gly","Ala","Leu","Met","Phe","Trp","Lys","Gln","Glu","Ser","Pro","Val","Ile","Cys","Tyr","His","Arg","Asn","Asp","Thr","Stop")
AA = data.frame(A,AAA)
data1 = merge(data1,AA, by.x = 'ancestral_aa',by.y='A'); data1$ancestral_aa = data1$AAA
data1 = merge(data1,AA, by.x = 'derived_aa',by.y='A'); data1$derived_aa = data1$AAA.y
data1$FromTo = paste(data1$ancestral_aa,data1$derived_aa,sep = ">")
VecOfNames = c("Tier2","chrom","position","FromTo","ancestral_aa","derived_aa")
data1 = data1[colnames(data1) %in% VecOfNames]
data1$FreqOfAASub = 1

FromTo = aggregate(data1$FreqOfAASub, by = list(data1$FromTo), FUN = sum)
names(FromTo) = c('FromAncestralToDerived','NumberOfEvents')

##### ADD DUMMY MATRIX WITH ZEROES 
AllAa1 = data.frame(unique(data$ancestral_aa)); nrow(AllAa1); names(AllAa1) = c('AA1')
AllAa2 = data.frame(unique(data$ancestral_aa)); nrow(AllAa2); names(AllAa2) = c('AA2')
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
exp = read.table("../../Body/2Derived/ExpectedAaSubstitutionDirection20AA.txt", header = TRUE, sep = " ")
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
summary(FromTo4$ExpectedMoreThanOne) # 0.5789  4.9191  7.3214     Inf 12.0074     Inf     167

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
wilcox.test(FromTo5$ExpectedMoreThanOne, mu = 1) # PAPER !!!! p-value = 7.153e-07

FromTo5$TotalSubst = FromTo5$NumberOfExpectedAaSubst + FromTo5$NumberOfUnexpectedAaSubst
Short = FromTo5[!is.na(FromTo5$ExpectedMoreThanOne),]
Short$DummyAhGh = 0
for (i in 1:nrow(Short))
{
  if (Short$NuclSubstLightChainNotation[i] == 'T>C') {Short$DummyAhGh[i] = 1}
  if (Short$NuclSubstLightChainNotation[i] == 'G>A') {Short$DummyAhGh[i] = 0}
}

## bias is higher in case of low Grantham!?
cor.test(Short$ExpectedMoreThanOne,Short$GranthamDistance, method = 'spearman') # a bit negative!!! 
cor.test(Short$TotalSubst,Short$GranthamDistance, method = 'spearman') # a bit negative
cor.test(Short$ExpectedMoreThanOne,Short$TotalSubst, method = 'spearman') 
## bias is higher in Ah>Gh than Ch>Th!!! why? Ah>Gh is more asymmetric on average???!!! may be yes!? check human global tree and mammalian average piechart
wilcox.test(Short[Short$NuclSubstLightChainNotation == 'G>A',]$ExpectedMoreThanOne,Short[Short$NuclSubstLightChainNotation == 'T>C',]$ExpectedMoreThanOne)
boxplot(Short[Short$NuclSubstLightChainNotation == 'G>A',]$ExpectedMoreThanOne,Short[Short$NuclSubstLightChainNotation == 'T>C',]$ExpectedMoreThanOne, names = c('Ch>Th','Ah>Gh'), ylab = 'expected shift')

Short = Short[,-c(2:6)]
write.table(Short,"../../Body/3Results/Alima04.AsymmetryInCancers.12genes.txt")

####################################
### 2B: for ND6 only
####################################
data1=data[grepl('MT-ND6', data$Annot),]

# from one letter to three letter code AA:
A = c("G","A","L","M","F","W","K","Q","E","S","P","V","I","C","Y","H","R","N","D","T","*")
AAA = c("Gly","Ala","Leu","Met","Phe","Trp","Lys","Gln","Glu","Ser","Pro","Val","Ile","Cys","Tyr","His","Arg","Asn","Asp","Thr","Stop")
AA = data.frame(A,AAA)
data1 = merge(data1,AA, by.x = 'ancestral_aa',by.y='A'); data1$ancestral_aa = data1$AAA
data1 = merge(data1,AA, by.x = 'derived_aa',by.y='A'); data1$derived_aa = data1$AAA.y
data1$FromTo = paste(data1$ancestral_aa,data1$derived_aa,sep = ">")
VecOfNames = c("Tier2","chrom","position","FromTo","ancestral_aa","derived_aa")
data1 = data1[colnames(data1) %in% VecOfNames]
data1$FreqOfAASub = 1

FromTo = aggregate(data1$FreqOfAASub, by = list(data1$FromTo), FUN = sum)
names(FromTo) = c('FromAncestralToDerived','NumberOfEvents')

##### ADD DUMMY MATRIX WITH ZEROES 
AllAa1 = data.frame(unique(data$ancestral_aa)); nrow(AllAa1); names(AllAa1) = c('AA1')
AllAa2 = data.frame(unique(data$ancestral_aa)); nrow(AllAa2); names(AllAa2) = c('AA2')
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
exp = read.table("../../Body/2Derived/ExpectedAaSubstitutionDirection20AA.txt", header = TRUE, sep = " ")
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
summary(FromTo4$ExpectedMoreThanOne) # 0.5789  4.9191  7.3214     Inf 12.0074     Inf     167

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

summary(FromTo5$ExpectedMoreThanOne)
wilcox.test(FromTo5$ExpectedMoreThanOne, mu = 1) # PAPER !!!! p-value = 7.153e-07

FromTo5$TotalSubst = FromTo5$NumberOfExpectedAaSubst + FromTo5$NumberOfUnexpectedAaSubst
Short = FromTo5[!is.na(FromTo5$ExpectedMoreThanOne),]
Short$DummyAhGh = 0
for (i in 1:nrow(Short))
{
  if (Short$NuclSubstLightChainNotation[i] == 'T>C') {Short$DummyAhGh[i] = 1}
  if (Short$NuclSubstLightChainNotation[i] == 'G>A') {Short$DummyAhGh[i] = 0}
}

## bias is higher in case of low Grantham!?
cor.test(Short$ExpectedMoreThanOne,Short$GranthamDistance, method = 'spearman') # a bit negative!!! 
cor.test(Short$TotalSubst,Short$GranthamDistance, method = 'spearman') # a bit negative
cor.test(Short$ExpectedMoreThanOne,Short$TotalSubst, method = 'spearman') 
## bias is higher in Ah>Gh than Ch>Th!!! why? Ah>Gh is more asymmetric on average???!!! may be yes!? check human global tree and mammalian average piechart
wilcox.test(Short[Short$NuclSubstLightChainNotation == 'G>A',]$ExpectedMoreThanOne,Short[Short$NuclSubstLightChainNotation == 'T>C',]$ExpectedMoreThanOne)

Short = Short[,-c(2:6)]
write.table(Short,"../../Body/3Results/Alima04.AsymmetryInCancers.ND6.txt")
dev.off()

