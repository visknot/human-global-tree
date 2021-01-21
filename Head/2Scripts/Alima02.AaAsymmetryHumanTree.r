rm(list=ls(all=TRUE)) 

pdf("../../Body/4Figures/Alima02.AaAsymmetryHumanTree.pdf")

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
data <- subset(data, synonymous == "non-synonymous" & derived_aa != "Ambiguous" & derived_aa != "Asn/Asp" & derived_aa != "Gln/Glu" & derived_aa != "Leu/Ile")  #убираю лишнее
table(data$derived_aa)

head(data)

#### From Leu to LeuCT and LeuTT
#### From Ser to SerTC and SerAG
for (i in 1:nrow(data))
{ # i = 1
  if (data$ancestral_aa[i] == 'Leu')
  {
    if (tolower(data$ancestral_codon[i]) == 'ctt' | tolower(data$ancestral_codon[i]) == 'ctc' | tolower(data$ancestral_codon[i]) == 'cta' | tolower(data$ancestral_codon[i]) == 'ctg') {data$ancestral_aa[i] = 'LeuCT'}
    if (tolower(data$ancestral_codon[i]) == 'tta' | tolower(data$ancestral_codon[i]) == 'ttg') {data$ancestral_aa[i] = 'LeuTT'}
  }
  if (data$ancestral_aa[i] == 'Ser')
  {
    if (tolower(data$ancestral_codon[i]) == 'tct' | tolower(data$ancestral_codon[i]) == 'tcc' | tolower(data$ancestral_codon[i]) == 'tca' | tolower(data$ancestral_codon[i]) == 'tcg') {data$ancestral_aa[i] = 'SerTC'}
    if (tolower(data$ancestral_codon[i]) == 'agt' | tolower(data$ancestral_codon[i]) == 'agc') {data$ancestral_aa[i] = 'SerAG'}
  }
  if (data$derived_aa[i] == 'Leu')
  {
    if (tolower(data$derived_codon[i]) == 'ctt' | tolower(data$derived_codon[i]) == 'ctc' | tolower(data$derived_codon[i]) == 'cta' | tolower(data$derived_codon[i]) == 'ctg') {data$derived_aa[i] = 'LeuCT'}
    if (tolower(data$derived_codon[i]) == 'tta' | tolower(data$derived_codon[i]) == 'ttg') {data$derived_aa[i] = 'LeuTT'}
  }
  if (data$derived_aa[i] == 'Ser')
  {
    if (tolower(data$derived_codon[i]) == 'tct' | tolower(data$derived_codon[i]) == 'tcc' | tolower(data$derived_codon[i]) == 'tca' | tolower(data$derived_codon[i]) == 'tcg') {data$derived_aa[i] = 'SerTC'}
    if (tolower(data$derived_codon[i]) == 'agt' | tolower(data$derived_codon[i]) == 'agc') {data$derived_aa[i] = 'SerAG'}
  }  
}
table(data$ancestral_aa)
table(data$derived_aa)

data$FromTo = paste(data$ancestral_aa,data$derived_aa,sep = '>')

#### read Victor Annotation:

mt = read.table("../../Body/1Raw/mtNucAnnotation_MergeRazerV6TurboMaxPro.csv", header = TRUE, sep = ";")
cor.test(mt$freq,mt$AfHomGnomad, method = 'spearman') # super positive as expected
data=merge(data,mt,by.x = 'position', by.y = 'pos')
summary(data$PhyloP)
summary(data$freq)
summary(data$AfHomGnomad)
table(data$gene_info)
VecOfHighlyConstrainedGenes = c('mRNA_ND1','mRNA_ND4L','mRNA_CYTB','mRNA_ND5','mRNA_ND3','mRNA_ND1')
VecOfLowConstrainedGenes = c('mRNA_COX1','mRNA_ND4','mRNA_ATP6','mRNA_COX3','mRNA_COX2','mRNA_ATP8')
table(data$tm) # tm should be more constrained!? I think yes! 
#VecOfGenesWithLongTimeBeingSingleStranded = c('mRNA_CYTB','mRNA_ND5','mRNA_ND4','mRNA_ND4L','mRNA_ND3') # ND3 ND4L ND4 ND5 CytB
#VecOfGenesWithShortTimeBeingSingleStranded = c('mRNA_COX1','mRNA_COX2','mRNA_ATP8','mRNA_ATP6','mRNA_COX3') # COX1 COX2 ATP8 ATP6 COX3
VecOfGenesWithLongTimeBeingSingleStranded = c('mRNA_CYTB') # ND3 ND4L ND4 ND5 CytB
VecOfGenesWithShortTimeBeingSingleStranded = c('mRNA_COX1') # COX1 COX2 ATP8 ATP6 COX3

for (set in 1:14)
{ # i =4
  if (set == 1) {data1 = data[data$gene_info != 'mRNA_ND6',]; outfile = "../../Body/3Results/Alima02.AaAsymmetry.HumanTree.12Genes.ExpectedVsObservedAaChanges.txt"} # 12Genes
  if (set == 2) {data1 = data[data$gene_info == 'mRNA_ND6',]; outfile = "../../Body/3Results/Alima02.AaAsymmetry.HumanTree.ND6.ExpectedVsObservedAaChanges.txt"} # Nd6
# NO if (set == 3) {data1 = data[data$gene_info  != 'mRNA_ND6' & data$PhyloP < median(data$PhyloP),]; outfile = "../../Body/3Results/Alima02.AaAsymmetry.HumanTree.12GenesLowPhyloP.ExpectedVsObservedAaChanges.txt"} # 12GenesLowPhyloP
# NO if (set == 4) {data1 = data[data$gene_info  != 'mRNA_ND6' & data$PhyloP >= median(data$PhyloP),]; outfile = "../../Body/3Results/Alima02.AaAsymmetry.HumanTree.12GenesHighPhyloP.ExpectedVsObservedAaChanges.txt"} # 12GenesHighPhyloP
# NO if (set == 3) {data1 = data[data$gene_info  != 'mRNA_ND6' & data$freq >= quantile(data$freq,0.9),]; outfile = "../../Body/3Results/Alima02.AaAsymmetry.HumanTree.12GenesLowConst.ExpectedVsObservedAaChanges.txt"} # 12GenesLowPhyloP
# NO if (set == 4) {data1 = data[data$gene_info  != 'mRNA_ND6' & data$freq <  quantile(data$freq,0.5),]; outfile = "../../Body/3Results/Alima02.AaAsymmetry.HumanTree.12GenesHighConst.ExpectedVsObservedAaChanges.txt"} # 12GenesHighPhyloP
# a bit  if (set == 3) {data1 = data[data$gene_info  != 'mRNA_ND6' & data$gene_info %in% VecOfLowConstrainedGenes,]; outfile = "../../Body/3Results/Alima02.AaAsymmetry.HumanTree.12GenesLowConst.ExpectedVsObservedAaChanges.txt"} # 12GenesLowPhyloP 
# a bit  if (set == 4) {data1 = data[data$gene_info  != 'mRNA_ND6' & data$gene_info %in% VecOfHighlyConstrainedGenes,]; outfile = "../../Body/3Results/Alima02.AaAsymmetry.HumanTree.12GenesHighConst.ExpectedVsObservedAaChanges.txt"} # 12GenesHighPhyloP
# NO  if (set == 3) {data1 = data[data$gene_info  != 'mRNA_ND6' & data$tm == 'notransm',]; outfile = "../../Body/3Results/Alima02.AaAsymmetry.HumanTree.12GenesLowConst.ExpectedVsObservedAaChanges.txt"} # 12GenesLowPhyloP 
# NO  if (set == 4) {data1 = data[data$gene_info  != 'mRNA_ND6' & data$tm == 'transm',]; outfile = "../../Body/3Results/Alima02.AaAsymmetry.HumanTree.12GenesHighConst.ExpectedVsObservedAaChanges.txt"} # 12GenesHighPhyloP
#  if (set == 3) {data1 = data[data$gene_info  != 'mRNA_ND6' & data$gene_info %in% VecOfGenesWithLongTimeBeingSingleStranded,]; outfile = "../../Body/3Results/Alima02.AaAsymmetry.HumanTree.12GenesLowConst.ExpectedVsObservedAaChanges.txt"} # 12GenesLowPhyloP 
#  if (set == 4) {data1 = data[data$gene_info  != 'mRNA_ND6' & data$gene_info %in% VecOfGenesWithShortTimeBeingSingleStranded  ,]; outfile = "../../Body/3Results/Alima02.AaAsymmetry.HumanTree.12GenesHighConst.ExpectedVsObservedAaChanges.txt"} # 12GenesHighPhyloP
  if (set == 3) {data1 = data[data$gene_info == 'mRNA_CYTB',]; outfile = "../../Body/3Results/Alima02.AaAsymmetry.HumanTree.CYTB.ExpectedVsObservedAaChanges.txt"} 
  if (set == 4) {data1 = data[data$gene_info == 'mRNA_ND5',]; outfile = "../../Body/3Results/Alima02.AaAsymmetry.HumanTree.ND5.ExpectedVsObservedAaChanges.txt"}
  if (set == 5) {data1 = data[data$gene_info == 'mRNA_ND4',]; outfile = "../../Body/3Results/Alima02.AaAsymmetry.HumanTree.ND4.ExpectedVsObservedAaChanges.txt"}
  if (set == 6) {data1 = data[data$gene_info == 'mRNA_ND3',]; outfile = "../../Body/3Results/Alima02.AaAsymmetry.HumanTree.ND3.ExpectedVsObservedAaChanges.txt"}
  if (set == 7) {data1 = data[data$gene_info == 'mRNA_ND2',]; outfile = "../../Body/3Results/Alima02.AaAsymmetry.HumanTree.ND2.ExpectedVsObservedAaChanges.txt"}
  if (set == 8) {data1 = data[data$gene_info == 'mRNA_ND1',]; outfile = "../../Body/3Results/Alima02.AaAsymmetry.HumanTree.ND1.ExpectedVsObservedAaChanges.txt"}
  if (set == 9) {data1 = data[data$gene_info == 'mRNA_ND4L',]; outfile = "../../Body/3Results/Alima02.AaAsymmetry.HumanTree.ND4L.ExpectedVsObservedAaChanges.txt"}
  if (set == 10) {data1 = data[data$gene_info == 'mRNA_COX1',]; outfile = "../../Body/3Results/Alima02.AaAsymmetry.HumanTree.COX1.ExpectedVsObservedAaChanges.txt"}
  if (set == 11) {data1 = data[data$gene_info == 'mRNA_COX2',]; outfile = "../../Body/3Results/Alima02.AaAsymmetry.HumanTree.COX2.ExpectedVsObservedAaChanges.txt"}
  if (set == 12) {data1 = data[data$gene_info == 'mRNA_COX3',]; outfile = "../../Body/3Results/Alima02.AaAsymmetry.HumanTree.COX3.ExpectedVsObservedAaChanges.txt"}
  if (set == 13) {data1 = data[data$gene_info == 'mRNA_ATP6',]; outfile = "../../Body/3Results/Alima02.AaAsymmetry.HumanTree.ATP6.ExpectedVsObservedAaChanges.txt"}
  if (set == 14) {data1 = data[data$gene_info == 'mRNA_ATP8',]; outfile = "../../Body/3Results/Alima02.AaAsymmetry.HumanTree.ATP8.ExpectedVsObservedAaChanges.txt"}
  
FromTo = data.frame(table(data1$FromTo))
names(FromTo) = c('FromAncestralToDerived', 'NumberOfEvents')
nrow(FromTo) # 236, but totally there are 21*21 possibilities
FromTo$FromAncestralToDerived = as.character(FromTo$FromAncestralToDerived)

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
wilcox.test(FromTo5$ExpectedMoreThanOne, mu = 1) # PAPER !!!! p-value = 7.47e-05

FromTo5$TotalSubst = FromTo5$NumberOfExpectedAaSubst + FromTo5$NumberOfUnexpectedAaSubst
Short = FromTo5[!is.na(FromTo5$ExpectedMoreThanOne),]
Short$DummyAhGh = 0
for (i in 1:nrow(Short))
{
  if (Short$NuclSubstLightChainNotation[i] == 'T>C') {Short$DummyAhGh[i] = 1}
  if (Short$NuclSubstLightChainNotation[i] == 'G>A') {Short$DummyAhGh[i] = 0}
}

#cor.test(Short$ExpectedMoreThanOne,Short$GranthamDistance, method = 'spearman') 
#cor.test(Short$TotalSubst,Short$GranthamDistance, method = 'spearman') # a bit negative: the lower Grantham, the higher number of substitutions. good
#cor.test(Short$ExpectedMoreThanOne,Short$TotalSubst, method = 'spearman') 
#summary(lm(Short$ExpectedMoreThanOne ~ Short$TotalSubst*Short$GranthamDistance))
#summary(lm(Short[Short$NuclSubstLightChainNotation == 'T>C',]$ExpectedMoreThanOne ~ Short[Short$NuclSubstLightChainNotation == 'T>C',]$TotalSubst + Short[Short$NuclSubstLightChainNotation == 'T>C',]$GranthamDistance))

## bias is higher in Ah>Gh than Ch>Th!!! why? Ah>Gh is more asymmetric on average???!!! may be yes!? check human global tree and mammalian average piechart
wilcox.test(Short[Short$NuclSubstLightChainNotation == 'G>A',]$ExpectedMoreThanOne,Short[Short$NuclSubstLightChainNotation == 'T>C',]$ExpectedMoreThanOne)
boxplot(Short[Short$NuclSubstLightChainNotation == 'G>A',]$ExpectedMoreThanOne,Short[Short$NuclSubstLightChainNotation == 'T>C',]$ExpectedMoreThanOne, names = c('Ch>Th','Ah>Gh'), ylab = 'expected shift')

Short = Short[,-c(2:6)]
write.table(Short,outfile)
}

dev.off()
