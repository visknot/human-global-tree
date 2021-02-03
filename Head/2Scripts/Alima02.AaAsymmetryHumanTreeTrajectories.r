rm(list=ls(all=TRUE)) 

##########################################################
############# 1: READ HUMAN TREE AND KEEP ONLY NONSYN MUTS
##########################################################

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

#### From Leu to LeuCT and LeuTT, from Ser to SerTC and SerAG etc SPLIT ALL NECESSARY CODONS.... 
#####LERA....

for (i in 1:nrow(data))
{ # i = 1
  if (data$ancestral_aa[i] == 'Leu')
  {
    if (tolower(data$ancestral_codon[i]) %in% c('ctt','ctc','cta','ctg')) {data$ancestral_aa[i] = 'LeuCT'}
  }
}
table(data$ancestral_aa)
table(data$derived_aa)
data$FromTo = paste(data$ancestral_aa,data$derived_aa,sep = '>')

##########################################################
############# 2: MAKE AA PREDICTIONS (LERA - IMPROVE IT, EXPAND IT?)
##########################################################

exp = read.table("../../Body/2Derived/ExpectedAaSubstitutionDirection20AA.txt", header = TRUE, sep = " ")
for (i in 1:nrow(exp)) {exp$AaPairId[i] = paste(sort(unlist(strsplit(exp$ExpectedAminoAcidSubstBias[i],'>'))),collapse = '>') }
exp$ExpDummyFromLoserToGainer = 1
exp$ChToThDummy = gsub('G>A',1,exp$NuclSubstLightChainNotation)
exp$ChToThDummy = gsub('T>C',0,exp$ChToThDummy)
exp=exp[,-c(2,5)]
names(exp)[2] = c('AaSub')
names(exp)

TO = gsub("(.*)>",'',exp$AaSub); FROM = gsub(">(.*)",'',exp$AaSub); 
unexp = paste(TO,FROM,sep = '>')

unexp = data.frame(unexp,exp$GranthamDistance,exp$NuclPosInCodon,exp$ChToThDummy)
unexp$ExpDummyFromLoserToGainer = 0
names(unexp)=c('AaSub','GranthamDistance','NuclPosInCodon','ChToThDummy','ExpDummyFromLoserToGainer')
AaPredictions = rbind(exp,unexp)
AaPredictions[is.na(AaPredictions$GranthamDistance),]$GranthamDistance <-300 # make from stop or to stop equals 300 (max)

AaPredictions$QuadrantOne = 0; AaPredictions$QuadrantTwo = 0; AaPredictions$QuadrantThree = 0;  AaPredictions$QuadrantFour = 0; 
for (i in 1:nrow(AaPredictions))
{ # i = 1
  if (grepl("Gln|His|Tyr|Cys|Trp|Arg",AaPredictions$AaSub[i])) {AaPredictions$QuadrantOne[i] = 1}
  if (grepl("Phe|Leu|Pro",AaPredictions$AaSub[i])) {AaPredictions$QuadrantTwo[i] = 1}
  if (grepl("Ile|Met|Val|Thr|Ala",AaPredictions$AaSub[i])) {AaPredictions$QuadrantThree[i] = 1}
  if (grepl("Asn|Lys|Asp|Glu|Gly",AaPredictions$AaSub[i])) {AaPredictions$QuadrantFour[i] = 1}
}


##########################################################
############# 3: TAKE CODON USAGE OF THE HUMAN MTDNA FROM VICTOR'S FILE (FOR NORMALIZATION OF THE NUMBER OF OBSERVED SUBSTITUTIONS) - LERA - DO IT FOR CODONS...?
##########################################################

MitAnno = read.table("../../Body/2Derived/mtNucAnnotation_MergeRazerV6TurboMaxPro_normMTcode_normND6_withGERP.AllTransitions.Tbss.txt")
MitAnno = MitAnno[MitAnno$AaSub!='',] # 7152
# THERE ARE 7152 POTENTIAL NONSYNONYMOUS SUBSTITUTIONS BECAUSE OF 4 TRANSITIONS



