rm(list=ls(all=TRUE)) 

### THIS FILE IS SIMILAR TO Alima04.AsymmetryInCancers.r
###########################################################################################################
##### 1. read expected variants (27) and make them bidirectional (all four transitions): AaPredictions (54)
###########################################################################################################

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

###########################################################################################################
##### 2. read extensively mutated file, take only nons and only driven by FOUR transitions, merge with AaPredictions and make base
###########################################################################################################

MitAnno = read.table("../../Body/2Derived/mtNucAnnotation_MergeRazerV6TurboMaxPro_normMTcode_normND6_withGERP.AllTransitions.Tbss.txt")
MitAnno = MitAnno[MitAnno$AaSub!='',] # 7152
MitAnno = merge(MitAnno,AaPredictions,by="AaSub") 
names(MitAnno)
for (i in 1:nrow(MitAnno))
{
  if (MitAnno$nuc[i] == 'G')   {MitAnno$NucSub[i] = 'G>A'}
  if (MitAnno$nuc[i] == 'A')   {MitAnno$NucSub[i] = 'A>G'}
  if (MitAnno$nuc[i] == 'T')   {MitAnno$NucSub[i] = 'T>C'}
  if (MitAnno$nuc[i] == 'C')   {MitAnno$NucSub[i] = 'C>T'}
}
MitAnno$position = MitAnno$pos
MitAnno$mutation = 0
MitAnno = MitAnno[names(MitAnno) %in% c('position','NucSub','AaSub','mutation')]

###########################################################################################################
##### 3. read human global tree and keep only nons mutations, driven by transitions
###########################################################################################################

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
data$AaSub =  paste(data$ancestral_aa,data$derived_aa,sep='>')
table(data$AaSub)
data$NucSub = paste(data$temp1,data$temp2,sep='>')
table(data$NucSub)
## save only transitions:
nrow(data) # 96216
data=data[data$NucSub %in% c('A>G','G>A','C>T','T>C'),]
nrow(data) # 88644

head(data)
######### how often mutations are the same (direction of aminoacids) as in RefSeq predictions?
DataShort = data[colnames(data) %in% c('ref_pos','NucSub','AaSub')]
DataShort = DataShort[order(DataShort$ref_pos),]  # there are many places with forward and backward substitutions for each site! - how I analyze it?

nrow(DataShort)
DataShort=merge(DataShort,MitAnno, by.x = 'ref_pos',by.y = 'position')
nrow(DataShort)

######### TO DO
### distribution of observed mutations by positions: rare (singletons) should be more biased (mutagenesis only), while common (high derived allele frequeincy) should be under selection
### we need to estimate frequency of a given variant (!!!) not the number of mutations towards it => so I need to use Helix data.. or Ask Kristina to derive more!!!

Pos = data.frame(table(data$ref_pos)); names(Pos) = c('Position','Freq')
Pos = Pos[order(Pos$Freq),]
summary(Pos$Freq) # 
