rm(list=ls(all=TRUE)) 

##################################################
#### 1. READ MITOMAP AND FILTER AA SUBSTITUTIONS TO WORK WITH: 
##################################################

data = read.table("../../Body/1Raw/PathogenicHumanMtDnaMutations/MitochondrialDeseasesUpdate.xlsx.txt", header = TRUE, sep = "\t")

## filter - only single nucleotide substitutions (like 'A-C')
table(data$Nucleotide.Change)
data$TempLength = 0
for (i in 1:nrow(data)) {data$TempLength[i] = length(unlist(strsplit(data$Nucleotide.Change[i],'')))}
data$Nucleotide.Change = gsub('-','>',data$Nucleotide.Change)
table(data$Nucleotide.Change)
nrow(data)
data = data[data$TempLength == 3,]
nrow(data)

## filter - take only transitions:
nrow(data)
data = data[data$Nucleotide.Change %in% c('A>G','C>T','G>A','T>C'),]
nrow(data)

## filter - only substitutions within protein-coding genes:
data = data[data$Amino.Acid.Change != 'noncoding',]
data = data[data$Amino.Acid.Change != 'HV1',]     # control region
data = data[!grepl('ATP8',data$Amino.Acid.Change),]  # overlap between ATP8 and ATP6 - not many, too complex to deeg
nrow(data)
data$Amino.Acid.Change = gsub('-','>',data$Amino.Acid.Change)

## read MitAnno and merge with data
MitAnno = read.table("../../Body/2Derived/mtNucAnnotation_MergeRazerV6TurboMaxPro_normMTcode_normND6_withGERP.AllTransitions.Tbss.txt")
MitAnno = MitAnno[MitAnno$AaSub!='',] # 7152
nrow(data) # 336
data = merge(data,MitAnno,by.x="Position",by.y='pos') # 7152 
nrow(data) # 306 - loose some variants with synonymous substitutions (first and second amino acids are the same)

table(data$Locus) # 26 mutations from Nd6 and the rest from TwelveGenes

### make matrix of AaSubstitutions for twelve genes:
FromTo12 = data.frame(table(data[data$Locus != 'MT-ND6',]$AaSub))
names(FromTo12) = c('AaSub', 'NumberOfEvents')

FromToNd6 = data.frame(table(data[data$Locus == 'MT-ND6',]$AaSub))
names(FromToNd6) = c('AaSub', 'NumberOfEvents')


##################################################
##### 2. read expected variants (27) and make them bidirectional (all four transitions): AaPredictions (54)
##################################################

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
for (i in 1:nrow(AaPredictions)) {AaPredictions$AaPairId[i] = paste(sort(unlist(strsplit(AaPredictions$AaSub[i],'>'))),collapse = '>') }

##################################################
##### 3. merge FromTo(12 and Nd6) with AaPredictions and write table 3Results/Alima02B.AaAsymmetryMitoMap.TwelveGenes.txt
##################################################

### FromTo12

FromTo12 = merge(FromTo12,AaPredictions, by = 'AaSub')
FromTo12Exp = FromTo12[FromTo12$ExpDummyFromLoserToGainer == 1,];
names(FromTo12Exp); FromTo12Exp = FromTo12Exp[,c(1,2,11)]; names(FromTo12Exp)[2]=c('NumberOfExpectedAaSubst')
FromTo12Unexp = FromTo12[FromTo12$ExpDummyFromLoserToGainer == 0,];
names(FromTo12Unexp); FromTo12Unexp = FromTo12Unexp[,c(2,11)]; names(FromTo12Unexp)[1]=c('NumberOfUnexpectedAaSubst')

MitoMap = merge(FromTo12Exp,FromTo12Unexp, by = 'AaPairId', all = TRUE)
MitoMap = MitoMap[!is.na(MitoMap$AaSub),]
MitoMap[is.na(MitoMap)]<-0

write.table(MitoMap,"../../Body/3Results/Alima02B.AaAsymmetryMitoMap.TwelveGenes.txt")

### FromToNd6

FromTo = merge(FromToNd6,AaPredictions, by = 'AaSub')
FromToExp = FromTo[FromTo$ExpDummyFromLoserToGainer == 1,];
names(FromToExp); FromToExp = FromToExp[,c(1,2,11)]; names(FromToExp)[2]=c('NumberOfExpectedAaSubst')
FromToUnexp = FromTo[FromTo$ExpDummyFromLoserToGainer == 0,];
names(FromToUnexp); FromToUnexp = FromToUnexp[,c(2,11)]; names(FromToUnexp)[1]=c('NumberOfUnexpectedAaSubst')

MitoMap = merge(FromToExp,FromToUnexp, by = 'AaPairId', all = TRUE)
MitoMap = MitoMap[!is.na(MitoMap$AaSub),]
MitoMap[is.na(MitoMap)]<-0

write.table(MitoMap,"../../Body/3Results/Alima02B.AaAsymmetryMitoMap.Nd6.txt")


