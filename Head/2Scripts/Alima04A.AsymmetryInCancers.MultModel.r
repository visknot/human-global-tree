rm(list=ls(all=TRUE))

##### 1. read expected variants and modify them to dataset AaPredictions

exp = read.table("../../Body/2Derived/ExpectedAaSubstitutionDirection20AA.txt", header = TRUE, sep = " ")
for (i in 1:nrow(exp)) {exp$AaPairId[i] = paste(sort(unlist(strsplit(exp$ExpectedAminoAcidSubstBias[i],'>'))),collapse = '>') }
exp$ExpDummy = 0
exp$AhToGhDummy = gsub('T>C',1,exp$NuclSubstLightChainNotation)
exp$AhToGhDummy = gsub('G>A',0,exp$AhToGhDummy)
exp=exp[,-c(2,5)]
names(exp)[2] = c('AASub')
names(exp)

TO = gsub("(.*)>",'',exp$AASub); FROM = gsub(">(.*)",'',exp$AASub); 
unexp = paste(TO,FROM,sep = '>')

unexp = data.frame(unexp,exp$GranthamDistance,exp$NuclPosInCodon,exp$AhToGhDummy)
unexp$ExpDummy = 1
names(unexp)=c('AASub','GranthamDistance','NuclPosInCodon','AhToGhDummy','ExpDummy')
AaPredictions = rbind(exp,unexp)

##### 2. read extensively mutated file, take only nons and only driven by two most common transitions and make base

ReturnSecond <- function(x) {unlist(strsplit(x, split = "\\|"))[2]}

mut = read.table("../../Body/1Raw/FromCovidProject/HumanMtDnaRefSeq.fasta.ExtensiveMut.vcf.ann", header = FALSE, sep = '\t')
mut$SubstType = apply(as.matrix(mut$V8),1,ReturnSecond)
table(mut$SubstType)
#initiator_codon_variant                            missense_variant      missense_variant&splice_region_variant 
#12                                       23410                                          30 
#splice_region_variant&stop_retained_variant    splice_region_variant&synonymous_variant                                  start_lost 
#5                                          39                                          58 
#stop_gained           stop_gained&splice_region_variant                                   stop_lost 
#1017                                           3                                         744 
#stop_lost&splice_region_variant                       stop_retained_variant                          synonymous_variant 
#36                                          92                                        8577 
#upstream_gene_variant 
#15684 

mut = mut[mut$SubstType == 'missense_variant',]
mut = mut[,c(2,4,5,8,9)]
names(mut)=c('position','Ref','Alt','Annotation','SubstType')
# save only nonsynonymous, overlapped with 27*2 above
mut$NucSub = paste(mut$Ref,mut$Alt,sep = '>')
table(mut$NucSub)
mut = mut[mut$NucSub %in% c('C>T','T>C','G>A','A>G'),]
# save only transitions

ReturnEleven <- function(x) 
{ # temp = 'p.Asp139Asn'
  temp = unlist(strsplit(x, split = "\\|"))[11]; 
  Aa1 = paste(unlist(strsplit(temp,split=''))[c(3,4,5)], collapse = '')
  Length = length(unlist(strsplit(temp,split='')));
  Aa2 = paste(unlist(strsplit(temp,split=''))[c(Length-2,Length-1,Length)], collapse = '')
  AaSubs = paste(Aa1,Aa2,sep='>')
  }
mut$AaSub = apply(as.matrix(mut$Annotation),1,ReturnEleven)
table(mut$AaSub)
MutShort = mut[names(mut) %in% c('position','NucSub','AaSub')]
MutShort$mutation = 0 # it means mutation which might happen

##### 3. 
data = read.table("../../Body/1Raw/CancerDataFromCampbell/mtDNA_snv_Oct2016.txt", header = TRUE, sep = '\t')
for (i in 1:nrow(data)) {data$AaChanges[i] = unlist(strsplit(data$Annot[i],','))[5]}
data$ancestral_aa = gsub("\\d(.*)",'',data$AaChanges)
data$derived_aa = gsub("(.*)\\d",'',data$AaChanges)
data=data[grepl('nsSNP', data$Annot),]
data$NucSub = paste(data$ref,data$var, sep = '>')
data$AASub = paste(data$ancestral_aa,data$derived_aa, sep = '>')
table(data$AASub)

# from one letter to three letter code AA:
A = c("G","A","L","M","F","W","K","Q","E","S","P","V","I","C","Y","H","R","N","D","T","*")
AAA = c("Gly","Ala","Leu","Met","Phe","Trp","Lys","Gln","Glu","Ser","Pro","Val","Ile","Cys","Tyr","His","Arg","Asn","Asp","Thr","Stop")
AA = data.frame(A,AAA)
data = merge(data,AA, by.x = 'ancestral_aa',by.y='A'); data$ancestral_aa = data$AAA
data = merge(data,AA, by.x = 'derived_aa',by.y='A'); data$derived_aa = data$AAA.y
data$AaSub = paste(data$ancestral_aa,data$derived_aa,sep = ">")
VecOfNames = c("position","NucSub","AaSub")
data = data[colnames(data) %in% VecOfNames]
names(data) # "position" "NucSub"   "AaSub" 
data$mutation = 1 # it means real mutation

### rbind MutShort and data and aggregate

data = rbind(data,MutShort)
agg = aggregate(data$mutation, by = list(data$position,data$NucSub,data$AaSub), FUN = sum)
names(agg) = c('position','NucSub','AaSub','mutation')
agg = agg[order(agg$position),]

#### КОЕ ГДЕ ПАРАША!!! (3308 например) => либо перегнать SnpEff, либо промутировать руками (на базе скрипта Вити)


### read mtDNA annotation by Victor:

library('Biostrings')

mt = read.table("../../Body/1Raw/mtNucAnnotation_MergeRazerV6TurboMaxPro_normMTcode_normND6_withGERP.csv", header = TRUE, sep = ';')
table(mt$role)
ProtCodGenes = c("mRNA_ATP6","mRNA_ATP6&COX3","mRNA_ATP8","mRNA_ATP8&ATP6","mRNA_COX1","mRNA_COX2","mRNA_COX3","mRNA_CYTB","mRNA_ND1","mRNA_ND2","mRNA_ND3","mRNA_ND4","mRNA_ND4L","mRNA_ND4L&ND4","mRNA_ND5","mRNA_ND6")
length(ProtCodGenes) # 13 genes + 3 overlaps
mt = mt[mt$role %in% ProtCodGenes,]

### for ND6 reverse complement codons
mt$RnaCodon = mt$codon
for (i in 1:nrow(mt))
{ # i = 12
  if (mt$role[i] == 'mRNA_ND6')
  {
  TempCodon = DNAString(x=mt$codon[i], start=1, nchar=NA);  
  mt$RnaCodon[i] = as.character(reverseComplement(TempCodon))
  }
}

### walk nucleotide by nucleotide, if it is G (Ch) or T (Ah) I mutate them accordingly and see - what happens (nons, syn and if nons - which aminoacids)

mt$MutatedCodon = ''
for (i in 1:nrow(mt))
{
  if (mt$nuc[i] == 'T' & mt$nnuc[i] == 1)  {mt$MutatedCodon[i] = paste('C',unlist(strsplit(mt$codon[i],''))[2],unlist(strsplit(mt$codon[i],''))[3], sep = '')}
  if (mt$nuc[i] == 'G' & mt$nnuc[i] == 1)  {mt$MutatedCodon[i] = paste('A',unlist(strsplit(mt$codon[i],''))[2],unlist(strsplit(mt$codon[i],''))[3], sep = '')}
  if (mt$nuc[i] == 'T' & mt$nnuc[i] == 2)  {mt$MutatedCodon[i] = paste(unlist(strsplit(mt$codon[i],''))[1],'C',unlist(strsplit(mt$codon[i],''))[3], sep = '')}
  if (mt$nuc[i] == 'G' & mt$nnuc[i] == 2)  {mt$MutatedCodon[i] = paste(unlist(strsplit(mt$codon[i],''))[1],'A',unlist(strsplit(mt$codon[i],''))[3], sep = '')}
}

MitCode = getGeneticCode(id_or_name2="2", full.search=FALSE, as.data.frame=FALSE)
mt$MutatedAa = ''
for (i in 1:nrow(mt))
{ # i = 1
  if (mt$MutatedCodon[i] != '') 
  {
  codon = DNAString(x=mt$MutatedCodon[i], start=1, nchar=NA); 
  AaOneLetter = AAString(translate(codon,MitCode)); 
  mt$MutatedAa[i] = as.character(AMINO_ACID_CODE[strsplit(as.character(AaOneLetter), NULL)[[1]]]);
  }
}
mt$MutatedAa[is.na(mt$MutatedAa)] <- 'Stop'


# if ND6 = reverse compliment each codon: reverse(data$codon[i]) !!!!!!!!!!!!!!!!!!!!!

### merge data with mt and rename Leu and Ser according to their codons [I need to revert codons in case of ND6 before this]
data = merge(data,mt,by.x = 'position', by.y = 'pos')
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

table(data$ancestral_aa)  # till now there are some 'L' and 'S' => WHY!!!!???? ND6 - parasha!!! victor!!!!
temp = data[data$ancestral_aa == 'L' | data$ancestral_aa == 'S',]

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
### preparation of a matrix 
####################################

##### 
nrow(data) # 3802
data = merge(data,AaPredictions, by.x = 'FromTo', by.y = 'AASub') # 
nrow(data) # 3540 == majority

for (i in 1:nrow(data)) {  if (data$role[i] == 'mRNA_ND6')  {data$NotNd6Dummy[i] = 0};   if (data$role[i] != 'mRNA_ND6')  {data$NotNd6Dummy[i] = 1}}
table(data$NotNd6Dummy)

for (i in 1:nrow(data))
{
if (data$position[i] > 5798) {data$Tbss[i] = (data$position[i]-5798)*2}
if (data$position[i] <= 5798) {data$Tbss[i] = (16569-5798 + data$position[i])}
}
summary(data$Tbss)

#################################
########### analysis
################################

table(data$ExpDummy) # 390 versus 3150!!!!!! logistic regression works better if units are rare -> recode!!!!

summary(glm(ExpDummy ~ GranthamDistance + AhToGhDummy +  Tbss, data = data[data$NotNd6Dummy == 1,]))
summary(glm(ExpDummy ~ GranthamDistance + Tbss, data = data[data$NotNd6Dummy == 1,]))
summary(glm(ExpDummy ~ Tbss, data = data[data$NotNd6Dummy == 1,]))

#### probably I need to consider EACH aminoacid in mt genome and code it as mutated or not and if mutated => which direction (expected or unexpected) and how many times...


for (i in 1:nrow(data))
{
if ()  
}


