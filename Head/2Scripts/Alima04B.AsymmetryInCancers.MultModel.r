rm(list=ls(all=TRUE))

##### 1. read expected variants (27) and make them bidirectional (all four transitions): AaPredictions (54)

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

##### 2. read extensively mutated file, take only nons and only driven by FOUR transitions and make base

MitAnno = read.table("../../Body/2Derived/mtNucAnnotation_MergeRazerV6TurboMaxPro_normMTcode_normND6_withGERP.AllTransitions.Tbss.txt")
MitAnno = MitAnno[MitAnno$AaSub!='',] # 7152
MitAnno = merge(MitAnno,AaPredictions,by="AaSub") 
# THERE ARE 7152 POTENTIAL NONSYNONYMOUS SUBSTITUTIONS BECAUSE OF 4 TRANSITIONS
#MitAnno = merge(MitAnno,AaPredictions,by="AaSub",all.x=TRUE) # 6504 => 7597-6504 (1093 substitutions were not in MitAnno - which ones)
#table(MitAnno$AaSub)
#MitAnnoNa = MitAnno[is.na(MitAnno$ExpDummy),] # zero

#mut = mut[mut$SubstType == 'missense_variant',]
#mut = mut[,c(2,4,5,8,9)]
#names(mut)=c('position','Ref','Alt','Annotation','SubstType')
# save only nonsynonymous, overlapped with 27*2 above
#mut$NucSub = paste(mut$Ref,mut$Alt,sep = '>')
#table(mut$NucSub)
#mut = mut[mut$NucSub %in% c('C>T','T>C','G>A','A>G'),]
# save only transitions

# ReturnEleven <- function(x) 
# { # temp = 'p.Asp139Asn'
#  temp = unlist(strsplit(x, split = "\\|"))[11]; 
#  Aa1 = paste(unlist(strsplit(temp,split=''))[c(3,4,5)], collapse = '')
#  Length = length(unlist(strsplit(temp,split='')));
#  Aa2 = paste(unlist(strsplit(temp,split=''))[c(Length-2,Length-1,Length)], collapse = '')
#  AaSubs = paste(Aa1,Aa2,sep='>')
#  }
# mut$AaSub = apply(as.matrix(mut$Annotation),1,ReturnEleven)
# table(mut$AaSub)
# MutShort = mut[names(mut) %in% c('position','NucSub','AaSub')]
# MutShort$mutation = 0 # it means mutation which might happen

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

##### 3. READ CANCER DATA AND REMOVE EVERYTHING EXCEPT 27*2 
data = read.table("../../Body/1Raw/CancerDataFromCampbell/mtDNA_snv_Oct2016.txt", header = TRUE, sep = '\t')
for (i in 1:nrow(data)) {data$AaChanges[i] = unlist(strsplit(data$Annot[i],','))[5]}
data$ancestral_aa = gsub("\\d(.*)",'',data$AaChanges)
data$derived_aa = gsub("(.*)\\d",'',data$AaChanges)
data=data[grepl('nsSNP', data$Annot),]
data$NucSub = paste(data$ref,data$var, sep = '>')
data$AaSub = paste(data$ancestral_aa,data$derived_aa, sep = '>')
table(data$AaSub)

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
nrow(data) # 3805
data = data[data$AaSub %in% AaPredictions$AaSub,] # delete all substitutions due to transversions and double mutations
nrow(data) # 3540
data = data[data$NucSub %in% c('C>T','T>C','A>G','G>A'),] # delete all substitutions due to transversions and double mutations
nrow(data) # 3528
# THERE ARE 3528 NONSYNONYMOUS SUBSTITUTIONS, DUE TO TRANSITIONS, OBSERVED IN CANCERS
length(unique(data$position)) # 1940
# THERE ARE 1940 NONSYNONYMOUS SUBSTITUTIONS, DUE TO TRANSITIONS, OBSERVED IN CANCERS

### rbind cancer with MitAnno and aggregate

data = rbind(data,MitAnno)
agg = aggregate(data$mutation, by = list(data$position,data$NucSub,data$AaSub), FUN = sum)
names(agg) = c('position','NucSub','AaSub','mutation')
agg = agg[order(agg$position),]
nrow(agg) # 7152
length(unique(agg$position)) # 7152
#DuplPositions = agg[duplicated(agg$position),]$position; length(DuplPositions)
#DuplPositions # 12766  !!!!!!!!!!!!!!!!!!!!!!!!!!! should be ZERO
#Problems = agg[agg$position %in% DuplPositions,]
#agg = agg[!agg$position %in% DuplPositions,]   # ZAPLATKA!!!!!!!!!!!!!

#### merge agg with AaPredictions
nrow(agg) # 71522
agg = merge(agg,AaPredictions, by = 'AaSub')
nrow(agg) # 71522

#### merge agg with Victor's annotation
MitAnno = read.table("../../Body/2Derived/mtNucAnnotation_MergeRazerV6TurboMaxPro_normMTcode_normND6_withGERP.AllTransitions.Tbss.txt")
nrow(agg) # 71522
agg = merge(agg,MitAnno, by.x = 'position', by.y = 'pos')
nrow(agg) # 71522

table(agg$role)
# derive DummyNd6!!!!!!!!!!!!!!!
agg$DummyNd6 = 0; 
for (i in 1:nrow(agg)) {if (agg$role[i] == 'mRNA_ND6') {agg$DummyNd6[i] = 1}}
names(agg)

################################
########### analysis
################################

summary(glm(mutation ~ ExpDummyFromLoserToGainer*DummyNd6 + ChToThDummy + scale(TimeBeingSingleStranded) + scale(GranthamDistance), data = agg, family = 'poisson'))
summary(glm(mutation ~ (ExpDummyFromLoserToGainer + ChToThDummy)*DummyNd6 + scale(TimeBeingSingleStranded) + scale(GranthamDistance), data = agg, family = 'poisson')) # final
#Call:
#  glm(formula = mutation ~ (ExpDummyFromLoserToGainer + ChToThDummy) * 
#        DummyNd6 + scale(TimeBeingSingleStranded) + scale(GranthamDistance), 
#      family = "poisson", data = agg)
#
#Deviance Residuals: 
#  Min       1Q   Median       3Q      Max  
#-2.2277  -0.6834  -0.4180  -0.2396  11.5692  
#
#Coefficients:
#  Estimate Std. Error z value Pr(>|z|)    
#(Intercept)                        -3.26965    0.06529 -50.080  < 2e-16 ***
#  ExpDummyFromLoserToGainer           2.70892    0.06217  43.576  < 2e-16 ***
#  ChToThDummy1                        1.11518    0.03550  31.410  < 2e-16 ***
#  DummyNd6                            3.64382    0.15356  23.728  < 2e-16 ***
#  scale(TimeBeingSingleStranded)      0.17637    0.01721  10.250  < 2e-16 ***
#  scale(GranthamDistance)            -0.06033    0.01564  -3.858 0.000114 ***
#  ExpDummyFromLoserToGainer:DummyNd6 -5.74148    0.29385 -19.539  < 2e-16 ***
#  ChToThDummy1:DummyNd6              -1.95682    0.19051 -10.272  < 2e-16 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for poisson family taken to be 1)
#
#Null deviance: 10720.5  on 7151  degrees of freedom
#Residual deviance:  6369.7  on 7144  degrees of freedom
#AIC: 10994
#
#Number of Fisher Scoring iterations: 6

summary(glm(mutation ~ ExpDummyFromLoserToGainer*DummyNd6 + scale(TimeBeingSingleStranded) + scale(GranthamDistance), data = agg, family = 'poisson'))

### Quadrant 3 is a bit more mutagenic... why? low Grantham? 
summary(glm(mutation ~ (ExpDummyFromLoserToGainer + ChToThDummy)*DummyNd6 + scale(TimeBeingSingleStranded) + scale(GranthamDistance) + QuadrantOne + QuadrantTwo + QuadrantThree + QuadrantFour, data = agg, family = 'poisson')) # final
summary(glm(mutation ~ (ExpDummyFromLoserToGainer + ChToThDummy)*DummyNd6 + scale(TimeBeingSingleStranded) + scale(GranthamDistance) + QuadrantThree + QuadrantFour, data = agg, family = 'poisson')) # final
summary(glm(mutation ~ (ExpDummyFromLoserToGainer + ChToThDummy)*DummyNd6 + scale(TimeBeingSingleStranded) + scale(GranthamDistance) + QuadrantThree, data = agg, family = 'poisson')) # final

