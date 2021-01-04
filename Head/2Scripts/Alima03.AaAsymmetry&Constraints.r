rm(list=ls(all=TRUE)) 

#indata <- read.csv("/home/alima/arrr/fulltreeCodons.csv", header = TRUE, sep = ";") #читаю файл 
data = read.csv("../../Body/2Derived/fulltreeCodons.csv", header = TRUE, sep = ";")

#### FILTER OUT: 
data <- subset(data, synonymous == "non-synonymous" & derived_aa != "Ambiguous" & derived_aa != "Asn/Asp" & derived_aa != "Gln/Glu" & derived_aa != "Leu/Ile" & gene_info != "mRNA_ND6")  #убираю лишнее
table(data$derived_aa)

## frequencies of changes in a given position versus AaSymmetry 

temp = data[(data$ancestral_aa == 'Leu' & data$derived_aa == 'Pro' ) | (data$ancestral_aa == 'Pro' & data$derived_aa == 'Leu') ,]
temp$AaSubst = paste(temp$ancestral_aa,temp$derived_aa, sep = '_'); nrow(temp)

EvolConstOfPositions = data.frame(table(temp$position)); 
names(EvolConstOfPositions)=c('position','HowManySubstitutions')
EvolConstOfPositions = EvolConstOfPositions[order(-EvolConstOfPositions$HowManySubstitutions),]
EvolConstOfPositionsVec1 = EvolConstOfPositions[EvolConstOfPositions$HowManySubstitutions <= quantile(EvolConstOfPositions$HowManySubstitutions,0.9),]$position
EvolConstOfPositionsVec2 = EvolConstOfPositions[EvolConstOfPositions$HowManySubstitutions > quantile(EvolConstOfPositions$HowManySubstitutions,0.9),]$position

length(EvolConstOfPositionsVec1)
length(EvolConstOfPositionsVec2)

temp1 = temp[temp$position %in% EvolConstOfPositionsVec1,]; nrow(temp1) # 158
temp2 = temp[temp$position %in% EvolConstOfPositionsVec2,]; nrow(temp2) # 2781

table(temp1$AaSubst)
table(temp2$AaSubst)



TempUnderSelection =   temp[temp$position %in% PositionsUnderSelection,]; nrow(TempUnderSelection) # 601
table(TempUnderSelection$AaSubst)   # 316 / 285  = 1.108;   491/460  = 1.06; 680/623 = 1.09; 87/71 = 1.22; 146/127 = 1/149; 
table(TempUnderMutagenesis$AaSubst) #  1230 / 1108 = 1.110; 1055/933 = 1.13; 866 / 770 = 1.12;  1459  /  1322 =  1.1; 1400/1266 = 1.105






## ADD FILTER OF BACKGROUND (THE SAME NEIGHBOR NUCLEOTIDES)

head(data)

##### DERIVE SUBSTITUTION MATRIX FromTo
data$FromTo = paste(data$ancestral_aa,data$derived_aa,sep = '>')
FromTo = data.frame(table(data$FromTo))
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

##### DERIVE UNIQUE IDENTIFIER FOR EACH PAIR OF AMINOACIDS:
for (i in 1:nrow(FromTo2)) {FromTo2$AaPairId[i] = paste(sort(unlist(strsplit(FromTo2$FromAncestralToDerived[i],'>'))),collapse = '>') }

##### from 420 to 210 rows and NumberOfEvents1To2
A = FromTo2[FromTo2$FromAncestralToDerived == FromTo2$AaPairId,]; A= A[c(1,2,5)]; names(A)=c('FromAncestralToDerived1','NumberOfEvents1','AaPairId')
B = FromTo2[FromTo2$FromAncestralToDerived != FromTo2$AaPairId,]; B= B[c(1,2,5)]; names(B)=c('FromAncestralToDerived2','NumberOfEvents2','AaPairId')
FromTo3 = merge(A,B)
FromTo3$NumberOfEvents1To2 = FromTo3$NumberOfEvents1 / FromTo3$NumberOfEvents2 

##### read expectations:
exp = read.table("../../Body/2Derived/ExpectedAaSubstitutionDirection.txt", header = TRUE, sep = " ")
for (i in 1:nrow(exp)) {exp$AaPairId[i] = paste(sort(unlist(strsplit(exp$ExpectedAminoAcidSubstBias[i],'>'))),collapse = '>') }
exp=exp[c(3,4)]
FromTo4 = merge(FromTo3,exp, by = 'AaPairId', all.x = TRUE)

##### 
FromTo4$ExpectedMoreThanOne=NA
for (i in 1:nrow(FromTo4))
{ # i = 17
  if (!is.na(FromTo4$ExpectedAminoAcidSubstBias[i]) & FromTo4$ExpectedAminoAcidSubstBias[i] == FromTo4$FromAncestralToDerived1[i]) {FromTo4$ExpectedMoreThanOne[i] = FromTo4$NumberOfEvents1[i] / FromTo4$NumberOfEvents2[i]}
  if (!is.na(FromTo4$ExpectedAminoAcidSubstBias[i]) & FromTo4$ExpectedAminoAcidSubstBias[i] == FromTo4$FromAncestralToDerived2[i]) {FromTo4$ExpectedMoreThanOne[i] = FromTo4$NumberOfEvents2[i] / FromTo4$NumberOfEvents1[i]}
}

FromTo4 = FromTo4[order(FromTo4$ExpectedMoreThanOne),]


summary(FromTo4$ExpectedMoreThanOne)
hist(FromTo4$ExpectedMoreThanOne,breaks = 100)
wilcox.test(FromTo4$ExpectedMoreThanOne, mu = 1)

FromTo4$TotalSubst = FromTo4$NumberOfEvents1 + FromTo4$NumberOfEvents2



write.table(FromTo4,"../../Body/3Results/Alima02.AaAsymmetry.ExpectedVsObservedAaChanges.txt")

###### controls and permutations  (how to kill signal and derive real p-value)?

