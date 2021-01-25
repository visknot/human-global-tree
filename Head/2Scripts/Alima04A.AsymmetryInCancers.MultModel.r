rm(list=ls(all=TRUE))

##### 1. read expected variants and modify a bit 

##### 2. read extensively mutated file, take only nons and only 


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
data$FromTo = paste(data$ancestral_aa,data$derived_aa,sep = ">")
VecOfNames = c("Tier2","chrom","position","FromTo","ancestral_aa","derived_aa")
data = data[colnames(data) %in% VecOfNames]

### read mtDNA annotation by Victor:
mt = read.table("../../Body/1Raw/mtNucAnnotation_MergeRazerV6TurboMaxPro_normMTcode_normND6_withGERP.csv", header = TRUE, sep = ';')

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

##### read expectations, add modify them to dataset AaPredictions
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


