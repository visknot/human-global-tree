rm(list=ls(all=TRUE)) 

data = read.table("../../Body/1Raw/ComparativeSpeciesTriosFromDima/Pars_Subs.txt", sep = ' ', header = TRUE)
# I need to reformat this table and make it similar to this one: 3Results/Alima04.AsymmetryInCancers.12genes.txt
data$FirstAa = gsub(">(.*)",'',data$Substitution)
data$SecondAa = gsub("(.*)>",'',data$Substitution)

# from one letter to three letter code AA:
A = c("G","A","L","M","F","W","K","Q","E","S","P","V","I","C","Y","H","R","N","D","T","*")
AAA = c("Gly","Ala","Leu","Met","Phe","Trp","Lys","Gln","Glu","Ser","Pro","Val","Ile","Cys","Tyr","His","Arg","Asn","Asp","Thr","Stop")
AA = data.frame(A,AAA)

data = merge(data,AA, by.x = 'FirstAa', by.y = 'A')
data = data[,-c(names(data) == 'FirstAa')]
names(data)[names(data) == 'AAA'] <- 'FirstAa'

data = merge(data,AA, by.x = 'SecondAa', by.y = 'A')
data = data[,-c(names(data) == 'SecondAa')]
names(data)[names(data) == 'AAA'] <- 'SecondAa'
data = data[,-c(names(data) == 'AA')]
data$AaSubst = paste(data$FirstAa,data$SecondAa,sep = '>')
names(data)
data=data[names(data) %in% c('AaSubst','Frequency','Count_Amino')]

#### expectations
exp = read.table("../../Body/2Derived/ExpectedAaSubstitutionDirection20AA.txt", header = TRUE, sep = " ")
first = gsub(">(.*)",'',exp$ExpectedAminoAcidSubstBias)
second = gsub("(.*)>",'',exp$ExpectedAminoAcidSubstBias)
exp$UnexpAaSubst = paste(second,first,sep= '>')
VecOfUnexpAaSubst = exp$UnexpAaSubst
VecOfExpAaSubst = exp$ExpectedAminoAcidSubstBias

DataExp = data[data$AaSubst %in% VecOfExpAaSubst,]
names(DataExp)=c('ExpFreq','ExpAncestor','ExpAaSubst')
for (i in 1:nrow(DataExp)) {DataExp$AaPairId[i] = paste(sort(unlist(strsplit(DataExp$ExpAaSubst[i],'>'))),collapse = '>') }

DataUnexp = data[data$AaSubst %in% VecOfUnexpAaSubst,]
names(DataUnexp)=c('UnexpFreq','UnexpAncestor','UnexpAaSubst')
for (i in 1:nrow(DataUnexp)) {DataUnexp$AaPairId[i] = paste(sort(unlist(strsplit(DataUnexp$UnexpAaSubst[i],'>'))),collapse = '>') }

data = merge(DataExp,DataUnexp, by = 'AaPairId', all = TRUE)
data[is.na(data)] <-0

#### run fisher test:
for (i in 1:nrow(data))
{
X = rbind(c(data$ExpFreq[i],data$ExpAncestor[i]),c(data$UnexpFreq[i],data$UnexpAncestor[i]))
fisher = fisher.test(X)
data$FisherP[i] = fisher$p.value
data$FisherOdds[i] = as.numeric(fisher$estimate)
data$FisherConf1[i] = fisher$conf.int[1]
data$FisherConf2[i] = fisher$conf.int[2]
}

data = data[order(data$FisherOdds),]
data = data[data$UnexpAaSubst != 0 & data$ExpAaSubst != 0,]

summary(data$FisherOdds)
wilcox.test(data$FisherOdds, mu = 1) # p = 0.79
hist(log2(data$FisherOdds), breaks = 20) # p = 0.79
abline(v=0, col = 'red', lwd = 2)

#data$ExpFreq[i]/data$ExpAncestor[i] # 2.15    2.15/1.66 = 1.29
#data$UnexpFreq[i]/data$UnexpAncestor[i] # 1.66

#### strange results - no trend, but really different fluctuations (in both directions). probably, something is wrong with normalization?
### https://www.ncbi.nlm.nih.gov/nuccore/NC_012920
CytbHumanRefSeq = "MTPMRKTNPLMKLINHSFIDLPTPSNISAWWNFGSLLGACLILQITTGLFLAMHYSPDASTAFSSIAHITRDVNYGWIIRYLHANGASMFFICLFLHIGRGLYYGSFLYSETWNIGIILLLATMATAFMGYVLPWGQMSFWGATVITNLLSAIPYIGTDLVQWIWGGYSVDSPTLTRFFTFHFILPFIIAALATLHLLFLHETGSNNPLGITSHSDKITFHPYYTIKDALGLLLFLLSLMTLTLFSPDLLGDPDNYTLANPLNTPPHIKPEWYFLFAYTILRSVPNKLGGVLALLLSILILAMIPILHMSKQQSMMFRPLSQSLYWLLAADLLILTWIGGQPVSYPFTIIGQVASVLYFTTILILMPTISLIENKMLKWA"
CytbAaFreqHumanRefSeq = data.frame(table(unlist(strsplit(CytbHumanRefSeq,split=''))))
names(CytbAaFreqHumanRefSeq)=c('AA','FreqRefSeq')

DataOnceMore = read.table("../../Body/1Raw/ComparativeSpeciesTriosFromDima/Pars_Subs.txt", sep = ' ', header = TRUE)
names(DataOnceMore)
DataOnceMore=unique(DataOnceMore[names(DataOnceMore) %in% c('AA','Count_Amino')])

Cytb = merge(CytbAaFreqHumanRefSeq,DataOnceMore, all = TRUE)  # Ile is absent !!! 
cor.test(Cytb$FreqRefSeq,Cytb$Count_Amino, method = 'spearman') # otherwise it is a good correlation
plot(Cytb$FreqRefSeq,Cytb$Count_Amino)

