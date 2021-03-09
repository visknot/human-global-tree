rm(list=ls(all=TRUE)) 

data = read.table("../../Body/1Raw/VictorSimulations/pseudogenomes1/generations.csv", sep = ';', header = TRUE, quote = '')

#### equilibrium was reached in 100 000 generations (stop codons!?)

#### go to last generation and switch to percentage of the codon usage:
data = data[data$generations == 1000000,]
data = data[,-c(1)]
sum(data[c(1:64)]) # 6400
data = data/6400 # 

#### from CodonUsage to AaUsage
DataT = data.frame(t(data))
DataT$Codons = row.names(DataT)
DataT$Aa = DataT$Codons

TranslateMitCodonsIntoThreeLetterAa<-function(x)
{
  if (x %in% c('TTT','TTC')) {return ("Phe")}
  if (x %in% c('TTA','TTG')) {return ("LeuTT")}
  if (x %in% c('CTT','CTC','CTA','CTG')) {return ("LeuCT")}
  if (x %in% c('ATT','ATC')) {return ("Ile")}
  if (x %in% c('ATA','ATG')) {return ("Met")}
  if (x %in% c('GTC','GTA','GTG','GTT')) {return ("Val")}
  
  if (x %in% c('TCT','TCC','TCA','TCG')) {return ("SerTC")}
  if (x %in% c('CCT','CCC','CCA','CCG')) {return ("Pro")}
  if (x %in% c('ACT','ACC','ACA','ACG')) {return ("Thr")}
  if (x %in% c('GCT','GCC','GCA','GCG')) {return ("Ala")}
  
  if (x %in% c('TAT','TAC')) {return ("Tyr")}
  if (x %in% c('TAA','TAG')) {return ("Stop")}
  if (x %in% c('CAT','CAC')) {return ("His")}
  if (x %in% c('CAA','CAG')) {return ("Gln")}
  if (x %in% c('AAT','AAC')) {return ("Asn")}
  if (x %in% c('AAA','AAG')) {return ("Lys")}
  if (x %in% c('GAT','GAC')) {return ("Asp")}
  if (x %in% c('GAA','GAG')) {return ("Glu")}
  
  if (x %in% c('TGT','TGC')) {return ("Cys")}
  if (x %in% c('TGA','TGG')) {return ("Trp")}
  if (x %in% c('CGT','CGC','CGA','CGG')) {return ("Arg")}
  if (x %in% c('AGT','AGC')) {return ("SerAG")}
  if (x %in% c('AGA','AGG')) {return ("Stop")}
  if (x %in% c('GGT','GGC','GGA','GGG')) {return ("Gly")}
}

DataT$Aa = apply(as.matrix(DataT$Codons),1,FUN = TranslateMitCodonsIntoThreeLetterAa)
names(DataT)[1]=c('CodonFreq')
agg = aggregate(DataT$CodonFreq, by = list(DataT$Aa), FUN = sum)
names(agg)=c('Aa','Freq')

VecOfGainers = c('Pro','Thr','His','Gln','Asn','Lys')
VecOfLosers  = c('Phe','Val','Gly','Cys','Trp')

Gainers = sum(agg[agg$Aa %in% VecOfGainers,]$Freq)
Losers = sum(agg[agg$Aa %in% VecOfLosers,]$Freq)
Intermediate = sum(agg$Freq) - Gainers - Losers
agg = agg[order(-agg$Freq),]
agg$AaType = 'InterMediate'
for (i in 1:nrow(agg)) 
{
  if (agg$Aa[i] %in% VecOfGainers) {agg$AaType[i] = 'gainer'}
  if (agg$Aa[i] %in% VecOfLosers) {agg$AaType[i] =  'loser'}
}

##### PREDICTIONS OF VALERIAN

# задача разбивается на 4 простых случая:
#1)  Если кодон состоит только из пиримидинов (типа TTT), то в финале в этом кодоне будет (в среднем) 0.999 T и 2.001 С.
#2)  Если кодон состоит только из пуринов (типа AAA), то в финале в этом кодоне будет (в среднем) 2.769 A и 0.231 G.
#3)  Более интересный случай – когда кодон изначально состоит из 2 пиримидинов и 1 пурина (скажем, TTA), мы после достаточно продолжительного периода свободных мутаций получим кодон, в среднем состоящий из: 0.666 T; 1.3334 C; 0.923 A; и 0.077 G.
#4)  Наконец, для кодонов, состоящих из 1 пиримидина и 2 пуринов (скажем, CGG), мы в конце должны получить кодон, в среднем состоящий из: 0.333 T; 0.667 C; 1.846 A; и 0.154 G.

## first class of codons - only Pyrimidines (C + T): 
Class = DataT[DataT$Codons %in% c('TTT','TTC','CTT','CTC','TCT','TCC','CCT','CCC'),]
for (i in 1:nrow(Class))
{ # i = 1
nucs = unlist(strsplit(Class$Codons[i],'')); Class$NumberOfT[i] = length(nucs[nucs == 'T'])
}
FreqOfT = sum(Class$CodonFreq*Class$NumberOfT)
FreqOfC = sum(Class$CodonFreq*(3-Class$NumberOfT))
FreqOfT # 0.1301563
FreqOfC # 0.2420312
FreqOfC/FreqOfT # 1.859544 ~ 2.001/0.999 = 2.003 

## second class of codons - only Pyrines (A + G): 
Class = DataT[DataT$Codons %in% c('AAA','AAG','GAA','GAG','AGA','AGG','GGA','GGG'),]
for (i in 1:nrow(Class))
{ # i = 1
  nucs = unlist(strsplit(Class$Codons[i],'')); Class$NumberOfA[i] = length(nucs[nucs == 'A'])
}
FreqOfA = sum(Class$CodonFreq*Class$NumberOfA)
FreqOfG = sum(Class$CodonFreq*(3-Class$NumberOfA))
FreqOfA # 0.3476562
FreqOfG # 0.0353125
FreqOfA/FreqOfG # 9.845133 ~ 2.769 / 0.231 = 11.987

##### read simulations of Valerian and merge them with DataT

Val = read.table("../../Body/1Raw/VictorSimulations/Valerian01.txt", header = FALSE, quote = '')
Val=Val[c(1,3)]
names(Val) = c('Codons','FreqVal')
Val$FreqVal = as.numeric(gsub(';','',Val$FreqVal))
Val$FreqVal = Val$FreqVal/sum(Val$FreqVal)

DataT = merge(DataT,Val, by = 'Codons')
sum(DataT$FreqVal); sum(DataT$CodonFreq)
cor.test(DataT$FreqVal,DataT$CodonFreq,method = 'spearman')
plot(DataT$FreqVal,DataT$CodonFreq, xlab = 'Valerian Codon Usage', ylab = 'Victor Codon Usage')
plot(log(DataT$FreqVal),log(DataT$CodonFreq), xlab = 'Log Valerian Codon Usage', ylab = 'Log Victor Codon Usage')

DataT


