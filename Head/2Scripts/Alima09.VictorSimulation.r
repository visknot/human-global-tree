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
agg = aggregate(DataT$X1001, by = list(DataT$Aa), FUN = sum)
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




