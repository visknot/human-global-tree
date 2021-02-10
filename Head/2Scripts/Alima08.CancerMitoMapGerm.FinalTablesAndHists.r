###################################
##### prepare final tables, stats and hists from all analyses
###################################

rm(list=ls(all=TRUE)) 
pdf("../../Body/4Figures/Alima08.CancerMitoMapGerm.FinalTablesAndHists.r.pdf")

####################
## CANCERS
####################  

###### A. READ EXTENSIVELY MUTATED (ALL 4 TRANSITIONS) HUMAN MITOCHONDRIAL GENOME, 
###### TAKE ONLY PROTEIN-CODING GENES EXCEPT ND6 AND TRANSLATE THEM TO 3LETTER CODE,
###### AND GET FINALLY FREQUENCY OF 23 AMINOACIDS (WITH TWO LEU AND SER) IN THESE 12 GENES: TwelveGenesAaFreq (WILL USE THEM FOR NORMALIZATION LATER)

MitAnno = read.table("../../Body/2Derived/mtNucAnnotation_MergeRazerV6TurboMaxPro_normMTcode_normND6_withGERP.AllTransitions.Tbss.txt")
MitAnno = MitAnno[!is.na(MitAnno$acid),] # 11341 POSITIONS ARE IN PROTEIN-CODING REGIONS
TwelveGenes = MitAnno[MitAnno$role != 'mRNA_ND6',]

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

TwelveGenes$AcidNew = apply(as.matrix(TwelveGenes$RnaCodon),1,FUN = TranslateMitCodonsIntoThreeLetterAa)
TwelveGenesAaFreq = data.frame(table(TwelveGenes$AcidNew))  # should not be simply Leu or Ser !!!!!!!!!!!!!!!!!!!!!!
names(TwelveGenesAaFreq)=c('Acid','Freq')
TwelveGenesAaFreq = TwelveGenesAaFreq[!TwelveGenesAaFreq$Acid %in% c('Ser','Leu'),]
TwelveGenesAaFreq$Freq = TwelveGenesAaFreq$Freq/3 # should be round!!!!!!!!!!!!
TwelveGenesAaFreq$Freq = round(TwelveGenesAaFreq$Freq,0)

###### B: READ CANCERS STATISTICS FOR 12 GENES, PREPARED BEFORE (NUMBER OF EXPECTED, NUMBER OF UNEXPECTED)
###### INTRODUCE ANNOTATION OF 2 LEU AND 2 SER
###### MERGE WITH FREQUENCIES OF ANCESTRAL AMINO-ACIDS FOR 12 GENES, CALCULATE RATES AND RATIO OF RATES:

AllExceptNd6 = read.table("../../Body/3Results/Alima04.AsymmetryInCancers.12genes.txt", header = TRUE)
names(AllExceptNd6)
AllExceptNd6 = AllExceptNd6[,c(2:8)]
AllExceptNd6$ExpectedAminoAcidSubstBias = gsub('Leu>Ser','LeuTT>SerTC',AllExceptNd6$ExpectedAminoAcidSubstBias)
AllExceptNd6$ExpectedAminoAcidSubstBias = gsub('Ser>Asn','SerAG>Asn',AllExceptNd6$ExpectedAminoAcidSubstBias)
AllExceptNd6$ExpectedAminoAcidSubstBias = gsub('Ser>Pro','SerTC>Pro',AllExceptNd6$ExpectedAminoAcidSubstBias)
AllExceptNd6$ExpectedAminoAcidSubstBias = gsub('Phe>Leu','Phe>LeuCT',AllExceptNd6$ExpectedAminoAcidSubstBias)
AllExceptNd6$ExpectedAminoAcidSubstBias = gsub('Leu>Pro','LeuCT>Pro',AllExceptNd6$ExpectedAminoAcidSubstBias)
AllExceptNd6$ExpectedAminoAcidSubstBias = gsub('Gly>Ser','Gly>SerAG',AllExceptNd6$ExpectedAminoAcidSubstBias)
AllExceptNd6$ExpectedAminoAcidSubstBias = gsub('Phe>Ser','Phe>SerTC',AllExceptNd6$ExpectedAminoAcidSubstBias)
AllExceptNd6$FROM = gsub(">(.*)",'',AllExceptNd6$ExpectedAminoAcidSubstBias)
AllExceptNd6$TO = gsub("(.*)>",'',AllExceptNd6$ExpectedAminoAcidSubstBias)
AllExceptNd6 = merge(AllExceptNd6,TwelveGenesAaFreq, by.x = 'FROM', by.y = 'Acid') #  AllExceptNd6$FROM, by.y = TwelveGenesAaFreq$Acid)
names(AllExceptNd6)[names(AllExceptNd6) == 'Freq'] <- 'FreqFrom'

AllExceptNd6 = merge(AllExceptNd6,TwelveGenesAaFreq, by.x = 'TO', by.y = 'Acid') #  AllExceptNd6$FROM, by.y = TwelveGenesAaFreq$Acid)
names(AllExceptNd6)[names(AllExceptNd6) == 'Freq'] <- 'FreqTo'

AllExceptNd6$ExpectedRate = AllExceptNd6$NumberOfExpectedAaSubst/AllExceptNd6$FreqFrom
AllExceptNd6$UnexpectedRate = AllExceptNd6$NumberOfUnexpectedAaSubst/AllExceptNd6$FreqTo
AllExceptNd6$LogRatioOfExpToUnexpRates = log2(AllExceptNd6$ExpectedRate / AllExceptNd6$UnexpectedRate)
AllExceptNd6 = AllExceptNd6[(c('ExpectedAminoAcidSubstBias','NumberOfExpectedAaSubst','FreqFrom','ExpectedRate','NumberOfUnexpectedAaSubst','FreqTo','UnexpectedRate','LogRatioOfExpToUnexpRates'))] # data <- data[c("A", "B", "C")]
AllExceptNd6 = AllExceptNd6[order(AllExceptNd6$LogRatioOfExpToUnexpRates),]
AllExceptNd6$gene = 'TwelveGenes'

###### C. READ EXTENSIVELY MUTATED (ALL 4 TRANSITIONS) HUMAN MITOCHONDRIAL GENOME, 
###### TAKE ONLY PROTEIN-CODING GENE ND6 AND TRANSLATE THEM TO 3LETTER CODE,
###### AND GET FINALLY FREQUENCY OF 23 AMINOACIDS (WITH TWO LEU AND SER) IN THESE 12 GENES: TwelveGenesAaFreq (WILL USE THEM FOR NORMALIZATION LATER)

MitAnno = read.table("../../Body/2Derived/mtNucAnnotation_MergeRazerV6TurboMaxPro_normMTcode_normND6_withGERP.AllTransitions.Tbss.txt")
MitAnno = MitAnno[!is.na(MitAnno$acid),] # 11341 POSITIONS ARE IN PROTEIN-CODING REGIONS
TwelveGenes = MitAnno[MitAnno$role == 'mRNA_ND6',]

TwelveGenes$AcidNew = apply(as.matrix(TwelveGenes$RnaCodon),1,FUN = TranslateMitCodonsIntoThreeLetterAa)
Nd6AaFreq = data.frame(table(TwelveGenes$AcidNew))  # should not be simply Leu or Ser !!!!!!!!!!!!!!!!!!!!!!
names(Nd6AaFreq)=c('Acid','Freq')
Nd6AaFreq = Nd6AaFreq[!Nd6AaFreq$Acid %in% c('Ser','Leu'),]
Nd6AaFreq$Freq = Nd6AaFreq$Freq/3 # should be round!!!!!!!!!!!!
Nd6AaFreq$Freq = round(Nd6AaFreq$Freq,0)  #  21 ROWS NOT 23!!! TWO AMINO ACIDS ARE ABSENT! WHICH? 

###### D: READ CANCERS STATISTICS FOR ND6 GENE, PREPARED BEFORE (NUMBER OF EXPECTED, NUMBER OF UNEXPECTED)
###### INTRODUCE ANNOTATION OF 2 LEU AND 2 SER
###### MERGE WITH FREQUENCIES OF ANCESTRAL AMINO-ACIDS FOR ND6, CALCULATE RATES AND RATIO OF RATES:

Nd6 = read.table("../../Body/3Results/Alima04.AsymmetryInCancers.ND6.txt", header = TRUE)
names(Nd6)
Nd6 = Nd6[,c(2:8)]
Nd6$ExpectedAminoAcidSubstBias = gsub('Leu>Ser','LeuTT>SerTC',Nd6$ExpectedAminoAcidSubstBias)
Nd6$ExpectedAminoAcidSubstBias = gsub('Ser>Asn','SerAG>Asn',Nd6$ExpectedAminoAcidSubstBias)
Nd6$ExpectedAminoAcidSubstBias = gsub('Ser>Pro','SerTC>Pro',Nd6$ExpectedAminoAcidSubstBias)
Nd6$ExpectedAminoAcidSubstBias = gsub('Phe>Leu','Phe>LeuCT',Nd6$ExpectedAminoAcidSubstBias)
Nd6$ExpectedAminoAcidSubstBias = gsub('Leu>Pro','LeuCT>Pro',Nd6$ExpectedAminoAcidSubstBias)
Nd6$ExpectedAminoAcidSubstBias = gsub('Gly>Ser','Gly>SerAG',Nd6$ExpectedAminoAcidSubstBias)
Nd6$ExpectedAminoAcidSubstBias = gsub('Phe>Ser','Phe>SerTC',Nd6$ExpectedAminoAcidSubstBias)
Nd6$FROM = gsub(">(.*)",'',Nd6$ExpectedAminoAcidSubstBias)
Nd6$TO = gsub("(.*)>",'',Nd6$ExpectedAminoAcidSubstBias)
Nd6 = merge(Nd6,Nd6AaFreq, by.x = 'FROM', by.y = 'Acid') #  Nd6$FROM, by.y = TwelveGenesAaFreq$Acid)
names(Nd6)[names(Nd6) == 'Freq'] <- 'FreqFrom'

Nd6 = merge(Nd6,Nd6AaFreq, by.x = 'TO', by.y = 'Acid') #  Nd6$FROM, by.y = TwelveGenesAaFreq$Acid)
names(Nd6)[names(Nd6) == 'Freq'] <- 'FreqTo'

Nd6$ExpectedRate = Nd6$NumberOfExpectedAaSubst/Nd6$FreqFrom
Nd6$UnexpectedRate = Nd6$NumberOfUnexpectedAaSubst/Nd6$FreqTo
Nd6$LogRatioOfExpToUnexpRates = log2(Nd6$ExpectedRate / Nd6$UnexpectedRate)
Nd6 = Nd6[(c('ExpectedAminoAcidSubstBias','NumberOfExpectedAaSubst','FreqFrom','ExpectedRate','NumberOfUnexpectedAaSubst','FreqTo','UnexpectedRate','LogRatioOfExpToUnexpRates'))] # data <- data[c("A", "B", "C")]
Nd6 = Nd6[order(-Nd6$LogRatioOfExpToUnexpRates),]
Nd6$gene = 'Nd6'

###### 
###### E: RBIND AllExceptNd6 and Nd6 and write final cancer file 'All'
###### 

Cancer = rbind(AllExceptNd6,Nd6)  # print it and make a table!!!
wilcox.test(Cancer[Cancer$gene == 'TwelveGenes',]$LogRatioOfExpToUnexpRates, mu = 0)     # 4.768e-07
wilcox.test(Cancer[Cancer$gene == 'Nd6',]$LogRatioOfExpToUnexpRates, mu = 0)             # 0.007813

MIN = round((min(Cancer$LogRatioOfExpToUnexpRates)-0.5),0)
MAX = round((max(Cancer$LogRatioOfExpToUnexpRates)+0.5),0)

breaks=seq(-10,10,0.5)
hist(Cancer[Cancer$gene == 'TwelveGenes',]$LogRatioOfExpToUnexpRates, xlim = c(MIN,MAX), ylim = c(0,6), col = rgb(1,0.1,0.1,0.3), xlab = 'log2(forward rate / backward rate)', main = 'somatic mutations from cancers', breaks = breaks); 
#hist(Cancer[Cancer$gene == 'Nd6',]$LogRatioOfExpToUnexpRates, xlim = c(MIN,MAX),ylim = c(0,6), col = rgb(0.1,0.1,0.1,0.3), xlab = 'log2(forward rate / backward rate)', main = 'somatic mutations from cancers', breaks = breaks);
abline(v=0, lt = 2, lwd = 2, col = 'black') # dev.off() print it as a main figure
#text(-4,6,'ND6', col = 'grey', cex = 2)
#text(6,6,'12 genes', col = 'pink', cex = 2)


breaks=seq(-10,10,0.5)
hist(Cancer[Cancer$gene == 'TwelveGenes',]$LogRatioOfExpToUnexpRates, xlim = c(MIN,MAX), ylim = c(0,6), col = rgb(1,0.1,0.1,0.3), xlab = '', main = '', breaks = breaks); par(new = TRUE)
hist(Cancer[Cancer$gene == 'Nd6',]$LogRatioOfExpToUnexpRates, xlim = c(MIN,MAX),ylim = c(0,6), col = rgb(0.1,0.1,0.1,0.3), xlab = 'log2(forward rate / backward rate)', main = 'somatic mutations from cancers', breaks = breaks);
abline(v=0, lt = 2, lwd = 2, col = 'black') # dev.off() print it as a main figure
text(-4,6,'ND6', col = 'grey', cex = 2)
text(6,6,'12 genes', col = 'pink', cex = 2)


Cancer$DataSet = 'Cancer'

#plot(Cancer[Cancer$gene == 'TwelveGenes',]$ExpectedRate,Cancer[Cancer$gene == 'TwelveGenes',]$UnexpectedRate)
#abline(a=0,b=1)

#########
## PATHOGENIC HUMAN MUTATIONS - 12 genes and Nd6 - all positions
########

###### A - will use numbers of amino acids from 12 genes, derived above (TwelveGenesAaFreq) 

###### B: READ MitoMap STATISTICS FOR 12 GENES, PREPARED BEFORE (NUMBER OF EXPECTED, NUMBER OF UNEXPECTED)
###### INTRODUCE ANNOTATION OF 2 LEU AND 2 SER
###### MERGE WITH FREQUENCIES OF ANCESTRAL AMINO-ACIDS FOR 12 GENES, CALCULATE RATES AND RATIO OF RATES:

AllExceptNd6 = read.table("../../Body/3Results/Alima02B.AaAsymmetryMitoMap.TwelveGenes.txt", header = TRUE)
AllExceptNd6$ExpectedAminoAcidSubstBias = AllExceptNd6$AaSub
AllExceptNd6$ExpectedAminoAcidSubstBias = gsub('Leu>Ser','LeuTT>SerTC',AllExceptNd6$ExpectedAminoAcidSubstBias)
AllExceptNd6$ExpectedAminoAcidSubstBias = gsub('Ser>Asn','SerAG>Asn',AllExceptNd6$ExpectedAminoAcidSubstBias)
AllExceptNd6$ExpectedAminoAcidSubstBias = gsub('Ser>Pro','SerTC>Pro',AllExceptNd6$ExpectedAminoAcidSubstBias)
AllExceptNd6$ExpectedAminoAcidSubstBias = gsub('Phe>Leu','Phe>LeuCT',AllExceptNd6$ExpectedAminoAcidSubstBias)
AllExceptNd6$ExpectedAminoAcidSubstBias = gsub('Leu>Pro','LeuCT>Pro',AllExceptNd6$ExpectedAminoAcidSubstBias)
AllExceptNd6$ExpectedAminoAcidSubstBias = gsub('Gly>Ser','Gly>SerAG',AllExceptNd6$ExpectedAminoAcidSubstBias)
AllExceptNd6$ExpectedAminoAcidSubstBias = gsub('Phe>Ser','Phe>SerTC',AllExceptNd6$ExpectedAminoAcidSubstBias)
AllExceptNd6$FROM = gsub(">(.*)",'',AllExceptNd6$ExpectedAminoAcidSubstBias)
AllExceptNd6$TO = gsub("(.*)>",'',AllExceptNd6$ExpectedAminoAcidSubstBias)
AllExceptNd6 = merge(AllExceptNd6,TwelveGenesAaFreq, by.x = 'FROM', by.y = 'Acid') #  AllExceptNd6$FROM, by.y = TwelveGenesAaFreq$Acid)
names(AllExceptNd6)[names(AllExceptNd6) == 'Freq'] <- 'FreqFrom'

AllExceptNd6 = merge(AllExceptNd6,TwelveGenesAaFreq, by.x = 'TO', by.y = 'Acid') #  AllExceptNd6$FROM, by.y = TwelveGenesAaFreq$Acid)
names(AllExceptNd6)[names(AllExceptNd6) == 'Freq'] <- 'FreqTo'

AllExceptNd6$ExpectedRate = AllExceptNd6$NumberOfExpectedAaSubst/AllExceptNd6$FreqFrom
AllExceptNd6$UnexpectedRate = AllExceptNd6$NumberOfUnexpectedAaSubst/AllExceptNd6$FreqTo
AllExceptNd6$LogRatioOfExpToUnexpRates = log2(AllExceptNd6$ExpectedRate / AllExceptNd6$UnexpectedRate)
AllExceptNd6 = AllExceptNd6[(c('ExpectedAminoAcidSubstBias','NumberOfExpectedAaSubst','FreqFrom','ExpectedRate','NumberOfUnexpectedAaSubst','FreqTo','UnexpectedRate','LogRatioOfExpToUnexpRates'))] # data <- data[c("A", "B", "C")]
AllExceptNd6 = AllExceptNd6[order(AllExceptNd6$LogRatioOfExpToUnexpRates),]
AllExceptNd6$gene = 'TwelveGenes'

###### C - will use numbers of amino acids from ND6, derived above (Nd6AaFreq) 

###### D read statistics about ND6 and merge with amino acid frequencies
# Nd6 = read.table("../../Body/3Results/Alima02.AaAsymmetry.HumanTree.ND6.ExpectedVsObservedAaChanges.txt", header = TRUE)
Nd6 = read.table("../../Body/3Results/Alima02B.AaAsymmetryMitoMap.Nd6.txt", header = TRUE)
Nd6$ExpectedAminoAcidSubstBias = Nd6$AaSub
Nd6$ExpectedAminoAcidSubstBias = gsub('Leu>Ser','LeuTT>SerTC',Nd6$ExpectedAminoAcidSubstBias)
Nd6$ExpectedAminoAcidSubstBias = gsub('Ser>Asn','SerAG>Asn',Nd6$ExpectedAminoAcidSubstBias)
Nd6$ExpectedAminoAcidSubstBias = gsub('Ser>Pro','SerTC>Pro',Nd6$ExpectedAminoAcidSubstBias)
Nd6$ExpectedAminoAcidSubstBias = gsub('Phe>Leu','Phe>LeuCT',Nd6$ExpectedAminoAcidSubstBias)
Nd6$ExpectedAminoAcidSubstBias = gsub('Leu>Pro','LeuCT>Pro',Nd6$ExpectedAminoAcidSubstBias)
Nd6$ExpectedAminoAcidSubstBias = gsub('Gly>Ser','Gly>SerAG',Nd6$ExpectedAminoAcidSubstBias)
Nd6$ExpectedAminoAcidSubstBias = gsub('Phe>Ser','Phe>SerTC',Nd6$ExpectedAminoAcidSubstBias)
Nd6$FROM = gsub(">(.*)",'',Nd6$ExpectedAminoAcidSubstBias)
Nd6$TO = gsub("(.*)>",'',Nd6$ExpectedAminoAcidSubstBias)
Nd6 = merge(Nd6,Nd6AaFreq, by.x = 'FROM', by.y = 'Acid') #  Nd6$FROM, by.y = TwelveGenesAaFreq$Acid)
names(Nd6)[names(Nd6) == 'Freq'] <- 'FreqFrom'

Nd6 = merge(Nd6,Nd6AaFreq, by.x = 'TO', by.y = 'Acid') #  Nd6$FROM, by.y = TwelveGenesAaFreq$Acid)
names(Nd6)[names(Nd6) == 'Freq'] <- 'FreqTo'

Nd6$ExpectedRate = Nd6$NumberOfExpectedAaSubst/Nd6$FreqFrom
Nd6$UnexpectedRate = Nd6$NumberOfUnexpectedAaSubst/Nd6$FreqTo
Nd6$LogRatioOfExpToUnexpRates = log2(Nd6$ExpectedRate / Nd6$UnexpectedRate)
Nd6 = Nd6[(c('ExpectedAminoAcidSubstBias','NumberOfExpectedAaSubst','FreqFrom','ExpectedRate','NumberOfUnexpectedAaSubst','FreqTo','UnexpectedRate','LogRatioOfExpToUnexpRates'))] # data <- data[c("A", "B", "C")]
Nd6 = Nd6[order(-Nd6$LogRatioOfExpToUnexpRates),]
Nd6$gene = 'Nd6'

###### 
###### E: RBIND AllExceptNd6 and Nd6 and write final MitoMap file 'MitoMap'
###### 

MitoMap = rbind(AllExceptNd6,Nd6)  # print it and make a table!!!
wilcox.test(MitoMap[MitoMap$gene == 'TwelveGenes',]$LogRatioOfExpToUnexpRates, mu = 0)     # 0.0001327
wilcox.test(MitoMap[MitoMap$gene == 'Nd6',]$LogRatioOfExpToUnexpRates, mu = 0)             # 0.1967

#MIN = round((min(MitoMap$LogRatioOfExpToUnexpRates)-0.5),0)
#MAX = round((max(MitoMap[MitoMap$LogRatioOfExpToUnexpRates < Inf,]$LogRatioOfExpToUnexpRates)+0.5),0)

hist(MitoMap[MitoMap$gene == 'TwelveGenes',]$LogRatioOfExpToUnexpRates, xlim = c(MIN,MAX), ylim = c(0,6), col = rgb(1,0.1,0.1,0.3), xlab = '', main = '', breaks = breaks); par(new = TRUE)
hist(MitoMap[MitoMap$gene == 'Nd6',]$LogRatioOfExpToUnexpRates, xlim = c(MIN,MAX),ylim = c(0,6), col = rgb(0.1,0.1,0.1,0.3), xlab = 'log2(forward rate / backward rate)', main = 'pathogenic mutations', breaks = breaks);
abline(v=0, lt = 2, lwd = 2, col = 'black') # dev.off() print it as a main figure
text(-4,6,'ND6', col = 'grey', cex = 2)
text(6,6,'12 genes', col = 'pink', cex = 2)


MitoMap$DataSet = 'MitoMap'

#########
## GERM-LINE SUBSTITUTIONS FROM THE HUMAN TREE 
#########

##### A - will use numbers of amino acids from 12 genes, derived above (TwelveGenesAaFreq) 

##### B READ GERM STATISTICS FOR 12 GENES, PREPARED BEFORE (NUMBER OF EXPECTED, NUMBER OF UNEXPECTED)
##### DON'T NEED TO INTRODUCE ANNOTATION OF 2 LEU AND 2 SER (IT IS ALREADY)
##### MERGE WITH FREQUENCIES OF ANCESTRAL AMINO-ACIDS FOR 12 GENES, CALCULATE RATES AND RATIO OF RATES

AllExceptNd6 = read.table("../../Body/3Results/Alima02.AaAsymmetry.HumanTree.12Genes.ExpectedVsObservedAaChanges.txt", header = TRUE)
names(AllExceptNd6)
AllExceptNd6 = AllExceptNd6[,c(2:8)]

AllExceptNd6$FROM = gsub(">(.*)",'',AllExceptNd6$ExpectedAminoAcidSubstBias)
AllExceptNd6$TO = gsub("(.*)>",'',AllExceptNd6$ExpectedAminoAcidSubstBias)
AllExceptNd6 = merge(AllExceptNd6,TwelveGenesAaFreq, by.x = 'FROM', by.y = 'Acid') #  AllExceptNd6$FROM, by.y = TwelveGenesAaFreq$Acid)
names(AllExceptNd6)[names(AllExceptNd6) == 'Freq'] <- 'FreqFrom'

AllExceptNd6 = merge(AllExceptNd6,TwelveGenesAaFreq, by.x = 'TO', by.y = 'Acid') #  AllExceptNd6$FROM, by.y = TwelveGenesAaFreq$Acid)
names(AllExceptNd6)[names(AllExceptNd6) == 'Freq'] <- 'FreqTo'

AllExceptNd6$ExpectedRate = AllExceptNd6$NumberOfExpectedAaSubst/AllExceptNd6$FreqFrom
AllExceptNd6$UnexpectedRate = AllExceptNd6$NumberOfUnexpectedAaSubst/AllExceptNd6$FreqTo
AllExceptNd6$LogRatioOfExpToUnexpRates = log2(AllExceptNd6$ExpectedRate / AllExceptNd6$UnexpectedRate)
AllExceptNd6 = AllExceptNd6[(c('ExpectedAminoAcidSubstBias','NumberOfExpectedAaSubst','FreqFrom','ExpectedRate','NumberOfUnexpectedAaSubst','FreqTo','UnexpectedRate','LogRatioOfExpToUnexpRates'))] # data <- data[c("A", "B", "C")]
AllExceptNd6 = AllExceptNd6[order(AllExceptNd6$LogRatioOfExpToUnexpRates),]
AllExceptNd6$gene = 'TwelveGenes'

###### C - will use numbers of amino acids from ND6, derived above (Nd6AaFreq) 

###### D read statistics about ND6 and merge with amino acid frequencies
Nd6 = read.table("../../Body/3Results/Alima02.AaAsymmetry.HumanTree.ND6.ExpectedVsObservedAaChanges.txt", header = TRUE)
names(Nd6)
Nd6 = Nd6[,c(2:8)]

Nd6$ExpectedAminoAcidSubstBias = gsub('Leu>Ser','LeuTT>SerTC',Nd6$ExpectedAminoAcidSubstBias)
Nd6$ExpectedAminoAcidSubstBias = gsub('Ser>Asn','SerAG>Asn',Nd6$ExpectedAminoAcidSubstBias)
Nd6$ExpectedAminoAcidSubstBias = gsub('Ser>Pro','SerTC>Pro',Nd6$ExpectedAminoAcidSubstBias)
Nd6$ExpectedAminoAcidSubstBias = gsub('Phe>Leu','Phe>LeuCT',Nd6$ExpectedAminoAcidSubstBias)
Nd6$ExpectedAminoAcidSubstBias = gsub('Leu>Pro','LeuCT>Pro',Nd6$ExpectedAminoAcidSubstBias)
Nd6$ExpectedAminoAcidSubstBias = gsub('Gly>Ser','Gly>SerAG',Nd6$ExpectedAminoAcidSubstBias)
Nd6$ExpectedAminoAcidSubstBias = gsub('Phe>Ser','Phe>SerTC',Nd6$ExpectedAminoAcidSubstBias)
Nd6$FROM = gsub(">(.*)",'',Nd6$ExpectedAminoAcidSubstBias)
Nd6$TO = gsub("(.*)>",'',Nd6$ExpectedAminoAcidSubstBias)
Nd6 = merge(Nd6,Nd6AaFreq, by.x = 'FROM', by.y = 'Acid') #  Nd6$FROM, by.y = TwelveGenesAaFreq$Acid)
names(Nd6)[names(Nd6) == 'Freq'] <- 'FreqFrom'

Nd6 = merge(Nd6,Nd6AaFreq, by.x = 'TO', by.y = 'Acid') #  Nd6$FROM, by.y = TwelveGenesAaFreq$Acid)
names(Nd6)[names(Nd6) == 'Freq'] <- 'FreqTo'

Nd6$ExpectedRate = Nd6$NumberOfExpectedAaSubst/Nd6$FreqFrom
Nd6$UnexpectedRate = Nd6$NumberOfUnexpectedAaSubst/Nd6$FreqTo
Nd6$LogRatioOfExpToUnexpRates = log2(Nd6$ExpectedRate / Nd6$UnexpectedRate)
Nd6 = Nd6[(c('ExpectedAminoAcidSubstBias','NumberOfExpectedAaSubst','FreqFrom','ExpectedRate','NumberOfUnexpectedAaSubst','FreqTo','UnexpectedRate','LogRatioOfExpToUnexpRates'))] # data <- data[c("A", "B", "C")]
Nd6 = Nd6[order(-Nd6$LogRatioOfExpToUnexpRates),]
Nd6$gene = 'Nd6'

###### 
###### E: RBIND AllExceptNd6 and Nd6 and write final file Germ
###### 

Germ = rbind(AllExceptNd6,Nd6)  # print it and make a table!!!
summary(Germ[Germ$gene == 'TwelveGenes',]$LogRatioOfExpToUnexpRates) # median and mean is a bit positive
wilcox.test(Germ[Germ$gene == 'TwelveGenes',]$LogRatioOfExpToUnexpRates, mu = 0)                              # 0.05625
wilcox.test(Germ[Germ$gene == 'TwelveGenes',]$LogRatioOfExpToUnexpRates, mu = 0, alternative = 'greater')     # 0.02813

summary(Germ[Germ$gene == 'Nd6',]$LogRatioOfExpToUnexpRates) # median and mean are negative
wilcox.test(Germ[Germ$gene == 'Nd6',]$LogRatioOfExpToUnexpRates, mu = 0)                                   # 0.09229
wilcox.test(Germ[Germ$gene == 'Nd6',]$LogRatioOfExpToUnexpRates, mu = 0, alternative = 'less')             # 0.04614

#MIN = round((min(Germ$LogRatioOfExpToUnexpRates)-0.5),0)
#MAX = round((max(Germ$LogRatioOfExpToUnexpRates)+0.5),0)

hist(Germ[Germ$gene == 'TwelveGenes',]$LogRatioOfExpToUnexpRates, xlim = c(MIN,MAX), ylim = c(0,6), col = rgb(1,0.1,0.1,0.3), xlab = '', main = '', breaks = breaks); par(new = TRUE)
hist(Germ[Germ$gene == 'Nd6',]$LogRatioOfExpToUnexpRates, xlim = c(MIN,MAX),ylim = c(0,6), col = rgb(0.1,0.1,0.1,0.3), xlab = 'log2(forward rate / backward rate)', main = 'germ-line substitutions in the human population', breaks = breaks);
abline(v=0, lt = 2, lwd = 2, col = 'black') # dev.off() print it as a main figure
text(-4,6,'ND6', col = 'grey', cex = 2)
text(6,6,'12 genes', col = 'pink', cex = 2)

Germ$DataSet = 'Germ'

# I can run fisher test for each substitution
#fisher.test(rbind(c(1077,183),c(981,48))) # first is expected data, second - unexpected
#fisher.test(rbind(c(168,183),c(160,63))) # first is expected data, second - unexpected
#fisher.test(rbind(c(502,22),c(502,124))) # first is expected data, second - unexpected
# think about permutation?
# plot(fisher.test(rbind(c(502,22),c(502,124))))

Final=rbind(Cancer,MitoMap,Germ)
Final$ExpectedRate = round(Final$ExpectedRate,3)
Final$UnexpectedRate = round(Final$UnexpectedRate,3)
Final$LogRatioOfExpToUnexpRates = round(Final$LogRatioOfExpToUnexpRates,3)
write.table(Final,'../../Body/3Results/Alima08.CancerMitoMapGerm.FinalTablesAndHists.r.txt') 
# copy and paste into LibreCalc => copy and special paste (rft)

dev.off()

############################################################################
######### OLD CODE FOR CONSTRAINTS, VERTEBRATE POLYMORPHISMS AND REFSEQS:
############################################################################

# #########
# ## GermConst - 12 genes: different constraints:
# ########
# 
# LowConst = read.table("../../Body/3Results/Alima02.AaAsymmetry.HumanTree.12GenesLowConst.ExpectedVsObservedAaChanges.txt", header = TRUE)
# names(LowConst)
# LowConst = LowConst[,c(2:8)]
# 
# HighConst = read.table("../../Body/3Results/Alima02.AaAsymmetry.HumanTree.12GenesHighConst.ExpectedVsObservedAaChanges.txt", header = TRUE)
# names(HighConst)
# HighConst = HighConst[,c(4,6,7,8)]
# names(HighConst)[2:4]=c('ExpectedMoreThanOne.HighConst','NumberOfExpectedAaSubst.HighConst','NumberOfUnexpectedAaSubst.HighConst')
# All = merge(LowConst,HighConst, all = TRUE)
# #All[is.na(All$NumberOfExpectedAaSubst.HighPh),]$NumberOfExpectedAaSubst.HighPh<-0
# #All[is.na(All$NumberOfUnexpectedAaSubst.HighPh),]$NumberOfUnexpectedAaSubst.Nd6<-0
# # All = All[order(-All$L),] - к черту
# #names(All)
# 
# All$Ratio = All$ExpectedMoreThanOne / All$ExpectedMoreThanOne.HighConst
# summary(All$Ratio) # a bit more than one (because LocConst to HighConst)
# wilcox.test(All$Ratio, mu = 1)
# wilcox.test(All$ExpectedMoreThanOne,All$ExpectedMoreThanOne.HighConst, paired = TRUE)
# # cor.test(All$Ratio,All$GranthamDistance, method = 'spearman')
# 
# #Final = All[,c(1,6,7,5,9,10,8)]
# #Final$ExpectedMoreThanOne = round(log2(Final$ExpectedMoreThanOne),2); summary(Final$ExpectedMoreThanOne)
# #Final$ExpectedMoreThanOne.HighPh = round(log2(Final$ExpectedMoreThanOne.HighPh),2); summary(Final$ExpectedMoreThanOne.HighPh)
# 
# #Final = Final[order(Final$ExpectedMoreThanOne),] 
# # dev.off()
# #breaks = seq(-1,1,0.05)
# #hist(Final$ExpectedMoreThanOne, xlim = c(-0.3,0.6), ylim = c(0,6), col = rgb(1,0.1,0.1,0.3), xlab = '', main = '', breaks = breaks); par(new = TRUE)
# #hist(Final$ExpectedMoreThanOne.HighPh, xlim = c(-0.3,0.6),ylim = c(0,6), col = rgb(0.1,0.1,0.1,0.3), xlab = 'log2(expected/unexpected)', main = '', breaks = breaks);
# #abline(v=0, lt = 2, lwd = 2, col = 'black')
# 
# #wilcox.test(Final$ExpectedMoreThanOne, mu = 0)     # 9.199e-05
# #wilcox.test(Final$ExpectedMoreThanOne.HighPh, mu = 0) # 0.0007255
# 
# GermConstr = Final
# 
# ########
# ## VertebratePolym (one type of files - all genes without ND6)
# ########
# 
# for (set in 1:6)
# { # set = 1
#   if (set == 1) {input = '../../Body/3Results/Alima05.AaAsymmetryVertebratePolymorphisms.r.AllClasses.txt'; class = 'AllClasses'}
#   if (set == 2) {input = '../../Body/3Results/Alima05.AaAsymmetryVertebratePolymorphisms.r.Actinopterygii.txt'; class = 'Actinopterygii'}
#   if (set == 3) {input = '../../Body/3Results/Alima05.AaAsymmetryVertebratePolymorphisms.r.Amphibia.txt'; class = 'Amphibia'}
#   if (set == 4) {input = '../../Body/3Results/Alima05.AaAsymmetryVertebratePolymorphisms.r.Reptilia.txt'; class = 'Reptilia'}
#   if (set == 5) {input = '../../Body/3Results/Alima05.AaAsymmetryVertebratePolymorphisms.r.Mammalia.txt'; class = 'Mammalia'}
#   if (set == 6) {input = '../../Body/3Results/Alima05.AaAsymmetryVertebratePolymorphisms.r.Aves.txt'; class = 'Aves'}    
#   
#   AllExceptNd6 = read.table(input, header = TRUE)
#   names(AllExceptNd6)
#   AllExceptNd6 = AllExceptNd6[,c(2:8)]
#   AllExceptNd6 = AllExceptNd6[]
#   
#   names(AllExceptNd6)
#   Final = AllExceptNd6[,c(3,6,7,5)]
#   Final$ExpectedMoreThanOne = round(log2(Final$ExpectedMoreThanOne),2); summary(Final$ExpectedMoreThanOne)
#   Final = Final[order(Final$ExpectedMoreThanOne),]
#   
#   wilcox.test(Final$ExpectedMoreThanOne, mu = 0)     # 6.333e-05
#   breaks = seq(-1,3,0.25)
#   hist(Final$ExpectedMoreThanOne, xlim = c(-1,3), ylim = c(0,6), col = rgb(1,0.1,0.1,0.3), xlab = '', main = class, breaks = breaks); # par(new = TRUE)
#   abline(v=0, lt = 2, lwd = 2, col = 'black')
#   
#   names(SuperFinal)
#   Final$NumberOfExpectedAaSubst.Nd6 = '';
#   Final$NumberOfUnexpectedAaSubst.Nd6 = '';
#   Final$ExpectedMoreThanOne.Nd6 = '';
#   Final$DataSet = class  # rbind in the loop!!!
#   
#   SuperFinal=rbind(SuperFinal,Final)
# }
# 
# ###### PRINT OUT EVERYTHING
# 
# 
# 
