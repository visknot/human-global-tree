### add to Victor's mtDNA annotation, all four possible transitions and corresponding AaSub as well as TimeBeingSingleStranded

rm(list=ls(all=TRUE))

# library('Biostrings')

#mt = read.table("../../Body/1Raw/mtNucAnnotation_MergeRazerV6TurboMaxPro_normMTcode_normND6_withGERP.csv", header = TRUE, sep = ';')
# mt = read.table("../../Body/1Raw/mtNucAnnotation_MergeRazerV6TurboMaxPro_normMTcode_normND6_withGERP_corrected280121.csv", header = TRUE, sep = ';')
mt = read.table("../../Body/1Raw/mtDNAannotationMerge_300121.KPmodif.csv", header = TRUE, sep = '\t')

table(mt$role)
ProtCodGenes = c("mRNA_ATP6","mRNA_ATP6&COX3","mRNA_ATP8","mRNA_ATP8&ATP6","mRNA_COX1","mRNA_COX2","mRNA_COX3","mRNA_CYTB","mRNA_ND1","mRNA_ND2","mRNA_ND3","mRNA_ND4","mRNA_ND4L","mRNA_ND4L&ND4","mRNA_ND5","mRNA_ND6")
ProtCodGenesExceptNd6 = c("mRNA_ATP6","mRNA_ATP6&COX3","mRNA_ATP8","mRNA_ATP8&ATP6","mRNA_COX1","mRNA_COX2","mRNA_COX3","mRNA_CYTB","mRNA_ND1","mRNA_ND2","mRNA_ND3","mRNA_ND4","mRNA_ND4L","mRNA_ND4L&ND4","mRNA_ND5")
length(ProtCodGenes) # 13 genes + 3 overlaps

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

### walk nucleotide by nucleotide within ProtCodGenes, make all possible transitions for each nucleotide: G(Ch) <=> A(Th) or T(Ah)<=>C(Gh) 

mt$MutatedCodon = ''
for (i in 1:nrow(mt))
{
  if (mt$role[i] %in% ProtCodGenesExceptNd6)
  {
  if (mt$nuc[i] == 'T' & mt$nnuc[i] == 1)  {mt$MutatedCodon[i] = paste('C',unlist(strsplit(mt$RnaCodon[i],''))[2],unlist(strsplit(mt$RnaCodon[i],''))[3], sep = '')}
  if (mt$nuc[i] == 'G' & mt$nnuc[i] == 1)  {mt$MutatedCodon[i] = paste('A',unlist(strsplit(mt$RnaCodon[i],''))[2],unlist(strsplit(mt$RnaCodon[i],''))[3], sep = '')}
  if (mt$nuc[i] == 'C' & mt$nnuc[i] == 1)  {mt$MutatedCodon[i] = paste('T',unlist(strsplit(mt$RnaCodon[i],''))[2],unlist(strsplit(mt$RnaCodon[i],''))[3], sep = '')}
  if (mt$nuc[i] == 'A' & mt$nnuc[i] == 1)  {mt$MutatedCodon[i] = paste('G',unlist(strsplit(mt$RnaCodon[i],''))[2],unlist(strsplit(mt$RnaCodon[i],''))[3], sep = '')}    

  if (mt$nuc[i] == 'T' & mt$nnuc[i] == 2)  {mt$MutatedCodon[i] = paste(unlist(strsplit(mt$RnaCodon[i],''))[1],'C',unlist(strsplit(mt$RnaCodon[i],''))[3], sep = '')}
  if (mt$nuc[i] == 'G' & mt$nnuc[i] == 2)  {mt$MutatedCodon[i] = paste(unlist(strsplit(mt$RnaCodon[i],''))[1],'A',unlist(strsplit(mt$RnaCodon[i],''))[3], sep = '')}
  if (mt$nuc[i] == 'C' & mt$nnuc[i] == 2)  {mt$MutatedCodon[i] = paste(unlist(strsplit(mt$RnaCodon[i],''))[1],'T',unlist(strsplit(mt$RnaCodon[i],''))[3], sep = '')}
  if (mt$nuc[i] == 'A' & mt$nnuc[i] == 2)  {mt$MutatedCodon[i] = paste(unlist(strsplit(mt$RnaCodon[i],''))[1],'G',unlist(strsplit(mt$RnaCodon[i],''))[3], sep = '')}

  if (mt$nuc[i] == 'T' & mt$nnuc[i] == 3)  {mt$MutatedCodon[i] = paste(unlist(strsplit(mt$RnaCodon[i],''))[1],unlist(strsplit(mt$RnaCodon[i],''))[2],'C', sep = '')}
  if (mt$nuc[i] == 'G' & mt$nnuc[i] == 3)  {mt$MutatedCodon[i] = paste(unlist(strsplit(mt$RnaCodon[i],''))[1],unlist(strsplit(mt$RnaCodon[i],''))[2],'A', sep = '')}
  if (mt$nuc[i] == 'C' & mt$nnuc[i] == 3)  {mt$MutatedCodon[i] = paste(unlist(strsplit(mt$RnaCodon[i],''))[1],unlist(strsplit(mt$RnaCodon[i],''))[2],'T', sep = '')}
  if (mt$nuc[i] == 'A' & mt$nnuc[i] == 3)  {mt$MutatedCodon[i] = paste(unlist(strsplit(mt$RnaCodon[i],''))[1],unlist(strsplit(mt$RnaCodon[i],''))[2],'G', sep = '')}
  }
  if (mt$role[i] == "mRNA_ND6") # the same as above but with complementarity rules
  {
  if (mt$nuc[i] == 'T' & mt$nnuc[i] == 1)  {mt$MutatedCodon[i] = paste('G',unlist(strsplit(mt$RnaCodon[i],''))[2],unlist(strsplit(mt$RnaCodon[i],''))[3], sep = '')}
  if (mt$nuc[i] == 'G' & mt$nnuc[i] == 1)  {mt$MutatedCodon[i] = paste('T',unlist(strsplit(mt$RnaCodon[i],''))[2],unlist(strsplit(mt$RnaCodon[i],''))[3], sep = '')}
  if (mt$nuc[i] == 'C' & mt$nnuc[i] == 1)  {mt$MutatedCodon[i] = paste('A',unlist(strsplit(mt$RnaCodon[i],''))[2],unlist(strsplit(mt$RnaCodon[i],''))[3], sep = '')}
  if (mt$nuc[i] == 'A' & mt$nnuc[i] == 1)  {mt$MutatedCodon[i] = paste('C',unlist(strsplit(mt$RnaCodon[i],''))[2],unlist(strsplit(mt$RnaCodon[i],''))[3], sep = '')}    
  
  if (mt$nuc[i] == 'T' & mt$nnuc[i] == 2)  {mt$MutatedCodon[i] = paste(unlist(strsplit(mt$RnaCodon[i],''))[1],'G',unlist(strsplit(mt$RnaCodon[i],''))[3], sep = '')}
  if (mt$nuc[i] == 'G' & mt$nnuc[i] == 2)  {mt$MutatedCodon[i] = paste(unlist(strsplit(mt$RnaCodon[i],''))[1],'T',unlist(strsplit(mt$RnaCodon[i],''))[3], sep = '')}
  if (mt$nuc[i] == 'C' & mt$nnuc[i] == 2)  {mt$MutatedCodon[i] = paste(unlist(strsplit(mt$RnaCodon[i],''))[1],'A',unlist(strsplit(mt$RnaCodon[i],''))[3], sep = '')}
  if (mt$nuc[i] == 'A' & mt$nnuc[i] == 2)  {mt$MutatedCodon[i] = paste(unlist(strsplit(mt$RnaCodon[i],''))[1],'C',unlist(strsplit(mt$RnaCodon[i],''))[3], sep = '')}
  
  if (mt$nuc[i] == 'T' & mt$nnuc[i] == 3)  {mt$MutatedCodon[i] = paste(unlist(strsplit(mt$RnaCodon[i],''))[1],unlist(strsplit(mt$RnaCodon[i],''))[2],'G', sep = '')}
  if (mt$nuc[i] == 'G' & mt$nnuc[i] == 3)  {mt$MutatedCodon[i] = paste(unlist(strsplit(mt$RnaCodon[i],''))[1],unlist(strsplit(mt$RnaCodon[i],''))[2],'T', sep = '')}
  if (mt$nuc[i] == 'C' & mt$nnuc[i] == 3)  {mt$MutatedCodon[i] = paste(unlist(strsplit(mt$RnaCodon[i],''))[1],unlist(strsplit(mt$RnaCodon[i],''))[2],'A', sep = '')}
  if (mt$nuc[i] == 'A' & mt$nnuc[i] == 3)  {mt$MutatedCodon[i] = paste(unlist(strsplit(mt$RnaCodon[i],''))[1],unlist(strsplit(mt$RnaCodon[i],''))[2],'C', sep = '')}
  }
}

#MitCode = getGeneticCode(id_or_name2="2", full.search=FALSE, as.data.frame=FALSE)

TranslateMitCodonsIntoThreeLetterAa<-function(x)
{
  if (x %in% c('TTT','TTC')) {return ("Phe")}
  if (x %in% c('TTA','TTG')) {return ("Leu")}
  if (x %in% c('CTT','CTC','CTA','CTG')) {return ("Leu")}
  if (x %in% c('ATT','ATC')) {return ("Ile")}
  if (x %in% c('ATA','ATG')) {return ("Met")}
  if (x %in% c('GTC','GTA','GTG','GTT')) {return ("Val")}
  
  if (x %in% c('TCT','TCC','TCA','TCG')) {return ("Ser")}
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
  if (x %in% c('AGT','AGC')) {return ("Ser")}
  if (x %in% c('AGA','AGG')) {return ("Stop")}
  if (x %in% c('GGT','GGC','GGA','GGG')) {return ("Gly")}
}
# TranslateMitCodonsIntoThreeLetterAa('TTC')
# MitCode = getGeneticCode("SGC1")
mt$MutatedAa = ''
mt$AcidTranslated = ''
for (i in 1:nrow(mt))
{ # i = 7599
  if (mt$MutatedCodon[i] != '') 
  {
  mt$MutatedAa[i] = TranslateMitCodonsIntoThreeLetterAa(mt$MutatedCodon[i])
  mt$AcidTranslated[i] = TranslateMitCodonsIntoThreeLetterAa(mt$RnaCodon[i])
  #codon = DNAString(x=mt$MutatedCodon[i], start=1, nchar=NA); 
  # AaOneLetter = AAString(translate(codon,genetic.code=MitCode)); 
  # mt$MutatedAa[i] = as.character(AMINO_ACID_CODE[strsplit(as.character(AaOneLetter), NULL)[[1]]]);
  }
}

# mt$acid = gsub("STOP",'Stop',mt$acid)
mt$AaSub = ''
for (i in 1:nrow(mt))
{ # i = 10500
  if (mt$MutatedAa[i] != '' & !is.na(mt$MutatedAa[i])) 
  {
    # AncestralAa = gsub("\\/(.*)",'',mt$acid[i])
    AncestralAa = mt$AcidTranslated[i]
    if (AncestralAa != mt$MutatedAa[i])
      {
      mt$AaSub[i] = paste(AncestralAa,as.character(mt$MutatedAa[i]),sep='>')
      }
  }
}
table(mt$AaSub)

mt$TimeBeingSingleStranded = 0
for (i in 1:nrow(mt))
{ # i = 10500
  if (mt$pos[i] >= 5798) {mt$TimeBeingSingleStranded[i] = (mt$pos[i]-5798)*2}
  if (mt$pos[i] <  5798) {mt$TimeBeingSingleStranded[i] = 16569 - 5798 - (5798 - mt$pos[i]) + mt$pos[i]} # path on major arc: 16569 - 5798 - (5798 - mt$pos[i]); path on minor arc = mt$pos[i]
}

write.table(mt,"../../Body/2Derived/mtNucAnnotation_MergeRazerV6TurboMaxPro_normMTcode_normND6_withGERP.AllTransitions.Tbss.txt")

