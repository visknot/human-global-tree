### add to Victor's mtDNA annotation, all four possible transitions and corresponding AaSub as well as TimeBeingSingleStranded

rm(list=ls(all=TRUE))

library('Biostrings')

mt = read.table("../../Body/1Raw/mtNucAnnotation_MergeRazerV6TurboMaxPro_normMTcode_normND6_withGERP.csv", header = TRUE, sep = ';')
table(mt$role)
ProtCodGenes = c("mRNA_ATP6","mRNA_ATP6&COX3","mRNA_ATP8","mRNA_ATP8&ATP6","mRNA_COX1","mRNA_COX2","mRNA_COX3","mRNA_CYTB","mRNA_ND1","mRNA_ND2","mRNA_ND3","mRNA_ND4","mRNA_ND4L","mRNA_ND4L&ND4","mRNA_ND5","mRNA_ND6")
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
  if (mt$role[i] %in% ProtCodGenes)
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
mt$MutatedAa[is.na(mt$MutatedAa)] <- 'STOP'

mt$AaSub = ''
for (i in 1:nrow(mt))
{ # i = 10500
  if (mt$MutatedAa[i] != '' & !is.na(mt$MutatedAa[i])) 
  {
    AncestralAa = gsub("\\/(.*)",'',mt$acid[i])
    if (AncestralAa != mt$MutatedAa[i])
      {
      mt$AaSub[i] = paste(AncestralAa,as.character(mt$MutatedAa[i]),sep='>')
      }
  }
}

mt$TimeBeingSingleStranded = 0
for (i in 1:nrow(mt))
{ # i = 10500
  if (mt$pos[i] >= 5798) {mt$TimeBeingSingleStranded[i] = (mt$pos[i]-5798)*2}
  if (mt$pos[i] <  5798) {mt$TimeBeingSingleStranded[i] = 16569 - 5798 - (5798 - mt$pos[i]) + mt$pos[i]} # path on major arc: 16569 - 5798 - (5798 - mt$pos[i]); path on minor arc = mt$pos[i]
}

write.table(mt,"../../Body/2Derived/'mtNucAnnotation_MergeRazerV6TurboMaxPro_normMTcode_normND6_withGERP.AllTransitions.Tbss.txt")

