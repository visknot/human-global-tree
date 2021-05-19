rm(list=ls(all=TRUE)) 

sm4b = read.table("../../Body/1Raw/DogCancer/Strakova2016SupplMat4b.txt", header = TRUE, quote = '')
table(sm4b$Substitution)
sm4b$SnpId = paste(sm4b$VariantPosition ,sm4b$Substitution,sep='.')
length(unique(sm4b$SnpId)) # 814

sm6a = read.table("../../Body/1Raw/DogCancer/Strakova2016SupplMat6a.txt", header = TRUE, quote = '', sep = '\t')
sm6a$SnpId = paste(sm6a$Position ,sm6a$BaseChangeRefToAlt,sep='.')
head(sm6a)
length(unique(sm6a$SnpId)) # 882

nrow(sm4b) # 928 > 814 ????? VAF !!!!!!!
mut = merge(sm4b,sm6a, by = 'SnpId')
nrow(mut) # 929 one more!!!!!

head(mut)
table(mut$Impact)

MutNons = mut[mut$Impact == 'Non-synonymous variant',]
table(MutNons$UniProtKBGeneName)
MutNons = MutNons[MutNons$UniProtKBGeneName != 'ND6',]

### observed:
table(MutNons$AminoAcidChange)
table(MutNons$AminoAcidChange)

### expected: Sasha Voronka, ideal table for RefSeq mtDNA dog!!!


