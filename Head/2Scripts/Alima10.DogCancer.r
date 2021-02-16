rm(list=ls(all=TRUE)) 

sm4b = read.table("../../Body/1Raw/DogCancer/Strakova2016SupplMat4b.txt", header = TRUE, quote = '')
sm4b$SnpId = paste(sm4b$VariantPosition ,sm4b$Substitution,sep='.')
length(unique(sm4b$SnpId))
sm6a = read.table("../../Body/1Raw/DogCancer/Strakova2016SupplMat6a.txt", header = TRUE, quote = '', sep = '\t')
sm6a$SnpId = paste(sm6a$Position ,sm6a$BaseChangeRefToAlt,sep='.')
length(unique(sm6a$SnpId))

nrow(sm4b) # 928
mut = merge(sm4b,sm6a, by = 'SnpId')
nrow(mut) # 929 one more!!!!

table(mut$Impact)
mut = mut[mut$Impact == 'Non-synonymous variant',]
table(mut$UniProtKBGeneName)
mut = mut[mut$UniProtKBGeneName != 'ND6',]

table(mut$AminoAcidChange)

