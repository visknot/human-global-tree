### PRODUCE MutSpec12 and MutSpec192 for different subsets of branches, genes and substitutions
rm(list=ls(all=TRUE))
library('dplyr')
library('gtools')

Mut = read.table('../../Body/2Derived/fulltreeCodons.csv', header = TRUE, sep = ';'); # "../../Body/3Results/fulltreeCodons.csv", sep=";"

##### DEFINE MAIN SUBSETS OF MUTATIONS:
SubsetBranches = c('AllBranches','TerminalBranches','InternalBranches')
SubsetGenes = c('withoutNd6')  # c('withNd6','WithoutNd6','Nd6Only')
SubsetSubstitutions = c('synonymous') #  SubsetSubstitutions = c('synonymous','non-synonymous','SynAndNons')

##### GENERAL FILTERS OF QUALITY:
###  protein-coding; only with normal Aa (21)
names(Mut)
table(Mut$note)
table(Mut$synonymous)
Pc= Mut[Mut$note == 'normal',] # filter out everything except protein-coding mutations:
table(Pc$ancestral_aa)
VecOfNormalAa = unique(Pc$ancestral_aa); length(VecOfNormalAa)
table(Pc$derived_aa) # Ambiguous, Asn/Asp, Gln/Glu, Leu/Ile => why all noisy AA only among derived??!!!!!!!!!!!!!!!! 
Pc = Pc[Pc$derived_aa %in% VecOfNormalAa,]

#### mutations are only a t g c;  Derive Subst and Context 
# Subst
ExtractThird<-function(x) {unlist(strsplit(x,''))[3]}
Pc$temp1 = apply(as.matrix(Pc$ancestor),1,FUN = ExtractThird); Pc = Pc[Pc$temp1 %in% c('A','T','G','C'),]
Pc$temp2 = apply(as.matrix(Pc$descendant),1,FUN = ExtractThird); Pc = Pc[Pc$temp2 %in% c('A','T','G','C'),]
Pc$Subst = paste(Pc$temp1,Pc$temp2,sep='')
table(Pc$Subst)

# FILTER FOR THE SAME BACKGROUND (Pos1Anc == Pos1Der; Pos5Anc == Pos5Der; ):
# very first (first) and very last (fifth) should be the same (important for codons - we will garantie, that in the codon there is only one substitution)
nrow(Pc) # 292532
ExtractFirst<-function(x) {unlist(strsplit(x,''))[1]}
Pc$Pos1Anc = apply(as.matrix(Pc$ancestor),1,FUN = ExtractFirst);   Pc = Pc[Pc$Pos1Anc %in% c('a','t','g','c'),]
Pc$Pos1Der = apply(as.matrix(Pc$descendant),1,FUN = ExtractFirst); Pc = Pc[Pc$Pos1Der %in% c('a','t','g','c'),]
Pc=Pc[Pc$Pos1Anc == Pc$Pos1Der,]; nrow(Pc)

ExtractFifth<-function(x) {unlist(strsplit(x,''))[5]}
Pc$Pos5Anc = apply(as.matrix(Pc$ancestor),1,FUN = ExtractFifth);   Pc = Pc[Pc$Pos5Anc %in% c('a','t','g','c'),]
Pc$Pos5Der = apply(as.matrix(Pc$descendant),1,FUN = ExtractFifth); Pc = Pc[Pc$Pos5Der %in% c('a','t','g','c'),]
Pc=Pc[Pc$Pos5Anc == Pc$Pos5Der,]; nrow(Pc)

# Context
ExtractSecond<-function(x) {unlist(strsplit(x,''))[2]}
Pc$Pos2Anc = apply(as.matrix(Pc$ancestor),1,FUN = ExtractSecond); Pc = Pc[Pc$Pos2Anc %in% c('a','t','g','c'),]
Pc$Pos2Der = apply(as.matrix(Pc$descendant),1,FUN = ExtractSecond); Pc = Pc[Pc$Pos2Der %in% c('a','t','g','c'),]
Pc=Pc[Pc$Pos2Anc == Pc$Pos2Der,]

ExtractFourth<-function(x) {unlist(strsplit(x,''))[4]}
Pc$Pos4Anc = apply(as.matrix(Pc$ancestor),1,FUN = ExtractFourth); Pc = Pc[Pc$Pos4Anc %in% c('a','t','g','c'),]
Pc$Pos4Der = apply(as.matrix(Pc$descendant),1,FUN = ExtractFourth); Pc = Pc[Pc$Pos4Der %in% c('a','t','g','c'),]
Pc=Pc[Pc$Pos4Anc == Pc$Pos4Der,]

Pc$Context = paste(toupper(Pc$Pos2Anc),Pc$temp1,toupper(Pc$Pos4Anc),sep='')

##### ADD HEAVY CHAIN NOTATION: SubstHeavy +  ContextHeavy (will use it in case of Nd6?):
library(seqinr) # install.packages('seqinr')
Comp <- function(x) {y = paste(comp(s2c(as.character(x))),collapse = ''); return(toupper(y))}
Pc$SubstHeavy = apply(as.matrix(Pc$Subst), 1, FUN = Comp)
RevComp <- function(x) {y = paste(rev(comp(s2c(as.character(x)))),collapse = ''); return(toupper(y))}
Pc$ContextHeavy = apply(as.matrix(Pc$Context), 1, FUN = RevComp)

PcBeforeFilters = Pc

### PREPARE ARTIFICIAL MATRIX WITH 192-COMPONENT CONTEXT:
Context = data.frame(permutations(4, 3, repeats.allowed = TRUE)); names(Context) = c('Before','Middle','After')
Subst = data.frame(permutations(4, 2, repeats.allowed = TRUE)); names(Subst) = c('ref','var')
Context = merge(Context,Subst)
Context[] <- lapply(Context, function(x) (gsub("1", "A", x))); Context[] <- lapply(Context, function(x) (gsub("2", "T", x))); 
Context[] <- lapply(Context, function(x) (gsub("3", "G", x))); Context[] <- lapply(Context, function(x) (gsub("4", "C", x)))
Context = Context[Context$ref != Context$var,];
Context = Context[Context$ref == Context$Middle,]; nrow(Context) # 192
Context$names = paste(Context$ref,Context$var,': ',Context$Before,Context$Middle,Context$After,sep = '')
Context$NumbObs = 0; Artificial192=select(Context,names,NumbObs)
Context = select(Context,names,NumbObs)

##### RUN ALL FILTERS FOR SUBSETS:
MutSpec12Fin = data.frame()
MutSpec192Fin = data.frame()

for (branches in SubsetBranches) # SubsetBranches = c('AllBranches','TerminalBranches','InternalBranches')
{
  if (branches == 'AllBranches')  {Pc = PcBeforeFilters;}
  if (branches == 'TerminalBranches') {Pc=PcBeforeFilters[grepl("_",PcBeforeFilters$second),]; Pc=Pc[grepl("[A-Za-z]",Pc$second),];}  # PC$second should have _ and letters like this 'AF_ET_0044' # both conditions are identical since terminal branch has both: "_" and [A-Za-z]
  if (branches == 'InternalBranches') {Pc=PcBeforeFilters[!grepl("_",PcBeforeFilters$second),]; Pc=Pc[!grepl("[A-Za-z]",Pc$second),];}  # PC$second should have _ and letters like this 'AF_ET_0044' # both conditions are identical since terminal branch has both: "_" and [A-Za-z]
  
for (genes in SubsetGenes) # SubsetGenes = c('withNd6','WithoutNd6','Nd6Only')
{    
  # if (genes == 'withNd6')  {}
  if (genes == 'WithoutNd6')  {Pc = Pc[Pc$gene_info != 'mRNA_ND6',];}    
  # if (genes == 'Nd6Only')  {Pc = Pc[Pc$gene_info == 'mRNA_ND6',];}    

for (subst in SubsetSubstitutions) # SubsetSubstitutions = c('synonymous','non-synonymous','SynAndNons')
{    
  if (subst == 'synonymous')  {Pc = Pc[Pc$synonymous == 'synonymous',]}
  #if (subst == 'non-synonymous')  {Pc = Pc[Pc$synonymous == 'non-synonymous',];}    
  #if (subst == 'SynAndNons')  {}    
    
######### MutSpec12
MutSpec12 = as.data.frame(table(Pc$Subst))
names(MutSpec12)=c('Subst','NumbObs')
MutSpec12 = MutSpec12[order(MutSpec12$Subst),]
MutSpec12$Branches = branches; 
MutSpec12$Genes = genes;
MutSpec12$SubstType = subst; 
MutSpec12Fin = rbind(MutSpec12Fin,MutSpec12)

######### MutSpec192
MutSpec192 = as.data.frame(table(Pc$Subst,Pc$Context)); 
MutSpec192 = MutSpec192[MutSpec192$Freq > 0,] 
names(MutSpec192)=c('Subst','Context','NumbObs')
MutSpec192$names = paste(MutSpec192$Subst,MutSpec192$Context,sep = ': ')
MutSpec192 = select(MutSpec192,names,NumbObs)

MutSpec192 = rbind(MutSpec192,Context)
MutSpec192 = aggregate(MutSpec192$NumbObs, by = list(MutSpec192$names),FUN = sum)
names(MutSpec192)=c('names','NumbObs')
nrow(MutSpec192) # 192!!!
MutSpec192 = MutSpec192[order(MutSpec192$names),]
MutSpec192$Branches = branches; MutSpec192$Genes = genes; MutSpec192$SubstType = subst; 
MutSpec192Fin = rbind(MutSpec192Fin,MutSpec192)

}}}

write.table(MutSpec12Fin,'../../Body/2Derived/02A.GlobalHumanTree.MutSpecsObserved192&12.MutSpecObs12.txt')
write.table(MutSpec192Fin,'../../Body/2Derived/02A.GlobalHumanTree.MutSpecsObserved192&12.MutSpecObs192.txt')
