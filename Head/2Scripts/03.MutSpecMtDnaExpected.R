#### GENERATE EXPECTED MutSpec, COMPARE it with observed and plot everything

rm(list=ls(all=TRUE))
library('dplyr')
library('gtools')

######### 1 GENERATE EXPECTED MUT SPECS 12 and 192 for different genes and substitutions: 
Mut = read.table("../../Body/1Raw/FromCovidProject/HumanMtDnaRefSeq.fasta.ExtensiveMut.vcf.ann")
names(Mut) = c('Chr','Position','Id','AncestralNuc','DerivedNuc','Qual','Filter','Info')
Mut$Subst = paste(Mut$AncestralNuc,Mut$DerivedNuc,sep = '')
Mut$ContextLeft = c(c('NA','NA','NA'),Mut$AncestralNuc[-c(nrow(Mut)-2,nrow(Mut)-1,nrow(Mut))])
Mut$ContextRight = c(Mut$AncestralNuc[-c(1,2,3)],c('NA','NA','NA'))
Mut$Context = paste(Mut$ContextLeft,Mut$AncestralNuc,Mut$ContextRight,sep = '')
Mut = Mut[Mut$ContextLeft !=  'NA' & Mut$ContextRight != 'NA',] 
extractGene<-function(x) {unlist(strsplit(x,"\\|"))[4];}; 
Mut$ProtCodGene = apply(as.matrix(Mut$Info),1,FUN=extractGene); table(Mut$ProtCodGene)
extractSubstType<-function(x) {unlist(strsplit(x,"\\|"))[2];} # synonymous_variant
Mut$SubstType = apply(as.matrix(Mut$Info),1,FUN=extractSubstType); table(Mut$SubstType)

######### 2  DEFINE MAIN SUBSETS OF MUTATIONS:
SubsetSubstitutions = c('synonymous','AllExceptStopGains') #  SubsetSubstitutions = c('synonymous','non-synonymous','SynAndNons','AllExceptStopGains')
SubsetGenes = c('withoutNd6')  # c('withNd6','withoutNd6','Nd6Only')

######### 3 PREPARE ARTIFICIAL MATRIX WITH 192-COMPONENT CONTEXT:
Context = data.frame(permutations(4, 3, repeats.allowed = TRUE)); names(Context) = c('Before','Middle','After')
Subst = data.frame(permutations(4, 2, repeats.allowed = TRUE)); names(Subst) = c('ref','var')
Context = merge(Context,Subst)
Context[] <- lapply(Context, function(x) (gsub("1", "A", x))); Context[] <- lapply(Context, function(x) (gsub("2", "T", x))); 
Context[] <- lapply(Context, function(x) (gsub("3", "G", x))); Context[] <- lapply(Context, function(x) (gsub("4", "C", x)))
Context = Context[Context$ref != Context$var,];
Context = Context[Context$ref == Context$Middle,]; nrow(Context) # 192
Context$names = paste(Context$ref,Context$var,': ',Context$Before,Context$Middle,Context$After,sep = '')
Context$NumbExp = 0; Artificial192=select(Context,names,NumbExp)
Context = select(Context,names,NumbExp)

##### RUN ALL FILTERS FOR SUBSETS:
PcBeforeFilters = Mut

MutSpec12Fin = data.frame()
MutSpec192Fin = data.frame()

#### inside the loops very important that internal loop is just one condition - withoutND6 - otherwise will be mistake!!!!! 
### loops are not good in filtering
for (subst in SubsetSubstitutions) # SubsetSubstitutions = c('synonymous','non-synonymous','SynAndNons','AllExceptStopGains')
{   # subst =  'synonymous' # subst = 'AllExceptStopGains'
  if (subst == 'synonymous')  {Pc = PcBeforeFilters[PcBeforeFilters$SubstType == 'synonymous_variant',]}
  if (subst == 'non-synonymous') {}
  if (subst == 'SynAndNons') {}
  if (subst == 'AllExceptStopGains')  {Pc = PcBeforeFilters[PcBeforeFilters$SubstType != 'stop_gained' & PcBeforeFilters$SubstType != 'stop_gained&splice_region_variant',]}
  
for (genes in SubsetGenes) # SubsetGenes = c('withNd6','withoutNd6','Nd6Only')
{    # genes = 'withoutNd6'
  # if (genes == 'withNd6')     {Pc = Pc}
  if (genes == 'withoutNd6')  {Pc = Pc[Pc$ProtCodGene != 'ND6',];} # table(PcBeforeFilters$ProtCodGene)
  # if (genes == 'Nd6Only')     {Pc = Pc[Pc$ProtCodGene == 'ND6',];}    
  
######### MutSpec12
MutSpec12 = as.data.frame(table(Pc$Subst))
names(MutSpec12)=c('Subst','NumbExp')
MutSpec12 = MutSpec12[order(MutSpec12$Subst),]
MutSpec12$Genes = genes;
MutSpec12$SubstType = subst; 
MutSpec12Fin = rbind(MutSpec12Fin,MutSpec12)
      
######### MutSpec192
MutSpec192 = as.data.frame(table(Pc$Subst,Pc$Context)); 
MutSpec192 = MutSpec192[MutSpec192$Freq > 0,] 
names(MutSpec192)=c('Subst','Context','NumbExp')
MutSpec192$names = paste(MutSpec192$Subst,MutSpec192$Context,sep = ': ')
MutSpec192 = select(MutSpec192,names,NumbExp)
      
MutSpec192 = rbind(MutSpec192,Context)
MutSpec192 = aggregate(MutSpec192$NumbExp, by = list(MutSpec192$names),FUN = sum)
names(MutSpec192)=c('names','NumbExp')
nrow(MutSpec192) # 192!!!
MutSpec192 = MutSpec192[order(MutSpec192$names),]
MutSpec192$Genes = genes; MutSpec192$SubstType = subst; 
MutSpec192Fin = rbind(MutSpec192Fin,MutSpec192)
      
}}

write.table(MutSpec12Fin,'../../Body/2Derived/03.MutSpecMtDnaExpected.MutSpecExp12.txt')
write.table(MutSpec192Fin,'../../Body/2Derived/03.MutSpecMtDnaExpected.MutSpecExp192.txt')

# old outputs: 
#write.table(MutSpec192,"../../Body/3Results/MutSpec192ObsToExpSynNoNd6.txt", row.names = FALSE, quote = FALSE, sep = '\t')
#write.table(MutSpec12,"../../Body/3Results/MutSpec12ObsToExpSynNoNd6.txt", row.names = FALSE, quote = FALSE, sep = '\t')

