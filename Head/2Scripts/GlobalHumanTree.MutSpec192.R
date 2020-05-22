rm(list=ls(all=TRUE))

pdf("../../Body/4Figures/GlobalHumanTree.MutSpec192.R.pdf", width = 40)

Mut = read.table('../../Body/2Derived/fulltreeCodons.csv', header = TRUE, sep = ';'); # "../../Body/3Results/fulltreeCodons.csv", sep=";"
names(Mut)
table(Mut$note)
table(Mut$synonymous)
Pc= Mut[Mut$note == 'normal',] # filter out everything except protein-coding mutations:
table(Pc$ancestral_aa)
VecOfNormalAa = unique(Pc$ancestral_aa); length(VecOfNormalAa)
table(Pc$derived_aa) # Ambiguous, Asn/Asp, Gln/Glu, Leu/Ile => why all noisy AA only among derived??!!!!!!!!!!!!!!!!
Pc = Pc[Pc$derived_aa %in% VecOfNormalAa,]

### Derive Subst and Context for all SYN + NONS (identical in ancestor and descendant) and in all mutations only a t g c  
ExtractThird<-function(x) {unlist(strsplit(x,''))[3]}
Pc$temp1 = apply(as.matrix(Pc$ancestor),1,FUN = ExtractThird); Pc = Pc[Pc$temp1 %in% c('A','T','G','C'),]
Pc$temp2 = apply(as.matrix(Pc$descendant),1,FUN = ExtractThird); Pc = Pc[Pc$temp2 %in% c('A','T','G','C'),]
Pc$Subst = paste(Pc$temp1,Pc$temp2,sep='')
table(Pc$Subst)

ExtractSecond<-function(x) {unlist(strsplit(x,''))[2]}
Pc$temp3 = apply(as.matrix(Pc$ancestor),1,FUN = ExtractSecond); Pc = Pc[Pc$temp3 %in% c('a','t','g','c'),]
Pc$temp4 = apply(as.matrix(Pc$descendant),1,FUN = ExtractSecond); Pc = Pc[Pc$temp4 %in% c('a','t','g','c'),]
Pc=Pc[Pc$temp3 == Pc$temp4,]

ExtractFourth<-function(x) {unlist(strsplit(x,''))[4]}
Pc$temp5 = apply(as.matrix(Pc$ancestor),1,FUN = ExtractFourth); Pc = Pc[Pc$temp5 %in% c('a','t','g','c'),]
Pc$temp6 = apply(as.matrix(Pc$descendant),1,FUN = ExtractFourth); Pc = Pc[Pc$temp6 %in% c('a','t','g','c'),]
Pc=Pc[Pc$temp5 == Pc$temp6,]

Pc$Context = paste(toupper(Pc$temp3),Pc$temp1,toupper(Pc$temp5),sep='')

#### SWITCH TO HEAVY CHAIN:

library(seqinr) # install.packages('seqinr')
#comp(seq) # gives complement
#rev(comp(seq)) # gi

Comp <- function(x) {y = paste(comp(s2c(as.character(x))),collapse = ''); return(toupper(y))}
Pc$SubstHeavy = apply(as.matrix(Pc$Subst), 1, FUN = Comp)

RevComp <- function(x) {y = paste(rev(comp(s2c(as.character(x)))),collapse = ''); return(toupper(y))}
Pc$ContextHeavy = apply(as.matrix(Pc$Context), 1, FUN = RevComp)

### create artificial 192 matrix with zeroes
VecOfSubst = c('AT','AG','AC','TA','TG','TC','CA','CG','CT','GA','GC','GT'); length(unique(VecOfSubst)) # 12
VecOfFirsts = c('A','T','G','C')
VecOfThirds = c('A','T','G','C')
VecOfContext = ''
Freq=0
flag = 0
for (i in 1:length(VecOfSubst))
{ # i = 1
  Subst = VecOfSubst[i]
  From = unlist(strsplit(VecOfSubst[i],''))[1]
  for (first in 1:4)
  { # first = 1
    for (third in 1:4)
    { # third = 1
      Context=paste(VecOfFirsts[first],From,VecOfThirds[third],sep='')
      NewLine = data.frame(Subst,Context,Freq)
      if (flag > 0) {Final = rbind(Final,NewLine)}
      if (flag == 0) {Final = NewLine; flag = 1}
    }
  }
}

# get real 192 MutSpec 4*12*4 

MutSpec192 = as.data.frame(table(Pc$SubstHeavy,Pc$ContextHeavy))
MutSpec192 = MutSpec192[MutSpec192$Freq > 0,] 
nrow(MutSpec192) # 121
names(MutSpec192)=c('Subst','Context','Freq')

# rbind both and keep zeroes

MutSpec192 = rbind(MutSpec192,Final)
MutSpec192 = aggregate(MutSpec192$Freq, by = list(MutSpec192$Subst,MutSpec192$Context),FUN = sum)
names(MutSpec192)=c('Subst','Context','Freq')
nrow(MutSpec192) # 192!!!

MutSpec192$Subst = as.character(MutSpec192$Subst)
MutSpec192$Context = as.character(MutSpec192$Context)
MutSpec192 = MutSpec192[order(MutSpec192$Subst,MutSpec192$Context),]
MutSpec192$names = paste(MutSpec192$Subst,': ',MutSpec192$Context,sep='')

nrow(Pc) # 298615
barplot(MutSpec192$Freq,names=MutSpec192$names,las = 2,cex.names = 0.4,col=c(rgb(1,0.0,0.0,0.5),rgb(0.0,1,0.0,0.5),rgb(0.0,0.0,1,0.5),rgb(0.0,1,1,0.5)), main = 'MutSpec192 (4*12*4) for all coding mtDNA substitutions (N=298615), heavy chain notation')

### ONLY SYN
nrow(Pc[Pc$synonymous == 'synonymous',]) # 201102
MutSpec192 = as.data.frame(table(Pc[Pc$synonymous == 'synonymous',]$SubstHeavy,Pc[Pc$synonymous == 'synonymous',]$ContextHeavy))
MutSpec192 = MutSpec192[MutSpec192$Freq > 0,] 
nrow(MutSpec192)
names(MutSpec192)=c('Subst','Context','Freq')
MutSpec192 = rbind(MutSpec192,Final)
MutSpec192 = aggregate(MutSpec192$Freq, by = list(MutSpec192$Subst,MutSpec192$Context),FUN = sum)
names(MutSpec192)=c('Subst','Context','Freq')
nrow(MutSpec192) # 192!!!

MutSpec192$Subst = as.character(MutSpec192$Subst)
MutSpec192$Context = as.character(MutSpec192$Context)
MutSpec192 = MutSpec192[order(MutSpec192$Subst,MutSpec192$Context),]
MutSpec192$names = paste(MutSpec192$Subst,': ',MutSpec192$Context,sep='')

barplot(MutSpec192$Freq,names=MutSpec192$names,las = 2,cex.names = 0.4,col=c(rgb(1,0.0,0.0,0.5),rgb(0.0,1,0.0,0.5),rgb(0.0,0.0,1,0.5),rgb(0.0,1,1,0.5)), main = 'MutSpec192 (4*12*4) for all synonymous mtDNA human substitutions (N=201102), heavy chain notation')

### ONLY NONS
nrow(Pc[Pc$synonymous == 'non-synonymous',]) # 97513
MutSpec192 = as.data.frame(table(Pc[Pc$synonymous == 'non-synonymous',]$SubstHeavy,Pc[Pc$synonymous == 'non-synonymous',]$ContextHeavy))
MutSpec192 = MutSpec192[MutSpec192$Freq > 0,] 
nrow(MutSpec192)
names(MutSpec192)=c('Subst','Context','Freq')
MutSpec192 = rbind(MutSpec192,Final)
MutSpec192 = aggregate(MutSpec192$Freq, by = list(MutSpec192$Subst,MutSpec192$Context),FUN = sum)
names(MutSpec192)=c('Subst','Context','Freq')
nrow(MutSpec192) # 192!!!

MutSpec192$Subst = as.character(MutSpec192$Subst)
MutSpec192$Context = as.character(MutSpec192$Context)
MutSpec192 = MutSpec192[order(MutSpec192$Subst,MutSpec192$Context),]
MutSpec192$names = paste(MutSpec192$Subst,': ',MutSpec192$Context,sep='')

barplot(MutSpec192$Freq,names=MutSpec192$names,las = 2,cex.names = 0.4,col=c(rgb(1,0.0,0.0,0.5),rgb(0.0,1,0.0,0.5),rgb(0.0,0.0,1,0.5),rgb(0.0,1,1,0.5)), main = 'MutSpec192 (4*12*4) for all non-synonymous mtDNA human substitutions (N=97513), heavy chain notation')

dev.off()

