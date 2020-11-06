rm(list=ls(all=TRUE)) 

#indata <- read.csv("/home/alima/arrr/fulltreeCodons.csv", header = TRUE, sep = ";") #читаю файл 
indata = read.csv("../../Body/2Derived/fulltreeCodons.csv", header = TRUE, sep = ";")

#### FILTER OUT: 
data <- subset(indata, synonymous == "non-synonymous" & derived_aa != "Ambiguous"& gene_info != "mRNA_ND6")  #убираю лишнее

## ADD FILTER OF BACKGROUND (THE SAME NEIGHBOR NUCLEOTIDES)


## DERIVE SUBSTITUTION MATRIX

data$FromTo = paste(data$ancestral_aa,data$derived_aa,sep = '>')

FromTo = data.frame(table(data$FromTo))
names(FromTo) = c('FromAncestralToDerived', 'NumberOfEvents')
nrow(FromTo) # 236, but totally there are 21*21 possibilities
AllAa = data.frame(unique(data$ancestral_aa)); nrow(AllAa) # including stop

# Make it into both directions: A1 A2 

#### we expect majority of ratios > 0
## 1st nucleotide, T>C
FromTo[FromTo$FromAncestralToDerived == 'Phe>Leu',]$NumberOfEvents / FromTo[FromTo$FromAncestralToDerived == 'Leu>Phe',]$NumberOfEvents
FromTo[FromTo$FromAncestralToDerived == 'Ser>Pro',]$NumberOfEvents / FromTo[FromTo$FromAncestralToDerived == 'Pro>Ser',]$NumberOfEvents
FromTo[FromTo$FromAncestralToDerived == 'Tyr>His',]$NumberOfEvents / FromTo[FromTo$FromAncestralToDerived == 'His>Tyr',]$NumberOfEvents
FromTo[FromTo$FromAncestralToDerived == 'Cys>Arg',]$NumberOfEvents / FromTo[FromTo$FromAncestralToDerived == 'Arg>Cys',]$NumberOfEvents
FromTo[FromTo$FromAncestralToDerived == 'Trp>Arg',]$NumberOfEvents / FromTo[FromTo$FromAncestralToDerived == 'Arg>Trp',]$NumberOfEvents
FromTo[FromTo$FromAncestralToDerived == 'Trp>Arg',]$NumberOfEvents / FromTo[FromTo$FromAncestralToDerived == 'Arg>Trp',]$NumberOfEvents
FromTo[FromTo$FromAncestralToDerived == 'Stop>Gln',]$NumberOfEvents / FromTo[FromTo$FromAncestralToDerived == 'Gln>Stop',]$NumberOfEvents 

## 2nd nucleotide, T>C = we expect ratios > 0:

#### make this dataset automatically:
FromTo$From = gsub(">.*",'',FromTo$FromAncestralToDerived)
FromTo$To = gsub(".*>",'',FromTo$FromAncestralToDerived)

# Merge with Alima file and obtain column: expected: 1 or 0

