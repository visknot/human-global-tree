rm(list=ls(all=TRUE))
########

library(readr)
df <- read_delim("~/mitoclub/human-global-tree/Body/1Raw/CancerDataFromCampbell/mtDNA_snv_Oct2016.txt", 
                                "\t", escape_double = FALSE, trim_ws = TRUE)


# ALL Total number of mutations (rows in table) - 7611
nrow(df)

########
# ALL - All mutations of a certain type
# ALL. A>G - 481
nrow(df[which(df$ref== "A" & df$var == "G"),])
# ALL. A>T - 72
nrow(df[which(df$ref== "A" & df$var == "T"),])
# ALL. A>C - 91
nrow(df[which(df$ref== "A" & df$var == "C"),])
# ALL. T>G - 27
nrow(df[which(df$ref== "T" & df$var == "G"),])
# ALL. T>C - 2120
nrow(df[which(df$ref== "T" & df$var == "C"),])
# ALL. T>A - 29
nrow(df[which(df$ref== "T" & df$var == "A"),])
# ALL. G>A - 3663
nrow(df[which(df$ref== "G" & df$var == "A"),])
# ALL. G>T - 49
nrow(df[which(df$ref== "G" & df$var == "T"),])
# ALL. G>C - 119
nrow(df[which(df$ref== "G" & df$var == "C"),])
# ALL. C>A - 167
nrow(df[which(df$ref== "C" & df$var == "A"),])
# ALL. C>T - 771
nrow(df[which(df$ref== "C" & df$var == "T"),])
# ALL. C>G - 22
nrow(df[which(df$ref== "C" & df$var == "G"),])

# #ALL Nodloop 6628
#> summary(nodloop$position)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#577    3841    7983    8132   12354   16023 
#ref var position  (16024..16569,1..576)
nodloop <- mtDNA[which(mtDNA$position > 576 & mtDNA$position < 16024),]
nrow(mtDNA[which(mtDNA$position > 576 & mtDNA$position < 16024),])

# NODLOOP. A>G - 350
nrow(nodloop[which(nodloop$ref== "A" & nodloop$var == "G"),])
# NODLOOP. A>T - 65
nrow(nodloop[which(nodloop$ref== "A" & nodloop$var == "T"),])
# NODLOOP. A>C - 80
nrow(nodloop[which(nodloop$ref== "A" & nodloop$var == "C"),])
# NODLOOP. T>G - 18
nrow(nodloop[which(nodloop$ref== "T" & nodloop$var == "G"),])
# NODLOOP. T>C - 1884
nrow(nodloop[which(nodloop$ref== "T" & nodloop$var == "C"),])
# NODLOOP. T>A - 23
nrow(nodloop[which(nodloop$ref== "T" & nodloop$var == "A"),])
# NODLOOP. G>A - 1452
nrow(nodloop[which(nodloop$ref== "G" & nodloop$var == "A"),])
# NODLOOP. G>T - 32
nrow(nodloop[which(nodloop$ref== "G" & nodloop$var == "T"),])
# NODLOOP. G>C - 112
nrow(nodloop[which(nodloop$ref== "G" & nodloop$var == "C"),])
# NODLOOP. C>A - 125
nrow(nodloop[which(nodloop$ref== "C" & nodloop$var == "A"),])
# NODLOOP. C>T - 471
nrow(nodloop[which(nodloop$ref== "C" & nodloop$var == "T"),])
# NODLOOP. C>G - 16
nrow(nodloop[which(nodloop$ref== "C" & nodloop$var == "G"),])


#Synonymall 874
synonymall = mtDNA[which(mtDNA$is_silent == 1),]
nrow(nodloop[which(nodloop$is_silent == 1),])

# SYNONYMALL. A>G - 91
nrow(synonymall[which(synonymall$ref== "A" & synonymall$var == "G"),])
# SYNONYMALL. A>T - 11
nrow(synonymall[which(synonymall$ref== "A" & synonymall$var == "T"),])
# SYNONYMALL. A>C - 12
nrow(synonymall[which(synonymall$ref== "A" & synonymall$var == "C"),])
# SYNONYMALL. T>G - 3
nrow(synonymall[which(synonymall$ref== "T" & synonymall$var == "G"),])
# SYNONYMALL. T>C - 262
nrow(synonymall[which(synonymall$ref== "T" & synonymall$var == "C"),])
# SYNONYMALL. T>A - 1
nrow(synonymall[which(synonymall$ref== "T" & synonymall$var == "A"),])
# SYNONYMALL. G>A - 255
nrow(synonymall[which(synonymall$ref== "G" & synonymall$var == "A"),])
# SYNONYMALL. G>T - 0
nrow(synonymall[which(synonymall$ref== "G" & synonymall$var == "T"),])
# SYNONYMALL. G>C - 4
nrow(synonymall[which(synonymall$ref== "G" & synonymall$var == "C"),])
# SYNONYMALL. C>A - 20
nrow(synonymall[which(synonymall$ref== "C" & synonymall$var == "A"),])
# SYNONYMALL. C>T - 211
nrow(synonymall[which(synonymall$ref== "C" & synonymall$var == "T"),])
# SYNONYMALL. C>G - 4
nrow(synonymall[which(synonymall$ref== "C" & synonymall$var == "G"),])



########
# Adding preceeding_nucl from expected
# Adding preceeding_nucl
datalist = c()
poslist = c()
count = 1
for (i in seq(4, nrow(df), by = 3)) {
  preceeding = i-1
  end = i + 2
  temp = unlist(strsplit(substr(as.character(df[preceeding,]), 1, 25), "\t"))
  preceeding_nucleotide = temp[4]
  print(temp)
  print(count)
  res = c(rep(preceeding_nucleotide, 3))
  datalist = c(datalist, res)
  poslist= c(poslist, c(rep(count, 3)))
  count = count+1
}
}
# Append to the top
last_nucl = c(rep("G", 3))
last_pos = c(rep(16567, 3))
datalist = c(last_nucl, datalist)
poslist = c(poslist, last_pos)

df_context = df
df_context['context'] = datalist
df_context['position'] = poslist

df_merged =  merge(df_context, df, by = "position")
df_merged['preceeding'] = df_merged$context.y

########
# ALL. A>G - 
nrow(df_merged[which(df_merged$ref== "A" & df_merged$var == "G" & df_merged$context == "A"),])/3
nrow(df_merged[which(df_merged$ref== "A" & df_merged$var == "G" & df_merged$context == "T"),])/3
nrow(df_merged[which(df_merged$ref== "A" & df_merged$var == "G" & df_merged$context == "G"),])/3
nrow(df_merged[which(df_merged$ref== "A" & df_merged$var == "G" & df_merged$context == "C"),])/3
# ALL. A>T - 5124
nrow(df_merged[which(df_merged$ref== "A" & df_merged$var == "T" & df_merged$context == "A"),])/3
nrow(df_merged[which(df_merged$ref== "A" & df_merged$var == "T" & df_merged$context == "T"),])/3
nrow(df_merged[which(df_merged$ref== "A" & df_merged$var == "T" & df_merged$context == "G"),])/3
nrow(df_merged[which(df_merged$ref== "A" & df_merged$var == "T" & df_merged$context == "C"),])/3
# ALL. A>C - 5124
nrow(df_merged[which(df_merged$ref== "A" & df_merged$var == "C" & df_merged$context == "A"),])/3
nrow(df_merged[which(df_merged$ref== "A" & df_merged$var == "C" & df_merged$context == "T"),])/3
nrow(df_merged[which(df_merged$ref== "A" & df_merged$var == "C" & df_merged$context == "G"),])/3
nrow(df_merged[which(df_merged$ref== "A" & df_merged$var == "C" & df_merged$context == "C"),])/3
# ALL. T>G - 4094
nrow(df_merged[which(df_merged$ref== "T" & df_merged$var == "G" & df_merged$context == "A"),])/3
nrow(df_merged[which(df_merged$ref== "T" & df_merged$var == "G" & df_merged$context == "T"),])/3
nrow(df_merged[which(df_merged$ref== "T" & df_merged$var == "G" & df_merged$context == "G"),])/3
nrow(df_merged[which(df_merged$ref== "T" & df_merged$var == "G" & df_merged$context == "C"),])/3
# ALL. T>C - 4094
nrow(df_merged[which(df_merged$ref== "T" & df_merged$var == "C" & df_merged$context == "A"),])/3
nrow(df_merged[which(df_merged$ref== "T" & df_merged$var == "C" & df_merged$context == "T"),])/3
nrow(df_merged[which(df_merged$ref== "T" & df_merged$var == "C" & df_merged$context == "G"),])/3
nrow(df_merged[which(df_merged$ref== "T" & df_merged$var == "C" & df_merged$context == "C"),])/3
# ALL. T>A - 4094
nrow(df_merged[which(df_merged$ref== "T" & df_merged$var == "A" & df_merged$context == "A"),])/3
nrow(df_merged[which(df_merged$ref== "T" & df_merged$var == "A" & df_merged$context == "T"),])/3
nrow(df_merged[which(df_merged$ref== "T" & df_merged$var == "A" & df_merged$context == "G"),])/3
nrow(df_merged[which(df_merged$ref== "T" & df_merged$var == "A" & df_merged$context == "C"),])/3
# ALL. G>A - 2169
nrow(df_merged[which(df_merged$ref== "G" & df_merged$var == "A" & df_merged$context == "A"),])/3
nrow(df_merged[which(df_merged$ref== "G" & df_merged$var == "A" & df_merged$context == "T"),])/3
nrow(df_merged[which(df_merged$ref== "G" & df_merged$var == "A" & df_merged$context == "G"),])/3
nrow(df_merged[which(df_merged$ref== "G" & df_merged$var == "A" & df_merged$context == "C"),])/3
# ALL. G>T - 2169
nrow(df_merged[which(df_merged$ref== "G" & df_merged$var == "T" & df_merged$context == "A"),])/3
nrow(df_merged[which(df_merged$ref== "G" & df_merged$var == "T" & df_merged$context == "T"),])/3
nrow(df_merged[which(df_merged$ref== "G" & df_merged$var == "T" & df_merged$context == "G"),])/3
nrow(df_merged[which(df_merged$ref== "G" & df_merged$var == "T" & df_merged$context == "C"),])/3
# ALL. G>C - 2169
nrow(df_merged[which(df_merged$ref== "G" & df_merged$var == "C" & df_merged$context == "A"),])/3
nrow(df_merged[which(df_merged$ref== "G" & df_merged$var == "C" & df_merged$context == "T"),])/3
nrow(df_merged[which(df_merged$ref== "G" & df_merged$var == "C" & df_merged$context == "G"),])/3
nrow(df_merged[which(df_merged$ref== "G" & df_merged$var == "C" & df_merged$context == "C"),])/3
# ALL. C>A - 5182
nrow(df_merged[which(df_merged$ref== "C" & df_merged$var == "A" & df_merged$context == "A"),])/3
nrow(df_merged[which(df_merged$ref== "C" & df_merged$var == "A" & df_merged$context == "T"),])/3
nrow(df_merged[which(df_merged$ref== "C" & df_merged$var == "A" & df_merged$context == "G"),])/3
nrow(df_merged[which(df_merged$ref== "C" & df_merged$var == "A" & df_merged$context == "C"),])/3
# ALL. C>T - 5182
nrow(df_merged[which(df_merged$ref== "C" & df_merged$var == "T" & df_merged$context == "A"),])/3
nrow(df_merged[which(df_merged$ref== "C" & df_merged$var == "T" & df_merged$context == "T"),])/3
nrow(df_merged[which(df_merged$ref== "C" & df_merged$var == "T" & df_merged$context == "G"),])/3
nrow(df_merged[which(df_merged$ref== "C" & df_merged$var == "T" & df_merged$context == "C"),])/3
# ALL. C>G - 5182
nrow(df_merged[which(df_merged$ref== "C" & df_merged$var == "G" & df_merged$context == "A"),])/3
nrow(df_merged[which(df_merged$ref== "C" & df_merged$var == "G" & df_merged$context == "T"),])/3
nrow(df_merged[which(df_merged$ref== "C" & df_merged$var == "G" & df_merged$context == "G"),])/3
nrow(df_merged[which(df_merged$ref== "C" & df_merged$var == "G" & df_merged$context == "C"),])/3


########
# NODLOOP - Everything outside of D-loop. D-loop coordinates from NC_012920 (16024..16569,1..576)
# Rows from df[1729,] to df[48072] 
nodloop <- df_merged[which(df_merged$position > 576 & df_merged$position < 16024),]
# NODLOOP. A>G - 4785
# NODLOOP. A>G - 
nrow(nodloop[which(nodloop$ref== "A" & nodloop$var == "G" & nodloop$context == "A"),])/3
nrow(nodloop[which(nodloop$ref== "A" & nodloop$var == "G" & nodloop$context == "T"),])/3
nrow(nodloop[which(nodloop$ref== "A" & nodloop$var == "G" & nodloop$context == "G"),])/3
nrow(nodloop[which(nodloop$ref== "A" & nodloop$var == "G" & nodloop$context == "C"),])/3
# NODLOOP. A>T - 5124
nrow(nodloop[which(nodloop$ref== "A" & nodloop$var == "T" & nodloop$context == "A"),])/3
nrow(nodloop[which(nodloop$ref== "A" & nodloop$var == "T" & nodloop$context == "T"),])/3
nrow(nodloop[which(nodloop$ref== "A" & nodloop$var == "T" & nodloop$context == "G"),])/3
nrow(nodloop[which(nodloop$ref== "A" & nodloop$var == "T" & nodloop$context == "C"),])/3
# NODLOOP. A>C - 5124
nrow(nodloop[which(nodloop$ref== "A" & nodloop$var == "C" & nodloop$context == "A"),])/3
nrow(nodloop[which(nodloop$ref== "A" & nodloop$var == "C" & nodloop$context == "T"),])/3
nrow(nodloop[which(nodloop$ref== "A" & nodloop$var == "C" & nodloop$context == "G"),])/3
nrow(nodloop[which(nodloop$ref== "A" & nodloop$var == "C" & nodloop$context == "C"),])/3
# NODLOOP. T>G - 4094
nrow(nodloop[which(nodloop$ref== "T" & nodloop$var == "G" & nodloop$context == "A"),])/3
nrow(nodloop[which(nodloop$ref== "T" & nodloop$var == "G" & nodloop$context == "T"),])/3
nrow(nodloop[which(nodloop$ref== "T" & nodloop$var == "G" & nodloop$context == "G"),])/3
nrow(nodloop[which(nodloop$ref== "T" & nodloop$var == "G" & nodloop$context == "C"),])/3
# NODLOOP. T>C - 4094
nrow(nodloop[which(nodloop$ref== "T" & nodloop$var == "C" & nodloop$context == "A"),])/3
nrow(nodloop[which(nodloop$ref== "T" & nodloop$var == "C" & nodloop$context == "T"),])/3
nrow(nodloop[which(nodloop$ref== "T" & nodloop$var == "C" & nodloop$context == "G"),])/3
nrow(nodloop[which(nodloop$ref== "T" & nodloop$var == "C" & nodloop$context == "C"),])/3
# NODLOOP. T>A - 4094
nrow(nodloop[which(nodloop$ref== "T" & nodloop$var == "A" & nodloop$context == "A"),])/3
nrow(nodloop[which(nodloop$ref== "T" & nodloop$var == "A" & nodloop$context == "T"),])/3
nrow(nodloop[which(nodloop$ref== "T" & nodloop$var == "A" & nodloop$context == "G"),])/3
nrow(nodloop[which(nodloop$ref== "T" & nodloop$var == "A" & nodloop$context == "C"),])/3
# NODLOOP. G>A - 2169
nrow(nodloop[which(nodloop$ref== "G" & nodloop$var == "A" & nodloop$context == "A"),])/3
nrow(nodloop[which(nodloop$ref== "G" & nodloop$var == "A" & nodloop$context == "T"),])/3
nrow(nodloop[which(nodloop$ref== "G" & nodloop$var == "A" & nodloop$context == "G"),])/3
nrow(nodloop[which(nodloop$ref== "G" & nodloop$var == "A" & nodloop$context == "C"),])/3
# NODLOOP. G>T - 2169
nrow(nodloop[which(nodloop$ref== "G" & nodloop$var == "T" & nodloop$context == "A"),])/3
nrow(nodloop[which(nodloop$ref== "G" & nodloop$var == "T" & nodloop$context == "T"),])/3
nrow(nodloop[which(nodloop$ref== "G" & nodloop$var == "T" & nodloop$context == "G"),])/3
nrow(nodloop[which(nodloop$ref== "G" & nodloop$var == "T" & nodloop$context == "C"),])/3
# NODLOOP. G>C - 2169
nrow(nodloop[which(nodloop$ref== "G" & nodloop$var == "C" & nodloop$context == "A"),])/3
nrow(nodloop[which(nodloop$ref== "G" & nodloop$var == "C" & nodloop$context == "T"),])/3
nrow(nodloop[which(nodloop$ref== "G" & nodloop$var == "C" & nodloop$context == "G"),])/3
nrow(nodloop[which(nodloop$ref== "G" & nodloop$var == "C" & nodloop$context == "C"),])/3
# NODLOOP. C>A - 5182
nrow(nodloop[which(nodloop$ref== "C" & nodloop$var == "A" & nodloop$context == "A"),])/3
nrow(nodloop[which(nodloop$ref== "C" & nodloop$var == "A" & nodloop$context == "T"),])/3
nrow(nodloop[which(nodloop$ref== "C" & nodloop$var == "A" & nodloop$context == "G"),])/3
nrow(nodloop[which(nodloop$ref== "C" & nodloop$var == "A" & nodloop$context == "C"),])/3
# NODLOOP. C>T - 5182
nrow(nodloop[which(nodloop$ref== "C" & nodloop$var == "T" & nodloop$context == "A"),])/3
nrow(nodloop[which(nodloop$ref== "C" & nodloop$var == "T" & nodloop$context == "T"),])/3
nrow(nodloop[which(nodloop$ref== "C" & nodloop$var == "T" & nodloop$context == "G"),])/3
nrow(nodloop[which(nodloop$ref== "C" & nodloop$var == "T" & nodloop$context == "C"),])/3
# NODLOOP. C>G - 5182
nrow(nodloop[which(nodloop$ref== "C" & nodloop$var == "G" & nodloop$context == "A"),])/3
nrow(nodloop[which(nodloop$ref== "C" & nodloop$var == "G" & nodloop$context == "T"),])/3
nrow(nodloop[which(nodloop$ref== "C" & nodloop$var == "G" & nodloop$context == "G"),])/3
nrow(nodloop[which(nodloop$ref== "C" & nodloop$var == "G" & nodloop$context == "C"),])/3


########
# SYNONYMALL - All synonymous variants - 8678
# SYNONYMALL. A>G - 
nrow(nodloop[which(nodloop$ref== "A" & nodloop$var == "G" & nodloop$context == "A" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "A" & nodloop$var == "G" & nodloop$context == "T" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "A" & nodloop$var == "G" & nodloop$context == "G" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "A" & nodloop$var == "G" & nodloop$context == "C" & nodloop$is_silent == 1),])/3
# SYNONYMALL. A>T - 5124
nrow(nodloop[which(nodloop$ref== "A" & nodloop$var == "T" & nodloop$context == "A" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "A" & nodloop$var == "T" & nodloop$context == "T" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "A" & nodloop$var == "T" & nodloop$context == "G" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "A" & nodloop$var == "T" & nodloop$context == "C" & nodloop$is_silent == 1),])/3
# SYNONYMALL. A>C - 5124
nrow(nodloop[which(nodloop$ref== "A" & nodloop$var == "C" & nodloop$context == "A" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "A" & nodloop$var == "C" & nodloop$context == "T" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "A" & nodloop$var == "C" & nodloop$context == "G" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "A" & nodloop$var == "C" & nodloop$context == "C" & nodloop$is_silent == 1),])/3
# SYNONYMALL. T>G - 4094
nrow(nodloop[which(nodloop$ref== "T" & nodloop$var == "G" & nodloop$context == "A" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "T" & nodloop$var == "G" & nodloop$context == "T" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "T" & nodloop$var == "G" & nodloop$context == "G" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "T" & nodloop$var == "G" & nodloop$context == "C" & nodloop$is_silent == 1),])/3
# SYNONYMALL. T>C - 4094
nrow(nodloop[which(nodloop$ref== "T" & nodloop$var == "C" & nodloop$context == "A" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "T" & nodloop$var == "C" & nodloop$context == "T" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "T" & nodloop$var == "C" & nodloop$context == "G" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "T" & nodloop$var == "C" & nodloop$context == "C" & nodloop$is_silent == 1),])/3
# SYNONYMALL. T>A - 4094
nrow(nodloop[which(nodloop$ref== "T" & nodloop$var == "A" & nodloop$context == "A" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "T" & nodloop$var == "A" & nodloop$context == "T" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "T" & nodloop$var == "A" & nodloop$context == "G" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "T" & nodloop$var == "A" & nodloop$context == "C" & nodloop$is_silent == 1),])/3
# SYNONYMALL. G>A - 2169
nrow(nodloop[which(nodloop$ref== "G" & nodloop$var == "A" & nodloop$context == "A" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "G" & nodloop$var == "A" & nodloop$context == "T" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "G" & nodloop$var == "A" & nodloop$context == "G" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "G" & nodloop$var == "A" & nodloop$context == "C" & nodloop$is_silent == 1),])/3
# SYNONYMALL. G>T - 2169
nrow(nodloop[which(nodloop$ref== "G" & nodloop$var == "T" & nodloop$context == "A" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "G" & nodloop$var == "T" & nodloop$context == "T" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "G" & nodloop$var == "T" & nodloop$context == "G" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "G" & nodloop$var == "T" & nodloop$context == "C" & nodloop$is_silent == 1),])/3
# SYNONYMALL. G>C - 2169
nrow(nodloop[which(nodloop$ref== "G" & nodloop$var == "C" & nodloop$context == "A" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "G" & nodloop$var == "C" & nodloop$context == "T" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "G" & nodloop$var == "C" & nodloop$context == "G" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "G" & nodloop$var == "C" & nodloop$context == "C" & nodloop$is_silent == 1),])/3
# SYNONYMALL. C>A - 5182
nrow(nodloop[which(nodloop$ref== "C" & nodloop$var == "A" & nodloop$context == "A" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "C" & nodloop$var == "A" & nodloop$context == "T" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "C" & nodloop$var == "A" & nodloop$context == "G" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "C" & nodloop$var == "A" & nodloop$context == "C" & nodloop$is_silent == 1),])/3
# SYNONYMALL. C>T - 5182
nrow(nodloop[which(nodloop$ref== "C" & nodloop$var == "T" & nodloop$context == "A" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "C" & nodloop$var == "T" & nodloop$context == "T" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "C" & nodloop$var == "T" & nodloop$context == "G" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "C" & nodloop$var == "T" & nodloop$context == "C" & nodloop$is_silent == 1),])/3
# SYNONYMALL. C>G - 5182
nrow(nodloop[which(nodloop$ref== "C" & nodloop$var == "G" & nodloop$context == "A" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "C" & nodloop$var == "G" & nodloop$context == "T" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "C" & nodloop$var == "G" & nodloop$context == "G" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "C" & nodloop$var == "G" & nodloop$context == "C" & nodloop$is_silent == 1),])/3