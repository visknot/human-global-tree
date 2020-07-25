rm(list=ls(all=TRUE))
########

library(readr)
df <- read_delim("~/mitoclub/human-global-tree/Body/1Raw/CancerDataFromCampbell/mtDNA_snv_Oct2016.txt",
                                "\t", escape_double = FALSE, trim_ws = TRUE)


# ALL. Total number of mutations (rows in table)
nrow(df)

########
# ALL - All mutations of a certain type
# ALL. A>G
nrow(df[which(df$ref== "A" & df$var == "G"),])
# ALL. A>T
nrow(df[which(df$ref== "A" & df$var == "T"),])
# ALL. A>C
nrow(df[which(df$ref== "A" & df$var == "C"),])
# ALL. T>G
nrow(df[which(df$ref== "T" & df$var == "G"),])
# ALL. T>C
nrow(df[which(df$ref== "T" & df$var == "C"),])
# ALL. T>A
nrow(df[which(df$ref== "T" & df$var == "A"),])
# ALL. G>A
nrow(df[which(df$ref== "G" & df$var == "A"),])
# ALL. G>T
nrow(df[which(df$ref== "G" & df$var == "T"),])
# ALL. G>C
nrow(df[which(df$ref== "G" & df$var == "C"),])
# ALL. C>A
nrow(df[which(df$ref== "C" & df$var == "A"),])
# ALL. C>T
nrow(df[which(df$ref== "C" & df$var == "T"),])
# ALL. C>G
nrow(df[which(df$ref== "C" & df$var == "G"),])

########
# NODLOOP
#> summary(nodloop$position)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#577    3841    7983    8132   12354   16023
# ref var position  (16024..16569,1..576)
nodloop <- mtDNA[which(mtDNA$position > 576 & mtDNA$position < 16024),]
nrow(mtDNA[which(mtDNA$position > 576 & mtDNA$position < 16024),])

# NODLOOP. A>G
nrow(nodloop[which(nodloop$ref== "A" & nodloop$var == "G"),])
# NODLOOP. A>T
nrow(nodloop[which(nodloop$ref== "A" & nodloop$var == "T"),])
# NODLOOP. A>C
nrow(nodloop[which(nodloop$ref== "A" & nodloop$var == "C"),])
# NODLOOP. T>G
nrow(nodloop[which(nodloop$ref== "T" & nodloop$var == "G"),])
# NODLOOP. T>C
nrow(nodloop[which(nodloop$ref== "T" & nodloop$var == "C"),])
# NODLOOP. T>A
nrow(nodloop[which(nodloop$ref== "T" & nodloop$var == "A"),])
# NODLOOP. G>A
nrow(nodloop[which(nodloop$ref== "G" & nodloop$var == "A"),])
# NODLOOP. G>T
nrow(nodloop[which(nodloop$ref== "G" & nodloop$var == "T"),])
# NODLOOP. G>C
nrow(nodloop[which(nodloop$ref== "G" & nodloop$var == "C"),])
# NODLOOP. C>A
nrow(nodloop[which(nodloop$ref== "C" & nodloop$var == "A"),])
# NODLOOP. C>T
nrow(nodloop[which(nodloop$ref== "C" & nodloop$var == "T"),])
# NODLOOP. C>G
nrow(nodloop[which(nodloop$ref== "C" & nodloop$var == "G"),])


########
# SYNONYMALL
synonymall = mtDNA[which(mtDNA$is_silent == 1),]
nrow(nodloop[which(nodloop$is_silent == 1),])

# SYNONYMALL. A>G
nrow(synonymall[which(synonymall$ref== "A" & synonymall$var == "G"),])
# SYNONYMALL. A>T
nrow(synonymall[which(synonymall$ref== "A" & synonymall$var == "T"),])
# SYNONYMALL. A>C
nrow(synonymall[which(synonymall$ref== "A" & synonymall$var == "C"),])
# SYNONYMALL. T>G
nrow(synonymall[which(synonymall$ref== "T" & synonymall$var == "G"),])
# SYNONYMALL. T>C
nrow(synonymall[which(synonymall$ref== "T" & synonymall$var == "C"),])
# SYNONYMALL. T>A
nrow(synonymall[which(synonymall$ref== "T" & synonymall$var == "A"),])
# SYNONYMALL. G>A
nrow(synonymall[which(synonymall$ref== "G" & synonymall$var == "A"),])
# SYNONYMALL. G>T
nrow(synonymall[which(synonymall$ref== "G" & synonymall$var == "T"),])
# SYNONYMALL. G>C
nrow(synonymall[which(synonymall$ref== "G" & synonymall$var == "C"),])
# SYNONYMALL. C>A
nrow(synonymall[which(synonymall$ref== "C" & synonymall$var == "A"),])
# SYNONYMALL. C>T
nrow(synonymall[which(synonymall$ref== "C" & synonymall$var == "T"),])
# SYNONYMALL. C>G
nrow(synonymall[which(synonymall$ref== "C" & synonymall$var == "G"),])



########
# Adding preceeding_nucl from Expected
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
# ALL. A>G
nrow(df_merged[which(df_merged$ref== "A" & df_merged$var == "G" & df_merged$context == "A"),])/3
nrow(df_merged[which(df_merged$ref== "A" & df_merged$var == "G" & df_merged$context == "T"),])/3
nrow(df_merged[which(df_merged$ref== "A" & df_merged$var == "G" & df_merged$context == "G"),])/3
nrow(df_merged[which(df_merged$ref== "A" & df_merged$var == "G" & df_merged$context == "C"),])/3
# ALL. A>T
nrow(df_merged[which(df_merged$ref== "A" & df_merged$var == "T" & df_merged$context == "A"),])/3
nrow(df_merged[which(df_merged$ref== "A" & df_merged$var == "T" & df_merged$context == "T"),])/3
nrow(df_merged[which(df_merged$ref== "A" & df_merged$var == "T" & df_merged$context == "G"),])/3
nrow(df_merged[which(df_merged$ref== "A" & df_merged$var == "T" & df_merged$context == "C"),])/3
# ALL. A>C
nrow(df_merged[which(df_merged$ref== "A" & df_merged$var == "C" & df_merged$context == "A"),])/3
nrow(df_merged[which(df_merged$ref== "A" & df_merged$var == "C" & df_merged$context == "T"),])/3
nrow(df_merged[which(df_merged$ref== "A" & df_merged$var == "C" & df_merged$context == "G"),])/3
nrow(df_merged[which(df_merged$ref== "A" & df_merged$var == "C" & df_merged$context == "C"),])/3
# ALL. T>G
nrow(df_merged[which(df_merged$ref== "T" & df_merged$var == "G" & df_merged$context == "A"),])/3
nrow(df_merged[which(df_merged$ref== "T" & df_merged$var == "G" & df_merged$context == "T"),])/3
nrow(df_merged[which(df_merged$ref== "T" & df_merged$var == "G" & df_merged$context == "G"),])/3
nrow(df_merged[which(df_merged$ref== "T" & df_merged$var == "G" & df_merged$context == "C"),])/3
# ALL. T>C
nrow(df_merged[which(df_merged$ref== "T" & df_merged$var == "C" & df_merged$context == "A"),])/3
nrow(df_merged[which(df_merged$ref== "T" & df_merged$var == "C" & df_merged$context == "T"),])/3
nrow(df_merged[which(df_merged$ref== "T" & df_merged$var == "C" & df_merged$context == "G"),])/3
nrow(df_merged[which(df_merged$ref== "T" & df_merged$var == "C" & df_merged$context == "C"),])/3
# ALL. T>A
nrow(df_merged[which(df_merged$ref== "T" & df_merged$var == "A" & df_merged$context == "A"),])/3
nrow(df_merged[which(df_merged$ref== "T" & df_merged$var == "A" & df_merged$context == "T"),])/3
nrow(df_merged[which(df_merged$ref== "T" & df_merged$var == "A" & df_merged$context == "G"),])/3
nrow(df_merged[which(df_merged$ref== "T" & df_merged$var == "A" & df_merged$context == "C"),])/3
# ALL. G>A
nrow(df_merged[which(df_merged$ref== "G" & df_merged$var == "A" & df_merged$context == "A"),])/3
nrow(df_merged[which(df_merged$ref== "G" & df_merged$var == "A" & df_merged$context == "T"),])/3
nrow(df_merged[which(df_merged$ref== "G" & df_merged$var == "A" & df_merged$context == "G"),])/3
nrow(df_merged[which(df_merged$ref== "G" & df_merged$var == "A" & df_merged$context == "C"),])/3
# ALL. G>T
nrow(df_merged[which(df_merged$ref== "G" & df_merged$var == "T" & df_merged$context == "A"),])/3
nrow(df_merged[which(df_merged$ref== "G" & df_merged$var == "T" & df_merged$context == "T"),])/3
nrow(df_merged[which(df_merged$ref== "G" & df_merged$var == "T" & df_merged$context == "G"),])/3
nrow(df_merged[which(df_merged$ref== "G" & df_merged$var == "T" & df_merged$context == "C"),])/3
# ALL. G>C
nrow(df_merged[which(df_merged$ref== "G" & df_merged$var == "C" & df_merged$context == "A"),])/3
nrow(df_merged[which(df_merged$ref== "G" & df_merged$var == "C" & df_merged$context == "T"),])/3
nrow(df_merged[which(df_merged$ref== "G" & df_merged$var == "C" & df_merged$context == "G"),])/3
nrow(df_merged[which(df_merged$ref== "G" & df_merged$var == "C" & df_merged$context == "C"),])/3
# ALL. C>A
nrow(df_merged[which(df_merged$ref== "C" & df_merged$var == "A" & df_merged$context == "A"),])/3
nrow(df_merged[which(df_merged$ref== "C" & df_merged$var == "A" & df_merged$context == "T"),])/3
nrow(df_merged[which(df_merged$ref== "C" & df_merged$var == "A" & df_merged$context == "G"),])/3
nrow(df_merged[which(df_merged$ref== "C" & df_merged$var == "A" & df_merged$context == "C"),])/3
# ALL. C>T
nrow(df_merged[which(df_merged$ref== "C" & df_merged$var == "T" & df_merged$context == "A"),])/3
nrow(df_merged[which(df_merged$ref== "C" & df_merged$var == "T" & df_merged$context == "T"),])/3
nrow(df_merged[which(df_merged$ref== "C" & df_merged$var == "T" & df_merged$context == "G"),])/3
nrow(df_merged[which(df_merged$ref== "C" & df_merged$var == "T" & df_merged$context == "C"),])/3
# ALL. C>G
nrow(df_merged[which(df_merged$ref== "C" & df_merged$var == "G" & df_merged$context == "A"),])/3
nrow(df_merged[which(df_merged$ref== "C" & df_merged$var == "G" & df_merged$context == "T"),])/3
nrow(df_merged[which(df_merged$ref== "C" & df_merged$var == "G" & df_merged$context == "G"),])/3
nrow(df_merged[which(df_merged$ref== "C" & df_merged$var == "G" & df_merged$context == "C"),])/3


########
# NODLOOP - Everything outside of D-loop. D-loop coordinates from NC_012920 (16024..16569,1..576)
# NODLOOP. A>G
nrow(nodloop[which(nodloop$ref== "A" & nodloop$var == "G" & nodloop$context == "A"),])/3
nrow(nodloop[which(nodloop$ref== "A" & nodloop$var == "G" & nodloop$context == "T"),])/3
nrow(nodloop[which(nodloop$ref== "A" & nodloop$var == "G" & nodloop$context == "G"),])/3
nrow(nodloop[which(nodloop$ref== "A" & nodloop$var == "G" & nodloop$context == "C"),])/3
# NODLOOP. A>T
nrow(nodloop[which(nodloop$ref== "A" & nodloop$var == "T" & nodloop$context == "A"),])/3
nrow(nodloop[which(nodloop$ref== "A" & nodloop$var == "T" & nodloop$context == "T"),])/3
nrow(nodloop[which(nodloop$ref== "A" & nodloop$var == "T" & nodloop$context == "G"),])/3
nrow(nodloop[which(nodloop$ref== "A" & nodloop$var == "T" & nodloop$context == "C"),])/3
# NODLOOP. A>C
nrow(nodloop[which(nodloop$ref== "A" & nodloop$var == "C" & nodloop$context == "A"),])/3
nrow(nodloop[which(nodloop$ref== "A" & nodloop$var == "C" & nodloop$context == "T"),])/3
nrow(nodloop[which(nodloop$ref== "A" & nodloop$var == "C" & nodloop$context == "G"),])/3
nrow(nodloop[which(nodloop$ref== "A" & nodloop$var == "C" & nodloop$context == "C"),])/3
# NODLOOP. T>G
nrow(nodloop[which(nodloop$ref== "T" & nodloop$var == "G" & nodloop$context == "A"),])/3
nrow(nodloop[which(nodloop$ref== "T" & nodloop$var == "G" & nodloop$context == "T"),])/3
nrow(nodloop[which(nodloop$ref== "T" & nodloop$var == "G" & nodloop$context == "G"),])/3
nrow(nodloop[which(nodloop$ref== "T" & nodloop$var == "G" & nodloop$context == "C"),])/3
# NODLOOP. T>C
nrow(nodloop[which(nodloop$ref== "T" & nodloop$var == "C" & nodloop$context == "A"),])/3
nrow(nodloop[which(nodloop$ref== "T" & nodloop$var == "C" & nodloop$context == "T"),])/3
nrow(nodloop[which(nodloop$ref== "T" & nodloop$var == "C" & nodloop$context == "G"),])/3
nrow(nodloop[which(nodloop$ref== "T" & nodloop$var == "C" & nodloop$context == "C"),])/3
# NODLOOP. T>A
nrow(nodloop[which(nodloop$ref== "T" & nodloop$var == "A" & nodloop$context == "A"),])/3
nrow(nodloop[which(nodloop$ref== "T" & nodloop$var == "A" & nodloop$context == "T"),])/3
nrow(nodloop[which(nodloop$ref== "T" & nodloop$var == "A" & nodloop$context == "G"),])/3
nrow(nodloop[which(nodloop$ref== "T" & nodloop$var == "A" & nodloop$context == "C"),])/3
# NODLOOP. G>A
nrow(nodloop[which(nodloop$ref== "G" & nodloop$var == "A" & nodloop$context == "A"),])/3
nrow(nodloop[which(nodloop$ref== "G" & nodloop$var == "A" & nodloop$context == "T"),])/3
nrow(nodloop[which(nodloop$ref== "G" & nodloop$var == "A" & nodloop$context == "G"),])/3
nrow(nodloop[which(nodloop$ref== "G" & nodloop$var == "A" & nodloop$context == "C"),])/3
# NODLOOP. G>T
nrow(nodloop[which(nodloop$ref== "G" & nodloop$var == "T" & nodloop$context == "A"),])/3
nrow(nodloop[which(nodloop$ref== "G" & nodloop$var == "T" & nodloop$context == "T"),])/3
nrow(nodloop[which(nodloop$ref== "G" & nodloop$var == "T" & nodloop$context == "G"),])/3
nrow(nodloop[which(nodloop$ref== "G" & nodloop$var == "T" & nodloop$context == "C"),])/3
# NODLOOP. G>C
nrow(nodloop[which(nodloop$ref== "G" & nodloop$var == "C" & nodloop$context == "A"),])/3
nrow(nodloop[which(nodloop$ref== "G" & nodloop$var == "C" & nodloop$context == "T"),])/3
nrow(nodloop[which(nodloop$ref== "G" & nodloop$var == "C" & nodloop$context == "G"),])/3
nrow(nodloop[which(nodloop$ref== "G" & nodloop$var == "C" & nodloop$context == "C"),])/3
# NODLOOP. C>A
nrow(nodloop[which(nodloop$ref== "C" & nodloop$var == "A" & nodloop$context == "A"),])/3
nrow(nodloop[which(nodloop$ref== "C" & nodloop$var == "A" & nodloop$context == "T"),])/3
nrow(nodloop[which(nodloop$ref== "C" & nodloop$var == "A" & nodloop$context == "G"),])/3
nrow(nodloop[which(nodloop$ref== "C" & nodloop$var == "A" & nodloop$context == "C"),])/3
# NODLOOP. C>T
nrow(nodloop[which(nodloop$ref== "C" & nodloop$var == "T" & nodloop$context == "A"),])/3
nrow(nodloop[which(nodloop$ref== "C" & nodloop$var == "T" & nodloop$context == "T"),])/3
nrow(nodloop[which(nodloop$ref== "C" & nodloop$var == "T" & nodloop$context == "G"),])/3
nrow(nodloop[which(nodloop$ref== "C" & nodloop$var == "T" & nodloop$context == "C"),])/3
# NODLOOP. C>G
nrow(nodloop[which(nodloop$ref== "C" & nodloop$var == "G" & nodloop$context == "A"),])/3
nrow(nodloop[which(nodloop$ref== "C" & nodloop$var == "G" & nodloop$context == "T"),])/3
nrow(nodloop[which(nodloop$ref== "C" & nodloop$var == "G" & nodloop$context == "G"),])/3
nrow(nodloop[which(nodloop$ref== "C" & nodloop$var == "G" & nodloop$context == "C"),])/3


########
# SYNONYMALL - All synonymous variants
# SYNONYMALL. A>G
nrow(nodloop[which(nodloop$ref== "A" & nodloop$var == "G" & nodloop$context == "A" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "A" & nodloop$var == "G" & nodloop$context == "T" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "A" & nodloop$var == "G" & nodloop$context == "G" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "A" & nodloop$var == "G" & nodloop$context == "C" & nodloop$is_silent == 1),])/3
# SYNONYMALL. A>T
nrow(nodloop[which(nodloop$ref== "A" & nodloop$var == "T" & nodloop$context == "A" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "A" & nodloop$var == "T" & nodloop$context == "T" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "A" & nodloop$var == "T" & nodloop$context == "G" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "A" & nodloop$var == "T" & nodloop$context == "C" & nodloop$is_silent == 1),])/3
# SYNONYMALL. A>C
nrow(nodloop[which(nodloop$ref== "A" & nodloop$var == "C" & nodloop$context == "A" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "A" & nodloop$var == "C" & nodloop$context == "T" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "A" & nodloop$var == "C" & nodloop$context == "G" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "A" & nodloop$var == "C" & nodloop$context == "C" & nodloop$is_silent == 1),])/3
# SYNONYMALL. T>G
nrow(nodloop[which(nodloop$ref== "T" & nodloop$var == "G" & nodloop$context == "A" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "T" & nodloop$var == "G" & nodloop$context == "T" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "T" & nodloop$var == "G" & nodloop$context == "G" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "T" & nodloop$var == "G" & nodloop$context == "C" & nodloop$is_silent == 1),])/3
# SYNONYMALL. T>C
nrow(nodloop[which(nodloop$ref== "T" & nodloop$var == "C" & nodloop$context == "A" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "T" & nodloop$var == "C" & nodloop$context == "T" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "T" & nodloop$var == "C" & nodloop$context == "G" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "T" & nodloop$var == "C" & nodloop$context == "C" & nodloop$is_silent == 1),])/3
# SYNONYMALL. T>A
nrow(nodloop[which(nodloop$ref== "T" & nodloop$var == "A" & nodloop$context == "A" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "T" & nodloop$var == "A" & nodloop$context == "T" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "T" & nodloop$var == "A" & nodloop$context == "G" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "T" & nodloop$var == "A" & nodloop$context == "C" & nodloop$is_silent == 1),])/3
# SYNONYMALL. G>A
nrow(nodloop[which(nodloop$ref== "G" & nodloop$var == "A" & nodloop$context == "A" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "G" & nodloop$var == "A" & nodloop$context == "T" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "G" & nodloop$var == "A" & nodloop$context == "G" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "G" & nodloop$var == "A" & nodloop$context == "C" & nodloop$is_silent == 1),])/3
# SYNONYMALL. G>T
nrow(nodloop[which(nodloop$ref== "G" & nodloop$var == "T" & nodloop$context == "A" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "G" & nodloop$var == "T" & nodloop$context == "T" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "G" & nodloop$var == "T" & nodloop$context == "G" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "G" & nodloop$var == "T" & nodloop$context == "C" & nodloop$is_silent == 1),])/3
# SYNONYMALL. G>C
nrow(nodloop[which(nodloop$ref== "G" & nodloop$var == "C" & nodloop$context == "A" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "G" & nodloop$var == "C" & nodloop$context == "T" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "G" & nodloop$var == "C" & nodloop$context == "G" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "G" & nodloop$var == "C" & nodloop$context == "C" & nodloop$is_silent == 1),])/3
# SYNONYMALL. C>A
nrow(nodloop[which(nodloop$ref== "C" & nodloop$var == "A" & nodloop$context == "A" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "C" & nodloop$var == "A" & nodloop$context == "T" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "C" & nodloop$var == "A" & nodloop$context == "G" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "C" & nodloop$var == "A" & nodloop$context == "C" & nodloop$is_silent == 1),])/3
# SYNONYMALL. C>T
nrow(nodloop[which(nodloop$ref== "C" & nodloop$var == "T" & nodloop$context == "A" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "C" & nodloop$var == "T" & nodloop$context == "T" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "C" & nodloop$var == "T" & nodloop$context == "G" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "C" & nodloop$var == "T" & nodloop$context == "C" & nodloop$is_silent == 1),])/3
# SYNONYMALL. C>G
nrow(nodloop[which(nodloop$ref== "C" & nodloop$var == "G" & nodloop$context == "A" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "C" & nodloop$var == "G" & nodloop$context == "T" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "C" & nodloop$var == "G" & nodloop$context == "G" & nodloop$is_silent == 1),])/3
nrow(nodloop[which(nodloop$ref== "C" & nodloop$var == "G" & nodloop$context == "C" & nodloop$is_silent == 1),])/3
