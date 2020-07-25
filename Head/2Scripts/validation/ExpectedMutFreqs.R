rm(list=ls(all=TRUE))
########

# Input data
df = read.table("/Users/xtinaushakova/mitoclub/human-global-tree/Body/1Raw/FromCovidProject/HumanMtDnaRefSeq.fasta.ExtensiveMut.vcf.ann", 
                sep='\n', 
                header = FALSE)

# Total number of mutations (rows in table) - 49707
nrow(df)

########
# ALL - All mutations of a certain type
# ALL. A>G - 5124
length(grep('A\tG', df$V1))
# ALL. A>T - 5124
length(grep('A\tT', df$V1))
# ALL. A>C - 5124
length(grep('A\tC', df$V1))
# ALL. T>G - 4094
length(grep('T\tG', df$V1))
# ALL. T>C - 4094
length(grep('T\tC', df$V1))
# ALL. T>A - 4094
length(grep('T\tA', df$V1))
# ALL. G>A - 2169
length(grep('G\tA', df$V1))
# ALL. G>T - 2169
length(grep('G\tT', df$V1))
# ALL. G>C - 2169
length(grep('G\tC', df$V1))
# ALL. C>A - 5182
length(grep('C\tA', df$V1))
# ALL. C>T - 5182
length(grep('C\tT', df$V1))
# ALL. C>G - 5182
length(grep('C\tG', df$V1))


########
# NODLOOP - Everything outside of D-loop. D-loop coordinates from NC_012920 (16024..16569,1..576)
# Rows from df[1729,] to df[48072] 
length(df$V1[1729:48072])
# NODLOOP. A>G - 4785
length(grep('A\tG', df$V1[1729:48072]))
# NODLOOP. A>T - 4785
length(grep('A\tT', df$V1[1729:48072]))
# NODLOOP. A>C - 4785
length(grep('A\tC', df$V1[1729:48072]))
# NODLOOP. T>G - 3835
length(grep('T\tG', df$V1[1729:48072]))
# NODLOOP. T>C - 3835
length(grep('T\tC', df$V1[1729:48072]))
# NODLOOP. T>A - 3835
length(grep('T\tA', df$V1[1729:48072]))
# NODLOOP. G>A - 2017
length(grep('G\tA', df$V1[1729:48072]))
# NODLOOP. G>T - 2017
length(grep('G\tT', df$V1[1729:48072]))
# NODLOOP. G>C - 2017
length(grep('G\tC', df$V1[1729:48072]))
# NODLOOP. C>A - 4811
length(grep('C\tA', df$V1[1729:48072]))
# NODLOOP. C>T - 4811
length(grep('C\tT', df$V1[1729:48072]))
# NODLOOP. C>G - 4811
length(grep('C\tG', df$V1[1729:48072]))


########
# SYNONYMALL - All synonymous variants - 8678
length(grep('synonymous_variant', df$V1))
# SYNONYMALL. A>G - 1154
length(grep('A\tG.*synonymous_variant', df$V1))
# SYNONYMALL. A>T - 979
length(grep('A\tT.*synonymous_variant', df$V1))
# SYNONYMALL. A>C - 968
length(grep('A\tC.*synonymous_variant', df$V1))
# SYNONYMALL. T>G - 293
length(grep('T\tG.*synonymous_variant', df$V1))
# SYNONYMALL. T>C - 665
length(grep('T\tC.*synonymous_variant', df$V1))
# SYNONYMALL. T>A - 404
length(grep('T\tA.*synonymous_variant', df$V1))
# SYNONYMALL. G>A - 149
length(grep('G\tA.*synonymous_variant', df$V1))
# SYNONYMALL. G>T - 103
length(grep('G\tT.*synonymous_variant', df$V1))
# SYNONYMALL. G>C - 101
length(grep('G\tC.*synonymous_variant', df$V1))
# SYNONYMALL. C>A - 1079
length(grep('C\tA.*synonymous_variant', df$V1))
# SYNONYMALL. C>T - 1928
length(grep('C\tT.*synonymous_variant', df$V1))
# SYNONYMALL. C>G - 855
length(grep('C\tG.*synonymous_variant', df$V1))


########
# SYN4FOLD - All 4-fold degenerate, 3 synonymous in one position in a row
count <- 0
datalist = list()
for (i in seq(1, nrow(df), by = 3)) {
  end = i + 2
  temp = as.data.frame(df[i:end,])
  res = grep('synonymous_variant', temp$`df[i:end, ]`)
  if (length(res) == 3) {
    print(length(res))
    count = count + 1
    datalist[[count]] <- temp
  }
}
df4fold = do.call(rbind, datalist)
########
# SYN4FOLD. Total 6144
length(grep('synonymous_variant', big_data$`df[i:end, ]`))
# SYNONYMALL. A>G - 802
length(grep('A\tG.*synonymous_variant', big_data$`df[i:end, ]`))
# SYNONYMALL. A>T - 802
length(grep('A\tT.*synonymous_variant', big_data$`df[i:end, ]`))
# SYNONYMALL. A>C - 802
length(grep('A\tC.*synonymous_variant', big_data$`df[i:end, ]`))
# SYNONYMALL. T>G - 290
length(grep('T\tG.*synonymous_variant', big_data$`df[i:end, ]`))
# SYNONYMALL. T>C - 290
length(grep('T\tC.*synonymous_variant', big_data$`df[i:end, ]`))
# SYNONYMALL. T>A - 290
length(grep('T\tA.*synonymous_variant', big_data$`df[i:end, ]`))
# SYNONYMALL. G>A - 101
length(grep('G\tA.*synonymous_variant', big_data$`df[i:end, ]`))
# SYNONYMALL. G>T - 101
length(grep('G\tT.*synonymous_variant', big_data$`df[i:end, ]`))
# SYNONYMALL. G>C - 101
length(grep('G\tC.*synonymous_variant', big_data$`df[i:end, ]`))
# SYNONYMALL. C>A - 855
length(grep('C\tA.*synonymous_variant', big_data$`df[i:end, ]`))
# SYNONYMALL. C>T - 855
length(grep('C\tT.*synonymous_variant', big_data$`df[i:end, ]`))
# SYNONYMALL. C>G - 855
length(grep('C\tG.*synonymous_variant', big_data$`df[i:end, ]`))

########
# Adding preceeding_nucl
datalist = c()
for (i in seq(4, nrow(df), by = 3)) {
  preceeding = i-1
  end = i + 2
  temp = unlist(strsplit(substr(as.character(df[preceeding,]), 1, 25), "\t"))
  preceeding_nucleotide = temp[4]
  print(temp)
  res = c(rep(preceeding_nucleotide, 3))
  datalist = c(datalist, res)
  }
}
# Append to the top
last_nucl = c(rep("G", 3))
datalist = c(last_nucl, datalist)

df_context = df
df_context['context'] = datalist


########
# ALL. A>G - 5124
length(grep('A\tG', df_context$V1[which(df_context$context == "A")]))
length(grep('A\tG', df_context$V1[which(df_context$context == "T")]))
length(grep('A\tG', df_context$V1[which(df_context$context == "G")]))
length(grep('A\tG', df_context$V1[which(df_context$context == "C")]))
# ALL. A>T - 5124
length(grep('A\tT', df_context$V1[which(df_context$context == "A")]))
length(grep('A\tT', df_context$V1[which(df_context$context == "T")]))
length(grep('A\tT', df_context$V1[which(df_context$context == "G")]))
length(grep('A\tT', df_context$V1[which(df_context$context == "C")]))
# ALL. A>C - 5124
length(grep('A\tC', df$V1))
length(grep('A\tC', df_context$V1[which(df_context$context == "A")]))
length(grep('A\tC', df_context$V1[which(df_context$context == "T")]))
length(grep('A\tC', df_context$V1[which(df_context$context == "G")]))
length(grep('A\tC', df_context$V1[which(df_context$context == "C")]))
# ALL. T>G - 4094
length(grep('T\tG', df$V1))
length(grep('T\tG', df_context$V1[which(df_context$context == "A")]))
length(grep('T\tG', df_context$V1[which(df_context$context == "T")]))
length(grep('T\tG', df_context$V1[which(df_context$context == "G")]))
length(grep('T\tG', df_context$V1[which(df_context$context == "C")]))
# ALL. T>C - 4094
length(grep('T\tC', df$V1))
length(grep('T\tC', df_context$V1[which(df_context$context == "A")]))
length(grep('T\tC', df_context$V1[which(df_context$context == "T")]))
length(grep('T\tC', df_context$V1[which(df_context$context == "G")]))
length(grep('T\tC', df_context$V1[which(df_context$context == "C")]))
# ALL. T>A - 4094
length(grep('T\tA', df$V1))
length(grep('T\tA', df_context$V1[which(df_context$context == "A")]))
length(grep('T\tA', df_context$V1[which(df_context$context == "T")]))
length(grep('T\tA', df_context$V1[which(df_context$context == "G")]))
length(grep('T\tA', df_context$V1[which(df_context$context == "C")]))
# ALL. G>A - 2169
length(grep('G\tA', df$V1))
length(grep('G\tA', df_context$V1[which(df_context$context == "A")]))
length(grep('G\tA', df_context$V1[which(df_context$context == "T")]))
length(grep('G\tA', df_context$V1[which(df_context$context == "G")]))
length(grep('G\tA', df_context$V1[which(df_context$context == "C")]))
# ALL. G>T - 2169
length(grep('G\tT', df$V1))
length(grep('G\tT', df_context$V1[which(df_context$context == "A")]))
length(grep('G\tT', df_context$V1[which(df_context$context == "T")]))
length(grep('G\tT', df_context$V1[which(df_context$context == "G")]))
length(grep('G\tT', df_context$V1[which(df_context$context == "C")]))
# ALL. G>C - 2169
length(grep('G\tC', df$V1))
length(grep('G\tC', df_context$V1[which(df_context$context == "A")]))
length(grep('G\tC', df_context$V1[which(df_context$context == "T")]))
length(grep('G\tC', df_context$V1[which(df_context$context == "G")]))
length(grep('G\tC', df_context$V1[which(df_context$context == "C")]))
# ALL. C>A - 5182
length(grep('C\tA', df$V1))
length(grep('C\tA', df_context$V1[which(df_context$context == "A")]))
length(grep('C\tA', df_context$V1[which(df_context$context == "T")]))
length(grep('C\tA', df_context$V1[which(df_context$context == "G")]))
length(grep('C\tA', df_context$V1[which(df_context$context == "C")]))
# ALL. C>T - 5182
length(grep('C\tT', df$V1))
length(grep('C\tT', df_context$V1[which(df_context$context == "A")]))
length(grep('C\tT', df_context$V1[which(df_context$context == "T")]))
length(grep('C\tT', df_context$V1[which(df_context$context == "G")]))
length(grep('C\tT', df_context$V1[which(df_context$context == "C")]))
# ALL. C>G - 5182
length(grep('C\tG', df$V1))
length(grep('C\tG', df_context$V1[which(df_context$context == "A")]))
length(grep('C\tG', df_context$V1[which(df_context$context == "T")]))
length(grep('C\tG', df_context$V1[which(df_context$context == "G")]))
length(grep('C\tG', df_context$V1[which(df_context$context == "C")]))


########
# NODLOOP - Everything outside of D-loop. D-loop coordinates from NC_012920 (16024..16569,1..576)
# Rows from df[1729,] to df[48072] 
# NODLOOP. A>G - 4785
length(grep('A\tG', df_context$V1[1729:48072][which(df_context$context == "A")]))
length(grep('A\tG', df_context$V1[1729:48072][which(df_context$context == "T")]))
length(grep('A\tG', df_context$V1[1729:48072][which(df_context$context == "G")]))
length(grep('A\tG', df_context$V1[1729:48072][which(df_context$context == "C")]))
# ALL. A>T - 5124
length(grep('A\tT', df_context$V1[1729:48072][which(df_context$context == "A")]))
length(grep('A\tT', df_context$V1[1729:48072][which(df_context$context == "T")]))
length(grep('A\tT', df_context$V1[1729:48072][which(df_context$context == "G")]))
length(grep('A\tT', df_context$V1[1729:48072][which(df_context$context == "C")]))
# ALL. A>C - 5124
length(grep('A\tC', df_context$V1[1729:48072]))
length(grep('A\tC', df_context$V1[1729:48072][which(df_context$context == "A")]))
length(grep('A\tC', df_context$V1[1729:48072][which(df_context$context == "T")]))
length(grep('A\tC', df_context$V1[1729:48072][which(df_context$context == "G")]))
length(grep('A\tC', df_context$V1[1729:48072][which(df_context$context == "C")]))
# ALL. T>G - 4094
length(grep('T\tG', df_context$V1[1729:48072]))
length(grep('T\tG', df_context$V1[1729:48072][which(df_context$context == "A")]))
length(grep('T\tG', df_context$V1[1729:48072][which(df_context$context == "T")]))
length(grep('T\tG', df_context$V1[1729:48072][which(df_context$context == "G")]))
length(grep('T\tG', df_context$V1[1729:48072][which(df_context$context == "C")]))
# ALL. T>C - 4094
length(grep('T\tC', df_context$V1[1729:48072]))
length(grep('T\tC', df_context$V1[1729:48072][which(df_context$context == "A")]))
length(grep('T\tC', df_context$V1[1729:48072][which(df_context$context == "T")]))
length(grep('T\tC', df_context$V1[1729:48072][which(df_context$context == "G")]))
length(grep('T\tC', df_context$V1[1729:48072][which(df_context$context == "C")]))
# ALL. T>A - 4094
length(grep('T\tA', df_context$V1[1729:48072]))
length(grep('T\tA', df_context$V1[1729:48072][which(df_context$context == "A")]))
length(grep('T\tA', df_context$V1[1729:48072][which(df_context$context == "T")]))
length(grep('T\tA', df_context$V1[1729:48072][which(df_context$context == "G")]))
length(grep('T\tA', df_context$V1[1729:48072][which(df_context$context == "C")]))
# ALL. G>A - 2169
length(grep('G\tA', df_context$V1[1729:48072]))
length(grep('G\tA', df_context$V1[1729:48072][which(df_context$context == "A")]))
length(grep('G\tA', df_context$V1[1729:48072][which(df_context$context == "T")]))
length(grep('G\tA', df_context$V1[1729:48072][which(df_context$context == "G")]))
length(grep('G\tA', df_context$V1[1729:48072][which(df_context$context == "C")]))
# ALL. G>T - 2169
length(grep('G\tT', df_context$V1[1729:48072]))
length(grep('G\tT', df_context$V1[1729:48072][which(df_context$context == "A")]))
length(grep('G\tT', df_context$V1[1729:48072][which(df_context$context == "T")]))
length(grep('G\tT', df_context$V1[1729:48072][which(df_context$context == "G")]))
length(grep('G\tT', df_context$V1[1729:48072][which(df_context$context == "C")]))
# ALL. G>C - 2169
length(grep('G\tC', df_context$V1[1729:48072]))
length(grep('G\tC', df_context$V1[1729:48072][which(df_context$context == "A")]))
length(grep('G\tC', df_context$V1[1729:48072][which(df_context$context == "T")]))
length(grep('G\tC', df_context$V1[1729:48072][which(df_context$context == "G")]))
length(grep('G\tC', df_context$V1[1729:48072][which(df_context$context == "C")]))
# ALL. C>A - 5182
length(grep('C\tA', df_context$V1[1729:48072]))
length(grep('C\tA', df_context$V1[1729:48072][which(df_context$context == "A")]))
length(grep('C\tA', df_context$V1[1729:48072][which(df_context$context == "T")]))
length(grep('C\tA', df_context$V1[1729:48072][which(df_context$context == "G")]))
length(grep('C\tA', df_context$V1[1729:48072][which(df_context$context == "C")]))
# ALL. C>T - 5182
length(grep('C\tT', df_context$V1[1729:48072]))
length(grep('C\tT', df_context$V1[1729:48072][which(df_context$context == "A")]))
length(grep('C\tT', df_context$V1[1729:48072][which(df_context$context == "T")]))
length(grep('C\tT', df_context$V1[1729:48072][which(df_context$context == "G")]))
length(grep('C\tT', df_context$V1[1729:48072][which(df_context$context == "C")]))
# ALL. C>G - 5182
length(grep('C\tG', df_context$V1[1729:48072]))
length(grep('C\tG', df_context$V1[1729:48072][which(df_context$context == "A")]))
length(grep('C\tG', df_context$V1[1729:48072][which(df_context$context == "T")]))
length(grep('C\tG', df_context$V1[1729:48072][which(df_context$context == "G")]))
length(grep('C\tG', df_context$V1[1729:48072][which(df_context$context == "C")]))


########
# SYNONYMALL - All synonymous variants - 8678
length(grep('synonymous_variant', df$V1))
# SYNONYMALL. A>G - 1154
length(grep('A\tG.*synonymous_variant', df_context$V1[which(df_context$context == "A")]))
length(grep('A\tG.*synonymous_variant', df_context$V1[which(df_context$context == "T")]))
length(grep('A\tG.*synonymous_variant', df_context$V1[which(df_context$context == "G")]))
length(grep('A\tG.*synonymous_variant', df_context$V1[which(df_context$context == "C")]))
# ALL. A>T - 5124
length(grep('A\tT.*synonymous_variant', df_context$V1[which(df_context$context == "A")]))
length(grep('A\tT.*synonymous_variant', df_context$V1[which(df_context$context == "T")]))
length(grep('A\tT.*synonymous_variant', df_context$V1[which(df_context$context == "G")]))
length(grep('A\tT.*synonymous_variant', df_context$V1[which(df_context$context == "C")]))
# ALL. A>C - 5124
length(grep('A\tC.*synonymous_variant', df$V1))
length(grep('A\tC.*synonymous_variant', df_context$V1[which(df_context$context == "A")]))
length(grep('A\tC.*synonymous_variant', df_context$V1[which(df_context$context == "T")]))
length(grep('A\tC.*synonymous_variant', df_context$V1[which(df_context$context == "G")]))
length(grep('A\tC.*synonymous_variant', df_context$V1[which(df_context$context == "C")]))
# ALL. T>G - 4094
length(grep('T\tG.*synonymous_variant', df$V1))
length(grep('T\tG.*synonymous_variant', df_context$V1[which(df_context$context == "A")]))
length(grep('T\tG.*synonymous_variant', df_context$V1[which(df_context$context == "T")]))
length(grep('T\tG.*synonymous_variant', df_context$V1[which(df_context$context == "G")]))
length(grep('T\tG.*synonymous_variant', df_context$V1[which(df_context$context == "C")]))
# ALL. T>C - 4094
length(grep('T\tC.*synonymous_variant', df$V1))
length(grep('T\tC.*synonymous_variant', df_context$V1[which(df_context$context == "A")]))
length(grep('T\tC.*synonymous_variant', df_context$V1[which(df_context$context == "T")]))
length(grep('T\tC.*synonymous_variant', df_context$V1[which(df_context$context == "G")]))
length(grep('T\tC.*synonymous_variant', df_context$V1[which(df_context$context == "C")]))
# ALL. T>A - 4094
length(grep('T\tA.*synonymous_variant', df$V1))
length(grep('T\tA.*synonymous_variant', df_context$V1[which(df_context$context == "A")]))
length(grep('T\tA.*synonymous_variant', df_context$V1[which(df_context$context == "T")]))
length(grep('T\tA.*synonymous_variant', df_context$V1[which(df_context$context == "G")]))
length(grep('T\tA.*synonymous_variant', df_context$V1[which(df_context$context == "C")]))
# ALL. G>A - 2169
length(grep('G\tA.*synonymous_variant', df$V1))
length(grep('G\tA.*synonymous_variant', df_context$V1[which(df_context$context == "A")]))
length(grep('G\tA.*synonymous_variant', df_context$V1[which(df_context$context == "T")]))
length(grep('G\tA.*synonymous_variant', df_context$V1[which(df_context$context == "G")]))
length(grep('G\tA.*synonymous_variant', df_context$V1[which(df_context$context == "C")]))
# ALL. G>T - 2169
length(grep('G\tT.*synonymous_variant', df$V1))
length(grep('G\tT.*synonymous_variant', df_context$V1[which(df_context$context == "A")]))
length(grep('G\tT.*synonymous_variant', df_context$V1[which(df_context$context == "T")]))
length(grep('G\tT.*synonymous_variant', df_context$V1[which(df_context$context == "G")]))
length(grep('G\tT.*synonymous_variant', df_context$V1[which(df_context$context == "C")]))
# ALL. G>C - 2169
length(grep('G\tC.*synonymous_variant', df$V1))
length(grep('G\tC.*synonymous_variant', df_context$V1[which(df_context$context == "A")]))
length(grep('G\tC.*synonymous_variant', df_context$V1[which(df_context$context == "T")]))
length(grep('G\tC.*synonymous_variant', df_context$V1[which(df_context$context == "G")]))
length(grep('G\tC.*synonymous_variant', df_context$V1[which(df_context$context == "C")]))
# ALL. C>A - 5182
length(grep('C\tA.*synonymous_variant', df$V1))
length(grep('C\tA.*synonymous_variant', df_context$V1[which(df_context$context == "A")]))
length(grep('C\tA.*synonymous_variant', df_context$V1[which(df_context$context == "T")]))
length(grep('C\tA.*synonymous_variant', df_context$V1[which(df_context$context == "G")]))
length(grep('C\tA.*synonymous_variant', df_context$V1[which(df_context$context == "C")]))
# ALL. C>T - 5182
length(grep('C\tT.*synonymous_variant', df$V1))
length(grep('C\tT.*synonymous_variant', df_context$V1[which(df_context$context == "A")]))
length(grep('C\tT.*synonymous_variant', df_context$V1[which(df_context$context == "T")]))
length(grep('C\tT.*synonymous_variant', df_context$V1[which(df_context$context == "G")]))
length(grep('C\tT.*synonymous_variant', df_context$V1[which(df_context$context == "C")]))
# ALL. C>G - 5182
length(grep('C\tG.*synonymous_variant', df$V1))
length(grep('C\tG.*synonymous_variant', df_context$V1[which(df_context$context == "A")]))
length(grep('C\tG.*synonymous_variant', df_context$V1[which(df_context$context == "T")]))
length(grep('C\tG.*synonymous_variant', df_context$V1[which(df_context$context == "G")]))
length(grep('C\tG.*synonymous_variant', df_context$V1[which(df_context$context == "C")]))




########
# SYN4FOLD - All 4-fold degenerate, 3 synonymous in one position in a row
count <- 0
datalist = list()
for (i in seq(1, nrow(df_context), by = 3)) {
  end = i + 2
  temp = as.data.frame(df_context[i:end,])
  res = grep('synonymous_variant', temp$V1)
  if (length(res) == 3) {
    count = count + 1
    datalist[[count]] <- temp
  }
}
df_context4fold = do.call(rbind, datalist)
########
# SYN4FOLD. Total 6144
nrow(df_context4fold)
# SYNONYMALL. A>G - 1154
# SYN4FOLD. A>G - 1154
length(grep('A\tG.*synonymous_variant', df_context4fold$V1[which(df_context4fold$context == "A")]))
length(grep('A\tG.*synonymous_variant', df_context4fold$V1[which(df_context4fold$context == "T")]))
length(grep('A\tG.*synonymous_variant', df_context4fold$V1[which(df_context4fold$context == "G")]))
length(grep('A\tG.*synonymous_variant', df_context4fold$V1[which(df_context4fold$context == "C")]))
# SYN4FOLD. A>T - 5124
length(grep('A\tT.*synonymous_variant', df_context4fold$V1[which(df_context4fold$context == "A")]))
length(grep('A\tT.*synonymous_variant', df_context4fold$V1[which(df_context4fold$context == "T")]))
length(grep('A\tT.*synonymous_variant', df_context4fold$V1[which(df_context4fold$context == "G")]))
length(grep('A\tT.*synonymous_variant', df_context4fold$V1[which(df_context4fold$context == "C")]))
# SYN4FOLD. A>C - 5124
length(grep('A\tC.*synonymous_variant', df$V1))
length(grep('A\tC.*synonymous_variant', df_context4fold$V1[which(df_context4fold$context == "A")]))
length(grep('A\tC.*synonymous_variant', df_context4fold$V1[which(df_context4fold$context == "T")]))
length(grep('A\tC.*synonymous_variant', df_context4fold$V1[which(df_context4fold$context == "G")]))
length(grep('A\tC.*synonymous_variant', df_context4fold$V1[which(df_context4fold$context == "C")]))
# SYN4FOLD. T>G - 4094
length(grep('T\tG.*synonymous_variant', df$V1))
length(grep('T\tG.*synonymous_variant', df_context4fold$V1[which(df_context4fold$context == "A")]))
length(grep('T\tG.*synonymous_variant', df_context4fold$V1[which(df_context4fold$context == "T")]))
length(grep('T\tG.*synonymous_variant', df_context4fold$V1[which(df_context4fold$context == "G")]))
length(grep('T\tG.*synonymous_variant', df_context4fold$V1[which(df_context4fold$context == "C")]))
# SYN4FOLD. T>C - 4094
length(grep('T\tC.*synonymous_variant', df$V1))
length(grep('T\tC.*synonymous_variant', df_context4fold$V1[which(df_context4fold$context == "A")]))
length(grep('T\tC.*synonymous_variant', df_context4fold$V1[which(df_context4fold$context == "T")]))
length(grep('T\tC.*synonymous_variant', df_context4fold$V1[which(df_context4fold$context == "G")]))
length(grep('T\tC.*synonymous_variant', df_context4fold$V1[which(df_context4fold$context == "C")]))
# SYN4FOLD. T>A - 4094
length(grep('T\tA.*synonymous_variant', df$V1))
length(grep('T\tA.*synonymous_variant', df_context4fold$V1[which(df_context4fold$context == "A")]))
length(grep('T\tA.*synonymous_variant', df_context4fold$V1[which(df_context4fold$context == "T")]))
length(grep('T\tA.*synonymous_variant', df_context4fold$V1[which(df_context4fold$context == "G")]))
length(grep('T\tA.*synonymous_variant', df_context4fold$V1[which(df_context4fold$context == "C")]))
# SYN4FOLD. G>A - 2169
length(grep('G\tA.*synonymous_variant', df$V1))
length(grep('G\tA.*synonymous_variant', df_context4fold$V1[which(df_context4fold$context == "A")]))
length(grep('G\tA.*synonymous_variant', df_context4fold$V1[which(df_context4fold$context == "T")]))
length(grep('G\tA.*synonymous_variant', df_context4fold$V1[which(df_context4fold$context == "G")]))
length(grep('G\tA.*synonymous_variant', df_context4fold$V1[which(df_context4fold$context == "C")]))
# SYN4FOLD. G>T - 2169
length(grep('G\tT.*synonymous_variant', df$V1))
length(grep('G\tT.*synonymous_variant', df_context4fold$V1[which(df_context4fold$context == "A")]))
length(grep('G\tT.*synonymous_variant', df_context4fold$V1[which(df_context4fold$context == "T")]))
length(grep('G\tT.*synonymous_variant', df_context4fold$V1[which(df_context4fold$context == "G")]))
length(grep('G\tT.*synonymous_variant', df_context4fold$V1[which(df_context4fold$context == "C")]))
# SYN4FOLD. G>C - 2169
length(grep('G\tC.*synonymous_variant', df$V1))
length(grep('G\tC.*synonymous_variant', df_context4fold$V1[which(df_context4fold$context == "A")]))
length(grep('G\tC.*synonymous_variant', df_context4fold$V1[which(df_context4fold$context == "T")]))
length(grep('G\tC.*synonymous_variant', df_context4fold$V1[which(df_context4fold$context == "G")]))
length(grep('G\tC.*synonymous_variant', df_context4fold$V1[which(df_context4fold$context == "C")]))
# SYN4FOLD. C>A - 5182
length(grep('C\tA.*synonymous_variant', df$V1))
length(grep('C\tA.*synonymous_variant', df_context4fold$V1[which(df_context4fold$context == "A")]))
length(grep('C\tA.*synonymous_variant', df_context4fold$V1[which(df_context4fold$context == "T")]))
length(grep('C\tA.*synonymous_variant', df_context4fold$V1[which(df_context4fold$context == "G")]))
length(grep('C\tA.*synonymous_variant', df_context4fold$V1[which(df_context4fold$context == "C")]))
# SYN4FOLD. C>T - 5182
length(grep('C\tT.*synonymous_variant', df$V1))
length(grep('C\tT.*synonymous_variant', df_context4fold$V1[which(df_context4fold$context == "A")]))
length(grep('C\tT.*synonymous_variant', df_context4fold$V1[which(df_context4fold$context == "T")]))
length(grep('C\tT.*synonymous_variant', df_context4fold$V1[which(df_context4fold$context == "G")]))
length(grep('C\tT.*synonymous_variant', df_context4fold$V1[which(df_context4fold$context == "C")]))
# SYN4FOLD. C>G - 5182
length(grep('C\tG.*synonymous_variant', df$V1))
length(grep('C\tG.*synonymous_variant', df_context4fold$V1[which(df_context4fold$context == "A")]))
length(grep('C\tG.*synonymous_variant', df_context4fold$V1[which(df_context4fold$context == "T")]))
length(grep('C\tG.*synonymous_variant', df_context4fold$V1[which(df_context4fold$context == "G")]))
length(grep('C\tG.*synonymous_variant', df_context4fold$V1[which(df_context4fold$context == "C")]))