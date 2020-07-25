
rm(list=ls(all=TRUE))
########

fulltreeCodons <- read_delim("~/mitoclub/human-global-tree/Body/1Raw/fulltreeCodons.csv", 
                             ";", escape_double = FALSE, trim_ws = TRUE)


nrow(fulltreeCodons) #577276 7😭0😩😔0076YL
######## CLEAN
# Identical context ancestor == descendant
distlist = c()
for (i in seq(1, nrow(fulltreeCodons))) {
  dist = adist(fulltreeCodons$ancestor[i], fulltreeCodons$descendant[i])[1]
  print(i)
  distlist = c(distlist, dist)
}
# Append to the top

fulltreeCodons['distance'] = distlist
identicontext = fulltreeCodons[which(fulltreeCodons$distance == 1),]
nrow(identicontext)

# add column with anc and derived nucl
# add column with preceeding
library(stringr)
identicontext['ancestral_nucl'] = str_extract(identicontext$ancestor, '[A-Z]')
identicontext['descendant_nucl'] = str_extract(identicontext$descendant, '[A-Z]')

# As for preceeding. The nucl of choice is always 3rd, so preceeding must be second
# summary(str_locate(identicontext$ancestor, '[A-Z]'))
identicontext['preceeding'] = substring(identicontext$ancestor, 2, 2)

# make sure it's only a t g c, filter out the rest.
# Before 531352
clean = identicontext[which(identicontext$ancestral_nucl == "A" | identicontext$ancestral_nucl == "T" | identicontext$ancestral_nucl == "G" | identicontext$ancestral_nucl == "C"),]
clean = clean[which(clean$descendant_nucl == "A" | clean$descendant_nucl == "T" | clean$descendant_nucl == "G" | clean$descendant_nucl == "C"),]
clean = clean[which(clean$preceeding == "a" | clean$preceeding == "t" | clean$preceeding == "g" | clean$preceeding == "c"),]
# After 502531

#############
df = clean
# ALL Total number of mutations (rows in table) - 502531
nrow(df)

########
# ALL - All mutations of a certain type
# ALL. A>G - 481
nrow(df[which(df$ancestral_nucl== "A" & df$descendant_nucl == "G"),])
# ALL. A>T - 72
nrow(df[which(df$ancestral_nucl== "A" & df$descendant_nucl == "T"),])
# ALL. A>C - 91
nrow(df[which(df$ancestral_nucl== "A" & df$descendant_nucl == "C"),])
# ALL. T>G - 27
nrow(df[which(df$ancestral_nucl== "T" & df$descendant_nucl == "G"),])
# ALL. T>C - 2120
nrow(df[which(df$ancestral_nucl== "T" & df$descendant_nucl == "C"),])
# ALL. T>A - 29
nrow(df[which(df$ancestral_nucl== "T" & df$descendant_nucl == "A"),])
# ALL. G>A - 3663
nrow(df[which(df$ancestral_nucl== "G" & df$descendant_nucl == "A"),])
# ALL. G>T - 49
nrow(df[which(df$ancestral_nucl== "G" & df$descendant_nucl == "T"),])
# ALL. G>C - 119
nrow(df[which(df$ancestral_nucl== "G" & df$descendant_nucl == "C"),])
# ALL. C>A - 167
nrow(df[which(df$ancestral_nucl== "C" & df$descendant_nucl == "A"),])
# ALL. C>T - 771
nrow(df[which(df$ancestral_nucl== "C" & df$descendant_nucl == "T"),])
# ALL. C>G - 22
nrow(df[which(df$ancestral_nucl== "C" & df$descendant_nucl == "G"),])

# #ALL Nodloop 6628
#> summary(nodloop$position)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#577    3841    7983    8132   12354   16023 
#ref var position  (16024..16569,1..576)
nodloop <- clean[which(clean$ref_pos > 576 & clean$ref_pos < 16024),]
nrow(nodloop) # 353515

# NODLOOP. A>G - 481
nrow(nodloop[which(nodloop$ancestral_nucl== "A" & nodloop$descendant_nucl == "G"),])
# NODLOOP. A>T - 72
nrow(nodloop[which(nodloop$ancestral_nucl== "A" & nodloop$descendant_nucl == "T"),])
# NODLOOP. A>C - 91
nrow(nodloop[which(nodloop$ancestral_nucl== "A" & nodloop$descendant_nucl == "C"),])
# NODLOOP. T>G - 27
nrow(nodloop[which(nodloop$ancestral_nucl== "T" & nodloop$descendant_nucl == "G"),])
# NODLOOP. T>C - 2120
nrow(nodloop[which(nodloop$ancestral_nucl== "T" & nodloop$descendant_nucl == "C"),])
# NODLOOP. T>A - 29
nrow(nodloop[which(nodloop$ancestral_nucl== "T" & nodloop$descendant_nucl == "A"),])
# NODLOOP. G>A - 3663
nrow(nodloop[which(nodloop$ancestral_nucl== "G" & nodloop$descendant_nucl == "A"),])
# NODLOOP. G>T - 49
nrow(nodloop[which(nodloop$ancestral_nucl== "G" & nodloop$descendant_nucl == "T"),])
# NODLOOP. G>C - 119
nrow(nodloop[which(nodloop$ancestral_nucl== "G" & nodloop$descendant_nucl == "C"),])
# NODLOOP. C>A - 167
nrow(nodloop[which(nodloop$ancestral_nucl== "C" & nodloop$descendant_nucl == "A"),])
# NODLOOP. C>T - 771
nrow(nodloop[which(nodloop$ancestral_nucl== "C" & nodloop$descendant_nucl == "T"),])
# NODLOOP. C>G - 22
nrow(nodloop[which(nodloop$ancestral_nucl== "C" & nodloop$descendant_nucl == "G"),])


#Synonymall 874
synonymall = df[which(df$synonymous == "synonymous"),]
nrow(synonymall)

# NODLOOP. A>G - 481
nrow(synonymall[which(synonymall$ancestral_nucl== "A" & synonymall$descendant_nucl == "G"),])
# NODLOOP. A>T - 72
nrow(synonymall[which(synonymall$ancestral_nucl== "A" & synonymall$descendant_nucl == "T"),])
# NODLOOP. A>C - 91
nrow(synonymall[which(synonymall$ancestral_nucl== "A" & synonymall$descendant_nucl == "C"),])
# NODLOOP. T>G - 27
nrow(synonymall[which(synonymall$ancestral_nucl== "T" & synonymall$descendant_nucl == "G"),])
# NODLOOP. T>C - 2120
nrow(synonymall[which(synonymall$ancestral_nucl== "T" & synonymall$descendant_nucl == "C"),])
# NODLOOP. T>A - 29
nrow(synonymall[which(synonymall$ancestral_nucl== "T" & synonymall$descendant_nucl == "A"),])
# NODLOOP. G>A - 3663
nrow(synonymall[which(synonymall$ancestral_nucl== "G" & synonymall$descendant_nucl == "A"),])
# NODLOOP. G>T - 49
nrow(synonymall[which(synonymall$ancestral_nucl== "G" & synonymall$descendant_nucl == "T"),])
# NODLOOP. G>C - 119
nrow(synonymall[which(synonymall$ancestral_nucl== "G" & synonymall$descendant_nucl == "C"),])
# NODLOOP. C>A - 167
nrow(synonymall[which(synonymall$ancestral_nucl== "C" & synonymall$descendant_nucl == "A"),])
# NODLOOP. C>T - 771
nrow(synonymall[which(synonymall$ancestral_nucl== "C" & synonymall$descendant_nucl == "T"),])
# NODLOOP. C>G - 22
nrow(synonymall[which(synonymall$ancestral_nucl== "C" & synonymall$descendant_nucl == "G"),])


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
