
rm(list=ls(all=TRUE))
########

fulltreeCodons <- read_delim("~/mitoclub/human-global-tree/Body/1Raw/fulltreeCodons.csv",
                             ";", escape_double = FALSE, trim_ws = TRUE)


nrow(fulltreeCodons) #577276
######## CLEAN
# Identical context ancestor == descendant
distlist = c()
for (i in seq(1, nrow(fulltreeCodons))) {
  dist = adist(fulltreeCodons$ancestor[i], fulltreeCodons$descendant[i])[1]
  print(i)
  distlist = c(distlist, dist)
}
# Append to the top

######## FILTER IDENTICAL CONTEXT
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
# ALL. A>G
nrow(df[which(df$ancestral_nucl== "A" & df$descendant_nucl == "G"),])
# ALL. A>T
nrow(df[which(df$ancestral_nucl== "A" & df$descendant_nucl == "T"),])
# ALL. A>C
nrow(df[which(df$ancestral_nucl== "A" & df$descendant_nucl == "C"),])
# ALL. T>G
nrow(df[which(df$ancestral_nucl== "T" & df$descendant_nucl == "G"),])
# ALL. T>C
nrow(df[which(df$ancestral_nucl== "T" & df$descendant_nucl == "C"),])
# ALL. T>A
nrow(df[which(df$ancestral_nucl== "T" & df$descendant_nucl == "A"),])
# ALL. G>A
nrow(df[which(df$ancestral_nucl== "G" & df$descendant_nucl == "A"),])
# ALL. G>T
nrow(df[which(df$ancestral_nucl== "G" & df$descendant_nucl == "T"),])
# ALL. G>C
nrow(df[which(df$ancestral_nucl== "G" & df$descendant_nucl == "C"),])
# ALL. C>A
nrow(df[which(df$ancestral_nucl== "C" & df$descendant_nucl == "A"),])
# ALL. C>T
nrow(df[which(df$ancestral_nucl== "C" & df$descendant_nucl == "T"),])
# ALL. C>G
nrow(df[which(df$ancestral_nucl== "C" & df$descendant_nucl == "G"),])

########
# NODLOOP
nodloop <- clean[which(clean$ref_pos > 576 & clean$ref_pos < 16024),]
nrow(nodloop) # 353515

# NODLOOP. A>G
nrow(nodloop[which(nodloop$ancestral_nucl== "A" & nodloop$descendant_nucl == "G"),])
# NODLOOP. A>T
nrow(nodloop[which(nodloop$ancestral_nucl== "A" & nodloop$descendant_nucl == "T"),])
# NODLOOP. A>C
nrow(nodloop[which(nodloop$ancestral_nucl== "A" & nodloop$descendant_nucl == "C"),])
# NODLOOP. T>G
nrow(nodloop[which(nodloop$ancestral_nucl== "T" & nodloop$descendant_nucl == "G"),])
# NODLOOP. T>C
nrow(nodloop[which(nodloop$ancestral_nucl== "T" & nodloop$descendant_nucl == "C"),])
# NODLOOP. T>A
nrow(nodloop[which(nodloop$ancestral_nucl== "T" & nodloop$descendant_nucl == "A"),])
# NODLOOP. G>A
nrow(nodloop[which(nodloop$ancestral_nucl== "G" & nodloop$descendant_nucl == "A"),])
# NODLOOP. G>T
nrow(nodloop[which(nodloop$ancestral_nucl== "G" & nodloop$descendant_nucl == "T"),])
# NODLOOP. G>C
nrow(nodloop[which(nodloop$ancestral_nucl== "G" & nodloop$descendant_nucl == "C"),])
# NODLOOP. C>A
nrow(nodloop[which(nodloop$ancestral_nucl== "C" & nodloop$descendant_nucl == "A"),])
# NODLOOP. C>T
nrow(nodloop[which(nodloop$ancestral_nucl== "C" & nodloop$descendant_nucl == "T"),])
# NODLOOP. C>G
nrow(nodloop[which(nodloop$ancestral_nucl== "C" & nodloop$descendant_nucl == "G"),])


########
# SYNONYMALL
synonymall = df[which(df$synonymous == "synonymous"),]
nrow(synonymall)

# SYNONYMALL. A>G
nrow(synonymall[which(synonymall$ancestral_nucl== "A" & synonymall$descendant_nucl == "G"),])
# SYNONYMALL. A>T
nrow(synonymall[which(synonymall$ancestral_nucl== "A" & synonymall$descendant_nucl == "T"),])
# SYNONYMALL. A>C
nrow(synonymall[which(synonymall$ancestral_nucl== "A" & synonymall$descendant_nucl == "C"),])
# SYNONYMALL. T>G
nrow(synonymall[which(synonymall$ancestral_nucl== "T" & synonymall$descendant_nucl == "G"),])
# SYNONYMALL. T>C
nrow(synonymall[which(synonymall$ancestral_nucl== "T" & synonymall$descendant_nucl == "C"),])
# SYNONYMALL. T>A
nrow(synonymall[which(synonymall$ancestral_nucl== "T" & synonymall$descendant_nucl == "A"),])
# SYNONYMALL. G>A
nrow(synonymall[which(synonymall$ancestral_nucl== "G" & synonymall$descendant_nucl == "A"),])
# SYNONYMALL. G>T
nrow(synonymall[which(synonymall$ancestral_nucl== "G" & synonymall$descendant_nucl == "T"),])
# SYNONYMALL. G>C
nrow(synonymall[which(synonymall$ancestral_nucl== "G" & synonymall$descendant_nucl == "C"),])
# SYNONYMALL. C>A
nrow(synonymall[which(synonymall$ancestral_nucl== "C" & synonymall$descendant_nucl == "A"),])
# SYNONYMALL. C>T
nrow(synonymall[which(synonymall$ancestral_nucl== "C" & synonymall$descendant_nucl == "T"),])
# SYNONYMALL. C>G
nrow(synonymall[which(synonymall$ancestral_nucl== "C" & synonymall$descendant_nucl == "G"),])


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
