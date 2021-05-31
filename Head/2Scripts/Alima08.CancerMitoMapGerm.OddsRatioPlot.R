library(ggplot2)

data = read.table('../../Body/3Results/Alima08.CancerMitoMapGerm.FinalTablesAndHists.r.txt',
                  header = TRUE)

# filter out ND6 
#data = data[data$gene == 'TwelveGenes',]

one_line = c()

for(i in 1:nrow(data)){ # i = 1
  temp = data[i, ]
  res = fisher.test(matrix(data = c(temp$NumberOfExpectedAaSubst, temp$FreqFrom, 
                                    temp$NumberOfUnexpectedAaSubst, temp$FreqTo), nrow=2))
  one_line = rbind(one_line, c(as.character(temp$ExpectedAminoAcidSubstBias),as.character(temp$gene),as.character(temp$DataSet), 
                               res$estimate, res$conf.int[1], res$conf.int[2]))

}

odds_results = as.data.frame(one_line)

names(odds_results) = c('Subs', 'gene','DataSet','estimate', 'lowCI', 'highCI')
odds_results$estimate = as.numeric(odds_results$estimate)
odds_results$DummyDataset = 0
for (i in 1:nrow(odds_results)) 
{
  if (odds_results$DataSet[i] == 'Cancer') {odds_results$DummyDataset[i] = 1 }
  if (odds_results$DataSet[i] == 'MitoMap') {odds_results$DummyDataset[i] = 2}
  if (odds_results$DataSet[i] == 'Germ') {odds_results$DummyDataset[i] = 3}
}

odds_results = odds_results[order(odds_results$DummyDataset,odds_results$estimate),]
odds_results1 = odds_results[odds_results$gene == 'TwelveGenes',]

# Create labels
boxLabels = as.character(odds_results1$Subs)

# Enter summary data. boxOdds are the odds ratios (calculated elsewhere), boxCILow is the lower bound of the CI, boxCIHigh is the upper bound.

df <- data.frame(
  yAxis = length(boxLabels):1,
  boxOdds = as.numeric(as.character(odds_results1$estimate)),
  boxCILow = as.numeric(as.character(odds_results1$lowCI)),
  boxCIHigh = as.numeric(as.character(odds_results1$highCI))
)


# Plot
pdf("../../Body/4Figures/Alima08.CancerMitoMapGerm.OddsRatioPlot.R.pdf", height = 10, width = 7)
CancerColor = 'orange'
p <- ggplot(df, aes(x = boxOdds, y = yAxis))
p + geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = boxCIHigh, xmin = boxCILow), size = .5, height = .2, color = "gray50") +
  geom_point(size = 3.5, color = CancerColor) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks = df$yAxis, labels = boxLabels) +
  scale_x_continuous(breaks = seq(0,7,1) ) +
  coord_trans(x = "log10") +
  ylab("") +
  xlab("Odds ratio (log scale)") +
  annotate(geom = "text", y =1.1, x = 3.5, label ="", size = 3.5, hjust = 0) + ggtitle("")

###### the same for ND6

odds_results1 = odds_results[odds_results$gene == 'Nd6',]

# Create labels
boxLabels = as.character(odds_results1$Subs)

# Enter summary data. boxOdds are the odds ratios (calculated elsewhere), boxCILow is the lower bound of the CI, boxCIHigh is the upper bound.

df <- data.frame(
  yAxis = length(boxLabels):1,
  boxOdds = as.numeric(as.character(odds_results1$estimate)),
  boxCILow = as.numeric(as.character(odds_results1$lowCI)),
  boxCIHigh = as.numeric(as.character(odds_results1$highCI))
)

# Plot
Color = 'red'
p <- ggplot(df, aes(x = boxOdds, y = yAxis))
p + geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = boxCIHigh, xmin = boxCILow), size = .5, height = .2, color = "gray50") +
  geom_point(size = 3.5, color = Color) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks = df$yAxis, labels = boxLabels) +
  scale_x_continuous(breaks = seq(0,7,1) ) +
  coord_trans(x = "log10") +
  ylab("") +
  xlab("Odds ratio (log scale)") +
  annotate(geom = "text", y =1.1, x = 3.5, label ="", size = 3.5, hjust = 0) + ggtitle("")


dev.off()