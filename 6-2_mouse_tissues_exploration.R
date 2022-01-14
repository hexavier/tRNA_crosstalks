library(ggplot2)

# Load data
sigdataset = read.delim("data/QuantM_mouse_tissues/interdependences_significant.tsv")
# Subset with at least 2/3 significant samples
data_2samp = sigdataset[(sigdataset$X.samples>=2),]

### Histogram ###
data_2samp$pair = apply(data_2samp,1,function(x) paste0(x[c("canon_var1","canon_var2")],collapse="-"))
paircounts = data.frame(table(data_2samp[,c("cond","pair")]))
paircounts$pair = factor(x = as.character(paircounts$pair), levels = as.character(unique(paircounts[order(-paircounts$Freq),"pair"])))
ggplot(paircounts,aes(x=pair, y=Freq)) +
  facet_wrap( ~ cond, ncol=1) +
  geom_bar(stat="identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color="black"),
        axis.text.y=element_text(colour="black"))

### Odds ratio distribution 58-Charged ###
subset_2526 = data_2samp[data_2samp$pair=="25-26",]
# Reformat
data_2526 = c()
for (r in rownames(subset_2526)){
  tempdata = data.frame(row.names = 1:subset_2526[r,"X.samples"])
  tempdata$reference = substr(subset_2526[r,"reference"],14,100)
  odds = sapply(strsplit(subset_2526[r,"odds"],"\\[|\\]| ")[[1]], function(x) if (x!=""){as.numeric(x)}else{NA})
  tempdata$odds = odds[!is.na(odds)]
  tempdata$avg_odds = subset_2526[r,"avg_odds"]
  tempdata$cond = subset_2526[r,"cond"]
  data_2526 = rbind(data_2526,tempdata)
}
data_2526$reference = factor(x = as.character(data_2526$reference), levels = as.character(unique(data_2526[order(data_2526$avg_odds),"reference"])))
# Plot
ggplot(data_2526,aes(x=reference, y=odds)) +
  facet_wrap( ~ cond, ncol=3) +
  geom_point(shape = 21, stroke=0.1, size=2,fill="#f1c232ff") +
  geom_hline(yintercept=0) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color="black"),
        axis.text.y=element_text(colour="black"))
