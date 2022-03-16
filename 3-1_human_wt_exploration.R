library(ggplot2)
library(reshape)
library(gplots)
library(patchwork)

# Load data
sigdataset = read.delim("data/hek293/dmtrnaseq/interdependences_significant.tsv")
# Subset wt data with al least 2 significant samples
wt_data = sigdataset[(sigdataset$cond=="dmseq-untreated")&(sigdataset$X.samples>1),]

### Histogram ###
wt_data$pair = apply(wt_data,1,function(x) paste0(x[c("canon_var1","canon_var2")],collapse="-"))
paircounts = data.frame(table(wt_data$pair))
paircounts$Var1 = factor(x = as.character(paircounts$Var1), levels = as.character(unique(paircounts[order(-paircounts$Freq),"Var1"])))
ggplot(paircounts[paircounts$Freq>2,],aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color="black"),
        axis.text.y=element_text(colour="black"))
ggsave("plots/histogram_human_wt_pairs.pdf",width=4,height=2.5)

### Odds ratio distribution 37-58 ###
subset_pair = wt_data[wt_data$pair=="37-58",]
# Reformat
pair_data = c()
for (r in rownames(subset_pair)){
  tempdata = data.frame(row.names = 1:subset_pair[r,"X.samples"])
  tempdata$reference = substr(subset_pair[r,"reference"],14,100)
  odds = sapply(strsplit(as.character(subset_pair[r,"odds"]),"\\[|\\]| ")[[1]], function(x) if (x!=""){as.numeric(x)}else{NA})
  tempdata$odds = odds[!is.na(odds)]
  tempdata$avg_odds = subset_pair[r,"avg_odds"]
  pair_data = rbind(pair_data,tempdata)
}
pair_data$reference = factor(x = as.character(pair_data$reference), levels = as.character(unique(pair_data[order(pair_data$avg_odds),"reference"])))
# Plot
p1 = ggplot(pair_data,aes(x=reference, y=odds)) +
  geom_point(shape = 21, stroke=0.1, size=2, fill="#f1c232ff") +
  geom_hline(yintercept=0) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color="black"),
        axis.text.y=element_text(colour="black")) +
  labs(title="37-58", x = NULL, y = "log2(OR)")

### Odds ratio distribution 26-58 ###
subset_pair = wt_data[wt_data$pair=="26-58",]
# Reformat
pair_data = c()
for (r in rownames(subset_pair)){
  tempdata = data.frame(row.names = 1:subset_pair[r,"X.samples"])
  tempdata$reference = substr(subset_pair[r,"reference"],14,100)
  odds = sapply(strsplit(as.character(subset_pair[r,"odds"]),"\\[|\\]| ")[[1]], function(x) if (x!=""){as.numeric(x)}else{NA})
  tempdata$odds = odds[!is.na(odds)]
  tempdata$avg_odds = subset_pair[r,"avg_odds"]
  pair_data = rbind(pair_data,tempdata)
}
pair_data$reference = factor(x = as.character(pair_data$reference), levels = as.character(unique(pair_data[order(pair_data$avg_odds),"reference"])))
# Plot
p2 = ggplot(pair_data,aes(x=reference, y=odds)) +
  geom_point(shape = 21, stroke=0.1, size=2, fill="#f1c232ff") +
  geom_hline(yintercept=0) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color="black"),
        axis.text.y=element_text(colour="black")) + 
  labs(title="26-58", x = NULL, y = "log2(OR)")
p1 | p2
ggsave("plots/distribution_2637-58_human_wt_pairs.pdf",width=9,height=3.2)

### Odds ratio distribution 9-Charged ###
subset_pair = wt_data[wt_data$pair=="9-Charged",]
# Reformat
pair_data = c()
for (r in rownames(subset_pair)){
  tempdata = data.frame(row.names = 1:subset_pair[r,"X.samples"])
  tempdata$reference = substr(subset_pair[r,"reference"],14,100)
  odds = sapply(strsplit(as.character(subset_pair[r,"odds"]),"\\[|\\]| ")[[1]], function(x) if (x!=""){as.numeric(x)}else{NA})
  tempdata$odds = odds[!is.na(odds)]
  tempdata$avg_odds = subset_pair[r,"avg_odds"]
  pair_data = rbind(pair_data,tempdata)
}
pair_data$reference = factor(x = as.character(pair_data$reference), levels = as.character(unique(pair_data[order(pair_data$avg_odds),"reference"])))
# Plot
ggplot(pair_data,aes(x=reference, y=odds)) +
  geom_point(shape = 21, stroke=0.1, size=2, fill="#f1c232ff") +
  geom_hline(yintercept=0) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color="black"),
        axis.text.y=element_text(colour="black"))
ggsave("plots/distribution_9-Charged_human_wt_pairs.pdf",width=7,height=3)

### Heatmap of all pairs ###
wt_data = sigdataset[(sigdataset$cond=="dmseq-untreated")&(sigdataset$X.samples>0),]
wt_data$pair = apply(wt_data,1,function(x) paste0(x[c("canon_var1","canon_var2")],collapse="-"))
heatdf = as.matrix(cast(wt_data, pair ~ reference, mean, value = 'avg_odds'))
colnames(heatdf) = substr(colnames(heatdf),14,100)
pdf(file="plots/heatmap_human_wt_pairs.pdf", width=17, height=17)
heatmap.2(heatdf, col=c("blue","red"), Rowv = F, Colv = F,trace="none",margins=c(8,6),
          colsep=1:nrow(heatdf),rowsep=1:nrow(heatdf),sepcolor = "black")
dev.off()
