library(ggplot2)
library(reshape)
library(gplots)
library(patchwork)

# Load data
sigdataset = read.delim("results/crosstalks_fragments_significant.tsv")
# Subset wt data with al least 2 significant samples
wt_data = sigdataset[(sigdataset$cond=="control_total"),]

### Histogram ###
wt_data$pair = apply(wt_data,1,function(x) paste0(x[c("canon_var1","var2")],collapse="-"))
paircounts = data.frame(table(wt_data$pair))
paircounts$Var1 = factor(x = as.character(paircounts$Var1), levels = as.character(unique(paircounts[order(-paircounts$Freq),"Var1"])))
ggplot(paircounts[paircounts$Freq>2,],aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color="black"),
        axis.text.y=element_text(colour="black"))

### Odds ratio distribution m22G26 ###
subset_pair = wt_data[grepl("26-",wt_data$pair),]
# Plot histo
paircounts = data.frame(table(subset_pair$pair))
paircounts$Var1 = factor(x = as.character(paircounts$Var1), levels = as.character(unique(paircounts[order(-paircounts$Freq),"Var1"])))
ggplot(paircounts,aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color="black"),
        axis.text.y=element_text(colour="black"))
ggsave("plots/histogram_fragments_wt_26pairs.pdf",width=2,height=2.5)
# Reformat
pair_data = c()
for (r in rownames(subset_pair)){
  tempdata = data.frame(row.names = 1:subset_pair[r,"X.samples"])
  tempdata$reference = substr(subset_pair[r,"reference"],14,100)
  odds = sapply(strsplit(as.character(subset_pair[r,"odds"]),"\\[|\\]| ")[[1]], function(x) if (x!=""){as.numeric(x)}else{NA})
  tempdata$odds = odds[!is.na(odds)]
  tempdata$avg_odds = subset_pair[r,"avg_odds"]
  tempdata$pair = subset_pair[r,"pair"]
  pair_data = rbind(pair_data,tempdata)
}
# Plot
p1 <- ggplot(pair_data[pair_data$pair=="26-30to39",],aes(x = factor(reference, level = as.character(unique(reference[order(avg_odds)]))), y=odds)) +
  geom_point(shape = 21, stroke=0.1, size=2, fill="#f1c232ff") +
  geom_hline(yintercept=0) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color="black"),
        axis.text.y=element_text(colour="black")) +
  labs(title="26-30to39", x = NULL, y = "log2(OR)")
p2 <- ggplot(pair_data[pair_data$pair=="26-40to49",],aes(x = factor(reference, level = as.character(unique(reference[order(avg_odds)]))), y=odds)) +
  geom_point(shape = 21, stroke=0.1, size=2, fill="#f1c232ff") +
  geom_hline(yintercept=0) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color="black"),
        axis.text.y=element_text(colour="black")) +
  labs(title="26-40to49", x = NULL, y = "log2(OR)")
p3 <- ggplot(pair_data[pair_data$pair=="26-50to59",],aes(x = factor(reference, level = as.character(unique(reference[order(avg_odds)]))), y=odds)) +
  geom_point(shape = 21, stroke=0.1, size=2, fill="#f1c232ff") +
  geom_hline(yintercept=0) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color="black"),
        axis.text.y=element_text(colour="black")) +
  labs(title="26-50to59", x = NULL, y = "log2(OR)")
p1 + p2 + p3 + plot_layout(widths = c(1.5,0.9,1.2))
ggsave("plots/distribution_26_fragments_wt_pairs.pdf",width=8,height=3)

### Heatmap of all pairs ###
wt_data = sigdataset[(sigdataset$cond=="control_total")&(sigdataset$X.samples>0),]
wt_data$pair = apply(wt_data,1,function(x) paste0(x[c("canon_var1","var2")],collapse="-"))
heatdf = as.matrix(cast(wt_data, pair ~ reference, mean, value = 'avg_odds'))
colnames(heatdf) = substr(colnames(heatdf),14,100)
heatmap.2(heatdf, col=c("blue","red"), Rowv = F, Colv = F,trace="none",margins=c(8,6),
          colsep=1:nrow(heatdf),rowsep=1:nrow(heatdf),sepcolor = "black")
