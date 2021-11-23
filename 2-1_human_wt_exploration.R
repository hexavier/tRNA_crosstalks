library(ggplot2)
library(reshape)
library(gplots)

# Load data
sigdataset = read.delim("data/Human_stress_polysome/interdependences_significant.tsv")
# Subset wt data with al least 2/3 significant samples
wt_data = sigdataset[(sigdataset$cond=="control_total")&(sigdataset$X.samples>1),]

### Histogram ###
wt_data$pair = apply(wt_data,1,function(x) paste0(x[c("canon_var1","canon_var2")],collapse="-"))
paircounts = data.frame(table(wt_data$pair))
paircounts$Var1 = factor(x = as.character(paircounts$Var1), levels = as.character(unique(paircounts[order(-paircounts$Freq),"Var1"])))
ggplot(paircounts[paircounts$Freq>2,],aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color="black"),
        axis.text.y=element_text(colour="black"))
ggsave("plots/histogram_human_wt_pairs.pdf",width=7,height=2.5)

### Odds ratio distribution 58-Charged ###
subset_Ch58 = wt_data[wt_data$pair=="58-Charged",]
# Reformat
Ch58_data = c()
for (r in rownames(subset_Ch58)){
  tempdata = data.frame(row.names = 1:subset_Ch58[r,"X.samples"])
  tempdata$reference = substr(subset_Ch58[r,"reference"],14,100)
  odds = sapply(strsplit(subset_Ch58[r,"odds"],"\\[|\\]| ")[[1]], function(x) if (x!=""){as.numeric(x)}else{NA})
  tempdata$odds = odds[!is.na(odds)]
  tempdata$avg_odds = subset_Ch58[r,"avg_odds"]
  Ch58_data = rbind(Ch58_data,tempdata)
}
Ch58_data$reference = factor(x = as.character(Ch58_data$reference), levels = as.character(unique(Ch58_data[order(Ch58_data$avg_odds),"reference"])))
# Plot
ggplot(Ch58_data,aes(x=reference, y=odds)) +
  geom_point(shape = 21, stroke=0.1, size=2,fill="#f1c232ff") +
  geom_hline(yintercept=0) +
  ylim(-1.05,1.05) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color="black"),
        axis.text.y=element_text(colour="black"))
ggsave("plots/distribution_Chrg58_human_wt_pairs.pdf",width=7,height=3)

### Heatmap of all pairs ###
wt_data = sigdataset[(sigdataset$cond=="control_total")&(sigdataset$X.samples>0),]
wt_data$pair = apply(wt_data,1,function(x) paste0(x[c("canon_var1","canon_var2")],collapse="-"))
heatdf = as.matrix(cast(wt_data, pair ~ reference, mean, value = 'avg_odds'))
colnames(heatdf) = substr(colnames(heatdf),14,100)
pdf(file="plots/heatmap_human_wt_pairs.pdf", width=25, height=25)
heatmap.2(heatdf, col=c("blue","red"), Rowv = F, Colv = F,trace="none",margins=c(8,6),
          colsep=1:nrow(heatdf),rowsep=1:nrow(heatdf),sepcolor = "black")
dev.off()