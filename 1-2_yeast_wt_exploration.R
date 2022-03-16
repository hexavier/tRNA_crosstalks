library(ggplot2)
library(reshape)
library(gplots)

# Load data
sigdataset = read.delim("data/Scer_mimseq/interdependences_significant.tsv")
# Subset wt data with al least 2/3 significant samples
wt_data = sigdataset[(sigdataset$cond=="WT")&(sigdataset$X.samples>1),]

### Histogram ###
wt_data$pair = apply(wt_data,1,function(x) paste0(x[c("canon_var1","canon_var2")],collapse="-"))
paircounts = data.frame(table(wt_data$pair))
paircounts$Var1 = factor(x = as.character(paircounts$Var1), levels = as.character(unique(paircounts[order(-paircounts$Freq),"Var1"])))
ggplot(paircounts,aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color="black"),
        axis.text.y=element_text(colour="black"))

### Heatmap of all pairs ###
heatdf = as.matrix(cast(wt_data, pair ~ reference, mean, value = 'avg_odds'))
colnames(heatdf) = substr(colnames(heatdf),26,100)
# Keep only known mod sites
heatdf = heatdf[sapply(rownames(heatdf),function(x) all(strsplit(x,"-")[[1]] %in% c(" 9","26","32","34","37","46","58","Charged"))),]
heatdf = heatdf[,colSums(!is.na(heatdf))>0]
# Plot
pdf(file="plots/heatmap_yeast_wt_pairs.pdf", width=5.5, height=6.5)
heatmap.2(heatdf, col=c("#304d63ff", "#f2d096ff"), Rowv = F, Colv = F, margins=c(9,6), trace="none",
          colsep=1:nrow(heatdf),rowsep=1:nrow(heatdf),sepcolor = "black")
dev.off()
