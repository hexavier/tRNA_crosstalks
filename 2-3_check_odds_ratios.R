library(reshape)
library(gplots)
library(RColorBrewer)

#### CROSSTALKS ####
# Load data
crosstalks = read.delim("data/simulated_reads/interdependences.tsv")

# Record input info
crosstalks$in26 = sapply(crosstalks$sample, function(x) as.numeric(gsub("-",".",substr(strsplit(as.character(x),"_")[[1]][1],5,9))))
crosstalks$in58 = sapply(crosstalks$sample, function(x) as.numeric(gsub("-",".",substr(strsplit(as.character(x),"_")[[1]][2],5,9))))
crosstalks$inodds = sapply(crosstalks$sample, function(x) as.numeric(gsub("-",".",substr(strsplit(as.character(x),"_")[[1]][3],5,9))))

# Mask non-significant OR
crosstalks[crosstalks$pval_corrected>0.05,"odds_ratio"] <- NA

# Plot
pdf(file="plots/simulation_odds_ratio.pdf", width=4, height=4)
for (o in unique(crosstalks$inodds)){
  # Make dataframe
  oddsdf = cast(crosstalks[crosstalks$inodds==o,], in26 ~ in58, sum, value="odds_ratio")
  rownames(oddsdf) = oddsdf$in26; oddsdf = oddsdf[,2:ncol(oddsdf)]
  # Compute errors
  diffdf = abs(log2(oddsdf) - log2(o))
  heatmap.2(as.matrix(diffdf),col=brewer.pal(9,"YlOrBr"),trace="none",Rowv=F,Colv=F,margins=c(5,5),na.color = "grey",
            breaks=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9), main=sprintf("log2(OR)=%.2f",log2(o)),ylab="ModA",xlab="ModB")
}
dev.off()