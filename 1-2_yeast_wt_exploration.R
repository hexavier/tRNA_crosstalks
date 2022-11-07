library(ggplot2)
library(reshape)
library(gplots)
library(eulerr)
library(reshape)

# Load data
sigdataset = read.delim("data/Scer_mimseq/interdependences_significant.tsv")
# Subset wt data with al least 2/3 significant samples
wtox1_data = sigdataset[(sigdataset$cond=="WT"),]
wtox0_data = sigdataset[(sigdataset$cond=="WTox0"),]

### Check overlap between WT and WTox0 ###
# Load mod data
zz=gzfile("data/Scer_mimseq/mismatchTable.csv.gz",'rt')
mod_dataset = read.delim(zz)
# WTox1
mod_dataset1 = mod_dataset[mod_dataset$condition=="WT",]
mod_dataset1 = cast(mod_dataset1, isodecoder*canon_pos*bam*condition ~ type, sum, value = 'proportion')
mod_dataset1$mismatches = rowSums(mod_dataset1[,c("A","T","G","C")],na.rm=T)
wtox1 = cast(mod_dataset1, isodecoder*canon_pos ~ condition, mean, value = 'mismatches')
wtox1 = wtox1[wtox1$WT>=0.05,] # keep modifications detected in 5% reads
# WTox0
mod_dataset2 = mod_dataset[mod_dataset$condition=="WTox0",]
mod_dataset2 = cast(mod_dataset2, isodecoder*canon_pos*bam*condition ~ type, sum, value = 'proportion')
mod_dataset2$mismatches = rowSums(mod_dataset2[,c("A","T","G","C")],na.rm=T)
wtox0 = cast(mod_dataset2, isodecoder*canon_pos ~ condition, mean, value = 'mismatches')
wtox0 = wtox0[wtox0$WTox0>=0.05,] # keep modifications detected in 5% reads
# Keep only crosstalk pairs for which modifications are detected in all methods
overlap = list(wtox1 = sprintf("%s|%s",wtox1$isodecoder,wtox1$canon_pos),
               wtox0 = sprintf("%s|%s",wtox0$isodecoder,wtox0$canon_pos))
fit <- euler(overlap)
plot(fit, fills=c("#F2D096","#8FB9AA"), edges=F, legend=F, labels=T, quantities=T,main="Modifications")
all2 = names(table(c(overlap$wtox1,overlap$wtox0))[table(c(overlap$wtox1,overlap$wtox0))==2])
idx = apply(wtox1_data,1,function(x) all(c(sprintf("%s|%s",x["reference"],x["canon_var1"]),sprintf("%s|%s",x["reference"],x["canon_var2"])) %in% all2))
wtox1_data = wtox1_data[idx,]
idx = apply(wtox0_data,1,function(x) all(c(sprintf("%s|%s",x["reference"],x["canon_var1"]),sprintf("%s|%s",x["reference"],x["canon_var2"])) %in% all2))
wtox0_data = wtox0_data[idx,]
# Discretize odds ratio
wtox1_data$bin_odds = wtox1_data$avg_odds>0
wtox0_data$bin_odds = wtox0_data$avg_odds>0
# Plot overlap
overlap = list(wtox1 = sprintf("%s|%s-%s|%s",wtox1_data$reference,wtox1_data$canon_var1,wtox1_data$canon_var2,wtox1_data$bin_odds),
               wtox0 = sprintf("%s|%s-%s|%s",wtox0_data$reference,wtox0_data$canon_var1,wtox0_data$canon_var2,wtox0_data$bin_odds))
fit <- euler(overlap)
plot(fit, fills=c("#F2D096","#8FB9AA"), edges=F, legend=F, labels=T, quantities=T, main="Crosstalks")

### MERGE CONTROLS ###
# As they highly overlap, merge the two controls
all_data = sigdataset[grepl("WT",sigdataset$cond),]
all_data$bin_odds = all_data$avg_odds>0
wt_data = cast(all_data, reference*canon_var1*canon_var2*bin_odds ~ cond, sum, value = "X.samples")
wt_data = as.data.frame(wt_data[rowSums(wt_data[,c("WT","WTox0")])>=2,])

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
heatdf = as.matrix(cast(wt_data, pair ~ reference, mean, value = 'bin_odds'))
colnames(heatdf) = substr(colnames(heatdf),26,100)
# Keep only known mod sites
heatdf = heatdf[sapply(rownames(heatdf),function(x) all(strsplit(x,"-")[[1]] %in% c(" 9","26","32","34","37","58","Charged"))),]
heatdf = heatdf[,colSums(!is.na(heatdf))>0]
# Plot
pdf(file="plots/heatmap_yeast_wt_pairs.pdf", width=6.5, height=5.5)
heatmap.2(heatdf, col=c("#304d63ff", "#f2d096ff"), Rowv = F, Colv = F, margins=c(8,6), trace="none",
          colsep=1:nrow(heatdf),rowsep=1:nrow(heatdf),sepcolor = "black")
dev.off()
