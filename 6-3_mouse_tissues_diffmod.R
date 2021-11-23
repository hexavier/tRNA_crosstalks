library(ggpubr)
library(ggplot2)
library(gplots)
library(reshape)

# Load data
mod_dataset = read.delim("data/QuantM_mouse_tissues/mismatchTable.csv")
dataset = cast(mod_dataset, isodecoder*canon_pos*bam*condition ~ type, sum, value = 'proportion')
dataset$mismatches = apply(dataset[,c("A","T","G","C")],1,sum,na.rm=T)
dataset$reference = substr(dataset$isodecoder,14,100)

# Differences between conditions #
dataset = dataset[dataset$canon_pos %in% c("58"),]
diffcond = compare_means(mismatches ~ condition, data = dataset, group.by = c("reference","canon_pos"),
                         method = "t.test", p.adjust.method="fdr")

# Compute odds ratio fold change
diffcond$delta = apply(diffcond,1,function(x) mean(dataset[(dataset$reference %in% x["reference"])&(dataset$condition %in% x["group2"])&(dataset$canon_pos %in% x["canon_pos"]),"mismatches"]) - mean(dataset[(dataset$reference %in% x["reference"])&(dataset$condition %in% x["group1"])&(dataset$canon_pos %in% x["canon_pos"]),"mismatches"]))
sigsubset = diffcond[order(diffcond$p)[1:80],]

# Plot
ggplot(dataset[(apply(dataset[,c("reference","canon_pos")],1,paste0,collapse="_") %in% apply(sigsubset[,c("reference","canon_pos")],1,paste0,collapse="_")),], aes(x=condition, y=mismatches)) + 
  facet_wrap( ~ reference*canon_pos, ncol=5, scales = "free") +
  geom_point(shape = 21, color="black", stroke=1, size=2, fill = "black") +
  stat_summary(fun.data=mean_sdl, geom="pointrange", color="red") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))
ggsave("plots/differences_mouse_tissues_58mods.pdf",width=15,height=12)

### Focus on 58 ###
data58 = dataset[dataset$canon_pos=="58",]

# Build heatmap matrix
mergedf = data.frame(row.names = unique(data58$reference))
for (c in unique(data58$bam)){
  idxname = paste0(c(strsplit(strsplit(c,"/")[[1]][2],"\\.")[[1]][1],unique(data58[data58$bam==c,"condition"])[1]),collapse="-")
  mergedf[,idxname] = sapply(rownames(mergedf),function(x) data58[(data58$reference %in% x)&(data58$bam %in% c),"mismatches"])
}
# Keep only those with at least 5% mismatch
mergedf = mergedf[rowMeans(mergedf)>0.05,]

# Plot
mergedf = mergedf[order(rowMeans(mergedf),decreasing = T),] # order based on avg
scaled = apply(mergedf,1,scale)
rownames(scaled) = colnames(mergedf)
# Heatmaps of % mod and scaled
hclust.ward <- function(x) hclust(x, method="ward.D2")
pdf(file="plots/mouse_tissues_mod58_heatmap.pdf", width=14, height=8)
heatmap.2(scaled, col=bluered, trace="none",margins=c(8,15), hclustfun=hclust.ward)
heatmap.2(t(mergedf), col=bluered, Rowv=F, Colv=F, trace="none",margins=c(8,15))
dev.off()
