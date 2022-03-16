library(ggpubr)
library(ggplot2)

# Load data
raw_dataset = read.delim("data/hek293/dmtrnaseq/interdependences_morethan2samples.tsv")

# Keep total
wt_dataset = raw_dataset[grepl("dmseq-untreated",raw_dataset$cond),]

# Reformat
dataset = c()
for (r in rownames(wt_dataset)){
  tempdata = data.frame(row.names = 1:wt_dataset[r,"X.samples"])
  tempdata$reference = substr(wt_dataset[r,"reference"],14,100)
  tempdata$pair = paste0(wt_dataset[r,c("canon_var1","canon_var2")],collapse="-")
  odds = sapply(strsplit(wt_dataset[r,"odds"],"\\[|\\]| ")[[1]], function(x) if (x!=""){as.numeric(x)}else{NA})
  tempdata$odds = odds[!is.na(odds)]
  pval = sapply(strsplit(wt_dataset[r,"pval_cor"],"\\[|\\]| ")[[1]], function(x) if (x!=""){as.numeric(x)}else{NA})
  tempdata$fisher = pval[!is.na(pval)]<0.05
  splitted = strsplit(wt_dataset[r,"reference"],"-")[[1]]
  tempdata$wobble = if(substr(splitted[length(splitted)-1],1,1) %in% c("A","T")){"AT"}else{"GC"}
  tempdata$TC = if(substr(splitted[length(splitted)-1],1,1) %in% "T"){"T"}else if(substr(splitted[length(splitted)-1],1,1) %in% "C"){"C"}else{"AG"}
  tempdata$aa = splitted[length(splitted)-2]
  dataset = rbind(dataset,tempdata)
}
# Remove mitochondrial
dataset = dataset[!grepl("mito",dataset$reference),]

### Differences between wobble positions ###
diffcond = compare_means(odds ~ wobble, data = dataset, group.by = c("pair","aa"),
                         method = "t.test", p.adjust.method="fdr")

# Plot per AA family
ggplot(dataset[(apply(dataset[,c("aa","pair")],1,paste0,collapse="_") %in% apply(diffcond[,c("aa","pair")],1,paste0,collapse="_")),], aes(x=pair, y=odds, fill=wobble)) + 
  facet_wrap( ~ aa, ncol=1, scales = "free") +
  geom_boxplot(lwd=0.5, notch=F, outlier.shape = NA , na.rm=T,alpha=0.2, 
               position = position_dodge(width=0.9)) +
  geom_jitter(na.rm=T , alpha=0.5, size=1, aes(color = wobble),position=position_jitterdodge(dodge.width=0.9)) +
  stat_compare_means(method = "t.test", label = "p.signif") +
  theme_classic()
ggsave("plots/differences_wt_wobble_AAfamilies.pdf",width=10,height=20)

# Plot all together
ggplot(dataset[(apply(dataset[,c("aa","pair")],1,paste0,collapse="_") %in% apply(diffcond[,c("aa","pair")],1,paste0,collapse="_")),], aes(x=pair, y=odds, fill=wobble)) + 
  geom_boxplot(lwd=0.5, notch=F, outlier.shape = NA , na.rm=T,alpha=0.2, 
               position = position_dodge(width=0.9)) +
  geom_jitter(na.rm=T , alpha=0.5, size=1, aes(color = wobble),position=position_jitterdodge(dodge.width=0.9)) +
  stat_compare_means(method = "t.test", label = "p.signif") +
  theme_classic()

subdataset = dataset[dataset$TC %in% c("C","T"),]
ggplot(subdataset[(apply(subdataset[,c("aa","pair")],1,paste0,collapse="_") %in% apply(diffcond[,c("aa","pair")],1,paste0,collapse="_")),], aes(x=pair, y=odds, fill=TC)) + 
  geom_boxplot(lwd=0.5, notch=F, outlier.shape = NA , na.rm=T,alpha=0.2, 
               position = position_dodge(width=0.9)) +
  geom_jitter(na.rm=T , alpha=0.5, size=1, aes(color = TC),position=position_jitterdodge(dodge.width=0.9)) +
  stat_compare_means(method = "t.test", label = "p.signif") +
  theme_classic()
