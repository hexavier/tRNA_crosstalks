library(ggpubr)
library(ggplot2)

# Load data
raw_dataset = read.delim("data/QuantM_mouse_tissues/interdependences_morethan2samples.tsv")

# Reformat
dataset = c()
for (r in rownames(raw_dataset)){
  tempdata = data.frame(row.names = 1:raw_dataset[r,"X.samples"])
  tempdata$reference = substr(raw_dataset[r,"reference"],14,100)
  tempdata$pair = paste0(raw_dataset[r,c("canon_var1","canon_var2")],collapse="-")
  tempdata$cond = raw_dataset[r,"cond"]
  odds = sapply(strsplit(raw_dataset[r,"odds"],"\\[|\\]| ")[[1]], function(x) if (x!=""){as.numeric(x)}else{NA})
  tempdata$odds = odds[!is.na(odds)]
  pval = sapply(strsplit(raw_dataset[r,"pval_cor"],"\\[|\\]| ")[[1]], function(x) if (x!=""){as.numeric(x)}else{NA})
  tempdata$fisher = pval[!is.na(pval)]<0.05
  splitted = strsplit(raw_dataset[r,"reference"],"-")[[1]]
  tempdata$wobble = if(substr(splitted[length(splitted)-1],1,1) %in% c("A","T")){"AT"}else{"GC"}
  tempdata$aa = splitted[length(splitted)-2]
  dataset = rbind(dataset,tempdata)
}
# Remove mitochondrial
dataset = dataset[!grepl("mito",dataset$reference),]

### Differences between conditions ###
diffcond = compare_means(odds ~ cond, data = dataset, group.by = c("pair","reference"),
                         method = "t.test", p.adjust.method="fdr")

# Compute odds ratio fold change
sigsubset = diffcond[order(diffcond$p)[1:20],]

# Plot
ggplot(dataset[(apply(dataset[,c("reference","pair")],1,paste0,collapse="_") %in% apply(sigsubset[,c("reference","pair")],1,paste0,collapse="_")),], aes(x=cond, y=odds)) + 
  facet_wrap( ~ reference*pair, ncol=5, scales = "free") +
  geom_point(shape = 21, color="black", stroke=1, size=2, aes(fill = fisher)) +
  scale_fill_manual(values=c("white", "black")) +
  stat_summary(fun.data=mean_sdl, geom="pointrange", color="red") +
  stat_compare_means(ref.group = "Heart", method = "t.test", label = "p.signif") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))
ggsave("plots/differences_mouse_tissues.pdf",width=16,height=10)

### Differences between isodecoders ###
# Focus on Glu-TTC
subsetGlu = dataset[dataset$reference %in% c("tRNA-Glu-TTC-1","tRNA-Glu-TTC-2"),]
subsetGlu$cond = factor(as.character(subsetGlu$cond), levels=c("Cerebellum","Cortex","MedullaOblongata","SpinalCord","Tibialis","Liver","Heart"))

# Plot per AA family
diffGlu = compare_means(odds ~ cond, data = subsetGlu, group.by = c("pair","reference"),
                         method = "anova", p.adjust.method="fdr")
posthoc = compare_means(odds ~ cond, data = subsetGlu[(subsetGlu$reference %in% diffGlu[diffGlu$p<0.05,"reference"])&(subsetGlu$pair %in% diffGlu[diffGlu$p<0.05,"pair"]),], group.by = c("pair","reference"),
                        method = "t.test", p.adjust.method="fdr")
ggplot(subsetGlu, aes(x=cond, y=odds)) + 
  facet_wrap( ~ reference*pair, ncol=6, scales = "free") +
  geom_point(shape = 21, color="black", stroke=1, size=2, aes(fill = fisher)) +
  geom_hline(yintercept=0) +
  scale_fill_manual(values=c("white", "black")) +
  stat_summary(fun.data=mean_sdl, geom="pointrange", color="red", alpha = 0.5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))
ggplot(subsetGlu, aes(x=cond, y=odds)) + 
  facet_wrap( ~ reference*pair, ncol=6, scales = "free") +
  geom_point(shape = 16, size=2) +
  geom_hline(yintercept=0) +
  scale_fill_manual(values=c("white", "black")) +
  stat_summary(fun.data=mean_sdl, geom="pointrange", color="red", alpha = 0.5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))
ggsave("plots/differences_mouse_tissues_GluTTC.pdf",width=11,height=2.3)
