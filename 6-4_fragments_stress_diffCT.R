library(ggpubr)
library(ggplot2)

# Load data
raw_dataset = read.delim("results/crosstalks_fragments_morethan2samples.tsv")

# Keep total
total_dataset = raw_dataset[grepl("_total",raw_dataset$cond),]

# Reformat
dataset = c()
for (r in rownames(total_dataset)){
  tempdata = data.frame(row.names = 1:total_dataset[r,"X.samples"])
  tempdata$reference = substr(total_dataset[r,"reference"],14,100)
  tempdata$pair = paste0(unlist(total_dataset[r,c("canon_var1","var2")]),collapse="-")
  tempdata$cond = strsplit(as.character(total_dataset[r,"cond"]),"_")[[1]][1]
  odds = sapply(strsplit(as.character(total_dataset[r,"odds"]),"\\[|\\]| ")[[1]], function(x) if (x!=""){as.numeric(x)}else{NA})
  tempdata$odds = odds[!is.na(odds)]
  pval = sapply(strsplit(as.character(total_dataset[r,"pval_cor"]),"\\[|\\]| ")[[1]], function(x) if (x!=""){as.numeric(x)}else{NA})
  tempdata$fisher = pval[!is.na(pval)]<0.05
  splitted = strsplit(as.character(total_dataset[r,"reference"]),"-")[[1]]
  # Record wobble only for 4-codon boxes
  tempdata$codbox = substr(as.character(splitted[length(splitted)-1]),2,3)
  tempdata$wobble = if(substr(splitted[length(splitted)-1],1,1) %in% c("A","T")){"AT"}else{"GC"}
  tempdata$aa = splitted[length(splitted)-2]
  dataset = rbind(dataset,tempdata)
}

### Differences between conditions ###
diffcond = compare_means(odds ~ cond, data = dataset, group.by = c("pair","reference"),
                         ref.group = "control", method = "t.test", p.adjust.method="fdr")

# Compute odds ratio fold change
diffcond$FC = apply(diffcond,1,function(x) mean(dataset[(dataset$reference==x["reference"])&(dataset$pair==x["pair"])&(dataset$cond==x["group2"]),"odds"]) - mean(dataset[(dataset$reference==x["reference"])&(dataset$pair==x["pair"])&(dataset$cond==x["group1"]),"odds"]))
sigsubset = diffcond[order(diffcond$p)[1:15],]

# Plot
ggplot(dataset[(apply(dataset[,c("reference","pair")],1,paste0,collapse="_") %in% apply(sigsubset[,c("reference","pair")],1,paste0,collapse="_")),], aes(x=cond, y=odds)) + 
  facet_wrap( ~ reference*pair, ncol=5, scales = "free") +
  geom_point(shape = 21, color="black", stroke=1, size=2, aes(fill = fisher)) +
  scale_fill_manual(values=c("white", "black")) +
  stat_summary(fun.data=mean_sdl, geom="pointrange", color="red") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1,color="black"),
        axis.text.y=element_text(colour="black"))

subsetLeu = dataset[grepl("tRNA-Leu-",dataset$reference)&(dataset$pair=="26-30to39"),]
diffLeu = compare_means(odds ~ cond, data = subsetLeu, group.by = c("reference"),
                        method = "anova", p.adjust.method="fdr")
posthoc = compare_means(odds ~ cond, data = subsetLeu[(subsetLeu$reference %in% diffLeu[diffLeu$p<0.05,"reference"]),], group.by = c("reference"),
                        method = "t.test", p.adjust.method="fdr")
ggplot(subsetLeu, aes(x=cond, y=odds)) + 
  facet_wrap( ~ reference, ncol=3, scales = "free") +
  geom_point(shape = 21, color="black", stroke=1, size=2, aes(fill = fisher)) +
  geom_hline(yintercept=0) +
  scale_fill_manual(values=c("white", "black")) +
  stat_summary(fun.data=mean_sdl, geom="pointrange", color="red") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1,color="black"),
        axis.text.y=element_text(colour="black"))
ggsave("plots/differencesCT_fragments_stress.pdf",width=6,height=2.5)


### Differences between wobble positions ###
diffcond = compare_means(odds ~ wobble, data = dataset, group.by = c("codbox","pair","aa","cond"),
                         method = "t.test", p.adjust.method="fdr")

# Focus on m22G26
leu26 = diffcond[sapply(diffcond$pair,function(x) strsplit(x,"-")[[1]][1])=="26",]

# Plot per AA family
ggplot(dataset[(apply(dataset[,c("aa","pair")],1,paste0,collapse="_") %in% apply(leu26[,c("aa","pair")],1,paste0,collapse="_")),], aes(x=cond, y=odds, fill=wobble)) + 
  facet_wrap( ~ codbox*aa*pair, ncol=6, scales = "free") +
  geom_boxplot(lwd=0.5, notch=F, outlier.shape = NA , na.rm=T,alpha=0.3, 
               position = position_dodge(width=0.9)) +
  scale_fill_manual(values=c("#ed8975ff", "#8fb9aaff")) +
  scale_color_manual(values=c("#ed8975ff", "#8fb9aaff")) +
  geom_jitter(na.rm=T , alpha=1, size=1, aes(color = wobble),position=position_jitterdodge(dodge.width=0.9)) +
  stat_compare_means(method = "t.test", label = "p.signif") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1,color="black"),
        axis.text.y=element_text(colour="black"))
ggsave("plots/differences_fragments_wobble_AAfamilies.pdf",width=12,height=10)
