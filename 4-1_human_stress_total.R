library(ggpubr)
library(ggplot2)

# Load data
raw_dataset = read.delim("data/Human_stress_polysome/interdependences_morethan2samples.tsv")

# Keep total
total_dataset = raw_dataset[grepl("_total",raw_dataset$cond),]

# Reformat
dataset = c()
for (r in rownames(total_dataset)){
  tempdata = data.frame(row.names = 1:total_dataset[r,"X.samples"])
  tempdata$reference = substr(total_dataset[r,"reference"],14,100)
  tempdata$pair = paste0(total_dataset[r,c("canon_var1","canon_var2")],collapse="-")
  tempdata$cond = strsplit(total_dataset[r,"cond"],"_")[[1]][1]
  odds = sapply(strsplit(total_dataset[r,"odds"],"\\[|\\]| ")[[1]], function(x) if (x!=""){as.numeric(x)}else{NA})
  tempdata$odds = odds[!is.na(odds)]
  pval = sapply(strsplit(total_dataset[r,"pval_cor"],"\\[|\\]| ")[[1]], function(x) if (x!=""){as.numeric(x)}else{NA})
  tempdata$fisher = pval[!is.na(pval)]<0.05
  splitted = strsplit(total_dataset[r,"reference"],"-")[[1]]
  # Record wobble only for 4-codon boxes
  tempdata$codbox = substr(splitted[length(splitted)-1],2,3)
  tempdata$wobble = if(substr(splitted[length(splitted)-1],1,1) %in% c("A","T")){"AT"}else{"GC"}
  tempdata$aa = splitted[length(splitted)-2]
  dataset = rbind(dataset,tempdata)
}
# Remove mitochondrial
dataset = dataset[!grepl("mito",dataset$reference),]

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
  stat_compare_means(ref.group = "control", method = "t.test", label = "p.signif") +
  theme_classic()
ggsave("plots/differences_total_stress.pdf",width=15,height=9)

# FOCUS ON CHARGING
chrgsubset = diffcond[sapply(diffcond$pair,function(x) strsplit(x,"-")[[1]][2])=="Charged",]
sigsubset = chrgsubset[order(chrgsubset$p)[1:15],]

# Plot
ggplot(dataset[(apply(dataset[,c("reference","pair")],1,paste0,collapse="_") %in% apply(sigsubset[,c("reference","pair")],1,paste0,collapse="_")),], aes(x=cond, y=odds)) + 
  facet_wrap( ~ reference*pair, ncol=5, scales = "free") +
  geom_point(shape = 21, color="black", stroke=1, size=2, aes(fill = fisher)) +
  scale_fill_manual(values=c("white", "black")) +
  stat_summary(fun.data=mean_sdl, geom="pointrange", color="red") +
  stat_compare_means(ref.group = "control", method = "t.test", label = "p.signif") +
  theme_classic()
ggsave("plots/differences_total_stress_charging.pdf",width=15,height=9)


### Differences between wobble positions ###
diffcond = compare_means(odds ~ wobble, data = dataset, group.by = c("codbox","pair","aa","cond"),
                         method = "t.test", p.adjust.method="fdr")

# Focus on charging
chrgsubset = diffcond[sapply(diffcond$pair,function(x) strsplit(x,"-")[[1]][2])=="Charged",]

# Plot per AA family
ggplot(dataset[(apply(dataset[,c("aa","pair")],1,paste0,collapse="_") %in% apply(chrgsubset[,c("aa","pair")],1,paste0,collapse="_")),], aes(x=cond, y=odds, fill=wobble)) + 
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
ggsave("plots/differences_total_wobble_AAfamilies.pdf",width=12,height=22)
