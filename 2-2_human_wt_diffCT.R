library(ggpubr)
library(ggplot2)

# Load data
raw_dataset = read.delim("data/Human_stress_polysome/interdependences_morethan2samples.tsv")

# Keep total
wt_dataset = raw_dataset[grepl("control_",raw_dataset$cond),]

# Reformat
dataset = c()
for (r in rownames(wt_dataset)){
  tempdata = data.frame(row.names = 1:wt_dataset[r,"X.samples"])
  tempdata$reference = substr(wt_dataset[r,"reference"],14,100)
  tempdata$pair = paste0(wt_dataset[r,c("canon_var1","canon_var2")],collapse="-")
  tempdata$cond = strsplit(wt_dataset[r,"cond"],"_")[[1]][2]
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
diffcond = compare_means(odds ~ wobble, data = dataset, group.by = c("pair","aa","cond"),
                         method = "t.test", p.adjust.method="fdr")

# Focus on charging
chrgsubset = diffcond[sapply(diffcond$pair,function(x) strsplit(x,"-")[[1]][2])=="Charged",]

# Plot per AA family
ggplot(dataset[(apply(dataset[,c("aa","pair")],1,paste0,collapse="_") %in% apply(chrgsubset[,c("aa","pair")],1,paste0,collapse="_")),], aes(x=cond, y=odds, fill=wobble)) + 
  facet_wrap( ~ aa*pair, ncol=6, scales = "free") +
  geom_boxplot(lwd=0.5, notch=F, outlier.shape = NA , na.rm=T,alpha=0.2, 
               position = position_dodge(width=0.9)) +
  geom_jitter(na.rm=T , alpha=0.5, size=1, aes(color = wobble),position=position_jitterdodge(dodge.width=0.9)) +
  stat_compare_means(method = "t.test", label = "p.signif") +
  theme_classic()
ggsave("plots/differences_wt_wobble_AAfamilies.pdf",width=16,height=20)

# Plot all together
ggplot(dataset[(apply(dataset[,c("aa","pair")],1,paste0,collapse="_") %in% apply(chrgsubset[,c("aa","pair")],1,paste0,collapse="_")),], aes(x=cond, y=odds, fill=wobble)) + 
  facet_wrap( ~ pair, ncol=6, scales = "free") +
  geom_boxplot(lwd=0.5, notch=F, outlier.shape = NA , na.rm=T,alpha=0.2, 
               position = position_dodge(width=0.9)) +
  geom_jitter(na.rm=T , alpha=0.5, size=1, aes(color = wobble),position=position_jitterdodge(dodge.width=0.9)) +
  stat_compare_means(method = "t.test", label = "p.signif") +
  theme_classic()
ggsave("plots/differences_wt_wobble.pdf",width=16,height=10)
subdataset = dataset[dataset$TC %in% c("C","T"),]
ggplot(subdataset[(apply(subdataset[,c("aa","pair")],1,paste0,collapse="_") %in% apply(chrgsubset[,c("aa","pair")],1,paste0,collapse="_")),], aes(x=cond, y=odds, fill=TC)) + 
  facet_wrap( ~ pair, ncol=6, scales = "free") +
  geom_boxplot(lwd=0.5, notch=F, outlier.shape = NA , na.rm=T,alpha=0.2, 
               position = position_dodge(width=0.9)) +
  geom_jitter(na.rm=T , alpha=0.5, size=1, aes(color = TC),position=position_jitterdodge(dodge.width=0.9)) +
  stat_compare_means(method = "t.test", label = "p.signif") +
  theme_classic()
ggsave("plots/differences_wt_TvsC.pdf",width=16,height=10)


### Differences between conditions ###
diffcond = compare_means(odds ~ cond, data = dataset, group.by = c("pair","reference"),
                         ref.group = "total", method = "t.test", p.adjust.method="fdr")
# Focus on charging
chrgsubset = diffcond[sapply(diffcond$pair,function(x) strsplit(x,"-")[[1]][2])=="Charged",]
sigsubset = chrgsubset[order(chrgsubset$p)[1:15],]

# Plot
ggplot(dataset[(apply(dataset[,c("reference","pair")],1,paste0,collapse="_") %in% apply(sigsubset[,c("reference","pair")],1,paste0,collapse="_")),], aes(x=cond, y=odds)) + 
  facet_wrap( ~ reference*pair, ncol=5, scales = "free") +
  geom_point(shape = 21, color="black", stroke=1, size=2, aes(fill = fisher)) +
  scale_fill_manual(values=c("white", "black")) +
  stat_summary(fun.data=mean_sdl, geom="pointrange", color="red") +
  stat_compare_means(ref.group = "total", method = "t.test", label = "p.signif") +
  theme_classic()
#ggsave("plots/differences_wt_totalVSpoly_charging.pdf",width=15,height=9)


### Differences between isodecoders ###
# Focus on 58-charging and total
chrgsubset = dataset[(dataset$pair=="58-Charged")&(dataset$cond=="total"),]

# Plot per AA family
ggplot(chrgsubset, aes(x=reference, y=odds)) + 
  facet_wrap( ~ aa*pair, ncol=6, scales = "free") +
  geom_point(shape = 21, color="black", stroke=1, size=2, aes(fill = fisher)) +
  scale_fill_manual(values=c("white", "black")) +
  stat_summary(fun.data=mean_sdl, geom="pointrange", color="red") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("plots/differences_wt_total_58Charged_AAfamilies.pdf",width=16,height=16)
