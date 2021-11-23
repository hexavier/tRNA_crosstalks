library(ggpubr)
library(ggplot2)

# Load data
raw_dataset = read.delim("data/Human_stress_polysome/interdependences_morethan2samples.tsv")

# Find unique pairs in control
uniqpairs = unique(raw_dataset[grepl("control_",raw_dataset$cond),c("reference","canon_var1","canon_var2","cond")])

# Reformat
dataset = c()
for (r in rownames(uniqpairs)){
  control = strsplit(uniqpairs[r,"cond"],"_")[[1]][2]
  tempdata = raw_dataset[raw_dataset$reference==uniqpairs[r,"reference"]&
                           raw_dataset$canon_var1==uniqpairs[r,"canon_var1"]&
                           raw_dataset$canon_var2==uniqpairs[r,"canon_var2"]&
                           grepl(control,raw_dataset$cond),]
  controldata = tempdata[grep("control_",tempdata$cond),]
  stressdata = tempdata[-grep("control_",tempdata$cond),]
  if (nrow(stressdata)>0){
    # Arrange data into new dataframe
    tomerge = data.frame(row.names = 1:nrow(stressdata))
    tomerge$reference = substr(uniqpairs[r,"reference"],14,100)
    tomerge$pair = paste0(uniqpairs[r,c("canon_var1","canon_var2")],collapse="-")
    splitted = strsplit(uniqpairs[r,"reference"],"-")[[1]]
    tomerge$wobble = if(substr(splitted[length(splitted)-1],1,1) %in% c("A","T")){"AT"}else{"GC"}
    tomerge$aa = splitted[length(splitted)-2]
    tomerge$cond = stressdata$cond
    # FC wrt control
    tomerge$FCodds = stressdata$avg_odds - controldata$avg_odds
    # Merge
    dataset = rbind(dataset,tomerge)
  }
}
# Remove mitochondrial
dataset = dataset[!grepl("mito",dataset$reference),]

# Focus on charging
chrgsubset = dataset[sapply(dataset$pair,function(x) strsplit(x,"-")[[1]][2])=="Charged",]

# Keep pairs with more than 3 conditions
filtered = c()
countpairs = table(apply(chrgsubset[,c("pair","aa")],1,function(x)paste0(x,collapse="_")))
for (p in names(countpairs[countpairs>=18])){
  splitted = strsplit(p,"_")[[1]]
  if (length(unique(chrgsubset[(chrgsubset$pair==splitted[1])&(chrgsubset$aa==splitted[2]),"wobble"]))>1){
    filtered = rbind(filtered,chrgsubset[(chrgsubset$pair==splitted[1])&(chrgsubset$aa==splitted[2]),]) 
  }
}


# Plot
ggplot(filtered, aes(x=cond, y=FCodds, fill=wobble)) + 
  facet_wrap( ~ pair*aa, ncol=5, scales = "free") +
  geom_boxplot(lwd=0.5, notch=F, outlier.shape = NA , na.rm=T,alpha=0.2, 
               position = position_dodge(width=0.9)) +
  geom_jitter(na.rm=T , alpha=0.5, size=1, aes(color = wobble),position=position_jitterdodge(dodge.width=0.9)) +
  stat_compare_means(method = "t.test", label = "p.signif") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))
ggsave("plots/differences_stress_overall_wobble_AAfamilies.pdf",width=15,height=20)

ggplot(filtered, aes(x=cond, y=FCodds, fill=wobble)) + 
  facet_wrap( ~ pair, ncol=5, scales = "free") +
  geom_boxplot(lwd=0.5, notch=F, outlier.shape = NA , na.rm=T,alpha=0.2, 
               position = position_dodge(width=0.9)) +
  geom_jitter(na.rm=T , alpha=0.5, size=1, aes(color = wobble),position=position_jitterdodge(dodge.width=0.9)) +
  stat_compare_means(method = "t.test", label = "p.signif") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))
ggsave("plots/differences_stress_overall_wobble.pdf",width=15,height=8)
