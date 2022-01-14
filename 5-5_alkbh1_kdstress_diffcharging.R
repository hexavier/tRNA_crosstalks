library(ggpubr)
library(ggplot2)
library(gplots)
library(patchwork)
library(reshape)
library(ggrepel)

# Load data
raw_dataset = read.delim("data/ALKBH1_KD_stress/CCAcounts.csv")

# Compute % charging
samples = unique(raw_dataset$sample)
dataset = c()
for (s in samples){
  refs = unique(raw_dataset[raw_dataset$sample==s,"gene"])
  tempdata = data.frame(row.names = 1:length(refs))
  tempdata$reference = refs
  tempdata$cond = raw_dataset[raw_dataset$sample==s,"condition"][1]
  for (r in rownames(tempdata)){
    subdata = raw_dataset[(raw_dataset$sample==s)&(raw_dataset$gene==tempdata[r,"reference"]),]
    cca = subdata[subdata$end=="CA","count"]
    cc = subdata[subdata$end=="CC","count"]
    tempdata[r,"charging"] = cca/(cca+cc)
  }
  dataset = rbind(dataset,tempdata)
}


### Differences between conditions ###
diffcond = compare_means(charging ~ cond, data = dataset, group.by = c("reference"),
                         ref.group = "minus-AsO2_control-shRNA", method = "t.test", p.adjust.method="fdr")

# Compute odds ratio fold change
diffcond$delta = apply(diffcond,1,function(x) mean(dataset[(dataset$reference==x["reference"])&(dataset$cond==x["group2"]),"charging"]) - mean(dataset[(dataset$reference==x["reference"])&(dataset$cond==x["group1"]),"charging"]))
sigsubset = diffcond[order(diffcond$p)[1:15],]

# Plot
ggplot(dataset[dataset$reference %in% sigsubset$reference,], aes(x=cond, y=charging)) + 
  facet_wrap( ~ reference, ncol=5, scales = "free") +
  geom_point(shape = 21, color="black", stroke=1, size=2, fill = "black") +
  stat_summary(fun.data=mean_sdl, geom="pointrange", color="red") +
  stat_compare_means(ref.group = "minus-AsO2_control-shRNA", method = "t.test", label = "p.signif") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))
ggsave("plots/differences_alkbh1_kdstress_charging.pdf",width=15,height=12)

# Heatmap
mergedf = data.frame(row.names = unique(dataset$reference))
for (c in unique(diffcond$group2)){
  mergedf[,c] = sapply(rownames(mergedf),function(x) as.numeric(diffcond[(diffcond$reference %in% x)&(diffcond$group2 %in% c),"delta"]))
}
mergedf[is.infinite(as.matrix(mergedf))] <- NA
mergedf = mergedf[rowSums(is.na(mergedf))==0,]
rownames(mergedf) = substr(rownames(mergedf),14,100)
# Plot
mergedf = mergedf[order(mergedf$`minus-AsO2_shRNAW`,decreasing = T),] # order based on KD
# Add info about prior knowledge
known58 = c("Ala-CGC","Ala-TGC","Asn-GTT","Cys-GCA","Gln-CTG","Glu-CTC","Glu-TTC","Gly-GCC","His-GTG","Leu-AAG","Lys-TTT","Tyr-GTA","Val-AAC","Val-CAC")
known58not = c("Phe-GAA","Val-TAC")
prior = sapply(rownames(mergedf),function(x) if("mito" %in% x){"grey"}else if(paste0(strsplit(x,"-")[[1]][2:3],collapse="-") %in% known58){"red"}else if(paste0(strsplit(x,"-")[[1]][2:3],collapse="-") %in% known58not){"blue"}else{"grey"})

pdf(file="plots/alkbh1_kdstress_charging_heatmap.pdf", width=5, height=25)
heatmap.2(as.matrix(mergedf), col=bluered, trace="none",Rowv=F,Colv=F,margins=c(22,8),RowSideColors=prior)
dev.off()

### Compare with OR ###
# Load modification data
zz=gzfile("data/ALKBH1_KD_stress/mismatchTable.csv.gz",'rt')   
mod_dataset = read.delim(zz)
mod_dataset = cast(mod_dataset, isodecoder*canon_pos*bam*condition ~ type, sum, value = 'proportion')
mod_dataset$mismatches = rowSums(mod_dataset[,c("A","T","G","C")],na.rm=T)
mod_dataset = mod_dataset[mod_dataset$canon_pos=="58",]
# Load interdependence data
sigdataset = read.delim("data/ALKBH1_KD_stress/interdependences_significant.tsv")

### WITHOUT STRESS ###
# Charging differences between conditions
subdataset = dataset[grepl("minus",dataset$cond),]
diffcond = compare_means(charging ~ cond, data = subdataset, group.by = c("reference"),
                         ref.group = "minus-AsO2_control-shRNA", method = "t.test", p.adjust.method="fdr")
diffcond$delta = apply(diffcond,1,function(x) mean(subdataset[(subdataset$reference==x["reference"])&(subdataset$cond==x["group2"]),"charging"]) - mean(subdataset[(subdataset$reference==x["reference"])&(subdataset$cond==x["group1"]),"charging"]))

# Subset wt 58-Charg data
wt_data = sigdataset[(sigdataset$cond=="minus-AsO2_control-shRNA")&(sigdataset$canon_var1=="58")&(sigdataset$canon_var2=="Charged")&(sigdataset$X.samples>0),]

# Modification differences between conditions
mod_subdataset = mod_dataset[grepl("minus",mod_dataset$cond),]
diffcond_mod = compare_means(mismatches ~ condition, data = mod_subdataset, group.by = c("isodecoder"),
                             ref.group = "minus-AsO2_control-shRNA", method = "t.test", p.adjust.method="fdr")
diffcond_mod$delta = apply(diffcond_mod,1,function(x) mean(mod_subdataset[(mod_subdataset$isodecoder %in% x["isodecoder"])&(mod_subdataset$condition %in% x["group2"]),"mismatches"]) - mean(mod_subdataset[(mod_subdataset$isodecoder %in% x["isodecoder"])&(mod_subdataset$condition %in% x["group1"]),"mismatches"]))


# Create heatdf
nonstressdf = data.frame(row.names = wt_data[order(wt_data$avg_odds),"reference"])
nonstressdf$reference = substr(rownames(nonstressdf),14,100)
nonstressdf$oddsWT = wt_data[order(wt_data$avg_odds),"avg_odds"]
for (c in unique(diffcond$group2)){
  nonstressdf[,sprintf("Chrg_%s",c)] = sapply(rownames(nonstressdf),function(x) as.numeric(diffcond[(diffcond$reference==x)&(diffcond$group2==c),"delta"]))
  nonstressdf[,sprintf("Mod_%s",c)] = sapply(rownames(nonstressdf),function(x) as.numeric(diffcond_mod[(diffcond_mod$isodecoder==x)&(diffcond_mod$group2==c),"delta"]))
}
nonstressdf = nonstressdf[rowSums(is.na(nonstressdf))==0,]

# Plot scatter
p1 = ggplot(nonstressdf[nonstressdf$`Mod_minus-AsO2_shRNAW`>0,], aes(x=oddsWT,y=`Chrg_minus-AsO2_shRNAW`,label=reference)) +
  geom_point(shape = 21, stroke=1, size=3, aes(fill = `Mod_minus-AsO2_shRNAW`, color=(oddsWT>0)==(`Chrg_minus-AsO2_shRNAW`>0))) +
  scale_fill_gradient2() +
  scale_color_manual(values=c("white","black")) +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_text_repel() +
  theme_classic() +
  labs(title ="No Stress", x ="log2(OR_WT)", y = "Delta_Chrg")

### WITH STRESS ###
# Charging differences between conditions
subdataset = dataset[grepl("plus",dataset$cond),]
diffcond = compare_means(charging ~ cond, data = subdataset, group.by = c("reference"),
                         ref.group = "plus-AsO2_control-shRNA", method = "t.test", p.adjust.method="fdr")
diffcond$delta = apply(diffcond,1,function(x) mean(subdataset[(subdataset$reference==x["reference"])&(subdataset$cond==x["group2"]),"charging"]) - mean(subdataset[(subdataset$reference==x["reference"])&(subdataset$cond==x["group1"]),"charging"]))

# Subset wt 58-Charg data
wt_data = sigdataset[(sigdataset$cond=="plus-AsO2_control-shRNA")&(sigdataset$canon_var1=="58")&(sigdataset$canon_var2=="Charged")&(sigdataset$X.samples>0),]

# Modification differences between conditions
mod_subdataset = mod_dataset[grepl("plus",mod_dataset$cond),]
diffcond_mod = compare_means(mismatches ~ condition, data = mod_subdataset, group.by = c("isodecoder"),
                             ref.group = "plus-AsO2_control-shRNA", method = "t.test", p.adjust.method="fdr")
diffcond_mod$delta = apply(diffcond_mod,1,function(x) mean(mod_subdataset[(mod_subdataset$isodecoder %in% x["isodecoder"])&(mod_subdataset$condition %in% x["group2"]),"mismatches"]) - mean(mod_subdataset[(mod_subdataset$isodecoder %in% x["isodecoder"])&(mod_subdataset$condition %in% x["group1"]),"mismatches"]))

# Create heatdf
stressdf = data.frame(row.names = wt_data[order(wt_data$avg_odds),"reference"])
stressdf$reference = substr(rownames(stressdf),14,100)
stressdf$oddsWT = wt_data[order(wt_data$avg_odds),"avg_odds"]
for (c in unique(diffcond$group2)){
  stressdf[,sprintf("Chrg_%s",c)] = sapply(rownames(stressdf),function(x) as.numeric(diffcond[(diffcond$reference==x)&(diffcond$group2==c),"delta"]))
  stressdf[,sprintf("Mod_%s",c)] = sapply(rownames(stressdf),function(x) as.numeric(diffcond_mod[(diffcond_mod$isodecoder==x)&(diffcond_mod$group2==c),"delta"]))
}
stressdf = stressdf[rowSums(is.na(stressdf))==0,]

# Plot scatter
p2 = ggplot(stressdf[stressdf$`Mod_plus-AsO2_shRNAV`>0,], aes(x=oddsWT,y=`Chrg_plus-AsO2_shRNAV`,label=reference)) +
  geom_point(shape = 21, stroke=1, size=3, aes(fill = `Mod_plus-AsO2_shRNAV`,color=(oddsWT>0)==(`Chrg_plus-AsO2_shRNAV`>0))) +
  scale_fill_gradient2() +
  scale_color_manual(values=c("white","black")) +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_text_repel() +
  theme_classic() +
  labs(title ="Stress", x ="log2(OR_WT)", y = "Delta_Chrg")

### Combine plots ###
p1 | p2

ggsave("plots/alkbh1_kdstress_charging_oddsWT_allsamples.pdf",width=20,height=8)
