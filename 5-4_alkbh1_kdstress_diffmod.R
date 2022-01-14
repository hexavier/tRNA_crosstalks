library(ggpubr)
library(ggplot2)
library(gplots)
library(reshape)

# Load data
zz=gzfile("data/ALKBH1_KD_stress/mismatchTable.csv.gz",'rt')   
mod_dataset = read.delim(zz)
dataset = cast(mod_dataset, isodecoder*canon_pos*bam*condition ~ type, sum, value = 'proportion')
dataset$mismatches = apply(dataset[,c("A","T","G","C")],1,sum,na.rm=T)
dataset$reference = substr(dataset$isodecoder,14,100)

# Differences between conditions #
dataset = dataset[dataset$canon_pos %in% c("58","e9","32","26","37"),]
diffcond = compare_means(mismatches ~ condition, data = dataset, group.by = c("reference","canon_pos"),
                         ref.group = "minus-AsO2_control-shRNA", method = "t.test", p.adjust.method="fdr")

# Compute odds ratio fold change
diffcond$delta = apply(diffcond,1,function(x) mean(dataset[(dataset$reference %in% x["reference"])&(dataset$condition %in% x["group2"])&(dataset$canon_pos %in% x["canon_pos"]),"mismatches"]) - mean(dataset[(dataset$reference %in% x["reference"])&(dataset$condition %in% x["group1"])&(dataset$canon_pos %in% x["canon_pos"]),"mismatches"]))
sigsubset = diffcond[order(diffcond$p)[1:15],]

# Plot
ggplot(dataset[(apply(dataset[,c("reference","canon_pos")],1,paste0,collapse="_") %in% apply(sigsubset[,c("reference","canon_pos")],1,paste0,collapse="_")),], aes(x=condition, y=mismatches)) + 
  facet_wrap( ~ reference*canon_pos, ncol=5, scales = "free") +
  geom_point(shape = 21, color="black", stroke=1, size=2, fill = "black") +
  stat_summary(fun.data=mean_sdl, geom="pointrange", color="red") +
  stat_compare_means(ref.group = "minus-AsO2_control-shRNA", method = "t.test", label = "p.signif") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))
ggsave("plots/differences_alkbh1_kdstress_mods.pdf",width=15,height=12)

### Focus on 58 ###
# ctrl VS KD - Without stress
mod_data = cast(dataset[grepl("minus",dataset$condition)&(dataset$canon_pos=="58"),],
                isodecoder*canon_pos ~ bam, mean, value = 'mismatches')
mod_data = mod_data[rowMeans(mod_data[,!(colnames(mod_data) %in% c("isodecoder","canon_pos"))],na.rm=T)>0.05,]
wt_samples = unique(dataset[dataset$condition=="minus-AsO2_control-shRNA","bam"])
avg_wt = rowMeans(mod_data[,wt_samples],na.rm=T)
mod_data = mod_data[,!(colnames(mod_data) %in% wt_samples)]
mod_data[,!(colnames(mod_data) %in% c("isodecoder","canon_pos"))] = apply(mod_data[,!(colnames(mod_data) %in% c("isodecoder","canon_pos"))],2,
                                                                          function(x) x-avg_wt)
# Build heatmap matrix
nonstressdf = data.frame(mod_data[,!(colnames(mod_data) %in% c("isodecoder","canon_pos"))],row.names = substr(mod_data$isodecoder,14,100))
# Plot
nonstressdf = nonstressdf[order(rowMeans(nonstressdf),decreasing = T),] # order based on KD
nonstressdf = nonstressdf[rowMeans(abs(nonstressdf))>0.01,]  # keep only isodecoders changing >1%
# Add info about prior knowledge
known58 = c("Ala-CGC","Ala-TGC","Asn-GTT","Cys-GCA","Gln-CTG","Glu-CTC","Glu-TTC","Gly-GCC","His-GTG","Leu-AAG","Lys-TTT","Tyr-GTA","Val-AAC","Val-CAC")
known58not = c("Phe-GAA","Val-TAC")
prior = sapply(rownames(nonstressdf),function(x) if("mito"==substr(x,1,4)){"grey"}else if(paste0(strsplit(x,"-")[[1]][2:3],collapse="-") %in% known58){"red"}else if(paste0(strsplit(x,"-")[[1]][2:3],collapse="-") %in% known58not){"blue"}else{"grey"})
# Plot
pdf(file="plots/alkbh1_kdstress_mod58_heatmap.pdf", width=8, height=5)
heatmap.2(t(nonstressdf), col=bluered, trace="none",Rowv=F,Colv=F,margins=c(10,10),ColSideColors=prior)


# ctrl VS KD - With stress
mod_data = cast(dataset[grepl("plus",dataset$condition)&(dataset$canon_pos=="58"),],
                isodecoder*canon_pos ~ bam, mean, value = 'mismatches')
mod_data = mod_data[rowMeans(mod_data[,!(colnames(mod_data) %in% c("isodecoder","canon_pos"))],na.rm=T)>0.05,]
wt_samples = unique(dataset[dataset$condition=="plus-AsO2_control-shRNA","bam"])
avg_wt = rowMeans(mod_data[,wt_samples],na.rm=T)
mod_data = mod_data[,!(colnames(mod_data) %in% wt_samples)]
mod_data[,!(colnames(mod_data) %in% c("isodecoder","canon_pos"))] = apply(mod_data[,!(colnames(mod_data) %in% c("isodecoder","canon_pos"))],2,
                                                                          function(x) x-avg_wt)
# Build heatmap matrix
stressdf = data.frame(mod_data[,!(colnames(mod_data) %in% c("isodecoder","canon_pos"))],row.names = substr(mod_data$isodecoder,14,100))
# Plot
stressdf = stressdf[order(rowMeans(stressdf),decreasing = T),] # order based on KD
stressdf = stressdf[rowMeans(abs(stressdf))>0.01,] # keep only isodecoders changing >1%
# Add info about prior knowledge
known58 = c("Ala-CGC","Ala-TGC","Asn-GTT","Cys-GCA","Gln-CTG","Glu-CTC","Glu-TTC","Gly-GCC","His-GTG","Leu-AAG","Lys-TTT","Tyr-GTA","Val-AAC","Val-CAC")
known58not = c("Phe-GAA","Val-TAC")
prior = sapply(rownames(stressdf),function(x) if("mito"==substr(x,1,4)){"grey"}else if(paste0(strsplit(x,"-")[[1]][2:3],collapse="-") %in% known58){"red"}else if(paste0(strsplit(x,"-")[[1]][2:3],collapse="-") %in% known58not){"blue"}else{"grey"})


heatmap.2(t(stressdf), col=bluered, trace="none",Rowv=F,Colv=F,margins=c(10,10),ColSideColors=prior)


# stress VS nonstress - Without KD
mod_data = cast(dataset[grepl("control-shRNA",dataset$condition)&(dataset$canon_pos=="58"),],
                isodecoder*canon_pos ~ bam, mean, value = 'mismatches')
mod_data = mod_data[rowMeans(mod_data[,!(colnames(mod_data) %in% c("isodecoder","canon_pos"))],na.rm=T)>0.05,]
wt_samples = unique(dataset[dataset$condition=="minus-AsO2_control-shRNA","bam"])
avg_wt = rowMeans(mod_data[,wt_samples],na.rm=T)
mod_data = mod_data[,!(colnames(mod_data) %in% wt_samples)]
mod_data[,!(colnames(mod_data) %in% c("isodecoder","canon_pos"))] = apply(mod_data[,!(colnames(mod_data) %in% c("isodecoder","canon_pos"))],2,
                                                                          function(x) x-avg_wt)
# Build heatmap matrix
stressdf = data.frame(mod_data[,!(colnames(mod_data) %in% c("isodecoder","canon_pos"))],row.names = substr(mod_data$isodecoder,14,100))
# Plot
stressdf = stressdf[order(rowMeans(stressdf),decreasing = T),] # order based on KD
stressdf = stressdf[rowMeans(abs(stressdf))>0.01,] # keep only isodecoders changing >1%
# Add info about prior knowledge
known58 = c("Ala-CGC","Ala-TGC","Asn-GTT","Cys-GCA","Gln-CTG","Glu-CTC","Glu-TTC","Gly-GCC","His-GTG","Leu-AAG","Lys-TTT","Tyr-GTA","Val-AAC","Val-CAC")
known58not = c("Phe-GAA","Val-TAC")
prior = sapply(rownames(stressdf),function(x) if("mito"==substr(x,1,4)){"grey"}else if(paste0(strsplit(x,"-")[[1]][2:3],collapse="-") %in% known58){"red"}else if(paste0(strsplit(x,"-")[[1]][2:3],collapse="-") %in% known58not){"blue"}else{"grey"})


heatmap.2(t(stressdf), col=bluered, trace="none",Rowv=F,Colv=F,margins=c(10,10),ColSideColors=prior)


# stress VS nonstress - With KD
mod_data = cast(dataset[(!grepl("control-shRNA",dataset$condition))&(dataset$canon_pos=="58"),],
                isodecoder*canon_pos ~ bam, mean, value = 'mismatches')
mod_data = mod_data[rowMeans(mod_data[,!(colnames(mod_data) %in% c("isodecoder","canon_pos"))],na.rm=T)>0.05,]
wt_samples = unique(dataset[dataset$condition=="minus-AsO2_shRNAW","bam"])
avg_wt = rowMeans(mod_data[,wt_samples],na.rm=T)
mod_data = mod_data[,!(colnames(mod_data) %in% wt_samples)]
mod_data[,!(colnames(mod_data) %in% c("isodecoder","canon_pos"))] = apply(mod_data[,!(colnames(mod_data) %in% c("isodecoder","canon_pos"))],2,
                                                                          function(x) x-avg_wt)
# Build heatmap matrix
stressdf = data.frame(mod_data[,!(colnames(mod_data) %in% c("isodecoder","canon_pos"))],row.names = substr(mod_data$isodecoder,14,100))
# Plot
stressdf = stressdf[order(rowMeans(stressdf),decreasing = T),] # order based on KD
stressdf = stressdf[rowMeans(abs(stressdf))>0.01,] # keep only isodecoders changing >1%
# Add info about prior knowledge
known58 = c("Ala-CGC","Ala-TGC","Asn-GTT","Cys-GCA","Gln-CTG","Glu-CTC","Glu-TTC","Gly-GCC","His-GTG","Leu-AAG","Lys-TTT","Tyr-GTA","Val-AAC","Val-CAC")
known58not = c("Phe-GAA","Val-TAC")
prior = sapply(rownames(stressdf),function(x) if("mito"==substr(x,1,4)){"grey"}else if(paste0(strsplit(x,"-")[[1]][2:3],collapse="-") %in% known58){"red"}else if(paste0(strsplit(x,"-")[[1]][2:3],collapse="-") %in% known58not){"blue"}else{"grey"})


heatmap.2(t(stressdf), col=bluered, trace="none",Rowv=F,Colv=F,margins=c(10,10),ColSideColors=prior)
dev.off()
