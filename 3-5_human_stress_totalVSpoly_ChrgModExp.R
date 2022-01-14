library(ggplot2)
library(patchwork)
library(reshape)
library(ggrepel)
library(ggpubr)
library(ggforce)

### CHARGING ###
# Load data
raw_dataset = read.delim("data/Human_stress_polysome/CCAcounts.csv")

# Compute % charging
samples = unique(raw_dataset$sample)
dataset = c()
for (s in samples){
  refs = unique(raw_dataset[raw_dataset$sample==s,"gene"])
  tempdata = data.frame(row.names = 1:length(refs))
  tempdata$reference = refs
  tempdata$cond = raw_dataset[raw_dataset$sample==s,"condition"][1]
  tempdata$sample = s
  for (r in rownames(tempdata)){
    subdata = raw_dataset[(raw_dataset$sample==s)&(raw_dataset$gene==tempdata[r,"reference"]),]
    cca = subdata[subdata$end=="CA","count"]
    cc = subdata[subdata$end=="CC","count"]
    tempdata[r,"charging"] = cca/(cca+cc)
  }
  dataset = rbind(dataset,tempdata)
}

### Subset tRNAs to analyze ###
# Isodecoders
isodecoders = unique(dataset[!grepl("mito",dataset$reference),"reference"])
# Pairs subset
modpos = sprintf("%s_58",isodecoders)
allpos = c(sprintf("%s|58",isodecoders),sprintf("%s|Charged",isodecoders),sprintf("%s|RPM",isodecoders))

### CHARGING ###
# Charging differences between conditions - Total
subdataset = dataset[grepl("total",dataset$cond),]
chrg_total = cast(subdataset, reference ~ sample, mean, value = 'charging')

# Charging differences between conditions - Poly
subdataset = dataset[grepl("poly",dataset$cond),]
chrg_poly = cast(subdataset, reference ~ sample, mean, value = 'charging')

### MODIFICATIONS ###
# Load modification data - Total
zz=gzfile("data/Human_stress_polysome/mismatchTable.csv.gz",'rt')   
mod_dataset = read.delim(zz)
mod_dataset = mod_dataset[grepl("total",mod_dataset$cond),]
mod_dataset = cast(mod_dataset, isodecoder*canon_pos*bam*condition ~ type, sum, value = 'proportion')
mod_dataset$mismatches = rowSums(mod_dataset[,c("A","T","G","C")],na.rm=T)
mod_dataset = mod_dataset[apply(as.data.frame(mod_dataset),1,function(x) paste0(x[c("isodecoder","canon_pos")],collapse="_")) %in% modpos,]
mod_total = cast(mod_dataset, isodecoder*canon_pos ~ bam, mean, value = 'mismatches')
# Load modification data - Poly
zz=gzfile("data/Human_stress_polysome/mismatchTable.csv.gz",'rt')   
mod_dataset = read.delim(zz)
mod_dataset = mod_dataset[grepl("poly",mod_dataset$cond),]
mod_dataset = cast(mod_dataset, isodecoder*canon_pos*bam*condition ~ type, sum, value = 'proportion')
mod_dataset$mismatches = rowSums(mod_dataset[,c("A","T","G","C")],na.rm=T)
mod_dataset = mod_dataset[apply(as.data.frame(mod_dataset),1,function(x) paste0(x[c("isodecoder","canon_pos")],collapse="_")) %in% modpos,]
mod_poly = cast(mod_dataset, isodecoder*canon_pos ~ bam, mean, value = 'mismatches')

### RPM ###
rpmdf = read.delim("data/Human_stress_polysome/Isodecoder_counts_DESEqNormalized.csv",sep="\t", row.names = 1)
# Total
rpmdf_total = rpmdf[,grepl("CW1_S1",colnames(rpmdf))]
# Poly
rpmdf_poly = rpmdf[,grepl("CW2_S2",colnames(rpmdf))]

### ODDS ###
# Load OR data
sigdataset = read.delim("data/Human_stress_polysome/interdependences_significant.tsv")
# Subset wt 58-Charg data
wt_data = sigdataset[(sigdataset$cond=="control_total")&(sigdataset$X.samples>0),]
wt_data$pairs = apply(wt_data,1,function(x) paste0(x[c("canon_var1","canon_var2")],collapse="-"))
wt_data = wt_data[wt_data$pairs=="58-Charged",]

# Create mergedf
conds = unique(mod_dataset[,c("bam","condition")])
mergedf = data.frame(row.names = 1:(length(allpos)*(ncol(chrg_total)-1)))
mergedf$reference = rep(allpos,(ncol(chrg_total)-1))
mergedf$pos = rep(sapply(allpos, function(x) strsplit(x,"\\|")[[1]][2]),(ncol(chrg_total)-1))
samplenames = as.character(sapply(colnames(chrg_total[,2:ncol(chrg_total)]), function(x) paste0(strsplit(x,"_|\\.")[[1]][5:7],collapse="_")))
mergedf$sample = c(sapply(samplenames,function(x) rep(x,length(allpos))))
for (idx in rownames(mergedf)){
  p = mergedf[idx,"reference"]
  s = mergedf[idx,"sample"]
  c = strsplit(as.character(conds[grepl(s,conds$bam),"condition"]),"_")[[1]][1]
  ref = strsplit(p,"\\|")[[1]][1]
  splitted = strsplit(ref,"-")[[1]]
  mergedf[idx,"cond"] = c
  mergedf[idx,"aa"] = splitted[length(splitted)-2]
  mergedf[idx,"wobble"] = if(substr(splitted[length(splitted)-1],1,1) %in% c("A","T")){"AT"}else{"GC"}
  ct = strsplit(p,"\\|")[[1]][2]
  if ((ct=="Charged")&(ref %in% chrg_poly$reference)&(ref %in% chrg_total$reference)){
    mergedf[idx,"diff"] = as.numeric(chrg_poly[(chrg_poly$reference==ref),grepl(s,colnames(chrg_poly))]) -
      as.numeric(chrg_total[(chrg_total$reference==ref),grepl(s,colnames(chrg_total))])
    mergedf[idx,"poly"] = as.numeric(chrg_poly[(chrg_poly$reference==ref),grepl(s,colnames(chrg_poly))])
    mergedf[idx,"total"] = as.numeric(chrg_total[(chrg_total$reference==ref),grepl(s,colnames(chrg_total))])
  }else if ((ct=="RPM")&(ref %in% rownames(rpmdf_poly))&(ref %in% rownames(rpmdf_total))){
    mergedf[idx,"diff"] = log2(as.numeric(rpmdf_poly[ref,grepl(s,colnames(rpmdf_poly))]) /
      as.numeric(rpmdf_total[ref,grepl(s,colnames(rpmdf_total))]))
    mergedf[idx,"poly"] = as.numeric(rpmdf_poly[ref,grepl(s,colnames(rpmdf_poly))])
    mergedf[idx,"total"] = as.numeric(rpmdf_total[ref,grepl(s,colnames(rpmdf_total))])
  }else if ((ref %in% mod_poly$isodecoder)&(ref %in% mod_total$isodecoder)){
    mergedf[idx,"diff"] = as.numeric(mod_poly[(mod_poly$canon_pos==ct)&(mod_poly$isodecoder==ref),grepl(s,colnames(mod_poly))]) -
      as.numeric(mod_total[(mod_total$canon_pos==ct)&(mod_total$isodecoder==ref),grepl(s,colnames(mod_total))])
    mergedf[idx,"poly"] = as.numeric(mod_poly[(mod_poly$canon_pos==ct)&(mod_poly$isodecoder==ref),grepl(s,colnames(mod_poly))])
    mergedf[idx,"total"] = as.numeric(mod_total[(mod_total$canon_pos==ct)&(mod_total$isodecoder==ref),grepl(s,colnames(mod_total))])
  }else{
    mergedf[idx,"diff"] = NA
    mergedf[idx,"poly"] = NA
    mergedf[idx,"total"] = NA
  }
  mergedf[idx,"odds"] = if(ref %in% wt_data$reference){wt_data[wt_data$reference==ref,"avg_odds"]}else{NA}
}
mergedf$reference = substr(mergedf$reference,14,100)
mergedf$reference = sapply(mergedf$reference,function(x) strsplit(x,"\\|")[[1]][1])

### Plot ###
ggplot(mergedf, aes(x=cond, y=diff, fill=wobble)) + 
  facet_wrap( ~ pos*aa, ncol=6, scales = "free") +
  geom_boxplot(lwd=0.5, notch=F, outlier.shape = NA , na.rm=T,alpha=0.3, 
               position = position_dodge(width=0.9)) +
  scale_fill_manual(values=c("#ed8975ff", "#8fb9aaff")) +
  scale_color_manual(values=c("#ed8975ff", "#8fb9aaff")) +
  geom_jitter(na.rm=T , alpha=1, size=1, aes(color = wobble),position=position_jitterdodge(dodge.width=0.9)) +
  stat_compare_means(method = "t.test", label = "p.signif") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1,color="black"),
        axis.text.y=element_text(colour="black"))
ggsave("plots/human_stress_totalVSpoly_wobble_AAfamilies.pdf",width=12,height=27)

# Plot scatter
ggplot(mergedf[(mergedf$pos=="58"),], aes(x=odds,y=diff,label=reference)) +
  facet_wrap( ~ cond, ncol=2, scales = "free") +
  geom_point(shape = 16, stroke=0.5, size=2, alpha=0.5, aes(color=wobble)) +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_text_repel() +
  theme_classic() +
  labs(x ="log2(ORwt)", y = "Delta 58 (poly-total)")
ggplot(mergedf[(mergedf$pos=="58")&(!is.na(mergedf$odds)),], aes(x=cond,y=diff,fill=odds>0)) +
  geom_violin() +
  stat_summary(fun=median, geom="point", shape=23, size=2, position = position_dodge(width=0.9)) +
  stat_compare_means(method = "wilcox.test", label = "p.signif") +
  scale_fill_manual(values=c("#304d63ff", "#f2d096ff")) +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  theme_classic() +
  labs(y = "Delta 58 (poly-total)")

ggsave("plots/human_stress_totalVSpolyMod_vs_oddsWT.pdf",width=8,height=3)

### Analyze ALKBH1 targets ###
# Based on previous paper
known58 = c("Ala-CGC","Ala-TGC","Asn-GTT","Cys-GCA","Gln-CTG","Glu-CTC","Glu-TTC","Gly-GCC","His-GTG","Leu-AAG","Lys-TTT","Tyr-GTA","Val-AAC","Val-CAC")
known58not = c("Phe-GAA","Val-TAC")
mergedf$alkbh1 = sapply(mergedf$reference,function(x) if("mito"==substr(x,1,4)){NA}else if(paste0(strsplit(x,"-")[[1]][2:3],collapse="-") %in% known58){"target"}else if(paste0(strsplit(x,"-")[[1]][2:3],collapse="-") %in% known58not){"nontarget"}else{NA})
# Plot
statdiff = compare_means(diff ~ alkbh1, data = mergedf[(mergedf$pos=="58"),], group.by = c("cond"),
                         method = "wilcox.test", ref.group = "target")
ggplot(mergedf[(mergedf$pos=="58")&(!is.na(mergedf$alkbh1)),], aes(x=cond,y=diff,fill=alkbh1)) +
  geom_violin() +
  stat_summary(fun=median, geom="point", shape=23, size=2, position = position_dodge(width=0.9)) +
  scale_fill_manual(values=c("#ed8975ff", "#8fb9aaff","grey")) +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  theme_classic() +
  labs(y = "Delta 58 (poly-total)")
# COmparison to polysome data from the old paper
knownup = c("Gly-GCC","His-GTG","Val-AAC","Val-CAC")
mergedf$alkbh1 = sapply(mergedf$reference,function(x) if(paste0(strsplit(x,"-")[[1]][2:3],collapse="-") %in% knownup){"knownup"}else{NA})
# Plot
statdiff = compare_means(diff ~ alkbh1, data = mergedf[(mergedf$pos=="58"),], group.by = c("cond"),
                         method = "wilcox.test")
ggplot(mergedf[(mergedf$pos=="58")&(!is.na(mergedf$alkbh1)),], aes(x=cond,y=diff,fill=alkbh1)) +
  geom_violin() +
  stat_summary(fun=median, geom="point", shape=23, size=2, position = position_dodge(width=0.9)) +
  scale_fill_manual(values=c("#ed8975ff", "#8fb9aaff","grey")) +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  theme_classic() +
  labs(y = "Delta 58 (poly-total)")
