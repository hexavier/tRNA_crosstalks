library(ggplot2)
library(patchwork)
library(reshape)
library(ggrepel)
library(ggpubr)

### CHARGING ###
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
  tempdata$sample = s
  for (r in rownames(tempdata)){
    subdata = raw_dataset[(raw_dataset$sample==s)&(raw_dataset$gene==tempdata[r,"reference"]),]
    cca = subdata[subdata$end=="CA","count"]
    cc = subdata[subdata$end=="CC","count"]
    tempdata[r,"charging"] = cca/(cca+cc)
  }
  dataset = rbind(dataset,tempdata)
}

##### TOTAL TRNA #####
### CHARGING ###
# Charging differences between conditions
chrg_data = cast(dataset, reference ~ sample, mean, value = 'charging')
wt_samples = unique(dataset[dataset$cond=="minus-AsO2_control-shRNA","sample"])
avg_wt = rowMeans(chrg_data[,wt_samples],na.rm=T)
chrg_data = chrg_data[,!(colnames(chrg_data) %in% wt_samples)]
chrg_data[,!(colnames(chrg_data) %in% "reference")] = apply(chrg_data[,!(colnames(chrg_data) %in% "reference")],2,
                                                            function(x) x-avg_wt)

### MODIFICATIONS ###
# Pairs subset
modpos = c("Homo_sapiens_tRNA-Glu-TTC-4_58",
           "Homo_sapiens_tRNA-Glu-TTC-1_58",
           "Homo_sapiens_tRNA-Glu-CTC-1_58")
allpos = c("Homo_sapiens_tRNA-Glu-TTC-4|58",
           "Homo_sapiens_tRNA-Glu-TTC-1|58",
           "Homo_sapiens_tRNA-Glu-CTC-1|58",
           "Homo_sapiens_tRNA-Glu-TTC-4|Charged",
           "Homo_sapiens_tRNA-Glu-TTC-1|Charged",
           "Homo_sapiens_tRNA-Glu-CTC-1|Charged",
           "Homo_sapiens_tRNA-Glu-TTC-4|RPM",
           "Homo_sapiens_tRNA-Glu-TTC-1|RPM",
           "Homo_sapiens_tRNA-Glu-CTC-1|RPM")
# Load modification data
zz=gzfile("data/ALKBH1_KD_stress/mismatchTable.csv.gz",'rt')   
mod_dataset = read.delim(zz)
mod_dataset = cast(mod_dataset, isodecoder*canon_pos*bam*condition ~ type, sum, value = 'proportion')
mod_dataset$mismatches = rowSums(mod_dataset[,c("A","T","G","C")],na.rm=T)
mod_dataset = mod_dataset[apply(as.data.frame(mod_dataset),1,function(x) paste0(x[c("isodecoder","canon_pos")],collapse="_")) %in% modpos,]
# Reformat dataset
mod_data = cast(mod_dataset, isodecoder*canon_pos ~ bam, mean, value = 'mismatches')
wt_samples = unique(mod_dataset[mod_dataset$condition=="minus-AsO2_control-shRNA","bam"])
avg_wt = rowMeans(mod_data[,wt_samples],na.rm=T)
mod_data = mod_data[,!(colnames(mod_data) %in% wt_samples)]
mod_data[,!(colnames(mod_data) %in% c("isodecoder","canon_pos"))] = apply(mod_data[,!(colnames(mod_data) %in% c("isodecoder","canon_pos"))],2,
                                                                          function(x) x-avg_wt)

### RPM ###
rpmdf = read.delim("data/ALKBH1_KD_stress/Isodecoder_counts_DESEqNormalized.csv",sep="\t", row.names = 1)
rpmdf = rpmdf[,-c(ncol(rpmdf))]
wt_samples = c("TP.CK.16S.CW05_S8_L1_bc1_TTCC.FINAL.unpaired_uniq","TP.CK.16S.CW05_S8_L1_bc2_TCTG.FINAL.unpaired_uniq")
avg_wt = rowMeans(rpmdf[,wt_samples],na.rm=T)
rpmdf = rpmdf[,!(colnames(rpmdf) %in% wt_samples)]
rpmdf = apply(rpmdf,2,function(x) log2(x/avg_wt))

# Create totaldf
conds = unique(mod_dataset[,c("bam","condition")])
totaldf = data.frame(row.names = 1:(length(allpos)*(ncol(chrg_data)-1)))
totaldf$reference = rep(allpos,(ncol(chrg_data)-1))
totaldf$pos = rep(sapply(allpos, function(x) strsplit(x,"\\|")[[1]][2]),(ncol(chrg_data)-1))
totaldf$sample = c(sapply(colnames(chrg_data[,2:ncol(chrg_data)]),function(x) rep(x,length(allpos))))
for (idx in rownames(totaldf)){
  p = totaldf[idx,"reference"]
  s = totaldf[idx,"sample"]
  c = conds[conds$bam==s,"condition"]
  ref = strsplit(p,"\\|")[[1]][1]
  splitted = strsplit(ref,"-")[[1]]
  totaldf[idx,"cond"] = c
  totaldf[idx,"aa"] = splitted[length(splitted)-2]
  ct = strsplit(p,"\\|")[[1]][2]
  if (ct=="Charged"){
    totaldf[idx,"perc"] = as.numeric(chrg_data[(chrg_data$reference==ref),s])
  }else if (ct=="RPM"){
    totaldf[idx,"perc"] = as.numeric(rpmdf[ref,grepl(substr(s,36,66),colnames(rpmdf))])
  }else{
    totaldf[idx,"perc"] = as.numeric(mod_data[(mod_data$canon_pos==ct)&(mod_data$isodecoder==ref),s])
  }
}
totaldf$reference = substr(totaldf$reference,14,100)
totaldf$reference = sapply(totaldf$reference,function(x) strsplit(x,"\\|")[[1]][1])

### Plot ###
diffcond = compare_means(perc ~ cond, data = totaldf, group.by = c("pos","reference"),
                         method = "t.test", p.adjust.method="fdr", ref.group = "plus-AsO2_shRNAV")
ggplot(totaldf, aes(x=pos, y=perc, color=cond, group=cond)) + 
  facet_wrap( ~ reference*pos, ncol=3, scales = "free") +
  geom_point(shape=16, na.rm=T, size=2, position=position_dodge(width=0.8)) +
  scale_color_manual(values=c("#ed8975ff", "#8fb9aaff","#304D63")) +
  geom_hline(yintercept=0) +
  theme_classic() +
  theme(axis.text.x = element_text(color="black"),
        axis.text.y=element_text(colour="black"))
ggsave("plots/alkbh1_kdstress_Glu.pdf",width=7,height=5)

