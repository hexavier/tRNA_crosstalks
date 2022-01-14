library(ggplot2)
library(patchwork)
library(reshape)
library(ggrepel)

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


##### WITHOUT STRESS #####
### CHARGING ###
# Charging differences between conditions
subdataset = dataset[grepl("minus",dataset$cond),]
chrg_data = cast(subdataset, reference ~ sample, mean, value = 'charging')
wt_samples = unique(subdataset[subdataset$cond=="minus-AsO2_control-shRNA","sample"])
avg_wt = rowMeans(chrg_data[,wt_samples],na.rm=T)
chrg_data = chrg_data[,!(colnames(chrg_data) %in% wt_samples)]
chrg_data[,!(colnames(chrg_data) %in% "reference")] = apply(chrg_data[,!(colnames(chrg_data) %in% "reference")],2,
                                                            function(x) x-avg_wt)

### ODDS RATIOS ###
# Load OR data
sigdataset = read.delim("data/ALKBH1_KD_stress/interdependences_significant.tsv")
# Subset wt 58-Charg data
wt_data = sigdataset[(sigdataset$cond=="minus-AsO2_control-shRNA")&((sigdataset$canon_var1=="58")|(sigdataset$canon_var2=="58"))&(sigdataset$X.samples>0),]
# Find positions crosstalking with m1A58
pos58 = unique(c(apply(wt_data,1,function(x) paste0(x[c("reference","canon_var1")],collapse="_")),
                   apply(wt_data,1,function(x) paste0(x[c("reference","canon_var2")],collapse="_"))))
pairs58 = apply(wt_data,1,function(x) paste0(c(x["reference"],paste0(x[c("canon_var1","canon_var2")],collapse="-")),collapse="|"))

### MODIFICATIONS ###
# Load modification data
zz=gzfile("data/ALKBH1_KD_stress/mismatchTable.csv.gz",'rt')   
mod_dataset = read.delim(zz)
mod_dataset = mod_dataset[grepl("minus",mod_dataset$cond),]
mod_dataset = cast(mod_dataset, isodecoder*canon_pos*bam*condition ~ type, sum, value = 'proportion')
mod_dataset$mismatches = rowSums(mod_dataset[,c("A","T","G","C")],na.rm=T)
mod_dataset = mod_dataset[apply(as.data.frame(mod_dataset),1,function(x) paste0(x[c("isodecoder","canon_pos")],collapse="_")) %in% pos58,]
# Reformat dataset
mod_data = cast(mod_dataset, isodecoder*canon_pos ~ bam, mean, value = 'mismatches')
wt_samples = unique(mod_dataset[mod_dataset$condition=="minus-AsO2_control-shRNA","bam"])
avg_wt = rowMeans(mod_data[,wt_samples],na.rm=T)
mod_data = mod_data[,!(colnames(mod_data) %in% wt_samples)]
mod_data[,!(colnames(mod_data) %in% c("isodecoder","canon_pos"))] = apply(mod_data[,!(colnames(mod_data) %in% c("isodecoder","canon_pos"))],2,
                                                                          function(x) x-avg_wt)

# Create heatdf
conds = unique(mod_dataset[,c("bam","condition")])
nonstressdf = data.frame(row.names = 1:(length(pairs58)*(ncol(chrg_data)-1)))
nonstressdf$reference = rep(pairs58,(ncol(chrg_data)-1))
nonstressdf$pair = rep(sapply(pairs58, function(x) strsplit(x,"\\|")[[1]][2]),(ncol(chrg_data)-1))
nonstressdf$sample = c(sapply(colnames(chrg_data[,2:ncol(chrg_data)]),function(x) rep(x,length(pairs58))))
nonstressdf$oddsWT = rep(wt_data$avg_odds,(ncol(chrg_data)-1))
for (idx in rownames(nonstressdf)){
  p = nonstressdf[idx,"reference"]
  s = nonstressdf[idx,"sample"]
  c = conds[conds$bam==s,"condition"]
  ref = strsplit(p,"\\|")[[1]][1]
  pair = strsplit(strsplit(p,"\\|")[[1]][2],"-")[[1]]
  ct = pair[pair!="58"]
  if (ct=="Charged"){
    nonstressdf[idx,sprintf("CT_%s",c)] = as.numeric(chrg_data[(chrg_data$reference==ref),s])
  }else{
    nonstressdf[idx,sprintf("CT_%s",c)] = as.numeric(mod_data[(mod_data$canon_pos==ct)&(mod_data$isodecoder==ref),s])
  }
  nonstressdf[idx,sprintf("Mod58_%s",c)] = as.numeric(mod_data[(mod_data$canon_pos=="58")&(mod_data$isodecoder==ref),s])
}
nonstressdf$reference = substr(nonstressdf$reference,14,100)
nonstressdf = nonstressdf[rowSums(is.na(nonstressdf))==0,]

##### WITH STRESS #####
### CHARGING ###
# Charging differences between conditions
subdataset = dataset[grepl("plus",dataset$cond),]
chrg_data = cast(subdataset, reference ~ sample, mean, value = 'charging')
wt_samples = unique(subdataset[subdataset$cond=="plus-AsO2_control-shRNA","sample"])
avg_wt = rowMeans(chrg_data[,wt_samples],na.rm=T)
chrg_data = chrg_data[,!(colnames(chrg_data) %in% wt_samples)]
chrg_data[,!(colnames(chrg_data) %in% "reference")] = apply(chrg_data[,!(colnames(chrg_data) %in% "reference")],2,
                                                            function(x) x-avg_wt)

### ODDS RATIOS ###
# Load OR data
sigdataset = read.delim("data/ALKBH1_KD_stress/interdependences_significant.tsv")
# Subset wt 58-Charg data
wt_data = sigdataset[(sigdataset$cond=="plus-AsO2_control-shRNA")&((sigdataset$canon_var1=="58")|(sigdataset$canon_var2=="58"))&(sigdataset$X.samples>0),]
# Find positions crosstalking with m1A58
pos58 = unique(c(apply(wt_data,1,function(x) paste0(x[c("reference","canon_var1")],collapse="_")),
                 apply(wt_data,1,function(x) paste0(x[c("reference","canon_var2")],collapse="_"))))
pairs58 = apply(wt_data,1,function(x) paste0(c(x["reference"],paste0(x[c("canon_var1","canon_var2")],collapse="-")),collapse="|"))

### MODIFICATIONS ###
# Load modification data
zz=gzfile("data/ALKBH1_KD_stress/mismatchTable.csv.gz",'rt')   
mod_dataset = read.delim(zz)
mod_dataset = mod_dataset[grepl("plus",mod_dataset$cond),]
mod_dataset = cast(mod_dataset, isodecoder*canon_pos*bam*condition ~ type, sum, value = 'proportion')
mod_dataset$mismatches = rowSums(mod_dataset[,c("A","T","G","C")],na.rm=T)
mod_dataset = mod_dataset[apply(as.data.frame(mod_dataset),1,function(x) paste0(x[c("isodecoder","canon_pos")],collapse="_")) %in% pos58,]
# Reformat dataset
mod_data = cast(mod_dataset, isodecoder*canon_pos ~ bam, mean, value = 'mismatches')
wt_samples = unique(mod_dataset[mod_dataset$condition=="plus-AsO2_control-shRNA","bam"])
avg_wt = rowMeans(mod_data[,wt_samples],na.rm=T)
mod_data = mod_data[,!(colnames(mod_data) %in% wt_samples)]
mod_data[,!(colnames(mod_data) %in% c("isodecoder","canon_pos"))] = apply(mod_data[,!(colnames(mod_data) %in% c("isodecoder","canon_pos"))],2,
                                                                          function(x) x-avg_wt)

# Create heatdf
conds = unique(mod_dataset[,c("bam","condition")])
stressdf = data.frame(row.names = 1:(length(pairs58)*(ncol(chrg_data)-1)))
stressdf$reference = rep(pairs58,(ncol(chrg_data)-1))
stressdf$pair = rep(sapply(pairs58, function(x) strsplit(x,"\\|")[[1]][2]),(ncol(chrg_data)-1))
stressdf$sample = c(sapply(colnames(chrg_data[,2:ncol(chrg_data)]),function(x) rep(x,length(pairs58))))
stressdf$oddsWT = rep(wt_data$avg_odds,(ncol(chrg_data)-1))
for (idx in rownames(stressdf)){
  p = stressdf[idx,"reference"]
  s = stressdf[idx,"sample"]
  c = conds[conds$bam==s,"condition"]
  ref = strsplit(p,"\\|")[[1]][1]
  pair = strsplit(strsplit(p,"\\|")[[1]][2],"-")[[1]]
  ct = pair[pair!="58"]
  if (ct=="Charged"){
    stressdf[idx,sprintf("CT_%s",c)] = as.numeric(chrg_data[(chrg_data$reference==ref),s])
  }else{
    stressdf[idx,sprintf("CT_%s",c)] = as.numeric(mod_data[(mod_data$canon_pos==ct)&(mod_data$isodecoder==ref),s])
  }
  stressdf[idx,sprintf("Mod58_%s",c)] = as.numeric(mod_data[(mod_data$canon_pos=="58")&(mod_data$isodecoder==ref),s])
}
stressdf$reference = substr(stressdf$reference,14,100)
stressdf = stressdf[rowSums(is.na(stressdf))==0,]

##### PLOT #####
# Minimum % difference in modification 58 to consider ther pair
thres = 0.03
# Plot isodecoders > thres% change m1A58
pval1 = binom.test(sum((nonstressdf[nonstressdf$`Mod58_minus-AsO2_shRNAW`>thres,"oddsWT"]>0)==(nonstressdf[nonstressdf$`Mod58_minus-AsO2_shRNAW`>thres,"CT_minus-AsO2_shRNAW"]>0)),
                   nrow(nonstressdf[nonstressdf$`Mod58_minus-AsO2_shRNAW`>thres,]), alternative = "greater")
p1 = ggplot(nonstressdf[nonstressdf$`Mod58_minus-AsO2_shRNAW`>thres,], aes(x=oddsWT,y=`CT_minus-AsO2_shRNAW`,label=reference)) +
  geom_point(shape = 21, stroke=1, size=3, aes(fill = `Mod58_minus-AsO2_shRNAW`, color=(oddsWT>0)==(`CT_minus-AsO2_shRNAW`>0))) +
  scale_fill_gradient2() +
  scale_color_manual(values=c("white","black")) +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_text_repel() +
  theme_classic() +
  labs(title ="No stress", subtitle = sprintf("Binomial Test Pval = %f",pval1$p.value), x ="log2(OR_WT)", y = "Delta_CT")

pval2 = binom.test(sum((stressdf[stressdf$`Mod58_plus-AsO2_shRNAV`>thres,"oddsWT"]>0)==(stressdf[stressdf$`Mod58_plus-AsO2_shRNAV`>thres,"CT_plus-AsO2_shRNAV"]>0)),
                   nrow(stressdf[stressdf$`Mod58_plus-AsO2_shRNAV`>thres,]), alternative = "greater")
p2 = ggplot(stressdf[stressdf$`Mod58_plus-AsO2_shRNAV`>thres,], aes(x=oddsWT,y=`CT_plus-AsO2_shRNAV`,label=reference)) +
  geom_point(shape = 21, stroke=1, size=3, aes(fill = `Mod58_plus-AsO2_shRNAV`, color=(oddsWT>0)==(`CT_plus-AsO2_shRNAV`>0))) +
  scale_fill_gradient2() +
  scale_color_manual(values=c("white","black")) +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_text_repel() +
  theme_classic() +
  labs(title ="Stress", subtitle = sprintf("Binomial Test Pval = %f",pval2$p.value), x ="log2(OR_WT)", y = "Delta_CT")

p1 | p2
ggsave(sprintf("plots/alkbh1_kdstress_allpairs_oddsWT_allsamples_delta%i.pdf",thres*100),width=20,height=8)

# Plot all
p1 = ggplot(nonstressdf, aes(x=`Mod58_minus-AsO2_shRNAW`,y=`CT_minus-AsO2_shRNAW`,label=reference)) +
  geom_point(shape = 21, stroke=1, size=3, color="black", aes(fill = oddsWT)) +
  scale_fill_gradient2(limits=c(-1,1), oob = scales::squish) +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  theme_classic() +
  labs(title ="No stress", x ="Delta_mod58", y = "Delta_CT")
p2 = ggplot(stressdf, aes(x=`Mod58_plus-AsO2_shRNAV`,y=`CT_plus-AsO2_shRNAV`,label=reference)) +
  geom_point(shape = 21, stroke=1, size=3, color="black", aes(fill = oddsWT)) +
  scale_fill_gradient2(limits=c(-1,1), oob = scales::squish) +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  theme_classic() +
  labs(title ="Stress", x ="Delta_mod58", y = "Delta_CT")
p1 | p2
#ggsave("plots/alkbh1_kdstress_allpairsVS58_allsamples.pdf",width=20,height=8)

# Without labels
# Plot isodecoders > thres% change m1A58
pval1 = binom.test(sum((nonstressdf[nonstressdf$`Mod58_minus-AsO2_shRNAW`>thres,"oddsWT"]>0)==(nonstressdf[nonstressdf$`Mod58_minus-AsO2_shRNAW`>thres,"CT_minus-AsO2_shRNAW"]>0)),
                   nrow(nonstressdf[nonstressdf$`Mod58_minus-AsO2_shRNAW`>thres,]), alternative = "greater")
p1 = ggplot(nonstressdf[nonstressdf$`Mod58_minus-AsO2_shRNAW`>thres,], aes(x=oddsWT,y=`CT_minus-AsO2_shRNAW`)) +
  geom_point(shape = 21, stroke=0.5, size=2, aes(fill = `Mod58_minus-AsO2_shRNAW`, color=(oddsWT>0)==(`CT_minus-AsO2_shRNAW`>0))) +
  scale_fill_gradient2() +
  scale_color_manual(values=c("white","black")) +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  theme_classic() +
  labs(title ="No stress", subtitle = sprintf("Binomial Test Pval = %f",pval1$p.value), x ="log2(OR_WT)", y = "Delta_CT")

pval2 = binom.test(sum((stressdf[stressdf$`Mod58_plus-AsO2_shRNAV`>thres,"oddsWT"]>0)==(stressdf[stressdf$`Mod58_plus-AsO2_shRNAV`>thres,"CT_plus-AsO2_shRNAV"]>0)),
                   nrow(stressdf[stressdf$`Mod58_plus-AsO2_shRNAV`>thres,]), alternative = "greater")
p2 = ggplot(stressdf[stressdf$`Mod58_plus-AsO2_shRNAV`>thres,], aes(x=oddsWT,y=`CT_plus-AsO2_shRNAV`)) +
  geom_point(shape = 21, stroke=0.5, size=2, aes(fill = `Mod58_plus-AsO2_shRNAV`, color=(oddsWT>0)==(`CT_plus-AsO2_shRNAV`>0))) +
  scale_fill_gradient2() +
  scale_color_manual(values=c("white","black")) +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  theme_classic() +
  labs(title ="Stress", subtitle = sprintf("Binomial Test Pval = %f",pval2$p.value), x ="log2(OR_WT)", y = "Delta_CT")

p1 | p2
ggsave(sprintf("plots/alkbh1_kdstress_allpairs_oddsWT_allsamples_delta%i_nolabels.pdf",thres*100),width=15,height=4)