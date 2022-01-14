library(ggplot2)
library(reshape)
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
  nonstressdf[idx,"cond"] = c
  ref = strsplit(p,"\\|")[[1]][1]
  pair = strsplit(strsplit(p,"\\|")[[1]][2],"-")[[1]]
  ct = pair[pair!="58"]
  if (ct=="Charged"){
    nonstressdf[idx,"CT"] = as.numeric(chrg_data[(chrg_data$reference==ref),s])
  }else{
    nonstressdf[idx,"CT"] = as.numeric(mod_data[(mod_data$canon_pos==ct)&(mod_data$isodecoder==ref),s])
  }
  nonstressdf[idx,"Mod58"] = as.numeric(mod_data[(mod_data$canon_pos=="58")&(mod_data$isodecoder==ref),s])
}
nonstressdf$reference = substr(nonstressdf$reference,14,100)

# Keep only top isodecoders
topiso = c("tRNA-Glu-TTC-4|58-Charged","tRNA-Glu-TTC-4|34-58",
           "tRNA-Glu-TTC-1|58-Charged","tRNA-Glu-TTC-1|34-58","tRNA-Glu-TTC-1|9-58",
           "mito_tRNA-Ser-TGA-1|58-Charged","mito_tRNA-Ser-TGA-1|32-58","mito_tRNA-Ser-TGA-1|37-58",
           "mito_tRNA-Lys-TTT-1|9-58")
topdata = nonstressdf[nonstressdf$reference %in% topiso,]
topdata$reference = sapply(topdata$reference,function(x) strsplit(x,"\\|")[[1]][1])
topdata$pair = factor(as.character(topdata$pair),levels = c("9-58","32-58","34-58","37-58","58-Charged"))
### Plot ###

ggplot(topdata, aes(x=pair, y=CT)) + 
  facet_wrap( ~ reference, ncol=4, scales = "free") +
  geom_point(shape = 16, stroke=1, size=2, fill = "black") +
  geom_hline(yintercept=0) +
  stat_summary(fun.data=mean_sdl, geom="pointrange", color="red", alpha=0.5) +
  theme_classic() +
  labs(x ="Position", y = "delta%")
ggsave("plots/alkbh1_kdstress_topiso.pdf",width=6,height=2)

# Keep only charging
topiso = c("tRNA-Glu-TTC-4|58-Charged",
           "tRNA-Glu-TTC-1|58-Charged",
           "tRNA-Glu-CTC-1|58-Charged",
           "tRNA-Glu-TTC-2|58-Charged")
topdata = nonstressdf[nonstressdf$reference %in% topiso,]
topdata$reference = sapply(topdata$reference,function(x) strsplit(x,"\\|")[[1]][1])
### Plot ###

ggplot(topdata, aes(x=reference, y=CT)) + 
  geom_point(shape = 16, stroke=1, size=2, fill = "black") +
  geom_hline(yintercept=0) +
  stat_summary(fun.data=mean_sdl, geom="pointrange", color="red", alpha=0.5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) +
  labs(x ="Position", y = "delta%")
ggsave("plots/alkbh1_kdstress_GluCharging.pdf",width=2,height=2)
