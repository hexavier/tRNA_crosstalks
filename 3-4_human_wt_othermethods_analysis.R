library(seqinr)
library(ggplot2)
library(reshape)
library(GGally)
library(ggrepel)
library(gplots)
library(RColorBrewer)

# QuantM: HEK293 T-Rex Flp-IN, Superscript IV at 55?C for 60 min
# mimseq: HEK293T, TGIRT 16h
# msrseq: HEK293T, overnight SSIV
# DM-tRNA-seq: HEK293T, TGIRT at 60?C for 60 min

# Load data
alldataset = read.delim("data/hek293/othermethods/interdependences.tsv")
sigdataset = read.delim("data/hek293/othermethods/interdependences_significant.tsv")
sigct = read.delim("results/crosstalks_modifications_overlap.tsv", sep="\t", row.names = 1)
sigct$n = rowSums(sigct)
sigct$name = sapply(rownames(sigct),function(x) strsplit(x,"\\|")[[1]][1])
sigct$mod1 = sapply(rownames(sigct),function(x) strsplit(strsplit(x,"\\|")[[1]][2],"-")[[1]][1])
sigct$mod2 = sapply(rownames(sigct),function(x) strsplit(strsplit(x,"\\|")[[1]][2],"-")[[1]][2])

### Are consensus CT having more extreme OR?
# Dataset of consesus with OR
dataset = c()
for (iso in rownames(sigct)){
  tempdf = data.frame(row.names = 1:4)
  tempdf$isodecoder = iso
  n=1
  for (meth in c("msrseq", "quantm", "mimseq", "dmseq-untreated")){
    tempdf[n, "number"] = sigct[iso,"n"]
    tempdf[n, "method"] = strsplit(meth,"-")[[1]][1]
    if (sigct[iso,strsplit(meth,"-")[[1]][1]]){
      tempdf[n,"OR"] = sigdataset[(sigdataset$reference %in% sigct[iso,"name"])&(sigdataset$canon_var1 %in% sigct[iso,"mod1"])&(sigdataset$canon_var2 %in% sigct[iso,"mod2"])&(sigdataset$cond %in% meth),"avg_odds"]
    }else{
      tempdf[n,"OR"] = NA
    }
    n=n+1
  }
  dataset = rbind(dataset,tempdf)
}

# Plot absolute OR
dataset$method = factor(as.character(dataset$method), levels=c("msrseq","quantm","mimseq","dmseq"))
ggplot(dataset, aes(x=number, y=abs(OR), color = method, group = number)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.3, alpha=0.5) +
  scale_color_manual(values=c("#F2D096","#8FB9AA","#ED8975","#304D63")) +
  theme(axis.text.x = element_text(color="black"),
        axis.text.y=element_text(colour="black")) +
  labs(x = "Consensus datasets", y = "Abs(OR)") +
  theme_classic()
ggsave("plots/othermethods_consensusVSor.pdf", height = 4, width = 4)

### Is depth and distance from 3' affecting methods?
# Load data
isoseq = read.fasta("data/hek293/othermethods/isodecoderTranscripts.fa")
mapping = list(msrseq=c("TP-CK-4S-CW1_S1_L1_bc1_TTCC","TP-CK-4S-CW1_S1_L1_bc2_TCTG","TP-CK-4S-CW1_S1_L1_bc3_TGGT"),
            quantm=c("SRR10587152","SRR10587153"),
            mimseq=c("SRR12026344","SRR12026345"),
            dmseq=c("SRR5398247","SRR5398248","SRR5398249"))

recover_depth <- function(textstr){
  splitted = strsplit(textstr, '\n')[[1]]
  cov = sum(sapply(splitted[2:5], function(x) as.numeric(strsplit(x," ")[[1]][length(strsplit(x," ")[[1]])])))
  return(cov)
}


# Dataset with all this raw data
dataset = c()
for (iso in rownames(sigct)){
  tempdf = data.frame(row.names = 1:10)
  tempdf$isodecoder = iso
  n=1
  for (meth in c("msrseq", "quantm", "mimseq", "dmseq")){
    subdf = alldataset[(alldataset$ref %in% sigct[iso,"name"])&(alldataset$canon_var1 %in% sigct[iso,"mod1"])&(alldataset$canon_var2 %in% sigct[iso,"mod2"])&(alldataset$sample %in% mapping[[meth]]),]
    for (s in rownames(subdf)){
      tempdf[n, "number"] = sigct[iso,"n"]
      tempdf[n, "method"] = meth
      tempdf[n, "moddist"] = length(isoseq[[subdf[1,"ref"]]]) - as.numeric(subdf[1,"var1"])
      if (sigct[iso,strsplit(meth,"-")[[1]][1]]){
        tempdf[n,"depth"] = recover_depth(subdf[s,"values"])
      }else{
        tempdf[n,"depth"] = NA
      }
      n=n+1
    }
  }
  dataset = rbind(dataset,tempdf)
}

# Plot number of reads vs moddist
dataset$method = factor(as.character(dataset$method), levels=c("msrseq","quantm","mimseq","dmseq"))
ggplot(dataset, aes(x=moddist, y=depth, color=method, alpha=number))+
  geom_point(na.rm = T) +
  scale_color_manual(values=c("#F2D096","#8FB9AA","#ED8975","#304D63")) +
  theme(axis.text.x = element_text(color="black"),
        axis.text.y=element_text(colour="black")) +
  labs(x = "Distance from 3'CCA", y = "Reads depth") +
  theme_classic()

### Sequencing depth
covdf = read.delim("data/hek293/othermethods/readme_files_cov.txt", header = F)
covdf = covdf[!(covdf$V2 %in% "dmseq-alkb"),]
covdf$V2 = factor(as.character(covdf$V2), levels=c("msrseq","quantm","mimseq","dmseq-untreated"))
ggplot(covdf,aes(x=V1, y=V3, fill=V2)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("#F2D096","#8FB9AA","#ED8975","#304D63")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color="black"),
        axis.text.y=element_text(colour="black")) +
  labs(title="Method read coverage", x = NULL, y = "# Uniquely mapped reads")
ggsave("plots/othermethods_cov.pdf", height = 8, width = 6)

### Check % modifications between methods
# Load data
zz=gzfile("data/hek293/othermethods/mismatchTable.csv.gz",'rt')
mod_dataset = read.delim(zz)
mod_dataset = cast(mod_dataset, isodecoder*canon_pos*bam*condition ~ type, sum, value = 'proportion')
mod_dataset$mismatches = rowSums(mod_dataset[,c("A","T","G","C")],na.rm=T)
methodsmod = cast(mod_dataset, isodecoder*canon_pos ~ condition, mean, value = 'mismatches')
methodsmod$label = sprintf("%s|%s",substr(methodsmod$isodecoder,14,100),methodsmod$canon_pos)
methodsmod = methodsmod[rowSums(is.na(methodsmod))==0,]
methodsmod = methodsmod[,c(-3)]
colnames(methodsmod) = sapply(colnames(methodsmod), function(x) strsplit(x,"-")[[1]][1])

# Plot
cleandf = methodsmod[rowSums(methodsmod[3:6])>0.05,]
ggpairs(cleandf[3:6], aes(alpha=0.2))+
  theme(axis.text.x = element_text(color="black"),
        axis.text.y=element_text(colour="black")) +
  labs(x = "% modifications", y = "% modifications") +
  theme_classic()
ggsave("plots/othermethods_mods_pairplot.pdf", height = 7, width = 8)

# Differences between dmseq and msrseq
methodsmod$delta = methodsmod$dmseq - methodsmod$msrseq
methodsmod$label = factor(x = as.character(methodsmod$label), levels = as.character(unique(methodsmod[order(-methodsmod$delta),"label"])))

ggplot(methodsmod[abs(methodsmod$delta)>0.5,],aes(x=label, y=delta)) +
  geom_bar(stat="identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color="black"),
        axis.text.y=element_text(colour="black")) +
  labs(title="DM-tRNA-seq - MSR-seq", x = NULL, y = "Delta %")

# Heatmap
heatdf = data.frame(methodsmod[,3:6],row.names = methodsmod$label)
heatdf = heatdf[apply(heatdf,1,sd)>0.2,]
Colors=brewer.pal(9,"YlOrRd")
pdf(file="plots/othermethods_mods_heatmap.pdf", width=18, height=4.5)
heatmap.2(t(heatdf), trace="none", margins = c(9,9), col=Colors,
          hclustfun = function(x) hclust(x, method="ward.D2"))
dev.off()
