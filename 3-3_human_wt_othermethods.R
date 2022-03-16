library(eulerr)
library(UpSetR)
library(reshape)

# QuantM: HEK293 T-Rex Flp-IN, Superscript IV at 55°C for 60 min
# mimseq: HEK293T, TGIRT 16h
# msrseq: HEK293T, overnight SSIV
# DM-tRNA-seq: HEK293T, TGIRT at 60◦C for 60 min

#### MODIFICATIONS ####
# Load data
zz=gzfile("data/hek293/othermethods/mismatchTable.csv.gz",'rt')
mod_dataset = read.delim(zz)
# MSR-seq
mod_dataset1 = mod_dataset[grepl("msrseq",mod_dataset$condition),]
mod_dataset1 = cast(mod_dataset1, isodecoder*canon_pos*bam*condition ~ type, sum, value = 'proportion')
mod_dataset1$mismatches = rowSums(mod_dataset1[,c("A","T","G","C")],na.rm=T)
msrseq = cast(mod_dataset1, isodecoder*canon_pos ~ condition, mean, value = 'mismatches')
msrseq = msrseq[msrseq$msrseq>=0.05,] # keep modifications detected in 5% reads
# Quantm
mod_dataset2 = mod_dataset[grepl("quantm",mod_dataset$condition),]
mod_dataset2 = cast(mod_dataset2, isodecoder*canon_pos*bam*condition ~ type, sum, value = 'proportion')
mod_dataset2$mismatches = rowSums(mod_dataset2[,c("A","T","G","C")],na.rm=T)
quantm = cast(mod_dataset2, isodecoder*canon_pos ~ condition, mean, value = 'mismatches')
quantm = quantm[quantm$quantm>=0.05,] # keep modifications detected in 5% reads
# Mimseq
mod_dataset3 = mod_dataset[grepl("mimseq",mod_dataset$condition),]
mod_dataset3 = cast(mod_dataset3, isodecoder*canon_pos*bam*condition ~ type, sum, value = 'proportion')
mod_dataset3$mismatches = rowSums(mod_dataset3[,c("A","T","G","C")],na.rm=T)
mimseq = cast(mod_dataset3, isodecoder*canon_pos ~ condition, mean, value = 'mismatches')
mimseq = mimseq[mimseq$mimseq>=0.05,] # keep modifications detected in 5% reads
# DM-tRNA-seq
mod_dataset4 = mod_dataset[grepl("dmseq-untreated",mod_dataset$condition),]
mod_dataset4 = cast(mod_dataset4, isodecoder*canon_pos*bam*condition ~ type, sum, value = 'proportion')
mod_dataset4$mismatches = rowSums(mod_dataset4[,c("A","T","G","C")],na.rm=T)
dmseq = cast(mod_dataset4, isodecoder*canon_pos ~ condition, mean, value = 'mismatches')
dmseq = dmseq[dmseq$dmseq>=0.05,] # keep modifications detected in 5% reads

# Create lists
overlap = list(msrseq = sprintf("%s|%s",msrseq$isodecoder,msrseq$canon_pos),
               quantm = sprintf("%s|%s",quantm$isodecoder,quantm$canon_pos),
               mimseq = sprintf("%s|%s",mimseq$isodecoder,mimseq$canon_pos),
               dmseq = sprintf("%s|%s",dmseq$isodecoder,dmseq$canon_pos))

# Get common mods
countdata = data.frame(m1AG9=unlist(lapply(overlap, function(x) sum(grepl("\\|9",x)))),
                       m22G26=unlist(lapply(overlap, function(x) sum(grepl("\\|26",x)))),
                       m1G37=unlist(lapply(overlap, function(x) sum(grepl("\\|37",x)))),
                       m1A58=unlist(lapply(overlap, function(x) sum(grepl("\\|58",x)))))

# Plot VENN diagrams
pdf("plots/venn_crosstalks_modifications_overlap.pdf", width=7, height=5)
## Plot diagram
fit <- euler(overlap)
print(plot(fit, fills=c("#F2D096","#8FB9AA","#ED8975","#304D63"), edges=F, legend=F, labels=T, quantities=T,main="Modifications"))
print(upset(fromList(overlap), order.by = "freq", point.size = 3, line.size = 1.2,
            mainbar.y.label = "Modifications Intersections", sets.x.label = "Detected modifications",
            text.scale = 1.3))

#### CROSSTALKS ####
# Load data
sigdataset = read.delim("data/hek293/othermethods/interdependences_significant.tsv")
msrseq = sigdataset[grepl("msrseq",sigdataset$cond),]
quantm = sigdataset[grepl("quantm",sigdataset$cond),]
mimseq = sigdataset[grepl("mimseq",sigdataset$cond),]
dmseq = sigdataset[grepl("dmseq-untreated",sigdataset$cond),]

# Keep only changes detected in all replicates
msrseq = msrseq[msrseq$X.samples>=3,]
quantm = quantm[quantm$X.samples>=2,]
mimseq = mimseq[mimseq$X.samples>=2,]
dmseq = dmseq[dmseq$X.samples>=3,]

# Discretize odds ratio
msrseq$avg_odds = msrseq$avg_odds>0
quantm$avg_odds = quantm$avg_odds>0
mimseq$avg_odds = mimseq$avg_odds>0
dmseq$avg_odds = dmseq$avg_odds>0

# Keep only crosstalk pairs for which modifications are detected in all methods
all4 = names(table(c(overlap$msrseq,overlap$quantm,overlap$mimseq,overlap$dmseq))[table(c(overlap$msrseq,overlap$quantm,overlap$mimseq,overlap$dmseq))==4])
idx = apply(msrseq,1,function(x) all(c(sprintf("%s|%s",x["reference"],x["canon_var1"]),sprintf("%s|%s",x["reference"],x["canon_var2"])) %in% all4))
msrseq = msrseq[idx,]
idx = apply(quantm,1,function(x) all(c(sprintf("%s|%s",x["reference"],x["canon_var1"]),sprintf("%s|%s",x["reference"],x["canon_var2"])) %in% all4))
quantm = quantm[idx,]
idx = apply(mimseq,1,function(x) all(c(sprintf("%s|%s",x["reference"],x["canon_var1"]),sprintf("%s|%s",x["reference"],x["canon_var2"])) %in% all4))
mimseq = mimseq[idx,]
idx = apply(dmseq,1,function(x) all(c(sprintf("%s|%s",x["reference"],x["canon_var1"]),sprintf("%s|%s",x["reference"],x["canon_var2"])) %in% all4))
dmseq = dmseq[idx,]

# Create lists
overlap = list(msrseq = sprintf("%s|%s-%s|%s",msrseq$reference,msrseq$canon_var1,msrseq$canon_var2,msrseq$avg_odds),
               quantm = sprintf("%s|%s-%s|%s",quantm$reference,quantm$canon_var1,quantm$canon_var2,quantm$avg_odds),
               mimseq = sprintf("%s|%s-%s|%s",mimseq$reference,mimseq$canon_var1,mimseq$canon_var2,mimseq$avg_odds),
               dmseq = sprintf("%s|%s-%s|%s",dmseq$reference,dmseq$canon_var1,dmseq$canon_var2,dmseq$avg_odds))

# Plot VENN diagrams
## Plot diagram
fit <- euler(overlap)
print(plot(fit, fills=c("#F2D096","#8FB9AA","#ED8975","#304D63"), edges=F, legend=F, labels=T, quantities=T, main="Crosstalks"))
print(upset(fromList(overlap), order.by = "freq", point.size = 3, line.size = 1.2,
            mainbar.y.label = "Crosstalks Intersections", sets.x.label = "Detected crosstalks",
            text.scale = 1.3))
dev.off()

# Create dataframe
outdf = data.frame(row.names = unique(c(overlap$msrseq,overlap$quantm,overlap$mimseq,overlap$dmseq)))
for (s in names(overlap)){
  outdf[,s] = rownames(outdf) %in% overlap[[s]]
}
rownames(outdf) = sapply(rownames(outdf), function(x) gsub(strsplit(x,"\\|")[[1]][3],if(strsplit(x,"\\|")[[1]][3]=="TRUE"){"OR>1"}else{"OR<1"},x))
write.table(outdf,"results/crosstalks_modifications_overlap.tsv", sep="\t", row.names = T, quote = F)

###############################################################################################
################################## PLOT HISTOGRAM #############################################
###############################################################################################

library(ggplot2)
library(patchwork)

#### MSRSEQ ####

# Load data
sigdataset = read.delim("data/hek293/othermethods/interdependences_significant.tsv")
# Subset wt data with 3 significant samples
wt_data = sigdataset[(sigdataset$cond=="msrseq")&(sigdataset$X.samples>2),]

### Histogram ###
wt_data$pair = apply(wt_data,1,function(x) paste0(x[c("canon_var1","canon_var2")],collapse="-"))
paircounts = data.frame(table(wt_data$pair))
paircounts$Var1 = factor(x = as.character(paircounts$Var1), levels = as.character(unique(paircounts[order(-paircounts$Freq),"Var1"])))
p1 = ggplot(paircounts[paircounts$Freq>1,],aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color="black"),
        axis.text.y=element_text(colour="black")) +
  labs(title="MSR-seq", x = NULL, y = "# Isodecoders")

#### MIMSEQ ####

# Subset wt data with 2 significant samples
wt_data = sigdataset[(sigdataset$cond=="mimseq")&(sigdataset$X.samples>1),]

### Histogram ###
wt_data$pair = apply(wt_data,1,function(x) paste0(x[c("canon_var1","canon_var2")],collapse="-"))
paircounts = data.frame(table(wt_data$pair))
paircounts$Var1 = factor(x = as.character(paircounts$Var1), levels = as.character(unique(paircounts[order(-paircounts$Freq),"Var1"])))
p2 = ggplot(paircounts[paircounts$Freq>1,],aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color="black"),
        axis.text.y=element_text(colour="black")) +
  labs(title="mim-tRNA-seq", x = NULL, y = "# Isodecoders")

#### QUANTM ####

# Subset wt data with 2 significant samples
wt_data = sigdataset[(sigdataset$cond=="quantm")&(sigdataset$X.samples>1),]

### Histogram ###
wt_data$pair = apply(wt_data,1,function(x) paste0(x[c("canon_var1","canon_var2")],collapse="-"))
paircounts = data.frame(table(wt_data$pair))
paircounts$Var1 = factor(x = as.character(paircounts$Var1), levels = as.character(unique(paircounts[order(-paircounts$Freq),"Var1"])))
p3 = ggplot(paircounts[paircounts$Freq>1,],aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color="black"),
        axis.text.y=element_text(colour="black")) +
  labs(title="QuantM-seq", x = NULL, y = "# Isodecoders")

#### DM-TRNA-SEQ ####

# Subset wt data with 2 significant samples
wt_data = sigdataset[(sigdataset$cond=="dmseq-untreated")&(sigdataset$X.samples>2),]

### Histogram ###
wt_data$pair = apply(wt_data,1,function(x) paste0(x[c("canon_var1","canon_var2")],collapse="-"))
paircounts = data.frame(table(wt_data$pair))
paircounts$Var1 = factor(x = as.character(paircounts$Var1), levels = as.character(unique(paircounts[order(-paircounts$Freq),"Var1"])))
p4 = ggplot(paircounts[paircounts$Freq>1,],aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color="black"),
        axis.text.y=element_text(colour="black")) +
  labs(title="DM-tRNA-seq", x = NULL, y = "# Isodecoders")

p1 / p2 / p3 / p4
ggsave("plots/histogram_human_wt_pairs_othermethods.pdf",width=7,height=6)
