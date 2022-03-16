library(ggpubr)
library(ggplot2)
library(ggrepel)
library(reshape)
library(gplots)

rscu <- function(data,codons){
  # Compute relative data
  aa = sapply(codons,function(x) substr(x,1,nchar(x)-3))
  uniqueaa = unique(aa)
  outdata = numeric(length=length(data))
  for (n in uniqueaa){
    idx = (aa %in% n)
    total = sum(data[idx],na.rm=T)
    outdata[idx] = (data[idx]*sum(idx))/total
    if (total %in% 0){
      outdata[idx] = 1.0
    }
  }
  return(outdata)
}


# Load data
genenames = read.table("data/Human_stress_polysome_mRNA/human_combine_name_dict.txt", header=T, row.names = 1)
filemapping = read.table("data/Human_stress_polysome_mRNA/mapping.txt", header=F)
# Load RCU
codons = read.csv("data/codons_table.tab", sep="\t", row.names = 1)
zz=gzfile("data/Human_stress_polysome_mRNA/human_combine_codon.tsv.gz",'rt')   
codus = read.table(zz, header = T, row.names = 1)

# Compute TE
TEdata = c()
for (f in rownames(filemapping)){
  # Load
  totaltemp = read.table(filemapping[f,"V1"],header=T,row.names = 1)
  totaltemp = totaltemp[rownames(totaltemp) %in% rownames(genenames)[genenames$feature=="cds"],,drop=F]
  polytemp = read.table(filemapping[f,"V2"],header=T,row.names = 1)
  polytemp = polytemp[rownames(polytemp) %in% rownames(genenames)[genenames$feature=="cds"],,drop=F]
  commongenes = rownames(totaltemp)[rownames(totaltemp) %in% rownames(polytemp)]
  # Save data
  outtemp = data.frame(row.names = 1:length(commongenes))
  outtemp$gene = commongenes
  outtemp$genename = genenames[commongenes,"gene_symbol"]
  outtemp$loc = genenames[commongenes,"chromosome"]
  kblen = rowSums(codus[commongenes,rownames(codons)[codons$AA!="Stop"]])*3/1000
  rpktotal = totaltemp[commongenes,"count"]/kblen
  rpkpoly = polytemp[commongenes,"count"]/kblen
  outtemp$tpm_total = (rpktotal/sum(rpktotal))*1000000
  outtemp$tpm_poly = (rpkpoly/sum(rpkpoly))*1000000
  outtemp$TE = log2(outtemp$tpm_poly/outtemp$tpm_total)
  outtemp$cond = filemapping[f,"V3"]
  TEdata = rbind(TEdata,outtemp)
}
# Exclude mito which have weird genetic code
TEdata = TEdata[!grepl("GRCh38:MT",TEdata$loc),]
# Keep only genes with at least TPM=10 in all conditions
overthres = table(TEdata[(TEdata$tpm_total>10)&(TEdata$tpm_poly>10),"gene"])
TEdata = TEdata[TEdata$gene %in% names(overthres)[overthres==12],]

# Find top TE genes
condTE = cast(TEdata, gene ~ cond, mean, value = "TE")
rownames(condTE) = condTE$gene; condTE = condTE[,-c(1)]
norm_condTE = apply(condTE,2,scale)
rownames(norm_condTE) = rownames(condTE); colnames(norm_condTE) = colnames(condTE)
norm_condTE[is.na(norm_condTE)] <-0

# Compute CU
codus = codus[rownames(codus) %in% TEdata$gene,rownames(codons)[codons$AA!="Stop"]]
colnames(codus) = paste0(codons[colnames(codus),"AA"],colnames(codus))
relcodus = t(apply(codus,1,rscu,colnames(codus))); colnames(relcodus) = colnames(codus)
weighCU = t(norm_condTE) %*% relcodus[rownames(norm_condTE),]
weighCU = t(sapply(rownames(weighCU), function(x) weighCU[x,]/sum(norm_condTE[,x]!=0)))
rownames(weighCU) = colnames(norm_condTE)

# Plot
CUdf = melt(weighCU[!(rownames(weighCU) %in% "control"),])
CUdf$control = apply(CUdf,1,function(x) weighCU["control",x["X2"]])
ggplot(CUdf, aes(x=control, y=value, label=X2)) +
  facet_wrap( ~ X1, ncol=2, scales = "free") +
  geom_point(shape = 16, size=2) +
  geom_abline(intercept = 0, slope=1) +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_text_repel() +
  theme_classic()

# Compute diffRCU
CUdf$diff = CUdf$value - CUdf$control
#CUdf$X2 = factor(as.character(CUdf$X2),levels = unique(CUdf[order(CUdf$diff,decreasing = T),"X2"]))
ggplot(CUdf, aes(x=X2, y=diff)) +
  facet_wrap( ~ X1, ncol=1, scales = "free") +
  geom_bar(stat="identity") +
  geom_hline(yintercept=0) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Heatmap
CUdf = CUdf[,!(colnames(CUdf) %in% "value")]
heatdf = as.matrix(cast(CUdf, X2 ~ X1, value="diff"))
heatdf = heatdf[!(rownames(heatdf) %in% c("MetATG","TrpTGG")),] # Keep only aa with more than 1 codon
hclust.ward <- function(x) hclust(x, method="ward.D2")
pdf(file="plots/codon_usage_stress_heatmap_tpmthres10.pdf", width=12, height=4.5)
heatmap.2(t(heatdf), col=bluered, trace="none",Rowv=F,Colv=T,margins=c(5,10),hclustfun=hclust.ward)
dev.off()

# Which AA family change the most upon stress
CUdf$aa = substr(CUdf$X2,1,3)
AAvar = melt(cast(CUdf, aa ~ X1, sd, value="diff"), id=c("aa"))
AAvar$aa = factor(as.character(AAvar$aa),levels = unique(AAvar[order(AAvar$value,decreasing = T),"aa"]))
ggplot(AAvar, aes(x=aa, y=value)) +
  facet_wrap( ~ X1, ncol=1, scales = "free") +
  geom_bar(stat="identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
