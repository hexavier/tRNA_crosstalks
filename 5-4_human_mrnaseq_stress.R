library(ggpubr)
library(ggplot2)
library(reshape)
library(gplots)

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
  totaltemp = read.table(as.character(filemapping[f,"V1"]),header=T,row.names = 1)
  totaltemp = totaltemp[rownames(totaltemp) %in% rownames(genenames)[genenames$feature=="cds"],,drop=F]
  polytemp = read.table(as.character(filemapping[f,"V2"]),header=T,row.names = 1)
  polytemp = polytemp[rownames(polytemp) %in% rownames(genenames)[genenames$feature=="cds"],,drop=F]
  commongenes = rownames(totaltemp)[rownames(totaltemp) %in% rownames(polytemp)]
  # Save data
  outtemp = data.frame(row.names = 1:length(commongenes))
  outtemp$gene = commongenes
  outtemp$file = as.character(filemapping[f,"V1"])
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

### Plot dot plot ###
# Structure data
polydata = cast(TEdata, gene ~ file, sum, value = "tpm_poly")
rownames(polydata) = polydata$gene; polydata = polydata[,-c(1)]
conddata = polydata[,!(colnames(polydata) %in% as.character(filemapping[filemapping$V3=="control","V1"]))]
avgwt = rowMeans(polydata[,(colnames(polydata) %in% as.character(filemapping[filemapping$V3=="control","V1"]))],na.rm=T)
FCdata = apply(conddata,2,function(x) log2(x/avgwt))
rownames(FCdata) = rownames(conddata); colnames(FCdata) = colnames(conddata)
FCdata[is.infinite(FCdata)] <- NA
# Subset
plotdf = melt(FCdata[rownames(FCdata) %in% c("ENST00000216489.8","ENST00000203001.7","ENST00000306108.10","ENST00000389749.9"),])
plotdf$cond = sapply(as.character(plotdf$X2), function(x) as.character(filemapping[filemapping$V1==x,"V3"]))
plotdf$X1 = factor(as.character(plotdf$X1), levels = c("ENST00000216489.8","ENST00000203001.7","ENST00000306108.10","ENST00000389749.9"))
# Plot
ggplot(plotdf, aes(x=cond, y=value, color=X1)) + 
  geom_boxplot(position=position_dodge(width=0.8), color="black", aes(fill=X1)) +
  geom_point(shape=16, na.rm=T, size=3, position=position_dodge(width=0.8), alpha=0.6) +
  scale_color_manual(values=c(scales::muted("red"), scales::muted("blue"),scales::muted("blue"),scales::muted("blue"))) +
  scale_fill_manual(values=c("white","white","white","white")) +
  geom_hline(yintercept=0) +
  theme_classic() +
  theme(axis.text.x = element_text(color="black"),
        axis.text.y=element_text(colour="black")) +
  ylab("FC(TPM)")

### Plot heatmap
polydata = cast(TEdata, gene ~ cond, mean, value = "tpm_poly")
rownames(polydata) = polydata$gene; polydata = polydata[,-c(1)]
plotdf = log2(polydata[,colnames(polydata)!="control"]/polydata[,colnames(polydata)=="control"])
plotdf = plotdf[c("ENST00000216489.8","ENST00000203001.7","ENST00000306108.10","ENST00000389749.9"),]

pdf(file="plots/human_mrnaseq_TPMenzymes.pdf", width=5, height=5.5)
heatmap.2(as.matrix(plotdf), col=bluered(20),trace="none",Rowv=F,Colv=F,margins=c(12,10),
          labRow = c("ALKBH1","TRMT6","TRMT61B","TRMT61A"))
dev.off()
