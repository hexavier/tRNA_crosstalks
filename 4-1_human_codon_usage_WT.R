library(ggpubr)
library(ggplot2)
library(ggrepel)
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

range01 <- function(x){(x-min(x, na.rm = T))/(max(x, na.rm = T)-min(x, na.rm = T))}

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
wtdata = TEdata[TEdata$cond=="control",]
wtdata = wtdata[!grepl("GRCh38:MT",wtdata$loc),]
# Keep only genes with at least TPM=10 in all conditions
overthres = table(wtdata[(wtdata$tpm_total>10)&(wtdata$tpm_poly>10),"gene"])
wtdata = wtdata[wtdata$gene %in% names(overthres)[overthres==3],]

# Find top TE genes
avgTE = sapply(unique(wtdata$gene),function(x) mean(wtdata[wtdata$gene %in% x,"TE"]))
norm_avgTE = scale(avgTE)

# Compute CU
codus = codus[rownames(codus) %in% wtdata$gene,rownames(codons)[codons$AA!="Stop"]]
colnames(codus) = paste0(codons[colnames(codus),"AA"],colnames(codus))
relcodus = t(apply(codus,1,rscu,colnames(codus))); colnames(relcodus) = colnames(codus)
weighCU = t(norm_avgTE) %*% relcodus[rownames(norm_avgTE),]
weighCU = weighCU/sum(norm_avgTE!=0)

topCU = colMeans(relcodus[names(avgTE[avgTE > 1]),])
bottomCU = colMeans(relcodus[names(avgTE[avgTE < (-1)]),])

# Plot
CUdf = data.frame(topCU,bottomCU); CUdf$codon = rownames(CUdf)
ggplot(CUdf, aes(x=topCU, y=bottomCU, label=codon)) +
  geom_point(shape = 16, size=2) +
  geom_abline(intercept = 0, slope=1) +
  geom_text_repel() +
  theme_classic()
CUdf = data.frame(codon=colnames(weighCU),CU=as.numeric(weighCU))
#CUdf$codon = factor(as.character(CUdf$codon),levels = CUdf[order(CUdf$CU,decreasing = T),"codon"])
ggplot(CUdf, aes(x=codon, y=CU)) +
  geom_bar(stat="identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# Heatmap
rownames(CUdf) = CUdf$codon; CUdf = CUdf[,-c(1), drop=F]
CUdf = CUdf[!(rownames(CUdf) %in% c("MetATG","TrpTGG")),, drop=F] # Keep only aa with more than 1 codon
hclust.ward <- function(x) hclust(x, method="ward.D2")
pdf(file="plots/codon_usage_control_heatmap_tpmthres10.pdf", width=12, height=4.5)
heatmap.2(t(cbind(CUdf,CUdf)), col=bluered, trace="none",Rowv=F,Colv=T,margins=c(5,10),hclustfun=hclust.ward)
dev.off()

# Which AA family change the most between codons
CUdf = data.frame(codon=colnames(weighCU),CU=as.numeric(weighCU))
CUdf$aa = substr(CUdf$codon,1,3)
AAvar = sapply(unique(CUdf$aa), function(x) sd(CUdf[CUdf$aa==x,"CU"]))
AAdf = data.frame(aa=names(AAvar),stdev=as.numeric(AAvar))
AAdf$aa = factor(as.character(AAdf$aa),levels = unique(AAdf[order(AAdf$stdev,decreasing = T),"aa"]))
ggplot(AAdf, aes(x=aa, y=stdev)) +
  geom_bar(stat="identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
