library(ggpubr)
library(ggplot2)
library(ggrepel)

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
  outtemp$rpm_total = (totaltemp[commongenes,"count"]/sum(totaltemp[commongenes,"count"]))*1000000
  outtemp$rpm_poly = (polytemp[commongenes,"count"]/sum(polytemp[commongenes,"count"]))*1000000
  outtemp$TE = log2(outtemp$rpm_poly/outtemp$rpm_total)
  outtemp$cond = filemapping[f,"V3"]
  TEdata = rbind(TEdata,outtemp)
}
# Exclude mito which have weird genetic code
TEdata = TEdata[!grepl("GRCh38:MT",TEdata$loc),]

# Find top TE genes
wtdata = TEdata[TEdata$cond=="control",]
avgTE = sapply(unique(wtdata$gene),function(x) mean(wtdata[wtdata$gene %in% x,"TE"]))

# Load RCU
codons = read.csv("data/codons_table.tab", sep="\t", row.names = 1)
codus = read.table("data/Human_stress_polysome_mRNA/human_combine_codon.tsv", header = T, row.names = 1)
codus = codus[rownames(codus) %in% wtdata$gene,rownames(codons)[codons$AA!="Stop"]]
colnames(codus) = paste0(codons[colnames(codus),"AA"],colnames(codus))
relcodus = t(apply(codus,1,rscu,colnames(codus))); colnames(relcodus) = colnames(codus)

# Compute CU
weighCU = t(scale(avgTE)) %*% relcodus[names(avgTE),]
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
