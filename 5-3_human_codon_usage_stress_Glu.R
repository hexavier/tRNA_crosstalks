library(ggpubr)
library(ggplot2)
library(reshape)

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
norm_condTE = apply(condTE,2,function(x) x)
rownames(norm_condTE) = rownames(condTE); colnames(norm_condTE) = colnames(condTE)

# Compute CU
codus = codus[rownames(codus) %in% TEdata$gene,rownames(codons)[codons$AA!="Stop"]]
colnames(codus) = paste0(codons[colnames(codus),"AA"],colnames(codus))
relcodus = t(apply(codus,1,rscu,colnames(codus))); colnames(relcodus) = colnames(codus)

# Build dataset
CUdf = melt(norm_condTE[,!(colnames(norm_condTE) %in% "control")])
CUdf$control = rep(norm_condTE[,"control"],sum(!(colnames(norm_condTE) %in% "control")))
CUdf$diff = CUdf$value - CUdf$control
CUdf = rbind(CUdf,CUdf) # concatenate as many times as Glu codons
CUdf$codon = c(rep("GluGAA",nrow(CUdf)/2),rep("GluGAG",nrow(CUdf)/2))
CUdf$RSCU = c(rep(relcodus[rownames(norm_condTE),"GluGAA"],3),rep(relcodus[rownames(norm_condTE),"GluGAG"],3))
CUdf$diffgroup = sapply(CUdf$diff, function(x) if(x > 1){"UP"}else if(x < -1){"DOWN"}else{NA})
CUdf$diffgroup = factor(as.character(CUdf$diffgroup),levels=c("UP","DOWN"))

# Plot
ggplot(CUdf[!is.na(CUdf$diffgroup),], aes(x=RSCU, fill=diffgroup)) +
  facet_wrap( ~ codon*X2, ncol=3) +
  geom_density(alpha=0.4, size=0.5, color="white") +
  scale_fill_manual(values=c(scales::muted("red"), scales::muted("blue"))) +
  theme_classic()
ggsave("plots/codon_usage_stress_GluDist.pdf",width=6,height=2.5)
