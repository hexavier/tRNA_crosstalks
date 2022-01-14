library(gplots)
library(reshape)

# Load data
raw_dataset = read.delim("data/Human_stress_polysome/interdependences_morethan2samples.tsv")

# Reformat
dataset = c()
for (r in rownames(raw_dataset)){
  tempdata = data.frame(row.names = 1:raw_dataset[r,"X.samples"])
  tempdata$reference = substr(raw_dataset[r,"reference"],14,100)
  tempdata$pair = paste0(raw_dataset[r,c("canon_var1","canon_var2")],collapse="-")
  tempdata$cond = strsplit(raw_dataset[r,"cond"],"_")[[1]][1]
  tempdata$fraction = strsplit(raw_dataset[r,"cond"],"_")[[1]][2]
  odds = sapply(strsplit(raw_dataset[r,"odds"],"\\[|\\]| ")[[1]], function(x) if (x!=""){as.numeric(x)}else{NA})
  tempdata$odds = odds[!is.na(odds)]
  pval = sapply(strsplit(raw_dataset[r,"pval_cor"],"\\[|\\]| ")[[1]], function(x) if (x!=""){as.numeric(x)}else{NA})
  tempdata$fisher = pval[!is.na(pval)]<0.05
  splitted = strsplit(raw_dataset[r,"reference"],"-")[[1]]
  tempdata$wobble = if(substr(splitted[length(splitted)-1],1,1) %in% c("A","T")){"AT"}else{"GC"}
  tempdata$TC = if(substr(splitted[length(splitted)-1],1,1) %in% "T"){"T"}else if(substr(splitted[length(splitted)-1],1,1) %in% "C"){"C"}else{"AG"}
  tempdata$aa = splitted[length(splitted)-2]
  dataset = rbind(dataset,tempdata)
}
# Remove mitochondrial
dataset = dataset[!grepl("mito",dataset$reference),]
# Keep only 58-Charging
dataset = dataset[dataset$pair=="58-Charged",]
# Separate poly and total
totaldata = dataset[dataset$fraction=="total",]
polydata = dataset[dataset$fraction=="poly",]

# Reformat dataset
totaldf = cast(totaldata, reference ~ cond, mean, value = 'odds')
rownames(totaldf) = totaldf$reference; totaldf = totaldf[,-c(1)]
polydf = cast(polydata, reference ~ cond, mean, value = 'odds')
rownames(polydf) = polydf$reference; polydf = polydf[,-c(1)]

# Keep only refs that have the same sign in total and poly
commonref = rownames(totaldf)[rownames(totaldf) %in% rownames(polydf)]
diffdf = sapply(commonref,function(x) if(all(na.omit((polydf[x,]>0)==(totaldf[x,]>0)))){abs(polydf[x,]) - abs(totaldf[x,])}else{rep(NA,ncol(totaldf))})
diffdf = as.matrix(apply(diffdf[,colSums(is.na(diffdf))<4],2,unlist))
rownames(diffdf) = colnames(totaldf)

# Heatmap
heatmap.2(diffdf, col=bluered, trace="none",Rowv=F,margins=c(10,10))

#### ONLY CONTROL ####
# Ctrl
totaldf = totaldf[,1,drop=F]; polydf = polydf[,1,drop=F]
# Keep only refs that have the same sign in total and poly
commonref = rownames(totaldf)[rownames(totaldf) %in% rownames(polydf)]
diffdf = matrix(sapply(commonref,function(x) if(all(na.omit((polydf[x,]>0)==(totaldf[x,]>0)))){abs(polydf[x,]) - abs(totaldf[x,])}else{rep(NA,ncol(totaldf))}))
colnames(diffdf) = colnames(totaldf); rownames(diffdf) = commonref
diffdf = diffdf[!is.na(diffdf),,drop=F]
diffdf = diffdf[order(diffdf),,drop=F]

# Heatmap
hclust.ward <- function(x) hclust(x, method="ward.D2")
pdf(file="plots/human_stress_totalVSpoly_odds_control_heatmap.pdf", width=12, height=4.5)
heatmap.2(t(cbind(diffdf,diffdf)), col=bluered, trace="none",Rowv=F,Colv=T,margins=c(7,10),hclustfun=hclust.ward)
dev.off()
