library(reshape)
library(gplots)
library(RColorBrewer)
# Function for making selection rectangles around selection cells
makeRects <- function(cells,ny,nx){
  coords = expand.grid(ny:1, 1:nx)[cells=="OR>1",]
  xl=coords[,2]-0.49
  yb=coords[,1]-0.49
  xr=coords[,2]+0.49
  yt=coords[,1]+0.49
  rect(xl,yb,xr,yt,border="#f2d096ff",lwd=2)
  coords = expand.grid(ny:1, 1:nx)[cells=="OR<1",]
  xl=coords[,2]-0.49
  yb=coords[,1]-0.49
  xr=coords[,2]+0.49
  yt=coords[,1]+0.49
  rect(xl,yb,xr,yt,border="#304d63ff",lwd=2)
}

### DIFFERANTIAL MODIFICATION AND CHARGING ANALYSIS ###
# Modification data
zz=gzfile("data/Scer_mimseq/mismatchTable.csv.gz",'rt')
mod_dataset = read.delim(zz)
mod_data = cast(mod_dataset, isodecoder*canon_pos*bam*condition*cov ~ type, sum, value = 'proportion')
mod_data$mismatches = apply(mod_data[,c("A","T","G","C")],1,sum,na.rm=T)
mod_data$nmod = mod_data$mismatches*mod_data$cov
mod_data$nunmod = mod_data$cov - mod_data$nmod
mod_data$reference = substr(mod_data$isodecoder,26,100)

# Charging data
# Load data
raw_dataset = read.delim("data/Scer_mimseq/CCAcounts.csv")
# Compute % charging
samples = unique(raw_dataset$sample)
chrg_data = c()
for (s in samples){
  refs = unique(raw_dataset[raw_dataset$sample==s,"gene"])
  tempdata = data.frame(row.names = 1:length(refs))
  tempdata$reference = refs
  tempdata$bam = s
  tempdata$cond = raw_dataset[raw_dataset$sample==s,"condition"][1]
  for (r in rownames(tempdata)){
    subdata = raw_dataset[(raw_dataset$sample==s)&(raw_dataset$gene==tempdata[r,"reference"]),]
    cca = subdata[subdata$end=="CA","count"]
    cc = subdata[subdata$end=="CC","count"]
    tempdata[r,"charging"] = cca/(cca+cc)
    tempdata[r,"cca"] = cca
    tempdata[r,"cc"] = cc
  }
  chrg_data = rbind(chrg_data,tempdata)
}
chrg_data$reference = substr(chrg_data$reference,26,100)
chrg_data = chrg_data[chrg_data$cond %in% c("WT","Trm7"),]

# Merge datasets
df1 = data.frame(mod_data[,c("bam","reference","condition","mismatches","nmod","nunmod","canon_pos")]); colnames(df1) = c("bam","reference","condition","perc","yes","no","pos")
df2 = data.frame(chrg_data[,c("bam","reference","cond","charging","cca","cc")]); df2$pos = "Charged" ; colnames(df2) = c("bam","reference","condition","perc","yes","no","pos")
dataset = rbind(df1,df2)
dataset$fullname = sprintf("%s|%s",dataset$reference,dataset$pos)

### ODDS RATIO TESTS ###
outdf = c()
for (cond in c("Trm1","Trm10","Trm7")){
  tempdf = data.frame(row.names = 1:length(unique(dataset$fullname)))
  tempdf$fullname = unique(dataset$fullname)
  tempdf$condition = cond
  n=1
  for (r in tempdf$fullname){
    ref = strsplit(r,"\\|")[[1]][1]; tempdf[n,"reference"] = ref
    pos = strsplit(r,"\\|")[[1]][2]; tempdf[n,"pos"] = pos
    # Counting
    if (cond=="Trm7"){
      wtyes = sum(dataset[(dataset$condition %in% "WT")&(dataset$reference %in% ref)&(dataset$pos %in% pos),"yes"], na.rm = T)
      wtno = sum(dataset[(dataset$condition %in% "WT")&(dataset$reference %in% ref)&(dataset$pos %in% pos),"no"], na.rm = T)
    }else{
      wtyes = sum(dataset[(dataset$condition %in% "WTox0")&(dataset$reference %in% ref)&(dataset$pos %in% pos),"yes"], na.rm = T)
      wtno = sum(dataset[(dataset$condition %in% "WTox0")&(dataset$reference %in% ref)&(dataset$pos %in% pos),"no"], na.rm = T)
    }
    if ((wtyes+wtno)>0){ # more than 5% mismatch
      condyes = sum(dataset[(dataset$condition %in% cond)&(dataset$reference %in% ref)&(dataset$pos %in% pos),"yes"], na.rm = T)
      condno = sum(dataset[(dataset$condition %in% cond)&(dataset$reference %in% ref)&(dataset$pos %in% pos),"no"], na.rm = T)
      tempdf[n,"log2OR"] = log2((condyes/wtyes)/(condno/wtno))
      tempdf[n,"deltaperc"] = (condyes/(condyes+condno)) - (wtyes/(wtyes+wtno))
      if ((wtyes/(wtyes+wtno))>0.05){
        tempdf[n,"pvalue"] = chisq.test(matrix(c(condyes,wtyes,condno,wtno),nrow=2,ncol=2))$p.value
      }
    }
    n=n+1
  }
  tempdf = tempdf[!is.na(tempdf$log2OR),]
  outdf = rbind(outdf,tempdf)
}
# FDR correction
outdf$p_cor = p.adjust(outdf$pvalue,method="fdr")

# Save
write.table(outdf, "results/yeast_mutants_diffmodchrg.tsv", sep = "\t", row.names = F, quote = F)
sigdf = outdf; sigdf[((sigdf$p_cor>0.05)|(abs(sigdf$deltaperc)<0.03)|(is.na(sigdf$p_cor))),"deltaperc"] <- NA

### PLOT RESULTS ###
# Load crosstalks
sigdataset = read.delim("data/Scer_mimseq/interdependences_significant.tsv",colClasses=c(canon_var1="character",canon_var2="character"))
sigdataset$reference = substr(sigdataset$reference,26,100)
all_data = sigdataset[grepl("WT",sigdataset$cond),]
all_data$bin_odds = all_data$avg_odds>0
wt_data = cast(all_data, reference*canon_var1*canon_var2*bin_odds ~ cond, sum, value = "X.samples")
wt_data = as.data.frame(wt_data[rowSums(wt_data[,c("WT","WTox0")])>=2,])

overlapdf = data.frame(row.names = c("Trm1","Trm10","Trm7"))
### ANALYZE TRM1
# Heatmap
subdf = sigdf[sigdf$condition=="Trm1",]
# Keep only isodecoders where 26 changes
subdf = subdf[subdf$reference %in% subdf[subdf$pos=="26","reference"],]
mutdf = cast(subdf,reference ~ pos, value="deltaperc")
rownames(mutdf) = mutdf$reference; mutdf = mutdf[,-c(1)]
mutdf = mutdf[,order(as.numeric(colnames(mutdf)))]
# Crosstalks
subdf$oddsWT = apply(subdf,1,function(x) if (x["pos"]=="26"){"mutant"}else if(sum((((wt_data$canon_var1=="26")&(wt_data$canon_var2==x["pos"]))|((wt_data$canon_var1==x["pos"])&(wt_data$canon_var2=="26")))&(grepl(x["reference"],wt_data$reference)))==0){NA}else if(wt_data[(((wt_data$canon_var1=="26")&(wt_data$canon_var2==x["pos"]))|((wt_data$canon_var1==x["pos"])&(wt_data$canon_var2=="26")))&(grepl(x["reference"],wt_data$reference)),"bin_odds"]){"OR>1"}else{"OR<1"})
ctdf = cast(subdf,reference ~ pos, value="oddsWT")
rownames(ctdf) = ctdf$reference; ctdf = ctdf[,-c(1)]
ctdf = ctdf[,order(as.numeric(colnames(ctdf)))]
# Plot heatmap
submutdf = mutdf[(rowSums(is.na(mutdf))<(ncol(mutdf)-1))|(rowSums(is.na(ctdf))<(ncol(ctdf)-1)),colnames(mutdf) %in% c("9","26","32","34","37","58","Charged")]
ctdf = ctdf[(rowSums(is.na(mutdf))<(ncol(mutdf)-1))|(rowSums(is.na(ctdf))<(ncol(ctdf)-1)),colnames(ctdf) %in% c("9","26","32","34","37","58","Charged")]
pdf(file="plots/heatmap_yeast_trm1.pdf", width=4.5, height=6)
Colors=colorRampPalette(brewer.pal(11,"RdBu"))(19); Colors[10] <- "grey"
heatmap.2(as.matrix(as.data.frame(submutdf)), trace="none", margins = c(4,12), col=Colors,Rowv = NULL, Colv = NULL,
          na.color = "grey", add.expr={makeRects(as.matrix(as.data.frame(ctdf)),nrow(ctdf),ncol(ctdf))},
          breaks=c(seq(-1, -0.03, length.out=10),seq(0.03, 1, length.out=10)),symkey=F)
dev.off()
# Save results in df
overlapdf["Trm1",c("True","False")] = c(sum(sum((ctdf=="OR>1")&(mutdf<0), na.rm=T),sum((ctdf=="OR<1")&(mutdf>0), na.rm=T)),
                                        sum(sum((ctdf=="OR>1")&(mutdf>0), na.rm=T),sum((ctdf=="OR<1")&(mutdf<0), na.rm=T)))

### ANALYZE TRM10
# Heatmap
subdf = sigdf[sigdf$condition=="Trm10",]
# Keep only isodecoders where 9 changes
subdf = subdf[subdf$reference %in% subdf[subdf$pos=="9","reference"],]
mutdf = cast(subdf,reference ~ pos, value="deltaperc")
rownames(mutdf) = mutdf$reference; mutdf = mutdf[,-c(1)]
mutdf = mutdf[,order(as.numeric(colnames(mutdf)))]
# Crosstalks
subdf$oddsWT = apply(subdf,1,function(x) if (x["pos"]=="9"){"mutant"}else if(sum((((wt_data$canon_var1=="9")&(wt_data$canon_var2==x["pos"]))|((wt_data$canon_var1==x["pos"])&(wt_data$canon_var2=="9")))&(grepl(x["reference"],wt_data$reference)))==0){NA}else if(wt_data[(((wt_data$canon_var1=="9")&(wt_data$canon_var2==x["pos"]))|((wt_data$canon_var1==x["pos"])&(wt_data$canon_var2=="9")))&(grepl(x["reference"],wt_data$reference)),"bin_odds"]){"OR>1"}else{"OR<1"})
ctdf = cast(subdf,reference ~ pos, value="oddsWT")
rownames(ctdf) = ctdf$reference; ctdf = ctdf[,-c(1)]
ctdf = ctdf[,order(as.numeric(colnames(ctdf)))]
# Plot heatmap
submutdf = mutdf[(rowSums(is.na(mutdf))<(ncol(mutdf)-1))|(rowSums(is.na(ctdf))<(ncol(ctdf)-1)),colnames(mutdf) %in% c("9","26","32","34","37","58","Charged")]
ctdf = ctdf[(rowSums(is.na(mutdf))<(ncol(mutdf)-1))|(rowSums(is.na(ctdf))<(ncol(ctdf)-1)),colnames(ctdf) %in% c("9","26","32","34","37","58","Charged")]
pdf(file="plots/heatmap_yeast_trm10.pdf", width=4.5, height=4)
Colors=colorRampPalette(brewer.pal(11,"RdBu"))(19); Colors[10] <- "grey"
heatmap.2(as.matrix(as.data.frame(submutdf)), trace="none", margins = c(4,12), col=Colors,Rowv = NULL, Colv = NULL,
          na.color = "grey", add.expr={makeRects(as.matrix(as.data.frame(ctdf)),nrow(ctdf),ncol(ctdf))},
          breaks=c(seq(-1, -0.03, length.out=10),seq(0.03, 1, length.out=10)),symkey=F)
dev.off()
# Save results in df
overlapdf["Trm10",c("True","False")] = c(sum(sum((ctdf=="OR>1")&(mutdf<0), na.rm=T),sum((ctdf=="OR<1")&(mutdf>0), na.rm=T)),
                                        sum(sum((ctdf=="OR>1")&(mutdf>0), na.rm=T),sum((ctdf=="OR<1")&(mutdf<0), na.rm=T)))

### ANALYZE TRM7
# Trm7 only producing yW in Phe tRNA
# Heatmap
subdf = sigdf[sigdf$condition=="Trm7",]
# Keep only isodecoders where Phe-37 changes
subdf = subdf[subdf$reference %in% subdf[(subdf$pos=="37")&(subdf$reference %in% c("tRNA-Phe-GAA-1","tRNA-Phe-GAA-2")),"reference"],]
mutdf = cast(subdf,reference ~ pos, value="deltaperc")
rownames(mutdf) = mutdf$reference; mutdf = mutdf[,-c(1)]
mutdf = mutdf[,order(as.numeric(colnames(mutdf)))]
# Crosstalks
subdf$oddsWT = apply(subdf,1,function(x) if (x["pos"]=="37"){"mutant"}else if(sum((((wt_data$canon_var1=="37")&(wt_data$canon_var2==x["pos"]))|((wt_data$canon_var1==x["pos"])&(wt_data$canon_var2=="37")))&(grepl(x["reference"],wt_data$reference)))==0){NA}else if(wt_data[(((wt_data$canon_var1=="37")&(wt_data$canon_var2==x["pos"]))|((wt_data$canon_var1==x["pos"])&(wt_data$canon_var2=="37")))&(grepl(x["reference"],wt_data$reference)),"bin_odds"]){"OR>1"}else{"OR<1"})
ctdf = cast(subdf,reference ~ pos, value="oddsWT")
rownames(ctdf) = ctdf$reference; ctdf = ctdf[,-c(1)]
ctdf = ctdf[,order(as.numeric(colnames(ctdf)))]
# Plot heatmap
submutdf = mutdf[(rowSums(is.na(mutdf))<(ncol(mutdf)-1))|(rowSums(is.na(ctdf))<(ncol(ctdf)-1)),colnames(mutdf) %in% c("9","26","32","34","37","58","Charged")]
ctdf = ctdf[(rowSums(is.na(mutdf))<(ncol(mutdf)-1))|(rowSums(is.na(ctdf))<(ncol(ctdf)-1)),colnames(ctdf) %in% c("9","26","32","34","37","58","Charged")]
pdf(file="plots/heatmap_yeast_trm7.pdf", width=4.5, height=4)
Colors=colorRampPalette(brewer.pal(11,"RdBu"))(19); Colors[10] <- "grey"
heatmap.2(as.matrix(as.data.frame(submutdf)), trace="none", margins = c(4,12), col=Colors,Rowv = NULL, Colv = NULL,
          na.color = "grey", add.expr={makeRects(as.matrix(as.data.frame(ctdf)),nrow(ctdf),ncol(ctdf))},
          breaks=c(seq(-1, -0.03, length.out=10),seq(0.03, 1, length.out=10)),symkey=F)
dev.off()
# Save results in df
overlapdf["Trm7",c("True","False")] = c(sum(sum((ctdf=="OR>1")&(mutdf<0), na.rm=T),sum((ctdf=="OR<1")&(mutdf>0), na.rm=T)),
                                        sum(sum((ctdf=="OR>1")&(mutdf>0), na.rm=T),sum((ctdf=="OR<1")&(mutdf<0), na.rm=T)))
pie(colSums(overlapdf), labels = colnames(overlapdf), main="Predictions")
