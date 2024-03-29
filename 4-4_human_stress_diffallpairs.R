library(ggplot2)
library(patchwork)
library(reshape)
library(ggrepel)

### CHARGING ###
# Load data
raw_dataset = read.delim("data/Human_stress_polysome/CCAcounts.csv")

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

##### TOTAL TRNA #####
### CHARGING ###
# Charging differences between conditions
subdataset = dataset[grepl("total",dataset$cond),]
chrg_data = cast(subdataset, reference ~ sample, mean, value = 'charging')
wt_samples = unique(subdataset[subdataset$cond=="control_total","sample"])
avg_wt = rowMeans(chrg_data[,wt_samples],na.rm=T)
chrg_data = chrg_data[,!(colnames(chrg_data) %in% wt_samples)]
chrg_data[,!(colnames(chrg_data) %in% "reference")] = apply(chrg_data[,!(colnames(chrg_data) %in% "reference")],2,
                                                            function(x) x-avg_wt)

### ODDS RATIOS ###
# Load OR data
sigdataset = read.delim("data/Human_stress_polysome/interdependences_significant.tsv")
# Subset wt 58-Charg data
wt_data = sigdataset[(sigdataset$cond=="control_total")&(sigdataset$X.samples>0),]
# Find positions crosstalking with m1A58
pos = unique(c(apply(wt_data,1,function(x) paste0(x[c("reference","canon_var1")],collapse="_")),
                   apply(wt_data,1,function(x) paste0(x[c("reference","canon_var2")],collapse="_"))))
pairs = apply(wt_data,1,function(x) paste0(c(x["reference"],paste0(x[c("canon_var1","canon_var2")],collapse="-")),collapse="|"))

### MODIFICATIONS ###
# Load modification data
zz=gzfile("data/Human_stress_polysome/mismatchTable.csv.gz",'rt')   
mod_dataset = read.delim(zz)
mod_dataset = mod_dataset[grepl("total",mod_dataset$cond),]
mod_dataset = cast(mod_dataset, isodecoder*canon_pos*bam*condition ~ type, sum, value = 'proportion')
mod_dataset$mismatches = rowSums(mod_dataset[,c("A","T","G","C")],na.rm=T)
mod_dataset = mod_dataset[apply(as.data.frame(mod_dataset),1,function(x) paste0(x[c("isodecoder","canon_pos")],collapse="_")) %in% pos,]
# Reformat dataset
mod_data = cast(mod_dataset, isodecoder*canon_pos ~ bam, mean, value = 'mismatches')
wt_samples = unique(mod_dataset[mod_dataset$condition=="control_total","bam"])
avg_wt = rowMeans(mod_data[,wt_samples],na.rm=T)
mod_data = mod_data[,!(colnames(mod_data) %in% wt_samples)]
mod_data[,!(colnames(mod_data) %in% c("isodecoder","canon_pos"))] = apply(mod_data[,!(colnames(mod_data) %in% c("isodecoder","canon_pos"))],2,
                                                                          function(x) x-avg_wt)

# Create heatdf
conds = unique(mod_dataset[,c("bam","condition")])
totaldf = data.frame(row.names = 1:(length(pairs)*(ncol(chrg_data)-1)))
totaldf$reference = rep(pairs,(ncol(chrg_data)-1))
totaldf$pair = rep(sapply(pairs, function(x) strsplit(x,"\\|")[[1]][2]),(ncol(chrg_data)-1))
totaldf$sample = c(sapply(colnames(chrg_data[,2:ncol(chrg_data)]),function(x) rep(x,length(pairs))))
totaldf$oddsWT = rep(wt_data$avg_odds,(ncol(chrg_data)-1))
for (idx in rownames(totaldf)){
  p = totaldf[idx,"reference"]
  s = totaldf[idx,"sample"]
  c = conds[conds$bam==s,"condition"]
  ref = strsplit(p,"\\|")[[1]][1]
  splitted = strsplit(ref,"-")[[1]]
  totaldf[idx,"cond"] = c
  totaldf[idx,"aa"] = splitted[length(splitted)-2]
  pair = strsplit(strsplit(p,"\\|")[[1]][2],"-")[[1]]
  n=1
  for (ct in pair){
    if (ct=="Charged"){
      totaldf[idx,sprintf("CT%i",n)] = as.numeric(chrg_data[(chrg_data$reference==ref),s])
    }else{
      totaldf[idx,sprintf("CT%i",n)] = as.numeric(mod_data[(mod_data$canon_pos==ct)&(mod_data$isodecoder==ref),s])
    }
    n = n+1
  }
}
totaldf$reference = substr(totaldf$reference,14,100)
totaldf = totaldf[rowSums(is.na(totaldf))==0,]
# Keep only pairs in which one of the member has delta% > thres
thres = 0.03
subsetdf = totaldf[apply(totaldf,1,function(x) any(abs(as.numeric(x["CT1"]))>thres,abs(as.numeric(x["CT2"]))>thres)),]


##### PLOT #####
# Plot all
countsdf = c()
annotations = c()
for (c in unique(subsetdf$cond)){
  tempdf = data.frame(row.names = 1:4)
  tempdf$xpos = c(-Inf,-Inf,Inf,Inf)
  tempdf$ypos = c(-Inf, Inf,-Inf,Inf)
  tempdf$hjustvar = c(-0.5,-0.5,1.5,1.5)
  tempdf$vjustvar = c(-1,2,-1,2)
  tempdf$cond = c
  nq1_ORp = sum((subsetdf[subsetdf$cond==c,"CT1"]>0)&(subsetdf[subsetdf$cond==c,"CT2"]>0)&(subsetdf[subsetdf$cond==c,"oddsWT"]>0))
  nq1_ORn = sum((subsetdf[subsetdf$cond==c,"CT1"]>0)&(subsetdf[subsetdf$cond==c,"CT2"]>0)&(subsetdf[subsetdf$cond==c,"oddsWT"]<0))
  nq2_ORp = sum((subsetdf[subsetdf$cond==c,"CT1"]<0)&(subsetdf[subsetdf$cond==c,"CT2"]>0)&(subsetdf[subsetdf$cond==c,"oddsWT"]>0))
  nq2_ORn = sum((subsetdf[subsetdf$cond==c,"CT1"]<0)&(subsetdf[subsetdf$cond==c,"CT2"]>0)&(subsetdf[subsetdf$cond==c,"oddsWT"]<0))
  nq3_ORp = sum((subsetdf[subsetdf$cond==c,"CT1"]<0)&(subsetdf[subsetdf$cond==c,"CT2"]<0)&(subsetdf[subsetdf$cond==c,"oddsWT"]>0))
  nq3_ORn = sum((subsetdf[subsetdf$cond==c,"CT1"]<0)&(subsetdf[subsetdf$cond==c,"CT2"]<0)&(subsetdf[subsetdf$cond==c,"oddsWT"]<0))
  nq4_ORp = sum((subsetdf[subsetdf$cond==c,"CT1"]>0)&(subsetdf[subsetdf$cond==c,"CT2"]<0)&(subsetdf[subsetdf$cond==c,"oddsWT"]>0))
  nq4_ORn = sum((subsetdf[subsetdf$cond==c,"CT1"]>0)&(subsetdf[subsetdf$cond==c,"CT2"]<0)&(subsetdf[subsetdf$cond==c,"oddsWT"]<0))
  tempdf$text = c(sprintf("OR>1 = %i\nOR<1 = %i",nq3_ORp,nq3_ORn),
                  sprintf("OR>1 = %i\nOR<1 = %i",nq2_ORp,nq2_ORn),
                  sprintf("OR>1 = %i\nOR<1 = %i",nq4_ORp,nq4_ORn),
                  sprintf("OR>1 = %i\nOR<1 = %i",nq1_ORp,nq1_ORn))
  annotations = rbind(annotations,tempdf)
  nexpected = nq1_ORp+nq3_ORp+nq2_ORn+nq4_ORn
  ntotal = nrow(subsetdf[subsetdf$cond==c,])
  pval = binom.test(nexpected, ntotal, alternative = "greater")
  print(sprintf("Binomial test p-value (%s): %f",c,pval$p.value))
  # Counts df
  tempcounts = data.frame(row.names = 1:8)
  tempcounts$quadrant = factor(c("UP-UP","DOWN-DOWN","UP-DOWN","DOWN-UP","UP-UP","DOWN-DOWN","UP-DOWN","DOWN-UP"),levels=c("UP-UP","DOWN-DOWN","UP-DOWN","DOWN-UP"))
  tempcounts$odds = c("OR>1","OR>1","OR>1","OR>1","OR<1","OR<1","OR<1","OR<1")
  tempcounts$counts = c(nq1_ORp,nq3_ORp,nq4_ORp,nq2_ORp,nq1_ORn,nq3_ORn,nq4_ORn,nq2_ORn)
  tempcounts$success = c(TRUE,TRUE,FALSE,FALSE,FALSE,FALSE,TRUE,TRUE)
  tempcounts$cond = c
  countsdf = rbind(countsdf,tempcounts)
}

ggplot(subsetdf, aes(x=CT1,y=CT2,label=reference)) +
  facet_wrap( ~ cond, ncol=3, scales = "free") +
  geom_point(shape = 16, stroke=0.5, size=1.5, alpha=0.5, aes(color = oddsWT)) +
  scale_color_gradient2(limits=c(-2,2), oob = scales::squish) +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=text),size=2.5) +
  theme_classic() +
  labs(title ="Total tRNA", x ="Delta_CT1", y = "Delta_CT2")

ggsave(sprintf("plots/human_stress_allpairs_oddsWT_allsamples_delta%i.pdf",thres*100),width=8,height=2.7)

# Barplot of counts
ggplot(countsdf, aes(fill=quadrant, y=counts, x=success)) + 
  facet_wrap( ~ cond, ncol=3) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values=c("#0057B7","#98C6FA","#FFDD00","#F9F0B8")) +
  theme_classic() +
  labs(x ="Success", y = "# pairs")
ggsave(sprintf("plots/human_stress_allpairs_oddsWT_allsamples_barplot_delta%i.pdf",thres*100),width=6,height=3)

# Test several thres
thresholds = c(0,0.01,0.02,0.03,0.04,0.05)
validation = data.frame(threshold=character(),condition=character(),nTRUE=numeric(),nFALSE=numeric(),
                        pvalue=numeric())
for (thr in thresholds){
  subsetdf = totaldf[apply(totaldf,1,function(x) any(abs(as.numeric(x["CT1"]))>thr,abs(as.numeric(x["CT2"]))>thr)),]
  for (c in unique(subsetdf$cond)){
    nq1_ORp = sum((subsetdf[subsetdf$cond==c,"CT1"]>0)&(subsetdf[subsetdf$cond==c,"CT2"]>0)&(subsetdf[subsetdf$cond==c,"oddsWT"]>0))
    nq1_ORn = sum((subsetdf[subsetdf$cond==c,"CT1"]>0)&(subsetdf[subsetdf$cond==c,"CT2"]>0)&(subsetdf[subsetdf$cond==c,"oddsWT"]<0))
    nq2_ORp = sum((subsetdf[subsetdf$cond==c,"CT1"]<0)&(subsetdf[subsetdf$cond==c,"CT2"]>0)&(subsetdf[subsetdf$cond==c,"oddsWT"]>0))
    nq2_ORn = sum((subsetdf[subsetdf$cond==c,"CT1"]<0)&(subsetdf[subsetdf$cond==c,"CT2"]>0)&(subsetdf[subsetdf$cond==c,"oddsWT"]<0))
    nq3_ORp = sum((subsetdf[subsetdf$cond==c,"CT1"]<0)&(subsetdf[subsetdf$cond==c,"CT2"]<0)&(subsetdf[subsetdf$cond==c,"oddsWT"]>0))
    nq3_ORn = sum((subsetdf[subsetdf$cond==c,"CT1"]<0)&(subsetdf[subsetdf$cond==c,"CT2"]<0)&(subsetdf[subsetdf$cond==c,"oddsWT"]<0))
    nq4_ORp = sum((subsetdf[subsetdf$cond==c,"CT1"]>0)&(subsetdf[subsetdf$cond==c,"CT2"]<0)&(subsetdf[subsetdf$cond==c,"oddsWT"]>0))
    nq4_ORn = sum((subsetdf[subsetdf$cond==c,"CT1"]>0)&(subsetdf[subsetdf$cond==c,"CT2"]<0)&(subsetdf[subsetdf$cond==c,"oddsWT"]<0))
    nexpected = nq1_ORp+nq3_ORp+nq2_ORn+nq4_ORn
    ntotal = nrow(subsetdf[subsetdf$cond==c,])
    pval = binom.test(nexpected, ntotal, alternative = "greater")
    validation[nrow(validation) + 1,] <- c(thr,c,nexpected,ntotal-nexpected,pval$p.value)
  }
}
# Try also taking corrected pval <0.05
difftab = read.delim("results/human_stress_diffmodchrg.tsv")
difftab = difftab[difftab$p_cor<0.05,]
for (c in unique(totaldf$cond)){
  difftabcond = difftab[difftab$condition==c,]
  subsetdf = totaldf[apply(totaldf,1,function(x) any(sprintf("%s|%s",strsplit(x["reference"],"\\|")[[1]][1],strsplit(x["pair"],"-")[[1]][1]) %in% difftabcond$fullname,
                                                     sprintf("%s|%s",strsplit(x["reference"],"\\|")[[1]][1],strsplit(x["pair"],"-")[[1]][2]) %in% difftabcond$fullname)),]
  nq1_ORp = sum((subsetdf[subsetdf$cond==c,"CT1"]>0)&(subsetdf[subsetdf$cond==c,"CT2"]>0)&(subsetdf[subsetdf$cond==c,"oddsWT"]>0))
  nq1_ORn = sum((subsetdf[subsetdf$cond==c,"CT1"]>0)&(subsetdf[subsetdf$cond==c,"CT2"]>0)&(subsetdf[subsetdf$cond==c,"oddsWT"]<0))
  nq2_ORp = sum((subsetdf[subsetdf$cond==c,"CT1"]<0)&(subsetdf[subsetdf$cond==c,"CT2"]>0)&(subsetdf[subsetdf$cond==c,"oddsWT"]>0))
  nq2_ORn = sum((subsetdf[subsetdf$cond==c,"CT1"]<0)&(subsetdf[subsetdf$cond==c,"CT2"]>0)&(subsetdf[subsetdf$cond==c,"oddsWT"]<0))
  nq3_ORp = sum((subsetdf[subsetdf$cond==c,"CT1"]<0)&(subsetdf[subsetdf$cond==c,"CT2"]<0)&(subsetdf[subsetdf$cond==c,"oddsWT"]>0))
  nq3_ORn = sum((subsetdf[subsetdf$cond==c,"CT1"]<0)&(subsetdf[subsetdf$cond==c,"CT2"]<0)&(subsetdf[subsetdf$cond==c,"oddsWT"]<0))
  nq4_ORp = sum((subsetdf[subsetdf$cond==c,"CT1"]>0)&(subsetdf[subsetdf$cond==c,"CT2"]<0)&(subsetdf[subsetdf$cond==c,"oddsWT"]>0))
  nq4_ORn = sum((subsetdf[subsetdf$cond==c,"CT1"]>0)&(subsetdf[subsetdf$cond==c,"CT2"]<0)&(subsetdf[subsetdf$cond==c,"oddsWT"]<0))
  nexpected = nq1_ORp+nq3_ORp+nq2_ORn+nq4_ORn
  ntotal = nrow(subsetdf[subsetdf$cond==c,])
  pval = binom.test(nexpected, ntotal, alternative = "greater")
  validation[nrow(validation) + 1,] <- c("FDR>0.05",c,nexpected,ntotal-nexpected,pval$p.value)
}
# Save
write.table(validation, "results/binomials_thesholds.tsv", sep = "\t", quote = F)

### Check subsets of AA  and pairs###
thres = 0.03
# AA
subsetdf = totaldf[apply(totaldf,1,function(x) any(abs(as.numeric(x["CT1"]))>thres,abs(as.numeric(x["CT2"]))>thres)),]
ggplot(subsetdf, aes(x=CT1,y=CT2,label=reference)) +
  facet_wrap( ~ cond*aa, ncol=3, scales = "free") +
  geom_point(shape = 16, stroke=0.5, size=2, alpha=0.8, aes(color = oddsWT)) +
  scale_color_gradient2(limits=c(-1,1), oob = scales::squish) +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  theme_classic() +
  labs(title ="Total tRNA", x ="Delta_CT1", y = "Delta_CT2")
ggsave(sprintf("plots/human_stress_allpairs_oddsWT_allsamples_delta%i_AAfamilies.pdf",thres*100),width=15,height=80, limitsize = F)

# Pairs with at least 10 points
paircounts = data.frame(table(subsetdf[,c("cond","pair")]))
pairs5 = apply(paircounts[paircounts$Freq>=5,c("cond","pair")],1,function(x) paste0(x,collapse="-"))
subsubset = subsetdf[apply(subsetdf[,c("cond","pair")],1,function(x) paste0(x,collapse="-")) %in% pairs5,]
ggplot(subsubset, aes(x=CT1,y=CT2,label=reference)) +
  facet_wrap( ~ cond*pair, ncol=3, scales = "free") +
  geom_point(shape = 16, stroke=0.5, size=2, alpha=0.8, aes(color = oddsWT)) +
  scale_color_gradient2(limits=c(-1,1), oob = scales::squish) +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  theme_classic() +
  labs(title ="Total tRNA", x ="Delta_CT1", y = "Delta_CT2")
ggsave(sprintf("plots/human_stress_allpairs_oddsWT_allsamples_delta%i_pairs.pdf",thres*100),width=15,height=100, limitsize = F)

# Check pairs and AA that are classified the best
checkdf = data.frame(by=character(),class=character(),cond=character(),nTRUE=numeric(),nFALSE=numeric(),percTRUE=numeric(),pval=numeric())
checkfunc <- function(tempdf){
  ntotal = nrow(tempdf)
  nq1_ORp = sum((tempdf["CT1"]>0)&(tempdf["CT2"]>0)&(tempdf["oddsWT"]>0))
  nq2_ORn = sum((tempdf["CT1"]<0)&(tempdf["CT2"]>0)&(tempdf["oddsWT"]<0))
  nq3_ORp = sum((tempdf["CT1"]<0)&(tempdf["CT2"]<0)&(tempdf["oddsWT"]>0))
  nq4_ORn = sum((tempdf["CT1"]>0)&(tempdf["CT2"]<0)&(tempdf["oddsWT"]<0))
  nexpected = nq1_ORp+nq3_ORp+nq2_ORn+nq4_ORn
  pval = binom.test(nexpected, ntotal)$p.value
  return(c(nexpected,ntotal-nexpected,nexpected/ntotal,pval))
}
for (c in unique(totaldf$cond)){
  for (aa in unique(subsetdf$aa)){
    tempdf = subsetdf[(subsetdf$aa==aa)&(subsetdf$cond==c),]
    if (nrow(tempdf)>5){
      checkdf[nrow(checkdf) + 1,] <- c("AA",aa,c,checkfunc(tempdf))
    }
  }
  for (pair in unique(subsetdf$pair)){
    tempdf = subsetdf[(subsetdf$pair==pair)&(subsetdf$cond==c),]
    if (nrow(tempdf)>5){
      checkdf[nrow(checkdf) + 1,] <- c("pair",pair,c,checkfunc(tempdf))
    }
  }
  for (odds in c(0,0.5,1,2)){
    tempdf = subsetdf[(abs(subsetdf$oddsWT)>odds)&(subsetdf$cond==c),]
    if (nrow(tempdf)>5){
      checkdf[nrow(checkdf) + 1,] <- c("OR",sprintf("Abs(OR)>%s",odds),c,checkfunc(tempdf))
    }
  }
}

# Plot
checkdf$percTRUE = as.numeric(checkdf$percTRUE)
checkdf$avgperc = apply(checkdf,1,function(x) mean(checkdf[(checkdf$class==x["class"]),"percTRUE"],na.rm=T))
checkdf$class = factor(as.character(checkdf$class),levels=unique(checkdf[order(checkdf$avgperc,decreasing = T),"class"]))
checkdf$total = as.numeric(checkdf$nTRUE) + as.numeric(checkdf$nFALSE)

ggplot(checkdf[checkdf$by=="OR",], aes(x=class,y=percTRUE,color=cond,group=cond,label=total))+
  geom_hline(yintercept=0.5, linetype="dashed") +
  geom_point(shape = 16, size=2) +
  geom_line() +
  geom_text_repel(size=3) +
  theme_classic()
ggsave(sprintf("plots/human_stress_allpairs_oddsWT_checkOR_delta%i.pdf",thres*100),width=6,height=2.5)

ggplot(checkdf[checkdf$by=="AA",], aes(x=class,y=percTRUE,color=cond,size=total))+
  geom_hline(yintercept=0.5, linetype="dashed") +
  geom_point(shape = 16,alpha=0.6) +
  theme_classic()
ggsave(sprintf("plots/human_stress_allpairs_oddsWT_checkAA_delta%i.pdf",thres*100),width=9,height=2.5)

ggplot(checkdf[checkdf$by=="pair",], aes(x=class,y=percTRUE,color=cond,size=total))+
  geom_hline(yintercept=0.5, linetype="dashed") +
  geom_point(shape = 16,alpha=0.6) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))
ggsave(sprintf("plots/human_stress_allpairs_oddsWT_checkpairs_delta%i.pdf",thres*100),width=9,height=2.5)
