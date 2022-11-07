library(reshape)

### LOAD MOFIFICATION AND CHARGING DATA ###
# Modification data
zz=gzfile("data/Human_stress_polysome/mismatchTable.csv.gz",'rt')
mod_dataset = read.delim(zz)
mod_data = cast(mod_dataset, isodecoder*canon_pos*bam*condition*cov ~ type, sum, value = 'proportion')
mod_data$mismatches = apply(mod_data[,c("A","T","G","C")],1,sum,na.rm=T)
mod_data$nmod = mod_data$mismatches*mod_data$cov
mod_data$nunmod = mod_data$cov - mod_data$nmod
mod_data$reference = substr(mod_data$isodecoder,14,100)

# Charging data
# Load data
raw_dataset = read.delim("data/Human_stress_polysome/CCAcounts.csv")
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
chrg_data$reference = substr(chrg_data$reference,14,100)

# Merge datasets
df1 = data.frame(mod_data[,c("bam","reference","condition","mismatches","nmod","nunmod","canon_pos")]); colnames(df1) = c("bam","reference","condition","perc","yes","no","pos")
df2 = data.frame(chrg_data[,c("bam","reference","cond","charging","cca","cc")]); df2$pos = "Charged" ; colnames(df2) = c("bam","reference","condition","perc","yes","no","pos")
dataset = rbind(df1,df2)
dataset$fullname = sprintf("%s|%s",dataset$reference,dataset$pos)

# Filter out poly conditions
dataset = dataset[(dataset$condition %in% c("control_total","heatshock_total","H2O2_total","NaAsO2_total")),]

### ODDS RATIO TESTS ###
outdf = c()
for (cond in c("heatshock_total","H2O2_total","NaAsO2_total")){
  tempdf = data.frame(row.names = 1:length(unique(dataset$fullname)))
  tempdf$fullname = unique(dataset$fullname)
  tempdf$condition = cond
  n=1
  for (r in tempdf$fullname){
    ref = strsplit(r,"\\|")[[1]][1]; tempdf[n,"reference"] = ref
    pos = strsplit(r,"\\|")[[1]][2]; tempdf[n,"pos"] = pos
    # Counting
    wtyes = sum(dataset[(dataset$condition %in% "control_total")&(dataset$reference %in% ref)&(dataset$pos %in% pos),"yes"])
    wtno = sum(dataset[(dataset$condition %in% "control_total")&(dataset$reference %in% ref)&(dataset$pos %in% pos),"no"])
    if (((wtyes+wtno)>0)&&(wtyes/(wtyes+wtno))>0.05){ # more than 5% mismatch
      condyes = sum(dataset[(dataset$condition %in% cond)&(dataset$reference %in% ref)&(dataset$pos %in% pos),"yes"])
      condno = sum(dataset[(dataset$condition %in% cond)&(dataset$reference %in% ref)&(dataset$pos %in% pos),"no"])
      tempdf[n,"log2OR"] = log2((condyes/wtyes)/(condno/wtno))
      tempdf[n,"pvalue"] = chisq.test(matrix(c(condyes,wtyes,condno,wtno),nrow=2,ncol=2))$p.value
    }
    n=n+1
  }
  tempdf = tempdf[!is.na(tempdf$log2OR),]
  outdf = rbind(outdf,tempdf)
}
# FDR correction
outdf$p_cor = p.adjust(outdf$pvalue,method="fdr")

# Save
write.table(outdf, "results/human_stress_diffmodchrg.tsv", sep = "\t", row.names = F, quote = F)
