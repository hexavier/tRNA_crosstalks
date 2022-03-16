library(reshape)

#### MODIFICATIONS ####
# Load data
zz=gzfile("data/simulated_reads/mismatchTable.csv.gz",'rt')
mod_dataset = read.delim(zz)
mod_dataset = cast(mod_dataset, isodecoder*canon_pos*bam*condition ~ type, sum, value = 'proportion')
mod_dataset$mismatches = rowSums(mod_dataset[,c("A","T","G","C")],na.rm=T)

# Check Phe-GAA-2
phegaa2 = mod_dataset[(mod_dataset$isodecoder=="Saccharomyces_cerevisiae_tRNA-Phe-GAA-2")&(mod_dataset$canon_pos %in% c("26","58")),]
phegaa2 = cast(phegaa2, condition ~ canon_pos, sum, value = 'mismatches')

# Record input info
phegaa2$in26 = sapply(phegaa2$condition, function(x) as.numeric(gsub("-",".",substr(strsplit(as.character(x),"_")[[1]][1],5,9))))
phegaa2$in58 = sapply(phegaa2$condition, function(x) as.numeric(gsub("-",".",substr(strsplit(as.character(x),"_")[[1]][2],5,9))))

# Plot
plot(phegaa2$`26`,phegaa2$in26)
plot(phegaa2$`58`,phegaa2$in58)
