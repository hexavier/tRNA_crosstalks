# -*- coding: utf-8 -*-
"""
Created on Wed Oct 13 11:01:53 2021

@author: user
"""
import pandas as pd
import numpy as np
from os import listdir
from os.path import isdir,isfile,join
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import multiprocessing
from functools import partial

# Fetch data files
dirpath = "/home/xhernandez/midway/Human_stress_polysome/2_mimseq_results/single_read_data"
thres = 0.05 # % of positions with mismatched to do test
samples = [f for f in listdir(dirpath) if isdir(join(dirpath, f))]
posinfo = pd.read_csv("/".join(dirpath.split("/")[:-1])+"/mods/mismatchTable.csv", sep="\t",
                      index_col=["isodecoder","pos"],usecols=["isodecoder","pos","canon_pos"],
                      dtype="category")
posinfo = posinfo[~posinfo.index.duplicated(keep='first')]

# Define function for parallelization
def analyze_1sample(s,dirpath,thres):
	outdf = pd.DataFrame(columns=["sample","ref","var1","var2","pval","odds_ratio","values"])
	n=0
	ref_files = [f for f in listdir(join(dirpath, s)) if isfile(join(dirpath, s, f))]
	for f in ref_files:
		# Load table
		ref = f.split(".")[0]
		ref = "-".join(ref.split("-")[:-1]) if not "chr" in ref else ref
		readsdf = pd.read_csv(join(dirpath, s, f),sep="\t",compression="gzip",index_col="READ",dtype="category")
		# Remove Charging and canonize numbering
		readsdf.drop("Charged", axis=1, inplace=True)
		# Exclude mito and low coverage genes
		if (ref in posinfo.index.get_level_values(0)) and "mito" not in ref:
			colnames = [posinfo.loc[(ref,n),"canon_pos"] for n in readsdf.columns]
			colmap = {posinfo.loc[(ref,s),"canon_pos"]:int(s) for s in readsdf.columns if posinfo.loc[(ref,s),"canon_pos"]!="-"}
			# Classify fragments
			tempdf = readsdf.iloc[:,np.argmax(colnames=="30"):]
			readsdf["fragment"] = ["30to39" if colmap["30"]<=(int(s)-1)<=colmap["39"] else \
						  "40to49" if colmap["40"]<=(int(s)-1)<=colmap["49"] else \
							  "50to59" if colmap["50"]<=(int(s)-1)<=colmap["59"] else \
								  "60+" if colmap["60"]<=(int(s)-1) else np.nan for s in tempdf.apply(lambda x: x.index[np.argmax(x.isna())], axis=1).values]
			readsdf.dropna(subset=["fragment"], inplace=True)
			for v1 in readsdf.columns[:-1]:
				for v2 in ["30to39","40to49","50to59"]:
					reads1 = readsdf[v1].dropna()
					reads2 = readsdf["fragment"].dropna()
					# Do test only if at least % positions are modified
					if reads1.shape[0]>0 and reads2.shape[0]>0:
						if sum(reads1!="0")/reads1.shape[0]>thres:
							# Keep only reads that contain both v1 and v2
							tempdf = readsdf[[v1,"fragment"]].copy()
							tempdf["fragment"] = [s if s in [v2,"60+"] else np.nan for s in tempdf["fragment"]]
							tempdf.dropna(how="any", inplace=True)
							# Build contingency table, var1 in rows and var 2 in columns
							counts = pd.concat([tempdf[v1]!="0",tempdf["fragment"]==v2], axis=1).value_counts()
							if counts.shape[0]==4:
								cont_tab = np.array([[counts.loc[True,True], counts.loc[True,False]],
													[counts.loc[False,True], counts.loc[False,False]]])
								oddsr, p = fisher_exact(cont_tab)
								outdf.loc[n] = [s,ref,v1,v2,p,oddsr,counts]
								n += 1
	return outdf

def pd_wrapper(dirpath,thres,samples):
    # Multiprocessed function
    pool = multiprocessing.Pool(8) # CHANGE DEPENDING ON THE ENVIRONMENT
    frozen_fun = partial(analyze_1sample, dirpath=dirpath, thres=thres)
    dfs = pool.map(func=frozen_fun, iterable=samples, chunksize=1)
    pool.close()
    pool.join()
    # Concatenate tables
    outdf = pd.concat(dfs, ignore_index=True)
    # Correct multiple comparisons
    outdf["pval_corrected"] = multipletests(outdf["pval"].values,method="fdr_bh")[1]
    # Include canonical positions
    outdf["canon_var1"] = ["Charged" if s[1]=="Charged" else posinfo.loc[(s[0],s[1]),"canon_pos"] if s[0] in posinfo.index.get_level_values(0) else "NA" for s in outdf[["ref","var1"]].to_numpy()]
    return outdf

if __name__ == '__main__':
    outdf = pd_wrapper(dirpath,thres,samples)
    outdf.to_csv("results/crosstalks_fragments.tsv",sep="\t",index=False)
