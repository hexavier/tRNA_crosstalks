# -*- coding: utf-8 -*-
"""
Created on Wed Oct 13 11:01:53 2021

@author: user
"""
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Fetch data files
files = {"tRNA-Leu-CAA-4":["/home/xhernandez/midway/Human_stress_polysome/2_mimseq_results/single_read_data/TP-CK-4S-CW1_S1_L1_bc1_TTCC/Homo_sapiens_tRNA-Leu-CAA-4-1.tsv.gz",
						   "/home/xhernandez/midway/Human_stress_polysome/2_mimseq_results/single_read_data/TP-CK-4S-CW1_S1_L1_bc2_TCTG/Homo_sapiens_tRNA-Leu-CAA-4-1.tsv.gz",
						   "/home/xhernandez/midway/Human_stress_polysome/2_mimseq_results/single_read_data/TP-CK-4S-CW1_S1_L1_bc3_TGGT/Homo_sapiens_tRNA-Leu-CAA-4-1.tsv.gz"],
		 "tRNA-Leu-CAA-1":["/home/xhernandez/midway/Human_stress_polysome/2_mimseq_results/single_read_data/TP-CK-4S-CW1_S1_L1_bc1_TTCC/Homo_sapiens_tRNA-Leu-CAA-1-1.tsv.gz",
						 "/home/xhernandez/midway/Human_stress_polysome/2_mimseq_results/single_read_data/TP-CK-4S-CW1_S1_L1_bc2_TCTG/Homo_sapiens_tRNA-Leu-CAA-1-1.tsv.gz",
						  "/home/xhernandez/midway/Human_stress_polysome/2_mimseq_results/single_read_data/TP-CK-4S-CW1_S1_L1_bc3_TGGT/Homo_sapiens_tRNA-Leu-CAA-1-1.tsv.gz"],
		 "tRNA-Leu-AAG-2":["/home/xhernandez/midway/Human_stress_polysome/2_mimseq_results/single_read_data/TP-CK-4S-CW1_S1_L1_bc1_TTCC/Homo_sapiens_tRNA-Leu-AAG-2-1.tsv.gz",
						 "/home/xhernandez/midway/Human_stress_polysome/2_mimseq_results/single_read_data/TP-CK-4S-CW1_S1_L1_bc2_TCTG/Homo_sapiens_tRNA-Leu-AAG-2-1.tsv.gz",
						 "/home/xhernandez/midway/Human_stress_polysome/2_mimseq_results/single_read_data/TP-CK-4S-CW1_S1_L1_bc3_TGGT/Homo_sapiens_tRNA-Leu-AAG-2-1.tsv.gz"],
		 "tRNA-Leu-TAA-1":["/home/xhernandez/midway/Human_stress_polysome/2_mimseq_results/single_read_data/TP-CK-4S-CW1_S1_L1_bc1_TTCC/Homo_sapiens_tRNA-Leu-TAA-1-1.tsv.gz",
						 "/home/xhernandez/midway/Human_stress_polysome/2_mimseq_results/single_read_data/TP-CK-4S-CW1_S1_L1_bc2_TCTG/Homo_sapiens_tRNA-Leu-TAA-1-1.tsv.gz",
						  "/home/xhernandez/midway/Human_stress_polysome/2_mimseq_results/single_read_data/TP-CK-4S-CW1_S1_L1_bc3_TGGT/Homo_sapiens_tRNA-Leu-TAA-1-1.tsv.gz"]}
posinfo = pd.read_csv("/home/xhernandez/midway/Human_stress_polysome/2_mimseq_results/mods/mismatchTable.csv", sep="\t",
                      index_col=["isodecoder","pos"],usecols=["isodecoder","pos","canon_pos"],
                      dtype="category")
posinfo = posinfo[~posinfo.index.duplicated(keep='first')]

#%% Analyze coverage saparating by m22G26
covfiles = {"tRNA-Leu-CAA-4":[],"tRNA-Leu-CAA-1":[],"tRNA-Leu-AAG-2":[],"tRNA-Leu-TAA-1":[]}
for t in files.keys():
	for f in files[t]:
		ref = f.split("/")[-1].split(".")[0]
		ref = "-".join(ref.split("-")[:-1]) if not "chr" in ref else ref
		readsdf = pd.read_csv(f,sep="\t",compression="gzip",index_col="READ",dtype="category")
		# Remove Charging and canonize numbering
		readsdf.drop("Charged", axis=1, inplace=True)
		# Classify fragments
		colnames = [posinfo.loc[(ref,n),"canon_pos"] for n in readsdf.columns]
		colmap = {posinfo.loc[(ref,s),"canon_pos"]:int(s) for s in readsdf.columns if posinfo.loc[(ref,s),"canon_pos"]!="-"}
		# Classify fragments
		tempdf = readsdf.iloc[:,np.argmax(colnames=="30"):]
		readsdf.columns = [posinfo.loc[(ref,n),"canon_pos"] for n in readsdf.columns]
		readsdf["fragment"] = ["30to39" if colmap["30"]<=(int(s)-1)<=colmap["39"] else \
					  "40to49" if colmap["40"]<=(int(s)-1)<=colmap["49"] else \
						  "50to59" if colmap["50"]<=(int(s)-1)<=colmap["59"] else \
							  "60+" if colmap["60"]<=(int(s)-1) else np.nan for s in tempdf.apply(lambda x: x.index[np.argmax(x.isna())], axis=1).values]
		readsdf.dropna(subset=["26","fragment"], how="any", inplace=True)
		# Create coverage dataset
		outdf = pd.DataFrame(index=["unmod","mod"], columns=readsdf.columns[:-1])
		for i in outdf.index:
			if i=="unmod":
				tempdf = readsdf.loc[readsdf["26"]=="0"]
			else:
				tempdf = readsdf.loc[readsdf["26"]!="0"]
			outdf.loc[i] = tempdf.iloc[:,:-1].notna().sum(axis=0)
		# Remove CCA
		outdf = outdf.iloc[:,:-3]
		# Scale by max per row
		outdf = outdf.divide(outdf.max(axis=1), axis=0)
		# Reshape
		longdf = outdf.melt(var_name="pos",ignore_index=False)
		longdf.reset_index(inplace=True)
		# Save
		covfiles[t].append(longdf)
	
#%% Plot
fig, ax = plt.subplots(nrows=len(covfiles), ncols=1, figsize=(10,8))
for r,t in enumerate(covfiles.keys()):
	for longdf in covfiles[t]:
		sns.lineplot(data=longdf,x="pos", y="value",hue="index",
			   palette=["#ed8975","#304d63"], linewidth=1.5,ax= ax[r])
		ax[r].set_ylabel("Reads coverage")
		ax[r].set_xlabel(t)
		ax[r].tick_params(axis='x', labelsize=6)
fig.tight_layout()
# save
fig.savefig('plots/Leu_fragments_coverage.svg', bbox_inches='tight')
