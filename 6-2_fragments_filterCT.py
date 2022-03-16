# -*- coding: utf-8 -*-
"""
Created on Tue Oct 19 12:00:16 2021

@author: user
"""
import pandas as pd
import numpy as np
import sys

# Load data
allinterdep = pd.read_csv("results/crosstalks_fragments.tsv", sep="\t", dtype={"canon_var1":"str","var2":"str"})
metadata = pd.read_csv("data/Human_stress_polysome/readme_mimseq.txt", sep="\t", names=["file","cond"])

#%% Analyze significant changes
# Filter out unsignificant changes
siginterdep = allinterdep.loc[allinterdep.pval_corrected < 0.05]
# Log2 odds ratio
siginterdep["log2_odds_ratio"] = np.log2(siginterdep["odds_ratio"])
# Merge samples
mergedoutdf = pd.DataFrame()
for ref in siginterdep.ref.unique():
    tempdf = siginterdep.loc[siginterdep.ref==ref]
    outdf = pd.DataFrame()
    for cond in metadata.cond.unique():
        samples_cond = [s.split("/")[-1].split(".")[0] for s in metadata.loc[metadata.cond==cond,"file"]]
        outdf_cond = pd.DataFrame(columns=["reference", "canon_var1", "var2", "cond", "#samples", "odds", "consistent"])
        tempdf_cond = tempdf.loc[[s in samples_cond for s in tempdf["sample"]]]
        detectedCT = [(s[0],s[1]) for s in tempdf_cond[["canon_var1","var2"]].drop_duplicates().values]
        # Fill table
        for n,pair in enumerate(detectedCT):
            tempdf_cond_pair = tempdf_cond.loc[np.logical_and(tempdf_cond.canon_var1==pair[0],tempdf_cond.var2==pair[1])]
            outdf_cond.loc[n] = [ref, pair[0], pair[1], cond, tempdf_cond_pair.shape[0],
                                 tempdf_cond_pair["log2_odds_ratio"].values,
                                 (tempdf_cond_pair["log2_odds_ratio"]>0).unique().size==1]
        # Merge with output table
        outdf = pd.concat([outdf,outdf_cond],axis=0,ignore_index=True)
    mergedoutdf = pd.concat([mergedoutdf,outdf],axis=0,ignore_index=True)
# Filter out unconsistent changes
filtereddf = mergedoutdf.loc[mergedoutdf["consistent"]].reset_index(drop=True)
filtereddf["avg_odds"] = [np.nanmean(s) for s in filtereddf["odds"]]
# Save results
filtereddf.to_csv("results/crosstalks_fragments_significant.tsv", sep="\t", index=False)

#%% Analyze all consistent changes
# Log2 odds ratio
allinterdep["log2_odds_ratio"] = np.log2(allinterdep["odds_ratio"])
# Merge samples
mergedoutdf = pd.DataFrame()
for ref in allinterdep.ref.unique():
    tempdf = allinterdep.loc[allinterdep.ref==ref]
    outdf = pd.DataFrame()
    for cond in metadata.cond.unique():
        samples_cond = [s.split("/")[-1].split(".")[0] for s in metadata.loc[metadata.cond==cond,"file"]]
        outdf_cond = pd.DataFrame(columns=["reference", "canon_var1", "var2", "cond", "#samples", "odds", "consistent","pval_cor"])
        tempdf_cond = tempdf.loc[[s in samples_cond for s in tempdf["sample"]]]
        detectedCT = [(s[0],s[1]) for s in tempdf_cond[["canon_var1","var2"]].drop_duplicates().values]
        # Fill table
        for n,pair in enumerate(detectedCT):
            tempdf_cond_pair = tempdf_cond.loc[np.logical_and(tempdf_cond.canon_var1==pair[0],tempdf_cond.var2==pair[1])]
            outdf_cond.loc[n] = [ref, pair[0], pair[1], cond, tempdf_cond_pair.shape[0],
                                 tempdf_cond_pair["log2_odds_ratio"].values,
                                 (tempdf_cond_pair["log2_odds_ratio"]>0).unique().size==1,
				 tempdf_cond_pair["pval_corrected"].values]
        # Merge with output table
        outdf = pd.concat([outdf,outdf_cond],axis=0,ignore_index=True)
    mergedoutdf = pd.concat([mergedoutdf,outdf],axis=0,ignore_index=True)
# Filter out detected in at least 2 samples
filtereddf = mergedoutdf.loc[mergedoutdf["#samples"]>1].reset_index(drop=True)
filtereddf["avg_odds"] = [np.nanmean(s) for s in filtereddf["odds"]]
# Save results
filtereddf.to_csv("results/crosstalks_fragments_morethan2samples.tsv", sep="\t", index=False)