#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import matplotlib.lines as mlines

plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"
plt.rcParams['pdf.use14corefonts'] = True

#%% Create dataframe with all features and targets
# Upload data
rawdata = pd.read_csv("data/ALKBH1_KD_stress/mismatchTable.csv",compression='gzip',sep="\t")
rawdata["bam"] = [s.split("/")[1].split(".")[0] for s in rawdata["bam"]]
moddata = rawdata.loc[rawdata.canon_pos!="-"].pivot(index=["isodecoder","canon_pos","bam","cov"], columns='type', values='proportion')
moddata["mismatches"] = moddata.sum(min_count=1, axis=1)

#%% Build matrix
# Rows are samples and columns are crosstalks
unique_mods = list(set([s[0]+"_"+str(s[1]) for s in moddata.index]))
uniquesamples = list(set([s[2] for s in moddata.index]))
modsdf = pd.DataFrame(np.nan, index=uniquesamples, columns = unique_mods)
for row in moddata.index:
    if row[3]>100:
        modsdf.loc[row[2],row[0]+"_"+str(row[1])] = moddata.loc[row,"mismatches"]
modsdf = modsdf.dropna(axis=1)
modsdf = modsdf.loc[:,modsdf.mean(axis=0)>0.05]
modsdf.columns = [s[13:] for s in modsdf.columns]

#%% PCA
pca = PCA(n_components=2)
principalComponents = pca.fit_transform(modsdf)
expl_var = pca.explained_variance_ratio_
features = pd.DataFrame(pca.components_.transpose(), columns=["PCA1","PCA2"], index = modsdf.columns)

principalDf = pd.DataFrame(principalComponents,columns=["PCA1","PCA2"],index=modsdf.index)

#%% Establish color labels
metadata = {"TP-CK-16S-CW05_S8_L1_bc1_TTCC":"minus-AsO2_control-shRNA",
            "TP-CK-16S-CW05_S8_L1_bc2_TCTG":"minus-AsO2_control-shRNA",
            "TP-CK-16S-CW05_S8_L1_bc3_TGGT":"plus-AsO2_control-shRNA",
            "TP-CK-16S-CW05_S8_L1_bc4_CTGA":"plus-AsO2_control-shRNA",
            "TP-CK-16S-CW05_S8_L1_bc5_CCAT":"minus-AsO2_shRNAW",
            "TP-CK-16S-CW05_S8_L1_bc6_CATC":"minus-AsO2_shRNAW",
            "TP-CK-16S-CW05_S8_L1_bc7_GTAG":"minus-AsO2_shRNAW",
            "TP-CK-16S-CW05_S8_L1_bc8_GGTA":"plus-AsO2_shRNAV",
            "TP-CK-16S-CW05_S8_L1_bc9_GACT":"plus-AsO2_shRNAV",
            "TP-CK-16S-CW05_S8_L1_bc10_ACCA":"plus-AsO2_shRNAV"}

#%% Plot pca
fig = plt.figure(figsize = (8,2.5))
ax = fig.add_subplot(1,1,1) 
ax.set_xlabel(str('Dim 1 (%1.2f%%)' % (expl_var[0]*100)), fontsize = 16)
ax.set_ylabel(str('Dim 2 (%1.2f%%)' % (expl_var[1]*100)), fontsize = 16)
ax.set_title('Dimension Reduction', fontsize = 20)

# Color based on tropism
labels = ["minus-AsO2_control-shRNA","plus-AsO2_control-shRNA","minus-AsO2_shRNAW","plus-AsO2_shRNAV"]
colors = ["#304D63","#8FB9AA","#ED8975","#F2D096"]
handles = [mlines.Line2D([], [], marker='o', color=c, alpha=0.5) for c in colors]

for label, color in zip(labels,colors):
    idx = list([s in [k for k,v in metadata.items() if v==label] for s in principalDf.index])
    ax.scatter(principalDf.loc[idx, 'PCA1'], principalDf.loc[idx, 'PCA2'], c = color, s = 50)
fig.legend(handles,labels,bbox_to_anchor=(1.02,0.8), loc="upper left")
fig.tight_layout()
ax.grid()
fig.savefig("plots/PCA_mods_alkbh1_kdstress.pdf", bbox_inches='tight')

# Plot features
ordfeat = features.sort_values(by="PCA1", ascending=True)
ax = pd.concat([ordfeat.iloc[:10],ordfeat.iloc[-10:]]).plot.bar(y="PCA1")
ax.figure.tight_layout()
ax.figure.savefig("plots/PCA_mods_features_alkbh1_kdstress.pdf")
