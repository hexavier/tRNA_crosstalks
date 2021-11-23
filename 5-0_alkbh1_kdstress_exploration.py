#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import matplotlib.lines as mlines

#%% Create dataframe with all features and targets
# Upload data
rawdata = pd.read_csv("data/ALKBH1_KD_stress/interdependences.tsv",sep="\t")

#%% Build matrix
# Rows are samples and columns are crosstalks
uniqueCT = rawdata[["ref","canon_var1","canon_var2"]].dropna().drop_duplicates()
uniquesamples = rawdata["sample"].drop_duplicates()
CTdf = pd.DataFrame(0, index=uniquesamples, columns = ["_".join(uniqueCT.loc[s].values) for s in uniqueCT.index])
for row in rawdata.dropna().index:
    CTdf.loc[rawdata.loc[row,"sample"],"_".join(rawdata.loc[row,["ref","canon_var1","canon_var2"]].values)] = np.log2(rawdata.loc[row,"odds_ratio"])

#%% PCA
pca = PCA(n_components=2)
principalComponents = pca.fit_transform(CTdf)
expl_var = pca.explained_variance_ratio_
features = pd.DataFrame(pca.components_.transpose(), columns=["PCA1","PCA2"], index = CTdf.columns)

principalDf = pd.DataFrame(principalComponents,columns=["PCA1","PCA2"],index=CTdf.index)

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
fig = plt.figure(figsize = (8,6))
ax = fig.add_subplot(1,1,1) 
ax.set_xlabel(str('Dim 1 (%1.2f%%)' % (expl_var[0]*100)), fontsize = 15)
ax.set_ylabel(str('Dim 2 (%1.2f%%)' % (expl_var[1]*100)), fontsize = 15)
ax.set_title('Dimension Reduction', fontsize = 20)

# Color based on tropism
labels = list(set(metadata.values()))
cmap = plt.cm.get_cmap('jet') # set color map
colors = cmap(np.arange(len(labels))/(len(labels)-1))
handles = [mlines.Line2D([], [], marker='o', color=c, alpha=0.5) for c in colors]

for label, color in zip(labels,colors):
    idx = list([s in [k for k,v in metadata.items() if v==label] for s in principalDf.index])
    ax.scatter(principalDf.loc[idx, 'PCA1'], principalDf.loc[idx, 'PCA2'], c = color, s = 50, alpha=0.5)
    for lab in principalDf.index[idx]:
       ax.annotate(metadata[lab],[principalDf.loc[lab, 'PCA1'],principalDf.loc[lab, 'PCA2']],fontsize=8)
fig.legend(handles,labels,bbox_to_anchor=(1.02,0.7), loc="upper left")
fig.tight_layout()
ax.grid()
#fig.savefig("PCA_human_stress.pdf")
