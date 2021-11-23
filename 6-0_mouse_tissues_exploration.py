#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import matplotlib.lines as mlines

#%% Create dataframe with all features and targets
# Upload data
rawdata = pd.read_csv("data/QuantM_mouse_tissues/interdependences.tsv",sep="\t")

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
metadata = {"SRR10587154":"Cerebellum",
            "SRR10587155":"Cerebellum",
            "SRR10587156":"Cerebellum",
            "SRR10587157":"Cortex",
            "SRR10587158":"Cortex",
            "SRR10587159":"Cortex",
            "SRR10587160":"MedullaOblongata",
            "SRR10587161":"MedullaOblongata",
            "SRR10587162":"MedullaOblongata",
            "SRR10587163":"SpinalCord",
            "SRR10587164":"SpinalCord",
            "SRR10587165":"SpinalCord",
            "SRR10587166":"Liver",
            "SRR10587167":"Liver",
            "SRR10587168":"Liver",
            "SRR10587169":"Heart",
            "SRR10587170":"Heart",
            "SRR10587171":"Heart",
            "SRR10587172":"Tibialis",
            "SRR10587173":"Tibialis",
            "SRR10587174":"Tibialis",}

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
fig.savefig("plots/PCA_mouse_tissues.pdf")
