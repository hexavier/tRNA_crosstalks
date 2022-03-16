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
rawdata = pd.read_csv("data/QuantM_mouse_tissues/mismatchTable.csv",compression="gzip",sep="\t")
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

#%% Compute variance of detected modifications
detectedmods = modsdf.T
detectedmods["std"] = detectedmods[modsdf.index].std(axis=1)
detectedmods["mean"] = detectedmods[modsdf.index].mean(axis=1)
detectedmods["CV"] = detectedmods["std"]/detectedmods["mean"]
detectedmods.sort_values(by="CV", ascending=False).iloc[:20].plot.bar(y="CV")

#%% PCA
pca = PCA(n_components=2)
principalComponents = pca.fit_transform(modsdf)
expl_var = pca.explained_variance_ratio_
features = pd.DataFrame(pca.components_.transpose(), columns=["PCA1","PCA2"], index = modsdf.columns)

principalDf = pd.DataFrame(principalComponents,columns=["PCA1","PCA2"],index=modsdf.index)

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
            "SRR10587174":"Tibialis"}

#%% Plot pca
fig = plt.figure(figsize = (8,2.5))
ax = fig.add_subplot(1,1,1) 
ax.set_xlabel(str('Dim 1 (%1.2f%%)' % (expl_var[0]*100)), fontsize = 16)
ax.set_ylabel(str('Dim 2 (%1.2f%%)' % (expl_var[1]*100)), fontsize = 16)
ax.set_title('Dimension Reduction', fontsize = 20)

# Color based on tropism
labels = ["Cerebellum","Cortex","MedullaOblongata","SpinalCord","Tibialis","Liver","Heart"]
colors = ["#304D63","#304D63","#304D63","#304D63","#8FB9AA","#ED8975","#F2D096"]
handles = [mlines.Line2D([], [], marker='o', color=c, alpha=0.5) for c in colors]

for label, color in zip(labels,colors):
    idx = list([s in [k for k,v in metadata.items() if v==label] for s in principalDf.index])
    ax.scatter(principalDf.loc[idx, 'PCA1'], principalDf.loc[idx, 'PCA2'], c = color, s = 50)
fig.legend(handles,labels,bbox_to_anchor=(1.02,0.8), loc="upper left")
fig.tight_layout()
ax.grid()
fig.savefig("plots/PCA_mods_mouse_tissues.pdf", bbox_inches='tight')

# Plot features
ordfeat = features.sort_values(by="PCA1", ascending=True)
ax = pd.concat([ordfeat.iloc[:10],ordfeat.iloc[-10:]]).plot.bar(y="PCA1")
ax.figure.tight_layout()
ax.figure.savefig("plots/PCA_mods_features_mouse_tissues.pdf")
