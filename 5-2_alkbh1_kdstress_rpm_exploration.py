#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.lines as mlines

plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"
plt.rcParams['pdf.use14corefonts'] = True

#%% Create dataframe with all features and targets
# Upload data
rawdata = pd.read_csv("data/ALKBH1_KD_stress/Isodecoder_counts_DESEqNormalized.csv",sep="\t", index_col=0)
rawdata.drop("size", axis=1, inplace=True)
rawdata.columns = [s.split(".")[0] for s in rawdata.columns]
rawdata.index = [s[13:] for s in rawdata.index]
pcadf = rawdata.T

#%% PCA
pca = PCA(n_components=2)
principalComponents = pca.fit_transform(StandardScaler().fit_transform(pcadf))
expl_var = pca.explained_variance_ratio_
features = pd.DataFrame(pca.components_.transpose(), columns=["PCA1","PCA2"], index = pcadf.columns)

principalDf = pd.DataFrame(principalComponents,columns=["PCA1","PCA2"],index=pcadf.index)

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
fig.savefig("plots/PCA_rpm_alkbh1_kdstress.pdf", bbox_inches='tight')

# Plot features
ordfeat = features.sort_values(by="PCA1", ascending=True)
ax = pd.concat([ordfeat.iloc[:10],ordfeat.iloc[-10:]]).plot.bar(y="PCA1")
ax.figure.tight_layout()
ax.figure.savefig("plots/PCA_rpm_features_alkbh1_kdstress.pdf")
