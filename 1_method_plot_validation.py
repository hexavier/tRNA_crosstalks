# -*- coding: utf-8 -*-
"""
Created on Thu Oct 28 13:31:56 2021

@author: xa_he
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
plt.rcParams['pdf.use14corefonts'] = True

# Load data
fulltab = pd.read_csv("data/Scer_mimseq/interdependences.tsv", sep="\t",
                      dtype={"canon_var1":"category","canon_var2":"category","var1":"category","var2":"category"})
samples = ["SRR12026331","SRR12026332","SRR12026333"]
gene = "Saccharomyces_cerevisiae_tRNA-Phe-GAA-2"
static_pos = "Charged"

# Subset data
idxsample = np.array([s in samples for s in fulltab["sample"]])
idxpos = np.array([any([p1==static_pos,p2==static_pos]) for p1,p2 in zip(fulltab["canon_var1"],fulltab["canon_var2"])])
subtab = fulltab.loc[np.logical_and(fulltab.ref==gene, idxsample & idxpos)]

def recover_conttab(textstr):
    columnsvals = [x for x in tempdf.loc[idx,"values"].split('\n')[0].split()]
    splitted = [x.split() for x in tempdf.loc[idx,"values"].split('\n')[1:-1]]
    # Reconstruct index
    indexvals = [(s[0],s[1]) if len(s)==3 else (splitted[n-1][0],s[0]) for n,s in enumerate(splitted)]
    myindex = pd.MultiIndex.from_tuples(indexvals, names=columnsvals)
    # Build dataframe
    outdf = pd.Series([s[-1] for s in splitted], index=myindex)
    return outdf
    
#%% Compute and plot relative numbers of modified vs unmodified reads
fig, ax = plt.subplots(1, figsize=(7,3))
nsamples = len(samples)
fullwidth = 0.40
width = fullwidth*0.85; space = fullwidth-width

for n,s in enumerate(samples):
    tempdf = subtab.loc[subtab["sample"]==s]
    ## CHARGED ##
    uniqpos = np.unique([p for s in tempdf[["canon_var1","canon_var2"]].to_numpy() for p in s if p!=static_pos])
    plotdf = pd.DataFrame(index = uniqpos, columns=["C_unmod","C_mod","U_unmod","U_mod"])
    for idx in tempdf.index:
        conttab = recover_conttab(tempdf.loc[idx,"values"])
        pos = tempdf.loc[idx,["canon_var1","canon_var2"]]; pos = pos.loc[pos!=static_pos][0]
        plotdf.loc[pos,"C_unmod"] = int(conttab.xs(("True", "False"), level=[static_pos,pos])[0])
        plotdf.loc[pos,"C_mod"] = int(conttab.xs(("True", "True"), level=[static_pos,pos])[0])
        plotdf.loc[pos,"U_unmod"] = int(conttab.xs(("False", "False"), level=[static_pos,pos])[0])
        plotdf.loc[pos,"U_mod"] = int(conttab.xs(("False", "True"), level=[static_pos,pos])[0])
    
    # Relative to total mod or unmod
    total_C = plotdf[["C_mod","C_unmod"]].sum(axis=1)
    plotdf[["C_mod","C_unmod"]] = plotdf[["C_mod","C_unmod"]].div(total_C,axis=0)
    total_U = plotdf[["U_mod","U_unmod"]].sum(axis=1)
    plotdf[["U_mod","U_unmod"]] = plotdf[["U_mod","U_unmod"]].div(total_U,axis=0)
    
    # Plot
    x = np.arange(plotdf.shape[0]); 
    ax.bar(x - fullwidth + fullwidth*(4*n+1)/(nsamples*2) + space/(2*nsamples), plotdf["C_mod"], width/nsamples, label= static_pos+'=1 + Mod', color= "#304D63", edgecolor="k", linewidth=0.5)
    ax.bar(x - fullwidth + fullwidth*(4*n+1)/(nsamples*2) + space/(2*nsamples), plotdf["C_unmod"], width/nsamples, label=static_pos+'=1 + Unmod', bottom=plotdf["C_mod"], color= "#8FB9AA", edgecolor="k", linewidth=0.5)
    ax.bar(x - fullwidth + fullwidth*(4*n+3)/(nsamples*2) - space/(2*nsamples), plotdf["U_mod"], width/nsamples, label=static_pos+'=0 + Mod', color= "#ED8975", edgecolor="k", linewidth=0.5)
    ax.bar(x - fullwidth + fullwidth*(4*n+3)/(nsamples*2) - space/(2*nsamples), plotdf["U_unmod"], width/nsamples, label=static_pos+'=0 + Unmod', bottom=plotdf["U_mod"], color= "#F2D096", edgecolor="k", linewidth=0.5)

# Plot
ax.set_xticks(x)
ax.set_xticklabels(plotdf.index)
ax.set(xlabel='Position', ylabel='% reads')
ax.legend(bbox_to_anchor=(1.04,1), loc="upper left")
ax.yaxis.grid()
ax.set_axisbelow(True)
# Legend
handles, labels = fig.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles))
ax.legend(by_label.values(), by_label.keys(),bbox_to_anchor=(1.04,1), loc="upper left")
fig.tight_layout()
#fig.savefig("PheGAA_58.pdf")

#%% Compute and plot relative numbers of charged vs uncharged reads
fig, ax = plt.subplots(1, figsize=(7,3))
nsamples = len(samples)
fullwidth = 0.40
width = fullwidth*0.85; space = fullwidth-width

for n,s in enumerate(samples):
    tempdf = subtab.loc[subtab["sample"]==s]
    ## CHARGED ##
    uniqpos = np.unique([p for s in tempdf[["canon_var1","canon_var2"]].to_numpy() for p in s if p!=static_pos])
    plotdf = pd.DataFrame(index = uniqpos, columns=["C_unmod","C_mod","U_unmod","U_mod"])
    for idx in tempdf.index:
        conttab = recover_conttab(tempdf.loc[idx,"values"])
        pos = tempdf.loc[idx,["canon_var1","canon_var2"]]; pos = pos.loc[pos!=static_pos][0]
        plotdf.loc[pos,"C_unmod"] = int(conttab.xs(("True", "False"), level=[static_pos,pos])[0])
        plotdf.loc[pos,"C_mod"] = int(conttab.xs(("True", "True"), level=[static_pos,pos])[0])
        plotdf.loc[pos,"U_unmod"] = int(conttab.xs(("False", "False"), level=[static_pos,pos])[0])
        plotdf.loc[pos,"U_mod"] = int(conttab.xs(("False", "True"), level=[static_pos,pos])[0])
    
    # Relative to total mod or unmod
    total_mod = plotdf[["C_mod","U_mod"]].sum(axis=1)
    plotdf[["C_mod","U_mod"]] = plotdf[["C_mod","U_mod"]].div(total_mod,axis=0)
    total_unmod = plotdf[["C_unmod","U_unmod"]].sum(axis=1)
    plotdf[["C_unmod","U_unmod"]] = plotdf[["C_unmod","U_unmod"]].div(total_unmod,axis=0)
    
    # Plot
    x = np.arange(plotdf.shape[0]); 
    ax.bar(x - fullwidth + fullwidth*(4*n+1)/(nsamples*2) + space/(2*nsamples), plotdf["C_unmod"], width/nsamples, label= static_pos+'=1 + Unmod', color= "#304D63", edgecolor="k", linewidth=0.5)
    ax.bar(x - fullwidth + fullwidth*(4*n+1)/(nsamples*2) + space/(2*nsamples), plotdf["U_unmod"], width/nsamples, label=static_pos+'=0 + Unmod', bottom=plotdf["C_unmod"], color= "#8FB9AA", edgecolor="k", linewidth=0.5)
    ax.bar(x - fullwidth + fullwidth*(4*n+3)/(nsamples*2) - space/(2*nsamples), plotdf["C_mod"], width/nsamples, label=static_pos+'=1 + Mod', color= "#ED8975", edgecolor="k", linewidth=0.5)
    ax.bar(x - fullwidth + fullwidth*(4*n+3)/(nsamples*2) - space/(2*nsamples), plotdf["U_mod"], width/nsamples, label=static_pos+'=0 + Mod', bottom=plotdf["C_mod"], color= "#F2D096", edgecolor="k", linewidth=0.5)

# Plot
ax.set_xticks(x)
ax.set_xticklabels(plotdf.index)
ax.set(xlabel='Position', ylabel='% reads')
ax.legend(bbox_to_anchor=(1.04,1), loc="upper left")
ax.yaxis.grid()
ax.set_axisbelow(True)
# Legend
handles, labels = fig.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles))
ax.legend(by_label.values(), by_label.keys(),bbox_to_anchor=(1.04,1), loc="upper left")
fig.tight_layout()
#fig.savefig("PheGAA_Charged.pdf")
