# -*- coding: utf-8 -*-
"""
Created on Tue Feb  1 14:58:32 2022

@author: xa_he
"""
import numpy as np
import os
import gzip

# Inputs (Scer Phe-GAA-2)
seq = "GCGGATTTAGCTCAGTTGGGAGAGCGCCAGACTGAAGATCTGGAGGtCCTGTGTTCGATCCACAGAATTCGCACCA"
mod1 = 26
mod2 = 58
reads_per_sample = 20000

# Iteration variables
p_mod1 = np.arange(0.075,1.0,0.1)
p_mod2 = np.arange(0.075,1.0,0.1)
odds = np.power(2,[-1.0,-0.5,-0.25,0.0,0.25,0.5,1.0])

os.mkdir("simulated_seqs")
for p1 in p_mod1:
    for p2 in p_mod2:
        for o in odds:
            fqout = gzip.open("simulated_seqs/modA%s_modB%s_odds%s.fq.gz" % (str("%.3f" % p1).replace(".","-"),str("%.3f" % p2).replace(".","-"),str("%.3f" % o).replace(".","-")),"wt")
            for n in range(reads_per_sample):
                tempseq = str(seq)
                # Compute joint probabilities based on marginals and OR
                if o!=1:
                    p1p2joint = (1+(p1+p2)*(o-1)-np.sqrt((1+(p1+p2)*(o-1))**2 + 4*o*(1-o)*p1*p2))/(2*(o-1))
                else:
                    p1p2joint = p1*p2
                # Add first modification and compute conditional prob of the second modification
                if np.random.choice([True,False],p=[p1,1.0-p1]):
                    tempseq = tempseq[:mod1-1] + "T" + tempseq[mod1:]
                    p2temp = p1p2joint/p1
                else:
                    p2temp = (p2 - p1p2joint)/(1-p1)
                # Add second modification
                if np.random.choice([True,False],p=[p2temp,1.0-p2temp]):
                    tempseq = tempseq[:mod2-1] + "G" + tempseq[mod2:]
                fqout.write(str("@Read%i\n" % n))
                fqout.write(str("%s\n" % tempseq))
                fqout.write("+\n"+"I"*len(tempseq)+"\n")
            fqout.close()