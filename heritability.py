#!/usr/bin/python
from __future__ import division
import numpy as np
import pandas as pd


def heritability(gwas_snps, ld_scores, n1, n2):
    np.seterr(invalid='ignore')
    m = len(gwas_snps)
    ld_scores = ld_scores.drop(['CHR', 'BP', 'CM', 'MAF'], axis=1, errors='ignore').reset_index(drop=True)
    merged = pd.merge(gwas_snps, ld_scores, on=['SNP'])
    ld_scores = merged.iloc[:,4]
    Z_x, Z_y = merged['Z_x'], merged['Z_y']
    
    Z_xx = Z_x ** 2
    Z_yy = Z_y ** 2

    nblock = 200
    blocksize = int(np.floor(m / nblock)) + 1
    q1_block = np.empty(nblock)
    q2_block = np.empty(nblock)

    for i in range(nblock):
        ind_start = blocksize * i 
        ind_end = min(m-1, blocksize * (i + 1))
        l = ld_scores.drop(range(ind_start, ind_end))
        z_xx = Z_xx.drop(range(ind_start, ind_end))
        z_yy = Z_yy.drop(range(ind_start, ind_end))
        m0 = len(l)
        q1_block[i] = (np.mean(z_xx) - 1) / np.mean(l) * m0 / n1
        q2_block[i] = (np.mean(z_yy) - 1) / np.mean(l) * m0 / n2
    h1 = np.mean(q1_block)
    h2 = np.mean(q2_block)
    var1 = np.cov(q1_block, bias=True) * (nblock - 1)
    var2 = np.cov(q2_block, bias=True) * (nblock - 1)
    return h1, h2, var1, var2