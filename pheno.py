#!/usr/bin/python
from __future__ import division
import numpy as np
import pandas as pd
from sklearn import linear_model


def pheno(gwas_snps, ld_scores, n1, n2, h1, h2):
    np.seterr(invalid='ignore')
    m = len(gwas_snps)
    ld_scores = ld_scores.drop(['CHR', 'BP', 'CM', 'MAF'], axis=1, errors='ignore').reset_index(drop=True)
    merged = pd.merge(gwas_snps, ld_scores, on=['SNP'])
    ld_scores = merged.iloc[:,4]
    Z_x, Z_y = merged['Z_x'], merged['Z_y']

    Z_xy = Z_x * Z_y

    Weight1 = 1 + n1 * h1 * ld_scores / m
    Weight2 = 1 + n2 * h2 * ld_scores / m
    Weight3 = np.mean(Z_xy) * ld_scores

    Weight = 1 / (Weight1 * Weight2 + Weight3^2)

    nblock = 200
    q_block = np.empty(nblock)
    for i, z_xy in enumerate(np.array_split(Z_xy, nblock))