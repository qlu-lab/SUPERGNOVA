#!/usr/bin/python
from __future__ import division
import numpy as np
import pandas as pd


def heritability(gwas_snps, ld_scores, n1, n2):
    Z_x = gwas_snps['Z_x']
    Z_y = gwas_snps['Z_y']
    m = len(gwas_snps)
    h1 = (np.mean(Z_x ** 2) - 1) / np.mean(ld_scores['L2']) * m / n1
    h2 = (np.mean(Z_y ** 2) - 1) / np.mean(ld_scores['L2']) * m / n2
    return h1, h2