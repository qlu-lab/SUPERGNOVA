#!/usr/bin/python
from __future__ import division, print_function
import multiprocessing
from subprocess import call
import numpy as np
import pandas as pd

import ld.ldscore as ld
import ld.parse as ps
import ld.ldmatrix as lm


def _supergnova(bfile, partition, thread, gwas_snps, n1, n2, h1, h2):
    snp_file, snp_obj = bfile+'.bim', ps.PlinkBIMFile
    ind_file, ind_obj = bfile+'.fam', ps.PlinkFAMFile
    array_file, array_obj = bfile+'.bed', ld.PlinkBEDFile




def calculate(bfile, partition, thread, gwas_snps, n1, n2, h1, h2):
    if thread is None:
        thread = multiprocessing.cpu_count()
        print('{C} CPUs are detected. Using {C} threads in computation  ... '.format(C=str(thread)))
    else:
        cpuNum = multiprocessing.cpu_count()
        thread = min(thread, cpuNum)
        print('{C} CPUs are detected. Using {N} threads in computation  ... '.format(C=str(cpuNum), N=str(thread)))

    df = None
    if '@' in bfile:
        all_dfs = []
        chrs = list(set(partition.ix[:,0]))
        for i in range(len(chrs)):
            cur_bfile = bfile.replace('@', str(chrs[i]))
            all_dfs.append(_supergnova(cur_bfile, partition, thread, gwas_snps, n1, n2, h1, h2))
            print('Computed local genetic covariance for chromosome {}'.format(chrs[i]))
        df = pd.concat(all_dfs)
    else:
        df = _supergnova(bfile, partition, thread, gwas_snps, n1, n2, h1, h2)
    
    return df