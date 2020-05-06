#!/usr/bin/python
from __future__ import division
import multiprocessing
from subprocess import call
import numpy as np
import pandas as pd


def calculate(bfile, partition, thread, gwas_snps, n1, n2, h1, h2):
    if thread is None:
        thread = multiprocessing.cpu_count()
        print('{C} CPUs are detected. Using {C} threads in computation  ... '.format(C=str(thread)))
    else:
        cpuNum = multiprocessing.cpu_count()
        thread = min(thread, cpuNum)
        print('{C} CPUs are detected. Using {N} threads in computation  ... '.format(C=str(cpuNum), N=str(thread)))

    bed = None
    if '@' in partition:
        bedfiles = []
        for i range(1, 23):
    else:
        bed = 

    df = None
    if '@' in bfile:
        all_dfs = []
        for i in range(1, 23):
            cur_bfile = bfile.replace('@', str(i))
            all_dfs.append(_supergnova(cur_bfile, gwas_snps))
            print('Computed LD scores for chromosome {}'.format(i))
        df = pd.concat(all_dfs)
    else:
        df = _supergnova(bfile, gwas_snps)

    return df