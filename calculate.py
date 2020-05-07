#!/usr/bin/python
from __future__ import division, print_function
import multiprocessing
from subprocess import call
import numpy as np
import pandas as pd

import ld.ldscore as ld
import ld.parse as ps
import ld.ldmatrix as lm
from ldsc_thin import __filter_bim__


def _supergnova(bfile, partition, thread, gwas_snps, n1, n2, h1, h2, pheno_corr, pheno_corr_var):
    m = len(gwas_snps)

    snp_file, snp_obj = bfile+'.bim', ps.PlinkBIMFile
    ind_file, ind_obj = bfile+'.fam', ps.PlinkFAMFile
    array_file, array_obj = bfile+'.bed', ld.PlinkBEDFile

    # read bim/snp
    array_snps = snp_obj(snp_file)
    chr_bfile = list(set(array_snps.df['CHR']))
    partition = partition[partition.ix[:,0].isin(chr_bfile)]
    gwas_snps = gwas_snps[gwas_snps.ix[:,0].isin(chr_bfile)]
    blockN = len(partition)
    # snp list
    annot_matrix, annot_colnames, keep_snps = None, None, None,
    n_annot = 1

    keep_snps = __filter_bim__(gwas_snps, array_snps)

    array_indivs = ind_obj(ind_file)
    n = len(array_indivs.IDList)
    keep_indivs = None

    ## reading genotype

    geno_array = array_obj(array_file, n, array_snps, keep_snps=keep_snps,
        keep_indivs=keep_indivs, mafMin=None)
    max_dist = 1
    coords = np.array(array_snps.df['CM'])[geno_array.kept_snps]

    pool = multiprocessing.Pool(processes = thread)
    queue = multiprocessing.Manager().Queue()
    for i in range(blockN):
        pool.apply_async(calBlockCorr, args=(i, queue))
    pool.close()
    pool.join()
    queue.put('STOP')


def calculate(bfile, partition, thread, gwas_snps, n1, n2, h1, h2, pheno_corr, pheno_corr_var):
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