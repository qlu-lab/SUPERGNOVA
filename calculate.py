#!/usr/bin/python
from __future__ import division, print_function
import multiprocessing
from subprocess import call
import numpy as np
import pandas as pd
import numpy.linalg as linalg

import ld.ldscore as ld
import ld.parse as ps
from ldsc_thin import __filter_bim__


def nearest_Corr(input_mat):
    d, v = linalg.eigh(input_mat)
    A = (v * np.maximum(d, 0)).dot(v.T)
    A = (A + A.T) / 2
    multiplier = 1 / np.sqrt(np.diag(A))
    A = A * multiplier
    A = (A.T * multiplier).T
    return A


def calLocalCov(i, partition, geno_array, coords, bps, gwas_snps, n1, n2, h1, h2, m, pheno_corr, pheno_corr_var, queue):
    CHR = partition.iloc[i, 0]
    START = partition.iloc[i, 1]
    END = partition.iloc[i, 2]

    idx = np.logical_and(np.logical_and(gwas_snps['CHR']==CHR, bps <= END), bps >= START)
    m0 = np.sum(idx)
    if m0 < 120:
        queue.put(None)
        return
    
    tmp_coords = coords[idx]

    block_gwas_snps = gwas_snps[idx]
    max_dist = 1
    block_left = ld.getBlockLefts(tmp_coords, max_dist)

    lN, blockLD = geno_array.ldCorrVarBlocks(block_left, idx)
    meanLD = np.mean(lN)
    local_LD = nearest_Corr(blockLD)

    d, v = linalg.eigh(local_LD)
    order = d.argsort()[::-1]
    d = d[order]
    v = v[:,order]
    if np.sum(d>0) < 120:
        queue.put(None)
        return
    
    sub_d = d[d>0]
    sub_v = v[:,d>0]

    tz1 = np.dot(sub_v.T, block_gwas_snps['Z_x'])
    tz2 = np.dot(sub_v.T, block_gwas_snps['Z_y'])
    y = tz1 * tz2

    wh1 = h1 * m0 / m
    wh2 = h2 * m0 / m
    Localh1 = (np.mean(block_gwas_snps['Z_x'] ** 2) - 1) / meanLD * m0 / n1
    Localh2 = (np.mean(block_gwas_snps['Z_y'] ** 2) - 1) / meanLD * m0 / n2
    Localrho = 

    queue.put(df)

def _supergnova(bfile, partition, thread, gwas_snps, n1, n2, h1, h2, pheno_corr, pheno_corr_var):
    m = len(gwas_snps)

    snp_file, snp_obj = bfile+'.bim', ps.PlinkBIMFile
    ind_file, ind_obj = bfile+'.fam', ps.PlinkFAMFile
    array_file, array_obj = bfile+'.bed', ld.PlinkBEDFile

    # read bim/snp
    array_snps = snp_obj(snp_file)
    chr_bfile = list(set(array_snps.df['CHR']))
    tmp_partition = partition[partition.iloc[:,0].isin(chr_bfile)]
    tmp_gwas_snps = gwas_snps[gwas_snps.iloc[:,0].isin(chr_bfile)]
    blockN = len(tmp_partition)
    # snp list
    annot_matrix, annot_colnames, keep_snps = None, None, None
    n_annot = 1

    keep_snps = __filter_bim__(tmp_gwas_snps, array_snps)

    array_indivs = ind_obj(ind_file)
    n = len(array_indivs.IDList)
    keep_indivs = None

    ## reading genotype

    geno_array = array_obj(array_file, n, array_snps, keep_snps=keep_snps,
        keep_indivs=keep_indivs, mafMin=None)
    coords = np.array(array_snps.df['CM'])[geno_array.kept_snps]
    bps = np.array(array_snps.df['BP'])[geno_array.kept_snps]

    ## Calculating local genetic covariance

    results = []
    pool = multiprocessing.Pool(processes = thread)
    queue = multiprocessing.Manager().Queue()
    for i in range(blockN):
        pool.apply_async(calLocalCov, args=(i, tmp_partition, geno_array, coords, 
            bps, tmp_gwas_snps, n1, n2, h1, h2, m, pheno_corr, pheno_corr_var, 
            queue))
    pool.close()
    pool.join()
    queue.put('STOP')


    for g_cov in iter(queue.get, 'STOP'):
        results.append(g_cov)
    df = pd.concat(results, ignore_index=True)
    return df

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
        chrs = list(set(partition.iloc[:,0]))
        for i in range(len(chrs)):
            cur_bfile = bfile.replace('@', str(chrs[i]))
            all_dfs.append(_supergnova(cur_bfile, partition, thread, gwas_snps, n1, n2, h1, h2, pheno_corr, pheno_corr_var))
            print('Computed local genetic covariance for chromosome {}'.format(chrs[i]))
        df = pd.concat(all_dfs, ignore_index=True)
    else:
        df = _supergnova(bfile, partition, thread, gwas_snps, n1, n2, h1, h2, pheno_corr, pheno_corr_var)
    
    return df