#!/usr/bin/python
from __future__ import division, print_function
import multiprocessing
from subprocess import call
import numpy as np
import pandas as pd
import numpy.linalg as linalg
from math import sqrt
import ld.ldscore as ld
import ld.parse as ps
from ldsc_thin import __filter_bim__
from scipy.stats import norm
from collections import OrderedDict


def nearest_Corr(input_mat):
    d, v = linalg.eigh(input_mat)
    A = (v * np.maximum(d, 0)).dot(v.T)
    A = (A + A.T) / 2
    multiplier = 1 / np.sqrt(np.diag(A))
    A = A * multiplier
    A = (A.T * multiplier).T
    return A


def calLocalCov(i, partition, geno_array, coords, bps, gwas_snps, ld_scores, n1, n2, pheno_corr, pheno_corr_var):
    m = len(gwas_snps)
    CHR = partition.iloc[i, 0]
    START = partition.iloc[i, 1]
    END = partition.iloc[i, 2]

    idx = np.logical_and(np.logical_and(gwas_snps['CHR']==CHR, bps <= END), bps >= START)
    m0 = np.sum(idx)
    if m0 < 120:
        df = pd.DataFrame(OrderedDict({"chr":[], "start":[], "end":[], "rho":[], "corr":[], "h2_1":[], "h2_2":[], "var":[], "p":[], "m":[]}))
        return df
    
    tmp_coords = coords[idx]

    block_gwas_snps = gwas_snps[idx]
    block_ld_scores = ld_scores[idx]
    max_dist = 0.03
    block_left = ld.getBlockLefts(tmp_coords, max_dist)

    lN, blockLD = geno_array.ldCorrVarBlocks(block_left, idx)
    lN = block_ld_scores["L2"]
    meanLD = np.mean(lN)
    local_LD = nearest_Corr(blockLD)

    d, v = linalg.eigh(local_LD)
    order = d.argsort()[::-1]
    d = d[order]
    v = v[:,order]
    if np.sum(d>0) < 120:
        df = pd.DataFrame(OrderedDict({"chr":[], "start":[], "end":[], "rho":[], "corr":[], "h2_1":[], "h2_2":[], "var":[], "p":[], "m":[]}))
        return df
    
    sub_d = d[d>0]
    sub_v = v[:,d>0]

    tz1 = np.dot(sub_v.T, block_gwas_snps['Z_x'])
    tz2 = np.dot(sub_v.T, block_gwas_snps['Z_y'])
    y = tz1 * tz2 - pheno_corr * sub_d

    Localh1 = (np.mean(block_gwas_snps['Z_x'] ** 2) - 1) / meanLD * m0 / n1
    Localh2 = (np.mean(block_gwas_snps['Z_y'] ** 2) - 1) / meanLD * m0 / n2

    Z_x = gwas_snps['Z_x']
    Z_y = gwas_snps['Z_y']

    h1 = (np.mean(Z_x ** 2) - 1) / np.mean(ld_scores['L2']) * m / n1
    h2 = (np.mean(Z_y ** 2) - 1) / np.mean(ld_scores['L2']) * m / n2

    wh1 = h1 * m0 / m
    wh2 = h2 * m0 / m
    #wh12 = np.max([Localh1, 0])
    #wh22 = np.max([Localh2, 0])
    #wh1 = (wh11 + wh12) / 2
    #wh2 = (wh21 + wh22) / 2
    Localrho = (np.sum(block_gwas_snps['Z_x'] * block_gwas_snps['Z_y']) - pheno_corr * m0) / meanLD / sqrt(n1 * n2)

    threshold = 1
    cur_d = sub_d[sub_d>threshold]
    cur_y = y[sub_d>threshold]
    cur_dsq = cur_d ** 2
    denominator = (wh1 * cur_d / m0 + 1 / n1) * (wh2 * cur_d / m0 + 1 / n2)
    cur_v1 = np.sum(cur_dsq / denominator)
    cur_v2 = np.sum(cur_y / sqrt(n1 * n2) / denominator)
    cur_v3 = np.sum(cur_y ** 2 / (n1 * n2) / (denominator * cur_dsq))

    emp_var = [(cur_v3 - (cur_v2 ** 2) / cur_v1) / (cur_v1 * (len(cur_d) - 1))]
    theo_var = [1 / cur_v1]

    for K in range(len(cur_d), len(sub_d)):
        eig = sub_d[K]
        tmp_y = y[K]
        cur_v1 += eig ** 2 / ((wh1 * eig / m0 + 1 / n1) * (wh2 * eig / m0 + 1 / n2))
        cur_v2 += tmp_y / sqrt(n1 * n2) / ((wh1 * eig / m0 + 1 / n1) * (wh2 * eig / m0 + 1 / n2))
        cur_v3 += tmp_y ** 2 / (n1 * n2) / ((wh1 * eig ** 2 / m0 + eig / n1) * (wh2 * eig ** 2 / m0 + eig / n2))
        emp_var.append((cur_v3 - (cur_v2 ** 2) / cur_v1) / (cur_v1 * K))
        theo_var.append(1 / cur_v1)
    
    max_emp_theo = np.maximum(emp_var, theo_var)
    min_idx = np.argmin(max_emp_theo)

    y = y[:(len(cur_d)+min_idx-1)]
    sub_d = sub_d[:(len(cur_d)+min_idx-1)]
    sub_dsq = sub_d ** 2

    var_rho = m0 ** 2 * min(max_emp_theo)
    q = (wh1 * sub_d / m0 + 1 / n1) * (wh2 * sub_d / m0 + 1 / n2)
    v4 = np.sum(sub_d/q)/np.sum(sub_dsq/q)
    var_phencorr = pheno_corr_var / (n1 * n2) * m0 ** 2 * v4 ** 2
    var_rho += var_phencorr

    se_rho = sqrt(var_rho)
    p_value = norm.sf(abs(Localrho / se_rho)) * 2

    if Localh1 < 0 or Localh2 < 0:
        corr = np.nan
    else:
        corr = Localrho / sqrt(Localh1 * Localh2)

    df = pd.DataFrame(OrderedDict({"chr":[CHR], "start":[START], "end":[END], "rho":[Localrho], "corr":[corr], "h2_1":[Localh1], "h2_2":[Localh2], "var":[var_rho], "p":[p_value], "m":[m0]}))

    return df

def _supergnova(bfile, partition, thread, gwas_snps, ld_scores, n1, n2, pheno_corr, pheno_corr_var):
    m = len(gwas_snps)

    snp_file, snp_obj = bfile+'.bim', ps.PlinkBIMFile
    ind_file, ind_obj = bfile+'.fam', ps.PlinkFAMFile
    array_file, array_obj = bfile+'.bed', ld.PlinkBEDFile

    # read bim/snp
    array_snps = snp_obj(snp_file)
    chr_bfile = list(set(array_snps.df['CHR']))
    tmp_partition = partition[partition.iloc[:,0].isin(chr_bfile)]
    tmp_gwas_snps = gwas_snps[gwas_snps.iloc[:,0].isin(chr_bfile)].reset_index(drop=True)
    tmp_ld_scores = ld_scores[ld_scores.iloc[:,0].isin(chr_bfile)].reset_index(drop=True)
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
    def collect_results(result):
        results.append(result)
    pool = multiprocessing.Pool(processes = thread)
    for i in range(blockN):
        pool.apply_async(calLocalCov, args=(i, tmp_partition, geno_array, coords, 
            bps, tmp_gwas_snps, tmp_ld_scores, n1, n2, pheno_corr, pheno_corr_var),
            callback=collect_results)
    pool.close()
    pool.join()
    df = pd.concat(results, ignore_index=True)
    #df = pd.DataFrame(results)
    #df.columns = ["chr", "start", "end", "rho", "corr", "h1", "h2", "var", "p", "m"]
    convert_dict = {"chr": int, "start": int, "end":int, "m":int}
    df = df.astype(convert_dict)
    return df

def calculate(bfile, partition, thread, gwas_snps, ld_scores, n1, n2, pheno_corr, pheno_corr_var):
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
            all_dfs.append(_supergnova(cur_bfile, partition, thread, gwas_snps, ld_scores, n1, n2, pheno_corr, pheno_corr_var))
            print('Computed local genetic covariance for chromosome {}'.format(chrs[i]))
        df = pd.concat(all_dfs, ignore_index=True)
    else:
        df = _supergnova(bfile, partition, thread, gwas_snps, ld_scores, n1, n2, pheno_corr, pheno_corr_var)
    
    return df
