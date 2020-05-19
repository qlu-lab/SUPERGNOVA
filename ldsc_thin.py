#!/usr/bin/env python

from __future__ import division, print_function
import ld.ldscore as ld
import ld.parse as ps
import numpy as np
import pandas as pd


def _remove_dtype(x):
    '''Removes dtype: float64 and dtype: int64 from pandas printouts'''
    x = str(x)
    x = x.replace('\ndtype: int64', '')
    x = x.replace('\ndtype: float64', '')
    return x


def loj_bim(filter_df, array):
    r = filter_df.columns[1]
    l = array.IDList.columns[0]
    merge_df = filter_df.iloc[:,[1]]
    merge_df.loc[:,'keep'] = True
    z = pd.merge(array.IDList, merge_df, how='left', left_on=l, right_on=r, sort=False)
    ii = z['keep'] == True
    return ii.to_numpy().nonzero()[0]


def __filter_bim__(filter_df, array):
    merged_list = loj_bim(filter_df, array)
    len_merged_list = len(merged_list)
    if len_merged_list > 0:
        c = 'After merging, {0} SNPs remain'
    else:
        error_msg = 'No SNPs retained for analysis'
        raise ValueError(error_msg)
    return merged_list


def subset_annot_file(a_df, GWAS_df, kept_cols):
    GWAS_df.loc[:,'idx'] = pd.Series(range(len(GWAS_df.SNP.values)))
    a_df = pd.merge(a_df, GWAS_df, on=['SNP'])
    a_df = a_df.sort_values(['idx'])
    a_df.drop('idx', axis=1, inplace=True)
    a_df.rename(columns={'CHR_x':'CHR', 'BP_x':'BP', 'CM_x':'CM'}, inplace=True)
    a_df = a_df.iloc[:,0:kept_cols]
    return a_df


def remove_brackets(x):
    return x.replace('[', '').replace(']', '').strip()


def _ldscore(bfile, gwas_snps):
    '''
    Wrapper function for estimating l1, l1^2, l2 and l4 (+ optionally standard errors) from
    reference panel genotypes.

    Annot format is
    chr snp bp cm <annotations>

    '''

    snp_file, snp_obj = bfile+'.bim', ps.PlinkBIMFile
    ind_file, ind_obj = bfile+'.fam', ps.PlinkFAMFile
    array_file, array_obj = bfile+'.bed', ld.PlinkBEDFile
    # read bim/snp
    array_snps = snp_obj(snp_file)
    # snp list
    m = len(array_snps.IDList)
    annot_matrix, annot_colnames, keep_snps = None, None, None,
    n_annot = 1

    keep_snps = __filter_bim__(gwas_snps, array_snps)


    # read fam
    array_indivs = ind_obj(ind_file)
    n = len(array_indivs.IDList)
    # read keep_indivs
    keep_indivs = None

    # read genotype array
    geno_array = array_obj(array_file, n, array_snps, keep_snps=keep_snps,
        keep_indivs=keep_indivs, mafMin=None)

    #determine block widths

    max_dist = 1
    coords = np.array(array_snps.df['CM'])[geno_array.kept_snps]

    block_left = ld.getBlockLefts(coords, max_dist)

    scale_suffix = ''

    lN = geno_array.ldScoreVarBlocks(block_left, 50, annot=annot_matrix)
    col_prefix = "L2"
        
    ldscore_colnames = [col_prefix+scale_suffix]

    # print .ldscore. Output columns: CHR, BP, RS, [LD Scores]
    new_colnames = geno_array.colnames + ldscore_colnames
    df = pd.DataFrame.from_records(np.c_[geno_array.df, lN])
    df.columns = new_colnames
    df.drop(['CM','MAF'], axis=1)

    # print LD Score summary
    pd.set_option('display.max_rows', 200)
    t = df.iloc[:,4:].describe()

    np.seterr(divide='ignore', invalid='ignore')  # print NaN instead of weird errors
    # print correlation matrix including all LD Scores and sample MAF

    np.seterr(divide='raise', invalid='raise')
    return df


def ldscore(bfile, gwas_snps):
    df = None
    if '@' in bfile:
        all_dfs = []
        for i in range(1, 23):
            cur_bfile = bfile.replace('@', str(i))
            all_dfs.append(_ldscore(cur_bfile, gwas_snps))
            print('Computed LD scores for chromosome {}'.format(i))
        df = pd.concat(all_dfs)
    else:
        df = _ldscore(bfile, gwas_snps)

    numeric = df._get_numeric_data()
    numeric[numeric < 0] = 0
    return df
