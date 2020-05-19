#!/usr/bin/env python
'''
Local genetic correlation estimation

SUPERGNOVA

Created on 2020-5-4

Happy birthday PKU!

@author: Yiliang
'''

import argparse, os.path, sys
from subprocess import call
import multiprocessing
import pandas as pd
import numpy as np
from prep import prep
from ldsc_thin import ldscore
from heritability import heritability
from pheno import pheno
from calculate import calculate


try:
    x = pd.DataFrame({'A': [1, 2, 3]})
    x.drop_duplicates(subset='A')
except TypeError:
    raise ImportError('SUPERGNOVA requires pandas version > 0.15.2')


pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
pd.set_option('precision', 4)
pd.set_option('max_colwidth',1000)
np.set_printoptions(linewidth=1000)
np.set_printoptions(precision=4)


# returns whether the parent directory of path exists
def parent_dir_exists(path):
    return os.path.exists(os.path.abspath(os.path.join(path, os.pardir)))

def pipeline(args):
    pd.options.mode.chained_assignment = None

    # Sanity check args
    if not parent_dir_exists(args.out):
        raise ValueError('--out flag points to an invalid path.')

    print('Preparing files for analysis...')
    gwas_snps, bed, N1, N2 = prep(args.bfile, args.partition, args.sumstats1, args.sumstats2, args.N1, args.N2)
    print('Calculating LD scores...')
    ld_scores = ldscore(args.bfile, gwas_snps)
    gwas_snps = gwas_snps[gwas_snps['SNP'].isin(ld_scores['SNP'])]
    print('{} SNPs included in our analysis...'.format(len(gwas_snps)))
    print('Calculating heritability...')
    h_1, h_2 = heritability(gwas_snps, ld_scores, N1, N2)
    print('The genome-wide heritability of the first trait is {}.\nThe genome-wide heritability of the second trait is {}.'.format(h_1, h_2))
    print('Calculating phenotypic correlation...')
    pheno_corr, pheno_corr_var = pheno(gwas_snps, ld_scores, N1, N2, h_1, h_2)
    print('Calculating local genetic covariance...')
    out = calculate(args.bfile, bed, args.thread, gwas_snps, ld_scores, N1, N2, pheno_corr, pheno_corr_var)
    out.to_csv(args.out, sep=' ', na_rep='NA', index=False)


parser = argparse.ArgumentParser()

parser.add_argument('sumstats1',
    help='The first sumstats file.')
parser.add_argument('sumstats2',
    help='The second sumstats file.')

parser.add_argument('--bfile', required=True, type=str,
    help='Prefix for Plink .bed/.bim/.fam file.')
parser.add_argument('--partition', required=True, type=str,
    help='Genome partition file in bed format')
parser.add_argument('--N1', type=int,
    help='N of the sumstats1 file. If not provided, this value will be inferred '
    'from the sumstats1 arg.')
parser.add_argument('--N2', type=int,
    help='N of the sumstats2 file. If not provided, this value will be inferred '
    'from the sumstats2 arg.')

parser.add_argument('--out', required=True, type=str,
    help='Location to output results.')
parser.add_argument('--thread', default= multiprocessing.cpu_count(), type=int,
    help='Thread numbers used for calculation. Default = CPU numbers.')

if __name__ == '__main__':
    pipeline(parser.parse_args())

