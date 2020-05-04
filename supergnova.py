#!/usr/bin/env python
import argparse, os.path, sys
from subprocess import call
import multiprocessing
import pandas as pd
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

# returns whether the parent directory of path exists
def parent_dir_exists(path):
    return os.path.exists(os.path.abspath(os.path.join(path, os.pardir)))

def pipeline(args):
    pd.options.mode.chained_assignment = None

    # Sanity check args
    if not parent_dir_exists(args.out):
        raise ValueError('--out flag points to an invalid path.')

    print('Preparing files for analysis...')
    gwas_snps, N1, N2 = prep(args.bfile, args.sumstats1, args.sumstats2, args.N1, args.N2)
    print('Calculating LD scores...')
    ld_scores = ldscore(args.bfile, gwas_snps)
    print('Calculating heritability...')
    h_1, h_2 = heritability(gwas_snps, ld_scores, N1, N2)
    print('The genome-wide heritability of the first trait is {}.\nThe genome-wide heritability of the second trait is {}.'.format(h1, h2))
    print('Calculating phenotypic correlation...')
    pheno_corr, pheno_corr_var = pheno
    print('Calculating local genetic covariance...')
    out = calculate()
    out.to_csv(args.out, sep=' ', na_rep='NA', index=False)
