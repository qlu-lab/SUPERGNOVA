#!/usr/bin/env python
import argparse, os.path, sys
from prep import prep
from ldsc_thin import ldscore
import pandas as pd
from calculate import calculate

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
    print('Performing eigen decomposition for LD matrix...')
    tilde, h1, h2, = ldeigen
    print('Calculating local genetic covariance...')
    out = calculate(gwas_snps, beds, N1, N2)
    out.to_csv(args.out, sep=' ', na_rep='NA', index=False)
