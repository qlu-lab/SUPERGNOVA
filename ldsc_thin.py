#!/usr/bin/env python
'''
(c) 2014 Brendan Bulik-Sullivan and Hilary Finucane

LDSC is a command line tool for estimating
    1. LD Score
    2. heritability / partitioned heritability
    3. genetic covariance / correlation

'''
from __future__ import division, print_function
import ldscore.ldscore as ld
import ldscore.parse as ps
import numpy as np
import pandas as pd

try:
    x = pd.DataFrame({'A': [1, 2, 3]})
    x.drop_duplicates(subset='A')
except TypeError:
    raise ImportError('LDSC requires pandas version > 0.15.2')

__version__ = '1.0.0'
MASTHEAD = "*********************************************************************\n"
MASTHEAD += "* LD Score Regression (LDSC)\n"
MASTHEAD += "* Version {V}\n".format(V=__version__)
MASTHEAD += "* (C) 2014-2015 Brendan Bulik-Sullivan and Hilary Finucane\n"
MASTHEAD += "* Broad Institute of MIT and Harvard / MIT Department of Mathematics\n"
MASTHEAD += "* GNU General Public License v3\n"
MASTHEAD += "*********************************************************************\n"
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
pd.set_option('precision', 4)
pd.set_option('max_colwidth',1000)
np.set_printoptions(linewidth=1000)
np.set_printoptions(precision=4)


def sec_to_str(t):
    '''Convert seconds to days:hours:minutes:seconds'''
    [d, h, m, s, n] = reduce(lambda ll, b : divmod(ll[0], b) + ll[1:], [(t, 1), 60, 60, 24])
    f = ''
    if d > 0:
        f += '{D}d:'.format(D=d)
    if h > 0:
        f += '{H}h:'.format(H=h)
    if m > 0:
        f += '{M}m:'.format(M=m)

    f += '{S}s'.format(S=s)
    return f


class Logger(object):
    '''
    Lightweight logging.
    TODO: replace with logging module

    '''
    def __init__(self, fh):
        self.log_fh = open(fh, 'wb')

    def log(self, msg):
        '''
        Print to log file.
        TODO: Also print to stdout if verbose flag is set.

        '''
        print(msg, file=self.log_fh)


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
    return np.nonzero(ii)[0]


def __filter_bim__(filter_df, array, log):
    merged_list = loj_bim(filter_df, array)
    len_merged_list = len(merged_list)
    if len_merged_list > 0:
        c = 'After merging, {0} SNPs remain'
        log.log(c.format(len_merged_list))
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


def _ldscore(bfile, annots, gwas_snps):
    '''
    Wrapper function for estimating l1, l1^2, l2 and l4 (+ optionally standard errors) from
    reference panel genotypes.

    Annot format is
    chr snp bp cm <annotations>

    '''
    log = Logger('ldsc.log')

    snp_file, snp_obj = bfile+'.bim', ps.PlinkBIMFile
    ind_file, ind_obj = bfile+'.fam', ps.PlinkFAMFile
    array_file, array_obj = bfile+'.bed', ld.PlinkBEDFile
    # read bim/snp
    array_snps = snp_obj(snp_file)
    # snp list
    m = len(array_snps.IDList)
    log.log('Read list of {m} SNPs from {f}'.format(m=m, f=snp_file))
    if annots is not None:  # read --annot
        try:
            annot = ps.AnnotFile(pd.concat(annots, ignore_index=True))
            n_annot, ma = len(annot.df.columns) - 4, len(annot.df)
            log.log("Read {A} annotations for {M} SNPs".format(A=n_annot, M=ma))
            annot_colnames = annot.df.columns[4:]
            keep_snps = None
            #take only annot SNPs in intersect
            kept_cols = len(annot.df.columns)
            annot.df = subset_annot_file(annot.df, gwas_snps, kept_cols)
            annot_matrix = np.array(annot.df.iloc[:,4:])
        except Exception:
            log.log('Error parsing .annot file')
            raise
    else:
        annot_matrix, annot_colnames, keep_snps = None, None, None,
        n_annot = 1

    keep_snps = __filter_bim__(gwas_snps, array_snps, log)


    # read fam
    array_indivs = ind_obj(ind_file)
    n = len(array_indivs.IDList)
    log.log('Read list of {n} individuals from {f}'.format(n=n, f=ind_file))
    # read keep_indivs
    keep_indivs = None

    # read genotype array
    log.log('Reading genotypes from {fname}'.format(fname=array_file))
    geno_array = array_obj(array_file, n, array_snps, log, keep_snps=keep_snps,
        keep_indivs=keep_indivs, mafMin=None)

    #determine block widths

    max_dist = 1
    coords = np.array(array_snps.df['CM'])[geno_array.kept_snps]

    block_left = ld.getBlockLefts(coords, max_dist)

    scale_suffix = ''

    lN = geno_array.ldScoreVarBlocks(block_left, 50, annot=annot_matrix)
    col_prefix = "L2";

    if n_annot == 1:
        ldscore_colnames = [col_prefix+scale_suffix]
    else:
        ldscore_colnames =  [y+col_prefix+scale_suffix for y in annot_colnames]

    # print .ldscore. Output columns: CHR, BP, RS, [LD Scores]
    new_colnames = geno_array.colnames + ldscore_colnames
    df = pd.DataFrame.from_records(np.c_[geno_array.df, lN])
    df.columns = new_colnames
    df.drop(['CM','MAF'], axis=1)

    # print LD Score summary
    pd.set_option('display.max_rows', 200)
    log.log('\nSummary of LD Scores')
    t = df.ix[:,4:].describe()
    log.log( t.ix[1:,:] )

    np.seterr(divide='ignore', invalid='ignore')  # print NaN instead of weird errors
    # print correlation matrix including all LD Scores and sample MAF
    log.log('')
    log.log('MAF/LD Score Correlation Matrix')
    log.log( df.ix[:,4:].corr() )

    # print condition number
    if n_annot > 1: # condition number of a column vector w/ nonzero var is trivially one
        log.log('\nLD Score Matrix Condition Number')
        cond_num = np.linalg.cond(df.ix[:,5:])
        log.log(remove_brackets(str(np.matrix(cond_num))))
        if cond_num > 10000:
            log.log('WARNING: ill-conditioned LD Score Matrix!')

    # summarize annot matrix if there is one
    if annot_matrix is not None:
        # covariance matrix
        x = pd.DataFrame(annot_matrix, columns=annot_colnames)
        log.log('\nAnnotation Correlation Matrix')
        log.log( x.corr() )

        # column sums
        log.log('\nAnnotation Matrix Column Sums')
        log.log(_remove_dtype(x.sum(axis=0)))

        # row sums
        log.log('\nSummary of Annotation Matrix Row Sums')
        row_sums = x.sum(axis=1).describe()
        log.log(_remove_dtype(row_sums))

    np.seterr(divide='raise', invalid='raise')
    return df


def ldscore(bfile, annots, gwas_snps, save_ld):
    df = None
    if '@' in bfile:
        all_dfs = []
        for i in range(1, 23):
            cur_bfile = bfile.replace('@', str(i))
            if annots is None:
                cur_annot = None
            elif len(annots) > 1:
                cur_annot = [annots[i - 1]]
            else:
                cur_annot = annots
            all_dfs.append(_ldscore(cur_bfile, cur_annot, gwas_snps))
            print('Computed LD scores for chromosome {}'.format(i))
        df = pd.concat(all_dfs)
    else:
        df = _ldscore(bfile, annots, gwas_snps)

    numeric = df._get_numeric_data()
    numeric[numeric < 0] = 0
    if save_ld is not None:
        file_name = save_ld + '.csv.gz'
        print('Saving computed LD scores to {}.'.format(file_name))
        df.to_csv(file_name, ' ', index=False, compression='gzip')
    return df
