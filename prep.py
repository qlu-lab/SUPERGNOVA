import os

import numpy as np
import pandas as pd

def allign_alleles(df):
    """Look for reversed alleles and inverts the z-score for one of them.

    Here, we take advantage of numpy's vectorized functions for performance.
    """
    d = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    alleles = []
    for colname in ['A1_x', 'A2_x', 'A1_y', 'A2_y']:
        tmp = np.empty(len(df[colname]), dtype=int)
        for k, v in d.items():
            tmp[np.array(df[colname]) == k] = v
        alleles.append(tmp)
    reversed_alleles = np.logical_and(alleles[0] == alleles[3],
        alleles[1] == alleles[2])
    reversed_strand_flip_alleles = np.logical_and(alleles[0] == 3 - alleles[3],
        alleles[1] == 3 - alleles[2])
    to_flip = np.logical_or(reversed_alleles, reversed_strand_flip_alleles)
    df['Z_y'] *= -2 * to_flip + 1


def matched_or_reversed(df):
    """Returns boolean array signifying whether rows have matched or reversed
    alleles.
    """
    d = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    a = []  # array of alleles
    for colname in ['A1_x', 'A2_x', 'A1_y', 'A2_y']:
        tmp = np.empty(len(df[colname]), dtype=int)
        for k, v in d.items():
            tmp[np.array(df[colname]) == k] = v
        a.append(tmp)
    matched_alleles = (((a[0] == a[2]) & (a[1] == a[3])) |
        ((a[0] == 3 - a[2]) | (a[1] == 3 - a[3])))
    reversed_alleles = (((a[0] == a[3]) & (a[1] == a[2])) |
        ((a[0] == 3 - a[0]) | (a[1] == 3 - a[2])))
    return matched_alleles | reversed_alleles


def get_files(file_name):
    if '@' in file_name:
        valid_files = []
        for i in range(1, 23):
            cur_file = file_name.replace('@', str(i))
            if os.path.isfile(cur_file):
                valid_files.append(cur_file)
            else:
                raise ValueError('No file matching {} for chr {}'.format(
                    file_name, i))
        return valid_files
    else:
        if os.path.isfile(file_name):
            return [file_name]
        else:
            ValueError('No files matching {}'.format(file_name))


def prep(bfile, annot, sumstats1, sumstats2):
    bim_files = get_files(bfile + '.bim')

    if annot is not None:
        annot_files = get_files(annot)
        len_b, len_a = len(bim_files), len(annot_files)
        if len_b != len_a and len_b > 1 and len_a > 1:
            raise ValueError("The number of bim files and annotation files " +
                             "should be the same.")

    # read in bim files
    bims = [pd.read_csv(f,
                        header=None,
                        names=['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2'],
                        delim_whitespace=True) for f in bim_files]
    bim = pd.concat(bims, ignore_index=True)

    # read in annotation files
    if annot is not None:
        annots = [pd.read_csv(f, delim_whitespace=True) for f in annot_files]
        for i, a in enumerate(annots):
            if len(a) != len(bims[i]):
                raise ValueError("Number of rows in bim and annotation " +
                                 "files are not equal for {} and {}".format(
                                    bim_files[i], annot_files[i]))
            annots[i] = pd.concat([bims[i][['CHR', 'SNP', 'BP', 'CM']],
                                   pd.Series(np.ones(len(a))),
                                   a], axis=1)
            annots[i].rename(columns={0: 'ALL_'}, inplace=True)
    else:
        annots = None

    dfs = [pd.read_csv(file, delim_whitespace=True)
        for file in [sumstats1, sumstats2]]

    # rename cols
    bim.rename(columns={'A1': 'A1_ref', 'A2': 'A2_ref'}, inplace=True)
    dfs[0].rename(columns={'A1': 'A1_x', 'A2': 'A2_x', 'N': 'N_x', 'Z': 'Z_x'},
        inplace=True)
    dfs[1].rename(columns={'A1': 'A1_y', 'A2': 'A2_y', 'N': 'N_y', 'Z': 'Z_y'},
        inplace=True)

    # take overlap between output and ref genotype files
    df = pd.merge(bim, dfs[1], on=['SNP']).merge(dfs[0], on=['SNP'])

    # flip sign of z-score for allele reversals
    allign_alleles(df)
    df = df[matched_or_reversed(df)]

    return (df[['CHR', 'SNP', 'Z_x', 'Z_y']],
            dfs[0]['N_x'].max(),
            dfs[1]['N_y'].max(),
            annots)
