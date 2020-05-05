# SUPERGNOVA

SUPERGNOVA (SUPER GeNetic cOVariance Analyzer) is a statistical framework to perform local genetic covariance analysis. SUPERGNOVA only needs GWAS summary data and a reference panel as input data. The preprint is available at biorxiv.

## Requirements

The software is developed and tested in Linux and Mac OS environments. The following softwares and packages are required:

1. **Python 2.7**
2. **numpy**
3. **scipy**
4. **pandas**
5. **sklearn**
6. **bitarray**

## Tutorial

Suppose you would like to calculate local genetic covariance between  Crohn's disease and ulcerative colitis. We'll need a few types of files:

- **Summary statistics files:** You can get your own GWAS summary statistics files for these two diseases [here](https://www.ibdgenetics.org). We assume that the files are in the standard format that ``ldsc`` understands. If not, make sure to run them through the included ``munge_sumstats.py`` file or use the one included in ``ldsc`` (see [here](https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation#reformatting-summary-statistics) for instructions).

- **Plink bfiles:** These are files .bed/.bim/.fam format. You can download some that we have prepared for you here. These files are from the 1000 Genomes Project, with rare variants (MAF < 5\%) filtered out.

- **Genome partition files**: These files should be in bed format. Please note that different population may have different genome partition. Here is an example dataset for European population.

More details about these supplied files can be found in here.

You may run the following command:

```
python2 supergnova.py data/CD.sumstats.gz data/UC.sumstats.gz \
--N1 27726 \
--N2 28738 \
--bfile data/bfiles/eur_chr@_SNPmaf5 \
--partition data/partition/eur_chr@.bed \
--out results.txt
```
## Explanation of Command-Line Arguments

- The first two arguments, `data/CD.sumstats.gz` and `data/UC.sumstats.gz`, denote the locations of the first and second summary statistics files. These files may be compressed using gzip, bz2, zip, xz, or not compressed at all. The program will infer the compression method if the files end with .gz, .bz2, .zip, xz, respectively. As previously mentioned, we assume that the files are in the standard format that `ldsc` understands.

- The `N1` and `N2` arguments (optional) denote the sample sizes of the summary statistics files. If they are not provided, they will be inferred from the summary statistics files.

- The `bfile` argument denotes the prefix of the `.bed/.bim/.fam` genotypic data files. Note the '@', which denotes a wildcard character that SUPERGNOVA will be replaced with 1-22. Alternatively, if you have one set of genotypic data files with 22 chromosome combined, you can just specify one bfile. We recommend you use the whole genome data as the reference panel although you may only be interested in the local genetic correlation in one specific region. 

- The `partition` argument denotes the genome partition file in bed format. Note the '@', which denotes a wildcard character that SUPERGNOVA will be replaced with 1-22. Alternatively, if you are only interested in some specific genomic regions, you can just specify one bed file that summarizes all the regions that you want to look into.

- The `out` flag denotes the file location for the results to be outputted to.

## Explanation of Output
The output will be a whitespace-delimited text file, with the rows corresponding to different annotations and the columns as such:

- `annot_name:` The name of the annotation.
- `rho:` The genetic covariance estimate.
- `rho_corrected:` The genetic covariance estimate with sample overlap correction.
- `se_rho:` The standard error of the estimate of `rho`.
- `pvalue:` The p-value from the statistical test for genetic covariance.
- `pvalue_corrected:` The p-value from the statistical test for genetic covariance with sample overlap correction.
- `corr`: The genetic correlation estimate.
- `corr_corrected`: The genetic correlation estimate with sample overlap correction.
- `h2_1`: The heritability estimate for the first trait.
- `h2_2`: The heritability estimate for the second trait.

NOTE: When functional annotations are present, the true heritability in each annotation category may be small. Although methods for estimating annotation-stratified heritability exist, they may provide unstable, in many cases negative heritability estimates, especially when a number of annotation categories are related to the repressed or non-functional genome. GNOVA ignores negative hertiability estimates, leaving the correlation estimates as 'NaN'. So, we recommend the users to focus on genetic covariance instead of genetic correlation when performing annotation-stratified analysis.


## Credits

Those using the SUPERGNOVA software should cite: Zhang, Y.L. et al. Local genetic correlation analysis reveals heterogeneous etiologic sharing of complex traits. 2020.

The LD score calculation  and the estimation of phenotypic covariance are adapted from `ldsc.py` in  `ldsc` and `ldsc_thin.py` in `GNOVA`. See [Bulik-Sullivan, B. et al. An Atlas of Genetic Correlations across Human Diseases and Traits. Nature Genetics, 2015.](https://www.nature.com/articles/ng.3406) and [Lu, Q.S. et al. A powerful approach to estimating annotation-stratified genetic covariance using GWAS summary statistics. The American Journal of Human Genetics, 2017.](https://www.cell.com/ajhg/fulltext/S0002-9297(17)30453-6)
