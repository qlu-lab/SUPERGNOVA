# SUPERGNOVA

SUPERGNOVA (SUPER GeNetic cOVariance Analyzer) is a statistical framework to perform local genetic covariance analysis. SUPERGNOVA only needs GWAS summary data and a reference panel as input data. The preprint is available at biorxiv.

## Requirements

The software is developed and tested in Linux and Mac OS environments. The following softwares and packages are required:

**Python 2.7:**
1. numpy
2. scipy
3. pandas
4. sklearn
5. bitarray
	
**R:**
1. data.table
2. optparse

## Tutorial

Suppose you would like to calculate local genetic covariance between two phenotypes. We'll need a few types of files:

- **Summary statistics files:** We assume that the files are in the standard format that ``ldsc`` understands. If not, make sure to run them through the included ``munge_sumstats.py`` file or use the one included in ``ldsc`` (see [here](https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation#reformatting-summary-statistics) for instructions).

- **Plink bfiles:** These are files .bed/.bim/.fam format. You can download some that we have prepared for you here. These files are from the 1000 Genomes Project, with rare variants (MAF < 5\%) filtered out.

- **Genome partition files**: These files should be in bed format, one for each chromosome. Please note that different population may have different genome partition. Here is an example data for European population.

More details about these supplied files can be found in here.

## Credits

Those using the SUPERGNOVA software should cite: Zhang, Y.L. et al. Local genetic correlation analysis reveals heterogeneous etiologic sharing of complex traits. 2020.

The LD score calculation  and the estimation of phenotypic covariance are adapted from `ldsc.py` in  `ldsc` and `ldsc_thin.py` in `GNOVA`. See [Bulik-Sullivan, B. et al. An Atlas of Genetic Correlations across Human Diseases and Traits. Nature Genetics, 2015.](https://www.nature.com/articles/ng.3406) and [Lu, Q.S. et al. A powerful approach to estimating annotation-stratified genetic covariance using GWAS summary statistics. The American Journal of Human Genetics, 2017.](https://www.cell.com/ajhg/fulltext/S0002-9297(17)30453-6)
