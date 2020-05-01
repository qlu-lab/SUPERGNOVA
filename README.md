# SUPERGNOVA

SUPERGNOVA (SUPER-GeNetic cOVariance Analyzer) is a statistical framework to perform local genetic covariance analysis. SUPERGNOVA only needs GWAS summary data and a reference panel as input data. The preprint is available at biorxiv.

## Requirements

The software is developed and tested in Linux and Mac OS environments. The following softwares and packages are required:

- **Python 2.7**
	1. numpy
	2. scipy
	3. pandas
	4. sklearn
	5. bitarray
	
- **R**:
	- **data.table**

## Tutorial

Suppose you would like to calculate local genetic covariance between phenotypes. We'll need a few types of files:

- **Summary statistics files:** We assume that the files are in the standard format that ldsc understands. 

## Credits

Those using the SUPERGNOVA software should cite: Zhang, Y.L. et al. Local genetic correlation analysis reveals heterogeneous etiologic sharing of complex traits. 2020.

The LD score calculation  and the estimation of phenotypic covariance are adapted from `ldsc`. See [Bulik-Sullivan, B. et al. An Atlas of Genetic Correlations across Human Diseases and Traits. Nature Genetics, 2015.](https://www.nature.com/articles/ng.3406)