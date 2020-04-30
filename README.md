# SUPERGNOVA

### Introduction

SUPERGNOVA (SUPER-GeNetic cOVariance Analyzer) is a statistical framework to perform local genetic covariance analysis. SUPERGNOVA only needs GWAS summary data and a reference panel as input data. The preprint is available at biorxiv.

![SUPERGNOVA workflow](https://github.com/qlu-lab/SUPERGNOVA/blob/master/Fig/Figure1.png)

### Prerequisites

The software is developed and tested in Linux and Mac OS environments. The following softwares and packages are required:

- **R**:
	- **data.table**
- **Python 2.7**
	- **numpy**

### Tutorial

### Credits

Those using the SUPERGNOVA software should cite: Local genetic correlation analysis reveals heterogeneous etiologic sharing of complex traits.

The LD score calculation is adapted from `ldsc`. See [Bulik-Sullivan, et al. LD Score Regression Distinguishes Confounding from Polygenicity in Genome-Wide Association Studies.
Nature Genetics, 2015.](http://www.nature.com/ng/journal/vaop/ncurrent/full/ng.3211.html)