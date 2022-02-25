[![R-CMD-check](https://github.com/Iuliana-Ionita-Laza/GeneScan3DKnock/workflows/R-CMD-check/badge.svg)](https://github.com/Iuliana-Ionita-Laza/GeneScan3DKnock/actions)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
# GeneScan3DKnock 
This is an R package for performing improved gene-based testing by integrating long-range chromatin interactions and knockoff statistics.

## Description
The package contain functions for the gene-based association tests (GeneScan3D) that integrate both common and rare genetic variation from promoter and enhancers for each gene, along with the knockoff-enhanced tests. The GeneScan3DKnock has two steps: Step 1. Knockoff generation using function GeneScan3D.KnockoffGeneration() and Step 2. Knockoff filter using function GeneScan3DKnock(), after obtaining the original and knockoff p-values for each gene. 

To deal with the case-control imbalance issue for binary traits, we apply the Saddlepoint approximation (SPA) of the gene-based tests to avoid the inflation of Type I error rate. 

Besides, we also optimize the knockoff generations for GeneScan3DKnock using Shrinkage leveraging (SL) algorithm. After the optimization, knockoffs generation of whole-genome UK biobank genotypes only take 11 CPU hours for 1,000 computational cores.

## Workflow
<img src="https://github.com/Iuliana-Ionita-Laza/GeneScan3DKnock/blob/master/docs/Workflow.jpg" width="800" height="900">

## Prerequisites
R (recommended version >= 3.6.0)

## Dependencies
GeneScan3DKnock depends on R packages SKAT, Matrix, MASS, WGScan, SPAtest, CompQuadForm, abind and irlba. Make sure to install those packages before installing GeneScan3DKnock.
    
## Installation
library(devtools) 

devtools::install_github("Iuliana-Ionita-Laza/GeneScan3DKnock")

The current version is 0.3 (August 29, 2021).

## Usage
Please see the GeneScan3DKnock <a href="https://github.com/Iuliana-Ionita-Laza/GeneScan3DKnock/blob/master/docs/GeneScan3DKnock_0.3.pdf"> **user manual** </a> for detailed usage of GeneScan3DKnock package. Please see the <a href="https://htmlpreview.github.io/?https://github.com/Iuliana-Ionita-Laza/GeneScan3DKnock/blob/master/docs/GeneScan3DKnock_vignette.html">**tutorial**</a> of using the GeneScan3DKnock package.

## Contact
If you have any questions about GeneScan3DKnock please contact

- <sm4857@cumc.columbia.edu>

If you want to submit a issue concerning the software please do so using the **GeneScan3DKnock** [Github repository](https://github.com/Iuliana-Ionita-Laza/GeneScan3DKnock/issues).

## Version history
Version 0.3 (Resampling option in functions can be false. Don't need to conduct resampling-based moment-matching when the sample size is large, especially for UK biobank-scale data; Add W.threshold in the GeneScan3DKnock() function; Apply SPA gene-based test for binary traits, to deal with imbalance case-control issues; The knockoff generations are optimized using shrinkage leveraging algorithm.)

## Citation
* The GeneScan3DKnock paper: Ma, S., Dalgleish, J., Lee, J., Wang, C., Liu, L., Gill, R., Buxbaum, J. D., Chung, W. K., Aschard, H., Silverman, E. K., Cho, M. H., He, Z. and Ionita-Laza, I. (2021). ["Powerful gene-based testing by integrating long-range chromatin interactions and knockoff genotypes".](https://doi.org/10.1073/pnas.2105191118) _Proceedings of the National Academy of Sciences of the United States of America_, **118**, e2105191118.

* Shrinkage leveraging algorithm for knockoff generation: He, Z., Guen, Y. L., Liu, L., Lee, J., Ma, S., Yang, A. C.,  Liu. X., Rutledge, J., Losada, P. M., Song, B., Belloy, M. E., Butler, R. R., Longo, F. M., Tang, H., Mormino, E. C., Wyss-Coray, T., Greicius, M. D. and Ionita-Laza, I. (2021) ["Genome-wide analysis of common and rare variants via multiple knockoffs at biobank scale, with an application to Alzheimer disease genetics".](https://doi.org/10.1016/j.ajhg.2021.10.009) _American Journal of Human Genetics_, **108**, 2336-2353.


## License
This software is licensed under GPL-3.

## Functional annotation scores
For functional annotation scores, we use the genome-wide functional annotations in 127 different cell types and tissues (GenoNet scores, https://www.nature.com/articles/s41467-018-07349-w). The precomputed GenoNet scores for human genome assembly GRCh37 (hg19) can be downloaded here: https://zenodo.org/record/3336209#.YhkNI-iZOUk. 

The functional annotation scores can help to increase the power of the gene-based test, but this is not mandatory. If one don't want to use any functional score, just let Z=NULL, and the package would automatically use the Beta(MAF;1,25) for rare variants and Beta(MAF;1,1) for common variants as the weights.
