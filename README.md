[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
# GeneScan3DKnock 
This is an R package for performing improved gene-based testing by integrating long-range chromatin interactions and knockoff statistics.

## Description
The package contain functions for the gene-based association tests (GeneScan3D) that integrate both common and rare genetic variation from promoter and enhancers for each gene, along with the knockoff-enhanced tests. The GeneScan3DKnock has two steps: Step 1. Knockoff generation using function GeneScan3D.KnockoffGeneration() and Step 2. Knockoff filter using function GeneScan3DKnock(), after obtaining the original and knockoff p-values for each gene. 

To deal with the case-control imbalance issue for binary traits, we apply the Saddlepoint approximation (SPA) of the gene-based tests to avoid the inflation of Type I error rate. 

Besides, we also optimize the knockoff generations for GeneScan3DKnock using Shrinkage leveraging (SL) algorithm. After the optimization, knockoffs generation of whole-genome UK biobank genotypes only take 11 CPU hours for 1,000 computational cores.

## Workflow
![Workflow](https://user-images.githubusercontent.com/57265092/99107266-8c690a80-25b3-11eb-8fe1-ceb388bffa38.jpg)

## Prerequisites
R (recommended version >= 3.6.0)

## Dependencies
GeneScan3DKnock depends on R packages SKAT, Matrix, MASS, WGScan, SPAtest, CompQuadForm, abind and irlba. Make sure to install those packages before installing GeneScan3DKnock.
    
## Installation
library(devtools) 

devtools::install_github("Iuliana-Ionita-Laza/GeneScan3DKnock")

## Usage
Please see the GeneScan3DKnock <a href="https://github.com/Iuliana-Ionita-Laza/GeneScan3DKnock/blob/master/GeneScan3DKnock_0.3.pdf"> **user manual** </a> for detailed usage of GeneScan3DKnock package. Please see the <a href="https://htmlpreview.github.io/?https://github.com/Iuliana-Ionita-Laza/GeneScan3DKnock/blob/master/vignettes/GeneScan3DKnock_vignette.html">**tutorial**</a> of using the GeneScan3DKnock package.

## Contact
If you have any questions about GeneScan3DKnock please contact

- <sm4857@cumc.columbia.edu>

If you want to submit a issue concerning the software please do so
using the **GeneScan3DKnock** [Github repository](https://github.com/Iuliana-Ionita-Laza/GeneScan3DKnock/issues).

## Version history
The current version is Version 0.3.

Version 0.3 (Resampling option in functions can be false. Don't need to conduct resampling-based moment-matching when the sample size is large, especially for UK biobank-scale data; Add W.threshold in the GeneScan3DKnock() function; Apply SPA gene-based test for binary traits, to deal with imbalance case-control issues; The knockoff generations are optimized using shrinkage leveraging algorithm.)

## Citation
The GeneScan3DKnock pre-print: Ma, S., Dalgleish, J. L ., Lee, J., Wang, C., Liu, L., Gill, R., Buxbaum, J. D., Chung, W., Aschard, H., Silverman, E. K., Cho, M. H., He, Z. and Ionita-Laza, I. "Powerful gene-based testing by integrating long-range chromatin interactions and knockoff genotypes". medRxiv 2021.07.14.21260405 (2021). doi: 10.1101/2021.07.14.21260405

The KnockoffScreen-AL (shrinkage algorithmic leveraging) pre-print: He, Z.\*, Guen, Y. L.\*, Liu, L., Lee, J., Ma, S., Yang, A. C., Liu. X., Rutledge, J., Losada, P. M., Song, B., Butler, R., Longo, F., Tang, H., Mormino, E., Wyss-Coray, T., Greicius, M. and Ionita-Laza, I. Genome-wide analysis of common and rare variants via multiple knockoffs at biobank-scale: methods and applications (2021).

## License
This software is licensed under GPL-3.
