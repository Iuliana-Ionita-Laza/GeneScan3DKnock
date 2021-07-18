# GeneScan3DKnock 
This is an R package for performing improved gene-based testing by integrating long-range chromatin interactions and knockoff statistics.

## Description
Functions for the gene-based association tests that integrate both common and rare genetic variation from promoter and enhancers for each gene, along with the knockoff-enhanced tests. The GeneScan3DKnock has two steps: Step 1. Knockoff generation using function GeneScan3D.KnockoffGeneration() and Step 2. Knockoff filter using function GeneScan3DKnock().

## Workflow
![Workflow](https://user-images.githubusercontent.com/57265092/99107266-8c690a80-25b3-11eb-8fe1-ceb388bffa38.jpg)

## Prerequisites
R (recommended version >= 3.5.0)

## Dependencies
GeneScan3DKnock imports R packages SKAT, Matrix, MASS, WGScan, SPAtest, CompQuadForm, KnockoffScreen and abind. Make sure to install those packages before installing GeneScan3DKnock.
    
## Installation
library(devtools) 

devtools::install_github("Iuliana-Ionita-Laza/GeneScan3DKnock")

## Usage
Please see the GeneScan3DKnock <a href="https://github.com/Iuliana-Ionita-Laza/GeneScan3DKnock/blob/master/GeneScan3DKnock_0.2.pdf"> **user manual** </a> for detailed usage of GeneScan3DKnock package. Please see the <a href="https://htmlpreview.github.io/?https://github.com/Iuliana-Ionita-Laza/GeneScan3DKnock/blob/master/GeneScan3DKnock_vignette.html">**tutorial**</a> of using the GeneScan3DKnock package.


## Authors and Maintainer: 
Authors: Shiyang Ma, James Dalgleish, Zihuai He, Iuliana Ionita-Laza

Maintainer: Shiyang Ma <sm4857@cumc.columbia.edu>

## Version
The current version is 0.2 (March 3, 2021).

## Citation
Ma, S., Dalgleish, J. L ., Lee, J., Wang, C., Liu, L., Gill, R., Buxbaum, J. D., Chung, W., Aschard, H., Silverman, E. K., Cho, M. H., He, Z. and Ionita-Laza, I. "Powerful gene-based testing by integrating long-range chromatin interactions and knockoff genotypes". medRxiv 2021.07.14.21260405 (2021）doi: 10.1101/2021.07.14.21260405

## License
This software is licensed under GPL-3.
