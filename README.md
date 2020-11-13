# GeneScan3DKnock 
(A unified framework for gene-based testing with joint analysis of coding and regulatory variation, and integration of knockoff statistics for causal gene identification)

## Description
GeneScan3DKnock is an R package for performing the gene-based testing with joint analysis of coding and regulatory variation, i.e., promoter and putative enhancers. Functions for the gene-based association tests that integrate both common and rare genetic variation from putative regulatory elements, including promoters and enhancers for each gene, along with the knockoff-enhanced tests.

## Prerequisites
R (recommended version >= 3.5.0)

## Dependencies
GeneScan3DKnock imports R packages SKAT, Matrix, MASS, WGScan, SPAtest, CompQuadForm and KnockoffScreen. Make sure to install those packages before using GeneScan3DKnock.

## Installation
library(devtools) 

devtools::install_github("Iuliana-Ionita-Laza/GeneScan3DKnock")

## Usage
Please see the GeneScan3DKnock user manual for detailed usage of GeneScan3DKnock package.

## Version
The current version is 0.1 (November 12, 2020).

## Citation

Shiyang Ma, James Dalgleish, Justin Lee, Chen Wang, Richard Gill, Edwin K. Silverman, Michael Cho, Zihuai He, Iuliana Ionita-Laza (2020+) "GeneScan3DKnock: A unified framework for
gene-based testing with joint analysis of coding and regulatory variation, and integration of knockoff statistics for causal gene identification". 

## License
This software is licensed under GPLv3.
