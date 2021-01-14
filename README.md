# GeneScan3DKnock 
(A unified framework for gene-based testing with joint analysis of coding and regulatory variation, and integration of knockoff statistics for causal gene identification)

## Description
GeneScan3DKnock is an R package for performing the gene-based testing with joint analysis of coding and regulatory variation, i.e., promoter and putative enhancers. Functions for the gene-based association tests that integrate both common and rare genetic variation from putative regulatory elements for each gene, along with the knockoff-enhanced tests.

## Workflow
![Workflow](https://user-images.githubusercontent.com/57265092/99107266-8c690a80-25b3-11eb-8fe1-ceb388bffa38.jpg)

## Prerequisites
R (recommended version >= 3.5.0)

## Dependencies
GeneScan3DKnock imports R packages SKAT, Matrix, MASS, WGScan, SPAtest, CompQuadForm and KnockoffScreen. Make sure to install those packages before installing GeneScan3DKnock.

## Installation
library(devtools) 

devtools::install_github("Iuliana-Ionita-Laza/GeneScan3DKnock")

## Usage
Please see the GeneScan3DKnock user manual for detailed usage of GeneScan3DKnock package. https://github.com/Iuliana-Ionita-Laza/GeneScan3DKnock/blob/master/GeneScan3DKnock_0.1.pdf

## Version
The current version is 0.1 (November 13, 2020).

## Citation

Shiyang Ma, James Dalgleish, Justin Lee, Chen Wang, Linxi Liu, Richard Gill, Wendy Chung, Hugues Aschard, Edwin K. Silverman, Michael H. Cho, Zihuai He, Iuliana Ionita-Laza, "A unified knockoff framework for gene-based testing with joint analysis of coding and regulatory variation", 2020

## License
This software is licensed under GPLv3.
