# GeneScan3DKnock 
(A unified framework for gene-based testing with joint analysis of coding and regulatory variation, and integration of knockoff statistics for causal gene identification)

## Description
Functions for the gene-based association tests that integrate both common and rare genetic variation from proximal and distal regulatory elements, including promoter and enhancers for each gene, along with the knockoff-enhanced tests. The GeneScan3DKnock has two steps: Step 1. Knockoff generation using AR_KnockoffGeneration() function and Step 2. Knockoff filter using GeneScan3DKnock().

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
