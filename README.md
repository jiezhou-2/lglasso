
<!-- README.md is generated from README.Rmd. Please edit that file -->

# lglasso

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/lglasso)](https://CRAN.R-project.org/package=lglasso)
<!-- badges: end -->

<div style="text-align: justify">

The previous version aimed to estimate the high-dimensional networks
from longitudinal data using Gaussian graphical models. This new update
adds two new features, which are

1)  Heterogeneous networks. Joint estimation of multiple networks from
    longitudinal data. The models in previous version assumed a
    stationary process for the longitudinal data. This might not be the
    case if the treatment has multiple stages in which each stage has
    its own correlation structure, or two cohorts which have their own
    correlation structure. In this version, the function *lglasso* is
    extended to accommodate these scenarios.

2)  Cross validation. This version added function *CVlglasso* to select
    tuning parameter using cross validation method, where MSE of testing
    sets is minimized. Note we conduct the cross validation on the
    subject level instead of individual data point level.

</div>

## Installation

You can install the development version of lglasso from
[GitHub](https://github.com/) with:

First, install the package remotes:

    install.packages("remotes")

Then install lglasso :

    remotes::install_github("jiezhou-2/lglasso", ref ="main") 

## How to use

Please click the following link for details [package
website](https://jiezhou-2.github.io/lglasso/).

**Reference**

\[1\] Friedman J et al (2019) Graphical Lasso: Estimation of Gaussian
Graphical Models, Version: 1.11.

\[2\] Danaher P et al. The joint graphical lasso for inverse covariance
estimation across multiple classes. J R Stat Soc Series B Stat Methodol.
2014 Mar;76(2):373-397. doi: 10.1111/rssb.12033. PMID: 24817823; PMCID:
PMC4012833.

\[3\] Zhou J et al. Identifying stationary microbial interaction
networks based on irregularly spaced longitudinal 16S rRNA gene
sequencing data. Front Microbiomes. 2024;3:1366948. doi:
10.3389/frmbi.2024.1366948. Epub 2024 Jun 2. PMID: 40687607; PMCID:
PMC12276884.
