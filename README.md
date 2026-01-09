
<!-- README.md is generated from README.Rmd. Please edit that file -->

# lglasso

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/lglasso)](https://CRAN.R-project.org/package=lglasso)
<!-- badges: end -->

<div style="text-align: justify">

The previous version implemented the algorithms proposed in *Zhou et
al. (2024)*, which aims to estimate the high-dimensional networks from
longitudinal data using Gaussian graphical models. Though the
overarching goal of the package is the same, this updated version add
two more features, which are

1)  Joint estimation of multiple networks from longitudinal data. The
    models in previous version (i.e., in the 2024 paper) assumed a
    stationary process for the longitudinal data. This might not be the
    case if the treatment has multiple stages, e.g, the antibody network
    before and after vaccination. In this version, the function
    *lglasso* is extended to accommodate such scenarios.

2)  Cross validation. This version added function *CVlglasso* to select
    tuning parameter using cross validation method, where MSE of testing
    sets is minimized. Note we conduct the cross validation on the
    subject level instead of individual data level.

</div>

## Installation

You can install the development version of lglasso from
[GitHub](https://github.com/) with:

First, install the package remotes:

    install.packages("remotes")

Then install lglasso :

    remotes::install_github("jiezhou-2/lglasso", ref ="main") 

## How to use

Please see [package website](https://jiezhou-2.github.io/lglasso/).

**Reference**

\[1\] Zhou J, Gui J, Viles WD, Chen H, Li S, Madan JC, Coker MO, Hoen
AG. Identifying stationary microbial interaction networks based on
irregularly spaced longitudinal 16S rRNA gene sequencing data. Front
Microbiomes. 2024;3:1366948. doi: 10.3389/frmbi.2024.1366948. Epub 2024
Jun 2. PMID: 40687607; PMCID: PMC12276884.

\[2\] Friedman J., Hastie T., Tibshirani R. (2019) Graphical Lasso:
Estimation of Gaussian Graphical Models, Version: 1.11.

\[3\] Matt Galloway (2025), CVglasso: Lasso Penalized Precision Matrix
Estimation, version 1.0

\[4\] Danaher P, Wang P, Witten DM. The joint graphical lasso for
inverse covariance estimation across multiple classes. J R Stat Soc
Series B Stat Methodol. 2014 Mar;76(2):373-397. doi: 10.1111/rssb.12033.
PMID: 24817823; PMCID: PMC4012833.
