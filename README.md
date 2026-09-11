
<!-- README.md is generated from README.Rmd. Please edit that file -->

# lglasso

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/lglasso)](https://CRAN.R-project.org/package=lglasso)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
<!-- badges: end -->

<div style="text-align: justify">

The previous version aimed to estimate a one-stage high-dimensional
network from longitudinal data using Gaussian graphical models. This
version added two new features, which are

1)  Estimate two-stage high-dimensional networks. Previous version
    assumed a stationary process that underlies the longitudinal data.
    The new version extended this assumption to include the scenarios
    where the longitudinal data cover two stages, e.g., one for
    pre-treatment, the other for post-treatment. Each stage has its own
    network structure.

2)  Tuning parameter selection. This version added functions for the
    selection of tuning parameters. First, the likelihood value is added
    to the output of the main function *lglasso* so that users can used
    it to compute AIC or BIC for the model selection. Second, cross
    validation is added to the package which can be used to select the
    tuning parameter as well. CV becomes very slow when the network is
    too big. I personally recommend to use likelihood-based method to
    select the tuning parameter. It should be pointed out that the cross
    validation is performed on the subject level instead of individual
    data point level.

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
