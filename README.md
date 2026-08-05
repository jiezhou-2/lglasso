
<!-- README.md is generated from README.Rmd. Please edit that file -->

# lglasso

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/lglasso)](https://CRAN.R-project.org/package=lglasso)
<!-- badges: end -->

<div style="text-align: justify">

The previous version aimed to estimate a high-dimensional network from
longitudinal data using Gaussian graphical models. This new update adds
two new features, which are

1)  Heterogeneous networks. Jointly estimate two networks from
    longitudinal data. Previous version assumed a stationary process
    that underlies the longitudinal data. The new version extended this
    assumption to include the cases where the longitudinal data cover
    two stages, one for pre-treatment, one for post-treatment. Each
    stage has its own correlation structure.

2)  Cross validation. This version added cross validation method to
    select tuning parameter. Note the cross validation is carried out on
    the subject level instead of individual data point level.

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
