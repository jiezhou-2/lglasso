
<!-- README.md is generated from README.Rmd. Please edit that file -->

# lglasso

<!-- badges: start -->
<!-- badges: end -->

<div style="text-align: justify">

The initial version of the package implemented the algorithms proposed
in Zhou et al. (2024), which aims to estimate the high-dimensional
networks from longitudinal data using Gassian grapical models. Though
the overarching goal of the package is the same as previous version,
i.e., explore the possible correlations between high-dimensional data to
improve the estimations of networks. this updated version does add
several essential extensions, which mainly include

1)  Estimate networks for given sub data sets. If we call the network
    estimated from the whole data set general network, then such sub
    data sets derived networks are called individal networks. They might
    correspond to the data for given tissue, or the data before
    vaccination etc. In these cases, a single general network is not
    enough to account for the difference between different tisses or
    befor and after vaccination;

2)  Compared to previous version, the current version can handle data
    from different tissues within one subject;

3)  This version provide more flexible modelings of the longitudinal
    data, and use more stable optimization algorithms. The package
    implements three network models which correspond to three different
    way how the correlated data are modeled.

4)  Besides the main function $lglasso$, this version also provide a
    function $cv.lglasso$ for selecting the optimal tuning parameter.

</div>

## Installation

You can install the development version of lglasso from
[GitHub](https://github.com/) with:

First, install the package remotes:

    install.packages("remotes")

Then install lglasso :

    remotes::install_github("jiezhou-2/lglasso", ref ="main") 

## How to use

Please se the vignette, [package
website](https://jiezhou-2.github.io/lglasso/).

\[1\] Zhou J, Gui J, Viles WD, Chen H, Li S, Madan JC, Coker MO, Hoen
AG. Identifying stationary microbial interaction networks based on
irregularly spaced longitudinal 16S rRNA gene sequencing data. Front
Microbiomes. 2024;3:1366948. doi: 10.3389/frmbi.2024.1366948. Epub 2024
Jun 2. PMID: 40687607; PMCID: PMC12276884.

\[2\] Friedman J., Hastie T., Tibshirani R. (2019) Graphical Lasso:
Estimation of Gaussian Graphical Models, Version: 1.11.
<https://CRAN.R-project.org/package=glasso>.
