
<!-- README.md is generated from README.Rmd. Please edit that file -->

# lglasso

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/lglasso)](https://CRAN.R-project.org/package=lglasso)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
<!-- badges: end -->

<div style="text-align: justify">

The goal of *lglasso* package is to estimate networks from longitudinal
high-dimensional data. It can be used to estimate either one-stage
models where a single network is underlying all the longitudinal
observations, or two-stage models where the networks before and after a
treatment (exposure) are different from each other. The one(two)-stage
model can further be classified to homogeneous model and heterogeneous
models. For details of the definitions of these models, please check the
refereces listed below.

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

\[1\] Zhou J et al. Identifying stationary microbial interaction
networks based on irregularly spaced longitudinal 16S rRNA gene
sequencing data. Front Microbiomes. 2024;3:1366948. doi:
10.3389/frmbi.2024.1366948. Epub 2024 Jun 2. PMID: 40687607; PMCID:
PMC12276884.
