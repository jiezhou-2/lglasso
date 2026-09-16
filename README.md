
<!-- README.md is generated from README.Rmd. Please edit that file -->

# lglasso

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/lglasso)](https://CRAN.R-project.org/package=lglasso)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
<!-- badges: end -->

<div style="text-align: justify">

R package *lglasso* is designed to estimate networks from longitudinal
high-dimensional data. It can be used to estimate either one-stage
models where a single network is underlying all the longitudinal
measurements, or two-stage models where the network before the treatment
is different from the one after the treatment (exposure). The
one(two)-stage model can further be classified to homogeneous model and
heterogeneous models. The difference between homogeneous and
heterogeneous models is that heterogeneous models contain individual
level random effects while homogeneous models treat all individuals as
i.i.d samples. For details of the definitions and its usage, please
check the reference and link below. If you have any questions, please
email *<chowstat@gmail.com>*. I will get back to you at my earliest
convenience.

</div>

## Installation

You can install the development version of lglasso in R from
[GitHub](https://github.com/) with:

First, install R package *remotes*:

    install.packages("remotes")

Then install *lglasso* :

    remotes::install_github("jiezhou-2/lglasso", ref ="main") 

## How to use

Here are some [examples](https://jiezhou-2.github.io/lglasso/) for the
usage of the package. If you are interested in the underlying
statistical model of *lglasso*, you could check the our manuscript here
[draft](https://drive.google.com/file/d/1oCNOGJODGfQITZJNyGDkGVMoaIgNU5on/view?usp=sharing),
which is an extension of the paper in the Reference.

**Reference**

\[1\] Zhou J et al. Identifying stationary microbial interaction
networks based on irregularly spaced longitudinal 16S rRNA gene
sequencing data. Front Microbiomes. 2024;3:1366948. doi:
10.3389/frmbi.2024.1366948. Epub 2024 Jun 2. PMID: 40687607; PMCID:
PMC12276884.
