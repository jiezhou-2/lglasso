
<!-- README.md is generated from README.Rmd. Please edit that file -->

# lglasso

<!-- badges: start -->
<!-- badges: end -->

The package lglasso  implements the algorithms proposed in Zhou et al. (2021), 
which aims to estimates the microbial interaction network (modeled by Gaussian graphical model) 
based on irregularly-spaced longitudinal abundance data. 
There are two models in Zhou et al. (2021)  which correspond two scenarios. 
The first model is homogeneous stationary Gaussian graphical model (SGGM) 
in which all the subjects share a same auto-correlation parameter $\tau$.
The second model is heterogeneous SGGM in which each subject possesses 
his/her own auto-correlation parameter $\tau_i$. 
For the details of the algorithms and the output of the packages, please see Zhou et al. (2021). The package also has been published on CRAN. 

## Installation

First, install the package remotes:

```
install.packages("remotes")
```

Then install lglasso :

```
remotes::install_github("jiezhou-2/lglasso", ref ="main") 
```

## How to use

 Please checkout the tutorial on the [package website](https://jiezhou-2.github.io/lglasso/).

