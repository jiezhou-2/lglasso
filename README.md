
<!-- README.md is generated from README.Rmd. Please edit that file -->

# lglasso

<!-- badges: start -->
<!-- badges: end -->

The initial version of the package implemented the algorithms proposed in Zhou et al. (2021), 
which can estimate the high-dimensional  networks  from longitudinal data using Gassian grapical models. 
The current version includes several important extensions and improvements compared to the previous ones. 
The overarching goal however is the same, i.e., utilize the possible correlations between high-dimensional 
data to improve the estimations of networks.  
Specifically, currently the package implements three network models which correspond to three 
different way how the correlated data are modeled. Furthermore, the package can also estimate 
the individual networks for each cohorts which is particular useful in clinical studies. 

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

