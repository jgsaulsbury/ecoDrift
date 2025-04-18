
<!-- README.md is generated from README.Rmd. Please edit that file -->

# turnoverNT

<!-- badges: start -->
<!-- badges: end -->

Contains functions for modeling change in fossil communities with
Hubbell’s neutral theory.

## Installation

You can install the development version of turnoverNT from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("jgsaulsbury/turnoverNT")
```

## Example

Simulate and visualize a timeseries under neutral theory with incomplete
sampling:

``` r
library(turnoverNT)
J <- 50000 #50,000 individuals
nsp <- 10 #10 species
tslength <- 5000 #run for 5000 timesteps
every <- 200 #sample every 200 timesteps...
ss <- 1000 #...and sample 1,000 individuals (with replacement) when you do
ages <- seq(0,tslength,every)
set.seed(1)
X <- simNT(startingabs=rep(J/nsp,nsp),ts=ages,ss=1000)
plot_spindles(X$simulation,X$times)
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />

Find the best-fit J for this simulated timeseries:

``` r
fitJ(occs=X$simulation,ages=X$times,CI=TRUE)
#> $loglik
#> [1] 513.0276
#> 
#> $J
#> [1] 42415.75
#> 
#> $CI
#> [1] 32343.11 55906.94
```
