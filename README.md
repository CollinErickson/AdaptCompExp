
<!-- README.md is generated from README.Rmd. Please edit that file -->

# AdaptCompExp

<!-- badges: start -->

[![Travis build
status](https://travis-ci.org/CollinErickson/AdaptCompExp.svg?branch=master)](https://travis-ci.org/CollinErickson/AdaptCompExp)
[![Codecov test
coverage](https://codecov.io/gh/CollinErickson/AdaptCompExp/branch/master/graph/badge.svg)](https://codecov.io/gh/CollinErickson/AdaptCompExp?branch=master)
<!-- badges: end -->

The goal of AdaptCompExp is to implement algorithms for adaptive
computer experiments.

## Installation

You can install the released version of AdaptCompExp from
[CRAN](https://CRAN.R-project.org) with:

    # install.packages("AdaptCompExp")

And the development version from [GitHub](https://github.com/) with:

    # install.packages("devtools")
    devtools::install_github("CollinErickson/AdaptCompExp")

## Example

``` r
library(AdaptCompExp)
```

``` r
a <- adapt.concept2.sFFLHD.R6$new(D=2,L=3,func=TestFunctions::gaussian1,obj="desirability",
                                  des_func=des_func_relmax, n0=20, take_until_maxpvar_below=.9,
                                  package="GauPro", design='sFFLHD', selection_method="max_des_red",
                                  alpha_des=1)
#> Registered S3 method overwritten by 'DoE.base':
#>   method           from       
#>   factorize.factor conf.design
a$run(2)
#> Starting iteration 1 at 2020-01-25 16:30:21 
#> no suitable  resolution IV or more  array found
#> Warning in DoE.base::oa.design(nruns = L^2, nfactors = D + 1, nlevels = L, :
#> resources were not sufficient for optimizing column selection
#> Error in deviance_fngr_joint(X, Z, K) : chol(): decomposition failed
#> Below nug.min, setting nug to nug.min and reoptimizing #82387
#> Starting iteration 2 at 2020-01-25 16:30:24
#> Below nug.min, setting nug to nug.min and reoptimizing #82387
```

<img src="man/figures/README-example exp1-1.png" width="100%" /><img src="man/figures/README-example exp1-2.png" width="100%" />
