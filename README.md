
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
#> Starting iteration 1 at 2020-03-01 20:39:43 
#> no suitable  resolution IV or more  array found
#> Warning in DoE.base::oa.design(nruns = L^2, nfactors = D + 1, nlevels = L, :
#> resources were not sufficient for optimizing column selection
#> Below nug.min, setting nug to nug.min and reoptimizing #82387
```

<img src="man/figures/README-example_exp1-1.png" width="100%" />

    #> Starting iteration 2 at 2020-03-01 20:39:46

<img src="man/figures/README-example_exp1-2.png" width="100%" />

## Example of a comparison

``` r
ca1 <- compare.adaptR6$new(func=TestFunctions::gaussian1, D=2, L=3,
                           batches=2, reps=2,
                           n0=6, obj="desirability",
                           selection_method=c('max_des', 'SMED'),
                           des_func=c('des_func_relmax', 'des_func_relmax')
)
#> [1] "TestFunctions::gaussian1"
#> [1] 2
#> [1] 3
#> [1] 3
#> [1] 2
#> [1] 2
#> [1] "desirability"
#> [1] 0
#> [1] 0
#> [1] 6
#> [1] "laGP"
#> [1] "max_des" "SMED"   
#> [1] "sFFLHD"
#> [1] "des_func_relmax" "des_func_relmax"
ca1$run_all(noplot = T)
#> Running 1, completed 0/4 Sun Mar 01 8:39:51 PM 2020
#> Starting iteration 1 at 2020-03-01 20:39:51 
#> no suitable  resolution IV or more  array found
#> Warning in DoE.base::oa.design(nruns = L^2, nfactors = D + 1, nlevels = L, :
#> resources were not sufficient for optimizing column selection
#> Starting iteration 2 at 2020-03-01 20:39:51 
#> Running 2, completed 1/4 Sun Mar 01 8:39:52 PM 2020
#> Starting iteration 1 at 2020-03-01 20:39:52 
#> no suitable  resolution IV or more  array found
#> Warning in DoE.base::oa.design(nruns = L^2, nfactors = D + 1, nlevels = L, :
#> resources were not sufficient for optimizing column selection
#> Starting iteration 2 at 2020-03-01 20:39:52 
#> Running 3, completed 2/4 Sun Mar 01 8:39:53 PM 2020
#> Starting iteration 1 at 2020-03-01 20:39:53 
#> no suitable  resolution IV or more  array found
#> Warning in DoE.base::oa.design(nruns = L^2, nfactors = D + 1, nlevels = L, :
#> resources were not sufficient for optimizing column selection
#> Starting iteration 2 at 2020-03-01 20:39:53 
#> Running 4, completed 3/4 Sun Mar 01 8:39:53 PM 2020
#> Starting iteration 1 at 2020-03-01 20:39:54 
#> no suitable  resolution IV or more  array found
#> Warning in DoE.base::oa.design(nruns = L^2, nfactors = D + 1, nlevels = L, :
#> resources were not sufficient for optimizing column selection
#> Starting iteration 2 at 2020-03-01 20:39:54
ca1$plot()
```

<img src="man/figures/README-compare_adapt_example_1-1.png" width="100%" /><img src="man/figures/README-compare_adapt_example_1-2.png" width="100%" /><img src="man/figures/README-compare_adapt_example_1-3.png" width="100%" />

    #> Warning: Removed 8 rows containing missing values (geom_path).
    #> Warning: Removed 4 rows containing missing values (geom_path).
    #> Warning: Removed 8 rows containing missing values (geom_point).

<img src="man/figures/README-compare_adapt_example_1-4.png" width="100%" />

## Code from papers

Below is some code that was used to create the results shown in papers.
