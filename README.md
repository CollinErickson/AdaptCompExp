
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
#> Starting iteration 1 at 2020-01-26 16:21:37 
#> no suitable  resolution IV or more  array found
#> Warning in DoE.base::oa.design(nruns = L^2, nfactors = D + 1, nlevels = L, :
#> resources were not sufficient for optimizing column selection
#> Below nug.min, setting nug to nug.min and reoptimizing #82387
```

<img src="man/figures/README-example exp1-1.png" width="100%" />

    #> Starting iteration 2 at 2020-01-26 16:21:39

<img src="man/figures/README-example exp1-2.png" width="100%" />

## Example of a comparison

``` r
ca1 <- compare.adaptR6$new(func=TestFunctions::gaussian1, D=2, L=3,
                           batches=2, reps=2,
                           n0=6, obj="desirability",
                           selection_method=c('max_des', 'SMED'),
                           des_func=c('des_func_relmax', 'des_func_relmax')
)
ca1$run_all(noplot = T)
#> Running 1, completed 0/4 Sun Jan 26 4:21:43 PM 2020
#> Starting iteration 1 at 2020-01-26 16:21:43 
#> no suitable  resolution IV or more  array found
#> Warning in DoE.base::oa.design(nruns = L^2, nfactors = D + 1, nlevels = L, :
#> resources were not sufficient for optimizing column selection
#> Starting iteration 2 at 2020-01-26 16:21:43
#> Warning in self$select_new_points_from_max_des_red(): uses_ALM, why was this a
#> browser spot?
#> Warning in int_werrors_red_func(add_points_indices = Xopts_inds_to_consider):
#> This was a browser spot too #92348

#> Warning in int_werrors_red_func(add_points_indices = Xopts_inds_to_consider):
#> This was a browser spot too #92348

#> Warning in int_werrors_red_func(add_points_indices = Xopts_inds_to_consider):
#> This was a browser spot too #92348
#> Running 2, completed 1/4 Sun Jan 26 4:21:44 PM 2020
#> Starting iteration 1 at 2020-01-26 16:21:44 
#> no suitable  resolution IV or more  array found
#> Warning in DoE.base::oa.design(nruns = L^2, nfactors = D + 1, nlevels = L, :
#> resources were not sufficient for optimizing column selection
#> Starting iteration 2 at 2020-01-26 16:21:44
#> Warning in self$select_new_points_from_max_des_red(): uses_ALM, why was this a
#> browser spot?

#> Warning in self$select_new_points_from_max_des_red(): This was a browser spot
#> too #92348

#> Warning in self$select_new_points_from_max_des_red(): This was a browser spot
#> too #92348

#> Warning in self$select_new_points_from_max_des_red(): This was a browser spot
#> too #92348
#> Running 3, completed 2/4 Sun Jan 26 4:21:44 PM 2020
#> Starting iteration 1 at 2020-01-26 16:21:45 
#> no suitable  resolution IV or more  array found
#> Warning in DoE.base::oa.design(nruns = L^2, nfactors = D + 1, nlevels = L, :
#> resources were not sufficient for optimizing column selection
#> Starting iteration 2 at 2020-01-26 16:21:45 
#> Running 4, completed 3/4 Sun Jan 26 4:21:45 PM 2020
#> Starting iteration 1 at 2020-01-26 16:21:45 
#> no suitable  resolution IV or more  array found
#> Warning in DoE.base::oa.design(nruns = L^2, nfactors = D + 1, nlevels = L, :
#> resources were not sufficient for optimizing column selection
#> Starting iteration 2 at 2020-01-26 16:21:45
#> Some values in $outdf are NaN, warning comes from making meanlogdf #92538
ca1$plot()
```

<img src="man/figures/README-compare adapt example 1-1.png" width="100%" /><img src="man/figures/README-compare adapt example 1-2.png" width="100%" /><img src="man/figures/README-compare adapt example 1-3.png" width="100%" />

    #> Warning: Removed 8 rows containing missing values (geom_path).
    #> Warning: Removed 4 rows containing missing values (geom_path).
    #> Warning: Removed 8 rows containing missing values (geom_point).

<img src="man/figures/README-compare adapt example 1-4.png" width="100%" />
