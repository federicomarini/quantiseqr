# quantiseqr

<!-- badges: start -->
[![R-CMD-check](https://github.com/federicomarini/quantiseqr/workflows/R-CMD-check/badge.svg)](https://github.com/federicomarini/quantiseqr/actions)
<!-- badges: end -->

The goal of `quantiseqr` is to streamline the quantification of tumor immune contexture from RNA-seq data.
It does so using the TIL10 siganture, designed and validated to identify 10 different cell types.

## Installation

You can install the development version of `quantiseqr` from GitHub with:

``` r
library("remotes")
remotes::install_github("federicomarini/quantiseqr", 
                        dependencies = TRUE, build_vignettes = TRUE)
```

## Usage overview

You can find the rendered version of the documentation of `quantiseqr` at the project website <https://federicomarini.github.io/quantiseqr>, created with `pkgdown`.


## Development

If you encounter a bug, have usage questions, or want to share ideas and functionality to make this package better, feel free to file an [issue](https://github.com/federicomarini/quantiseqr/issues).

## Code of Conduct

Please note that the quantiseqr project is released with a [Contributor Code of Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html). By contributing to this project, you agree to abide by its terms.

