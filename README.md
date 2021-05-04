# quantiseqr

<!-- badges: start -->
[![R-CMD-check](https://github.com/federicomarini/quantiseqr/workflows/R-CMD-check/badge.svg)](https://github.com/federicomarini/quantiseqr/actions)
<!-- badges: end -->

The goal of `quantiseqr` is to streamline the quantification of tumor immune contexture from RNA-seq data using quanTIseq. 
quanTIseq is deconvolution method based on an RNA-seq-derived signature matrix (i.e., *TIL10*) designed and validated to quantify 10 different immune cell types.
Please see [(Finotello F et al., Genome Medicine 2019)](https://doi.org/10.1186/s13073-019-0638-6) for additional details on quanTIseq.

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

## Citation

If you use this package in your work, please cite the original quanTIseq study:

Finotello F, Mayer C, Plattner C, Laschober G, Rieder D, Hackl H, Krogsdam A, Loncova Z, Posch W, Wilflingseder D, Sopper S, Ijsselsteijn M, Brouwer TP, Johnson D, Xu Y, Wang Y, Sanders ME, Estrada MV, Ericsson-Gonzalez P, Charoentong P, Balko J, de Miranda NFDCC, Trajanoski Z. Molecular and pharmacological modulators of the tumor immune contexture revealed by deconvolution of RNA-seq data. Genome Medicine, 2019. 11(1):34. https://doi.org/10.1186/s13073-019-0638-6

