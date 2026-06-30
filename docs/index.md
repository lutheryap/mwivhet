# mwivhet: Inference with Many Weak Instruments and Heterogeneity

This R package implements inference in a linear instrumental variable
regression model with many potentially weak instruments, in the presence
of heterogeneous treatment effects, using procedures from [Yap
(2025)](https://www.arxiv.org/abs/2408.11193).

See the vignette `mwivhet` for a description of the package (available
through
[`vignette("mwivhet")`](https://lutheryap.github.io/mwivhet/articles/mwivhet.md)
once installed), and the [package
website](https://lutheryap.github.io/mwivhet) for full documentation.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
if (!requireNamespace("remotes")) {
  install.packages("remotes")
}
remotes::install_github("lutheryap/mwivhet")
```

## Citation

If you use mwivhet, please cite:
<https://doi.org/10.48550/arXiv.2408.11193>
