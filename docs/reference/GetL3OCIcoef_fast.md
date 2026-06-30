# Compute Quadratic Coefficients for CI (Stratified Asymmetric)

Calculates the coefficients \\(a, b, c)\\ for the confidence interval
quadratic inequality in a general design that includes both grouping
(stratification) and covariates. This function is optimized for designs
where projection matrices are block-diagonal with respect to `groupW`
but potentially asymmetric within blocks due to covariate adjustments.

## Usage

``` r
GetL3OCIcoef_fast(
  df,
  groupW,
  groupQ,
  X,
  Y,
  MX,
  MY,
  q = qnorm(0.975)^2,
  noisy = FALSE
)
```

## Arguments

- df:

  Data frame. Contains the observable variables and grouping indicators.

- groupW:

  Column name (unquoted). The covariate stratification variable.

- groupQ:

  Column name (unquoted). The instrument grouping variable.

- X:

  Column name (unquoted). The endogenous regressor.

- Y:

  Column name (unquoted). The outcome variable.

- MX:

  Column name (unquoted). Leverage-adjusted regressor (\\M X\\).

- MY:

  Column name (unquoted). Leverage-adjusted outcome (\\M Y\\).

- q:

  Numeric scalar. The critical value for test inversion (typically
  \\1.96^2\\). Defaults to `qnorm(.975)^2`.

- noisy:

  Logical. If `TRUE`, prints progress. Defaults to `FALSE`.

## Value

Numeric vector of length 3 containing `c(a, b, c)`.

## Details

This function performs a loop:

1.  Splits data by `groupW`. Computes local projection matrices \\P_Q\\
    and \\P_W\\ for each stratum.

2.  Calculates the full UJIVE weighting matrix \\G = U(P_Q) - U(P_W)\\
    locally. Since \\G\\ is asymmetric, it computes all 5 variance
    components.

3.  Accumulates all polynomial coefficients for \\\hat{V}(\beta)\\
    simultaneously using extensive pre-calculation of vector products.

## References

Yap, L. (2025). "Inference with Many Weak Instruments and
Heterogeneity". Working Paper.
