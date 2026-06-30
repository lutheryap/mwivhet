# Compute Quadratic Coefficients for CI (No Covariates)

Calculates the coefficients \\(a, b, c)\\ of the quadratic inequality
\\a\beta^2 + b\beta + c \leq 0\\ used to construct confidence intervals
for \\\beta\\ in the grouped design setting (no covariates). This
function inverts the UJIVE/LIML score test using optimized variance
estimators for block-diagonal designs.

## Usage

``` r
GetCIcoef_nocov(df, groupZ, X, Y, MX, MY, q = qnorm(0.975)^2, noisy = FALSE)
```

## Arguments

- df:

  Data frame. Contains the observable variables \\X, Y\\ and their
  projections.

- groupZ:

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

  Numeric scalar. Critical value for the test inversion (e.g.,
  \\1.96^2\\). Defaults to `qnorm(.975)^2` (approx. 3.84) for a 95
  percent confidence interval.

- noisy:

  Logical. If `TRUE`, prints progress dots during calculation. Defaults
  to `FALSE`.

## Value

Numeric vector of length 3 containing `c(a, b, c)`.

## Details

This function performs the same logic as
[`GetCIcoef`](https://lutheryap.github.io/mwivhet/reference/GetCIcoef.md)
but is specifically appropriate for designs where instruments form
mutually exclusive groups and no global covariates link them.

The returned coefficients define the confidence set: \$\$\\ \beta :
(P\_{XX}^2 - q C_2)\beta^2 + (-2 P\_{XY} P\_{XX} - q C_1)\beta +
(P\_{XY}^2 - q C_0) \leq 0 \\\$\$

## References

Yap, L. (2025). "Inference with Many Weak Instruments and
Heterogeneity". Working Paper.
