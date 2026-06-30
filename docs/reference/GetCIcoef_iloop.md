# Compute Quadratic Coefficients for CI (General Design)

Computes the coefficients \\(a, b, c)\\ of the quadratic inequality
\\a\beta^2 + b\beta + c \leq 0\\ used to construct confidence intervals
for \\\beta\\ in general instrumental variable designs. This function
supports asymmetric weighting matrices (e.g., UJIVE with continuous
covariates) by using the fully generalized `_iloop` variance estimators.

## Usage

``` r
GetCIcoef_iloop(
  df,
  P,
  G,
  X,
  Y,
  MX,
  MY,
  Z,
  W,
  q = qnorm(0.975)^2,
  noisy = FALSE
)
```

## Arguments

- df:

  Data frame. Contains the observable variables \\X, Y\\ and their
  projections.

- P:

  Matrix of dimension n x n. The full projection matrix \\P\\.

- G:

  Matrix of dimension n x n. The UJIVE weighting matrix \\G\\.

- X:

  Column name (unquoted). The endogenous regressor.

- Y:

  Column name (unquoted). The outcome variable.

- MX:

  Column name (unquoted). Leverage-adjusted regressor (\\M X\\).

- MY:

  Column name (unquoted). Leverage-adjusted outcome (\\M Y\\).

- Z:

  Matrix of instruments.

- W:

  Matrix of covariates.

- q:

  Numeric scalar. Critical value for the test inversion (typically
  \\1.96^2\\). Defaults to `qnorm(.975)^2` (approx. 3.84) for a 95
  percent confidence interval.

- noisy:

  Logical. If `TRUE`, prints progress dots during calculation. Defaults
  to `FALSE`.

## Value

Numeric vector of length 3 containing `c(a, b, c)`.

## Details

This is the most general version of the confidence interval coefficient
calculator. It does not assume symmetry of the weighting matrix \\G\\.
Consequently, it computes all five distinct variance components (\\A_1\\
through \\A_5\\) for each term in the polynomial expansion of
\\\hat{V}(\beta)\\.

The returned coefficients define the curvature and position of the
confidence set parabola:

The function explicitly constructs the diagonal adjustments for the
UJIVE signal calculation using the matrices \\Z\\ and \\W\\.

## References

Yap, L. (2025). "Inference with Many Weak Instruments and
Heterogeneity". Working Paper.
