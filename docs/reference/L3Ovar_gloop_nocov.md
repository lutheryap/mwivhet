# Compute L3O Variance for Score Statistic (Grouped Design, No Covariates)

Computes the Leave-Three-Out (L3O) variance estimator
(\\\hat{V}\_{LM}\\) for the score statistic in a grouped instrument
design without covariates. This function leverages the block-diagonal
structure of the projection matrix to compute all variance components
locally within each group.

## Usage

``` r
L3Ovar_gloop_nocov(df, group, X, e, MX, Me)
```

## Arguments

- df:

  Data frame. Contains the variables used in estimation.

- group:

  Column name (unquoted). The instrument grouping variable.

- X:

  Column name (unquoted). The endogenous regressor.

- e:

  Column name (unquoted). Residuals under the null hypothesis (\\Y -
  X\beta_0\\).

- MX:

  Column name (unquoted). Leverage-adjusted regressor (\\M_i X_i\\).

- Me:

  Column name (unquoted). Leverage-adjusted residual (\\M_i e_i\\).

## Value

Scalar. The estimated variance \\\hat{V}\_{LM}\\.

## Details

This function implements the variance estimator \\\hat{V}\_{LM}\\ for
the specific case where \\G = P\\ (symmetric weights) and the design
matrix is block-diagonal (no covariates).

Instead of operating on \\N \times N\\ matrices, the function iterates
through groups. Within each group subset, the projection matrix is dense
but small. The function computes: \$\$\hat{V} = \sum_g (A\_{1g} +
2A\_{2g} + A\_{3g} - A\_{4g} - A\_{5g})\$\$

**Components:**

- **A1, A2, A3**: Signal, covariance, and error variance terms
  calculated using the symmetric weighting matrix \\G=P\\.

- **A4, A5**: Bias correction terms utilizing squared weights
  \\G\_{ij}^2\\.
