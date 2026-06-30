# Compute L3O Variance for Score Statistic (Grouped Design)

Computes the Leave-Three-Out (L3O) variance estimator
(\\\hat{V}\_{LM}\\) for the score statistic in a grouped instrument
design with covariates. This function simultaneously calculates all five
variance components by exploiting the block-diagonal symmetry of the
projection matrices inherent to discrete instruments.

## Usage

``` r
L3Ovar_gloop_cov(df, group, groupW, X, e, MX, Me)
```

## Arguments

- df:

  Data frame. Contains the variables used in estimation.

- group:

  Column name (unquoted). The instrument grouping variable.

- groupW:

  Column name (unquoted). The covariate stratification variable.

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
the Limited Information Maximum Likelihood (LIML) or UJIVE estimator. It
computes: \$\$\hat{V} = A_1 + 2A_2 + A_3 - A_4 - A_5\$\$

Because the design implies symmetric weighting matrices (\\G\_{ij} =
G\_{ji}\\), the function simplifies the asymmetric components.

**Components:**

- **A1, A2, A3**: Variance and covariance of the score statistic.
  Calculated via 5 sub-components each.

- **A4, A5**: Bias correction terms for "own-observation" variance
  contributions. Calculated via 4 sub-components each.
