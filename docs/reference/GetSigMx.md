# Compute Covariance Matrix of UJIVE Quadratic Forms

Estimates the joint variance-covariance matrix of the three core
quadratic forms used in UJIVE/LIML estimation: \\Y'GY\\, \\X'GY\\, and
\\X'GX\\. This function is highly optimized for large-scale datasets by
leveraging block-diagonal geometries and scalar algebra.

## Usage

``` r
GetSigMx(df, groupW, group, X, Y, MX, MY, noisy = FALSE)
```

## Arguments

- df:

  Data frame. Contains the variables used in estimation.

- groupW:

  Column name (unquoted). The covariate stratification variable.

- group:

  Column name (unquoted). The instrument grouping variable.

- X:

  Column name (unquoted). The endogenous regressor.

- Y:

  Column name (unquoted). The outcome variable.

- MX:

  Column name (unquoted). Leverage-adjusted regressor (\\M X\\).

- MY:

  Column name (unquoted). Leverage-adjusted outcome (\\M Y\\).

- noisy:

  Logical. If `TRUE`, prints progress during variance component
  calculation. Defaults to `FALSE`.

## Value

Numeric vector of length 6. Contains
`c(sig11, sig22, sig33, sig12, sig23, sig13)`.

## Details

The function estimates the covariance components for the vector of
quadratic forms: \$\$\Psi = \[Y'GY, \quad X'GY, \quad X'GX\]^T\$\$

**Algorithmic Implementation:** To achieve high performance, the
function processes the data using a two-level nested loop over covariate
strata (`groupW`) and instrument groups (`group`). The heavy \\N \times
N\\ geometry matrices (such as the projection matrix \\P\\, the
leverage-adjusted weight matrix \\G\\, and the residual variance matrix
\\D_2\\) are pre-computed exactly once per stratum.

Within each instrument group, the function relies on localized
extraction and exact scalar algebra sub-routines (via \\A_1\\ and
\\A_4\\ component helpers) to evaluate the leave-three-out (L3O)
variance interactions. This approach minimizes computational complexity
to \\O(N)\\ at the group level and completely eliminates redundant
matrix allocations and inversions.

The returned vector contains the unique elements of the symmetric
covariance matrix \\\Sigma\_\Psi\\:

- `sig11`: \\Var(Y'GY)\\

- `sig22`: \\Var(X'GY)\\

- `sig33`: \\Var(X'GX)\\

- `sig12`: \\Cov(Y'GY, X'GY)\\

- `sig23`: \\Cov(X'GY, X'GX)\\

- `sig13`: \\Cov(Y'GY, X'GX)\\

## References

Yap, L. (2025). "Inference with Many Weak Instruments and
Heterogeneity". Working Paper.
