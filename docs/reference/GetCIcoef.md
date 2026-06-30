# Compute Quadratic Coefficients for Confidence Sets (Grouped Data)

Estimates the coefficients \\a\\, \\b\\, and \\c\\ for the quadratic
inequality \\a\beta^2 + b\beta + c \le 0\\, which defines the
\\1-\alpha\\ confidence set for the structural parameter \\\beta\\. This
function is highly optimized for large-scale datasets, relying on
block-diagonal geometries and scalar algebra.

## Usage

``` r
GetCIcoef(df, groupW, group, X, Y, MX, MY, q = qnorm(0.975)^2, noisy = FALSE)
```

## Arguments

- df:

  Data frame. Contains the observable variables and their projections.

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

- q:

  Numeric scalar. Critical value for the test statistic inversion (e.g.,
  \\\chi^2\_{1, 1-\alpha}\\). Defaults to `qnorm(.975)^2` (approx. 3.84)
  for a 95 percent confidence interval.

- noisy:

  Logical. If `TRUE`, prints progress dots during calculation. Defaults
  to `FALSE`.

## Value

Numeric vector of length 3: `c(a, b, c)`.

## Details

The confidence set is constructed by inverting a test statistic based on
the quadratic form \\Q(\beta) = (\mathbf{Y} - \beta \mathbf{X})' G
(\mathbf{Y} - \beta \mathbf{X})\\. The coefficients are derived from the
variance estimator \\\hat{V}(\beta)\\ of this quadratic form, decomposed
into interactions between the outcome and the regressor.

**Algorithmic Implementation:** To achieve high performance, the
function processes the data using a two-level nested loop over covariate
strata (`groupW`) and instrument groups (`group`). The heavy \\N \times
N\\ geometry matrices (such as the projection matrices \\P\\, \\M\\, and
the leverage-adjusted weight matrix \\G\\) are pre-computed exactly once
per stratum.

Furthermore, the target inner products (\\P\_{XY}\\ and \\P\_{XX}\\) are
aggregated inline during the stratum loop to avoid redundant passes over
the data. Within each instrument group, unique Leave-Three-Out (L3O)
variance interactions (\\A_1\\ and \\A_4\\ terms) are evaluated using
exact scalar algebra helpers. This approach minimizes computational
complexity to \\O(N)\\ at the group level and completely eliminates
redundant matrix allocations.

The returned coefficients correspond to: \$\$a = P\_{XX}^2 - q \cdot
C_2\$\$ \$\$b = -2 P\_{XY} P\_{XX} - q \cdot C_1\$\$ \$\$c = P\_{XY}^2 -
q \cdot C_0\$\$

Where \\P\_{XY}\\ and \\P\_{XX}\\ are the UJIVE estimators for the
cross-products, and \\C_0, C_1, C_2\\ are the variance components
compiled via the L3O adjustment framework.

## References

Yap, L. (2025). "Inference with Many Weak Instruments and
Heterogeneity". Working Paper.
