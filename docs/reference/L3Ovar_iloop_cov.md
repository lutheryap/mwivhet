# Leave-Three-Out Variance Estimator for LM Statistic

Calculates the consistent variance estimator (\\\hat{V}\_{LM}\\) for the
Lagrange Multiplier (LM) test statistic using the "Leave-Three-Out"
(L3O) adjustment.

## Usage

``` r
L3Ovar_iloop_cov(X, e, P, G, noisy = FALSE)
```

## Arguments

- X:

  Numeric vector of length n. The endogenous variable.

- e:

  Numeric vector of length n. Residuals under the null hypothesis (\\e =
  Y - X\beta_0\\).

- P:

  Matrix of dimension n x n. Projection matrix of the instruments (and
  potentially covariates). Corresponds to matrix \\P\\ or \\H_Q\\ in the
  paper.

- G:

  Matrix of dimension n x n. Weighting matrix used in the JIVE/UJIVE
  estimator. For standard JIVE, G is equal to P. For UJIVE with
  covariates, G is the adjusted matrix defined in Section 3.1.

- noisy:

  Logical. If `TRUE`, print progress dots during the loop. Defaults to
  `FALSE`.

## Value

Scalar. The estimated variance \\\hat{V}\_{LM}\\.

## Details

This variance estimator is robust to both many weak instruments and
heterogeneous treatment effects. It corrects for biases in variance
estimation that arise when reduced-form coefficients are not
consistently estimable.

The function computes the variance estimator defined in Equation (9) of
Yap (2025): \$\$\hat{V}\_{LM} = A_1 + A_2 + A_3 + A_4 + A_5\$\$

It iterates through each observation \\i\\ to compute the necessary
adjustments (Leave-Three-Out determinants \\D\_{ijk}\\) and aggregates
the components using optimized matrix operations to handle the double
sums over \\j\\ and \\k\\.

Specifically:

- **A1, A2, A3** capture the core variance components involving
  interactions between the instruments, endogenous variable, and
  residuals.

- **A4, A5** are correction terms that account for the variability from
  estimating the reduced-form coefficients (which cannot be treated as
  fixed in the many-instrument setting).

The calculation relies on determinants \\D\_{ij}\\ and \\D\_{ijk}\\
derived from the annihilator matrix \\M = I - P\\ to ensure the
estimator is unbiased.

## References

Yap, L. (2025). "Inference with Many Weak Instruments and
Heterogeneity".
