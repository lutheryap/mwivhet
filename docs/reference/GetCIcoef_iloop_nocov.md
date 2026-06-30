# Compute Quadratic Coefficients for CI (No Covariates, General P)

Calculates the coefficients \\(a, b, c)\\ for the confidence interval
quadratic inequality in the "No Covariates" setting (\\G = P\\), but for
**general symmetric projection matrices**. This function performs a
highly optimized single-pass loop to compute all polynomial coefficients
of the variance estimator \\\hat{V}(\beta)\\ simultaneously.

## Usage

``` r
GetCIcoef_iloop_nocov(X, Y, P, q = qnorm(0.975)^2, noisy = FALSE)
```

## Arguments

- X:

  Numeric vector of length n. The endogenous regressor.

- Y:

  Numeric vector of length n. The outcome variable.

- P:

  Matrix of dimension n x n. The symmetric projection matrix.

- q:

  Numeric scalar. The critical value for the test inversion (typically
  \\1.96^2\\). Defaults to `qnorm(0.975)^2`.

- noisy:

  Logical. If `TRUE`, prints progress through the N loops. Defaults to
  `FALSE`.

## Value

Numeric vector of length 3 containing `c(a, b, c)`.

## Details

This function is the solver for generic symmetric designs. It inverts
the test statistic: \$\$\frac{(P\_{XY} - \beta P\_{XX})^2}{C_2 \beta^2 +
C_1 \beta + C_0} \leq q\$\$

**Optimization:** Rather than calling variance estimators multiple
times, it decomposes the variance formula into geometric components
(depending only on \\P\\) and data components (\\X, Y\\). It iterates
through observations \\i\\ once, accumulating the weighted contributions
for \\C_0\\, \\C_1\\, and \\C_2\\ in parallel.

## References

Yap, L. (2025). "Inference with Many Weak Instruments and
Heterogeneity". Working Paper.
