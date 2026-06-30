# Compute UJIVE Signal Component (Stratified Design)

Calculates the UJIVE "signal" or cross-product term \\X' G e\\ for a
design where instruments are nested within discrete covariate strata
(e.g., Judges within Years). This function iterates through covariate
groups to compute the quadratic form locally, handling the centering of
instruments within each block.

## Usage

``` r
GetLM(df, X, e, groupW, group, noisy = FALSE)
```

## Arguments

- df:

  Data frame. Contains the observable variables and grouping indicators.

- X:

  Column name (unquoted). The first variable (e.g., endogenous
  regressor).

- e:

  Column name (unquoted). The second variable (e.g., outcome or
  residual).

- groupW:

  Column name (unquoted). The covariate stratification variable (defines
  blocks).

- group:

  Column name (unquoted). The instrument grouping variable (defines
  treatments).

- noisy:

  Logical. If `TRUE`, prints progress of the stratum iteration. Defaults
  to `FALSE`.

## Value

Numeric scalar. The sum of stratum-specific quadratic forms.

## Details

This function implements the estimator for the signal component \\S\\ in
a stratified design: \$\$S = \sum\_{s} e_s' \[ U(P\_{Q,s}) - U(P\_{W,s})
\] X_s\$\$

Within each stratum \\s\\:

- \\P\_{Q,s}\\ is the projection onto the instrument groups.

- \\P\_{W,s}\\ is the projection onto the stratum intercept (local
  mean).

- \\U(P)\\ denotes the projection matrix with its diagonal elements
  removed, rescaled by the inverse of the annihilator diagonal \\(1 -
  P\_{ii})^{-1}\\.

This corresponds to the numerator terms (\\P\_{XY}, P\_{XX}\\) for test
inversion in designs with discrete controls.

**Algorithmic Implementation:** To ensure high performance and low
memory overhead on large datasets, this function computes the projection
transformation \\G_s X_s\\ using strictly \\O(N)\\ vector operations. It
computes group means via optimized aggregation, derives the diagonal
leverage adjustments inline, and directly applies the leave-one-out
transformation without ever constructing the dense \\N \times N\\
projection matrices.

## References

Yap, L. (2025). "Inference with Many Weak Instruments and
Heterogeneity". Working Paper.
