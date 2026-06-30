# Compute UJIVE Signal Component (General Design)

Calculates the UJIVE "signal" or cross-product term \\X' G Y\\ for a
general design with covariates. The weighting matrix \\G\\ is defined as
\\U(P_Q) - U(P_W)\\, where \\U(P)\\ represents the projection matrix
with diagonal elements set to zero and rows rescaled by
\\1/(1-P\_{ii})\\.

## Usage

``` r
GetLM_WQ(df, IdPW, IdPQ, dPW, dPQ, W, Q, X, Y)
```

## Arguments

- df:

  Data frame. Contains the observable variables.

- IdPW:

  Numeric vector. The annihilator diagonals \\1 - P\_{W,ii}\\.

- IdPQ:

  Numeric vector. The annihilator diagonals \\1 - P\_{Q,ii}\\.

- dPW:

  Numeric vector. The projection diagonals \\P\_{W,ii}\\.

- dPQ:

  Numeric vector. The projection diagonals \\P\_{Q,ii}\\.

- W:

  Matrix. The covariate matrix.

- Q:

  Matrix. The combined instrument and covariate matrix \\\[Z, W\]\\.

- X:

  Column name (unquoted). The first variable (e.g., regressor).

- Y:

  Column name (unquoted). The second variable (e.g., outcome).

## Value

Numeric scalar. The computed cross-product.

## Details

This function computes the scalar: \$\$S = X' \[ (I-D\_{P_Q})^{-1}(P_Q -
D\_{P_Q}) - (I-D\_{P_W})^{-1}(P_W - D\_{P_W}) \] Y\$\$

It uses an efficient matrix algebra expansion that avoids constructing
the \\N \times N\\ projection matrices explicitly. The computation
complexity is linear in \\N\\ (given pre-computed basis matrices),
making it suitable for large datasets where \\G\\ is dense.

The result corresponds to the cross-term \\P\_{XY}\\ (if \\X \neq Y\\)
or the self-term \\P\_{XX}\\ (if \\X = Y\\) used in the confidence
interval inequality.

## References

Yap, L. (2025). "Inference with Many Weak Instruments and
Heterogeneity". Working Paper.
