# Compute UJIVE Signal Component (No Covariates)

Calculates the UJIVE "signal" or cross-product term \\X' G e\\ for the
grouped instrument design. This function exploits the block-diagonal
structure of the projection matrix implied by mutually exclusive
instrument groups to compute the quadratic form via efficient group-wise
summation.

## Usage

``` r
GetLM_nocov(df, X, e, groupZ, noisy = FALSE)
```

## Arguments

- df:

  Data frame. Contains the observable variables and grouping indicator.

- X:

  Column name (unquoted). The first variable (e.g., endogenous
  regressor).

- e:

  Column name (unquoted). The second variable (e.g., outcome or
  residual).

- groupZ:

  Column name (unquoted). The instrument grouping variable.

- noisy:

  Logical. If `TRUE`, prints progress of the group iteration. Defaults
  to `FALSE`.

## Value

Numeric scalar. The total sum of group-specific quadratic forms.

## Details

This function computes the scalar: \$\$S = \sum\_{g=1}^J e_g' \[
(I-D\_{P_g})^{-1}(P_g - D\_{P_g}) \] X_g\$\$ where \\P_g\\ is the
projection matrix onto the intercept for group \\g\\ (i.e., the group
mean).

This calculation corresponds to the numerator components (\\P\_{XY},
P\_{XX}\\) of the score statistic variance estimator in designs without
covariates.

## References

Yap, L. (2025). "Inference with Many Weak Instruments and
Heterogeneity". Working Paper.
