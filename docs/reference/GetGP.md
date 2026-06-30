# Construct UJIVE Weighting and Projection Matrices

Generates the weighting matrix \\G\\ and projection matrix \\P\\
required for the Unbiased Jackknife Instrumental Variables Estimator
(UJIVE). It constructs these matrices based on group/cluster indicators
for instruments (\\Z\\) and stratification covariates (\\W\\).

## Usage

``` r
GetGP(group, groupW, n)
```

## Arguments

- group:

  Integer vector. Primary grouping variable used to define the
  instrument structure (e.g., judge or examiner IDs).

- groupW:

  Integer vector. Grouping variable for stratification or covariates
  (e.g., time periods or court locations).

- n:

  Integer. The total sample size (number of observations).

## Value

A list containing four matrices:

- `G`: The \\N \times N\\ UJIVE weighting matrix.

- `P`: The \\N \times N\\ full projection matrix on \\\[Z, W\]\\.

- `Z`: The \\N \times K\\ matrix of instrument indicators.

- `W`: The \\N \times L\\ matrix of covariate indicators.

## Details

The function constructs the design matrices for instruments (\\Z\\) and
covariates (\\W\\) based on the provided grouping vectors. It creates
dummy variable matrices where \\Z\_{ij} = 1\\ if observation \\i\\
belongs to group \\j\\.

It calculates the standard projection matrices: \$\$P_Z =
Z(Z'Z)^{-1}Z'\$\$ \$\$P_W = W(W'W)^{-1}W'\$\$ \$\$P = P\_{\[Z,W\]} \quad
(\text{Projection onto both } Z \text{ and } W)\$\$

The UJIVE weighting matrix \\G\\ is then computed as the difference
between the "leave-one-out" adjusted projection of the full set and the
covariates: \$\$G = D_P^{-1}(P - \text{diag}(P)) - D_W^{-1}(P_W -
\text{diag}(P_W))\$\$ where \\D_P\\ and \\D_W\\ are diagonal matrices
containing the annihilator diagonals (\\1 - P\_{ii}\\) for the
respective projections.

Note: This function assumes that the resulting design matrices are full
rank.

## References

Yap, L. (2025). "Inference with Many Weak Instruments and
Heterogeneity". Working Paper.
