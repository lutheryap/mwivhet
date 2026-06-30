# Compute Quadratic Roots for Confidence Sets

Computes the roots of the quadratic equation \\a\beta^2 + b\beta + c =
0\\. These roots serve as the boundaries for the confidence sets
constructed by inverting the UJIVE/LIML score test.

## Usage

``` r
GetCIvals(CIcoef)
```

## Arguments

- CIcoef:

  Numeric vector of length 3. The coefficients \\(a, b, c)\\ of the
  quadratic inequality, typically obtained from
  [`GetCIcoef`](https://lutheryap.github.io/mwivhet/reference/GetCIcoef.md).

## Value

Numeric vector of length 2. Contains the two roots.

## Details

This function applies the standard quadratic formula: \$\$\beta\_{1,2} =
\frac{-b \pm \sqrt{b^2 - 4ac}}{2a}\$\$

**Note:** This function does not check the sign of the discriminant. It
is intended to be called by a wrapper function (like
[`GetCItypebd`](https://lutheryap.github.io/mwivhet/reference/GetCItypebd.md))
that first verifies the existence of real roots (i.e., \\b^2 - 4ac \ge
0\\).

Depending on the sign of \\a\\ (convexity), these values represent
either the endpoints of a bounded confidence interval or the inner
boundaries of a disjoint ("donut") confidence set.

## References

Yap, L. (2025). "Inference with Many Weak Instruments and
Heterogeneity". Working Paper.
