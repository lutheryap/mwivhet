# Classify and Compute Confidence Interval Bounds

Solves the quadratic inequality \\a\beta^2 + b\beta + c \leq 0\\ derived
from the score test inversion to determine the topology and boundaries
of the confidence set for \\\beta\\.

## Usage

``` r
GetCItypebd(CIcoef)
```

## Arguments

- CIcoef:

  Numeric vector of length 3. The coefficients \\(a, b, c)\\ obtained
  from
  [`GetCIcoef`](https://lutheryap.github.io/mwivhet/reference/GetCIcoef.md).

## Value

Numeric vector of length 3. Format: `c(CItype, LowerBound, UpperBound)`.

## Details

The confidence set is defined as \\\\ \beta : a\beta^2 + b\beta + c \leq
0 \\\\. Depending on the sign of \\a\\ and the discriminant \\\Delta =
b^2 - 4ac\\, this set can take one of four forms:

- **Type 1: Bounded Interval** (\\a \ge 0, \Delta \ge 0\\). The parabola
  opens upward with real roots. The CI is the closed interval
  \\\[\beta\_{min}, \beta\_{max}\]\\.

- **Type 2: Disjoint Union (Donut)** (\\a \< 0, \Delta \ge 0\\). The
  parabola opens downward with real roots. The CI is the union of two
  infinite rays: \\(-\infty, \beta\_{min}\] \cup \[\beta\_{max},
  \infty)\\.

- **Type 3: Real Line** (\\a \< 0, \Delta \< 0\\). The parabola is
  always negative. The CI includes the entire real line. Bounds are
  returned as placeholders `c(-100, 100)`.

- **Type 4: Empty Set** (\\a \ge 0, \Delta \< 0\\). The parabola is
  always positive. The confidence set is empty, implying the model is
  rejected at the specified significance level for all \\\beta\\.
