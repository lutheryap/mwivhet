# Generate Integer Group Indices

Converts a grouping variable into a vector of contiguous integer indices
ranging from 1 to \\G\\, where \\G\\ is the number of unique groups.
This is a helper function used to facilitate matrix indexing for
group-level variance calculations.

## Usage

``` r
Getgroupindex(ds, group)
```

## Arguments

- ds:

  Data frame. The dataset containing the grouping variable.

- group:

  Name of the grouping variable (passed as an unquoted symbol) within
  `ds`. This variable defines the clusters (e.g., judges, examiners, or
  time periods).

## Value

An integer vector of length \\n\\, containing values in \\\\1, \dots,
G\\\\.

## Details

The function extracts the specified `group` column from the data frame
using non-standard evaluation. It then maps each unique value of the
group to an integer index corresponding to its order of appearance in
`unique(group)`.

This transformation is functionally equivalent to
`as.numeric(as.factor(group))`, but ensures explicit handling within the
provided data frame context.

## Examples

``` r
data <- data.frame(judge_id = c("JudgeA", "JudgeB", "JudgeA", "JudgeC"))
Getgroupindex(data, judge_id)
#> [1] 1 2 1 3
```
