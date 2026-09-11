# Compare the decomposition across several models

Runs
[`decompose_cindex()`](https://Lainsm.github.io/cindexdecomp/reference/decompose_cindex.md)
for each of several risk scores measured on the same cohort, and
assembles the results into one table with bootstrap standard deviations.
This is the input consumed by
[`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html) to
draw the dumbbell comparison.

## Usage

``` r
compare_decompositions(
  risks,
  time,
  status,
  weights = weights_harrell(),
  higher_is_riskier = TRUE,
  n_boot = 1000,
  conf_level = 0.95
)

# S3 method for class 'cindex_comparison'
print(x, ...)

# S3 method for class 'cindex_comparison'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- risks:

  A **named** list, data frame or matrix of risk-score vectors/columns,
  each the same length as `time` and `status`. The names (or, for a
  matrix, the column names) label the models.

- time, status:

  Numeric vectors describing the shared cohort.

- weights:

  A `cindex_weights` object. Defaults to
  [`weights_harrell()`](https://Lainsm.github.io/cindexdecomp/reference/cindex_weights.md).

- higher_is_riskier:

  Logical; see
  [`decompose_cindex()`](https://Lainsm.github.io/cindexdecomp/reference/decompose_cindex.md).
  Applied to every model.

- n_boot:

  Bootstrap replicates per model. Defaults to 1000, since the comparison
  table exists to carry uncertainty.

- conf_level:

  Confidence level stored on each fit, and used for the
  `ci_*_lo`/`ci_*_hi` percentile bounds below.

- x:

  A `cindex_comparison` object.

- ...:

  Ignored.

- row.names:

  Ignored.

- optional:

  Ignored.

## Value

An object of class `cindex_comparison`. Its `$table` element has one row
per model, with columns `model`, `ci_ee`, `ci_ec`, `global_c`, `gap`,
`sd_ee`, `sd_ec`, `sd_global`, `sd_gap`, `n_ee`, `n_ec`, `w_ee`, `w_ec`,
`ci_ee_lo`, `ci_ee_hi`, `ci_ec_lo` and `ci_ec_hi`. `w_ee`/ `w_ec` are
the weighted totals behind `ci_ee`/`ci_ec`, so
`global_c == (w_ee * ci_ee + w_ec * ci_ec) / (w_ee + w_ec)` is
verifiable for every row. `ci_ee_lo`/`ci_ee_hi` and `ci_ec_lo`/
`ci_ec_hi` are `conf_level` percentile bootstrap bounds (`NA` when
`n_boot = 0`); `sd_*` columns are kept alongside them for backward
compatibility.

## Examples

``` r
set.seed(42)
n <- 200
time <- rexp(n, rate = 0.1)
status <- rbinom(n, 1, 0.6)
risks <- list(model_a = rnorm(n), model_b = rnorm(n))
compare_decompositions(risks, time, status, n_boot = 0)
#> C-Index Decomposition Comparison
#> Weighting: Harrell | 2 model(s) | n_boot = 0
#> 
#>    model   C_ee   C_ec C_global      gap
#>  model_a 0.4855 0.4545   0.4754 -0.03099
#>  model_b 0.4740 0.5578   0.5014  0.08385
#> 
#> No bootstrap: sd_* columns are NA. Set n_boot > 0 for uncertainty.
```
