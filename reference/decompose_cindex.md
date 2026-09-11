# Decompose a concordance index into Event-Event and Event-Censored parts

Splits a concordance index into the two ranking tasks it pools: pairs
where both subjects had the event (`C_ee`), and pairs where one had the
event and one was censored (`C_ec`). As censoring rises the event-event
pairs vanish and `C_ec` dominates the global score, so a model close to
chance on true events can still report a stable global C-index.

The decomposition applies to any concordance index of the form
`C = sum(w_ij * c_ij) / sum(w_ij)`. The weighting is chosen with
[`weights_harrell()`](https://Lainsm.github.io/cindexdecomp/reference/cindex_weights.md),
[`weights_uno()`](https://Lainsm.github.io/cindexdecomp/reference/cindex_weights.md),
[`weights_truncated()`](https://Lainsm.github.io/cindexdecomp/reference/cindex_weights.md)
or
[`weights_custom()`](https://Lainsm.github.io/cindexdecomp/reference/cindex_weights.md),
and the identity `C = (W_ee * C_ee + W_ec * C_ec) / (W_ee + W_ec)` holds
for each.

## Usage

``` r
decompose_cindex(x, ...)

# Default S3 method
decompose_cindex(
  x,
  status,
  risk,
  weights = weights_harrell(),
  higher_is_riskier = TRUE,
  n_boot = 0,
  conf_level = 0.95,
  ...
)

# S3 method for class 'formula'
decompose_cindex(
  x,
  data = parent.frame(),
  weights = weights_harrell(),
  higher_is_riskier = TRUE,
  n_boot = 0,
  conf_level = 0.95,
  ...
)

# S3 method for class 'cindex_decomp'
print(x, ...)

# S3 method for class 'cindex_decomp'
summary(object, ...)

# S3 method for class 'cindex_decomp'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  A `cindex_decomp` object.

- ...:

  Ignored.

- status:

  Numeric vector of event indicators (1 = event, 0 = censored). Default
  method only.

- risk:

  Numeric vector of risk scores. Default method only.

- weights:

  A `cindex_weights` object. Defaults to
  [`weights_harrell()`](https://Lainsm.github.io/cindexdecomp/reference/cindex_weights.md).

- higher_is_riskier:

  Logical. `TRUE` (default) if larger values of `risk` indicate higher
  hazard. If `FALSE`, `risk` is negated once before any pair counting,
  so every component stays mutually consistent. Values below 0.5 are
  reported as they are: a worse-than-chance model is a finding, not an
  error.

- n_boot:

  Number of bootstrap replicates. `0` (default) skips resampling, in
  which case [`confint()`](https://rdrr.io/r/stats/confint.html) is
  unavailable. Replicates resample subjects, not pairs.

- conf_level:

  Confidence level stored on the object and used as the default for
  [`confint()`](https://rdrr.io/r/stats/confint.html).

- data:

  A data frame in which to evaluate the formula. Formula method only.

- object:

  A `cindex_decomp` object.

- row.names:

  Ignored.

- optional:

  Ignored.

## Value

An object of class `cindex_decomp`.

## Examples

``` r
set.seed(42)
time <- rexp(200, rate = 0.1)
status <- rbinom(200, 1, 0.6)
risk <- rnorm(200)

decompose_cindex(time, status, risk)
#> C-Index Decomposition
#> Weighting: Harrell | n = 200, events = 131 (65.5%), censoring 34.5%
#> 
#>                     C-index       pairs    share
#>   Event-Event        0.4855       8,515    67.3%
#>   Event-Censored     0.4545       4,132    32.7%
#>   ----------------------------------------------
#>   Global             0.4754      12,647   100.0%
#> 
#>   Masking gap (C_ec - C_ee): -0.031
#>   ! Event-Event concordance is near chance (0.50).
#>   ! Global C-index is below 0.50. If the risk score is oriented
#>     so that higher means lower hazard, set higher_is_riskier = FALSE.
#>   Pair counts describe composition, not precision. Set n_boot > 0 for intervals.
decompose_cindex(time, status, risk, weights = weights_uno())
#> C-Index Decomposition
#> Weighting: Uno (IPCW) | n = 200, events = 131 (65.5%), censoring 34.5%
#> 
#>                     C-index       pairs    share
#>   Event-Event        0.4956       8,515    67.3%
#>   Event-Censored     0.4520       4,132    32.7%
#>   ----------------------------------------------
#>   Global             0.4818      12,647   100.0%
#> 
#>   Masking gap (C_ec - C_ee): -0.044
#>   ! Event-Event concordance is near chance (0.50).
#>   ! Global C-index is below 0.50. If the risk score is oriented
#>     so that higher means lower hazard, set higher_is_riskier = FALSE.
#>   Pair counts describe composition, not precision. Set n_boot > 0 for intervals.
```
