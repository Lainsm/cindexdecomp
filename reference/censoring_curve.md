# Trace the decomposition as administrative censoring increases

Sweeps a series of administrative cut-off times over the cohort. At each
cut-off, subjects whose event falls after the cut-off are treated as
censored there, and the decomposition is recomputed. This traces how the
global C-index behaves as censoring rises, while `C_ee` – the harder
ranking task – is tracked separately.

The cut-off is applied to the **whole cohort**, so the reported
censoring rate is the cohort's actual censoring rate at that cut-off.

## Usage

``` r
censoring_curve(
  time,
  status,
  risk,
  weights = weights_harrell(),
  higher_is_riskier = TRUE,
  n_thresholds = 20,
  probs = seq(0.05, 0.7, length.out = n_thresholds),
  min_pairs = NULL
)

# S3 method for class 'cindex_curve'
print(x, ...)

# S3 method for class 'cindex_curve'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- time, status, risk:

  Numeric vectors of equal length.

- weights:

  A `cindex_weights` object. Defaults to
  [`weights_harrell()`](https://lainsm.github.io/cindexdecomp/reference/cindex_weights.md).

- higher_is_riskier:

  Logical; see
  [`decompose_cindex()`](https://lainsm.github.io/cindexdecomp/reference/decompose_cindex.md).

- n_thresholds:

  Number of cut-offs to evaluate. Ignored if `probs` is supplied.

- probs:

  Quantiles of the observed event times at which to place the cut-offs.

- min_pairs:

  Minimum number of event-event pairs for an estimate to be considered
  precise. Rows below it are **flagged** in the `low_precision` column,
  never removed and never set to `NA`. `NULL` (default) selects
  `max(50, round(n_events^1.5 * 0.1))`.

- x:

  A `cindex_curve` object.

- ...:

  Ignored.

- row.names:

  Ignored.

- optional:

  Ignored.

## Value

An object of class `cindex_curve`. Its `$data` element has one row per
threshold, with columns `threshold`, `censoring`, `C_ee`, `C_ec`,
`C_global`, `gap`, `N_ee`, `N_ec`, `low_precision`, `W_ee` and `W_ec`.
`W_ee`/`W_ec` are the weighted totals behind `C_ee`/`C_ec`
(`C_ee = S_ee / W_ee`, etc.), so the decomposition identity is
verifiable at every row:
`C_global == (W_ee * C_ee + W_ec * C_ec) / (W_ee + W_ec)`. `N_ee`/
`N_ec` are pair *counts* and only equal `W_ee`/`W_ec` under a weighting
that assigns every comparable pair weight 1 (Harrell's); under
[`weights_uno()`](https://lainsm.github.io/cindexdecomp/reference/cindex_weights.md)
or another non-unit weighting the identity does not reconstruct from
counts alone.

## Examples

``` r
set.seed(42)
time <- rexp(300, rate = 0.1)
status <- rbinom(300, 1, 0.6)
risk <- rnorm(300)
censoring_curve(time, status, risk, n_thresholds = 8)
#> Censoring Curve
#> Weighting: Harrell | n = 300, events = 183 | 8 thresholds
#> Censoring range: 57.3% to 96.7%
#> 
#>  censoring   C_ee   C_ec C_global N_ee low_precision
#>      57.3% 0.5336 0.5415   0.5389 8128         FALSE
#>      63.0% 0.5276 0.5456   0.5408 6105         FALSE
#>      68.7% 0.5134 0.5522   0.5439 4371         FALSE
#>      74.3% 0.5465 0.5302   0.5328 2926         FALSE
#>      80.0% 0.5486 0.5383   0.5395 1770         FALSE
#>      85.7% 0.5249 0.5678   0.5644  903         FALSE
#>      91.0% 0.5413 0.5608   0.5598  351         FALSE
#>      96.7% 0.4667 0.5424   0.5413   45          TRUE
#> 
#> 1 row(s) flagged low_precision (N_ee < 248). Estimates are shown, not removed.
```
