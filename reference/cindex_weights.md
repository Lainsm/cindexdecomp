# Pair weightings for concordance decomposition

A concordance index of the form `C = sum(w_ij * c_ij) / sum(w_ij)` is
fully determined by its pair weighting `w_ij`. These constructors supply
that weighting to
[`decompose_cindex()`](https://Lainsm.github.io/cindexdecomp/reference/decompose_cindex.md).
`weights_harrell()` gives Harrell's C, `weights_uno()` gives Uno's
inverse-probability-of-censoring-weighted C, and `weights_custom()`
accepts any user-supplied rule.

## Usage

``` r
weights_harrell()

weights_uno(tau = NULL)

weights_truncated(tau)

weights_custom(fn, name = "Custom")
```

## Arguments

- tau:

  Truncation time. For `weights_uno()`, pairs whose earlier event time
  exceeds `tau` receive weight zero; `NULL` (default) means no
  truncation. For `weights_truncated()`, `tau` is required.

- fn:

  A function of two arguments `(t, G)` returning a single non-negative
  finite number. `t` is the event time of the earlier member of the
  pair; `G` is a function returning the Kaplan-Meier censoring survival
  probability at a given time.

- name:

  A length-one character label used in printed output and plot
  annotations.

## Value

An object of class `cindex_weights`.

## Examples

``` r
weights_harrell()
#> <cindex_weights> Harrell 
weights_uno()
#> <cindex_weights> Uno (IPCW) 
weights_truncated(tau = 365)
#> <cindex_weights> Truncated (tau = 365) 
weights_custom(function(t, G) 1 / G(t), name = "IPCW (power 1)")
#> <cindex_weights> IPCW (power 1) 
```
