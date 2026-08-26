# cindexdecomp 0.2.0

## Breaking changes

* `simulate_censoring()` is renamed `censoring_curve()`.
* `plot_simulation()` and `plot_decomposition()` are replaced by
  `autoplot()` methods.
* `decompose_cindex()` returns a `cindex_decomp` object rather than a
  plain list. Existing element access (`$CI_ee`) becomes `$C_ee`.

## Bug fixes

* Event-censored pairs at tied times were dropped. The size of the error
  grows with tie density. Survival times are normally recorded in whole
  days, so ties are common: on a realistic day-scale fixture (300
  subjects, ~211 distinct times) the bug loses 53 comparable pairs and
  biases the estimate by +9.0e-05 relative to `survival::concordance()`;
  on an earlier, more heavily tied fixture (~31 distinct times for the
  same 300 subjects) the loss reached 1,062 pairs and the bias +0.0007.
* `simulate_censoring()` discarded every already-censored subject before
  computing, then reported the censoring rate among events only. A
  cut-off whose true cohort censoring was 51% could be reported as 20%.
* Global C-index values below 0.50 were inverted, reporting a 0.42 model
  as 0.58 and breaking the decomposition identity. Orientation is now
  declared once via `higher_is_riskier`.
* `C_ec` values below 0.50 were silently set to `NA`. Estimates are no
  longer discarded because of their magnitude; `min_pairs` sets a
  `low_precision` flag instead.
* Plots no longer request the "Arial" font family, which is absent on many
  Linux systems, and no longer impose fixed axis limits that hid
  worse-than-chance models.

## New features

* Pluggable pair weightings: `weights_harrell()`, `weights_uno()`,
  `weights_truncated()`, `weights_custom()`.
* Bootstrap inference over subjects, with `confint()` covering `C_ee`,
  `C_ec`, `C_global` and the masking gap.
* `compare_decompositions()` for multi-model comparison.
* `theme_cindex()`, light by default with an opt-in dark variant.
