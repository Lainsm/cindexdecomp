# cindexdecomp

> C-Index Decomposition for Survival Models

[![R-CMD-check](https://github.com/Lainsm/cindexdecomp/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Lainsm/cindexdecomp/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)]()

## Overview

Harrell's concordance index pools two different ranking tasks:

- **Event-Event pairs (`C_ee`)** — both subjects had the event. The harder task.
- **Event-Censored pairs (`C_ec`)** — one had the event, one was censored later.

As censoring rises the event-event pairs vanish and `C_ec` dominates the
pooled score, so a model close to chance on true events can still report a
stable global C-index. `cindexdecomp` makes that visible.

## Installation

```r
# install.packages("devtools")
devtools::install_github("Lainsm/cindexdecomp")
```

## Quick start

```r
library(cindexdecomp)
library(survival)

lung2 <- lung[complete.cases(lung), ]
cox <- coxph(Surv(time, status - 1) ~ age + sex + ph.ecog, data = lung2)

fit <- decompose_cindex(
  lung2$time,
  lung2$status - 1,            # lung codes 1/2; the package expects 0/1
  predict(cox, type = "lp"),
  n_boot = 1000
)
fit
confint(fit)
```

## The identity

Any concordance index that averages over comparable pairs,

```
C_w = sum(w_ij * c_ij) / sum(w_ij)
```

decomposes exactly:

```
C_w = (W_ee * C_ee + W_ec * C_ec) / (W_ee + W_ec)
```

where `W` denotes sums of weights. Setting `w_ij = 1` recovers Harrell's C.

## Weightings

| Constructor | Estimator |
|---|---|
| `weights_harrell()` | Harrell's C |
| `weights_uno()` | Uno's IPCW C |
| `weights_truncated(tau)` | C truncated at `tau` |
| `weights_custom(fn, name)` | any user-supplied rule |

The pooled result agrees with `survival::concordance()` to machine
precision for Harrell's weighting under every tie pattern, and for Uno's
weighting when event times are distinct. Under tied times, Uno's
pairwise estimator diverges from `survival`'s counting-process form by
around 1e-04 -- immaterial in practice, and explained in the vignette
(`vignette("cindexdecomp")`).

## Functions

| Function | Description |
|---|---|
| `decompose_cindex()` | Decomposes into `C_ee` and `C_ec` |
| `confint()` | Bootstrap intervals, including for the masking gap |
| `censoring_curve()` | Traces the decomposition as censoring rises |
| `compare_decompositions()` | Runs several models over one cohort |
| `autoplot()` | Plots any of the above |
| `theme_cindex()` | Light (default) and dark plot themes |

## Model agnostic

Any model that produces a risk score works:

```r
decompose_cindex(time, status, predict(cox_fit, type = "lp"))
decompose_cindex(time, status, predict(xgb_fit, xgb.DMatrix(X)))
decompose_cindex(time, status, predict(rsf_fit, newdata = test)$predicted)
```

## Limitations

- Assumes non-informative censoring.
- `C_ee` is imprecise when few event-event pairs exist. Use `confint()`;
  pair counts describe composition, not precision.
- Evaluates discrimination only, not calibration.
- Model-based estimators such as Gönen–Heller are out of scope: they
  integrate over an assumed model rather than counting observed pairs, so
  the event-event partition is undefined for them.

## Citation

```
Lainsbury MJ (2026). cindexdecomp: C-Index Decomposition for Survival
Models. R package version 0.2.0.
https://github.com/Lainsm/cindexdecomp
```

## References

Harrell FE et al. (1982) Evaluating the yield of medical tests.
*JAMA*, 247(18), 2543-2546.

Uno H et al. (2011) On the C-statistics for evaluating overall adequacy of
risk prediction procedures with censored survival data.
*Statistics in Medicine*, 30(10), 1105-1117.
