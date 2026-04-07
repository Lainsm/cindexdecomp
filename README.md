# cindexdecomp

> C-Index Decomposition for Survival Models

[![R package](https://img.shields.io/badge/R-package-blue)]()
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)]()

## Overview

Harrell's Concordance Index is used in over 80% of survival model 
publications. In highly censored datasets it creates a false illusion 
of performance by pooling two fundamentally different ranking tasks:

- **Event-Event pairs (C_ee):** Both patients experience the event. This is the harder task to rank.
- **Event-Censored pairs (C_ec):** One patient experiences the event, 
 one is censored.

As censoring rises, Event-Event pairs vanish and C_ec dominates the 
global score. A model near random chance on true events (C_ee ≈ 0.50) 
can still report a misleadingly stable Global C-Index.

`cindexdecomp` makes this visible.

## Installation
```r
# Install from GitHub
devtools::install_github("Lainsm/cindexdecomp")
```

## Quick Start
```r
library(cindexdecomp)
library(survival)

# Prepare data
lung_clean <- lung[complete.cases(lung), ]
fit  <- coxph(Surv(time, status - 1) ~ age + sex + ph.ecog,
              data = lung_clean)
risk <- predict(fit, type = "lp")

# Decompose C-Index
result <- decompose_cindex(
  time   = lung_clean$time,
  status = lung_clean$status - 1,
  risk   = risk
)
print(result)

# Run Time Machine simulation
sim <- simulate_censoring(
  time   = lung_clean$time,
  status = lung_clean$status - 1,
  risk   = risk
)

# Visualise
plot_simulation(sim)
```

## Functions

| Function | Description |
|---|---|
| `decompose_cindex()` | Decomposes C-Index into C_ee and C_ec |
| `simulate_censoring()` | Time Machine simulation |
| `plot_simulation()` | Illusion Zone line chart |
| `plot_decomposition()` | Dumbbell chart across models |

## The Math

$$C_{global} = \frac{N_{ee} \cdot C_{ee} + N_{ec} \cdot C_{ec}}{N_{ee} + N_{ec}}$$

$$CI_{ee} = \frac{N^+_{ee} + \frac{1}{2}N^=_{ee}}{N^+_{ee} + N^-_{ee} + N^=_{ee}}$$

$$CI_{ec} = \frac{N^+_{ec} + \frac{1}{2}N^=_{ec}}{N^+_{ec} + N^-_{ec} + N^=_{ec}}$$

## Model Agnostic

`cindexdecomp` works with any model that outputs a risk score:
```r
# Cox PH
cox_risk <- predict(cox_fit, type = "lp")
decompose_cindex(time, status, cox_risk)

# XGBoost
xgb_risk <- predict(xgb_fit, xgb.DMatrix(X))
decompose_cindex(time, status, xgb_risk)

# Random Survival Forest
rsf_risk <- predict(rsf_fit, newdata = test)$predicted
decompose_cindex(time, status, rsf_risk)
```

## Limitations

- Assumes non-informative censoring meaning that the censoring is not dependant of the event of interest.
- C_ee becomes extremely unstable when N_ee is small.
- Evaluates discrimination only, not calibration.
- The simulation freezes risk scores from pre-trained model.

## Citation
```
Lainsbury MJ (2025). cindexdecomp: C-Index Decomposition for 
Survival Models. R package version 0.1.0.
https://github.com/Lainsm/cindexdecomp
```

## References

Harrell FE et al. (1982) Evaluating the yield of medical tests.
*JAMA*, 247(18), 2543-2546.

Uno H et al. (2011) On the C-statistics for evaluating overall
adequacy of risk prediction procedures with censored survival data.
*Statistics in Medicine*, 30(10), 1105-1117.
