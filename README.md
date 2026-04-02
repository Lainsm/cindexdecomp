# Package Documentation: cindexdecomp

## Overview

`cindexdecomp` is an R package that decomposes Harrell's Concordance 
Index into its Event-Event (C_ee) and Event-Censored (C_ec) components. 
It exposes performance masking in highly censored survival data through 
decomposition and simulation.


## Installation
```r
# Install from GitHub
devtools::install_github("lainsm/cindexdecomp")
```

## Core Functions  

### 1. `decompose_cindex()`

Decomposes Harrell's C-Index into Event-Event and Event-Censored 
components.

**Arguments:**
- `time`: Numeric vector of survival times
- `status`: Numeric vector of event (1 = event, 0 = censored)
- `risk`: Numeric vector of predicted risk scores

**Returns:**
- `CI_ee`: Event-Event concordance index
- `CI_ec`: Event-Censored concordance index
- `N_ee`: Number of Event-Event comparable pairs
- `N_ec`: Number of Event-Censored comparable pairs

**Mathematical basis:**

$$CI_{ee} = \frac{N^+_{ee} + \frac{1}{2}N^=_{ee}}{N^+_{ee} + N^-_{ee} + N^=_{ee}}$$

$$CI_{ec} = \frac{N^+_{ec} + \frac{1}{2}N^=_{ec}}{N^+_{ec} + N^-_{ec} + N^=_{ec}}$$

$$C_{global} = \frac{N_{ee} \cdot CI_{ee} + N_{ec} \cdot CI_{ec}}{N_{ee} + N_{ec}}$$

**Example:**
```r
library(cindexdecomp)
library(survival)

# Fit Cox model on lung dataset
lung_clean <- lung[complete.cases(lung), ]
fit  <- coxph(Surv(time, status - 1) ~ age + sex + ph.ecog,
              data = lung_clean)
risk <- predict(fit, type = "lp")

# Decompose
result <- decompose_cindex(
  time   = lung_clean$time,
  status = lung_clean$status - 1,
  risk   = risk
)

print(result)
```

---

### 2. `simulate_censoring()`

Implements a "Time Machine" simulation. It freezes risk scores and 
applies artificial censoring thresholds iteratively to expose the 
Censoring Illusion.

**Arguments:**
- `time`: Numeric vector of survival times
- `status`: Numeric vector of event (1 = event, 0 = censored)
- `risk`: Numeric vector of predicted risk scores
- `n_thresholds`: Number of censoring points (default 20)
- `min_pairs`: Minimum EE pairs required. If nothing is selected it has a default auto-calculated minumum threshold

**Returns:**

A data frame with columns:
- `threshold`: Applied time threshold
- `censoring`: Resulting censoring rate
- `global_c`: Global C-Index score
- `ci_ee`: Event-Event C-Index score
- `ci_ec`: Event-Censored C-Index score
- `n_ee`: Number of EE pairs
- `n_ec`: Number of EC pairs

**Design note:**

The simulation is restricted to patients with confirmed events 
(ground-truth cohort). Risk scores are frozen from the pre-trained 
model you input. Any change in C-Index is therefore attributable purely 
to censoring, not model performance.

**Example:**
```r
sim <- simulate_censoring(
  time   = lung_clean$time,
  status = lung_clean$status - 1,
  risk   = risk
)

head(sim)
```

---

### 3. `plot_simulation()`

Visualises the Time Machine simulation as a line chart with the 
Censoring Illusion Zone shaded between Global C and C_ee.

**Arguments:**
- `sim_df`: Output from `simulate_censoring()`
- `col_ee`: Colour for C_ee line (default IBM magenta `#DC267F`)
- `col_ec`: Colour for C_ec line (default IBM blue `#648FFF`)
- `background`: Plot background (default Dark/Grey `#282a36`)

**Returns:** A ggplot2 object

**Example:**
```r
plot_simulation(sim)
```

---

### 4. `plot_decomposition()`

Visualises C-Index decomposition across one or multiple models built as a 
dumbbell chart with uncertainty bands.

**Arguments:**
- `results_df`: Data frame with columns: model, ci_ee, ci_ec, 
  global_c, sd_ee, sd_ec, sd_global
- `col_ee`: Colour for C_ee (default IBM magenta `#DC267F`)
- `col_ec`: Colour for C_ec (default IBM blue `#648FFF`)
- `background`: Plot background (default Dark/Grey `#282a36`)

**Returns:** A ggplot2 object

**Example:**
```r
plot_decomposition(results_df)
```

---

## Complete Workflow
```r
library(cindexdecomp)
library(survival)

# 1. Prepare data
lung_clean <- lung[complete.cases(lung), ]

# 2. Train model
fit  <- coxph(Surv(time, status - 1) ~ age + sex + ph.ecog,
              data = lung_clean)
risk <- predict(fit, type = "lp")

# 3. Decompose C-Index
result <- decompose_cindex(
  time   = lung_clean$time,
  status = lung_clean$status - 1,
  risk   = risk
)
print(result)

# 4. Run simulation
sim <- simulate_censoring(
  time   = lung_clean$time,
  status = lung_clean$status - 1,
  risk   = risk
)

# 5. Plot simulation
plot_simulation(sim)

# 6. Plot decomposition across models
plot_decomposition(results_df)
```

---

## Mathematical Background

Harrell's C-Index pools two fundamentally different tasks:

**Event-Event pairs (C_ee):** Both patients experience the event. 
The model must correctly rank which patient fails first. This is 
the clinically meaningful task.

**Event-Censored pairs (C_ec):** One patient experiences the event, 
one is censored. The model must rank the event patient as higher risk 
than the censored patient before their censoring time. This is 
inherently easier.

As censoring increases:
- N_ee → 0
- N_ec dominates
- C_global → C_ec
- True model performance (C_ee) is masked

## Limitations

- C_ee becomes unstable when N_ee is small — use `min_pairs` to 
  control this
- Does not account for informative censoring
- Simulation assumes censoring is non-informative
- The C-index decompostion checks the nature of the score and not the model 
- Results depend on quality of risk scores from the user's model

## Citation

If you use `cindexdecomp` in your research please cite:

> Lainsbury MJ (2025). cindexdecomp: C-Index Decomposition for 
> Survival Models. R package version 0.1.0.
> https://github.com/lainsm/cindexdecomp

## References

Harrell FE et al. (1982) Evaluating the yield of medical tests. 
*JAMA*, 247(18), 2543-2546.

Uno H et al. (2011) On the C-statistics for evaluating overall 
adequacy of risk prediction procedures with censored survival data. 
*Statistics in Medicine*, 30(10), 1105-1117.