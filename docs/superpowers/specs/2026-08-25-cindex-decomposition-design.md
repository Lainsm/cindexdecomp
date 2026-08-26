# cindexdecomp: Design Specification

**Date:** 2026-08-25
**Status:** Approved design, pending implementation plan
**Package version at time of writing:** 0.1.0

---

## 1. Purpose and publication path

`cindexdecomp` decomposes concordance indices for survival models into
Event-Event (`C_ee`) and Event-Censored (`C_ec`) components, exposing
performance masking in highly censored data.

The method is not novel; the implementation is. No R package currently
provides this decomposition. The package is therefore the scholarly
contribution.

**Two-paper sequence:**

1. **The R Journal** — a software paper. Contribution: the implementation
   fills a gap in the R ecosystem. Requires a CRAN package.
2. **A methods journal** (e.g. BMC Medical Research Methodology,
   Diagnostic and Prognostic Research) — an empirical paper on how often
   published prognostic models are masked. Cites paper 1.

These are distinct contributions (an implementation vs. a finding), so the
sequence is legitimate rather than salami-slicing. The package must be able
to *power* paper 2: multi-model comparison, batch runs, reproducibility.

**Quality bar:** CRAN-grade. `R CMD check` clean on three platforms, full
test suite, vignette, stable API.

---

## 2. Statistical foundations

### 2.1 The generalised decomposition

The decomposition is a **partition of comparable pairs**, not a property of
Harrell's C specifically. It applies to any concordance estimator of the
form:

```
C_w = sum_ij w_ij * c_ij / sum_ij w_ij
```

where the sum runs over comparable pairs `(i, j)`, `c_ij` in `{0, 0.5, 1}`
is the concordance indicator, and `w_ij` is a pair weight.

Partitioning comparable pairs into event-event and event-censored sets
gives the identity:

```
C_w = (W_ee * C_ee^w + W_ec * C_ec^w) / (W_ee + W_ec)
```

where `W` denotes **sums of weights**, not counts. This generalises the
README's current count-based identity, which is the `w_ij = 1` special
case.

Weightings covered:

| Estimator      | Pair weight `w_ij`     |
|----------------|------------------------|
| Harrell's C    | `1`                    |
| Uno's C        | `1 / Ghat(T_i)^2`      |
| Truncated C(t) | `1{T_i <= tau}`        |
| user-supplied  | any function of `(T_i, Ghat)` |

`Ghat` is the Kaplan-Meier estimate of the censoring distribution.

### 2.2 Explicit scope boundary

Model-based concordance estimators — notably **Gönen–Heller** — are **out
of scope by construction**, not by omission. They integrate over an assumed
model rather than counting observed comparable pairs, so the event-event /
event-censored partition is undefined for them. The paper states this as a
result: it defines the class of estimators to which the decomposition
applies.

### 2.3 Tie convention

Three cases, matching `survival::concordance()` exactly:

| Case                              | Rule                                    | Current code |
|-----------------------------------|-----------------------------------------|--------------|
| Two events at the same time       | not comparable, excluded                | correct      |
| Event and censoring at same time  | comparable; event precedes censoring    | **BUG: dropped** |
| Two identical risk scores         | comparable, half credit                 | correct      |

**The bug.** `R/decompose_cindex.R:51` uses `which(time > time[i])`, which
drops event-censored pairs at tied times. Measured on n=300 with times
rounded to whole days: **1,062 pairs lost**, pooled estimate biased
`+0.00071` relative to `survival::concordance()`.

**Fix:**

```r
j <- which(time > time[i] | (time == time[i] & status == 0))
```

Verified to restore agreement to `+0.00e+00`.

This is not an edge case. Real survival data is recorded in whole days
(`survival::lung` included), so ties are the normal case.

### 2.4 Uno's C under tied event times

`cindexdecomp` implements Uno's **pairwise** estimator as published (Uno et
al. 2011): every comparable pair is weighted by `1 / Ghat(T_i)^2`, where
`T_i` is the event time of the earlier member.

`survival::concordance(timewt = "n/G2")` uses a counting-process form that
applies its weight per event **time** rather than per pair. The two
formulations coincide exactly when event times are distinct, and diverge
once times are tied:

| Data | Harrell agreement | Uno agreement |
|---|---|---|
| Continuous times | `0.00e+00` | `-6.66e-16` |
| Day-scale, rounded | `0.00e+00` | `1.2e-4` |
| `survival::lung` (real) | `0.00e+00` | `2.7e-4` |

Two consequences, both deliberate:

1. **Per-pair weighting is required, not preferred.** A weight attached to
   an event *time* cannot be partitioned into event-event and
   event-censored pair sets, so the decomposition this package exists to
   compute could not be defined under `survival`'s formulation. The
   convention is forced by the package's purpose.
2. **The discrepancy is immaterial in size but must be documented.** At
   `~1e-4` it is two orders of magnitude below the C-index's own standard
   error (typically `0.02`-`0.03`). It is nonetheless stated plainly in
   `weights_uno()`'s documentation and in the vignette, because a user who
   cross-checks against `survival` on tied data will otherwise think they
   have found a bug.

Evaluating `Ghat` at `t-` rather than `t` was tested and does **not**
reconcile the two: it widens the gap from `3.4e-4` to `5.5e-4`.

Harrell's weighting is unaffected and matches exactly under all tie
patterns, so the package's central correctness claim stands unqualified.

---

## 3. Architecture

### 3.1 Public API

```r
# Core generic, two methods
decompose_cindex(Surv(time, status) ~ risk, data = df,
                 weights = weights_harrell(),
                 higher_is_riskier = TRUE, n_boot = 0)
decompose_cindex(time, status, risk, weights = weights_uno())

# Weight constructors (the extension point)
weights_harrell()                 # w = 1
weights_uno(tau = NULL)           # w = 1 / Ghat(T_i)^2
weights_truncated(tau)            # w = 1{T_i <= tau}
weights_custom(fn, name)          # w = fn(T_i, Ghat)

# S3 methods on class "cindex_decomp"
print()  summary()  confint()  as.data.frame()  autoplot()

# Multi-model comparison -> class "cindex_comparison"
compare_decompositions(risks, time, status,
                       weights = weights_harrell(), n_boot = 1000)

# Censoring sweep -> class "cindex_curve"
censoring_curve(time, status, risk, weights = weights_harrell(),
                n_thresholds = 20, min_pairs = NULL)
```

**Argument contracts** (pinned here because each has a plausible wrong
reading):

- `risks` in `compare_decompositions()` is a **named list or data frame** of
  risk-score vectors, each the same length as `time` and `status`. The names
  become the model labels in the printed table and the plot.
- `tau` in `weights_uno()` is a **truncation time**: pairs with `T_i > tau`
  receive weight zero. `NULL` (default) means no truncation. `Ghat` is
  always estimated from the full data regardless of `tau`.
- `fn` in `weights_custom()` must be `function(t, G)` returning a single
  **non-negative** scalar. `t` is the event time of the earlier member of
  the pair; `G` is a function returning the Kaplan-Meier censoring survival
  probability at a given time. Negative or non-finite weights are an error.
- `higher_is_riskier` declares the orientation of the risk score once, for
  the whole call. See Section 5.1.
- `n_boot = 0` (the default for `decompose_cindex()`) skips resampling.
  Calling `confint()` on a result computed with `n_boot = 0` is an error
  with a message naming the argument to set, not a silent `NA`.

**Naming constraint:** the generic must NOT be called `decompose()` —
`stats::decompose()` already exists (time-series decomposition) and masking
it triggers `R CMD check` complaints. Retain `decompose_cindex()`.

**Weight constructor names** are prefixed `weights_*` to avoid collisions
(`truncated` exists in `truncdist`; `uno` and `custom` are too generic for
the global namespace).

A weight constructor returns a list with two elements:

```r
weights_uno <- function(tau = NULL) {
  structure(list(name = "Uno (IPCW)", fn = function(t, G) 1 / G(t)^2),
            class = "cindex_weights")
}
```

Constructors must format their own `name` for display. The auto-generated
label `"Truncated tau=4.11833385100334"` observed during prototyping is
unacceptable; names are rounded and formatted at construction.

### 3.2 Renames

The package is pre-CRAN with no users, so renames are free. No deprecation
shims.

| Current              | New                        |
|----------------------|----------------------------|
| `simulate_censoring` | `censoring_curve`          |
| `plot_simulation`    | `autoplot.cindex_curve`    |
| `plot_decomposition` | `autoplot.cindex_comparison` |
| (missing)            | `compare_decompositions`   |

`simulate_censoring` is renamed because it does not simulate censoring from
a model; it sweeps administrative cutoffs.

### 3.3 The orphaned function

`plot_decomposition()` is currently **unreachable**. Its help page names
`monte_carlo_decompose()` as the source of its required
`ci_ee/ci_ec/global_c/sd_ee/sd_ec/sd_global` columns. That function does not
exist in the package. `compare_decompositions()` is that missing supplier.

### 3.4 File layout

```
R/  cindexdecomp-package.R   package doc, imports
    weights.R                weight constructors + validation
    engine.R                 internal pair-counting core (not exported)
    decompose.R              generic + formula/default methods
    methods.R                print/summary/confint/as.data.frame
    compare.R                compare_decompositions()
    censoring_curve.R        threshold sweep
    autoplot.R               ggplot2 methods
    utils.R                  input validation helpers
```

The engine takes vectors plus a weight function and returns weight-sums and
concordance sums. It knows nothing about formulas, models, or plotting.
This is what makes `weights_custom()` a genuine extension point.

### 3.5 Implementation language

**Pure R. No compiled code.** Measured performance of the existing loop:

| n      | time  |
|--------|-------|
| 500    | 0.013s |
| 2,000  | 0.049s |
| 5,000  | 0.234s |
| 10,000 | 0.795s |

A 1000-replicate bootstrap at n=5,000 costs ~4 minutes. Acceptable. Avoiding
Rcpp removes the compilation toolchain requirement, simplifies CRAN
submission, and eliminates binary-build failures on CRAN's platforms. If the
methods paper needs more speed, Rcpp can go behind the same API later.

### 3.6 Dependencies

Drop `tidyr` (unused). Drop `dplyr` and `magrittr` if the only uses are
`bind_rows` and a single pipe — both replaceable with base R. Retain
`survival`, `ggplot2`, `scales`.

---

## 4. Inference

### 4.1 Method: bootstrap over subjects

`confint()` on a `cindex_decomp` object returns percentile intervals for
`C_ee`, `C_ec`, `C_global`, and the **gap** `C_ec - C_ee`.

Rationale for bootstrap over an analytic variance:

- Works for **every** weighting including `weights_custom()`. An analytic
  variance would need re-derivation per weighting, defeating the framework.
- Gives a valid interval for the gap, which is the paper's headline
  quantity. The gap cannot be obtained by subtracting two marginal
  intervals because the components are correlated.
- Measured cost is acceptable (Section 3.5). `n_boot = 0` skips it.

### 4.2 Resample subjects, not pairs

Pairs are dependent — each subject appears in hundreds of them. Measured on
n=300 (16,471 event-event pairs):

| Method                              | SE(C_ee) | CI width |
|-------------------------------------|----------|----------|
| Pairs treated as independent        | 0.00381  | 0.0149   |
| Subject-level bootstrap             | 0.02174  | 0.0852   |

The pair-based SE is **5.7x too small**. Documented and pinned by a test.

### 4.3 The gap requires a paired bootstrap

Measured correlation between `C_ee` and `C_ec` across bootstrap replicates:
**+0.384**. Assuming independence overstates `SE(gap)` by **26%**
(0.03614 vs 0.02871). Because the bootstrap resamples subjects, both
components move together within each replicate and the correlation is
handled automatically.

### 4.4 Reporting constraint

`print()` shows the compact table (components, pair counts, shares, global,
masking gap) plus any warnings. `summary()` adds the confidence intervals,
the weighting's name, and the bootstrap settings used.

`print()` must show `N_ee` **alongside** its confidence interval, never
alone. A pair count of 16,471 from 300 patients looks like a large sample
and is not one. The pair count describes composition; the CI describes
precision. Conflating them is the exact error the package exists to prevent.

---

## 5. Behaviours to remove

Both are in `R/simulate_censoring.R` and both suppress unfavourable results.

### 5.1 The sign flip (line 80-81)

```r
global_c <- ifelse(global_c_raw < 0.5, 1 - global_c_raw, global_c_raw)
```

A model scoring 0.42 is reported as 0.58. Worse-than-chance performance is a
genuine finding, often the most interesting one. This also breaks the
package's own identity: only the global value is flipped, so it ceases to
equal the weighted mean of its components, and the plot shows a global line
contradicting the components beneath it.

**Replacement:** the explicit `higher_is_riskier` argument (default `TRUE`),
applied once to the whole call rather than inferred per threshold from the
answer. When `FALSE`, the risk score is negated before any pair counting, so
the identity holds and all three quantities stay mutually consistent. Values
below 0.5 are then reported as they are.

### 5.2 Value-based suppression (line 88-89)

```r
ci_ec <- ifelse(decomp$N_ec < min_pairs | decomp$CI_ec < 0.5, NA, decomp$CI_ec)
```

The `N_ec < min_pairs` half is defensible (imprecision). The `CI_ec < 0.5`
half deletes results **because of their value** — precisely the points that
demonstrate the package's thesis.

**Replacement:** always report the estimate; attach pair count and CI width;
let `autoplot()` visually de-emphasise low-precision points rather than
deleting them. `min_pairs` becomes a flag, not a filter.

---

## 6. Censoring curve

### 6.1 Keep the whole cohort

`R/simulate_censoring.R:47` begins `event_idx <- which(status == 1)`,
discarding every already-censored observation before anything else happens.
The resulting curve describes the event-only subset, not the user's cohort,
while the x-axis is labelled "Artificial Censoring Rate".

Measured on n=500 with 39% baseline censoring:

| tau | reported censoring (current) | actual cohort censoring | understatement |
|-----|------------------------------|-------------------------|----------------|
| 1.1 | 80%                          | 88%                     | 8 pp           |
| 2.7 | 60%                          | 76%                     | 16 pp          |
| 4.9 | 40%                          | 64%                     | 24 pp          |
| 8.6 | 20%                          | 51%                     | **31 pp**      |

Every point is plotted at the correct height and the wrong horizontal
position. `C_ec` and `C_global` also differ materially (at tau=8.6: global
0.6433 vs 0.7067).

`C_ee` is **identical** under both approaches — event-event pairs only
involve observed events, so the core quantity was never wrong.

**Fix:** apply the cutoff to everyone.

```r
sim_time   <- pmin(time, tau)
sim_status <- status * as.numeric(time <= tau)
```

The corrected version demonstrates the thesis more cleanly: as censoring
rises 51% -> 88%, `C_global` moves 0.002 (0.7067 -> 0.7086, flat) while
`C_ee` falls 0.053 (0.5967 -> 0.5432, toward chance).

### 6.2 Documentation fix

`min_pairs` is documented as "Default is 5000"; the actual default is
`NULL`. Correct the documentation to match the auto-selection behaviour.

---

## 7. Plotting

All plots become `autoplot()` methods on their result classes, so a result
can never be paired with the wrong plot. (Currently
`plot_simulation(decompose_result)` fails with
`"incorrect number of dimensions"`.)

Three CRAN-blocking issues:

1. **`family = "Arial"` is hardcoded 14 times** across the two plot files.
   Arial does not exist on most Linux systems including CRAN check machines;
   this emits font warnings on every call and can fail checks. Default to
   the system font.
2. **The Dracula dark theme becomes opt-in.** `#282a36` backgrounds suit a
   talk and fail a journal figure printed on white paper. Default to a clean
   light theme; ship the dark variant as `theme_cindex(dark = TRUE)`.
3. **Hardcoded axis limits removed.** `scale_x_continuous(limits =
   c(0.48, 0.74))` in `plot_decomposition()` silently drops any model
   outside that window, including worse-than-chance models.

---

## 8. Testing

The suite deleted in commit `b77fb30` returns. Organising principle:
**`survival::concordance()` is ground truth.**

1. **Identity tests.** Under Harrell's weighting, the pooled decomposition
   equals `survival::concordance()` to machine precision (`1e-12`) across
   continuous times, tied times, tied risks, and both tied. Measured at
   `0.00e+00` in every case including real `survival::lung` data. This layer
   would have caught the 1,062 missing pairs. Uno's weighting agrees exactly
   only when event times are distinct -- see Section 2.4.
2. **Weighting tests.** `weights_harrell()` matches `survival`'s default at
   `1e-12` under every tie pattern. `weights_uno()` matches
   `timewt = "n/G2"` at `1e-10` when event times are distinct, and to `1e-3`
   when they are tied (Section 2.4).
3. **Property tests.** Negating `risk` gives `1 - C`; the weighted identity
   holds for arbitrary custom weights; `N_ee + N_ec` equals `survival`'s
   comparable-pair count.
4. **Regression tests.** No sign flip; no value-based suppression; censoring
   curve reports whole-cohort rates.

---

## 9. Infrastructure

- GitHub Actions: `R CMD check` on Linux, macOS, Windows.
- Test coverage reporting.
- `pkgdown` site.
- `NEWS.md`.
- **`.Rbuildignore` currently excludes `LICENSE.md` while `DESCRIPTION`
  declares `MIT + file LICENSE`.** That combination fails `R CMD check`. Fix
  by shipping the licence file.

---

## 10. Vignette

The vignette is the backbone of the R Journal paper:

1. The decomposition and its identity.
2. The tie convention and why it matters.
3. The weighting framework, with the three-estimator table below.
4. A worked `survival::lung` example.
5. The censoring curve.
6. The Uno tie convention of Section 2.4, stated plainly so a user
   cross-checking against `survival` on tied data is not surprised.

**The central table.** One engine, three estimators, n=400, 44% censoring:

| weighting     | `C_global` | `C_ee` | `C_ec` |
|---------------|------------|--------|--------|
| Harrell       | 0.6936     | 0.6235 | 0.7736 |
| Uno (IPCW)    | 0.6825     | 0.6171 | 0.7523 |
| Truncated     | 0.7035     | 0.6302 | 0.7914 |

`C_global` varies with the estimator chosen (0.68-0.70), and every value is
defensible and publishable. `C_ee` sits at ~0.62 regardless, roughly 0.15
below `C_ec` in every row.

**The masking is invariant to the choice of estimator.** Uno's IPCW
correction moves the headline number by 0.011 and leaves the composition
problem entirely invisible. This is the package's answer to the
predictable reviewer objection ("doesn't Uno's C already handle censoring?")
and it is an empirical claim, not a rhetorical one.

---

## 11. Non-goals

- Gönen–Heller and other model-based estimators (Section 2.2).
- Calibration assessment. The package evaluates discrimination only.
- Model fitting. The package consumes risk scores from any source.
- Competing risks. Single-event survival only in v1.
- Informative censoring. Non-informative censoring is assumed throughout
  and stated as a limitation.

---

## 12. Sequencing

The package ships first; the methods paper cites it. Nothing in this spec
depends on the methods paper's analyses. `compare_decompositions()` is the
component that will carry that work when it begins.
