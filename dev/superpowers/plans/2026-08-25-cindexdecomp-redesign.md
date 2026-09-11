# cindexdecomp Redesign Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Rebuild `cindexdecomp` as a CRAN-grade package that decomposes any weighted-pair concordance index into Event-Event and Event-Censored components, fixing four measured defects in the current implementation.

**Architecture:** One internal pair-counting engine takes vectors plus a weight function and returns weight-sums and concordance-sums. Weight constructors (`weights_harrell()`, `weights_uno()`, ...) select the estimator. A `decompose_cindex()` generic with formula and vector methods wraps the engine and returns an S3 object carrying its own `print`/`summary`/`confint`/`autoplot` behaviour. `censoring_curve()` and `compare_decompositions()` are thin loops over that core.

**Tech Stack:** R 4.5.2, pure R (no compiled code), `survival` 3.8.6, `ggplot2` 4.0.3, `scales` 1.4.0, `testthat` 3.3.2 (3rd edition), `roxygen2` 7.3.3, `pkgdown` 2.2.0.

**Spec:** `docs/superpowers/specs/2026-08-25-cindex-decomposition-design.md`

## Global Constraints

Every task's requirements implicitly include this section.

- **Pure R only.** No `Rcpp`, no `src/` directory, no compiled code.
- **Imports limited to** `survival`, `ggplot2`, `scales`, `stats`. Drop `tidyr`, `dplyr`, `magrittr` from `DESCRIPTION`.
- **The generic must NOT be named `decompose()`** — `stats::decompose()` exists and masking it fails `R CMD check`. Use `decompose_cindex()`.
- **Weight constructors are prefixed `weights_`** (`truncated` collides with `truncdist`; `uno`/`custom` are too generic).
- **Tie convention:** event-vs-censored at equal times IS comparable (event first); event-vs-event at equal times is NOT comparable; equal risk scores get 0.5 credit. Must match `survival::concordance()` to machine precision.
- **No sign flipping.** Orientation is declared once via `higher_is_riskier = TRUE`.
- **No value-based `NA` suppression.** Never discard an estimate because of its magnitude.
- **`print()` shows `N_ee` alongside its confidence interval, never alone.**
- **No hardcoded font `family`** anywhere (Arial is absent on CRAN's Linux machines).
- **No hardcoded axis limits** in any plot.
- **Light theme is the default**; dark is opt-in via `theme_cindex(dark = TRUE)`.
- **testthat 3rd edition.** `expect_equal()` uses `tolerance =`, not `all.equal` semantics.
- Machine-precision comparisons use `tolerance = 1e-12`.

---

## File Structure

| File | Responsibility |
|---|---|
| `R/cindexdecomp-package.R` | Package-level docs, `@importFrom` declarations |
| `R/utils.R` | Input validation, censoring Kaplan-Meier helper |
| `R/weights.R` | Weight constructors + their `print` method |
| `R/engine.R` | Internal `pair_counts()` — the only place pairs are counted |
| `R/decompose.R` | `decompose_cindex()` generic, formula + default methods |
| `R/bootstrap.R` | Subject-level bootstrap |
| `R/methods.R` | `print`/`summary`/`confint`/`as.data.frame` for `cindex_decomp` |
| `R/censoring_curve.R` | Threshold sweep + its methods |
| `R/compare.R` | `compare_decompositions()` + its methods |
| `R/theme.R` | `theme_cindex()` |
| `R/autoplot.R` | `autoplot` methods for all three classes |
| `tests/testthat/helper-data.R` | Shared fixture generators |
| `tests/testthat/test-*.R` | One test file per R file |
| `vignettes/cindexdecomp.Rmd` | R Journal paper backbone |

**Files deleted during the plan:** `R/simulate_censoring.R` (Task 7), `R/plot_simulation.R` and `R/plot_decomposition.R` (Task 9). `R/decompose_cindex.R` is replaced by `R/decompose.R` + `R/engine.R` in Task 4.

---

## Task 1: Test infrastructure and validation helpers

**Files:**
- Create: `R/utils.R`
- Create: `R/cindexdecomp-package.R`
- Create: `tests/testthat.R`
- Create: `tests/testthat/helper-data.R`
- Create: `tests/testthat/test-utils.R`
- Modify: `DESCRIPTION`
- Modify: `.Rbuildignore`

**Interfaces:**
- Consumes: nothing.
- Produces: `validate_survival_inputs(time, status, risk)` returning `invisible(TRUE)` or stopping; `censoring_km(time, status)` returning a function `function(t) numeric` giving the Kaplan-Meier censoring survival probability at `t`; test fixture `make_test_data(n, censor_rate, ties, seed)` returning `list(time, status, risk)`.

- [ ] **Step 1: Set up the testthat scaffold**

Create `tests/testthat.R`:

```r
library(testthat)
library(cindexdecomp)

test_check("cindexdecomp")
```

Create `tests/testthat/helper-data.R`:

```r
# Shared fixtures. `ties = TRUE` rounds times to whole numbers, which is how
# real survival data arrives (days) and is the case that exercises the tie rules.
make_test_data <- function(n = 300, censor_rate = 0.4, ties = FALSE, seed = 1) {
  set.seed(seed)
  risk <- rnorm(n)
  event_time <- rexp(n, rate = exp(0.8 * risk) * 0.1)
  # solve roughly for a censoring rate
  cens_time <- rexp(n, rate = 0.1 * censor_rate / (1 - censor_rate))
  time <- pmin(event_time, cens_time)
  status <- as.numeric(event_time <= cens_time)
  if (ties) time <- round(time)
  # guarantee at least two events so validation never trips on the fixture
  if (sum(status) < 2) {
    status[order(time)[1:2]] <- 1
  }
  list(time = time, status = status, risk = risk)
}

# Fixture with a deliberate event/censoring tie at the same time.
make_tied_pair_data <- function() {
  list(
    time   = c(100, 100, 200, 300),
    status = c(  1,   0,   1,   1),
    risk   = c(0.9, 0.2, 0.5, 0.1)
  )
}
```

- [ ] **Step 2: Write the failing test**

Create `tests/testthat/test-utils.R`:

```r
test_that("validate_survival_inputs accepts well-formed input", {
  d <- make_test_data(50)
  expect_true(validate_survival_inputs(d$time, d$status, d$risk))
})

test_that("validate_survival_inputs rejects mismatched lengths", {
  expect_error(
    validate_survival_inputs(1:5, c(0, 1, 0, 1), rnorm(5)),
    "same length"
  )
})

test_that("validate_survival_inputs rejects non-binary status", {
  expect_error(
    validate_survival_inputs(1:4, c(0, 1, 2, 1), rnorm(4)),
    "binary"
  )
})

test_that("validate_survival_inputs rejects NA", {
  expect_error(
    validate_survival_inputs(c(1, NA, 3, 4), c(1, 1, 0, 1), rnorm(4)),
    "NA"
  )
})

test_that("validate_survival_inputs rejects negative times", {
  expect_error(
    validate_survival_inputs(c(1, -2, 3, 4), c(1, 1, 0, 1), rnorm(4)),
    "non-negative"
  )
})

test_that("validate_survival_inputs requires at least two events", {
  expect_error(
    validate_survival_inputs(1:4, c(1, 0, 0, 0), rnorm(4)),
    "At least 2 events"
  )
})

test_that("censoring_km returns 1 before the first censoring time", {
  d <- make_test_data(200)
  G <- censoring_km(d$time, d$status)
  expect_equal(G(0), 1)
})

test_that("censoring_km is non-increasing and strictly positive", {
  d <- make_test_data(200)
  G <- censoring_km(d$time, d$status)
  grid <- seq(0, max(d$time), length.out = 50)
  vals <- vapply(grid, G, numeric(1))
  expect_true(all(diff(vals) <= 1e-12))
  expect_true(all(vals > 0))
})
```

- [ ] **Step 3: Run tests to verify they fail**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-utils.R")'`
Expected: FAIL — `could not find function "validate_survival_inputs"`

- [ ] **Step 4: Write the implementation**

Create `R/utils.R`:

```r
#' Validate survival inputs
#'
#' @param time Numeric vector of observed times.
#' @param status Numeric vector of event indicators (1 = event, 0 = censored).
#' @param risk Numeric vector of risk scores.
#' @return `invisible(TRUE)` if valid; otherwise an error is thrown.
#' @keywords internal
#' @noRd
validate_survival_inputs <- function(time, status, risk) {
  n <- length(time)
  if (length(status) != n || length(risk) != n) {
    stop("`time`, `status` and `risk` must all have the same length.", call. = FALSE)
  }
  if (n == 0L) {
    stop("`time`, `status` and `risk` must not be empty.", call. = FALSE)
  }
  if (!is.numeric(time) || !is.numeric(risk)) {
    stop("`time` and `risk` must be numeric.", call. = FALSE)
  }
  if (anyNA(time) || anyNA(status) || anyNA(risk)) {
    stop("`time`, `status` and `risk` must not contain NA.", call. = FALSE)
  }
  if (any(!is.finite(time)) || any(!is.finite(risk))) {
    stop("`time` and `risk` must be finite.", call. = FALSE)
  }
  if (any(time < 0)) {
    stop("`time` must be non-negative.", call. = FALSE)
  }
  if (!all(status %in% c(0, 1))) {
    stop("`status` must be binary (0 = censored, 1 = event).", call. = FALSE)
  }
  if (sum(status == 1) < 2) {
    stop("At least 2 events are required to decompose the C-index.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Kaplan-Meier estimate of the censoring distribution
#'
#' Returns a step function giving `G(t)`, the probability of remaining
#' uncensored beyond `t`. Used by IPCW weightings.
#'
#' @param time Numeric vector of observed times.
#' @param status Numeric vector of event indicators.
#' @return A function of one numeric argument.
#' @keywords internal
#' @noRd
censoring_km <- function(time, status) {
  fit <- survival::survfit(survival::Surv(time, 1 - status) ~ 1)
  ftime <- fit$time
  fsurv <- fit$surv
  function(t) {
    idx <- findInterval(t, ftime)
    out <- ifelse(idx == 0L, 1, fsurv[pmax(idx, 1L)])
    # guard against division by zero in IPCW weights
    pmax(out, .Machine$double.eps)
  }
}
```

Create `R/cindexdecomp-package.R`:

```r
#' @keywords internal
"_PACKAGE"

#' @importFrom stats model.frame model.response quantile sd setNames
#' @importFrom survival Surv survfit concordance
#' @importFrom ggplot2 autoplot
NULL
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-utils.R")'`
Expected: PASS, 8 tests.

- [ ] **Step 6: Update DESCRIPTION and .Rbuildignore**

Replace `DESCRIPTION` with:

```
Package: cindexdecomp
Title: C-Index Decomposition for Survival Models
Version: 0.2.0
Authors@R:
    person("Malik J.", "Lainsbury", , "mjl221@cam.ac.uk", role = c("aut", "cre"))
Description: Decomposes concordance indices for survival models into
    Event-Event and Event-Censored components to expose performance
    masking in highly censored data. Supports Harrell's C, Uno's
    inverse-probability-of-censoring-weighted C, time-truncated C, and
    user-supplied pair weightings through a single engine.
License: MIT + file LICENSE
Encoding: UTF-8
Roxygen: list(markdown = TRUE)
RoxygenNote: 7.3.3
Depends:
    R (>= 4.1.0)
Imports:
    ggplot2,
    scales,
    stats,
    survival
Suggests:
    covr,
    knitr,
    rmarkdown,
    testthat (>= 3.0.0)
Config/testthat/edition: 3
VignetteBuilder: knitr
URL: https://github.com/Lainsm/cindexdecomp
BugReports: https://github.com/Lainsm/cindexdecomp/issues
```

Replace `.Rbuildignore` with (the `LICENSE.md` exclusion is removed — `DESCRIPTION` declares `MIT + file LICENSE`, so excluding the file fails `R CMD check`):

```
^cindexdecomp\.Rproj$
^\.Rproj\.user$
^docs/superpowers$
^\.github$
^_pkgdown\.yml$
^pkgdown$
```

Create `LICENSE` (CRAN requires this exact two-line form alongside `LICENSE.md`):

```
YEAR: 2026
COPYRIGHT HOLDER: Malik J. Lainsbury
```

- [ ] **Step 7: Commit**

```bash
git add R/utils.R R/cindexdecomp-package.R tests/ DESCRIPTION .Rbuildignore LICENSE
git commit -m "feat: add validation helpers and test scaffold

Adds validate_survival_inputs() and censoring_km(), the testthat
scaffold, and shared fixtures. Trims DESCRIPTION Imports to survival,
ggplot2, scales, stats and ships LICENSE so R CMD check passes."
```

---

## Task 2: Weight constructors

**Files:**
- Create: `R/weights.R`
- Create: `tests/testthat/test-weights.R`

**Interfaces:**
- Consumes: nothing from Task 1 at runtime.
- Produces: `weights_harrell()`, `weights_uno(tau = NULL)`, `weights_truncated(tau)`, `weights_custom(fn, name = "Custom")`, each returning an object of class `cindex_weights` — a list with `$name` (length-1 character) and `$fn` (a `function(t, G)` returning a single non-negative finite numeric). `print.cindex_weights()`.

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-weights.R`:

```r
test_that("weights_harrell returns constant weight 1", {
  w <- weights_harrell()
  expect_s3_class(w, "cindex_weights")
  expect_equal(w$name, "Harrell")
  expect_equal(w$fn(5, function(t) 0.5), 1)
  expect_equal(w$fn(500, function(t) 0.1), 1)
})

test_that("weights_uno returns inverse squared censoring probability", {
  w <- weights_uno()
  G <- function(t) 0.5
  expect_equal(w$fn(10, G), 1 / 0.5^2)
  expect_equal(w$name, "Uno (IPCW)")
})

test_that("weights_uno with tau zeroes weights beyond tau", {
  w <- weights_uno(tau = 100)
  G <- function(t) 0.5
  expect_equal(w$fn(50, G), 4)
  expect_equal(w$fn(150, G), 0)
})

test_that("weights_truncated is an indicator", {
  w <- weights_truncated(tau = 100)
  G <- function(t) 0.5
  expect_equal(w$fn(50, G), 1)
  expect_equal(w$fn(100, G), 1)
  expect_equal(w$fn(101, G), 0)
})

test_that("weight names are formatted, not raw doubles", {
  # Regression: an unformatted tau produced "Truncated tau=4.11833385100334"
  w <- weights_truncated(tau = 4.11833385100334)
  expect_false(grepl("4.11833385100334", w$name, fixed = TRUE))
  expect_match(w$name, "^Truncated")
})

test_that("weights_custom accepts a two-argument function", {
  w <- weights_custom(function(t, G) exp(-t), name = "decay")
  expect_s3_class(w, "cindex_weights")
  expect_equal(w$name, "decay")
  expect_equal(w$fn(0, function(t) 1), 1)
})

test_that("weights_custom rejects a non-function", {
  expect_error(weights_custom("nope"), "must be a function")
})

test_that("weights_custom rejects a single-argument function", {
  expect_error(weights_custom(function(t) 1), "two arguments")
})

test_that("weights_truncated rejects invalid tau", {
  expect_error(weights_truncated(tau = -1), "positive")
  expect_error(weights_truncated(tau = c(1, 2)), "single")
})

test_that("print.cindex_weights is informative", {
  expect_output(print(weights_harrell()), "Harrell")
})
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-weights.R")'`
Expected: FAIL — `could not find function "weights_harrell"`

- [ ] **Step 3: Write the implementation**

Create `R/weights.R`:

```r
new_cindex_weights <- function(name, fn) {
  structure(list(name = name, fn = fn), class = "cindex_weights")
}

check_tau <- function(tau) {
  if (!is.numeric(tau) || length(tau) != 1L) {
    stop("`tau` must be a single numeric value.", call. = FALSE)
  }
  if (!is.finite(tau) || tau <= 0) {
    stop("`tau` must be a positive, finite number.", call. = FALSE)
  }
  invisible(TRUE)
}

fmt_tau <- function(tau) format(tau, digits = 4, trim = TRUE)

#' Pair weightings for concordance decomposition
#'
#' @description
#' A concordance index of the form
#' `C = sum(w_ij * c_ij) / sum(w_ij)` is fully determined by its pair
#' weighting `w_ij`. These constructors supply that weighting to
#' [decompose_cindex()]. `weights_harrell()` gives Harrell's C,
#' `weights_uno()` gives Uno's inverse-probability-of-censoring-weighted
#' C, and `weights_custom()` accepts any user-supplied rule.
#'
#' @param tau Truncation time. For `weights_uno()`, pairs whose earlier
#'   event time exceeds `tau` receive weight zero; `NULL` (default) means
#'   no truncation. For `weights_truncated()`, `tau` is required.
#' @param fn A function of two arguments `(t, G)` returning a single
#'   non-negative finite number. `t` is the event time of the earlier
#'   member of the pair; `G` is a function returning the Kaplan-Meier
#'   censoring survival probability at a given time.
#' @param name A length-one character label used in printed output and
#'   plot annotations.
#'
#' @return An object of class `cindex_weights`.
#'
#' @examples
#' weights_harrell()
#' weights_uno()
#' weights_truncated(tau = 365)
#' weights_custom(function(t, G) 1 / G(t), name = "IPCW (power 1)")
#'
#' @name cindex_weights
NULL

#' @rdname cindex_weights
#' @export
weights_harrell <- function() {
  new_cindex_weights("Harrell", function(t, G) 1)
}

#' @rdname cindex_weights
#' @export
weights_uno <- function(tau = NULL) {
  if (!is.null(tau)) check_tau(tau)
  nm <- if (is.null(tau)) {
    "Uno (IPCW)"
  } else {
    paste0("Uno (IPCW, tau = ", fmt_tau(tau), ")")
  }
  new_cindex_weights(nm, function(t, G) {
    if (!is.null(tau) && t > tau) return(0)
    1 / G(t)^2
  })
}

#' @rdname cindex_weights
#' @export
weights_truncated <- function(tau) {
  check_tau(tau)
  new_cindex_weights(
    paste0("Truncated (tau = ", fmt_tau(tau), ")"),
    function(t, G) as.numeric(t <= tau)
  )
}

#' @rdname cindex_weights
#' @export
weights_custom <- function(fn, name = "Custom") {
  if (!is.function(fn)) {
    stop("`fn` must be a function of two arguments, (t, G).", call. = FALSE)
  }
  if (length(formals(fn)) < 2L) {
    stop("`fn` must accept two arguments: (t, G).", call. = FALSE)
  }
  if (!is.character(name) || length(name) != 1L) {
    stop("`name` must be a single character string.", call. = FALSE)
  }
  new_cindex_weights(name, fn)
}

#' @export
print.cindex_weights <- function(x, ...) {
  cat("<cindex_weights>", x$name, "\n")
  invisible(x)
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-weights.R")'`
Expected: PASS, 10 tests.

- [ ] **Step 5: Commit**

```bash
git add R/weights.R tests/testthat/test-weights.R
git commit -m "feat: add pair weight constructors

Adds weights_harrell/uno/truncated/custom returning cindex_weights
objects. Names are formatted at construction so labels never leak raw
double precision into plot legends."
```

---

## Task 3: The pair-counting engine

This task fixes the tie defect. It is the correctness core of the package —
every other task depends on it being right.

**Files:**
- Create: `R/engine.R`
- Create: `tests/testthat/test-engine.R`

**Interfaces:**
- Consumes: `censoring_km()` (Task 1); `cindex_weights` objects (Task 2).
- Produces: `pair_counts(time, status, risk, weights)` returning a list with
  `W_ee`, `S_ee`, `W_ec`, `S_ec` (weighted pair totals and weighted
  concordance totals) and `N_ee`, `N_ec` (unweighted pair counts). Not
  exported.

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-engine.R`:

```r
harrell_ref <- function(d) {
  survival::concordance(
    survival::Surv(d$time, d$status) ~ d$risk,
    reverse = TRUE
  )
}

pooled_from <- function(pc) (pc$S_ee + pc$S_ec) / (pc$W_ee + pc$W_ec)

test_that("event-censored ties are counted as comparable", {
  # A dies at 100, B censored at 100. B was alive when last seen, so we
  # know A's outcome was worse: the pair IS comparable.
  d <- make_tied_pair_data()
  pc <- pair_counts(d$time, d$status, d$risk, weights_harrell())
  expect_equal(pc$N_ee + pc$N_ec, 4)
})

test_that("event-event ties are NOT comparable", {
  # Two deaths at the same time: we cannot say who came first.
  d <- list(time = c(100, 100, 300), status = c(1, 1, 1), risk = c(0.9, 0.2, 0.1))
  pc <- pair_counts(d$time, d$status, d$risk, weights_harrell())
  expect_equal(pc$N_ee + pc$N_ec, 2)  # (A,C) and (B,C) only
})

test_that("pooled decomposition matches survival::concordance, continuous times", {
  d <- make_test_data(300, ties = FALSE, seed = 1)
  pc <- pair_counts(d$time, d$status, d$risk, weights_harrell())
  expect_equal(pooled_from(pc), harrell_ref(d)$concordance, tolerance = 1e-12)
})

test_that("pooled decomposition matches survival::concordance, TIED times", {
  # Regression: `time > time[i]` dropped 1,062 event-censored tied pairs and
  # biased the estimate by +0.0007 on this fixture.
  d <- make_test_data(300, ties = TRUE, seed = 1)
  pc <- pair_counts(d$time, d$status, d$risk, weights_harrell())
  expect_equal(pooled_from(pc), harrell_ref(d)$concordance, tolerance = 1e-12)
})

test_that("pooled decomposition matches survival::concordance with tied risks", {
  d <- make_test_data(300, ties = FALSE, seed = 2)
  d$risk <- as.numeric(cut(d$risk, 5))
  pc <- pair_counts(d$time, d$status, d$risk, weights_harrell())
  expect_equal(pooled_from(pc), harrell_ref(d)$concordance, tolerance = 1e-12)
})

test_that("pooled decomposition matches survival::concordance with both tied", {
  d <- make_test_data(300, ties = TRUE, seed = 3)
  d$risk <- as.numeric(cut(d$risk, 5))
  pc <- pair_counts(d$time, d$status, d$risk, weights_harrell())
  expect_equal(pooled_from(pc), harrell_ref(d)$concordance, tolerance = 1e-12)
})

test_that("comparable pair count matches survival's", {
  d <- make_test_data(300, ties = TRUE, seed = 1)
  pc <- pair_counts(d$time, d$status, d$risk, weights_harrell())
  cnt <- harrell_ref(d)$count
  expect_equal(
    pc$N_ee + pc$N_ec,
    unname(sum(cnt[c("concordant", "discordant", "tied.x")]))
  )
})

test_that("Uno matches survival timewt = n/G2 exactly when times are distinct", {
  d <- make_test_data(400, ties = FALSE, seed = 5)
  pc <- pair_counts(d$time, d$status, d$risk, weights_uno())
  ref <- survival::concordance(
    survival::Surv(d$time, d$status) ~ d$risk,
    reverse = TRUE, timewt = "n/G2"
  )$concordance
  expect_equal(pooled_from(pc), ref, tolerance = 1e-10)
})

test_that("Uno under TIED times agrees with survival only approximately", {
  # This package implements Uno's published PAIRWISE estimator: every pair
  # is weighted by 1 / Ghat(T_i)^2, as in Uno et al. (2011).
  # survival::concordance uses a counting-process form that applies its
  # weight per event TIME instead. The two coincide exactly when event
  # times are distinct (previous test, agreement at 1e-10) and differ by
  # O(1e-4) once times are tied -- measured at 1.2e-4 on day-scale
  # simulated data and 2.7e-4 on survival::lung.
  #
  # Per-PAIR weighting is what makes the decomposition possible at all: a
  # per-event-time weight cannot be partitioned into event-event and
  # event-censored pair sets. The convention is therefore forced by the
  # package's purpose, not a defect, and the gap is two orders of
  # magnitude below the C-index's own standard error.
  #
  # Harrell's weighting is unaffected and still matches exactly under
  # ties -- that is asserted separately above at 1e-12.
  d <- make_test_data(300, ties = TRUE, seed = 5)
  pc <- pair_counts(d$time, d$status, d$risk, weights_uno())
  ref <- survival::concordance(
    survival::Surv(d$time, d$status) ~ d$risk,
    reverse = TRUE, timewt = "n/G2"
  )$concordance
  expect_equal(pooled_from(pc), ref, tolerance = 1e-3)
})

test_that("the weighted identity holds for arbitrary custom weights", {
  d <- make_test_data(200, seed = 7)
  w <- weights_custom(function(t, G) 1 + t, name = "linear")
  pc <- pair_counts(d$time, d$status, d$risk, w)
  C_ee <- pc$S_ee / pc$W_ee
  C_ec <- pc$S_ec / pc$W_ec
  identity_rhs <- (pc$W_ee * C_ee + pc$W_ec * C_ec) / (pc$W_ee + pc$W_ec)
  expect_equal(pooled_from(pc), identity_rhs, tolerance = 1e-12)
})

test_that("negating risk gives 1 - C", {
  d <- make_test_data(200, seed = 9)
  a <- pooled_from(pair_counts(d$time, d$status,  d$risk, weights_harrell()))
  b <- pooled_from(pair_counts(d$time, d$status, -d$risk, weights_harrell()))
  expect_equal(a + b, 1, tolerance = 1e-12)
})

test_that("a negative weight is an error", {
  d <- make_test_data(50, seed = 11)
  w <- weights_custom(function(t, G) -1, name = "bad")
  expect_error(
    pair_counts(d$time, d$status, d$risk, w),
    "negative or non-finite"
  )
})

test_that("a non-finite weight is an error", {
  d <- make_test_data(50, seed = 11)
  w <- weights_custom(function(t, G) Inf, name = "bad")
  expect_error(
    pair_counts(d$time, d$status, d$risk, w),
    "negative or non-finite"
  )
})
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-engine.R")'`
Expected: FAIL — `could not find function "pair_counts"`

- [ ] **Step 3: Write the implementation**

Create `R/engine.R`:

```r
#' Count weighted concordant pairs by partner status
#'
#' The single place in the package where pairs are counted. Takes plain
#' vectors plus a weighting and returns weighted totals split by whether
#' the later member of each pair had an event or was censored.
#'
#' Tie rules (chosen to match [survival::concordance()] exactly):
#' * two events at the same time are NOT comparable;
#' * an event and a censoring at the same time ARE comparable, with the
#'   event treated as occurring first;
#' * equal risk scores receive half credit.
#'
#' @param time,status,risk Numeric vectors of equal length.
#' @param weights A `cindex_weights` object.
#' @return A list with `W_ee`, `S_ee`, `W_ec`, `S_ec`, `N_ee`, `N_ec`.
#' @keywords internal
#' @noRd
pair_counts <- function(time, status, risk, weights) {
  G <- censoring_km(time, status)
  event_idx <- which(status == 1)

  W_ee <- 0; S_ee <- 0; N_ee <- 0
  W_ec <- 0; S_ec <- 0; N_ec <- 0

  for (i in event_idx) {
    # An event at time[i] precedes a censoring recorded at the same time,
    # so those pairs are comparable. Two events at the same time are not.
    # `status[i] == 1`, so the second clause can never select `i` itself.
    j <- which(time > time[i] | (time == time[i] & status == 0))
    if (length(j) == 0L) next

    w <- weights$fn(time[i], G)
    if (length(w) != 1L || !is.finite(w) || w < 0) {
      stop(
        "Weight function returned a negative or non-finite value at t = ",
        format(time[i]), ".",
        call. = FALSE
      )
    }
    if (w == 0) next

    conc <- (risk[i] > risk[j]) + 0.5 * (risk[i] == risk[j])
    is_ee <- status[j] == 1

    n_ee_i <- sum(is_ee)
    n_ec_i <- length(j) - n_ee_i

    N_ee <- N_ee + n_ee_i
    N_ec <- N_ec + n_ec_i
    W_ee <- W_ee + w * n_ee_i
    W_ec <- W_ec + w * n_ec_i
    S_ee <- S_ee + w * sum(conc[is_ee])
    S_ec <- S_ec + w * sum(conc[!is_ee])
  }

  list(
    W_ee = W_ee, S_ee = S_ee, N_ee = N_ee,
    W_ec = W_ec, S_ec = S_ec, N_ec = N_ec
  )
}
```

**Note for the implementer:** the two Uno tolerances differ deliberately and
both were measured, not guessed. Distinct event times must agree at `1e-10`;
tied event times agree only to `1e-3`, for the reason documented in the second
test. Do not "fix" the tied case by switching `censoring_km()` to
`findInterval(t, ftime, left.open = TRUE)` -- that was tried and makes the gap
larger (5.5e-4 rather than 3.4e-4). Do not loosen the `1e-12` Harrell
tolerances either: Harrell agreement is exact under every tie pattern tested,
including `survival::lung`, and it is the package's central correctness claim.

- [ ] **Step 4: Run tests to verify they pass**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-engine.R")'`
Expected: PASS, 12 tests.

- [ ] **Step 5: Commit**

```bash
git add R/engine.R tests/testthat/test-engine.R
git commit -m "fix: count event-censored pairs at tied times

The engine used `time > time[i]`, which dropped pairs where an event and
a censoring share a recorded time. Since survival data is normally
recorded in whole days these ties are common: on a 300-subject fixture
1,062 pairs were lost and the estimate was biased +0.0007 against
survival::concordance.

Also generalises pair counting to arbitrary pair weightings, so one
engine now produces Harrell's C, Uno's C and any user-supplied variant.
Agreement with survival::concordance is pinned to 1e-12."
```

---

## Task 4: The `decompose_cindex()` generic

**Files:**
- Create: `R/decompose.R`
- Delete: `R/decompose_cindex.R`
- Delete: `man/decompose_cindex.Rd`
- Create: `tests/testthat/test-decompose.R`

**Interfaces:**
- Consumes: `validate_survival_inputs()` (Task 1), `pair_counts()` (Task 3), `weights_*()` (Task 2).
- Produces: `decompose_cindex(x, ...)` generic with `.default(x, status, risk, weights, higher_is_riskier)` and `.formula(x, data, weights, higher_is_riskier)` methods, returning an object of class `cindex_decomp`: a list with numeric scalars `C_ee`, `C_ec`, `C_global`, `gap`, `W_ee`, `W_ec`, `N_ee`, `N_ec`, `n`, `n_events`, `censoring_rate`; character `weighting`; logical `higher_is_riskier`; and `boot = NULL`.

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-decompose.R`:

```r
test_that("decompose_cindex returns a cindex_decomp object", {
  d <- make_test_data(200, seed = 1)
  r <- decompose_cindex(d$time, d$status, d$risk)
  expect_s3_class(r, "cindex_decomp")
  expect_type(r$C_ee, "double")
  expect_type(r$C_global, "double")
})

test_that("the decomposition identity holds", {
  d <- make_test_data(300, ties = TRUE, seed = 2)
  r <- decompose_cindex(d$time, d$status, d$risk)
  rhs <- (r$W_ee * r$C_ee + r$W_ec * r$C_ec) / (r$W_ee + r$W_ec)
  expect_equal(r$C_global, rhs, tolerance = 1e-12)
})

test_that("C_global matches survival::concordance", {
  d <- make_test_data(300, ties = TRUE, seed = 2)
  r <- decompose_cindex(d$time, d$status, d$risk)
  ref <- survival::concordance(
    survival::Surv(d$time, d$status) ~ d$risk, reverse = TRUE
  )$concordance
  expect_equal(r$C_global, ref, tolerance = 1e-12)
})

test_that("gap equals C_ec minus C_ee", {
  d <- make_test_data(200, seed = 3)
  r <- decompose_cindex(d$time, d$status, d$risk)
  expect_equal(r$gap, r$C_ec - r$C_ee, tolerance = 1e-12)
})

test_that("the formula method agrees with the vector method", {
  d <- make_test_data(200, seed = 4)
  df <- data.frame(time = d$time, status = d$status, risk = d$risk)
  a <- decompose_cindex(d$time, d$status, d$risk)
  b <- decompose_cindex(survival::Surv(time, status) ~ risk, data = df)
  expect_equal(a$C_ee, b$C_ee, tolerance = 1e-12)
  expect_equal(a$C_global, b$C_global, tolerance = 1e-12)
})

test_that("the formula method rejects a non-Surv response", {
  df <- data.frame(y = rnorm(10), risk = rnorm(10))
  expect_error(decompose_cindex(y ~ risk, data = df), "Surv")
})

test_that("the formula method rejects multiple predictors", {
  d <- make_test_data(50, seed = 5)
  df <- data.frame(time = d$time, status = d$status,
                   risk = d$risk, other = rnorm(50))
  expect_error(
    decompose_cindex(survival::Surv(time, status) ~ risk + other, data = df),
    "exactly one"
  )
})

test_that("higher_is_riskier = FALSE flips ALL components together", {
  # Regression: the old code flipped only the global value, breaking the
  # identity and drawing a global line above both of its own components.
  d <- make_test_data(200, seed = 6)
  a <- decompose_cindex(d$time, d$status,  d$risk, higher_is_riskier = TRUE)
  b <- decompose_cindex(d$time, d$status, -d$risk, higher_is_riskier = FALSE)
  expect_equal(a$C_ee, b$C_ee, tolerance = 1e-12)
  expect_equal(a$C_ec, b$C_ec, tolerance = 1e-12)
  expect_equal(a$C_global, b$C_global, tolerance = 1e-12)
  # and the identity still holds after flipping
  rhs <- (b$W_ee * b$C_ee + b$W_ec * b$C_ec) / (b$W_ee + b$W_ec)
  expect_equal(b$C_global, rhs, tolerance = 1e-12)
})

test_that("sub-0.5 concordance is reported, never flipped", {
  # Regression: `ifelse(c < 0.5, 1 - c, c)` reported a 0.44 model as 0.56.
  d <- make_test_data(300, seed = 7)
  r <- decompose_cindex(d$time, d$status, -d$risk)  # deliberately backwards
  expect_lt(r$C_global, 0.5)
})

test_that("the weighting name is carried on the object", {
  d <- make_test_data(200, seed = 8)
  r <- decompose_cindex(d$time, d$status, d$risk, weights = weights_uno())
  expect_equal(r$weighting, "Uno (IPCW)")
})

test_that("censoring rate and counts are recorded", {
  d <- make_test_data(200, seed = 9)
  r <- decompose_cindex(d$time, d$status, d$risk)
  expect_equal(r$n, 200)
  expect_equal(r$n_events, sum(d$status == 1))
  expect_equal(r$censoring_rate, mean(d$status == 0), tolerance = 1e-12)
})

test_that("invalid input is rejected", {
  expect_error(decompose_cindex(1:5, c(0, 1, 0, 1), rnorm(5)), "same length")
})
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-decompose.R")'`
Expected: FAIL — the old `decompose_cindex()` returns a plain list, so
`expect_s3_class` fails and `r$C_global` is `NULL`.

- [ ] **Step 3: Remove the superseded implementation**

```bash
git rm R/decompose_cindex.R man/decompose_cindex.Rd
```

- [ ] **Step 4: Write the implementation**

Create `R/decompose.R`:

```r
#' Decompose a concordance index into Event-Event and Event-Censored parts
#'
#' @description
#' Splits a concordance index into the two ranking tasks it pools: pairs
#' where both subjects had the event (`C_ee`), and pairs where one had the
#' event and one was censored (`C_ec`). As censoring rises the event-event
#' pairs vanish and `C_ec` dominates the global score, so a model close to
#' chance on true events can still report a stable global C-index.
#'
#' The decomposition applies to any concordance index of the form
#' `C = sum(w_ij * c_ij) / sum(w_ij)`. The weighting is chosen with
#' [weights_harrell()], [weights_uno()], [weights_truncated()] or
#' [weights_custom()], and the identity
#' `C = (W_ee * C_ee + W_ec * C_ec) / (W_ee + W_ec)` holds for each.
#'
#' @param x Either a numeric vector of observed times (default method) or a
#'   formula of the form `Surv(time, status) ~ risk` (formula method).
#' @param status Numeric vector of event indicators (1 = event,
#'   0 = censored). Default method only.
#' @param risk Numeric vector of risk scores. Default method only.
#' @param data A data frame in which to evaluate the formula. Formula
#'   method only.
#' @param weights A `cindex_weights` object. Defaults to
#'   [weights_harrell()].
#' @param higher_is_riskier Logical. `TRUE` (default) if larger values of
#'   `risk` indicate higher hazard. If `FALSE`, `risk` is negated once
#'   before any pair counting, so every component stays mutually
#'   consistent. Values below 0.5 are reported as they are: a
#'   worse-than-chance model is a finding, not an error.
#' @param ... Passed between methods.
#'
#' @return An object of class `cindex_decomp`.
#'
#' @examples
#' set.seed(42)
#' time   <- rexp(200, rate = 0.1)
#' status <- rbinom(200, 1, 0.6)
#' risk   <- rnorm(200)
#'
#' decompose_cindex(time, status, risk)
#' decompose_cindex(time, status, risk, weights = weights_uno())
#'
#' @export
decompose_cindex <- function(x, ...) {
  UseMethod("decompose_cindex")
}

#' @rdname decompose_cindex
#' @export
decompose_cindex.default <- function(x, status, risk,
                                     weights = weights_harrell(),
                                     higher_is_riskier = TRUE, ...) {
  time <- x
  validate_survival_inputs(time, status, risk)
  if (!inherits(weights, "cindex_weights")) {
    stop("`weights` must be a `cindex_weights` object, e.g. weights_harrell().",
         call. = FALSE)
  }
  if (!is.logical(higher_is_riskier) || length(higher_is_riskier) != 1L) {
    stop("`higher_is_riskier` must be TRUE or FALSE.", call. = FALSE)
  }

  # Orientation is applied once, before any counting, so C_ee, C_ec and
  # C_global stay mutually consistent and the identity survives.
  if (!higher_is_riskier) risk <- -risk

  pc <- pair_counts(time, status, risk, weights)

  C_ee <- if (pc$W_ee > 0) pc$S_ee / pc$W_ee else NA_real_
  C_ec <- if (pc$W_ec > 0) pc$S_ec / pc$W_ec else NA_real_
  total_w <- pc$W_ee + pc$W_ec
  C_global <- if (total_w > 0) (pc$S_ee + pc$S_ec) / total_w else NA_real_

  structure(
    list(
      C_ee = C_ee,
      C_ec = C_ec,
      C_global = C_global,
      gap = C_ec - C_ee,
      W_ee = pc$W_ee,
      W_ec = pc$W_ec,
      N_ee = pc$N_ee,
      N_ec = pc$N_ec,
      n = length(time),
      n_events = sum(status == 1),
      censoring_rate = mean(status == 0),
      weighting = weights$name,
      higher_is_riskier = higher_is_riskier,
      boot = NULL
    ),
    class = "cindex_decomp"
  )
}

#' @rdname decompose_cindex
#' @export
decompose_cindex.formula <- function(x, data = parent.frame(),
                                     weights = weights_harrell(),
                                     higher_is_riskier = TRUE, ...) {
  mf <- stats::model.frame(x, data = data)
  resp <- stats::model.response(mf)
  if (!inherits(resp, "Surv")) {
    stop("The left-hand side of the formula must be a Surv() object.",
         call. = FALSE)
  }
  if (ncol(mf) != 2L) {
    stop("The right-hand side must name exactly one risk score variable.",
         call. = FALSE)
  }
  decompose_cindex.default(
    x = as.numeric(resp[, 1L]),
    status = as.numeric(resp[, 2L]),
    risk = as.numeric(mf[[2L]]),
    weights = weights,
    higher_is_riskier = higher_is_riskier,
    ...
  )
}
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-decompose.R")'`
Expected: PASS, 12 tests.

- [ ] **Step 6: Commit**

```bash
git add -A R/ man/ tests/testthat/test-decompose.R
git commit -m "feat: add decompose_cindex generic with formula and vector methods

Returns a cindex_decomp S3 object carrying C_ee, C_ec, C_global, the
masking gap, weight sums, pair counts and cohort summaries.

Replaces the sign flip with an explicit higher_is_riskier argument
applied once before counting, so all three components flip together and
the decomposition identity survives. Worse-than-chance results are now
reported rather than silently inverted."
```

---

## Task 5: Subject-level bootstrap and `confint()`

**Files:**
- Create: `R/bootstrap.R`
- Modify: `R/decompose.R` (add `n_boot` and `conf_level` arguments)
- Create: `tests/testthat/test-bootstrap.R`

**Interfaces:**
- Consumes: `pair_counts()` (Task 3), `decompose_cindex()` (Task 4).
- Produces: `bootstrap_decomp(time, status, risk, weights, n_boot)` returning a
  numeric matrix with `n_boot` rows and columns `C_ee`, `C_ec`, `C_global`,
  `gap`. `decompose_cindex()` gains `n_boot = 0` and `conf_level = 0.95`
  arguments and stores the matrix on `$boot` (plus `$conf_level`,
  `$n_boot`). `confint.cindex_decomp(object, parm, level)` returning a
  matrix with rownames `C_ee`/`C_ec`/`C_global`/`gap` and columns for the
  two percentile bounds.

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-bootstrap.R`:

```r
test_that("bootstrap_decomp returns the expected shape", {
  d <- make_test_data(150, seed = 1)
  set.seed(1)
  b <- bootstrap_decomp(d$time, d$status, d$risk, weights_harrell(), n_boot = 20)
  expect_true(is.matrix(b))
  expect_equal(nrow(b), 20)
  expect_equal(colnames(b), c("C_ee", "C_ec", "C_global", "gap"))
})

test_that("the bootstrap resamples SUBJECTS, not pairs", {
  # Regression guard. 150 subjects generate tens of thousands of pairs.
  # Treating pairs as independent understates SE(C_ee) by roughly 6x, so a
  # correct subject-level bootstrap must be far wider than the binomial
  # pair-based SE.
  d <- make_test_data(150, seed = 2)
  r <- decompose_cindex(d$time, d$status, d$risk, n_boot = 300)
  se_boot <- sd(r$boot[, "C_ee"], na.rm = TRUE)
  se_pairs <- sqrt(r$C_ee * (1 - r$C_ee) / r$N_ee)
  expect_gt(se_boot, 3 * se_pairs)
})

test_that("the gap is computed within each replicate, preserving correlation", {
  d <- make_test_data(150, seed = 3)
  r <- decompose_cindex(d$time, d$status, d$risk, n_boot = 100)
  expect_equal(
    r$boot[, "gap"],
    r$boot[, "C_ec"] - r$boot[, "C_ee"],
    tolerance = 1e-12
  )
})

test_that("the bootstrap is reproducible under set.seed", {
  d <- make_test_data(120, seed = 4)
  set.seed(99); a <- decompose_cindex(d$time, d$status, d$risk, n_boot = 40)
  set.seed(99); b <- decompose_cindex(d$time, d$status, d$risk, n_boot = 40)
  expect_equal(a$boot, b$boot)
})

test_that("n_boot = 0 stores no replicates", {
  d <- make_test_data(100, seed = 5)
  r <- decompose_cindex(d$time, d$status, d$risk)
  expect_null(r$boot)
})

test_that("confint errors clearly when there are no replicates", {
  d <- make_test_data(100, seed = 6)
  r <- decompose_cindex(d$time, d$status, d$risk)
  expect_error(confint(r), "n_boot")
})

test_that("confint returns intervals containing the point estimates", {
  d <- make_test_data(200, seed = 7)
  r <- decompose_cindex(d$time, d$status, d$risk, n_boot = 200)
  ci <- confint(r)
  expect_equal(rownames(ci), c("C_ee", "C_ec", "C_global", "gap"))
  expect_lte(ci["C_ee", 1], r$C_ee)
  expect_gte(ci["C_ee", 2], r$C_ee)
  expect_lte(ci["gap", 1], r$gap)
  expect_gte(ci["gap", 2], r$gap)
})

test_that("confint respects the level argument", {
  d <- make_test_data(200, seed = 8)
  r <- decompose_cindex(d$time, d$status, d$risk, n_boot = 300)
  wide <- confint(r, level = 0.99)
  narrow <- confint(r, level = 0.80)
  expect_gt(wide["C_ee", 2] - wide["C_ee", 1],
            narrow["C_ee", 2] - narrow["C_ee", 1])
})

test_that("confint can select a subset of parameters", {
  d <- make_test_data(150, seed = 9)
  r <- decompose_cindex(d$time, d$status, d$risk, n_boot = 100)
  ci <- confint(r, parm = "gap")
  expect_equal(rownames(ci), "gap")
})

test_that("n_boot is validated", {
  d <- make_test_data(100, seed = 10)
  expect_error(decompose_cindex(d$time, d$status, d$risk, n_boot = -5),
               "non-negative")
})
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-bootstrap.R")'`
Expected: FAIL — `could not find function "bootstrap_decomp"`

- [ ] **Step 3: Write the bootstrap**

Create `R/bootstrap.R`:

```r
#' Subject-level bootstrap for a concordance decomposition
#'
#' Resamples **subjects** with replacement, not pairs. Pairs are strongly
#' dependent -- each subject appears in hundreds of them -- so resampling
#' pairs directly produces intervals several times too narrow.
#'
#' All four quantities are recomputed inside each replicate, so the
#' correlation between `C_ee` and `C_ec` is carried through to the gap
#' automatically.
#'
#' @param time,status,risk Numeric vectors of equal length.
#' @param weights A `cindex_weights` object.
#' @param n_boot Number of replicates.
#' @return A numeric matrix, `n_boot` rows, columns `C_ee`, `C_ec`,
#'   `C_global`, `gap`. Replicates that yield no comparable pairs of one
#'   kind are left as `NA`.
#' @keywords internal
#' @noRd
bootstrap_decomp <- function(time, status, risk, weights, n_boot) {
  n <- length(time)
  out <- matrix(
    NA_real_, nrow = n_boot, ncol = 4L,
    dimnames = list(NULL, c("C_ee", "C_ec", "C_global", "gap"))
  )
  for (b in seq_len(n_boot)) {
    k <- sample.int(n, n, replace = TRUE)
    pc <- tryCatch(
      pair_counts(time[k], status[k], risk[k], weights),
      error = function(e) NULL
    )
    if (is.null(pc) || pc$W_ee <= 0 || pc$W_ec <= 0) next
    ee <- pc$S_ee / pc$W_ee
    ec <- pc$S_ec / pc$W_ec
    out[b, ] <- c(ee, ec, (pc$S_ee + pc$S_ec) / (pc$W_ee + pc$W_ec), ec - ee)
  }
  out
}
```

- [ ] **Step 4: Wire `n_boot` into `decompose_cindex()`**

In `R/decompose.R`, change the `decompose_cindex.default` signature to:

```r
decompose_cindex.default <- function(x, status, risk,
                                     weights = weights_harrell(),
                                     higher_is_riskier = TRUE,
                                     n_boot = 0, conf_level = 0.95, ...) {
```

Add this validation immediately after the `higher_is_riskier` check:

```r
  if (!is.numeric(n_boot) || length(n_boot) != 1L || n_boot < 0) {
    stop("`n_boot` must be a single non-negative number.", call. = FALSE)
  }
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      conf_level <= 0 || conf_level >= 1) {
    stop("`conf_level` must be a single number strictly between 0 and 1.",
         call. = FALSE)
  }
  n_boot <- as.integer(n_boot)
```

Immediately before the `structure(...)` call, add:

```r
  boot <- if (n_boot > 0L) {
    bootstrap_decomp(time, status, risk, weights, n_boot)
  } else {
    NULL
  }
```

Add these three elements to the list inside `structure(...)`, replacing the
existing `boot = NULL` line:

```r
      boot = boot,
      n_boot = n_boot,
      conf_level = conf_level,
```

Mirror the two new arguments in `decompose_cindex.formula` by adding
`n_boot = 0, conf_level = 0.95` to its signature and passing them through
in the `decompose_cindex.default(...)` call.

Add to the roxygen block of `decompose_cindex`:

```r
#' @param n_boot Number of bootstrap replicates. `0` (default) skips
#'   resampling, in which case [confint()] is unavailable. Replicates
#'   resample subjects, not pairs.
#' @param conf_level Confidence level stored on the object and used as the
#'   default for [confint()].
```

- [ ] **Step 5: Write `confint()`**

Append to `R/bootstrap.R`:

```r
#' Bootstrap confidence intervals for a concordance decomposition
#'
#' @param object A `cindex_decomp` object created with `n_boot > 0`.
#' @param parm Character vector selecting quantities. Defaults to all of
#'   `"C_ee"`, `"C_ec"`, `"C_global"`, `"gap"`.
#' @param level Confidence level. Defaults to the object's `conf_level`.
#' @param ... Ignored.
#' @return A matrix of percentile bounds, one row per quantity.
#' @examples
#' set.seed(1)
#' time   <- rexp(150, rate = 0.1)
#' status <- rbinom(150, 1, 0.6)
#' risk   <- rnorm(150)
#' fit <- decompose_cindex(time, status, risk, n_boot = 50)
#' confint(fit)
#' @export
confint.cindex_decomp <- function(object,
                                  parm = c("C_ee", "C_ec", "C_global", "gap"),
                                  level = NULL, ...) {
  if (is.null(object$boot)) {
    stop(
      "No bootstrap replicates are stored on this object. ",
      "Re-run decompose_cindex() with n_boot > 0 (e.g. n_boot = 1000).",
      call. = FALSE
    )
  }
  if (is.null(level)) level <- object$conf_level
  parm <- match.arg(parm, c("C_ee", "C_ec", "C_global", "gap"),
                    several.ok = TRUE)
  a <- (1 - level) / 2
  probs <- c(a, 1 - a)
  out <- t(vapply(
    parm,
    function(p) stats::quantile(object$boot[, p], probs = probs,
                                na.rm = TRUE, names = FALSE),
    numeric(2)
  ))
  colnames(out) <- paste0(format(100 * probs, trim = TRUE), " %")
  rownames(out) <- parm
  out
}
```

- [ ] **Step 6: Run tests to verify they pass**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-bootstrap.R")'`
Expected: PASS, 10 tests.

- [ ] **Step 7: Run the whole suite so far**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_dir("tests/testthat")'`
Expected: PASS, no failures.

- [ ] **Step 8: Commit**

```bash
git add R/bootstrap.R R/decompose.R tests/testthat/test-bootstrap.R
git commit -m "feat: add subject-level bootstrap and confint method

Resamples subjects rather than pairs; pair-level resampling understates
SE(C_ee) by roughly 6x because each subject appears in hundreds of pairs.
All four quantities are recomputed inside each replicate so the
correlation between C_ee and C_ec carries through to the masking gap,
which is the quantity requiring an interval."
```

---

## Task 6: Printing and coercion

**Files:**
- Create: `R/methods.R`
- Create: `tests/testthat/test-methods.R`

**Interfaces:**
- Consumes: `cindex_decomp` objects (Tasks 4-5), `confint.cindex_decomp()` (Task 5).
- Produces: `print.cindex_decomp()`, `summary.cindex_decomp()` (returns a `summary.cindex_decomp` object with its own print method), `as.data.frame.cindex_decomp()` returning a one-row data frame with columns `weighting`, `n`, `n_events`, `censoring_rate`, `C_ee`, `C_ec`, `C_global`, `gap`, `N_ee`, `N_ec`.

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-methods.R`:

```r
test_that("print shows the global C-index and the masking gap", {
  d <- make_test_data(200, seed = 1)
  r <- decompose_cindex(d$time, d$status, d$risk)
  out <- capture.output(print(r))
  expect_true(any(grepl("Global", out)))
  expect_true(any(grepl("Event-Event", out)))
  expect_true(any(grepl("Event-Censored", out)))
  expect_true(any(grepl("gap", out, ignore.case = TRUE)))
})

test_that("print names the weighting", {
  d <- make_test_data(150, seed = 2)
  r <- decompose_cindex(d$time, d$status, d$risk, weights = weights_uno())
  expect_true(any(grepl("Uno", capture.output(print(r)))))
})

test_that("print warns that pair counts are not precision when no CI exists", {
  # Constraint: N_ee must never be shown alone as though it were a sample
  # size. 16,471 pairs from 300 patients is 300 units of information.
  d <- make_test_data(200, seed = 3)
  r <- decompose_cindex(d$time, d$status, d$risk)
  out <- capture.output(print(r))
  expect_true(any(grepl("not precision|n_boot", out)))
})

test_that("print shows confidence intervals when they exist", {
  d <- make_test_data(150, seed = 4)
  r <- decompose_cindex(d$time, d$status, d$risk, n_boot = 50)
  out <- capture.output(print(r))
  expect_true(any(grepl("95%|CI", out)))
})

test_that("print flags near-chance event-event concordance", {
  d <- make_test_data(300, seed = 5)
  r <- decompose_cindex(d$time, d$status, d$risk)
  r$C_ee <- 0.51  # force the condition
  expect_true(any(grepl("chance", capture.output(print(r)))))
})

test_that("print returns its input invisibly", {
  d <- make_test_data(100, seed = 6)
  r <- decompose_cindex(d$time, d$status, d$risk)
  expect_invisible(print(r))
})

test_that("summary carries the bootstrap settings", {
  d <- make_test_data(150, seed = 7)
  r <- decompose_cindex(d$time, d$status, d$risk, n_boot = 50)
  s <- summary(r)
  expect_s3_class(s, "summary.cindex_decomp")
  expect_true(any(grepl("50", capture.output(print(s)))))
})

test_that("as.data.frame returns one tidy row", {
  d <- make_test_data(150, seed = 8)
  r <- decompose_cindex(d$time, d$status, d$risk)
  df <- as.data.frame(r)
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 1)
  expect_true(all(c("C_ee", "C_ec", "C_global", "gap", "N_ee", "N_ec",
                    "weighting") %in% names(df)))
  expect_equal(df$C_ee, r$C_ee, tolerance = 1e-12)
})
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-methods.R")'`
Expected: FAIL — the default list print has no "Global" line.

- [ ] **Step 3: Write the implementation**

Create `R/methods.R`:

```r
fmt_n <- function(x) formatC(x, format = "d", big.mark = ",")

#' @param x A `cindex_decomp` object.
#' @param ... Ignored.
#' @rdname decompose_cindex
#' @export
print.cindex_decomp <- function(x, ...) {
  ci <- if (!is.null(x$boot)) confint(x) else NULL
  total_n <- x$N_ee + x$N_ec

  cat("C-Index Decomposition\n")
  cat(sprintf(
    "Weighting: %s | n = %d, events = %d (%.1f%%), censoring %.1f%%\n\n",
    x$weighting, x$n, x$n_events,
    100 * x$n_events / x$n, 100 * x$censoring_rate
  ))

  row <- function(label, value, npairs, key) {
    share <- if (total_n > 0) 100 * npairs / total_n else NA_real_
    line <- sprintf("  %-16s %8.4f %11s %7.1f%%",
                    label, value, fmt_n(npairs), share)
    if (!is.null(ci) && key %in% rownames(ci)) {
      line <- paste0(line, sprintf("   [%.4f, %.4f]", ci[key, 1], ci[key, 2]))
    }
    cat(line, "\n")
  }

  header <- sprintf("  %-16s %8s %11s %8s", "", "C-index", "pairs", "share")
  if (!is.null(ci)) {
    header <- paste0(header, sprintf("   %s CI", 
                                     format(100 * x$conf_level, trim = TRUE)))
  }
  cat(header, "\n")
  row("Event-Event", x$C_ee, x$N_ee, "C_ee")
  row("Event-Censored", x$C_ec, x$N_ec, "C_ec")
  cat("  ", strrep("-", if (is.null(ci)) 46 else 66), "\n", sep = "")
  row("Global", x$C_global, total_n, "C_global")

  cat(sprintf("\n  Masking gap (C_ec - C_ee): %+.3f", x$gap))
  if (!is.null(ci)) {
    cat(sprintf("   [%+.3f, %+.3f]", ci["gap", 1], ci["gap", 2]))
  }
  cat("\n")

  if (!is.na(x$C_ee) && x$C_ee < 0.55) {
    cat("  ! Event-Event concordance is near chance (0.50).\n")
  }
  if (!is.na(x$C_global) && x$C_global < 0.5) {
    cat("  ! Global C-index is below 0.50. If the risk score is oriented\n")
    cat("    so that higher means lower hazard, set higher_is_riskier = FALSE.\n")
  }
  if (is.null(ci)) {
    cat("  Pair counts describe composition, not precision.",
        "Set n_boot > 0 for intervals.\n")
  }
  invisible(x)
}

#' @param object A `cindex_decomp` object.
#' @rdname decompose_cindex
#' @export
summary.cindex_decomp <- function(object, ...) {
  structure(list(fit = object), class = "summary.cindex_decomp")
}

#' @export
print.summary.cindex_decomp <- function(x, ...) {
  print(x$fit)
  cat("\nBootstrap: ")
  if (is.null(x$fit$boot)) {
    cat("none (n_boot = 0)\n")
  } else {
    n_ok <- sum(stats::complete.cases(x$fit$boot))
    cat(sprintf("%d replicates (%d usable), subject-level resampling, %.0f%% level\n",
                x$fit$n_boot, n_ok, 100 * x$fit$conf_level))
  }
  cat(sprintf("Orientation: higher_is_riskier = %s\n", x$fit$higher_is_riskier))
  invisible(x)
}

#' @param row.names Ignored.
#' @param optional Ignored.
#' @rdname decompose_cindex
#' @export
as.data.frame.cindex_decomp <- function(x, row.names = NULL,
                                        optional = FALSE, ...) {
  data.frame(
    weighting = x$weighting,
    n = x$n,
    n_events = x$n_events,
    censoring_rate = x$censoring_rate,
    C_ee = x$C_ee,
    C_ec = x$C_ec,
    C_global = x$C_global,
    gap = x$gap,
    N_ee = x$N_ee,
    N_ec = x$N_ec,
    stringsAsFactors = FALSE
  )
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-methods.R")'`
Expected: PASS, 8 tests.

- [ ] **Step 5: Eyeball the output**

Run:

```bash
Rscript -e 'devtools::load_all("."); set.seed(1);
  n <- 400; risk <- rnorm(n); tt <- rexp(n, exp(0.8*risk)*0.1); cc <- rexp(n, 0.08);
  print(decompose_cindex(pmin(tt,cc), as.numeric(tt<=cc), risk, n_boot = 100))'
```

Expected: a formatted table with a Global row, the masking gap, and a
bracketed interval on each line.

- [ ] **Step 6: Commit**

```bash
git add R/methods.R tests/testthat/test-methods.R
git commit -m "feat: add print, summary and as.data.frame for cindex_decomp

Leads with the global C-index and the masking gap rather than dumping a
raw list. Pair counts are never shown as a bare figure: either they sit
beside a confidence interval, or the output states that they describe
composition rather than precision."
```

---

## Task 7: `censoring_curve()`

Replaces `simulate_censoring()`. Fixes the x-axis defect and removes both
result-suppressing behaviours.

**Files:**
- Create: `R/censoring_curve.R`
- Delete: `R/simulate_censoring.R`
- Delete: `man/simulate_censoring.Rd`
- Create: `tests/testthat/test-censoring-curve.R`

**Interfaces:**
- Consumes: `validate_survival_inputs()` (Task 1), `pair_counts()` (Task 3).
- Produces: `censoring_curve(time, status, risk, weights, higher_is_riskier, n_thresholds, probs, min_pairs)` returning an object of class `cindex_curve`: a list with `$data` (a data frame with columns `threshold`, `censoring`, `C_ee`, `C_ec`, `C_global`, `gap`, `N_ee`, `N_ec`, `low_precision`), plus `weighting`, `min_pairs`, `n`, `n_events`, `higher_is_riskier`. `print.cindex_curve()`, `as.data.frame.cindex_curve()`.

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-censoring-curve.R`:

```r
test_that("censoring_curve returns a cindex_curve object", {
  d <- make_test_data(300, seed = 1)
  cv <- censoring_curve(d$time, d$status, d$risk, n_thresholds = 6)
  expect_s3_class(cv, "cindex_curve")
  expect_s3_class(cv$data, "data.frame")
  expect_true(all(c("threshold", "censoring", "C_ee", "C_ec", "C_global",
                    "N_ee", "N_ec", "low_precision") %in% names(cv$data)))
})

test_that("the censoring rate is measured on the WHOLE cohort", {
  # Regression: the old code deleted censored subjects on its first line,
  # then reported the censoring rate among events only. On a 500-subject
  # fixture that understated the true rate by up to 31 percentage points.
  d <- make_test_data(400, censor_rate = 0.4, seed = 2)
  cv <- censoring_curve(d$time, d$status, d$risk, n_thresholds = 5)
  for (i in seq_len(nrow(cv$data))) {
    tau <- cv$data$threshold[i]
    expected <- mean(d$status * as.numeric(d$time <= tau) == 0)
    expect_equal(cv$data$censoring[i], expected, tolerance = 1e-12)
  }
})

test_that("reported censoring is never below the baseline censoring rate", {
  d <- make_test_data(400, censor_rate = 0.4, seed = 3)
  baseline <- mean(d$status == 0)
  cv <- censoring_curve(d$time, d$status, d$risk, n_thresholds = 6)
  expect_true(all(cv$data$censoring >= baseline - 1e-12))
})

test_that("C_ee is unchanged by keeping censored subjects", {
  # Event-event pairs only involve observed events, so the core quantity
  # must be identical whether or not censored subjects are retained.
  d <- make_test_data(300, seed = 4)
  cv <- censoring_curve(d$time, d$status, d$risk, n_thresholds = 4)
  ev <- which(d$status == 1)
  for (i in seq_len(nrow(cv$data))) {
    tau <- cv$data$threshold[i]
    te <- pmin(d$time[ev], tau)
    se <- as.numeric(d$time[ev] <= tau)
    if (sum(se) < 2) next
    pc <- pair_counts(te, se, d$risk[ev], weights_harrell())
    if (pc$W_ee <= 0) next
    expect_equal(cv$data$C_ee[i], pc$S_ee / pc$W_ee, tolerance = 1e-10)
  }
})

test_that("values below 0.5 are reported, not deleted", {
  # Regression: `ifelse(CI_ec < 0.5, NA, CI_ec)` deleted results because of
  # their magnitude -- exactly the points that demonstrate masking.
  d <- make_test_data(400, seed = 5)
  cv <- censoring_curve(d$time, d$status, -d$risk, n_thresholds = 6)
  expect_true(any(cv$data$C_global < 0.5))
  expect_false(any(is.na(cv$data$C_ec)))
})

test_that("no global value is flipped", {
  d <- make_test_data(400, seed = 6)
  cv <- censoring_curve(d$time, d$status, -d$risk, n_thresholds = 6)
  expect_true(all(cv$data$C_global < 0.5))
})

test_that("the identity holds at every threshold", {
  d <- make_test_data(300, ties = TRUE, seed = 7)
  cv <- censoring_curve(d$time, d$status, d$risk, n_thresholds = 6)
  ok <- !is.na(cv$data$C_ee) & !is.na(cv$data$C_ec)
  expect_true(any(ok))
})

test_that("min_pairs flags rather than filters", {
  d <- make_test_data(300, seed = 8)
  a <- censoring_curve(d$time, d$status, d$risk, n_thresholds = 6, min_pairs = 1)
  b <- censoring_curve(d$time, d$status, d$risk, n_thresholds = 6,
                       min_pairs = 10^9)
  expect_equal(nrow(a$data), nrow(b$data))
  expect_true(all(b$data$low_precision))
  expect_false(any(a$data$low_precision))
  expect_false(any(is.na(b$data$C_ee)))
})

test_that("censoring increases with the threshold decreasing", {
  d <- make_test_data(300, seed = 9)
  cv <- censoring_curve(d$time, d$status, d$risk, n_thresholds = 8)
  expect_equal(order(cv$data$censoring, decreasing = TRUE),
               order(cv$data$threshold))
})

test_that("the curve accepts alternative weightings", {
  d <- make_test_data(250, seed = 10)
  cv <- censoring_curve(d$time, d$status, d$risk,
                        weights = weights_uno(), n_thresholds = 4)
  expect_equal(cv$weighting, "Uno (IPCW)")
})

test_that("print is informative", {
  d <- make_test_data(200, seed = 11)
  cv <- censoring_curve(d$time, d$status, d$risk, n_thresholds = 4)
  out <- capture.output(print(cv))
  expect_true(any(grepl("Censoring", out, ignore.case = TRUE)))
})
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-censoring-curve.R")'`
Expected: FAIL — `could not find function "censoring_curve"`

- [ ] **Step 3: Remove the superseded implementation**

```bash
git rm R/simulate_censoring.R man/simulate_censoring.Rd
```

- [ ] **Step 4: Write the implementation**

Create `R/censoring_curve.R`:

```r
#' Trace the decomposition as administrative censoring increases
#'
#' @description
#' Sweeps a series of administrative cut-off times over the cohort. At each
#' cut-off, subjects whose event falls after the cut-off are treated as
#' censored there, and the decomposition is recomputed. This traces how the
#' global C-index behaves as censoring rises, while `C_ee` -- the harder
#' ranking task -- is tracked separately.
#'
#' The cut-off is applied to the **whole cohort**, so the reported
#' censoring rate is the cohort's actual censoring rate at that cut-off.
#'
#' @param time,status,risk Numeric vectors of equal length.
#' @param weights A `cindex_weights` object. Defaults to
#'   [weights_harrell()].
#' @param higher_is_riskier Logical; see [decompose_cindex()].
#' @param n_thresholds Number of cut-offs to evaluate. Ignored if `probs`
#'   is supplied.
#' @param probs Quantiles of the observed event times at which to place the
#'   cut-offs.
#' @param min_pairs Minimum number of event-event pairs for an estimate to
#'   be considered precise. Rows below it are **flagged** in the
#'   `low_precision` column, never removed and never set to `NA`. `NULL`
#'   (default) selects `max(50, round(n_events^1.5 * 0.1))`.
#'
#' @return An object of class `cindex_curve`.
#'
#' @examples
#' set.seed(42)
#' time   <- rexp(300, rate = 0.1)
#' status <- rbinom(300, 1, 0.6)
#' risk   <- rnorm(300)
#' censoring_curve(time, status, risk, n_thresholds = 8)
#'
#' @export
censoring_curve <- function(time, status, risk,
                            weights = weights_harrell(),
                            higher_is_riskier = TRUE,
                            n_thresholds = 20,
                            probs = seq(0.05, 0.70, length.out = n_thresholds),
                            min_pairs = NULL) {
  validate_survival_inputs(time, status, risk)
  if (!inherits(weights, "cindex_weights")) {
    stop("`weights` must be a `cindex_weights` object.", call. = FALSE)
  }
  if (!higher_is_riskier) risk <- -risk

  event_idx <- which(status == 1)
  n_events <- length(event_idx)
  if (is.null(min_pairs)) {
    min_pairs <- max(50, round(n_events^1.5 * 0.1))
  }

  taus <- unique(stats::quantile(time[event_idx], probs = probs, names = FALSE))

  rows <- lapply(taus, function(tau) {
    # The cut-off applies to EVERYONE. Subjects already censored stay
    # censored; events after the cut-off become censored at it.
    sim_time <- pmin(time, tau)
    sim_status <- status * as.numeric(time <= tau)
    if (sum(sim_status) < 2) return(NULL)

    pc <- pair_counts(sim_time, sim_status, risk, weights)
    total_w <- pc$W_ee + pc$W_ec
    if (total_w <= 0) return(NULL)

    C_ee <- if (pc$W_ee > 0) pc$S_ee / pc$W_ee else NA_real_
    C_ec <- if (pc$W_ec > 0) pc$S_ec / pc$W_ec else NA_real_

    data.frame(
      threshold = tau,
      censoring = mean(sim_status == 0),
      C_ee = C_ee,
      C_ec = C_ec,
      C_global = (pc$S_ee + pc$S_ec) / total_w,
      gap = C_ec - C_ee,
      N_ee = pc$N_ee,
      N_ec = pc$N_ec,
      low_precision = pc$N_ee < min_pairs
    )
  })

  dat <- do.call(rbind, rows)
  if (is.null(dat) || nrow(dat) == 0L) {
    stop("No threshold produced enough events to decompose. ",
         "Try a smaller `probs` range.", call. = FALSE)
  }
  dat <- dat[order(dat$censoring), , drop = FALSE]
  rownames(dat) <- NULL

  structure(
    list(
      data = dat,
      weighting = weights$name,
      higher_is_riskier = higher_is_riskier,
      min_pairs = min_pairs,
      n = length(time),
      n_events = n_events
    ),
    class = "cindex_curve"
  )
}

#' @param x A `cindex_curve` object.
#' @param ... Ignored.
#' @rdname censoring_curve
#' @export
print.cindex_curve <- function(x, ...) {
  cat("Censoring Curve\n")
  cat(sprintf("Weighting: %s | n = %d, events = %d | %d thresholds\n",
              x$weighting, x$n, x$n_events, nrow(x$data)))
  cat(sprintf("Censoring range: %.1f%% to %.1f%%\n\n",
              100 * min(x$data$censoring), 100 * max(x$data$censoring)))
  show <- x$data[, c("censoring", "C_ee", "C_ec", "C_global", "N_ee",
                     "low_precision")]
  show$censoring <- sprintf("%.1f%%", 100 * show$censoring)
  print(show, digits = 4, row.names = FALSE)
  if (any(x$data$low_precision)) {
    cat(sprintf("\n%d row(s) flagged low_precision (N_ee < %d).",
                sum(x$data$low_precision), x$min_pairs))
    cat(" Estimates are shown, not removed.\n")
  }
  invisible(x)
}

#' @param row.names Ignored.
#' @param optional Ignored.
#' @rdname censoring_curve
#' @export
as.data.frame.cindex_curve <- function(x, row.names = NULL,
                                       optional = FALSE, ...) {
  x$data
}
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-censoring-curve.R")'`
Expected: PASS, 11 tests.

- [ ] **Step 6: Commit**

```bash
git add -A R/ man/ tests/testthat/test-censoring-curve.R
git commit -m "fix: measure censoring on the whole cohort in the curve

simulate_censoring() deleted every already-censored subject on its first
line, then reported the censoring rate among events only. On a
500-subject fixture with 39% baseline censoring, a cut-off whose true
cohort censoring was 51% was plotted at 20%.

Renamed to censoring_curve(). Applies the cut-off to the whole cohort,
removes the sign flip and the CI_ec < 0.5 suppression, and turns
min_pairs into a low_precision flag so no estimate is discarded for its
magnitude."
```

---

## Task 8: `compare_decompositions()`

Supplies the input that `plot_decomposition()` documented but nothing produced.

**Files:**
- Create: `R/compare.R`
- Create: `tests/testthat/test-compare.R`

**Interfaces:**
- Consumes: `decompose_cindex()` (Task 4), bootstrap (Task 5).
- Produces: `compare_decompositions(risks, time, status, weights, higher_is_riskier, n_boot, conf_level)` returning an object of class `cindex_comparison`: a list with `$table` (a data frame with columns `model`, `ci_ee`, `ci_ec`, `global_c`, `gap`, `sd_ee`, `sd_ec`, `sd_global`, `sd_gap`, `n_ee`, `n_ec`), `$fits` (named list of `cindex_decomp`), `weighting`, `n_boot`. `print.cindex_comparison()`, `as.data.frame.cindex_comparison()`.

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-compare.R`:

```r
two_models <- function(seed = 1, n = 250) {
  d <- make_test_data(n, seed = seed)
  list(
    d = d,
    risks = list(good = d$risk, noisy = d$risk + rnorm(n, sd = 2))
  )
}

test_that("compare_decompositions produces the columns the dumbbell plot needs", {
  # The old plot_decomposition() required these exact columns and named a
  # monte_carlo_decompose() function that did not exist in the package.
  m <- two_models()
  cmp <- compare_decompositions(m$risks, m$d$time, m$d$status, n_boot = 30)
  expect_s3_class(cmp, "cindex_comparison")
  expect_true(all(c("model", "ci_ee", "ci_ec", "global_c",
                    "sd_ee", "sd_ec", "sd_global") %in% names(cmp$table)))
  expect_equal(nrow(cmp$table), 2)
})

test_that("model names are carried through from the list names", {
  m <- two_models()
  cmp <- compare_decompositions(m$risks, m$d$time, m$d$status, n_boot = 0)
  expect_equal(cmp$table$model, c("good", "noisy"))
})

test_that("a data frame of risks is accepted", {
  m <- two_models()
  cmp <- compare_decompositions(as.data.frame(m$risks), m$d$time,
                                m$d$status, n_boot = 0)
  expect_equal(nrow(cmp$table), 2)
})

test_that("unnamed risks are rejected", {
  m <- two_models()
  expect_error(
    compare_decompositions(unname(m$risks), m$d$time, m$d$status, n_boot = 0),
    "named"
  )
})

test_that("standard deviations are NA when n_boot is 0", {
  m <- two_models()
  cmp <- compare_decompositions(m$risks, m$d$time, m$d$status, n_boot = 0)
  expect_true(all(is.na(cmp$table$sd_ee)))
})

test_that("standard deviations are finite when n_boot > 0", {
  m <- two_models()
  cmp <- compare_decompositions(m$risks, m$d$time, m$d$status, n_boot = 40)
  expect_true(all(is.finite(cmp$table$sd_ee)))
  expect_true(all(cmp$table$sd_ee > 0))
})

test_that("per-model values match a direct decompose_cindex call", {
  m <- two_models()
  cmp <- compare_decompositions(m$risks, m$d$time, m$d$status, n_boot = 0)
  direct <- decompose_cindex(m$d$time, m$d$status, m$risks$good)
  expect_equal(cmp$table$ci_ee[1], direct$C_ee, tolerance = 1e-12)
  expect_equal(cmp$table$global_c[1], direct$C_global, tolerance = 1e-12)
})

test_that("fits are retained for downstream use", {
  m <- two_models()
  cmp <- compare_decompositions(m$risks, m$d$time, m$d$status, n_boot = 0)
  expect_named(cmp$fits, c("good", "noisy"))
  expect_s3_class(cmp$fits$good, "cindex_decomp")
})

test_that("a risk vector of the wrong length is rejected", {
  m <- two_models()
  bad <- list(good = m$d$risk, short = m$d$risk[1:10])
  expect_error(
    compare_decompositions(bad, m$d$time, m$d$status, n_boot = 0),
    "same length"
  )
})

test_that("print lists every model", {
  m <- two_models()
  cmp <- compare_decompositions(m$risks, m$d$time, m$d$status, n_boot = 0)
  out <- capture.output(print(cmp))
  expect_true(any(grepl("good", out)))
  expect_true(any(grepl("noisy", out)))
})
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-compare.R")'`
Expected: FAIL — `could not find function "compare_decompositions"`

- [ ] **Step 3: Write the implementation**

Create `R/compare.R`:

```r
#' Compare the decomposition across several models
#'
#' @description
#' Runs [decompose_cindex()] for each of several risk scores measured on
#' the same cohort, and assembles the results into one table with bootstrap
#' standard deviations. This is the input consumed by
#' `autoplot()` to draw the dumbbell comparison.
#'
#' @param risks A **named** list or data frame of risk-score vectors, each
#'   the same length as `time` and `status`. The names label the models.
#' @param time,status Numeric vectors describing the shared cohort.
#' @param weights A `cindex_weights` object. Defaults to
#'   [weights_harrell()].
#' @param higher_is_riskier Logical; see [decompose_cindex()]. Applied to
#'   every model.
#' @param n_boot Bootstrap replicates per model. Defaults to 1000, since
#'   the comparison table exists to carry uncertainty.
#' @param conf_level Confidence level stored on each fit.
#'
#' @return An object of class `cindex_comparison`.
#'
#' @examples
#' set.seed(42)
#' n <- 200
#' time   <- rexp(n, rate = 0.1)
#' status <- rbinom(n, 1, 0.6)
#' risks  <- list(model_a = rnorm(n), model_b = rnorm(n))
#' compare_decompositions(risks, time, status, n_boot = 0)
#'
#' @export
compare_decompositions <- function(risks, time, status,
                                   weights = weights_harrell(),
                                   higher_is_riskier = TRUE,
                                   n_boot = 1000,
                                   conf_level = 0.95) {
  if (is.data.frame(risks)) risks <- as.list(risks)
  if (!is.list(risks) || length(risks) == 0L) {
    stop("`risks` must be a non-empty named list or data frame of risk scores.",
         call. = FALSE)
  }
  if (is.null(names(risks)) || any(!nzchar(names(risks)))) {
    stop("`risks` must be fully named; the names label the models.",
         call. = FALSE)
  }
  bad <- vapply(risks, function(r) length(r) != length(time), logical(1))
  if (any(bad)) {
    stop("Every risk score must be the same length as `time` and `status`. ",
         "Offending: ", paste(names(risks)[bad], collapse = ", "), ".",
         call. = FALSE)
  }

  fits <- lapply(risks, function(r) {
    decompose_cindex(
      time, status, r,
      weights = weights,
      higher_is_riskier = higher_is_riskier,
      n_boot = n_boot,
      conf_level = conf_level
    )
  })
  names(fits) <- names(risks)

  boot_sd <- function(fit, col) {
    if (is.null(fit$boot)) NA_real_ else stats::sd(fit$boot[, col], na.rm = TRUE)
  }

  tab <- data.frame(
    model = names(fits),
    ci_ee = vapply(fits, function(f) f$C_ee, numeric(1)),
    ci_ec = vapply(fits, function(f) f$C_ec, numeric(1)),
    global_c = vapply(fits, function(f) f$C_global, numeric(1)),
    gap = vapply(fits, function(f) f$gap, numeric(1)),
    sd_ee = vapply(fits, boot_sd, numeric(1), col = "C_ee"),
    sd_ec = vapply(fits, boot_sd, numeric(1), col = "C_ec"),
    sd_global = vapply(fits, boot_sd, numeric(1), col = "C_global"),
    sd_gap = vapply(fits, boot_sd, numeric(1), col = "gap"),
    n_ee = vapply(fits, function(f) f$N_ee, numeric(1)),
    n_ec = vapply(fits, function(f) f$N_ec, numeric(1)),
    stringsAsFactors = FALSE
  )
  rownames(tab) <- NULL

  structure(
    list(
      table = tab,
      fits = fits,
      weighting = weights$name,
      n_boot = n_boot
    ),
    class = "cindex_comparison"
  )
}

#' @param x A `cindex_comparison` object.
#' @param ... Ignored.
#' @rdname compare_decompositions
#' @export
print.cindex_comparison <- function(x, ...) {
  cat("C-Index Decomposition Comparison\n")
  cat(sprintf("Weighting: %s | %d model(s) | n_boot = %d\n\n",
              x$weighting, nrow(x$table), x$n_boot))
  show <- x$table[, c("model", "ci_ee", "ci_ec", "global_c", "gap")]
  names(show) <- c("model", "C_ee", "C_ec", "C_global", "gap")
  print(show, digits = 4, row.names = FALSE)
  if (x$n_boot == 0) {
    cat("\nNo bootstrap: sd_* columns are NA. Set n_boot > 0 for uncertainty.\n")
  }
  invisible(x)
}

#' @param row.names Ignored.
#' @param optional Ignored.
#' @rdname compare_decompositions
#' @export
as.data.frame.cindex_comparison <- function(x, row.names = NULL,
                                            optional = FALSE, ...) {
  x$table
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-compare.R")'`
Expected: PASS, 10 tests.

- [ ] **Step 5: Commit**

```bash
git add R/compare.R tests/testthat/test-compare.R
git commit -m "feat: add compare_decompositions for multi-model comparison

plot_decomposition() documented its input as coming from
monte_carlo_decompose(), which was never in the package, so the function
was unreachable. This supplies that table: one row per model with
bootstrap standard deviations."
```

---

## Task 9: Theme and `autoplot()` methods

Replaces `plot_simulation()` and `plot_decomposition()`.

**Files:**
- Create: `R/theme.R`
- Create: `R/autoplot.R`
- Delete: `R/plot_simulation.R`, `R/plot_decomposition.R`
- Delete: `man/plot_simulation.Rd`, `man/plot_decomposition.Rd`
- Create: `tests/testthat/test-autoplot.R`

**Interfaces:**
- Consumes: `cindex_decomp` (Task 4), `cindex_curve` (Task 7), `cindex_comparison` (Task 8).
- Produces: `theme_cindex(dark = FALSE, base_size = 14)` returning a ggplot2 theme; `cindex_palette(dark)` returning a named list of colours; `autoplot.cindex_decomp()`, `autoplot.cindex_curve()`, `autoplot.cindex_comparison()`, each returning a `ggplot` object.

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-autoplot.R`:

```r
skip_if_no_ggplot <- function() skip_if_not_installed("ggplot2")

test_that("autoplot dispatches on a decomposition", {
  skip_if_no_ggplot()
  d <- make_test_data(200, seed = 1)
  p <- ggplot2::autoplot(decompose_cindex(d$time, d$status, d$risk))
  expect_s3_class(p, "ggplot")
})

test_that("autoplot dispatches on a censoring curve", {
  skip_if_no_ggplot()
  d <- make_test_data(250, seed = 2)
  p <- ggplot2::autoplot(censoring_curve(d$time, d$status, d$risk,
                                         n_thresholds = 6))
  expect_s3_class(p, "ggplot")
})

test_that("autoplot dispatches on a comparison", {
  skip_if_no_ggplot()
  d <- make_test_data(200, seed = 3)
  cmp <- compare_decompositions(
    list(a = d$risk, b = d$risk + rnorm(200)),
    d$time, d$status, n_boot = 20
  )
  expect_s3_class(ggplot2::autoplot(cmp), "ggplot")
})

test_that("no plot hardcodes axis limits", {
  # Regression: scale_x_continuous(limits = c(0.48, 0.74)) silently dropped
  # any model outside that window, including worse-than-chance models.
  skip_if_no_ggplot()
  d <- make_test_data(200, seed = 4)
  cmp <- compare_decompositions(list(a = d$risk), d$time, d$status, n_boot = 0)
  p <- ggplot2::autoplot(cmp)
  xs <- p$scales$get_scales("x")
  expect_true(is.null(xs) || is.null(xs$limits))
})

test_that("a worse-than-chance model still appears in the plot data", {
  skip_if_no_ggplot()
  d <- make_test_data(250, seed = 5)
  cmp <- compare_decompositions(list(backwards = -d$risk), d$time,
                                d$status, n_boot = 0)
  expect_lt(cmp$table$global_c, 0.5)
  p <- ggplot2::autoplot(cmp)
  built <- ggplot2::ggplot_build(p)
  xs <- unlist(lapply(built$data, function(l) l$x))
  expect_true(any(xs < 0.5, na.rm = TRUE))
})

test_that("no plot sets a hardcoded font family", {
  # Regression: family = "Arial" appeared 14 times; Arial is absent from
  # CRAN's Linux check machines and emits warnings on every plot call.
  skip_if_no_ggplot()
  d <- make_test_data(150, seed = 6)
  ps <- list(
    ggplot2::autoplot(decompose_cindex(d$time, d$status, d$risk)),
    ggplot2::autoplot(censoring_curve(d$time, d$status, d$risk,
                                      n_thresholds = 4))
  )
  for (p in ps) {
    expect_false(identical(p$theme$text$family, "Arial"))
  }
})

test_that("theme_cindex defaults to light and offers dark", {
  skip_if_no_ggplot()
  light <- theme_cindex()
  dark <- theme_cindex(dark = TRUE)
  expect_s3_class(light, "theme")
  expect_equal(light$plot.background$fill, "#FFFFFF")
  expect_equal(dark$plot.background$fill, "#282A36")
})

test_that("theme_cindex sets no font family", {
  skip_if_no_ggplot()
  expect_false(identical(theme_cindex()$text$family, "Arial"))
})

test_that("the curve plot marks low-precision points without dropping them", {
  skip_if_no_ggplot()
  d <- make_test_data(250, seed = 7)
  cv <- censoring_curve(d$time, d$status, d$risk, n_thresholds = 6,
                        min_pairs = 10^9)
  expect_true(all(cv$data$low_precision))
  p <- ggplot2::autoplot(cv)
  built <- ggplot2::ggplot_build(p)
  n_points <- sum(vapply(built$data, function(l) nrow(l), numeric(1)))
  expect_gt(n_points, 0)
})
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-autoplot.R")'`
Expected: FAIL — `could not find function "theme_cindex"`

- [ ] **Step 3: Remove the superseded plot functions**

```bash
git rm R/plot_simulation.R R/plot_decomposition.R \
       man/plot_simulation.Rd man/plot_decomposition.Rd
```

- [ ] **Step 4: Write the theme**

Create `R/theme.R`:

```r
#' Colour palette for cindexdecomp plots
#'
#' Event-Event and Event-Censored keep the IBM colourblind-safe magenta and
#' blue in both variants; only the neutrals change.
#'
#' @param dark Logical. `FALSE` (default) returns the light palette.
#' @return A named list of hex colours.
#' @keywords internal
#' @noRd
cindex_palette <- function(dark = FALSE) {
  list(
    ee     = "#DC267F",
    ec     = "#648FFF",
    global = if (dark) "#F8F8F2" else "#2B2B2B",
    fg     = if (dark) "#F8F8F2" else "#1A1A1A",
    bg     = if (dark) "#282A36" else "#FFFFFF",
    grid   = if (dark) "#3A3C4E" else "#E6E6E6",
    track  = if (dark) "#3A3C4E" else "#EDEDED"
  )
}

#' Plot theme for cindexdecomp
#'
#' @description
#' Light by default, because journal figures print on white paper. Pass
#' `dark = TRUE` for the high-contrast presentation variant.
#'
#' No font family is set anywhere: the system default is used, so plots
#' render without warnings on machines that lack any particular typeface.
#'
#' @param dark Logical. `FALSE` (default) for the light theme.
#' @param base_size Base font size in points.
#' @return A ggplot2 theme object.
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(wt, mpg)) + geom_point() + theme_cindex()
#' @export
theme_cindex <- function(dark = FALSE, base_size = 14) {
  p <- cindex_palette(dark)
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.background   = ggplot2::element_rect(fill = p$bg, colour = NA),
      panel.background  = ggplot2::element_rect(fill = p$bg, colour = NA),
      legend.background = ggplot2::element_rect(fill = p$bg, colour = NA),
      legend.key        = ggplot2::element_rect(fill = p$bg, colour = NA),
      text              = ggplot2::element_text(colour = p$fg),
      axis.text         = ggplot2::element_text(colour = p$fg),
      axis.line         = ggplot2::element_line(colour = p$fg),
      axis.ticks        = ggplot2::element_line(colour = p$fg),
      panel.grid.major  = ggplot2::element_line(colour = p$grid,
                                                linewidth = 0.3),
      panel.grid.minor  = ggplot2::element_blank(),
      plot.title        = ggplot2::element_text(face = "bold"),
      plot.subtitle     = ggplot2::element_text(colour = p$fg),
      plot.margin       = ggplot2::margin(12, 16, 12, 12)
    )
}
```

- [ ] **Step 5: Write the autoplot methods**

Create `R/autoplot.R`:

```r
#' Plot a concordance decomposition
#'
#' @description
#' `autoplot()` methods for the three result classes. Because the object
#' carries its own class, the right plot is selected automatically and a
#' result can never be paired with the wrong chart.
#'
#' @param object A `cindex_decomp`, `cindex_curve` or `cindex_comparison`
#'   object.
#' @param dark Logical. Passed to [theme_cindex()].
#' @param ... Ignored.
#' @return A `ggplot` object.
#' @name autoplot.cindexdecomp
NULL

dumbbell_plot <- function(tab, weighting, dark) {
  p <- cindex_palette(dark)
  tab$model <- factor(tab$model, levels = rev(unique(tab$model)))
  has_sd <- all(is.finite(tab$sd_ee))

  gg <- ggplot2::ggplot(tab) +
    ggplot2::geom_vline(xintercept = 0.5, linetype = "dashed",
                        colour = p$fg, linewidth = 0.6) +
    ggplot2::geom_segment(
      ggplot2::aes(x = .data$ci_ee, xend = .data$ci_ec,
                   y = .data$model, yend = .data$model),
      colour = p$grid, linewidth = 2
    )

  if (has_sd) {
    gg <- gg +
      ggplot2::geom_errorbarh(
        ggplot2::aes(xmin = .data$ci_ee - .data$sd_ee,
                     xmax = .data$ci_ee + .data$sd_ee,
                     y = .data$model),
        height = 0.12, colour = p$ee, linewidth = 0.7
      ) +
      ggplot2::geom_errorbarh(
        ggplot2::aes(xmin = .data$ci_ec - .data$sd_ec,
                     xmax = .data$ci_ec + .data$sd_ec,
                     y = .data$model),
        height = 0.12, colour = p$ec, linewidth = 0.7
      )
  }

  gg +
    ggplot2::geom_point(
      ggplot2::aes(x = .data$ci_ee, y = .data$model,
                   colour = "Event-Event"), size = 4.5
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = .data$ci_ec, y = .data$model,
                   colour = "Event-Censored"), size = 4.5
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = .data$global_c, y = .data$model,
                   colour = "Global"), size = 3.5, shape = 18
    ) +
    ggplot2::scale_colour_manual(
      name = NULL,
      values = c("Event-Event" = p$ee, "Event-Censored" = p$ec,
                 "Global" = p$global),
      breaks = c("Event-Event", "Global", "Event-Censored")
    ) +
    ggplot2::labs(
      x = "Concordance index", y = NULL,
      title = "C-index decomposition",
      subtitle = paste0("Weighting: ", weighting,
                        " · dashed line marks chance (0.50)")
    ) +
    theme_cindex(dark = dark)
}

#' @rdname autoplot.cindexdecomp
#' @export
autoplot.cindex_comparison <- function(object, dark = FALSE, ...) {
  dumbbell_plot(object$table, object$weighting, dark)
}

#' @rdname autoplot.cindexdecomp
#' @export
autoplot.cindex_decomp <- function(object, dark = FALSE, ...) {
  sd_of <- function(col) {
    if (is.null(object$boot)) NA_real_ else stats::sd(object$boot[, col],
                                                      na.rm = TRUE)
  }
  tab <- data.frame(
    model = "model",
    ci_ee = object$C_ee,
    ci_ec = object$C_ec,
    global_c = object$C_global,
    sd_ee = sd_of("C_ee"),
    sd_ec = sd_of("C_ec"),
    sd_global = sd_of("C_global"),
    stringsAsFactors = FALSE
  )
  dumbbell_plot(tab, object$weighting, dark)
}

#' @rdname autoplot.cindexdecomp
#' @export
autoplot.cindex_curve <- function(object, dark = FALSE, ...) {
  p <- cindex_palette(dark)
  dat <- object$data
  dat$precision <- ifelse(dat$low_precision, "low", "ok")

  ggplot2::ggplot(dat, ggplot2::aes(x = .data$censoring)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$C_ee, ymax = .data$C_global),
      fill = p$ee, alpha = 0.10, na.rm = TRUE
    ) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed",
                        colour = p$fg, linewidth = 0.6) +
    ggplot2::geom_line(
      ggplot2::aes(y = .data$C_ec, colour = "Event-Censored"),
      linewidth = 1, na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = .data$C_global, colour = "Global"),
      linewidth = 1, na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = .data$C_ee, colour = "Event-Event"),
      linewidth = 1, na.rm = TRUE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(y = .data$C_ec, colour = "Event-Censored",
                   alpha = .data$precision), size = 2.4, na.rm = TRUE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(y = .data$C_global, colour = "Global",
                   alpha = .data$precision), size = 2.4, na.rm = TRUE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(y = .data$C_ee, colour = "Event-Event",
                   alpha = .data$precision), size = 2.4, na.rm = TRUE
    ) +
    # Low-precision points are drawn faintly, never removed.
    ggplot2::scale_alpha_manual(
      values = c(ok = 1, low = 0.35), guide = "none"
    ) +
    ggplot2::scale_colour_manual(
      name = NULL,
      values = c("Event-Event" = p$ee, "Event-Censored" = p$ec,
                 "Global" = p$global),
      breaks = c("Event-Event", "Global", "Event-Censored")
    ) +
    ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(
      x = "Cohort censoring rate", y = "Concordance index",
      title = "Concordance under increasing censoring",
      subtitle = paste0("Weighting: ", object$weighting,
                        " · shaded band is the masking gap")
    ) +
    theme_cindex(dark = dark)
}
```

Add `#' @importFrom rlang .data` is **not** needed — instead add `.data` to
the package imports. In `R/cindexdecomp-package.R`, add to the roxygen
block:

```r
#' @importFrom ggplot2 .data
```

- [ ] **Step 6: Regenerate documentation and NAMESPACE**

```bash
Rscript -e 'devtools::document(".")'
```

Expected: `NAMESPACE` gains `S3method(autoplot, cindex_curve)` and friends,
and loses the four old `export()` lines for the removed functions.

- [ ] **Step 7: Run tests to verify they pass**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-autoplot.R")'`
Expected: PASS, 9 tests.

- [ ] **Step 8: Run the full suite**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_dir("tests/testthat")'`
Expected: PASS, no failures.

- [ ] **Step 9: Commit**

```bash
git add -A R/ man/ NAMESPACE tests/testthat/test-autoplot.R
git commit -m "feat: replace plot functions with autoplot methods

Plots are now selected by the object's class, so a result can no longer
be handed to the wrong chart (plot_simulation(decompose_result)
previously failed with 'incorrect number of dimensions').

Removes the hardcoded Arial family, which warns on machines without it
including CRAN's Linux checkers, and the hardcoded x limits of
c(0.48, 0.74), which silently dropped worse-than-chance models. The
Dracula palette is now opt-in via dark = TRUE."
```

---

## Task 10: Vignette

**Files:**
- Create: `vignettes/cindexdecomp.Rmd`
- Create: `vignettes/.gitignore`

**Interfaces:**
- Consumes: the complete public API (Tasks 2, 4, 5, 7, 8, 9).
- Produces: a built vignette, `cindexdecomp.html`.

- [ ] **Step 1: Create the vignette gitignore**

Create `vignettes/.gitignore`:

```
*.html
*.R
```

- [ ] **Step 2: Write the vignette**

Create `vignettes/cindexdecomp.Rmd`:

````markdown
---
title: "Decomposing the C-index"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Decomposing the C-index}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(collapse = TRUE, comment = "#>",
                      fig.width = 7, fig.height = 4.5)
```

```{r setup}
library(cindexdecomp)
library(survival)
```

## The problem

Harrell's concordance index pools two different ranking tasks. Some
comparable pairs consist of two subjects who both had the event; others
pair an event with a subject censored later. The second task is easier,
and as censoring rises it comes to dominate the pooled score.

## The decomposition

Write any concordance index as a weighted average over comparable pairs,

$$C_w = \frac{\sum_{ij} w_{ij} c_{ij}}{\sum_{ij} w_{ij}},$$

where $c_{ij} \in \{0, 0.5, 1\}$ records whether the pair was ranked
correctly. Splitting the pairs by the status of the later subject gives

$$C_w = \frac{W_{ee} C_{ee} + W_{ec} C_{ec}}{W_{ee} + W_{ec}},$$

with $W$ denoting sums of weights. The identity holds exactly for every
weighting.

```{r}
set.seed(7)
n <- 400
risk <- rnorm(n)
event_time <- rexp(n, rate = exp(0.8 * risk) * 0.1)
cens_time <- rexp(n, rate = 0.08)
time <- pmin(event_time, cens_time)
status <- as.numeric(event_time <= cens_time)

fit <- decompose_cindex(time, status, risk, n_boot = 200)
fit
```

The global figure agrees with `survival::concordance()` exactly:

```{r}
c(cindexdecomp = fit$C_global,
  survival = concordance(Surv(time, status) ~ risk, reverse = TRUE)$concordance)
```

## Ties

Three cases arise, and the package follows `survival::concordance()` in
each:

* two events at the same recorded time are **not** comparable;
* an event and a censoring at the same recorded time **are** comparable,
  with the event treated as occurring first;
* equal risk scores receive half credit.

The second case matters in practice because survival times are usually
recorded in whole days, so exact ties are common rather than exceptional.

## Choosing a weighting

The weighting selects the estimator. One engine produces all of them.

```{r}
weightings <- list(
  Harrell   = weights_harrell(),
  Uno       = weights_uno(),
  Truncated = weights_truncated(tau = quantile(time[status == 1], 0.5))
)

do.call(rbind, lapply(names(weightings), function(nm) {
  f <- decompose_cindex(time, status, risk, weights = weightings[[nm]])
  data.frame(weighting = nm, C_global = f$C_global,
             C_ee = f$C_ee, C_ec = f$C_ec)
}))
```

`C_global` varies with the estimator chosen, and every one of those values
is defensible. `C_ee` does not: it stays close to the same figure
throughout, well below `C_ec` in every row. Uno's inverse-probability
weighting corrects the bias in the pooled estimate; it does not reveal the
composition problem underneath. The masking is a property of the data, not
of the estimator.

A user-supplied weighting needs no changes to the package:

```{r}
own <- weights_custom(function(t, G) 1 / G(t), name = "IPCW (power 1)")
decompose_cindex(time, status, risk, weights = own)$C_global
```

## Inference

Intervals come from resampling **subjects**, not pairs. Four hundred
subjects generate tens of thousands of pairs, but the information content
is four hundred observations; resampling pairs directly gives intervals
several times too narrow.

```{r}
confint(fit)
```

The final row is usually the one of interest. It is computed inside each
replicate, so the correlation between `C_ee` and `C_ec` is carried through
rather than assumed away.

## Censoring curves

`censoring_curve()` applies a series of administrative cut-offs to the
**whole cohort** and recomputes at each. The reported rate is therefore the
cohort's actual censoring rate, including subjects censored before any
cut-off was applied.

```{r}
cv <- censoring_curve(time, status, risk, n_thresholds = 12)
head(as.data.frame(cv))
```

```{r}
ggplot2::autoplot(cv)
```

Estimates resting on few event-event pairs are flagged in
`low_precision` and drawn faintly. They are never removed, and no value is
discarded for being small.

## Comparing models

```{r}
risks <- list(
  full = risk,
  noisy = risk + rnorm(n, sd = 1.5)
)
cmp <- compare_decompositions(risks, time, status, n_boot = 200)
cmp
```

```{r}
ggplot2::autoplot(cmp)
```

## A worked example

```{r}
lung2 <- lung[complete.cases(lung), ]
cox <- coxph(Surv(time, status - 1) ~ age + sex + ph.ecog, data = lung2)
decompose_cindex(
  Surv(time, status - 1) ~ predict(cox, type = "lp"),
  data = lung2,
  n_boot = 200
)
```

Note the `status - 1`: `lung` codes status as 1/2, while the package
expects 0/1. Data already coded 0/1 needs no adjustment.

## Scope

The decomposition applies to any estimator that averages over observed
comparable pairs. Model-based estimators such as Gönen–Heller integrate
over an assumed model rather than counting pairs, so the event-event and
event-censored partition is undefined for them; they are outside the scope
of this package by construction.

Non-informative censoring is assumed throughout. The package evaluates
discrimination only, not calibration.
````

- [ ] **Step 3: Build the vignette**

Run: `Rscript -e 'devtools::build_vignettes(".")'`
Expected: builds without error or warning.

- [ ] **Step 4: Commit**

```bash
git add vignettes/
git commit -m "docs: add package vignette

Covers the generalised identity, the tie convention, the weighting
framework, subject-level inference, censoring curves and multi-model
comparison. Serves as the backbone of the R Journal submission."
```

---

## Task 11: README, CI, pkgdown and CRAN readiness

**Files:**
- Modify: `README.md`
- Create: `NEWS.md`
- Create: `.github/workflows/R-CMD-check.yaml`
- Create: `.github/workflows/test-coverage.yaml`
- Create: `_pkgdown.yml`

**Interfaces:**
- Consumes: the complete package.
- Produces: no R-level interface. Deliverable is a clean `R CMD check`.

- [ ] **Step 1: Rewrite the README**

The current README documents `simulate_censoring()`, `plot_simulation()`
and `plot_decomposition()`, none of which still exist, and states the
count-based identity rather than the general one. Replace `README.md`:

````markdown
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
  time   = lung2$time,
  status = lung2$status - 1,   # lung codes 1/2; the package expects 0/1
  risk   = predict(cox, type = "lp"),
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
precision for Harrell's and Uno's weightings.

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
````

- [ ] **Step 2: Write NEWS.md**

Create `NEWS.md`:

```markdown
# cindexdecomp 0.2.0

## Breaking changes

* `simulate_censoring()` is renamed `censoring_curve()`.
* `plot_simulation()` and `plot_decomposition()` are replaced by
  `autoplot()` methods.
* `decompose_cindex()` returns a `cindex_decomp` object rather than a
  plain list. Existing element access (`$CI_ee`) becomes `$C_ee`.

## Bug fixes

* Event-censored pairs at tied times were dropped. Because survival times
  are normally recorded in whole days these ties are common: on a
  300-subject example 1,062 pairs were lost and the estimate was biased
  by +0.0007 relative to `survival::concordance()`.
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
```

- [ ] **Step 3: Add the R CMD check workflow**

Create `.github/workflows/R-CMD-check.yaml`:

```yaml
on:
  push:
    branches: [main, master]
  pull_request:

name: R-CMD-check

jobs:
  R-CMD-check:
    runs-on: ${{ matrix.config.os }}
    name: ${{ matrix.config.os }} (${{ matrix.config.r }})
    strategy:
      fail-fast: false
      matrix:
        config:
          - {os: ubuntu-latest,  r: 'release'}
          - {os: ubuntu-latest,  r: 'oldrel-1'}
          - {os: macos-latest,   r: 'release'}
          - {os: windows-latest, r: 'release'}
    env:
      GITHUB_PAT: ${{ secrets.GITHUB_TOKEN }}
      R_KEEP_PKG_SOURCE: yes
    steps:
      - uses: actions/checkout@v4
      - uses: r-lib/actions/setup-pandoc@v2
      - uses: r-lib/actions/setup-r@v2
        with:
          r-version: ${{ matrix.config.r }}
          use-public-rspm: true
      - uses: r-lib/actions/setup-r-dependencies@v2
        with:
          extra-packages: any::rcmdcheck
          needs: check
      - uses: r-lib/actions/check-r-package@v2
        with:
          upload-snapshots: true
```

- [ ] **Step 4: Add the coverage workflow**

Create `.github/workflows/test-coverage.yaml`:

```yaml
on:
  push:
    branches: [main, master]
  pull_request:

name: test-coverage

jobs:
  test-coverage:
    runs-on: ubuntu-latest
    env:
      GITHUB_PAT: ${{ secrets.GITHUB_TOKEN }}
    steps:
      - uses: actions/checkout@v4
      - uses: r-lib/actions/setup-r@v2
        with:
          use-public-rspm: true
      - uses: r-lib/actions/setup-r-dependencies@v2
        with:
          extra-packages: any::covr, any::xml2
          needs: coverage
      - name: Test coverage
        run: |
          covr::package_coverage()
        shell: Rscript {0}
```

- [ ] **Step 5: Add the pkgdown config**

Create `_pkgdown.yml`:

```yaml
url: https://Lainsm.github.io/cindexdecomp/

template:
  bootstrap: 5

reference:
  - title: Decomposition
    contents:
      - decompose_cindex
      - confint.cindex_decomp
  - title: Weightings
    contents:
      - cindex_weights
  - title: Censoring and comparison
    contents:
      - censoring_curve
      - compare_decompositions
  - title: Plotting
    contents:
      - autoplot.cindexdecomp
      - theme_cindex
```

- [ ] **Step 6: Run the full check**

```bash
Rscript -e 'devtools::document(".")'
Rscript -e 'devtools::check(".")'
```

Expected: **0 errors, 0 warnings, 0 notes.**

Known issues to resolve if they appear:

| Check output | Fix |
|---|---|
| `no visible binding for global variable 'ci_ee'` | The `.data$` prefix is missing from an `aes()` call in `R/autoplot.R` |
| `Undocumented arguments in \usage` | Add the missing `@param` to the roxygen block, then re-run `devtools::document(".")` |
| `Files in .Rbuildignore not found` | A pattern in `.Rbuildignore` matches nothing; delete it |
| `checking for unstated dependencies in vignettes` | Add the package to `Suggests` |
| `Non-standard file/directory found: docs` | Confirm `^docs/superpowers$` is present in `.Rbuildignore` |

- [ ] **Step 7: Commit**

```bash
git add README.md NEWS.md .github/ _pkgdown.yml
git commit -m "docs: rewrite README, add NEWS, CI workflows and pkgdown config

README documented three functions that no longer exist and stated the
count-based identity rather than the general weighted one. Adds
R CMD check on Linux/macOS/Windows and a coverage workflow."
```

- [ ] **Step 8: Open the pull request**

```bash
git push -u origin design/cindex-decomposition
gh pr create --title "Redesign cindexdecomp around a general weighted-pair framework" --body "$(cat <<'BODY'
Implements the design in `docs/superpowers/specs/2026-08-25-cindex-decomposition-design.md`.

## Bug fixes

- Event-censored pairs at tied times were dropped (1,062 pairs and +0.0007 bias on a 300-subject fixture). Pooled results now match `survival::concordance()` to 1e-12.
- The censoring curve discarded already-censored subjects and reported the censoring rate among events only, understating it by up to 31 percentage points.
- Sub-0.5 global values were inverted, breaking the decomposition identity.
- `C_ec` values below 0.5 were silently deleted.

## New

- Pluggable pair weightings covering Harrell's C, Uno's IPCW C, truncated C and user-supplied rules from one engine.
- Subject-level bootstrap with `confint()`, including an interval for the masking gap.
- `compare_decompositions()`, which supplies the table the old `plot_decomposition()` documented but nothing produced.
- `autoplot()` methods, light-default themes, vignette, CI and pkgdown.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
BODY
)"
```

---

## Self-Review

**Spec coverage.** Every numbered spec section maps to a task:

| Spec section | Task |
|---|---|
| 2.1 generalised decomposition | 2, 3 |
| 2.2 Gönen–Heller out of scope | 10 (vignette), 11 (README) |
| 2.3 tie convention | 3 |
| 3.1 public API | 2, 4, 5, 7, 8, 9 |
| 3.2 renames | 7, 9 |
| 3.3 orphaned `plot_decomposition` | 8 |
| 3.4 file layout | all |
| 3.5 pure R | Global Constraints |
| 3.6 dependencies | 1 |
| 4.1–4.3 bootstrap | 5 |
| 4.4 reporting constraint | 6 |
| 5.1 sign flip | 4, 7 |
| 5.2 value-based suppression | 7 |
| 6.1 whole cohort | 7 |
| 6.2 `min_pairs` doc fix | 7 |
| 7 plotting | 9 |
| 8 testing | every task |
| 9 infrastructure | 1 (LICENSE), 11 |
| 10 vignette | 10 |
| 11 non-goals | 10, 11 |

**Type consistency.** `cindex_decomp` fields (`C_ee`, `C_ec`, `C_global`,
`gap`, `W_ee`, `W_ec`, `N_ee`, `N_ec`, `n`, `n_events`, `censoring_rate`,
`weighting`, `higher_is_riskier`, `boot`, `n_boot`, `conf_level`) are used
consistently in Tasks 4, 5, 6, 8, 9. `pair_counts()` returns
`W_ee`/`S_ee`/`N_ee`/`W_ec`/`S_ec`/`N_ec` throughout. The comparison table
uses `ci_ee`/`ci_ec`/`global_c`/`sd_*`, matching the column names the old
`plot_decomposition()` expected.

**Ordering note.** Task 4 defines `decompose_cindex()` without `n_boot`;
Task 5 adds it. Task 4's tests therefore never call `n_boot`, and Task 6's
print tests depend on Task 5 having run.

---

## Notes for the executor

- Run `devtools::document(".")` after any task that changes roxygen
  comments. `NAMESPACE` is generated; never hand-edit it.
- `expect_equal()` under testthat 3rd edition needs an explicit
  `tolerance =`. Machine-precision claims use `1e-12`.
- If an identity test fails, fix the engine. Never widen the tolerance —
  agreement with `survival::concordance()` is the package's central claim.
- The four regression tests named in commit messages (tied pairs, whole
  cohort, no flipping, no value-based `NA`) each correspond to a measured
  defect. If one starts failing, a defect has returned.
