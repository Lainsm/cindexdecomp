# cindexdecomp (development version)

## Packaging and tooling

No user-facing behaviour changed in this section; it records the build and
documentation infrastructure so a future reader knows what is generated and
what is hand-written.

* The landing page is now generated: edit `README.Rmd` and re-render with
  `devtools::build_readme()`. Its examples execute against `survival::lung`,
  so the printed decompositions, intervals and both figures are real output
  rather than transcribed. The chunk options seed the bootstrap, so
  re-rendering does not churn the diff.
* Added a `pkgdown` site (`_pkgdown.yml`), deployed to GitHub Pages by a new
  `pkgdown` workflow. The site builds to `docs/`, which is now git-ignored;
  the design and plan documents that used to live there moved to
  `dev/superpowers/`.
* Added `inst/CITATION`, so `citation("cindexdecomp")` gives the intended
  entry rather than an auto-generated one.
* Added a `spelling` check (`tests/spelling.R` plus `inst/WORDLIST`) and set
  `Language: en-GB`.
* Added a `lint` workflow and an `.lintr.R` config. The config documents why
  `object_name_linter` is off (the package deliberately mirrors the
  estimator's notation -- `C_ee`, `W_ec`, `G(t)`) and why
  `indentation_linter` is off (`styler` is the formatter of record and the
  two tools disagree by two spaces on wrapped `if` conditions).
* The whole package is now formatted with `styler::style_pkg()` and lints
  clean under `lintr::lint_package()`.
* The `test-coverage` workflow computed coverage but never uploaded it; it
  now reports to Codecov. The upload is deliberately non-fatal until a
  `CODECOV_TOKEN` secret exists, so an unconfigured Codecov does not red the
  build.
* The documentation URL is now `https://lainsm.github.io/cindexdecomp`, in
  lower case. GitHub Pages serves user sites from the lower-cased account
  name and only reaches the capitalised form through a 301, which CRAN's URL
  check reports.
* CI actions pinned to current versions: `actions/checkout@v6`,
  `codecov/codecov-action` v7 and `github-pages-deploy-action` v4.8.0, both
  pinned by commit SHA as r-lib now does.
* `tests/testthat/test-censoring-curve.R` is renamed
  `test-censoring_curve.R` to mirror `R/censoring_curve.R`, the theme tests
  moved from `test-autoplot.R` into their own `test-theme.R`, and the shared
  `skip_if_no_ggplot()` helper moved into `helper-ggplot.R` so any test file
  can reach it.

# cindexdecomp 0.2.0

## Breaking changes

* The first argument of `decompose_cindex()` is now `x` rather than
  `time`, as required for S3 dispatch. Calls written as
  `decompose_cindex(time = t, status = s, risk = r)` fail with a
  self-explanatory error pointing at the fix (call positionally, use
  `x = `, or use the formula method) rather than a bare
  `argument "x" is missing` dispatch failure.
* `decompose_cindex()` returns a `cindex_decomp` object rather than a
  plain list. Existing element access (`$CI_ee`) becomes `$C_ee`.
* `simulate_censoring()` is renamed `censoring_curve()` and now returns
  a `cindex_curve` object, not a plain data frame. Results are in
  `$data` (or via `as.data.frame()`), with columns renamed:
  `ci_ee`/`ci_ec`/`global_c` become `C_ee`/`C_ec`/`C_global`,
  `n_ee`/`n_ec` become `N_ee`/`N_ec`, and `gap` and `low_precision` are
  new.
* `plot_simulation()` and `plot_decomposition()` are replaced by
  `autoplot()` methods.

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
  longer discarded because of their magnitude; `censoring_curve()`'s
  `min_pairs` argument flags rows as `low_precision` instead of deleting
  them (`decompose_cindex()` has no equivalent threshold).
* Plots no longer request the "Arial" font family, which is absent on many
  Linux systems, and no longer impose fixed axis limits that hid
  worse-than-chance models.
* `bootstrap_decomp()` discarded an entire replicate (`C_ee`, `C_ec`,
  `C_global` and the gap all set to `NA`) whenever either side had zero
  weight, instead of NA'ing only the empty side as the point estimate
  already does. Under small samples with heavy censoring this silently
  threw away usable single-sided draws and could shrink the effective
  bootstrap size.
* The dumbbell `autoplot()` suppressed error bars for every model in a
  `compare_decompositions()` plot if even one model had an unstable
  bootstrap confidence interval. Bars are now shown per model.
* `autoplot()` on a `cindex_curve` could draw its shaded band inverted at
  a threshold where the masking gap was negative.
* `autoplot()` on a single `cindex_decomp` no longer labels its row with
  the literal placeholder `"model"`.
* `validate_survival_inputs()` now rejects non-numeric `status` (e.g.
  `status = c("0", "1")`), which previously passed silently.
* `censoring_curve()`'s `n_thresholds` argument is now validated; an
  invalid value errors immediately instead of surfacing later as a
  confusing message from the internal quantile calculation.
* `n_boot = Inf` now gives the package's own validation error instead of a
  confusing internal one.
* A malformed custom weight function (`weights_custom()`) returning the
  wrong length now gets its own error message, distinct from "negative or
  non-finite".

## New features

* Pluggable pair weightings: `weights_harrell()`, `weights_uno()`,
  `weights_truncated()`, `weights_custom()`.
* Bootstrap inference over subjects, with `confint()` covering `C_ee`,
  `C_ec`, `C_global` and the masking gap.
* `compare_decompositions()` for multi-model comparison.
* `theme_cindex()`, light by default with an opt-in dark variant.
