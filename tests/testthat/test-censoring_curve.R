test_that("censoring_curve returns a cindex_curve object", {
  d <- make_test_data(300, seed = 1)
  cv <- censoring_curve(d$time, d$status, d$risk, n_thresholds = 6)
  expect_s3_class(cv, "cindex_curve")
  expect_s3_class(cv$data, "data.frame")
  expect_true(all(c(
    "threshold", "censoring", "C_ee", "C_ec", "C_global",
    "N_ee", "N_ec", "low_precision", "W_ee",
    "W_ec"
  ) %in% names(cv$data)))
})

test_that("higher_is_riskier is validated, not silently coerced", {
  # Regression: `if (!higher_is_riskier) risk <- -risk` treated any
  # non-logical, non-zero value (e.g. 0) as if it were FALSE via `!0`,
  # silently negating risk instead of erroring. decompose_cindex() already
  # rejects this; censoring_curve() must match.
  d <- make_test_data(200, seed = 20)
  expect_error(
    censoring_curve(d$time, d$status, d$risk,
      higher_is_riskier = 0,
      n_thresholds = 4
    ),
    "TRUE or FALSE"
  )
  expect_error(
    censoring_curve(d$time, d$status, d$risk,
      higher_is_riskier = c(TRUE, FALSE), n_thresholds = 4
    ),
    "TRUE or FALSE"
  )
})

test_that("n_thresholds is validated", {
  # Regression: unlike decompose_cindex()'s explicit higher_is_riskier
  # check, n_thresholds had no validation at all -- a bad value only
  # surfaced later as a cryptic error from seq()'s `length.out`.
  d <- make_test_data(200, seed = 20)
  expect_error(
    censoring_curve(d$time, d$status, d$risk, n_thresholds = 0),
    "n_thresholds"
  )
  expect_error(
    censoring_curve(d$time, d$status, d$risk, n_thresholds = -3),
    "n_thresholds"
  )
  expect_error(
    censoring_curve(d$time, d$status, d$risk, n_thresholds = 2.5),
    "n_thresholds"
  )
})

test_that("n_thresholds is NOT validated when probs is supplied directly", {
  # n_thresholds is documented as "ignored if probs is supplied" -- an
  # invalid n_thresholds must not block a call that never uses it.
  d <- make_test_data(200, seed = 20)
  cv <- censoring_curve(d$time, d$status, d$risk,
    n_thresholds = -1,
    probs = c(0.1, 0.3, 0.5)
  )
  expect_s3_class(cv, "cindex_curve")
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
  # Regression: this test used to assert only `any(ok)` -- true even if the
  # identity never actually holds -- so it verified nothing about the
  # decomposition itself. Now it checks
  # C_global == (W_ee * C_ee + W_ec * C_ec) / (W_ee + W_ec) row by row,
  # using the curve's own stored W_ee/W_ec (weighted totals), not N_ee/N_ec
  # (pair counts). Under weights_uno() the two diverge, so reconstructing
  # from counts would fail here even though the true identity holds.
  d <- make_test_data(300, ties = TRUE, seed = 7)
  cv <- censoring_curve(d$time, d$status, d$risk,
    weights = weights_uno(),
    n_thresholds = 6
  )
  ok <- !is.na(cv$data$C_ee) & !is.na(cv$data$C_ec)
  expect_true(any(ok))
  for (i in which(ok)) {
    total_w <- cv$data$W_ee[i] + cv$data$W_ec[i]
    rhs <- (cv$data$W_ee[i] * cv$data$C_ee[i] +
      cv$data$W_ec[i] * cv$data$C_ec[i]) / total_w
    expect_equal(cv$data$C_global[i], rhs, tolerance = 1e-12)
  }
})

test_that("W_ee/W_ec equal N_ee/N_ec under a unit weighting but not under IPCW", {
  d <- make_test_data(300, ties = TRUE, seed = 22)
  harrell <- censoring_curve(d$time, d$status, d$risk, n_thresholds = 5)
  expect_equal(harrell$data$W_ee, as.numeric(harrell$data$N_ee),
    tolerance = 1e-12
  )
  expect_equal(harrell$data$W_ec, as.numeric(harrell$data$N_ec),
    tolerance = 1e-12
  )

  uno <- censoring_curve(d$time, d$status, d$risk,
    weights = weights_uno(),
    n_thresholds = 5
  )
  expect_false(isTRUE(all.equal(uno$data$W_ee, as.numeric(uno$data$N_ee))))
})

test_that("min_pairs flags rather than filters", {
  d <- make_test_data(300, seed = 8)
  a <- censoring_curve(
    d$time, d$status, d$risk,
    n_thresholds = 6, min_pairs = 1
  )
  b <- censoring_curve(d$time, d$status, d$risk,
    n_thresholds = 6,
    min_pairs = 10^9
  )
  expect_equal(nrow(a$data), nrow(b$data))
  expect_true(all(b$data$low_precision))
  expect_false(any(a$data$low_precision))
  expect_false(any(is.na(b$data$C_ee)))
})

test_that("censoring increases with the threshold decreasing", {
  d <- make_test_data(300, seed = 9)
  cv <- censoring_curve(d$time, d$status, d$risk, n_thresholds = 8)
  expect_equal(
    order(cv$data$censoring, decreasing = TRUE),
    order(cv$data$threshold)
  )
})

test_that("the curve accepts alternative weightings", {
  d <- make_test_data(250, seed = 10)
  cv <- censoring_curve(d$time, d$status, d$risk,
    weights = weights_uno(), n_thresholds = 4
  )
  expect_equal(cv$weighting, "Uno (IPCW)")
})

test_that("print is informative", {
  d <- make_test_data(200, seed = 11)
  cv <- censoring_curve(d$time, d$status, d$risk, n_thresholds = 4)
  out <- capture.output(print(cv))
  expect_true(any(grepl("Censoring", out, ignore.case = TRUE)))
})
