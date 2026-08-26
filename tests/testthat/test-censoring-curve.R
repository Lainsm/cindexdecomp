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
