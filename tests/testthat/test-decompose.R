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
    survival::Surv(d$time, d$status) ~ d$risk,
    reverse = TRUE
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
  df <- data.frame(
    time = d$time, status = d$status,
    risk = d$risk, other = rnorm(50)
  )
  expect_error(
    decompose_cindex(survival::Surv(time, status) ~ risk + other, data = df),
    "exactly one"
  )
})

test_that("higher_is_riskier = FALSE flips ALL components together", {
  # Regression: the old code flipped only the global value, breaking the
  # identity and drawing a global line above both of its own components.
  d <- make_test_data(200, seed = 6)
  a <- decompose_cindex(d$time, d$status, d$risk, higher_is_riskier = TRUE)
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
  r <- decompose_cindex(d$time, d$status, -d$risk) # deliberately backwards
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

test_that("zero comparable pairs gives NA, never NaN, and doesn't error", {
  # Two events tied at the same time, no censored subjects: the tie rule
  # makes two-events-at-once NOT comparable, and there's nothing else to
  # compare against, so there are zero pairs of either kind. Pooling
  # 0/0 must read as "undefined" (NA), not as NaN from a bare division.
  r <- decompose_cindex(c(5, 5), c(1, 1), c(0.1, 0.9))
  expect_equal(r$N_ee, 0)
  expect_equal(r$N_ec, 0)
  expect_true(is.na(r$C_ee) && !is.nan(r$C_ee))
  expect_true(is.na(r$C_ec) && !is.nan(r$C_ec))
  expect_true(is.na(r$C_global) && !is.nan(r$C_global))
  expect_true(is.na(r$gap) && !is.nan(r$gap))
  expect_output(print(r), "NA")
})

test_that("the formula method rejects NA the same way the vector method does", {
  # Regression: stats::model.frame()'s default na.action is na.omit, so the
  # formula method silently dropped NA rows (changing the analysis
  # population with no message) while the vector method, via
  # validate_survival_inputs(), errored on the same data. Both entry
  # points must now agree.
  set.seed(30)
  n <- 100
  d <- make_test_data(n, seed = 30)
  risk_na <- d$risk
  risk_na[c(3, 47, 90)] <- NA
  df <- data.frame(time = d$time, status = d$status, risk = risk_na)

  expect_error(
    decompose_cindex(d$time, d$status, risk_na),
    "NA"
  )
  expect_error(
    decompose_cindex(survival::Surv(time, status) ~ risk, data = df),
    "missing"
  )
})

test_that("time = gives a self-explanatory error, not a missing-argument one", {
  # Regression: decompose_cindex(time = t, status = s, risk = r) failed
  # with `argument "x" is missing, with no default`, because S3 dispatch
  # forces the generic's first parameter to be `x`. censoring_curve()'s
  # equivalent call style works fine, so this asymmetry needs an
  # explanatory error rather than a bare dispatch failure.
  d <- make_test_data(50, seed = 31)
  expect_error(
    decompose_cindex(time = d$time, status = d$status, risk = d$risk),
    "not `time`"
  )
  # x = still works (it's the correct, S3-required spelling).
  expect_s3_class(
    decompose_cindex(x = d$time, status = d$status, risk = d$risk),
    "cindex_decomp"
  )
})
