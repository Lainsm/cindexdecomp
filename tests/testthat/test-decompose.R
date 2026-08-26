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
