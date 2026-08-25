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
