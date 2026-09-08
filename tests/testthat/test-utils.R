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

test_that("validate_weights_and_orientation rejects a non-cindex_weights object", {
  expect_error(validate_weights_and_orientation(list(), TRUE), "cindex_weights")
})

test_that("validate_weights_and_orientation rejects a non-scalar-logical higher_is_riskier", {
  expect_error(validate_weights_and_orientation(weights_harrell(), 0), "TRUE or FALSE")
  expect_error(validate_weights_and_orientation(weights_harrell(), c(TRUE, TRUE)),
               "TRUE or FALSE")
})

test_that("decompose_cindex and censoring_curve give the IDENTICAL weights error", {
  # Regression: the two entry points hand-duplicated this validation and
  # their messages had already drifted (one named weights_harrell() as an
  # example, the other didn't). Pinning byte-identical text so a future
  # edit to one call site can't silently leave the other behind.
  d <- make_test_data(50, seed = 1)
  e1 <- tryCatch(
    decompose_cindex(d$time, d$status, d$risk, weights = list()),
    error = function(e) conditionMessage(e)
  )
  e2 <- tryCatch(
    censoring_curve(d$time, d$status, d$risk, weights = list()),
    error = function(e) conditionMessage(e)
  )
  expect_equal(e1, e2)
})
