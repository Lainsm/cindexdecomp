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
  expect_equal(w$fn(100, G), 1 / 0.5^2) # t == tau is INSIDE the window
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
