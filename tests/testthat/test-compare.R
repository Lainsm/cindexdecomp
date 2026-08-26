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
